# =============================================================================
# CONTROLLED INTERRUPTED TIME SERIES (CITS) ANALYSIS
# IRA Insulin Copay Cap Impact
# =============================================================================
# CITS is ideal when:
# 1. Pre-existing trends are present (✓ you have this)
# 2. You have a clear intervention point (✓ January 2023)
# 3. You have a control group (✓ Commercial insurance)
#
# Key advantage: Explicitly models trends, asks if trend CHANGED at IRA
# =============================================================================

library(tidyverse)
library(fixest)
library(modelsummary)
library(ggplot2)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. LOAD AND PREPARE DATA
# -----------------------------------------------------------------------------

cat("=== LOADING PANEL DATA FOR CITS ===\n\n")

monthly_panel <- read_csv("analytical_monthly_panel.csv") %>%
  mutate(
    year_month = as.Date(year_month),
    treatment = as.numeric(treatment),
    post = as.numeric(year_month >= as.Date("2023-01-01")),
    
    # CRITICAL CITS VARIABLES:
    # 1. Time (continuous, runs through entire period)
    time = as.numeric(year_month - min(year_month)) / 30,  # Months since start
    
    # 2. Time since intervention (0 before, increases after)
    time_since_ira = pmax(0, as.numeric(year_month - as.Date("2023-01-01")) / 30),
    
    # 3. Interactions
    treat_time = treatment * time,              # MA-specific pre-trend
    treat_post = treatment * post,              # MA level shift at IRA
    treat_time_since = treatment * time_since_ira,  # MA post-IRA trend change
    
    # Calendar month for seasonality
    calendar_month = factor(month(year_month), levels = 1:12, labels = month.abb)
  )

cat("Prepared", nrow(monthly_panel), "patient-months\n")
cat("Time range:", min(monthly_panel$year_month), "to", max(monthly_panel$year_month), "\n\n")

# -----------------------------------------------------------------------------
# 2. CITS MODEL SPECIFICATION
# -----------------------------------------------------------------------------

cat("=== RUNNING CITS MODELS ===\n\n")

# The CITS model equation:
# Y_it = β0 + β1*time + β2*treatment + β3*post + β4*time_since_ira
#      + β5*treat_time + β6*treat_post + β7*treat_time_since
#      + controls + month FE + patient FE
#
# Key parameters:
# β1 = Commercial time trend (slope before IRA)
# β2 = Baseline difference between MA and Commercial
# β3 = Level shift at IRA for Commercial
# β4 = Commercial trend change after IRA
# β5 = DIFFERENCE in pre-IRA trends (MA vs Commercial)
# β6 = DIFFERENCE in level shifts (DiD estimator)
# β7 = DIFFERENCE in trend changes (DiDiD estimator)

run_cits_model <- function(data, outcome, controls) {
  
  formula_str <- paste0(
    outcome, " ~ ",
    "time + ",                      # Overall time trend
    "treatment + ",                 # Treatment group indicator
    "post + ",                      # Post-IRA indicator
    "time_since_ira + ",           # Post-IRA time trend
    "treat_time + ",               # MA-specific pre-trend
    "treat_post + ",               # MA-specific level shift (DiD)
    "treat_time_since + ",         # MA-specific trend change (DiDiD)
    "calendar_month + ",           # Seasonality
    paste(controls, collapse = " + "),
    " | pat_id"                    # Patient FE
  )
  
  model <- feols(as.formula(formula_str), data = data, cluster = ~pat_id)
  
  # Extract key coefficients
  coef_summary <- tidy(model) %>%
    filter(term %in% c("time", "post", "time_since_ira",
                       "treat_time", "treat_post", "treat_time_since")) %>%
    mutate(
      outcome = outcome,
      ci_lower = estimate - 1.96 * std.error,
      ci_upper = estimate + 1.96 * std.error,
      sig = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  return(list(model = model, coefs = coef_summary))
}

# Define outcomes and controls
primary_outcomes <- c("insulin_oop_per_supply", "insulin_copay_per_supply")
adherence_outcomes <- c("insulin_standardized_supplies", "insulin_gap")
falsification_outcomes <- c("lipid_copay", "metformin_copay")

all_outcomes <- c(primary_outcomes, adherence_outcomes, falsification_outcomes)
controls <- c("n_charlson", "dcsi_total", "n_unique_drugs_excl_insulin")

# Run CITS for all outcomes
cits_results <- map(all_outcomes, function(outcome) {
  cat("Running CITS for:", outcome, "\n")
  run_cits_model(monthly_panel, outcome, controls)
})
names(cits_results) <- all_outcomes

# Combine coefficients
cits_coefs <- map_df(cits_results, ~.x$coefs)
write_csv(cits_coefs, "cits_all_coefficients.csv")

# -----------------------------------------------------------------------------
# 3. INTERPRETATION FUNCTION
# -----------------------------------------------------------------------------

interpret_cits <- function(outcome_name) {
  
  coefs <- cits_coefs %>% filter(outcome == outcome_name)
  
  # Extract key parameters
  comm_pre_trend <- coefs %>% filter(term == "time") %>% pull(estimate)
  comm_level_shift <- coefs %>% filter(term == "post") %>% pull(estimate)
  comm_trend_change <- coefs %>% filter(term == "time_since_ira") %>% pull(estimate)
  
  ma_diff_pre_trend <- coefs %>% filter(term == "treat_time") %>% pull(estimate)
  ma_diff_level_shift <- coefs %>% filter(term == "treat_post") %>% pull(estimate)
  ma_diff_trend_change <- coefs %>% filter(term == "treat_time_since") %>% pull(estimate)
  
  # Calculate MA totals
  ma_pre_trend <- comm_pre_trend + ma_diff_pre_trend
  ma_level_shift <- comm_level_shift + ma_diff_level_shift
  ma_trend_change <- comm_trend_change + ma_diff_trend_change
  ma_post_trend <- ma_pre_trend + ma_trend_change
  
  cat("\n========================================\n")
  cat("CITS INTERPRETATION:", outcome_name, "\n")
  cat("========================================\n\n")
  
  cat("COMMERCIAL (Control Group):\n")
  cat("---------------------------\n")
  cat("  Pre-IRA trend:     ", sprintf("%+.3f", comm_pre_trend), "per month\n")
  cat("  Level shift at IRA:", sprintf("%+.3f", comm_level_shift), "\n")
  cat("  Trend change:      ", sprintf("%+.3f", comm_trend_change), "per month\n")
  cat("  Post-IRA trend:    ", sprintf("%+.3f", comm_pre_trend + comm_trend_change), "per month\n\n")
  
  cat("MEDICARE ADVANTAGE (Treatment Group):\n")
  cat("--------------------------------------\n")
  cat("  Pre-IRA trend:     ", sprintf("%+.3f", ma_pre_trend), "per month\n")
  cat("  Level shift at IRA:", sprintf("%+.3f", ma_level_shift), "\n")
  cat("  Trend change:      ", sprintf("%+.3f", ma_trend_change), "per month\n")
  cat("  Post-IRA trend:    ", sprintf("%+.3f", ma_post_trend), "per month\n\n")
  
  cat("DIFFERENCE-IN-DIFFERENCES (DiD):\n")
  cat("---------------------------------\n")
  cat("  Pre-trend difference: ", sprintf("%+.3f", ma_diff_pre_trend), 
      coefs %>% filter(term == "treat_time") %>% pull(sig), "\n")
  cat("  Level shift (DiD):    ", sprintf("%+.3f", ma_diff_level_shift),
      coefs %>% filter(term == "treat_post") %>% pull(sig), "\n")
  cat("  Trend change (DiDiD): ", sprintf("%+.3f", ma_diff_trend_change),
      coefs %>% filter(term == "treat_time_since") %>% pull(sig), "\n\n")
  
  cat("INTERPRETATION:\n")
  cat("---------------\n")
  
  # Pre-trend
  if(abs(ma_diff_pre_trend) > 0.1 && 
     coefs %>% filter(term == "treat_time") %>% pull(p.value) < 0.05) {
    cat("⚠ MA had DIFFERENT pre-trend than Commercial (parallel trends violated)\n")
  } else {
    cat("✓ Similar pre-trends between groups\n")
  }
  
  # Level shift
  if(coefs %>% filter(term == "treat_post") %>% pull(p.value) < 0.05) {
    if(ma_diff_level_shift < 0) {
      cat("✓ IRA caused IMMEDIATE REDUCTION in MA (relative to Commercial)\n")
    } else {
      cat("⚠ IRA caused IMMEDIATE INCREASE in MA (paradoxical!)\n")
    }
  } else {
    cat("○ No significant immediate effect at IRA\n")
  }
  
  # Trend change
  if(coefs %>% filter(term == "treat_time_since") %>% pull(p.value) < 0.05) {
    if(ma_diff_trend_change < 0) {
      cat("✓ MA trend became MORE NEGATIVE after IRA (accelerated decline)\n")
    } else {
      cat("○ MA trend became LESS NEGATIVE after IRA (slowed decline)\n")
    }
  } else {
    cat("○ No significant change in trend after IRA\n")
  }
  
  cat("\n")
  
  # Return for plotting
  return(list(
    comm_pre = comm_pre_trend,
    comm_post = comm_pre_trend + comm_trend_change,
    comm_level = comm_level_shift,
    ma_pre = ma_pre_trend,
    ma_post = ma_post_trend,
    ma_level = ma_level_shift
  ))
}

# Interpret primary outcomes
params_oop <- interpret_cits("insulin_oop_per_supply")
params_copay <- interpret_cits("insulin_copay_per_supply")

# -----------------------------------------------------------------------------
# 4. VISUALIZATION: OBSERVED DATA + FITTED LINES
# -----------------------------------------------------------------------------

cat("\n=== CREATING CITS VISUALIZATIONS ===\n\n")

create_cits_plot <- function(outcome_name, ylab) {
  
  # Calculate observed means
  plot_data <- monthly_panel %>%
    group_by(year_month, treatment) %>%
    summarise(
      observed = mean(!!sym(outcome_name), na.rm = TRUE),
      se = sd(!!sym(outcome_name), na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      group_label = ifelse(treatment == 1, 
                           "Medicare Advantage", 
                           "Commercial"),
      time = as.numeric(year_month - min(year_month)) / 30
    )
  
  # Get CITS coefficients for this outcome
  coefs <- cits_coefs %>% filter(outcome == outcome_name)
  
  # Extract parameters
  b_time <- coefs %>% filter(term == "time") %>% pull(estimate)
  b_post <- coefs %>% filter(term == "post") %>% pull(estimate)
  b_time_since <- coefs %>% filter(term == "time_since_ira") %>% pull(estimate)
  b_treat_time <- coefs %>% filter(term == "treat_time") %>% pull(estimate)
  b_treat_post <- coefs %>% filter(term == "treat_post") %>% pull(estimate)
  b_treat_time_since <- coefs %>% filter(term == "treat_time_since") %>% pull(estimate)
  
  # Time at IRA
  ira_time <- as.numeric(as.Date("2023-01-01") - min(monthly_panel$year_month)) / 30
  
  # Create fitted lines
  fitted_lines <- expand_grid(
    year_month = seq(min(monthly_panel$year_month),
                     max(monthly_panel$year_month),
                     by = "month"),
    treatment = c(0, 1)
  ) %>%
    mutate(
      time = as.numeric(year_month - min(year_month)) / 30,
      post = as.numeric(year_month >= as.Date("2023-01-01")),
      time_since_ira = pmax(0, time - ira_time),
      
      # Fitted values from CITS model
      fitted = b_time * time +
        b_post * post +
        b_time_since * time_since_ira +
        treatment * (b_treat_time * time +
                       b_treat_post * post +
                       b_treat_time_since * time_since_ira),
      
      group_label = ifelse(treatment == 1, 
                           "Medicare Advantage", 
                           "Commercial"),
      period = ifelse(year_month < as.Date("2023-01-01"), "Pre-IRA", "Post-IRA")
    )
  
  # Create counterfactual (what if IRA never happened?)
  counterfactual <- fitted_lines %>%
    filter(treatment == 1, year_month >= as.Date("2023-01-01")) %>%
    mutate(
      # Remove IRA effects: no level shift, no trend change
      counterfactual = b_time * time +
        (b_time + b_treat_time) * time,  # Just continue pre-trend
      type = "Counterfactual\n(No IRA)"
    )
  
  # Base observed mean at Dec 2022 for counterfactual
  base_ma <- plot_data %>%
    filter(treatment == 1, year_month == as.Date("2022-12-01")) %>%
    pull(observed)
  
  base_fitted <- fitted_lines %>%
    filter(treatment == 1, year_month == as.Date("2022-12-01")) %>%
    pull(fitted)
  
  # Adjust counterfactual to observed level
  counterfactual <- counterfactual %>%
    mutate(counterfactual = counterfactual - base_fitted + base_ma)
  
  # Create plot
  p <- ggplot() +
    
    # IRA line
    geom_vline(xintercept = as.Date("2023-01-01"), 
               linetype = "dashed", color = "red", linewidth = 1) +
    
    # Observed data points
    geom_point(data = plot_data,
               aes(x = year_month, y = observed, color = group_label),
               size = 2.5, alpha = 0.7) +
    
    geom_errorbar(data = plot_data,
                  aes(x = year_month, y = observed, 
                      ymin = observed - 1.96*se,
                      ymax = observed + 1.96*se,
                      color = group_label),
                  width = 10, alpha = 0.4) +
    
    # Fitted lines (segmented by period)
    geom_line(data = fitted_lines %>% filter(period == "Pre-IRA"),
              aes(x = year_month, y = fitted, color = group_label,
                  group = interaction(group_label, period)),
              linewidth = 1.5, linetype = "solid") +
    
    geom_line(data = fitted_lines %>% filter(period == "Post-IRA"),
              aes(x = year_month, y = fitted, color = group_label,
                  group = interaction(group_label, period)),
              linewidth = 1.5, linetype = "solid") +
    
    # Counterfactual line (dashed)
    geom_line(data = counterfactual,
              aes(x = year_month, y = counterfactual),
              color = "#2c7bb6", linewidth = 1.2, linetype = "dashed", alpha = 0.8) +
    
    # Shaded region showing effect
    geom_ribbon(data = counterfactual %>%
                  left_join(fitted_lines %>% 
                              filter(treatment == 1, year_month >= as.Date("2023-01-01")) %>%
                              select(year_month, fitted),
                            by = "year_month"),
                aes(x = year_month,
                    ymin = pmin(counterfactual, fitted),
                    ymax = pmax(counterfactual, fitted)),
                fill = "purple", alpha = 0.2) +
    
    # Colors
    scale_color_manual(
      values = c("Medicare Advantage" = "#2c7bb6",
                 "Commercial" = "#d7191c"),
      name = "Group"
    ) +
    
    # Labels
    labs(
      title = paste("CITS Analysis:", ylab),
      subtitle = paste0(
        "Solid lines = CITS fitted trends | Dashed line = MA counterfactual (no IRA)\n",
        "Purple area = IRA effect (observed vs counterfactual)"
      ),
      x = "Month",
      y = ylab,
      caption = paste0(
        "Points = observed monthly means (95% CI). Red line = IRA (Jan 2023).\n",
        "CITS model controls for pre-existing trends and seasonality.\n",
        "Counterfactual shows what would have happened if MA continued pre-IRA trend."
      )
    ) +
    
    # Theme
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 11, color = "gray30", lineheight = 1.2),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      plot.caption = element_text(size = 9, hjust = 0, color = "gray50", lineheight = 1.2)
    )
  
  return(p)
}

# Create plots for primary outcomes
plot_oop <- create_cits_plot("insulin_oop_per_supply", "OOP per 30-Day Supply ($)")
plot_copay <- create_cits_plot("insulin_copay_per_supply", "Copay per 30-Day Supply ($)")

# Save
ggsave("cits_insulin_oop.png", plot_oop, width = 14, height = 9, dpi = 300)
ggsave("cits_insulin_oop.pdf", plot_oop, width = 14, height = 9)

ggsave("cits_insulin_copay.png", plot_copay, width = 14, height = 9, dpi = 300)
ggsave("cits_insulin_copay.pdf", plot_copay, width = 14, height = 9)

# Combined
combined_cits <- plot_oop / plot_copay +
  plot_annotation(
    title = "Controlled Interrupted Time Series: IRA Insulin Copay Cap",
    subtitle = "Addresses pre-existing trends by modeling time explicitly",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

ggsave("cits_combined.png", combined_cits, width = 14, height = 16, dpi = 300)
ggsave("cits_combined.pdf", combined_cits, width = 14, height = 16)

# -----------------------------------------------------------------------------
# 5. REGRESSION TABLES
# -----------------------------------------------------------------------------

cat("\n=== CREATING REGRESSION TABLES ===\n\n")

# Primary outcomes table
primary_models <- list(
  "OOP per Supply" = cits_results$insulin_oop_per_supply$model,
  "Copay per Supply" = cits_results$insulin_copay_per_supply$model
)

modelsummary(
  primary_models,
  stars = c('*' = 0.05, '**' = 0.01, '***' = 0.001),
  coef_map = c(
    "time" = "Time (Commercial pre-trend)",
    "post" = "Post-IRA (Commercial level shift)",
    "time_since_ira" = "Time×Post (Commercial trend change)",
    "treat_time" = "MA × Time (pre-trend difference)",
    "treat_post" = "MA × Post (DiD: level shift)",
    "treat_time_since" = "MA × Time×Post (DiDiD: trend change)"
  ),
  gof_map = c("nobs", "r.squared", "adj.r.squared"),
  notes = c(
    "CITS model with patient fixed effects and month dummies.",
    "Standard errors clustered at patient level.",
    "treat_post = DiD estimator (immediate effect)",
    "treat_time_since = DiDiD estimator (change in trend)"
  ),
  output = "cits_primary_outcomes.docx"
)

# All outcomes table
all_models <- list(
  "OOP/Supply" = cits_results$insulin_oop_per_supply$model,
  "Copay/Supply" = cits_results$insulin_copay_per_supply$model,
  "Supplies" = cits_results$insulin_standardized_supplies$model,
  "Gap" = cits_results$insulin_gap$model
)

modelsummary(
  all_models,
  stars = c('*' = 0.05, '**' = 0.01, '***' = 0.001),
  coef_map = c(
    "treat_time" = "MA Pre-Trend Difference",
    "treat_post" = "MA Level Shift (DiD)",
    "treat_time_since" = "MA Trend Change (DiDiD)"
  ),
  gof_map = c("nobs", "r.squared"),
  notes = c(
    "CITS models with patient FE and month FE.",
    "Only key treatment interactions shown.",
    "See full specification for Commercial parameters."
  ),
  output = "cits_all_outcomes.docx"
)

# Falsification tests
falsification_models <- list(
  "Lipid Copay" = cits_results$lipid_copay$model,
  "Metformin Copay" = cits_results$metformin_copay$model
)

modelsummary(
  falsification_models,
  stars = c('*' = 0.05, '**' = 0.01, '***' = 0.001),
  coef_map = c(
    "treat_time" = "MA Pre-Trend Difference",
    "treat_post" = "MA Level Shift",
    "treat_time_since" = "MA Trend Change"
  ),
  gof_map = c("nobs", "r.squared"),
  notes = c(
    "Falsification tests: Non-insulin medications.",
    "Should show no effects if design is valid."
  ),
  output = "cits_falsification.docx"
)

# -----------------------------------------------------------------------------
# 6. SUMMARY AND INTERPRETATION GUIDE
# -----------------------------------------------------------------------------

cat("\n\n=============================================================\n")
cat("CITS ANALYSIS COMPLETE\n")
cat("=============================================================\n\n")

cat("KEY PARAMETERS TO REPORT:\n")
cat("--------------------------\n")
cat("1. treat_time (β5): MA pre-trend different from Commercial?\n")
cat("   → If significant: Parallel trends violated (expected!)\n\n")

cat("2. treat_post (β6): Immediate effect at IRA (DiD)\n")
cat("   → Main policy effect controlling for trends\n\n")

cat("3. treat_time_since (β7): Change in MA trend after IRA (DiDiD)\n")
cat("   → Did decline accelerate/slow after IRA?\n\n")

cat("INTERPRETATION FRAMEWORK:\n")
cat("-------------------------\n")
cat("β5 < 0, significant: MA was declining faster pre-IRA ✓\n")
cat("β6 = 0, NS: No immediate level shift\n")
cat("β7 = 0, NS: No change in trend\n")
cat("→ Conclusion: IRA had no additional effect beyond pre-existing trends\n\n")

cat("ADVANTAGES OF CITS OVER EVENT STUDY:\n")
cat("-------------------------------------\n")
cat("✓ Explicitly models pre-trends (not just controls for them)\n")
cat("✓ Provides clear counterfactual (dashed line in plots)\n")
cat("✓ Tests TWO mechanisms: level shift AND trend change\n")
cat("✓ Visual interpretation is intuitive\n")
cat("✓ Standard in policy evaluation literature\n\n")

cat("FILES CREATED:\n")
cat("--------------\n")
cat("Data:\n")
cat("  - cits_all_coefficients.csv\n\n")
cat("Plots:\n")
cat("  - cits_insulin_oop.png/pdf\n")
cat("  - cits_insulin_copay.png/pdf\n")
cat("  - cits_combined.png/pdf\n\n")
cat("Tables:\n")
cat("  - cits_primary_outcomes.docx\n")
cat("  - cits_all_outcomes.docx\n")
cat("  - cits_falsification.docx\n\n")

cat("FOR YOUR DISSERTATION:\n")
cat("----------------------\n")
cat("1. Use CITS as PRIMARY analysis (not event study)\n")
cat("2. Show CITS plot in main results (clearest visual)\n")
cat("3. Report β6 (treat_post) as main effect estimate\n")
cat("4. Report β7 (treat_time_since) as trend change\n")
cat("5. Acknowledge β5 (treat_time) shows pre-trend violation\n")
cat("6. Frame: 'CITS accounts for anticipatory market changes'\n\n")

cat("NEXT STEPS:\n")
cat("-----------\n")
cat("1. Run this code\n")
cat("2. Examine CITS plots (much clearer than event study!)\n")
cat("3. Check if β6 and β7 are significant\n")
cat("4. Compare Commercial vs MA trends visually\n")
cat("5. Use counterfactual to show 'what if no IRA'\n\n")