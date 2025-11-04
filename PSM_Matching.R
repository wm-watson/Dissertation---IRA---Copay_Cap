# =============================================================================
# TREND-ADJUSTED DID MODELS - IRA INSULIN COPAY CAP
# =============================================================================
# Addresses parallel trends violations by allowing group-specific time trends
# 
# Key features:
# 1. Group-specific pre-IRA trends (treatment:months_since_start)
# 2. Group-specific post-IRA trend changes (treatment:months_since_ira)
# 3. Month-by-cohort FE for seasonality (calendar_month*treatment)
# 4. Patient fixed effects
#
# Interpretation:
# - treat_post = Immediate level shift at IRA (net of trends)
# - treatment:months_since_ira = Change in trend slope post-IRA
# =============================================================================

library(tidyverse)
library(fixest)
library(broom)
library(modelsummary)

# -----------------------------------------------------------------------------
# 1. LOAD DATA
# -----------------------------------------------------------------------------

cat("=== LOADING PANEL DATA FOR TREND-ADJUSTED ANALYSIS ===\n\n")

monthly_panel <- read_csv("analytical_monthly_panel.csv")

# Clean and prepare
monthly_panel <- monthly_panel %>%
  mutate(
    # Replace Inf/-Inf with NA
    across(where(is.numeric), ~ifelse(is.infinite(.), NA, .)),
    
    # Replace NA with 0 for cost outcomes
    across(c(insulin_copay, insulin_oop, insulin_copay_per_supply, insulin_oop_per_supply,
             metformin_copay, lipid_copay), ~replace_na(., 0)),
    
    # Ensure proper types
    year_month = as.Date(year_month),
    treatment = as.numeric(treatment),
    post = as.numeric(post_ira),
    treat_post = treatment * post,
    
    # CRITICAL: Create time trend variables
    # Time since panel start (continuous, spans both years)
    months_since_start = as.numeric(year_month - min(year_month)) / 30,
    
    # Time since IRA (0 before IRA, increases after)
    months_since_ira = pmax(0, as.numeric(year_month - as.Date("2023-01-01")) / 30),
    
    # Calendar month factor for seasonality
    calendar_month = factor(month(year_month), levels = 1:12, labels = month.abb)
  )

cat("Panel prepared:\n")
cat("  Patients:", n_distinct(monthly_panel$pat_id), "\n")
cat("  Patient-months:", nrow(monthly_panel), "\n")
cat("  Time range:", min(monthly_panel$year_month), "to", max(monthly_panel$year_month), "\n\n")

# -----------------------------------------------------------------------------
# 2. DEFINE OUTCOMES
# -----------------------------------------------------------------------------

primary_outcomes <- c("insulin_oop_per_supply", "insulin_copay_per_supply")
adherence_outcomes <- c("insulin_standardized_supplies", "insulin_gap")
secondary_outcomes <- c("insulin_oop", "insulin_copay")
falsification_outcomes <- c("lipid_copay", "metformin_copay")

all_outcomes <- c(primary_outcomes, adherence_outcomes, 
                  secondary_outcomes, falsification_outcomes)

binned_controls <- c("charlson_bin", "dcsi_bin", "polypharm_bin")

# -----------------------------------------------------------------------------
# 3. TREND-ADJUSTED DID FUNCTION
# -----------------------------------------------------------------------------

run_trend_adjusted_did <- function(data, outcome, controls) {
  # Model specification:
  # Outcome = β0 + β1*treat_post                          [LEVEL SHIFT at IRA]
  #         + β2*treatment*months_since_start             [MA pre-trend]
  #         + β3*treatment*months_since_ira               [MA post-trend CHANGE]
  #         + β4*calendar_month*treatment                 [Seasonality by group]
  #         + controls + patient FE
  #
  # Interpretation:
  # β1 = Immediate effect of IRA (level shift), controlling for trends
  # β2 = MA was declining by β2 per month before IRA
  # β3 = MA trend changed by β3 per month after IRA
  #      - If β3 < 0: Decline accelerated post-IRA
  #      - If β3 > 0: Decline slowed/reversed post-IRA
  #      - If β3 ≈ 0: Trend unchanged, only level shift
  
  formula_str <- paste0(
    outcome, " ~ treat_post + ",
    "treatment:months_since_start + ",       # Pre-IRA MA trend
    "treatment:months_since_ira + ",         # Post-IRA MA trend change
    "calendar_month*treatment + ",           # Seasonality
    paste(controls, collapse = " + "),
    " | pat_id"                              # Patient FE
  )
  
  model <- feols(as.formula(formula_str), data = data, cluster = ~pat_id)
  
  # Extract key coefficients
  coef_summary <- tidy(model) %>%
    filter(term %in% c("treat_post", 
                       "treatment:months_since_start",
                       "treatment:months_since_ira")) %>%
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

# -----------------------------------------------------------------------------
# 4. RUN TREND-ADJUSTED MODELS FOR ALL OUTCOMES
# -----------------------------------------------------------------------------

cat("=== RUNNING TREND-ADJUSTED DID MODELS ===\n\n")

trend_adjusted_results <- map(all_outcomes, function(outcome) {
  cat("Processing:", outcome, "\n")
  run_trend_adjusted_did(monthly_panel, outcome, binned_controls)
})
names(trend_adjusted_results) <- all_outcomes

# Combine coefficients
trend_adjusted_coefs <- map_df(trend_adjusted_results, ~.x$coefs)

write_csv(trend_adjusted_coefs, "trend_adjusted_did_coefficients.csv")

# -----------------------------------------------------------------------------
# 5. COMPARISON: SIMPLE DID VS TREND-ADJUSTED DID
# -----------------------------------------------------------------------------

cat("\n=== COMPARING SIMPLE DID VS TREND-ADJUSTED DID ===\n\n")

# Run simple DiD for comparison (no trend controls)
run_simple_did_panel <- function(data, outcome, controls) {
  formula_str <- paste0(
    outcome, " ~ treat_post + calendar_month*treatment + ",
    paste(controls, collapse = " + "),
    " | pat_id"
  )
  
  model <- feols(as.formula(formula_str), data = data, cluster = ~pat_id)
  
  did_estimate <- tidy(model) %>%
    filter(term == "treat_post") %>%
    mutate(outcome = outcome)
  
  return(list(model = model, did = did_estimate))
}

simple_results <- map(all_outcomes, ~run_simple_did_panel(monthly_panel, .x, binned_controls))
names(simple_results) <- all_outcomes

simple_coefs <- map_df(simple_results, ~.x$did) %>%
  mutate(
    specification = "Simple DiD",
    ci_lower = estimate - 1.96 * std.error,
    ci_upper = estimate + 1.96 * std.error
  )

# Extract just treat_post from trend-adjusted
trend_adjusted_treat_post <- trend_adjusted_coefs %>%
  filter(term == "treat_post") %>%
  mutate(specification = "Trend-Adjusted DiD")

# Combine for comparison
comparison <- bind_rows(
  simple_coefs %>% select(specification, outcome, estimate, std.error, ci_lower, ci_upper, p.value),
  trend_adjusted_treat_post %>% select(specification, outcome, estimate, std.error, ci_lower, ci_upper, p.value)
) %>%
  mutate(
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

write_csv(comparison, "comparison_simple_vs_trend_adjusted.csv")

# Print comparison for primary outcomes
cat("\nComparison for Primary Outcomes:\n")
print(comparison %>% 
        filter(outcome %in% primary_outcomes) %>%
        select(specification, outcome, estimate, std.error, p.value, sig))

# -----------------------------------------------------------------------------
# 6. INTERPRETATION HELPER FUNCTION
# -----------------------------------------------------------------------------

interpret_trend_adjusted_results <- function(outcome_name) {
  results <- trend_adjusted_coefs %>% filter(outcome == outcome_name)
  
  level_shift <- results %>% filter(term == "treat_post")
  pre_trend <- results %>% filter(term == "treatment:months_since_start")
  post_trend_change <- results %>% filter(term == "treatment:months_since_ira")
  
  cat("\n=== INTERPRETATION:", outcome_name, "===\n")
  cat("\n1. PRE-IRA TREND (MA vs Commercial):\n")
  if(nrow(pre_trend) > 0) {
    cat("   MA was changing by", round(pre_trend$estimate, 3), "per month BEFORE IRA\n")
    cat("   (p =", format.pval(pre_trend$p.value, digits = 3), pre_trend$sig, ")\n")
    if(pre_trend$p.value < 0.05) {
      if(pre_trend$estimate < 0) {
        cat("   → MA costs were DECLINING before IRA (parallel trends violated)\n")
      } else {
        cat("   → MA costs were INCREASING before IRA (parallel trends violated)\n")
      }
    } else {
      cat("   → No significant pre-trend (parallel trends OK)\n")
    }
  }
  
  cat("\n2. IMMEDIATE IRA EFFECT (level shift):\n")
  if(nrow(level_shift) > 0) {
    cat("   IRA caused immediate change of", round(level_shift$estimate, 3), "\n")
    cat("   (p =", format.pval(level_shift$p.value, digits = 3), level_shift$sig, ")\n")
    cat("   95% CI: [", round(level_shift$ci_lower, 3), ",", round(level_shift$ci_upper, 3), "]\n")
  }
  
  cat("\n3. POST-IRA TREND CHANGE:\n")
  if(nrow(post_trend_change) > 0) {
    cat("   MA trend changed by", round(post_trend_change$estimate, 3), "per month after IRA\n")
    cat("   (p =", format.pval(post_trend_change$p.value, digits = 3), post_trend_change$sig, ")\n")
    if(post_trend_change$p.value < 0.05) {
      if(pre_trend$estimate < 0 & post_trend_change$estimate < 0) {
        cat("   → Pre-existing decline ACCELERATED after IRA\n")
      } else if(pre_trend$estimate < 0 & post_trend_change$estimate > 0) {
        cat("   → Pre-existing decline SLOWED/REVERSED after IRA\n")
      } else if(abs(post_trend_change$estimate) < 0.1) {
        cat("   → Trend essentially unchanged, IRA caused level shift only\n")
      }
    } else {
      cat("   → No significant change in trend slope\n")
    }
  }
  
  cat("\n4. TOTAL EFFECT (10 months post-IRA):\n")
  if(nrow(level_shift) > 0 & nrow(post_trend_change) > 0) {
    # Effect = Level shift + (Trend change × 10 months)
    total_effect <- level_shift$estimate + (post_trend_change$estimate * 10)
    cat("   ", round(total_effect, 3), "by October 2023\n")
    cat("   (Level shift + 10 months of trend change)\n")
  }
  cat("\n")
}

# Interpret primary outcomes
interpret_trend_adjusted_results("insulin_oop_per_supply")
interpret_trend_adjusted_results("insulin_copay_per_supply")

# -----------------------------------------------------------------------------
# 7. VISUALIZATION: FITTED VALUES FROM BOTH MODELS (FIXED)
# -----------------------------------------------------------------------------

cat("\n=== CREATING COMPARISON VISUALIZATION ===\n\n")

# Function to create counterfactual predictions (FIXED)
create_trend_comparison_data <- function(outcome_name) {
  
  # Get actual data and calculate group means by month
  actual_means <- monthly_panel %>%
    group_by(year_month, treatment) %>%
    summarise(
      actual_mean = mean(!!sym(outcome_name), na.rm = TRUE),
      se = sd(!!sym(outcome_name), na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      group_label = ifelse(treatment == 1, "Medicare Advantage (65+)", "Commercial (18-64)"),
      months_since_start = as.numeric(year_month - min(year_month)) / 30,
      months_since_ira = pmax(0, as.numeric(year_month - as.Date("2023-01-01")) / 30),
      post = as.numeric(year_month >= as.Date("2023-01-01")),
      treat_post = treatment * post
    )
  
  # Get model coefficients
  simple_model <- simple_results[[outcome_name]]$model
  trend_model <- trend_adjusted_results[[outcome_name]]$model
  
  # Extract key coefficients from simple model
  simple_coefs <- tidy(simple_model) %>%
    filter(term == "treat_post") %>%
    pull(estimate)
  
  # Extract key coefficients from trend-adjusted model
  trend_coefs_df <- tidy(trend_model)
  
  trend_treat_post <- trend_coefs_df %>%
    filter(term == "treat_post") %>%
    pull(estimate)
  
  trend_pre_trend <- trend_coefs_df %>%
    filter(term == "treatment:months_since_start") %>%
    pull(estimate)
  
  if(length(trend_pre_trend) == 0) trend_pre_trend <- 0
  
  trend_post_change <- trend_coefs_df %>%
    filter(term == "treatment:months_since_ira") %>%
    pull(estimate)
  
  if(length(trend_post_change) == 0) trend_post_change <- 0
  
  # Calculate fitted values based on model structure
  plot_data <- actual_means %>%
    mutate(
      # Simple DiD fitted: Just add treat_post effect to treatment group post-period
      simple_fitted = case_when(
        treatment == 1 & post == 1 ~ actual_mean - (simple_coefs * treat_post),
        TRUE ~ actual_mean
      ),
      
      # Trend-adjusted fitted: Remove trends and add back model predictions
      trend_fitted = case_when(
        treatment == 1 ~ actual_mean - 
          (trend_pre_trend * months_since_start) -
          (trend_treat_post * post) -
          (trend_post_change * months_since_ira),
        TRUE ~ actual_mean
      )
    )
  
  return(plot_data)
}

# Create plot for primary outcome
plot_data_oop <- create_trend_comparison_data("insulin_oop_per_supply")

# Create visualization
plot_trend_comparison <- ggplot(plot_data_oop, 
                                aes(x = year_month, color = group_label)) +
  geom_vline(xintercept = as.Date("2023-01-01"), linetype = "dashed", 
             color = "red", linewidth = 1, alpha = 0.7) +
  
  # Actual observed data (points with error bars)
  geom_point(aes(y = actual_mean), size = 2, alpha = 0.6) +
  geom_errorbar(aes(y = actual_mean, ymin = actual_mean - 1.96*se, 
                    ymax = actual_mean + 1.96*se), 
                width = 10, alpha = 0.3) +
  
  # Simple DiD trend line
  geom_line(aes(y = actual_mean, linetype = "Observed Data"), 
            linewidth = 1.2, alpha = 0.8) +
  
  scale_color_manual(
    values = c("Medicare Advantage (65+)" = "#2c7bb6",
               "Commercial (18-64)" = "#d7191c"),
    name = "Group"
  ) +
  
  scale_linetype_manual(
    values = c("Observed Data" = "solid"),
    name = NULL
  ) +
  
  labs(
    title = "Trend-Adjusted DiD: Insulin OOP per 30-Day Supply",
    subtitle = paste0(
      "MA Pre-IRA Trend: -$0.37/month (p=0.026)\n",
      "IRA Level Shift: +$3.17 (p=0.12, NS)\n",
      "Simple DiD Effect: -$1.78 vs Trend-Adjusted: +$3.17"
    ),
    x = "Month",
    y = "OOP per 30-Day Supply ($)",
    caption = "Points = monthly means with 95% CI. Red line = IRA (Jan 2023).\nPre-existing MA decline explains most of simple DiD effect."
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid.minor = element_blank(),
    plot.caption = element_text(size = 9, hjust = 0, color = "gray50")
  )

ggsave("trend_comparison_observed_data.png", plot_trend_comparison,
       width = 12, height = 8, dpi = 300)

cat("Comparison plot saved: trend_comparison_observed_data.png\n")

# -----------------------------------------------------------------------------
# ALTERNATIVE: Create separate counterfactual plot
# -----------------------------------------------------------------------------

# Calculate counterfactuals manually
create_counterfactual_plot <- function(outcome_name) {
  
  # Get treatment group means
  ma_data <- monthly_panel %>%
    filter(treatment == 1) %>%
    group_by(year_month) %>%
    summarise(
      observed = mean(!!sym(outcome_name), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      months_since_start = as.numeric(year_month - min(year_month)) / 30,
      months_since_ira = pmax(0, as.numeric(year_month - as.Date("2023-01-01")) / 30),
      post = as.numeric(year_month >= as.Date("2023-01-01"))
    )
  
  # Get trend-adjusted coefficients
  trend_coefs_df <- tidy(trend_adjusted_results[[outcome_name]]$model)
  
  pre_trend <- trend_coefs_df %>%
    filter(term == "treatment:months_since_start") %>%
    pull(estimate)
  if(length(pre_trend) == 0) pre_trend <- 0
  
  level_shift <- trend_coefs_df %>%
    filter(term == "treat_post") %>%
    pull(estimate)
  if(length(level_shift) == 0) level_shift <- 0
  
  post_trend_change <- trend_coefs_df %>%
    filter(term == "treatment:months_since_ira") %>%
    pull(estimate)
  if(length(post_trend_change) == 0) post_trend_change <- 0
  
  # Create counterfactual: What if IRA never happened?
  ma_data <- ma_data %>%
    mutate(
      # Counterfactual = remove IRA effects, keep only pre-trend
      counterfactual = observed - 
        (level_shift * post) - 
        (post_trend_change * months_since_ira),
      
      # Effect = observed - counterfactual
      effect = observed - counterfactual
    )
  
  # Plot
  ggplot(ma_data, aes(x = year_month)) +
    geom_vline(xintercept = as.Date("2023-01-01"), linetype = "dashed", 
               color = "red", linewidth = 1) +
    
    # Observed
    geom_line(aes(y = observed, color = "Observed MA"), linewidth = 1.2) +
    geom_point(aes(y = observed, color = "Observed MA"), size = 2) +
    
    # Counterfactual
    geom_line(aes(y = counterfactual, color = "Counterfactual (no IRA)"), 
              linewidth = 1.2, linetype = "dashed") +
    geom_point(aes(y = counterfactual, color = "Counterfactual (no IRA)"), 
               size = 2, shape = 1) +
    
    # Shaded area showing effect
    geom_ribbon(aes(ymin = pmin(observed, counterfactual),
                    ymax = pmax(observed, counterfactual)),
                alpha = 0.2, fill = "purple") +
    
    scale_color_manual(
      values = c("Observed MA" = "#2c7bb6",
                 "Counterfactual (no IRA)" = "#d7191c"),
      name = NULL
    ) +
    
    labs(
      title = "Trend-Adjusted Counterfactual Analysis",
      subtitle = paste0(
        outcome_name, 
        "\nPurple area = IRA effect after controlling for pre-existing trend"
      ),
      x = "Month",
      y = "Mean Value",
      caption = "Counterfactual assumes MA would have continued pre-IRA trend.\nShaded area = treatment effect net of trend."
    ) +
    
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

# Create counterfactual plots for primary outcomes
plot_cf_oop <- create_counterfactual_plot("insulin_oop_per_supply")
plot_cf_copay <- create_counterfactual_plot("insulin_copay_per_supply")

ggsave("counterfactual_oop.png", plot_cf_oop, width = 12, height = 8, dpi = 300)
ggsave("counterfactual_copay.png", plot_cf_copay, width = 12, height = 8, dpi = 300)

cat("Counterfactual plots saved\n")

# -----------------------------------------------------------------------------
# 8. REGRESSION TABLES
# -----------------------------------------------------------------------------

cat("\n=== CREATING REGRESSION TABLES ===\n\n")

# Table 1: Primary outcomes - Full specification comparison
primary_models_comparison <- list(
  "OOP (Simple)" = simple_results$insulin_oop_per_supply$model,
  "OOP (Trend-Adj)" = trend_adjusted_results$insulin_oop_per_supply$model,
  "Copay (Simple)" = simple_results$insulin_copay_per_supply$model,
  "Copay (Trend-Adj)" = trend_adjusted_results$insulin_copay_per_supply$model
)

modelsummary(
  primary_models_comparison,
  stars = c('*' = 0.05, '**' = 0.01, '***' = 0.001),
  coef_map = c(
    "treat_post" = "Treatment × Post-IRA",
    "treatment:months_since_start" = "MA × Pre-IRA Trend",
    "treatment:months_since_ira" = "MA × Post-IRA Trend Change"
  ),
  gof_map = c("nobs", "r.squared", "adj.r.squared"),
  notes = c(
    "All models include patient fixed effects and month-by-cohort fixed effects.",
    "Standard errors clustered at patient level.",
    "Trend-adjusted models allow MA to have different time trends than Commercial.",
    "Month dummy coefficients omitted for brevity."
  ),
  output = "trend_adjusted_PRIMARY_comparison_table.docx"
)

# Table 2: All outcomes - Trend-adjusted only
all_trend_models <- list(
  "OOP/Supply" = trend_adjusted_results$insulin_oop_per_supply$model,
  "Copay/Supply" = trend_adjusted_results$insulin_copay_per_supply$model,
  "Supplies" = trend_adjusted_results$insulin_standardized_supplies$model,
  "Gap" = trend_adjusted_results$insulin_gap$model
)

modelsummary(
  all_trend_models,
  stars = c('*' = 0.05, '**' = 0.01, '***' = 0.001),
  coef_map = c(
    "treat_post" = "Treatment × Post-IRA (Level Shift)",
    "treatment:months_since_start" = "MA × Pre-IRA Trend ($/month)",
    "treatment:months_since_ira" = "MA × Post-IRA Trend Change ($/month)"
  ),
  gof_map = c("nobs", "r.squared"),
  notes = c(
    "All models include patient FE and month-by-cohort FE.",
    "Clustered SEs at patient level.",
    "Pre-IRA trend shows MA was already declining before policy.",
    "Post-IRA trend change shows if decline accelerated/slowed after policy."
  ),
  output = "trend_adjusted_ALL_outcomes_table.docx"
)

# Table 3: Falsification tests
falsification_models <- list(
  "Lipid (Simple)" = simple_results$lipid_copay$model,
  "Lipid (Trend-Adj)" = trend_adjusted_results$lipid_copay$model,
  "Metformin (Simple)" = simple_results$metformin_copay$model,
  "Metformin (Trend-Adj)" = trend_adjusted_results$metformin_copay$model
)

modelsummary(
  falsification_models,
  stars = c('*' = 0.05, '**' = 0.01, '***' = 0.001),
  coef_map = c(
    "treat_post" = "Treatment × Post-IRA",
    "treatment:months_since_start" = "MA × Pre-IRA Trend",
    "treatment:months_since_ira" = "MA × Post-IRA Trend Change"
  ),
  gof_map = c("nobs", "r.squared"),
  notes = c(
    "Falsification tests: Non-insulin medications should show no effect.",
    "All models include patient FE and month-by-cohort FE.",
    "Finding effects here would invalidate identification strategy."
  ),
  output = "trend_adjusted_FALSIFICATION_table.docx"
)

# -----------------------------------------------------------------------------
# 9. SUMMARY STATISTICS
# -----------------------------------------------------------------------------

cat("\n\n=== TREND-ADJUSTED DID ANALYSIS SUMMARY ===\n\n")

cat("KEY FINDINGS:\n\n")

# Primary outcome summary
oop_results <- trend_adjusted_coefs %>% filter(outcome == "insulin_oop_per_supply")

cat("INSULIN OOP PER 30-DAY SUPPLY:\n")
cat("  Pre-IRA MA trend:", 
    round(oop_results$estimate[oop_results$term == "treatment:months_since_start"], 3),
    "per month", 
    oop_results$sig[oop_results$term == "treatment:months_since_start"], "\n")
cat("  IRA level shift:", 
    round(oop_results$estimate[oop_results$term == "treat_post"], 3),
    oop_results$sig[oop_results$term == "treat_post"], "\n")
cat("  Post-IRA trend change:", 
    round(oop_results$estimate[oop_results$term == "treatment:months_since_ira"], 3),
    "per month",
    oop_results$sig[oop_results$term == "treatment:months_since_ira"], "\n\n")

# Comparison
comparison_oop <- comparison %>% filter(outcome == "insulin_oop_per_supply")
cat("EFFECT ESTIMATE COMPARISON:\n")
cat("  Simple DiD:        ", round(comparison_oop$estimate[comparison_oop$specification == "Simple DiD"], 3), "\n")
cat("  Trend-Adjusted:    ", round(comparison_oop$estimate[comparison_oop$specification == "Trend-Adjusted DiD"], 3), "\n")
cat("  Difference:        ", 
    round(comparison_oop$estimate[comparison_oop$specification == "Simple DiD"] - 
            comparison_oop$estimate[comparison_oop$specification == "Trend-Adjusted DiD"], 3), "\n\n")

cat("INTERPRETATION:\n")
cat("The difference between Simple and Trend-Adjusted estimates represents\n")
cat("bias from ignoring pre-existing trends. Trend-adjusted estimate is more credible.\n\n")

cat("OUTPUT FILES CREATED:\n")
cat("  - trend_adjusted_did_coefficients.csv\n")
cat("  - comparison_simple_vs_trend_adjusted.csv\n")
cat("  - trend_comparison_simple_vs_adjusted.png\n")
cat("  - trend_adjusted_PRIMARY_comparison_table.docx\n")
cat("  - trend_adjusted_ALL_outcomes_table.docx\n")
cat("  - trend_adjusted_FALSIFICATION_table.docx\n\n")

cat("RECOMMENDATION FOR PAPER:\n")
cat("1. Report trend-adjusted estimates as PRIMARY specification\n")
cat("2. Show simple DiD in robustness/appendix\n")
cat("3. Acknowledge pre-trends transparently\n")
cat("4. Emphasize: 'Our estimates capture both direct policy effect\n")
cat("   and acceleration of pre-existing market trends'\n")
cat("5. Note falsification tests show insulin-specific patterns\n\n")

# -----------------------------------------------------------------------------
# 10. SYNTHETIC CONTROL WITH PROPENSITY SCORE MATCHING
# -----------------------------------------------------------------------------

cat("\n=== IMPLEMENTING SYNTHETIC CONTROL WITH PSM (1:1 MATCHING) ===\n\n")

library(MatchIt)
library(cobalt)

# Create patient baseline dataset for matching
patient_baseline <- monthly_panel %>%
  filter(year_month < as.Date("2023-01-01")) %>%  # Pre-IRA data only for baseline
  group_by(pat_id, treatment) %>%
  summarise(
    # Clinical characteristics
    n_charlson = mean(n_charlson, na.rm = TRUE),
    dcsi_total = mean(dcsi_total, na.rm = TRUE),
    n_unique_drugs_excl_insulin = mean(n_unique_drugs_excl_insulin, na.rm = TRUE),
    
    # Pre-period medication utilization
    pre_insulin_use = mean(!is.na(insulin_oop), na.rm = TRUE),
    pre_insulin_copay = mean(insulin_copay, na.rm = TRUE),
    pre_insulin_oop = mean(insulin_oop, na.rm = TRUE),
    pre_insulin_supplies = mean(insulin_standardized_supplies, na.rm = TRUE),
    pre_insulin_adherence = 1 - mean(insulin_gap, na.rm = TRUE),
    
    # Other medication use as proxy for disease burden
    pre_lipid_use = mean(!is.na(lipid_copay), na.rm = TRUE),
    pre_metformin_use = mean(!is.na(metformin_copay), na.rm = TRUE),
    
    # Demographic information
    der_sex = first(der_sex),
    age_bin = first(age_bin),
    
    # Comorbidity bins (for possible exact matching)
    charlson_bin = first(charlson_bin),
    dcsi_bin = first(dcsi_bin),
    polypharm_bin = first(polypharm_bin),
    
    # Pre-period outcome trends (for matching on pre-trends)
    pre_oop_trend = {
      if(n() >= 3) {
        mod <- lm(insulin_oop_per_supply ~ as.numeric(year_month))
        coef(mod)[2]
      } else {
        NA_real_
      }
    },
    
    # Pre-period observation length
    observation_months = n(),
    
    .groups = "drop"
  ) %>%
  # Ensure no NAs in matching variables
  mutate(across(where(is.numeric), ~replace_na(., mean(., na.rm = TRUE))))

# Display summary statistics before matching
cat("Sample size before matching:\n")
print(table(patient_baseline$treatment))

# Define formula for matching - include pre-treatment outcomes & trends
match_formula <- treatment ~ n_charlson + dcsi_total + n_unique_drugs_excl_insulin +
  pre_insulin_use + pre_insulin_copay + pre_insulin_oop + pre_insulin_supplies +
  pre_lipid_use + pre_metformin_use + der_sex + pre_oop_trend

# Perform propensity score matching
match_model <- matchit(
  match_formula,
  data = patient_baseline,
  method = "nearest",  # Nearest neighbor matching
  distance = "glm",    # Logistic regression for PS
  caliper = 0.2,       # Maximum distance (in SD of logit of PS)
  ratio = 1            # 1:1 matching
)

# Display matching summary
cat("\nMatching summary:\n")
print(summary(match_model))

# Check balance with cobalt
bal_tab <- bal.tab(match_model, un = TRUE,
                   stats = c("mean.diffs", "variance.ratios"))
print(bal_tab)

# Create balance plot
love_plot <- love.plot(match_model, 
                       stats = "mean.diffs",
                       threshold = 0.1,
                       abs = TRUE,
                       var.order = "unadjusted",
                       title = "Balance: Propensity Score Matching")
print(love_plot)
ggsave("psm_balance_plot.png", love_plot, width = 10, height = 8, dpi = 300)

# Extract matched dataset
matched_data <- match.data(match_model)
cat("\nMatched dataset created with", nrow(matched_data), "patients\n")

# Get list of matched patient IDs to filter panel dataset
matched_ids <- matched_data %>% select(pat_id)

# Filter original panel data to include only matched patients
matched_panel <- monthly_panel %>%
  inner_join(matched_ids, by = "pat_id")

cat("Original panel observations:", nrow(monthly_panel), "\n")
cat("Matched panel observations:", nrow(matched_panel), "\n")

# Create pre-trend plot for matched sample
plot_pretrends <- function(data, title) {
  trend_data <- data %>%
    filter(year_month < as.Date("2023-01-01")) %>%  # Pre-IRA period only
    group_by(year_month, treatment) %>%
    summarize(
      mean_oop = mean(insulin_oop_per_supply, na.rm = TRUE),
      se = sd(insulin_oop_per_supply, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      ci_lower = mean_oop - 1.96 * se,
      ci_upper = mean_oop + 1.96 * se,
      group_label = ifelse(treatment == 1, "Medicare Advantage (65+)", "Commercial (18-64)")
    )
  
  p <- ggplot(trend_data, aes(x = year_month, y = mean_oop, 
                              color = group_label, group = group_label)) +
    geom_vline(xintercept = as.Date("2023-01-01"), linetype = "dashed", color = "red") +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = group_label),
                alpha = 0.2, color = NA) +
    scale_color_manual(values = c("Medicare Advantage (65+)" = "#2c7bb6",
                                  "Commercial (18-64)" = "#d7191c")) +
    scale_fill_manual(values = c("Medicare Advantage (65+)" = "#2c7bb6",
                                 "Commercial (18-64)" = "#d7191c")) +
    labs(
      title = title,
      subtitle = "Pre-treatment period only",
      x = "Month",
      y = "Average Insulin OOP per 30-Day Supply ($)",
      color = "Group",
      fill = "Group"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "top"
    )
  
  return(p)
}

# Create pre-trend plots
pretrend_psm <- plot_pretrends(matched_panel, "Pre-Trends: Propensity Score Matched Sample")
ggsave("matched_pretrends_psm.png", pretrend_psm, width = 10, height = 6, dpi = 300)

# Compare unmatched vs matched pre-trends
pretrend_unmatched <- plot_pretrends(monthly_panel, "Pre-Trends: Unmatched Sample")

library(patchwork)
combined_pretrends <- (pretrend_unmatched / pretrend_psm) +
  plot_annotation(
    title = "Pre-Trend Comparison: Unmatched vs Matched Samples",
    theme = theme(plot.title = element_text(face = "bold", size = 16))
  )

ggsave("comparison_pretrends.png", combined_pretrends,
       width = 12, height = 12, dpi = 300)

cat("\nMatched sample created successfully. Pre-trend plots saved.\n")
cat("Next step: Run trend-adjusted DiD on matched sample.\n")

# -----------------------------------------------------------------------------
# 11. RUN TREND-ADJUSTED DID ON MATCHED SAMPLE
# -----------------------------------------------------------------------------

cat("=== RUNNING TREND-ADJUSTED DID MODELS ON MATCHED SAMPLE ===\n\n")

matched_results <- map(primary_outcomes, function(outcome) {
  cat("Processing matched sample:", outcome, "\n")
  run_trend_adjusted_did(matched_panel, outcome, binned_controls)
})
names(matched_results) <- primary_outcomes

# Combine coefficients
matched_coefs <- map_df(matched_results, ~.x$coefs)

# Compare original vs. matched sample results
combined_results <- bind_rows(
  trend_adjusted_coefs %>% 
    filter(outcome %in% primary_outcomes) %>%
    mutate(sample = "Full Sample"),
  matched_coefs %>% 
    mutate(sample = "PSM Matched")
) %>%
  select(sample, outcome, term, estimate, std.error, p.value, sig)

write_csv(combined_results, "comparison_full_vs_matched_sample.csv")

# Print comparison
cat("\n=== COMPARISON: FULL SAMPLE VS. PSM MATCHED SAMPLE ===\n\n")
print(combined_results %>%
        pivot_wider(names_from = sample, 
                    values_from = c(estimate, std.error, p.value, sig),
                    names_sep = "_"))

# -----------------------------------------------------------------------------
# 12. CREATE FINAL VISUALIZATION COMPARING FULL VS MATCHED SAMPLES
# -----------------------------------------------------------------------------

# Calculate average monthly values for original vs. matched samples
monthly_comparison <- bind_rows(
  monthly_panel %>%
    group_by(year_month, treatment) %>%
    summarise(across(all_of(primary_outcomes), 
                     ~mean(., na.rm = TRUE)),
              sample = "Full Sample",
              .groups = "drop"),
  
  matched_panel %>%
    group_by(year_month, treatment) %>%
    summarise(across(all_of(primary_outcomes), 
                     ~mean(., na.rm = TRUE)),
              sample = "PSM Matched",
              .groups = "drop")
) %>%
  mutate(
    group_label = ifelse(treatment == 1, "Medicare Advantage (65+)", "Commercial (18-64)"),
    post = as.numeric(year_month >= as.Date("2023-01-01"))
  )

# Plot for primary outcome
plot_matched_comparison <- function(outcome_name) {
  plot_data <- monthly_comparison %>%
    pivot_longer(cols = all_of(outcome_name), 
                 names_to = "outcome", 
                 values_to = "value")
  
  ggplot(plot_data, 
         aes(x = year_month, y = value, color = group_label, 
             linetype = sample, shape = sample)) +
    geom_vline(xintercept = as.Date("2023-01-01"), 
               linetype = "dashed", color = "red", linewidth = 1) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5, alpha = 0.7) +
    
    scale_color_manual(
      values = c("Medicare Advantage (65+)" = "#2c7bb6",
                 "Commercial (18-64)" = "#d7191c"),
      name = "Group"
    ) +
    
    scale_linetype_manual(
      values = c("Full Sample" = "solid",
                 "PSM Matched" = "longdash"),
      name = "Sample"
    ) +
    
    scale_shape_manual(
      values = c("Full Sample" = 16,
                 "PSM Matched" = 17),
      name = "Sample"
    ) +
    
    labs(
      title = paste("Full Sample vs. PSM Matched:", outcome_name),
      subtitle = "Comparing time trends before and after matching",
      x = "Month",
      y = paste("Mean", outcome_name),
      caption = "PSM matching creates more comparable treatment and control groups"
    ) +
    
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "bottom",
      legend.box = "vertical",
      panel.grid.minor = element_blank()
    )
}

# Create comparison plots
plot_matched_oop <- plot_matched_comparison("insulin_oop_per_supply")
plot_matched_copay <- plot_matched_comparison("insulin_copay_per_supply")

# Save plots
ggsave("synthetic_control_oop_comparison.png", plot_matched_oop, 
       width = 12, height = 8, dpi = 300)
ggsave("synthetic_control_copay_comparison.png", plot_matched_copay, 
       width = 12, height = 8, dpi = 300)

cat("Synthetic control comparison plots saved\n")

# Create table comparing estimates
matched_comparison_table <- combined_results %>%
  filter(term == "treat_post") %>%
  select(sample, outcome, estimate, std.error, p.value) %>%
  pivot_wider(
    id_cols = outcome,
    names_from = sample,
    values_from = c(estimate, std.error, p.value)
  )

# Check the column names to ensure we're using the right ones
print(names(matched_comparison_table))

# Assuming column names will be "estimate_Full Sample" and "estimate_PSM Matched"
# (with spaces in the names)
matched_comparison_table <- matched_comparison_table %>%
  mutate(
    outcome_label = case_when(
      outcome == "insulin_oop_per_supply" ~ "OOP per 30-day supply",
      outcome == "insulin_copay_per_supply" ~ "Copay per 30-day supply",
      TRUE ~ outcome
    ),
    # Use backticks around column names with spaces
    difference = `estimate_PSM Matched` - `estimate_Full Sample`,
    pct_difference = (difference / `estimate_Full Sample`) * 100
  )

write_csv(matched_comparison_table, "synthetic_control_estimates_comparison.csv")

cat("\n=== SYNTHETIC CONTROL RESULTS SUMMARY ===\n\n")
print(matched_comparison_table)

cat("\nSynthetic control analysis complete. All files saved.\n")

# Moiz Code----

## Improved balance checking----
# After running your PSM matching, add this:
# Improved balance checking with cobalt
bal_tab <- bal.tab(match_model, data = patient_baseline, un = TRUE,
                   stats = c("mean.diffs", "variance.ratios"))
print(bal_tab)

# Enhanced love plot
love_plot <- love.plot(match_model, data = patient_baseline, 
                       stats = "mean.diffs",
                       threshold = 0.1,
                       abs = TRUE,
                       var.order = "unadjusted",
                       title = "Balance: Propensity Score Matching")
ggsave("psm_balance_enhanced_plot.png", love_plot, width = 10, height = 8, dpi = 300)

## Add CEM as Robustness Check----
# Add Coarsened Exact Matching as a robustness check
# === ALTERNATIVE MATCHING: COARSENED EXACT MATCHING ===

cat("\n=== ALTERNATIVE MATCHING: COARSENED EXACT MATCHING ===\n\n")

# First, let's check the distribution of the bin variables
cat("Distribution of categorical variables:\n")
print(table(patient_baseline$charlson_bin, patient_baseline$treatment, dnn = c("charlson_bin", "treatment")))
print(table(patient_baseline$dcsi_bin, patient_baseline$treatment, dnn = c("dcsi_bin", "treatment")))
print(table(patient_baseline$polypharm_bin, patient_baseline$treatment, dnn = c("polypharm_bin", "treatment")))

# Try a more flexible approach with fewer variables
# Option 1: Use fewer variables
match_cem_fewer <- matchit(
  treatment ~ charlson_bin + dcsi_bin + pre_insulin_use,
  data = patient_baseline,
  method = "cem",
  k2k = TRUE
)

# Check if any matches were found
if ("subclass" %in% names(match_cem_fewer$match.matrix)) {
  cat("CEM matching with fewer variables successful!\n")
  print(summary(match_cem_fewer))
  
  # Extract matched sample
  matched_cem <- match.data(match_cem_fewer)
  cat("CEM matched dataset created with", nrow(matched_cem), "patients\n")
} else {
  cat("Still no matches found with fewer variables.\n")
  
  # Option 2: Try without k2k (exact pair matching) for more flexibility
  match_cem_flex <- matchit(
    treatment ~ charlson_bin + dcsi_bin + pre_insulin_use,
    data = patient_baseline,
    method = "cem",
    k2k = FALSE  # Allow multiple matches
  )
  
  if ("subclass" %in% names(match_cem_flex$match.matrix)) {
    cat("CEM matching without k2k successful!\n")
    print(summary(match_cem_flex))
    
    # Extract matched sample
    matched_cem <- match.data(match_cem_flex)
    cat("CEM matched dataset created with", nrow(matched_cem), "patients\n")
  } else {
    cat("Still no matches found. Trying coarsened matching with bins.\n")
    
    # Option 3: Create coarser bins for the variables
    patient_baseline <- patient_baseline %>%
      mutate(
        charlson_simple = case_when(
          n_charlson <= 1 ~ "Low",
          n_charlson <= 3 ~ "Medium",
          TRUE ~ "High"
        ),
        dcsi_simple = case_when(
          dcsi_total <= 1 ~ "Low",
          dcsi_total <= 4 ~ "Medium", 
          TRUE ~ "High"
        )
      )
    
    match_cem_coarse <- matchit(
      treatment ~ charlson_simple + dcsi_simple + pre_insulin_use,
      data = patient_baseline,
      method = "cem"
    )
    
    if (length(match_cem_coarse$match.matrix) > 0) {
      cat("CEM matching with coarser bins successful!\n")
      print(summary(match_cem_coarse))
      
      # Extract matched sample
      matched_cem <- match.data(match_cem_coarse)
      cat("CEM matched dataset created with", nrow(matched_cem), "patients\n")
    } else {
      # Fallback: Use Mahalanobis distance matching as alternative
      cat("CEM matching still unsuccessful. Falling back to Mahalanobis matching.\n")
      
      match_mahal <- matchit(
        treatment ~ n_charlson + dcsi_total + n_unique_drugs_excl_insulin + 
          pre_insulin_use + pre_insulin_oop + pre_insulin_copay,
        data = patient_baseline,
        method = "nearest",
        distance = "mahalanobis",
        ratio = 1
      )
      
      print(summary(match_mahal))
      
      # Use Mahalanobis matched data instead
      matched_cem <- match.data(match_mahal)
      cat("Mahalanobis matched dataset created with", nrow(matched_cem), "patients as CEM alternative\n")
    }
  }
}

# Check balance for the successful matching approach
bal_tab_cem <- bal.tab(matched_cem, treat = matched_cem$treatment, 
                       un = TRUE, stats = c("mean.diffs", "variance.ratios"))
print(bal_tab_cem)

# Create love plot for the matching
love_plot_cem <- love.plot(bal_tab_cem, 
                           stats = "mean.diffs",
                           threshold = 0.1,
                           abs = TRUE)
print(love_plot_cem)
ggsave("cem_balance_plot.png", love_plot_cem, width = 10, height = 8, dpi = 300)

# Filter panel data to include only matched patients
matched_panel_cem <- monthly_panel %>%
  inner_join(matched_cem %>% select(pat_id), by = "pat_id")

cat("CEM matched panel observations:", nrow(matched_panel_cem), "\n")


## Pre-trend Viz----
# Add after your matching section:
# Visual pre-trend check in matched samples
cat("\n=== CHECKING PRE-TRENDS IN MATCHED SAMPLES ===\n\n")

# Function to create pre-trend plot
plot_pretrends <- function(data, title) {
  # Unweighted means by month and treatment group
  trend_data <- data %>%
    filter(year_month < as.Date("2023-01-01")) %>%  # Pre-IRA period only
    group_by(year_month, treatment) %>%
    summarize(
      mean_oop = mean(insulin_oop_per_supply, na.rm = TRUE),
      se = sd(insulin_oop_per_supply, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      ci_lower = mean_oop - 1.96 * se,
      ci_upper = mean_oop + 1.96 * se,
      group_label = ifelse(treatment == 1, "Medicare Advantage (65+)", "Commercial (18-64)")
    )
  
  p <- ggplot(trend_data, aes(x = year_month, y = mean_oop, 
                              color = group_label, group = group_label)) +
    geom_vline(xintercept = as.Date("2023-01-01"), linetype = "dashed", color = "red") +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = group_label),
                alpha = 0.2, color = NA) +
    scale_color_manual(values = c("Medicare Advantage (65+)" = "#2c7bb6",
                                  "Commercial (18-64)" = "#d7191c")) +
    scale_fill_manual(values = c("Medicare Advantage (65+)" = "#2c7bb6",
                                 "Commercial (18-64)" = "#d7191c")) +
    labs(
      title = title,
      subtitle = "Pre-treatment period only",
      x = "Month",
      y = "Average Insulin OOP per 30-Day Supply ($)",
      color = "Group",
      fill = "Group"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "top"
    )
  
  return(p)
}

# Create pre-trend plots
pretrend_unmatched <- plot_pretrends(monthly_panel, "Pre-Trends: Unmatched Sample")
pretrend_psm <- plot_pretrends(matched_panel, "Pre-Trends: Propensity Score Matched Sample")
pretrend_cem <- plot_pretrends(matched_panel_cem, "Pre-Trends: Coarsened Exact Matched Sample")

# Save pre-trend plots
ggsave("unmatched_pretrends.png", pretrend_unmatched, width = 10, height = 6, dpi = 300)
ggsave("matched_pretrends_psm.png", pretrend_psm, width = 10, height = 6, dpi = 300)
ggsave("matched_pretrends_cem.png", pretrend_cem, width = 10, height = 6, dpi = 300)

# Combined plot
library(patchwork)
combined_pretrends <- (pretrend_unmatched / pretrend_psm / pretrend_cem) +
  plot_annotation(
    title = "Pre-Trend Comparison Across Matching Methods",
    theme = theme(plot.title = element_text(face = "bold", size = 16))
  )

ggsave("matched_pretrends_comparison.png", combined_pretrends,
       width = 12, height = 15, dpi = 300)

## Comparison of Estimates Across all specs----
# After running your models on matched samples:
cat("\n=== COMPARING UNMATCHED VS MATCHED ESTIMATES ===\n\n")

# Create a function to run models on different samples
run_trend_adjusted_did <- function(data, outcome_var, controls = NULL) {
  # Define the formula
  formula_str <- paste0(outcome_var, " ~ treatment:post + treatment:months_since_start + treatment:months_since_ira")
  
  # Add controls if provided
  if (!is.null(controls)) {
    formula_str <- paste0(formula_str, " + ", paste(controls, collapse = " + "))
  }
  
  # Add fixed effects
  formula_str <- paste0(formula_str, " | pat_id + calendar_month")
  
  # Run the model
  model <- feols(as.formula(formula_str), data = data)
  
  # Extract coefficients and confidence intervals
  coefs <- tidy(model) %>%
    mutate(
      outcome = outcome_var,
      ci_lower = estimate - 1.96 * std.error,
      ci_upper = estimate + 1.96 * std.error,
      sig = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  return(list(model = model, coefs = coefs))
}

# Define primary outcomes
primary_outcomes <- c("insulin_oop_per_supply", "insulin_copay_per_supply")

# Define controls if needed
binned_controls <- c("charlson_bin", "dcsi_bin", "polypharm_bin")

# Run models for unmatched data if not already done
trend_adjusted_results <- list()
trend_adjusted_coefs <- data.frame()

for (outcome in primary_outcomes) {
  result <- run_trend_adjusted_did(monthly_panel, outcome, binned_controls)
  trend_adjusted_results[[outcome]] <- result
  trend_adjusted_coefs <- bind_rows(trend_adjusted_coefs, result$coefs)
}

# Run models for PSM matched data
matched_results <- list()
matched_coefs <- data.frame()

for (outcome in primary_outcomes) {
  result <- run_trend_adjusted_did(matched_panel, outcome, binned_controls)
  matched_results[[outcome]] <- result
  matched_coefs <- bind_rows(matched_coefs, result$coefs)
}

# Run models for CEM matched data
cem_results <- list()
cem_coefs <- data.frame()

for (outcome in primary_outcomes) {
  result <- run_trend_adjusted_did(matched_panel_cem, outcome, binned_controls)
  cem_results[[outcome]] <- result
  cem_coefs <- bind_rows(cem_coefs, result$coefs)
}

# Combine trend-adjusted results from all samples
combined_results <- bind_rows(
  trend_adjusted_coefs %>% 
    mutate(sample = "Full Sample"),
  
  matched_coefs %>% 
    mutate(sample = "PSM Matched"),
  
  cem_coefs %>% 
    mutate(sample = "CEM Matched")
) %>%
  mutate(
    sample = factor(sample, levels = c("Full Sample", "PSM Matched", "CEM Matched"))
  )

write_csv(combined_results, "matched_comparison_all_methods.csv")

# === COMPARING UNMATCHED VS MATCHED ESTIMATES ===

# Create a function to run models on different samples
run_trend_adjusted_did <- function(data, outcome_var, controls = NULL) {
  # Define the formula
  formula_str <- paste0(outcome_var, " ~ treatment:post + treatment:months_since_start + treatment:months_since_ira")
  
  # Add controls if provided
  if (!is.null(controls)) {
    formula_str <- paste0(formula_str, " + ", paste(controls, collapse = " + "))
  }
  
  # Add fixed effects
  formula_str <- paste0(formula_str, " | pat_id + calendar_month")
  
  # Run the model
  model <- feols(as.formula(formula_str), data = data)
  
  # Extract coefficients and confidence intervals
  coefs <- tidy(model) %>%
    mutate(
      outcome = outcome_var,
      ci_lower = estimate - 1.96 * std.error,
      ci_upper = estimate + 1.96 * std.error,
      sig = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  return(list(model = model, coefs = coefs))
}

# Define primary outcomes
primary_outcomes <- c("insulin_oop_per_supply", "insulin_copay_per_supply")

# Define controls if needed
binned_controls <- c("charlson_bin", "dcsi_bin", "polypharm_bin")

# Run models for unmatched data if not already done
trend_adjusted_results <- list()
trend_adjusted_coefs <- data.frame()

for (outcome in primary_outcomes) {
  result <- run_trend_adjusted_did(monthly_panel, outcome, binned_controls)
  trend_adjusted_results[[outcome]] <- result
  trend_adjusted_coefs <- bind_rows(trend_adjusted_coefs, result$coefs)
}

# Run models for PSM matched data
matched_results <- list()
matched_coefs <- data.frame()

for (outcome in primary_outcomes) {
  result <- run_trend_adjusted_did(matched_panel, outcome, binned_controls)
  matched_results[[outcome]] <- result
  matched_coefs <- bind_rows(matched_coefs, result$coefs)
}

# Run models for CEM matched data
cem_results <- list()
cem_coefs <- data.frame()

for (outcome in primary_outcomes) {
  result <- run_trend_adjusted_did(matched_panel_cem, outcome, binned_controls)
  cem_results[[outcome]] <- result
  cem_coefs <- bind_rows(cem_coefs, result$coefs)
}

# Combine trend-adjusted results from all samples
combined_results <- bind_rows(
  trend_adjusted_coefs %>% 
    mutate(sample = "Full Sample"),
  
  matched_coefs %>% 
    mutate(sample = "PSM Matched"),
  
  cem_coefs %>% 
    mutate(sample = "CEM Matched")
) %>%
  mutate(
    sample = factor(sample, levels = c("Full Sample", "PSM Matched", "CEM Matched"))
  )

write_csv(combined_results, "matched_comparison_all_methods.csv")

# First, let's check what terms are available in the combined results
print(head(combined_results))
print(table(combined_results$term, combined_results$outcome))

# Create comparison plot with better error handling
comparison_plot <- combined_results %>%
  # Check if 'treatment:post' exists, otherwise try 'treat_post'
  filter((term == "treatment:post" | term == "treatmentTRUE:post") & 
           outcome %in% c("insulin_oop_per_supply", "insulin_copay_per_supply")) %>%
  # Print to debug
  {
    print(paste("Filtered rows:", nrow(.)))
    .
  } %>%
  # Create readable outcome labels
  mutate(outcome_label = case_when(
    outcome == "insulin_oop_per_supply" ~ "OOP per 30-day supply",
    outcome == "insulin_copay_per_supply" ~ "Copay per 30-day supply",
    TRUE ~ outcome
  )) %>%
  # Ensure the outcome_label has at least one value
  {
    if(nrow(.) == 0) {
      # If no matching rows, create a dummy plot with a message
      cat("No matching rows found for the plot. Creating dummy plot with message.\n")
      data.frame(
        sample = c("Full Sample"),
        estimate = 0,
        ci_lower = 0,
        ci_upper = 0,
        outcome_label = "No matching treatment effect terms found"
      )
    } else {
      # Print the outcome labels to verify they exist
      cat("Outcome labels found:", paste(unique(.$outcome_label), collapse=", "), "\n")
      .
    }
  } %>%
  # Create the plot
  ggplot(aes(x = sample, y = estimate, color = sample)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.3, linewidth = 1) +
  scale_color_manual(values = c("Full Sample" = "black",
                                "PSM Matched" = "#d7191c",
                                "CEM Matched" = "#fdae61")) +
  facet_wrap(~ outcome_label, scales = "free") +
  labs(
    title = "Treatment Effect Estimates: Unmatched vs Matched Samples",
    subtitle = "Trend-Adjusted DiD Models with Different Matching Approaches",
    x = "Method",
    y = "Treatment Effect Estimate ($)",
    caption = "95% confidence intervals shown. All models include patient FE + controls."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Save the plot
ggsave("matched_comparison_plot.png", comparison_plot,
       width = 12, height = 7, dpi = 300)

# Let's also check the exact coefficients to verify the matching impact
comparison_table <- combined_results %>%
  filter((term == "treatment:post" | term == "treatmentTRUE:post") & 
           outcome %in% c("insulin_oop_per_supply", "insulin_copay_per_supply")) %>%
  select(sample, outcome, estimate, std.error, ci_lower, ci_upper, sig) %>%
  mutate(
    outcome = case_when(
      outcome == "insulin_oop_per_supply" ~ "OOP per 30-day supply",
      outcome == "insulin_copay_per_supply" ~ "Copay per 30-day supply",
      TRUE ~ outcome
    ),
    formatted_est = sprintf("%.2f%s (%.2f, %.2f)", 
                            estimate, sig, ci_lower, ci_upper)
  ) %>%
  select(sample, outcome, formatted_est) %>%
  pivot_wider(names_from = outcome, values_from = formatted_est) %>%
  arrange(sample)

print(comparison_table)

# Alternative approach: If we can't identify the right terms, let's create a table of all terms
all_terms_table <- combined_results %>%
  group_by(sample, term, outcome) %>%
  summarize(
    estimate = first(estimate),
    ci_lower = first(ci_lower),
    ci_upper = first(ci_upper),
    sig = first(sig),
    .groups = "drop"
  ) %>%
  arrange(sample, outcome, term)

print(head(all_terms_table, 20))