# =============================================================================
# COMPLETE MONTHLY ANALYSIS - ALL COHORTS (WITH INTERACTED MONTH EFFECTS)
# =============================================================================
# Analysis for:
# - Original: 18-64 Commercial vs 65+ Medicare Advantage
# - Restrictive 1: 62-64 vs 65-67
# - Restrictive 2: 54-64 vs 65-75
# - Restrictive 3: 54-64 vs 65+
# - Restrictive 4: 62-64 vs 65+
#
# UPDATED: Added calendar month*treatment interactions (month-by-cohort FE)
# to allow MA and Commercial to have different seasonal patterns
# =============================================================================

library(tidyverse)
library(fixest)
library(broom)
library(modelsummary)
library(lmtest)
library(sandwich)

# -----------------------------------------------------------------------------
# 1. LOAD ALL COHORTS
# -----------------------------------------------------------------------------

cat("=== LOADING ALL COHORTS ===\n\n")

# ORIGINAL COHORTS (18-64 vs 65+)
monthly_cross_2022 <- read_csv("analytical_monthly_cross_2022.csv")
monthly_cross_2023 <- read_csv("analytical_monthly_cross_2023.csv")
monthly_panel <- read_csv("analytical_monthly_panel.csv")
monthly_panel_balanced <- read_csv("analytical_monthly_panel_balanced.csv")

# RESTRICTIVE COHORT 1: 62-64 vs 65-67
restrictive_62_64_vs_65_67_2022 <- read_csv("analytical_monthly_restrictive_2022.csv")
restrictive_62_64_vs_65_67_2023 <- read_csv("analytical_monthly_restrictive_2023.csv")

# RESTRICTIVE COHORT 2: 54-64 vs 65-75
restrictive_54_64_vs_65_75_2022 <- read_csv("analytical_monthly_restrictive_54_64_vs_65_75_2022.csv")
restrictive_54_64_vs_65_75_2023 <- read_csv("analytical_monthly_restrictive_54_64_vs_65_75_2023.csv")

# RESTRICTIVE COHORT 3: 54-64 vs 65+
restrictive_54_64_vs_65_plus_2022 <- read_csv("analytical_monthly_restrictive_54_64_vs_65_plus_2022.csv")
restrictive_54_64_vs_65_plus_2023 <- read_csv("analytical_monthly_restrictive_54_64_vs_65_plus_2023.csv")

# RESTRICTIVE COHORT 4: 62-64 vs 65+
restrictive_62_64_vs_65_plus_2022 <- read_csv("analytical_monthly_restrictive_62_64_vs_65_plus_2022.csv")
restrictive_62_64_vs_65_plus_2023 <- read_csv("analytical_monthly_restrictive_62_64_vs_65_plus_2023.csv")

# -----------------------------------------------------------------------------
# 2. DATA CLEANING FUNCTION (WITH MONTH DUMMIES)
# -----------------------------------------------------------------------------

clean_data <- function(data) {
  data %>%
    mutate(
      # Replace Inf/-Inf with NA
      across(where(is.numeric), ~ifelse(is.infinite(.), NA, .)),
      
      # For per-supply metrics: keep NA (indicates no fills that month)
      # For total costs: replace NA with 0
      across(c(insulin_copay, insulin_deductible, insulin_coinsamt, insulin_cobamt, 
               insulin_oop, total_copay, total_oop,
               starts_with("copay_"), starts_with("oop_"),
               metformin_copay, metformin_oop, lipid_copay, lipid_oop), 
             ~replace_na(., 0)),
      
      # For count/supply outcomes: replace NA with 0
      across(c(insulin_n_fills, insulin_quantity, insulin_standardized_supplies,
               total_fills, starts_with("n_fills_"), starts_with("total_quan_"),
               starts_with("supplies_"),
               metformin_n_fills, lipid_n_fills),
             ~replace_na(., 0)),
      
      # For adherence flags
      insulin_gap = replace_na(insulin_gap, TRUE),
      fully_covered = replace_na(fully_covered, FALSE),
      estimated_days_covered = replace_na(estimated_days_covered, 0),
      
      # CRITICAL: Add calendar month as factor for seasonality controls
      calendar_month = factor(month(year_month), levels = 1:12, labels = month.abb)
    )
}

# Combine and clean all cohorts
all_cohorts <- list(
  # Original
  "original_18_64_vs_65_plus" = bind_rows(monthly_cross_2022, monthly_cross_2023),
  
  # Restrictive cohorts
  "restrictive_62_64_vs_65_67" = bind_rows(restrictive_62_64_vs_65_67_2022, restrictive_62_64_vs_65_67_2023),
  "restrictive_54_64_vs_65_75" = bind_rows(restrictive_54_64_vs_65_75_2022, restrictive_54_64_vs_65_75_2023),
  "restrictive_54_64_vs_65_plus" = bind_rows(restrictive_54_64_vs_65_plus_2022, restrictive_54_64_vs_65_plus_2023),
  "restrictive_62_64_vs_65_plus" = bind_rows(restrictive_62_64_vs_65_plus_2022, restrictive_62_64_vs_65_plus_2023)
)

# Clean and prepare all cohorts
all_cohorts <- map(all_cohorts, function(data) {
  clean_data(data) %>%
    mutate(
      treatment = as.numeric(treatment),
      post = as.numeric(post_ira),
      treat_post = treatment * post
    )
})

# Also prepare panel data (original only)
monthly_panel <- monthly_panel %>%
  clean_data() %>%
  mutate(
    treatment = as.numeric(treatment),
    post = as.numeric(post_ira),
    treat_post = treatment * post
  )

monthly_panel_balanced <- monthly_panel_balanced %>%
  clean_data() %>%
  mutate(
    treatment = as.numeric(treatment),
    post = as.numeric(post_ira),
    treat_post = treatment * post
  )

# -----------------------------------------------------------------------------
# 3. DEFINE OUTCOMES AND CONTROLS
# -----------------------------------------------------------------------------

# PRIMARY OUTCOMES (per 30-day supply)
primary_cost_outcomes <- c(
  "insulin_oop_per_supply",
  "insulin_copay_per_supply",
  "insulin_deductible",
  "insulin_coinsamt",
  "insulin_cobamt"
)

# ADHERENCE OUTCOMES
adherence_outcomes <- c(
  "insulin_standardized_supplies",
  "insulin_gap",
  "fully_covered",
  "estimated_days_covered"
)

# SECONDARY OUTCOMES (monthly totals)
secondary_cost_outcomes <- c(
  "insulin_copay",
  "insulin_oop"
)

# UTILIZATION OUTCOMES
utilization_outcomes <- c(
  "insulin_n_fills",
  "insulin_quantity"
)

# FALSIFICATION OUTCOMES
falsification_outcomes <- c(
  "metformin_copay", "metformin_oop", "metformin_n_fills",
  "lipid_copay", "lipid_oop", "lipid_n_fills"
)

# Binned controls (time-varying for panel FE)
binned_controls <- c("charlson_bin", "dcsi_bin", "polypharm_bin")

# All outcomes for analysis
all_outcomes <- c(primary_cost_outcomes, adherence_outcomes, 
                  secondary_cost_outcomes, utilization_outcomes)

# -----------------------------------------------------------------------------
# 4. RUN DID MODELS FOR ALL COHORTS (CROSS-SECTIONAL WITH INTERACTED MONTHS)
# -----------------------------------------------------------------------------

cat("\n=== RUNNING DID MODELS FOR ALL COHORTS (WITH INTERACTED MONTH EFFECTS) ===\n\n")

run_monthly_did_interacted <- function(data, outcome, controls) {
  # Month-by-cohort FE: calendar_month*treatment
  # Allows MA and Commercial to have different seasonal patterns
  formula_str <- paste0(outcome, " ~ treatment + post + treat_post + calendar_month*treatment + ", 
                        paste(controls, collapse = " + "))
  
  model <- lm(as.formula(formula_str), data = data)
  
  did_estimate <- tidy(model) %>%
    filter(term == "treat_post") %>%
    mutate(outcome = outcome)
  
  return(list(model = model, did = did_estimate))
}

# Run DID for all cohorts and all outcomes
all_did_results_interacted <- map(names(all_cohorts), function(cohort_name) {
  cat("Processing cohort:", cohort_name, "\n")
  
  cohort_data <- all_cohorts[[cohort_name]]
  
  # Add age_bin and der_sex to controls for cross-sectional models
  cross_sectional_controls <- c("age_bin", "der_sex", binned_controls)
  
  # Run DID for all outcomes
  cohort_results <- map(all_outcomes, ~run_monthly_did_interacted(cohort_data, .x, cross_sectional_controls))
  names(cohort_results) <- all_outcomes
  
  # Extract DID estimates
  cohort_estimates <- map_df(cohort_results, ~.x$did) %>%
    mutate(
      cohort = cohort_name,
      ci_lower = estimate - 1.96 * std.error,
      ci_upper = estimate + 1.96 * std.error,
      sig = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        TRUE ~ ""
      ),
      outcome_type = case_when(
        outcome %in% primary_cost_outcomes ~ "Primary Cost (per supply)",
        outcome %in% adherence_outcomes ~ "Adherence",
        outcome %in% secondary_cost_outcomes ~ "Secondary Cost (monthly)",
        outcome %in% utilization_outcomes ~ "Utilization",
        TRUE ~ "Other"
      )
    )
  
  return(list(models = cohort_results, estimates = cohort_estimates))
})

names(all_did_results_interacted) <- names(all_cohorts)

# Combine all estimates
combined_did_estimates_interacted <- map_df(all_did_results_interacted, ~.x$estimates)

write_csv(combined_did_estimates_interacted, "monthly_all_cohorts_did_estimates_interacted_months.csv")

# -----------------------------------------------------------------------------
# 5. PANEL FE MODELS WITH INTERACTED MONTH EFFECTS (ORIGINAL 22-MONTH PANEL)
# -----------------------------------------------------------------------------

cat("\n=== RUNNING PANEL FE MODELS WITH INTERACTED MONTH EFFECTS (22 MONTHS) ===\n\n")

run_panel_fe_with_interacted_months <- function(data, outcome, controls) {
  # Month-by-cohort FE: Allows MA and Commercial to have different seasonal patterns
  # This addresses potential violations of parallel trends due to:
  # - Different formulary reset timing
  # - Part D donut hole effects (MA only)
  # - Different end-of-year behaviors
  
  formula_str <- paste0(outcome, " ~ treat_post + calendar_month*treatment + ", 
                        paste(controls, collapse = " + "),
                        " | pat_id")
  
  model <- feols(as.formula(formula_str), data = data, cluster = ~pat_id)
  
  did_estimate <- tidy(model) %>%
    filter(term == "treat_post") %>%
    mutate(outcome = outcome)
  
  return(list(model = model, did = did_estimate))
}

# Run for all outcomes
panel_monthly_fe_interacted <- map(all_outcomes, ~run_panel_fe_with_interacted_months(monthly_panel, .x, binned_controls))
names(panel_monthly_fe_interacted) <- all_outcomes

panel_monthly_estimates_interacted <- map_df(panel_monthly_fe_interacted, ~.x$did) %>%
  mutate(
    cohort = "panel_22_months_interacted",
    ci_lower = estimate - 1.96 * std.error,
    ci_upper = estimate + 1.96 * std.error,
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    outcome_type = case_when(
      outcome %in% primary_cost_outcomes ~ "Primary Cost (per supply)",
      outcome %in% adherence_outcomes ~ "Adherence",
      outcome %in% secondary_cost_outcomes ~ "Secondary Cost (monthly)",
      outcome %in% utilization_outcomes ~ "Utilization",
      TRUE ~ "Other"
    )
  )

write_csv(panel_monthly_estimates_interacted, "monthly_panel_fe_estimates_interacted_months.csv")

# -----------------------------------------------------------------------------
# 6. COMPARE SIMPLE VS INTERACTED MONTH SPECIFICATIONS
# -----------------------------------------------------------------------------

cat("\n=== COMPARING MONTH SPECIFICATIONS ===\n\n")

# Run simple month dummies for comparison
run_panel_fe_with_simple_months <- function(data, outcome, controls) {
  formula_str <- paste0(outcome, " ~ treat_post + calendar_month + ", 
                        paste(controls, collapse = " + "),
                        " | pat_id")
  
  model <- feols(as.formula(formula_str), data = data, cluster = ~pat_id)
  
  did_estimate <- tidy(model) %>%
    filter(term == "treat_post") %>%
    mutate(outcome = outcome)
  
  return(list(model = model, did = did_estimate))
}

panel_monthly_fe_simple <- map(all_outcomes, ~run_panel_fe_with_simple_months(monthly_panel, .x, binned_controls))
names(panel_monthly_fe_simple) <- all_outcomes

panel_monthly_estimates_simple <- map_df(panel_monthly_fe_simple, ~.x$did) %>%
  mutate(
    cohort = "panel_22_months_simple",
    ci_lower = estimate - 1.96 * std.error,
    ci_upper = estimate + 1.96 * std.error,
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

write_csv(panel_monthly_estimates_simple, "monthly_panel_fe_estimates_simple_months.csv")

# Create comparison table
comparison_specs <- bind_rows(
  panel_monthly_estimates_simple %>% 
    filter(outcome %in% c("insulin_oop_per_supply", "insulin_copay_per_supply",
                          "insulin_standardized_supplies", "insulin_gap")) %>%
    mutate(specification = "Simple month dummies"),
  panel_monthly_estimates_interacted %>% 
    filter(outcome %in% c("insulin_oop_per_supply", "insulin_copay_per_supply",
                          "insulin_standardized_supplies", "insulin_gap")) %>%
    mutate(specification = "Interacted month-by-cohort FE")
) %>%
  select(specification, outcome, estimate, std.error, ci_lower, ci_upper, p.value, sig)

print(comparison_specs)

write_csv(comparison_specs, "monthly_panel_specification_comparison.csv")

cat("\nInterpretation:\n")
cat("- Simple month dummies: Assumes MA and Commercial have same seasonal patterns\n")
cat("- Interacted (month-by-cohort): Allows different seasonal patterns by insurance type\n")
cat("- If estimates are similar: Seasonality not driving results\n")
cat("- If estimates differ substantially: Seasonal heterogeneity matters\n\n")

# -----------------------------------------------------------------------------
# 7. EXTRACT GROUP-SPECIFIC SEASONALITY PATTERNS
# -----------------------------------------------------------------------------

cat("\n=== EXTRACTING GROUP-SPECIFIC MONTH EFFECTS ===\n\n")

extract_group_month_effects <- function(model) {
  coefs <- tidy(model)
  
  # Base month effects (Commercial group, reference = January)
  commercial_months <- coefs %>%
    filter(str_detect(term, "^calendar_month[A-Z]") & !str_detect(term, ":")) %>%
    mutate(
      month = str_remove(term, "calendar_month"),
      group = "Commercial",
      estimate_adj = estimate
    ) %>%
    select(month, group, estimate = estimate_adj, std.error, p.value)
  
  # Add January (reference month = 0)
  commercial_months <- bind_rows(
    tibble(month = "Jan", group = "Commercial", estimate = 0, std.error = 0, p.value = NA),
    commercial_months
  )
  
  # Interaction effects (difference for MA)
  ma_interactions <- coefs %>%
    filter(str_detect(term, "calendar_month.*:treatment$")) %>%
    mutate(
      month = str_extract(term, "(?<=calendar_month)[A-Z][a-z]{2}"),
      interaction = estimate
    ) %>%
    select(month, interaction, std.error_int = std.error)
  
  # Base treatment effect (MA vs Commercial in January)
  treatment_base <- coefs %>%
    filter(term == "treatment") %>%
    pull(estimate)
  
  # MA month effects = Commercial month effect + interaction
  ma_months <- commercial_months %>%
    left_join(ma_interactions, by = "month") %>%
    mutate(
      group = "Medicare Advantage",
      estimate = estimate + coalesce(interaction, 0),
      std.error = sqrt(std.error^2 + coalesce(std.error_int^2, 0))
    ) %>%
    select(month, group, estimate, std.error, p.value)
  
  # Combine
  all_months <- bind_rows(commercial_months, ma_months) %>%
    mutate(
      month = factor(month, levels = c("Jan", month.abb[2:12])),
      ci_lower = estimate - 1.96 * std.error,
      ci_upper = estimate + 1.96 * std.error
    )
  
  return(all_months)
}

# Extract for primary outcome
if("insulin_oop_per_supply" %in% names(panel_monthly_fe_interacted)) {
  month_effects_by_group <- extract_group_month_effects(
    panel_monthly_fe_interacted$insulin_oop_per_supply$model
  )
  
  write_csv(month_effects_by_group, "month_effects_by_insurance_group.csv")
  
  # Visualize
  plot_seasonality_by_group <- ggplot(month_effects_by_group, 
                                      aes(x = month, y = estimate, color = group, group = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, linewidth = 1) +
    scale_color_manual(values = c("Commercial" = "#d7191c", 
                                  "Medicare Advantage" = "#2c7bb6")) +
    labs(
      title = "Group-Specific Seasonal Patterns in Insulin OOP Costs",
      subtitle = "Month-by-cohort fixed effects reveal different seasonal patterns",
      x = "Month",
      y = "Effect on OOP per 30-Day Supply ($)",
      color = "Insurance Type",
      caption = "January is reference month (estimate = 0). 95% confidence intervals shown.\nInteracted model allows MA and Commercial to have different seasonality."
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "top",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.caption = element_text(size = 9, hjust = 0)
    )
  
  ggsave("seasonality_by_insurance_group.png", plot_seasonality_by_group, 
         width = 12, height = 8, dpi = 300)
  
  cat("Group-specific seasonality plot saved: seasonality_by_insurance_group.png\n")
  
  # Test if seasonal patterns differ significantly
  model <- panel_monthly_fe_interacted$insulin_oop_per_supply$model
  interaction_terms <- names(coef(model))[str_detect(names(coef(model)), "calendar_month.*:treatment")]
  
  if(length(interaction_terms) > 0) {
    joint_test <- tryCatch({
      wald(model, interaction_terms)
    }, error = function(e) NULL)
    
    if(!is.null(joint_test)) {
      cat("\n=== TEST: DO SEASONAL PATTERNS DIFFER BY GROUP? ===\n")
      cat("Joint F-test of all month × treatment interactions:\n")
      cat("F-statistic:", round(joint_test$stat, 3), "\n")
      cat("p-value:", format.pval(joint_test$p, digits = 3), "\n")
      if(joint_test$p < 0.05) {
        cat("CONCLUSION: Seasonal patterns DIFFER significantly between MA and Commercial (p < 0.05)\n")
        cat("→ Interacted month FE specification is appropriate\n\n")
      } else {
        cat("CONCLUSION: Seasonal patterns do NOT differ significantly (p >= 0.05)\n")
        cat("→ Simple month dummies may be sufficient\n\n")
      }
    }
  }
}

# -----------------------------------------------------------------------------
# 8. CAP ANALYSIS FOR ALL COHORTS
# -----------------------------------------------------------------------------

cat("\n=== CALCULATING % EXCEEDING $35 PER-SUPPLY CAP (ALL COHORTS) ===\n\n")

cap_analysis_all <- map_df(names(all_cohorts), function(cohort_name) {
  all_cohorts[[cohort_name]] %>%
    filter(insulin_n_fills > 0) %>%
    group_by(year, treatment) %>%
    summarise(
      n_patient_months = n(),
      n_exceeding_cap = sum(exceeds_cap, na.rm = TRUE),
      pct_exceeding_cap = mean(exceeds_cap, na.rm = TRUE) * 100,
      mean_oop_per_supply = mean(insulin_oop_per_supply[insulin_n_fills > 0], na.rm = TRUE),
      median_oop_per_supply = median(insulin_oop_per_supply[insulin_n_fills > 0], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(cohort = cohort_name)
})

write_csv(cap_analysis_all, "monthly_cap_analysis_all_cohorts.csv")

# DID for cap exceeding (all cohorts) - with interacted months
cap_did_all_interacted <- map_df(names(all_cohorts), function(cohort_name) {
  cap_data <- all_cohorts[[cohort_name]] %>%
    filter(insulin_n_fills > 0) %>%
    mutate(exceeds_cap_numeric = as.numeric(exceeds_cap))
  
  cap_model <- lm(
    exceeds_cap_numeric ~ treatment + post + treat_post + calendar_month*treatment + age_bin + der_sex,
    data = cap_data
  )
  
  tidy(cap_model) %>%
    filter(term == "treat_post") %>%
    mutate(
      cohort = cohort_name,
      ci_lower = estimate - 1.96 * std.error,
      ci_upper = estimate + 1.96 * std.error,
      pct_point_change = estimate * 100,
      sig = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
})

write_csv(cap_did_all_interacted, "monthly_cap_did_all_cohorts_interacted_months.csv")

# -----------------------------------------------------------------------------
# 9. CREATE COMPARISON TABLES
# -----------------------------------------------------------------------------

cat("\n=== CREATING COMPARISON TABLES ACROSS COHORTS ===\n\n")

# Primary outcomes comparison
primary_comparison <- combined_did_estimates_interacted %>%
  filter(outcome %in% c("insulin_oop_per_supply", "insulin_copay_per_supply", 
                        "insulin_standardized_supplies", "insulin_gap")) %>%
  select(cohort, outcome, estimate, std.error, ci_lower, ci_upper, p.value, sig) %>%
  arrange(outcome, cohort)

write_csv(primary_comparison, "monthly_primary_outcomes_comparison_interacted_months.csv")

# Create wide format for easy comparison
primary_comparison_wide <- primary_comparison %>%
  select(cohort, outcome, estimate, sig) %>%
  mutate(estimate_sig = paste0(round(estimate, 3), sig)) %>%
  select(-estimate, -sig) %>%
  pivot_wider(names_from = cohort, values_from = estimate_sig)

write_csv(primary_comparison_wide, "monthly_primary_outcomes_comparison_wide_interacted_months.csv")

# -----------------------------------------------------------------------------
# 10. REGRESSION TABLES FOR ORIGINAL COHORT
# -----------------------------------------------------------------------------

cat("\n=== CREATING REGRESSION TABLES (ORIGINAL COHORT) ===\n\n")

# PRIMARY TABLE: Per-supply costs and adherence (Panel FE with interacted months)
primary_models <- list(
  "OOP per Supply" = panel_monthly_fe_interacted$insulin_oop_per_supply$model,
  "Copay per Supply" = panel_monthly_fe_interacted$insulin_copay_per_supply$model,
  "Supplies (Adherence)" = panel_monthly_fe_interacted$insulin_standardized_supplies$model,
  "Treatment Gap" = panel_monthly_fe_interacted$insulin_gap$model
)

modelsummary(
  primary_models,
  stars = c('*' = 0.05, '**' = 0.01, '***' = 0.001),
  coef_rename = c("treat_post" = "Treatment × Post-IRA"),
  coef_omit = "calendar_month",  # Omit month dummies from table for readability
  notes = c(
    "Patient fixed effects and month-by-cohort fixed effects included.",
    "Month-by-cohort FE allows MA and Commercial to have different seasonal patterns.",
    "Standard errors clustered at patient level.",
    "Per-supply outcomes standardized to 30-day equivalents.",
    "Month dummy coefficients omitted for brevity (see separate seasonality analysis)."
  ),
  gof_map = c("nobs", "r.squared"),
  output = "monthly_panel_PRIMARY_table_interacted_months.docx"
)

# COMPARISON TABLE: Simple vs Interacted specification
spec_comparison_models <- list(
  "Simple Months" = panel_monthly_fe_simple$insulin_oop_per_supply$model,
  "Interacted Months" = panel_monthly_fe_interacted$insulin_oop_per_supply$model
)

modelsummary(
  spec_comparison_models,
  stars = c('*' = 0.05, '**' = 0.01, '***' = 0.001),
  coef_rename = c("treat_post" = "Treatment × Post-IRA"),
  coef_omit = "calendar_month",
  notes = c(
    "Dependent variable: Insulin OOP per 30-day supply",
    "Simple: Common month dummies for both groups",
    "Interacted: Separate seasonal patterns by insurance type",
    "Patient fixed effects included in both specifications."
  ),
  gof_map = c("nobs", "r.squared"),
  output = "monthly_panel_specification_comparison_table.docx"
)

# CROSS-SECTIONAL TABLE: Compare original vs restrictive cohorts
key_outcome_models <- list(
  "Original (18-64 vs 65+)" = all_did_results_interacted$original_18_64_vs_65_plus$models$insulin_oop_per_supply$model,
  "Restrictive (62-64 vs 65-67)" = all_did_results_interacted$restrictive_62_64_vs_65_67$models$insulin_oop_per_supply$model,
  "Restrictive (54-64 vs 65-75)" = all_did_results_interacted$restrictive_54_64_vs_65_75$models$insulin_oop_per_supply$model,
  "Restrictive (54-64 vs 65+)" = all_did_results_interacted$restrictive_54_64_vs_65_plus$models$insulin_oop_per_supply$model,
  "Restrictive (62-64 vs 65+)" = all_did_results_interacted$restrictive_62_64_vs_65_plus$models$insulin_oop_per_supply$model
)

modelsummary(
  key_outcome_models,
  stars = c('*' = 0.05, '**' = 0.01, '***' = 0.001),
  coef_rename = c("treat_post" = "Treatment × Post-IRA"),
  coef_omit = "calendar_month",
  notes = c(
    "Dependent variable: Insulin OOP per 30-day supply",
    "All models include age, sex, comorbidity controls, and month-by-cohort FE.",
    "Different age specifications test robustness to control group definition.",
    "Month dummy coefficients omitted for brevity."
  ),
  gof_map = c("nobs", "r.squared"),
  output = "monthly_cohort_comparison_table_interacted_months.docx"
)

# -----------------------------------------------------------------------------
# 11. COMPREHENSIVE SUMMARY
# -----------------------------------------------------------------------------

cat("\n\n=== COMPLETE ANALYSIS SUMMARY (WITH INTERACTED MONTH EFFECTS) ===\n\n")

cat("COHORTS ANALYZED:\n")
cat("1. Original: 18-64 Commercial vs 65+ Medicare Advantage\n")
cat("2. Restrictive: 62-64 vs 65-67\n")
cat("3. Restrictive: 54-64 vs 65-75\n")
cat("4. Restrictive: 54-64 vs 65+\n")
cat("5. Restrictive: 62-64 vs 65+\n\n")

cat("PRIMARY OUTCOMES (all cohorts with interacted month-by-cohort FE):\n")
primary_summary <- combined_did_estimates_interacted %>%
  filter(outcome == "insulin_oop_per_supply") %>%
  select(cohort, estimate, ci_lower, ci_upper, p.value, sig)
print(primary_summary)

cat("\n\nPANEL FE RESULTS (original cohort only - 22 months):\n")
cat("Comparison of specifications:\n")
panel_summary <- comparison_specs %>%
  filter(outcome == "insulin_oop_per_supply") %>%
  select(specification, estimate, ci_lower, ci_upper, p.value, sig)
print(panel_summary)

cat("\n\n% EXCEEDING CAP BY COHORT:\n")
cap_summary <- cap_analysis_all %>%
  filter(year == 2023, treatment == 1) %>%
  select(cohort, pct_exceeding_cap, mean_oop_per_supply)
print(cap_summary)

cat("\n\nOUTPUT FILES CREATED:\n")
cat("COMPREHENSIVE RESULTS (WITH INTERACTED MONTH-BY-COHORT FE):\n")
cat("- monthly_all_cohorts_did_estimates_interacted_months.csv (all cohorts, all outcomes)\n")
cat("- monthly_primary_outcomes_comparison_interacted_months.csv (key outcomes compared)\n")
cat("- monthly_primary_outcomes_comparison_wide_interacted_months.csv (wide format)\n")
cat("- monthly_cap_analysis_all_cohorts.csv (cap compliance all cohorts)\n")
cat("- monthly_cap_did_all_cohorts_interacted_months.csv (cap DID estimates)\n\n")

cat("SPECIFICATION COMPARISON:\n")
cat("- monthly_panel_fe_estimates_simple_months.csv (simple month dummies)\n")
cat("- monthly_panel_fe_estimates_interacted_months.csv (month-by-cohort FE)\n")
cat("- monthly_panel_specification_comparison.csv (side-by-side comparison)\n\n")

cat("SEASONALITY ANALYSIS:\n")
cat("- month_effects_by_insurance_group.csv (group-specific seasonal patterns)\n")
cat("- seasonality_by_insurance_group.png (visualization)\n\n")

cat("REGRESSION TABLES:\n")
cat("- monthly_panel_PRIMARY_table_interacted_months.docx (main results)\n")
cat("- monthly_panel_specification_comparison_table.docx (simple vs interacted)\n")
cat("- monthly_cohort_comparison_table_interacted_months.docx (cross-cohort comparison)\n\n")

cat("KEY METHODOLOGICAL UPDATE:\n")
cat("✓ Calendar month × treatment interactions (month-by-cohort FE)\n")
cat("✓ Allows MA and Commercial to have DIFFERENT seasonal patterns\n")
cat("✓ Controls for Part D donut hole, formulary timing differences\n")
cat("✓ Addresses potential parallel trends violations from seasonal heterogeneity\n")
cat("✓ F-test conducted to formally test if seasonal patterns differ\n\n")

cat("INTERPRETATION:\n")
cat("- Compare simple vs interacted specs: If similar, seasonality not driving results\n")
cat("- Examine seasonality plot: Do MA/Commercial have different patterns?\n")
cat("- Check F-test p-value: Significant = interacted spec is necessary\n")
cat("- Compare estimates across cohorts to assess robustness\n")
cat("- Panel FE (original) provides within-patient estimates\n")
cat("- All use corrected per-supply outcomes (30-day standardized)\n\n")

cat("NEXT STEPS FOR YOUR CHAIR:\n")
cat("1. Check seasonality_by_insurance_group.png - do patterns differ?\n")
cat("2. Look at F-test result - is interaction significant?\n")
cat("3. Compare simple vs interacted estimates - does it matter?\n")
cat("4. If patterns don't differ much: simple month dummies sufficient\n")
cat("5. If patterns differ significantly: interacted spec is necessary\n\n")