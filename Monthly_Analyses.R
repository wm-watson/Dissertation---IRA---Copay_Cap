# =============================================================================
# COMPLETE MONTHLY ANALYSIS - ALL COHORTS (UPDATED WITH MONTH DUMMIES)
# =============================================================================
# Analysis for:
# - Original: 18-64 Commercial vs 65+ Medicare Advantage
# - Restrictive 1: 62-64 vs 65-67
# - Restrictive 2: 54-64 vs 65-75
# - Restrictive 3: 54-64 vs 65+
# - Restrictive 4: 62-64 vs 65+
#
# UPDATED: Added calendar month dummies to capture seasonality/cyclicality
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
# 2. DATA CLEANING FUNCTION (UPDATED WITH MONTH DUMMIES)
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
# 4. RUN DID MODELS FOR ALL COHORTS (CROSS-SECTIONAL WITH MONTH DUMMIES)
# -----------------------------------------------------------------------------

cat("\n=== RUNNING DID MODELS FOR ALL COHORTS (WITH MONTH DUMMIES) ===\n\n")

run_monthly_did <- function(data, outcome, controls) {
  # Add calendar_month to controls for cross-sectional models
  formula_str <- paste0(outcome, " ~ treatment + post + treat_post + calendar_month + ", 
                        paste(controls, collapse = " + "))
  
  model <- lm(as.formula(formula_str), data = data)
  
  did_estimate <- tidy(model) %>%
    filter(term == "treat_post") %>%
    mutate(outcome = outcome)
  
  return(list(model = model, did = did_estimate))
}

# Run DID for all cohorts and all outcomes
all_did_results <- map(names(all_cohorts), function(cohort_name) {
  cat("Processing cohort:", cohort_name, "\n")
  
  cohort_data <- all_cohorts[[cohort_name]]
  
  # Add age_bin and der_sex to controls for cross-sectional models
  cross_sectional_controls <- c("age_bin", "der_sex", binned_controls)
  
  # Run DID for all outcomes
  cohort_results <- map(all_outcomes, ~run_monthly_did(cohort_data, .x, cross_sectional_controls))
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

names(all_did_results) <- names(all_cohorts)

# Combine all estimates
combined_did_estimates <- map_df(all_did_results, ~.x$estimates)

write_csv(combined_did_estimates, "monthly_all_cohorts_did_estimates_with_month_dummies.csv")

# -----------------------------------------------------------------------------
# 5. PANEL FE MODELS WITH MONTH DUMMIES (ORIGINAL 22-MONTH PANEL ONLY)
# -----------------------------------------------------------------------------

cat("\n=== RUNNING PANEL FE MODELS WITH MONTH DUMMIES (22 MONTHS - ORIGINAL COHORT) ===\n\n")

run_panel_fe_with_months <- function(data, outcome, controls) {
  # Two-way FE with calendar month dummies
  # Patient FE absorbs time-invariant characteristics (age_bin, der_sex)
  # Year-month FE would absorb calendar_month, so we use calendar_month instead
  
  formula_str <- paste0(outcome, " ~ treat_post + calendar_month + ", 
                        paste(controls, collapse = " + "),
                        " | pat_id")
  
  model <- feols(as.formula(formula_str), data = data, cluster = ~pat_id)
  
  did_estimate <- tidy(model) %>%
    filter(term == "treat_post") %>%
    mutate(outcome = outcome)
  
  return(list(model = model, did = did_estimate))
}

# Run for all outcomes
panel_monthly_fe <- map(all_outcomes, ~run_panel_fe_with_months(monthly_panel, .x, binned_controls))
names(panel_monthly_fe) <- all_outcomes

panel_monthly_estimates <- map_df(panel_monthly_fe, ~.x$did) %>%
  mutate(
    cohort = "panel_22_months",
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

write_csv(panel_monthly_estimates, "monthly_panel_fe_estimates_with_month_dummies.csv")

# -----------------------------------------------------------------------------
# 6. EXTRACT AND SAVE MONTH COEFFICIENTS (SEASONALITY ANALYSIS)
# -----------------------------------------------------------------------------

cat("\n=== EXTRACTING MONTH COEFFICIENTS (SEASONALITY PATTERNS) ===\n\n")

# Extract month coefficients from primary outcome models
extract_month_coefs <- function(model) {
  tidy(model) %>%
    filter(str_detect(term, "calendar_month")) %>%
    mutate(
      month = str_remove(term, "calendar_month"),
      month = factor(month, levels = month.abb)
    )
}

# Panel FE: Extract month effects for primary outcomes
month_effects_panel <- map_df(c("insulin_oop_per_supply", "insulin_copay_per_supply",
                                "insulin_standardized_supplies", "insulin_n_fills"), 
                              function(outcome) {
                                extract_month_coefs(panel_monthly_fe[[outcome]]$model) %>%
                                  mutate(outcome = outcome, model = "panel_fe")
                              })

write_csv(month_effects_panel, "month_seasonality_effects_panel.csv")

# Cross-sectional: Extract month effects for original cohort
month_effects_cross <- map_df(c("insulin_oop_per_supply", "insulin_copay_per_supply",
                                "insulin_standardized_supplies", "insulin_n_fills"), 
                              function(outcome) {
                                extract_month_coefs(all_did_results$original_18_64_vs_65_plus$models[[outcome]]$model) %>%
                                  mutate(outcome = outcome, model = "cross_sectional")
                              })

write_csv(month_effects_cross, "month_seasonality_effects_cross_sectional.csv")

# Combine for comparison
month_effects_combined <- bind_rows(month_effects_panel, month_effects_cross)

# Plot seasonality patterns
library(ggplot2)

plot_seasonality <- month_effects_combined %>%
  filter(outcome == "insulin_oop_per_supply") %>%
  ggplot(aes(x = month, y = estimate, color = model, group = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = estimate - 1.96 * std.error, 
                    ymax = estimate + 1.96 * std.error), width = 0.2) +
  scale_color_manual(values = c("panel_fe" = "#2c7bb6", "cross_sectional" = "#d7191c"),
                     labels = c("Panel FE", "Cross-sectional")) +
  labs(
    title = "Seasonality in Insulin Out-of-Pocket Costs",
    subtitle = "Calendar month effects (relative to January)",
    x = "Month",
    y = "Coefficient Estimate ($)",
    color = "Model",
    caption = "Note: 95% confidence intervals shown. January is the reference month."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("seasonality_insulin_oop.png", plot_seasonality, width = 12, height = 8, dpi = 300)

cat("\nSeasonality plot saved: seasonality_insulin_oop.png\n")

# -----------------------------------------------------------------------------
# 7. CAP ANALYSIS FOR ALL COHORTS
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

# DID for cap exceeding (all cohorts) - with month dummies
cap_did_all <- map_df(names(all_cohorts), function(cohort_name) {
  cap_data <- all_cohorts[[cohort_name]] %>%
    filter(insulin_n_fills > 0) %>%
    mutate(exceeds_cap_numeric = as.numeric(exceeds_cap))
  
  cap_model <- lm(
    exceeds_cap_numeric ~ treatment + post + treat_post + calendar_month + age_bin + der_sex,
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

write_csv(cap_did_all, "monthly_cap_did_all_cohorts_with_month_dummies.csv")

# -----------------------------------------------------------------------------
# 8. CREATE COMPARISON TABLES
# -----------------------------------------------------------------------------

cat("\n=== CREATING COMPARISON TABLES ACROSS COHORTS ===\n\n")

# Primary outcomes comparison
primary_comparison <- combined_did_estimates %>%
  filter(outcome %in% c("insulin_oop_per_supply", "insulin_copay_per_supply", 
                        "insulin_standardized_supplies", "insulin_gap")) %>%
  select(cohort, outcome, estimate, std.error, ci_lower, ci_upper, p.value, sig) %>%
  arrange(outcome, cohort)

write_csv(primary_comparison, "monthly_primary_outcomes_comparison_with_month_dummies.csv")

# Create wide format for easy comparison
primary_comparison_wide <- primary_comparison %>%
  select(cohort, outcome, estimate, sig) %>%
  mutate(estimate_sig = paste0(round(estimate, 3), sig)) %>%
  select(-estimate, -sig) %>%
  pivot_wider(names_from = cohort, values_from = estimate_sig)

write_csv(primary_comparison_wide, "monthly_primary_outcomes_comparison_wide_with_month_dummies.csv")

# -----------------------------------------------------------------------------
# 9. REGRESSION TABLES FOR ORIGINAL COHORT
# -----------------------------------------------------------------------------

cat("\n=== CREATING REGRESSION TABLES (ORIGINAL COHORT) ===\n\n")

# PRIMARY TABLE: Per-supply costs and adherence (Panel FE with month dummies)
primary_models <- list(
  "OOP per Supply" = panel_monthly_fe$insulin_oop_per_supply$model,
  "Copay per Supply" = panel_monthly_fe$insulin_copay_per_supply$model,
  "Supplies (Adherence)" = panel_monthly_fe$insulin_standardized_supplies$model,
  "Treatment Gap" = panel_monthly_fe$insulin_gap$model
)

modelsummary(
  primary_models,
  stars = c('*' = 0.05, '**' = 0.01, '***' = 0.001),
  coef_rename = c("treat_post" = "Treatment Ã Post-IRA"),
  coef_omit = "calendar_month",  # Omit month dummies from table for readability
  notes = c(
    "Patient fixed effects and calendar month dummies included.",
    "Standard errors clustered at patient level.",
    "Per-supply outcomes standardized to 30-day equivalents.",
    "Month dummy coefficients omitted for brevity (see separate seasonality analysis)."
  ),
  gof_map = c("nobs", "r.squared"),
  output = "monthly_panel_PRIMARY_table_with_month_dummies.docx"
)

# CROSS-SECTIONAL TABLE: Compare original vs restrictive cohorts
# Extract key models from each cohort for main outcome
key_outcome_models <- list(
  "Original (18-64 vs 65+)" = all_did_results$original_18_64_vs_65_plus$models$insulin_oop_per_supply$model,
  "Restrictive (62-64 vs 65-67)" = all_did_results$restrictive_62_64_vs_65_67$models$insulin_oop_per_supply$model,
  "Restrictive (54-64 vs 65-75)" = all_did_results$restrictive_54_64_vs_65_75$models$insulin_oop_per_supply$model,
  "Restrictive (54-64 vs 65+)" = all_did_results$restrictive_54_64_vs_65_plus$models$insulin_oop_per_supply$model,
  "Restrictive (62-64 vs 65+)" = all_did_results$restrictive_62_64_vs_65_plus$models$insulin_oop_per_supply$model
)

modelsummary(
  key_outcome_models,
  stars = c('*' = 0.05, '**' = 0.01, '***' = 0.001),
  coef_rename = c("treat_post" = "Treatment Ã Post-IRA"),
  coef_omit = "calendar_month",  # Omit month dummies from table
  notes = c(
    "Dependent variable: Insulin OOP per 30-day supply",
    "All models include age, sex, comorbidity controls, and calendar month dummies.",
    "Different age specifications test robustness to control group definition.",
    "Month dummy coefficients omitted for brevity."
  ),
  gof_map = c("nobs", "r.squared"),
  output = "monthly_cohort_comparison_table_with_month_dummies.docx"
)

# -----------------------------------------------------------------------------
# 10. COMPREHENSIVE SUMMARY
# -----------------------------------------------------------------------------

cat("\n\n=== COMPLETE ANALYSIS SUMMARY (WITH MONTH DUMMIES) ===\n\n")

cat("COHORTS ANALYZED:\n")
cat("1. Original: 18-64 Commercial vs 65+ Medicare Advantage\n")
cat("2. Restrictive: 62-64 vs 65-67\n")
cat("3. Restrictive: 54-64 vs 65-75\n")
cat("4. Restrictive: 54-64 vs 65+\n")
cat("5. Restrictive: 62-64 vs 65+\n\n")

cat("PRIMARY OUTCOMES (all cohorts with month dummies):\n")
primary_summary <- combined_did_estimates %>%
  filter(outcome == "insulin_oop_per_supply") %>%
  select(cohort, estimate, ci_lower, ci_upper, p.value, sig)
print(primary_summary)

cat("\n\nPANEL FE RESULTS (original cohort only - 22 months with month dummies):\n")
panel_summary <- panel_monthly_estimates %>%
  filter(outcome %in% c("insulin_oop_per_supply", "insulin_copay_per_supply",
                        "insulin_standardized_supplies", "insulin_gap")) %>%
  select(outcome, estimate, ci_lower, ci_upper, p.value, sig)
print(panel_summary)

cat("\n\n% EXCEEDING CAP BY COHORT:\n")
cap_summary <- cap_analysis_all %>%
  filter(year == 2023, treatment == 1) %>%
  select(cohort, pct_exceeding_cap, mean_oop_per_supply)
print(cap_summary)

cat("\n\nOUTPUT FILES CREATED:\n")
cat("COMPREHENSIVE RESULTS (WITH MONTH DUMMIES):\n")
cat("- monthly_all_cohorts_did_estimates_with_month_dummies.csv (all cohorts, all outcomes)\n")
cat("- monthly_primary_outcomes_comparison_with_month_dummies.csv (key outcomes compared)\n")
cat("- monthly_primary_outcomes_comparison_wide_with_month_dummies.csv (wide format)\n")
cat("- monthly_cap_analysis_all_cohorts.csv (cap compliance all cohorts)\n")
cat("- monthly_cap_did_all_cohorts_with_month_dummies.csv (cap DID estimates)\n\n")

cat("SEASONALITY ANALYSIS (NEW):\n")
cat("- month_seasonality_effects_panel.csv (month effects from panel FE)\n")
cat("- month_seasonality_effects_cross_sectional.csv (month effects from cross-sectional)\n")
cat("- seasonality_insulin_oop.png (visualization of seasonal patterns)\n\n")

cat("ORIGINAL COHORT SPECIFIC:\n")
cat("- monthly_panel_fe_estimates_with_month_dummies.csv (22-month panel FE)\n")
cat("- monthly_panel_PRIMARY_table_with_month_dummies.docx (formatted regression table)\n")
cat("- monthly_cohort_comparison_table_with_month_dummies.docx (cross-cohort comparison)\n\n")

cat("KEY METHODOLOGICAL UPDATE:\n")
cat("â Calendar month dummies included to capture seasonality/cyclicality\n")
cat("â Controls for patterns like January deductible resets, summer utilization drops\n")
cat("â Two-way fixed effects: patient FE + month dummies (not year-month FE)\n")
cat("â Seasonality coefficients extracted and visualized separately\n\n")

cat("INTERPRETATION:\n")
cat("- Compare estimates across cohorts to assess robustness\n")
cat("- Examine month coefficients to understand seasonal patterns\n")
cat("- Panel FE (original) provides within-patient estimates\n")
cat("- All use corrected per-supply outcomes (30-day standardized)\n")
cat("- Month dummies control for cyclical healthcare spending patterns\n\n")