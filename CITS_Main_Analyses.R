# =============================================================================
# CITS ANALYSIS: ALL COHORTS + HETEROGENEITY - FIXED
# =============================================================================

library(tidyverse)
library(fixest)
library(modelsummary)
library(ggplot2)
library(patchwork)

cat("=======================================================\n")
cat("CITS ANALYSIS: MULTIPLE COHORTS + HETEROGENEITY\n")
cat("=======================================================\n\n")

# -----------------------------------------------------------------------------
# 1. LOAD AND COMBINE DATA FOR ALL COHORTS
# -----------------------------------------------------------------------------

cat("=== LOADING AND COMBINING COHORT DATA ===\n\n")

# Function to load and combine 2022 + 2023 monthly data
load_cohort_panel <- function(cohort_name) {
  file_2022 <- paste0("analytical_monthly_", cohort_name, "_2022.csv")
  file_2023 <- paste0("analytical_monthly_", cohort_name, "_2023.csv")
  
  if(!file.exists(file_2022) | !file.exists(file_2023)) {
    return(NULL)
  }
  
  data_2022 <- read_csv(file_2022, show_col_types = FALSE)
  data_2023 <- read_csv(file_2023, show_col_types = FALSE)
  
  combined <- bind_rows(data_2022, data_2023) %>%
    mutate(
      year_month = as.Date(year_month),
      treatment = as.numeric(treatment),
      post_ira = as.numeric(year_month >= as.Date("2023-01-01"))
    )
  
  return(combined)
}

# Load all cohorts
cohort_data_list <- list(
  
  # PRIMARY: Main analysis (uses the full panel)
  primary = read_csv("analytical_monthly_panel.csv", show_col_types = FALSE) %>%
    mutate(
      year_month = as.Date(year_month),
      treatment = as.numeric(treatment),
      post_ira = as.numeric(year_month >= as.Date("2023-01-01"))
    ),
  
  # CROSS-SECTIONAL (combine 2022 + 2023)
  cross_sectional = bind_rows(
    read_csv("analytical_monthly_cross_2022.csv", show_col_types = FALSE),
    read_csv("analytical_monthly_cross_2023.csv", show_col_types = FALSE)
  ) %>%
    mutate(
      year_month = as.Date(year_month),
      treatment = as.numeric(treatment),
      post_ira = as.numeric(year_month >= as.Date("2023-01-01"))
    ),
  
  # RESTRICTIVE COHORTS
  restrictive_62_64_vs_65_plus = load_cohort_panel("restrictive_62_64_vs_65_plus"),
  restrictive_54_64_vs_65_75 = load_cohort_panel("restrictive_54_64_vs_65_75"),
  restrictive_54_64_vs_65_plus = load_cohort_panel("restrictive_54_64_vs_65_plus"),
  restrictive = load_cohort_panel("restrictive")
)

# Remove NULL entries
cohort_data_list <- compact(cohort_data_list)

# Print summary
cat("Loaded cohorts:\n")
for(cohort_name in names(cohort_data_list)) {
  data <- cohort_data_list[[cohort_name]]
  cat(sprintf("  %-35s: %8d obs | %6d patients\n", 
              cohort_name, 
              nrow(data), 
              n_distinct(data$pat_id)))
}
cat("\n")

# -----------------------------------------------------------------------------
# 2. DEFINE COHORT SPECIFICATIONS
# -----------------------------------------------------------------------------

cohort_specs <- list(
  
  primary = list(
    name = "Primary (18-64 vs 65+)",
    description = "Main analysis: Commercial 18-64 vs MA 65+",
    data = cohort_data_list$primary,
    treatment_var = "treatment"
  ),
  
  cross_sectional = list(
    name = "Cross-Sectional",
    description = "Separate 2022 and 2023 cohorts combined",
    data = cohort_data_list$cross_sectional,
    treatment_var = "treatment"
  ),
  
  restrictive_62_64_vs_65_plus = list(
    name = "Restrictive (62-64 vs 65+)",
    description = "Ages 62-64 vs 65+",
    data = cohort_data_list$restrictive_62_64_vs_65_plus,
    treatment_var = "treatment"
  ),
  
  restrictive_54_64_vs_65_75 = list(
    name = "Restrictive (54-64 vs 65-75)",
    description = "Ages 54-64 vs 65-75",
    data = cohort_data_list$restrictive_54_64_vs_65_75,
    treatment_var = "treatment"
  ),
  
  restrictive_54_64_vs_65_plus = list(
    name = "Restrictive (54-64 vs 65+)",
    description = "Ages 54-64 vs 65+",
    data = cohort_data_list$restrictive_54_64_vs_65_plus,
    treatment_var = "treatment"
  ),
  
  restrictive = list(
    name = "Restrictive (62-64 vs 65-67)",
    description = "Narrow window: 62-64 vs 65-67",
    data = cohort_data_list$restrictive,
    treatment_var = "treatment"
  )
)

# Remove specs where data is NULL
cohort_specs <- cohort_specs[!sapply(cohort_specs, function(x) is.null(x$data))]

cat("Valid cohorts for analysis:\n")
walk(names(cohort_specs), ~cat(paste0("  - ", cohort_specs[[.x]]$name, "\n")))
cat("\n")

# -----------------------------------------------------------------------------
# 3. DEFINE HETEROGENEITY SUBGROUPS
# -----------------------------------------------------------------------------

cat("=== DEFINING HETEROGENEITY SUBGROUPS ===\n\n")

primary_data <- cohort_specs$primary$data

heterogeneity_specs <- list(
  
  male = list(
    name = "Male",
    filter_expr = expression(der_sex == "M")
  ),
  
  female = list(
    name = "Female", 
    filter_expr = expression(der_sex == "F")
  ),
  
  low_comorbidity = list(
    name = "Low Comorbidity",
    filter_expr = expression(n_charlson <= 1)
  ),
  
  high_comorbidity = list(
    name = "High Comorbidity",
    filter_expr = expression(n_charlson >= 2)
  ),
  
  low_dcsi = list(
    name = "Low Diabetes Severity",
    filter_expr = expression(dcsi_total <= 2)
  ),
  
  high_dcsi = list(
    name = "High Diabetes Severity",
    filter_expr = expression(dcsi_total >= 3)
  ),
  
  low_polypharm = list(
    name = "Low Polypharmacy",
    filter_expr = expression(n_unique_drugs_excl_insulin <= 5)
  ),
  
  high_polypharm = list(
    name = "High Polypharmacy",
    filter_expr = expression(n_unique_drugs_excl_insulin >= 10)
  )
)

# Check feasibility
cat("Checking feasibility of heterogeneity subgroups...\n")
heterogeneity_specs_valid <- list()

for(subgroup_name in names(heterogeneity_specs)) {
  spec <- heterogeneity_specs[[subgroup_name]]
  
  tryCatch({
    test_data <- primary_data %>%
      filter(eval(spec$filter_expr))
    
    if(nrow(test_data) > 0) {
      n_treat <- sum(test_data$treatment == 1, na.rm = TRUE)
      n_control <- sum(test_data$treatment == 0, na.rm = TRUE)
      
      if(n_treat > 100 && n_control > 100) {
        cat("  ✓", subgroup_name, ": Treat:", n_treat, "| Control:", n_control, "\n")
        heterogeneity_specs_valid[[subgroup_name]] <- spec
      } else {
        cat("  ✗", subgroup_name, ": Insufficient sample\n")
      }
    }
  }, error = function(e) {
    cat("  ✗", subgroup_name, ": Error -", e$message, "\n")
  })
}

heterogeneity_specs <- heterogeneity_specs_valid
cat("\n")

# -----------------------------------------------------------------------------
# 4. CITS MODEL FUNCTION - COMPLETELY REWRITTEN
# -----------------------------------------------------------------------------

run_cits_model_robust <- function(data, outcome, controls, treatment_var = "treatment") {
  
  # Check if outcome exists
  if(!outcome %in% names(data)) {
    return(NULL)
  }
  
  # Prepare CITS variables
  data <- data %>%
    mutate(
      time = as.numeric(year_month - min(year_month, na.rm = TRUE)) / 30,
      post = as.numeric(year_month >= as.Date("2023-01-01")),
      time_since_ira = pmax(0, as.numeric(year_month - as.Date("2023-01-01")) / 30),
      calendar_month = factor(month(year_month), levels = 1:12, labels = month.abb),
      treat_time = !!sym(treatment_var) * time,
      treat_post = !!sym(treatment_var) * post,
      treat_time_since = !!sym(treatment_var) * time_since_ira
    )
  
  # Remove missing/infinite values for outcome
  data_clean <- data %>%
    filter(!is.na(!!sym(outcome)), 
           !is.infinite(!!sym(outcome)))
  
  # Check if we have enough data
  if(nrow(data_clean) < 100) {
    return(NULL)
  }
  
  # Check variance
  outcome_variance <- var(data_clean[[outcome]], na.rm = TRUE)
  if(is.na(outcome_variance) || outcome_variance == 0) {
    return(NULL)
  }
  
  # Build formula
  formula_str <- paste0(
    outcome, " ~ ",
    "time + ",
    treatment_var, " + ",
    "post + ",
    "time_since_ira + ",
    "treat_time + ",
    "treat_post + ",
    "treat_time_since + ",
    "calendar_month + ",
    paste(controls, collapse = " + "),
    " | pat_id"
  )
  
  # Run model
  model <- tryCatch({
    feols(as.formula(formula_str), data = data_clean, cluster = ~pat_id)
  }, error = function(e) {
    return(NULL)
  })
  
  if(is.null(model)) return(NULL)
  
  # Extract coefficients
  coef_values <- coef(model)
  se_values <- se(model)
  t_stats <- coef_values / se_values
  p_values <- 2 * pt(abs(t_stats), df = model$nobs - length(coef_values), lower.tail = FALSE)
  
  coef_summary <- tibble(
    term = names(coef_values),
    estimate = coef_values,
    std.error = se_values,
    statistic = t_stats,
    p.value = p_values
  ) %>%
    filter(term %in% c("time", treatment_var, "post", "time_since_ira",
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

# -----------------------------------------------------------------------------
# 5. RUN CITS FOR ALL COHORTS
# -----------------------------------------------------------------------------

cat("=== RUNNING CITS FOR ALL COHORTS ===\n\n")

primary_outcomes <- c("insulin_oop_per_supply", "insulin_copay_per_supply")
adherence_outcomes <- c("insulin_standardized_supplies", "insulin_gap")
all_outcomes <- c(primary_outcomes, adherence_outcomes)
controls <- c("n_charlson", "dcsi_total", "n_unique_drugs_excl_insulin")

all_cohort_results <- list()

for(cohort_name in names(cohort_specs)) {
  
  cat("\n--- COHORT:", cohort_specs[[cohort_name]]$name, "---\n")
  
  cohort_data <- cohort_specs[[cohort_name]]$data
  
  cat("  N observations:", nrow(cohort_data), "\n")
  cat("  N patients:", n_distinct(cohort_data$pat_id), "\n")
  
  # Check available outcomes
  available_outcomes <- all_outcomes[all_outcomes %in% names(cohort_data)]
  cat("  Available outcomes:", paste(available_outcomes, collapse = ", "), "\n")
  
  cohort_results <- list()
  
  for(outcome in available_outcomes) {
    cat("    Running:", outcome, "...")
    
    result <- run_cits_model_robust(
      data = cohort_data,
      outcome = outcome,
      controls = controls,
      treatment_var = cohort_specs[[cohort_name]]$treatment_var
    )
    
    if(!is.null(result)) {
      result$coefs <- result$coefs %>%
        mutate(cohort = cohort_name)
      cohort_results[[outcome]] <- result
      cat(" ✓\n")
    } else {
      cat(" ✗ Failed\n")
    }
  }
  
  all_cohort_results[[cohort_name]] <- cohort_results
}

# Combine coefficients
cohort_coefs <- map_df(all_cohort_results, function(cohort) {
  map_df(cohort, ~.x$coefs)
})

write_csv(cohort_coefs, "cits_all_cohorts_coefficients.csv")
cat("\n✓ Saved: cits_all_cohorts_coefficients.csv\n\n")

# -----------------------------------------------------------------------------
# 6. RUN CITS FOR HETEROGENEITY SUBGROUPS
# -----------------------------------------------------------------------------

cat("=== RUNNING HETEROGENEITY ANALYSES ===\n\n")

heterogeneity_results <- list()

for(subgroup_name in names(heterogeneity_specs)) {
  
  cat("\n--- SUBGROUP:", heterogeneity_specs[[subgroup_name]]$name, "---\n")
  
  subgroup_data <- primary_data %>%
    filter(eval(heterogeneity_specs[[subgroup_name]]$filter_expr))
  
  cat("  N observations:", nrow(subgroup_data), "\n")
  cat("  N patients:", n_distinct(subgroup_data$pat_id), "\n")
  
  subgroup_results <- list()
  
  for(outcome in all_outcomes) {
    cat("    Running:", outcome, "...")
    
    result <- run_cits_model_robust(
      data = subgroup_data,
      outcome = outcome,
      controls = controls,
      treatment_var = "treatment"
    )
    
    if(!is.null(result)) {
      result$coefs <- result$coefs %>%
        mutate(subgroup = subgroup_name)
      subgroup_results[[outcome]] <- result
      cat(" ✓\n")
    } else {
      cat(" ✗ Failed\n")
    }
  }
  
  heterogeneity_results[[subgroup_name]] <- subgroup_results
}

# Combine coefficients
heterogeneity_coefs <- map_df(heterogeneity_results, function(subgroup) {
  map_df(subgroup, ~.x$coefs)
})

write_csv(heterogeneity_coefs, "cits_heterogeneity_coefficients.csv")
cat("\n✓ Saved: cits_heterogeneity_coefficients.csv\n\n")

# -----------------------------------------------------------------------------
# 7. CREATE SUMMARY TABLES (FIXED)
# -----------------------------------------------------------------------------

cat("=== CREATING SUMMARY TABLES ===\n\n")

# Cohort comparison
oop_comparison <- cohort_coefs %>%
  filter(outcome == "insulin_oop_per_supply",
         term %in% c("treat_time", "treat_post", "treat_time_since")) %>%
  select(cohort, term, estimate, std.error, p.value, sig) %>%
  pivot_wider(
    names_from = term,
    values_from = c(estimate, std.error, p.value, sig)
  )

write_csv(oop_comparison, "cits_cohort_comparison_table.csv")

# Heterogeneity comparison
hetero_comparison <- heterogeneity_coefs %>%
  filter(outcome == "insulin_oop_per_supply", term == "treat_post") %>%
  select(subgroup, estimate, std.error, p.value, sig, ci_lower, ci_upper) %>%
  arrange(estimate)

write_csv(hetero_comparison, "cits_heterogeneity_comparison_table.csv")

cat("✓ Summary tables created\n\n")

# Print key results
cat("KEY RESULTS - OOP LEVEL SHIFTS:\n\n")
cat("BY COHORT:\n")
print(oop_comparison %>% select(cohort, estimate_treat_post, sig_treat_post), n = Inf)
cat("\n")

cat("BY SUBGROUP:\n")
print(hetero_comparison, n = Inf)
cat("\n")

# -----------------------------------------------------------------------------
# 8. CREATE PLOTS
# -----------------------------------------------------------------------------

cat("=== CREATING VISUALIZATIONS ===\n\n")

# Cohort labels
cohort_labels <- c(
  "primary" = "Primary\n18-64 vs 65+",
  "cross_sectional" = "Cross-Sect\n2022+2023",
  "restrictive_62_64_vs_65_plus" = "Restrict\n62-64 vs 65+",
  "restrictive_54_64_vs_65_75" = "Restrict\n54-64 vs 65-75",
  "restrictive_54_64_vs_65_plus" = "Restrict\n54-64 vs 65+",
  "restrictive" = "Restrict\n62-64 vs 65-67"
)

# PLOT 1: Cohort comparison
plot_cohorts <- cohort_coefs %>%
  filter(term == "treat_post") %>%
  mutate(
    outcome_label = case_when(
      outcome == "insulin_oop_per_supply" ~ "OOP",
      outcome == "insulin_copay_per_supply" ~ "Copay",
      outcome == "insulin_standardized_supplies" ~ "Fills",
      outcome == "insulin_gap" ~ "Gaps"
    ),
    cohort_label = cohort_labels[cohort]
  ) %>%
  filter(!is.na(cohort_label)) %>%
  ggplot(aes(x = cohort_label, y = estimate, color = outcome_label, group = outcome_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                width = 0.25, position = position_dodge(width = 0.6)) +
  geom_text(aes(label = sig, y = ci_upper), 
            vjust = -0.5, size = 3.5,
            position = position_dodge(width = 0.6)) +
  scale_color_manual(
    values = c("OOP" = "#2c7bb6", "Copay" = "#d7191c",
               "Fills" = "#abdda4", "Gaps" = "#fdae61")
  ) +
  labs(
    title = "CITS: Treatment Effects Across Cohorts",
    subtitle = "Immediate IRA effect (treat_post)",
    x = NULL,
    y = "Effect Size",
    color = "Outcome"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave("cits_cohort_comparison.png", plot_cohorts, width = 12, height = 7, dpi = 300)

# PLOT 2: Heterogeneity
plot_hetero <- heterogeneity_coefs %>%
  filter(term == "treat_post", 
         outcome %in% c("insulin_oop_per_supply", "insulin_standardized_supplies")) %>%
  mutate(
    outcome_label = ifelse(outcome == "insulin_oop_per_supply", "OOP ($)", "Fills"),
    subgroup_label = factor(subgroup, levels = names(heterogeneity_specs),
                            labels = sapply(heterogeneity_specs, function(x) x$name))
  ) %>%
  ggplot(aes(x = subgroup_label, y = estimate, color = outcome_label, group = outcome_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                width = 0.25, position = position_dodge(width = 0.6)) +
  geom_text(aes(label = sig, y = ci_upper),
            vjust = -0.5, size = 3.5,
            position = position_dodge(width = 0.6)) +
  scale_color_manual(values = c("OOP ($)" = "#2c7bb6", "Fills" = "#abdda4")) +
  labs(
    title = "CITS: Heterogeneous Treatment Effects",
    subtitle = "Primary cohort by patient characteristics",
    x = NULL,
    y = "Effect Size",
    color = "Outcome"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave("cits_heterogeneity_comparison.png", plot_hetero, width = 12, height = 7, dpi = 300)

cat("✓ Plots saved\n\n")

# -----------------------------------------------------------------------------
# 9. FINAL SUMMARY
# -----------------------------------------------------------------------------

cat("=======================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=======================================================\n\n")

cat("COHORTS ANALYZED:\n")
walk(names(cohort_specs), ~cat(paste0("  - ", cohort_specs[[.x]]$name, "\n")))
cat("\n")

cat("HETEROGENEITY SUBGROUPS:\n")
walk(names(heterogeneity_specs), ~cat(paste0("  - ", heterogeneity_specs[[.x]]$name, "\n")))
cat("\n")

cat("FILES CREATED:\n")
cat("  - cits_all_cohorts_coefficients.csv\n")
cat("  - cits_heterogeneity_coefficients.csv\n")
cat("  - cits_cohort_comparison_table.csv\n")
cat("  - cits_heterogeneity_comparison_table.csv\n")
cat("  - cits_cohort_comparison.png\n")
cat("  - cits_heterogeneity_comparison.png\n\n")

# Display key findings
cat("KEY FINDINGS SUMMARY:\n")
cat("---------------------\n\n")

cat("1. PRIMARY COHORT (OOP per supply):\n")
primary_oop <- cohort_coefs %>%
  filter(cohort == "primary", outcome == "insulin_oop_per_supply",
         term %in% c("treat_time", "treat_post", "treat_time_since"))

if(nrow(primary_oop) > 0) {
  for(i in 1:nrow(primary_oop)) {
    cat(sprintf("   %s: %.3f (SE=%.3f) %s\n", 
                primary_oop$term[i], 
                primary_oop$estimate[i], 
                primary_oop$std.error[i],
                primary_oop$sig[i]))
  }
}
cat("\n")

cat("2. HETEROGENEITY (treat_post for OOP):\n")
hetero_summary <- heterogeneity_coefs %>%
  filter(outcome == "insulin_oop_per_supply", term == "treat_post") %>%
  arrange(estimate) %>%
  select(subgroup, estimate, std.error, sig)

if(nrow(hetero_summary) > 0) {
  for(i in 1:nrow(hetero_summary)) {
    cat(sprintf("   %-25s: %7.3f (SE=%.3f) %s\n",
                hetero_summary$subgroup[i],
                hetero_summary$estimate[i],
                hetero_summary$std.error[i],
                hetero_summary$sig[i]))
  }
}
cat("\n")

cat("DONE!\n")