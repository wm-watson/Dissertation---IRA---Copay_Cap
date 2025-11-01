# =============================================================================
# EVENT STUDY ANALYSIS - ALL COHORTS (UPDATED WITH MONTH DUMMIES)
# =============================================================================
# Runs event studies for:
# - Original: 18-64 Commercial vs 65+ Medicare Advantage
# - Restrictive 1: 62-64 vs 65-67
# - Restrictive 2: 54-64 vs 65-75
# - Restrictive 3: 54-64 vs 65+
# - Restrictive 4: 62-64 vs 65+
#
# UPDATED: Added calendar month dummies to control for seasonality/cyclicality
# =============================================================================

library(tidyverse)
library(fixest)
library(ggplot2)
library(patchwork)
library(lubridate)

# -----------------------------------------------------------------------------
# 1. LOAD ALL PANEL DATASETS
# -----------------------------------------------------------------------------

cat("=== LOADING PANEL DATASETS ===\n\n")

# ORIGINAL PANEL (18-64 vs 65+)
monthly_panel_original <- read_csv("analytical_monthly_panel.csv")

# Note: Restrictive cohorts don't have pre-built panel datasets
# We need to create panels from the cross-sectional data
# Load cross-sectional restrictive cohorts
restrictive_62_64_vs_65_67_2022 <- read_csv("analytical_monthly_restrictive_2022.csv")
restrictive_62_64_vs_65_67_2023 <- read_csv("analytical_monthly_restrictive_2023.csv")

restrictive_54_64_vs_65_75_2022 <- read_csv("analytical_monthly_restrictive_54_64_vs_65_75_2022.csv")
restrictive_54_64_vs_65_75_2023 <- read_csv("analytical_monthly_restrictive_54_64_vs_65_75_2023.csv")

restrictive_54_64_vs_65_plus_2022 <- read_csv("analytical_monthly_restrictive_54_64_vs_65_plus_2022.csv")
restrictive_54_64_vs_65_plus_2023 <- read_csv("analytical_monthly_restrictive_54_64_vs_65_plus_2023.csv")

restrictive_62_64_vs_65_plus_2022 <- read_csv("analytical_monthly_restrictive_62_64_vs_65_plus_2022.csv")
restrictive_62_64_vs_65_plus_2023 <- read_csv("analytical_monthly_restrictive_62_64_vs_65_plus_2023.csv")

# Create pseudo-panels for restrictive cohorts (patients present in both years)
create_pseudo_panel <- function(data_2022, data_2023) {
  # Find patients present in both years
  patients_2022 <- unique(data_2022$pat_id)
  patients_2023 <- unique(data_2023$pat_id)
  panel_patients <- intersect(patients_2022, patients_2023)
  
  # Filter to panel patients and combine
  panel <- bind_rows(
    data_2022 %>% filter(pat_id %in% panel_patients),
    data_2023 %>% filter(pat_id %in% panel_patients)
  ) %>%
    arrange(pat_id, year_month)
  
  return(panel)
}

# Create panels for each restrictive cohort
monthly_panel_62_64_vs_65_67 <- create_pseudo_panel(
  restrictive_62_64_vs_65_67_2022, restrictive_62_64_vs_65_67_2023
)

monthly_panel_54_64_vs_65_75 <- create_pseudo_panel(
  restrictive_54_64_vs_65_75_2022, restrictive_54_64_vs_65_75_2023
)

monthly_panel_54_64_vs_65_plus <- create_pseudo_panel(
  restrictive_54_64_vs_65_plus_2022, restrictive_54_64_vs_65_plus_2023
)

monthly_panel_62_64_vs_65_plus <- create_pseudo_panel(
  restrictive_62_64_vs_65_plus_2022, restrictive_62_64_vs_65_plus_2023
)

cat("Panel datasets created:\n")
cat("- Original: ", n_distinct(monthly_panel_original$pat_id), "patients\n")
cat("- 62-64 vs 65-67: ", n_distinct(monthly_panel_62_64_vs_65_67$pat_id), "patients\n")
cat("- 54-64 vs 65-75: ", n_distinct(monthly_panel_54_64_vs_65_75$pat_id), "patients\n")
cat("- 54-64 vs 65+: ", n_distinct(monthly_panel_54_64_vs_65_plus$pat_id), "patients\n")
cat("- 62-64 vs 65+: ", n_distinct(monthly_panel_62_64_vs_65_plus$pat_id), "patients\n\n")

# -----------------------------------------------------------------------------
# 2. PREPARE DATA WITH RELATIVE TIME AND MONTH DUMMIES
# -----------------------------------------------------------------------------

prepare_panel_data <- function(data) {
  data %>%
    mutate(
      year_month = as.Date(year_month),
      relative_month = interval(as.Date("2023-01-01"), year_month) %/% months(1),
      treatment = as.numeric(treatment),
      # CRITICAL: Add calendar month as factor for seasonality controls
      calendar_month = factor(month(year_month), levels = 1:12, labels = month.abb)
    )
}

all_panels <- list(
  "original_18_64_vs_65_plus" = prepare_panel_data(monthly_panel_original),
  "restrictive_62_64_vs_65_67" = prepare_panel_data(monthly_panel_62_64_vs_65_67),
  "restrictive_54_64_vs_65_75" = prepare_panel_data(monthly_panel_54_64_vs_65_75),
  "restrictive_54_64_vs_65_plus" = prepare_panel_data(monthly_panel_54_64_vs_65_plus),
  "restrictive_62_64_vs_65_plus" = prepare_panel_data(monthly_panel_62_64_vs_65_plus)
)

# -----------------------------------------------------------------------------
# 3. DEFINE EVENT STUDY FUNCTIONS (UPDATED WITH MONTH DUMMIES)
# -----------------------------------------------------------------------------

run_event_study <- function(data, outcome, omit_month = -1) {
  
  # Event study with patient FE + calendar month dummies
  # This controls for:
  # 1. Patient-specific time-invariant characteristics (via pat_id FE)
  # 2. Seasonal/cyclical patterns (via calendar_month dummies)
  # 3. Time-varying confounders (via continuous controls)
  
  formula_str <- paste0(
    outcome, " ~ ",
    paste0("i(relative_month, treatment, ref = ", omit_month, ")"),
    " + calendar_month",  # ADDED: Month dummies for seasonality
    " + n_charlson + dcsi_total + n_unique_drugs_excl_insulin",
    " | pat_id"  # Patient FE (absorbs age_bin, der_sex)
  )
  
  model <- feols(
    as.formula(formula_str),
    data = data,
    cluster = ~pat_id
  )
  
  # Extract coefficients for event time
  coefs <- tidy(model) %>%
    filter(str_detect(term, "relative_month")) %>%
    mutate(
      relative_month = case_when(
        str_detect(term, "::-") ~ as.numeric(str_extract(term, "(?<=::)-\\d+")),
        str_detect(term, "::") ~ as.numeric(str_extract(term, "(?<=::)\\d+")),
        TRUE ~ NA_real_
      )
    ) %>%
    filter(!is.na(relative_month)) %>%
    select(relative_month, estimate, std.error, p.value) %>%
    mutate(
      ci_lower = estimate - 1.96 * std.error,
      ci_upper = estimate + 1.96 * std.error
    )
  
  return(list(model = model, coefs = coefs))
}

# Extract month coefficients (seasonality pattern)
extract_month_effects <- function(model) {
  tidy(model) %>%
    filter(str_detect(term, "calendar_month")) %>%
    mutate(
      month = str_remove(term, "calendar_month"),
      month = factor(month, levels = month.abb)
    ) %>%
    select(month, estimate, std.error, p.value) %>%
    mutate(
      ci_lower = estimate - 1.96 * std.error,
      ci_upper = estimate + 1.96 * std.error
    )
}

create_event_plot <- function(coefs, title, ylab, cohort_label = NULL) {
  subtitle <- if(!is.null(cohort_label)) {
    paste0("Cohort: ", cohort_label)
  } else {
    NULL
  }
  
  ggplot(coefs, aes(x = relative_month, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "solid", color = "red", linewidth = 0.8) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
    geom_line() +
    scale_x_continuous(
      breaks = seq(-12, 9, 3),
      labels = c("Jan\n2022", "Apr\n2022", "Jul\n2022", "Oct\n2022", 
                 "Jan\n2023", "Apr\n2023", "Jul\n2023", "Oct\n2023")
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Month",
      y = ylab,
      caption = "Reference: December 2022 (month -1). 95% CI shown. Red line = IRA implementation.\nMonth dummies included to control for seasonality."
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      panel.grid.minor = element_blank(),
      plot.caption = element_text(size = 8, hjust = 0)
    )
}

create_seasonality_plot <- function(month_effects, title) {
  ggplot(month_effects, aes(x = month, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(size = 3, color = "#2c7bb6") +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, linewidth = 1) +
    geom_line(aes(group = 1), linewidth = 1, color = "#2c7bb6") +
    labs(
      title = title,
      subtitle = "Calendar month effects (January is reference)",
      x = "Month",
      y = "Effect ($)",
      caption = "95% confidence intervals shown. Captures seasonal patterns in healthcare spending."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.caption = element_text(size = 8, hjust = 0)
    )
}

# -----------------------------------------------------------------------------
# 4. RUN EVENT STUDIES FOR ALL COHORTS
# -----------------------------------------------------------------------------

cat("=== RUNNING EVENT STUDIES FOR ALL COHORTS (WITH MONTH DUMMIES) ===\n\n")

# Define outcomes to analyze
primary_outcomes <- c("insulin_oop_per_supply", "insulin_copay_per_supply")
adherence_outcomes <- c("insulin_standardized_supplies", "insulin_gap")
secondary_outcomes <- c("insulin_oop", "insulin_copay")  # Old for comparison
falsification_outcomes <- c("lipid_copay", "metformin_copay")

all_outcomes <- c(primary_outcomes, adherence_outcomes, 
                  secondary_outcomes, falsification_outcomes)

# Run event studies for all cohorts and outcomes
all_event_studies <- map(names(all_panels), function(cohort_name) {
  cat("Processing cohort:", cohort_name, "\n")
  
  cohort_data <- all_panels[[cohort_name]]
  
  # Run event study for each outcome
  cohort_results <- map(all_outcomes, function(outcome) {
    tryCatch({
      run_event_study(cohort_data, outcome)
    }, error = function(e) {
      cat("  Error in", outcome, ":", e$message, "\n")
      return(NULL)
    })
  })
  names(cohort_results) <- all_outcomes
  
  # Remove failed results
  cohort_results <- cohort_results[!sapply(cohort_results, is.null)]
  
  return(cohort_results)
})

names(all_event_studies) <- names(all_panels)

# -----------------------------------------------------------------------------
# 5. EXTRACT AND SAVE SEASONALITY PATTERNS
# -----------------------------------------------------------------------------

cat("\n=== EXTRACTING SEASONALITY PATTERNS ===\n\n")

# Extract month effects for original cohort, primary outcome
if("original_18_64_vs_65_plus" %in% names(all_event_studies)) {
  if("insulin_oop_per_supply" %in% names(all_event_studies$original_18_64_vs_65_plus)) {
    
    model_oop <- all_event_studies$original_18_64_vs_65_plus$insulin_oop_per_supply$model
    month_effects_oop <- extract_month_effects(model_oop)
    
    write_csv(month_effects_oop, "event_study_seasonality_insulin_oop.csv")
    
    # Create seasonality plot
    plot_seasonality <- create_seasonality_plot(
      month_effects_oop,
      "Seasonality in Insulin Out-of-Pocket Costs"
    )
    
    ggsave("event_study_seasonality_pattern.png", plot_seasonality, 
           width = 10, height = 6, dpi = 300)
    
    cat("Seasonality pattern saved: event_study_seasonality_pattern.png\n")
  }
}

# Extract month effects for all outcomes (original cohort only)
if("original_18_64_vs_65_plus" %in% names(all_event_studies)) {
  month_effects_all <- map_df(names(all_event_studies$original_18_64_vs_65_plus), function(outcome) {
    model <- all_event_studies$original_18_64_vs_65_plus[[outcome]]$model
    if(!is.null(model)) {
      extract_month_effects(model) %>%
        mutate(outcome = outcome)
    }
  })
  
  write_csv(month_effects_all, "event_study_seasonality_all_outcomes.csv")
}

# -----------------------------------------------------------------------------
# 6. SAVE ALL COEFFICIENTS
# -----------------------------------------------------------------------------

cat("\n=== SAVING COEFFICIENTS ===\n\n")

# Save coefficients for each cohort and outcome
for(cohort_name in names(all_event_studies)) {
  for(outcome in names(all_event_studies[[cohort_name]])) {
    coefs <- all_event_studies[[cohort_name]][[outcome]]$coefs
    if(!is.null(coefs) && nrow(coefs) > 0) {
      coefs <- coefs %>% mutate(cohort = cohort_name, outcome = outcome)
      filename <- paste0("event_study_", cohort_name, "_", outcome, "_coefs_with_month_dummies.csv")
      write_csv(coefs, filename)
    }
  }
}

# Create combined coefficients file
all_coefs_combined <- map_df(names(all_event_studies), function(cohort_name) {
  map_df(names(all_event_studies[[cohort_name]]), function(outcome) {
    coefs <- all_event_studies[[cohort_name]][[outcome]]$coefs
    if(!is.null(coefs) && nrow(coefs) > 0) {
      coefs %>% mutate(cohort = cohort_name, outcome = outcome)
    }
  })
})

write_csv(all_coefs_combined, "event_study_all_cohorts_all_outcomes_coefs_with_month_dummies.csv")

# -----------------------------------------------------------------------------
# 7. CREATE PLOTS FOR PRIMARY OUTCOME (ALL COHORTS)
# -----------------------------------------------------------------------------

cat("\n=== CREATING COMPARISON PLOTS ===\n\n")

# Create comparison plot for insulin_oop_per_supply across all cohorts
cohort_labels <- c(
  "original_18_64_vs_65_plus" = "Original (18-64 vs 65+)",
  "restrictive_62_64_vs_65_67" = "Restrictive (62-64 vs 65-67)",
  "restrictive_54_64_vs_65_75" = "Restrictive (54-64 vs 65-75)",
  "restrictive_54_64_vs_65_plus" = "Restrictive (54-64 vs 65+)",
  "restrictive_62_64_vs_65_plus" = "Restrictive (62-64 vs 65+)"
)

# Individual plots for each cohort (primary outcome)
oop_plots <- map(names(all_event_studies), function(cohort_name) {
  if("insulin_oop_per_supply" %in% names(all_event_studies[[cohort_name]])) {
    coefs <- all_event_studies[[cohort_name]]$insulin_oop_per_supply$coefs
    create_event_plot(
      coefs,
      "Insulin OOP per 30-Day Supply",
      "Effect ($)",
      cohort_labels[cohort_name]
    )
  }
})
names(oop_plots) <- names(all_event_studies)
oop_plots <- oop_plots[!sapply(oop_plots, is.null)]

# Combined 5-panel plot
if(length(oop_plots) == 5) {
  combined_oop <- wrap_plots(oop_plots, ncol = 2) +
    plot_annotation(
      title = "Event Study: Insulin OOP per 30-Day Supply Across Cohorts",
      subtitle = "Comparing estimates across different age specifications (with month dummies)",
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  
  ggsave("event_study_oop_per_supply_all_cohorts_with_month_dummies.png", combined_oop, 
         width = 16, height = 20, dpi = 300)
  ggsave("event_study_oop_per_supply_all_cohorts_with_month_dummies.pdf", combined_oop, 
         width = 16, height = 20)
}

# -----------------------------------------------------------------------------
# 8. CREATE DETAILED PLOTS FOR ORIGINAL COHORT
# -----------------------------------------------------------------------------

cat("\n=== CREATING DETAILED PLOTS FOR ORIGINAL COHORT ===\n\n")

if("original_18_64_vs_65_plus" %in% names(all_event_studies)) {
  
  original_results <- all_event_studies$original_18_64_vs_65_plus
  
  # Primary outcomes
  plot_oop <- create_event_plot(
    original_results$insulin_oop_per_supply$coefs,
    "A. Insulin OOP per 30-Day Supply (PRIMARY)",
    "Effect ($)"
  )
  
  plot_copay <- create_event_plot(
    original_results$insulin_copay_per_supply$coefs,
    "B. Insulin Copay per 30-Day Supply",
    "Effect ($)"
  )
  
  # Adherence outcomes
  plot_supplies <- create_event_plot(
    original_results$insulin_standardized_supplies$coefs,
    "C. Insulin Supplies (30-Day Equivalents)",
    "Effect (Supplies)"
  )
  
  plot_gaps <- create_event_plot(
    original_results$insulin_gap$coefs,
    "D. Probability of Treatment Gap",
    "Effect (pp)"
  )
  
  # Falsification
  plot_lipid <- create_event_plot(
    original_results$lipid_copay$coefs,
    "E. Lipid Drug Copay (Falsification)",
    "Effect ($)"
  )
  
  plot_metformin <- create_event_plot(
    original_results$metformin_copay$coefs,
    "F. Metformin Copay (Falsification)",
    "Effect ($)"
  )
  
  # Save individual plots
  ggsave("event_study_original_oop_per_supply_with_month_dummies.png", plot_oop, 
         width = 10, height = 6, dpi = 300)
  ggsave("event_study_original_supplies_with_month_dummies.png", plot_supplies, 
         width = 10, height = 6, dpi = 300)
  ggsave("event_study_original_gaps_with_month_dummies.png", plot_gaps, 
         width = 10, height = 6, dpi = 300)
  
  # Combined main figure (4-panel)
  combined_main <- (plot_oop | plot_copay) / (plot_supplies | plot_gaps) +
    plot_annotation(
      title = "Event Study Analysis: IRA Insulin Copay Cap Impact (Original Cohort)",
      subtitle = "Treatment Effects by Month (MA 65+ vs Commercial 18-64) - Controls for Seasonality",
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  
  ggsave("event_study_original_combined_main_with_month_dummies.png", combined_main, 
         width = 16, height = 12, dpi = 300)
  ggsave("event_study_original_combined_main_with_month_dummies.pdf", combined_main, 
         width = 16, height = 12)
  
  # Falsification tests (2-panel)
  combined_false <- (plot_lipid | plot_metformin) +
    plot_annotation(
      title = "Falsification Tests: Non-Insulin Medications (Original Cohort)",
      subtitle = "Should show no systematic effects if identification is valid"
    )
  
  ggsave("event_study_original_falsification_with_month_dummies.png", combined_false, 
         width = 16, height = 6, dpi = 300)
  ggsave("event_study_original_falsification_with_month_dummies.pdf", combined_false, 
         width = 16, height = 6)
}

# -----------------------------------------------------------------------------
# 9. CREATE OVERLAY PLOT (ALL COHORTS, PRIMARY OUTCOME)
# -----------------------------------------------------------------------------

cat("\n=== CREATING OVERLAY COMPARISON PLOT ===\n\n")

# Combine coefficients for insulin_oop_per_supply from all cohorts
oop_coefs_all <- all_coefs_combined %>%
  filter(outcome == "insulin_oop_per_supply") %>%
  mutate(cohort_label = cohort_labels[cohort])

# Overlay plot
overlay_plot <- ggplot(oop_coefs_all, aes(x = relative_month, y = estimate, 
                                          color = cohort_label, group = cohort_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "solid", color = "red", linewidth = 0.8) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, alpha = 0.5) +
  scale_x_continuous(
    breaks = seq(-12, 9, 3),
    labels = c("Jan\n2022", "Apr\n2022", "Jul\n2022", "Oct\n2022", 
               "Jan\n2023", "Apr\n2023", "Jul\n2023", "Oct\n2023")
  ) +
  scale_color_brewer(palette = "Set1", name = "Cohort") +
  labs(
    title = "Event Study Comparison: Insulin OOP per 30-Day Supply",
    subtitle = "Robustness check across different age specifications (with month dummies)",
    x = "Month",
    y = "Effect ($)",
    caption = "Reference: December 2022. 95% CI shown. Red line = IRA implementation.\nMonth dummies control for seasonality."
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.caption = element_text(size = 8, hjust = 0)
  )

ggsave("event_study_overlay_all_cohorts_with_month_dummies.png", overlay_plot, 
       width = 12, height = 8, dpi = 300)
ggsave("event_study_overlay_all_cohorts_with_month_dummies.pdf", overlay_plot, 
       width = 12, height = 8)

# -----------------------------------------------------------------------------
# 10. PRE-TREND TEST (F-TEST)
# -----------------------------------------------------------------------------

cat("\n=== TESTING FOR PRE-TRENDS ===\n\n")

test_pretrends <- function(model) {
  # Extract pre-period coefficients
  coefs <- tidy(model) %>%
    filter(str_detect(term, "relative_month")) %>%
    mutate(
      relative_month = case_when(
        str_detect(term, "::-") ~ as.numeric(str_extract(term, "(?<=::)-\\d+")),
        str_detect(term, "::") ~ as.numeric(str_extract(term, "(?<=::)\\d+")),
        TRUE ~ NA_real_
      )
    ) %>%
    filter(!is.na(relative_month), relative_month < 0)
  
  # F-test: all pre-period coefficients = 0
  if(nrow(coefs) > 0) {
    # Get coefficient names
    pre_terms <- paste0("relative_month::", coefs$relative_month, ":treatment")
    
    # Test joint significance
    test_result <- tryCatch({
      wald(model, pre_terms)
    }, error = function(e) {
      return(NULL)
    })
    
    if(!is.null(test_result)) {
      return(tibble(
        f_stat = test_result$stat,
        p_value = test_result$p,
        df1 = test_result$df1,
        df2 = test_result$df2
      ))
    }
  }
  return(NULL)
}

# Test pre-trends for original cohort
if("original_18_64_vs_65_plus" %in% names(all_event_studies)) {
  pretrend_tests <- map_df(names(all_event_studies$original_18_64_vs_65_plus), function(outcome) {
    model <- all_event_studies$original_18_64_vs_65_plus[[outcome]]$model
    test_result <- test_pretrends(model)
    if(!is.null(test_result)) {
      test_result %>% mutate(outcome = outcome)
    }
  })
  
  write_csv(pretrend_tests, "event_study_pretrend_tests_original.csv")
  
  cat("\nPre-trend test results (Original cohort):\n")
  print(pretrend_tests %>% select(outcome, f_stat, p_value))
}

# -----------------------------------------------------------------------------
# 11. SUMMARY OUTPUT
# -----------------------------------------------------------------------------

cat("\n\n=== EVENT STUDY ANALYSIS COMPLETE (WITH MONTH DUMMIES) ===\n\n")

cat("COHORTS ANALYZED:\n")
for(cohort_name in names(all_panels)) {
  n_patients <- n_distinct(all_panels[[cohort_name]]$pat_id)
  cat("  ", cohort_labels[cohort_name], ": ", n_patients, " patients\n", sep = "")
}

cat("\nOUTCOMES ANALYZED:\n")
cat("  Primary: insulin_oop_per_supply, insulin_copay_per_supply\n")
cat("  Adherence: insulin_standardized_supplies, insulin_gap\n")
cat("  Secondary: insulin_oop, insulin_copay (monthly totals)\n")
cat("  Falsification: lipid_copay, metformin_copay\n\n")

cat("KEY METHODOLOGICAL UPDATE:\n")
cat("â Calendar month dummies included in all event study models\n")
cat("â Controls for January deductible resets, summer utilization drops, etc.\n")
cat("â Seasonality patterns extracted and visualized separately\n")
cat("â Pre-trend F-tests conducted for parallel trends assumption\n\n")

cat("FILES CREATED:\n")
cat("COEFFICIENTS:\n")
cat("  - event_study_all_cohorts_all_outcomes_coefs_with_month_dummies.csv\n")
cat("  - Individual files: event_study_[cohort]_[outcome]_coefs_with_month_dummies.csv\n\n")

cat("SEASONALITY ANALYSIS:\n")
cat("  - event_study_seasonality_insulin_oop.csv (month effects, primary outcome)\n")
cat("  - event_study_seasonality_all_outcomes.csv (month effects, all outcomes)\n")
cat("  - event_study_seasonality_pattern.png (visualization)\n\n")

cat("PRE-TREND TESTS:\n")
cat("  - event_study_pretrend_tests_original.csv (F-tests for parallel trends)\n\n")

cat("PLOTS:\n")
cat("  - event_study_oop_per_supply_all_cohorts_with_month_dummies.* (5-panel)\n")
cat("  - event_study_overlay_all_cohorts_with_month_dummies.* (overlay)\n")
cat("  - event_study_original_combined_main_with_month_dummies.* (4-panel main)\n")
cat("  - event_study_original_falsification_with_month_dummies.* (falsification)\n")
cat("  - event_study_seasonality_pattern.png (seasonal patterns)\n")
cat("  - Individual outcome plots for original cohort\n\n")

cat("FOR YOUR PAPER:\n")
cat("  - Main text: event_study_original_combined_main_with_month_dummies.*\n")
cat("  - Robustness: event_study_overlay_all_cohorts_with_month_dummies.*\n")
cat("  - Appendix: event_study_original_falsification_with_month_dummies.*\n")
cat("  - Supplementary: event_study_seasonality_pattern.png\n\n")

cat("INTERPRETATION:\n")
cat("  - Compare point estimates across cohorts for robustness\n")
cat("  - Check pre-trend F-tests: p-value > 0.05 supports parallel trends\n")
cat("  - Examine seasonality plot to understand cyclical patterns\n")
cat("  - Month dummies ensure treatment effects aren't confounded by seasons\n")
cat("  - Similar patterns across cohorts strengthen causal interpretation\n\n")