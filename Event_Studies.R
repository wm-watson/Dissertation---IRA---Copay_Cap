# =============================================================================
# EVENT STUDY ANALYSIS - ALL COHORTS (WITH INTERACTED MONTH EFFECTS)
# =============================================================================
# Runs event studies for:
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
# 3. DEFINE EVENT STUDY FUNCTIONS (WITH INTERACTED MONTH EFFECTS)
# -----------------------------------------------------------------------------

# SIMPLE SPECIFICATION: Common month dummies
run_event_study_simple <- function(data, outcome, omit_month = -1) {
  # Event study with patient FE + simple calendar month dummies
  
  formula_str <- paste0(
    outcome, " ~ ",
    paste0("i(relative_month, treatment, ref = ", omit_month, ")"),
    " + calendar_month",  # Simple month dummies
    " + n_charlson + dcsi_total + n_unique_drugs_excl_insulin",
    " | pat_id"
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

# INTERACTED SPECIFICATION: Month-by-cohort FE
run_event_study_interacted <- function(data, outcome, omit_month = -1) {
  # Event study with patient FE + month-by-cohort FE
  # Allows MA and Commercial to have different seasonal patterns
  
  formula_str <- paste0(
    outcome, " ~ ",
    paste0("i(relative_month, treatment, ref = ", omit_month, ")"),
    " + calendar_month*treatment",  # INTERACTED: Month-by-cohort FE
    " + n_charlson + dcsi_total + n_unique_drugs_excl_insulin",
    " | pat_id"
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

create_event_plot <- function(coefs, title, ylab, cohort_label = NULL, specification = NULL) {
  subtitle_parts <- c()
  if(!is.null(cohort_label)) subtitle_parts <- c(subtitle_parts, paste0("Cohort: ", cohort_label))
  if(!is.null(specification)) subtitle_parts <- c(subtitle_parts, specification)
  subtitle <- if(length(subtitle_parts) > 0) paste(subtitle_parts, collapse = " | ") else NULL
  
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
      caption = "Reference: December 2022 (month -1). 95% CI shown. Red line = IRA implementation."
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
# 4. RUN EVENT STUDIES FOR ALL COHORTS (BOTH SPECIFICATIONS)
# -----------------------------------------------------------------------------

cat("=== RUNNING EVENT STUDIES FOR ALL COHORTS (BOTH SPECIFICATIONS) ===\n\n")

# Define outcomes to analyze
primary_outcomes <- c("insulin_oop_per_supply", "insulin_copay_per_supply")
adherence_outcomes <- c("insulin_standardized_supplies", "insulin_gap")
secondary_outcomes <- c("insulin_oop", "insulin_copay")  # Old for comparison
falsification_outcomes <- c("lipid_copay", "metformin_copay")

all_outcomes <- c(primary_outcomes, adherence_outcomes, 
                  secondary_outcomes, falsification_outcomes)

# Run BOTH specifications for all cohorts and outcomes
all_event_studies_simple <- map(names(all_panels), function(cohort_name) {
  cat("Processing cohort (SIMPLE):", cohort_name, "\n")
  
  cohort_data <- all_panels[[cohort_name]]
  
  # Run event study for each outcome
  cohort_results <- map(all_outcomes, function(outcome) {
    tryCatch({
      run_event_study_simple(cohort_data, outcome)
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

names(all_event_studies_simple) <- names(all_panels)

all_event_studies_interacted <- map(names(all_panels), function(cohort_name) {
  cat("Processing cohort (INTERACTED):", cohort_name, "\n")
  
  cohort_data <- all_panels[[cohort_name]]
  
  # Run event study for each outcome
  cohort_results <- map(all_outcomes, function(outcome) {
    tryCatch({
      run_event_study_interacted(cohort_data, outcome)
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

names(all_event_studies_interacted) <- names(all_panels)

# Use interacted as primary specification
all_event_studies <- all_event_studies_interacted

# -----------------------------------------------------------------------------
# 5. EXTRACT AND SAVE SEASONALITY PATTERNS
# -----------------------------------------------------------------------------

cat("\n=== EXTRACTING SEASONALITY PATTERNS ===\n\n")

# Extract month effects for original cohort, primary outcome (simple spec)
if("original_18_64_vs_65_plus" %in% names(all_event_studies_simple)) {
  if("insulin_oop_per_supply" %in% names(all_event_studies_simple$original_18_64_vs_65_plus)) {
    
    model_oop_simple <- all_event_studies_simple$original_18_64_vs_65_plus$insulin_oop_per_supply$model
    month_effects_oop_simple <- extract_month_effects(model_oop_simple)
    
    write_csv(month_effects_oop_simple, "event_study_seasonality_insulin_oop_simple.csv")
    
    # Create seasonality plot
    plot_seasonality_simple <- create_seasonality_plot(
      month_effects_oop_simple,
      "Seasonality in Insulin OOP Costs (Simple Month Dummies)"
    )
    
    ggsave("event_study_seasonality_pattern_simple.png", plot_seasonality_simple, 
           width = 10, height = 6, dpi = 300)
    
    cat("Simple seasonality pattern saved: event_study_seasonality_pattern_simple.png\n")
  }
}

# Extract month effects for all outcomes (simple spec, original cohort)
if("original_18_64_vs_65_plus" %in% names(all_event_studies_simple)) {
  month_effects_all_simple <- map_df(names(all_event_studies_simple$original_18_64_vs_65_plus), function(outcome) {
    model <- all_event_studies_simple$original_18_64_vs_65_plus[[outcome]]$model
    if(!is.null(model)) {
      extract_month_effects(model) %>%
        mutate(outcome = outcome)
    }
  })
  
  write_csv(month_effects_all_simple, "event_study_seasonality_all_outcomes_simple.csv")
}

# Note: Interacted spec doesn't have simple month effects to extract (they're group-specific)

# -----------------------------------------------------------------------------
# 6. SAVE ALL COEFFICIENTS
# -----------------------------------------------------------------------------

cat("\n=== SAVING COEFFICIENTS ===\n\n")

# Save coefficients for each cohort and outcome (INTERACTED spec)
for(cohort_name in names(all_event_studies)) {
  for(outcome in names(all_event_studies[[cohort_name]])) {
    coefs <- all_event_studies[[cohort_name]][[outcome]]$coefs
    if(!is.null(coefs) && nrow(coefs) > 0) {
      coefs <- coefs %>% mutate(cohort = cohort_name, outcome = outcome)
      filename <- paste0("event_study_", cohort_name, "_", outcome, "_coefs_interacted.csv")
      write_csv(coefs, filename)
    }
  }
}

# Create combined coefficients file (INTERACTED)
all_coefs_combined_interacted <- map_df(names(all_event_studies), function(cohort_name) {
  map_df(names(all_event_studies[[cohort_name]]), function(outcome) {
    coefs <- all_event_studies[[cohort_name]][[outcome]]$coefs
    if(!is.null(coefs) && nrow(coefs) > 0) {
      coefs %>% mutate(cohort = cohort_name, outcome = outcome)
    }
  })
})

write_csv(all_coefs_combined_interacted, "event_study_all_cohorts_all_outcomes_coefs_interacted.csv")

# Also save SIMPLE spec coefficients for comparison
all_coefs_combined_simple <- map_df(names(all_event_studies_simple), function(cohort_name) {
  map_df(names(all_event_studies_simple[[cohort_name]]), function(outcome) {
    coefs <- all_event_studies_simple[[cohort_name]][[outcome]]$coefs
    if(!is.null(coefs) && nrow(coefs) > 0) {
      coefs %>% mutate(cohort = cohort_name, outcome = outcome)
    }
  })
})

write_csv(all_coefs_combined_simple, "event_study_all_cohorts_all_outcomes_coefs_simple.csv")

# -----------------------------------------------------------------------------
# 7. CREATE PLOTS FOR PRIMARY OUTCOME (ALL COHORTS) - INTERACTED SPEC
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

# Individual plots for each cohort (primary outcome, interacted spec)
oop_plots <- map(names(all_event_studies), function(cohort_name) {
  if("insulin_oop_per_supply" %in% names(all_event_studies[[cohort_name]])) {
    coefs <- all_event_studies[[cohort_name]]$insulin_oop_per_supply$coefs
    create_event_plot(
      coefs,
      "Insulin OOP per 30-Day Supply",
      "Effect ($)",
      cohort_labels[cohort_name],
      "Interacted month-by-cohort FE"
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
      subtitle = "Interacted month-by-cohort FE specification (allows different seasonal patterns)",
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  
  ggsave("event_study_oop_per_supply_all_cohorts_interacted.png", combined_oop, 
         width = 16, height = 20, dpi = 300)
  ggsave("event_study_oop_per_supply_all_cohorts_interacted.pdf", combined_oop, 
         width = 16, height = 20)
}

# -----------------------------------------------------------------------------
# 8. CREATE DETAILED PLOTS FOR ORIGINAL COHORT (INTERACTED SPEC)
# -----------------------------------------------------------------------------

cat("\n=== CREATING DETAILED PLOTS FOR ORIGINAL COHORT ===\n\n")

if("original_18_64_vs_65_plus" %in% names(all_event_studies)) {
  
  original_results <- all_event_studies$original_18_64_vs_65_plus
  
  # Primary outcomes
  plot_oop <- create_event_plot(
    original_results$insulin_oop_per_supply$coefs,
    "A. Insulin OOP per 30-Day Supply (PRIMARY)",
    "Effect ($)",
    specification = "Interacted month-by-cohort FE"
  )
  
  plot_copay <- create_event_plot(
    original_results$insulin_copay_per_supply$coefs,
    "B. Insulin Copay per 30-Day Supply",
    "Effect ($)",
    specification = "Interacted month-by-cohort FE"
  )
  
  # Adherence outcomes
  plot_supplies <- create_event_plot(
    original_results$insulin_standardized_supplies$coefs,
    "C. Insulin Supplies (30-Day Equivalents)",
    "Effect (Supplies)",
    specification = "Interacted month-by-cohort FE"
  )
  
  plot_gaps <- create_event_plot(
    original_results$insulin_gap$coefs,
    "D. Probability of Treatment Gap",
    "Effect (pp)",
    specification = "Interacted month-by-cohort FE"
  )
  
  # Falsification
  plot_lipid <- create_event_plot(
    original_results$lipid_copay$coefs,
    "E. Lipid Drug Copay (Falsification)",
    "Effect ($)",
    specification = "Interacted month-by-cohort FE"
  )
  
  plot_metformin <- create_event_plot(
    original_results$metformin_copay$coefs,
    "F. Metformin Copay (Falsification)",
    "Effect ($)",
    specification = "Interacted month-by-cohort FE"
  )
  
  # Save individual plots
  ggsave("event_study_original_oop_per_supply_interacted.png", plot_oop, 
         width = 10, height = 6, dpi = 300)
  ggsave("event_study_original_supplies_interacted.png", plot_supplies, 
         width = 10, height = 6, dpi = 300)
  ggsave("event_study_original_gaps_interacted.png", plot_gaps, 
         width = 10, height = 6, dpi = 300)
  
  # Combined main figure (4-panel)
  combined_main <- (plot_oop | plot_copay) / (plot_supplies | plot_gaps) +
    plot_annotation(
      title = "Event Study Analysis: IRA Insulin Copay Cap Impact (Original Cohort)",
      subtitle = "Treatment Effects by Month (MA 65+ vs Commercial 18-64) - Interacted Month-by-Cohort FE",
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  
  ggsave("event_study_original_combined_main_interacted.png", combined_main, 
         width = 16, height = 12, dpi = 300)
  ggsave("event_study_original_combined_main_interacted.pdf", combined_main, 
         width = 16, height = 12)
  
  # Falsification tests (2-panel)
  combined_false <- (plot_lipid | plot_metformin) +
    plot_annotation(
      title = "Falsification Tests: Non-Insulin Medications (Original Cohort)",
      subtitle = "Should show no systematic effects if identification is valid"
    )
  
  ggsave("event_study_original_falsification_interacted.png", combined_false, 
         width = 16, height = 6, dpi = 300)
  ggsave("event_study_original_falsification_interacted.pdf", combined_false, 
         width = 16, height = 6)
}

# -----------------------------------------------------------------------------
# 9. SPECIFICATION COMPARISON: SIMPLE VS INTERACTED (ORIGINAL COHORT)
# -----------------------------------------------------------------------------

cat("\n=== CREATING SPECIFICATION COMPARISON PLOTS ===\n\n")

if("original_18_64_vs_65_plus" %in% names(all_event_studies)) {
  
  # Get coefficients from both specs
  coefs_simple <- all_event_studies_simple$original_18_64_vs_65_plus$insulin_oop_per_supply$coefs %>%
    mutate(specification = "Simple month dummies")
  
  coefs_interacted <- all_event_studies$original_18_64_vs_65_plus$insulin_oop_per_supply$coefs %>%
    mutate(specification = "Interacted month-by-cohort FE")
  
  coefs_comparison <- bind_rows(coefs_simple, coefs_interacted)
  
  # Comparison plot
  plot_spec_comparison <- ggplot(coefs_comparison, 
                                 aes(x = relative_month, y = estimate, 
                                     color = specification, group = specification)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "solid", color = "red", linewidth = 0.8) +
    geom_point(size = 2, position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2,
                  position = position_dodge(width = 0.3)) +
    geom_line(linewidth = 1, position = position_dodge(width = 0.3)) +
    scale_color_manual(values = c("Simple month dummies" = "#d7191c",
                                  "Interacted month-by-cohort FE" = "#2c7bb6")) +
    scale_x_continuous(
      breaks = seq(-12, 9, 3),
      labels = c("Jan\n2022", "Apr\n2022", "Jul\n2022", "Oct\n2022", 
                 "Jan\n2023", "Apr\n2023", "Jul\n2023", "Oct\n2023")
    ) +
    labs(
      title = "Event Study: Specification Comparison",
      subtitle = "Testing sensitivity to seasonal pattern assumptions (Original Cohort)",
      x = "Month",
      y = "Effect on OOP per 30-Day Supply ($)",
      color = "Specification",
      caption = "Reference: December 2022. Red line = IRA implementation.\nInteracted spec allows MA and Commercial to have different seasonality."
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "top",
      panel.grid.minor = element_blank(),
      plot.caption = element_text(size = 9, hjust = 0)
    )
  
  ggsave("event_study_specification_comparison.png", plot_spec_comparison,
         width = 12, height = 8, dpi = 300)
  ggsave("event_study_specification_comparison.pdf", plot_spec_comparison,
         width = 12, height = 8)
  
  cat("\nSpecification comparison plot saved: event_study_specification_comparison.png\n")
  
  # Save comparison coefficients
  write_csv(coefs_comparison, "event_study_specification_comparison_coefs.csv")
}

# -----------------------------------------------------------------------------
# 10. CREATE OVERLAY PLOT (ALL COHORTS, PRIMARY OUTCOME, INTERACTED SPEC)
# -----------------------------------------------------------------------------

cat("\n=== CREATING OVERLAY COMPARISON PLOT ===\n\n")

# Combine coefficients for insulin_oop_per_supply from all cohorts (interacted)
oop_coefs_all <- all_coefs_combined_interacted %>%
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
    subtitle = "Robustness check across age specifications (Interacted month-by-cohort FE)",
    x = "Month",
    y = "Effect ($)",
    caption = "Reference: December 2022. 95% CI shown. Red line = IRA implementation.\nInteracted month-by-cohort FE allows different seasonal patterns by insurance type."
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

ggsave("event_study_overlay_all_cohorts_interacted.png", overlay_plot, 
       width = 12, height = 8, dpi = 300)
ggsave("event_study_overlay_all_cohorts_interacted.pdf", overlay_plot, 
       width = 12, height = 8)

# -----------------------------------------------------------------------------
# 11. PRE-TREND TEST (F-TEST)
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

# Test pre-trends for original cohort (BOTH specifications)
if("original_18_64_vs_65_plus" %in% names(all_event_studies)) {
  
  # Simple spec
  pretrend_tests_simple <- map_df(names(all_event_studies_simple$original_18_64_vs_65_plus), function(outcome) {
    model <- all_event_studies_simple$original_18_64_vs_65_plus[[outcome]]$model
    test_result <- test_pretrends(model)
    if(!is.null(test_result)) {
      test_result %>% mutate(outcome = outcome, specification = "Simple")
    }
  })
  
  # Interacted spec
  pretrend_tests_interacted <- map_df(names(all_event_studies$original_18_64_vs_65_plus), function(outcome) {
    model <- all_event_studies$original_18_64_vs_65_plus[[outcome]]$model
    test_result <- test_pretrends(model)
    if(!is.null(test_result)) {
      test_result %>% mutate(outcome = outcome, specification = "Interacted")
    }
  })
  
  # Combine
  pretrend_tests_all <- bind_rows(pretrend_tests_simple, pretrend_tests_interacted)
  
  write_csv(pretrend_tests_all, "event_study_pretrend_tests_both_specs.csv")
  
  cat("\nPre-trend test results (Original cohort - both specifications):\n")
  print(pretrend_tests_all %>% select(specification, outcome, f_stat, p_value))
  
  cat("\nInterpretation:\n")
  cat("- p-value > 0.05: Parallel trends assumption holds (good!)\n")
  cat("- p-value < 0.05: Pre-trends differ significantly (violation of parallel trends)\n")
  cat("- Compare simple vs interacted: Does allowing different seasonality help?\n\n")
}

# -----------------------------------------------------------------------------
# 12. SUMMARY OUTPUT
# -----------------------------------------------------------------------------

cat("\n\n=== EVENT STUDY ANALYSIS COMPLETE (WITH INTERACTED MONTH EFFECTS) ===\n\n")

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
cat("✓ Calendar month × treatment interactions (month-by-cohort FE)\n")
cat("✓ Allows MA and Commercial to have DIFFERENT seasonal patterns\n")
cat("✓ Controls for Part D donut hole, formulary timing differences\n")
cat("✓ Addresses potential parallel trends violations from seasonal heterogeneity\n")
cat("✓ Pre-trend F-tests conducted for both simple and interacted specs\n\n")

cat("FILES CREATED:\n")
cat("COEFFICIENTS:\n")
cat("  - event_study_all_cohorts_all_outcomes_coefs_interacted.csv (INTERACTED spec)\n")
cat("  - event_study_all_cohorts_all_outcomes_coefs_simple.csv (SIMPLE spec)\n")
cat("  - event_study_specification_comparison_coefs.csv (side-by-side)\n")
cat("  - Individual files: event_study_[cohort]_[outcome]_coefs_interacted.csv\n\n")

cat("SEASONALITY ANALYSIS:\n")
cat("  - event_study_seasonality_insulin_oop_simple.csv (month effects, simple spec)\n")
cat("  - event_study_seasonality_all_outcomes_simple.csv (all outcomes, simple spec)\n")
cat("  - event_study_seasonality_pattern_simple.png (visualization)\n\n")

cat("PRE-TREND TESTS:\n")
cat("  - event_study_pretrend_tests_both_specs.csv (F-tests for parallel trends)\n\n")

cat("PLOTS:\n")
cat("  INTERACTED SPEC (PRIMARY):\n")
cat("  - event_study_oop_per_supply_all_cohorts_interacted.* (5-panel)\n")
cat("  - event_study_overlay_all_cohorts_interacted.* (overlay)\n")
cat("  - event_study_original_combined_main_interacted.* (4-panel main)\n")
cat("  - event_study_original_falsification_interacted.* (falsification)\n\n")
cat("  COMPARISON:\n")
cat("  - event_study_specification_comparison.* (simple vs interacted)\n\n")
cat("  SEASONALITY:\n")
cat("  - event_study_seasonality_pattern_simple.png (simple spec only)\n\n")

cat("FOR YOUR PAPER:\n")
cat("  - Main text: event_study_original_combined_main_interacted.*\n")
cat("  - Robustness: event_study_overlay_all_cohorts_interacted.*\n")
cat("  - Appendix: event_study_original_falsification_interacted.*\n")
cat("  - Specification check: event_study_specification_comparison.*\n\n")

cat("INTERPRETATION:\n")
cat("  - Compare simple vs interacted event studies\n")
cat("  - Check pre-trend F-tests: Does interacted spec improve parallel trends?\n")
cat("  - Examine specification comparison plot: Are estimates robust?\n")
cat("  - If estimates similar: Seasonality not driving results\n")
cat("  - If estimates differ: Seasonal heterogeneity matters\n")
cat("  - Similar patterns across cohorts strengthen causal interpretation\n\n")

cat("NEXT STEPS:\n")
cat("  1. Compare pre-trend tests: simple vs interacted (which passes better?)\n")
cat("  2. Look at specification comparison plot: stable estimates?\n")
cat("  3. If interacted improves pre-trends: use as primary specification\n")
cat("  4. Report both specs as robustness check\n\n")