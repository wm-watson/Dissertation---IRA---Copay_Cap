library(tidyverse)
library(haven)
library(lubridate)

# =============================================================================
# MONTHLY ANALYTICAL DATASETS - IRA INSULIN COPAY CAP
# WITH COMPREHENSIVE TRACKING FOR WRITE-UP
# =============================================================================

# -----------------------------------------------------------------------------
# TRACKING FRAMEWORK
# -----------------------------------------------------------------------------

# Create tracking list to store all counts
tracking <- list()

track_step <- function(step, description, data = NULL, n_patients = NULL, 
                       n_patient_months = NULL, n_claims = NULL) {
  entry <- list(
    step = step,
    description = description,
    n_patients = n_patients %||% (if(!is.null(data)) n_distinct(data$pat_id) else NA),
    n_patient_months = n_patient_months %||% (if(!is.null(data) && "year_month" %in% names(data)) nrow(data) else NA),
    n_claims = n_claims
  )
  
  tracking[[length(tracking) + 1]] <<- entry
  
  cat(sprintf("\n[Step %d] %s\n", step, description))
  if (!is.na(entry$n_patients)) cat(sprintf("  Patients: %s\n", format(entry$n_patients, big.mark = ",")))
  if (!is.na(entry$n_patient_months)) cat(sprintf("  Patient-months: %s\n", format(entry$n_patient_months, big.mark = ",")))
  if (!is.na(entry$n_claims)) cat(sprintf("  Claims: %s\n", format(entry$n_claims, big.mark = ",")))
}

cat("\n=============================================================================\n")
cat("MONTHLY ANALYTICAL DATASET CONSTRUCTION WITH TRACKING\n")
cat("=============================================================================\n")

# -----------------------------------------------------------------------------
# 1. LOAD DATA WITH TRACKING
# -----------------------------------------------------------------------------

cat("\n--- STEP 1: LOADING DATA ---\n")

# Load cohorts
cohort_cross_raw <- read_sas("cohort_cross_sectional.sas7bdat")
cohort_panel_raw <- read_sas("cohort_panel_22.sas7bdat")

track_step(1.1, "Cross-sectional cohorts loaded (before sex filter)",
           n_patients = n_distinct(c(
             cohort_cross_raw$pat_id[cohort_cross_raw$in_cohort_2022 == 1],
             cohort_cross_raw$pat_id[cohort_cross_raw$in_cohort_2023 == 1]
           )))

track_step(1.2, "Panel cohort loaded (before sex filter)",
           data = cohort_panel_raw %>% filter(in_cohort_panel == 1))

# Filter out unknown sex
cohort_cross <- cohort_cross_raw %>% filter(der_sex != "U")
cohort_panel <- cohort_panel_raw %>% filter(der_sex != "U")

n_excluded_sex_cross <- n_distinct(c(
  cohort_cross_raw$pat_id[cohort_cross_raw$in_cohort_2022 == 1],
  cohort_cross_raw$pat_id[cohort_cross_raw$in_cohort_2023 == 1]
)) - n_distinct(c(
  cohort_cross$pat_id[cohort_cross$in_cohort_2022 == 1],
  cohort_cross$pat_id[cohort_cross$in_cohort_2023 == 1]
))

n_excluded_sex_panel <- n_distinct(cohort_panel_raw$pat_id[cohort_panel_raw$in_cohort_panel == 1]) -
  n_distinct(cohort_panel$pat_id[cohort_panel$in_cohort_panel == 1])

cat("\n=== SEX EXCLUSION ===\n")
cat(sprintf("Cross-sectional: Excluded %d patients with unknown sex (der_sex='U')\n", n_excluded_sex_cross))
cat(sprintf("Panel: Excluded %d patients with unknown sex (der_sex='U')\n", n_excluded_sex_panel))

track_step(1.3, "After excluding unknown sex - Cross-sectional 2022",
           data = cohort_cross %>% filter(in_cohort_2022 == 1))

track_step(1.4, "After excluding unknown sex - Cross-sectional 2023",
           data = cohort_cross %>% filter(in_cohort_2023 == 1))

track_step(1.5, "After excluding unknown sex - Panel 22 months",
           data = cohort_panel %>% filter(in_cohort_panel == 1))

# Load claims
rx_2022_raw <- read_sas("master_rx_2022.sas7bdat")
rx_2023_raw <- read_sas("master_rx_2023.sas7bdat")

track_step(1.6, "Pharmacy claims 2022 loaded (before reversed claim filter)",
           n_claims = nrow(rx_2022_raw))

track_step(1.7, "Pharmacy claims 2023 loaded (before reversed claim filter)",
           n_claims = nrow(rx_2023_raw))

# Filter reversed claims
rx_2022 <- rx_2022_raw %>% filter(!(copay == 0 & is.na(quan)))
rx_2023 <- rx_2023_raw %>% filter(!(copay == 0 & is.na(quan)))

cat(sprintf("\nReversed claims excluded 2022: %s\n", 
            format(nrow(rx_2022_raw) - nrow(rx_2022), big.mark = ",")))
cat(sprintf("Reversed claims excluded 2023: %s\n", 
            format(nrow(rx_2023_raw) - nrow(rx_2023), big.mark = ",")))

track_step(1.8, "Pharmacy claims 2022 after removing reversed claims",
           n_claims = nrow(rx_2022))

track_step(1.9, "Pharmacy claims 2023 after removing reversed claims",
           n_claims = nrow(rx_2023))

# Load medical claims
outpat_2022 <- read_sas("master_outpat_2022.sas7bdat")
outpat_2023 <- read_sas("master_outpat_2023.sas7bdat")
inpat_2022 <- read_sas("master_inpat_2022.sas7bdat")
inpat_2023 <- read_sas("master_inpat_2023.sas7bdat")

track_step(1.10, "Medical claims loaded",
           n_claims = nrow(outpat_2022) + nrow(outpat_2023) + 
             nrow(inpat_2022) + nrow(inpat_2023))

# -----------------------------------------------------------------------------
# 2. ADD MONTH IDENTIFIER
# -----------------------------------------------------------------------------

cat("\n--- STEP 2: ADDING TIME IDENTIFIERS ---\n")

add_month_id <- function(data) {
  data %>%
    mutate(
      from_dt = as.Date(from_dt),
      year_month = floor_date(from_dt, "month"),
      year = year(from_dt),
      month = month(from_dt)
    )
}

rx_2022 <- add_month_id(rx_2022)
rx_2023 <- add_month_id(rx_2023)

track_step(2.1, "Time identifiers added to pharmacy claims")

# -----------------------------------------------------------------------------
# 3. DEDUPLICATE AND AGGREGATE PHARMACY CLAIMS
# -----------------------------------------------------------------------------

cat("\n--- STEP 3: DEDUPLICATING AND AGGREGATING PHARMACY CLAIMS ---\n")

aggregate_to_patient_month <- function(rx_data, year) {
  
  # Count before deduplication
  n_before_dedup <- nrow(rx_data)
  
  # STEP 1: Deduplicate at fill level
  rx_deduplicated <- rx_data %>%
    group_by(pat_id, from_dt, ndc, prscbr_id) %>%
    summarise(
      year_month = first(year_month),
      gpi14 = first(gpi14),
      generic_name = first(generic_name),
      product_name = first(product_name),
      copay = max(copay, na.rm = TRUE),
      deductible = max(deductible, na.rm = TRUE),
      coinsamt = max(coinsamt, na.rm = TRUE),
      cobamt = max(cobamt, na.rm = TRUE),
      paid = max(paid, na.rm = TRUE),
      allowed = max(allowed, na.rm = TRUE),
      quan = ifelse(all(is.na(quan)), 0, max(quan, na.rm = TRUE)),
      dayssup = ifelse(all(is.na(dayssup)), 30, max(dayssup, na.rm = TRUE)),
      .groups = "drop"
    )
  
  n_after_dedup <- nrow(rx_deduplicated)
  n_removed_dedup <- n_before_dedup - n_after_dedup
  
  cat(sprintf("\n  %d claims -> %d unique fills (removed %d duplicate lines)\n",
              n_before_dedup, n_after_dedup, n_removed_dedup))
  
  # STEP 2: Add drug classifications AND 30-day standardization
  rx_classified <- rx_deduplicated %>%
    mutate(
      dayssup_clean = case_when(
        is.na(dayssup) | dayssup == 0 ~ 30,
        dayssup > 365 ~ 365,
        dayssup < 1 ~ 1,
        TRUE ~ dayssup
      ),
      
      # Standardize to 30-day supply units
      standardized_supplies = pmax(1, dayssup_clean / 30),
      
      # Calculate OOP per standardized 30-day supply
      oop = copay + deductible + coinsamt + cobamt,
      copay_per_supply = copay / standardized_supplies,
      deductible_per_supply = deductible / standardized_supplies,
      coinsamt_per_supply = coinsamt / standardized_supplies,
      cobamt_per_supply = cobamt / standardized_supplies,
      oop_per_supply = oop / standardized_supplies,
      
      # Normalize names
      product_upper = toupper(product_name),
      generic_upper = toupper(generic_name),
      
      # Insulin classification
      is_insulin_gpi = str_detect(gpi14, "^2710"),
      is_insulin_product = str_detect(product_upper, 
                                      "INSULIN|HUMALOG|NOVOLOG|LANTUS|LEVEMIR|TRESIBA|TOUJEO|BASAGLAR|AFREZZA|SEMGLEE|FIASP|APIDRA|ADMELOG|LYUMJEV|XULTOPHY|SOLIQUA|EXUBERA"),
      is_insulin_generic = str_detect(generic_upper, 
                                      "INSULIN|LISPRO|ASPART|GLARGINE|DETEMIR|DEGLUDEC|GLULISINE"),
      is_insulin = (is_insulin_gpi | is_insulin_product | is_insulin_generic),
      
      # Concentrated insulins
      is_concentrated = case_when(
        !is_insulin ~ FALSE,
        str_detect(product_upper, "U-500|U500|U-300|U300|TOUJEO|U-200|U200") ~ TRUE,
        TRUE ~ FALSE
      ),
      
      # Insulin TYPE classification
      insulin_type = case_when(
        !is_insulin ~ NA_character_,
        str_detect(product_upper, "HUMALOG|ADMELOG|LISPRO|NOVOLOG|FIASP|APIDRA") |
          str_detect(generic_upper, "LISPRO|ASPART|GLULISINE") |
          str_sub(gpi14, 1, 8) %in% c("27104005", "27104002") ~ "rapid_acting",
        str_detect(product_upper, "HUMULIN R|NOVOLIN R|REGULAR|VELOSULIN|MYXREDLIN") |
          str_detect(generic_upper, "REGULAR") |
          str_sub(gpi14, 1, 8) == "27104010" ~ "short_acting",
        str_detect(product_upper, "HUMULIN N|NOVOLIN N|NPH|INSULATARD") |
          str_detect(generic_upper, "NPH|ISOPHANE") |
          str_sub(gpi14, 1, 8) == "27104020" ~ "intermediate_acting",
        str_detect(product_upper, "LANTUS|BASAGLAR|TOUJEO|SEMGLEE|LEVEMIR|TRESIBA") |
          str_detect(generic_upper, "GLARGINE|DETEMIR|DEGLUDEC") |
          str_sub(gpi14, 1, 8) %in% c("27104003", "27104006", "27104007") ~ "long_acting",
        str_detect(product_upper, "70/30|75/25|50/50|MIX") |
          str_detect(generic_upper, "PROTAMINE") |
          str_sub(gpi14, 1, 8) %in% c("27104080", "27104090") ~ "mixed",
        str_detect(product_upper, "AFREZZA|EXUBERA") |
          str_detect(generic_upper, "INHALED") ~ "inhaled",
        str_detect(product_upper, "XULTOPHY|SOLIQUA") |
          str_detect(generic_upper, "LIRAGLUTIDE|LIXISENATIDE") ~ "combination_glp1",
        str_detect(product_upper, "ILETIN|PORK|BEEF|BOVINE") |
          str_detect(generic_upper, "PORK|BEEF|BOVINE") ~ "animal",
        TRUE ~ "other"
      ),
      
      # Cap compliance
      exceeds_cap_per_supply = (oop_per_supply > 35),
      
      # Other classifications
      is_metformin = str_detect(gpi14, "^2730"),
      is_lipid = str_detect(gpi14, "^39"),
      gpi6 = str_sub(gpi14, 1, 6)
    )
  
  # Count insulin fills
  n_insulin_fills <- sum(rx_classified$is_insulin)
  pct_insulin <- round(n_insulin_fills / nrow(rx_classified) * 100, 1)
  
  cat(sprintf("  Insulin fills identified: %s (%.1f%% of total)\n",
              format(n_insulin_fills, big.mark = ","), pct_insulin))
  
  # STEP 3: Aggregate to patient-month
  patient_month <- rx_classified %>%
    group_by(pat_id, year_month) %>%
    summarise(
      n_unique_drugs = n_distinct(gpi6),
      total_fills = n(),
      total_copay = sum(copay, na.rm = TRUE),
      total_deductible = sum(deductible, na.rm = TRUE),
      total_coinsamt = sum(coinsamt, na.rm = TRUE),
      total_cobamt = sum(cobamt, na.rm = TRUE),
      total_oop = sum(oop, na.rm = TRUE),
      total_paid = sum(paid, na.rm = TRUE),
      total_allowed = sum(allowed, na.rm = TRUE),
      
      # Standardized supply metrics
      total_standardized_supplies = sum(standardized_supplies, na.rm = TRUE),
      insulin_standardized_supplies = sum(standardized_supplies[is_insulin], na.rm = TRUE),
      
      # Insulin metrics
      insulin_n_fills = sum(is_insulin),
      insulin_copay = sum(copay[is_insulin], na.rm = TRUE),
      insulin_deductible = sum(deductible[is_insulin], na.rm = TRUE),
      insulin_coinsamt = sum(coinsamt[is_insulin], na.rm = TRUE),
      insulin_cobamt = sum(cobamt[is_insulin], na.rm = TRUE),
      insulin_oop = sum(oop[is_insulin], na.rm = TRUE),
      insulin_paid = sum(paid[is_insulin], na.rm = TRUE),
      insulin_allowed = sum(allowed[is_insulin], na.rm = TRUE),
      insulin_quantity = sum(quan[is_insulin], na.rm = TRUE),
      insulin_n_types = n_distinct(insulin_type[is_insulin & !is.na(insulin_type)]),
      
      # Per-supply costs
      insulin_copay_per_supply = weighted.mean(copay_per_supply[is_insulin], 
                                               standardized_supplies[is_insulin], na.rm = TRUE),
      insulin_oop_per_supply = weighted.mean(oop_per_supply[is_insulin], 
                                             standardized_supplies[is_insulin], na.rm = TRUE),
      
      # Cap compliance
      pct_insulin_supplies_exceeding_cap = mean(exceeds_cap_per_supply[is_insulin], na.rm = TRUE),
      any_insulin_supply_exceeds_cap = any(exceeds_cap_per_supply[is_insulin], na.rm = TRUE),
      
      # Concentrated
      n_concentrated_fills = sum(is_insulin & is_concentrated, na.rm = TRUE),
      
      # By type
      n_fills_rapid_acting = sum(is_insulin & insulin_type == "rapid_acting", na.rm = TRUE),
      n_fills_short_acting = sum(is_insulin & insulin_type == "short_acting", na.rm = TRUE),
      n_fills_intermediate_acting = sum(is_insulin & insulin_type == "intermediate_acting", na.rm = TRUE),
      n_fills_long_acting = sum(is_insulin & insulin_type == "long_acting", na.rm = TRUE),
      n_fills_mixed = sum(is_insulin & insulin_type == "mixed", na.rm = TRUE),
      n_fills_inhaled = sum(is_insulin & insulin_type == "inhaled", na.rm = TRUE),
      n_fills_combination = sum(is_insulin & insulin_type == "combination_glp1", na.rm = TRUE),
      
      supplies_rapid_acting = sum(standardized_supplies[is_insulin & insulin_type == "rapid_acting"], na.rm = TRUE),
      supplies_long_acting = sum(standardized_supplies[is_insulin & insulin_type == "long_acting"], na.rm = TRUE),
      
      copay_rapid_acting = sum(copay[is_insulin & insulin_type == "rapid_acting"], na.rm = TRUE),
      copay_long_acting = sum(copay[is_insulin & insulin_type == "long_acting"], na.rm = TRUE),
      
      oop_rapid_acting = sum(oop[is_insulin & insulin_type == "rapid_acting"], na.rm = TRUE),
      oop_long_acting = sum(oop[is_insulin & insulin_type == "long_acting"], na.rm = TRUE),
      
      # Falsification drugs
      metformin_copay = sum(copay[is_metformin], na.rm = TRUE),
      metformin_oop = sum(oop[is_metformin], na.rm = TRUE),
      metformin_n_fills = sum(is_metformin),
      
      lipid_copay = sum(copay[is_lipid], na.rm = TRUE),
      lipid_oop = sum(oop[is_lipid], na.rm = TRUE),
      lipid_n_fills = sum(is_lipid),
      
      .groups = "drop"
    ) %>%
    mutate(
      exceeds_cap = (insulin_oop_per_supply > 35) & !is.na(insulin_oop_per_supply),
      insulin_copay_per_supply = ifelse(is.nan(insulin_copay_per_supply), 0, insulin_copay_per_supply),
      insulin_oop_per_supply = ifelse(is.nan(insulin_oop_per_supply), 0, insulin_oop_per_supply),
      pct_insulin_supplies_exceeding_cap = ifelse(is.nan(pct_insulin_supplies_exceeding_cap), 0, pct_insulin_supplies_exceeding_cap),
      
      # Adherence metrics
      insulin_gap = (insulin_n_fills == 0),
      fully_covered = (insulin_standardized_supplies >= 1.0) & (insulin_n_fills > 0),
      estimated_days_covered = pmin(insulin_standardized_supplies * 30, 30),
      supply_adequacy_ratio = insulin_standardized_supplies / 1.0,
      n_unique_drugs_excl_insulin = n_unique_drugs - (insulin_n_fills > 0)
    )
  
  return(list(
    data = patient_month,
    n_fills_before_dedup = n_before_dedup,
    n_fills_after_dedup = n_after_dedup,
    n_insulin_fills = n_insulin_fills
  ))
}

cat("\nAggregating 2022 pharmacy claims to patient-month level...\n")
rx_month_2022_result <- aggregate_to_patient_month(rx_2022, 2022)
rx_month_2022 <- rx_month_2022_result$data

track_step(3.1, "Pharmacy claims 2022 aggregated to patient-month",
           data = rx_month_2022,
           n_claims = rx_month_2022_result$n_fills_after_dedup)

cat("\nAggregating 2023 pharmacy claims to patient-month level...\n")
rx_month_2023_result <- aggregate_to_patient_month(rx_2023, 2023)
rx_month_2023 <- rx_month_2023_result$data

track_step(3.2, "Pharmacy claims 2023 aggregated to patient-month",
           data = rx_month_2023,
           n_claims = rx_month_2023_result$n_fills_after_dedup)

# Insulin classification validation
cat("\n=== INSULIN CLASSIFICATION VALIDATION ===\n")
insulin_validation_2022 <- rx_2022 %>%
  mutate(
    product_upper = toupper(product_name),
    gpi_insulin = str_detect(gpi14, "^2710"),
    name_insulin = str_detect(product_upper, 
                              "INSULIN|HUMALOG|NOVOLOG|LANTUS|LEVEMIR|TRESIBA|TOUJEO|BASAGLAR|AFREZZA|SEMGLEE|FIASP|APIDRA|ADMELOG|LYUMJEV|XULTOPHY|SOLIQUA")
  ) %>%
  summarise(
    gpi_only = sum(gpi_insulin & !name_insulin, na.rm = TRUE),
    name_only = sum(!gpi_insulin & name_insulin, na.rm = TRUE),
    both = sum(gpi_insulin & name_insulin, na.rm = TRUE),
    total_insulin = sum(gpi_insulin | name_insulin, na.rm = TRUE)
  )

cat("2022 Classification method overlap:\n")
cat(sprintf("  GPI only (2710*): %s claims\n", format(insulin_validation_2022$gpi_only, big.mark = ",")))
cat(sprintf("  Name only: %s claims\n", format(insulin_validation_2022$name_only, big.mark = ",")))
cat(sprintf("  Both methods: %s claims\n", format(insulin_validation_2022$both, big.mark = ",")))
cat(sprintf("  Total insulin: %s claims\n", format(insulin_validation_2022$total_insulin, big.mark = ",")))

if (insulin_validation_2022$name_only > 0) {
  cat(sprintf("\nNOTE: %s insulin claims captured by name but not by GPI.\n", 
              format(insulin_validation_2022$name_only, big.mark = ",")))
  cat("These may include: inhaled insulin (Afrezza), GLP-1 combinations (Xultophy/Soliqua),\n")
  cat("or products with non-standard GPI codes.\n")
}

# -----------------------------------------------------------------------------
# 4. CREATE COMPLETE PATIENT-MONTH PANELS
# -----------------------------------------------------------------------------

cat("\n--- STEP 4: CREATING COMPLETE PATIENT-MONTH PANELS ---\n")

create_complete_patient_months <- function(cohort_data, cohort_flag_col, 
                                           rx_month_2022, rx_month_2023, year) {
  cohort_patients <- cohort_data %>%
    filter(!!sym(cohort_flag_col) == 1) %>%
    pull(pat_id)
  
  if (year == 2022) {
    month_grid <- expand_grid(
      pat_id = cohort_patients,
      year_month = seq(as.Date("2022-01-01"), as.Date("2022-12-01"), by = "month")
    )
    rx_month <- rx_month_2022
  } else {
    month_grid <- expand_grid(
      pat_id = cohort_patients,
      year_month = seq(as.Date("2023-01-01"), as.Date("2023-10-01"), by = "month")
    )
    rx_month <- rx_month_2023
  }
  
  complete_months <- month_grid %>%
    left_join(rx_month, by = c("pat_id", "year_month")) %>%
    mutate(across(where(is.numeric), ~replace_na(., 0)))
  
  return(complete_months)
}

complete_2022 <- create_complete_patient_months(cohort_cross, "in_cohort_2022", 
                                                rx_month_2022, rx_month_2023, 2022)
complete_2023 <- create_complete_patient_months(cohort_cross, "in_cohort_2023", 
                                                rx_month_2022, rx_month_2023, 2023)

track_step(4.1, "Complete patient-month panel 2022 (12 months)",
           data = complete_2022)

track_step(4.2, "Complete patient-month panel 2023 (10 months)",
           data = complete_2023)

# -----------------------------------------------------------------------------
# 5. CALCULATE COMORBIDITIES
# -----------------------------------------------------------------------------

cat("\n--- STEP 5: CALCULATING COMORBIDITIES ---\n")

calculate_comorbidities <- function(outpat, inpat, year) {
  
  # Extract diagnoses
  outpat_diags <- outpat %>%
    select(pat_id, from_dt, diag1, diag2, diag3, diag4, diag5, diag6, 
           diag7, diag8, diag9, diag10, diag11, diag12) %>%
    pivot_longer(cols = starts_with("diag"), values_to = "diag") %>%
    filter(!is.na(diag), diag != "") %>%
    select(pat_id, diag)
  
  inpat_diags <- inpat %>%
    select(pat_id, from_dt, diag_admit, diag1, diag2, diag3, diag4, diag5, 
           diag6, diag7, diag8, diag9, diag10, diag11, diag12) %>%
    pivot_longer(cols = c(diag_admit, starts_with("diag")), values_to = "diag") %>%
    filter(!is.na(diag), diag != "") %>%
    select(pat_id, diag)
  
  all_diags <- bind_rows(outpat_diags, inpat_diags) %>% distinct()
  
  cat(sprintf("  %d: Total diagnosis codes: %s from %s patients\n",
              year, 
              format(nrow(all_diags), big.mark = ","),
              format(n_distinct(all_diags$pat_id), big.mark = ",")))
  
  # Charlson comorbidities
  comorbid <- all_diags %>%
    mutate(
      mi = str_detect(diag, "^(I21|I22|I252)"),
      chf = str_detect(diag, "^(I099|I110|I130|I132|I255|I420|I425|I426|I427|I428|I429|I43|I50|P290)"),
      pvd = str_detect(diag, "^(I70|I71|I731|I738|I739|I771|I790|I792|K551|K558|K559|Z958|Z959)"),
      cevd = str_detect(diag, "^(G45|G46|H340|I60|I61|I62|I63|I64|I65|I66|I67|I68|I69)"),
      dementia = str_detect(diag, "^(F00|F01|F02|F03|F051|G30|G311)"),
      copd = str_detect(diag, "^(I278|I279|J40|J41|J42|J43|J44|J45|J46|J47|J60|J61|J62|J63|J64|J65|J66|J67|J684|J701|J703)"),
      rheum = str_detect(diag, "^(M05|M06|M315|M32|M33|M34|M351|M353|M360)"),
      pud = str_detect(diag, "^(K25|K26|K27|K28)"),
      mild_liver = str_detect(diag, "^(B18|K700|K701|K702|K703|K709|K713|K714|K715|K717|K73|K74|K760|K762|K763|K764|K768|K769|Z944)"),
      hemiplegia = str_detect(diag, "^(G041|G114|G801|G802|G81|G82|G830|G831|G832|G833|G834|G839)"),
      ckd = str_detect(diag, "^(I120|I131|N032|N033|N034|N035|N036|N037|N052|N053|N054|N055|N056|N057|N18|N19|N250|Z490|Z491|Z492|Z940|Z992)"),
      cancer = str_detect(diag, "^C0[0-2]|^C3[0-4]|^C37|^C38|^C39|^C4[0-1]|^C43|^C4[5-9]|^C5[0-8]|^C6[0-9]|^C7[0-6]|^C81|^C82|^C83|^C84|^C85|^C88|^C90|^C91|^C92|^C93|^C94|^C95|^C96|^C97"),
      severe_liver = str_detect(diag, "^(I850|I859|I864|I982|K704|K711|K721|K729|K765|K766|K767)"),
      metastatic = str_detect(diag, "^(C7[7-9]|C80)"),
      hiv = str_detect(diag, "^(B20|B21|B22|B24)")
    ) %>%
    group_by(pat_id) %>%
    summarise(across(mi:hiv, ~as.integer(any(.))), .groups = "drop") %>%
    mutate(
      n_charlson = mi + chf + pvd + cevd + dementia + copd + rheum + pud + 
        mild_liver + hemiplegia + ckd + cancer + severe_liver + 
        metastatic + hiv
    )
  
  # DCSI
  dcsi <- all_diags %>%
    mutate(
      retinopathy = case_when(
        str_detect(diag, "^(E1032|E1033|E1134|E1135|E1136|E1137|E1139)") ~ 2,
        str_detect(diag, "^(H360)") ~ 1,
        TRUE ~ 0
      ),
      nephropathy = case_when(
        str_detect(diag, "^(E1022|E1122|N18)") ~ 2,
        str_detect(diag, "^(E1021|E1121)") ~ 1,
        TRUE ~ 0
      ),
      neuropathy = case_when(
        str_detect(diag, "^(E1042|E1043|E1044|E1142|E1143|E1144)") ~ 2,
        str_detect(diag, "^(E1040|E1041|E1140|E1141)") ~ 1,
        TRUE ~ 0
      ),
      cerebrovascular = case_when(
        str_detect(diag, "^(I6[0-9])") ~ 2,
        str_detect(diag, "^(G45)") ~ 1,
        TRUE ~ 0
      ),
      cardiovascular = case_when(
        str_detect(diag, "^(I2[1-2]|I50)") ~ 2,
        str_detect(diag, "^(I20|I24|I25)") ~ 1,
        TRUE ~ 0
      ),
      pvd_dcsi = case_when(
        str_detect(diag, "^(E1051|E1052|E1151|E1152|I70|I73)") ~ 2,
        str_detect(diag, "^(E1059|E1159)") ~ 1,
        TRUE ~ 0
      ),
      metabolic = case_when(
        str_detect(diag, "^(E1010|E1110|E1165|E1169)") ~ 2,
        str_detect(diag, "^(E1000|E1001|E1100|E1101)") ~ 1,
        TRUE ~ 0
      )
    ) %>%
    group_by(pat_id) %>%
    summarise(
      retinopathy = max(retinopathy),
      nephropathy = max(nephropathy),
      neuropathy = max(neuropathy),
      cerebrovascular = max(cerebrovascular),
      cardiovascular = max(cardiovascular),
      pvd_dcsi = max(pvd_dcsi),
      metabolic = max(metabolic),
      .groups = "drop"
    ) %>%
    mutate(
      dcsi_total = retinopathy + nephropathy + neuropathy + cerebrovascular + 
        cardiovascular + pvd_dcsi + metabolic
    )
  
  comorbid_final <- comorbid %>%
    left_join(dcsi, by = "pat_id") %>%
    mutate(across(starts_with("retinopathy"):dcsi_total, ~replace_na(., 0)))
  
  return(comorbid_final)
}

comorbid_2022 <- calculate_comorbidities(outpat_2022, inpat_2022, 2022)
comorbid_2023 <- calculate_comorbidities(outpat_2023, inpat_2023, 2023)

track_step(5.1, "Comorbidities calculated for 2022",
           data = comorbid_2022)

track_step(5.2, "Comorbidities calculated for 2023",
           data = comorbid_2023)

# -----------------------------------------------------------------------------
# 6. CREATE MONTHLY ANALYTICAL DATASETS
# -----------------------------------------------------------------------------

cat("\n--- STEP 6: CREATING MONTHLY ANALYTICAL DATASETS ---\n")

create_monthly_analytical <- function(cohort_data, cohort_flag_col, year, 
                                      age_col, treatment_col, payer_col,
                                      comorbid, restrictive = FALSE) {
  
  cohort_patients <- cohort_data %>%
    filter(!!sym(cohort_flag_col) == 1) %>%
    pull(pat_id)
  
  if (year == 2022) {
    month_grid <- expand_grid(
      pat_id = cohort_patients,
      year_month = seq(as.Date("2022-01-01"), as.Date("2022-12-01"), by = "month")
    )
    rx_month <- rx_month_2022
  } else {
    month_grid <- expand_grid(
      pat_id = cohort_patients,
      year_month = seq(as.Date("2023-01-01"), as.Date("2023-10-01"), by = "month")
    )
    rx_month <- rx_month_2023
  }
  
  complete_months <- month_grid %>%
    left_join(rx_month, by = c("pat_id", "year_month")) %>%
    mutate(across(where(is.numeric), ~replace_na(., 0)))
  
  patient_info <- cohort_data %>%
    filter(!!sym(cohort_flag_col) == 1) %>%
    select(pat_id, der_sex, pat_state, 
           age = !!sym(age_col), 
           treatment = !!sym(treatment_col),
           payer = !!sym(payer_col))
  
  analytical_monthly <- complete_months %>%
    left_join(patient_info, by = "pat_id") %>%
    left_join(comorbid, by = "pat_id") %>%
    mutate(
      year = year(year_month),
      month = month(year_month),
      post_ira = (year == 2023),
      months_since_ira = ifelse(year == 2023, month, 0),
      charlson_bin = cut(n_charlson, breaks = c(-1, 0, 1, 2, Inf),
                         labels = c("0", "1", "2", "3+"), right = TRUE),
      dcsi_bin = cut(dcsi_total, breaks = c(-1, 0, 2, 4, Inf),
                     labels = c("0", "1-2", "3-4", "5+"), right = TRUE),
      polypharm_bin = cut(n_unique_drugs_excl_insulin, 
                          breaks = c(-1, 5, 10, 15, Inf),
                          labels = c("0-5", "6-10", "11-15", "16+"), right = TRUE)
    )
  
  if (restrictive) {
    analytical_monthly <- analytical_monthly %>%
      mutate(
        age_bin = case_when(
          treatment == 1 & age == 65 ~ "65",
          treatment == 1 & age == 66 ~ "66",
          treatment == 1 & age == 67 ~ "67",
          treatment == 1 & age >= 68 ~ "68+",
          treatment == 0 & age == 62 ~ "62",
          treatment == 0 & age == 63 ~ "63",
          treatment == 0 & age == 64 ~ "64",
          treatment == 0 & age >= 65 ~ "65+",
          TRUE ~ NA_character_
        ),
        age_bin = factor(age_bin)
      )
  } else {
    analytical_monthly <- analytical_monthly %>%
      mutate(
        age_bin = case_when(
          treatment == 1 & age >= 65 & age <= 69 ~ "65-69",
          treatment == 1 & age >= 70 & age <= 74 ~ "70-74",
          treatment == 1 & age >= 75 & age <= 79 ~ "75-79",
          treatment == 1 & age >= 80 & age <= 84 ~ "80-84",
          treatment == 1 & age >= 85 ~ "85+",
          treatment == 0 & age >= 18 & age <= 29 ~ "18-29",
          treatment == 0 & age >= 30 & age <= 39 ~ "30-39",
          treatment == 0 & age >= 40 & age <= 49 ~ "40-49",
          treatment == 0 & age >= 50 & age <= 59 ~ "50-59",
          treatment == 0 & age >= 60 & age <= 64 ~ "60-64",
          treatment == 0 & age >= 65 ~ "65+",
          TRUE ~ NA_character_
        ),
        age_bin = factor(age_bin)
      )
  }
  
  return(analytical_monthly)
}

# Create cross-sectional datasets
cat("\nCreating cross-sectional 2022 monthly dataset...\n")
monthly_cross_2022 <- create_monthly_analytical(
  cohort_cross, "in_cohort_2022", 2022, "age_2022", "treatment_2022", "payer_2022",
  comorbid_2022, restrictive = FALSE
)

track_step(6.1, "Monthly analytical cross-sectional 2022",
           data = monthly_cross_2022)

cat("\nCreating cross-sectional 2023 monthly dataset...\n")
monthly_cross_2023 <- create_monthly_analytical(
  cohort_cross, "in_cohort_2023", 2023, "age_2023", "treatment_2023", "payer_2023",
  comorbid_2023, restrictive = FALSE
)

track_step(6.2, "Monthly analytical cross-sectional 2023",
           data = monthly_cross_2023)

# Create restrictive datasets
cat("\nCreating restrictive 2022 monthly dataset (62-64 vs 65-67)...\n")
monthly_restrictive_2022 <- cohort_cross %>%
  filter(in_cohort_2022 == 1,
         ((treatment_2022 == 0 & age_2022 >= 62 & age_2022 <= 64) |
            (treatment_2022 == 1 & age_2022 >= 65 & age_2022 <= 67))) %>%
  {create_monthly_analytical(., "in_cohort_2022", 2022, "age_2022", 
                             "treatment_2022", "payer_2022", comorbid_2022, 
                             restrictive = TRUE)}

track_step(6.3, "Monthly analytical restrictive 2022 (62-64 vs 65-67)",
           data = monthly_restrictive_2022)

cat("\nCreating restrictive 2023 monthly dataset (62-64 vs 65-67)...\n")
monthly_restrictive_2023 <- cohort_cross %>%
  filter(in_cohort_2023 == 1,
         ((treatment_2023 == 0 & age_2023 >= 62 & age_2023 <= 64) |
            (treatment_2023 == 1 & age_2023 >= 65 & age_2023 <= 67))) %>%
  {create_monthly_analytical(., "in_cohort_2023", 2023, "age_2023", 
                             "treatment_2023", "payer_2023", comorbid_2023, 
                             restrictive = TRUE)}

track_step(6.4, "Monthly analytical restrictive 2023 (62-64 vs 65-67)",
           data = monthly_restrictive_2023)

# Create panel datasets
cat("\nCreating panel monthly dataset (22 months)...\n")
panel_patients <- cohort_panel %>% filter(in_cohort_panel == 1) %>% pull(pat_id)

monthly_panel <- bind_rows(
  monthly_cross_2022 %>% filter(pat_id %in% panel_patients),
  monthly_cross_2023 %>% filter(pat_id %in% panel_patients)
) %>%
  arrange(pat_id, year_month)

track_step(6.5, "Monthly analytical panel (22 months: 12 in 2022 + 10 in 2023)",
           data = monthly_panel)

cat("\nCreating balanced panel monthly dataset (20 months)...\n")
monthly_panel_balanced <- bind_rows(
  monthly_cross_2022 %>% 
    filter(pat_id %in% panel_patients,
           month(year_month) <= 10),
  monthly_cross_2023 %>% 
    filter(pat_id %in% panel_patients)
) %>%
  arrange(pat_id, year_month)

track_step(6.6, "Monthly analytical balanced panel (20 months: Jan-Oct both years)",
           data = monthly_panel_balanced)

# Additional restrictive cohorts
cat("\n--- CREATING ADDITIONAL RESTRICTIVE COHORTS ---\n")

# 54-64 vs 65-75
cat("\nCohort: 54-64 vs 65-75...\n")
cohort_cross_54_64_vs_65_75_2022 <- cohort_cross %>%
  filter(in_cohort_2022 == 1,
         ((treatment_2022 == 0 & age_2022 >= 54 & age_2022 <= 64) |
            (treatment_2022 == 1 & age_2022 >= 65 & age_2022 <= 75)))

monthly_restrictive_54_64_vs_65_75_2022 <- create_monthly_analytical(
  cohort_cross_54_64_vs_65_75_2022, "in_cohort_2022", 2022, "age_2022", 
  "treatment_2022", "payer_2022", comorbid_2022, restrictive = FALSE
)

track_step(6.7, "Restrictive 54-64 vs 65-75 (2022)",
           data = monthly_restrictive_54_64_vs_65_75_2022)

cohort_cross_54_64_vs_65_75_2023 <- cohort_cross %>%
  filter(in_cohort_2023 == 1,
         ((treatment_2023 == 0 & age_2023 >= 54 & age_2023 <= 64) |
            (treatment_2023 == 1 & age_2023 >= 65 & age_2023 <= 75)))

monthly_restrictive_54_64_vs_65_75_2023 <- create_monthly_analytical(
  cohort_cross_54_64_vs_65_75_2023, "in_cohort_2023", 2023, "age_2023", 
  "treatment_2023", "payer_2023", comorbid_2023, restrictive = FALSE
)

track_step(6.8, "Restrictive 54-64 vs 65-75 (2023)",
           data = monthly_restrictive_54_64_vs_65_75_2023)

# 54-64 vs 65+
cat("\nCohort: 54-64 vs 65+...\n")
cohort_cross_54_64_vs_65_plus_2022 <- cohort_cross %>%
  filter(in_cohort_2022 == 1,
         ((treatment_2022 == 0 & age_2022 >= 54 & age_2022 <= 64) |
            (treatment_2022 == 1 & age_2022 >= 65)))

monthly_restrictive_54_64_vs_65_plus_2022 <- create_monthly_analytical(
  cohort_cross_54_64_vs_65_plus_2022, "in_cohort_2022", 2022, "age_2022", 
  "treatment_2022", "payer_2022", comorbid_2022, restrictive = FALSE
)

track_step(6.9, "Restrictive 54-64 vs 65+ (2022)",
           data = monthly_restrictive_54_64_vs_65_plus_2022)

cohort_cross_54_64_vs_65_plus_2023 <- cohort_cross %>%
  filter(in_cohort_2023 == 1,
         ((treatment_2023 == 0 & age_2023 >= 54 & age_2023 <= 64) |
            (treatment_2023 == 1 & age_2023 >= 65)))

monthly_restrictive_54_64_vs_65_plus_2023 <- create_monthly_analytical(
  cohort_cross_54_64_vs_65_plus_2023, "in_cohort_2023", 2023, "age_2023", 
  "treatment_2023", "payer_2023", comorbid_2023, restrictive = FALSE
)

track_step(6.10, "Restrictive 54-64 vs 65+ (2023)",
           data = monthly_restrictive_54_64_vs_65_plus_2023)

# 62-64 vs 65+
cat("\nCohort: 62-64 vs 65+...\n")
cohort_cross_62_64_vs_65_plus_2022 <- cohort_cross %>%
  filter(in_cohort_2022 == 1,
         ((treatment_2022 == 0 & age_2022 >= 62 & age_2022 <= 64) |
            (treatment_2022 == 1 & age_2022 >= 65)))

monthly_restrictive_62_64_vs_65_plus_2022 <- create_monthly_analytical(
  cohort_cross_62_64_vs_65_plus_2022, "in_cohort_2022", 2022, "age_2022", 
  "treatment_2022", "payer_2022", comorbid_2022, restrictive = FALSE
)

track_step(6.11, "Restrictive 62-64 vs 65+ (2022)",
           data = monthly_restrictive_62_64_vs_65_plus_2022)

cohort_cross_62_64_vs_65_plus_2023 <- cohort_cross %>%
  filter(in_cohort_2023 == 1,
         ((treatment_2023 == 0 & age_2023 >= 62 & age_2023 <= 64) |
            (treatment_2023 == 1 & age_2023 >= 65)))

monthly_restrictive_62_64_vs_65_plus_2023 <- create_monthly_analytical(
  cohort_cross_62_64_vs_65_plus_2023, "in_cohort_2023", 2023, "age_2023", 
  "treatment_2023", "payer_2023", comorbid_2023, restrictive = FALSE
)

track_step(6.12, "Restrictive 62-64 vs 65+ (2023)",
           data = monthly_restrictive_62_64_vs_65_plus_2023)

# -----------------------------------------------------------------------------
# 7. EXPORT DATASETS
# -----------------------------------------------------------------------------

cat("\n--- STEP 7: EXPORTING DATASETS ---\n")

write_csv(monthly_cross_2022, "analytical_monthly_cross_2022.csv")
write_csv(monthly_cross_2023, "analytical_monthly_cross_2023.csv")
write_csv(monthly_restrictive_2022, "analytical_monthly_restrictive_2022.csv")
write_csv(monthly_restrictive_2023, "analytical_monthly_restrictive_2023.csv")
write_csv(monthly_panel, "analytical_monthly_panel.csv")
write_csv(monthly_panel_balanced, "analytical_monthly_panel_balanced.csv")
write_csv(monthly_restrictive_54_64_vs_65_75_2022, "analytical_monthly_restrictive_54_64_vs_65_75_2022.csv")
write_csv(monthly_restrictive_54_64_vs_65_75_2023, "analytical_monthly_restrictive_54_64_vs_65_75_2023.csv")
write_csv(monthly_restrictive_54_64_vs_65_plus_2022, "analytical_monthly_restrictive_54_64_vs_65_plus_2022.csv")
write_csv(monthly_restrictive_54_64_vs_65_plus_2023, "analytical_monthly_restrictive_54_64_vs_65_plus_2023.csv")
write_csv(monthly_restrictive_62_64_vs_65_plus_2022, "analytical_monthly_restrictive_62_64_vs_65_plus_2022.csv")
write_csv(monthly_restrictive_62_64_vs_65_plus_2023, "analytical_monthly_restrictive_62_64_vs_65_plus_2023.csv")

cat("\nAll datasets exported successfully.\n")

# -----------------------------------------------------------------------------
# 8. CREATE TRACKING SUMMARY AND EXPORT
# -----------------------------------------------------------------------------

cat("\n--- STEP 8: CREATING TRACKING SUMMARY ---\n")

tracking_df <- bind_rows(tracking) %>%
  mutate(
    n_patients = format(n_patients, big.mark = ","),
    n_patient_months = format(n_patient_months, big.mark = ","),
    n_claims = format(n_claims, big.mark = ",")
  )

write_csv(tracking_df, "analytical_dataset_construction_tracking.csv")

cat("\n=============================================================================\n")
cat("MONTHLY ANALYTICAL DATASET CONSTRUCTION COMPLETE\n")
cat("=============================================================================\n\n")

print(tracking_df, n = Inf)

cat("\n\nOUTPUT FILES CREATED:\n")
cat("  1. analytical_monthly_cross_2022.csv\n")
cat("  2. analytical_monthly_cross_2023.csv\n")
cat("  3. analytical_monthly_restrictive_2022.csv (62-64 vs 65-67)\n")
cat("  4. analytical_monthly_restrictive_2023.csv (62-64 vs 65-67)\n")
cat("  5. analytical_monthly_panel.csv (22 months)\n")
cat("  6. analytical_monthly_panel_balanced.csv (20 months)\n")
cat("  7-12. Additional restrictive cohort files\n")
cat("  13. analytical_dataset_construction_tracking.csv (THIS FILE FOR WRITE-UP)\n\n")

# Final summary statistics
cat("=== FINAL SUMMARY ===\n\n")

summary_stats <- tribble(
  ~Cohort, ~N_Patients, ~N_Patient_Months, ~Months_Per_Patient,
  "Cross-sectional 2022", n_distinct(monthly_cross_2022$pat_id), 
  nrow(monthly_cross_2022), 12,
  "Cross-sectional 2023", n_distinct(monthly_cross_2023$pat_id), 
  nrow(monthly_cross_2023), 10,
  "Restrictive 2022 (62-64 vs 65-67)", n_distinct(monthly_restrictive_2022$pat_id),
  nrow(monthly_restrictive_2022), 12,
  "Restrictive 2023 (62-64 vs 65-67)", n_distinct(monthly_restrictive_2023$pat_id),
  nrow(monthly_restrictive_2023), 10,
  "Panel (22 months)", n_distinct(monthly_panel$pat_id),
  nrow(monthly_panel), 22,
  "Panel Balanced (20 months)", n_distinct(monthly_panel_balanced$pat_id),
  nrow(monthly_panel_balanced), 20
)

print(summary_stats)

write_csv(summary_stats, "analytical_dataset_summary.csv")

cat("\n\nKEY POINTS FOR WRITE-UP:\n")
cat("1. Patients with unknown sex (der_sex='U') excluded from all cohorts\n")
cat("2. Reversed pharmacy claims (copay=0 & quan=NA) removed\n")
cat("3. Claims deduplicated at pat_id + from_dt + ndc + prscbr_id level\n")
cat("4. All costs standardized to 30-day supply equivalents per IRA regulations\n")
cat("5. Primary outcome: insulin_oop_per_supply (OOP per 30-day supply)\n")
cat("6. Missing months filled with zeros for complete patient-month panels\n")
cat("7. 2023 data limited to January-October due to enrollment data availability\n\n")

cat("Use 'analytical_dataset_construction_tracking.csv' for detailed Methods section.\n")
cat("=============================================================================\n")