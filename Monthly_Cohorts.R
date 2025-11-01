library(tidyverse)
library(haven)
library(lubridate)

# =============================================================================
# MONTHLY ANALYTICAL DATASETS - IRA INSULIN COPAY CAP (UPDATED)
# =============================================================================
# Changes:
# 1. Added deduplication at pat_id + from_dt + ndc + prscbr_id level
# 2. Control age 18-64 (not 62-64) - though age bins still flexible
# 3. Fixed panel to be 22 months (not 24)
# 4. CRITICAL: Added 30-day supply standardization per IRA regulations
#    - IRA cap is $35 per "month's supply" = per 30 days
#    - 90-day fill = 3 Ã 30-day supplies = $105 max is compliant
#    - New primary outcomes: insulin_oop_per_supply, insulin_copay_per_supply
# 5. NEW: Exclude patients with der_sex = "U" (unknown sex)
#
# ENROLLMENT CRITERIA:
# - Cross-sectional 2022: 12 months continuous enrollment (Jan-Dec 2022)
# - Cross-sectional 2023: 10 months continuous enrollment (Jan-Oct 2023)
# - Panel 22-month: 22+ consecutive months spanning 2022-2023
# - Balanced panel: 20 months (Jan-Oct in both years)
# =============================================================================

# -----------------------------------------------------------------------------
# 1. LOAD DATA AND EXCLUDE UNKNOWN SEX
# -----------------------------------------------------------------------------

# Cohort files - FILTER OUT der_sex = "U"
cohort_cross <- read_sas("cohort_cross_sectional.sas7bdat") %>%
  filter(der_sex != "U")

cohort_panel <- read_sas("cohort_panel_22.sas7bdat") %>%
  filter(der_sex != "U")

cat("\n=== SEX EXCLUSION SUMMARY ===\n")
cat("Patients with der_sex = 'U' have been excluded from all cohorts\n")
cat("Remaining cross-sectional patients:", n_distinct(cohort_cross$pat_id), "\n")
cat("Remaining panel patients:", n_distinct(cohort_panel$pat_id), "\n\n")

# Pharmacy claims (already processed to remove reversed claims)
rx_2022 <- read_sas("master_rx_2022.sas7bdat") %>%
  filter(!(copay == 0 & is.na(quan)))

rx_2023 <- read_sas("master_rx_2023.sas7bdat") %>%
  filter(!(copay == 0 & is.na(quan)))

# Medical claims
outpat_2022 <- read_sas("master_outpat_2022.sas7bdat")
outpat_2023 <- read_sas("master_outpat_2023.sas7bdat")
inpat_2022 <- read_sas("master_inpat_2022.sas7bdat")
inpat_2023 <- read_sas("master_inpat_2023.sas7bdat")

# -----------------------------------------------------------------------------
# 2. CREATE MONTH IDENTIFIER FROM DATES
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# 3. DEDUPLICATE AND AGGREGATE PHARMACY CLAIMS TO PATIENT-MONTH LEVEL
# -----------------------------------------------------------------------------

aggregate_to_patient_month <- function(rx_data, year) {
  
  # STEP 1: Deduplicate at fill level (pat_id + from_dt + ndc + prscbr_id)
  # This prevents double-counting from multiple claim lines for same fill
  rx_deduplicated <- rx_data %>%
    group_by(pat_id, from_dt, ndc, prscbr_id) %>%
    summarise(
      year_month = first(year_month),
      gpi14 = first(gpi14),
      generic_name = first(generic_name),
      product_name = first(product_name),
      copay = max(copay, na.rm = TRUE),          # Take max across lines
      deductible = max(deductible, na.rm = TRUE),
      coinsamt = max(coinsamt, na.rm = TRUE),
      cobamt = max(cobamt, na.rm = TRUE),
      paid = max(paid, na.rm = TRUE),
      allowed = max(allowed, na.rm = TRUE),
      quan = ifelse(all(is.na(quan)), 0, max(quan, na.rm = TRUE)),
      dayssup = ifelse(all(is.na(dayssup)), 30, max(dayssup, na.rm = TRUE)),  # Default to 30 if missing
      .groups = "drop"
    )
  
  # STEP 2: Add drug classifications AND 30-day standardization at fill level
  rx_classified <- rx_deduplicated %>%
    mutate(
      # Clean days supply
      dayssup_clean = case_when(
        is.na(dayssup) | dayssup == 0 ~ 30,  # Assume 30 if missing/zero
        dayssup > 365 ~ 365,                 # Cap at 1 year
        dayssup < 1 ~ 1,                     # Minimum 1 day
        TRUE ~ dayssup
      ),
      
      # CRITICAL: Standardize to 30-day supply units per IRA regulations
      standardized_supplies = pmax(1, dayssup_clean / 30),
      
      # Calculate OOP per standardized 30-day supply
      oop = copay + deductible + coinsamt + cobamt,
      copay_per_supply = copay / standardized_supplies,
      deductible_per_supply = deductible / standardized_supplies,
      coinsamt_per_supply = coinsamt / standardized_supplies,
      cobamt_per_supply = cobamt / standardized_supplies,
      oop_per_supply = oop / standardized_supplies,
      
      # Normalize names for matching
      product_upper = toupper(product_name),
      generic_upper = toupper(generic_name),
      
      # =====================================================================
      # COMPREHENSIVE INSULIN CLASSIFICATION (IRA scope)
      # =====================================================================
      
      # Primary: GPI-based identification
      is_insulin_gpi = str_detect(gpi14, "^2710"),
      
      # Secondary: Product/generic name (catches inhaled, combinations)
      is_insulin_product = str_detect(product_upper, 
                                      "INSULIN|HUMALOG|NOVOLOG|LANTUS|LEVEMIR|TRESIBA|TOUJEO|BASAGLAR|AFREZZA|SEMGLEE|FIASP|APIDRA|ADMELOG|LYUMJEV|XULTOPHY|SOLIQUA|EXUBERA"),
      
      is_insulin_generic = str_detect(generic_upper, 
                                      "INSULIN|LISPRO|ASPART|GLARGINE|DETEMIR|DEGLUDEC|GLULISINE"),
      
      # Final insulin flag (comprehensive)
      is_insulin = (is_insulin_gpi | is_insulin_product | is_insulin_generic),
      
      # Concentrated insulins (U-200, U-300, U-500)
      is_concentrated = case_when(
        !is_insulin ~ FALSE,
        str_detect(product_upper, "U-500|U500|U-300|U300|TOUJEO") ~ TRUE,
        str_detect(product_upper, "U-200|U200") ~ TRUE,
        TRUE ~ FALSE
      ),
      
      # Insulin TYPE classification (detailed)
      insulin_type = case_when(
        !is_insulin ~ NA_character_,
        
        # RAPID-ACTING
        str_detect(product_upper, "HUMALOG|ADMELOG|LISPRO") |
          str_detect(generic_upper, "LISPRO") |
          str_sub(gpi14, 1, 8) %in% c("27104005", "27104002") ~ "rapid_acting",
        
        str_detect(product_upper, "NOVOLOG|FIASP") |
          str_detect(generic_upper, "ASPART") ~ "rapid_acting",
        
        str_detect(product_upper, "APIDRA") |
          str_detect(generic_upper, "GLULISINE") ~ "rapid_acting",
        
        # SHORT-ACTING / REGULAR
        str_detect(product_upper, "HUMULIN R|NOVOLIN R|REGULAR|VELOSULIN|MYXREDLIN") |
          str_detect(generic_upper, "REGULAR") |
          str_sub(gpi14, 1, 8) == "27104010" ~ "short_acting",
        
        # INTERMEDIATE-ACTING / NPH
        str_detect(product_upper, "HUMULIN N|NOVOLIN N|NPH|INSULATARD") |
          str_detect(generic_upper, "NPH|ISOPHANE") |
          str_sub(gpi14, 1, 8) == "27104020" ~ "intermediate_acting",
        
        # LONG-ACTING
        str_detect(product_upper, "LANTUS|BASAGLAR|TOUJEO|SEMGLEE") |
          str_detect(generic_upper, "GLARGINE") |
          str_sub(gpi14, 1, 8) %in% c("27104003", "27104006") ~ "long_acting",
        
        str_detect(product_upper, "LEVEMIR") |
          str_detect(generic_upper, "DETEMIR") |
          str_sub(gpi14, 1, 8) == "27104007" ~ "long_acting",
        
        str_detect(product_upper, "TRESIBA") |
          str_detect(generic_upper, "DEGLUDEC") ~ "long_acting",
        
        # MIXED / COMBINATION
        str_detect(product_upper, "70/30|75/25|50/50|MIX") |
          str_detect(generic_upper, "PROTAMINE") |
          str_sub(gpi14, 1, 8) %in% c("27104080", "27104090") ~ "mixed",
        
        # INHALED
        str_detect(product_upper, "AFREZZA|EXUBERA") |
          str_detect(generic_upper, "INHALED") ~ "inhaled",
        
        # COMBINATION WITH GLP-1
        str_detect(product_upper, "XULTOPHY|SOLIQUA") |
          str_detect(generic_upper, "LIRAGLUTIDE|LIXISENATIDE") ~ "combination_glp1",
        
        # ANIMAL INSULINS (legacy)
        str_detect(product_upper, "ILETIN|PORK|BEEF|BOVINE") |
          str_detect(generic_upper, "PORK|BEEF|BOVINE") ~ "animal",
        
        TRUE ~ "other"
      ),
      
      # Flag if this fill's per-supply cost exceeds $35 cap
      exceeds_cap_per_supply = (oop_per_supply > 35),
      
      # Other drug classifications (falsification tests)
      is_metformin = str_detect(gpi14, "^2730"),
      is_lipid = str_detect(gpi14, "^39"),
      
      # GPI-6 for unique drug counting
      gpi6 = str_sub(gpi14, 1, 6)
    )
  
  # STEP 3: Aggregate to patient-month level with standardized supplies
  patient_month <- rx_classified %>%
    group_by(pat_id, year_month) %>%
    summarise(
      # Overall metrics
      n_unique_drugs = n_distinct(gpi6),
      total_fills = n(),
      total_copay = sum(copay, na.rm = TRUE),
      total_deductible = sum(deductible, na.rm = TRUE),
      total_coinsamt = sum(coinsamt, na.rm = TRUE),
      total_cobamt = sum(cobamt, na.rm = TRUE),
      total_oop = sum(oop, na.rm = TRUE),
      total_paid = sum(paid, na.rm = TRUE),
      total_allowed = sum(allowed, na.rm = TRUE),
      
      # Standardized supply metrics (NEW - accounts for 90-day fills)
      total_standardized_supplies = sum(standardized_supplies, na.rm = TRUE),
      insulin_standardized_supplies = sum(standardized_supplies[is_insulin], na.rm = TRUE),
      
      # Insulin metrics (total costs in month)
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
      
      # CRITICAL NEW: Per-supply costs (standardized to 30-day equivalents)
      # This is what the IRA $35 cap actually applies to
      insulin_copay_per_supply = weighted.mean(copay_per_supply[is_insulin], 
                                               standardized_supplies[is_insulin], na.rm = TRUE),
      insulin_oop_per_supply = weighted.mean(oop_per_supply[is_insulin], 
                                             standardized_supplies[is_insulin], na.rm = TRUE),
      
      # Cap compliance metrics (CORRECTED)
      pct_insulin_supplies_exceeding_cap = mean(exceeds_cap_per_supply[is_insulin], na.rm = TRUE),
      any_insulin_supply_exceeds_cap = any(exceeds_cap_per_supply[is_insulin], na.rm = TRUE),
      
      # Adherence metrics (NEW)
      # Count of concentrated insulin fills
      n_concentrated_fills = sum(is_insulin & is_concentrated, na.rm = TRUE),
      
      # Insulin by type (total costs)
      n_fills_rapid_acting = sum(is_insulin & insulin_type == "rapid_acting", na.rm = TRUE),
      n_fills_short_acting = sum(is_insulin & insulin_type == "short_acting", na.rm = TRUE),
      n_fills_intermediate_acting = sum(is_insulin & insulin_type == "intermediate_acting", na.rm = TRUE),
      n_fills_long_acting = sum(is_insulin & insulin_type == "long_acting", na.rm = TRUE),
      n_fills_mixed = sum(is_insulin & insulin_type == "mixed", na.rm = TRUE),
      n_fills_inhaled = sum(is_insulin & insulin_type == "inhaled", na.rm = TRUE),
      n_fills_combination = sum(is_insulin & insulin_type == "combination_glp1", na.rm = TRUE),
      
      supplies_rapid_acting = sum(standardized_supplies[is_insulin & insulin_type == "rapid_acting"], na.rm = TRUE),
      supplies_short_acting = sum(standardized_supplies[is_insulin & insulin_type == "short_acting"], na.rm = TRUE),
      supplies_intermediate_acting = sum(standardized_supplies[is_insulin & insulin_type == "intermediate_acting"], na.rm = TRUE),
      supplies_long_acting = sum(standardized_supplies[is_insulin & insulin_type == "long_acting"], na.rm = TRUE),
      supplies_mixed = sum(standardized_supplies[is_insulin & insulin_type == "mixed"], na.rm = TRUE),
      supplies_inhaled = sum(standardized_supplies[is_insulin & insulin_type == "inhaled"], na.rm = TRUE),
      supplies_combination = sum(standardized_supplies[is_insulin & insulin_type == "combination_glp1"], na.rm = TRUE),
      
      copay_rapid_acting = sum(copay[is_insulin & insulin_type == "rapid_acting"], na.rm = TRUE),
      copay_short_acting = sum(copay[is_insulin & insulin_type == "short_acting"], na.rm = TRUE),
      copay_intermediate_acting = sum(copay[is_insulin & insulin_type == "intermediate_acting"], na.rm = TRUE),
      copay_long_acting = sum(copay[is_insulin & insulin_type == "long_acting"], na.rm = TRUE),
      copay_mixed = sum(copay[is_insulin & insulin_type == "mixed"], na.rm = TRUE),
      
      oop_rapid_acting = sum(oop[is_insulin & insulin_type == "rapid_acting"], na.rm = TRUE),
      oop_short_acting = sum(oop[is_insulin & insulin_type == "short_acting"], na.rm = TRUE),
      oop_intermediate_acting = sum(oop[is_insulin & insulin_type == "intermediate_acting"], na.rm = TRUE),
      oop_long_acting = sum(oop[is_insulin & insulin_type == "long_acting"], na.rm = TRUE),
      oop_mixed = sum(oop[is_insulin & insulin_type == "mixed"], na.rm = TRUE),
      
      total_quan_rapid_acting = sum(quan[is_insulin & insulin_type == "rapid_acting"], na.rm = TRUE),
      total_quan_short_acting = sum(quan[is_insulin & insulin_type == "short_acting"], na.rm = TRUE),
      total_quan_intermediate_acting = sum(quan[is_insulin & insulin_type == "intermediate_acting"], na.rm = TRUE),
      total_quan_long_acting = sum(quan[is_insulin & insulin_type == "long_acting"], na.rm = TRUE),
      total_quan_mixed = sum(quan[is_insulin & insulin_type == "mixed"], na.rm = TRUE),
      
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
      # CORRECTED: Flag if AVERAGE per-supply OOP in this month exceeds $35
      # This is the proper interpretation of the IRA cap
      exceeds_cap = (insulin_oop_per_supply > 35) & !is.na(insulin_oop_per_supply),
      
      # Replace NaN with 0 for per-supply metrics (occurs when no insulin fills)
      insulin_copay_per_supply = ifelse(is.nan(insulin_copay_per_supply), 0, insulin_copay_per_supply),
      insulin_oop_per_supply = ifelse(is.nan(insulin_oop_per_supply), 0, insulin_oop_per_supply),
      pct_insulin_supplies_exceeding_cap = ifelse(is.nan(pct_insulin_supplies_exceeding_cap), 0, pct_insulin_supplies_exceeding_cap),
      
      # ADHERENCE METRICS (NEW)
      # 1. Binary: Did patient fill ANY insulin this month?
      insulin_gap = (insulin_n_fills == 0),
      
      # 2. Binary: Did patient have enough insulin for full month?
      fully_covered = (insulin_standardized_supplies >= 1.0) & (insulin_n_fills > 0),
      
      # 3. Estimated days covered (capped at 30)
      estimated_days_covered = pmin(insulin_standardized_supplies * 30, 30),
      
      # 4. Supply adequacy ratio (>1 = oversupply/stockpiling, <1 = undersupply)
      supply_adequacy_ratio = insulin_standardized_supplies / 1.0,
      
      # Exclude non-insulin drugs from unique drug count
      n_unique_drugs_excl_insulin = n_unique_drugs - (insulin_n_fills > 0)
    )
  
  return(patient_month)
}

rx_month_2022 <- aggregate_to_patient_month(rx_2022, 2022)
rx_month_2023 <- aggregate_to_patient_month(rx_2023, 2023)

# Validation: Check insulin classification
cat("\n=== INSULIN CLASSIFICATION VALIDATION (2022) ===\n")

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

cat("Classification method overlap:\n")
cat("  GPI only (27104):       ", insulin_validation_2022$gpi_only, "claims\n")
cat("  Name only (no GPI):     ", insulin_validation_2022$name_only, "claims\n")
cat("  Both methods agree:     ", insulin_validation_2022$both, "claims\n")
cat("  Total insulin claims:   ", insulin_validation_2022$total_insulin, "claims\n")

if (insulin_validation_2022$name_only > 0) {
  cat("\nWARNING:", insulin_validation_2022$name_only, 
      "insulin claims captured by name but missed by GPI.\n")
  cat("These may be: inhaled insulin (Afrezza), GLP-1 combos (Xultophy/Soliqua),\n")
  cat("or other products with non-standard GPI codes.\n")
}

cat("\n")

# -----------------------------------------------------------------------------
# 4. CREATE COMPLETE PATIENT-MONTH PANEL (FILL IN MISSING MONTHS)
# -----------------------------------------------------------------------------

create_complete_patient_months <- function(cohort_data, cohort_flag_col, rx_month_2022, rx_month_2023, year) {
  # Get patients in cohort
  cohort_patients <- cohort_data %>%
    filter(!!sym(cohort_flag_col) == 1) %>%
    pull(pat_id)
  
  # Create complete month grid for year
  # NOTE: 2023 enrollment data only available Jan-Oct
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
  
  # Join and fill missing months with zeros
  complete_months <- month_grid %>%
    left_join(rx_month, by = c("pat_id", "year_month")) %>%
    mutate(across(where(is.numeric), ~replace_na(., 0)))
  
  return(complete_months)
}

complete_2022 <- create_complete_patient_months(cohort_cross, "in_cohort_2022", 
                                                rx_month_2022, rx_month_2023, 2022)
complete_2023 <- create_complete_patient_months(cohort_cross, "in_cohort_2023", 
                                                rx_month_2022, rx_month_2023, 2023)

# -----------------------------------------------------------------------------
# 5. CALCULATE COMORBIDITIES (YEARLY)
# -----------------------------------------------------------------------------

calculate_comorbidities <- function(outpat, inpat) {
  
  # Extract all diagnosis codes
  # NOTE: Exclude diagprc_ind (it's a numeric indicator, not a diagnosis code)
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
  
  all_diags <- bind_rows(outpat_diags, inpat_diags) %>%
    distinct()
  
  # Charlson comorbidities (excluding diabetes)
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
      # diabetes excluded per documentation
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
  
  # DCSI (Diabetes Complications Severity Index)
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
  
  # Combine
  comorbid_final <- comorbid %>%
    left_join(dcsi, by = "pat_id") %>%
    mutate(across(starts_with("retinopathy"):dcsi_total, ~replace_na(., 0)))
  
  return(comorbid_final)
}

comorbid_2022 <- calculate_comorbidities(outpat_2022, inpat_2022)
comorbid_2023 <- calculate_comorbidities(outpat_2023, inpat_2023)

# -----------------------------------------------------------------------------
# 6. CREATE MONTHLY ANALYTICAL DATASETS
# -----------------------------------------------------------------------------

create_monthly_analytical <- function(cohort_data, cohort_flag_col, year, 
                                      age_col, treatment_col, payer_col,
                                      comorbid, restrictive = FALSE) {
  
  # CRITICAL FIX: Create complete patient-months from the FILTERED cohort_data
  # not from global complete_2022/complete_2023
  cohort_patients <- cohort_data %>%
    filter(!!sym(cohort_flag_col) == 1) %>%
    pull(pat_id)
  
  # Create month grid for THIS cohort's patients
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
  
  # Join and fill missing months with zeros
  complete_months <- month_grid %>%
    left_join(rx_month, by = c("pat_id", "year_month")) %>%
    mutate(across(where(is.numeric), ~replace_na(., 0)))
  
  # Get patient info from the FILTERED cohort_data
  patient_info <- cohort_data %>%
    filter(!!sym(cohort_flag_col) == 1) %>%
    select(pat_id, der_sex, pat_state, 
           age = !!sym(age_col), 
           treatment = !!sym(treatment_col),
           payer = !!sym(payer_col))
  
  # Join everything together
  analytical_monthly <- complete_months %>%
    left_join(patient_info, by = "pat_id") %>%
    left_join(comorbid, by = "pat_id") %>%
    mutate(
      year = year(year_month),
      month = month(year_month),
      post_ira = (year == 2023),
      months_since_ira = ifelse(year == 2023, month, 0),
      # Create binned controls
      charlson_bin = cut(n_charlson, breaks = c(-1, 0, 1, 2, Inf),
                         labels = c("0", "1", "2", "3+"), right = TRUE),
      dcsi_bin = cut(dcsi_total, breaks = c(-1, 0, 2, 4, Inf),
                     labels = c("0", "1-2", "3-4", "5+"), right = TRUE),
      polypharm_bin = cut(n_unique_drugs_excl_insulin, 
                          breaks = c(-1, 5, 10, 15, Inf),
                          labels = c("0-5", "6-10", "11-15", "16+"), right = TRUE)
    )
  
  # Create age bins based on treatment and restrictive specification
  if (restrictive) {
    # Restrictive: different bins for treatment vs control
    analytical_monthly <- analytical_monthly %>%
      mutate(
        age_bin = case_when(
          # Treatment group (65-67)
          treatment == 1 & age == 65 ~ "65",
          treatment == 1 & age == 66 ~ "66",
          treatment == 1 & age == 67 ~ "67",
          treatment == 1 & age >= 68 ~ "68+",
          # Control group (62-64)
          treatment == 0 & age == 62 ~ "62",
          treatment == 0 & age == 63 ~ "63",
          treatment == 0 & age == 64 ~ "64",
          treatment == 0 & age >= 65 ~ "65+",
          TRUE ~ NA_character_
        ),
        age_bin = factor(age_bin)
      )
  } else {
    # Standard: different bins for treatment vs control
    analytical_monthly <- analytical_monthly %>%
      mutate(
        age_bin = case_when(
          # Treatment group (65+)
          treatment == 1 & age >= 65 & age <= 69 ~ "65-69",
          treatment == 1 & age >= 70 & age <= 74 ~ "70-74",
          treatment == 1 & age >= 75 & age <= 79 ~ "75-79",
          treatment == 1 & age >= 80 & age <= 84 ~ "80-84",
          treatment == 1 & age >= 85 ~ "85+",
          # Control group (18-64)
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

# Create monthly datasets for cross-sectional cohorts
monthly_cross_2022 <- create_monthly_analytical(
  cohort_cross, "in_cohort_2022", 2022, "age_2022", "treatment_2022", "payer_2022",
  comorbid_2022, restrictive = FALSE
)

monthly_cross_2023 <- create_monthly_analytical(
  cohort_cross, "in_cohort_2023", 2023, "age_2023", "treatment_2023", "payer_2023",
  comorbid_2023, restrictive = FALSE
)

# Create monthly datasets for restrictive cohorts (62-64 vs 65-67)
monthly_restrictive_2022 <- cohort_cross %>%
  filter(in_cohort_2022 == 1,
         ((treatment_2022 == 0 & age_2022 >= 62 & age_2022 <= 64) |
            (treatment_2022 == 1 & age_2022 >= 65 & age_2022 <= 67))) %>%
  {create_monthly_analytical(., "in_cohort_2022", 2022, "age_2022", 
                             "treatment_2022", "payer_2022", comorbid_2022, 
                             restrictive = TRUE)}

monthly_restrictive_2023 <- cohort_cross %>%
  filter(in_cohort_2023 == 1,
         ((treatment_2023 == 0 & age_2023 >= 62 & age_2023 <= 64) |
            (treatment_2023 == 1 & age_2023 >= 65 & age_2023 <= 67))) %>%
  {create_monthly_analytical(., "in_cohort_2023", 2023, "age_2023", 
                             "treatment_2023", "payer_2023", comorbid_2023, 
                             restrictive = TRUE)}

# Create panel dataset (22 months per patient: 12 in 2022 + 10 in 2023)
panel_patients <- cohort_panel %>% filter(in_cohort_panel == 1) %>% pull(pat_id)

monthly_panel <- bind_rows(
  monthly_cross_2022 %>% filter(pat_id %in% panel_patients),
  monthly_cross_2023 %>% filter(pat_id %in% panel_patients)
) %>%
  arrange(pat_id, year_month)

# Create BALANCED panel (Jan-Oct only for both years = 20 months)
monthly_panel_balanced <- bind_rows(
  monthly_cross_2022 %>% 
    filter(pat_id %in% panel_patients,
           month(year_month) <= 10),
  monthly_cross_2023 %>% 
    filter(pat_id %in% panel_patients)
) %>%
  arrange(pat_id, year_month)

# -----------------------------------------------------------------------------
# 7. EXPORT MONTHLY DATASETS
# -----------------------------------------------------------------------------

write_csv(monthly_cross_2022, "analytical_monthly_cross_2022.csv")
write_csv(monthly_cross_2023, "analytical_monthly_cross_2023.csv")
write_csv(monthly_restrictive_2022, "analytical_monthly_restrictive_2022.csv")
write_csv(monthly_restrictive_2023, "analytical_monthly_restrictive_2023.csv")
write_csv(monthly_panel, "analytical_monthly_panel.csv")
write_csv(monthly_panel_balanced, "analytical_monthly_panel_balanced.csv")

# -----------------------------------------------------------------------------
# 8. SUMMARY STATISTICS
# -----------------------------------------------------------------------------

cat("\n=== MONTHLY ANALYTICAL DATASETS CREATED ===\n\n")
cat("NOTE: Patients with der_sex = 'U' excluded from all cohorts\n")
cat("NOTE: 2023 data includes January-October only (enrollment data limitation)\n\n")

cat("Cross-sectional 2022:", nrow(monthly_cross_2022), "patient-months\n")
cat("  - Unique patients:", n_distinct(monthly_cross_2022$pat_id), "\n")
cat("  - Months: 12 (Jan-Dec 2022)\n\n")

cat("Cross-sectional 2023:", nrow(monthly_cross_2023), "patient-months\n")
cat("  - Unique patients:", n_distinct(monthly_cross_2023$pat_id), "\n")
cat("  - Months: 10 (Jan-Oct 2023)\n\n")

cat("Panel (2022-2023):", nrow(monthly_panel), "patient-months\n")
cat("  - Unique patients:", n_distinct(monthly_panel$pat_id), "\n")
cat("  - Months: 22 per patient (12 in 2022 + 10 in 2023)\n\n")

cat("Panel BALANCED (Jan-Oct both years):", nrow(monthly_panel_balanced), "patient-months\n")
cat("  - Unique patients:", n_distinct(monthly_panel_balanced$pat_id), "\n")
cat("  - Months: 20 per patient (10 in 2022 + 10 in 2023)\n\n")

cat("Restrictive 2022:", nrow(monthly_restrictive_2022), "patient-months\n")
cat("Restrictive 2023:", nrow(monthly_restrictive_2023), "patient-months\n\n")

# Check months exceeding $35 cap (CORRECTED with standardization)
cap_summary <- bind_rows(
  monthly_cross_2022 %>% mutate(cohort = "2022"),
  monthly_cross_2023 %>% mutate(cohort = "2023")
) %>%
  filter(insulin_n_fills > 0) %>%  # Only months with insulin fills
  group_by(cohort, treatment) %>%
  summarise(
    n_patient_months = n(),
    # OLD measure (for comparison)
    pct_exceeding_cap_old = mean(insulin_copay > 35) * 100,
    # NEW measure (CORRECTED - per 30-day supply)
    pct_exceeding_cap_new = mean(exceeds_cap, na.rm = TRUE) * 100,
    mean_insulin_copay = mean(insulin_copay),
    median_insulin_copay = median(insulin_copay),
    # NEW: Per-supply metrics
    mean_oop_per_supply = mean(insulin_oop_per_supply, na.rm = TRUE),
    median_oop_per_supply = median(insulin_oop_per_supply, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== INSULIN COPAY CAP COMPLIANCE (CORRECTED WITH 30-DAY STANDARDIZATION) ===\n")
cat("(Among patient-months with insulin fills)\n\n")
print(cap_summary)

# Diagnostic: Check days supply distribution
days_supply_summary <- bind_rows(
  rx_month_2022 %>% filter(insulin_n_fills > 0) %>% mutate(year = 2022),
  rx_month_2023 %>% filter(insulin_n_fills > 0) %>% mutate(year = 2023)
) %>%
  group_by(year) %>%
  summarise(
    mean_supplies_per_month = mean(insulin_standardized_supplies, na.rm = TRUE),
    pct_1_supply = mean(insulin_standardized_supplies <= 1.2, na.rm = TRUE) * 100,
    pct_2_supply = mean(insulin_standardized_supplies > 1.2 & insulin_standardized_supplies <= 2.5, na.rm = TRUE) * 100,
    pct_3_supply = mean(insulin_standardized_supplies > 2.5, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cat("\n\n=== 30-DAY SUPPLY STANDARDIZATION DIAGNOSTIC ===\n")
cat("Distribution of standardized insulin supplies per patient-month:\n\n")
print(days_supply_summary)
cat("\nInterpretation:\n")
cat("- 1 supply = 30-day fills (most common)\n")
cat("- 2 supplies = 60-day fills\n")
cat("- 3+ supplies = 90-day fills\n")
cat("Note: IRA cap is $35 PER supply, so 90-day fill legally costs up to $105\n")

cat("\n\nOUTPUT FILES:\n")
cat("- analytical_monthly_cross_2022.csv (12 months)\n")
cat("- analytical_monthly_cross_2023.csv (10 months)\n")
cat("- analytical_monthly_restrictive_2022.csv (12 months)\n")
cat("- analytical_monthly_restrictive_2023.csv (10 months)\n")
cat("- analytical_monthly_panel.csv (22 months: Jan-Dec 2022 + Jan-Oct 2023)\n")
cat("- analytical_monthly_panel_balanced.csv (20 months: Jan-Oct 2022 + Jan-Oct 2023)\n")

cat("\n\nKEY VARIABLES:\n")
cat("PRIMARY OUTCOMES (NEW - 30-day standardized):\n")
cat("- insulin_oop_per_supply: Average OOP per 30-day supply (PRIMARY for cap analysis)\n")
cat("- insulin_copay_per_supply: Average copay per 30-day supply\n")
cat("- exceeds_cap: Binary flag if avg per-supply OOP > $35 (CORRECTED definition)\n")
cat("- pct_insulin_supplies_exceeding_cap: % of supplies with OOP > $35\n")
cat("- insulin_standardized_supplies: Total 30-day equivalents filled in month\n\n")
cat("SECONDARY OUTCOMES (total monthly costs):\n")
cat("- year_month: Month identifier (first day of month)\n")
cat("- post_ira: Binary indicator (1 = 2023)\n")
cat("- months_since_ira: Months since January 2023 (0 for 2022)\n")
cat("- insulin_copay, insulin_oop, insulin_n_fills: Monthly totals\n")
cat("- All insulin type-specific metrics at monthly level\n")

cat("\n\nKEY UPDATES FROM PREVIOUS VERSION:\n")
cat("1. Deduplication at pat_id + from_dt + ndc + prscbr_id level\n")
cat("2. Panel correctly shows 22 months (not 24)\n")
cat("3. Control age range 18-64 reflected in age bins\n")
cat("4. CRITICAL: 30-day supply standardization per IRA regulations\n")
cat("   - Properly accounts for 90-day fills (3 Ã $35 = $105 max)\n")
cat("   - New primary outcomes: *_per_supply metrics\n")
cat("   - Corrected exceeds_cap definition\n")
cat("5. NEW: Excluded patients with der_sex = 'U' (unknown sex)\n")
cat("\n\nIMPORTANT: Use insulin_oop_per_supply as PRIMARY outcome for IRA analysis.\n")
cat("This properly measures compliance with the $35 per-month's-supply cap.\n")

# =============================================================================
# ADDITIONAL RESTRICTIVE COHORTS (CORRECTED)
# =============================================================================

cat("\n=== CREATING ADDITIONAL RESTRICTIVE COHORTS (CORRECTED) ===\n\n")

# -----------------------------------------------------------------------------
# Cohort 1: 62-64 Commercial vs 65-67 Medicare Advantage
# (Already created above as monthly_restrictive_2022/2023)
# -----------------------------------------------------------------------------

cat("Cohort 1 (62-64 vs 65-67): Already created as monthly_restrictive files\n")

# -----------------------------------------------------------------------------
# Cohort 2: 54-64 Commercial vs 65-75 Medicare Advantage
# -----------------------------------------------------------------------------

cat("\nCohort 2: 54-64 Commercial vs 65-75 Medicare Advantage\n")

# CRITICAL FIX: Create filtered cohort subset FIRST
cohort_cross_54_64_vs_65_75_2022 <- cohort_cross %>%
  filter(in_cohort_2022 == 1,
         ((treatment_2022 == 0 & age_2022 >= 54 & age_2022 <= 64) |
            (treatment_2022 == 1 & age_2022 >= 65 & age_2022 <= 75)))

monthly_restrictive_54_64_vs_65_75_2022 <- create_monthly_analytical(
  cohort_cross_54_64_vs_65_75_2022, "in_cohort_2022", 2022, "age_2022", 
  "treatment_2022", "payer_2022", comorbid_2022, restrictive = FALSE
)

cohort_cross_54_64_vs_65_75_2023 <- cohort_cross %>%
  filter(in_cohort_2023 == 1,
         ((treatment_2023 == 0 & age_2023 >= 54 & age_2023 <= 64) |
            (treatment_2023 == 1 & age_2023 >= 65 & age_2023 <= 75)))

monthly_restrictive_54_64_vs_65_75_2023 <- create_monthly_analytical(
  cohort_cross_54_64_vs_65_75_2023, "in_cohort_2023", 2023, "age_2023", 
  "treatment_2023", "payer_2023", comorbid_2023, restrictive = FALSE
)

write_csv(monthly_restrictive_54_64_vs_65_75_2022, "analytical_monthly_restrictive_54_64_vs_65_75_2022.csv")
write_csv(monthly_restrictive_54_64_vs_65_75_2023, "analytical_monthly_restrictive_54_64_vs_65_75_2023.csv")

cat("  2022:", nrow(monthly_restrictive_54_64_vs_65_75_2022), "patient-months,",
    n_distinct(monthly_restrictive_54_64_vs_65_75_2022$pat_id), "unique patients\n")
cat("  2023:", nrow(monthly_restrictive_54_64_vs_65_75_2023), "patient-months,",
    n_distinct(monthly_restrictive_54_64_vs_65_75_2023$pat_id), "unique patients\n")

# -----------------------------------------------------------------------------
# Cohort 3: 54-64 Commercial vs 65+ Medicare Advantage
# -----------------------------------------------------------------------------

cat("\nCohort 3: 54-64 Commercial vs 65+ Medicare Advantage\n")

cohort_cross_54_64_vs_65_plus_2022 <- cohort_cross %>%
  filter(in_cohort_2022 == 1,
         ((treatment_2022 == 0 & age_2022 >= 54 & age_2022 <= 64) |
            (treatment_2022 == 1 & age_2022 >= 65)))

monthly_restrictive_54_64_vs_65_plus_2022 <- create_monthly_analytical(
  cohort_cross_54_64_vs_65_plus_2022, "in_cohort_2022", 2022, "age_2022", 
  "treatment_2022", "payer_2022", comorbid_2022, restrictive = FALSE
)

cohort_cross_54_64_vs_65_plus_2023 <- cohort_cross %>%
  filter(in_cohort_2023 == 1,
         ((treatment_2023 == 0 & age_2023 >= 54 & age_2023 <= 64) |
            (treatment_2023 == 1 & age_2023 >= 65)))

monthly_restrictive_54_64_vs_65_plus_2023 <- create_monthly_analytical(
  cohort_cross_54_64_vs_65_plus_2023, "in_cohort_2023", 2023, "age_2023", 
  "treatment_2023", "payer_2023", comorbid_2023, restrictive = FALSE
)

write_csv(monthly_restrictive_54_64_vs_65_plus_2022, "analytical_monthly_restrictive_54_64_vs_65_plus_2022.csv")
write_csv(monthly_restrictive_54_64_vs_65_plus_2023, "analytical_monthly_restrictive_54_64_vs_65_plus_2023.csv")

cat("  2022:", nrow(monthly_restrictive_54_64_vs_65_plus_2022), "patient-months,",
    n_distinct(monthly_restrictive_54_64_vs_65_plus_2022$pat_id), "unique patients\n")
cat("  2023:", nrow(monthly_restrictive_54_64_vs_65_plus_2023), "patient-months,",
    n_distinct(monthly_restrictive_54_64_vs_65_plus_2023$pat_id), "unique patients\n")

# -----------------------------------------------------------------------------
# Cohort 4: 62-64 Commercial vs 65+ Medicare Advantage
# -----------------------------------------------------------------------------

cat("\nCohort 4: 62-64 Commercial vs 65+ Medicare Advantage\n")

cohort_cross_62_64_vs_65_plus_2022 <- cohort_cross %>%
  filter(in_cohort_2022 == 1,
         ((treatment_2022 == 0 & age_2022 >= 62 & age_2022 <= 64) |
            (treatment_2022 == 1 & age_2022 >= 65)))

monthly_restrictive_62_64_vs_65_plus_2022 <- create_monthly_analytical(
  cohort_cross_62_64_vs_65_plus_2022, "in_cohort_2022", 2022, "age_2022", 
  "treatment_2022", "payer_2022", comorbid_2022, restrictive = FALSE
)

cohort_cross_62_64_vs_65_plus_2023 <- cohort_cross %>%
  filter(in_cohort_2023 == 1,
         ((treatment_2023 == 0 & age_2023 >= 62 & age_2023 <= 64) |
            (treatment_2023 == 1 & age_2023 >= 65)))

monthly_restrictive_62_64_vs_65_plus_2023 <- create_monthly_analytical(
  cohort_cross_62_64_vs_65_plus_2023, "in_cohort_2023", 2023, "age_2023", 
  "treatment_2023", "payer_2023", comorbid_2023, restrictive = FALSE
)

write_csv(monthly_restrictive_62_64_vs_65_plus_2022, "analytical_monthly_restrictive_62_64_vs_65_plus_2022.csv")
write_csv(monthly_restrictive_62_64_vs_65_plus_2023, "analytical_monthly_restrictive_62_64_vs_65_plus_2023.csv")

cat("  2022:", nrow(monthly_restrictive_62_64_vs_65_plus_2022), "patient-months,",
    n_distinct(monthly_restrictive_62_64_vs_65_plus_2022$pat_id), "unique patients\n")
cat("  2023:", nrow(monthly_restrictive_62_64_vs_65_plus_2023), "patient-months,",
    n_distinct(monthly_restrictive_62_64_vs_65_plus_2023$pat_id), "unique patients\n")

# -----------------------------------------------------------------------------
# Summary of All Restrictive Cohorts
# -----------------------------------------------------------------------------

cat("\n=== SUMMARY OF ALL RESTRICTIVE COHORTS ===\n\n")

cat("Files created:\n")
cat("1. analytical_monthly_restrictive_2022.csv (62-64 vs 65-67)\n")
cat("2. analytical_monthly_restrictive_2023.csv (62-64 vs 65-67)\n")
cat("3. analytical_monthly_restrictive_54_64_vs_65_75_2022.csv\n")
cat("4. analytical_monthly_restrictive_54_64_vs_65_75_2023.csv\n")
cat("5. analytical_monthly_restrictive_54_64_vs_65_plus_2022.csv\n")
cat("6. analytical_monthly_restrictive_54_64_vs_65_plus_2023.csv\n")
cat("7. analytical_monthly_restrictive_62_64_vs_65_plus_2022.csv\n")
cat("8. analytical_monthly_restrictive_62_64_vs_65_plus_2023.csv\n\n")

cat("All cohorts compare:\n")
cat("  Treatment = Medicare Advantage (age varies by cohort)\n")
cat("  Control = Commercial (age varies by cohort)\n")
cat("  All patients have T2D and insulin use in their respective year\n")
cat("  Patients with der_sex = 'U' excluded from all cohorts\n\n")