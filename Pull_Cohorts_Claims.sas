/*=============================================================================*/
  /* IRA INSULIN COPAY CAP ANALYSIS - EXTRACT ALL CLAIMS FOR ALL COHORTS       */
  /*                                                                             */
  /* Pulls ALL claims (medical + pharmacy) for patients in any cohort:          */
  /*   - Cross-sectional 2022                                                   */
  /*   - Cross-sectional 2023                                                   */
  /*   - Panel 22 months                                                        */
  /*                                                                             */
  /* Output: Master claims files with cohort membership flags                   */
  /*=============================================================================*/
    
    libname inpat "R:/IQVIA PharMetrics Plus (2024)/Claims - Medical - Inpatient";
  libname outpat "R:/IQVIA PharMetrics Plus (2024)/Claims - Medical - Outpatient";
  libname pharm "R:/IQVIA PharMetrics Plus (2024)/Claims - Retail Pharmacy";
  libname output "R:/IQVIA PharMetrics Plus (2024) Members/WatsonWilliam/New_copay_project";
  
  options compress=yes sortsize=8G sumsize=8G memsize=32G threads cpucount=6 fullstimer;
  
  %let start_time = %sysfunc(datetime());
  
  %put =============================================================================;
  %put EXTRACTING ALL CLAIMS FOR ALL COHORTS;
  %put =============================================================================;
  %put Start time: %sysfunc(putn(&start_time, datetime20.));
  %put;
  
  /*=============================================================================*/
    /* STEP 1: CREATE MASTER PATIENT LIST WITH ALL COHORT FLAGS                  */
    /*=============================================================================*/
      %put =============================================================================;
      %put STEP 1: CREATING MASTER PATIENT LIST;
      %put =============================================================================;
      
      /* Read cross-sectional cohorts */
        data work.cross_sectional;
      set output.cohort_cross_sectional;
      keep pat_id in_cohort_2022 treatment_2022 payer_2022 age_2022 age_group_2022
      in_cohort_2023 treatment_2023 payer_2023 age_2023 age_group_2023
      der_sex pat_state has_t2d_2022 has_insulin_2022 has_t2d_2023 has_insulin_2023;
      run;
      
      /* Read panel cohort */
        data work.panel;
      set output.cohort_panel_22;
      keep pat_id in_cohort_panel treatment_panel payer_panel age_2023 age_group_2023
      der_sex pat_state has_t2d_2022 has_insulin_2022;
      run;
      
      /* Merge all cohorts - full outer join to get everyone */
        proc sql;
      create table work.master_patient_list as
      select 
      coalesce(c.pat_id, p.pat_id) as pat_id,
      
      /* Demographics - take from whichever is not missing */
        coalesce(c.der_sex, p.der_sex) as der_sex,
      coalesce(c.pat_state, p.pat_state) as pat_state,
      
      /* Cross-sectional 2022 flags */
        coalesce(c.in_cohort_2022, 0) as in_cohort_2022,
      c.treatment_2022,
      c.payer_2022,
      c.age_2022,
      c.age_group_2022,
      
      /* Cross-sectional 2023 flags */
        coalesce(c.in_cohort_2023, 0) as in_cohort_2023,
      c.treatment_2023,
      c.payer_2023,
      
      /* Panel flags */
        coalesce(p.in_cohort_panel, 0) as in_cohort_panel,
      p.treatment_panel,
      p.payer_panel,
      
      /* Age 2023 - prefer cross-sectional if available */
        coalesce(c.age_2023, p.age_2023) as age_2023,
      coalesce(c.age_group_2023, p.age_group_2023) as age_group_2023,
      
      /* Clinical flags */
        coalesce(c.has_t2d_2022, p.has_t2d_2022, 0) as has_t2d_2022,
      coalesce(c.has_insulin_2022, p.has_insulin_2022, 0) as has_insulin_2022,
      coalesce(c.has_t2d_2023, 0) as has_t2d_2023,
      coalesce(c.has_insulin_2023, 0) as has_insulin_2023,
      
      /* Flag for which cohorts patient is in */
        case when c.in_cohort_2022=1 then 1 else 0 end as any_2022,
      case when c.in_cohort_2023=1 then 1 else 0 end as any_2023,
      case when p.in_cohort_panel=1 then 1 else 0 end as any_panel
      
      from work.cross_sectional as c
      full join work.panel as p
      on c.pat_id = p.pat_id;
      quit;
      
      /* Save master patient list */
        data output.master_patient_list;
      set work.master_patient_list;
      run;
      
      /* Get summary counts */
        proc sql noprint;
      select count(distinct pat_id) into :n_total from output.master_patient_list;
      select count(distinct pat_id) into :n_2022 from output.master_patient_list where in_cohort_2022=1;
      select count(distinct pat_id) into :n_2023 from output.master_patient_list where in_cohort_2023=1;
      select count(distinct pat_id) into :n_panel from output.master_patient_list where in_cohort_panel=1;
      quit;
      
      %put;
      %put NOTE: Master patient list created;
      %put NOTE: Total unique patients across all cohorts: &n_total;
      %put NOTE:   - In 2022 cohort: &n_2022;
      %put NOTE:   - In 2023 cohort: &n_2023;
      %put NOTE:   - In panel cohort: &n_panel;
      %put;
      
      /* Summary of cohort overlap */
        title 'COHORT OVERLAP SUMMARY';
      proc freq data=output.master_patient_list;
      tables in_cohort_2022 * in_cohort_2023 * in_cohort_panel / list nocum nopercent;
      run;
      
      /*=============================================================================*/
        /* STEP 2: EXTRACT 2022 CLAIMS                                               */
        /*=============================================================================*/
        %put =============================================================================;
      %put STEP 2: EXTRACTING 2022 CLAIMS;
      %put =============================================================================;
      
      /* Pharmacy 2022 */
        %put   Extracting 2022 pharmacy claims...;
      data output.master_rx_2022;
      if _n_ = 1 then do;
      declare hash h(dataset: "output.master_patient_list(keep=pat_id)");
      h.definekey('pat_id');
      h.definedone();
      end;
      
      set pharm.clm_rx_rtl_22;
      
      if h.check() = 0;
      
      /* Add insulin flag for convenience */
        insulin_claim = (substr(gpi14, 1, 4) = '2710');
      run;
      
      proc sql noprint;
      select count(*) format=comma15. into :rx_2022_n from output.master_rx_2022;
      select count(distinct pat_id) format=comma12. into :rx_2022_pts from output.master_rx_2022;
      quit;
      
      %put NOTE: 2022 Pharmacy - &rx_2022_n claims, &rx_2022_pts patients;
      
      /* Outpatient 2022 */
        %put   Extracting 2022 outpatient claims...;
      data output.master_outpat_2022;
      if _n_ = 1 then do;
      declare hash h(dataset: "output.master_patient_list(keep=pat_id)");
      h.definekey('pat_id');
      h.definedone();
      end;
      
      set outpat.clm_outpat_22;
      
      if h.check() = 0;
      run;
      
      proc sql noprint;
      select count(*) format=comma15. into :op_2022_n from output.master_outpat_2022;
      select count(distinct pat_id) format=comma12. into :op_2022_pts from output.master_outpat_2022;
      quit;
      
      %put NOTE: 2022 Outpatient - &op_2022_n claims, &op_2022_pts patients;
      
      /* Inpatient 2022 */
        %put   Extracting 2022 inpatient claims...;
      data output.master_inpat_2022;
      if _n_ = 1 then do;
      declare hash h(dataset: "output.master_patient_list(keep=pat_id)");
      h.definekey('pat_id');
      h.definedone();
      end;
      
      set inpat.clm_inpat_22;
      
      if h.check() = 0;
      run;
      
      proc sql noprint;
      select count(*) format=comma15. into :ip_2022_n from output.master_inpat_2022;
      select count(distinct pat_id) format=comma12. into :ip_2022_pts from output.master_inpat_2022;
      quit;
      
      %put NOTE: 2022 Inpatient - &ip_2022_n claims, &ip_2022_pts patients;
      
      /*=============================================================================*/
        /* STEP 3: EXTRACT 2023 CLAIMS                                               */
        /*=============================================================================*/
        %put =============================================================================;
      %put STEP 3: EXTRACTING 2023 CLAIMS;
      %put =============================================================================;
      
      /* Pharmacy 2023 */
        %put   Extracting 2023 pharmacy claims...;
      data output.master_rx_2023;
      if _n_ = 1 then do;
      declare hash h(dataset: "output.master_patient_list(keep=pat_id)");
      h.definekey('pat_id');
      h.definedone();
      end;
      
      set pharm.clm_rx_rtl_23;
      
      if h.check() = 0;
      
      /* Add insulin flag for convenience */
        insulin_claim = (substr(gpi14, 1, 4) = '2710');
      run;
      
      proc sql noprint;
      select count(*) format=comma15. into :rx_2023_n from output.master_rx_2023;
      select count(distinct pat_id) format=comma12. into :rx_2023_pts from output.master_rx_2023;
      quit;
      
      %put NOTE: 2023 Pharmacy - &rx_2023_n claims, &rx_2023_pts patients;
      
      /* Outpatient 2023 */
        %put   Extracting 2023 outpatient claims...;
      data output.master_outpat_2023;
      if _n_ = 1 then do;
      declare hash h(dataset: "output.master_patient_list(keep=pat_id)");
      h.definekey('pat_id');
      h.definedone();
      end;
      
      set outpat.clm_outpat_23;
      
      if h.check() = 0;
      run;
      
      proc sql noprint;
      select count(*) format=comma15. into :op_2023_n from output.master_outpat_2023;
      select count(distinct pat_id) format=comma12. into :op_2023_pts from output.master_outpat_2023;
      quit;
      
      %put NOTE: 2023 Outpatient - &op_2023_n claims, &op_2023_pts patients;
      
      /* Inpatient 2023 */
        %put   Extracting 2023 inpatient claims...;
      data output.master_inpat_2023;
      if _n_ = 1 then do;
      declare hash h(dataset: "output.master_patient_list(keep=pat_id)");
      h.definekey('pat_id');
      h.definedone();
      end;
      
      set inpat.clm_inpat_23;
      
      if h.check() = 0;
      run;
      
      proc sql noprint;
      select count(*) format=comma15. into :ip_2023_n from output.master_inpat_2023;
      select count(distinct pat_id) format=comma12. into :ip_2023_pts from output.master_inpat_2023;
      quit;
      
      %put NOTE: 2023 Inpatient - &ip_2023_n claims, &ip_2023_pts patients;
      
      /*=============================================================================*/
        /* STEP 4: CREATE SUMMARY TABLE OF EXTRACTED CLAIMS                          */
        /*=============================================================================*/
        %put =============================================================================;
      %put STEP 4: CREATING CLAIMS EXTRACTION SUMMARY;
      %put =============================================================================;
      
      data output.claims_extraction_summary;
      length claim_type $20 year 8 n_claims 8 n_patients 8;
      
      claim_type = 'Pharmacy'; year = 2022; 
      n_claims = &rx_2022_n; n_patients = &rx_2022_pts; output;
      
      claim_type = 'Outpatient'; year = 2022; 
      n_claims = &op_2022_n; n_patients = &op_2022_pts; output;
      
      claim_type = 'Inpatient'; year = 2022; 
      n_claims = &ip_2022_n; n_patients = &ip_2022_pts; output;
      
      claim_type = 'Pharmacy'; year = 2023; 
      n_claims = &rx_2023_n; n_patients = &rx_2023_pts; output;
      
      claim_type = 'Outpatient'; year = 2023; 
      n_claims = &op_2023_n; n_patients = &op_2023_pts; output;
      
      claim_type = 'Inpatient'; year = 2023; 
      n_claims = &ip_2023_n; n_patients = &ip_2023_pts; output;
      
      format n_claims comma15. n_patients comma12.;
      run;
      
      title 'CLAIMS EXTRACTION SUMMARY';
      proc print data=output.claims_extraction_summary noobs;
      run;
      
      proc export data=output.claims_extraction_summary
      outfile="R:/IQVIA PharMetrics Plus (2024) Members/WatsonWilliam/New_copay_project/claims_extraction_summary.csv"
      dbms=csv replace;
      run;
      
      /*=============================================================================*/
        /* COMPLETION SUMMARY                                                         */
        /*=============================================================================*/
        
        %let end_time = %sysfunc(datetime());
      %let elapsed = %sysevalf(&end_time - &start_time);
      
      %put;
      %put =============================================================================;
      %put CLAIMS EXTRACTION COMPLETE;
      %put =============================================================================;
      %put Elapsed time: %sysfunc(putn(&elapsed, time12.2));
      %put;
      %put =============================================================================;
      %put FILES CREATED:;
      %put =============================================================================;
      %put;
      %put MASTER PATIENT LIST:;
      %put   output.master_patient_list - All patients with cohort flags;
      %put     Variables: pat_id, der_sex, pat_state, age_2022, age_2023;
      %put                in_cohort_2022, treatment_2022, payer_2022, age_group_2022;
      %put                in_cohort_2023, treatment_2023, payer_2023, age_group_2023;
      %put                in_cohort_panel, treatment_panel, payer_panel;
      %put                has_t2d_2022, has_insulin_2022, has_t2d_2023, has_insulin_2023;
      %put;
      %put 2022 CLAIMS FILES (ALL CLAIMS FOR COHORT PATIENTS):;
      %put   output.master_rx_2022 - Pharmacy claims (&rx_2022_n claims);
      %put   output.master_outpat_2022 - Outpatient claims (&op_2022_n claims);
      %put   output.master_inpat_2022 - Inpatient claims (&ip_2022_n claims);
      %put;
      %put 2023 CLAIMS FILES (ALL CLAIMS FOR COHORT PATIENTS):;
      %put   output.master_rx_2023 - Pharmacy claims (&rx_2023_n claims);
      %put   output.master_outpat_2023 - Outpatient claims (&op_2023_n claims);
      %put   output.master_inpat_2023 - Inpatient claims (&ip_2023_n claims);
      %put;
      %put SUMMARY FILES:;
      %put   output.claims_extraction_summary - Claims counts by type/year;
      %put;
      %put =============================================================================;
      %put HOW TO USE THESE FILES:;
      %put =============================================================================;
      %put;
      %put 1. Merge claims with master_patient_list to add cohort flags:;
      %put;
      %put    proc sql;;
      %put        create table work.rx_2022_analysis as;
      %put        select c.*, p.in_cohort_2022, p.treatment_2022, p.age_2022;
      %put        from output.master_rx_2022 as c;
      %put        inner join output.master_patient_list as p;
      %put            on c.pat_id = p.pat_id;
      %put    quit;;
      %put;
      %put 2. Filter claims by cohort and date range:;
      %put;
      %put    /* For 2023 cross-sectional (Jan-Oct only) */;
      %put    data work.rx_2023_jan_oct;
      %put        set work.rx_2023_analysis;
      %put        where in_cohort_2023=1 and month(from_dt) <= 10;
      %put    run;;
      %put;
      %put 3. Calculate insulin out-of-pocket costs:;
      %put;
      %put    data work.insulin_oop;
      %put        set work.rx_analysis;
      %put        where insulin_claim=1;
      %put        oop_cost = sum(copay, coinsamt, deductible);
      %put    run;;
      %put;
      %put =============================================================================;
      %put;
      
      /* Clean up work library */
        proc datasets library=work nolist kill;
      quit;
