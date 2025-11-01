/*=============================================================================*/
/* IRA INSULIN COPAY CAP ANALYSIS - PANEL 22 MONTHS COHORT                   */
/* Only showing CHANGES from cross-sectional version                          */
/*=============================================================================*/
    
  libname inpat "R:/IQVIA PharMetrics Plus (2024)/Claims - Medical - Inpatient";
  libname outpat "R:/IQVIA PharMetrics Plus (2024)/Claims - Medical - Outpatient";
  libname pharm "R:/IQVIA PharMetrics Plus (2024)/Claims - Retail Pharmacy";
  libname enroll "R:/IQVIA PharMetrics Plus (2024)/Enrollment";
  libname output "R:/IQVIA PharMetrics Plus (2024) Members/WatsonWilliam/New_copay_project";
  
  options compress=yes sortsize=8G sumsize=8G memsize=32G threads cpucount=6 fullstimer;
  
  /*=============================================================================*/
    /* STEP 1: DEMOGRAPHICS - CHANGE: Only need 2023 age for classification      */
    /*=============================================================================*/
      data work.demographics;
    set enroll.enroll_synth;
    
    /* Only need 2023 age for panel classification */
      if der_yob = 0 then age_2023 = 86;
      else if der_yob > 0 then age_2023 = 2023 - der_yob;
      else delete;
      
      if age_2023 >= 18 and age_2023 <= 110;
      
      if pat_state not in ('CA','CO','CT','DE','IL','ME','NH','NM',
                           'NY','UT','VA','WA','WV');
      
      /* Age group for 2023 only */
        length age_group_2023 $10;
      if age_2023 >= 18 and age_2023 <= 61 then age_group_2023 = '18-61';
      else if age_2023 >= 62 and age_2023 <= 64 then age_group_2023 = '62-64';
      else if age_2023 >= 65 then age_group_2023 = '65+';
      
      keep pat_id der_sex pat_state der_yob age_2023 age_group_2023;
      run;
      
      /*=============================================================================*/
        /* STEP 2: ENROLLMENT - CHANGE: Extract 22-month period                      */
        /*=============================================================================*/
          proc sql;
        create table work.enrollment_raw as
        select e.pat_id,
        max(case when e.string_type = 'pay_type' then e.string_value else '' end) as pay_string,
        max(case when e.string_type = 'ben_rx' then e.string_value else '' end) as rx_string
        from enroll.enroll2 as e
        inner join work.demographics as d
        on e.pat_id = d.pat_id
        where e.string_type in ('pay_type', 'ben_rx')
        group by e.pat_id;
        quit;
        
        /* CHANGE: Extract 22 months (positions 253-274) */
          data work.enrollment;
        set work.enrollment_raw;
        
        /* 22 months: Jan 2022 (pos 253) through Oct 2023 (pos 274) */
          pay_22mo = substr(pay_string, 253, 22);
        rx_22mo = substr(rx_string, 253, 22);
        
        keep pat_id pay_22mo rx_22mo;
        run;
        
        /*=============================================================================*/
          /* STEP 3: ENROLLMENT ELIGIBILITY - CHANGE: Check 22 consecutive months      */
          /*=============================================================================*/
          data work.enrollment_eligible;
        set work.enrollment;
        
        length eligible_panel 8 payer_panel $15;
        eligible_panel = 0;
        payer_panel = '';
        
        /* Parse 22 months */
          array pay{22} $1 _temporary_;
        array rx{22} $1 _temporary_;
        
        do i = 1 to 22;
        pay{i} = substr(pay_22mo, i, 1);
        rx{i} = substr(rx_22mo, i, 1);
        end;
        
        /* Determine primary payer across 22 months */
          n_c = 0; n_ma = 0;
        do i = 1 to 22;
        if pay{i} = 'C' then n_c + 1;
        else if pay{i} in ('A','R') then n_ma + 1;
        end;
        
        if n_c > n_ma then payer_panel = 'C';
        else if n_ma > n_c then payer_panel = 'MA';
        
        /* Check for 22 consecutive months with payer AND Rx */
          if payer_panel ne '' then do;
        consec = 0;
        max_consec = 0;
        
        do i = 1 to 22;
        if payer_panel = 'C' then do;
        if pay{i} = 'C' and rx{i} = 'Y' then do;
        consec + 1;
        max_consec = max(max_consec, consec);
        end;
        else consec = 0;
        end;
        else if payer_panel = 'MA' then do;
        if pay{i} in ('A','R') and rx{i} = 'Y' then do;
        consec + 1;
        max_consec = max(max_consec, consec);
        end;
        else consec = 0;
        end;
        end;
        
        if max_consec >= 22 then eligible_panel = 1;
        end;
        
        keep pat_id eligible_panel payer_panel;
        drop i n_c n_ma consec max_consec;
        run;
        
        /*=============================================================================*/
          /* STEP 4: APPLY AGE FILTERS - CHANGE: Use 2023 age, single cohort flag      */
          /*=============================================================================*/
          proc sql;
        create table work.cohort_enrollment as
        select d.*, e.*
          from work.demographics as d
        inner join work.enrollment_eligible as e
        on d.pat_id = e.pat_id;
        quit;
        
        data work.cohort_eligible;
        set work.cohort_enrollment;
        
        /* CHANGE: Single panel cohort flag */
          length in_cohort_panel 8 treatment_panel 8;
        in_cohort_panel = 0;
        treatment_panel = .;
        
        if eligible_panel = 1 then do;
        /* Commercial 18-64 (Control) - use 2023 age */
          if payer_panel = 'C' and age_2023 >= 18 and age_2023 <= 64 then do;
        in_cohort_panel = 1;
        treatment_panel = 0;
        end;
        /* Medicare Advantage 65+ (Treatment) - use 2023 age */
          else if payer_panel = 'MA' and age_2023 >= 65 then do;
        in_cohort_panel = 1;
        treatment_panel = 1;
        end;
        end;
        
        if in_cohort_panel = 1;
        run;
        
        /*=============================================================================*/
          /* STEP 5: T2D IDENTIFICATION - CHANGE: Only need 2022 baseline              */
          /*=============================================================================*/
          proc sql;
        create table work.pts_panel as
        select distinct pat_id from work.cohort_eligible where in_cohort_panel = 1;
        quit;
        
        /* CHANGE: Only check 2022 for baseline T2D */
          data work.op_2022;
        if _n_ = 1 then do;
        declare hash h(dataset: "work.pts_panel");
        h.definekey('pat_id');
        h.definedone();
        end;
        
        set outpat.clm_outpat_22;
        if h.check() = 0;
        
        array diags{12} $ diag1-diag12;
        has_t2d = 0;
        do i = 1 to 12;
        if substr(diags{i}, 1, 3) = 'E11' then do;
        has_t2d = 1;
        leave;
        end;
        end;
        
        if has_t2d = 1;
        keep pat_id;
        run;
        
        data work.ip_2022;
        if _n_ = 1 then do;
        declare hash h(dataset: "work.pts_panel");
        h.definekey('pat_id');
        h.definedone();
        end;
        
        set inpat.clm_inpat_22;
        if h.check() = 0;
        
        if substr(diag_admit, 1, 3) = 'E11' then has_t2d = 1;
        else do;
        array diags{12} $ diag1-diag12;
        has_t2d = 0;
        do i = 1 to 12;
        if substr(diags{i}, 1, 3) = 'E11' then do;
        has_t2d = 1;
        leave;
        end;
        end;
        end;
        
        if has_t2d = 1;
        keep pat_id;
        run;
        
        data work.t2d_2022;
        set work.op_2022 work.ip_2022;
        run;
        
        proc sort data=work.t2d_2022 nodupkey;
        by pat_id;
        run;
        
        proc datasets library=work nolist;
        delete op_2022 ip_2022;
        quit;
        
        /*=============================================================================*/
          /* STEP 6: INSULIN IDENTIFICATION - CHANGE: Only need 2022 baseline          */
          /*=============================================================================*/
          data work.insulin_2022;
        if _n_ = 1 then do;
        declare hash h(dataset: "work.t2d_2022");
        h.definekey('pat_id');
        h.definedone();
        end;
        
        set pharm.clm_rx_rtl_22;
        if h.check() = 0;
        
        if substr(gpi14, 1, 4) = '2710';
        keep pat_id;
        run;
        
        proc sort data=work.insulin_2022 nodupkey;
        by pat_id;
        run;
        
        /*=============================================================================*/
          /* STEP 7: FINAL COHORT - CHANGE: Single set of flags                        */
          /*=============================================================================*/
          proc sql;
        create table work.cohort_with_clinical as
        select c.*,
        case when t.pat_id is not null then 1 else 0 end as has_t2d_2022,
        case when i.pat_id is not null then 1 else 0 end as has_insulin_2022
        from work.cohort_eligible as c
        left join work.t2d_2022 as t
        on c.pat_id = t.pat_id
        left join work.insulin_2022 as i
        on c.pat_id = i.pat_id;
        quit;
        
        /* CHANGE: Only check 2022 baseline for T2D and insulin */
          data work.cohort_final;
        set work.cohort_with_clinical;
        
        if in_cohort_panel = 1 and (has_t2d_2022 = 0 or has_insulin_2022 = 0) then do;
        in_cohort_panel = 0;
        treatment_panel = .;
        end;
        
        if in_cohort_panel = 1;
        run;
        
        data output.cohort_panel_22;
        set work.cohort_final;
        run;
        
        /*=============================================================================*/
          /* STEP 8: SUMMARY STATISTICS - CHANGE: Single cohort summary                */
          /*=============================================================================*/
          title 'PANEL 22 MONTHS COHORT (Jan 2022 - Oct 2023)';
        proc freq data=output.cohort_panel_22;
        tables treatment_panel * payer_panel / nocol norow nopercent missing;
        run;
        
        proc freq data=output.cohort_panel_22;
        tables treatment_panel * age_group_2023 / nocol norow nopercent;
        run;
        
        proc means data=output.cohort_panel_22 n mean std min p50 max;
        class treatment_panel;
        var age_2023;
        run;
        
        proc sql;
        create table output.panel_summary as
        select 
        'Panel 22 Months' as cohort length=30,
        count(*) as n_total format=comma12.,
        sum(case when treatment_panel=0 then 1 else 0 end) as n_control format=comma12.,
        sum(case when treatment_panel=1 then 1 else 0 end) as n_treatment format=comma12.,
        sum(case when treatment_panel=0 and age_group_2023='62-64' then 1 else 0 end) 
        as n_control_62_64 format=comma12. label='Control 62-64'
        from output.cohort_panel_22;
        quit;
        
        proc print data=output.panel_summary noobs;
        run;
        
        proc export data=output.panel_summary
        outfile="R:/IQVIA PharMetrics Plus (2024) Members/WatsonWilliam/New_copay_project/panel_summary.csv"
        dbms=csv replace;
        run;
        
        proc export data=output.cohort_panel_22
        outfile="R:/IQVIA PharMetrics Plus (2024) Members/WatsonWilliam/New_copay_project/cohort_panel_22.csv"
        dbms=csv replace;
        run;
        
        proc datasets library=work nolist kill;
        quit;
