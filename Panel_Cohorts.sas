/*=============================================================================*/
/* IRA INSULIN COPAY CAP ANALYSIS - PANEL 22 MONTHS COHORT                   */
/* WITH DETAILED INCLUSION/EXCLUSION TRACKING                                 */
/*=============================================================================*/

libname inpat "R:/IQVIA PharMetrics Plus (2024)/Claims - Medical - Inpatient";
libname outpat "R:/IQVIA PharMetrics Plus (2024)/Claims - Medical - Outpatient";
libname pharm "R:/IQVIA PharMetrics Plus (2024)/Claims - Retail Pharmacy";
libname enroll "R:/IQVIA PharMetrics Plus (2024)/Enrollment";
libname output "R:/IQVIA PharMetrics Plus (2024) Members/WatsonWilliam/New_copay_project";

options compress=yes sortsize=8G sumsize=8G memsize=32G threads cpucount=6 fullstimer;

%let start_time = %sysfunc(datetime());

/*=============================================================================*/
/* CREATE TRACKING DATASET FOR CONSORT FLOWCHART                             */
/*=============================================================================*/
data work.consort_tracking_panel;
    length step 8 description $200 n_patients 8 n_excluded 8;
    stop;
run;

%macro track_step(step=, desc=, dataset=);
    proc sql noprint;
        select count(distinct pat_id) into :n_current trimmed
        from &dataset;
    quit;
    
    data work.consort_tracking_panel;
        set work.consort_tracking_panel end=eof;
        output;
        if eof then do;
            step = &step;
            description = "&desc";
            n_patients = &n_current;
            if step > 1 then n_excluded = lag(n_patients) - n_patients;
            else n_excluded = .;
            output;
        end;
    run;
    
    %put NOTE: Step &step: &desc - N=&n_current;
%mend;

%put =============================================================================;
%put PANEL 22 MONTHS COHORT WITH COMPREHENSIVE TRACKING;
%put Analysis Period: January 2022 - October 2023 (22 months);
%put =============================================================================;

/*=============================================================================*/
/* STEP 0: COUNT TOTAL UNIVERSE                                              */
/*=============================================================================*/
%put =============================================================================;
%put STEP 0: TOTAL DATABASE UNIVERSE;
%put =============================================================================;

proc sql noprint;
    select count(distinct pat_id) into :n_universe trimmed
    from enroll.enroll_synth;
    
    select count(*) into :n_total_rx_22 trimmed
    from pharm.clm_rx_rtl_22;
    
    select count(*) into :n_total_rx_23 trimmed
    from pharm.clm_rx_rtl_23;
    
    select count(*) into :n_total_op_22 trimmed
    from outpat.clm_outpat_22;
    
    select count(*) into :n_total_op_23 trimmed
    from outpat.clm_outpat_23;
    
    select count(*) into :n_total_ip_22 trimmed
    from inpat.clm_inpat_22;
    
    select count(*) into :n_total_ip_23 trimmed
    from inpat.clm_inpat_23;
quit;

%put;
%put DATABASE UNIVERSE:;
%put   Total unique patients: &n_universe;
%put   Total pharmacy claims 2022: &n_total_rx_22;
%put   Total pharmacy claims 2023: &n_total_rx_23;
%put   Total outpatient claims 2022: &n_total_op_22;
%put   Total outpatient claims 2023: &n_total_op_23;
%put   Total inpatient claims 2022: &n_total_ip_22;
%put   Total inpatient claims 2023: &n_total_ip_23;
%put;

/*=============================================================================*/
/* STEP 1: DEMOGRAPHICS WITH TRACKING                                         */
/*=============================================================================*/
%put =============================================================================;
%put STEP 1: APPLYING BASE DEMOGRAPHIC FILTERS;
%put =============================================================================;

proc sql noprint;
    select count(distinct pat_id) into :n_before_demo trimmed
    from enroll.enroll_synth;
quit;
%put NOTE: Patients before demographic filters: &n_before_demo;

data work.demographics;
    set enroll.enroll_synth;
    
    /* Calculate 2023 age for panel classification */
    if der_yob = 0 then age_2023 = 86;
    else if der_yob > 0 then age_2023 = 2023 - der_yob;
    else delete;
    
    /* Flag exclusions */
    length exclude_age 8 exclude_state 8;
    exclude_age = 0;
    exclude_state = 0;
    
    /* Age filter */
    if age_2023 < 18 or age_2023 > 110 then exclude_age = 1;
    
    /* State filter */
    if pat_state in ('CA','CO','CT','DE','IL','ME','NH','NM',
                     'NY','UT','VA','WA','WV') then exclude_state = 1;
    
    /* Apply filters */
    if exclude_age = 0 and exclude_state = 0;
    
    /* Age group for 2023 */
    length age_group_2023 $10;
    if age_2023 >= 18 and age_2023 <= 61 then age_group_2023 = '18-61';
    else if age_2023 >= 62 and age_2023 <= 64 then age_group_2023 = '62-64';
    else if age_2023 >= 65 then age_group_2023 = '65+';
    
    keep pat_id der_sex pat_state der_yob age_2023 age_group_2023;
run;

proc sql noprint;
    select count(distinct pat_id) into :n_after_demo trimmed
    from work.demographics;
    
    %let n_excl_age_state = %sysevalf(&n_before_demo - &n_after_demo);
quit;

%put NOTE: After demographic filters: &n_after_demo;
%put NOTE: Excluded (age <18/state with copay cap): &n_excl_age_state;

%track_step(step=1, desc=After demographic filters (age 18+ in 2023 and no copay cap state), 
            dataset=work.demographics);

/*=============================================================================*/
/* STEP 2: EXTRACT ENROLLMENT STRINGS FOR 22 MONTHS                          */
/*=============================================================================*/
%put =============================================================================;
%put STEP 2: EXTRACTING 22-MONTH ENROLLMENT STRINGS;
%put =============================================================================;

proc sql;
    create table work.enrollment_raw as
    select e.pat_id,
           max(case when e.string_type = 'pay_type' then e.string_value else '' end) as pay_string,
           max(case when e.string_type = 'ben_rx' then e.string_value else '' end) as rx_string,
           max(case when e.string_type = 'mcob_type' then e.string_value else '' end) as mcob_string,
           max(case when e.string_type = 'pcob_type' then e.string_value else '' end) as pcob_string
    from enroll.enroll2 as e
    inner join work.demographics as d
        on e.pat_id = d.pat_id
    where e.string_type in ('pay_type', 'ben_rx', 'mcob_type', 'pcob_type')
    group by e.pat_id;
quit;

%put NOTE: Enrollment strings merged - N=%trim(%left(&sqlobs));

/* Extract 22 months: Jan 2022 (pos 253) through Oct 2023 (pos 274) */
data work.enrollment;
    set work.enrollment_raw;
    
    /* 22 months spanning 2022-2023 */
    pay_22mo = substr(pay_string, 253, 22);
    rx_22mo = substr(rx_string, 253, 22);
    mcob_22mo = substr(mcob_string, 253, 22);
    pcob_22mo = substr(pcob_string, 253, 22);
    
    keep pat_id pay_22mo rx_22mo mcob_22mo pcob_22mo;
run;

/*=============================================================================*/
/* STEP 3: CALCULATE ENROLLMENT ELIGIBILITY WITH TRACKING                    */
/*=============================================================================*/
%put =============================================================================;
%put STEP 3: CALCULATING 22-MONTH ENROLLMENT ELIGIBILITY;
%put =============================================================================;

data work.enrollment_eligible;
    set work.enrollment;
    
    length eligible_panel 8 payer_panel $15 has_medicare_cob 8;
    length reason_inelig $100;
    eligible_panel = 0;
    payer_panel = '';
    has_medicare_cob = 0;
    reason_inelig = '';
    
    /* Parse 22 months */
    array pay{22} $1 _temporary_;
    array rx{22} $1 _temporary_;
    array mcob{22} $1 _temporary_;
    array pcob{22} $1 _temporary_;
    
    do i = 1 to 22;
        pay{i} = substr(pay_22mo, i, 1);
        rx{i} = substr(rx_22mo, i, 1);
        mcob{i} = substr(mcob_22mo, i, 1);
        pcob{i} = substr(pcob_22mo, i, 1);
        
        /* Check for Medicare COB at any point */
        if mcob{i} = 'M' or pcob{i} = 'M' then has_medicare_cob = 1;
    end;
    
    /* Determine primary payer across 22 months */
    n_c = 0; n_ma = 0;
    do i = 1 to 22;
        if pay{i} = 'C' then n_c + 1;
        else if pay{i} in ('A','R') then n_ma + 1;
    end;
    
    if n_c > n_ma then payer_panel = 'C';
    else if n_ma > n_c then payer_panel = 'MA';
    else if n_c = 0 and n_ma = 0 then reason_inelig = 'No payer coverage';
    else reason_inelig = 'Mixed payer (no majority)';
    
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
        else reason_inelig = catx(': ', reason_inelig, 
                                   cats('Max consecutive=', max_consec, ' months'));
    end;
    
    keep pat_id eligible_panel payer_panel has_medicare_cob reason_inelig;
    drop i n_c n_ma consec max_consec;
run;

/* Count enrollment eligibility */
proc sql noprint;
    select count(distinct pat_id) into :n_elig_panel trimmed
    from work.enrollment_eligible where eligible_panel=1;
    
    select count(distinct pat_id) into :n_not_elig_panel trimmed
    from work.enrollment_eligible where eligible_panel=0;
quit;

%put;
%put PANEL ENROLLMENT ELIGIBILITY:;
%put   Eligible (22+ consecutive months): &n_elig_panel;
%put   Not eligible: &n_not_elig_panel;
%put;

/* Reasons for ineligibility */
title 'Panel 22 Months - Enrollment Ineligibility Reasons';
proc freq data=work.enrollment_eligible;
    where eligible_panel = 0;
    tables reason_inelig / nocum;
run;
title;

/*=============================================================================*/
/* STEP 4: APPLY AGE AND COB FILTERS                                         */
/*=============================================================================*/
%put =============================================================================;
%put STEP 4: APPLYING AGE AND COORDINATION OF BENEFITS FILTERS;
%put =============================================================================;

proc sql;
    create table work.cohort_enrollment as
    select d.*, e.*
    from work.demographics as d
    inner join work.enrollment_eligible as e
        on d.pat_id = e.pat_id;
quit;

data work.cohort_eligible;
    set work.cohort_enrollment;
    
    /* Track reasons for exclusion */
    length excl_reason $100;
    length in_cohort_panel 8 treatment_panel 8;
    in_cohort_panel = 0;
    treatment_panel = .;
    excl_reason = '';
    
    if eligible_panel = 1 then do;
        /* Commercial 18-64 (Control) - use 2023 age */
        if payer_panel = 'C' then do;
            if age_2023 < 18 or age_2023 > 64 then 
                excl_reason = 'Commercial: Age not 18-64';
            else if has_medicare_cob = 1 then 
                excl_reason = 'Commercial: Has Medicare COB';
            else do;
                in_cohort_panel = 1;
                treatment_panel = 0;
            end;
        end;
        /* Medicare Advantage 65+ (Treatment) - use 2023 age */
        else if payer_panel = 'MA' then do;
            if age_2023 < 65 then 
                excl_reason = 'MA: Age <65';
            else do;
                in_cohort_panel = 1;
                treatment_panel = 1;
            end;
        end;
        else excl_reason = 'No valid payer type';
    end;
    else excl_reason = reason_inelig;
    
    /* Keep only if eligible for cohort */
    if in_cohort_panel = 1;
run;

/* Track counts */
proc sql noprint;
    select count(distinct pat_id) into :n_after_age_cob trimmed
    from work.cohort_eligible;
    
    select count(distinct pat_id) into :n_control trimmed
    from work.cohort_eligible where treatment_panel=0;
    
    select count(distinct pat_id) into :n_treatment trimmed
    from work.cohort_eligible where treatment_panel=1;
quit;

%put;
%put AFTER AGE AND COB FILTERS:;
%put   Patients in panel cohort: &n_after_age_cob;
%put   Control (Commercial 18-64): &n_control;
%put   Treatment (MA 65+): &n_treatment;
%put;

%track_step(step=2, desc=After enrollment (22 consecutive months with RX) and age/COB filters, 
            dataset=work.cohort_eligible);

/* Exclusion reasons */
title 'Panel 22 Months - Exclusion Reasons (Among Enrolled)';
proc freq data=work.cohort_enrollment;
    where eligible_panel=1 and in_cohort_panel ne 1;
    tables excl_reason / nocum;
run;
title;

/*=============================================================================*/
/* STEP 5: IDENTIFY T2D PATIENTS (2022 BASELINE) WITH TRACKING               */
/*=============================================================================*/
%put =============================================================================;
%put STEP 5: IDENTIFYING TYPE 2 DIABETES PATIENTS (2022 BASELINE);
%put =============================================================================;

/* Create patient list */
proc sql;
    create table work.pts_panel as
    select distinct pat_id from work.cohort_eligible where in_cohort_panel = 1;
quit;

%put NOTE: Panel patient list created - N=%trim(%left(&sqlobs));

/* Outpatient 2022 - count claims */
data work.op_2022_claims;
    if _n_ = 1 then do;
        declare hash h(dataset: "work.pts_panel");
        h.definekey('pat_id');
        h.definedone();
    end;
    
    set outpat.clm_outpat_22;
    
    if h.check() = 0;
    
    /* Check for E11* */
    array diags{12} $ diag1-diag12;
    has_t2d = 0;
    do i = 1 to 12;
        if substr(diags{i}, 1, 3) = 'E11' then do;
            has_t2d = 1;
            leave;
        end;
    end;
    
    keep pat_id claimno linenum has_t2d;
run;

proc sql noprint;
    select count(*) into :n_op_claims_panel trimmed
    from work.op_2022_claims;
    
    select count(*) into :n_op_t2d_claims_panel trimmed
    from work.op_2022_claims where has_t2d=1;
quit;

/* Get unique patients */
proc sql;
    create table work.op_2022 as
    select distinct pat_id
    from work.op_2022_claims
    where has_t2d = 1;
quit;

/* Inpatient 2022 - count claims */
data work.ip_2022_claims;
    if _n_ = 1 then do;
        declare hash h(dataset: "work.pts_panel");
        h.definekey('pat_id');
        h.definedone();
    end;
    
    set inpat.clm_inpat_22;
    
    if h.check() = 0;
    
    /* Check diag_admit and diag1-diag12 */
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
    
    keep pat_id claimno linenum has_t2d;
run;

proc sql noprint;
    select count(*) into :n_ip_claims_panel trimmed
    from work.ip_2022_claims;
    
    select count(*) into :n_ip_t2d_claims_panel trimmed
    from work.ip_2022_claims where has_t2d=1;
quit;

/* Get unique patients */
proc sql;
    create table work.ip_2022 as
    select distinct pat_id
    from work.ip_2022_claims
    where has_t2d = 1;
quit;

/* Combine */
data work.t2d_2022;
    set work.op_2022 work.ip_2022;
run;

proc sort data=work.t2d_2022 nodupkey;
    by pat_id;
run;

%put;
%put T2D IDENTIFICATION - 2022 BASELINE:;
%put   Outpatient claims searched: &n_op_claims_panel;
%put   Outpatient claims with E11*: &n_op_t2d_claims_panel;
%put   Inpatient claims searched: &n_ip_claims_panel;
%put   Inpatient claims with E11*: &n_ip_t2d_claims_panel;
%put   Unique patients with T2D: %trim(%left(&sysnobs));
%put;

%track_step(step=3, desc=After T2D diagnosis requirement in 2022 (E11*), 
            dataset=work.t2d_2022);

/* Clean up claim datasets */
proc datasets library=work nolist;
    delete op_2022_claims ip_2022_claims op_2022 ip_2022;
quit;

/*=============================================================================*/
/* STEP 6: IDENTIFY INSULIN USERS (2022 BASELINE) WITH TRACKING              */
/*=============================================================================*/
%put =============================================================================;
%put STEP 6: IDENTIFYING INSULIN USERS (2022 BASELINE);
%put =============================================================================;

data work.insulin_2022_claims;
    if _n_ = 1 then do;
        declare hash h(dataset: "work.t2d_2022");
        h.definekey('pat_id');
        h.definedone();
    end;
    
    set pharm.clm_rx_rtl_22;
    
    if h.check() = 0;
    
    /* Check for insulin */
    is_insulin = (substr(gpi14, 1, 4) = '2710');
    
    keep pat_id claimno linenum is_insulin;
run;

proc sql noprint;
    select count(*) into :n_rx_claims_panel trimmed
    from work.insulin_2022_claims;
    
    select count(*) into :n_insulin_claims_panel trimmed
    from work.insulin_2022_claims where is_insulin=1;
quit;

/* Get unique patients */
proc sql;
    create table work.insulin_2022 as
    select distinct pat_id
    from work.insulin_2022_claims
    where is_insulin = 1;
quit;

%put;
%put INSULIN IDENTIFICATION - 2022 BASELINE:;
%put   Pharmacy claims searched: &n_rx_claims_panel;
%put   Insulin claims (GPI 2710*): &n_insulin_claims_panel;
%put   Unique patients with insulin: %trim(%left(&sysnobs));
%put;

%track_step(step=4, desc=Final panel cohort (after insulin requirement in 2022), 
            dataset=work.insulin_2022);

/* Clean up */
proc datasets library=work nolist;
    delete insulin_2022_claims;
quit;

/*=============================================================================*/
/* STEP 7: CREATE FINAL COHORT WITH TRACKING                                 */
/*=============================================================================*/
%put =============================================================================;
%put STEP 7: CREATING FINAL PANEL COHORT;
%put =============================================================================;

/* Merge T2D and insulin flags */
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

/* Apply T2D and insulin requirements with tracking */
data work.cohort_final;
    set work.cohort_with_clinical;
    
    /* Update exclusion reasons */
    if in_cohort_panel = 1 then do;
        if has_t2d_2022 = 0 then do;
            in_cohort_panel = 0;
            treatment_panel = .;
            excl_reason = 'No T2D diagnosis in 2022';
        end;
        else if has_insulin_2022 = 0 then do;
            in_cohort_panel = 0;
            treatment_panel = .;
            excl_reason = 'No insulin use in 2022';
        end;
    end;
    
    /* Keep only if in cohort */
    if in_cohort_panel = 1;
run;

/* Final counts */
proc sql noprint;
    select count(distinct pat_id) into :n_final_panel trimmed
    from work.cohort_final;
    
    select count(distinct pat_id) into :n_final_control trimmed
    from work.cohort_final where treatment_panel=0;
    
    select count(distinct pat_id) into :n_final_treatment trimmed
    from work.cohort_final where treatment_panel=1;
quit;

%put;
%put FINAL PANEL COHORT COUNTS:;
%put   Total: &n_final_panel;
%put   Control (Commercial 18-64): &n_final_control;
%put   Treatment (MA 65+): &n_final_treatment;
%put;

/* Reasons for exclusion at T2D/insulin stage */
title 'Panel 22 Months - Exclusions at T2D/Insulin Stage';
proc freq data=work.cohort_with_clinical;
    where eligible_panel=1 and payer_panel ne '' and 
          ((payer_panel='C' and age_2023 between 18 and 64 and has_medicare_cob=0) or
           (payer_panel='MA' and age_2023>=65));
    tables has_t2d_2022*has_insulin_2022 / nocol norow nopercent missing;
run;
title;

/* Save final cohort */
data output.cohort_panel_22;
    set work.cohort_final;
run;

/*=============================================================================*/
/* STEP 8: COMPREHENSIVE SUMMARY STATISTICS                                  */
/*=============================================================================*/
%put =============================================================================;
%put STEP 8: COMPREHENSIVE SUMMARY STATISTICS;
%put =============================================================================;

/* Print consort tracking */
title 'CONSORT FLOW - PANEL 22 MONTHS COHORT CONSTRUCTION';
proc print data=work.consort_tracking_panel noobs label;
    var step description n_patients n_excluded;
    format n_patients n_excluded comma12.;
    label step = 'Step'
          description = 'Description'
          n_patients = 'N Patients'
          n_excluded = 'N Excluded';
run;

/* Export consort tracking */
proc export data=work.consort_tracking_panel
    outfile="R:/IQVIA PharMetrics Plus (2024) Members/WatsonWilliam/New_copay_project/consort_flowchart_panel_22.csv"
    dbms=csv replace;
run;

/* Detailed breakdowns */
title 'PANEL 22 MONTHS COHORT (Jan 2022 - Oct 2023) - FINAL';
proc freq data=output.cohort_panel_22;
    tables treatment_panel * payer_panel / nocol norow nopercent missing;
    tables treatment_panel * age_group_2023 / nocol norow nopercent;
    tables treatment_panel * der_sex / nocol norow nopercent;
run;

proc means data=output.cohort_panel_22 n mean std min p25 p50 p75 max;
    class treatment_panel;
    var age_2023;
    title2 'Age Distribution (Age in 2023)';
run;

/* Comprehensive summary table */
proc sql;
    create table output.panel_summary_detailed as
    
    /* Overall */
    select 
        'Panel 22 Months' as cohort length=30,
        'Overall' as group length=20,
        count(*) as n_total format=comma12.,
        sum(case when treatment_panel=0 then 1 else 0 end) as n_control format=comma12.,
        sum(case when treatment_panel=1 then 1 else 0 end) as n_treatment format=comma12.,
        round(mean(case when treatment_panel=0 then age_2023 else . end), 0.1) as age_control_mean format=5.1,
        round(mean(case when treatment_panel=1 then age_2023 else . end), 0.1) as age_treatment_mean format=5.1,
        sum(case when treatment_panel=0 and der_sex='M' then 1 else 0 end) as n_male_control format=comma12.,
        sum(case when treatment_panel=1 and der_sex='M' then 1 else 0 end) as n_male_treatment format=comma12.
    from output.cohort_panel_22
    
    union all
    
    /* By Age Group */
    select 
        'Panel 22 Months' as cohort,
        age_group_2023 as group,
        count(*) as n_total format=comma12.,
        sum(case when treatment_panel=0 then 1 else 0 end) as n_control format=comma12.,
        sum(case when treatment_panel=1 then 1 else 0 end) as n_treatment format=comma12.,
        round(mean(case when treatment_panel=0 then age_2023 else . end), 0.1) as age_control_mean format=5.1,
        round(mean(case when treatment_panel=1 then age_2023 else . end), 0.1) as age_treatment_mean format=5.1,
        sum(case when treatment_panel=0 and der_sex='M' then 1 else 0 end) as n_male_control format=comma12.,
        sum(case when treatment_panel=1 and der_sex='M' then 1 else 0 end) as n_male_treatment format=comma12.
    from output.cohort_panel_22
    group by age_group_2023
    
    order by cohort, group;
quit;

title 'DETAILED PANEL COHORT SUMMARY';
proc print data=output.panel_summary_detailed noobs;
run;
title;

/* Export summary files */
proc export data=output.panel_summary_detailed
    outfile="R:/IQVIA PharMetrics Plus (2024) Members/WatsonWilliam/New_copay_project/panel_summary_detailed.csv"
    dbms=csv replace;
run;

proc export data=output.cohort_panel_22
    outfile="R:/IQVIA PharMetrics Plus (2024) Members/WatsonWilliam/New_copay_project/cohort_panel_22.csv"
    dbms=csv replace;
run;

/*=============================================================================*/
/* COMPLETION WITH COMPREHENSIVE LOG                                         */
/*=============================================================================*/

%let end_time = %sysfunc(datetime());
%let elapsed = %sysevalf(&end_time - &start_time);

%put;
%put =============================================================================;
%put PANEL 22 MONTHS COHORT COMPLETE - COMPREHENSIVE TRACKING;
%put =============================================================================;
%put;
%put ELAPSED TIME: %sysfunc(putn(&elapsed, time12.2));
%put;
%put DATABASE UNIVERSE:;
%put   Total patients in database: &n_universe;
%put   Total RX claims 2022: &n_total_rx_22;
%put   Total RX claims 2023: &n_total_rx_23;
%put   Total Outpatient claims 2022: &n_total_op_22;
%put   Total Outpatient claims 2023: &n_total_op_23;
%put   Total Inpatient claims 2022: &n_total_ip_22;
%put   Total Inpatient claims 2023: &n_total_ip_23;
%put;
%put FINAL PANEL COHORT SUMMARY:;
%put   Panel 22 Months Total: &n_final_panel;
%put     - Control (Commercial 18-64 in 2023): &n_final_control;
%put     - Treatment (MA 65+ in 2023): &n_final_treatment;
%put;
%put CLINICAL CRITERIA (2022 BASELINE):;
%put   T2D diagnosis requirement: E11* in medical claims;
%put   Insulin use requirement: GPI 2710* in pharmacy claims;
%put;
%put CLAIMS EXAMINED:;
%put   Outpatient claims (T2D): &n_op_claims_panel (&n_op_t2d_claims_panel with E11*);
%put   Inpatient claims (T2D): &n_ip_claims_panel (&n_ip_t2d_claims_panel with E11*);
%put   Pharmacy claims (insulin): &n_rx_claims_panel (&n_insulin_claims_panel with GPI 2710*);
%put;
%put OUTPUT FILES CREATED:;
%put   1. output.cohort_panel_22 - Master panel cohort file;
%put   2. output.panel_summary_detailed - Detailed summary statistics;
%put   3. consort_flowchart_panel_22.csv - Step-by-step inclusion/exclusion;
%put   4. panel_summary_detailed.csv - Detailed summary export;
%put   5. cohort_panel_22.csv - Full cohort export;
%put;
%put =============================================================================;

/* Clean up work library */
proc datasets library=work nolist kill;
quit;
