/*=============================================================================*/
/* IRA INSULIN COPAY CAP ANALYSIS - CROSS-SECTIONAL COHORTS                  */
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
data work.consort_tracking;
    length step 8 year $4 description $200 n_patients 8 n_excluded 8;
    stop;
run;

%macro track_step(step=, year=, desc=, dataset=);
    proc sql noprint;
        select count(distinct pat_id) into :n_current trimmed
        from &dataset;
    quit;
    
    data work.consort_tracking;
        set work.consort_tracking end=eof;
        output;
        if eof then do;
            step = &step;
            year = "&year";
            description = "&desc";
            n_patients = &n_current;
            if step > 1 then n_excluded = lag(n_patients) - n_patients;
            else n_excluded = .;
            output;
        end;
    run;
    
    %put NOTE: Step &step (&year): &desc - N=&n_current;
%mend;

%put =============================================================================;
%put CROSS-SECTIONAL COHORTS WITH TRACKING - 2022 AND 2023;
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
/* STEP 1: LOAD BASE DEMOGRAPHICS WITH TRACKING                              */
/*=============================================================================*/
%put =============================================================================;
%put STEP 1: APPLYING BASE DEMOGRAPHIC FILTERS;
%put =============================================================================;

/* Count before any filters */
proc sql noprint;
    select count(distinct pat_id) into :n_before_demo trimmed
    from enroll.enroll_synth;
quit;
%put NOTE: Patients before demographic filters: &n_before_demo;

data work.demographics;
    set enroll.enroll_synth;
    
    /* Calculate ages */
    if der_yob = 0 then do;
        age_2022 = 86;
        age_2023 = 86;
    end;
    else if der_yob > 0 then do;
        age_2022 = 2022 - der_yob;
        age_2023 = 2023 - der_yob;
    end;
    else delete;
    
    /* Flag exclusions for tracking */
    length exclude_age 8 exclude_state 8;
    exclude_age = 0;
    exclude_state = 0;
    
    /* Age filter */
    if age_2022 < 18 and age_2023 < 18 then exclude_age = 1;
    if age_2022 > 110 or age_2023 > 110 then exclude_age = 1;
    
    /* State filter */
    if pat_state in ('CA','CO','CT','DE','IL','ME','NH','NM',
                     'NY','UT','VA','WA','WV') then exclude_state = 1;
    
    /* Apply filters */
    if exclude_age = 0 and exclude_state = 0;
    
    /* Create age groups */
    length age_group_2022 $10 age_group_2023 $10;
    if age_2022 >= 18 and age_2022 <= 61 then age_group_2022 = '18-61';
    else if age_2022 >= 62 and age_2022 <= 64 then age_group_2022 = '62-64';
    else if age_2022 >= 65 then age_group_2022 = '65+';
    
    if age_2023 >= 18 and age_2023 <= 61 then age_group_2023 = '18-61';
    else if age_2023 >= 62 and age_2023 <= 64 then age_group_2023 = '62-64';
    else if age_2023 >= 65 then age_group_2023 = '65+';
    
    keep pat_id der_sex pat_state der_yob age_2022 age_2023 
         age_group_2022 age_group_2023;
run;

/* Track exclusions */
proc sql noprint;
    select count(distinct pat_id) into :n_after_demo trimmed
    from work.demographics;
    
    %let n_excl_age = %sysevalf(&n_before_demo - &n_after_demo);
quit;

%put NOTE: After demographic filters: &n_after_demo;
%put NOTE: Excluded (age <18 both years or state with copay cap): &n_excl_age;

%track_step(step=1, year=Both, desc=After demographic filters (age 18+ in either year and no copay cap state), 
            dataset=work.demographics);

/*=============================================================================*/
/* STEP 2: EXTRACT ENROLLMENT STRINGS                                        */
/*=============================================================================*/
%put =============================================================================;
%put STEP 2: EXTRACTING ENROLLMENT STRINGS;
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

/* Extract 2022 and 2023 periods */
data work.enrollment;
    set work.enrollment_raw;
    
    /* 2022: positions 253-264 (12 months) */
    pay_2022 = substr(pay_string, 253, 12);
    rx_2022 = substr(rx_string, 253, 12);
    mcob_2022 = substr(mcob_string, 253, 12);
    pcob_2022 = substr(pcob_string, 253, 12);
    
    /* 2023: positions 265-276 */
    pay_2023 = substr(pay_string, 265, 10);
    rx_2023 = substr(rx_string, 265, 10);
    mcob_2023 = substr(mcob_string, 265, 10);
    pcob_2023 = substr(pcob_string, 265, 10);
    
    keep pat_id pay_2022 rx_2022 mcob_2022 pcob_2022
         pay_2023 rx_2023 mcob_2023 pcob_2023;
run;

/*=============================================================================*/
/* STEP 3: CALCULATE ENROLLMENT ELIGIBILITY WITH TRACKING                    */
/*=============================================================================*/
%put =============================================================================;
%put STEP 3: CALCULATING ENROLLMENT ELIGIBILITY;
%put =============================================================================;

data work.enrollment_eligible;
    set work.enrollment;
    
    /*=========================================================================*/
    /* 2022 ELIGIBILITY                                                       */
    /*=========================================================================*/
    length eligible_2022 8 payer_2022 $15 has_medicare_cob_2022 8;
    length reason_inelig_2022 $100;
    eligible_2022 = 0;
    payer_2022 = '';
    has_medicare_cob_2022 = 0;
    reason_inelig_2022 = '';
    
    /* Parse strings */
    array pay22{12} $1 _temporary_;
    array rx22{12} $1 _temporary_;
    array mcob22{12} $1 _temporary_;
    array pcob22{12} $1 _temporary_;
    
    do i = 1 to 12;
        pay22{i} = substr(pay_2022, i, 1);
        rx22{i} = substr(rx_2022, i, 1);
        mcob22{i} = substr(mcob_2022, i, 1);
        pcob22{i} = substr(pcob_2022, i, 1);
        
        if mcob22{i} = 'M' or pcob22{i} = 'M' then has_medicare_cob_2022 = 1;
    end;
    
    /* Determine primary payer */
    n_c = 0; n_ma = 0;
    do i = 1 to 12;
        if pay22{i} = 'C' then n_c + 1;
        else if pay22{i} in ('A','R') then n_ma + 1;
    end;
    
    if n_c > n_ma then payer_2022 = 'C';
    else if n_ma > n_c then payer_2022 = 'MA';
    else if n_c = 0 and n_ma = 0 then reason_inelig_2022 = 'No payer coverage';
    else reason_inelig_2022 = 'Mixed payer (no majority)';
    
    /* Check for 12 consecutive months */
    if payer_2022 ne '' then do;
        consec = 0;
        max_consec = 0;
        
        do i = 1 to 12;
            if payer_2022 = 'C' then do;
                if pay22{i} = 'C' and rx22{i} = 'Y' then do;
                    consec + 1;
                    max_consec = max(max_consec, consec);
                end;
                else consec = 0;
            end;
            else if payer_2022 = 'MA' then do;
                if pay22{i} in ('A','R') and rx22{i} = 'Y' then do;
                    consec + 1;
                    max_consec = max(max_consec, consec);
                end;
                else consec = 0;
            end;
        end;
        
        if max_consec >= 12 then eligible_2022 = 1;
        else reason_inelig_2022 = catx(': ', reason_inelig_2022, 
                                       cats('Max consecutive=', max_consec, ' months'));
    end;
    
    /*=========================================================================*/
    /* 2023 ELIGIBILITY                                                       */
    /*=========================================================================*/
    length eligible_2023 8 payer_2023 $15 has_medicare_cob_2023 8;
    length reason_inelig_2023 $100;
    eligible_2023 = 0;
    payer_2023 = '';
    has_medicare_cob_2023 = 0;
    reason_inelig_2023 = '';
    
    /* Parse strings */
    array pay23{10} $1 _temporary_;
    array rx23{10} $1 _temporary_;
    array mcob23{10} $1 _temporary_;
    array pcob23{10} $1 _temporary_;
    
    do i = 1 to 10;
        pay23{i} = substr(pay_2023, i, 1);
        rx23{i} = substr(rx_2023, i, 1);
        mcob23{i} = substr(mcob_2023, i, 1);
        pcob23{i} = substr(pcob_2023, i, 1);
        
        if mcob23{i} = 'M' or pcob23{i} = 'M' then has_medicare_cob_2023 = 1;
    end;
    
    /* Determine primary payer */
    n_c = 0; n_ma = 0;
    do i = 1 to 10;
        if pay23{i} = 'C' then n_c + 1;
        else if pay23{i} in ('A','R') then n_ma + 1;
    end;
    
    if n_c > n_ma then payer_2023 = 'C';
    else if n_ma > n_c then payer_2023 = 'MA';
    else if n_c = 0 and n_ma = 0 then reason_inelig_2023 = 'No payer coverage';
    else reason_inelig_2023 = 'Mixed payer (no majority)';
    
    /* Check for 10 consecutive months */
    if payer_2023 ne '' then do;
        consec = 0;
        max_consec = 0;
        
        do i = 1 to 10;
            if payer_2023 = 'C' then do;
                if pay23{i} = 'C' and rx23{i} = 'Y' then do;
                    consec + 1;
                    max_consec = max(max_consec, consec);
                end;
                else consec = 0;
            end;
            else if payer_2023 = 'MA' then do;
                if pay23{i} in ('A','R') and rx23{i} = 'Y' then do;
                    consec + 1;
                    max_consec = max(max_consec, consec);
                end;
                else consec = 0;
            end;
        end;
        
        if max_consec >= 10 then eligible_2023 = 1;
        else reason_inelig_2023 = catx(': ', reason_inelig_2023,
                                       cats('Max consecutive=', max_consec, ' months'));
    end;
    
    keep pat_id eligible_2022 payer_2022 has_medicare_cob_2022 reason_inelig_2022
         eligible_2023 payer_2023 has_medicare_cob_2023 reason_inelig_2023;
    drop i n_c n_ma consec max_consec;
run;

/* Count enrollment eligibility */
proc sql noprint;
    select count(distinct pat_id) into :n_elig_2022 trimmed
    from work.enrollment_eligible where eligible_2022=1;
    
    select count(distinct pat_id) into :n_elig_2023 trimmed
    from work.enrollment_eligible where eligible_2023=1;
    
    select count(distinct pat_id) into :n_not_elig_2022 trimmed
    from work.enrollment_eligible where eligible_2022=0;
    
    select count(distinct pat_id) into :n_not_elig_2023 trimmed
    from work.enrollment_eligible where eligible_2023=0;
quit;

%put;
%put ENROLLMENT ELIGIBILITY:;
%put   2022 - Eligible: &n_elig_2022;
%put   2022 - Not eligible: &n_not_elig_2022;
%put   2023 - Eligible: &n_elig_2023;
%put   2023 - Not eligible: &n_not_elig_2023;
%put;

/* Reasons for ineligibility */
title '2022 Enrollment Ineligibility Reasons';
proc freq data=work.enrollment_eligible;
    where eligible_2022 = 0;
    tables reason_inelig_2022 / nocum;
run;

title '2023 Enrollment Ineligibility Reasons';
proc freq data=work.enrollment_eligible;
    where eligible_2023 = 0;
    tables reason_inelig_2023 / nocum;
run;
title;

/*=============================================================================*/
/* STEP 4: MERGE WITH DEMOGRAPHICS AND APPLY AGE/COB FILTERS                 */
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
    length excl_reason_2022 $100 excl_reason_2023 $100;
    
    /*=========================================================================*/
    /* 2022 COHORT                                                            */
    /*=========================================================================*/
    length in_cohort_2022 8 treatment_2022 8;
    in_cohort_2022 = 0;
    treatment_2022 = .;
    excl_reason_2022 = '';
    
    if eligible_2022 = 1 then do;
        /* Commercial 18-64 (Control) */
        if payer_2022 = 'C' then do;
            if age_2022 < 18 or age_2022 > 64 then 
                excl_reason_2022 = 'Commercial: Age not 18-64';
            else if has_medicare_cob_2022 = 1 then 
                excl_reason_2022 = 'Commercial: Has Medicare COB';
            else do;
                in_cohort_2022 = 1;
                treatment_2022 = 0;
            end;
        end;
        /* Medicare Advantage 65+ (Treatment) */
        else if payer_2022 = 'MA' then do;
            if age_2022 < 65 then 
                excl_reason_2022 = 'MA: Age <65';
            else do;
                in_cohort_2022 = 1;
                treatment_2022 = 1;
            end;
        end;
        else excl_reason_2022 = 'No valid payer type';
    end;
    else excl_reason_2022 = reason_inelig_2022;
    
    /*=========================================================================*/
    /* 2023 COHORT                                                            */
    /*=========================================================================*/
    length in_cohort_2023 8 treatment_2023 8;
    in_cohort_2023 = 0;
    treatment_2023 = .;
    excl_reason_2023 = '';
    
    if eligible_2023 = 1 then do;
        /* Commercial 18-64 (Control) */
        if payer_2023 = 'C' then do;
            if age_2023 < 18 or age_2023 > 64 then 
                excl_reason_2023 = 'Commercial: Age not 18-64';
            else if has_medicare_cob_2023 = 1 then 
                excl_reason_2023 = 'Commercial: Has Medicare COB';
            else do;
                in_cohort_2023 = 1;
                treatment_2023 = 0;
            end;
        end;
        /* Medicare Advantage 65+ (Treatment) */
        else if payer_2023 = 'MA' then do;
            if age_2023 < 65 then 
                excl_reason_2023 = 'MA: Age <65';
            else do;
                in_cohort_2023 = 1;
                treatment_2023 = 1;
            end;
        end;
        else excl_reason_2023 = 'No valid payer type';
    end;
    else excl_reason_2023 = reason_inelig_2023;
    
    /* Keep only if eligible for at least one cohort */
    if in_cohort_2022 = 1 or in_cohort_2023 = 1;
run;

/* Track exclusions */
proc sql noprint;
    select count(distinct pat_id) into :n_after_age_cob trimmed
    from work.cohort_eligible;
    
    select count(distinct pat_id) into :n_in_2022 trimmed
    from work.cohort_eligible where in_cohort_2022=1;
    
    select count(distinct pat_id) into :n_in_2023 trimmed
    from work.cohort_eligible where in_cohort_2023=1;
quit;

%put;
%put AFTER AGE AND COB FILTERS:;
%put   Patients in either cohort: &n_after_age_cob;
%put   Patients in 2022 cohort: &n_in_2022;
%put   Patients in 2023 cohort: &n_in_2023;
%put;

%track_step(step=2, year=2022, desc=After enrollment (12 consecutive months with RX), 
            dataset=work.cohort_eligible(where=(in_cohort_2022=1)));

%track_step(step=2, year=2023, desc=After enrollment (10 consecutive months with RX), 
            dataset=work.cohort_eligible(where=(in_cohort_2023=1)));

/* Exclusion reasons */
title '2022 Exclusion Reasons (Among Enrolled)';
proc freq data=work.cohort_enrollment;
    where eligible_2022=1 and in_cohort_2022 ne 1;
    tables excl_reason_2022 / nocum;
run;

title '2023 Exclusion Reasons (Among Enrolled)';
proc freq data=work.cohort_enrollment;
    where eligible_2023=1 and in_cohort_2023 ne 1;
    tables excl_reason_2023 / nocum;
run;
title;

/*=============================================================================*/
/* STEP 5: IDENTIFY T2D PATIENTS WITH CLAIM TRACKING                         */
/*=============================================================================*/
%put =============================================================================;
%put STEP 5: IDENTIFYING TYPE 2 DIABETES PATIENTS;
%put =============================================================================;

/* Create patient lists */
proc sql;
    create table work.pts_2022 as
    select distinct pat_id from work.cohort_eligible where in_cohort_2022 = 1;
    
    create table work.pts_2023 as
    select distinct pat_id from work.cohort_eligible where in_cohort_2023 = 1;
quit;

/* Macro to find T2D with claim tracking */
%macro find_t2d(year=);
    
    %put   Finding T2D diagnoses in &year...;
    
    /* Outpatient - count claims */
    data work.op_&year._claims;
        if _n_ = 1 then do;
            declare hash h(dataset: "work.pts_&year");
            h.definekey('pat_id');
            h.definedone();
        end;
        
        %if &year = 2022 %then %do;
            set outpat.clm_outpat_22;
        %end;
        %else %do;
            set outpat.clm_outpat_23;
        %end;
        
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
        select count(*) into :n_op_claims_&year trimmed
        from work.op_&year._claims;
        
        select count(*) into :n_op_t2d_claims_&year trimmed
        from work.op_&year._claims where has_t2d=1;
    quit;
    
    /* Get unique patients */
    proc sql;
        create table work.op_&year as
        select distinct pat_id
        from work.op_&year._claims
        where has_t2d = 1;
    quit;
    
    /* Inpatient - count claims */
    data work.ip_&year._claims;
        if _n_ = 1 then do;
            declare hash h(dataset: "work.pts_&year");
            h.definekey('pat_id');
            h.definedone();
        end;
        
        %if &year = 2022 %then %do;
            set inpat.clm_inpat_22;
        %end;
        %else %do;
            set inpat.clm_inpat_23;
        %end;
        
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
        select count(*) into :n_ip_claims_&year trimmed
        from work.ip_&year._claims;
        
        select count(*) into :n_ip_t2d_claims_&year trimmed
        from work.ip_&year._claims where has_t2d=1;
    quit;
    
    /* Get unique patients */
    proc sql;
        create table work.ip_&year as
        select distinct pat_id
        from work.ip_&year._claims
        where has_t2d = 1;
    quit;
    
    /* Combine */
    data work.t2d_&year;
        set work.op_&year work.ip_&year;
    run;
    
    proc sort data=work.t2d_&year nodupkey;
        by pat_id;
    run;
    
    %put;
    %put T2D IDENTIFICATION - &year:;
    %put   Outpatient claims searched: &&n_op_claims_&year;
    %put   Outpatient claims with E11*: &&n_op_t2d_claims_&year;
    %put   Inpatient claims searched: &&n_ip_claims_&year;
    %put   Inpatient claims with E11*: &&n_ip_t2d_claims_&year;
    %put   Unique patients with T2D: %trim(%left(&sysnobs));
    %put;
    
    proc datasets library=work nolist;
        delete op_&year._claims ip_&year._claims;
    quit;
    
%mend;

%find_t2d(year=2022);
%find_t2d(year=2023);

%track_step(step=3, year=2022, desc=After T2D diagnosis requirement (E11*), 
            dataset=work.t2d_2022);

%track_step(step=3, year=2023, desc=After T2D diagnosis requirement (E11*), 
            dataset=work.t2d_2023);

/*=============================================================================*/
/* STEP 6: IDENTIFY INSULIN USERS WITH CLAIM TRACKING                        */
/*=============================================================================*/
%put =============================================================================;
%put STEP 6: IDENTIFYING INSULIN USERS;
%put =============================================================================;

/* Macro to find insulin users with tracking */
%macro find_insulin(year=);
    
    %put   Finding insulin users in &year...;
    
    data work.insulin_&year._claims;
        if _n_ = 1 then do;
            declare hash h(dataset: "work.t2d_&year");
            h.definekey('pat_id');
            h.definedone();
        end;
        
        %if &year = 2022 %then %do;
            set pharm.clm_rx_rtl_22;
        %end;
        %else %do;
            set pharm.clm_rx_rtl_23;
        %end;
        
        if h.check() = 0;
        
        /* Check for insulin */
        is_insulin = (substr(gpi14, 1, 4) = '2710');
        
        keep pat_id claimno linenum is_insulin;
    run;
    
    proc sql noprint;
        select count(*) into :n_rx_claims_&year trimmed
        from work.insulin_&year._claims;
        
        select count(*) into :n_insulin_claims_&year trimmed
        from work.insulin_&year._claims where is_insulin=1;
    quit;
    
    /* Get unique patients */
    proc sql;
        create table work.insulin_&year as
        select distinct pat_id
        from work.insulin_&year._claims
        where is_insulin = 1;
    quit;
    
    %put;
    %put INSULIN IDENTIFICATION - &year:;
    %put   Pharmacy claims searched: &&n_rx_claims_&year;
    %put   Insulin claims (GPI 2710*): &&n_insulin_claims_&year;
    %put   Unique patients with insulin: %trim(%left(&sysnobs));
    %put;
    
    proc datasets library=work nolist;
        delete insulin_&year._claims;
    quit;
    
%mend;

%find_insulin(year=2022);
%find_insulin(year=2023);

%track_step(step=4, year=2022, desc=Final cohort (after insulin requirement), 
            dataset=work.insulin_2022);

%track_step(step=4, year=2023, desc=Final cohort (after insulin requirement), 
            dataset=work.insulin_2023);

/*=============================================================================*/
/* STEP 7: CREATE FINAL COHORTS                                              */
/*=============================================================================*/
%put =============================================================================;
%put STEP 7: CREATING FINAL COHORTS;
%put =============================================================================;

/* Merge T2D and insulin flags */
proc sql;
    create table work.cohort_with_clinical as
    select c.*,
           case when t22.pat_id is not null then 1 else 0 end as has_t2d_2022,
           case when i22.pat_id is not null then 1 else 0 end as has_insulin_2022,
           case when t23.pat_id is not null then 1 else 0 end as has_t2d_2023,
           case when i23.pat_id is not null then 1 else 0 end as has_insulin_2023
    from work.cohort_eligible as c
    left join work.t2d_2022 as t22
        on c.pat_id = t22.pat_id
    left join work.insulin_2022 as i22
        on c.pat_id = i22.pat_id
    left join work.t2d_2023 as t23
        on c.pat_id = t23.pat_id
    left join work.insulin_2023 as i23
        on c.pat_id = i23.pat_id;
quit;

/* Apply T2D and insulin requirements with tracking */
data work.cohort_final;
    set work.cohort_with_clinical;
    
    /* Update exclusion reasons */
    if in_cohort_2022 = 1 then do;
        if has_t2d_2022 = 0 then do;
            in_cohort_2022 = 0;
            treatment_2022 = .;
            excl_reason_2022 = 'No T2D diagnosis';
        end;
        else if has_insulin_2022 = 0 then do;
            in_cohort_2022 = 0;
            treatment_2022 = .;
            excl_reason_2022 = 'No insulin use';
        end;
    end;
    
    if in_cohort_2023 = 1 then do;
        if has_t2d_2023 = 0 then do;
            in_cohort_2023 = 0;
            treatment_2023 = .;
            excl_reason_2023 = 'No T2D diagnosis';
        end;
        else if has_insulin_2023 = 0 then do;
            in_cohort_2023 = 0;
            treatment_2023 = .;
            excl_reason_2023 = 'No insulin use';
        end;
    end;
    
    /* Keep only if in at least one cohort */
    if in_cohort_2022 = 1 or in_cohort_2023 = 1;
run;

/* Final counts */
proc sql noprint;
    select count(distinct pat_id) into :n_final_2022 trimmed
    from work.cohort_final where in_cohort_2022=1;
    
    select count(distinct pat_id) into :n_final_2023 trimmed
    from work.cohort_final where in_cohort_2023=1;
    
    select count(distinct pat_id) into :n_final_control_2022 trimmed
    from work.cohort_final where in_cohort_2022=1 and treatment_2022=0;
    
    select count(distinct pat_id) into :n_final_treatment_2022 trimmed
    from work.cohort_final where in_cohort_2022=1 and treatment_2022=1;
    
    select count(distinct pat_id) into :n_final_control_2023 trimmed
    from work.cohort_final where in_cohort_2023=1 and treatment_2023=0;
    
    select count(distinct pat_id) into :n_final_treatment_2023 trimmed
    from work.cohort_final where in_cohort_2023=1 and treatment_2023=1;
quit;

%put;
%put FINAL COHORT COUNTS:;
%put   2022 Total: &n_final_2022;
%put   2022 Control (Commercial 18-64): &n_final_control_2022;
%put   2022 Treatment (MA 65+): &n_final_treatment_2022;
%put   2023 Total: &n_final_2023;
%put   2023 Control (Commercial 18-64): &n_final_control_2023;
%put   2023 Treatment (MA 65+): &n_final_treatment_2023;
%put;

/* Reasons for exclusion at T2D/insulin stage */
title '2022 Exclusions at T2D/Insulin Stage';
proc freq data=work.cohort_with_clinical;
    where eligible_2022=1 and payer_2022 ne '' and 
          ((payer_2022='C' and age_2022 between 18 and 64 and has_medicare_cob_2022=0) or
           (payer_2022='MA' and age_2022>=65));
    tables has_t2d_2022*has_insulin_2022 / nocol norow nopercent missing;
run;

title '2023 Exclusions at T2D/Insulin Stage';
proc freq data=work.cohort_with_clinical;
    where eligible_2023=1 and payer_2023 ne '' and 
          ((payer_2023='C' and age_2023 between 18 and 64 and has_medicare_cob_2023=0) or
           (payer_2023='MA' and age_2023>=65));
    tables has_t2d_2023*has_insulin_2023 / nocol norow nopercent missing;
run;
title;

/* Save final cohorts */
data output.cohort_cross_sectional;
    set work.cohort_final;
run;

/*=============================================================================*/
/* STEP 8: COMPREHENSIVE SUMMARY STATISTICS                                  */
/*=============================================================================*/
%put =============================================================================;
%put STEP 8: COMPREHENSIVE SUMMARY STATISTICS;
%put =============================================================================;

/* Print consort tracking */
title 'CONSORT FLOW - COHORT CONSTRUCTION';
proc print data=work.consort_tracking noobs label;
    var step year description n_patients n_excluded;
    format n_patients n_excluded comma12.;
    label step = 'Step'
          year = 'Year'
          description = 'Description'
          n_patients = 'N Patients'
          n_excluded = 'N Excluded';
run;

/* Export consort tracking */
proc export data=work.consort_tracking
    outfile="R:/IQVIA PharMetrics Plus (2024) Members/WatsonWilliam/New_copay_project/consort_flowchart.csv"
    dbms=csv replace;
run;

/* Detailed breakdowns */
title 'CROSS-SECTIONAL 2022 COHORT - FINAL';
proc freq data=output.cohort_cross_sectional;
    where in_cohort_2022 = 1;
    tables treatment_2022 * payer_2022 / nocol norow nopercent missing;
    tables treatment_2022 * age_group_2022 / nocol norow nopercent;
    tables treatment_2022 * der_sex / nocol norow nopercent;
run;

proc means data=output.cohort_cross_sectional n mean std min p25 p50 p75 max;
    where in_cohort_2022 = 1;
    class treatment_2022;
    var age_2022;
run;

title 'CROSS-SECTIONAL 2023 COHORT - FINAL';
proc freq data=output.cohort_cross_sectional;
    where in_cohort_2023 = 1;
    tables treatment_2023 * payer_2023 / nocol norow nopercent missing;
    tables treatment_2023 * age_group_2023 / nocol norow nopercent;
    tables treatment_2023 * der_sex / nocol norow nopercent;
run;

proc means data=output.cohort_cross_sectional n mean std min p25 p50 p75 max;
    where in_cohort_2023 = 1;
    class treatment_2023;
    var age_2023;
run;

/* Comprehensive summary table */
proc sql;
    create table output.cohort_summary_detailed as
    
    /* 2022 Cohort */
    select 
        '2022 Cross-Sectional' as cohort length=30,
        'Overall' as group length=20,
        sum(in_cohort_2022) as n_total format=comma12.,
        sum(case when treatment_2022=0 then 1 else 0 end) as n_control format=comma12.,
        sum(case when treatment_2022=1 then 1 else 0 end) as n_treatment format=comma12.,
        round(mean(case when treatment_2022=0 then age_2022 else . end), 0.1) as age_control_mean format=5.1,
        round(mean(case when treatment_2022=1 then age_2022 else . end), 0.1) as age_treatment_mean format=5.1,
        sum(case when treatment_2022=0 and der_sex='M' then 1 else 0 end) as n_male_control format=comma12.,
        sum(case when treatment_2022=1 and der_sex='M' then 1 else 0 end) as n_male_treatment format=comma12.
    from output.cohort_cross_sectional
    where in_cohort_2022=1
    
    union all
    
    /* 2022 By Age Group */
    select 
        '2022 Cross-Sectional' as cohort,
        age_group_2022 as group,
        count(*) as n_total format=comma12.,
        sum(case when treatment_2022=0 then 1 else 0 end) as n_control format=comma12.,
        sum(case when treatment_2022=1 then 1 else 0 end) as n_treatment format=comma12.,
        round(mean(case when treatment_2022=0 then age_2022 else . end), 0.1) as age_control_mean format=5.1,
        round(mean(case when treatment_2022=1 then age_2022 else . end), 0.1) as age_treatment_mean format=5.1,
        sum(case when treatment_2022=0 and der_sex='M' then 1 else 0 end) as n_male_control format=comma12.,
        sum(case when treatment_2022=1 and der_sex='M' then 1 else 0 end) as n_male_treatment format=comma12.
    from output.cohort_cross_sectional
    where in_cohort_2022=1
    group by age_group_2022
    
    union all
    
    /* 2023 Cohort */
    select 
        '2023 Cross-Sectional' as cohort,
        'Overall' as group,
        sum(in_cohort_2023) as n_total format=comma12.,
        sum(case when treatment_2023=0 then 1 else 0 end) as n_control format=comma12.,
        sum(case when treatment_2023=1 then 1 else 0 end) as n_treatment format=comma12.,
        round(mean(case when treatment_2023=0 then age_2023 else . end), 0.1) as age_control_mean format=5.1,
        round(mean(case when treatment_2023=1 then age_2023 else . end), 0.1) as age_treatment_mean format=5.1,
        sum(case when treatment_2023=0 and der_sex='M' then 1 else 0 end) as n_male_control format=comma12.,
        sum(case when treatment_2023=1 and der_sex='M' then 1 else 0 end) as n_male_treatment format=comma12.
    from output.cohort_cross_sectional
    where in_cohort_2023=1
    
    union all
    
    /* 2023 By Age Group */
    select 
        '2023 Cross-Sectional' as cohort,
        age_group_2023 as group,
        count(*) as n_total format=comma12.,
        sum(case when treatment_2023=0 then 1 else 0 end) as n_control format=comma12.,
        sum(case when treatment_2023=1 then 1 else 0 end) as n_treatment format=comma12.,
        round(mean(case when treatment_2023=0 then age_2023 else . end), 0.1) as age_control_mean format=5.1,
        round(mean(case when treatment_2023=1 then age_2023 else . end), 0.1) as age_treatment_mean format=5.1,
        sum(case when treatment_2023=0 and der_sex='M' then 1 else 0 end) as n_male_control format=comma12.,
        sum(case when treatment_2023=1 and der_sex='M' then 1 else 0 end) as n_male_treatment format=comma12.
    from output.cohort_cross_sectional
    where in_cohort_2023=1
    group by age_group_2023
    
    order by cohort, group;
quit;

title 'DETAILED COHORT SUMMARY';
proc print data=output.cohort_summary_detailed noobs;
run;
title;

/* Export all summary files */
proc export data=output.cohort_summary_detailed
    outfile="R:/IQVIA PharMetrics Plus (2024) Members/WatsonWilliam/New_copay_project/cohort_summary_detailed.csv"
    dbms=csv replace;
run;

proc export data=output.cohort_cross_sectional
    outfile="R:/IQVIA PharMetrics Plus (2024) Members/WatsonWilliam/New_copay_project/cohort_cross_sectional.csv"
    dbms=csv replace;
run;

/*=============================================================================*/
/* COMPLETION WITH COMPREHENSIVE LOG                                         */
/*=============================================================================*/

%let end_time = %sysfunc(datetime());
%let elapsed = %sysevalf(&end_time - &start_time);

%put;
%put =============================================================================;
%put CROSS-SECTIONAL COHORTS COMPLETE - COMPREHENSIVE TRACKING;
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
%put FINAL COHORT SUMMARY:;
%put   2022 Cohort Total: &n_final_2022;
%put     - Control (Commercial 18-64): &n_final_control_2022;
%put     - Treatment (MA 65+): &n_final_treatment_2022;
%put   2023 Cohort Total: &n_final_2023;
%put     - Control (Commercial 18-64): &n_final_control_2023;
%put     - Treatment (MA 65+): &n_final_treatment_2023;
%put;
%put OUTPUT FILES CREATED:;
%put   1. output.cohort_cross_sectional - Master cohort file;
%put   2. output.cohort_summary_detailed - Detailed summary statistics;
%put   3. consort_flowchart.csv - Step-by-step inclusion/exclusion tracking;
%put   4. cohort_summary_detailed.csv - Detailed summary export;
%put   5. cohort_cross_sectional.csv - Full cohort export;
%put;
%put =============================================================================;

/* Clean up */
proc datasets library=work nolist kill;
quit;
