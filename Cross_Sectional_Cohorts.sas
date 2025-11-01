/*=============================================================================*/
/* IRA INSULIN COPAY CAP ANALYSIS - CROSS-SECTIONAL COHORTS                  */
/* CORRECTED VERSION - Uses pay_type + ben_rx (NOT mstr_enroll_cd)           */
/*=============================================================================*/

libname inpat "R:/IQVIA PharMetrics Plus (2024)/Claims - Medical - Inpatient";
libname outpat "R:/IQVIA PharMetrics Plus (2024)/Claims - Medical - Outpatient";
libname pharm "R:/IQVIA PharMetrics Plus (2024)/Claims - Retail Pharmacy";
libname enroll "R:/IQVIA PharMetrics Plus (2024)/Enrollment";
libname output "R:/IQVIA PharMetrics Plus (2024) Members/WatsonWilliam/New_copay_project";

options compress=yes sortsize=8G sumsize=8G memsize=32G threads cpucount=6 fullstimer;

%let start_time = %sysfunc(datetime());

%put =============================================================================;
%put CROSS-SECTIONAL COHORTS - 2022 AND 2023;
%put =============================================================================;

/*=============================================================================*/
/* STEP 1: LOAD BASE DEMOGRAPHICS                                            */
/*=============================================================================*/
data work.demographics;
    set enroll.enroll_synth;
    
    /* Calculate ages - CORRECTED top-coding */
    if der_yob = 0 then do;
        age_2022 = 86;  /* CORRECTED: 0 = 86+ */
        age_2023 = 86;
    end;
    else if der_yob > 0 then do;
        age_2022 = 2022 - der_yob;
        age_2023 = 2023 - der_yob;
    end;
    else delete;
    
    /* Base eligibility */
    if (age_2022 >= 18 or age_2023 >= 18) and age_2022 <= 110 and age_2023 <= 110;
    
    /* Exclude states with insulin copay caps */
    if pat_state not in ('CA','CO','CT','DE','IL','ME','NH','NM',
                         'NY','UT','VA','WA','WV');
    
    /* Create age groups */
    length age_group_2022 $10 age_group_2023 $10;
    if age_2022 >= 18 and age_2022 <= 61 then age_group_2022 = '18-61';
    else if age_2022 >= 62 and age_2022 <= 64 then age_group_2022 = '62-64';
    else if age_2022 >= 65 then age_group_2022 = '65+';
    
    if age_2023 >= 18 and age_2023 <= 61 then age_group_2023 = '18-61';
    else if age_2023 >= 62 and age_2023 <= 64 then age_group_2023 = '62-64';
    else if age_2023 >= 65 then age_group_2023 = '65+';
    
    keep pat_id der_sex pat_state der_yob age_2022 age_2023 age_group_2022 age_group_2023;
run;

%put NOTE: Demographics loaded - N=%trim(%left(&sysnobs));

/*=============================================================================*/
/* STEP 2: EXTRACT ENROLLMENT STRINGS (CORRECTED)                            */
/*=============================================================================*/
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
    
    /* 2023: positions 265-276 (only first 10 have data) */
    pay_2023 = substr(pay_string, 265, 10);
    rx_2023 = substr(rx_string, 265, 10);
    mcob_2023 = substr(mcob_string, 265, 10);
    pcob_2023 = substr(pcob_string, 265, 10);
    
    keep pat_id pay_2022 rx_2022 mcob_2022 pcob_2022
         pay_2023 rx_2023 mcob_2023 pcob_2023;
run;

%put NOTE: Enrollment periods extracted;

/*=============================================================================*/
/* STEP 3: CALCULATE ENROLLMENT ELIGIBILITY                                  */
/*=============================================================================*/
data work.enrollment_eligible;
    set work.enrollment;
    
    /*=========================================================================*/
    /* 2022 ELIGIBILITY (12 consecutive months)                              */
    /*=========================================================================*/
    length eligible_2022 8 payer_2022 $15 has_medicare_cob_2022 8;
    eligible_2022 = 0;
    payer_2022 = '';
    has_medicare_cob_2022 = 0;
    
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
        
        /* Check for Medicare COB */
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
    
    /* Check for 12 consecutive months with payer AND Rx coverage */
    if payer_2022 ne '' then do;
        consec = 0;
        max_consec = 0;
        
        do i = 1 to 12;
            /* Commercial needs Rx, MA doesn't (Part D is optional) */
            if payer_2022 = 'C' then do;
                if pay22{i} = 'C' and rx22{i} = 'Y' then do;
                    consec + 1;
                    max_consec = max(max_consec, consec);
                end;
                else consec = 0;
            end;
            else if payer_2022 = 'MA' then do;
                /* MA: Check if they have Rx coverage */
                if pay22{i} in ('A','R') and rx22{i} = 'Y' then do;
                    consec + 1;
                    max_consec = max(max_consec, consec);
                end;
                else consec = 0;
            end;
        end;
        
        if max_consec >= 12 then eligible_2022 = 1;
    end;
    
    /*=========================================================================*/
    /* 2023 ELIGIBILITY (10 consecutive months)                              */
    /*=========================================================================*/
    length eligible_2023 8 payer_2023 $15 has_medicare_cob_2023 8;
    eligible_2023 = 0;
    payer_2023 = '';
    has_medicare_cob_2023 = 0;
    
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
        
        /* Check for Medicare COB */
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
    
    /* Check for 10 consecutive months with payer AND Rx coverage */
    if payer_2023 ne '' then do;
        consec = 0;
        max_consec = 0;
        
        do i = 1 to 10;
            /* Commercial needs Rx, MA needs Rx too (Part D) */
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
    end;
    
    keep pat_id eligible_2022 payer_2022 has_medicare_cob_2022
         eligible_2023 payer_2023 has_medicare_cob_2023;
    drop i n_c n_ma consec max_consec;
run;

%put NOTE: Enrollment eligibility calculated;

/* Add diagnostics */
proc sql noprint;
    select count(*) into :n_elig_2022 from work.enrollment_eligible where eligible_2022=1;
    select count(*) into :n_elig_2023 from work.enrollment_eligible where eligible_2023=1;
quit;

%put NOTE: Patients eligible for 2022 enrollment: &n_elig_2022;
%put NOTE: Patients eligible for 2023 enrollment: &n_elig_2023;

/*=============================================================================*/
/* STEP 4: MERGE WITH DEMOGRAPHICS AND APPLY AGE/COB FILTERS                 */
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
    
    /* 2022 COHORT */
    length in_cohort_2022 8 treatment_2022 8;
    in_cohort_2022 = 0;
    treatment_2022 = .;
    
    if eligible_2022 = 1 then do;
        /* Commercial 18-64 (Control) */
        if payer_2022 = 'C' and age_2022 >= 18 and age_2022 <= 64 then do;
            if has_medicare_cob_2022 = 0 then do;
                in_cohort_2022 = 1;
                treatment_2022 = 0;
            end;
        end;
        /* Medicare Advantage 65+ (Treatment) */
        else if payer_2022 = 'MA' and age_2022 >= 65 then do;
            in_cohort_2022 = 1;
            treatment_2022 = 1;
        end;
    end;
    
    /* 2023 COHORT */
    length in_cohort_2023 8 treatment_2023 8;
    in_cohort_2023 = 0;
    treatment_2023 = .;
    
    if eligible_2023 = 1 then do;
        /* Commercial 18-64 (Control) */
        if payer_2023 = 'C' and age_2023 >= 18 and age_2023 <= 64 then do;
            if has_medicare_cob_2023 = 0 then do;
                in_cohort_2023 = 1;
                treatment_2023 = 0;
            end;
        end;
        /* Medicare Advantage 65+ (Treatment) */
        else if payer_2023 = 'MA' and age_2023 >= 65 then do;
            in_cohort_2023 = 1;
            treatment_2023 = 1;
        end;
    end;
    
    /* Keep only if eligible for at least one cohort */
    if in_cohort_2022 = 1 or in_cohort_2023 = 1;
run;

%put NOTE: Age and COB filters applied - N=%trim(%left(&sysnobs));

proc freq data=work.cohort_eligible;
    tables in_cohort_2022 in_cohort_2023;
    title 'Patients eligible by year (before T2D/insulin filters)';
run;

/*=============================================================================*/
/* STEP 5: IDENTIFY T2D PATIENTS                                             */
/*=============================================================================*/
%put =============================================================================;
%put STEP 5: IDENTIFYING TYPE 2 DIABETES PATIENTS;
%put =============================================================================;

/* Create patient lists for each year */
proc sql;
    create table work.pts_2022 as
    select distinct pat_id from work.cohort_eligible where in_cohort_2022 = 1;
    
    create table work.pts_2023 as
    select distinct pat_id from work.cohort_eligible where in_cohort_2023 = 1;
quit;

%put NOTE: 2022 patient list: N=%trim(%left(&sqlobs));

/* Macro to find T2D in medical claims */
%macro find_t2d(year=);
    
    %put   Finding T2D diagnoses in &year...;
    
    /* Outpatient */
    data work.op_&year;
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
        
        /* Check for E11* in any diagnosis field */
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
    
    /* Inpatient */
    data work.ip_&year;
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
        
        if has_t2d = 1;
        keep pat_id;
    run;
    
    /* Combine and deduplicate */
    data work.t2d_&year;
        set work.op_&year work.ip_&year;
    run;
    
    proc sort data=work.t2d_&year nodupkey;
        by pat_id;
    run;
    
    %put NOTE: T2D patients in &year: N=%trim(%left(&sysnobs));
    
    proc datasets library=work nolist;
        delete op_&year ip_&year;
    quit;
    
%mend;

%find_t2d(year=2022);
%find_t2d(year=2023);

/*=============================================================================*/
/* STEP 6: IDENTIFY INSULIN USERS                                            */
/*=============================================================================*/
%put =============================================================================;
%put STEP 6: IDENTIFYING INSULIN USERS;
%put =============================================================================;

/* Macro to find insulin users */
%macro find_insulin(year=);
    
    %put   Finding insulin users in &year...;
    
    data work.insulin_&year;
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
        
        /* Check for insulin (GPI 2710*) */
        if substr(gpi14, 1, 4) = '2710';
        
        keep pat_id;
    run;
    
    proc sort data=work.insulin_&year nodupkey;
        by pat_id;
    run;
    
    %put NOTE: Insulin users in &year: N=%trim(%left(&sysnobs));
    
%mend;

%find_insulin(year=2022);
%find_insulin(year=2023);

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

/* Apply T2D and insulin requirements */
data work.cohort_final;
    set work.cohort_with_clinical;
    
    /* 2022 cohort: Need T2D AND insulin in 2022 */
    if in_cohort_2022 = 1 and (has_t2d_2022 = 0 or has_insulin_2022 = 0) then do;
        in_cohort_2022 = 0;
        treatment_2022 = .;
    end;
    
    /* 2023 cohort: Need T2D AND insulin in 2023 */
    if in_cohort_2023 = 1 and (has_t2d_2023 = 0 or has_insulin_2023 = 0) then do;
        in_cohort_2023 = 0;
        treatment_2023 = .;
    end;
    
    /* Keep only if in at least one cohort */
    if in_cohort_2022 = 1 or in_cohort_2023 = 1;
run;

%put NOTE: Final cohorts created - N=%trim(%left(&sysnobs));

/* Save final cohorts */
data output.cohort_cross_sectional;
    set work.cohort_final;
run;

/*=============================================================================*/
/* STEP 8: SUMMARY STATISTICS                                                */
/*=============================================================================*/
%put =============================================================================;
%put STEP 8: COHORT SUMMARY STATISTICS;
%put =============================================================================;

title 'CROSS-SECTIONAL 2022 COHORT';
proc freq data=output.cohort_cross_sectional;
    where in_cohort_2022 = 1;
    tables treatment_2022 * payer_2022 / nocol norow nopercent missing;
run;

proc freq data=output.cohort_cross_sectional;
    where in_cohort_2022 = 1;
    tables treatment_2022 * age_group_2022 / nocol norow nopercent;
run;

proc means data=output.cohort_cross_sectional n mean std min p50 max;
    where in_cohort_2022 = 1;
    class treatment_2022;
    var age_2022;
run;

title 'CROSS-SECTIONAL 2023 COHORT';
proc freq data=output.cohort_cross_sectional;
    where in_cohort_2023 = 1;
    tables treatment_2023 * payer_2023 / nocol norow nopercent missing;
run;

proc freq data=output.cohort_cross_sectional;
    where in_cohort_2023 = 1;
    tables treatment_2023 * age_group_2023 / nocol norow nopercent;
run;

proc means data=output.cohort_cross_sectional n mean std min p50 max;
    where in_cohort_2023 = 1;
    class treatment_2023;
    var age_2023;
run;

/* Summary table */
proc sql;
    create table output.cohort_summary as
    select 
        '2022 Cross-Sectional' as cohort length=30,
        sum(in_cohort_2022) as n_total format=comma12.,
        sum(case when treatment_2022=0 then 1 else 0 end) as n_control format=comma12.,
        sum(case when treatment_2022=1 then 1 else 0 end) as n_treatment format=comma12.,
        sum(case when treatment_2022=0 and age_group_2022='62-64' then 1 else 0 end) 
            as n_control_62_64 format=comma12. label='Control 62-64'
    from output.cohort_cross_sectional
    
    union all
    
    select 
        '2023 Cross-Sectional' as cohort,
        sum(in_cohort_2023) as n_total format=comma12.,
        sum(case when treatment_2023=0 then 1 else 0 end) as n_control format=comma12.,
        sum(case when treatment_2023=1 then 1 else 0 end) as n_treatment format=comma12.,
        sum(case when treatment_2023=0 and age_group_2023='62-64' then 1 else 0 end) 
            as n_control_62_64 format=comma12.
    from output.cohort_cross_sectional;
quit;

title 'COHORT SUMMARY';
proc print data=output.cohort_summary noobs;
run;

/* Export */
proc export data=output.cohort_summary
    outfile="R:/IQVIA PharMetrics Plus (2024) Members/WatsonWilliam/New_copay_project/cohort_summary_cross_sectional.csv"
    dbms=csv replace;
run;

proc export data=output.cohort_cross_sectional
    outfile="R:/IQVIA PharMetrics Plus (2024) Members/WatsonWilliam/New_copay_project/cohort_cross_sectional.csv"
    dbms=csv replace;
run;

/*=============================================================================*/
/* COMPLETION                                                                 */
/*=============================================================================*/

%let end_time = %sysfunc(datetime());
%let elapsed = %sysevalf(&end_time - &start_time);

%put;
%put =============================================================================;
%put CROSS-SECTIONAL COHORTS COMPLETE;
%put =============================================================================;
%put Elapsed time: %sysfunc(putn(&elapsed, time12.2));
%put;
%put OUTPUT FILES CREATED:;
%put   output.cohort_cross_sectional - Master cohort file;
%put   output.cohort_summary - Summary statistics;
%put;

/* Clean up work library */
proc datasets library=work nolist kill;
quit;
