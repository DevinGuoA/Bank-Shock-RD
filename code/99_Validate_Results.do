version 18.0
clear all
set more off
set varabbrev off

if `"$HW2_ROOT"' == "" {
    local launch = subinstr(`"`c(pwd)'"', "\\", "/", .)
    capture confirm file `"`launch'/tmp/results/Question_1_Results.json"'
    if !_rc global HW2_ROOT `"`launch'"'
    else {
        capture confirm file `"`launch'/HW2/tmp/results/Question_1_Results.json"'
        if !_rc global HW2_ROOT `"`launch'/HW2"'
        else {
            capture confirm file `"`launch'/../tmp/results/Question_1_Results.json"'
            if !_rc {
                quietly cd `"`launch'/.."'
                global HW2_ROOT `"`c(pwd)'"'
                quietly cd `"`launch'"'
            }
        }
    }
}
capture confirm file `"$HW2_ROOT/tmp/results/Question_1_Results.json"'
if _rc exit 601
capture mkdir `"$HW2_ROOT/tmp/results/stata"'

* Confirm that every core Stata artifact was produced before comparing values.
local required_artifacts ///
    data/processed/stata/Loans_Cleaned_Stata.csv ///
    data/processed/stata/Loans_Cleaned_Stata.dta ///
    data/processed/stata/Banks_Cleaned_Stata.csv ///
    data/processed/stata/Banks_Cleaned_Stata.dta ///
    data/processed/stata/Firms_Cleaned_Stata.csv ///
    data/processed/stata/Firms_Cleaned_Stata.dta ///
    output/stata/Figure_1_Distributions_Stata.png ///
    output/stata/Figure_1_2_Banks_Per_Firm_Stata.png ///
    output/stata/Figure_2_1_RD_Fits_Stata.png ///
    output/stata/Figure_2_2_Asset_Density_Stata.png ///
    output/stata/Figure_2_3_Bank_Characteristics_Stata.png ///
    output/stata/Figure_4_Magnitude_Comparison_Stata.png ///
    output/stata/Table_1_1_Summary_Stata.csv ///
    output/stata/Table_1_2_Sample_Stata.csv ///
    output/stata/Table_2_1_RD_Windows_Stata.csv ///
    output/stata/Table_2_2_RD_Robust_Stata.csv ///
    output/stata/Table_2_3_Validity_Checks_Stata.csv ///
    output/stata/Table_3_1_Loan_Regressions_Stata.csv ///
    output/stata/Table_3_2_Credit_Quality_Sorting_Stata.csv ///
    output/stata/Table_4_1_Firm_Regressions_Stata.csv ///
    tmp/results/stata/Question_1_Results_Stata.csv ///
    tmp/results/stata/Question_2_Results_Stata.csv ///
    tmp/results/stata/Question_3_Results_Stata.csv ///
    tmp/results/stata/Question_4_Results_Stata.csv
tempfile artifactcheck
tempname artifactpost
postfile `artifactpost' str100 artifact byte exists using `artifactcheck', replace
local missing_artifacts = 0
foreach relative of local required_artifacts {
    capture confirm file `"$HW2_ROOT/`relative'"'
    local exists = (_rc==0)
    post `artifactpost' ("`relative'") (`exists')
    if !`exists' local ++missing_artifacts
}
postclose `artifactpost'
if `missing_artifacts' {
    display as error "`missing_artifacts' required Stata artifacts are missing."
    exit 601
}

* JSON parsing only: all estimators above are native Stata/Mata.  The Python
* standard library converts the already-saved Python JSONs into a long baseline.
python:
import csv
import json
import os
from itertools import zip_longest
from sfi import Macro

root = Macro.getGlobal("HW2_ROOT")
out = os.path.join(root, "tmp", "results", "stata", "Python_JSON_Baseline.csv")
rows = []
processed_rows = []

def flatten(value, prefix, question):
    if isinstance(value, bool):
        return
    if isinstance(value, (int, float)):
        rows.append((question, prefix, value))
    elif isinstance(value, dict):
        for key, item in value.items():
            name = f"{prefix}.{key}" if prefix else str(key)
            flatten(item, name, question)
    elif isinstance(value, list):
        for index, item in enumerate(value):
            flatten(item, f"{prefix}[{index}]", question)

for question in range(1, 5):
    path = os.path.join(root, "tmp", "results", f"Question_{question}_Results.json")
    with open(path, "r", encoding="utf-8") as handle:
        flatten(json.load(handle), "", question)

with open(out, "w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle)
    writer.writerow(["question", "metric", "python_value"])
    writer.writerows(rows)

# Cell-by-cell checks for the three cleaned datasets created in Question 1.
processed_out = os.path.join(root, "tmp", "results", "stata", "Processed_Data_Validation_Stata.csv")

for stem in ("Loans_Cleaned", "Banks_Cleaned", "Firms_Cleaned"):
    python_path = os.path.join(root, "data", "processed", f"{stem}.csv")
    stata_path = os.path.join(root, "data", "processed", "stata",
                              f"{stem}_Stata.csv")
    python_n = stata_n = cells_over = 0
    max_diff = 0.0
    with open(python_path, newline="", encoding="utf-8") as p_handle, \
         open(stata_path, newline="", encoding="utf-8") as s_handle:
        p_reader = csv.DictReader(p_handle)
        s_reader = csv.DictReader(s_handle)
        columns_match = int(p_reader.fieldnames == s_reader.fieldnames)
        for p_row, s_row in zip_longest(p_reader, s_reader):
            if p_row is not None:
                python_n += 1
            if s_row is not None:
                stata_n += 1
            if p_row is None or s_row is None or not columns_match:
                continue
            for column in p_reader.fieldnames:
                difference = abs(float(p_row[column])-float(s_row[column]))
                max_diff = max(max_diff, difference)
                cells_over += int(difference > 1e-12)
    passed = int(columns_match and python_n == stata_n and cells_over == 0)
    processed_rows.append((stem, python_n, stata_n, columns_match, max_diff,
                           cells_over, 1e-12, passed))

with open(processed_out, "w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle)
    writer.writerow(["dataset", "python_rows", "stata_rows", "columns_match",
                     "max_abs_diff", "cells_over_tolerance", "tolerance", "pass"])
    writer.writerows(processed_rows)
end

tempfile allstata
forvalues q=1/4 {
    capture confirm file `"$HW2_ROOT/tmp/results/stata/Question_`q'_Results_Stata.csv"'
    if _rc {
        display as error "Missing Stata result CSV for Question `q'."
        exit 601
    }
    import delimited using ///
        `"$HW2_ROOT/tmp/results/stata/Question_`q'_Results_Stata.csv"', ///
        varnames(1) case(preserve) stringcols(1) asdouble clear
    generate byte question = `q'
    rename value stata_value
    keep question metric stata_value
    if `q' == 1 save `allstata', replace
    else {
        append using `allstata'
        save `allstata', replace
    }
}

import delimited using `"$HW2_ROOT/tmp/results/stata/Python_JSON_Baseline.csv"', ///
    varnames(1) case(preserve) stringcols(2) asdouble clear
merge 1:1 question metric using `allstata'
generate str18 merge_status = cond(_merge==1,"Python only", ///
    cond(_merge==2,"Stata only","matched"))

* Every numeric JSON leaf (excluding JSON booleans) must have a Stata analogue.
quietly count if _merge!=3
if r(N) {
    sort question metric
    format python_value stata_value %24.17g
    export delimited using ///
        `"$HW2_ROOT/tmp/results/stata/Python_Stata_Validation_Stata.csv"', replace
    list question metric merge_status if _merge!=3, noobs abbreviate(60)
    display as error "Validation mapping is incomplete."
    exit 459
}
drop _merge merge_status

generate double abs_diff = abs(stata_value-python_value)
generate double tolerance = 1e-9
replace tolerance = .002 if question==2 & ///
    inlist(metric,"density_test.h_left","density_test.h_right")
replace tolerance = 1e-5 if question==2 & ///
    inlist(metric,"density_test.estimate_left","density_test.estimate_right", ///
    "density_test.estimate_difference","density_test.standard_error")
replace tolerance = 5e-4 if question==2 & ///
    inlist(metric,"density_test.z_statistic","density_test.p_value")
generate byte pass = abs_diff<=tolerance
generate str16 tolerance_basis = "strict"
replace tolerance_basis = "density port" if question==2 & ///
    strpos(metric,"density_test.")==1 & tolerance>1e-9
format python_value stata_value abs_diff tolerance %24.17g
sort question metric
export delimited using ///
    `"$HW2_ROOT/tmp/results/stata/Python_Stata_Validation_Stata.csv"', replace

preserve
    collapse (count) checks=python_value (sum) failures=pass (max) max_abs_diff=abs_diff, ///
        by(question)
    replace failures = checks-failures
    format max_abs_diff %24.17g
    export delimited using ///
        `"$HW2_ROOT/tmp/results/stata/Validation_Summary_Stata.csv"', replace
restore

quietly count if !pass
local failures = r(N)
quietly count
local checks = r(N)
if `failures' {
    list question metric stata_value python_value abs_diff tolerance if !pass, ///
        noobs abbreviate(60)
    display as error "`failures' of `checks' Python/Stata checks failed."
    exit 459
}

display as result "All `checks' Python/Stata numerical checks passed."
display as text "Density-only tolerances reflect rddensity Stata 3.0 versus Python 2.4.6; all reported three-decimal values coincide."

preserve
    import delimited using ///
        `"$HW2_ROOT/tmp/results/stata/Processed_Data_Validation_Stata.csv"', ///
        varnames(1) case(preserve) stringcols(1) asdouble clear
    quietly count if pass!=1
    if r(N) {
        list, noobs abbreviate(40)
        display as error "Cleaned-data comparison failed."
        exit 459
    }
    quietly count
    local datasets = r(N)
restore
display as result "All `datasets' cleaned datasets match Python cell by cell."

preserve
    use `artifactcheck', clear
    export delimited using ///
        `"$HW2_ROOT/tmp/results/stata/Artifact_Manifest_Stata.csv"', replace
    quietly count
    local artifacts = r(N)
restore
display as result "All `artifacts' required Stata artifacts are present."
