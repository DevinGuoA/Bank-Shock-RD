version 18.0
clear all
set more off
set varabbrev off

* Locate HW2 from the usual launch points: project root, HW2 root, or code/.
if `"$HW2_ROOT"' == "" {
    local launch = subinstr(`"`c(pwd)'"', "\\", "/", .)
    capture confirm file `"`launch'/data/PhdFinance_FNCE7020_HW2_Data.csv"'
    if !_rc global HW2_ROOT `"`launch'"'
    else {
        capture confirm file `"`launch'/HW2/data/PhdFinance_FNCE7020_HW2_Data.csv"'
        if !_rc global HW2_ROOT `"`launch'/HW2"'
        else {
            capture confirm file `"`launch'/../data/PhdFinance_FNCE7020_HW2_Data.csv"'
            if !_rc {
                quietly cd `"`launch'/.."'
                global HW2_ROOT `"`c(pwd)'"'
                quietly cd `"`launch'"'
            }
        }
    }
}

capture confirm file `"$HW2_ROOT/data/PhdFinance_FNCE7020_HW2_Data.csv"'
if _rc {
    display as error "Could not locate HW2. Run from the project root, HW2, or HW2/code."
    exit 601
}

capture mkdir `"$HW2_ROOT/output/stata"'
capture mkdir `"$HW2_ROOT/data/processed/stata"'
capture mkdir `"$HW2_ROOT/tmp/results/stata"'
capture mkdir `"$HW2_ROOT/tmp/stata_ado_hw2_20260527"'

capture log close _all
log using `"$HW2_ROOT/tmp/results/stata/00_Main_Stata.log"', text replace name(hw2main)

foreach q in 1 2 3 4 {
    local script : display %02.0f `q'
    display as text _newline "Running `script'_Question_`q'.do"
    do `"$HW2_ROOT/code/`script'_Question_`q'.do"'
}

display as text _newline "Running 99_Validate_Results.do"
do `"$HW2_ROOT/code/99_Validate_Results.do"'

display as result _newline "All HW2 Stata analyses and Python/Stata checks completed."
log close hw2main
