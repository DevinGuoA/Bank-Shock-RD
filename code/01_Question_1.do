version 18.0
clear all
set more off
set varabbrev off

* Standalone root discovery (00_Main.do sets this global in normal use).
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
if _rc exit 601
capture mkdir `"$HW2_ROOT/output/stata"'
capture mkdir `"$HW2_ROOT/data/processed/stata"'
capture mkdir `"$HW2_ROOT/tmp/results/stata"'

import delimited using `"$HW2_ROOT/data/PhdFinance_FNCE7020_HW2_Data.csv"', ///
    varnames(1) case(preserve) asdouble clear

local sourcevars LoanID FirmID BankID Char1_BankAssets Char2_BankSupervised ///
    Treatment_BankShock Char3_PreCreditShare Char4_FirmCreditQuality ///
    Outcome_LoanGrowth Outcome_FirmCreditGrowth Outcome_FirmInvestmentGrowth

* Structural checks: the source is a unique firm-bank relationship file.
egen long __rowmiss = rowmiss(`sourcevars')
quietly summarize __rowmiss, meanonly
scalar q1_missing = r(sum)
assert q1_missing == 0
drop __rowmiss
isid LoanID
isid FirmID BankID
bysort BankID: assert Char1_BankAssets == Char1_BankAssets[1]
bysort BankID: assert Char2_BankSupervised == Char2_BankSupervised[1]
bysort BankID: assert Treatment_BankShock == Treatment_BankShock[1]
bysort FirmID: assert Char4_FirmCreditQuality == Char4_FirmCreditQuality[1]
bysort FirmID: assert Outcome_FirmCreditGrowth == Outcome_FirmCreditGrowth[1]
bysort FirmID: assert Outcome_FirmInvestmentGrowth == Outcome_FirmInvestmentGrowth[1]
assert inlist(Char2_BankSupervised, 0, 1)
assert inlist(Treatment_BankShock, 0, 1)
assert Char2_BankSupervised == (Char1_BankAssets >= 30)
assert inrange(Char3_PreCreditShare, 0, 1)

quietly count
local loans = r(N)
quietly levelsof BankID, local(__banks)
local banks_n : word count `__banks'
quietly levelsof FirmID, local(__firms)
local firms_n : word count `__firms'
quietly count if Char3_PreCreditShare == 0
local zero_shares = r(N)

bysort FirmID: generate long LoanCount = _N
bysort FirmID: egen double __sharesum = total(Char3_PreCreditShare)
assert __sharesum > 0
egen byte __firmtag = tag(FirmID)
quietly count if __firmtag & LoanCount >= 2
local firms_multi = r(N)
quietly count if __firmtag & LoanCount == 1
local firms_one = r(N)
quietly summarize LoanCount if __firmtag, meanonly
local max_banks = r(max)
quietly count if __firmtag & abs(__sharesum - 1) > 1e-12
local share_round_firms = r(N)
generate double __share_error = abs(__sharesum - 1)
quietly summarize __share_error if __firmtag, meanonly
local max_share_error = r(max)
quietly summarize __sharesum if __firmtag, meanonly
local min_share_sum = r(min)
local max_share_sum = r(max)

generate double PreCreditShareNormalized = Char3_PreCreditShare / __sharesum
bysort FirmID: egen double __normsum = total(PreCreditShareNormalized)
generate double __normerror = abs(__normsum - 1)
quietly summarize __normerror if __firmtag, meanonly
local max_norm_error = r(max)

quietly summarize Char2_BankSupervised, meanonly
local supervised_loan_share = r(mean)
quietly summarize Treatment_BankShock, meanonly
local shock_loan_share = r(mean)

sort LoanID
preserve
    keep `sourcevars' PreCreditShareNormalized LoanCount
    order `sourcevars' PreCreditShareNormalized LoanCount
    export delimited using `"$HW2_ROOT/data/processed/stata/Loans_Cleaned_Stata.csv"', replace
    save `"$HW2_ROOT/data/processed/stata/Loans_Cleaned_Stata.dta"', replace
restore

tempfile loansdata banksdata firmsdata bankcounts q1results
save `loansdata', replace

preserve
    collapse (firstnm) Char1_BankAssets Char2_BankSupervised Treatment_BankShock ///
        (mean) Mean_LoanGrowth=Outcome_LoanGrowth ///
        Mean_FirmCreditQuality=Char4_FirmCreditQuality ///
        (count) LoanCount=LoanID, by(BankID)
    quietly summarize Char2_BankSupervised, meanonly
    local supervised_bank_share = r(mean)
    quietly summarize Treatment_BankShock, meanonly
    local shock_bank_share = r(mean)
    quietly count if Char2_BankSupervised == 1
    local supervised_bank_count = r(N)
    quietly count if Treatment_BankShock == 1
    local shock_bank_count = r(N)
    save `banksdata', replace
    export delimited using `"$HW2_ROOT/data/processed/stata/Banks_Cleaned_Stata.csv"', replace
    save `"$HW2_ROOT/data/processed/stata/Banks_Cleaned_Stata.dta"', replace
restore

generate double __rawshock = Char3_PreCreditShare * Treatment_BankShock
generate double __normshock = PreCreditShareNormalized * Treatment_BankShock
preserve
    collapse (firstnm) Char4_FirmCreditQuality Outcome_FirmCreditGrowth ///
        Outcome_FirmInvestmentGrowth LoanCount ///
        (sum) Exposure_raw=__rawshock Exposure=__normshock, by(FirmID)
    save `firmsdata', replace
    export delimited using `"$HW2_ROOT/data/processed/stata/Firms_Cleaned_Stata.csv"', replace
    save `"$HW2_ROOT/data/processed/stata/Firms_Cleaned_Stata.dta"', replace
restore

* Pandas describe() uses the type-7 quantile definition.  Implement it directly.
capture mata: mata drop hw2_type7()
mata:
real rowvector hw2_type7(string scalar v)
{
    real colvector x
    real rowvector probs, out
    real scalar n, j, h, lo, frac
    x = sort(st_data(., v), 1)
    x = select(x, x :< .)
    n = rows(x)
    probs = (0.25, 0.50, 0.75)
    out = J(1, cols(probs), .)
    for (j=1; j<=cols(probs); j++) {
        h = (n-1)*probs[j] + 1
        lo = floor(h)
        frac = h-lo
        if (lo>=n) out[j] = x[n]
        else out[j] = x[lo] + frac*(x[lo+1]-x[lo])
    }
    return(out)
}
end

tempname q1
postfile `q1' str100 metric double value using `q1results', replace
post `q1' ("audit.loans") (`loans')
post `q1' ("audit.banks") (`banks_n')
post `q1' ("audit.firms") (`firms_n')
post `q1' ("audit.firms_two_or_more_banks") (`firms_multi')
post `q1' ("audit.firms_one_bank") (`firms_one')
post `q1' ("audit.max_banks_per_firm") (`max_banks')
post `q1' ("audit.missing_cells") (q1_missing)
post `q1' ("audit.duplicate_loan_ids") (0)
post `q1' ("audit.duplicate_firm_bank_pairs") (0)
post `q1' ("audit.supervision_rule_violations") (0)
post `q1' ("audit.zero_credit_shares") (`zero_shares')
post `q1' ("audit.firms_with_share_rounding_error") (`share_round_firms')
post `q1' ("audit.max_credit_share_sum_error") (`max_share_error')
post `q1' ("audit.min_credit_share_sum") (`min_share_sum')
post `q1' ("audit.max_credit_share_sum") (`max_share_sum')
post `q1' ("audit.max_normalized_share_sum_error") (`max_norm_error')
foreach v of local sourcevars {
    post `q1' ("audit.missing_by_variable.`v'") (0)
}
foreach v in Char1_BankAssets Char2_BankSupervised Treatment_BankShock {
    post `q1' ("audit.inconsistent_bank_variables.`v'") (0)
}
foreach v in Char4_FirmCreditQuality Outcome_FirmCreditGrowth Outcome_FirmInvestmentGrowth {
    post `q1' ("audit.inconsistent_firm_variables.`v'") (0)
}

foreach v of local sourcevars {
    quietly summarize `v'
    local n = r(N)
    local mean = r(mean)
    local sd = r(sd)
    local min = r(min)
    local max = r(max)
    mata: st_matrix("__q", hw2_type7("`v'"))
    matrix __q = __q
    post `q1' ("summary_loan_level.`v'.mean") (`mean')
    post `q1' ("summary_loan_level.`v'.median") (__q[1,2])
    post `q1' ("summary_loan_level.`v'.min") (`min')
    post `q1' ("summary_loan_level.`v'.max") (`max')
    post `q1' ("summary_loan_level.`v'.q25") (__q[1,1])
    post `q1' ("summary_loan_level.`v'.q75") (__q[1,3])
    post `q1' ("summary_loan_level.`v'.std") (`sd')
    post `q1' ("summary_loan_level.`v'.count") (`n')
}

use `firmsdata', clear
quietly summarize Exposure
local exp_n = r(N)
local exp_mean = r(mean)
local exp_sd = r(sd)
local exp_min = r(min)
local exp_max = r(max)
mata: st_matrix("__q", hw2_type7("Exposure"))
post `q1' ("summary_exposure_firm_level.Exposure.mean") (`exp_mean')
post `q1' ("summary_exposure_firm_level.Exposure.median") (__q[1,2])
post `q1' ("summary_exposure_firm_level.Exposure.min") (`exp_min')
post `q1' ("summary_exposure_firm_level.Exposure.max") (`exp_max')
post `q1' ("summary_exposure_firm_level.Exposure.q25") (__q[1,1])
post `q1' ("summary_exposure_firm_level.Exposure.q75") (__q[1,3])
post `q1' ("summary_exposure_firm_level.Exposure.std") (`exp_sd')
post `q1' ("summary_exposure_firm_level.Exposure.count") (`exp_n')

quietly summarize LoanCount
local bpf_n = r(N)
local bpf_mean = r(mean)
local bpf_sd = r(sd)
local bpf_min = r(min)
local bpf_max = r(max)
mata: st_matrix("__q", hw2_type7("LoanCount"))
post `q1' ("banks_per_firm_summary.count") (`bpf_n')
post `q1' ("banks_per_firm_summary.mean") (`bpf_mean')
post `q1' ("banks_per_firm_summary.std") (`bpf_sd')
post `q1' ("banks_per_firm_summary.min") (`bpf_min')
post `q1' ("banks_per_firm_summary.25%") (__q[1,1])
post `q1' ("banks_per_firm_summary.50%") (__q[1,2])
post `q1' ("banks_per_firm_summary.75%") (__q[1,3])
post `q1' ("banks_per_firm_summary.max") (`bpf_max')

preserve
    contract LoanCount, freq(firms)
    sort LoanCount
    save `bankcounts', replace
    quietly count
    forvalues i=1/`r(N)' {
        local k = LoanCount[`i']
        local f = firms[`i']
        post `q1' ("banks_per_firm_counts.`k'") (`f')
    }
restore

post `q1' ("indicator_shares_loan_rows.Char2_BankSupervised") (`supervised_loan_share')
post `q1' ("indicator_shares_loan_rows.Treatment_BankShock") (`shock_loan_share')
post `q1' ("indicator_shares_unique_banks.Char2_BankSupervised") (`supervised_bank_share')
post `q1' ("indicator_shares_unique_banks.Treatment_BankShock") (`shock_bank_share')
post `q1' ("indicator_counts_unique_banks.Char2_BankSupervised") (`supervised_bank_count')
post `q1' ("indicator_counts_unique_banks.Treatment_BankShock") (`shock_bank_count')

use `banksdata', clear
quietly summarize Char1_BankAssets
local ba_n=r(N)
local ba_mean=r(mean)
local ba_sd=r(sd)
local ba_min=r(min)
local ba_max=r(max)
mata: st_matrix("__q", hw2_type7("Char1_BankAssets"))
post `q1' ("bank_assets_unique_bank_summary.count") (`ba_n')
post `q1' ("bank_assets_unique_bank_summary.mean") (`ba_mean')
post `q1' ("bank_assets_unique_bank_summary.std") (`ba_sd')
post `q1' ("bank_assets_unique_bank_summary.min") (`ba_min')
post `q1' ("bank_assets_unique_bank_summary.25%") (__q[1,1])
post `q1' ("bank_assets_unique_bank_summary.50%") (__q[1,2])
post `q1' ("bank_assets_unique_bank_summary.75%") (__q[1,3])
post `q1' ("bank_assets_unique_bank_summary.max") (`ba_max')
quietly summarize Char1_BankAssets, detail
local ba_skew = sqrt(r(N)*(r(N)-1))/(r(N)-2)*r(skewness)
post `q1' ("bank_assets_skewness") (`ba_skew')

use `loansdata', clear
quietly summarize Outcome_LoanGrowth, detail
local lg_skew = sqrt(r(N)*(r(N)-1))/(r(N)-2)*r(skewness)
post `q1' ("loan_growth_skewness") (`lg_skew')
* Exposure extrema require firm sums, already saved in firmsdata.
use `firmsdata', clear
quietly summarize Exposure_raw, meanonly
post `q1' ("raw_exposure_min") (r(min))
post `q1' ("raw_exposure_max") (r(max))
quietly summarize Exposure, meanonly
post `q1' ("normalized_exposure_min") (r(min))
post `q1' ("normalized_exposure_max") (r(max))
generate double __expdiff = abs(Exposure_raw-Exposure)
quietly summarize __expdiff, meanonly
post `q1' ("max_exposure_normalization_change") (r(max))
postclose `q1'

use `q1results', clear
format value %24.17g
sort metric
export delimited using `"$HW2_ROOT/tmp/results/stata/Question_1_Results_Stata.csv"', replace
preserve
    keep if strpos(metric, "summary_loan_level.") == 1 | ///
        strpos(metric, "summary_exposure_firm_level.") == 1
    export delimited using `"$HW2_ROOT/output/stata/Table_1_1_Summary_Stata.csv"', replace
restore
preserve
    keep if strpos(metric, "audit.") == 1 | strpos(metric, "indicator_") == 1
    export delimited using `"$HW2_ROOT/output/stata/Table_1_2_Sample_Stata.csv"', replace
restore

* Stata-native figures, kept separate from the Python PNGs.
use `loansdata', clear
histogram Outcome_LoanGrowth, density color(navy%80) lcolor(white) ///
    title("A. Loan growth (`loans' loans)") xtitle("Change in log lending") ///
    ytitle("Density") name(q1_lhist, replace)
graph hbox Outcome_LoanGrowth, box(1, color(navy%35)) ///
    ytitle("Change in log lending") name(q1_lbox, replace)
use `banksdata', clear
histogram Char1_BankAssets, density color(navy%80) lcolor(white) xline(30, lpattern(dash)) ///
    title("B. Bank assets (`banks_n' unique banks)") xtitle("Bank assets") ///
    ytitle("Density") name(q1_ahist, replace)
graph hbox Char1_BankAssets, box(1, color(navy%35)) ///
    ytitle("Bank assets") name(q1_abox, replace)
graph combine q1_lhist q1_ahist q1_lbox q1_abox, cols(2) ///
    graphregion(color(white)) name(q1_distributions, replace)
graph export `"$HW2_ROOT/output/stata/Figure_1_Distributions_Stata.png"', ///
    width(3120) replace

use `bankcounts', clear
graph bar (asis) firms, over(LoanCount) blabel(bar, format(%12.0fc)) ///
    bar(1, color(navy)) ytitle("Number of firms") ///
    title("Number of Banks per Firm") graphregion(color(white)) ///
    name(q1_bankcounts, replace)
graph export `"$HW2_ROOT/output/stata/Figure_1_2_Banks_Per_Firm_Stata.png"', ///
    width(2700) replace
graph drop _all

display as result "Question 1 Stata outputs saved."
