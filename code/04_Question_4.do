version 18.0
clear all
set more off
set varabbrev off

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
capture mkdir `"$HW2_ROOT/tmp/results/stata"'

import delimited using `"$HW2_ROOT/data/PhdFinance_FNCE7020_HW2_Data.csv"', ///
    varnames(1) case(preserve) asdouble clear
egen long __rowmiss = rowmiss(FirmID BankID Char3_PreCreditShare ///
    Treatment_BankShock Char4_FirmCreditQuality Outcome_LoanGrowth ///
    Outcome_FirmCreditGrowth Outcome_FirmInvestmentGrowth)
assert __rowmiss == 0
drop __rowmiss
isid FirmID BankID
bysort FirmID: assert Char4_FirmCreditQuality == Char4_FirmCreditQuality[1]
bysort FirmID: assert Outcome_FirmCreditGrowth == Outcome_FirmCreditGrowth[1]
bysort FirmID: assert Outcome_FirmInvestmentGrowth == Outcome_FirmInvestmentGrowth[1]
assert Char3_PreCreditShare >= 0

bysort FirmID: egen double __sharesum = total(Char3_PreCreditShare)
assert __sharesum > 0
generate double __shareerror = abs(__sharesum-1)
quietly summarize __shareerror, meanonly
local max_share_error = r(max)
assert `max_share_error' < .0001
generate double NormalizedShare = Char3_PreCreditShare/__sharesum
generate double WeightedShock = NormalizedShare*Treatment_BankShock
generate double RawWeightedShock = Char3_PreCreditShare*Treatment_BankShock
bysort FirmID: generate long LoanCount = _N

tempfile loans firms q4results effects
save `loans', replace
collapse (firstnm) Char4_FirmCreditQuality Outcome_FirmCreditGrowth ///
    Outcome_FirmInvestmentGrowth LoanCount ///
    (sum) Exposure=WeightedShock RawExposure=RawWeightedShock, by(FirmID)
assert inrange(Exposure,-1e-12,1+1e-12)
save `firms', replace

tempname q4
postfile `q4' str125 metric double value using `q4results', replace

local names credit_unadjusted credit_adjusted investment_unadjusted investment_adjusted
local outcomes Outcome_FirmCreditGrowth Outcome_FirmCreditGrowth ///
    Outcome_FirmInvestmentGrowth Outcome_FirmInvestmentGrowth
local controls 0 1 0 1

forvalues i=1/4 {
    local name : word `i' of `names'
    local outcome : word `i' of `outcomes'
    local control : word `i' of `controls'
    local rhs Exposure
    local terms exposure
    if `control' {
        local rhs `rhs' Char4_FirmCreditQuality
        local terms `terms' quality
    }
    quietly regress `outcome' `rhs', vce(hc3)
    local n = e(N)
    local rdf = e(df_r)
    local r2 = e(r2)
    post `q4' ("models.`name'.n") (`n')
    post `q4' ("models.`name'.residual_df") (`rdf')
    post `q4' ("models.`name'.r_squared") (`r2')

    local j = 0
    foreach term of local terms {
        local ++j
        local v Exposure
        if "`term'" == "quality" local v Char4_FirmCreditQuality
        local b = _b[`v']
        local se = _se[`v']
        local t = `b'/`se'
        local p = 2*ttail(`rdf',abs(`t'))
        local crit = invttail(`rdf',.025)
        post `q4' ("models.`name'.coefficients.`term'.estimate") (`b')
        post `q4' ("models.`name'.coefficients.`term'.std_error") (`se')
        post `q4' ("models.`name'.coefficients.`term'.t_statistic") (`t')
        post `q4' ("models.`name'.coefficients.`term'.p_value") (`p')
        post `q4' ("models.`name'.coefficients.`term'.ci_lower") (`b'-`crit'*`se')
        post `q4' ("models.`name'.coefficients.`term'.ci_upper") (`b'+`crit'*`se')
        if "`term'" == "exposure" & "`name'" == "credit_adjusted" ///
            local credit_effect = `b'
        if "`term'" == "exposure" & "`name'" == "investment_adjusted" ///
            local investment_effect = `b'
    }
    local b = _b[_cons]
    local se = _se[_cons]
    local t = `b'/`se'
    local p = 2*ttail(`rdf',abs(`t'))
    local crit = invttail(`rdf',.025)
    post `q4' ("models.`name'.coefficients.constant.estimate") (`b')
    post `q4' ("models.`name'.coefficients.constant.std_error") (`se')
    post `q4' ("models.`name'.coefficients.constant.t_statistic") (`t')
    post `q4' ("models.`name'.coefficients.constant.p_value") (`p')
    post `q4' ("models.`name'.coefficients.constant.ci_lower") (`b'-`crit'*`se')
    post `q4' ("models.`name'.coefficients.constant.ci_upper") (`b'+`crit'*`se')

    local rawrhs RawExposure
    if `control' local rawrhs `rawrhs' Char4_FirmCreditQuality
    quietly regress `outcome' `rawrhs'
    post `q4' ("models.`name'.raw_exposure_coefficient") (_b[RawExposure])
}

quietly summarize Exposure
local exp_n = r(N)
local exp_mean = r(mean)
local exp_sd = r(sd)
local exp_min = r(min)
local exp_max = r(max)
capture mata: mata drop hw2_q50()
mata:
real scalar hw2_q50(string scalar v)
{
    real colvector x
    real scalar n, h, lo, frac
    x=sort(st_data(.,v),1)
    x=select(x,x:<.)
    n=rows(x)
    h=(n-1)*.5+1
    lo=floor(h)
    frac=h-lo
    if (lo>=n) return(x[n])
    return(x[lo]+frac*(x[lo+1]-x[lo]))
}
end
mata: st_numscalar("__expmedian",hw2_q50("Exposure"))
quietly count if Exposure==0
local zero_exp = r(N)
quietly count if abs(Exposure-1)<=1e-12
local full_exp = r(N)
generate double __expdiff = abs(Exposure-RawExposure)
quietly summarize __expdiff, meanonly
local max_exp_diff = r(max)
quietly correlate Exposure Char4_FirmCreditQuality
matrix __C = r(C)
local exp_quality_corr = __C[1,2]

post `q4' ("exposure_summary.mean") (`exp_mean')
post `q4' ("exposure_summary.median") (scalar(__expmedian))
post `q4' ("exposure_summary.minimum") (`exp_min')
post `q4' ("exposure_summary.maximum") (`exp_max')
post `q4' ("exposure_summary.std_dev") (`exp_sd')
post `q4' ("exposure_summary.zero_exposure_firms") (`zero_exp')
post `q4' ("exposure_summary.full_exposure_firms") (`full_exp')
post `q4' ("exposure_summary.maximum_share_sum_error") (`max_share_error')
post `q4' ("exposure_summary.maximum_raw_exposure_difference") (`max_exp_diff')
post `q4' ("exposure_quality_correlation") (`exp_quality_corr')

* Recompute the multiple-bank within-firm slope independently.
use `loans', clear
keep if LoanCount>=2
bysort FirmID: egen double __xbar = mean(Treatment_BankShock)
bysort FirmID: egen double __ybar = mean(Outcome_LoanGrowth)
generate double __xdm = Treatment_BankShock-__xbar
generate double __ydm = Outcome_LoanGrowth-__ybar
generate double __xy = __xdm*__ydm
generate double __xx = __xdm^2
quietly summarize __xy, meanonly
local sumxy = r(sum)
quietly summarize __xx, meanonly
local sumxx = r(sum)
local loan_effect = `sumxy'/`sumxx'
post `q4' ("loan_firm_fe_effect") (`loan_effect')
post `q4' ("firm_credit_adjusted_effect") (`credit_effect')
post `q4' ("firm_investment_adjusted_effect") (`investment_effect')
post `q4' ("credit_to_loan_magnitude_ratio") (abs(`credit_effect'/`loan_effect'))
post `q4' ("investment_to_credit_magnitude_ratio") (abs(`investment_effect'/`credit_effect'))
postclose `q4'

use `q4results', clear
format value %24.17g
sort metric
export delimited using `"$HW2_ROOT/tmp/results/stata/Question_4_Results_Stata.csv"', replace
preserve
    keep if strpos(metric,"models.")==1
    export delimited using `"$HW2_ROOT/output/stata/Table_4_1_Firm_Regressions_Stata.csv"', replace
restore

* Magnitude comparison figure.
clear
set obs 3
generate byte order = _n
generate double effect = .
replace effect = `loan_effect' in 1
replace effect = `credit_effect' in 2
replace effect = `investment_effect' in 3
generate str38 specification = ""
replace specification = "Loan growth: firm fixed effects" in 1
replace specification = "Total credit growth: quality adjusted" in 2
replace specification = "Investment growth: quality adjusted" in 3
graph hbar (asis) effect, over(specification, sort(order)) ///
    blabel(bar, format(%7.3f) position(inside) color(white)) bar(1,color(navy)) ///
    yline(0,lcolor(gs6)) ytitle("Estimated change in original outcome units") ///
    title("From Bank Lending to Firm Investment") graphregion(color(white)) ///
    name(q4_magnitude, replace)
graph export `"$HW2_ROOT/output/stata/Figure_4_Magnitude_Comparison_Stata.png"', ///
    width(2550) replace
graph drop _all

display as result "Question 4 Stata outputs saved."
