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
capture mkdir `"$HW2_ROOT/tmp/stata_ado_hw2_20260527"'

* Use project-local, current RD packages.  The globally installed rdrobust on
* some machines is a 2017 release and does not implement masspoints(adjust).
sysdir set PLUS `"$HW2_ROOT/tmp/stata_ado_hw2_20260527"'
if !fileexists(`"$HW2_ROOT/tmp/stata_ado_hw2_20260527/r/rdrobust.ado"') {
    noisily net install rdrobust, ///
        from(https://raw.githubusercontent.com/rdpackages/rdrobust/d3823d73ff3098d5e87736be9bf0d90f35d6eafa/stata) replace
}
if !fileexists(`"$HW2_ROOT/tmp/stata_ado_hw2_20260527/r/rddensity.ado"') {
    noisily net install rddensity, ///
        from(https://raw.githubusercontent.com/rdpackages/rddensity/790a2fd973351c67b5204c8ac4b621d71a2bd75a/stata) replace
}
discard
capture which rdrobust
if _rc {
    display as error "rdrobust is required; see the project-local installation block above."
    exit 111
}
capture which rddensity
if _rc {
    display as error "rddensity is required; see the project-local installation block above."
    exit 111
}

import delimited using `"$HW2_ROOT/data/PhdFinance_FNCE7020_HW2_Data.csv"', ///
    varnames(1) case(preserve) asdouble clear
egen long __rowmiss = rowmiss(_all)
assert __rowmiss == 0
drop __rowmiss
isid FirmID BankID
bysort BankID: assert Char1_BankAssets == Char1_BankAssets[1]
bysort BankID: assert Char2_BankSupervised == Char2_BankSupervised[1]
bysort BankID: assert Treatment_BankShock == Treatment_BankShock[1]
quietly count
local loan_n = r(N)

collapse (firstnm) assets=Char1_BankAssets supervised=Char2_BankSupervised ///
    shock=Treatment_BankShock ///
    (mean) loan_growth=Outcome_LoanGrowth borrower_quality=Char4_FirmCreditQuality ///
    (count) loan_count=LoanID, by(BankID)
assert supervised == (assets >= 30)
generate double centered_assets = assets - 30
generate double sx = supervised * centered_assets
generate double shock_percent = 100 * shock
quietly count
local bank_n = r(N)
tempfile banks q2results balancebins densitybins
save `banks', replace

tempname q2
postfile `q2' str110 metric double value using `q2results', replace
post `q2' ("loan_observations") (`loan_n')
post `q2' ("bank_observations") (`bank_n')
quietly count if assets < 30
post `q2' ("bank_left") (r(N))
quietly count if assets >= 30
post `q2' ("bank_right") (r(N))
quietly count if assets == 30
post `q2' ("bank_assets_at_cutoff") (r(N))
quietly levelsof assets, local(__assetlevels)
local asset_unique : word count `__assetlevels'
post `q2' ("bank_assets_unique") (`asset_unique')
quietly summarize assets if assets < 30, meanonly
post `q2' ("nearest_assets_left") (r(max))
quietly summarize assets if assets >= 30, meanonly
post `q2' ("nearest_assets_right") (r(min))

local wi = 0
foreach bw in 3 5 8 {
    quietly regress loan_growth supervised centered_assets sx ///
        if abs(centered_assets) <= `bw', vce(hc3)
    local b = _b[supervised]
    local se = _se[supervised]
    local df = e(df_r)
    local t = `b'/`se'
    local p = 2*ttail(`df', abs(`t'))
    local crit = invttail(`df', .025)
    quietly count if abs(centered_assets)<=`bw' & supervised==0
    local nl = r(N)
    quietly count if abs(centered_assets)<=`bw' & supervised==1
    local nr = r(N)
    post `q2' ("windows[`wi'].bandwidth") (`bw')
    post `q2' ("windows[`wi'].estimate") (`b')
    post `q2' ("windows[`wi'].standard_error") (`se')
    post `q2' ("windows[`wi'].t_statistic") (`t')
    post `q2' ("windows[`wi'].p_value") (`p')
    post `q2' ("windows[`wi'].ci_lower") (`b'-`crit'*`se')
    post `q2' ("windows[`wi'].ci_upper") (`b'+`crit'*`se')
    post `q2' ("windows[`wi'].n_left") (`nl')
    post `q2' ("windows[`wi'].n_right") (`nr')
    post `q2' ("windows[`wi'].degrees_of_freedom") (`df')
    post `q2' ("windows[`wi'].coefficients[0]") (_b[_cons])
    post `q2' ("windows[`wi'].coefficients[1]") (_b[supervised])
    post `q2' ("windows[`wi'].coefficients[2]") (_b[centered_assets])
    post `q2' ("windows[`wi'].coefficients[3]") (_b[sx])
    local ++wi
}

* Automatic local-linear RD, matching Python rdrobust 2.0.0 options.
quietly rdrobust loan_growth assets, c(30) p(1) q(2) kernel(triangular) ///
    bwselect(mserd) vce(nn 3) masspoints(adjust) level(95)
post `q2' ("rdrobust.conventional_estimate") (e(tau_cl))
post `q2' ("rdrobust.conventional_standard_error") (e(se_tau_cl))
post `q2' ("rdrobust.conventional_p_value") (e(pv_cl))
post `q2' ("rdrobust.conventional_ci[0]") (e(ci_l_cl))
post `q2' ("rdrobust.conventional_ci[1]") (e(ci_r_cl))
post `q2' ("rdrobust.bias_corrected_estimate") (e(tau_bc))
post `q2' ("rdrobust.robust_standard_error") (e(se_tau_rb))
post `q2' ("rdrobust.robust_z_statistic") (e(tau_bc)/e(se_tau_rb))
post `q2' ("rdrobust.robust_p_value") (e(pv_rb))
post `q2' ("rdrobust.robust_ci[0]") (e(ci_l_rb))
post `q2' ("rdrobust.robust_ci[1]") (e(ci_r_rb))
post `q2' ("rdrobust.h_left") (e(h_l))
post `q2' ("rdrobust.h_right") (e(h_r))
post `q2' ("rdrobust.b_left") (e(b_l))
post `q2' ("rdrobust.b_right") (e(b_r))
post `q2' ("rdrobust.n_h_left") (e(N_h_l))
post `q2' ("rdrobust.n_h_right") (e(N_h_r))
post `q2' ("rdrobust.n_b_left") (e(N_b_l))
post `q2' ("rdrobust.n_b_right") (e(N_b_r))

* Limited balance diagnostics within [25,35].
local bi = 0
foreach v in shock_percent borrower_quality loan_count {
    quietly regress `v' supervised centered_assets sx ///
        if abs(centered_assets) <= 5, vce(hc3)
    local b = _b[supervised]
    local se = _se[supervised]
    local df = e(df_r)
    local t = `b'/`se'
    local p = 2*ttail(`df', abs(`t'))
    local crit = invttail(`df', .025)
    quietly count if abs(centered_assets)<=5 & supervised==0
    local nl = r(N)
    quietly count if abs(centered_assets)<=5 & supervised==1
    local nr = r(N)
    post `q2' ("balance_diagnostics[`bi'].estimate") (`b')
    post `q2' ("balance_diagnostics[`bi'].standard_error") (`se')
    post `q2' ("balance_diagnostics[`bi'].t_statistic") (`t')
    post `q2' ("balance_diagnostics[`bi'].p_value") (`p')
    post `q2' ("balance_diagnostics[`bi'].ci_lower") (`b'-`crit'*`se')
    post `q2' ("balance_diagnostics[`bi'].ci_upper") (`b'+`crit'*`se')
    post `q2' ("balance_diagnostics[`bi'].n_left") (`nl')
    post `q2' ("balance_diagnostics[`bi'].n_right") (`nr')
    local ++bi
}

* rddensity 3.0 double precision.  Its port differs slightly from Python
* rddensity 2.4.6 in bandwidth optimization; validation documents a narrow,
* density-only tolerance while retaining the fully automatic specification.
quietly rddensity assets, c(30) p(2) q(3) fitselect(unrestricted) ///
    kernel(triangular) vce(jackknife) bwselect(comb) precision(double)
post `q2' ("density_test.estimate_left") (e(f_ql))
post `q2' ("density_test.estimate_right") (e(f_qr))
post `q2' ("density_test.estimate_difference") (e(f_qr)-e(f_ql))
post `q2' ("density_test.standard_error") (e(se_q))
post `q2' ("density_test.z_statistic") (e(T_q))
post `q2' ("density_test.p_value") (e(pv_q))
post `q2' ("density_test.h_left") (e(h_l))
post `q2' ("density_test.h_right") (e(h_r))
post `q2' ("density_test.n_left") (e(N_h_l))
post `q2' ("density_test.n_right") (e(N_h_r))
post `q2' ("density_test.p") (e(p))
post `q2' ("density_test.q") (e(q))
postclose `q2'

use `q2results', clear
format value %24.17g
sort metric
export delimited using `"$HW2_ROOT/tmp/results/stata/Question_2_Results_Stata.csv"', replace
preserve
    keep if strpos(metric, "windows[") == 1
    export delimited using `"$HW2_ROOT/output/stata/Table_2_1_RD_Windows_Stata.csv"', replace
restore
preserve
    keep if strpos(metric, "rdrobust.") == 1
    export delimited using `"$HW2_ROOT/output/stata/Table_2_2_RD_Robust_Stata.csv"', replace
restore
preserve
    keep if strpos(metric, "balance_diagnostics[") == 1 | ///
        strpos(metric, "density_test.") == 1
    export delimited using `"$HW2_ROOT/output/stata/Table_2_3_Validity_Checks_Stata.csv"', replace
restore

* RD fits at the three prespecified windows.
use `banks', clear
local gi = 0
foreach bw in 3 5 8 {
    local ++gi
    twoway ///
        (scatter loan_growth assets if abs(centered_assets)<=`bw' & supervised==0, ///
            mcolor(navy%60) msize(small)) ///
        (scatter loan_growth assets if abs(centered_assets)<=`bw' & supervised==1, ///
            mcolor(orange_red%60) msize(small)) ///
        (lfit loan_growth assets if abs(centered_assets)<=`bw' & supervised==0, ///
            lcolor(navy) lwidth(medthick)) ///
        (lfit loan_growth assets if abs(centered_assets)<=`bw' & supervised==1, ///
            lcolor(orange_red) lwidth(medthick)), ///
        xline(30, lpattern(dot) lcolor(gs8)) legend(off) ///
        title("Asset window: 30 +/- `bw'") xtitle("Bank assets") ///
        ytitle("Mean change in log lending") graphregion(color(white)) ///
        name(q2_rd`gi', replace)
}
graph combine q2_rd1 q2_rd2 q2_rd3, cols(3) ycommon ///
    graphregion(color(white)) name(q2_rdall, replace)
graph export `"$HW2_ROOT/output/stata/Figure_2_1_RD_Fits_Stata.png"', ///
    width(3180) replace

* Equal-width density bars use all 320 banks in the denominator, as in Python.
use `banks', clear
generate int __bin = floor((assets-10)/2)
keep if inrange(__bin,0,19)
collapse (count) count=BankID, by(__bin)
generate double center = 11 + 2*__bin
generate double density = count/(`bank_n'*2)
save `densitybins', replace
twoway (bar density center if center<30, barwidth(1.85) color(navy%65)) ///
    (bar density center if center>=30, barwidth(1.85) color(orange_red%65)), ///
    xline(30, lpattern(dot) lcolor(gs6)) xscale(range(10 50)) ///
    xlabel(10(5)50) xtitle("Bank assets") ytitle("Density") legend(off) ///
    graphregion(color(white)) name(q2_density, replace)
graph export `"$HW2_ROOT/output/stata/Figure_2_2_Asset_Density_Stata.png"', ///
    width(2340) replace

* Binned descriptive points and full-sample local-linear fits for proxies.
use `banks', clear
keep if abs(centered_assets)<=5
generate byte asset_bin = min(floor(assets-25),9)
foreach v in shock_percent borrower_quality loan_count {
    quietly regress `v' supervised centered_assets sx, vce(hc3)
    generate double fit_`v' = _b[_cons] + _b[supervised]*supervised + ///
        _b[centered_assets]*centered_assets + _b[sx]*sx
}
collapse (mean) assets shock_percent borrower_quality loan_count ///
    fit_shock_percent fit_borrower_quality fit_loan_count, ///
    by(supervised asset_bin)
save `balancebins', replace
local gi = 0
foreach v in shock_percent borrower_quality loan_count {
    local ++gi
    local ttl "Bank funding shock"
    local yttl "Shocked banks (%)"
    if "`v'" == "borrower_quality" {
        local ttl "Mean borrower credit quality"
        local yttl "Credit quality (index)"
    }
    if "`v'" == "loan_count" {
        local ttl "Number of loan relationships"
        local yttl "Loan relationships"
    }
    twoway ///
        (scatter `v' assets if supervised==0, mcolor(navy) msize(medsmall)) ///
        (scatter `v' assets if supervised==1, mcolor(orange_red) msize(medsmall)) ///
        (line fit_`v' assets if supervised==0, sort lcolor(navy) lwidth(medthick)) ///
        (line fit_`v' assets if supervised==1, sort lcolor(orange_red) lwidth(medthick)), ///
        xline(30, lpattern(dot) lcolor(gs8)) legend(off) title("`ttl'") ///
        xtitle("Bank assets") ytitle("`yttl'") graphregion(color(white)) ///
        name(q2_bal`gi', replace)
}
graph combine q2_bal1 q2_bal2 q2_bal3, cols(3) ///
    graphregion(color(white)) name(q2_balall, replace)
graph export `"$HW2_ROOT/output/stata/Figure_2_3_Bank_Characteristics_Stata.png"', ///
    width(3180) replace
graph drop _all

display as result "Question 2 Stata outputs saved."
