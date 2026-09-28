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
egen long __rowmiss = rowmiss(LoanID FirmID BankID Treatment_BankShock ///
    Char4_FirmCreditQuality Outcome_LoanGrowth)
assert __rowmiss == 0
drop __rowmiss
isid FirmID BankID
assert inlist(Treatment_BankShock,0,1)
bysort FirmID: generate long LoanCount = _N
generate double one = 1
bysort FirmID: egen double __shock_min = min(Treatment_BankShock)
bysort FirmID: egen double __shock_max = max(Treatment_BankShock)
egen byte __firmtag = tag(FirmID)
quietly count if __firmtag & LoanCount>=2
local multi_firms = r(N)
quietly count if __firmtag & LoanCount>=2 & __shock_min<__shock_max
local switching_firms = r(N)
quietly count if __firmtag & LoanCount>=2 & __shock_min==__shock_max
local nonswitching_firms = r(N)
quietly count if LoanCount>=2 & __shock_min<__shock_max
local switching_loans = r(N)
tempfile alldata q3results
save `alldata', replace

* Exact covariance used in the Python script:
* V_bank + V_firm - V_(bank,firm), with component-specific cluster factors
* and (N-1)/(N-K).  K includes absorbed firm effects in the FE model.
capture mata: mata drop hw2_twoway_fit()
mata:
void hw2_twoway_fit(string scalar yvar, string scalar xvars,
    string scalar bankvar, string scalar firmvar, real scalar fullrank)
{
    real colvector y, u, bank, firm, ob, of
    real matrix X, bread, scores, sb, sf, mb, mf, mp, V, Vb
    real matrix ib, iff
    real colvector b, yc
    real scalar n, gb, gf, finite, rss, tssc, tssu

    y = st_data(., yvar)
    X = st_data(., tokens(xvars))
    bank = st_data(., bankvar)
    firm = st_data(., firmvar)
    n = rows(y)
    b = qrsolve(X, y)
    u = y-X*b
    bread = invsym(quadcross(X,X))
    scores = X:*u

    ob = order(bank,1)
    ib = panelsetup(bank[ob],1)
    sb = panelsum(scores[ob,.],ib)
    gb = rows(ib)
    of = order(firm,1)
    iff = panelsetup(firm[of],1)
    sf = panelsum(scores[of,.],iff)
    gf = rows(iff)

    finite = (n-1)/(n-fullrank)
    mb = gb/(gb-1)*finite*quadcross(sb,sb)
    mf = gf/(gf-1)*finite*quadcross(sf,sf)
    mp = n/(n-1)*finite*quadcross(scores,scores)
    V = bread*(mb+mf-mp)*bread
    Vb = bread*mb*bread

    rss = quadcross(u,u)
    yc = y :- mean(y)
    tssc = quadcross(yc,yc)
    tssu = quadcross(y,y)
    st_matrix("__hw2_b", b')
    st_matrix("__hw2_V", V)
    st_matrix("__hw2_Vbank", Vb)
    st_numscalar("__hw2_n",n)
    st_numscalar("__hw2_gb",gb)
    st_numscalar("__hw2_gf",gf)
    st_numscalar("__hw2_rss",rss)
    st_numscalar("__hw2_tssc",tssc)
    st_numscalar("__hw2_tssu",tssu)
}
end

capture program drop hw2_post_twoway
program define hw2_post_twoway, rclass
    version 18.0
    syntax, MODEL(string) OUTCOME(name) XVARS(string) TERMS(string) ///
        FULLRank(real) POSTName(name) [WITHIN]
    mata: hw2_twoway_fit("`outcome'", "`xvars'", "BankID", "FirmID", `fullrank')
    matrix __b = __hw2_b
    matrix __V = __hw2_V
    matrix __Vb = __hw2_Vbank
    local n = scalar(__hw2_n)
    local gb = scalar(__hw2_gb)
    local gf = scalar(__hw2_gf)
    local df = min(`gb',`gf')-1
    local rdf = `n'-`fullrank'
    local denom = scalar(__hw2_tssc)
    if `"`within'"' != "" local denom = scalar(__hw2_tssu)
    local r2 = 1-scalar(__hw2_rss)/`denom'
    post `postname' ("models.`model'.n") (`n')
    post `postname' ("models.`model'.firms") (`gf')
    post `postname' ("models.`model'.banks") (`gb')
    post `postname' ("models.`model'.full_rank") (`fullrank')
    post `postname' ("models.`model'.residual_df") (`rdf')
    post `postname' ("models.`model'.inference_df") (`df')
    post `postname' ("models.`model'.r_squared") (`r2')
    local j = 0
    foreach term of local terms {
        local ++j
        local b = __b[1,`j']
        local se = sqrt(__V[`j',`j'])
        local t = `b'/`se'
        local p = 2*ttail(`df',abs(`t'))
        local crit = invttail(`df',.025)
        local seb = sqrt(__Vb[`j',`j'])
        post `postname' ("models.`model'.coefficients.`term'.estimate") (`b')
        post `postname' ("models.`model'.coefficients.`term'.std_error") (`se')
        post `postname' ("models.`model'.coefficients.`term'.t_statistic") (`t')
        post `postname' ("models.`model'.coefficients.`term'.p_value") (`p')
        post `postname' ("models.`model'.coefficients.`term'.ci_lower") (`b'-`crit'*`se')
        post `postname' ("models.`model'.coefficients.`term'.ci_upper") (`b'+`crit'*`se')
        post `postname' ("models.`model'.coefficients.`term'.bank_cluster_std_error") (`seb')
    }
    return matrix b = __b
    return matrix V = __V
    return scalar r2 = `r2'
end

tempname q3
postfile `q3' str120 metric double value using `q3results', replace

use `alldata', clear
hw2_post_twoway, model(pooled_all) outcome(Outcome_LoanGrowth) ///
    xvars("Treatment_BankShock one") terms("shock constant") ///
    fullrank(2) postname(`q3')
matrix __b_pooled = r(b)

use `alldata', clear
keep if LoanCount>=2
hw2_post_twoway, model(pooled_multi) outcome(Outcome_LoanGrowth) ///
    xvars("Treatment_BankShock one") terms("shock constant") ///
    fullrank(2) postname(`q3')

use `alldata', clear
hw2_post_twoway, model(pooled_quality) outcome(Outcome_LoanGrowth) ///
    xvars("Treatment_BankShock Char4_FirmCreditQuality one") ///
    terms("shock quality constant") fullrank(3) postname(`q3')
matrix __b_quality = r(b)

use `alldata', clear
keep if LoanCount>=2
bysort FirmID: egen double __ybar = mean(Outcome_LoanGrowth)
bysort FirmID: egen double __xbar = mean(Treatment_BankShock)
generate double __ydm = Outcome_LoanGrowth-__ybar
generate double __xdm = Treatment_BankShock-__xbar
local fe_rank = `multi_firms'+1
hw2_post_twoway, model(firm_fe) outcome(__ydm) xvars("__xdm") ///
    terms("shock") fullrank(`fe_rank') postname(`q3') within
matrix __b_fe = r(b)

use `alldata', clear
hw2_post_twoway, model(quality_sorting) outcome(Char4_FirmCreditQuality) ///
    xvars("Treatment_BankShock one") terms("shock constant") ///
    fullrank(2) postname(`q3')
matrix __b_sort = r(b)

quietly summarize Char4_FirmCreditQuality if Treatment_BankShock==0, meanonly
local quality0 = r(mean)
local n0 = r(N)
quietly summarize Char4_FirmCreditQuality if Treatment_BankShock==1, meanonly
local quality1 = r(mean)
local n1 = r(N)
local qdiff = `quality1'-`quality0'
local ovb = __b_quality[1,2]*`qdiff'
local observed = __b_pooled[1,1]-__b_quality[1,1]
assert reldif(`ovb',`observed') < 1e-12
assert reldif(`qdiff',__b_sort[1,1]) < 1e-12

post `q3' ("multi_bank_firms") (`multi_firms')
post `q3' ("firms_with_within_shock_variation") (`switching_firms')
post `q3' ("firms_without_within_shock_variation") (`nonswitching_firms')
post `q3' ("loans_in_firms_with_within_shock_variation") (`switching_loans')
post `q3' ("quality_no_shock") (`quality0')
post `q3' ("quality_shock") (`quality1')
post `q3' ("quality_difference") (`qdiff')
post `q3' ("omitted_quality_bias") (`ovb')
post `q3' ("pooled_coefficient_change_with_quality") (`observed')
postclose `q3'

use `q3results', clear
format value %24.17g
sort metric
export delimited using `"$HW2_ROOT/tmp/results/stata/Question_3_Results_Stata.csv"', replace
preserve
    keep if strpos(metric,"models.pooled_all.")==1 | ///
        strpos(metric,"models.pooled_multi.")==1 | ///
        strpos(metric,"models.pooled_quality.")==1 | ///
        strpos(metric,"models.firm_fe.")==1
    export delimited using `"$HW2_ROOT/output/stata/Table_3_1_Loan_Regressions_Stata.csv"', replace
restore
preserve
    keep if strpos(metric,"models.quality_sorting.")==1 | ///
        inlist(metric,"quality_no_shock","quality_shock","quality_difference")
    export delimited using `"$HW2_ROOT/output/stata/Table_3_2_Credit_Quality_Sorting_Stata.csv"', replace
restore

display as result "Question 3 Stata outputs saved."
