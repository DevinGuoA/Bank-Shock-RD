### Imports and Paths ###

from pathlib import Path
import json

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "output"
RESULTS = ROOT / "tmp" / "results"
OUTPUT.mkdir(parents=True, exist_ok=True)
RESULTS.mkdir(parents=True, exist_ok=True)


### Load Data ###

data = pd.read_csv(ROOT / "data" / "PhdFinance_FNCE7020_HW2_Data.csv")
required = [
    "LoanID", "FirmID", "BankID", "Treatment_BankShock",
    "Char4_FirmCreditQuality", "Outcome_LoanGrowth",
]
assert data[required].notna().all().all()
assert not data.duplicated(["FirmID", "BankID"]).any()
assert data.groupby("BankID")["Treatment_BankShock"].nunique().eq(1).all()
assert data.groupby("FirmID")["Char4_FirmCreditQuality"].nunique().eq(1).all()

# A loan is one firm and bank relationship.
loan_count = data.groupby("FirmID")["LoanID"].transform("size")
multi = data.loc[loan_count.ge(2)].copy()
shock_counts = multi.groupby("FirmID")["Treatment_BankShock"].nunique()
switching_firms = shock_counts.index[shock_counts.eq(2)]


### Loan Regressions ###

model_definitions = [
    ("pooled_all", data, False, False, "Outcome_LoanGrowth"),
    ("pooled_multi", multi, False, False, "Outcome_LoanGrowth"),
    ("pooled_quality", data, True, False, "Outcome_LoanGrowth"),
    ("firm_fe", multi, False, True, "Outcome_LoanGrowth"),
    ("quality_sorting", data, False, False, "Char4_FirmCreditQuality"),
]
models = {}
for name, sample, include_quality, absorb_firm, outcome in model_definitions:
    y = sample[outcome].to_numpy(dtype=float)
    shock = sample["Treatment_BankShock"].to_numpy(dtype=float)
    n = len(sample)
    bank_codes, banks = pd.factorize(sample["BankID"], sort=True)
    firm_codes, firms = pd.factorize(sample["FirmID"], sort=True)

    if absorb_firm:
        # Demeaning removes all firm constants, including credit quality.
        y = y - sample.groupby("FirmID")[outcome].transform("mean").to_numpy()
        shock = shock - sample.groupby("FirmID")["Treatment_BankShock"].transform("mean").to_numpy()
        x = shock[:, None]
        names = ["shock"]
        # Full model rank includes one intercept for every firm.
        full_rank = len(firms) + 1
    else:
        columns = [shock]
        names = ["shock"]
        if include_quality:
            columns.append(sample["Char4_FirmCreditQuality"].to_numpy(dtype=float))
            names.append("quality")
        columns.append(np.ones(n))
        names.append("constant")
        x = np.column_stack(columns)
        full_rank = x.shape[1]

    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    residual = y - np.einsum("ij,j->i", x, beta)
    bread = np.linalg.inv(np.einsum("ij,ik->jk", x, x))
    scores = x * residual[:, None]

    # Two way clustering allows dependence within banks and within firms.
    bank_scores = np.zeros((len(banks), x.shape[1]))
    firm_scores = np.zeros((len(firms), x.shape[1]))
    np.add.at(bank_scores, bank_codes, scores)
    np.add.at(firm_scores, firm_codes, scores)
    finite_sample = (n - 1) / (n - full_rank)
    bank_meat = len(banks) / (len(banks) - 1) * finite_sample * np.einsum("ij,ik->jk", bank_scores, bank_scores)
    firm_meat = len(firms) / (len(firms) - 1) * finite_sample * np.einsum("ij,ik->jk", firm_scores, firm_scores)
    # Each firm and bank intersection contains exactly one loan.
    pair_meat = n / (n - 1) * finite_sample * np.einsum("ij,ik->jk", scores, scores)
    covariance = np.einsum("ij,jk,kl->il", bread, bank_meat + firm_meat - pair_meat, bread)
    bank_covariance = np.einsum("ij,jk,kl->il", bread, bank_meat, bread)
    assert np.all(np.diag(covariance) > 0)
    se = np.sqrt(np.diag(covariance))
    t_stat = beta / se
    inference_df = min(len(banks), len(firms)) - 1
    p_value = 2 * stats.t.sf(np.abs(t_stat), inference_df)
    critical = stats.t.ppf(0.975, inference_df)
    denominator = np.sum(y ** 2) if absorb_firm else np.sum((y - y.mean()) ** 2)

    models[name] = {
        "n": n,
        "firms": len(firms),
        "banks": len(banks),
        "full_rank": full_rank,
        "residual_df": n - full_rank,
        "inference_df": inference_df,
        "r_squared": float(1 - np.sum(residual ** 2) / denominator),
        "r_squared_type": "within" if absorb_firm else "overall",
        "firm_fixed_effects": absorb_firm,
        "coefficients": {
            term: {
                "estimate": float(beta[j]),
                "std_error": float(se[j]),
                "t_statistic": float(t_stat[j]),
                "p_value": float(p_value[j]),
                "ci_lower": float(beta[j] - critical * se[j]),
                "ci_upper": float(beta[j] + critical * se[j]),
                "bank_cluster_std_error": float(np.sqrt(bank_covariance[j, j])),
            }
            for j, term in enumerate(names)
        },
    }


### Sorting and Omitted Variable Check ###

quality_means = data.groupby("Treatment_BankShock")["Char4_FirmCreditQuality"].mean()
quality_counts = data.groupby("Treatment_BankShock").size()
quality_difference = float(quality_means.loc[1] - quality_means.loc[0])
ovb = (
    models["pooled_quality"]["coefficients"]["quality"]["estimate"]
    * quality_difference
)
observed_change = (
    models["pooled_all"]["coefficients"]["shock"]["estimate"]
    - models["pooled_quality"]["coefficients"]["shock"]["estimate"]
)
assert np.isclose(ovb, observed_change, atol=1e-12)
assert np.isclose(
    quality_difference,
    models["quality_sorting"]["coefficients"]["shock"]["estimate"],
    atol=1e-12,
)


### Regression Table ###

reported = ["pooled_all", "pooled_multi", "pooled_quality", "firm_fe"]
lines = [
    r"\begin{table}[H]",
    r"\centering",
    r"\caption{Funding Shocks and Loan Growth}",
    r"\label{tab:q3_loan_regressions}",
    r"\small",
    r"\renewcommand{\arraystretch}{1.10}",
    r"\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}lrrrr@{}}",
    r"\toprule",
    r" & (1) & (2) & (3) & (4) \\",
    r" & Pooled & Pooled & Pooled & Firm FE \\",
    r"Sample & All & Multiple banks & All & Multiple banks \\",
    r"\midrule",
]
for term, label in [("shock", "Bank shock"), ("quality", "Credit quality"), ("constant", "Constant")]:
    estimates = []
    errors = []
    for name in reported:
        coefficient = models[name]["coefficients"].get(term)
        if coefficient is None:
            estimates.append("")
            errors.append("")
            continue
        p = coefficient["p_value"]
        stars = "***" if p < 0.01 else "**" if p < 0.05 else "*" if p < 0.10 else ""
        estimates.append(f"{coefficient['estimate']:.3f}" + (rf"$^{{{stars}}}$" if stars else ""))
        error_text = "$<0.001$" if coefficient["std_error"] < 0.001 else f"{coefficient['std_error']:.3f}"
        errors.append(f"({error_text})")
    lines.append(label + " & " + " & ".join(estimates) + r" \\")
    lines.append(" & " + " & ".join(errors) + r" \\")
lines.append(r"\midrule")
for field, label in [("n", "Loans"), ("firms", "Firms"), ("banks", "Banks")]:
    lines.append(label + " & " + " & ".join(f"{models[name][field]:,}" for name in reported) + r" \\")
lines.extend([
    r"Firm fixed effects & No & No & No & Yes \\",
    r"$R^2$ & " + " & ".join(f"{models[name]['r_squared']:.3f}" for name in reported) + r" \\",
    r"\bottomrule",
    r"\end{tabular*}",
    r"\par\smallskip",
    r"\begin{minipage}{\linewidth}\footnotesize",
    r"\textit{Notes:} Specifications are $Y_{fb}=\alpha+\beta Shock_b+\varepsilon_{fb}$ in columns (1) and (2), with $\gamma Quality_f$ added in (3), and $Y_{fb}=\alpha_f+\beta Shock_b+\varepsilon_{fb}$ in (4). The dependent variable is the change in log lending. Multiple banks means at least two banks per firm. Credit quality is absorbed by firm effects. Parentheses report standard errors clustered by bank and firm. Each cluster component uses its cluster count correction and $(N-1)/(N-K)$, including absorbed firm effects in $K$. Tests use 319 degrees of freedom. The last $R^2$ is within firm. Positive standard errors below 0.001 are shown as $<0.001$; saved numerical results retain full precision. $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$.",
    r"\end{minipage}",
    r"\end{table}",
])
(OUTPUT / "Table_3_1_Loan_Regressions.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")


### Credit Quality Table ###

sorting = models["quality_sorting"]["coefficients"]["shock"]
p_text = "$<0.001$" if sorting["p_value"] < 0.001 else f"{sorting['p_value']:.3f}"
lines = [
    r"\begin{table}[H]",
    r"\centering",
    r"\caption{Borrower Credit Quality by Bank Shock Status}",
    r"\label{tab:q3_quality_sorting}",
    r"\small",
    r"\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}lrrrrr@{}}",
    r"\toprule",
    r" & No shock & Shock & Difference & Std. error & $p$ value \\",
    r"\midrule",
    "Credit quality & " + f"{quality_means.loc[0]:.3f} & {quality_means.loc[1]:.3f} & {quality_difference:.3f} & {sorting['std_error']:.3f} & {p_text}" + r" \\",
    "Loans & " + f"{quality_counts.loc[0]:,} & {quality_counts.loc[1]:,} & & &" + r" \\",
    r"\bottomrule",
    r"\end{tabular*}",
    r"\par\smallskip",
    r"\begin{minipage}{\linewidth}\footnotesize",
    r"\textit{Notes:} Means weight each observed loan equally, matching the pooled regression sample. Difference is Shock minus No shock. The standard error is clustered by bank and firm in a regression of credit quality on the shock indicator and a constant. The two sided test uses 319 degrees of freedom. A firm borrowing from both bank types can enter both means.",
    r"\end{minipage}",
    r"\end{table}",
]
(OUTPUT / "Table_3_2_Credit_Quality_Sorting.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")


### Save Numerical Results ###

results = {
    "models": models,
    "multi_bank_firms": int(len(shock_counts)),
    "firms_with_within_shock_variation": int(len(switching_firms)),
    "firms_without_within_shock_variation": int(shock_counts.eq(1).sum()),
    "loans_in_firms_with_within_shock_variation": int(multi["FirmID"].isin(switching_firms).sum()),
    "quality_no_shock": float(quality_means.loc[0]),
    "quality_shock": float(quality_means.loc[1]),
    "quality_difference": quality_difference,
    "omitted_quality_bias": ovb,
    "pooled_coefficient_change_with_quality": observed_change,
}
(RESULTS / "Question_3_Results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
print("Question 3: loan regressions and credit quality comparison saved.")
