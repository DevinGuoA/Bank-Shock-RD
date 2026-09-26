### Imports and Paths ###

from pathlib import Path
import json
import os

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "tmp" / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(ROOT / "tmp" / "cache"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


OUTPUT = ROOT / "output"
RESULTS = ROOT / "tmp" / "results"
OUTPUT.mkdir(parents=True, exist_ok=True)
RESULTS.mkdir(parents=True, exist_ok=True)


### Construct Firm Exposure ###

data = pd.read_csv(ROOT / "data" / "PhdFinance_FNCE7020_HW2_Data.csv")
firm_variables = [
    "Char4_FirmCreditQuality",
    "Outcome_FirmCreditGrowth",
    "Outcome_FirmInvestmentGrowth",
]
required = firm_variables + [
    "FirmID", "BankID", "Char3_PreCreditShare",
    "Treatment_BankShock", "Outcome_LoanGrowth",
]
assert data[required].notna().all().all()
assert data.groupby("FirmID")[firm_variables].nunique().eq(1).all().all()
assert not data.duplicated(["FirmID", "BankID"]).any()
assert data["Char3_PreCreditShare"].ge(0).all()

# Normalize only the small discrepancy from five decimal source rounding.
share_sum = data.groupby("FirmID")["Char3_PreCreditShare"].transform("sum")
assert share_sum.gt(0).all()
assert np.max(np.abs(share_sum - 1)) < 0.0001
data["NormalizedShare"] = data["Char3_PreCreditShare"] / share_sum
data["WeightedShock"] = data["NormalizedShare"] * data["Treatment_BankShock"]
data["RawWeightedShock"] = data["Char3_PreCreditShare"] * data["Treatment_BankShock"]

# Each firm has one observation, irrespective of its number of lenders.
firms = data.groupby("FirmID")[firm_variables].first()
firms["Exposure"] = data.groupby("FirmID")["WeightedShock"].sum()
firms["RawExposure"] = data.groupby("FirmID")["RawWeightedShock"].sum()
assert firms["Exposure"].between(-1e-12, 1 + 1e-12).all()


### Firm Regressions ###

model_definitions = [
    ("credit_unadjusted", "Outcome_FirmCreditGrowth", False),
    ("credit_adjusted", "Outcome_FirmCreditGrowth", True),
    ("investment_unadjusted", "Outcome_FirmInvestmentGrowth", False),
    ("investment_adjusted", "Outcome_FirmInvestmentGrowth", True),
]
models = {}
for name, outcome, include_quality in model_definitions:
    y = firms[outcome].to_numpy(dtype=float)
    columns = [firms["Exposure"].to_numpy(dtype=float)]
    names = ["exposure"]
    if include_quality:
        columns.append(firms["Char4_FirmCreditQuality"].to_numpy(dtype=float))
        names.append("quality")
    columns.append(np.ones(len(firms)))
    names.append("constant")
    x = np.column_stack(columns)
    n, k = x.shape
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    residual = y - np.einsum("ij,j->i", x, beta)
    bread = np.linalg.inv(np.einsum("ij,ik->jk", x, x))

    # HC3 adjusts each squared residual for its observation leverage.
    leverage = np.einsum("ij,jk,ik->i", x, bread, x)
    adjusted_scores = x * (residual / (1 - leverage))[:, None]
    meat = np.einsum("ij,ik->jk", adjusted_scores, adjusted_scores)
    covariance = np.einsum("ij,jk,kl->il", bread, meat, bread)
    se = np.sqrt(np.diag(covariance))
    t_stat = beta / se
    p_value = 2 * stats.t.sf(np.abs(t_stat), n - k)
    critical = stats.t.ppf(0.975, n - k)

    # Raw shares provide a sensitivity check without changing the sample.
    x_raw = x.copy()
    x_raw[:, 0] = firms["RawExposure"].to_numpy(dtype=float)
    beta_raw = np.linalg.lstsq(x_raw, y, rcond=None)[0]
    models[name] = {
        "outcome": outcome,
        "n": n,
        "residual_df": n - k,
        "r_squared": float(1 - np.sum(residual ** 2) / np.sum((y - y.mean()) ** 2)),
        "credit_quality_control": include_quality,
        "raw_exposure_coefficient": float(beta_raw[0]),
        "coefficients": {
            term: {
                "estimate": float(beta[j]),
                "std_error": float(se[j]),
                "t_statistic": float(t_stat[j]),
                "p_value": float(p_value[j]),
                "ci_lower": float(beta[j] - critical * se[j]),
                "ci_upper": float(beta[j] + critical * se[j]),
            }
            for j, term in enumerate(names)
        },
    }


### Regression Table ###

reported = [name for name, outcome, include_quality in model_definitions]
lines = [
    r"\begin{table}[H]",
    r"\centering",
    r"\caption{Shock Exposure, Total Credit, and Investment}",
    r"\label{tab:q4_firm_regressions}",
    r"\small",
    r"\renewcommand{\arraystretch}{1.10}",
    r"\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}lrrrr@{}}",
    r"\toprule",
    r" & \multicolumn{2}{c}{Total credit growth} & \multicolumn{2}{c}{Investment growth} \\",
    r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}",
    r" & (1) & (2) & (3) & (4) \\",
    r"\midrule",
]
for term, label in [("exposure", "Shock exposure"), ("quality", "Credit quality"), ("constant", "Constant")]:
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
        estimate_text = f"{coefficient['estimate']:.3f}"
        if estimate_text == "-0.000":
            estimate_text = "0.000"
        estimates.append(estimate_text + (rf"$^{{{stars}}}$" if stars else ""))
        error_text = "$<0.001$" if coefficient["std_error"] < 0.001 else f"{coefficient['std_error']:.3f}"
        errors.append(f"({error_text})")
    lines.append(label + " & " + " & ".join(estimates) + r" \\")
    lines.append(" & " + " & ".join(errors) + r" \\")
lines.extend([
    r"\midrule",
    "Firms & " + " & ".join(f"{models[name]['n']:,}" for name in reported) + r" \\",
    r"Credit quality control & No & Yes & No & Yes \\",
    r"$R^2$ & " + " & ".join(f"{models[name]['r_squared']:.3f}" for name in reported) + r" \\",
    r"\bottomrule",
    r"\end{tabular*}",
    r"\par\smallskip",
    r"\begin{minipage}{\linewidth}\footnotesize",
    r"\textit{Notes:} Specifications are $Y_f=\alpha+\beta Exposure_f+u_f$ in columns (1) and (3), with $\gamma Quality_f$ added in (2) and (4). Outcomes are total credit growth in (1) and (2) and investment growth in (3) and (4). Exposure is the preperiod borrowing share supplied by shocked banks, with weights normalized within firm to sum to one. Firms are equally weighted. Parentheses report HC3 heteroskedasticity robust standard errors; $t$ tests use $N-K$ degrees of freedom. Inference assumes independence across firms. Positive standard errors below 0.001 are displayed as $<0.001$. $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$.",
    r"\end{minipage}",
    r"\end{table}",
])
(OUTPUT / "Table_4_1_Firm_Regressions.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")


### Compare Loan and Firm Magnitudes ###

# Recompute the within firm slope so this script also runs independently.
loan_count = data.groupby("FirmID")["BankID"].transform("size")
multi = data.loc[loan_count.ge(2)]
within_shock = (
    multi["Treatment_BankShock"]
    - multi.groupby("FirmID")["Treatment_BankShock"].transform("mean")
).to_numpy(dtype=float)
within_growth = (
    multi["Outcome_LoanGrowth"]
    - multi.groupby("FirmID")["Outcome_LoanGrowth"].transform("mean")
).to_numpy(dtype=float)
loan_effect = float(np.sum(within_shock * within_growth) / np.sum(within_shock ** 2))
credit_effect = models["credit_adjusted"]["coefficients"]["exposure"]["estimate"]
investment_effect = models["investment_adjusted"]["coefficients"]["exposure"]["estimate"]

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.spines.left": False,
    "axes.titlepad": 14,
    "savefig.facecolor": "white",
})
effects = [loan_effect, credit_effect, investment_effect]
labels = ["Loan growth\nFirm fixed effects", "Total credit growth\nQuality adjusted", "Investment growth\nQuality adjusted"]
fig, ax = plt.subplots(figsize=(8.5, 3.6))
ax.barh(range(3), effects, height=0.54, color=["#203B55", "#43788E", "#86ADB5"])
ax.set_yticks(range(3), labels)
ax.invert_yaxis()
ax.set_xlim(min(effects) * 1.18, 0.007)
ax.axvline(0, color="#404040", linewidth=0.8)
ax.tick_params(axis="y", length=0, pad=12)
ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.3f"))
ax.set_xlabel("Estimated change in the original outcome units")
ax.set_title("From Bank Lending to Firm Investment", loc="left", fontweight="bold")
ax.grid(axis="x", alpha=0.15)
ax.set_axisbelow(True)
for position, effect in enumerate(effects):
    ax.text(effect - 0.004, position, f"{effect:.3f}", ha="right", va="center", fontsize=11)
fig.tight_layout()
fig.savefig(OUTPUT / "Figure_4_Magnitude_Comparison.png", dpi=300, bbox_inches="tight")
plt.close(fig)


### Save Numerical Results ###

results = {
    "models": models,
    "exposure_summary": {
        "mean": float(firms["Exposure"].mean()),
        "median": float(firms["Exposure"].median()),
        "minimum": float(firms["Exposure"].min()),
        "maximum": float(firms["Exposure"].max()),
        "std_dev": float(firms["Exposure"].std()),
        "zero_exposure_firms": int(firms["Exposure"].eq(0).sum()),
        "full_exposure_firms": int(np.isclose(firms["Exposure"], 1, atol=1e-12, rtol=0).sum()),
        "maximum_share_sum_error": float(np.max(np.abs(share_sum - 1))),
        "maximum_raw_exposure_difference": float(np.max(np.abs(firms["Exposure"] - firms["RawExposure"]))),
    },
    "loan_firm_fe_effect": loan_effect,
    "firm_credit_adjusted_effect": credit_effect,
    "firm_investment_adjusted_effect": investment_effect,
    "credit_to_loan_magnitude_ratio": float(abs(credit_effect / loan_effect)),
    "investment_to_credit_magnitude_ratio": float(abs(investment_effect / credit_effect)),
    "exposure_quality_correlation": float(firms[["Exposure", "Char4_FirmCreditQuality"]].corr().iloc[0, 1]),
    "std_error_method": "HC3 with Student t inference using N minus K degrees of freedom",
}
(RESULTS / "Question_4_Results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
print("Question 4: firm regressions and magnitude comparison saved.")
