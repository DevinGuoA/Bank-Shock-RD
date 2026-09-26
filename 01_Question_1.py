### Setup ###

from pathlib import Path
import json
import os

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "tmp" / "matplotlib"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUTPUT = ROOT / "output"
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "tmp" / "results"
for folder in [OUTPUT, PROCESSED, RESULTS]:
    folder.mkdir(parents=True, exist_ok=True)

data = pd.read_csv(ROOT / "data" / "PhdFinance_FNCE7020_HW2_Data.csv")
source_columns = data.columns.tolist()


### Data Checks ###

# One row is a firm-bank lending relationship, not a unique firm or bank.
bank_columns = [
    "Char1_BankAssets", "Char2_BankSupervised", "Treatment_BankShock"
]
firm_columns = [
    "Char4_FirmCreditQuality", "Outcome_FirmCreditGrowth",
    "Outcome_FirmInvestmentGrowth"
]
bank_variation = data.groupby("BankID")[bank_columns].nunique(dropna=False)
firm_variation = data.groupby("FirmID")[firm_columns].nunique(dropna=False)
bank_inconsistency = (bank_variation > 1).sum()
firm_inconsistency = (firm_variation > 1).sum()
loan_count = data.groupby("FirmID").size()
bank_count = data.groupby("FirmID")["BankID"].nunique()
share_sum = data.groupby("FirmID")["Char3_PreCreditShare"].sum()
threshold_mismatch = (
    data["Char2_BankSupervised"] != (data["Char1_BankAssets"] >= 30)
)

audit = {
    "loans": int(len(data)),
    "banks": int(data["BankID"].nunique()),
    "firms": int(data["FirmID"].nunique()),
    "firms_two_or_more_banks": int((loan_count >= 2).sum()),
    "firms_one_bank": int((loan_count == 1).sum()),
    "max_banks_per_firm": int(loan_count.max()),
    "missing_cells": int(data.isna().sum().sum()),
    "missing_by_variable": {c: int(data[c].isna().sum()) for c in source_columns},
    "duplicate_loan_ids": int(data["LoanID"].duplicated().sum()),
    "duplicate_firm_bank_pairs": int(data.duplicated(["FirmID", "BankID"]).sum()),
    "inconsistent_bank_variables": {c: int(bank_inconsistency[c]) for c in bank_columns},
    "inconsistent_firm_variables": {c: int(firm_inconsistency[c]) for c in firm_columns},
    "supervision_rule_violations": int(threshold_mismatch.sum()),
    "zero_credit_shares": int(data["Char3_PreCreditShare"].eq(0).sum()),
    "firms_with_share_rounding_error": int((share_sum.sub(1).abs() > 1e-12).sum()),
    "max_credit_share_sum_error": float(share_sum.sub(1).abs().max()),
    "min_credit_share_sum": float(share_sum.min()),
    "max_credit_share_sum": float(share_sum.max()),
}

# Stop on structural errors rather than silently deleting observations.
assert audit["missing_cells"] == 0
assert np.isfinite(data.to_numpy()).all()
assert audit["duplicate_loan_ids"] == 0
assert audit["duplicate_firm_bank_pairs"] == 0
assert loan_count.equals(bank_count)
assert (bank_inconsistency == 0).all() and (firm_inconsistency == 0).all()
assert data["Char2_BankSupervised"].isin([0, 1]).all()
assert data["Treatment_BankShock"].isin([0, 1]).all()
assert audit["supervision_rule_violations"] == 0
assert data["Char3_PreCreditShare"].between(0, 1).all()
assert (share_sum > 0).all()
assert audit["max_credit_share_sum_error"] < 0.000031


### Cleaned Data for Reuse ###

# Preserve rounded source shares and add weights summing to one within firm.
data["PreCreditShareNormalized"] = (
    data["Char3_PreCreditShare"] / data["FirmID"].map(share_sum)
)
data["LoanCount"] = data["FirmID"].map(loan_count)
data.to_csv(PROCESSED / "Loans_Cleaned.csv", index=False)

# Each bank receives one row, with an unweighted average over its loans.
banks = data.groupby("BankID", as_index=False).agg(
    Char1_BankAssets=("Char1_BankAssets", "first"),
    Char2_BankSupervised=("Char2_BankSupervised", "first"),
    Treatment_BankShock=("Treatment_BankShock", "first"),
    Mean_LoanGrowth=("Outcome_LoanGrowth", "mean"),
    Mean_FirmCreditQuality=("Char4_FirmCreditQuality", "mean"),
    LoanCount=("LoanID", "size"),
)
banks.to_csv(PROCESSED / "Banks_Cleaned.csv", index=False)

# Exposure is the pre-period-credit-weighted fraction of shocked lenders.
raw_exposure = (
    data["Char3_PreCreditShare"] * data["Treatment_BankShock"]
).groupby(data["FirmID"]).sum()
normalized_exposure = (
    data["PreCreditShareNormalized"] * data["Treatment_BankShock"]
).groupby(data["FirmID"]).sum()
firms = data.groupby("FirmID", as_index=False)[firm_columns].first()
firms["LoanCount"] = firms["FirmID"].map(loan_count)
firms["Exposure_raw"] = firms["FirmID"].map(raw_exposure)
firms["Exposure"] = firms["FirmID"].map(normalized_exposure)
firms.to_csv(PROCESSED / "Firms_Cleaned.csv", index=False)


### Summary Statistics ###

# Table 1 uses all source rows, so bank and firm attributes repeat by loan.
summary = data[source_columns].describe(percentiles=[0.25, 0.5, 0.75]).T
summary = summary.rename(columns={"50%": "median", "25%": "q25", "75%": "q75"})
summary = summary[["mean", "median", "min", "max", "q25", "q75", "std", "count"]]
# IDs remain in saved diagnostics but not in the displayed economic summaries.
table_summary = summary.drop(index=["LoanID", "FirmID", "BankID"])
# Exposure is summarized once per firm, matching the Question 4 sample.
exposure_summary = firms[["Exposure"]].describe(percentiles=[0.25, 0.5, 0.75]).T
exposure_summary = exposure_summary.rename(
    columns={"50%": "median", "25%": "q25", "75%": "q75"}
)[summary.columns]
bank_count_distribution = bank_count.value_counts().sort_index()
indicator_columns = ["Char2_BankSupervised", "Treatment_BankShock"]
loan_shares = data[indicator_columns].mean()
bank_shares = banks[indicator_columns].mean()

labels = {
    "LoanID": "Loan ID",
    "FirmID": "Firm ID",
    "BankID": "Bank ID",
    "Char1_BankAssets": "Bank assets",
    "Char2_BankSupervised": "Supervised",
    "Treatment_BankShock": "Bank shock",
    "Char3_PreCreditShare": "Pre credit share",
    "Char4_FirmCreditQuality": "Firm credit quality",
    "Outcome_LoanGrowth": "Loan growth",
    "Outcome_FirmCreditGrowth": "Firm credit growth",
    "Outcome_FirmInvestmentGrowth": "Investment growth",
    "Exposure": "Shock exposure",
}


### Tables ###

lines = [
    r"\begin{table}[H]",
    r"\centering",
    r"\caption{Summary Statistics}",
    r"\label{tab:q1_summary}",
    r"\begingroup",
    r"\footnotesize",
    r"\setlength{\tabcolsep}{2pt}",
    r"\renewcommand{\arraystretch}{1.15}",
    r"\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lrrrrrrrr@{}}",
    r"\toprule",
    r"Variable & Mean & Median & Min & Max & P25 & P75 & SD & $N$ \\",
    r"\midrule",
]
for column, row in pd.concat([table_summary, exposure_summary]).iterrows():
    cells = [f"{row[key]:,.3f}" for key in ["mean", "median", "min", "max", "q25", "q75", "std"]]
    cells.append(f"{int(row['count']):,}")
    lines.append(labels[column] + " & " + " & ".join(cells) + r" \\")
lines.extend([
    r"\bottomrule",
    r"\end{tabular*}",
    r"\par\smallskip",
    r"\begin{minipage}{\textwidth}\footnotesize",
    r"\textit{Notes:} Source variables weight loans equally, so firm and bank characteristics repeat across loans. Exposure counts each firm once. "
    + f"Nonmissing IDs identify {audit['loans']:,} unique loans, {audit['firms']:,} firms, and {audit['banks']:,} banks. "
    + r"$Exposure_f=\sum_b\widetilde w_{fb}Shock_b$, with $\widetilde w_{fb}=w_{fb}/\sum_jw_{fj}$, where $w_{fb}$ is the firm's preperiod borrowing share from bank $b$ and $Shock_b$ indicates a funding shock. Normalized shares sum to one within firm. SD is the sample standard deviation. P25 and P75 use linear interpolation. No observations are trimmed.",
    r"\end{minipage}",
    r"\endgroup",
    r"\end{table}",
])
(OUTPUT / "Table_1_1_Summary.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")

lines = [
    r"\begin{table}[H]",
    r"\centering",
    r"\begin{minipage}{0.85\textwidth}",
    r"\centering",
    r"\caption{Sample Composition, Supervision, and Funding Shocks}",
    r"\label{tab:q1_sample}",
    r"\begingroup\small",
    r"\renewcommand{\arraystretch}{1.15}",
    r"\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}lrr@{}}",
    r"\toprule",
    r"\multicolumn{3}{l}{\textit{Panel A: Sample counts}} \\",
    r"\midrule",
]
for label, key in [
    ("Loans", "loans"), ("Banks", "banks"), ("Firms", "firms"),
]:
    lines.append(label + r" & \multicolumn{2}{r}{" + f"{audit[key]:,}" + r"} \\")
# Indented subgroups partition all unique firms, not loan relationships.
for label, key in [
    ("Firms borrowing from at least two banks", "firms_two_or_more_banks"),
    ("Firms borrowing from one bank", "firms_one_bank"),
]:
    firm_percentage = 100 * audit[key] / audit["firms"]
    lines.append(
        r"\hspace*{1.5em}" + label + r" & \multicolumn{2}{r}{"
        + f"{audit[key]:,} ({firm_percentage:.3f}" + r"\%)} \\")
lines.extend([
    r"\midrule",
    r"\multicolumn{3}{l}{\textit{Panel B: Supervision and funding shock percentages}} \\",
    r"\midrule",
    r"Bank status & Loan relationships (\%) & Unique banks (\%) \\",
])
for column in indicator_columns:
    lines.append(
        labels[column] + f" & {100 * loan_shares[column]:.3f} & {100 * bank_shares[column]:.3f}" + r" \\")
lines.extend([
    r"\bottomrule",
    r"\end{tabular*}",
    r"\par\smallskip",
    r"\begin{minipage}{\linewidth}\footnotesize",
    r"\textit{Notes:} A loan is a unique firm and bank relationship. Percentages in Panel A use all 41,000 firms as the denominator. In Panel B, the first percentage column uses all 161,704 loan relationships as the denominator and reports the share involving a bank with the stated status. The second uses all 320 unique banks, counting each bank once. Supervised means assets of at least 30; Bank shock means the bank experienced the funding shock. Each status is considered separately.",
    r"\end{minipage}",
    r"\endgroup",
    r"\end{minipage}",
    r"\end{table}",
])
(OUTPUT / "Table_1_2_Sample.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")


### Distribution Figure ###

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titleweight": "semibold", "axes.titlesize": 11,
    "figure.facecolor": "white", "savefig.facecolor": "white",
})
fig, axes = plt.subplots(
    2, 2, figsize=(10.4, 5.5), sharex="col",
    gridspec_kw={"height_ratios": [4, 1], "hspace": 0.09, "wspace": 0.25},
)
plot_series = [data["Outcome_LoanGrowth"], banks["Char1_BankAssets"]]
plot_titles = [
    f"A. Loan growth ({len(data):,} loans)",
    f"B. Bank assets ({len(banks):,} unique banks)",
]
plot_labels = ["Change in log lending", "Bank assets"]
for index, series in enumerate(plot_series):
    axes[0, index].hist(series, bins="fd", density=True, color="#2B5870", alpha=0.85, edgecolor="white", linewidth=0.5)
    axes[0, index].set_title(plot_titles[index], loc="left", pad=10)
    axes[0, index].set_ylabel("Density")
    axes[0, index].grid(axis="y", color="#E5E5E5", linewidth=0.6)
    axes[0, index].set_axisbelow(True)
    axes[1, index].boxplot(
        series, vert=False, widths=0.45, patch_artist=True,
        boxprops={"facecolor": "#C7D6DF", "edgecolor": "#2B5870"},
        medianprops={"color": "#203D4A", "linewidth": 1.5},
        whiskerprops={"color": "#2B5870"}, capprops={"color": "#2B5870"},
        flierprops={"marker": ".", "markersize": 3, "markerfacecolor": "#2B5870", "markeredgecolor": "none", "alpha": 0.5},
    )
    axes[1, index].set_yticks([])
    axes[1, index].spines["left"].set_visible(False)
    axes[1, index].set_xlabel(plot_labels[index])
axes[0, 1].axvline(30, color="#A64D3E", linestyle="--", linewidth=1.3, label="Supervision cutoff")
axes[0, 1].legend(frameon=False, fontsize=9)
fig.subplots_adjust(left=0.075, right=0.985, top=0.91, bottom=0.12)
fig.savefig(OUTPUT / "Figure_1_Distributions.png", dpi=300, bbox_inches="tight")
plt.close(fig)


### Number of Banks per Firm ###

# Count each firm once, regardless of its number of loan relationships.
fig, ax = plt.subplots(figsize=(9.0, 3.4))
bars = ax.bar(
    bank_count_distribution.index, bank_count_distribution.values,
    width=0.65, color="#2B5870",
)
ax.bar_label(
    bars, labels=[f"{count:,}" for count in bank_count_distribution.values],
    padding=4, fontsize=10,
)
ax.set_xlabel("Number of banks per firm")
ax.set_ylabel("Number of firms")
ax.set_xticks(bank_count_distribution.index)
ax.set_ylim(0, bank_count_distribution.max() * 1.18)
ax.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter("{x:,.0f}"))
ax.grid(axis="y", color="#E5E5E5", linewidth=0.6)
ax.set_axisbelow(True)
fig.subplots_adjust(left=0.10, right=0.985, top=0.95, bottom=0.18)
fig.savefig(OUTPUT / "Figure_1_2_Banks_Per_Firm.png", dpi=300, bbox_inches="tight")
plt.close(fig)


### Numerical Results ###

audit["max_normalized_share_sum_error"] = float(
    data.groupby("FirmID")["PreCreditShareNormalized"].sum().sub(1).abs().max()
)
results = {
    "audit": audit,
    "summary_loan_level": json.loads(summary.to_json(orient="index")),
    "summary_exposure_firm_level": json.loads(exposure_summary.to_json(orient="index")),
    "banks_per_firm_counts": {str(k): int(v) for k, v in bank_count_distribution.items()},
    "banks_per_firm_summary": bank_count.describe().to_dict(),
    "indicator_shares_loan_rows": {c: float(loan_shares[c]) for c in indicator_columns},
    "indicator_shares_unique_banks": {c: float(bank_shares[c]) for c in indicator_columns},
    "indicator_counts_unique_banks": {c: int(banks[c].sum()) for c in indicator_columns},
    "bank_assets_unique_bank_summary": banks["Char1_BankAssets"].describe().to_dict(),
    "loan_growth_skewness": float(data["Outcome_LoanGrowth"].skew()),
    "bank_assets_skewness": float(banks["Char1_BankAssets"].skew()),
    "raw_exposure_min": float(raw_exposure.min()),
    "raw_exposure_max": float(raw_exposure.max()),
    "normalized_exposure_min": float(normalized_exposure.min()),
    "normalized_exposure_max": float(normalized_exposure.max()),
    "max_exposure_normalization_change": float((raw_exposure - normalized_exposure).abs().max()),
}
(RESULTS / "Question_1_Results.json").write_text(
    json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8"
)
print("Question 1: summary tables, two figures, and three cleaned datasets saved.")
