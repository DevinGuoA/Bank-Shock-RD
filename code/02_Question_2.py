### Setup ###

import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "output"
RESULTS = ROOT / "tmp" / "results"
OUTPUT.mkdir(exist_ok=True)
RESULTS.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "tmp" / "matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from rddensity import rddensity
from rdrobust import rdrobust

plt.rcParams.update({
    "font.family": "DejaVu Serif",
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 150,
    "savefig.dpi": 300,
})
BLUE = "#245780"
ORANGE = "#B56B32"

### Bank Level Sample ###

data = pd.read_csv(ROOT / "data" / "PhdFinance_FNCE7020_HW2_Data.csv")
bank_columns = ["Char1_BankAssets", "Char2_BankSupervised", "Treatment_BankShock"]
assert data.notna().all().all(), "Resolve missing observations before estimation."
assert data.groupby("BankID")[bank_columns].nunique().eq(1).all().all()
assert not data.duplicated(["BankID", "FirmID"]).any()

# Each bank receives equal weight after averaging its loan observations.
banks = data.groupby("BankID", sort=True).agg(
    assets=("Char1_BankAssets", "first"),
    supervised=("Char2_BankSupervised", "first"),
    shock=("Treatment_BankShock", "first"),
    loan_growth=("Outcome_LoanGrowth", "mean"),
    borrower_quality=("Char4_FirmCreditQuality", "mean"),
    loan_count=("LoanID", "size"),
).reset_index()
assert banks["supervised"].eq((banks["assets"] >= 30).astype(int)).all()

# Assets use the assignment's artificial units. Growth remains in log units.
banks["centered_assets"] = banks["assets"] - 30
banks["shock_percent"] = 100 * banks["shock"]

### Local Linear RD ###

window_results = []
window_models = {}
for bandwidth in [3, 5, 8]:
    local = banks.loc[banks["centered_assets"].abs() <= bandwidth].copy()
    # Columns: intercept, supervision, centered assets, supervision by assets.
    design = np.column_stack([
        np.ones(len(local)),
        local["supervised"],
        local["centered_assets"],
        local["supervised"] * local["centered_assets"],
    ])
    # Uniform kernel, separate slopes, HC3 standard errors, finite sample t tests.
    fit = sm.OLS(local["loan_growth"].to_numpy(), design).fit(
        cov_type="HC3", use_t=True
    )
    window_models[bandwidth] = fit
    window_results.append({
        "bandwidth": bandwidth,
        "estimate": float(fit.params[1]),
        "standard_error": float(fit.bse[1]),
        "t_statistic": float(fit.tvalues[1]),
        "p_value": float(fit.pvalues[1]),
        "ci_lower": float(fit.conf_int()[1, 0]),
        "ci_upper": float(fit.conf_int()[1, 1]),
        "n_left": int((local["supervised"] == 0).sum()),
        "n_right": int((local["supervised"] == 1).sum()),
        "degrees_of_freedom": float(fit.df_resid),
        "coefficients": fit.params.tolist(),
    })

### Automatic Bandwidth and Robust Bias Correction ###

robust_fit = rdrobust(
    banks["loan_growth"].to_numpy(),
    banks["assets"].to_numpy(),
    c=30, p=1, q=2, kernel="tri", bwselect="mserd", vce="nn", nnmatch=3,
    masspoints="adjust", level=95,
)
robust_results = {
    "conventional_estimate": float(robust_fit.coef.loc["Conventional"].iloc[0]),
    "conventional_standard_error": float(robust_fit.se.loc["Conventional"].iloc[0]),
    "conventional_p_value": float(robust_fit.pv.loc["Conventional"].iloc[0]),
    "conventional_ci": robust_fit.ci.loc["Conventional"].tolist(),
    "bias_corrected_estimate": float(robust_fit.coef.loc["Robust"].iloc[0]),
    "robust_standard_error": float(robust_fit.se.loc["Robust"].iloc[0]),
    "robust_z_statistic": float(robust_fit.z.loc["Robust"].iloc[0]),
    "robust_p_value": float(robust_fit.pv.loc["Robust"].iloc[0]),
    "robust_ci": robust_fit.ci.loc["Robust"].tolist(),
    "h_left": float(robust_fit.bws.loc["h", "left"]),
    "h_right": float(robust_fit.bws.loc["h", "right"]),
    "b_left": float(robust_fit.bws.loc["b", "left"]),
    "b_right": float(robust_fit.bws.loc["b", "right"]),
    "n_h_left": int(robust_fit.N_h[0]),
    "n_h_right": int(robust_fit.N_h[1]),
    "n_b_left": int(robust_fit.N_b[0]),
    "n_b_right": int(robust_fit.N_b[1]),
    "kernel": "triangular",
    "bandwidth_selection": "mserd",
    "variance_estimator": "nearest neighbor, nnmatch=3",
}

### Limited Balance Diagnostics ###

local = banks.loc[banks["centered_assets"].abs() <= 5].copy()
design = np.column_stack([
    np.ones(len(local)), local["supervised"], local["centered_assets"],
    local["supervised"] * local["centered_assets"],
])
# These proxies are not verified pre supervision bank characteristics.
balance_results = []
balance_models = {}
for variable, label, units in [
    ("shock_percent", "Bank funding shock", "percentage points"),
    ("borrower_quality", "Mean borrower credit quality", "index units"),
    ("loan_count", "Number of loan relationships", "relationships"),
]:
    fit = sm.OLS(local[variable].to_numpy(), design).fit(cov_type="HC3", use_t=True)
    balance_models[variable] = fit
    balance_results.append({
        "variable": variable, "label": label, "units": units,
        "estimate": float(fit.params[1]),
        "standard_error": float(fit.bse[1]),
        "t_statistic": float(fit.tvalues[1]),
        "p_value": float(fit.pvalues[1]),
        "ci_lower": float(fit.conf_int()[1, 0]),
        "ci_upper": float(fit.conf_int()[1, 1]),
        "n_left": int((local["supervised"] == 0).sum()),
        "n_right": int((local["supervised"] == 1).sum()),
    })

### Bank Characteristic Figure ###

# Bins summarize points only. Fits retain all equally weighted banks in [25, 35].
plot_sample = local.copy()
# One asset unit per bin, with 30 on the right and 35 in the final bin.
plot_sample["asset_bin"] = np.minimum(
    np.floor(plot_sample["assets"] - 25).astype(int), 9
)
bin_means = plot_sample.groupby(["supervised", "asset_bin"], sort=True).agg(
    assets=("assets", "mean"),
    shock_percent=("shock_percent", "mean"),
    borrower_quality=("borrower_quality", "mean"),
    loan_count=("loan_count", "mean"),
    n_banks=("BankID", "size"),
).reset_index()
assert bin_means["n_banks"].sum() == len(plot_sample)

figure, axes = plt.subplots(1, 3, figsize=(10.6, 3.9))
panel_titles = [
    "A. Bank funding shock", "B. Mean borrower credit quality",
    "C. Number of loan relationships",
]
axis_labels = ["Shocked banks (%)", "Credit quality (index)", "Loan relationships"]
for axis, result, title, ylabel in zip(axes, balance_results, panel_titles, axis_labels):
    variable = result["variable"]
    coefficient = np.asarray(balance_models[variable].params)
    for status, color in [(0, BLUE), (1, ORANGE)]:
        selected = bin_means.loc[bin_means["supervised"] == status]
        axis.scatter(selected["assets"], selected[variable], s=38,
                     color=color, alpha=0.85, edgecolor="white", linewidth=0.5,
                     zorder=3)
        grid = np.linspace(-5 if status == 0 else 0,
                           0 if status == 0 else 5, 100)
        line = (coefficient[0] + coefficient[1] * status
                + coefficient[2] * grid + coefficient[3] * status * grid)
        axis.plot(grid + 30, line, color=color, linewidth=2)
        axis.scatter([30], [coefficient[0] + coefficient[1] * status],
                     s=35, facecolor="white", edgecolor=color, zorder=4)
    axis.axvline(30, color="#777777", linestyle=":", linewidth=1)
    axis.set_xlim(25, 35)
    axis.set_xticks([25, 27.5, 30, 32.5, 35])
    axis.set_title(title, fontsize=10)
    axis.set_xlabel("Bank assets")
    axis.set_ylabel(ylabel)
    if variable == "shock_percent":
        axis.set_ylim(0, 105)
    else:
        axis.margins(y=0.22)
    axis.text(0.97, 0.97,
              f"Jump = {result['estimate']:.3f}\n$p$ = {result['p_value']:.3f}",
              transform=axis.transAxes, ha="right", va="top", fontsize=8.5,
              bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"})
    axis.grid(axis="y", color="#dddddd", linewidth=0.5)
figure.tight_layout(pad=1.2, w_pad=1.3)
figure.savefig(OUTPUT / "Figure_2_3_Bank_Characteristics.png", bbox_inches="tight")
plt.close(figure)

### Density Continuity ###

# Test bank assets once per bank, not once per loan.
density_fit = rddensity(
    banks["assets"].to_numpy(), c=30, p=2, q=3,
    fitselect="unrestricted", kernel="triangular", vce="jackknife",
    bwselect="comb", massPoints=True, bino_flag=False,
)
density_results = {
    "estimate_left": float(density_fit.hat["left"]),
    "estimate_right": float(density_fit.hat["right"]),
    "estimate_difference": float(density_fit.hat["diff"]),
    "standard_error": float(density_fit.sd_jk["diff"]),
    "z_statistic": float(density_fit.test["t_jk"]),
    "p_value": float(density_fit.test["p_jk"]),
    "h_left": float(density_fit.h["left"]),
    "h_right": float(density_fit.h["right"]),
    "n_left": int(density_fit.n["eff_left"]),
    "n_right": int(density_fit.n["eff_right"]),
    "kernel": "triangular", "bandwidth_selection": "comb",
    "variance_estimator": "jackknife", "p": 2, "q": 3,
}

### RD Figure ###

figure, axes = plt.subplots(1, 3, figsize=(10.6, 3.8), sharey=True)
for axis, result in zip(axes, window_results):
    bandwidth = result["bandwidth"]
    local = banks.loc[banks["centered_assets"].abs() <= bandwidth]
    coefficient = np.array(result["coefficients"])
    for status, color in [(0, BLUE), (1, ORANGE)]:
        selected = local.loc[local["supervised"] == status]
        axis.scatter(selected["assets"], selected["loan_growth"],
                     s=19, color=color, alpha=0.6, edgecolor="none")
        grid = np.linspace(-bandwidth if status == 0 else 0,
                           0 if status == 0 else bandwidth, 100)
        line = (coefficient[0] + coefficient[1] * status
                + coefficient[2] * grid + coefficient[3] * status * grid)
        axis.plot(grid + 30, line, color=color, linewidth=2)
        axis.scatter([30], [coefficient[0] + coefficient[1] * status],
                     s=30, facecolor="white", edgecolor=color, zorder=4)
    axis.axvline(30, color="#777777", linestyle=":", linewidth=1)
    axis.set_title(f"Asset window: $30 \\pm {bandwidth}$")
    axis.set_xlabel("Bank assets (source units)")
    axis.text(0.97, 0.97, f"RD = {result['estimate']:.3f}\n"
              f"Banks: {result['n_left']} / {result['n_right']}",
              transform=axis.transAxes, ha="right", va="top", fontsize=8.5,
              bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"})
    axis.grid(axis="y", color="#dddddd", linewidth=0.5)
axes[0].set_ylabel("Mean change in log lending")
figure.tight_layout(pad=1.2, w_pad=1.0)
figure.savefig(OUTPUT / "Figure_2_1_RD_Fits.png", bbox_inches="tight")
plt.close(figure)

### Density Figure ###

figure, axis = plt.subplots(figsize=(7.8, 3.8))
# Common bin widths and the full bank count preserve density comparability.
edges = np.arange(10, 52, 2)
counts, _ = np.histogram(banks["assets"], bins=edges)
centers = (edges[:-1] + edges[1:]) / 2
colors = [BLUE if center < 30 else ORANGE for center in centers]
axis.bar(centers, counts / (len(banks) * 2), width=1.85,
         color=colors, alpha=0.65, edgecolor="white", linewidth=0.5)
axis.axvline(30, color="#555555", linestyle=":", linewidth=1.4)
axis.set_xlim(10, 50)
axis.set_xlabel("Bank assets")
axis.set_ylabel("Density")
axis.grid(axis="y", color="#dddddd", linewidth=0.5)
axis.text(0.98, 0.96,
          f"Density continuity test: p = {density_results['p_value']:.3f}\n"
          "Bin width: 2 asset units",
          transform=axis.transAxes, ha="right", va="top", fontsize=9,
          bbox={"facecolor": "white", "alpha": 0.9, "edgecolor": "none"})
figure.tight_layout()
figure.savefig(OUTPUT / "Figure_2_2_Asset_Density.png", bbox_inches="tight")
plt.close(figure)

### Window Regression Table ###

table = [
    r"\begin{table}[H]", r"\centering", r"\small",
    r"\caption{Local Linear Effects of Bank Supervision}",
    r"\label{tab:q2_windows}",
    r"\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lccc}",
    r"\toprule", r" & (1) & (2) & (3) \\",
    r"Asset window & $[27,33]$ & $[25,35]$ & $[22,38]$ \\",
    r"\midrule",
]
coefficient_cells = []
for result in window_results:
    p_value = result["p_value"]
    stars = "***" if p_value < 0.01 else "**" if p_value < 0.05 else "*" if p_value < 0.10 else ""
    coefficient_cells.append(f"${result['estimate']:.3f}^{{{stars}}}$")
table.append("Supervised & " + " & ".join(coefficient_cells) + r" \\")
table.append(" & " + " & ".join(f"({r['standard_error']:.3f})" for r in window_results) + r" \\")
table.append("$t$ statistic & " + " & ".join(f"{r['t_statistic']:.3f}" for r in window_results) + r" \\")
table.append("$p$ value & " + " & ".join("$<0.001$" if r["p_value"] < .001 else f"{r['p_value']:.3f}" for r in window_results) + r" \\")
table.append(r"\midrule")
for label, key in [("Banks below cutoff", "n_left"), ("Banks above cutoff", "n_right")]:
    table.append(label + " & " + " & ".join(f"{r[key]:,}" for r in window_results) + r" \\")
table.extend([
    r"Separate asset slopes & Yes & Yes & Yes \\",
    r"\bottomrule", r"\end{tabular*}",
    r"\begin{minipage}{\textwidth}\footnotesize",
    r"\textit{Notes:} Specification: $\overline{Y}_b=\alpha+\tau S_b+\beta x_b+\gamma S_bx_b+u_b$, where $x_b=\mathrm{Assets}_b-30$ and $S_b=1\{\mathrm{Assets}_b\geq30\}$. The dependent variable is the bank mean change in log lending, kept in source units. Each bank has equal weight within the stated window, equivalent to a uniform kernel. HC3 standard errors are in parentheses. Tests and confidence intervals use a $t$ distribution with $N-4$ degrees of freedom. These manual estimates do not apply bias correction. $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$.",
    r"\end{minipage}", r"\end{table}",
])
(OUTPUT / "Table_2_1_RD_Windows.tex").write_text("\n".join(table) + "\n", encoding="utf-8")

### Robust RD Table ###

table = [
    r"\begin{table}[H]", r"\centering", r"\small",
    r"\caption{Automatic Bandwidth Selection and Robust RD Inference}",
    r"\label{tab:q2_robust}",
    r"\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lcc}",
    r"\toprule", r" & Conventional & Robust bias corrected \\", r"\midrule",
]
for label, key_c, key_r in [
    ("RD estimate", "conventional_estimate", "bias_corrected_estimate"),
    ("Standard error", "conventional_standard_error", "robust_standard_error"),
    ("$p$ value", "conventional_p_value", "robust_p_value"),
]:
    cells = []
    for key in [key_c, key_r]:
        value = robust_results[key]
        formatted = f"{value:.3f}"
        if label == "RD estimate":
            p_key = "conventional_p_value" if key == key_c else "robust_p_value"
            p_value = robust_results[p_key]
            stars = "***" if p_value < .01 else "**" if p_value < .05 else "*" if p_value < .10 else ""
            formatted = f"${formatted}^{{{stars}}}$"
        elif label == "$p$ value" and value < .001:
            formatted = "$<0.001$"
        cells.append(formatted)
    table.append(label + " & " + " & ".join(cells) + r" \\")
table.append("95\\% confidence interval & " + " & ".join(
    f"$[{robust_results[key][0]:.3f},\ {robust_results[key][1]:.3f}]$"
    for key in ["conventional_ci", "robust_ci"]) + r" \\")
table.extend([
    r"\midrule", r" & Below cutoff & Above cutoff \\",
    f"Estimation bandwidth $h$ & {robust_results['h_left']:.3f} & {robust_results['h_right']:.3f}" + r" \\",
    f"Bias bandwidth $b$ & {robust_results['b_left']:.3f} & {robust_results['b_right']:.3f}" + r" \\",
    f"Banks within $h$ & {robust_results['n_h_left']} & {robust_results['n_h_right']}" + r" \\",
    f"Banks within $b$ & {robust_results['n_b_left']} & {robust_results['n_b_right']}" + r" \\",
    r"\bottomrule", r"\end{tabular*}",
    r"\begin{minipage}{\textwidth}\footnotesize",
    r"\textit{Notes:} The estimand is the jump in bank mean log loan growth at assets of 30. Estimates and standard errors use the source log units. The \texttt{rdrobust} specification uses separate local linear fits ($p=1$), local quadratic bias estimation ($q=2$), a triangular kernel, the common MSE optimal bandwidth selector \texttt{mserd}, and nearest neighbor variance estimation with three neighbors. The right column reports the bias corrected estimate with its robust standard error and normal based confidence interval. $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$.",
    r"\end{minipage}", r"\end{table}",
])
(OUTPUT / "Table_2_2_RD_Robust.tex").write_text("\n".join(table) + "\n", encoding="utf-8")

### Validity Diagnostics Table ###

table = [
    r"\begin{table}[H]", r"\centering", r"\small",
    r"\caption{Limited Balance and Density Diagnostics}",
    r"\label{tab:q2_validity}",
    r"\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lrrr}",
    r"\toprule", r"Diagnostic & Jump & Std. error & $p$ value \\",
    r"\midrule", r"\multicolumn{4}{l}{\textit{Panel A. Available bank and borrower proxies}} \\",
]
for result in balance_results:
    p_formatted = "$<0.001$" if result["p_value"] < .001 else f"{result['p_value']:.3f}"
    table.append(f"{result['label']} & {result['estimate']:.3f} & "
                 f"{result['standard_error']:.3f} & {p_formatted}" + r" \\")
table.extend([
    r"\midrule", r"\multicolumn{4}{l}{\textit{Panel B. Bank asset density}} \\",
    f"Density at the cutoff & {density_results['estimate_difference']:.3f} & "
    f"{density_results['standard_error']:.3f} & {density_results['p_value']:.3f}" + r" \\",
    r"\bottomrule", r"\end{tabular*}",
    r"\begin{minipage}{\textwidth}\footnotesize",
    r"\textit{Notes:} Panel A uses the local linear specification in Table~\ref{tab:q2_windows} within assets $[25,35]$, with HC3 standard errors and $t$ tests. The shock jump is in percentage points, credit quality is in index units, and relationships are counts. These variables are not verified pre supervision covariates: shock timing relative to supervision is unspecified, and borrower composition may respond to supervision. Panel B uses the robust \texttt{rddensity} test with local polynomial orders $p=2$, $q=3$, a triangular kernel, automatic \texttt{comb} bandwidths, and jackknife variance estimation. The density jump is measured per asset unit. All jumps are right minus left. Individual $p$ values are unadjusted. Nonrejection does not establish RD validity.",
    r"\end{minipage}", r"\end{table}",
])
(OUTPUT / "Table_2_3_Validity_Checks.tex").write_text("\n".join(table) + "\n", encoding="utf-8")

### Save Numerical Results ###

results = {
    "loan_observations": len(data),
    "bank_observations": len(banks),
    "bank_left": int((banks["assets"] < 30).sum()),
    "bank_right": int((banks["assets"] >= 30).sum()),
    "bank_assets_at_cutoff": int((banks["assets"] == 30).sum()),
    "bank_assets_unique": int(banks["assets"].nunique()),
    "nearest_assets_left": float(banks.loc[banks["assets"] < 30, "assets"].max()),
    "nearest_assets_right": float(banks.loc[banks["assets"] >= 30, "assets"].min()),
    "windows": window_results,
    "rdrobust": robust_results,
    "balance_diagnostics": balance_results,
    "density_test": density_results,
    "units": {"assets": "artificial source units", "growth": "change in log lending"},
}
(RESULTS / "Question_2_Results.json").write_text(
    json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8"
)
print("Question 2: bank RD tables, figures, and numerical results generated.")
