Python analysis of regulatory supervision, bank funding shocks, and their transmission to firms. The project includes a sharp regression discontinuity design, within firm loan regressions, and firm credit and investment regressions. Written answers are maintained separately in LaTeX.

## Project structure

| Directory | Contents |
| --- | --- |
| `code/` | Four question scripts, main runner, dependencies, and LaTeX sources |
| `data/` | Original course CSV, never overwritten by the scripts |
| `data/processed/` | Cleaned loan, bank, and firm datasets for HW3 |
| `output/` | PNG figures, TeX tables, and the separately compiled report |
| `README/` | Assignment PDF and reproduction guides in Markdown, HTML, and PDF |
| `tmp/` | Full precision JSON results, caches, and document build files |

This Markdown file can be placed at a GitHub repository root. All commands below are run from the project root, not from the `README/` directory. No HW2 repository URL has been supplied; the report leaves a GitHub placeholder.

## Data and observation levels

The simulated dataset is supplied with the course assignment:

```text
data/PhdFinance_FNCE7020_HW2_Data.csv
README/HW2_FNCE7020_Fall2026.pdf
```

There are 161,704 loan relationships, 41,000 firms, and 320 banks. A row is one firm and bank pair. There are no missing source values, duplicate loan identifiers, or duplicate firm and bank pairs. A total of 38,826 firms borrow from at least two banks. Bank attributes are constant within bank, and firm attributes are constant within firm.

| Variable | Meaning |
| --- | --- |
| `LoanID` | Unique loan relationship identifier |
| `FirmID` | Firm identifier |
| `BankID` | Bank identifier |
| `Char1_BankAssets` | Bank assets in arbitrary source units |
| `Char2_BankSupervised` | One if bank assets are at least 30 |
| `Treatment_BankShock` | One if the bank experiences the funding shock |
| `Char3_PreCreditShare` | Firm's preperiod borrowing share from this bank |
| `Char4_FirmCreditQuality` | Firm credit quality, with higher values indicating stronger quality |
| `Outcome_LoanGrowth` | Change in log lending from the bank to the firm |
| `Outcome_FirmCreditGrowth` | Change in the firm's total borrowing |
| `Outcome_FirmInvestmentGrowth` | Change in the firm's investment |

Only loan growth is explicitly defined as a log change in the assignment. Firm outcomes retain their source growth units. ID statistics are mechanical. Summaries over loan rows weight repeated bank and firm attributes by their numbers of loans; bank shares and the asset histogram instead use one row per bank.

Source CSV SHA256:

```text
638ed7359180c0a595f0d320b9b17c8bc7285724d977f670762abfea9f1148bb
```

## Install and reproduce

The analysis was tested with Python 3.9.6. Use Python 3.9 through 3.12 with these pinned package versions. A virtual environment is recommended:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r code/requirements.txt
python code/00_Main.py
```

On Windows, activate with `.venv\Scripts\activate` instead. `requirements.txt` is needed to recreate the tested direct dependencies, including `rdrobust` and `rddensity`; keep it with the project. It pins direct dependencies, not every transitive package or system library.

The main script runs Questions 1 through 4 and stops if any script fails. It uses the active Python interpreter. Each script can also run independently, for example:

```sh
python code/02_Question_2.py
```

Paths are resolved relative to the script location, so the analysis does not rely on a particular working directory. Scripts use the full source sample without random subsampling. Running again overwrites generated analysis outputs and processed datasets, but never the original CSV.

## What each script does

| Script | Analysis |
| --- | --- |
| `00_Main.py` | Runs the four analysis scripts only |
| `01_Question_1.py` | Audits the source, saves cleaned datasets, summarizes all variables, and plots distributions |
| `02_Question_2.py` | Estimates the bank level RD at assets of 30, varies the window, runs `rdrobust`, and checks balance and density |
| `03_Question_3.py` | Estimates pooled and firm fixed effect loan regressions and documents credit quality sorting |
| `04_Question_4.py` | Constructs firm exposure and estimates four firm outcome regressions |

Python generates analysis tables, figures, and machine readable results only. It does not write answer fragments or create or compile the report or README documents.

## Analytical choices and limitations

1. Source columns and rows are preserved. There is no imputation, trimming, or winsorization. Three recorded precredit shares equal zero and are retained. Within firm, rounded source shares differ from one by at most 0.00003. The primary exposure uses shares divided by their firm sum. Raw exposure is retained, and Question 4 verifies that the raw share formula yields essentially identical coefficients.
2. Question 2 averages loan growth within bank and then weights banks equally. The manual local linear regressions use separate slopes, uniform windows of 3, 5, and 8, and HC3 standard errors. The automatic specification uses `rdrobust`, a triangular kernel, local linear fitting, quadratic bias correction, and MSE bandwidth selection. Its robust bias corrected confidence interval is the primary inference.
3. The asset density check uses `rddensity`. The available balance checks use shock status, mean borrower credit quality, and loan count. Their status as predetermined bank characteristics is not established by the assignment. These are limited diagnostics, not a substitute for genuinely predetermined bank balance sheet characteristics. Failure to reject a discontinuity does not prove identification.
4. Question 3 reports pooled results for all loans, a pooled comparison restricted to multiple bank firms, a pooled model controlling for quality, and firm fixed effects. The fixed effect is estimated by within firm demeaning, without constructing thousands of dummy columns. Standard errors are clustered by bank and firm. The finite sample correction includes absorbed firm effects, and tests use 319 degrees of freedom. Firms with no within firm variation in shock status do not identify the shock coefficient.
5. Question 4 uses one observation per firm, equal firm weights, and HC3 standard errors. HC3 does not account for dependence across firms sharing lenders; the reported significance is conditional on independent firm errors. Credit quality is a useful observable control, not proof of conditional random exposure or exclusion of other bank related channels.
6. All displayed statistics use three decimal places, with integer counts and very small p values shown as `<0.001`. Regression tables use significance stars and standard errors in parentheses. Full precision values are saved in `tmp/results/Question_*_Results.json`.

## Cleaned datasets for HW3

`01_Question_1.py` creates:

- `data/processed/Loans_Cleaned.csv`: all source columns, `PreCreditShareNormalized`, and firm `LoanCount`.
- `data/processed/Banks_Cleaned.csv`: one row per bank, bank assets, supervision and shock indicators, mean loan growth, mean borrower quality, and loan count.
- `data/processed/Firms_Cleaned.csv`: one row per firm, credit quality, both firm outcomes, loan count, `Exposure_raw`, and normalized `Exposure`.

Keep the original CSV alongside these files so that every transformation remains reproducible.

## Compile the written report separately

Written answers are in `code/HW2_Concise_Report.tex`. Numerical statements are edited manually, so update them if the data or model changes. A TeX installation with `latexmk` and `pdflatex` is required. After running the analysis:

```sh
mkdir -p tmp/documents
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=tmp/documents code/HW2_Concise_Report.tex
cp tmp/documents/HW2_Concise_Report.pdf output/HW2_Concise_Report.pdf
```

The report starts each question on a new page and references the generated tables and figures. Its title identifies Devin Guo as the author.

The README PDF is maintained in `code/README.tex` and compiled separately:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=tmp/documents code/README.tex
cp tmp/documents/README.pdf README/README.pdf
```

`README/README.html` is the browser friendly version. Neither documentation build is part of `00_Main.py`.
