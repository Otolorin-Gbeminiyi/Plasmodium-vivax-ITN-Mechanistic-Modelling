# Mechanistic modelling of *Plasmodium vivax* transmission and ITN effects

**R code accompanying a mechanistic modelling analysis of insecticide-treated nets (ITNs) in Ixcán, Guatemala.**

## Overview

This repository contains the original, unmodified analytical R script for Chapter 6 of the thesis work by Gbeminiyi R. Otolorin, intended to support a mathematical modelling manuscript. The script implements a deterministic human SEIT–mosquito SEI system with a relapse-eligible carrier state (`H_h`, represented as a subset of susceptible humans), a beta-binomial observation likelihood and joint calibration to two trial arms.

Key components implemented in the supplied script include:

- Preparation and validation of the control and intervention screening datasets, with October observations treated as passive-surveillance records and excluded from the primary likelihood.
- Human–mosquito differential equations, initial-state construction, burn-in and ODE solution using `deSolve`.
- Joint beta-binomial maximum-likelihood calibration of the baseline contact parameter (`b_m`), the ITN multiplier (`q_itn`) and observation precision (`phi_obs`), using 40 optimisation starts for the primary model.
- A logit-transformed ITN multiplier within the numerical parameter bounds of **0 < q_itn < 2**, permitting a protective effect (`q_itn < 1`), no effect (`q_itn = 1`) or increased effective contact (`q_itn > 1`). The direction of effect is not hard-coded as protective.
- Profile-likelihood assessment for `q_itn`, calibration diagnostics and model-specific reproduction-number calculations.
- One-way sensitivity refits, effective-contact and ITN-coverage scenarios, and five-year model-based projections.
- A biological-parameter sensitivity ensemble comprising 200 intended joint ±10% draws of 11 fixed biological parameters, each refitted using 10 starting values. These empirical percentiles are **not** a conventional confidence interval or Bayesian posterior.

## Repository contents

| Path | Description |
|---|---|
| `scripts/Chapter6_Mechanistic_Modelling.R` | Complete analytical script; content preserved byte-for-byte from the author-supplied file |
| `data/README.md` | Required input fields and data-access limitations |
| `figures/README.md` | Guidance on publishing selected figures |
| `results/README.md` | Guidance on generated result tables and diagnostics |
| `R-environment.md` | Package requirements and recording exact software versions |
| `CITATION.cff` | Draft software citation metadata; check contributor details before publication |
| `.gitignore` | Common R temporary files and locally generated/private data outputs excluded from normal Git tracking |

No underlying trial observations or fitted outputs are included in this package.

## Requirements

- R and preferably RStudio for interactive file selection.
- Required CRAN packages: `deSolve`, `dplyr`, `ggplot2`, `lubridate`, `scales` and `patchwork`.
- Two **authorised** CSV files: one control-arm and one intervention-arm dataset.

The script attempts to install missing packages. Its exact R and package versions have **not** been established for this deposit; record them from the environment in which the analysis was validated. See [R-environment.md](R-environment.md).

## Input data

At the beginning of the script, `file.choose()` opens two dialogs **in this order**:

1. Select the control-arm CSV.
2. Select the intervention-arm CSV.

Each input file must have the following case-sensitive columns:

| Column | Description |
|---|---|
| `Date_of_the_survey` | Survey date (use unambiguous ISO `YYYY-MM-DD` dates if possible) |
| `number_tested` | Number screened on that record/date |
| `number_positive` | Number of positive screening results; must not exceed `number_tested` |

The script optionally reads `cum_positive`, but it is not required. Data are aggregated by date unless `ROUND_COLUMN` is configured. Primary fitting diagnostics enforce the year 2001. By default, all October observations are omitted from the primary calibration likelihood and retained descriptively. See [data/README.md](data/README.md).

**Data are not included** because access or redistribution permissions must be confirmed with the trial data custodians and relevant institutions. A deidentified dataset is not automatically approved for unrestricted public release.

## How to run

1. Download or clone the repository and open RStudio.
2. Set the working directory to the repository root; for example, open an RStudio Project there, or use RStudio's Session > Set Working Directory > Choose Directory.
3. Ensure you can access the two original input CSV files with the required field names.
4. From an interactive R session run:

   ```r
   source("scripts/Chapter6_Mechanistic_Modelling.R", echo = FALSE)
   ```

5. Select the **control** CSV when prompted, followed by the **intervention** CSV.
6. Allow the run to complete. The primary multistart calibration, approximately 200 profile grid evaluations, one-way sensitivity refits and 200 biological-parameter uncertainty refits can require substantial computational time depending on the computer.

The R script uses the current working directory for many generated CSV and PNG files, and places uncertainty ensemble outputs in `Chapter6_200_refit_uncertainty_qITN_0_to_2/`. The script records session information in `runtime_session_info.txt` after the core analyses. The random seeds are defined in the script.

**Important:** This source code has been *packaged* for GitHub, not executed or independently validated in this environment. A run with the authorised original inputs is required to verify that all current analyses reproduce the manuscript outputs.

## Main script sections

| Section | Purpose |
|---|---|
| 6.1–6.2 | Model configuration and screening data preparation |
| 6.3–6.4 | Mechanistic model and observation likelihood |
| 6.5–6.6 | Numerical optimisation and calibration diagnostics |
| 6.7 | ITN multiplier profile likelihood |
| 6.8–6.9 | Optional hypnozoite sensitivity and parameter reporting |
| 6.10–6.12 | Transmission potential, long-run projections and coverage scenarios |
| 6.13–6.15 | One-way sensitivity, final checks and runtime diagnostics |
| 6.16 | Biological-parameter uncertainty: 200 intended refits |

## Reproducibility and limitations

Results are conditional on the specified mechanistic equations, fixed biological parameters, data ascertainment, likelihood, intervention structure and optimisation settings. The optional hypnozoite-parameter sensitivity block is disabled in the supplied script (`RUN_HYPNOZOITE_SENSITIVITY <- FALSE`), whereas the primary ITN profile, limited one-way sensitivity analysis and 200-refit ensemble are enabled.

Full independent numerical reproduction is dependent on authorised access to the screening input files and recording the software versions used to generate the final results.

## Code and data availability

The analytic R source is available in this repository. Original trial data are **not publicly deposited here**. Requests for access should be directed to the relevant data custodian, subject to applicable ethics, institutional and data-sharing requirements. Do not publish raw inputs, checkpoint files, machine-specific paths or generated tables containing non-public observations without permission.

## Citation and licence

A draft `CITATION.cff` file is supplied. Confirm contributor names and the final manuscript title, then add the repository URL, release information and archival DOI after creating a permanent release (for example through Zenodo).

**No software licence is supplied.** Confirm code ownership and institutional or collaborator approvals before choosing a licence or making the repository public. Without a licence, others generally lack permission to reuse or modify the code beyond applicable legal exceptions.
