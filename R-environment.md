# R environment and package dependencies

The uploaded code uses the following **direct** CRAN package dependencies:

| Package | Role in the script |
|---|---|
| `deSolve` | Numerical integration of the differential equations with `lsoda` |
| `dplyr` | Data manipulation and summaries |
| `ggplot2` | Visualisation and output figures |
| `lubridate` | Survey date parsing and calendar calculations |
| `scales` | Percentage axis formatting |
| `patchwork` | Composite visualisation panels |

The script also uses R's built-in `stats`, `utils`, `tools` and base functions.

## Record the validated environment

**Do not invent package version numbers or a `renv.lock` file.** Run the analysis on the computer where the final results were obtained, and capture:

```r
R.version.string
sapply(c("deSolve", "dplyr", "ggplot2", "lubridate", "scales", "patchwork"),
       function(p) as.character(packageVersion(p)))
sessionInfo()
```

The script already writes `runtime_session_info.txt` during execution. Review the file for machine-specific information before publication.

## Optional exact dependency lockfile

If you wish to use `renv`, initialise it **in a working copy of the repository** and snapshot the environment containing the validated packages:

```r
install.packages("renv")
renv::init()
renv::snapshot()
```

Confirm that the changes do not affect a verified modelling environment. Commit `renv.lock` once generated. This package does not contain a fabricated `renv.lock` because the source file alone cannot establish the versions of R and its libraries.
