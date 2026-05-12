# GitHub repository 

Code and scripts to reproduce the simulations and analyses from our paper "A comparative overview of win ratio and joint frailty models for
recurrent event endpoints with applications in oncology and cardiology".

## Repository layout
- `Simulations/Generation.R` - Monte Carlo data generation across several scenarios; uses `helpers/dataJFM_fast.cpp` (and the optional `dataJFM_fast_heterog.cpp`) via `Rcpp`.
- `Simulations/Estimation.R` - performance assessment for the win ratio (`WR` package) and the joint frailty model (`frailtypack`), producing the tables reported in the paper. 
- `Simulations/Description.R` - quick descriptive summaries of the generated datasets.
- `Simulations/helpers/` - data generator (`genCpp.R`), C++ back end, and estimation helpers (`estimationWR.R`, `estimationJFM.R`).
- `master_script.R` - master script to facilitate the script execution.
- `Analyses/` - code for the real-data applications:
  - `HFaction analysis.R` - HF-ACTION heart failure trial.
  - `Readmission analysis.R` - colorectal cancer readmission data from `frailtypack`.
  - `Sample size based on HF-action.R` - analytic and simulation-based sample size calculations for both approaches.
  - `hfaction.csv` - modified HF-ACTION dataset 
- `WR_1.0.tar.gz` / `WR_1.0.pdf` - source bundle and manual for our modified version of the WR package.

## Requirements
- R (tested on version 4.5.2).
- Compiler toolchain for `Rcpp`.
- R packages: `WR`, `frailtypack`, `dplyr`, `knitr`, `kableExtra`, `tidyr`, `MASS`, `Rcpp`, `here`.
- Install from CRAN where available:
  ```{r}
  install.packages(c(
    "Rcpp", "dplyr", "knitr", "kableExtra",
    "tidyr", "MASS", "frailtypack", "here"
  ))
  install.packages("WR") # does not include the modification we've made to the package
  ```
- If you want to use our modified version of `WR`, you can download `WR_1.0.tar.gz` file and run the following commands:
  ```{r}
  install.packages(c("cubature", "gumbel", "survival")) # dependancies of the WR package
  install.packages("WR_1.0.tar.gz", repos = NULL, type = "source") # adapt to your file path
  ```
  Alternatively, if you don't want to install the package directly from R, run the following commands:
  ```{r}
  install.packages(c("cubature", "gumbel", "survival")) # dependancies of the WR package
  install.packages("remotes") # to run the following command
  remotes::install_url("https://github.com/adrienorue/comparison_wrRec_jfm/raw/main/WR_1.0.tar.gz")
  ```
  If installation fails, please verify that "Rtools" is installed on your system.



## Reproducing the simulations

For ease of use, `master_script.R` provides an ordered workflow for running the simulation and analysis scripts: `Simulations/Generation.R`, `Simulations/Description.R`, `Simulations/Estimation.R`, `Analyses/HFaction analysis.R`, `Analyses/Readmission analysis.R`, and `Analyses/Sample size based on HF-action.R`.

We recommend running the master script section by section rather than all at once, especially because some steps are time-consuming. You can also run each script separately to reproduce a specific scenario or analysis.

By default, the master script generates datasets for all simulation scenarios, then runs the description and estimation scripts for scenario 1 only. The simulation-based WR sample-size calculation is skipped by default because `run_wr_simulation <- FALSE`.

To run the description and estimation for scenario 2, set `simulationDataPath <- "datasets scenario 2"` and update `wrTrue` to use the scenario 2 values: `1.2901` for the unstratified analysis and `1.2904` for the stratified analysis.

Precisely, the workflow consists of the following steps:

1) Generate datasets  
   - Open `Simulations/Generation.R` and run the script. It creates folders such as `datasets scenario 1`, `datasets scenario 2`, ... (500 datasets each) plus `*BIG` versions with `1e5` subjects for "asymptotic" win ratio values.  
   - Scenarios vary the frailty variance (`theta`), baseline rates (`lambda0`, `r0`), and baseline shape (exponential vs. log-logistic). Uncomment the heterogeneous generator (`dataJFM_fast_heterog.cpp`) for the bonus scenario if desired (and comment out `dataJFM_fast.cpp`).

2) Describe simulated datasets  
   - Run `Simulations/Description.R` to summarise mortality, recurrence counts, and censoring across all saved `.rds` files for a scenario.

3) Estimate performance  
   - Run `Simulations/Estimation.R` after generation completes. Set / adjust in the script:
     - Win ratio: dataset path and the pre-computed "true" (asymptotic) WR .
     - Joint frailty: dataset path, baseline hazard, or `alpha_power` to explore other scenarios.
   - Output is printed as markdown tables (bias, SE, coverage, power, number of ties/pairs).

## Applications
- `Analyses/HFaction analysis.R`: load the prepared HF-ACTION dataset, fits spline-based joint frailty models (unadjusted and age-adjusted), and computes last-event-assisted WR with and without stratification.
- `Analyses/Readmission analysis.R`: same workflow for the readmission dataset bundled with `frailtypack`, including multiple stratified WR variants.
- `Analyses/Sample size based on HF-action.R`: Schoenfeld-style event calculations, WR sample size via `WRSS()`, joint frailty sample size via `JFM.ssize()`, and a WR power curve simulated with `sz_lwr()`.

Note: if you, somehow, encounter a likelihood computation problem, please try to rerun the code.

## Data generator notes
- `Simulations/helpers/dataJFM_fast.cpp` simulates recurrent events and a terminal event under either exponential or log-logistic baselines and gamma frailty (mean 1, variance `theta`).
- `Simulations/helpers/genCpp.R` wraps the C++ generator and saves datasets with seeds that match the file names.

## How to cite
Please cite the companion paper. If you reuse or adapt the code, a link to this repository is appreciated.
