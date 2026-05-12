library(here)

# ------------------------------------------------------------
# ========== Simulations

# generates scenarios' datasets ; may be long
source(here("Simulations", "Generation.R"), print.eval = TRUE)


# For the "wrTrue" vector (lines 28 to 31), you need the real WR values for the "BIG" datasets.
# - Values are presented in the paper, Table 4, and in the code below (lines 34-35)
# - If you want to reproduce them, please uncomment the code below and run it. 
#   Caution: it may take a while to run.

# library(WR) # please install our optimized version of WR for this exact purpose
# source(here("Simulations", "helpers", "estimationWR.R"))
# trueWR_mc_run <- runMC_WR(
#     seedMax = 50,
#     seedStart = 1,
#     dataPath = here("Simulations", "datasets scenario 1BIG"), # adapt the path
#     id = "id",
#     time = "time",
#     status = "status",
#     trt = "z1",
#     # strata = "z2", # uncomment for stratified ; comment for unstratified
#     trueWR = 1
# )


# change the dataset path according to folders created by Generation.R,
# that is "datasets scenario X", for X = 1,2,3,4,5 or 6
# Values for each scenario are presented in the paper, Table 4. For ease of use,
# here are the values for scenario 1,2,3,4, 5 and 6 respectively:
# - unstratified : 1.2253, 1.2901, 1.1912, 1.1962, 1.3465, 1.2514
# - stratified   : 1.2254, 1.2904, 1.1912, 1.1962, 1.3470, 1.2515
wrTrue <- c(
    unstratified = 1.2253,
    stratified = 1.2254
)
simulationDataPath <- "datasets scenario 1"
source(here("Simulations", "Description.R"), print.eval = TRUE)
source(here("Simulations", "Estimation.R"), print.eval = TRUE)



# ------------------------------------------------------------
# ========== Applications

source(here("Analyses", "Readmission analysis.R"), print.eval = TRUE)
source(here("Analyses", "HFaction analysis.R"), print.eval = TRUE)

#if TRUE, perform sample-size calculation using the simulation-based approach for the win ratio.
run_wr_simulation <- FALSE # if TRUE, will be long to run
source(here("Analyses", "Sample size based on HF-action.R"), print.eval = TRUE)
