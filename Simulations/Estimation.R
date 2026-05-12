# ------------------------------------------------------------
# Set current script location as working directory
# ------------------------------------------------------------
library(here)
setwd(here("Simulations"))


# ------------------------------------------------------------
# Win ratio simulations
# ------------------------------------------------------------
library(WR)
source("helpers/estimationWR.R")

mc_wr_unstratified <- runMC_WR(
  seedMax = 500,
  seedStart = 1,
  dataPath = simulationDataPath,
  id = "id",
  time = "time",
  status = "status",
  trt = "z1",
  trueWR = unname(wrTrue["unstratified"])
)

mc_wr_stratified <- runMC_WR(
  seedMax = 500,
  seedStart = 1,
  dataPath = simulationDataPath,
  id = "id",
  time = "time",
  status = "status",
  trt = "z1",
  strata = "z2",
  trueWR = unname(wrTrue["stratified"])
)

wr_summary <- rbind(
  cbind(Analysis = "Unstratified", mc_wr_unstratified$a),
  cbind(Analysis = "Stratified", mc_wr_stratified$a)
)
row.names(wr_summary) <- NULL

wr_ties <- rbind(
  cbind(Analysis = "Unstratified", mc_wr_unstratified$b),
  cbind(Analysis = "Stratified", mc_wr_stratified$b)
)
row.names(wr_ties) <- NULL


md1 <- knitr::kable(
  wr_summary,
  digits = 5,
  format = "markdown",
  caption = "Win ratio Monte-Carlo simulation",
  row.names = FALSE
)

md2 <- knitr::kable(
  wr_ties,
  digits = 4,
  format = "markdown",
  caption = "Ties summary",
  row.names = FALSE
)

cat(md1, sep = "\n")
cat("\n\n")
cat(md2, sep = "\n")


# ------------------------------------------------------------
# Joint frailty model simulations
# ------------------------------------------------------------
source("helpers/estimationJFM.R")

mc <- runMC(
  seedMax = 500L,
  seedStart = 1L,
  hazard = "Weibull",
  dataPath = simulationDataPath,
  alpha_power = 0.05,
  use_robust_test = TRUE,
  count_fail_as_no_reject = TRUE
)

sprintf(
  "Successful fits: %d out of %d\n\nGlobal Wald (2-df) power at alpha = %.3f:\n  Power = %.3f | MCSE = %.3f | 95%% MC CI = (%.3f, %.3f)",
  mc$nSuccess,
  mc$nSuccess + mc$nFail,
  mc$alpha,
  mc$power,
  mc$mcse,
  mc$ci95[1],
  mc$ci95[2]
)

kable(
  mc$summary,
  digits = 4,
  caption = "Monte-Carlo performance of the joint frailty model"
)
