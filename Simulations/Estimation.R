# ------------------------------------------------------------
# Set current script location as working directory
# ------------------------------------------------------------
library(rstudioapi)
a <- getActiveDocumentContext()
setwd(dirname(a$path))
remove(a)


# ------------------------------------------------------------
# Win ratio simulations
# ------------------------------------------------------------
library(WR)
source("helpers/estimationWR.R")

mc_wr <- runMC_WR(
  seedMax = 500,
  seedStart = 1,
  dataPath = "datasets scenario 1",
  id = "id",
  time = "time",
  status = "status",
  trt = "z1",
  # strata = "z2", # uncomment to use stratified LWR
  trueWR = 1.2515 # use the computed value from the 50 big scenario-specific datasets
)


md1 <- knitr::kable(
  mc_wr$a,
  digits = 5,
  format = "markdown",
  caption = "Win ratio Monte-Carlo simulation",
  row.names = FALSE
)

md2 <- knitr::kable(
  mc_wr$b,
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
  dataPath = "datasets scenario 2",
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
