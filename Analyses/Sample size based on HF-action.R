script_path <- rstudioapi::getSourceEditorContext()$path
setwd(dirname(script_path))
source("timing_helpers.R")
rm(script_path)

# "Schoenfeld" ----------------------------------------------------------------
nbEvent <- function(alpha, power, p, HR) {
  z_alpha <- qnorm(1 - alpha / 2)
  z_beta <- qnorm(power)
  d <- (z_alpha + z_beta)^2 / (log(HR)^2 * p * (1 - p))
  return(ceiling(d))
}

lambdaC <- -log(1 - 0.3)
lambdaT <- -1 / 2 * log(1 - 0.89 * (1 - exp(-2 * lambdaC)))
HRest <- lambdaT / lambdaC
nbEvent80P <- nbEvent(alpha = 0.05, power = 0.8, p = 0.5, HR = HRest)
nbEvent90P <- nbEvent(alpha = 0.05, power = 0.9, p = 0.5, HR = HRest)

probaEventTreat <- 1 - (exp(-lambdaT) - exp(-4 * lambdaT)) / (3 * lambdaT)
probaEventControl <- 1 - (exp(-lambdaC) - exp(-4 * lambdaC)) / (3 * lambdaC)
meanProbaEvent <- (probaEventTreat + probaEventControl) / 2

nbEvent80P / meanProbaEvent # round => 2132
nbEvent90P / meanProbaEvent # round + nearest odd => 2856


# Standard win ratio sample size ----------------------------------------------
library(WR)

wr_sample_size <- time_expr({
  dat <- read.csv("hfaction.csv", sep = ",")
  pilot <- dat[dat$trt_ab == 0, ] # control group
  id <- pilot$patid
  time <- pilot$time / 12 # months -> years
  status <- pilot$status
  # # Baseline parameters for the Gumbel-Hougaard copula
  gum <- gumbel.est(id, time, status)
  lambda_D <- gum$lambda_D
  lambda_H <- gum$lambda_H
  kappa <- gum$kappa
  tau <- 4 # 4 years of follow-up max
  tau_b <- 3 # 3 years accrual
  lambda_L <- 0.05 # loss to follow-up rate
  bparam <- base(lambda_D, lambda_H, kappa, tau_b, tau, lambda_L)

  n80 <- WRSS(
    xi = log(c(0.9, 0.8)),
    bparam = bparam,
    q = 0.5,
    alpha = 0.05,
    power = 0.8
  )$n

  n90 <- WRSS(
    xi = log(c(0.9, 0.8)),
    bparam = bparam,
    q = 0.5,
    alpha = 0.05,
    power = 0.9
  )$n

  list(n80 = n80, n90 = n90)
})

wr_sample_size$value$n80
wr_sample_size$value$n90

# Joint frailty model sample size ---------------------------------------------
library(frailtypack)

# # Assumptions
# # # - Accrual duration: 3 years
# # # - End of study: 4 years (3 years of accrual + 1 additional minimum year of follow-up)
# # # - Allocation ratio experimental:control of 1:1
# # # - Mean number of 3 recurrent events per patient in the control group (~Poisson distribution)
# # # - Median time to first hospitalization in the control group: 3.6 months
# # # - Median time to death in the control group: 28 months
# # # - The hazard increases over time (shapes > 1),
# # #   faster for death (shape = 2) than for hospitalization (shape = 1.5)
# # # - Hazard ratio for first hospitalization: 0.8
# # # - Hazard ratio for death: 0.9
# # # - Between-event correlation assumed (frailty variance = 1)
# # # - Recurrent events are positively associated with the terminal event (alpha = 1)
# # # - Bilateral joint test, with a two-sided alpha = 0.05

jfm_sample_size <- time_expr({
  res80 <- JFM.ssize(
    power = 0.8,
    ni = 3,
    ni.type = "Pois",
    Acc.Dur = 3,
    FUP = 4,
    FUP.type = "UpToEnd",
    medianR.H0 = 3.6 / 12,
    medianD.H0 = 28 / 12,
    betaR.H0 = log(1),
    betaR.HA = log(0.8),
    betaD.H0 = log(1),
    betaD.HA = log(0.9),
    shapeR.W = 1,
    shapeD.W = 1,
    theta = 1,
    alpha = 1,
    ratio = 1,
    samples.mc = 1e5,
    seed = 42,
    timescale = "calendar",
    betaTest.type = "joint",
    statistic = "Wald",
    typeIerror = 0.05,
    test.type = "2-sided"
  )

  res90 <- JFM.ssize(
    power = 0.9,
    ni = 3,
    ni.type = "Pois",
    Acc.Dur = 3,
    FUP = 4,
    FUP.type = "UpToEnd",
    medianR.H0 = 3.6 / 12,
    medianD.H0 = 28 / 12,
    betaR.H0 = log(1),
    betaR.HA = log(0.8),
    betaD.H0 = log(1),
    betaD.HA = log(0.9),
    shapeR.W = 1,
    shapeD.W = 1,
    theta = 1,
    alpha = 1,
    ratio = 1,
    samples.mc = 1e5,
    seed = 42,
    timescale = "calendar",
    betaTest.type = "joint",
    statistic = "Wald",
    typeIerror = 0.05,
    test.type = "2-sided"
  )

  list(res80 = res80, res90 = res90)
})

jfm_sample_size$value$res80
jfm_sample_size$value$res90

print_timing("Sample size - Win ratio (80% + 90%)", wr_sample_size$elapsed)
print_timing("Sample size - Joint frailty model (80% + 90%)", jfm_sample_size$elapsed)


# Win ratio sample size through simulations -----------------------------------

# # Same assumptions as above
library(WR)

run_wr_simulation <- FALSE

if (run_wr_simulation) {
  res <- sz_lwr(
    power = 0.80,
    type1Error = 0.05,
    weibShapeRec = 1,
    weibScaleRec = (3.6 / 12) / log(2),
    weibShapeTerm = 1,
    weibScaleTerm = (28 / 12) / log(2),
    HR_recurrent = 0.70,
    HR_terminal = 0.80,
    niType = "poisson",
    ni = 3,
    theta = 0.5,
    alpha = 1.0,
    fupType = "uptoend",
    FUP = 4,
    accrualDuration = 3,
    sample_size_min = 1000,
    sample_size_max = 1400,
    sample_size_step = 100,
    n_sim = 2000,
    baseSeed = 42,
    output_plot_file = sprintf(
      "power_curve_%s.png",
      gsub("[ :\\-]", "_", round(Sys.time(), 0))
    )
  )
}
