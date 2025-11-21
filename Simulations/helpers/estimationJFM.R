# ------------------------------------------------------------
# Libraries
# ------------------------------------------------------------
library(frailtypack)
library(knitr)
library(dplyr)
library(MASS)


# ------------------------------------------------------------
# Numeric utilities
# ------------------------------------------------------------
safe_mat_inv <- function(M) {
  inv_try <- try(chol2inv(chol(M)), silent = TRUE)
  if (inherits(inv_try, "try-error")) MASS::ginv(M) else inv_try
}

wald_global_from_fit <- function(fit, use_robust = TRUE) {
  Vfull <- if (use_robust && !is.null(fit$varHIH)) fit$varHIH else fit$varH
  idx_v <- c(2L, 4L) #beta1 & betastar positions in varH/varHIH
  b <- c(fit$coef[1], fit$coef[3])

  Vsub <- Vfull[idx_v, idx_v, drop = FALSE]
  Vinv <- safe_mat_inv(Vsub)

  W <- as.numeric(crossprod(b, Vinv %*% b))
  p <- pchisq(W, df = 2, lower.tail = FALSE)
  list(W = W, df = 2L, p = p, robust = use_robust && !is.null(fit$varHIH))
}

# ------------------------------------------------------------
# Model fitting
# ------------------------------------------------------------
fitModel <- function(
  longData,
  hazard = "Weibull",
  n.knots = 7,
  kappa = c(1e-3, 1e-3),
  use_robust_test = TRUE
) {
  fit <- NULL
  nb_gl_sequence <- c(50, 32, 20)

  base_args <- list(
    formula = Surv(tStart, tStop, recurrence) ~ z1 +
      z2 +
      cluster(id) +
      terminal(death),
    formula.terminalEvent = ~z1,
    recurrentAG = TRUE,
    hazard = hazard,
    data = longData,
    Alpha = "None",
    init.B = c(log(0.7), log(0.9), log(0.8)),
    init.Theta = 0.5,
    maxit = 100
  )

  if (hazard == "Weibull") {
    for (gl_val in nb_gl_sequence) {
      fit_args <- c(base_args, list(nb.gl = gl_val))
      fit_attempt <- try(do.call("frailtyPenal", fit_args), silent = TRUE)
      if (!inherits(fit_attempt, "try-error") && fit_attempt$istop == 1) {
        fit <- fit_attempt
        break
      }
    }
  } else if (hazard == "Splines") {
    current_kappa <- kappa
    while (all(current_kappa <= 1e12)) {
      for (gl_val in nb_gl_sequence) {
        fit_args <- c(
          base_args,
          list(n.knots = n.knots, kappa = current_kappa, nb.gl = gl_val)
        )
        fit_attempt <- try(do.call("frailtyPenal", fit_args), silent = TRUE)
        if (!inherits(fit_attempt, "try-error") && fit_attempt$istop == 1) {
          fit <- fit_attempt
          break
        }
      }
      if (!is.null(fit)) {
        break
      }
      current_kappa <- current_kappa * 10
    }
  }

  if (is.null(fit)) {
    return(NULL)
  }

  # --- SEs based on model-based varH
  seVec <- sqrt(diag(fit$varH))
  names(seVec) <- c("sqrtTheta", "beta1", "beta2", "betaStar")
  seTheta <- abs(2 * sqrt(fit$theta)) * seVec["sqrtTheta"] #delta method

  # --- 2-df global Wald test for H0: beta1 = betaStar = 0
  gtest <- wald_global_from_fit(fit, use_robust = use_robust_test)

  list(
    estimate = c(
      beta1 = fit$coef[1],
      beta2 = fit$coef[2],
      betaStar = fit$coef[3],
      theta = fit$theta
    ),
    se = c(
      beta1 = seVec["beta1"],
      beta2 = seVec["beta2"],
      betaStar = seVec["betaStar"],
      theta = seTheta
    ),
    global = gtest # includes W, df = 2, p-value and whether robust variance was used
  )
}

# ------------------------------------------------------------
# Monte-Carlo wrapper
# ------------------------------------------------------------
runMC <- function(
  seedMax = 500,
  seedStart = 1L,
  hazard = "Weibull",
  dataPath,
  alpha_power = 0.05,
  use_robust_test = TRUE,
  count_fail_as_no_reject = TRUE
) {
  results <- vector("list", seedMax)

  # True parameter values from first dataset
  first_dataset_path <- file.path(
    dataPath,
    paste0("dataset_seed_", seedStart, ".rds")
  )
  true <- unlist(readRDS(first_dataset_path)$truth)
  true <- true[c("beta1", "beta2", "betaStar", "theta")] # drop alpha

  for (r in seq_len(seedMax)) {
    current_seed <- seedStart + r - 1L
    dataset_path <- file.path(
      dataPath,
      paste0("dataset_seed_", current_seed, ".rds")
    )
    sim_data <- readRDS(dataset_path)
    fit <- tryCatch(
      fitModel(
        sim_data$longData,
        hazard = hazard,
        use_robust_test = use_robust_test
      ),
      error = function(e) NULL
    )
    results[[r]] <- fit
  }

  ok <- vapply(results, Negate(is.null), logical(1))
  n_success <- sum(ok)
  n_fail <- seedMax - n_success
  if (n_success == 0) {
    stop(
      "No successful model fits. Check dataset files or model specifications."
    )
  }

  # --- Estimation performance
  estMat <- do.call(rbind, lapply(results[ok], `[[`, "estimate"))
  seMat <- do.call(rbind, lapply(results[ok], `[[`, "se"))

  means <- colMeans(estMat)
  absBias <- means - true
  relBias <- absBias / true
  empSE <- apply(estMat, 2, sd)
  meanSE <- colMeans(seMat)

  trueMat <- matrix(true, nrow(estMat), length(true), byrow = TRUE)
  cover <- colMeans(
    ((estMat - qnorm(0.975) * seMat) <= trueMat) &
      ((estMat + qnorm(0.975) * seMat) >= trueMat),
    na.rm = TRUE
  )

  summary_df <- data.frame(
    `True Value` = true,
    `Mean estimate` = means,
    `Absolute bias` = absBias,
    `Relative bias` = relBias * 100,
    `Empirical SE` = empSE,
    `Asymptotic SE` = meanSE,
    `Coverage rate` = cover
  )
  row.names(summary_df) <- names(true)

  # --- Global Wald p-values and power
  pvals <- vapply(
    results,
    function(x) if (is.null(x)) NA_real_ else x$global$p,
    numeric(1)
  )

  if (count_fail_as_no_reject) {
    # Treat failed fits as "do not reject"
    reject <- pvals < alpha_power
    reject[is.na(reject)] <- FALSE
    B <- length(reject)
  } else {
    # Compute over successful fits only
    good <- !is.na(pvals)
    reject <- (pvals[good] < alpha_power)
    B <- length(reject)
  }

  power <- mean(reject)
  mcse <- sqrt(power * (1 - power) / B)
  ci95 <- binom.test(sum(reject), B, p = alpha_power)$conf.int

  list(
    summary = summary_df,
    nSuccess = n_success,
    nFail = n_fail,
    alpha = alpha_power,
    power = power,
    mcse = mcse,
    ci95 = ci95,
    pvals = pvals
  )
}
