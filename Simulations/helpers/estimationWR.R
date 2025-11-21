unwrap_dataset <- function(obj) {
    if (is.data.frame(obj)) {
        return(obj)
    }

    if (is.list(obj)) {
        if (!is.null(obj$longData) && is.data.frame(obj$longData)) {
            return(obj$longData)
        }
    }
    stop("Cannot find a data.frame ($longData) inside the .rds object.")
}

fitWR <- function(
    dat,
    id = "id",
    time = "time",
    status = "status",
    trt = "z1",
    strata = NULL,
    ...
) {
    if (!is.data.frame(dat)) {
        stop("`dat` must be a data.frame-like object.")
    }

    obj <- WRrec(
        ID = dat[[id]],
        time = dat[[time]],
        status = dat[[status]],
        trt = dat[[trt]],
        strata = if (!is.null(strata)) dat[[strata]] else NULL,
        ...
    )

    logWR <- unname(obj$log.WR)
    se_log <- unname(obj$se)
    WR_est <- exp(logWR)
    se_WR <- WR_est * se_log
    theta <- unname(obj$theta)

    trt_vec <- dat[[trt]]
    id_vec <- dat[[id]]

    if (is.null(strata)) {
        n1 <- length(unique(id_vec[trt_vec == 1]))
        n0 <- length(unique(id_vec[trt_vec == 0]))
        Pairs <- as.numeric(n1 * n0)
    } else {
        strata_vec <- dat[[strata]]
        strata_levels <- unique(strata_vec)
        Pairs <- 0
        for (lev in strata_levels) {
            n1_lev <- length(unique(id_vec[trt_vec == 1 & strata_vec == lev]))
            n0_lev <- length(unique(id_vec[trt_vec == 0 & strata_vec == lev]))
            Pairs <- as.numeric(Pairs + n1_lev * n0_lev)
        }
    }

    tie_prop <- 1 - sum(theta)
    ties_cnt <- round(Pairs * max(0, tie_prop))

    c(
        logWR = logWR,
        WR = WR_est,
        invWR = 1 / WR_est,
        se_log = se_log,
        se_WR = se_WR,
        pval = unname(obj$pval),
        pairs = Pairs,
        ties = ties_cnt,
        ties_prop = tie_prop
    )
}

runMC_WR <- function(
    seedMax = 100L,
    seedStart = 1L,
    dataPath = ".",
    id = "id",
    time = "time",
    status = "status",
    trt = "z1",
    strata = NULL,
    filePattern = "dataset_seed_%d.rds",
    trueWR,
    verbose = FALSE,
    ...
) {
    R <- seedMax
    WR_vec <- rep(NA_real_, R)
    inv_vec <- rep(NA_real_, R)
    se_vec <- rep(NA_real_, R)
    cov_vec <- rep(NA, R)
    pairs <- rep(NA_real_, R)
    ties <- rep(NA_real_, R)
    ties_p <- rep(NA_real_, R)
    signif <- rep(NA_real_, R)

    for (i in seq_len(R)) {
        seed_i <- seedStart + i - 1L
        file_i <- file.path(dataPath, sprintf(filePattern, seed_i))

        if (!file.exists(file_i)) {
            warning("File not found: ", file_i)
            next
        }

        raw_i <- readRDS(file_i)
        dat_i <- tryCatch(unwrap_dataset(raw_i), error = function(e) {
            if (isTRUE(verbose)) {
                message(sprintf(
                    "Seed %d – unwrap failed: %s",
                    seed_i,
                    e$message
                ))
            }
            NULL
        })
        if (is.null(dat_i)) {
            next
        }

        res <- tryCatch(
            fitWR(
                dat_i,
                id = id,
                time = time,
                status = status,
                trt = trt,
                strata = strata,
                ...
            ),
            error = function(e) {
                if (isTRUE(verbose)) {
                    message(sprintf(
                        "Seed %d – WRrec failed: %s",
                        seed_i,
                        e$message
                    ))
                }
                NULL
            }
        )

        if (is.null(res)) {
            next
        }

        WR_vec[i] <- res["WR"]
        inv_vec[i] <- res["invWR"]
        se_vec[i] <- res["se_WR"]

        ci_low <- exp(res["logWR"] - abs(qnorm(0.025)) * res["se_log"])
        ci_high <- exp(res["logWR"] + abs(qnorm(0.025)) * res["se_log"])
        cov_vec[i] <- (trueWR >= ci_low) & (trueWR <= ci_high)

        signif[i] <- (res["pval"] < 0.05)

        pairs[i] <- res["pairs"]
        ties[i] <- res["ties"]
        ties_p[i] <- res["ties_prop"]
    }

    ok <- !is.na(WR_vec)
    if (!any(ok)) {
        stop("No successful replicates were obtained.")
    }

    ties_df <- data.frame(
        `Pairs mean` = as.numeric(mean(pairs[ok], na.rm = TRUE)),
        `Pairs min` = as.numeric(min(pairs[ok], na.rm = TRUE)),
        `Pairs max` = as.numeric(max(pairs[ok], na.rm = TRUE)),
        `Ties mean` = as.numeric(mean(ties[ok], na.rm = TRUE)),
        `Ties min` = as.numeric(min(ties[ok], na.rm = TRUE)),
        `Ties max` = as.numeric(max(ties[ok], na.rm = TRUE)),
        `Tie proportion mean (%)` = mean(ties_p[ok], na.rm = TRUE) * 100,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    summary_df <- data.frame(
        `True WR` = trueWR,
        `Mean WR` = mean(WR_vec[ok]),
        `Abs. bias` = mean(WR_vec[ok] - trueWR),
        `Rel. bias (%)` = mean((WR_vec[ok] - trueWR) / trueWR) * 100,
        `Mean 1/WR` = mean(inv_vec[ok]),
        `Empirical SE` = sd(WR_vec[ok]),
        `Asymptotic SE` = mean(se_vec[ok]),
        `Coverage rate (%)` = mean(cov_vec[ok]) * 100,
        `Power (%)` = mean(signif[ok]) * 100,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    attr(summary_df, "nSuccess") <- sum(ok)
    attr(summary_df, "nFail") <- R - sum(ok)

    return(list(a = summary_df, b = ties_df))
}
