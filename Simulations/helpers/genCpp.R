dataJFM <- function(
  nSubjects = 400,
  theta = 0.5,
  alpha = 1,
  fixedCensor = 3.0,

  # exponential baselines rates
  r0 = 1.5, # recurrent
  lambda0 = 0.5, # death
  betaGap = c(log(0.7), log(0.9)),
  betaDeathStar = log(0.8),

  effect = "ph", # ph=multiplicative, ah=additive
  baseline = "exponential",

  # only used if baseline == "loglogistic"
  shapeRec = 1.3,
  scaleRec = 2.0,
  shapeDeath = 1.2,
  scaleDeath = 3.0,
  seed = 42
) {
  effect <- match.arg(effect, c("ph", "ah"))
  baseline <- match.arg(baseline, c("exponential", "loglogistic"))

  out <- dataJFM_fast(
    nSubjects,
    theta,
    alpha,
    fixedCensor,
    lambda0,
    r0,
    as.numeric(betaGap),
    betaDeathStar,
    if (is.null(seed)) NA_integer_ else as.integer(seed),
    baseline_type = switch(baseline, exponential = 0L, loglogistic = 1L),
    effect_type = switch(effect, ph = 0L, ah = 1L),
    shapeRec,
    scaleRec,
    shapeDeath,
    scaleDeath
  )

  longData <- as.data.frame(out$longData)
  longData <- longData[order(longData$id, longData$time), ]

  # compute tStart, tStop, and status in base R
  longData <- within(longData, {
    tStart <- ave(time, id, FUN = function(x) c(0, head(x, -1)))
    tStop <- time
    status <- ifelse(
      death == 0 & recurrence == 0,
      0L,
      ifelse(
        death == 1 & recurrence == 0,
        1L,
        ifelse(death == 0 & recurrence == 1, 2L, NA_integer_)
      )
    )
  })

  list(longData = longData, subjectData = out$subjectData, truth = out$truth)
}

genDatas <- function(nbDatasets, outputDir, ...) {
  if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = TRUE)
  }

  for (i in 1:nbDatasets) {
    dataSim <- dataJFM(seed = i, ...)
    saveRDS(dataSim, file = paste0(outputDir, "dataset_seed_", i, ".rds"))
  }
}
