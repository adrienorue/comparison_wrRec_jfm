library(WR)
library(frailtypack)
library(dplyr)

# Data ---------------------------------------------------------------
data(readmission)
df <- readmission

df <- df %>%
    mutate(
        status = case_when(
            death == 1 ~ 1L,
            event == 1 ~ 2L,
            event == 0 & death == 0 ~ 0L
        )
    )

df$chemo <- ifelse(df$chemo == "Treated", 1L, 0L)

df <- df %>%
    group_by(id) %>%
    mutate(
        # for stratified LWR
        strata = first(charlson),
        strata2 = last(charlson)
    ) %>%
    ungroup()
# View(df)

# Description ----------------------------------------------------------

# Summary of the number of recurrent events (min - max - mean - median per individual)
df %>%
    group_by(id) %>%
    summarise(n_events = sum(event == 1)) %>%
    summarise(
        min = min(n_events),
        max = max(n_events),
        mean = mean(n_events),
        median = median(n_events)
    )

#    min   max  mean median
#  <int> <int> <dbl>  <int>
#      0    22  1.14      1

# How many died
df %>%
    group_by(id) %>%
    summarise(death = max(death)) %>%
    summarise(n_death = sum(death == 1), n_alive = sum(death == 0))
#  n_death n_alive
#    <int>   <int>
#      109     294

# How many died after at least one rehospitalization
df %>%
    group_by(id) %>%
    summarise(n_events = sum(event == 1), death = max(death)) %>%
    filter(n_events > 0) %>%
    summarise(n_death = sum(death == 1), n_alive = sum(death == 0))
#   n_death n_alive
#     <int>   <int>
#        73     131

# Follow-up time (min - max - mean - median)
df %>%
    group_by(id) %>%
    summarise(lastTime = max(t.stop) / 365.25) %>%
    summarise(
        median = median(lastTime),
        mean = mean(lastTime),
        min = min(lastTime),
        max = max(lastTime)
    )
#  median  mean     min   max
#   <dbl> <dbl>   <dbl> <dbl>
#    3.09  2.81 0.00274  5.96

# Joint frailty models -------------------------------------------------
NB_GL <- 50

# # Unadjusted models
redSFM <- frailtyPenal(
    Surv(t.start, t.stop, event) ~ cluster(id) + as.factor(chemo),
    hazard = "Splines-per",
    nb.gl = NB_GL,
    cross.validation = TRUE,
    kappa = 1e4,
    n.knots = 6,
    recurrentAG = TRUE,
    data = df
)

redCox <- frailtyPenal(
    Surv(t.stop, death) ~ as.factor(chemo),
    hazard = "Splines-per",
    nb.gl = NB_GL,
    cross.validation = TRUE,
    kappa = 1e4,
    n.knots = 6,
    data = df[df$event == 0, ]
)

initBetas <- unname(c(redSFM$coef, redCox$coef))
initTheta <- redSFM$theta
initKappas <- c(redSFM$kappa, redCox$kappa)

fitJFM <- frailtyPenal(
    Surv(t.start, t.stop, event) ~ cluster(id) +
        as.factor(chemo) +
        terminal(death),
    formula.terminalEvent = ~ as.factor(chemo),
    hazard = "Splines-per",
    nb.gl = NB_GL,
    # Alpha = "None",
    kappa = initKappas,
    recurrentAG = TRUE,
    n.knots = 6, # quartiles
    init.B = initBetas,
    init.Alpha = 1,
    data = df
)

# Recurrences:
# -------------
#            coef exp(coef) SE coef (H) SE coef (HIH)        z       p
# chemo1 -0.25275  0.776662    0.163126      0.209654 -1.54941 0.12128

# Terminal event:
# ----------------
#            coef exp(coef) SE coef (H) SE coef (HIH)       z        p
# chemo1 0.477422   1.61191    0.271976      0.298052 1.75538 0.079195

#  Frailty parameters:
#    theta (variance of Frailties, w): 1.37657 (SE (H): 0.120022 ) p = < 1e-16
#    alpha (w^alpha for terminal event): 1.28281 (SE (H): 0.183086 ) p = 2.4414e-12

pnorm(-abs(-0.25275 / 0.209654)) * 2
pnorm(-abs(0.477422 / 0.298052)) * 2
exp(-0.25275 + c(qnorm(0.025), qnorm(0.975)) * 0.209654)
exp(0.477422 + c(qnorm(0.025), qnorm(0.975)) * 0.298052)


# # Adjusted models
redSFM_adjust <- frailtyPenal(
    Surv(t.start, t.stop, event) ~ cluster(id) +
        as.factor(chemo) +
        sex +
        charlson +
        dukes,
    hazard = "Splines-per",
    nb.gl = NB_GL,
    cross.validation = TRUE,
    recurrentAG = TRUE,
    kappa = 1e5,
    n.knots = 6,
    data = df
)

redCox_adjust <- frailtyPenal(
    Surv(t.stop, death) ~ as.factor(chemo) + sex + charlson + dukes,
    hazard = "Splines-per",
    nb.gl = NB_GL,
    cross.validation = TRUE,
    kappa = 1e5,
    n.knots = 6,
    data = df[df$event == 0, ]
)

initBetas <- unname(c(redSFM_adjust$coef, redCox_adjust$coef))
initTheta <- redSFM_adjust$theta
initKappas <- c(redSFM_adjust$kappa, redCox_adjust$kappa)

fitJFM_adjust <- frailtyPenal(
    Surv(t.start, t.stop, event) ~ cluster(id) +
        as.factor(chemo) +
        sex +
        charlson +
        dukes +
        terminal(death),
    formula.terminalEvent = ~ as.factor(chemo) + sex + charlson + dukes,
    recurrentAG = TRUE,
    hazard = "Splines-per",
    nb.gl = NB_GL,
    kappa = initKappas,
    n.knots = 6, # quartiles
    init.B = initBetas,
    init.Alpha = 1,
    data = df
)


# Recurrences:
# -------------
#                  coef exp(coef) SE coef (H) SE coef (HIH)        z          p
# chemo1      -0.153123  0.858025    0.164076      0.171210 -0.93324 3.5070e-01
# sexFemale   -0.635488  0.529677    0.156084      0.166447 -4.07146 4.6720e-05
# charlson1-2  0.487969  1.629005    0.284269      0.430219  1.71658 8.6056e-02
# charlson3    0.581851  1.789348    0.148098      0.188607  3.92882 8.5365e-05
# dukesC       0.358838  1.431664    0.182580      0.210085  1.96537 4.9371e-02
# dukesD       1.540986  4.669192    0.229083      0.256687  6.72677 1.7347e-11

#            chisq df global p
# charlson 18.3822  2 1.02e-04
# dukes    49.1121  2 2.16e-11

# Terminal event:
# ----------------
#                  coef exp(coef) SE coef (H) SE coef (HIH)         z          p
# chemo1       1.063744  2.897197    0.248877      0.296035  4.274179 1.9184e-05
# sexFemale   -0.248823  0.779718    0.224189      0.251862 -1.109879 2.6705e-01
# charlson1-2  0.482610  1.620298    0.625671      0.594918  0.771348 4.4050e-01
# charlson3    1.282153  3.604391    0.250719      0.298436  5.113895 3.1558e-07
# dukesC       1.279572  3.595102    0.357803      0.409117  3.576189 3.4864e-04
# dukesD       3.086509 21.900482    0.403609      0.450012  7.647282 2.0539e-14

#            chisq df global p
# charlson 26.1529  2 2.09e-06
# dukes    67.7511  2 1.89e-15

#  Frailty parameters:
#    theta (variance of Frailties, w): 1.07163 (SE (H): 0.114609 ) p = < 1e-16
#    alpha (w^alpha for terminal event): 0.632122 (SE (H): 0.15691 ) p = 5.6117e-05

# Recurrences
pnorm(-abs(-0.153123 / 0.171210)) * 2
pnorm(-abs(-0.635488 / 0.166447)) * 2
pnorm(-abs(0.487969 / 0.430219)) * 2
pnorm(-abs(0.581851 / 0.188607)) * 2
pnorm(-abs(0.358838 / 0.210085)) * 2
pnorm(-abs(1.540986 / 0.256687)) * 2
exp(-0.153123 + c(qnorm(0.025), qnorm(0.975)) * 0.171210)
exp(-0.635488 + c(qnorm(0.025), qnorm(0.975)) * 0.166447)
exp(0.487969 + c(qnorm(0.025), qnorm(0.975)) * 0.430219)
exp(0.581851 + c(qnorm(0.025), qnorm(0.975)) * 0.188607)
exp(0.358838 + c(qnorm(0.025), qnorm(0.975)) * 0.210085)
exp(1.540986 + c(qnorm(0.025), qnorm(0.975)) * 0.256687)

# Terminal
pnorm(-abs(1.063744 / 0.296035)) * 2
pnorm(-abs(-0.248823 / 0.251862)) * 2
pnorm(-abs(0.482610 / 0.594918)) * 2
pnorm(-abs(1.282153 / 0.298436)) * 2
pnorm(-abs(1.279572 / 0.409117)) * 2
pnorm(-abs(3.086509 / 0.450012)) * 2
exp(1.063744 + c(qnorm(0.025), qnorm(0.975)) * 0.296035)
exp(-0.248823 + c(qnorm(0.025), qnorm(0.975)) * 0.251862)
exp(0.482610 + c(qnorm(0.025), qnorm(0.975)) * 0.594918)
exp(1.282153 + c(qnorm(0.025), qnorm(0.975)) * 0.298436)
exp(1.279572 + c(qnorm(0.025), qnorm(0.975)) * 0.409117)
exp(3.086509 + c(qnorm(0.025), qnorm(0.975)) * 0.450012)


# Win ratio -------------------------------------------------

WRrec(
    ID = as.numeric(df$id),
    time = as.numeric(df$t.stop),
    status = as.numeric(df$status),
    trt = as.numeric(df$chemo)
)
#             N Rec. Event Death Med. Follow-up
# Control   186        282    51         1258.5
# Treatment 217        176    58         1054.0

# Analysis of last-event-assisted WR (LWR; recommended), first-event-assisted WR (FWR), and naive WR (NWR):
#     Win prob Loss prob WR (95% CI)*        p-value
# LWR 38.5%    39.3%     0.98 (0.754, 1.27)  0.878
# -----
# Total number of pairs:  40362

WRrec(
    ID = as.numeric(df$id),
    time = as.numeric(df$t.stop),
    status = as.numeric(df$status),
    trt = as.numeric(df$chemo),
    strata = as.numeric(df$sex)
)
#             N Rec. Event Death Med. Follow-up
# Control   186        282    51         1258.5
# Treatment 217        176    58         1054.0

# Analysis of last-event-assisted WR (LWR):
#     Win prob Loss prob WR (95% CI)*        p-value
# LWR 38.8%    39.4%     0.983 (0.756, 1.28) 0.897
# -----
# Total number of pairs:  20694

WRrec(
    ID = as.numeric(df$id),
    time = as.numeric(df$t.stop),
    status = as.numeric(df$status),
    trt = as.numeric(df$chemo),
    strata = as.numeric(df$strata) # -> first(charlson)
)
#             N Rec. Event Death Med. Follow-up
# Control   186        282    51         1258.5
# Treatment 217        176    58         1054.0

# Analysis of last-event-assisted WR (LWR):
#     Win prob Loss prob WR (95% CI)*        p-value
# LWR 34.3%    43.4%     0.789 (0.607, 1.03) 0.0761
# -----
# Total number of pairs:  22023

WRrec(
    ID = as.numeric(df$id),
    time = as.numeric(df$t.stop),
    status = as.numeric(df$status),
    trt = as.numeric(df$chemo),
    strata = as.numeric(df$dukes)
)
#             N Rec. Event Death Med. Follow-up
# Control   186        282    51         1258.5
# Treatment 217        176    58         1054.0

# Analysis of last-event-assisted WR (LWR):
#     Win prob Loss prob WR (95% CI)*        p-value
# LWR 32.6%    43.1%     0.755 (0.565, 1.01) 0.0579
# -----
# Total number of pairs:  12519

# # dukes x charlson
df$stratificationVar1 <- as.numeric(interaction(df$dukes, df$strata))
WRrec(
    ID = as.numeric(df$id),
    time = as.numeric(df$t.stop),
    status = as.numeric(df$status),
    trt = as.numeric(df$chemo),
    strata = as.numeric(df$stratificationVar1)
)
#             N Rec. Event Death Med. Follow-up
# Control   186        282    51         1258.5
# Treatment 217        176    58         1054.0

# Analysis of last-event-assisted WR (LWR):
#     Win prob Loss prob WR (95% CI)*        p-value
# LWR 33.2%    42.4%     0.782 (0.582, 1.05) 0.104
# -----
# Total number of pairs:  7975

# Number of pairs
df_unique <- df %>%
    group_by(id) %>%
    slice(1) %>%
    ungroup()
t0 <- table(df_unique$chemo)
sum(t0[1] * t0[2])

t1 <- table(df_unique$sex, df_unique$chemo)
sum(t1[, 1] * t1[, 2])

t2 <- table(df_unique$strata, df_unique$chemo)
sum(t2[, 1] * t2[, 2])

t3 <- table(df_unique$dukes, df_unique$chemo)
sum(t3[, 1] * t3[, 2])

t4 <- table(df_unique$strata, df_unique$dukes, df_unique$chemo)
sum(t4[,, 1] * t4[,, 2])
