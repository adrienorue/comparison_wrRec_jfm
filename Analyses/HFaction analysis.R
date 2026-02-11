script_path <- rstudioapi::getSourceEditorContext()$path
setwd(dirname(script_path))
source("timing_helpers.R")
rm(script_path)
timing <- list(jfm = 0, wr = 0)

library(WR)
library(frailtypack)
library(dplyr)

# Data ---------------------------------------------------------------
# df <- hfaction_cpx9

# df <- df %>%
#   group_by(patid) %>%
#   mutate(
#     # mainly for frailtypack::frailtyPenal()
#     t.start = lag(time, default = 0),
#     t.stop = time,
#     gap = t.stop - t.start,
#     event = as.numeric(status == 2),
#     death = as.numeric(status == 1),
#     trt_ab = factor(trt_ab)
#   ) %>%
#   ungroup()

# df$trt_ab <- as.numeric(df$trt_ab) - 1
# # Manage potentially problematics observations
# # # df[df$patid %in% c("HFACT00662", "HFACT01359"), ]
# # # - HFACT00662 has an event at the same time as the censoring time
# # # - HFACT01359 has an event at time = 0
# df <- df %>%
#   group_by(patid) %>%
#   mutate(
#     time = ifelse(gap == 0, time + 1e-3, time),
#     t.stop = time,
#     gap = t.stop - t.start
#   ) %>%
#   ungroup()

# write.csv(df, "hfaction.csv", row.names = FALSE)

df <- read.csv("hfaction.csv", sep = ",")


# Description ----------------------------------------------------------

# Summary of the number of recurrent events
df %>%
  group_by(patid) %>%
  summarise(n_events = sum(event == 1)) %>%
  summarise(
    min = min(n_events),
    max = max(n_events),
    mean = mean(n_events),
    median = median(n_events)
  )
#   min   max  mean median
# <int> <int> <dbl>  <dbl>
#     0    26  2.40      1

# How many died ?
df %>%
  group_by(patid) %>%
  summarise(death = max(death)) %>%
  summarise(n_death = sum(death == 1), n_alive = sum(death == 0))
#  n_death n_alive
#    <int>   <int>
#      93     333

# How many died after at least one rehospitalization
df %>%
  group_by(patid) %>%
  summarise(n_events = sum(event == 1), death = max(death)) %>%
  filter(n_events > 0) %>%
  summarise(n_death = sum(death == 1), n_alive = sum(death == 0))
#   n_death n_alive
#     <int>   <int>
#        82     233

# Follow-up time (min - max - mean - median)
df %>%
  group_by(patid) %>%
  summarise(lastTime = max(t.stop)) %>%
  summarise(
    median = median(lastTime),
    mean = mean(lastTime),
    min = min(lastTime),
    max = max(lastTime)
  )
# median  mean   min   max
#  <dbl> <dbl> <dbl> <dbl>
#   28.0  28.6 0.328  52.8

# Joint frailty models -------------------------------------------------
NB_GL <- 50

# # Unadjusted models
tmp <- time_expr(frailtyPenal(
  Surv(t.start, t.stop, event) ~ cluster(patid) + trt_ab,
  hazard = "Splines-per",
  nb.gl = NB_GL,
  cross.validation = TRUE,
  kappa = 1e4,
  n.knots = 6,
  recurrentAG = TRUE,
  data = df
))
redSFM <- tmp$value
timing$jfm <- timing$jfm + tmp$elapsed

tmp <- time_expr(frailtyPenal(
  Surv(t.stop, death) ~ trt_ab,
  hazard = "Splines-per",
  nb.gl = NB_GL,
  cross.validation = TRUE,
  kappa = 1e4,
  n.knots = 6,
  data = df[df$event == 0, ]
))
redCox <- tmp$value
timing$jfm <- timing$jfm + tmp$elapsed

initBetas <- unname(c(redSFM$coef, redCox$coef))
initTheta <- redSFM$theta
initKappas <- c(redSFM$kappa, redCox$kappa)

tmp <- time_expr(frailtyPenal(
  Surv(t.start, t.stop, event) ~ cluster(patid) + trt_ab + terminal(death),
  formula.terminalEvent = ~trt_ab,
  hazard = "Splines-per",
  nb.gl = NB_GL,
  kappa = initKappas,
  recurrentAG = TRUE,
  n.knots = 6,
  init.B = initBetas,
  init.Theta = initTheta,
  init.Alpha = 1,
  data = df
))
fitJFM <- tmp$value
timing$jfm <- timing$jfm + tmp$elapsed

# Recurrences:
# -------------
#             coef exp(coef) SE coef (H) SE coef (HIH)        z        p
# trt_ab -0.253996  0.775695    0.120925      0.134265 -1.89175 0.058524

# Terminal event:
# ----------------
#             coef exp(coef) SE coef (H) SE coef (HIH)        z       p
# trt_ab -0.510773  0.600031    0.258358      0.262681 -1.94446 0.05184

#  Frailty parameters:
#    theta (variance of Frailties, w): 0.998177 (SE (HIH): 0.0802848 ) p = < 1e-16
#    alpha (w^alpha for terminal event): 1.36958 (SE (HIH): 0.206061 ) p = 3.0016e-11

# p-values and 95% CI using se(HIH)
pnorm(-abs(-0.253996 / 0.134265)) * 2
pnorm(-abs(-0.510774 / 0.262681)) * 2

exp(-0.253996 + c(qnorm(0.025), qnorm(0.975)) * 0.134265)
exp(-0.510774 + c(qnorm(0.025), qnorm(0.975)) * 0.262681)


# # Adjusted models
tmp <- time_expr(frailtyPenal(
  Surv(t.start, t.stop, event) ~ cluster(patid) + trt_ab + age60,
  hazard = "Splines-per",
  nb.gl = NB_GL,
  cross.validation = TRUE,
  kappa = 1e5,
  n.knots = 6,
  recurrentAG = TRUE,
  data = df
))
redSFM_adjust <- tmp$value
timing$jfm <- timing$jfm + tmp$elapsed

tmp <- time_expr(frailtyPenal(
  Surv(t.stop, death) ~ trt_ab + age60,
  hazard = "Splines-per",
  nb.gl = NB_GL,
  cross.validation = TRUE,
  kappa = 1e5,
  n.knots = 6,
  data = df[df$event == 0, ]
))
redCox_adjust <- tmp$value
timing$jfm <- timing$jfm + tmp$elapsed

initBetas_adjust <- unname(c(redSFM_adjust$coef, redCox_adjust$coef))
initTheta_adjust <- redSFM_adjust$theta
initKappas_adjust <- c(redSFM_adjust$kappa, redCox_adjust$kappa)

tmp <- time_expr(frailtyPenal(
  Surv(t.start, t.stop, event) ~ cluster(patid) +
    trt_ab +
    age60 +
    terminal(death),
  formula.terminalEvent = ~ trt_ab + age60,
  hazard = "Splines-per",
  nb.gl = NB_GL,
  kappa = initKappas_adjust,
  recurrentAG = TRUE,
  n.knots = 6,
  init.B = initBetas_adjust,
  init.Theta = initTheta_adjust,
  init.Alpha = 1.4,
  data = df
))
fitJFM_adjust <- tmp$value
timing$jfm <- timing$jfm + tmp$elapsed


# Recurrences:
# -------------
#             coef exp(coef) SE coef (H) SE coef (HIH)        z        p
# trt_ab -0.282816  0.753659    0.120818      0.129896 -2.17724 0.029462
# age60  -0.327295  0.720871    0.124361      0.128447 -2.54810 0.010831

# Terminal event:
# ----------------
#             coef exp(coef) SE coef (H) SE coef (HIH)        z        p
# trt_ab -0.520146  0.594434    0.271562      0.272083 -1.91172 0.055912
# age60   0.386695  1.472108    0.267207      0.274293  1.40979 0.158600

#  Frailty parameters:
#    theta (variance of Frailties, w): 0.974917 (SE (HIH): 0.0803426 ) p = < 1e-16
#    alpha (w^alpha for terminal event): 1.57187 (SE (HIH): 0.233596 ) p = 1.7081e-11

pnorm(-abs(-0.282828 / 0.129896)) * 2
pnorm(-abs(-0.327306 / 0.128443)) * 2
pnorm(-abs(-0.520181 / 0.272079)) * 2
pnorm(-abs(0.386671 / 0.274251)) * 2

exp(-0.282828 + c(qnorm(0.025), qnorm(0.975)) * 0.129896)
exp(-0.327306 + c(qnorm(0.025), qnorm(0.975)) * 0.128443)
exp(-0.520181 + c(qnorm(0.025), qnorm(0.975)) * 0.272079)
exp(0.386671 + c(qnorm(0.025), qnorm(0.975)) * 0.274251)


# Win ratio -------------------------------------------------

tmp <- time_expr(WRrec(
  ID = df$patid,
  time = df$time,
  status = df$status,
  trt = df$trt_ab,
  naive = TRUE
))
timing$wr <- timing$wr + tmp$elapsed

length(unique(df[df$trt_ab == 1, "patid"])) *
  length(unique(df[df$trt_ab == 0, "patid"]))

# N Rec. Event Death Med. Follow-up
# Control   221        571    57       28.62295
# Treatment 205        451    36       27.57377

# Analysis of last-event-assisted WR (LWR; recommended), first-event-assisted WR (FWR), and naive WR (NWR):
#     Win prob Loss prob WR (95% CI)*      p-value
# LWR 50.3%    38.5%     1.31 (1.04, 1.64) 0.0233
# FWR 50.3%    38.5%      1.3 (1.04, 1.64) 0.0238
# NWR 46.9%    35.2%     1.33 (1.04,  1.7) 0.0229
# -----
# Total number of pairs:  45,305

tmp <- time_expr(WRrec(
  ID = df$patid,
  time = df$time,
  status = df$status,
  trt = df$trt_ab,
  strata = df$age60,
  naive = TRUE
))
timing$wr <- timing$wr + tmp$elapsed

length(unique(df[df$trt_ab == 1 & df$age60 == 1, "patid"])) *
  length(unique(df[df$trt_ab == 0 & df$age60 == 1, "patid"])) +
  length(unique(df[df$trt_ab == 1 & df$age60 == 0, "patid"])) *
    length(unique(df[df$trt_ab == 0 & df$age60 == 0, "patid"]))

#             N Rec. Event Death Med. Follow-up
# Control   221        571    57       28.62295
# Treatment 205        451    36       27.57377

# Analysis of last-event-assisted WR (LWR; recommended), first-event-assisted WR (FWR), and naive WR (NWR):
#     Win prob Loss prob WR (95% CI)*      p-value
# LWR 50.4%    38.2%     1.32 (1.05, 1.66) 0.0189
# FWR 50.4%    38.3%     1.32 (1.04, 1.66) 0.0202
# NWR   47%      35%     1.34 (1.05, 1.72) 0.0193
# -----
# Total number of pairs:  23,239

df |>
  group_by(patid) |>
  slice(1) |>
  ungroup() |>
  select(age60, trt_ab) |>
  table()

print_timing("HF-ACTION - JFM (all frailtyPenal fits)", timing$jfm)
print_timing("HF-ACTION - Win ratio", timing$wr)
