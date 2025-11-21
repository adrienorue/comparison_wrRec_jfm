# ------------------------------------------------------------
# Set current script location as working directory
# ------------------------------------------------------------
library(rstudioapi)
a <- getActiveDocumentContext()
setwd(dirname(a$path))
remove(a)

source("helpers/genCpp.R")

Rcpp::sourceCpp(
    "helpers/dataJFM_fast.cpp",
    verbose = FALSE
)

# # UNCOMMENT IF YOU WANT TO USE THE "HETEROGENEOUS" VERSION
# # TO INCREASE STRATIFIED LWR POWER

# Rcpp::sourceCpp(
#   "helpers/dataJFM_fast_heterog.cpp", # -> -0.1 if Z2 = 0, +0.1 if Z2 = 1
#   verbose = FALSE
# )

# Datasets generation simulations
# ---- Scénario 1 : alpha = 1, theta = 0.5, r0=1.5, lambda0 = 1/2 ----

# for the MC simulations
genDatas(
    nbDatasets = 500,
    outputDir = "datasets scenario 1/"
)

# for the "asymptotic" LWR value
genDatas(
    nbDatasets = 50,
    outputDir = "datasets scenario 1BIG/",
    nSubjects = 1e5
)

# ---- Scénario 2 : theta = 0.01 ----

genDatas(
    nbDatasets = 500,
    outputDir = "datasets scenario 2/",
    theta = 0.01
)

genDatas(
    nbDatasets = 50,
    outputDir = "datasets scenario 2BIG/",
    nSubjects = 1e5,
    theta = 0.01
)


# ---- Scénario 3 : theta = 1 ----

genDatas(
    nbDatasets = 500,
    outputDir = "datasets scenario 3/",
    theta = 1
)

genDatas(
    nbDatasets = 50,
    outputDir = "datasets scenario 3BIG/",
    nSubjects = 1e5,
    theta = 1
)

# ---- Scénario 4 : r0 = 1/5, lambda0 = 2 ----

genDatas(
    nbDatasets = 500,
    outputDir = "datasets scenario 4/",
    lambda0 = 2,
    r0 = 1 / 5
)

genDatas(
    nbDatasets = 50,
    outputDir = "datasets scenario 4BIG/",
    nSubjects = 1e5,
    lambda0 = 2,
    r0 = 1 / 5
)

# ---- Scénario 5 : r0 = 2, lambda0 = 1/7 ----

genDatas(
    nbDatasets = 500,
    outputDir = "datasets scenario 5/",
    lambda0 = 1 / 7,
    r0 = 2
)

genDatas(
    nbDatasets = 50,
    outputDir = "datasets scenario 5BIG/",
    nSubjects = 1e5,
    lambda0 = 1 / 7,
    r0 = 2
)

# ---- Scénarios 6 :  ----

# # ---- Scénario 6/ log-logistic ----

genDatas(
    nbDatasets = 500,
    outputDir = "datasets scenario 6b/",
    baseline = "loglogistic",
    shapeRec = 1.4,
    scaleRec = 2.0,
    shapeDeath = 1.2,
    scaleDeath = 3.0
)

genDatas(
    nbDatasets = 50,
    outputDir = "datasets scenario 6bBIG/",
    nSubjects = 1e5,
    baseline = "loglogistic",
    shapeRec = 1.4,
    scaleRec = 2.0,
    shapeDeath = 1.2,
    scaleDeath = 3.0
)

# # ---- Scénario bonus : scénario 5 & rate of -0.1 if Z2=0 / +0.1 if Z2 = 1  ----

# needs "dataJFM_fast_heterog.cpp"
genDatas(
    nbDatasets = 500,
    outputDir = "datasets scenario bonus/",
    lambda0 = 1 / 7,
    r0 = 2
)
