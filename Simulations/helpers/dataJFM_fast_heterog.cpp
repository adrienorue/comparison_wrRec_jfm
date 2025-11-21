// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

enum Baseline
{
  EXPONENTIAL = 0,
  LOGLOGISTIC = 1
};
enum Effect
{
  PH = 0,
  AH = 1
};

inline double runif01() { return R::runif(0.0, 1.0); }
inline double pos(double x) { return (x > 1e-12) ? x : 1e-12; }

// ---------- Proportional hazards samplers ----------
inline double sample_PH_exp(double rate0, double accel)
{
  double u = runif01();
  return -std::log(u) / pos(rate0 * accel);
}
inline double sample_PH_llogis(double scale, double shape, double accel)
{
  // t = scale * (u^{-1/accel} - 1)^{1/shape}
  double u = runif01();
  return scale * std::pow(std::pow(u, -1.0 / accel) - 1.0, 1.0 / shape);
}

// ---------- Additive hazards sampler (exponential only) ----------
inline double sample_AH_exp(double rate0, double add)
{
  // add = W + eta or alpha*W + eta
  double u = runif01();
  return -std::log(u) / pos(rate0 + add);
}

// [[Rcpp::export]]
List dataJFM_fast(
    int nSubjects,
    double theta,
    double alpha,
    double fixedCensor,
    double lambda0,
    double r0,
    NumericVector betaGap,
    double betaDeathStar,
    Nullable<int> seed,
    int baseline_type,
    int effect_type,
    double shapeRec,
    double scaleRec,
    double shapeDeath,
    double scaleDeath)
{
  if (seed.isNotNull())
  {
    Environment base("package:base");
    Function set_seed = base["set.seed"];
    set_seed(as<int>(seed));
  }
  RNGScope rngScope;

  // Disallow AH + log-logistic
  if (effect_type == AH && baseline_type == LOGLOGISTIC)
  {
    stop("Additive hazards is implemented only with exponential baseline. "
         "Use baseline_type = 0 when effect_type = 1.");
  }

  // Basic checks for log-logistic params (only when used)
  if (baseline_type == LOGLOGISTIC)
  {
    if (NumericVector::is_na(shapeRec) || NumericVector::is_na(scaleRec) ||
        NumericVector::is_na(shapeDeath) || NumericVector::is_na(scaleDeath))
      stop("Provide finite log-logistic shape/scale parameters.");
    if (shapeRec <= 0.0 || scaleRec <= 0.0 || shapeDeath <= 0.0 || scaleDeath <= 0.0)
      stop("Log-logistic shape/scale must be positive.");
  }
  if (betaGap.size() < 2)
    stop("betaGap must have length >= 2 (z1, z2).");

  // 1) Subject-level draws
  NumericVector frailty = rgamma(nSubjects, 1.0 / theta, theta); // mean=1, var=theta
  NumericVector z1 = rbinom(nSubjects, 1, 0.5);
  NumericVector z2 = rbinom(nSubjects, 1, 0.5);

  // 2) Allocate long data buffers
  std::vector<int> id;
  std::vector<double> time;
  std::vector<int> recurrence;
  std::vector<int> death;
  std::vector<int> z1Long;
  std::vector<int> z2Long;

  id.reserve(4 * nSubjects);
  time.reserve(4 * nSubjects);
  recurrence.reserve(4 * nSubjects);
  death.reserve(4 * nSubjects);
  z1Long.reserve(4 * nSubjects);
  z2Long.reserve(4 * nSubjects);

  // Subject-level summaries
  NumericVector eventEndTime(nSubjects);
  IntegerVector died(nSubjects);
  NumericVector censorTime(nSubjects, fixedCensor);

  // 3) Simulate per subject
  for (int i = 0; i < nSubjects; ++i)
  {
    const double W = frailty[i];
    const int zz1 = static_cast<int>(z1[i]);
    const int zz2 = static_cast<int>(z2[i]);

    const double etaR = betaGap[0] * zz1 + betaGap[1] * zz2; // recurrent
    const double etaD = betaDeathStar * zz1;                 // death

    // Baseline risk change for z2=1
    const double delta = 0.1;
    const bool z2_one = (zz2 == 1);
    const double rateR0 = z2_one ? (r0 + delta) : r0 - delta;           // EXP baselines
    const double rateD0 = z2_one ? (lambda0 + delta) : lambda0 - delta; // EXP baselines
    const double baseMult = z2_one ? (1.0 + delta) : 1.0;               // INUTILE, pas utilise (loglogistic)

    // Death time
    double Tdeath;
    if (effect_type == PH)
    {
      double accelD = std::pow(W, alpha) * std::exp(etaD);
      if (baseline_type == EXPONENTIAL)
      {
        Tdeath = sample_PH_exp(rateD0, accelD);
      }
      else
      { // LOGLOGISTIC
        Tdeath = sample_PH_llogis(scaleDeath, shapeDeath, accelD * baseMult);
      }
    }
    else
    { // AH + EXPONENTIAL
      double addD = alpha * W + etaD;
      // add +0.5 to baseline rate if z2=1
      Tdeath = sample_AH_exp(rateD0, addD);
    }

    const double stopTime = std::min(fixedCensor, Tdeath);
    eventEndTime[i] = stopTime;
    died[i] = static_cast<int>(Tdeath < fixedCensor);

    // Recurrent events until stopTime
    double t = 0.0;
    while (true)
    {
      double gap;
      if (effect_type == PH)
      {
        double accelR = W * std::exp(etaR);
        if (baseline_type == EXPONENTIAL)
        {
          gap = sample_PH_exp(rateR0, accelR);
        }
        else
        { // LOGLOGISTIC
          gap = sample_PH_llogis(scaleRec, shapeRec, accelR * baseMult);
        }
      }
      else
      { // AH + EXPONENTIAL
        double addR = W + etaR;
        gap = sample_AH_exp(rateR0, addR);
      }

      double tnext = t + gap;
      if (tnext < stopTime)
      {
        id.push_back(i + 1);
        time.push_back(tnext);
        recurrence.push_back(1);
        death.push_back(0);
        z1Long.push_back(zz1);
        z2Long.push_back(zz2);
        t = tnext;
      }
      else
      {
        id.push_back(i + 1);
        time.push_back(stopTime);
        recurrence.push_back(0);
        death.push_back(static_cast<int>(Tdeath < fixedCensor));
        z1Long.push_back(zz1);
        z2Long.push_back(zz2);
        break;
      }
    }
  }

  // 4) Return
  return List::create(
      _["longData"] = DataFrame::create(
          _["id"] = wrap(id),
          _["time"] = wrap(time),
          _["recurrence"] = wrap(recurrence),
          _["death"] = wrap(death),
          _["z1"] = wrap(z1Long),
          _["z2"] = wrap(z2Long)),
      _["subjectData"] = DataFrame::create(
          _["id"] = seq_len(nSubjects),
          _["omega"] = frailty,
          _["z1"] = as<IntegerVector>(z1),
          _["z2"] = as<IntegerVector>(z2),
          _["deathTime"] = eventEndTime, // min(Tdeath, C)
          _["censorTime"] = censorTime,
          _["eventEndTime"] = eventEndTime,
          _["died"] = died),
      _["truth"] = List::create(
          _["beta1"] = betaGap[0],
          _["beta2"] = betaGap[1],
          _["betaStar"] = betaDeathStar,
          _["alpha"] = alpha,
          _["theta"] = theta,
          _["r0"] = r0,
          _["lambda0"] = lambda0,
          _["effect_type"] = effect_type,     // 0=PH, 1=AH
          _["baseline_type"] = baseline_type, // 0=exp, 1=loglogistic
          _["shapeRec"] = shapeRec,
          _["scaleRec"] = scaleRec,
          _["shapeDeath"] = shapeDeath,
          _["scaleDeath"] = scaleDeath));
}
