#' Assurance (expected power) for LWYY
#'
#' Function to compute assurance given fixed sample size N and follow-up tau, under input assumptions. Assumes a log-normal prior distribution.
#' Here assurance is taken to mean probability to achieve statistical significance (p-value < alp).
#'
#' If working with a blinded estimate of frailty variance (i.e. misspecified model), it is recommended to use frailty.type = "blind" and lam.type = "pool".
#' In which case the frailty variance (i.e. model with treatment effect) is derived using thta and the quantile drawn from the prior distribution of the log-rate ratio.
#' If working with an estimate of frailty variance, should use frailty.type = "unblind" instead.
#'
#' Function assumes a rate ratio < 1 is favourable to treatment.
#' @param N Sample size.
#' @param bta1 log-transform of rate ratio.
#' @param bta1_sd assumed standard deviation of log(rate ratio)
#' @param thta Variance of frailty parameter. If frailty.type = "blind", assumes thta derives from pooled model; if frailty.type = "unblind" assumes thta is from correctly specified model. Default "unblind".
#' @param tau Expected follow-up time.
#' @param lam Event rate. If lam.type = "pool", assumes lam is pooled rate; if lam.type = "base", assumes lam is baseline control event rate. Default "base".
#' @param alp Two-sided alpha-level.
#' @param ar Allocation ratio (Number control / Total).
#' @param ns Maximum number of subintervals (if method = "integration") or Number of draws from prior distribution (if method = "montecarlo").
#' @param method Whether to numerically solve ("integration") or estimate by random draws ("montecarlo"). Defaults to numerical.
#' @param frailty.type Indicates whether frailty variance is based on blinded information ("blind") or unblinded ("unblind"). Default "unblind".
#' @param lam.type Indicates whether event rate is based on control rate ("base") or pooled rate ("pool"). Default "base".
#' @param thtawarning If TRUE indicates how many estimates of theta were negative before setting to 0. Default FALSE.
#' @return The assurance given the input assumptions.
#' @examples
#'
#' assurance(N = 500, bta1 = log(0.8), bta1_sd = 1, thta = 2, tau = 1, lam = 1.1, alp = 0.05,
#'           ns = 100000)
#'
#'
#' if (require("dplyr") & require("tidyr")) {
#'
#'   assumptions = tibble(alp = 0.05) %>%
#'   crossing(
#'     N = c(500, 1000),
#'     RR = c(0.6, 0.7, 0.8),
#'     bta1_sd = 1,
#'     thta = c(2, 3, 4),
#'     tau = c(0.8,0.9, 1.0),
#'     lam0 = c(3, 3.5)
#'   ) %>%
#'     mutate(pow = pow(N = N, bta1 = log(RR), thta = thta, tau = tau, lam = lam0, alp = alp)) %>%
#'     mutate(
#'       assurance_in_blind = assurance(N = N, bta1 = log(RR), bta1_sd = bta1_sd, thta = thta,
#'                                      tau = tau, lam = lam0, alp = alp, ns = 1000,
#'                                      frailty.type = "blind", lam.type = "pool")
#'     ) %>%
#'     mutate(
#'       assurance_mc_blind = assurance(N = N, bta1 = log(RR), bta1_sd = bta1_sd, thta = thta,
#'                                      tau = tau, lam = lam0, alp = alp, ns = 1000,
#'                                      method = "monte",
#'                                      frailty.type = "blind", lam.type = "pool")
#'     ) %>%
#'     mutate(
#'       assurance_in_unblind = assurance(N = N, bta1 = log(RR), bta1_sd = bta1_sd, thta = thta,
#'                                        tau = tau, lam = lam0, alp = alp, ns = 1000)
#'     ) %>%
#'     mutate(
#'       assurance_mc_unblind = assurance(N = N, bta1 = log(RR), bta1_sd = bta1_sd, thta = thta,
#'                                        tau = tau, lam = lam0, alp = alp, ns = 1000)
#'     )
#'
#'   assumptions %>% data.frame()
#'
#' }
#' @export
assurance = function(N, bta1, bta1_sd, thta, tau, lam, alp = 0.05, ar = 0.5, ns = 1000, method = c("integration", "montecarlo"), frailty.type = c("unblind", "blind"), lam.type = c("base", "pool"), thtawarning = FALSE) {
  method = match.arg(method)
  frailty.type = match.arg(frailty.type)
  lam.type = match.arg(lam.type)

  if (method == "integration") {
    usefunc = assurance_int
  } else if (method == "montecarlo") {
    usefunc = assurance_mc
  }

  pts = usefunc(N = N, bta1 = bta1, bta1_sd = bta1_sd, thta = thta, tau = tau, lam = lam, alp = alp, ar = ar, ns = ns, frailty.type = frailty.type, lam.type = lam.type, thtawarning = thtawarning)

  return(pts)
}

assurance_int = function(N, bta1, bta1_sd, thta, tau, lam, alp = 0.05, ar = 0.5, ns = 1000, frailty.type = c("unblind", "blind"), lam.type = c("base", "pool"), thtawarning = FALSE) {

  X = data.frame(bta1, bta1_sd, N, thta, tau, lam, alp, ar)
  pts = vapply(split(X, seq(nrow(X))), function(x) {
    power = function(p) {
      q = stats::qnorm(p, mean = x$bta1, sd = x$bta1_sd)
      if (frailty.type == "blind") {
        thta_ = thtap2thta(bta1 = q, thtap = x$thta, ar = x$ar, thtawarning = thtawarning)
      } else if (frailty.type == "unblind") {
        thta_ = x$thta
      }
      if (lam.type == "pool") {
        lam0_ = x$lam / (x$ar + (1-x$ar) * exp(q))
      } else if (lam.type == "base") {
        lam0_ = x$lam
      }
      stats::pnorm(zpow(N = x$N, bta1 = q, thta = thta_, tau = x$tau, lam0 = lam0_, alp = x$alp, ar = x$ar))
    }

    stats::integrate(power, lower = 0, upper = 1, subdivisions = ns)$value

  }, numeric(1), USE.NAMES = FALSE)


  return(pts)
}

assurance_mc = function(N, bta1, bta1_sd, thta, tau, lam, alp = 0.05, ar = 0.5, ns = 1000, frailty.type = c("unblind", "blind"), lam.type = c("base", "pool"), thtawarning = FALSE) {

  X = data.frame(bta1, bta1_sd, N, thta, tau, lam, alp, ar)
  pts = vapply(split(X, seq(nrow(X))), function(x) {
    qs = stats::rnorm(ns, mean = x$bta1, sd = x$bta1_sd)
    if (frailty.type == "blind") {
      thta_ = thtap2thta(bta1 = qs, thtap = x$thta, ar = x$ar, thtawarning = thtawarning)
    } else if (frailty.type == "unblind") {
      thta_ = x$thta
    }
    if (lam.type == "pool") {
      lam0_ = x$lam / (x$ar + (1-x$ar) * exp(qs))
    } else if (lam.type == "base") {
      lam0_ = x$lam
    }
    mean(stats::pnorm(zpow(N = x$N, bta1 = qs, thta = thta_, tau = x$tau, lam0 = lam0_, alp = x$alp, ar = x$ar)))
  }, numeric(1), USE.NAMES = FALSE)

  return(pts)
}

