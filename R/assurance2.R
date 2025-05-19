#' Assurance (expected power) for LWYY
#'
#' Function to compute assurance given fixed sample size N and number of events, under input assumptions. Assumes a log-normal prior distribution.
#' Here assurance is taken to mean probability to achieve statistical significance (p-value < alp).
#'
#' If working with a blinded estimate of frailty variance (i.e. misspecified model), it is recommended to use frailty.type = "blind".
#' In which case the frailty variance (i.e. model with treatment effect) is derived using thta and the quantile drawn from the prior distribution of the log-rate ratio.
#' If working with an estimate of frailty variance, should use frailty.type = "unblind" instead.
#'
#' Function assumes a rate ratio < 1 is favourable to treatment.
#' @param N Sample size.
#' @param bta1 log-transform of rate ratio.
#' @param bta1_sd assumed standard deviation of log(rate ratio)
#' @param thta Variance of frailty parameter. If frailty.type = "blind", assumes thta derives from pooled model; if frailty.type = "unblind" assumes thta is from correctly specified model. Default "unblind".
#' @param L Number of events
#' @param alp Two-sided alpha-level.
#' @param ar Allocation ratio (Number control / Total).
#' @param ns Maximum number of subintervals (if method = "integration") or Number of draws from prior distribution (if method = "montecarlo").
#' @param method Whether to numerically solve ("integration") or estimate by random draws ("montecarlo"). Defaults to numerical.
#' @param frailty.type Indicates whether frailty variance is based on blinded information ("blind") or unblinded ("unblind"). Default "unblind".
#' @param thtawarning If TRUE indicates how many estimates of theta were negative before setting to 0. Default FALSE.
#' @return The assurance given the input assumptions.
#' @examples
#'
#' assurance2(N = 500, bta1 = log(0.8), bta1_sd = 1, thta = 2, L = 1000, alp = 0.05,
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
#'     L = c(500, 1000, 1500)
#'   ) %>%
#'     mutate(pow = pow2(N = N, bta1 = log(RR), thta = thta, L = L, alp = alp)) %>%
#'     mutate(
#'       assurance_in_blind = assurance2(N = N, bta1 = log(RR), bta1_sd = bta1_sd, thta = thta,
#'                                      L = L, alp = alp, ns = 1000,
#'                                      frailty.type = "blind")
#'     ) %>%
#'     mutate(
#'       assurance_mc_blind = assurance2(N = N, bta1 = log(RR), bta1_sd = bta1_sd, thta = thta,
#'                                      L = L, alp = alp, ns = 1000,
#'                                      method = "monte",
#'                                      frailty.type = "blind")
#'     ) %>%
#'     mutate(
#'       assurance_in_unblind = assurance2(N = N, bta1 = log(RR), bta1_sd = bta1_sd, thta = thta,
#'                                        L = L, alp = alp, ns = 1000)
#'     ) %>%
#'     mutate(
#'       assurance_mc_unblind = assurance2(N = N, bta1 = log(RR), bta1_sd = bta1_sd, thta = thta,
#'                                        L = L, alp = alp, ns = 1000)
#'     )
#'
#'   assumptions %>% data.frame()
#'
#' }
#' @export
assurance2 = function(N, bta1, bta1_sd, thta, L, alp = 0.05, ar = 0.5, ns = 1000, method = c("integration", "montecarlo"), frailty.type = c("unblind", "blind"), thtawarning = FALSE) {
  method = match.arg(method)
  frailty.type = match.arg(frailty.type)

  if (method == "integration") {
    usefunc = assurance2_int
  } else if (method == "montecarlo") {
    usefunc = assurance2_mc
  }

  pts = usefunc(N = N, bta1 = bta1, bta1_sd = bta1_sd, thta = thta, L = L, alp = alp, ar = ar, ns = ns, frailty.type = frailty.type, thtawarning = thtawarning)

  return(pts)
}

assurance2_int = function(N, bta1, bta1_sd, thta, L, alp = 0.05, ar = 0.5, ns = 1000, frailty.type = c("unblind", "blind"), thtawarning = FALSE) {

  X = data.frame(bta1, bta1_sd, N, thta, L, alp, ar)
  pts = vapply(split(X, seq(nrow(X))), function(x) {
    power = function(p) {
      q = stats::qnorm(p, mean = x$bta1, sd = x$bta1_sd)
      if (frailty.type == "blind") {
        thta_ = thtap2thta(bta1 = q, thtap = x$thta, ar = x$ar, thtawarning = thtawarning)
      } else if (frailty.type == "unblind") {
        thta_ = x$thta
      }
      stats::pnorm(zpow2(N = x$N, bta1 = q, thta = thta_, L = x$L, alp = x$alp, ar = x$ar))
    }

    stats::integrate(power, lower = 0, upper = 1, subdivisions = ns)$value

  }, numeric(1), USE.NAMES = FALSE)


  return(pts)
}

assurance2_mc = function(N, bta1, bta1_sd, thta, L, alp = 0.05, ar = 0.5, ns = 1000, frailty.type = c("unblind", "blind"), thtawarning = FALSE) {

  X = data.frame(bta1, bta1_sd, N, thta, L, alp, ar)
  pts = vapply(split(X, seq(nrow(X))), function(x) {
    qs = stats::rnorm(ns, mean = x$bta1, sd = x$bta1_sd)
    if (frailty.type == "blind") {
      thta_ = thtap2thta(bta1 = qs, thtap = x$thta, ar = x$ar, thtawarning = thtawarning)
    } else if (frailty.type == "unblind") {
      thta_ = x$thta
    }
    mean(stats::pnorm(zpow2(N = x$N, bta1 = qs, thta = thta_, L = x$L, alp = x$alp, ar = x$ar)))
  }, numeric(1), USE.NAMES = FALSE)

  return(pts)
}

