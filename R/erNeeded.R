#' Baseline event rate needed (fixed sample size)
#'
#' Function to compute baseline (control) event rate needed to achieve given power at one-sided Type I control level alp/2.
#'
#' Assumes rate ratio < 1 is favourable to treatment.
#' A negative estimated event rate indicates no event rate is sufficient under the input assumptions.
#' @param N Sample size.
#' @param bta1 log-transform of rate ratio.
#' @param thta Variance of frailty parameter.
#' @param tau Expected follow-up time.
#' @param alp Two-sided alpha-level.
#' @param pow Target power.
#' @param ar Allocation ratio (Number control / Total).
#' @param lam0warning If TRUE indicates how many estimates of lam0 were negative before setting to Inf. Default FALSE.
#' @return The baseline event rate needed to achieve target power at one-sided Type I control level alpha/2, given the input assumptions.
#' @examples
#'
#' erNeeded(N = 500, bta1 = log(0.6), thta = 2, tau = 1, alp = 0.05, pow = 0.8)
#' erNeeded(N = 500, bta1 = log(0.6), thta = 3, tau = 1, alp = 0.05, pow = 0.8)
#'
#'
#' if (require("dplyr") & require("tidyr")) {
#'
#'   assumptions = tibble(alp = 0.05) %>%
#'   crossing(
#'     tau = c(0.8,0.9, 1.0),
#'     RR = c(0.6, 0.7, 0.8),
#'     thta = c(2, 3, 4),
#'     pow = 0.8,
#'     N = c(500, 1000)
#'   ) %>%
#'     mutate(er = erNeeded(N = N, bta1 = log(RR), thta = thta, tau = tau, alp, pow))
#'
#'   assumptions %>% data.frame()
#'
#' }
#'
#' @export
erNeeded = function(N, bta1, thta, tau, alp = 0.05, pow = 0.8, ar = 0.5, lam0warning = FALSE) {

  zalp = stats::qnorm(1-alp/2, 0, 1)
  zpow = stats::qnorm(pow, 0, 1)

  num = (ar + (1-ar)*exp(bta1))

  den = (N * (ar + (1-ar)*exp(bta1))^2 * ar * (1-ar) * (bta1^2)/(zalp + zpow)^2 - thta * (ar + (1-ar)*exp(2*bta1))) * tau

  lam0 = num / den

  if (sum(lam0 < 0) > 0) {
    nlz = sum(lam0 < 0)
    lam0 = ifelse(lam0 < 0, Inf, lam0)
    if (lam0warning) warning(paste0(nlz, " out of ", length(lam0), " estimates of lam0 were negative, setting these to Inf."))
  }

  return(lam0)
}

