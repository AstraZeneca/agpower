#' Power for LWYY (fixed sample size)
#'
#' Function to compute power at one-sided Type I control level alp/2.
#'
#' Note: Approximation breaks down in no event scenario. For example, pow(N=1000, bta1 = log(1), thta = 1, lam0 = 0, tau = 1, alp = 0.05) returns a power of 0.025
#' @param N Sample size.
#' @param bta1 log-transform of rate ratio.
#' @param thta Variance of frailty parameter.
#' @param tau Expected follow-up time.
#' @param lam0 Baseline rate for control.
#' @param alp Two-sided alpha-level.
#' @param ar Allocation ratio (Number control / Total)
#' @return The power given the input assumptions.
#' @examples
#'
#' pow(N = 500, bta1 = log(0.8), thta = 2, tau = 1, lam0 = 1.1, alp = 0.05)
#' pow(N = 500, bta1 = log(0.8), thta = 3, tau = 1, lam0 = 1.1, alp = 0.05)
#'
#' if (require("dplyr") & require("tidyr")) {
#'
#'   assumptions = tibble(alp = 0.05) %>%
#'   crossing(
#'     tau = c(0.8,0.9, 1.0),
#'     RR = c(0.6, 0.7, 0.8),
#'     lam0 = c(3, 3.5),
#'     thta = c(2, 3, 4),
#'     N = c(500, 1000)
#'   ) %>%
#'     mutate(pow = pow(N = N, bta1 = log(RR), thta = thta, tau = tau, lam0 = lam0, alp = alp))
#'
#'   assumptions %>% data.frame()
#'
#' }
#'
#' @export
pow = function(N, bta1, thta, tau, lam0, alp = 0.05, ar = 0.5) {

  zpow = zpow(N = N, bta1 = bta1, thta = thta, tau = tau, lam0 = lam0, alp = alp, ar = ar)

  pow = stats::pnorm(zpow)

  return(pow)
}

zpow = function(N, bta1, thta, tau, lam0, alp = 0.05, ar = 0.5) {

  L = N * tau * lam0 * (ar + (1-ar) * exp(bta1))

  eta = 1 + thta * tau * (lam0) * (ar + (1-ar) * exp(2*bta1)) / (ar + (1-ar) * exp(bta1))

  zalp = stats::qnorm(1-alp/2)

  zpow = -(bta1) * sqrt(L *(ar*(1-ar)) / eta) - zalp

  return(zpow)
}
