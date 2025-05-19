#' Function to compute sample size needed to achieve target power
#'
#' Computes sample size needed to achieve target power at one-sided Type I control level alp/2.
#' @param bta1 log-transform of rate ratio.
#' @param thta Variance of frailty parameter.
#' @param tau Expected follow-up time.
#' @param lam0 Baseline rate for control.
#' @param alp Two-sided alpha-level.
#' @param pow Target power.
#' @param ar Allocation ratio (Number control / Total)
#' @return The sample size required given the input assumptions to achieve target power.
#' @examples
#'
#' nNeeded(bta1 = log(0.8), thta = 1, tau = 0.8, lam0 = 3.5, alp = 0.05, pow = 0.8)
#'
#' if (require("dplyr") & require("tidyr")) {
#'
#'   assumptions = tibble(alp = 0.05) %>%
#'   crossing(
#'     tau = c(0.8,0.9, 1.0),
#'     RR = c(0.6, 0.7, 0.8),
#'     lam0 = c(3, 3.5),
#'     thta = c(2, 3, 4),
#'     pow = 0.8
#'   ) %>%
#'     mutate(N = nNeeded(bta1 = log(RR), thta = thta, tau = tau, lam0 = lam0, alp = alp, pow = pow))
#'
#'   assumptions %>% data.frame()
#'
#' }
#'
#' @export
nNeeded = function(bta1, thta, tau, lam0, alp = 0.05, pow = 0.8, ar = 0.5) {

  L = eventsNeeded(bta1 = bta1, thta = thta, tau = tau, lam0 = lam0, alp = alp, pow = pow, ar = ar)

  n = L / (ar + (1-ar) * exp(bta1)) / (tau * (lam0))

  return(n)
}
