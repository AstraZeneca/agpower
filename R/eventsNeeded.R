#' Number of events needed
#'
#' Function to compute number of events (L) needed to achieve power (pow) at one-sided Type I control level alpha/2 (alp/2).
#'
#' Assumes rate ratio < 1 is favourable to treatment.
#' @param bta1 log-transform of rate ratio.
#' @param thta Variance of frailty parameter.
#' @param tau Expected follow-up time.
#' @param lam0 Baseline rate for control.
#' @param alp Two-sided alpha-level.
#' @param pow Target power.
#' @param ar Allocation ratio (Number control / Total)
#' @return The number of events (L) needed.
#' @examples
#'
#' eventsNeeded(bta1 = log(0.8), thta = 1, tau = 0.8, lam0 = 3.5, alp = 0.05, pow = 0.8)
#'
#' if (require("dplyr") & require("tidyr")) {
#'
#'   assumptions = tibble(alp = 0.05) %>%
#'     crossing(
#'       tau = c(0.8,0.9, 1.0),
#'       RR = c(0.6, 0.7, 0.8),
#'       lam0 = c(3, 3.5),
#'       thta = c(2, 3, 4),
#'       pow = 0.8
#'     ) %>%
#'     mutate(
#'       L = eventsNeeded(bta1 = log(RR), thta = thta, tau = tau,
#'                        lam0 = lam0, alp = alp, pow = pow)
#'     )
#'
#'   assumptions %>% data.frame()
#'
#' }
#'
#' @export
eventsNeeded = function(bta1, thta, tau, lam0, alp = 0.05, pow = 0.8, ar = 0.5) {

  eta = 1 + thta * tau * (lam0) * (ar + (1-ar) * exp(2*bta1)) / (ar + (1-ar) * exp(bta1))

  zalp = stats::qnorm(1-alp/2, 0, 1)
  zpow = stats::qnorm(pow, 0, 1)

  L = 1/(1-ar)/ar * eta * ((zalp + zpow)/bta1)^2

  return(L)
}
