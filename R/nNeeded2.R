#' Function to compute sample size needed to achieve target power
#'
#' Computes sample size needed to achieve target power at one-sided Type I control level alp/2.
#' A negative estimated sample size indicates no sample size is sufficient under the input assumptions.
#' @param bta1 log-transform of rate ratio.
#' @param thta Variance of frailty parameter.
#' @param L Number of events
#' @param alp Two-sided alpha-level.
#' @param pow Target power.
#' @param ar Allocation ratio (Number control / Total)
#' @return The sample size required given the input assumptions to achieve target power.
#' @examples
#'
#' nNeeded2(bta1 = log(0.8), thta = 1, L = 1000, alp = 0.05, pow = 0.8)
#'
#' if (require("dplyr") & require("tidyr")) {
#'
#'   assumptions = tibble(alp = 0.05) %>%
#'   crossing(
#'     L = c(500, 1000, 1500),
#'     RR = c(0.6, 0.7, 0.8),
#'     thta = c(2, 3, 4),
#'     pow = 0.8
#'   ) %>%
#'     mutate(N = nNeeded2(bta1 = log(RR), thta = thta, L = L, alp = alp, pow = pow))
#'
#'   assumptions %>% data.frame()
#'
#' }
#'
#' @export
nNeeded2 = function(bta1, thta, L, alp = 0.05, pow = 0.8, ar = 0.5) {

  zalp = stats::qnorm(1-alp/2)
  zpow = stats::qnorm(pow)

  num = thta * L * (ar + (1-ar)*exp(2*bta1))
  den = ar * (1-ar) * L * (bta1 / (zalp + zpow))^2 - 1

  n = num / den

  return(n)
}
