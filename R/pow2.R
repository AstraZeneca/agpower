#' Power for LWYY (fixed sample size)
#'
#' Function to compute power at one-sided Type I control level alp/2.
#'
#' @param N Sample size.
#' @param bta1 log-transform of rate ratio.
#' @param thta Variance of frailty parameter.
#' @param L Number of events.
#' @param alp Two-sided alpha-level.
#' @param ar Allocation ratio (Number control / Total)
#' @return The power given the input assumptions.
#' @examples
#'
#' pow2(N = 500, bta1 = log(0.8), thta = 2, L = 1000, alp = 0.05)
#' pow2(N = 500, bta1 = log(0.8), thta = 3, L = 1000, alp = 0.05)
#'
#' if (require("dplyr") & require("tidyr")) {
#'
#'   assumptions = tibble(alp = 0.05) %>%
#'   crossing(
#'     thta = c(1),
#'     RR = c(0.6, 0.7, 0.8),
#'     L = c(1000, 1500, 2000),
#'     N = c(500, 1000)
#'   ) %>%
#'     mutate(pow = pow2(N = N, bta1 = log(RR), thta = thta, L = L, alp = alp))
#'
#'   assumptions %>% data.frame()
#'
#' }
#'
#' @export
pow2 = function(N, bta1, thta, L, alp = 0.05, ar = 0.5) {

  zpow = zpow2(N = N, bta1 = bta1, thta = thta, L = L, alp = alp, ar = ar)

  pow = stats::pnorm(zpow)

  return(pow)
}

zpow2 = function(N, bta1, thta, L, alp = 0.05, ar = 0.5) {

  c1 = ar + (1-ar) * exp(bta1)
  dnum = ar * (1-ar) * N * L
  dden = N * c1^2 + thta * L * (ar + (1-ar) * exp(2*bta1))
  c2 = -bta1 * c1 * sqrt(dnum/dden)

  zalp = stats::qnorm(1-alp/2)

  zpow = c2 - zalp

  return(zpow)
}
