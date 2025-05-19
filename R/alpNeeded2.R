#' Function to compute alpha needed (fixed sample size)
#'
#' Function to compute two-sided alpha needed to achieve target power given a rate ratio. Useful for computing probability to achieve hurdles.
#'
#' This function computes the two-sided alpha alp.
#' Function assumes a rate ratio < 1 is favourable to treatment.
#' @param N Sample size.
#' @param bta1 log-transform of rate ratio.
#' @param thta Variance of frailty parameter.
#' @param L Number of events.
#' @param pow Target power.
#' @param ar Allocation ratio (Number control / Total).
#' @return The two-sided alpha level.
#' @examples
#'
#' # alpha needed to achieve multiple powers given rate ratio (and other input).
#' alpNeeded2(N = 1000, bta1 = log(0.8), thta = 2, L = 1000, pow = c( .7, .8))
#'
#' # alpha needed for many inputs
#' if (require("dplyr") & require("tidyr")) {
#'
#   assumptions = tibble(RR = 0.8) %>%
#     crossing(
#       thta = c(2, 3, 4),
#       pow = c(0.7, 0.8),
#       N = c(500, 1000),
#       L = 1000
#     ) %>%
#     mutate(
#       alp = alpNeeded2(N = N, bta1 = log(RR), thta = thta, L = L, pow = pow)
#     )
# #  assumptions %>% data.frame()
#'
#' }
#'
#'
#' @export
alpNeeded2 = function(N, bta1, thta, L, pow = 0.8, ar = 0.5) {

  zalp = zalp2(N = N, bta1 = bta1, thta = thta, L = L, pow = pow, ar = ar)

  alp = (1 - stats::pnorm(zalp)) * 2

  return(alp)
}

zalp2 = function(N, bta1, thta, L, pow, ar = 0.5) {

  c1 = ar + (1-ar) * exp(bta1)
  dnum = ar * (1-ar) * N * L
  dden = N * c1^2 + thta * L * (ar + (1-ar) * exp(2*bta1))
  c2 = -bta1 * c1 * sqrt(dnum/dden)

  zpow = stats::qnorm(pow)

  zalp = c2 - zpow

  return(zalp)
}





