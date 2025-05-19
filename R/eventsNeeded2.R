#' Number events needed (fixed sample size)
#'
#' Function to compute number events (L) needed to achieve given power with fixed sample size N.
#'
#' Assumes rate ratio < 1 is favourable to treatment.
#' A negative estimated number of events indicates no number of events is sufficient under the input assumptions.
#' @param N Sample size.
#' @param bta1 log-transform of rate ratio.
#' @param thta Variance of frailty parameter.
#' @param alp Two-sided alpha-level.
#' @param pow Target power.
#' @param ar Allocation ratio (Number control / Total).
#' @param Lmessage If TRUE indicates how many estimates of L were negative before setting to Inf. Default FALSE.
#' @return The number of events (L) needed to achieve target power at one-sided Type I control level alpha/2, given the input assumptions.
#' @examples
#'
#' eventsNeeded2(N = 500, bta1 = log(0.8), thta = 2, alp = 0.05, pow = 0.8)
#' eventsNeeded2(N = 500, bta1 = log(0.8), thta = 3, alp = 0.05, pow = 0.8)
#'
#' if (require("dplyr") & require("tidyr")) {
#'
#'   assumptions = tibble(alp = 0.05) %>%
#'   crossing(
#'     RR = c(0.6, 0.7, 0.8),
#'     thta = c(2, 3, 4),
#'     pow = 0.8,
#'     N = c(500, 1000)
#'   ) %>%
#'     mutate(L = eventsNeeded2(N = N, bta1 = log(RR), thta = thta, alp, pow))
#'
#'   assumptions %>% data.frame()
#'
#' }
#'
#' @export
eventsNeeded2 = function(N, bta1, thta, alp = 0.05, pow = 0.8, ar = 0.5, Lmessage = FALSE) {

  zalp = stats::qnorm(1-alp/2, 0, 1)
  zpow = stats::qnorm(pow, 0, 1)

  num = N * (ar + (1-ar)*exp(bta1))^2

  den = num * ar * (1-ar) * (bta1^2)/(zalp + zpow)^2 - thta * (ar + (1-ar)*exp(2*bta1))

  L = num / den

  if (sum(L < 0) > 0) {
    nlz = sum(L < 0)
    L = ifelse(L < 0, Inf, L)
    if (Lmessage) warning(paste0(nlz, " out of ", length(L), " estimates of L were negative, setting these to Inf."))
  }

  return(L)
}
