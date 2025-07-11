#' Function to compute log rate ratio needed (fixed sample size)
#'
#' Function to compute log rate ratio needed to achieve target power at one-sided Type I control level alp/2. Useful to compute critical value (set pow = 0.5).
#'
#' This function computes the log rate ratio bta1 as the root of the equation pow(bta1, N, thta_, tau, lam0_, alp, ar) - pow = 0, where thta_ and lam0_ also depend on bta1 if using estimates from blinded analyses.
#' If frailty.type = "blind": thta_ = thtap2thta(thta, bta1); otherwise if frailty = "unblind": thta_ = thta. If lam.type = "pool" then lam0_ = lam * 2 / (1 + exp(bta1)); otherwise if lam.type = "base": lam0_ = lam.
#' Function assumes a rate ratio < 1 is favourable to treatment.
#' @param N Sample size.
#' @param thta Variance of frailty parameter. If frailty.type = "blind", assumes thta derives from pooled model; if frailty.type = "unblind" assumes thta is from correctly specified model. Default "unblind".
#' @param tau Expected follow-up time.
#' @param lam Event rate. If lam.type = "pool", assumes lam is pooled rate; if lam.type = "base", assumes lam is baseline control event rate. Default "base".
#' @param alp Two-sided alpha-level.
#' @param pow Target power.
#' @param ar Allocation ratio (Number control / Total)
#' @param frailty.type Indicates whether frailty variance is based on blinded information ("blind") or unblinded ("unblind"). Default "unblind".
#' @param lam.type Indicates whether event rate is based on control rate ("base") or pooled rate ("pool"). Default "base".
#' @param interval Initial search interval for bta1.
#' @return The log rate ratio.
#' @examples
#'
#' # Based on unblinded estimates
#' btaNeeded(N = 1000, thta = 2, tau = 1, lam = 1.1, alp = c(0.01, 0.05), pow = c(.5))
#' exp(btaNeeded(N = 1000, thta = 2, tau = 1, lam = 1.1, alp = c(0.01, 0.05), pow = c(.5)))
#'
#' # Based on blinded estimates
#' btaNeeded(N = 1000, thta = 2, tau = 1, lam = 0.7, alp = c(0.01, 0.05), pow = c( .5),
#'   frailty.type = "bl", lam.type = "po")
#' exp(btaNeeded(N = 1000, thta = 2, tau = 1, lam = 0.7, alp = c(0.01, 0.05), pow = c( .5),
#'                frailty.type = "bl", lam.type = "po"))
#'
#' # Based on blinded estimates
#' if (require("dplyr") & require("tidyr")) {
#'
#'   assumptions = tibble(alp = 0.05) %>%
#'     crossing(
#'       thta = c(2, 3, 4),
#'       lam = 1.1,
#'       pow = c(0.5, 0.8),
#'       N = c(500, 1000),
#'       tau = 1
#'     ) %>%
#'     mutate(
#'       bta1 = btaNeeded(N = N, thta = thta, tau = tau, lam = lam, alp = alp, pow = pow),
#'       RR = exp(bta1)
#'     )
#
#'   assumptions %>% data.frame()
#'
#' }
#'
#' @export
btaNeeded = function(N, thta, tau, lam, alp = 0.05, pow = 0.8, ar = 0.5, frailty.type = c("unblind", "blind"), lam.type = c("base", "pool"), interval = c(log(0.5),log(1))) {
  frailty.type = match.arg(frailty.type)
  lam.type = match.arg(lam.type)

  fdt = data.frame(N, tau, lam, thta, alp, pow, ar, frailty.type, lam.type)
  rn = 1:nrow(fdt)
  fdt$rn = rn

  oout = vapply(split(fdt, list(rn)), function(x) {
    N=x$N; tau=x$tau; lam=x$lam; thta=x$thta; alp=x$alp; pow=x$pow; frailty.type = x$frailty.type; lam.type = x$lam.type; ar = x$ar
    zerofunc_solvebta1 = function(bta1) {
      zerofunc(N=x$N, bta1 = bta1, thta=x$thta, tau=x$tau, lam=x$lam, alp=x$alp, pow=x$pow, ar = x$ar, frailty.type = x$frailty.type, lam.type = x$lam.type)
    }
    y = stats::uniroot(f = zerofunc_solvebta1, interval = interval, extendInt = "yes")$root
    return(y)
  }, numeric(1), USE.NAMES = FALSE)

  return(oout)
}

zerofunc = function(N, bta1, thta, tau, lam, alp, pow, ar, frailty.type, lam.type) {

  if (frailty.type == "blind") {
    thta_ = thtap2thta(bta1 = bta1, thtap = thta, ar = ar)
  } else if (frailty.type == "unblind") {
    thta_ = thta
  }
  if (lam.type == "pool") {
    lam0_ = lam / (ar + (1-ar) * exp(bta1))
  } else if (lam.type == "base") {
    lam0_ = lam
  }

  oout = pow(N = N, bta1 = bta1, thta = thta_, tau = tau, lam0 = lam0_, alp = alp, ar = ar) - pow

  return(oout)
}


