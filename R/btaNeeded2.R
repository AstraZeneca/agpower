#' Function to compute log rate ratio needed (fixed sample size)
#'
#' Function to compute log rate ratio needed to achieve target power at one-sided Type I control level alp/2. Useful to compute critical value (set pow = 0.5).
#'
#' This function computes the log rate ratio bta1 as the root of the equation pow2(bta1, N, thta_, L, alp, ar) - pow = 0, where thta_ depends on bta1 if using estimates from blinded analyses.
#' If frailty.type = "blind": thta_ = thtap2thta(thta, bta1); otherwise if frailty = "unblind": thta_ = thta.
#' Function assumes a rate ratio < 1 is favourable to treatment.
#' @param N Sample size.
#' @param thta Variance of frailty parameter. If frailty.type = "blind", assumes thta derives from pooled model; if frailty.type = "unblind" assumes thta is from correctly specified model. Default "unblind".
#' @param L Number of events
#' @param alp Two-sided alpha-level.
#' @param pow Target power.
#' @param ar Allocation ratio (Number control / Total)
#' @param frailty.type Indicates whether frailty variance is based on blinded information ("blind") or unblinded ("unblind"). Default "unblind".
#' @param interval Initial search interval for bta1.
#' @return The log rate ratio.
#' @examples
#'
#' # Based on unblinded estimates
# btaNeeded2(N = 1000, thta = 2, L = 1000, alp = c(0.01, 0.05), pow = c(.5))
# exp(btaNeeded2(N = 1000, thta = 2, L = 1000, alp = c(0.01, 0.05), pow = c(.5)))
#'
#' # Based on blinded estimates
#' btaNeeded2(N = 1000, thta = 2, L = 1000, alp = c(0.01, 0.05), pow = c(.5), frailty.type = "bl")
#' exp(btaNeeded2(N = 1000, thta = 2, L = 1000, alp = c(0.01, 0.05), pow = c(.5), frailty.type = "bl"))
#'
#'
#' @export
btaNeeded2 = function(N, thta, L, alp = 0.05, pow = 0.8, ar = 0.5, frailty.type = c("unblind", "blind"), interval = c(log(0.5),log(1))) {
  frailty.type = match.arg(frailty.type)

  fdt = data.frame(N, L, thta, alp, pow, ar, frailty.type)
  rn = 1:nrow(fdt)
  fdt$rn = rn

  oout = vapply(split(fdt, list(rn)), function(x) {
    N=x$N; L=x$L; thta=x$thta; alp=x$alp; pow=x$pow; frailty.type = x$frailty.type; ar = x$ar
    zerofunc2_solvebta1 = function(bta1) {
      zerofunc2(N=x$N, bta1 = bta1, thta=x$thta, L=x$L, alp=x$alp, pow=x$pow, ar = x$ar, frailty.type = x$frailty.type)
    }
    y = stats::uniroot(f = zerofunc2_solvebta1, interval = interval, extendInt = "yes")$root
    return(y)
  }, numeric(1), USE.NAMES = FALSE)

  return(oout)
}

zerofunc2 = function(N, bta1, thta, L, alp, pow, ar, frailty.type) {

  if (frailty.type == "blind") {
    thta_ = thtap2thta(bta1 = bta1, thtap = thta, ar = ar)
  } else if (frailty.type == "unblind") {
    thta_ = thta
  }

  oout = pow2(N = N, bta1 = bta1, thta = thta_, L = L, alp = alp, ar = ar) - pow

  return(oout)
}


