#' Obtain theta from pooled theta (model without treatment).
#'
#' Function transforms estimate of blinded pooled theta (from misspecified model without adjusting for treatment) to estimate of theta for model adjusting for treatment.
#'
#' @param thtap Estimate of the variance of the frailty parameter under misspecification (no adjustment for treatment).
#' @param bta1 log-transform of assumed rate ratio.
#' @param ar Allocation ratio (Number in control / Total).
#' @param thtawarning If TRUE indicates how many estimates of theta were negative before setting to 0. Default FALSE.
#' @return An estimate of the variance of the frailty parameter when the model includes treatment.
#' @details
#' This function assumes a recurrent event distribution with exponential baseline function and frailty parameter distribution with mean 1 and variance \eqn{\theta},
#' and derives from expectations of the variance under misspecification.
#'
#' Specifically, the following relationship between the frailty variance for the misspecified model (\eqn{\theta_p}) and correctly specified model (\eqn{\theta}) is used
#' \deqn{\theta_p = Var(Z \cdot \exp{\beta}) = c^2 - 2pc - 2(1-p)\exp{\beta}c + (p + (1-p)\exp{2\beta})(\theta + 1)},
#' where \eqn{c = E(Z\exp{X\beta})} and \eqn{p = ar}.
#'
#' @examples
#'
#' thtap2thta(2, log(c(0.8, 0.7, 0.6)))
#'
#' @export
thtap2thta = function(bta1, thtap, ar = 0.5, thtawarning = FALSE) {

  c = ar + (1-ar)*exp(bta1)

  thta = (thtap + 1) * c^2 / (ar + (1-ar)*exp(2*bta1)) - 1

  thta[which(is.infinite(bta1))] = 0
  thta[which(is.nan(thta))] = 0

  if (sum(thta < 0) > 0) {
    nlz = sum(thta < 0)
    thta = pmax(thta, 0)
    if (thtawarning) warning(paste0(nlz, " out of ", length(thta), " estimates of theta were negative, setting these to zero."))
  }

  return(thta)

}
