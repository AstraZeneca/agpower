#' Power calculations for one-sided two-sample test of LWYY model parameter.
#'
#' Function to compute power at one-sided Type I control level `alp/2` or
#' determine parameters to target given power to test the hypotheses
#' `H0: RR = 1` vs `HA: RR < 1` for the LWYY (Andersen-Gill with robust standard
#' errors) model.
#'
#' @param N Sample size.
#' @param RR,bta1 Rate ratio, log-transform of rate ratio. Provide at most one of RR and bta1.
#' @param thta Variance of frailty parameter.
#' @param L Number of events.
#' @param tau Expected follow-up time.
#' @param lam Event rate. If lam.type = "base", it is the event rate for control. If lam.type = "pool", it is the pooled rate.
#' @param alp Two-sided alpha-level.
#' @param pow Target power.
#' @param ar Allocation ratio (Number control / Total)
#' @param frailty.type Indicates whether frailty variance is based on blinded information ("blind") or unblinded ("unblind"). Default "unblind".
#' @param lam.type Indicates whether event rate is based on control rate ("base") or pooled rate ("pool"). Ignored if lam = NULL. Default "base".
#' @return An object of class "power.lwyytest", a list of arguments, including those computed, with `method` and note `elements`.
#' @details
#' An object of class "power.lwyytest" has the following components:
#'
#' \itemize{
#'   \item `method`: Description of the test.
#'   \item `N`: Total sample size. `N = Nc + Nt`, where `Nc = ar * N` is the number on control and `Nt = (1-ar) * N` is the number on treatment.
#'   \item `RR`, `RRCV`: Rate ratio powered for and critical value for rate ratio.
#'   \item `bta1`, `bta1CV`: The log of `RR`, `RRCV`.
#'   \item `thta`: The variance of the frailty parameter.
#'   \item `thtap`: If frailty.type = "blind", the input variance from a blinded source.
#'   \item `L`: Number of events. Note that all else being equal, an increase (decrease) in `N` requires a decrease (increase) in `L`.
#'   \item `tau`, `lam0`: The mean follow-up time and control event rate. Note if `tau = 1`, then `lam0` may be interpreted as mean cumulative intensity on control.
#'   \item `lamp`: If lam.type = "pool" and argument `lam` is not null, the input rate from a blinded source.
#'   \item `alp`: `alp/2` is the one-sided Type I control level of the test `H0: RR = 1` vs `HA: RR < 1`.
#'   \item `pow`: Power of the test.
#'   \item `ar`: Allocation to control ratio `Nc/N`.
#'   \item `note`: Set of key notes about the assessment.
#' }
#'
#'
#'
#'
#' @examples
#'
#' x = power.lwyy.test(N = 1000, RR = 0.8, thta = 1, L = 1000, tau = 0.9, alp = 0.05, ar = 0.5)
#' print(x)
#'
#' x = power.lwyy.test(N = 1000, RR = 0.8, thta = 1, tau = 0.9, lam = 1.23, alp = 0.05, ar = 0.5)
#' print(x)
#'
#' x = power.lwyy.test(N = 1000, RR = 0.8, thta = 1, tau = NULL, lam = 1.23, alp = 0.05, ar = 0.5)
#' print(x, digits = 3)
#'
#' @export
power.lwyy.test = function(N = NULL, RR = NULL, bta1 = NULL, thta, L = NULL, tau = NULL, lam = NULL, alp = 0.05, pow = NULL, ar = 0.5, frailty.type = c("unblind", "blind"), lam.type = c("base", "pool")) {
  frailty.type = match.arg(frailty.type)
  lam.type = match.arg(lam.type)

  note = c(
    "N is the total sample size N = Nc + Nt, where Nc = ar * N, Nt = (1-ar) * N.",
    "Test of H0: RR = 1 vs HA: RR < 1. alp/2 is one-sided alpha level."
  )

  # Check bta1/RR
  if (sum(vapply(list(bta1, RR), is.null, NA)) == 0) {
    stop("Provide at most one of (RR, bta1). The other is derived through RR = exp(bta1)")
  } else if (is.null(RR) & !is.null(bta1)) {
    RR = exp(bta1)
  } else if (!is.null(RR) & is.null(bta1)) {
    bta1 = log(RR)
  }

  if (sum(vapply(list(N, bta1, L, lam), is.null, NA)) %in% c(0, 3, 4)) {
    stop("At least one, and no more than two, of (N, bta1, L, lam) must be null.")
  }
  if (sum(vapply(list(N, bta1, L, lam), is.null, NA)) == 1) {
    if (sum(vapply(list(alp, pow), is.null, NA)) %in% c(0, 2)) {
      stop("Three of  (N, bta1, L, lam) were provided. Therefore exactly one of (alp, pow) must be set to null.")
    }
  }
  if (sum(vapply(list(N, bta1, L, lam), is.null, NA)) == 2) {
    if (sum(vapply(list(alp, pow), is.null, NA)) != 0) {
      stop("Two of  (N, bta1, L, lam) were provided. Therefore exactly none of (alp, pow) must be set to null.")
    }
  }

  # Check input values
  if (!is.null(N)) {
    if (!is.numeric(N)) stop("N must be numeric.")
    if (any(N <= 0)) stop("N must be positive.")
  }
  if (!is.null(RR)) {
    if (!is.numeric(RR)) stop("RR and bta1 must be numeric.")
    if (any(RR <= 0)) stop("RR = exp(bta1) must be positive.")
    if (any(RR >= 1)) stop("At least one RR > 1 (bta1 > 0). Should be RR < 1 (bta1 < 0).")
  }
  if (!is.null(thta)) {
    if (!is.numeric(thta)) stop("thta must be numeric.")
    if (any(thta < 0)) stop("thta must be non-negative.")
  }
  if (!is.null(L)) {
    if (!is.numeric(L)) stop("L must be numeric.")
    if (any(L <= 0)) stop("L must be positive.")
  }
  if (!is.null(tau)) {
    if (!is.numeric(tau)) stop("tau must be numeric.")
    if (any(tau <= 0)) stop("tau must be positive.")
  } else {
    tau = 1
    note = c(note, ("Mean follow-up time not provided. Setting tau = 1.\nlam0 may be interpreted as mean cumulative intensity on control."))
  }
  if (!is.null(lam)) {
    if (!is.numeric(lam)) stop("lam must be numeric.")
    if (any(lam <= 0)) stop("lam must be positive.")
  }
  if (!is.null(alp)) {
    if (!is.numeric(alp)) stop("alp must be numeric.")
    if (any(alp < 0 | alp > 1)) stop("alp must fall in the unit interval [0,1].")
  }
  if (!is.null(pow)) {
    if (!is.numeric(pow)) stop("alp must be numeric.")
    if (any(pow < 0 | pow > 1)) stop("pow must fall in the unit interval [0,1].")
  }
  if (!is.null(ar)) {
    if (!is.numeric(ar)) stop("alp must be numeric.")
    if (any(ar <= 0 | ar >= 1)) stop("ar must fall in the open unit interval (0,1).")
  }

  # Calculate missing values
  if (sum(vapply(list(N, bta1, L, lam), is.null, NA)) == 1) {
    if (!is.null(lam) & lam.type == "pool") {
      lamp = lam
      lam0 = lamp / (ar + (1-ar) * exp(bta1))
    } else {
      lamp = NULL
      lam0 = lam
    }

    if (is.null(N)) {
      N = L / (tau * lam0 * (ar + (1-ar) * exp(bta1)))
    }
    if (is.null(bta1)) {
      bta1 = log((L/(N*tau*lam0) - ar) / (1 - ar))
    }
    if (is.null(L)) {
      L = N * tau * lam0 * (ar + (1-ar) * exp(bta1))
    }
    if (is.null(lam0)) {
      lam0 = L / (N * tau * (ar + (1-ar) * exp(bta1)))
    }

    if (frailty.type == "blind") {
      thtap = thta
      thta = thtap2thta(bta1, thta)
    } else {
      thtap = NULL
      thta = thta
    }

    if (is.null(alp)) {
      alp = alpNeeded2(N = N, bta1 = bta1, thta = thta, L = L, pow = pow, ar = ar)
    }
    if (is.null(pow)) {
      pow = pow2(N = N, bta1 = bta1, thta = thta, L = L, alp = alp, ar = ar)
    }
  }

  if (sum(vapply(list(N, bta1, L, lam), is.null, NA)) == 2) {
    if (is.null(bta1)) {
      if (is.null(N)) {
        bta1 = btaNeeded3(thta = thta, L = L, tau = tau, lam = lam, alp = alp, pow = pow, ar = ar, frailty.type = frailty.type, lam.type = lam.type)
        if (lam.type == "pool") {
          lamp = lam
          lam0 = lamp / (ar + (1-ar) * exp(bta1))
        } else {
          lamp = NULL
          lam0 = lam
        }
        N = L / ((ar + (1-ar) * exp(bta1)) * tau * (lam0))
      } else
      if (is.null(L)) {
        bta1 = btaNeeded(N = N, thta = thta, tau = tau, lam = lam, alp = alp, pow = pow, ar = ar, frailty.type = frailty.type, lam.type = lam.type)
        if (lam.type == "pool") {
          lamp = lam
          lam0 = lamp / (ar + (1-ar) * exp(bta1))
        } else {
          lamp = NULL
          lam0 = lam
        }
        L = N * tau * lam0 * (ar + (1-ar) * exp(bta1))
      } else
      if (is.null(lam)) {
        bta1 = btaNeeded2(N = N, thta = thta, L = L, alp = alp, pow = pow, ar = ar, frailty.type = frailty.type)
        lamp = NULL
        lam0 = L / (N * tau * (ar + (1-ar) * exp(bta1)))
      }

      if (frailty.type == "blind") {
        thtap = thta
        thta = thtap2thta(bta1, thta)
      } else {
        thtap = NULL
        thta = thta
      }
    } else {
      if (frailty.type == "blind") {
        thtap = thta
        thta = thtap2thta(bta1, thta)
      } else {
        thtap = NULL
        thta = thta
      }
      if (!is.null(lam) & lam.type == "pool") {
        lamp = lam
        lam0 = lamp / (ar + (1-ar) * exp(bta1))
      } else {
        lamp = NULL
        lam0 = lam
      }

      if (is.null(N) & is.null(L)) {
        L = eventsNeeded(bta1 = bta1, thta = thta, tau = tau, lam0 = lam0, alp = alp, pow = pow, ar = ar)
        N = L / ((ar + (1-ar) * exp(bta1)) * tau * (lam0))
      } else
      if (is.null(N) & is.null(lam)) {
        N = nNeeded2(bta1 = bta1, thta = thta, L = L, alp = alp, pow = pow, ar = ar)
        lam0 = L / (N * tau * (ar + (1-ar) * exp(bta1)))
      } else
      if (is.null(L) & is.null(lam)) {
        L = eventsNeeded2(N = N, bta1 = bta1, thta = thta, alp = alp, pow = pow, ar = ar)
        lam0 = L / (N * tau * (ar + (1-ar) * exp(bta1)))
      }
    }
  }

  if (is.null(RR)) {
    RR = exp(bta1)
  }

  # Derive critical values
  bta1_critval = RR_critval = NA
  if (!any(is.nan(L) | is.infinite(L)) & !any(is.nan(N) | is.infinite(N))) {
    bta1_critval = btaNeeded2(N = N, thta = thta, L = L, alp = alp, pow = 0.5, ar = ar)
    RR_critval = exp(bta1_critval)
  }

  structure(
    list(
      method = "Two sample LWYY rate ratio power calculation.",
      N = N,
      RR = RR, RRCV = RR_critval,
      bta1 = bta1, bta1CV = bta1_critval,
      thta = thta, thtap = thtap,
      L = L, tau = tau,
      lam0 = lam0, lamp = lamp,
      alp = alp, pow = pow, ar = ar,
      note = note
    ),
    class = "power.lwyytest"
  )



}

#' @export
print.power.lwyytest = function(x, digits = getOption("digits"), exclude = c("bta1", "bta1CV"), ...) {
  x <- x[!sapply(x,is.null)]
  cat("\n    ", x$method, "\n\n")
  note = x$note
  x[c("method", "note")] = NULL
  y = x
  y[exclude] = NULL
  cat(paste(format(names(y), width = 15L, justify = "right"),
            format(y, digits = digits), sep = " = "), sep = "\n")
  if (!is.null(note))
    cat("\n", "NOTES: ", paste0("\n - ", note), "\n\n", sep = "")
  else cat("\n")
  invisible(x)

}


#' @export
as.data.frame.power.lwyytest = function(x, ...) {
  method = x$method
  note = x$note

  x[c("method", "note")] = NULL
  x = x[!sapply(x,is.null)]

  class(x) = "list"
  x = as.data.frame(x, ...)
}



