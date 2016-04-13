#' Tidy coxme
#'
#' Tidy method for coxme objects
#'
#' @param x coxme object
#' @param exponentiate boolean
#' @param ... other params
#' @export
#' @import coxme
#' @import broom

tidy.coxme <- function(x, exponentiate = TRUE, ...){
  beta <- x$coefficients
  nvar <- length(beta)
  nfrail <- nrow(x$var) - nvar
  nn <- c("estimate", "exp()", "std.error", "statistic", "p.value")
  se <- sqrt(diag(x$var)[nfrail + 1:nvar])
  z <- qnorm((1 + 0.95)/2, 0, 1)
  ret <- data.frame(
    "term" = names(beta),
    "estimate" = beta,
    "std.error" = se,
    "statistic" = round(beta/se, 2),
    "p.value" = round(signif(1 - pchisq((beta/se)^2, 1), 2), 3),
    "conf.low" =  beta - z * se,
    "conf.high" =  beta + z * se
  )
  if (exponentiate) {
    ret$estimate <- exp(ret$estimate)
    ret$conf.low <- exp(ret$conf.low)
    ret$conf.high <- exp(ret$conf.high)
  }
  rownames(ret) <- c(1:nrow(ret))
  ret
}

#' Glance coxme
#'
#' Glance method for coxme objects
#'
#' @param x coxme object
#' @param ... other params
#' @export


glance.coxme <- function(x, ...){
  loglik <- x$loglik + c(0, 0, x$penalty)
  chi1 <- 2 * diff(x$loglik[c(1, 2)])
  chi1 <- 2 * diff(loglik[1:2])
  chi2 <- 2 * diff(loglik[c(1, 3)])
  ret <- rbind(c(round(chi1, 2), round(x$df[1], 2), signif(1 -
                                                             pchisq(chi1, x$df[1]), 5), round(chi1 - 2 * x$df[1],
                                                                                              2), round(chi1 - log(x$n[1]) * x$df[1], 2)), c(round(chi2,
                                                                                                                                                   2), round(x$df[2], 2), signif(1 - pchisq(chi2, x$df[2]),
                                                                                                                                                                                 5), round(chi2 - 2 * x$df[2], 2), round(chi2 - log(x$n[1]) *
                                                                                                                                                                                                                           x$df[2], 2)))
  dimnames(ret) <- list(c("Integrated loglik", " Penalized loglik"),
                        c("Chisq", "df", "p", "AIC", "BIC"))
  # finish_glance(ret, x)
  ret
}