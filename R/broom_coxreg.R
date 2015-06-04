#' Tidy coxreg
#'
#' Tidy method for coxme objects
#'
#' @param x coxreg object
#' @param exponentiate boolean
#' @param ... other params
#' @export
#' @import eha
#' @examples
#' 
#' require(eha)
#' require(broom)
#' data(oldmort)
#' fit <- coxreg(Surv(enter, exit, event) ~ civ + sex + birthdate, data = oldmort)
#' tidy(fit)

tidy.coxreg <- function(x, exponentiate = TRUE, ...){
  beta <- x$coefficients
  nvar <- length(beta)
  nn <- c("estimate", "exp()", "std.error", "statistic", "p.value")
  se <- sqrt(diag(x$var)[1:nvar])
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

#' Glance coxreg
#'
#' Glance method for coxreg objects
#'
#' @param x coxme object
#' @param ... other params
#' @export

glance.coxreg <- function(x, ...){
  df <- ifelse(is.null(x$df), 
    sum(!is.na(x$coef)),
    x$df
  )
  logtest <- -2 * (x$loglik[1] - x$loglik[2])
  ret <- list(
    n = x$n, 
    nevent = x$events, 
    ttr = x$ttr,
    df = df,
    max.loglik = x$loklik[2],
    lt.test = logtest,
    concordance = x$concordance[1], 
    std.error.concordance = x$concordance[2],
    wald.p = 1 - pchisq(logtest, df)
  )
  ret <- as.data.frame(compact(ret))
  finish_glance(ret, x)
}


