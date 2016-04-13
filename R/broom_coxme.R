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
#' @examples 
#' library(broom)
#' library(coxme)
#' fit <- coxme(Surv(y, uncens) ~ trt + (1|center), eortc)
#' tidy(fit)
#' 

tidy.coxme <- function(x, exponentiate = FALSE, ...){
  beta <- x$coefficients
  nvar <- length(beta)
  nfrail <- nrow(x$var) - nvar
  nn <- c("estimate", "exp()", "std.error", "statistic", "p.value")
  se <- sqrt(diag(as.matrix(x$var))[nfrail + 1:nvar])
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
#' @examples
#' library(broom)
#' library(coxme)
#' fit <- coxme(Surv(y, uncens) ~ trt + (1|center), eortc)
#' glance(fit)
#' 


glance.coxme <- function(x, ...){
  loglik <- x$loglik + c(0, 0, x$penalty)
  chi1 <- 2 * diff(loglik[1:2]) 
  chi2 <- 2 * diff(loglik[c(1,3)])
  temp0 <- as.list(c(
    x$n[2], x$n[1],
    round(chi1, 2), round(x$df[1], 2),
    as.numeric(loglik[3]),
    signif(1 - pchisq(chi1, x$df[1]), 5),
    round(chi1 - 2 * x$df[1], 2),
    round(chi1 - log(x$n[1]) * x$df[1], 2)))
  names(temp0) <- c("n", "events" , "Chisq", "df", "logLik", "p", "AIC", "BIC")
  
  ## random effects
  random <- VarCorr(x)
  nrow <-  sapply(random, 
                  function(x) if (is.matrix(x)) nrow(x) else length(x))
  maxcol <-max(sapply(random,
                      function(x) if (is.matrix(x)) 1 + ncol(x) else 2))
  temp1 <- matrix(NA, nrow=sum(nrow), ncol=maxcol)
  indx <- 0
  for (term in  random) {
    if (is.matrix(term)) {
      k <- nrow(term)
      nc <- ncol(term)
      for (j in 1:k) {
        temp1[j + indx, 1] <- sqrt(term[j,j])
        temp1[j + indx, 2] <- term[j,j]
        if (nc > j) {
          indx2 <- (j + 1):nc
          temp1[j + indx, 1 + indx2] <- term[j, indx2]
        }
      }
    }
    else {
      k <- length(term)
      temp1[1:k + indx, 1] <- sqrt(term)
      temp1[1:k + indx, 2] <- term
    }
    indx <- indx + k
  }
  
  indx <- cumsum(c(1, nrow))   # starting row of each effect
  temp3 <- rep("", nrow(temp1))
  temp3[indx[-length(indx)]] <- names(random)
  xname <- unlist(lapply(random, 
                         function(x) if (is.matrix(x)) dimnames(x)[[1]] else names(x)))
 
  b <- as.vector(temp1)
  names(b) <- paste(rep(c("random_variance", "random_sd"), each = length(b)/2), temp3, sep = "_")
  
  # groups
  grps <- unlist(lapply(x$frail, length))
  names(grps) <- paste0("random_n_", names(grps))
  
  ret <- as.list(c(temp0, grps, b))
  ret <- as.data.frame(ret)
  broom:::finish_glance(ret, x)
}