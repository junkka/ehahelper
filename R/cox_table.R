#' Make a Results Table From Cox Regression
#' 
#' @param cox A cox.ph regression
#' @export

cox_table <- function(cox)
{
  beta <- coef(cox)
  se   <- sqrt(diag(cox$var))
  p    <- 1 - pchisq((beta/se)^2, 1)
  CI   <- round(confint(cox), 3)
  
  # Bind columns together, and select desired rows
  res <- cbind(beta, se = exp(beta), CI, p)
  res <- formatC(res, format="f", digits=3)
  res <- data.frame(res)
  colnames(res) <- c("beta", "se","2,5%","97,5%","p")
  return(res)
}


#' Transform coxph objects to long format table
#'
#' For nested models
#'
#' @param ... coxph objects or list of objects
#' @param cilevel confidence interval
#' @author Johan Junkka \email{johan.junkka@@gmail.com}
#' @export
#' @examples
#' library(survival)
#' data('cancer', package = 'survival')
#' cancer$sex <- factor(cancer$sex, labels = c('male', 'female'))
#' model1 <- coxph(Surv(time, status) ~ age + sex  + ph.ecog + strata(inst), data=cancer)
#' model2 <- coxph(Surv(time, status) ~ age + strata(inst), data=cancer)
#' coxph_to_long(model1, model2)
#' 



coxph_to_long <- function(..., cilevel=0.95) {
  x <- list(...)
  if (any(as.logical(lapply(x, inherits, "list")))) {
    if (!inherits(x[[1]], 'list'))
      stop('No coxph or list of coxph')
    x <- x[[1]]
  } else  {
    names(x) <- obj_names(...)
  }

  dat <- ldply(x, .fun=function(y) coxph_df(y, cilevel))
  colnames(dat)[1] <- 'model'

  # get diff between unique model var and all var by model
  model_diff <- ddply(dat, .(model), .fun = function(a){
    var_names <- as.character(levels(dat$var))
    res <- length(levels(dat$var)) - sum(unique(as.character(a$var)) %in% var_names)
  })
  
  # get var where model = model with min V1
  full_model  <- model_diff$model[model_diff$V1 == min(model_diff$V1)][1]
  nm          <- unique(as.character(dat$levs[dat$model == full_model]))
  # Order levs by largest model
  dat$levs    <- ordered(dat$levs, levels=rev(nm))

  # Add summary stat by model
  attr(dat, "summary") <- llply(x, .fun = coxph_sum)
  class(dat)           <- c("cox_tb", class(dat))

  return(dat)
}

obj_names <- function(...) {
   x <- deparse(substitute(list(...)))
   x <- str_replace(str_replace(x, '^list\\(', ''), '\\)', '')
   x <- str_trim(unlist(str_split(x, ',')))
   return(x)
}

coxph_df <- function(x, cilevel = 0.9) {
  # by assign get xlevels
  ldply(names(x$assign), .fun = function(var) {

    s       <- summary(x)$coefficients
    beta    <- exp(coef(x))
    pvalue  <- s[ ,ncol(s)]
    loci    <- exp(confint(x, level = cilevel))[ ,1]
    hici    <- exp(confint(x, level = cilevel))[ ,2]
    n_coef  <- x$assign[[var]]

    if (str_detect(var, ':')){
      # is interaction
      vars <- str_split(var, ':')[[1]]
      lev1 <- x$xlevels[[vars[1]]]
      lev1 <- lev1[2:length(lev1)]
      lev2 <- x$xlevels[[vars[2]]]
      levs <- unlist(lapply(lev1, function(a) paste(a, lev2, sep = ' * ')))
    } else {
      levs <- x$xlevels[[var]]
      if (is.null(levs)) 
        levs <- var
    }
    if (length(levs) > 1) {
      data.frame(
        var     = rep(var, length(levs)),
        levs    = levs,
        haz     = c(1, beta[n_coef]),
        p       = c(NA, pvalue[n_coef]),
        ci_low  = c(NA, loci[n_coef]),
        ci_high = c(NA, hici[n_coef]),
        signif  = ifelse(c(NA, pvalue[n_coef]) < (1 - cilevel), TRUE, FALSE)
      )
    } else {
      data.frame(
        var     = rep(var, length(levs)),
        levs    = levs,
        haz     = beta[n_coef],
        p       = pvalue[n_coef],
        ci_low  = loci[n_coef],
        ci_high = hici[n_coef],
        signif  = ifelse(pvalue[n_coef] < (1 - cilevel), TRUE, FALSE)
      )
    }
  })
}

#' Summary stats for coxph object
#'
#' Description
#'
#' @param x coxph object
#' @export


coxph_sum <- function(x) {
  coef       <- x$coef
  df         <- ifelse(is.null(x$df), sum(!is.na(coef)), round(sum(x$df),2)) 
  getp <- function(x) 1 - pchisq(x, df)
  data.frame(
    var = c("Events", "Observations", "AIC", "Likelihood ratio test", "Logrank score", "Wald test"),
    value = c(
      as.numeric(x$nevent),
      as.numeric(x$n),
      extractAIC(x)[2],
      (-2 * (x$loglik[1] - x$loglik[2])),
      as.numeric(x$score),
      as.numeric(x$wald.test)
      ),
    p = c(rep(NA, 3), getp(-2 * (x$loglik[1] - x$loglik[2])), getp(x$score), getp(x$wald))
  )
}
