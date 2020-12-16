#' Predict relative risk
#' 
#' predict centered relative risk from coxph and coxme objects
#' 
#' @param model coxph or coxme object
#' @param ... others
#' @export


rr_pred <- function(model, ...) UseMethod("rr_pred", model)


#' Predict relative risk
#' 
#' predict centered relative risk from coxph object
#' 
#' @param model coxph object
#' @param newdata data.frame 
#' @param center center variable
#' @param conf.level The confidence level to use for the confidence interval.
#' @export
#' @importFrom stringr str_detect
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
#' @examples 
#' 
#' library(survival)
#' fit <- coxph(Surv(time, status) ~ age + ph.ecog + strata(inst), lung)
#' newdata <- expand.grid(age = 50, ph.ecog = 0:3, inst = 11)
#' 
#' rr_pred(fit, newdata, center = c(ph.ecog = 0))
#' 

rr_pred.coxph <- function(model, newdata, center, conf.level = .95){
  
  
  if(conf.level >= 1 | conf.level <= 0)
    stop("conf.level must be between 0 and 1")
  
  terms <- Terms <- delete.response(terms(model))
  
  mm <- stats::model.frame(terms,newdata, na.action = na.pass, xlev = model$xlevels)
  mm2 <- model.matrix(model, mm)
  mm3 <- mm2
  # center
  mm3[,str_detect(colnames(mm2), names(center)) ] <- center
  mm_c <- mm2 - mm3
  
  # extract coefficients
  coef <- coef(model)
  # calculate predicted
  lp <- mm_c %*% coef 
  # approximate standard error
  se <- sqrt(diag(mm_c %*% vcov(model) %*% t(mm_c)))
  
  alpha <- 1-conf.level
  crit <- -qnorm(alpha/2)
  
  nd <- as_tibble(newdata) 
  mutate(nd,
      lp = as.vector(lp),
      se = se,
      rr = exp(lp),
      rr_l = exp(lp-(crit*se)),
      rr_h = exp(lp+(crit*se))
    )
}



#' Predict relative risk
#' 
#' predict centered relative risk from  coxme object
#' 
#' @param model coxme object
#' @param newdata data.frame 
#' @param center center variable
#' @param conf.level The confidence level to use for the confidence interval.
#' @export
#' @importFrom stringr str_detect
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
#' @examples 
#' 
#' library(coxme)
#' fit <- coxme(Surv(time, status) ~ age + ph.ecog + (1 | inst), data = lung)
#' newdata <- expand.grid(age = 50, ph.ecog = 0:3, inst = 11)
#' rr_pred(fit, newdata, center = c(ph.ecog = 0))

rr_pred.coxme <- function(model, 
                     newdata = NULL, 
                     center, 
                     conf.level = .95){
  
  if(conf.level >= 1 | conf.level <= 0)
    stop("conf.level must be between 0 and 1")
  
  
  n <- model$n[2]
  Terms <- delete.response(terms(model))
  has_strata <- !is.null(attr(Terms, "specials")$strata) 
  
  if (has_strata) 
    has_strata <- ifelse(length(attr(Terms, "specials")$strata) == 0, FALSE, has_strata)
  has_newdata  <- !is.null(newdata)
  
  
  # if strata update terms
  if (has_strata){
    strata_terms <- untangle.specials(Terms, "strata")
    Terms2 <- Terms[-strata_terms$terms]
  } else {
    Terms2 <- Terms
  }
  
  # Extract model.matrix
  mf <- survival:::model.frame.coxph(model)
  m <- model.frame(Terms, newdata)
  mm <- model.matrix(Terms2, m)
  mm <- mm[ ,-1]
  
  # Center mm
  mm3 <- mm
  mm3[ ,str_detect(colnames(mm3), names(center))] <- center
  mm_c <- mm - mm3
  
  # Extract coefficients
  coef <- fixed.effects(model)
  
  # centered predicted
  if (length(coef) == 1){
    lp <- mm_c * coef
  } else {
    lp <- (mm_c %*% coef)
  }
  
  # Approximate standard error
  se <- sqrt(diag(mm_c %*% vcov(model) %*% t(mm_c)))

  alpha <- 1-conf.level
  crit <- -qnorm(alpha/2)
  
  nd <- as_tibble(newdata) 
  mutate(
    nd,
    lp = as.vector(lp),
    se = se,
    rr = exp(lp),
    rr_l = exp(lp-(crit*se)),
    rr_h = exp(lp+(crit*se))
  )
}
