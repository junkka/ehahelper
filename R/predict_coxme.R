#' Predict coxme 
#' 
#' Prediction method for coxme objects
#' 
#' @param object coxme object
#' @param newdata data.frame new data set
#' @param type type of prediction, linear predictor ("lp") or relative risk ("risk")
#' @param strata_ref logical, use strata as reference
#' @export
#' @import coxme
#' @examples 
#' library(coxme)
#' fit <- coxme(Surv(y, uncens) ~ trt + (1|center), eortc)
#' predict_coxme(fit)
#' 

predict_coxme <- function(object, 
                          newdata = NULL, 
                          type = c("lp", "risk"), 
                          strata_ref = TRUE){
  
  type <- match.arg(type)
  n <- object$n[2]
  Terms <- delete.response(terms(object))
  has_strata <- !is.null(attr(Terms, "specials")$strata) 
  if (has_strata) 
    has_strata <- ifelse(length(attr(Terms, "specials")$strata) == 0, FALSE, has_strata)
  has_newdata  <- !is.null(newdata)
  
  
  coef <- fixed.effects(object)
  mf <- stats::model.frame(object)
  if (has_newdata){
    m <- model.frame(Terms, newdata)
  } else {
    m <- mf
  }
  
  # if strata update terms
  if (has_strata){
    strata_terms <- untangle.specials(Terms, "strata")
    Terms2 <- Terms[-stemp$terms]
  } else {
    Terms2 <- Terms
  }
  
  if (has_newdata){
    mm <- model.matrix(Terms2, m)
    mm <- mm[ ,-1]
  }
  
  # has strata and reference is strata
  # calculate strata means
  if (has_strata & strata_ref){
    # Full model matrix
    x <- model.matrix(Terms, data = mf)
    
    oldstrat <- mf[[stemp$vars]]
    xmeans <- rowsum(x, oldstrat)/as.vector(table(oldstrat))
  }
  
  if (!has_newdata){
    # cols <- rownames(attr(Terms, "factor"))[-stemp$terms]
    # extract all cols in x which matches Terms
    mm <- model.matrix(Terms2, data =mf)[ ,-1]
    m <- mf
  }
  
  if (has_strata & strata_ref){
    newstrat <- m[[stemp$vars]]
    mm <- mm - xmeans[match(newstrat, row.names(xmeans)), colnames(mm)]
  } else {
    mm <- mm - rep(object$means, each = nrow(m))
  }
  
  if (length(coef) == 1){
    pred <- mm * coef
  } else {
    pred <- (mm %*% coef)
  }
  
  
  se <- sqrt(rowSums((mm %*% vcov(object)) * mm))
  if (type == "risk"){
    pred <- exp(pred)
    se <- se * sqrt(pred)
  }
  
  list(fit = pred, se.fit = se)
}
