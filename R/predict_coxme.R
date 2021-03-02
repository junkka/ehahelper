#' Predict coxme 
#' 
#' Prediction method for coxme objects
#' 
#' @param object coxme object.
#' @param newdata data.frame new data set.
#' @param type type of prediction, linear predictor ("lp") or relative risk ("risk").
#' @param se.fit if TRUE, pointwise standard errors are produced for the predictions.
#' @param strata_ref logical, use strata as reference.
#' @export
#' @import coxme
#' @examples 
#' library(coxme)
#' data(eortc)
#' fit <- coxme(Surv(y, uncens) ~ trt + (1|center), eortc)
#' predict_coxme(fit)
#' 

predict_coxme <- function(object, 
                          newdata = NULL, 
                          type = c("lp", "risk"), 
                          se.fit = FALSE,
                          strata_ref = TRUE){
  
  if (!inherits(object, 'coxme'))
    stop("Primary argument much be a coxme object")
  
  type <- match.arg(type)
  n <- object$n[2]
  Terms <- delete.response(terms(object))
  has_strata <- !is.null(attr(Terms, "specials")$strata) 
  if (has_strata) 
    has_strata <- ifelse(length(attr(Terms, "specials")$strata) == 0, FALSE, has_strata)
  has_newdata  <- !is.null(newdata)
  
  if (!se.fit & type == "lp" & !has_newdata) return(object$linear.predictor)
  
  coef <- fixed.effects(object)
  mf <- survival:::model.frame.coxph(object)
  
  # boot.ci
  
  
  
  if (has_newdata){
    m <- model.frame(Terms, newdata)
  } else {
    m <- mf
  }
  
  # if strata update terms
  if (has_strata){
    strata_terms <- untangle.specials(Terms, "strata")
    Terms2 <- Terms[-strata_terms$terms]
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
    
    oldstrat <- mf[[strata_terms$vars]]
    xmeans <- rowsum(x, oldstrat)/as.vector(table(oldstrat))
  }
  
  if (!has_newdata){
    # extract all cols in x which matches Terms
    mm <- model.matrix(Terms2, data =mf)[ ,-1]
    m <- mf
  }
  
  if (has_strata & strata_ref){
    newstrat <- m[[strata_terms$vars]]
    mm <- mm - xmeans[match(newstrat, row.names(xmeans)), colnames(mm)]
  } else {
    mm <- mm - rep(object$means, each = nrow(m))
  }
  
  # if (!has_newdata & !has_strata){
  #   pred <- object$linear.predictor
  # }
  if (length(coef) == 1){
    pred <- mm * coef
  } else {
    pred <- (mm %*% coef)
  }
  
  
  if (se.fit) se <- sqrt(rowSums((mm %*% vcov(object)) * mm))
  if (type == "risk"){
    pred <- exp(pred)
    if (se.fit) se <- se * sqrt(pred)
  }
  if (se.fit) list(fit = pred, se.fit = se)
  else pred
}
