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
