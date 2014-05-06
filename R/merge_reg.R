#' Merge Cox Regessions to One Table
#' 
#' @param full A model with all relavent variables or cector 
#'   f covariate names
#' @param list_res a list with all regressions to be merged
#' @export
merge_reg <- function(full, list_res){
  # If full is coxph extract to covariates
  if (inherits(full, "coxph")) {
    # Set a maximum number of covariates from largest model
    df <- data.frame(n = c(1:length(full$coefficients)))
    df$vars <- rownames(cbind(full$coefficients))
  } else {
    df <- data.frame(n = c(1:length(full)))
    df$vars <- full
  }
  
  n <- nrow(df)
  df <- rbind(df, c(n+1, "Test stat"),c(n+2, "LRT"), c(n+3,"Logrank"), c(n+4, "p"), c(n+5, "n"), c(n+16, "events"))
  df$n <- as.numeric(df$n)
  # for each model add column with beta and p
  for(i in 1:length(list_res)){
    beta <- coef(list_res[[i]])
    df2 <- data.frame(rownames(cbind(list_res[[i]]$coefficients)))
    se   <- sqrt(diag(list_res[[i]]$var))
    df2$se   <- as.numeric(formatC(exp(beta), format="f", digits=3))
    df2$p    <- 1 - pchisq((beta/se)^2, 1)
    df2$p    <- as.numeric(formatC(df2$p, format="f", digits=3))
    
    # add significans if p < 0.05
    df2$se <- ifelse(df2$p <= 0.05, paste(df2$se, "*"), ifelse(df2$p <= 0.1, paste(df2$se, "."), df2$se ) )
    df2$p <- NULL
    
    colnames(df2) <- c("vars", paste(names(list_res)[i]))
    df2$vars <- as.character(df2$vars)
    # Add test stat
    df2 <- rbind(df2, 
                 c("LRT", as.numeric( formatC( (list_res[[i]]$loglik[2]-list_res[[i]]$loglik[1])*2), format="f", digits=1)),
                 c("Logrank", as.numeric( formatC( list_res[[i]]$score), format="f", digits=1)),
                 c("n", as.numeric( formatC( list_res[[i]]$n), format="f", digits=1)),
                 c("events", as.numeric( formatC( list_res[[i]]$nevent), format="f", digits=1)))
    
    
    df <- merge(df, df2, by="vars", all.x=T, sort=F)    
  }
  
  return(df[order(df$n),c(1,3:ncol(df))])
}

