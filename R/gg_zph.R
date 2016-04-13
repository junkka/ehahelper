#' cox.zph as data.frame
#'
#' Extract scaled Schoenfeld residuals by variable for a cox.zph object
#' @param x a cox.zph object
#' @param val index of variables to extract
#' @export
#' @examples
#' library(survival)
#' fit <- coxph(Surv(time, status) ~ age + sex + disease, data = kidney)
#' fit_zph <- cox.zph(fit, transform = "identity")
#' zph_df <- as.data.frame(fit_zph)
#'


as.data.frame.cox.zph <- function(x, val = NA){
  if (all(is.na(val)))
    val <- 1:nrow(x$var)
  
  time <- c()
  beta_t <- numeric()
  var <- character()
  p_value <- numeric()
  
  for (i in val){
    time = c(time, x[i]$x)
    beta_t = c(beta_t, as.numeric(x[i]$y[ ,1]))
    var = c(var, rep(as.character(attr(x$y, "dimnames")[[2]][i]), length(x[i]$x)))
    p_value = c(p_value, rep(round(x[i]$table[ ,3], 3), length(x[i]$x)))
  }
  data.frame(time, beta_t, var, p_value,
    stringsAsFactors = FALSE
  )
}

#' ggplot cox.zph
#'
#' A ggplot2 implementation of plotting scaled Schoenfeld residuals
#'   from a cox.zph objects
#'
#' @param x a cox.zph object
#' @param val index of variables to plot
#' @param log logical, if TRUE log transform time
#' @param alpha geom_point alpha value
#' @import ggplot2
#' @export
#' @examples
#' library(survival)
#' fit <- coxph(Surv(time, status) ~ age + sex + disease, data = kidney)
#' fit_zph <- cox.zph(fit, transform = "identity")
#' gg_zph(fit_zph)
#' gg_zph(fit_zph, log = TRUE)
#'

gg_zph <- function(x, val = NA, log = FALSE, alpha = 0.2){
  x <- as.data.frame(x, val)
  
  x_lab = paste(x$transform, "time")
  
  x$var <- paste0(x$var, " (", x$p_value, ")")
  p <- ggplot(x, aes(time, beta_t)) + geom_point(alpha = alpha) +
    stat_smooth() +
    xlab(x_lab)
  
  if(log)
    p <- p + scale_x_log10()
  
  if (all(is.na(val)) | length(val) > 1)
    p <- p +facet_wrap(~var, scales = "free")
  else
    p <- p + labs(title = x$var)
  p
}