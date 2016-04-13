#' cox.zph as data.frame
#'
#' Extract scaled Schoenfeld residuals by variable for a cox.zph object
#' @param x a cox.zph object
#' @param val index of variables to extract
#' @export
#' @examples
#' library(survival)
#' surv_object <- coxph(Surv(time, status) ~ x, data = aml)
#' fit_zph <- cox.zph(surv_object, transform = "identity")
#' zph_df <- as.data.frame(fit_zph)
#'


as.data.frame.cox.zph <- function(x, val = NA){
  if (all(is.na(val)))
    val <- 1:nrow(x$var)
  
  res <- plyr::ldply(val, function(a){
    data.frame(
      time = x[a]$x,
      beta_t = as.numeric(x[a]$y[ ,1]),
      var = as.character(attr(x[a]$y, "dimnames")[[2]]),
      p_value = round(x[a]$table[ ,3], 3),
      stringsAsFactors = FALSE
    )
  })
  res
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
#' surv_object <- coxph(Surv(time, status) ~ x, data = aml)
#' fit_zph <- cox.zph(surv_object, transform = "identity")
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