#' Survfit for ggplot2
#'
#' ggsurv creates a data.frame of a coxph or a survfit object for 
#'   visualization with ggplot2
#' @param x a coxph or survfit object
#' @param lims not implemented
#' @importFrom purrr map_dfr
#' @importFrom dplyr bind_rows
#' @export
#' @return data.frame with time, strata, std_err, upper, lower, n_risk, n_event, n_centor, 
#'   and surv. Optionaly if object is coxph, it returns cumhaz and if the object is 
#' stratified, strata.
#' @author Johan Junkka
#' @examples
#' library(survival)
#' library(ggplot2)
#' surv_object <- coxph(Surv(time, status) ~ strata(x), data = aml)
#' ggplot(ggsurv(surv_object), aes(time, surv, color=strata)) + geom_step()
#'
ggsurv <- function(x, lims = NULL) {
  if (!inherits(x, c("coxph","survfit"))) stop("x is not a coxph or survfit object")
  if (inherits(x, "coxph")) x <- survival::survfit(x)

  t_lim <- range(x$time)
  
  dat <- data.frame(
    time      = x$time,
    std_err   = x$std.err,
    upper     = x$upper,
    lower     = x$lower,
    n_risk    = x$n.risk,
    n_event   = x$n.event,
    n_censor  = x$n.censor,
    surv      = x$surv
  )
  
  
  # If coxph add cumhaz
  if(!is.null(x$cumhaz)) dat$cumhaz <- x$cumhaz
  
  # for each in x$strata rep and add to vector
  if (!is.null(x$strata)) {
    
    strata <- c()
    for (i in 1:length(x$strata)) {
      # Extract only label of factor
      s_name <- gsub("^.*?=","",names(x$strata[i]))
      # Add to strata vector
      strata <- c(strata, rep(s_name, as.numeric(x$strata[i])))
      
    }
    dat$strata <- as.factor(strata)
    
    f <- function(st){
      data.frame(
        time      = 0,
        std_err   = 0,
        upper     = 0,
        lower     = 0,
        n_risk    = max(dat$n_risk[dat$strata == st]),
        n_event   = 0,
        n_censor  = 0,
        surv      = 1,
        strata    = as.character(st),
        cumhaz    = 0
      )
    }
    
    dat <- bind_rows(dat, map_dfr(names(x$strata), f) )
    
    return(dat)
  } else {
    dat <- bind_rows(dat, data.frame(
      time      = 0,
      std_err   = 0,
      upper     = 0,
      lower     = 0,
      n_risk    = max(dat$n_risk),
      n_event   = 0,
      n_censor  = 0,
      surv      = 1,
      cumhaz    = 0
    ))
  }
  
  return(dat)
}
