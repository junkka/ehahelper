context('Merge reg')
# tests

test_that("merge_reg", {
  library(survival)
  data('cancer')
  fit_sml <- coxph(Surv(time, status) ~ age + strata(inst), data=cancer)
  fit_lrg <- coxph(Surv(time, status) ~ age + factor(sex)  + ph.ecog + strata(inst), data=cancer)
  res_tbl <- merge_reg(fit_lrg, list(small = fit_sml, larg=fit_lrg ))

  vars  <- as.character(c("age","factor(sex)2","ph.ecog","Test stat","LRT","Logrank","p","n","events"))
  small <- as.character(c("1.022 *",NA,NA,NA,"4.986","4.868",NA,"227","164"))
  large <- as.character(c("1.01","0.578 *","1.817 *",NA,"32.03","31.9",NA,"226","163"))
  exp_res <- data.frame(vars = vars, small=small, large=large, stringsAsFactors=FALSE)

  expect_that(res_tbl, is_equivalent_to(exp_res))

  covars    <- c("age","factor(sex)2","ph.ecog")
  res_tbl2  <- merge_reg(covars, list(small = fit_sml, larg=fit_lrg ))
  expect_that(res_tbl2, is_equivalent_to(exp_res))

})