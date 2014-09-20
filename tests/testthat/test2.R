context('coxph to long')

# This function will create a long format data.frame. 
# Each variable has one row with model, var, level, haz, ci_low, ci_high, p, signif
# | model  | var | level  |   haz | etc...
# |:-------|:----|:-------|------:|
# | model1 | sex | male   |    NA |
# | model1 | sex | female | 0.987 |
# | model1 | age | age    | 0.905 |
# | model2 | age | age    | 1.111 |


test_that('coxph_to_long', {
  library(survival)
  data('cancer')
  cancer$sex <- factor(cancer$sex, labels = c('male', 'female'))
  model1 <- coxph(Surv(time, status) ~ age + sex  + ph.ecog + strata(inst), data=cancer)
  model2 <- coxph(Surv(time, status) ~ age + strata(inst), data=cancer)
  res1 <- coxph_to_long(model1)
  res2 <- coxph_to_long(model1, model2)
  
  expect_that(dim(res1), equals(c(4, 8)))
  expect_that(dim(res2), equals(c(5, 8)))
  # expext_that(unique(res1$model), equals('model1'))
  # expect_that(unique(res1$var), is_equivalent_to(c('sex', 'age', 'ph.ecog')))

})