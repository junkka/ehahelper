Event history analysis helper package
=====================================

Helper functions for event history analysis with the survival package.

Install
-------

Install from GitHub repository.

``` r
library(devtools)
install_github('junkka/ehahelper')
```

Functions
---------

#### ggsurv

Make a data.frame of a `survfit` or `coxph` object for visualization with ggplot2.

``` r
library(ehahelper)
library(survival)
library(ggplot2)
surv_object <- coxph(Surv(time, status) ~ strata(x), data = aml)
ggplot(ggsurv(surv_object), aes(time, surv, color=strata)) + geom_step()
```

![](README_files/figure-markdown_github/ggsurv-1.png)

#### gg\_zph

ggplot2 impementation for visualizing scaled Schoenfield residuals form a cox.zph object.

``` r
bladder1 <- bladder[bladder$enum < 5, ] 
fit <- coxph(Surv(stop, event) ~ (rx + size + number) * strata(enum) + 
          cluster(id), bladder1)
x <- cox.zph(fit, transform = "identity")
gg_zph(x)
```

![](README_files/figure-markdown_github/gg_zph-1.png)

``` r
gg_zph(x, log = TRUE)
```

![](README_files/figure-markdown_github/gg_zph-2.png)
