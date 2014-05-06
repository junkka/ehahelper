# Event history analysis helper package

Helper functions for event history analysis with the survival package.

## Functions

### `ggsurv`

Make a data.frame of a survfit or coxph object for visualization with ggplot2.

**Example:**

```{r}
surv_object <- coxph(Surv(time, status) ~ strata(x), data = aml)
ggplot(ggsurv(surv_object), aes(time, surv, color=strata)) + geom_step()
```

### `split_rows`

Split rows at a split point if the split point falls between start and stop time.

**Example:**

```{r}
split_rows(mgus1, 500, "id", "start","stop")
```

### `age_split`

Split event history data into a equal sequences such as age groups. A row is split if the star and stop time stretches over a split point. Each row is assigned a new start and end time, a group and a duration within that group.

**Example:**

```{r}
mgus1$status <- ifelse(mgus1$event == "death", 1, 0)
x <- age_split(mgus1, "start","stop","id","status",500)
```

### `merge_reg`

Merge results from nested coxph regressions

### `cox_table`

Create a simple results table from a coxph regression.
