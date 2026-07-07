# Marginal (age- and sex-adjusted) per-feature hazard ratios: each domain feature
# fit on its own, Surv(time, status) ~ age + sex + feature, on its own complete-case
# sample. These are the correctly-signed univariate effect sizes (unlike the joint
# model, where summed scores compete with their own components). Plotted by
# code/supplementary_forest.R.
#
# Run from the repository root:  Rscript code/marginal_hr.R

source("code/_config.R")
suppressMessages({
  library(tidyverse)
  library(survival)
})

survival <- read_rds(out_path("processed", "survival.rds")) %>% na.omit %>% column_to_rownames("eid")

rows <- list()
for (d in c("cognitive", "psychological", "locomotion", "vitality", "sensory")) {
  dom <- read_rds(out_path("processed", "domains", paste0(d, ".rds")))
  dom <- dom[, colSums(!is.na(dom)) > cfg$feature_cutoff, drop = FALSE]
  dom <- dom[, !colnames(dom) %in% domain_feature_exclude, drop = FALSE]
  ids <- intersect(rownames(dom), rownames(survival))
  for (f in colnames(dom)) {
    df <- data.frame(survival[ids, c("time", "status", "age", "sex")], x = dom[ids, f])
    df <- df[complete.cases(df), ]
    if (nrow(df) < 1000 || sd(df$x) == 0) next
    sm <- summary(coxph(Surv(time, status) ~ age + sex + x, data = df))$coefficients
    b <- sm["x", "coef"]; se <- sm["x", "se(coef)"]; sdx <- sd(df$x)
    rows[[paste(d, f)]] <- tibble(domain = str_to_sentence(d), feature = f, n = nrow(df),
      HR_perSD = exp(b * sdx), HR_lo = exp((b - 1.96 * se) * sdx), HR_hi = exp((b + 1.96 * se) * sdx),
      std_beta = b * sdx, p = sm["x", "Pr(>|z|)"])
  }
}
write_csv(bind_rows(rows), out_path("processed", "marginal_hr.csv"))
