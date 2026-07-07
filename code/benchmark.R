# Benchmark functional age against established mortality-risk measures (R1 #1).
#
# Every measure is fit as Cox(Surv ~ age + sex + measure) with 5-fold CV, folds
# stratified by age x sex, scored by held-out Harrell's C. Measures compared:
# the five functional-age domains, UK Biobank-refit PhenoAge and KDM, the
# Charlson comorbidity index and the Williams frailty index. Each is evaluated
# on its own complete-case cohort, so cohort sizes differ across rows. We also
# report a nested likelihood-ratio test against age + sex, and the correlation
# of each domain's functional age with PhenoAge and KDM.
#
# Run from the repository root:  Rscript code/benchmark.R

source("code/_config.R")
suppressMessages({
  library(tidyverse)
  library(survival)
  library(BioAge)
})

seed <- cfg$seed
domains <- c("cognitive", "psychological", "locomotion", "vitality", "sensory")

# PhenoAge / KDM inputs: the nine biomarkers, in SI units, keyed by field id.
biomarkers <- c("albumin_gL", "lymph", "mcv", "glucose_mmol",
                "rdw", "creat_umol", "lncrp", "alp", "wbc")
biomarker_fields <- c(albumin_gL = "30600", lymph = "30180", mcv = "30040",
                      glucose_mmol = "30740", rdw = "30070", creat_umol = "30700",
                      alp = "30610", wbc = "30000", age = "21003")

cidx <- function(time, status, lp) {
  concordance(Surv(time, status) ~ lp, reverse = TRUE)$concordance
}

make_folds <- function(age, sex, k = 5) {
  strat <- interaction(cut(age, 5), sex, drop = TRUE)
  set.seed(seed)
  fold <- integer(length(age))
  for (lv in levels(strat)) {
    idx <- sample(which(strat == lv))
    fold[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  fold
}

# Held-out C for age + sex + rhs, averaged over folds.
cv_cindex <- function(dat, rhs) {
  fold <- make_folds(dat$age, dat$sex)
  f <- as.formula(paste("Surv(time, status) ~ age + sex +", rhs))
  cs <- sapply(seq_len(5), function(i) {
    tr <- dat[fold != i, ]; te <- dat[fold == i, ]
    cidx(te$time, te$status, predict(coxph(f, tr), te, type = "lp"))
  })
  mean(cs)
}

# KDM is refit per fold (its own regression on age) and projected to held out.
cv_cindex_kdm <- function(dat) {
  fold <- make_folds(dat$age, dat$sex)
  cs <- sapply(seq_len(5), function(i) {
    tr <- dat[fold != i, ]; te <- dat[fold == i, ]
    fit <- kdm_calc(tr, biomarkers = biomarkers)
    tr$kdm <- fit$data$kdm
    te$kdm <- kdm_calc(te, biomarkers = biomarkers, fit = fit$fit)$data$kdm
    m <- coxph(Surv(time, status) ~ age + sex + kdm, tr)
    cidx(te$time, te$status, predict(m, te, type = "lp"))
  })
  mean(cs)
}

lr_test <- function(dat, rhs) {
  m0 <- coxph(Surv(time, status) ~ age + sex, dat)
  m1 <- coxph(as.formula(paste("Surv(time, status) ~ age + sex +", rhs)), dat)
  chisq <- 2 * (m1$loglik[2] - m0$loglik[2])
  df <- length(m1$coefficients) - length(m0$coefficients)
  tibble(chisq = chisq, df = df, p = pchisq(chisq, df, lower.tail = FALSE))
}

# ---- cohorts ----------------------------------------------------------------
survival <- read_rds(out_path("processed", "survival.rds")) %>%
  filter(!is.na(eid)) %>%
  mutate(eid = as.character(eid))

# Biomarker cohort for PhenoAge and KDM.
bio <- as.data.frame(read_rds(in_path("biochemistry")))
biomarker_data <- tibble(eid = rownames(bio))
for (nm in names(biomarker_fields)) biomarker_data[[nm]] <- bio[[biomarker_fields[[nm]]]]
biomarker_data$lncrp <- log1p(bio[["30710"]] / 10)
biomarker_data <- biomarker_data %>%
  inner_join(survival %>% select(eid, time, status, sex), by = "eid") %>%
  filter(complete.cases(across(all_of(c(biomarkers, "age", "time", "status", "sex")))))

# Out-of-fold functional age for a domain (the fitted mortality model's held-out
# risk), used for the discrimination row and the clock correlations.
domain_prediction <- function(d) {
  m <- read_rds(out_path("models", "clinical", paste0(d, ".rds")))$pred[[1]]
  tibble(eid = as.character(m$eid), risk = m$risk, age = m$age) %>%
    inner_join(survival %>% select(eid, time, status), by = "eid") %>%
    filter(complete.cases(.))
}

# One complete-case cohort per domain, from that domain's clinical features
# (used for the nested LR test).
domain_cohort <- function(d) {
  dom <- read_rds(out_path("processed", "domains", paste0(d, ".rds")))
  dom <- dom[, colSums(!is.na(dom)) > cfg$feature_cutoff, drop = FALSE]
  dom <- as.data.frame(dom)
  dom$eid <- as.character(rownames(dom))
  survival %>%
    select(eid, time, status, age, sex) %>%
    inner_join(dom, by = "eid") %>%
    filter(complete.cases(.))
}

# Derived comorbidity / frailty scores (one row per participant, see prep scripts).
charlson <- read_rds(out_path("processed", "charlson.rds")) %>%
  mutate(eid = as.character(eid)) %>%
  inner_join(survival %>% select(eid, time, status, age, sex), by = "eid") %>%
  filter(complete.cases(.))
frailty <- read_rds(out_path("processed", "frailty_index.rds")) %>%
  mutate(eid = as.character(eid)) %>%
  select(eid, frailty_index) %>%
  inner_join(survival %>% select(eid, time, status, age, sex), by = "eid") %>%
  filter(complete.cases(.))

# ---- C-index table ----------------------------------------------------------
rows <- list()
add_row <- function(label, dat, C) {
  rows[[label]] <<- tibble(measure = label, n = nrow(dat),
                           deaths = sum(dat$status), C = round(C, 3))
}

for (d in domains) {
  dat <- domain_prediction(d)
  add_row(paste0("Functional age: ", d), dat, cidx(dat$time, dat$status, dat$risk))
}
add_row("age + sex", biomarker_data, cv_cindex(biomarker_data, "1"))
add_row("PhenoAge (UKBB-refit)", biomarker_data,
        cv_cindex(biomarker_data, paste(biomarkers, collapse = " + ")))
add_row("KDM (UKBB-refit)", biomarker_data, cv_cindex_kdm(biomarker_data))
add_row("Charlson comorbidity index", charlson, cv_cindex(charlson, "charlson"))
add_row("Williams frailty index", frailty, cv_cindex(frailty, "frailty_index"))

cindex_table <- bind_rows(rows)
cat("\n=== Mortality discrimination (5-fold CV, held-out Harrell's C) ===\n")
print(as.data.frame(cindex_table), row.names = FALSE)
write_csv(cindex_table, out_path("benchmark", "cindex.csv"))

# ---- likelihood-ratio test vs age + sex -------------------------------------
lr <- bind_rows(
  lapply(domains, function(d) {
    dat <- domain_cohort(d)
    feat <- setdiff(colnames(dat), c("eid", "time", "status", "age", "sex"))
    lr_test(dat, paste(sprintf("`%s`", feat), collapse = " + ")) %>%
      mutate(measure = paste0("Functional age: ", d), .before = 1)
  })
) %>%
  bind_rows(
    lr_test(biomarker_data, paste(biomarkers, collapse = " + ")) %>%
      mutate(measure = "PhenoAge (UKBB-refit)", .before = 1),
    lr_test(charlson, "charlson") %>% mutate(measure = "Charlson comorbidity index", .before = 1),
    lr_test(frailty, "frailty_index") %>% mutate(measure = "Williams frailty index", .before = 1)
  )
cat("\n=== Nested LR test vs age + sex (full cohort) ===\n")
print(as.data.frame(lr), row.names = FALSE)
write_csv(lr, out_path("benchmark", "lrtest.csv"))
