# Per-domain Cox proportional-hazards mortality models on the intrinsic-capacity
# clinical features. Features are kept if measured in more than `feature_cutoff`
# participants, with the follow-up-only and help-seeking items in
# `domain_feature_exclude` removed. Discrimination is the held-out Harrell's C
# over age x sex-stratified 5-fold cross-validation; the out-of-fold linear
# predictors are the functional-age model for each domain.
#
# Run from the repository root:  Rscript code/01_mortality_models.R

library(tidyverse)
library(survival)
library(caret)

source("code/_config.R")

cutoff <- cfg$feature_cutoff
survival <- read_rds(out_path("processed", "survival.rds")) %>% na.omit %>% column_to_rownames("eid")

for (d in c("cognitive", "psychological", "locomotion", "vitality", "sensory")) {
  message(d)
  domain <- read_rds(out_path("processed", "domains", paste0(d, ".rds")))
  domain <- domain[, colSums(!is.na(domain)) > cutoff]
  domain <- domain[, !colnames(domain) %in% domain_feature_exclude, drop = FALSE]
  domain <- domain %>% na.omit

  s <- survival[intersect(rownames(domain), rownames(survival)), ]
  domain <- domain[rownames(s), ]
  features <- cbind(s %>% dplyr::select(time, status, age, sex), domain)

  strat <- interaction(cut(features$age, breaks = 5), features$sex, drop = TRUE)
  set.seed(cfg$seed)
  folds <- createFolds(strat, k = 5)
  c_index <- numeric(5)
  pred <- tibble()
  for (i in seq_along(folds)) {
    train <- features[-folds[[i]], ]
    test  <- features[folds[[i]], ]
    fit <- coxph(Surv(time, status) ~ ., data = train)
    risk <- predict(fit, newdata = test, type = "lp")
    pred <- rbind(pred, tibble(eid = rownames(test), fold = i, risk = risk, age = test$age, sex = test$sex))
    c_index[i] <- concordance(Surv(test$time, test$status) ~ risk, reverse = TRUE)$concordance
  }

  fit <- coxph(Surv(time, status) ~ ., data = features)
  results <- tibble(domain = d, cutoff = cutoff,
                    controls = sum(s$status == 0), cases = sum(s$status == 1),
                    nfeatures = ncol(features) - 2, nsamples = nrow(features),
                    c = list(c_index), coef = list(coef(fit) %>% enframe), pred = list(pred))
  write_rds(results, out_path("models", "clinical", paste0(d, ".rds")))
}
