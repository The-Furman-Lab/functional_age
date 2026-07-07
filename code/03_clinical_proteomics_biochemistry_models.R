# Combined mortality model for one domain: Cox-LASSO on the clinical domain
# features (forced in, penalty 0) plus the proteins and biochemistry markers
# selected by the domain's proteomic and biochemistry surrogates (penalised).
# Out-of-fold linear predictors give the discrimination of clinical + both omics.
#
# Run from the repository root, with the domain id 1..5:
#   Rscript code/03_clinical_proteomics_biochemistry_models.R <id>

library(glmnet)
library(survival)
library(tidyverse)
source("code/_config.R")

d <- c("cognitive", "psychological", "locomotion", "vitality", "sensory")[[as.numeric(commandArgs(trailingOnly = TRUE)[1])]]

load_selected <- function(kind) {
  do.call(rbind, lapply(list.files(out_path("models", kind), full.names = TRUE),
                        function(x) read_rds(x) %>% mutate(domain = gsub(".rds", "", basename(x))))) %>%
    dplyr::select(domain, coef) %>% unnest(coef) %>% filter(gene != "(Intercept)")
}
prot <- load_selected("proteomics")
bioc <- load_selected("biochemistry")

survival <- read_rds(out_path("processed", "survival.rds")) %>% na.omit() %>% column_to_rownames("eid")
proteomics <- read_rds(in_path("proteomics"))
biochemistry <- read_rds(in_path("biochemistry"))

domain <- read_rds(out_path("processed", "domains", paste0(d, ".rds")))
domain <- domain[, colSums(!is.na(domain)) > cfg$feature_cutoff]
domain <- domain[, !colnames(domain) %in% domain_feature_exclude, drop = FALSE]
domain <- domain %>% na.omit()

survival <- survival[intersect(rownames(domain), rownames(survival)), ]
survival <- survival[intersect(rownames(proteomics), rownames(survival)), ]
survival <- survival[intersect(rownames(biochemistry), rownames(survival)), ]
domain <- domain[rownames(survival), ]

proteins <- prot$gene[prot$domain == str_to_sentence(d)]
markers  <- bioc$gene[bioc$domain == str_to_sentence(d)]
features <- cbind(survival %>% dplyr::select(time, status, age, sex), domain,
                  proteomics[rownames(survival), proteins], biochemistry[rownames(survival), markers])

x <- as.matrix(features[, 3:ncol(features)])
p_factors <- rep(1, ncol(x))
p_factors[1:(ncol(domain) + 2)] <- 0
set.seed(cfg$seed)
cv_fit <- cv.glmnet(x, Surv(features$time, features$status), family = "cox", alpha = 1,
                    keep = TRUE, penalty.factor = p_factors)

pred <- cv_fit$fit.preval[, cv_fit$lambda == cv_fit$lambda.1se] %>% enframe %>% set_names("sample", "prediction")
coef <- coef(cv_fit, s = "lambda.1se") %>% as.matrix() %>% data.frame %>% rownames_to_column("feature") %>%
  set_names("feature", "coefficient") %>% filter(coefficient != 0)

write_rds(tibble(domain = d, coef = list(coef), pred = list(pred)),
          out_path("models", "clinical_proteomics_biochemistry", paste0(d, ".rds")))
