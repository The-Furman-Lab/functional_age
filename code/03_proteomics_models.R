# Proteomic surrogate of functional age: LASSO (Olink proteins -> functional age)
# for one domain. Functional age is the age-rescaled out-of-fold linear predictor
# of the domain's clinical mortality model.
#
# Run from the repository root, with the domain id 1..5:
#   Rscript code/03_proteomics_models.R <id>

library(glmnet)
library(tidyverse)
source("code/_config.R")

d <- c("Cognitive", "Locomotion", "Psychological", "Sensory", "Vitality")[[as.numeric(commandArgs(trailingOnly = TRUE)[1])]]

age_transformation <- function(age, xb) {
  ((xb - mean(xb)) / sd(xb)) * sd(age) + mean(age)
}

results <- do.call(rbind, lapply(list.files(out_path("models", "clinical"), full.names = TRUE), read_rds))
results$domain <- str_to_sentence(results$domain)
results <- results %>% filter(domain == d) %>% dplyr::select(domain, pred) %>% unnest(pred) %>%
  group_by(domain) %>% mutate(functional_age = age_transformation(age, risk)) %>% na.omit()
colnames(results)[2] <- "sample"

proteomics <- read_rds(in_path("proteomics"))
data <- results %>% filter(sample %in% rownames(proteomics))
proteomics <- proteomics[data$sample, ]

lambdas <- 10^seq(-5, 5, by = 0.1)
set.seed(cfg$seed)
age_model <- cv.glmnet(x = as.matrix(proteomics), y = data$functional_age,
                       type.measure = "mae", family = "gaussian", alpha = 1,
                       keep = TRUE, lambda = lambdas, nfolds = 10)

preds <- tibble(sample = rownames(proteomics), functional_age = data$functional_age,
                predicted_functional_age = age_model$fit.preval[, age_model$lambda == age_model$lambda.1se]) %>% na.omit()
coef <- coef(age_model, s = age_model$lambda.1se) %>%
  {tibble(gene = rownames(.), coefficient = as.matrix(.)[, 1])} %>% filter(coefficient != 0)

write_rds(tibble(domain = d, coef = list(coef), preds = list(preds)),
          out_path("models", "proteomics", paste0(d, ".rds")))
