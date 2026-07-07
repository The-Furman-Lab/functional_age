# Blood-biochemistry surrogate of age for one domain: LASSO (biochemistry markers
# -> chronological age), evaluated against functional age. Companion to the
# proteomic surrogate; used as the biochemistry feature set in the combined
# mortality models.
#
# Run from the repository root, with the domain id 1..5:
#   Rscript code/03_biochemistry_models.R <id>

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

biochemistry <- read_rds(in_path("biochemistry"))[, 3:38]
biochemistry <- biochemistry[, !colnames(biochemistry) %in% "30770"]
data <- results %>% filter(sample %in% rownames(biochemistry))
biochemistry <- biochemistry[data$sample, ]

lambdas <- 10^seq(-5, 5, by = 0.1)
set.seed(cfg$seed)
age_model <- cv.glmnet(x = as.matrix(biochemistry), y = data$age,
                       type.measure = "mae", family = "gaussian", alpha = 1,
                       keep = TRUE, lambda = lambdas, nfolds = 10)

preds <- tibble(sample = rownames(biochemistry), functional_age = data$functional_age,
                predicted_functional_age = age_model$fit.preval[, age_model$lambda == age_model$lambda.1se]) %>% na.omit()
coef <- coef(age_model, s = age_model$lambda.1se) %>%
  {tibble(gene = rownames(.), coefficient = as.matrix(.)[, 1])} %>% filter(coefficient != 0)

write_rds(tibble(domain = d, coef = list(coef), preds = list(preds)),
          out_path("models", "biochemistry", paste0(d, ".rds")))
