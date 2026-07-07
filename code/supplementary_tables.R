# Supplementary tables of cohort coverage: per-domain sample size and deaths
# (Table 1), and per-feature sample size, deaths and UK Biobank field id
# (Table 2), for the features that enter each functional-age model.
#
# Run from the repository root:  Rscript code/supplementary_tables.R

source("code/_config.R")
suppressMessages({
  library(tidyverse)
  library(reshape2)
})

survival <- read_rds(out_path("processed", "survival.rds"))
fields <- read_tsv(in_path("field_metadata"), show_col_types = FALSE) %>% mutate(field_id = as.character(field_id))

per_domain <- tibble()
per_feature <- tibble()
for (d in c("cognitive", "psychological", "locomotion", "vitality", "sensory")) {
  domain <- read_rds(out_path("processed", "domains", paste0(d, ".rds")))
  deaths <- reshape2::melt(domain) %>% na.omit %>% left_join(survival[, 1:2], by = c("Var1" = "eid")) %>%
    group_by(Var2) %>% summarise(cases = sum(status, na.rm = TRUE)) %>% set_names("feature", "cases")
  samples <- colSums(!is.na(domain)) %>% enframe %>% set_names("feature", "n") %>% left_join(deaths, by = "feature")

  domain <- domain[, colSums(!is.na(domain)) > cfg$feature_cutoff]
  domain <- domain[, !colnames(domain) %in% domain_feature_exclude, drop = FALSE]
  feats <- colnames(domain)
  # model cohort: complete domain features intersected with the survival table
  ids <- intersect(rownames(na.omit(domain)), survival$eid)
  n_samples <- length(ids)
  n_deaths <- survival %>% filter(eid %in% ids) %>% pull(status) %>% sum
  per_domain <- rbind(per_domain, tibble(domain = str_to_sentence(d), n = n_samples, cases = n_deaths))
  per_feature <- rbind(per_feature, tibble(domain = str_to_sentence(d), feature = feats) %>% left_join(samples, by = "feature"))
}

per_feature <- per_feature %>% left_join(fields[, 1:2] %>% set_names("field_id", "feature"), by = "feature") %>%
  group_by(domain, feature) %>% summarise(field_id = paste0(field_id, collapse = ", "), n = unique(n), cases = unique(cases), .groups = "drop")
per_feature$field_id[per_feature$feature == "Hand grip strength"] <- "46, 47"

per_domain <- per_domain %>% set_names("IC domain", "Sample size", "Deaths")
per_feature <- per_feature %>% set_names("IC domain", "Description", "UKB Field", "Sample size", "Deaths")

write_csv(per_domain, out_path("figures", "Supplementary_Table1_cohorts.csv"))
write_csv(per_feature, out_path("figures", "Supplementary_Table2_features.csv"))
