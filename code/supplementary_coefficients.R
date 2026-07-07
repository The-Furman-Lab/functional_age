# Supplementary coefficient table: every non-zero model coefficient across all
# data types (clinical, proteomics, biochemistry, and their combinations) for each
# functional-age domain. Biochemistry field codes are mapped to readable names via
# the UK Biobank crosswalk.
#
# Run from the repository root:  Rscript code/supplementary_coefficients.R

source("code/_config.R")
suppressMessages({
  library(tidyverse)
  library(readxl)
})

crosswalk <- read_tsv(in_path("field_metadata"), show_col_types = FALSE) %>%
  dplyr::select(field_id, title) %>% set_names("code", "name") %>% mutate(code = as.character(code))

load_coef <- function(kind) {
  do.call(rbind, lapply(list.files(out_path("models", kind), full.names = TRUE), read_rds)) %>%
    dplyr::select(domain, coef) %>% unnest(coef) %>% set_names("Domain", "Feature", "Coefficient")
}
map_codes <- function(df) {
  df %>% left_join(crosswalk, by = c("Feature" = "code")) %>%
    mutate(Feature = coalesce(name, Feature)) %>% dplyr::select(Domain, Feature, Coefficient)
}

coefficients <- bind_rows(
  load_coef("clinical") %>% mutate(Feature = str_to_sentence(str_replace_all(Feature, "`", ""))) %>% mutate(Data = "Clinical"),
  map_codes(load_coef("biochemistry")) %>% mutate(Data = "Biochemistry"),
  load_coef("proteomics") %>% mutate(Data = "Proteomics"),
  load_coef("clinical_proteomics") %>% mutate(Data = "Clinical + Proteomics"),
  map_codes(load_coef("clinical_biochemistry")) %>% mutate(Data = "Clinical + Biochemistry"),
  map_codes(load_coef("clinical_proteomics_biochemistry")) %>% mutate(Data = "Clinical + Proteomics + Biochemistry")
) %>%
  filter(Feature != "(Intercept)") %>%
  mutate(Domain = str_to_sentence(Domain),
         Feature = recode(Feature, age = "Age", sex = "Sex", Sex = "Sex")) %>%
  dplyr::select(Data, Domain, Feature, Coefficient)

write_csv(coefficients, out_path("figures", "Supplementary_Coefficients.csv"))
cat(sprintf("wrote %d coefficients across %d data types x %d domains\n",
            nrow(coefficients), n_distinct(coefficients$Data), n_distinct(coefficients$Domain)))
