# Charlson comorbidity index for UK Biobank participants (Charlson et al.,
# J Chronic Dis 1987), from hospital-episode ICD-10 diagnoses (field 41270) with
# their diagnosis dates (field 41280). Only diagnoses recorded on or before the
# assessment date (field 53) are counted, so the score reflects comorbidity at
# baseline and cannot borrow information from later follow-up. Codes are mapped
# with the Quan ICD-10 algorithm and weighted per Charlson; participants with no
# qualifying diagnosis score zero.
#
# Run from the repository root:  Rscript code/charlson.R

source("code/_config.R")
suppressMessages({library(data.table); library(comorbidity)})

code_cols <- sprintf("41270-0.%d", 0:258)
date_cols <- sprintf("41280-0.%d", 0:258)
d <- fread(in_path("phenotype_table"),
           select = c("eid", "53-0.0", code_cols, date_cols),
           colClasses = list(character = code_cols), showProgress = FALSE)

assessment <- d[, .(eid, assess = as.IDate(`53-0.0`))]

# Pair each diagnosis with its date, then keep those recorded by the baseline visit.
long <- melt(d, id.vars = "eid",
             measure.vars = list(code_cols, date_cols),
             value.name = c("code", "date"))
long <- long[!is.na(code) & code != ""]
long[, date := as.IDate(date)]
long <- assessment[long, on = "eid"]
long <- long[!is.na(date) & date <= assess, .(eid, code)]

cm <- comorbidity(x = long, id = "eid", code = "code",
                  map = "charlson_icd10_quan", assign0 = FALSE)
scored <- data.table(eid = cm$eid,
                     charlson = score(cm, weights = "charlson", assign0 = FALSE))

res <- merge(data.table(eid = d$eid), scored, by = "eid", all.x = TRUE)
res[is.na(charlson), charlson := 0]

cat(sprintf("Charlson index: %d participants, mean=%.3f, %% with score>0 = %.1f%%\n",
            nrow(res), mean(res$charlson), 100 * mean(res$charlson > 0)))
saveRDS(res, out_path("processed", "charlson.rds"))
