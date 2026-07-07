# Functional-age acceleration versus diseases first diagnosed during follow-up.
#
# For the disease-association panel we count only strictly incident cases: a
# first hospital diagnosis after the baseline assessment (field 53) that was not
# already present as self-report or an earlier hospital diagnosis. Conditions
# prevalent at baseline are excluded, so the association reflects functional
# decline that precedes diagnosis rather than existing disease. Per domain and
# disease we report the mean acceleration (years) over incident cases, a
# one-sample t-test against zero, and the BH false-discovery rate, plus the lead
# time from assessment to first diagnosis. The disease taxonomy is the Tian et
# al. map (tian_fields.xlsx sheet 7); cancer is added from ICD-10 C00-C97
# (excluding C44), ICD-9 140-208, and self-report (field 20001).
#
# Run from the repository root:  Rscript code/incident_disease.R

source("code/_config.R")
suppressMessages({library(data.table); library(readxl); library(stringr)})

domains <- c("cognitive", "psychological", "locomotion", "vitality", "sensory")
age_transformation <- function(age, xb) ((xb - mean(xb)) / sd(xb)) * sd(age) + mean(age)

# ---- functional-age acceleration on the 5-domain-complete cohort ------------
pred <- rbindlist(lapply(domains, function(d) {
  p <- as.data.table(readRDS(out_path("models", "clinical", paste0(d, ".rds")))$pred[[1]])
  p[, `:=`(domain = str_to_sentence(d), eid = as.character(eid))]
  p
}))
common <- pred[, .N, by = eid][N == length(domains), eid]
acc <- pred[eid %in% common]
acc[, functional_age := age_transformation(age, risk), by = domain]
acc[, functional_acc := resid(lm(functional_age ~ age + sex)), by = domain]
acc <- acc[, .(domain, eid, functional_acc)]

# ---- disease taxonomy and one read of the wide table ------------------------
map <- as.data.table(read_xlsx(in_path("tian_fields"), sheet = 7))
map[, field_id := str_trim(as.character(field_id))]
sr_map    <- map[grepl("Self Report", type), .(field_id, disease)]
icd9_map  <- map[type == "ICD 9",  .(field_id, disease)]
icd10_map <- map[type == "ICD 10", .(field_id, disease)]

c10 <- sprintf("41270-0.%d", 0:258); d10 <- sprintf("41280-0.%d", 0:258)
c9  <- sprintf("41271-0.%d", 0:46);  d9  <- sprintf("41281-0.%d", 0:46)
sr  <- sprintf("20002-0.%d", 0:33);  src <- sprintf("20001-0.%d", 0:5)
w <- fread(in_path("phenotype_table"),
           select = c("eid", "53-0.0", sr, src, c10, d10, c9, d9),
           colClasses = list(character = c("53-0.0", sr, src, c10, d10, c9, d9)),
           showProgress = FALSE)
w[, eid := as.character(eid)]
setnames(w, "53-0.0", "assess"); w[, assess := as.IDate(assess)]

# prevalent from self-report: mapped conditions (20002) and any cancer (20001)
sr_long <- melt(w[, c("eid", sr), with = FALSE], id.vars = "eid", value.name = "code")[code != "" & !is.na(code)]
sr_long[, code := str_trim(code)]
prev_sr <- unique(merge(sr_long[, .(eid, code)], sr_map, by.x = "code", by.y = "field_id")[, .(eid, disease)])
cancer_sr <- unique(melt(w[, c("eid", src), with = FALSE], id.vars = "eid", value.name = "code")[code != "" & !is.na(code), .(eid)])
prev_sr <- unique(rbind(prev_sr, cancer_sr[, .(eid, disease = "Cancer")]))

# dated hospital diagnoses, code paired to date by array position
hes_long <- function(code_cols, date_cols, code_map, system) {
  cd <- melt(w[, c("eid", code_cols), with = FALSE], id.vars = "eid", variable.name = "col", value.name = "code")
  cd[, k := sub(".*-0\\.", "", col)]; cd <- cd[code != "" & !is.na(code), .(eid, k, code = str_trim(code))]
  dt <- melt(w[, c("eid", date_cols), with = FALSE], id.vars = "eid", variable.name = "col", value.name = "dt")
  dt[, k := sub(".*-0\\.", "", col)]; dt <- dt[, .(eid, k, dt = as.IDate(dt))]
  x <- merge(cd, dt, by = c("eid", "k"))
  x <- merge(x[, .(eid, code, dt)], w[, .(eid, assess)], by = "eid")
  tian <- merge(x, code_map, by.x = "code", by.y = "field_id")[, .(eid, disease, dt, assess)]
  if (system == "icd10") {
    canc <- x[grepl("^C[0-9]{2}", code) & suppressWarnings(as.integer(substr(code, 2, 3))) <= 97 & substr(code, 1, 3) != "C44"]
  } else {
    canc <- x[{ v <- suppressWarnings(as.integer(substr(code, 1, 3))); !is.na(v) & v %in% 140:208 }]
  }
  rbind(tian, canc[, .(eid, disease = "Cancer", dt, assess)])
}
hes <- rbind(hes_long(c10, d10, icd10_map, "icd10"), hes_long(c9, d9, icd9_map, "icd9"))
hes <- hes[!is.na(dt) & dt > as.IDate("1910-01-01")]

# prevalent vs strict incident, with first-incident date
agg <- hes[, .(prev_hes = any(dt <= assess), inc_hes = any(dt > assess),
               first_inc = if (any(dt > assess)) min(dt[dt > assess]) else as.IDate(NA),
               assess = assess[1]), by = .(eid, disease)]
prev <- unique(rbind(prev_sr, agg[prev_hes == TRUE, .(eid, disease)]))[, is_prev := TRUE]
inc <- agg[inc_hes == TRUE, .(eid, disease, first_inc, assess)]
inc <- merge(inc, prev, by = c("eid", "disease"), all.x = TRUE)[is.na(is_prev)]
inc[, lead_y := as.numeric(first_inc - assess) / 365.25]
inc <- inc[lead_y > 0]

# ---- association in years, per domain and disease ---------------------------
inc_cohort <- inc[eid %in% common]
assoc <- merge(acc, inc_cohort[, .(eid, disease)], by = "eid", allow.cartesian = TRUE)
assoc <- assoc[, {
  mu <- mean(functional_acc)
  p <- if (.N >= 3 && sd(functional_acc) > 0) t.test(functional_acc, mu = 0)$p.value else NA_real_
  .(n = .N, mean_years = mu, p = p)
}, by = .(domain, disease)]
assoc[, fdr := p.adjust(p, "BH")][, sig := !is.na(fdr) & fdr < 0.05]
setorder(assoc, -mean_years)

cat(sprintf("incident pairs (cohort): %d ; diseases: %d ; median lead time: %.1f years\n",
            nrow(inc_cohort), uniqueN(inc_cohort$disease), median(inc_cohort$lead_y)))
fwrite(assoc[, .(domain, disease, n, mean_years = round(mean_years, 3), p, fdr, sig)],
       out_path("processed", "incident_disease_years.csv"))
