# Williams frailty index for UK Biobank participants (Williams et al.,
# J Gerontol A 2019): 49 self-report deficits at the baseline visit. Each item
# scores 0 (absent) to 1 (present), graded for a few ordinal items; the index is
# the mean over non-missing items, kept for participants with at least 40/49.
#
# Run from the repository root:  Rscript code/frailty_index.R

source("code/_config.R")
suppressMessages({library(data.table); library(tidyverse)})

# Baseline (instance 0) fields: single-answer items and the multi-select fields.
singles <- c("2247-0.0","2316-0.0","2178-0.0","2080-0.0","1200-0.0","2050-0.0",
             "1970-0.0","2020-0.0","1930-0.0","2188-0.0","2296-0.0","2463-0.0",
             "2443-0.0","2453-0.0","2335-0.0","134-0.0")
ms <- list("6148"=0:4,"6150"=0:3,"6152"=0:4,"6149"=0:5,"6159"=0:6,"20002"=0:33)
ms_cols <- unlist(lapply(names(ms), function(f) sprintf("%s-0.%d", f, ms[[f]])))

d <- fread(in_path("phenotype_table"), select = c("eid", singles, ms_cols),
           showProgress = FALSE)

# Negative sentinels (don't-know / prefer-not-to-answer) become missing.
neg_na <- function(x) { x[x %in% c(-1, -3)] <- NA; x }
recode_map <- function(col, map) as.numeric(map[as.character(neg_na(d[[col]]))])
recode_bin <- function(col) as.numeric(neg_na(d[[col]]) == 1)

ms_mat <- function(f) as.matrix(d[, sprintf("%s-0.%d", f, ms[[f]]), with = FALSE])
ms_assessed <- function(M) { M[M %in% c(-1, -3)] <- NA; rowSums(!is.na(M)) > 0 }
ms_item <- function(M, codes, assessed) {
  hit <- Reduce(`|`, lapply(codes, function(t) rowSums(M == t, na.rm = TRUE) > 0))
  ifelse(hit, 1, ifelse(assessed, 0, NA_real_))
}

M6148 <- ms_mat("6148"); a6148 <- ms_assessed(M6148)
M6150 <- ms_mat("6150"); a6150 <- ms_assessed(M6150)
M6152 <- ms_mat("6152"); a6152 <- ms_assessed(M6152)
M6149 <- ms_mat("6149"); a6149 <- ms_assessed(M6149)
M6159 <- ms_mat("6159"); a6159 <- ms_assessed(M6159)
M20002 <- ms_mat("20002")                       # self-reported conditions: absence = 0
sr <- function(code) as.numeric(rowSums(M20002 == code, na.rm = TRUE) > 0)

FI <- data.table(eid = d$eid)
FI[, i01_glaucoma  := ms_item(M6148, 2, a6148)]
FI[, i02_cataract  := ms_item(M6148, 4, a6148)]
FI[, i03_hearing   := { x <- neg_na(d[["2247-0.0"]]); ifelse(is.na(x), NA, as.numeric(x == 1 | x == 99)) }]
FI[, i04_migraine  := sr(1265)]
FI[, i05_dental    := ms_item(M6149, 1:6, a6149)]
FI[, i06_srh       := recode_map("2178-0.0", c("1"=0,"2"=0.25,"3"=0.5,"4"=1))]
FI[, i07_fatigue   := recode_map("2080-0.0", c("1"=0,"2"=0.25,"3"=0.5,"4"=1))]
FI[, i08_sleep     := recode_map("1200-0.0", c("1"=0,"2"=0.5,"3"=1))]
FI[, i09_depmood   := recode_map("2050-0.0", c("1"=0,"2"=0.5,"3"=0.75,"4"=1))]
FI[, i10_nervous   := recode_bin("1970-0.0")]
FI[, i11_anxiety   := sr(1287)]
FI[, i12_lonely    := recode_bin("2020-0.0")]
FI[, i13_misery    := recode_bin("1930-0.0")]
FI[, i14_longstand := recode_bin("2188-0.0")]
FI[, i15_falls     := recode_map("2296-0.0", c("1"=0,"2"=0.5,"3"=1))]
FI[, i16_fracture  := recode_bin("2463-0.0")]
FI[, i17_diabetes  := recode_bin("2443-0.0")]
FI[, i18_mi        := ms_item(M6150, 1, a6150)]
FI[, i19_angina    := ms_item(M6150, 2, a6150)]
FI[, i20_stroke    := ms_item(M6150, 3, a6150)]
FI[, i21_highbp    := ms_item(M6150, 4, a6150)]
FI[, i22_hypothyr  := sr(1226)]
FI[, i23_dvt       := ms_item(M6152, 5, a6152)]
FI[, i24_highchol  := sr(1473)]
FI[, i25_wheeze    := recode_bin("2316-0.0")]
FI[, i26_pneumonia := sr(1398)]
FI[, i27_emphysema := ms_item(M6152, 7, a6152)]
FI[, i28_asthma    := ms_item(M6152, 8, a6152)]
FI[, i29_ra        := sr(1464)]
FI[, i30_oa        := sr(1465)]
FI[, i31_gout      := sr(1466)]
FI[, i32_osteopor  := sr(1309)]
FI[, i33_hayfever  := ms_item(M6152, 9, a6152)]
FI[, i34_psoriasis := sr(1453)]
FI[, i35_cancer    := recode_bin("2453-0.0")]
FI[, i36_multicanc := { x <- d[["134-0.0"]]; ifelse(is.na(x), NA, as.numeric(x >= 2)) }]
FI[, i37_chestpain := recode_bin("2335-0.0")]
FI[, i38_headneck  := ms_item(M6159, c(1, 3), a6159)]
FI[, i39_backpain  := ms_item(M6159, 4, a6159)]
FI[, i40_stomach   := ms_item(M6159, 5, a6159)]
FI[, i41_hippain   := ms_item(M6159, 6, a6159)]
FI[, i42_kneepain  := ms_item(M6159, 7, a6159)]
FI[, i43_allover   := ms_item(M6159, 8, a6159)]
FI[, i44_facial    := ms_item(M6159, 2, a6159)]
FI[, i45_sciatica  := sr(1476)]
FI[, i46_reflux    := sr(1138)]
FI[, i47_hiatus    := sr(1474)]
FI[, i48_gallstone := sr(1162)]
FI[, i49_divertic  := sr(1458)]

items <- grep("^i[0-9]", names(FI), value = TRUE)
stopifnot(length(items) == 49)
IM <- as.matrix(FI[, items, with = FALSE])
n_items <- rowSums(!is.na(IM))
FI[, n_items := n_items]
FI[, frailty_index := ifelse(n_items >= 40, rowSums(IM, na.rm = TRUE) / n_items, NA_real_)]

cat(sprintf("frailty index: %d participants (>=40/49 items), mean=%.4f sd=%.4f\n",
            sum(!is.na(FI$frailty_index)), mean(FI$frailty_index, na.rm = TRUE),
            sd(FI$frailty_index, na.rm = TRUE)))
saveRDS(FI[, .(eid, n_items, frailty_index)], out_path("processed", "frailty_index.rds"))
