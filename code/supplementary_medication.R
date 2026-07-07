# Supplementary Figure 1 — functional-age acceleration difference between regular
# users and non-users of each self-reported medication and dietary supplement
# (touchscreen fields; see medication in 00_tidy_domains.R). Per medication and
# domain we report the mean acceleration in users minus non-users, a Wilcoxon
# test and BH false-discovery rate; cells passing FDR < 0.05 are outlined.
#
# Run from the repository root:  Rscript code/supplementary_medication.R

source("code/_config.R")
suppressMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
})

medication <- read_rds(out_path("processed", "medication.rds")) %>% mutate(eid = as.character(eid))
acc <- read_rds(out_path("processed", "age.rds")) %>% ungroup() %>%
  transmute(domain, eid = as.character(eid), functional_acc)
medication <- medication %>% inner_join(acc, by = "eid")

# keep medications with enough users and non-users to compare
keep <- medication %>% count(meaning, domain, group) %>%
  pivot_wider(names_from = group, values_from = n, values_fill = 0) %>%
  filter(cases > 100, controls > 100)
medication <- medication %>% semi_join(keep, by = c("meaning", "domain"))

res <- medication %>% group_by(meaning, domain) %>%
  summarise(mean_diff = mean(functional_acc[group == "cases"]) - mean(functional_acc[group == "controls"]),
            p = wilcox.test(functional_acc[group == "cases"], functional_acc[group == "controls"])$p.value,
            .groups = "drop")
res$fdr <- p.adjust(res$p, "BH")

diff_mat <- reshape2::acast(res, meaning ~ domain, value.var = "mean_diff")
fdr_mat  <- reshape2::acast(res, meaning ~ domain, value.var = "fdr")

col_fun <- colorRamp2(seq(-1.5, 1.5, 0.3), rev(brewer.pal(11, "RdBu")))
ht <- Heatmap(diff_mat, border = TRUE, col = col_fun,
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.2f", diff_mat[i, j]), x, y, gp = gpar(fontsize = 8))
    if (fdr_mat[i, j] < 0.05)
      grid.rect(x, y, w * 0.8, h * 0.8, gp = gpar(col = "black", fill = NA, lwd = 1))
  },
  cluster_columns = TRUE, cluster_rows = TRUE,
  clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
  row_names_side = "left", column_names_gp = gpar(fontsize = 10),
  show_row_dend = FALSE, show_column_dend = FALSE, show_heatmap_legend = FALSE)

pdf(out_path("figures", "Sup_Figure1.pdf"), width = 5.5, height = 8.35)
draw(ht, padding = unit(c(1, 1, 1, 25), "mm"))
draw(Legend(col_fun = col_fun, title = "Age\nacceleration\ndifference",
            direction = "vertical", legend_height = unit(4, "cm"), title_position = "topcenter"),
     x = unit(0.95, "npc"), y = unit(0.5, "npc"), just = c("right", "center"))
dev.off()
