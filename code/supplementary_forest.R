# Supplementary Figure 3 — forest of marginal (age- and sex-adjusted) hazard ratios
# per 1 SD for each domain feature, from marginal_hr.R. Age and sex are omitted from
# the panels (age HR/SD ~2.35 would dominate the scale). Direction is carried by
# point colour, so only reverse-coded / opaque-unit features are relabelled.
#
# Run from the repository root:  Rscript code/supplementary_forest.R
# (requires code/marginal_hr.R to have been run first)

source("code/_config.R")
suppressMessages(library(tidyverse))

relabel <- c(
  "Health satisfaction"                               = "Poorer health satisfaction",
  "Financial situation satisfaction"                  = "Poorer financial satisfaction",
  "Happiness"                                         = "Poorer happiness",
  "Speech-reception-threshold (SRT) estimate (right)" = "Worse hearing (SRT, right)",
  "Speech-reception-threshold (SRT) estimate (left)"  = "Worse hearing (SRT, left)")

feat <- read_csv(out_path("processed", "marginal_hr.csv"), show_col_types = FALSE) %>%
  mutate(term_disp = ifelse(feature %in% names(relabel), relabel[feature], feature),
         lab = str_wrap(str_trunc(term_disp, 46), 34),
         dir = ifelse(HR_perSD > 1, "higher value → higher risk", "higher value → lower risk"),
         abseff = abs(std_beta),
         domain = factor(domain, levels = c("Vitality", "Locomotion", "Psychological", "Cognitive", "Sensory"))) %>%
  arrange(domain, abseff)
feat$lab <- factor(feat$lab, levels = unique(feat$lab))

p <- ggplot(feat, aes(x = HR_perSD, y = lab, color = dir)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey55") +
  geom_errorbarh(aes(xmin = HR_lo, xmax = HR_hi), height = 0.28, linewidth = 0.5) +
  geom_point(size = 2) +
  facet_grid(rows = vars(domain), scales = "free_y", space = "free_y") +
  scale_x_log10(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4, 1.6)) +
  scale_color_manual(values = c("higher value → higher risk" = "#B2182B",
                                "higher value → lower risk" = "#2166AC"), name = NULL) +
  labs(x = "Marginal hazard ratio per 1 SD (age + sex adjusted, log scale)", y = NULL,
       title = "Marginal association of intrinsic capacity features with all-cause mortality") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5),
        strip.text.y = element_text(angle = 0, face = "bold"),
        panel.grid.minor = element_blank())

ggsave(out_path("figures", "Sup_Figure3.pdf"), p, width = 9.5, height = max(12, nrow(feat) * 0.23 + 2))
