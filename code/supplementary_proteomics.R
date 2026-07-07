# Supplementary Figure 2 — how the proteomic surrogate weights distribute across
# the five functional-age domains: (a) number of proteins selected in 1..5 domains,
# and (b) their absolute correlation with chronological age, by number of domains.
#
# Run from the repository root:  Rscript code/supplementary_proteomics.R

source("code/_config.R")
suppressMessages({
  library(tidyverse)
  library(ggpubr)
  library(cowplot)
})

data <- do.call(rbind, lapply(list.files(out_path("models", "proteomics"), full.names = TRUE),
                              function(x) read_rds(x) %>% mutate(domain = gsub(".rds", "", basename(x))))) %>%
  dplyr::select(domain, coef) %>% unnest(coef) %>% filter(gene != "(Intercept)")

freq <- table(data$gene) %>% enframe %>% set_names("gene", "freq")
counts <- freq %>% count(freq, name = "value") %>% mutate(freq = as.numeric(freq))
p1 <- ggplot(counts, aes(factor(freq), value, fill = factor(freq))) +
  geom_col(color = "black", width = 0.7) +
  geom_text(aes(label = value), vjust = -0.5) +
  theme_pubr(border = TRUE) +
  scale_fill_manual(values = colors$green[2:6]) +
  labs(x = "Associated domains", y = "Number of proteins") +
  theme(legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

proteomics <- read_rds(in_path("proteomics"))
age <- read_rds(out_path("processed", "age.rds")) %>% ungroup %>%
  dplyr::select(eid, age) %>% unique %>% filter(eid %in% rownames(proteomics))
proteomics <- proteomics[age$eid, ]
genes <- freq %>% left_join(apply(proteomics, 2, function(x) cor(x, age$age)) %>% enframe %>%
                              set_names("gene", "r"), by = "gene")
p2 <- ggplot(genes %>% mutate(freq = as.character(freq)), aes(freq, abs(r), color = freq)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.2, width = 0.3) +
  stat_summary(fun = median, geom = "text", aes(label = after_stat(round(y, 3))), vjust = -0.25, color = "red") +
  theme_pubr(border = TRUE) +
  scale_color_manual(values = colors$green[2:6]) +
  labs(x = "Associated domains", y = "Absolute age correlation") +
  theme(legend.position = "none")

pdf(out_path("figures", "Sup_Figure2.pdf"), width = 8.30, height = 3.65)
print(plot_grid(p1, p2, nrow = 1, labels = c("a", "b")))
dev.off()
