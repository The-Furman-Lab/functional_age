library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(purrr)
library(ggpubr)
library(ggpmisc)
library(cowplot)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(survival)
library(dplyr)
library(ggplot2)
library(ggrepel)

source("code/_config.R")

colors <- list(
  tone = c("#f6f2ed","#dedbca","#c4c0a5","#a59f83","#878264","#5e5948"),
  grey = c("#e5e4e8","#c5c9d6","#959fb3","#6e788c","#425468","#1b2942"),
  olive = c("#f2edb2","#dbdc63","#c4c400","#95a008","#637314","#304215"),
  green = c("#d7e5c5","#9fc978","#5db342","#41912f","#1d6e29","#0e3716"),
  teal = c("#c9e4ef","#96ced3","#48bcbc","#00959f","#006479","#003547"),
  blue = c("#c5e4fb","#9bc9e8","#5495ce","#006eae","#01478c","#002259"),
  purple = c("#e9d3e7","#d1a9ce","#b678b3","#a44990","#792373","#430b4d"),
  red = c("#f5cfc9","#e9a0a4","#db6463","#c5373c","#9c241b","#730c0d"),
  orange = c("#fbdcbc","#f9bd7b","#f29741","#e96700","#b34a00","#832a00"),
  yellow = c("#ffedc1","#f7dc86","#e8c54d","#ca9a23","#9b730a","#685409"),
  skin = c("#f6e5d3","#dcbc9f","#bc9778","#906852","#734d3d","#422a17"))
tone <- 3

data <- do.call("rbind", lapply(list.files(out_path("models", "proteomics"), full.names = TRUE), function(x) read_rds(x) %>% mutate(domain = gsub(".rds","",basename(x)))))
n_features <- data %>% dplyr::select(domain,coef) %>% unnest(coef) %>% group_by(domain) %>% summarise(n = n()-1)
data <- data %>% dplyr::select(domain,preds) %>% unnest(preds) 
data <- data %>% left_join(data %>% group_by(domain) %>% summarise(r = cor(functional_age, predicted_functional_age), mae = mean(abs(functional_age - predicted_functional_age))))
data <- data %>% left_join(n_features)

a <- ggplot(data, aes(functional_age, predicted_functional_age, fill = domain, color = domain))+
  geom_point(shape = 21, alpha = 0.2, size = 1.5, color = "white")+
  stat_correlation(label.x = 0.95, label.y = 0.17, color = "black", parse = FALSE, aes(label = paste("Features = ",paste0(n))))+
  stat_correlation(label.x = 0.95, label.y = 0.10, color = "black", parse = FALSE, aes(label = paste("R = ",paste0(round(r,2)))))+
  stat_correlation(label.x = 0.95, label.y = 0.03, color = "black", parse = FALSE, aes(label = paste("MAE = ",paste0(round(mae,2)))))+
  scale_fill_manual(values = c(colors$red[tone],colors$blue[tone],colors$teal[tone],colors$orange[tone],colors$yellow[tone]))+
  scale_color_manual(values = c(colors$red[5],colors$blue[5],colors$teal[5],colors$orange[5],colors$yellow[5]))+
  geom_smooth(method = "lm")+
  geom_abline(slope = 1, lty = 3)+
  facet_wrap(.~domain, nrow = 1)+
  theme_pubr(border = TRUE)+
  labs(x = "Functional age", y = "Predicted Functional age")+
  theme(legend.position = "none")

#figure 2b
data <- do.call("rbind", lapply(list.files(out_path("models", "proteomics"), full.names = TRUE), function(x) read_rds(x) %>% mutate(domain = gsub(".rds","",basename(x)))))
data <- data %>% dplyr::select(domain,coef) %>% unnest(coef) %>% filter(gene!="(Intercept)")
keep <- table(data$gene) %>% enframe %>% filter(value==1) %>% pull(name)

####
proteomics <- read_rds(in_path("proteomics"))
results <- read_rds(out_path("processed", "age.rds")) %>% ungroup %>% dplyr::select(eid,age) %>% unique %>% filter(eid%in%rownames(proteomics))
proteomics <- proteomics[results$eid,]
data <- data %>% left_join(apply(proteomics, 2, function(x) cor(x,results$age)) %>% enframe %>% set_names("gene","r"))

# 1. Prepare the data for plotting
# We identify unique proteins based on your 'keep' list
callouts <- c("NTRK2")   # shared proteins named in the text; labelled in grey, in their max-weight domain
data_plot <- data %>%
  mutate(status = ifelse(gene %in% keep, "Unique", "Shared")) %>%
  group_by(gene) %>%
  mutate(gene_max = abs(coefficient) == max(abs(coefficient))) %>%
  group_by(domain) %>%
  mutate(label = ifelse(
    (status == "Unique" &
       abs(coefficient) %in% sort(abs(coefficient[status == "Unique"]), decreasing = TRUE)[1:10]) |
      (gene %in% callouts & gene_max),
    gene, NA)) %>%
  ungroup()

b <- ggplot(data_plot, aes(x = r, y = coefficient, color = domain)) +
  geom_hline(yintercept = 0, lty = 3, color = "grey60") +
  geom_vline(xintercept = 0, lty = 3, color = "grey60") +
  geom_point(data = filter(data_plot, status == "Shared"), color = "#e5e4e8", size = 1) +
  geom_point(data = filter(data_plot, status == "Unique"), color = "black", size = 1) +
  scale_color_manual(values = c(colors$red[5],colors$blue[5],colors$teal[5],colors$orange[5],colors$yellow[5]))+
  geom_text_repel(data = filter(data_plot, status == "Unique"), aes(label = label),
                  size = 2.8, fontface = "bold",
                  max.overlaps = 100, box.padding = 0.5, point.padding = 0.3,
                  segment.color = "grey40", min.segment.length = 0, na.rm = TRUE) +
  geom_text_repel(data = filter(data_plot, status == "Shared"), aes(label = label),
                  color = "grey55", size = 2.8, fontface = "bold",
                  max.overlaps = 100, box.padding = 0.5, point.padding = 0.3,
                  segment.color = "grey70", min.segment.length = 0, na.rm = TRUE) +
  facet_wrap(~domain, scales = "free", nrow = 1) +
  theme_pubr(border = TRUE)+
  theme(legend.position = "none")+
  labs(x = "Correlation with age (R)", y = "Lasso Coefficient")

#enrichment
data <- do.call("rbind", lapply(list.files(out_path("models", "proteomics"), full.names = TRUE), function(x) read_rds(x) %>% mutate(domain = gsub(".rds","",basename(x)))))
data <- data %>% dplyr::select(domain,coef) %>% unnest(coef) %>% filter(gene!="(Intercept)")

ora_results_kegg <- data %>%
  filter(abs(coefficient) > 0) %>%
  group_by(domain) %>%
  group_split() %>%
  set_names(map_chr(., ~unique(.x$domain))) %>%
  map(~{
    gene_conv <- suppressWarnings(
      bitr(.x$gene, 
           fromType = "SYMBOL", 
           toType   = "ENTREZID", 
           OrgDb    = org.Hs.eg.db,
           drop     = TRUE)
    )
    enrichKEGG(
      gene          = gene_conv$ENTREZID,
      organism      = 'hsa',
      pAdjustMethod = "BH",
      minGSSize     = 10,
      maxGSSize     = 500,
      pvalueCutoff  = 1, 
      qvalueCutoff  = 1
    )
  })

combined_results <- map_dfr(ora_results_kegg, ~ .x@result, .id = "Domain") %>% filter(p.adjust < 0.05)
combined_results$score <- -log10(combined_results$pvalue)
combined_results$genelist <- sapply(combined_results$geneID, function(x) strsplit(x, split = "\\/"))

mat <- reshape2::acast(combined_results, Domain~Description, value.var = "score")
mat[is.na(mat)] <- 0
dim(mat)

col_fun <- colorRamp2(
  breaks = c(0, 2, 5, 8, 10), 
  colors = c(
    "white",    # Zero (Background)
    "#f7dc86",  # Yellow (Low signal)
    "#f29741",  # Orange (Medium signal)
    "#c5373c",  # Red (High signal)
    "#730c0d"   # Dark Red (Extreme signal)
  )
)

descriptions <- colSums(mat!=0, na.rm = TRUE) %>% enframe %>% filter(value==5) %>% pull(name)
positions <- sapply(descriptions, function(x) which(colnames(mat)==x))

gene_labels <- sapply(descriptions, function(desc) {
  genes <- table(combined_results %>% 
                   filter(Description == desc) %>% 
                   pull(genelist) %>% list %>% unlist) %>% 
    enframe %>% filter(value == 5) %>% pull(name)
  symbols <- bitr(genes, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL
  paste(symbols, collapse = "\n")
})

c(strsplit(gene_labels[["JAK-STAT signaling pathway"]], split = "\n")[[1]],
  strsplit(gene_labels[["Ras signaling pathway"]], split = "\n")[[1]],
  strsplit(gene_labels[["MAPK signaling pathway"]], split = "\n")[[1]],
  strsplit(gene_labels[["PI3K-Akt signaling pathway"]], split = "\n")[[1]]) %>% table %>% enframe %>% arrange(desc(value))

table(lapply(gene_labels, function(x) strsplit(x, split = "\n")) %>% unlist %>% as.vector) %>% enframe %>% arrange(desc(value))

my_14_cols <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", # 1-5
  "#A65628", "#F781BF", "#1B9E77", "#D95F02", "#7570B3", # 6-10
  "#E7298A", "#66A61E", "#E6AB02", "#000000", "#882255", # 11-15
  "#117733", "#004488", "#543005"                       # 16-18
)

ta = columnAnnotation(
  genes = anno_mark(
    at = positions, 
    labels = gene_labels,
    side = "top",
    # Keep your original styling, just add the color vector
    labels_gp = gpar(
      fontsize = 7, 
      fontface = "italic", 
      col = my_14_cols        # <--- Added unique colors to text
    ),
    labels_rot = 0,        
    link_height = unit(10, "mm"),
    padding = -2,
    extend = 0.1,
    link_width = unit(0, "mm"),
    # Add matching colors to the connector lines
    link_gp = gpar(
      col = my_14_cols,       # <--- Added unique colors to lines
      lwd = 1                 # Keeps lines clean and precise
    )
  )
)

ht <- Heatmap(mat,
              top_annotation = ta,  # <--- Add this line
              rect_gp = gpar(col = "black", lwd = 1, lty = 3),
              border = TRUE,
              height = unit(30, "mm"),
              column_title = "",
              cluster_columns = TRUE,
              row_dend_side = "right",
              column_names_gp = gpar(fontsize = 10),
              clustering_method_rows = "ward.D2",
              clustering_method_columns = "ward.D2",
              cluster_rows = TRUE,
              show_row_names = TRUE,
              row_names_side = "left",
              row_title_side = "left",
              show_column_names = TRUE,
              show_row_dend = FALSE,
              show_heatmap_legend = TRUE,
              show_column_dend = FALSE,
              heatmap_legend_param = list(
                direction = "vertical", 
                title = "-log(p)", 
                title_position = "topcenter", 
                legend_height = unit(4, "cm")
              ),
              col = col_fun
)

ordered <- descriptions[c(4,14,18,2,8,17,11,12,6,9,3,10,16,7,5,15,13,1)]

c <- grid.grabExpr({
  draw(ht, padding = unit(c(35, 1, 1, 1), "mm"))
  decorate_annotation("genes", {
    n_points <- length(positions)
    grid.points(
      x = c(1:11,13:19), 
      y = rep(0.01, 18),
      pch = 21,                          
      size = unit(3, "mm"), 
      gp = gpar(
        fill = my_14_cols[sapply(ordered, function(x) which(names(gene_labels)==x))],               
        col = "white",
        lwd = 0.5
      )
    )
  })
})

#figure d
survival <- read_rds(out_path("processed", "survival.rds"))
all_cindex <- tibble()

#clinical
results <- do.call("rbind", lapply(list.files(out_path("models", "clinical"), full.names = TRUE), function(x) read_rds(x)))
results$domain <- str_to_sentence(results$domain)
results <- results %>% dplyr::select(domain,pred) %>% unnest(pred)
results <- results %>% dplyr::select(domain,eid,risk)
colnames(results)[2] <- "eid"
colnames(results)[3] <- "prediction"
results$group <- "Clinical"
results <- results %>% left_join(survival %>% mutate(eid = as.character(eid)) %>% dplyr::select(eid,time,status))
all_cindex <- rbind(all_cindex,results)

#clinical + proteomics
results <- do.call("rbind", lapply(list.files(out_path("models", "clinical_proteomics"), full.names = TRUE), function(x) read_rds(x)))
results$domain <- str_to_sentence(results$domain)
results <- results %>% dplyr::select(domain,pred) %>% unnest(pred)
colnames(results)[2] <- "eid"
colnames(results)[3] <- "prediction"
results$group <- "Clinical + Proteomics"
results <- results %>% left_join(survival %>% mutate(eid = as.character(eid)) %>% dplyr::select(eid,time,status))
all_cindex <- rbind(all_cindex,results)

#clinical + biochemistry
results <- do.call("rbind", lapply(list.files(out_path("models", "clinical_biochemistry"), full.names = TRUE), function(x) read_rds(x)))
results$domain <- str_to_sentence(results$domain)
results <- results %>% dplyr::select(domain,pred) %>% unnest(pred)
colnames(results)[2] <- "eid"
colnames(results)[3] <- "prediction"
results$group <- "Clinical + Biochemistry"
results <- results %>% left_join(survival %>% mutate(eid = as.character(eid)) %>% dplyr::select(eid,time,status))
all_cindex <- rbind(all_cindex,results)

#clinical + proteomics + biochemistry
results <- do.call("rbind", lapply(list.files(out_path("models", "clinical_proteomics_biochemistry"), full.names = TRUE), function(x) read_rds(x)))
results$domain <- str_to_sentence(results$domain)
results <- results %>% dplyr::select(domain,pred) %>% unnest(pred)
colnames(results)[2] <- "eid"
colnames(results)[3] <- "prediction"
results$group <- "Clinical + Biochemistry + Proteomics"
results <- results %>% left_join(survival %>% mutate(eid = as.character(eid)) %>% dplyr::select(eid,time,status))
all_cindex <- rbind(all_cindex,results)

all_cindex <- all_cindex %>%
  group_by(domain) %>%
  mutate(n_groups_in_domain = n_distinct(group)) %>%
  group_by(domain, eid) %>%
  filter(n_distinct(group) == n_groups_in_domain) %>%
  dplyr::select(-n_groups_in_domain) %>%
  ungroup()

domain_results <- all_cindex %>% group_by(domain, group) %>% summarise(n = n(), c_index = 1-concordance(Surv(time, status) ~ prediction)$concordance) %>%ungroup()
domain_results %>% group_by(group) %>% summarise(mean = mean(c_index))

domain_results$group <- factor(domain_results$group, levels = c("Clinical","Clinical + Biochemistry","Clinical + Proteomics","Clinical + Biochemistry + Proteomics"))
domain_results$domain <- paste0(domain_results$domain, "\n(n = ",domain_results$n,")")

d <- ggplot(domain_results, aes(x = domain, y = c_index, fill = group)) +
  geom_col(
    position = position_dodge(width = 0.8), 
    width = 0.7, 
    color = "black", 
    linewidth = 0.2  
  ) +
  geom_text(
    aes(label = round(c_index, 3)), 
    position = position_dodge(width = 0.8), 
    vjust = -0.5, 
    size = 4, 
  ) +
  coord_cartesian(ylim = c(0.5, 0.85)) +
  labs(
    x = "Intrinsic Capacity domain",
    y = "C-index",
    fill = "Data type"
  ) +
  scale_fill_manual(values = c("#cccccc", "#9ecae1", "#4292c6", "#084594")) + 
  theme_pubr()+
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank()
  )

#plot all
pdf(file = out_path("figures", "Figure2.pdf"), width = 12, height = 18)
plot_grid(a,b,c,d, nrow = 4, rel_heights = c(0.7,0.65,2.2,1), labels = c("a","b","c","d"))
dev.off()


 