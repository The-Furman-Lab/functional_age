library(ggpmisc)
library(ggpubr)
library(tidyverse)
library(data.table)
library(survival)
library(caret)  
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)

source("code/_config.R")
tone <- 3

age_transformation <- function(age,xb){
  m.age=mean(age)
  sd.age=sd(age)
  m.cox=mean(xb)
  sd.cox=sd(xb)
  Y0 <- xb
  Y=(Y0-m.cox)/sd.cox
  bioage <- as.numeric((Y*sd.age)+m.age)
  return(bioage)  
}

calc_aa <- function(x,age,sex){
  return(resid(lm(x~age+sex)))
}

results <- do.call("rbind", lapply(list.files(out_path("models", "clinical"), full.names = TRUE), function(x) read_rds(x)))
results$domain <- str_to_sentence(results$domain)
results <- results %>% dplyr::select(domain,pred) %>% unnest(pred)
common <- table(results$eid) %>% enframe %>% filter(value==5) %>% pull(name)
results <- results %>% filter(eid%in%common)
results <- results %>% group_by(domain) %>% mutate(functional_age = age_transformation(age,risk)) %>% na.omit
results <- results %>% group_by(domain) %>% mutate(functional_acc = resid(lm(functional_age~age))) %>% na.omit
results_fixed <- results

#figure a
results <- results_fixed %>% group_by(domain) %>% mutate(r = cor(age,risk)) 
add_variables <- do.call("rbind", lapply(list.files(out_path("models", "clinical"), full.names = TRUE), function(x) read_rds(x)))
add_variables <- add_variables %>% dplyr::select(domain,c,nfeatures) %>% unnest(c) %>%
  mutate(nfeatures = nfeatures - 2) %>%   # count intrinsic-capacity features, excluding the age and sex covariates
  group_by(domain,nfeatures) %>% summarise(c = mean(c))
add_variables$domain <- str_to_sentence(add_variables$domain)
results <- results %>% left_join(add_variables)

a <- ggplot(results, aes(age, functional_age, fill = domain, color = domain))+
  geom_point(shape = 21, alpha = 0.2, size = 1.5, color = "white")+
  stat_correlation(label.x = 0.95, label.y = 0.17, color = "black", parse = FALSE, aes(label = paste("Features = ",paste0(nfeatures))))+
  stat_correlation(label.x = 0.95, label.y = 0.10, color = "black", parse = FALSE, aes(label = paste("C-index = ",paste0(round(c,2)))))+
  stat_correlation(label.x = 0.95, label.y = 0.03, color = "black", parse = FALSE, aes(label = paste("R = ",paste0(round(r,2)))))+
  scale_fill_manual(values = c(colors$red[tone],colors$blue[tone],colors$teal[tone],colors$orange[tone],colors$yellow[tone]))+
  scale_color_manual(values = c(colors$red[5],colors$blue[5],colors$teal[5],colors$orange[5],colors$yellow[5]))+
  geom_smooth(method = "lm")+
  geom_abline(slope = 1, lty = 3)+
  facet_wrap(.~domain, nrow = 1)+
  theme_pubr(border = TRUE)+
  labs(x = "Age", y = "Functional age")+
  theme(legend.position = "none")

#figure b
results_sex <- results_fixed %>% group_by(domain,sex) %>%
  summarise(m = mean(functional_acc), se = sd(functional_acc)/sqrt(n()), df = n()-1, .groups = "drop") %>%
  transmute(domain, sex = ifelse(sex==0,"Female","Male"), mean = m,
            lower = m - qt(0.995, df)*se, upper = m + qt(0.995, df)*se)
results_sex$domain <- factor(results_sex$domain, levels = results_sex %>% group_by(domain) %>% summarise(sum = sum(abs(mean))) %>% arrange(sum) %>% pull(domain))

b <- ggplot(results_sex, aes(x = domain, y = mean, fill = sex, group = sex)) +
  geom_bar(stat = "identity", position = "identity", width = 0.6, alpha = 0.7) +  
  geom_errorbar(aes(ymin = lower, ymax = upper), position = "identity", width = 0.2) + 
  theme_pubr(border = TRUE) +
  scale_fill_manual(values = c(colors$orange[[2]], colors$green[[2]])) +
  geom_hline(yintercept = 0, lty = 3) +
  coord_flip() +
  ylim(c(-3, 3)) +
  labs(x = "Domain", y = "Mean functional-age acceleration", fill = "Sex") +
  geom_text(aes(label = round(mean, 2), 
                y = mean, 
                hjust = ifelse(mean > 0, -0.2, 1.2)), 
            position = "identity")

#figure c
results <- do.call("rbind", lapply(list.files(out_path("models", "clinical"), full.names = TRUE), function(x) read_rds(x)))
results$domain <- str_to_sentence(results$domain)
results <- results %>% dplyr::select(domain,pred) %>% unnest(pred)
common <- table(results$eid) %>% enframe %>% filter(value==5) %>% pull(name)
results <- results %>% filter(eid%in%common)
results <- results %>% group_by(domain) %>% mutate(functional_age = age_transformation(age,risk)) %>% na.omit
results <- results %>% group_by(domain) %>% mutate(functional_acc = resid(lm(functional_age~age+sex))) %>% na.omit
write_rds(results, file = out_path("processed", "age.rds"))
mat <- reshape2::acast(results, eid~domain, value.var = "functional_acc")
cor_mat <- cor(mat)

col_fun <- colorRamp2(seq(-1,1,0.2), c(colors$blue[5:1],"white",colors$red[1:5]))
ht <- Heatmap(cor_mat,
              border = TRUE,
              rect_gp = gpar(col = "black", lwd = 1, lty = 3),
              cell_fun = function(j, i, x, y, w, h, col) {
                grid.text(round(cor_mat,2)[i, j], x, y, gp = gpar(fontsize = 12))
              },
              cluster_columns = TRUE,
              row_dend_side = "right",
              row_names_side = "left",
              cluster_rows = TRUE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              show_row_dend = TRUE,
              show_column_dend = TRUE,
              #show_heatmap_legend = F,
              heatmap_legend_param = list(direction = "vertical", title = "R", title_position = "topcenter", legend_height = unit(4, "cm")),
              col = col_fun)

c <- grid.grabExpr({
  draw(ht, padding = unit(c(1, 1, 1, 7.5), "mm"))  # Draw heatmap only
})

## panel d — diseases first diagnosed during follow-up (incident only), in years
mat <- read_csv(out_path("processed", "incident_disease_years.csv"), show_col_types = FALSE)
mat$disease <- factor(mat$disease, levels = mat %>% group_by(disease) %>% summarise(mean = mean(mean_years)) %>% arrange(desc(mean)) %>% pull(disease) %>% rev)
mat <- mat %>% left_join(mat %>% group_by(domain) %>% summarise(mean_domain = mean(mean_years)), by = "domain")

d <- ggplot(mat, aes(disease, mean_years, fill = domain, color = domain)) +
  geom_hline(yintercept = 0, color = "grey70", lwd = 0.3) +
  geom_hline(aes(yintercept = mean_domain, color = domain), lty = 3, lwd = 0.5) +
  geom_point(data = mat %>% filter(!sig), shape = 16, size = 3, show.legend = TRUE, alpha = 0.7) +
  geom_point(data = mat %>% filter(sig), shape = 21, size = 3, stroke = 1, color = "black", show.legend = TRUE, alpha = 0.7) +
  coord_flip() +
  theme_pubr() +
  theme(legend.position = "right") +
  scale_fill_manual(values = c(colors$red[tone],colors$blue[tone],colors$teal[tone],colors$orange[tone],colors$yellow[tone]))+
  scale_color_manual(values = c(colors$red[tone],colors$blue[tone],colors$teal[tone],colors$orange[tone],colors$yellow[tone]))+
  labs(x = "Disease (incident during follow-up)", y = "Mean functional-age acceleration (years)", fill = "Domain") +
  guides(size = "none", stroke = "none", color = "none", fill = guide_legend(override.aes = list(stroke = 0)))

#panel e
all_fields <- read_tsv(in_path("field_metadata")) %>% mutate(field_id = as.character(field_id))
subset <- read_rds(out_path("processed", "lifestyle.rds"))
subset <- subset[, !colnames(subset) %in% lifestyle_exclude]
results <- read_rds(out_path("processed", "age.rds"))

correlations <- tibble()
for (col in 2:ncol(subset)) {
  temp <- results %>% dplyr::select(domain,eid,functional_acc,age,sex) %>% left_join(subset[,c(1,col)]) %>% na.omit 
  field_id <- colnames(temp)[6]
  temp <- temp %>% set_names("domain","eid","functional_acc","age","sex","value") %>% group_by(domain) %>% summarise(r = cor.test(functional_acc,value)$estimate) %>% mutate(field_id = field_id)
  correlations <- rbind(correlations,temp)
}
correlations <- correlations %>% left_join(all_fields)
correlations$title[correlations$field_id=="Number of age-related diseases"] <- "Number of age-related diseases"
correlations$title[correlations$field_id=="Healthspan"] <- "Healthspan"
correlations$domain <- str_to_sentence(correlations$domain)

col_fun <- colorRamp2(seq(-0.5,0.5,0.1),brewer.pal(11, "RdBu") %>% rev)

mat <- reshape2::acast(correlations, title~domain, value.var = "r")
mat <- mat[rowSums(abs(mat)>=0.05) %>% enframe %>% filter(value>0) %>% pull(name),]
mat <- mat[rowMeans(mat) %>% enframe %>% arrange(desc(value)) %>% pull(name),]
mat <- t(mat)
colnames(mat) <- ifelse(colnames(mat) %in% names(lifestyle_labels), lifestyle_labels[colnames(mat)], colnames(mat))

sq_legend = Legend(
  labels = c("> 0.3", "0.2 - 0.3", "0.1 - 0.2", "0.05 - 0.1"),
  title = "|r|",
  type = "lines",        # Changed from "points" to "lines"
  legend_gp = gpar(
    col = "black", 
    lwd = c(4, 3, 2, 1)  # Matches your cell_fun logic
  ),
  grid_width = unit(6, "mm") # Adjusts the length of the legend lines
)

ht <- Heatmap(mat,
              border = TRUE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                if (abs(mat[i,j])>0.3) {
                  grid.rect(x = x, y = y, width = width*0.8, height = height*0.8, gp = gpar(col = "black", fill = NA, lty = 1, lwd = 4))
                } else if (abs(mat[i,j])>0.2) {
                  grid.rect(x = x, y = y, width = width*0.8, height = height*0.8, gp = gpar(col = "black", fill = NA, lty = 1, lwd = 3))
                } else if (abs(mat[i,j])>0.1) {
                  grid.rect(x = x, y = y, width = width*0.8, height = height*0.8, gp = gpar(col = "black", fill = NA, lty = 1, lwd = 2))
                } else if (abs(mat[i,j])>0.05) {
                  grid.rect(x = x, y = y, width = width*0.8, height = height*0.8, gp = gpar(col = "black", fill = NA, lty = 1, lwd = 1))
                }
              },
              #row_title = "Lifestyle factor",
              column_title = "",
              cluster_columns = FALSE,
              row_dend_side = "right",
              column_names_gp = gpar(fontsize = 10),
              #column_names_rot = 75,
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
              heatmap_legend_param = list(direction = "vertical", title = "R", title_position = "topcenter", legend_height = unit(4, "cm")),
              col = col_fun)

e <- grid.grabExpr({
  draw(ht, padding = unit(c(52, 1, 5, 5), "mm"),
  annotation_legend_list = list(sq_legend),
  merge_legend = TRUE) 
})

library(cowplot)
pdf(file = out_path("figures", "Figure1.pdf"), width = 15, height = 13.4)
plot_grid(plot_grid(a, nrow = 1, labels = c("a")),
          plot_grid(b,c,d, nrow = 1, labels = c("b","c","d"), rel_widths = c(1,1.1,1.4)), 
          plot_grid(e, nrow = 1, labels = c("e")), 
          nrow = 3, rel_heights = c(1,1.1,1.6))
dev.off()


