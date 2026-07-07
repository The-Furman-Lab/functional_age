# Build the survival table: follow-up time, vital status, sex and exact age.
# in : mortality (date/age/cause of death), phenotype table (assessment date,
#      year/month of birth, sex), withdrawals.
# out: output/<visit_mode>/processed/survival.rds
library(glmnet)
library(tidyverse)
library(data.table)
library(survival)
library(impute)
library(MatchIt)
library(flexsurv)

source("code/_config.R")

mortality <- read_csv(in_path("mortality"))
mortality <- mortality[,c("eid","p40000_i0","p40007_i0","p40001_i0")] %>% set_names("eid","40000","40007","40001") #eid, date_death, age_death, cause_death
mortality$`40001` <- sapply(mortality$`40001`, function(x) strsplit(x, split = " ")[[1]][1]) %>% as.character()
mortality$`40001` <- gsub("\\.","",mortality$`40001`)
mortality$`40001`[(!is.na(mortality$`40000`)&is.na(mortality$`40001`))] <- "XXXX" #unknown causes of death

survival_param <- c("53-0.0","34-0.0","52-0.0","31-0.0")
selected <- c("eid",survival_param)
subset <- fread(in_path("phenotype_table"), select = selected) %>% data.frame(check.names = FALSE) 
colnames(subset) <- gsub("-0.0","",colnames(subset)) #remove instance number
subset <- mortality %>% left_join(subset, by = c("eid")) #combine mortality and features

subset$follow_time_dead <- as.numeric(difftime(subset$`40000`, subset$`53`, units = "days"))/365.25 # date_death - date_assessment  = follow_up_time
subset$follow_time_dead[is.na(subset$follow_time_dead)] <- 0

date_censor <- max(subset$`40000`, na.rm = TRUE) #max date
subset$follow_time_alive <- as.numeric(difftime(date_censor, subset$`53`, units = "days"))/365.25 # date_censor - date_assessment = follow_up_time
subset$status <- ifelse(subset$follow_time_dead==0, 0, 1)

subset$time <- ifelse(subset$follow_time_dead==0, subset$follow_time_alive, subset$follow_time_dead)
subset$cod <- subset$`40001` # add cause of death
subset$sex <- subset$`31`
subset$age <- as.numeric(difftime(subset$`53`, paste0(subset$`34`,"-",sprintf("%02d", subset$`52`),"-","01"), units = "days")/365.25) #exact age

remove <- read_csv(in_path("withdrawals"), col_names = FALSE)
survival <- subset %>% dplyr::select(eid,status,time,sex,age) %>% na.omit %>% filter(!eid%in%remove$X1)
write_rds(survival, file = out_path("processed", "survival.rds"))
