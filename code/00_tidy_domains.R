library(glmnet)
library(tidyverse)
library(data.table)
library(survival)

source("code/_config.R")
args = commandArgs(trailingOnly=TRUE)
id <- args[1] %>% as.numeric()

column_names <- fread(in_path("phenotype_table"),  nrows = 1) %>% data.frame(check.names = FALSE) %>% colnames()
all_fields <- read_tsv(in_path("field_metadata"))
min_cutoff <- 10000

#cognition
if (id==1) {
  selected <- c("20023", #reaction time
                "4282","20240", #numeric memory
                "20016","20191", #fluid intelligence
                "6350","20157", #trial making
                "6373","20760", #matrix pattern completition
                "23324","20159", #symbol digit substitution
                "20139", #broken letter recognition
                "26302", #picture vocabulary
                "21004", #tower rearranging
                "20197", #paired associate learning
                "20018", #prospective memory
                "399","20132") #pairs matching
  selected <- intersect(all_fields %>% filter(num_participants>=min_cutoff) %>% pull(field_id),selected)
  fields <- lapply(selected, function(x) column_names[grepl(paste0("^",x,"-"), column_names)])
  names(fields) <- selected
  fields <- fields %>% enframe %>% unnest(value) %>% set_names("field_id","exact_field") %>% left_join(all_fields[,1:2] %>% mutate(field_id = as.character(field_id))) %>% dplyr::select(field_id,title,exact_field)
  subset <- fread(in_path("phenotype_table"), select = c("eid",fields$exact_field)) %>% data.frame(check.names = FALSE) %>% mutate(eid = as.character(eid)) 
  subset <- subset %>% mutate(across(matches("4282-"), ~ ifelse(. == -1, NA, .)))
  subset <- subset %>% mutate(across(matches("20016-"), ~ ifelse(. == -1, NA, .)))
  subset <- subset %>% mutate(across(matches("20191-"), ~ ifelse(. == -1, NA, .)))
  subset <- subset %>% mutate(across(matches("6350-"), ~ ifelse(. == 0, NA, .)))
  subset <- subset %>% mutate(across(matches("20018-"), ~ ifelse(. == 2, 0, .)))
  data <- reshape2::melt(subset %>% column_to_rownames(var = "eid") %>% as.matrix) %>% na.omit %>% set_names("eid","exact_field","value") %>% left_join(fields[,c("exact_field","title")])
  data$instance <- sapply(data$exact_field, function(x) strsplit(x, split = "\\-|\\.")[[1]][2])
  data$array <- sapply(data$exact_field, function(x) strsplit(x, split = "\\-|\\.")[[1]][3])
  data <- data %>% group_by(eid,title,instance) %>% summarise(value = median(value, na.rm = TRUE))
  data <- data %>% group_by(eid,title) %>% summarise(value = median(value, na.rm = TRUE))
  cognitive <- reshape2::acast(data, eid~title, value.var = "value") 
  write_rds(cognitive, file = out_path("processed", "domains", "cognitive.rds"))
}

#vitality
if (id==2) {
  selected <- c("3062","3063","3064",#spirometry
                "46","47",#grip strength
                "30020",#hemoglobin, biomarker
                "30770",#IGF-1, biomaker
                "21001","48",#anthropometry
                "2306")#weight change
  selected <- intersect(all_fields %>% filter(num_participants>=min_cutoff) %>% pull(field_id),selected)
  fields <- lapply(selected, function(x) column_names[grepl(paste0("^",x,"-"), column_names)])
  names(fields) <- selected
  fields <- fields %>% enframe %>% unnest(value) %>% set_names("field_id","exact_field") %>% left_join(all_fields[,1:2] %>% mutate(field_id = as.character(field_id))) %>% dplyr::select(field_id,title,exact_field)
  subset <- fread(in_path("phenotype_table"), select = c("eid",fields$exact_field)) %>% data.frame(check.names = FALSE) %>% mutate(eid = as.character(eid)) 
  subset <- subset %>% mutate(across(matches("2306-"), ~ ifelse(. %in%c(-1,-3), NA, .)))
  subset <- subset %>% mutate(across(matches("2306-"), ~ ifelse(. == 0, 1, .)))
  subset <- subset %>% mutate(across(matches("2306-"), ~ ifelse(. == 3, 0, .)))
  data <- reshape2::melt(subset %>% column_to_rownames(var = "eid") %>% as.matrix) %>% na.omit %>% set_names("eid","exact_field","value") %>% left_join(fields[,c("exact_field","title")])
  data$instance <- sapply(data$exact_field, function(x) strsplit(x, split = "\\-|\\.")[[1]][2])
  data$array <- sapply(data$exact_field, function(x) strsplit(x, split = "\\-|\\.")[[1]][3])
  data <- data %>% group_by(eid,title,instance) %>% summarise(value = median(value, na.rm = TRUE))
  data <- data %>% group_by(eid,title) %>% summarise(value = median(value, na.rm = TRUE))
  data$title[data$title=="Hand grip strength (left)"] <- "Hand grip strength"
  data$title[data$title=="Hand grip strength (right)"] <- "Hand grip strength"
  data <- data %>% group_by(eid,title) %>% summarise(value = median(value, na.rm = TRUE))
  vitality <- reshape2::acast(data, eid~title, value.var = "value") 
  write_rds(vitality, file = out_path("processed", "domains", "vitality.rds"))
}

#locomotion
if (id==3) {
  selected <- c()
  #selected <- c("2296","2634","1021","894","3647","914","874","981","2624","1011","3637","943","971","884","904","864","924","22035","22036","22032","22038","22039","22037","22040","22033","22034")
  selected <- c("22035","22036","22032","22038","22039","22037","22040","22033","22034",
                "2296",#falls
                #"2634","1021","894","3647","914","874","981","884","904","864",#physical activity "1001" small sample size
                #"2624","1011","3637","943","971",#frequency activity last 4 weeks "991" small sample size
                #"1090","1080","1070",#time spend sitting
                "924") #walking pace
  selected <- intersect(all_fields %>% filter(num_participants>=min_cutoff) %>% pull(field_id),selected)
  fields <- lapply(selected, function(x) column_names[grepl(paste0("^",x,"-"), column_names)])
  names(fields) <- selected
  fields <- fields %>% enframe %>% unnest(value) %>% set_names("field_id","exact_field") %>% left_join(all_fields[,1:2] %>% mutate(field_id = as.character(field_id))) %>% dplyr::select(field_id,title,exact_field)
  subset <- fread(in_path("phenotype_table"), select = c("eid",fields$exact_field)) %>% data.frame(check.names = FALSE) %>% mutate(eid = as.character(eid)) 
  subset <- subset %>% mutate(across(matches("864-"), ~ ifelse(. == -2, 0, .)))
  subset <- subset %>% mutate(across(matches("1070-"), ~ ifelse(. == -10, 0, .)))
  subset <- subset %>% mutate(across(matches("1080-"), ~ ifelse(. == -10, 0, .)))
  subset <- subset %>% mutate(across(matches("1090-"), ~ ifelse(. == -10, 0, .)))
  subset[subset==-1] <- NA
  subset[subset==-3] <- NA
  subset[subset==-7] <- NA
  data <- reshape2::melt(subset %>% column_to_rownames(var = "eid") %>% as.matrix) %>% na.omit %>% set_names("eid","exact_field","value") %>% left_join(fields[,c("exact_field","title")])
  data$instance <- sapply(data$exact_field, function(x) strsplit(x, split = "\\-|\\.")[[1]][2])
  data$array <- sapply(data$exact_field, function(x) strsplit(x, split = "\\-|\\.")[[1]][3])
  data <- data %>% group_by(eid,title,instance) %>% summarise(value = median(value, na.rm = TRUE))
  data <- data %>% group_by(eid,title) %>% summarise(value = median(value, na.rm = TRUE))
  locomotion <- reshape2::acast(data, eid~title, value.var = "value") 
  write_rds(locomotion, file = out_path("processed", "domains", "locomotion.rds"))
}

#sensory
if (id==4) {
  selected <- c("20021","20019","2247","3393", #hearing
                "5201","5208","6148","2207") #audition
  selected <- intersect(all_fields %>% filter(num_participants>=min_cutoff) %>% pull(field_id),selected)
  fields <- lapply(selected, function(x) column_names[grepl(paste0("^",x,"-"), column_names)])
  names(fields) <- selected
  fields <- fields %>% enframe %>% unnest(value) %>% set_names("field_id","exact_field") %>% left_join(all_fields[,1:2] %>% mutate(field_id = as.character(field_id))) %>% dplyr::select(field_id,title,exact_field)
  subset <- fread(in_path("phenotype_table"), select = c("eid",fields$exact_field)) %>% data.frame(check.names = FALSE) %>% mutate(eid = as.character(eid)) 
  subset <- subset %>% mutate(across(matches("2247-"), ~ ifelse(. == 99, 1, .)))
  subset <- subset %>% mutate(across(matches("2247-"), ~ ifelse(. %in%c(-1,-3), NA, .)))
  subset <- subset %>% mutate(across(matches("3393-"), ~ ifelse(. == -3, NA, .)))
  subset <- subset %>% mutate(across(matches("6148-"), ~ ifelse(. %in%c(-1,-3), NA, .)))
  subset <- subset %>% mutate(across(matches("6148-"), ~ ifelse(. %in%c(1:6), 1, .)))
  subset <- subset %>% mutate(across(matches("6148-"), ~ ifelse(. == -7, 0, .)))
  data <- reshape2::melt(subset %>% column_to_rownames(var = "eid") %>% as.matrix) %>% na.omit %>% set_names("eid","exact_field","value") %>% left_join(fields[,c("exact_field","title")])
  data$instance <- sapply(data$exact_field, function(x) strsplit(x, split = "\\-|\\.")[[1]][2])
  data$array <- sapply(data$exact_field, function(x) strsplit(x, split = "\\-|\\.")[[1]][3])
  data <- data %>% group_by(eid,title,instance) %>% summarise(value = median(value, na.rm = TRUE))
  data <- data %>% group_by(eid,title) %>% summarise(value = median(value, na.rm = TRUE))
  sensory <- reshape2::acast(data, eid~title, value.var = "value") 
  write_rds(sensory, file = out_path("processed", "domains", "sensory.rds"))
}

#psychological
if (id==5) {
  selected <- all_fields %>% filter(main_category=="100060") %>% pull(field_id)
  selected <- intersect(all_fields %>% filter(num_participants>=min_cutoff) %>% pull(field_id),selected)
  fields <- lapply(selected, function(x) column_names[grepl(paste0("^",x,"-"), column_names)])
  names(fields) <- selected
  fields <- fields %>% enframe %>% unnest(value) %>% set_names("field_id","exact_field") %>% left_join(all_fields[,1:2] %>% mutate(field_id = as.character(field_id))) %>% dplyr::select(field_id,title,exact_field)
  subset <- fread(in_path("phenotype_table"), select = c("eid",fields$exact_field)) %>% data.frame(check.names = FALSE) %>% mutate(eid = as.character(eid)) 
  subset[subset == -1 | subset == -3] <- NA
  data <- reshape2::melt(subset %>% column_to_rownames(var = "eid") %>% as.matrix) %>% na.omit %>% set_names("eid","exact_field","value") %>% left_join(fields[,c("exact_field","title")])
  data$instance <- sapply(data$exact_field, function(x) strsplit(x, split = "\\-|\\.")[[1]][2])
  data$array <- sapply(data$exact_field, function(x) strsplit(x, split = "\\-|\\.")[[1]][3])
  data <- data %>% group_by(eid,title,instance) %>% summarise(value = median(value, na.rm = TRUE))
  data <- data %>% group_by(eid,title) %>% summarise(value = median(value, na.rm = TRUE))
  psychological <- reshape2::acast(data, eid~title, value.var = "value") 
  write_rds(psychological, file = out_path("processed", "domains", "psychological.rds"))
}

##diseases
column_names <- fread(in_path("phenotype_table"),  nrows = 1) %>% data.frame(check.names = FALSE) %>% colnames()
fields <- readxl::read_xlsx(in_path("tian_fields"), sheet = 7) %>% mutate(field_id = as.character(field_id))

#self reported
self <- fields[grepl("UK Biobank Self Report", fields$type),]
subset <- fread(in_path("phenotype_table"), select = c("eid",column_names[grepl("^20002-", column_names)])) %>%
  data.frame(check.names = FALSE) %>% mutate(eid = as.character(eid))
age_disease <- fread(in_path("phenotype_table"), select = c("eid",column_names[grepl("^20009-", column_names)])) %>%
  data.frame(check.names = FALSE) %>% mutate(eid = as.character(eid))
age_disease <- reshape2::melt(age_disease, id.vars = "eid") %>% na.omit %>% dplyr::select(eid,variable,value) %>% set_names("eid","col","age") %>%
  mutate(col = gsub("20009-0\\.", "", col))
self <- reshape2::melt(subset, id.vars = "eid") %>% na.omit %>% dplyr::select(eid,variable,value) %>% filter(value%in%self$field_id) %>%
  set_names("eid","col","field_id") %>%
  mutate(col = gsub("20002-0\\.", "", col)) %>%
  mutate(field_id = as.character(field_id)) %>% left_join(self) %>% left_join(age_disease)
healthspan <- self %>% filter(age>0) %>% group_by(eid) %>% summarise(healthspan = min(age))
write_rds(healthspan, file = out_path("processed", "healthspan.rds"))
self <- self %>% dplyr::select(eid,disease) %>% unique

#icd-9
icd9 <- fields %>% filter(type=="ICD 9")
subset <- fread(in_path("phenotype_table"), select = c("eid",column_names[grepl("^41271-", column_names)])) %>%
  data.frame(check.names = FALSE) %>% mutate(eid = as.character(eid))
subset[subset==""] <- NA
icd9 <- reshape2::melt(subset, id.vars = "eid") %>% na.omit %>% dplyr::select(eid,value) %>% filter(value%in%icd9$field_id) %>%
  set_names("eid","field_id") %>% mutate(field_id = as.character(field_id)) %>% left_join(icd9)
icd9 <- icd9 %>% dplyr::select(eid,disease) %>% unique

#icd-10
icd10 <- fields %>% filter(type=="ICD 10")
subset <- fread(in_path("phenotype_table"), select = c("eid",column_names[grepl("^41270-", column_names)])) %>%
  data.frame(check.names = FALSE) %>% mutate(eid = as.character(eid))
subset[subset==""] <- NA
icd10 <- reshape2::melt(subset, id.vars = "eid") %>% na.omit %>% dplyr::select(eid,value) %>% filter(value%in%icd10$field_id) %>%
  set_names("eid","field_id") %>% mutate(field_id = as.character(field_id)) %>% left_join(icd10)
icd10 <- icd10 %>% dplyr::select(eid,disease) %>% unique

diseases <- rbind(self,icd9,icd10)
write_rds(diseases, file = out_path("processed", "diseases.rds"))

multimorbidity <- diseases %>% dplyr::select(eid,disease) %>%  unique %>% group_by(eid) %>% summarise(n = length(unique(disease)))
write_rds(multimorbidity, file = out_path("processed", "multimorbidity.rds"))

#lifestyle
all_fields <- read_tsv(in_path("field_metadata")) %>% mutate(field_id = as.character(field_id))
column_names <- fread(in_path("phenotype_table"),  nrows = 1) %>% data.frame(check.names = FALSE) %>% colnames()

exclude <- do.call("rbind", lapply(list.files(out_path("models", "clinical"), full.names = TRUE), function(x) read_rds(x)))
exclude <- exclude %>% dplyr::select(domain,coef) %>% unnest(coef)
exclude <- all_fields$field_id[all_fields$title%in%unique(gsub("`","",exclude$name))]
exclude <- c(exclude,all_fields %>% filter(main_category%in%c("54","100054")) %>% pull(field_id))
fields <- readxl::read_xlsx(in_path("tian_fields"), sheet = 6, skip = 2) %>% filter(!`Field ID`%in%c("6138-0.0","1080-0.0","6142-0.0","6143-0.0","6144-0.0"))
fields <- fields[grepl("-0.0", fields$`Field ID`),1:2] %>% set_names("field_id","title")
fields$field_id[fields$field_id=="189-0.0"] <- "22189-0.0"
fields <- rbind(fields, tibble(field_id="1807-0.0",title="Father's age at death"))
fields <- rbind(fields, tibble(field_id="3526-0.0",title="Mother's age at death"))
fields <- fields %>% filter(!field_id%in%paste0(exclude,"-0.0"))

#check
#qualifications
array1 <- fread(in_path("phenotype_table"), select = c("eid",column_names[grepl("^6138-0",column_names)])) %>% data.frame(check.names = FALSE) %>% mutate(eid = as.character(eid))
array1[array1<0] <- NA
literacy_map <- c(
  "1" = 6, # Degree
  "6" = 5, # Professional
  "5" = 4, # NVQ/HND/HNC
  "2" = 3, # A Levels
  "3" = 2, # GCSE/O Level
  "4" = 1  # CSE
)
df_final <- array1 %>%
  mutate(across(starts_with("6138"), ~ literacy_map[as.character(.)])) %>%
  mutate(max_literacy = do.call(pmax, c(dplyr::select(., starts_with("6138")), na.rm = TRUE)))
array1 <- tibble(eid = array1$eid, `6138` = df_final$max_literacy)
array1$`6138`[is.infinite(array1$`6138`)] <- NA

#all others
subset <- fread(in_path("phenotype_table"), select = c("eid",unique(fields$field_id))) %>%
  data.frame(check.names = FALSE) %>% mutate(eid = as.character(eid))

subset <- subset[,c("eid",table(fields$field_id) %>% enframe %>% filter(value<=2) %>% pull(name))]
colnames(subset) <- gsub("-0.0","",colnames(subset))
subset[subset<0] <- NA
subset <- subset[,!colnames(subset)%in%exclude]
subset$`1031` <- 7-subset$`1031`
subset$`1210`[subset$`1210`==2] <- 0
subset$`1249` <- 4-subset$`1249`
subset <- subset %>% mutate(`1239` = recode(`1239`, `1` = 2, `2` = 1, `0` = 0))
subset <- subset %>% mutate(`1687` = recode(`1687`, `3` = 2, `2` = 3, `1` = 1))
subset <- subset %>% mutate(`1697` = recode(`1697`, `3` = 2, `2` = 3, `1` = 1))
subset <- subset %>% mutate(`1757` = recode(`1757`, `3` = 2, `2` = 3, `1` = 1))
subset$`2724`[subset$`2724`==2] <- NA
subset$`2724`[subset$`2724`==3] <- NA
subset$`26427` <- NULL
subset$`26410` <- NULL
subset$`26426` <- NULL
subset <- subset %>% left_join(read_rds(out_path("processed", "multimorbidity.rds")) %>% set_names("eid","Number of age-related diseases"))
subset$`Number of age-related diseases`[is.na(subset$`Number of age-related diseases`)] <- 0
subset <- subset %>% left_join(read_rds(out_path("processed", "healthspan.rds")) %>% set_names("eid","Healthspan"))
subset <- subset %>% left_join(array1)
write_rds(subset, file = out_path("processed", "lifestyle.rds"))

#medication
fields <- c("6177","6153", "6154", "10004", "10005", "6155", "10007", "6179", "10723", "10854")
fields <- lapply(fields, function(x) paste0(x, paste0("-0.",0:10))) %>% unlist
subset <- fread(in_path("phenotype_table"), select = c("eid",fields)) %>%
  data.frame(check.names = FALSE) %>% mutate(eid = as.character(eid))

codes <- rbind(tibble(field = "6177", code = c(1:3, -7), meaning = c("Cholesterol lowering medication","Blood pressure medication","Insulin","None of the above")),
            tibble(field = "6153", code = c(1:5, -7), meaning = c("Cholesterol lowering medication","Blood pressure medication","Insulin","Hormone replacement therapy","Oral contraceptive pill or minipill","None of the above")),
            tibble(field = "6154", code = c(1:6, -7), meaning = c("Aspirin","Ibuprofen","Paracetamol","Ranitidine","Omeprazole","Laxatives","None of the above")),
            tibble(field = "10004", code = c(1:5, -7), meaning = c("Aspirin","Ibuprofen","Paracetamol","Codeine","Ranitidine","None of the above")),
            tibble(field = "10005", code = c(1:4, -7), meaning = c("Omeprazole","Laxatives","Nicotine","Antihistamines","None of the above")),
            tibble(field = "6155", code = c(1:7, -7), meaning = c("Vitamin A","Vitamin B","Vitamin C","Vitamin D","Vitamin E","Folic acid" ,"Multivitamins","None of the above")),
            tibble(field = "10007", code = c(1:6, -7), meaning = c("Evening primrose oil","Fish oil","Garlic","Ginkgo","Glucosamine","Other supplements","None of the above")),
            tibble(field = "6179", code = c(1:6, -7), meaning = c("Fish oil","Glucosamine","Calcium","Zinc","Iron","Selenium", "None of the above")),
            tibble(field = "10723", code = c(1:6, -7), meaning = c("Vitamin A","Vitamin B","Vitamin C","Vitamin D","Vitamin E","Folic acid","None of the above")),
            tibble(field = "10854", code = c(1:6, -7), meaning = c("Iron","Zinc","Calcium","Selenium","Multivitamins","Multivitamins","None of the above")))

medication <- reshape2::melt(subset) %>% na.omit %>% filter(!value%in%c("-1","-3")) %>% set_names("eid","field","code") %>% mutate(field = gsub("-0.0|-0.1|-0.2|-0.3|-0.4|-0.5|-0.6|-0.7|-0.8|-0.9|-0.10","",field)) %>% left_join(codes) %>% na.omit
medication <- medication %>% filter(meaning!="None of the above") %>% group_by(meaning, field) %>% summarise(cases = list(unique(eid))) %>% 
  left_join(medication %>% filter(meaning=="None of the above") %>% group_by(field) %>% summarise(controls = list(unique(eid))))

medication <- medication %>%
  pivot_longer(
    cols = c(cases, controls),
    names_to = "group",
    values_to = "eid"
  ) %>%
  unnest(eid) %>% dplyr::select(meaning,group,eid) %>% unique

medication <- medication %>%
  group_by(meaning, eid) %>%
  filter(!any(group == "cases") | !any(group == "controls")) %>%
  ungroup()

write_rds(medication, file = out_path("processed", "medication.rds"))

#drugs 
coding <- read_tsv(file.path(in_path("coding_dir"), "coding4.tsv"))
fields <- c("20003")
fields <- paste0(fields, "-0.",0:48)
subset <- fread(in_path("phenotype_table"), select = c("eid",fields)) %>%
  data.frame(check.names = FALSE) %>% mutate(eid = as.character(eid))
drugs <- reshape2::melt(subset) %>% na.omit %>% dplyr::select(eid,value) %>% set_names("eid","coding") %>% left_join(coding) %>% na.omit
no_drugs <- subset$eid[!subset$eid%in%drugs$eid]
drugs <- drugs %>% group_by(meaning) %>% summarise(cases = list(eid), n = n()) %>% filter(n>1000) %>% filter(meaning!="Free-text entry, unable to be coded")
drugs$n <- NULL
drugs <- drugs %>% mutate(controls = list(no_drugs))
write_rds(drugs, file = out_path("processed", "drugs.rds"))
