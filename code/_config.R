# Central configuration for the IC pipeline.
#
# All input/output locations and run parameters live in config.yml at the repo
# root, so the analysis scripts carry no machine-specific paths. Run the scripts
# from the repository root, e.g.  Rscript code/00_tidy_domains.R 1
#
# in_path(key)  -> an input location, looked up from config.yml
# out_path(...) -> an output location under output_root/<visit_mode>/...,
#                  creating the parent directory if needed
#
# The large phenotype table can be kept outside the repo by setting the
# environment variable IC_PHENOTYPE_TABLE, without editing config.yml.

library(yaml)

cfg <- yaml::read_yaml(Sys.getenv("IC_CONFIG", "config.yml"))$default

.pheno <- Sys.getenv("IC_PHENOTYPE_TABLE", "")
if (nzchar(.pheno)) cfg$phenotype_table <- .pheno

in_path <- function(key) {
  p <- cfg[[key]]
  if (is.null(p)) stop("unknown config key: ", key)
  p
}

out_path <- function(...) {
  p <- file.path(cfg$output_root, cfg$visit_mode, ...)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}

colors <- list(
  tone   = c("#f6f2ed","#dedbca","#c4c0a5","#a59f83","#878264","#5e5948"),
  grey   = c("#e5e4e8","#c5c9d6","#959fb3","#6e788c","#425468","#1b2942"),
  olive  = c("#f2edb2","#dbdc63","#c4c400","#95a008","#637314","#304215"),
  green  = c("#d7e5c5","#9fc978","#5db342","#41912f","#1d6e29","#0e3716"),
  teal   = c("#c9e4ef","#96ced3","#48bcbc","#00959f","#006479","#003547"),
  blue   = c("#c5e4fb","#9bc9e8","#5495ce","#006eae","#01478c","#002259"),
  purple = c("#e9d3e7","#d1a9ce","#b678b3","#a44990","#792373","#430b4d"),
  red    = c("#f5cfc9","#e9a0a4","#db6463","#c5373c","#9c241b","#730c0d"),
  orange = c("#fbdcbc","#f9bd7b","#f29741","#e96700","#b34a00","#832a00"),
  yellow = c("#ffedc1","#f7dc86","#e8c54d","#ca9a23","#9b730a","#685409"),
  skin   = c("#f6e5d3","#dcbc9f","#bc9778","#906852","#734d3d","#422a17"))

# Domain features excluded from the mortality models: cognitive tests administered
# only at a follow-up visit (introduce healthy-volunteer selection bias) and
# psychological help-seeking items (outcome-related rather than functional).
domain_feature_exclude <- c(
  "Duration to complete alphanumeric path (trail #2)",
  "Number of symbol digit matches made correctly",
  "Seen a psychiatrist for nerves, anxiety, tension or depression",
  "Seen doctor (GP) for nerves, anxiety, tension or depression")

# Lifestyle factors excluded from Fig 1e, and display labels for the retained ones.
lifestyle_exclude <- c("2188","1120","1110","1458","1408","1478","1548","2237","24508","1050","1060","24501","24504",
                       "Number of age-related diseases","Healthspan")
lifestyle_labels <- c(
  "Pack years adult smoking as proportion of life span exposed to smoking" = "Smoking pack-years (% lifespan)",
  "Pack years of smoking"                                    = "Smoking pack-years",
  "Townsend deprivation index at recruitment"                = "Townsend deprivation index",
  "Someone to take to doctor when needed as a child"         = "Someone to take to doctor as child",
  "Morning/evening person (chronotype)"                      = "Evening person (Chronotype)",
  "Nitrogen dioxide air pollution; 2010"                     = "NO2 air pollution",
  "Nitrogen oxides air pollution; 2010"                      = "NOx air pollution",
  "Particulate matter air pollution (pm2.5); 2010"           = "PM2.5 air pollution",
  "Particulate matter air pollution (pm2.5) absorbance; 2010"= "PM2.5 absorbance",
  "Particulate matter air pollution (pm10); 2010"            = "PM10 air pollution",
  "Job involves heavy manual or physical work"               = "Job: heavy physical",
  "Job involves mainly walking or standing"                  = "Job: walking/standing",
  "Job involves shift work"                                  = "Job: shift work",
  "Job involves night shift work"                            = "Job: night shift",
  "Age started smoking in former smokers"                    = "Age started smoking (former)",
  "Age started smoking in current smokers"                   = "Age started smoking (current)",
  "Age completed full time education"                        = "Age completed education",
  "Age when periods started (menarche)"                      = "Age at menarche",
  "Age at menopause (last menstrual period)"                 = "Age at menopause",
  "Age at first live birth"                                  = "Age at first birth",
  "Ever had stillbirth, spontaneous miscarriage or termination" = "Pregnancy loss (ever)",
  "Natural environment percentage, buffer 1000m"             = "Natural environment % within 1000m",
  "Natural environment percentage, buffer 300m"              = "Natural environment % within 300m")
