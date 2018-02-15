cat("===========================================================================================
                META-ANALYSIS ATHERO-EXPRESS METHYLATION STUDIES 450K 1 & 2
            -- Linear & Logistic regressions of Smoking vs. Plaque Phenotypes --
    Version:      v2.8
    
    Last update:  2018-01-31
    Written by:   Sander W. van der Laan (s.w.vanderlaan-2@umcutrecht.nl);
    Marten A. Siemelink
    
    Description:  Script to do meta-analysis of Athero-Express Methylation 
    Study 450K 1 (2013) & 2 (2016) EWAS on smoking.
    Based on DNAmArray: 
    - https://molepi.github.io/DNAmArray_workflow/index.html
    - https://github.com/molepi/DNAmArray
    
    Requirements: R version 3.4.1 (2017-06-30) -- 'Single Candle', Linux CentOS7, Mac OS X El Capitan+
    
    ===========================================================================================")
cat("\n===========================================================================================")
cat("CLEAR THE BOARD")
rm(list = ls())

cat("\n===========================================================================================")
cat("GENERAL R SETUP")
### FUNCTION TO INSTALL PACKAGES
### This function will automatically check in both CRAN and Bioconductor. This is 
### a function found by Sander W. van der Laan online from @Samir: 
### http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
### 
cat("\n* Creating funxtion to install and load packages...")
install.packages.auto <- function(x) { 
  x <- as.character(substitute(x)) 
  if (isTRUE(x %in% .packages(all.available = TRUE))) { 
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else { 
    # Update installed packages - this may mean a full upgrade of R, which in turn
    # may not be warrented. 
    #update.packages(ask = FALSE) 
    eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE, repos = \"https://cloud.r-project.org/\")", x)))
  }
  if (isTRUE(x %in% .packages(all.available = TRUE))) { 
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else {
    source("http://bioconductor.org/biocLite.R")
    # Update installed packages - this may mean a full upgrade of R, which in turn
    # may not be warrented.
    #biocLite(character(), ask = FALSE) 
    eval(parse(text = sprintf("biocLite(\"%s\")", x)))
    eval(parse(text = sprintf("require(\"%s\")", x)))
  }
}
# In this case I'm keeping track of the various packages, as versions and 
# actual loading of the libraries gave issues before.
cat("\n* General packages...\n")
# for survival analyses
install.packages.auto("survival")
install.packages.auto("survminer")
# for general statistics
install.packages.auto("Hmisc")
install.packages.auto("openxlsx")
install.packages.auto("devtools")
install.packages.auto("dplyr")
install.packages.auto("data.table")
install.packages.auto("tableone")
install.packages.auto("haven")
# for methylation/rna data
install.packages.auto("RMySQL")
install.packages.auto("GenomicFeatures")
install.packages.auto("bumphunter")
install.packages.auto("minfi")
install.packages.auto("SummarizedExperiment")
install.packages.auto("IlluminaHumanMethylation450kmanifest")
install.packages.auto("IlluminaHumanMethylation450kanno.ilmn12.hg19")
install.packages.auto("FDb.InfiniumMethylation.hg19")
install.packages.auto("TxDb.Hsapiens.UCSC.hg19.knownGene")
install.packages.auto("org.Hs.eg.db")
install.packages.auto("AnnotationDbi")
# for plotting
install.packages.auto("pheatmap")
install.packages.auto("qqman")
install.packages.auto("forestplot")
# for meta-analysis
install.packages.auto("meta")
install.packages.auto("bacon")

# The actual DNAmArray package
cat("\n* DNAmArray package...\n")
# Also refer to: 
# - https://molepi.github.io/DNAmArray_workflow/index.html
# - https://github.com/molepi/DNAmArray
# - https://github.com/bbmri-nl/BBMRIomics
library(devtools)
install_github("molepi/DNAmArray", force = FALSE)
library(DNAmArray)
install_github("molepi/omicsPrint", ref = "R3.4", force = FALSE)
library(omicsPrint)
install_github("bbmri-nl/BBMRIomics", subdir = "BBMRIomics", force = FALSE)
library(BBMRIomics)

### Create datestamp
Today = format(as.Date(as.POSIXlt(Sys.time())), "%Y%m%d")
Today.Report = format(as.Date(as.POSIXlt(Sys.time())), "%A, %B %d, %Y")

###	UtrechtSciencePark Colours Scheme
###
### Website to convert HEX to RGB: http://hex.colorrrs.com.
### For some functions you should divide these numbers by 255.
###
###	No.	Color				    HEX		  RGB							          CMYK					          CHR		MAF/INFO
### --------------------------------------------------------------------------------------------------------------------
###	1	  yellow				  #FBB820 (251,184,32)				      (0,26.69,87.25,1.57) 	  =>	1 		or 1.0 > INFO
###	2	  gold				    #F59D10 (245,157,16)				      (0,35.92,93.47,3.92) 	  =>	2		
###	3	  salmon				  #E55738 (229,87,56) 				      (0,62.01,75.55,10.2) 	  =>	3 		or 0.05 < MAF < 0.2 or 0.4 < INFO < 0.6
###	4	  darkpink			  #DB003F ((219,0,63)					      (0,100,71.23,14.12) 	  =>	4		
###	5	  lightpink			  #E35493 (227,84,147)				      (0,63,35.24,10.98) 		  =>	5 		or 0.8 < INFO < 1.0
###	6	  pink			  	  #D5267B (213,38,123)				      (0,82.16,42.25,16.47) 	=>	6		
###	7	  hardpink			  #CC0071 (204,0,113)					      (0,0,0,0) 	=>	7		
###	8	  lightpurple		  #A8448A (168,68,138)				      (0,0,0,0) 	=>	8		
###	9	  purple				  #9A3480 (154,52,128)				      (0,0,0,0) 	=>	9		
###	10	lavendel		  	#8D5B9A (141,91,154)				      (0,0,0,0) 	=>	10		
###	11	bluepurple			#705296 (112,82,150)				      (0,0,0,0) 	=>	11		
###	12	purpleblue			#686AA9 (104,106,169)				      (0,0,0,0) 	=>	12		
###	13	lightpurpleblue	#6173AD (97,115,173/101,120,180)	(0,0,0,0) 	=>	13		
###	14	seablue				  #4C81BF (76,129,191)				      (0,0,0,0) 	=>	14		
###	15	skyblue				  #2F8BC9 (47,139,201)				      (0,0,0,0) 	=>	15		
###	16	azurblue			  #1290D9 (18,144,217)				      (0,0,0,0) 	=>	16		 or 0.01 < MAF < 0.05 or 0.2 < INFO < 0.4
###	17	lightazurblue	  #1396D8 (19,150,216)				      (0,0,0,0) 	=>	17		
###	18	greenblue			  #15A6C1 (21,166,193)				      (0,0,0,0) 	=>	18		
###	19	seaweedgreen	  #5EB17F (94,177,127)				      (0,0,0,0) 	=>	19		
###	20	yellowgreen		  #86B833 (134,184,51)				      (0,0,0,0) 	=>	20		
###	21	lightmossgreen  #C5D220 (197,210,32)				      (0,0,0,0) 	=>	21		
###	22	mossgreen			  #9FC228 (159,194,40)				      (0,0,0,0) 	=>	22		or MAF > 0.20 or 0.6 < INFO < 0.8
###	23	lightgreen		  #78B113 (120,177,19)				      (0,0,0,0) 	=>	23/X
###	24	green				    #49A01D (73,160,29)					      (0,0,0,0) 	=>	24/Y
###	25	grey				    #595A5C (89,90,92)					      (0,0,0,0) 	=>	25/XY	or MAF < 0.01 or 0.0 < INFO < 0.2
###	26	lightgrey			  #A2A3A4	(162,163,164)				      (0,0,0,0) 	=> 	26/MT
### 
### ADDITIONAL COLORS
### 27	midgrey				  #D7D8D7
### 28	very lightgrey	#ECECEC
### 29	white				    #FFFFFF
### 30	black				    #000000
### --------------------------------------------------------------------------------------------------------------------

uithof_color = c("#FBB820","#F59D10","#E55738","#DB003F","#E35493","#D5267B",
                 "#CC0071","#A8448A","#9A3480","#8D5B9A","#705296","#686AA9",
                 "#6173AD","#4C81BF","#2F8BC9","#1290D9","#1396D8","#15A6C1",
                 "#5EB17F","#86B833","#C5D220","#9FC228","#78B113","#49A01D",
                 "#595A5C","#A2A3A4", "#D7D8D7", "#ECECEC", "#FFFFFF", "#000000")

uithof_color_legend = c("#FBB820", "#F59D10", "#E55738", "#DB003F", "#E35493",
                        "#D5267B", "#CC0071", "#A8448A", "#9A3480", "#8D5B9A",
                        "#705296", "#686AA9", "#6173AD", "#4C81BF", "#2F8BC9",
                        "#1290D9", "#1396D8", "#15A6C1", "#5EB17F", "#86B833",
                        "#C5D220", "#9FC228", "#78B113", "#49A01D", "#595A5C",
                        "#A2A3A4", "#D7D8D7", "#ECECEC", "#FFFFFF", "#000000")

### ----------------------------------------------------------------------------

cat("===========================================================================================")
cat("\nSETUP ANALYSIS")
# Assess where we are
getwd()
# Set locations
### Operating System Version

### Mac Pro
# ROOT_loc = "/Volumes/EliteProQx2Media"

### MacBook
ROOT_loc = "/Users/swvanderlaan"

### SOME VARIABLES WE NEED DOWN THE LINE
PROJECTDATASET = "AEMS450KMETA"
PROJECTNAME = "metasmoke"
SUBPROJECTNAME1 = "AEMS450K1"
SUBPROJECTNAME2 = "AEMS450K2"
EWAS_trait = "SmokerCurrent" # Phenotype

INP_AE_loc = paste0(ROOT_loc, "/PLINK/_AE_Originals")
INP_AEMS450K1_loc = paste0(INP_AE_loc, "/AEMS450K1")
INP_AEMS450K2_loc = paste0(INP_AE_loc, "/AEMS450K2")

EPIGENETICS_loc = paste0(ROOT_loc, "/PLINK/analyses/epigenetics")
RES_AEMS450K1_loc = paste0(EPIGENETICS_loc, "/AEMS450K1")
RES_AEMS450K2_loc = paste0(EPIGENETICS_loc, "/AEMS450K2")

ifelse(!dir.exists(file.path(EPIGENETICS_loc, "/",PROJECTDATASET)), 
       dir.create(file.path(EPIGENETICS_loc, "/",PROJECTDATASET)), 
       FALSE)
INP_loc = paste0(EPIGENETICS_loc, "/",PROJECTDATASET)

cat("\nCreate a new analysis directory...")
ifelse(!dir.exists(file.path(INP_loc, "/",PROJECTNAME)), 
       dir.create(file.path(INP_loc, "/",PROJECTNAME)), 
       FALSE)
ANALYSIS_loc = paste0(INP_loc,"/",PROJECTNAME)
ifelse(!dir.exists(file.path(ANALYSIS_loc, "/PLOTS")), 
       dir.create(file.path(ANALYSIS_loc, "/PLOTS")), 
       FALSE)
PLOT_loc = paste0(ANALYSIS_loc,"/PLOTS")
ifelse(!dir.exists(file.path(PLOT_loc, "/COX")), 
       dir.create(file.path(PLOT_loc, "/COX")), 
       FALSE)
COX_loc = paste0(PLOT_loc,"/COX")
ifelse(!dir.exists(file.path(PLOT_loc, "/QC")), 
       dir.create(file.path(PLOT_loc, "/QC")), 
       FALSE)
QC_loc = paste0(PLOT_loc,"/QC")
ifelse(!dir.exists(file.path(ANALYSIS_loc, "/OUTPUT")), 
       dir.create(file.path(ANALYSIS_loc, "/OUTPUT")), 
       FALSE)
OUT_loc = paste0(ANALYSIS_loc, "/OUTPUT")

cat("===========================================================================================")
cat("\nLOAD ATHERO-EXPRESS METHYLATION STUDY DATASETS")
setwd(INP_loc)
list.files()

cat("\n  - loading AE database...")
load(paste0(OUT_loc,"/20180207.AEData_2016.update.DF.CEA.TFvar.RData"))

AEData_2016.update.DF.CEA.SampleSize = nrow(AEData_2016.update.DF.CEA.TFvar)

cat("===========================================================================================")
cat("SMOKING vs PLAQUE PHENOTYPES")

# Function to grep data from glm()/lm()
GLM.CON <- function(fit, DATASET, x_name, y){
  cat("Analyzing in dataset '", DATASET ,"' the association of '", x_name ,"' with '", y ,"' .\n")
  if (nrow(summary(fit)$coefficients) == 1) {
    output = c(DATASET, x_name, y, rep(NA,8))
    cat("Model not fitted; probably singular.\n")
  }else {
    cat("Collecting data.\n\n")
    effectsize = summary(fit)$coefficients[2,1]
    SE = summary(fit)$coefficients[2,2]
    OReffect = exp(summary(fit)$coefficients[2,1])
    CI_low = exp(effectsize - 1.96 * SE)
    CI_up = exp(effectsize + 1.96 * SE)
    tvalue = summary(fit)$coefficients[2,3]
    pvalue = summary(fit)$coefficients[2,4]
    R = summary(fit)$r.squared
    R.adj = summary(fit)$adj.r.squared
    sample_size = nrow(model.frame(fit))
    AE_N = AEData_2016.update.DF.CEA.SampleSize
    Perc_Miss = 100 - ((sample_size * 100)/AE_N)
    
    output = c(DATASET, x_name, y, effectsize, SE, OReffect, CI_low, CI_up, tvalue, pvalue, R, R.adj, AE_N, sample_size, Perc_Miss)
    cat("We have collected the following:\n")
    cat("Dataset...................:", DATASET, "\n")
    cat("Score.....................:", x_name, "\n")
    cat("Trait.....................:", y, "\n")
    cat("Effect size...............:", round(effectsize, 6), "\n")
    cat("Standard error............:", round(SE, 6), "\n")
    cat("Odds ratio (effect size)..:", round(OReffect, 3), "\n")
    cat("Lower 95% CI..............:", round(CI_low, 3), "\n")
    cat("Upper 95% CI..............:", round(CI_up, 3), "\n")
    cat("T-value...................:", round(tvalue, 6), "\n")
    cat("P-value...................:", signif(pvalue, 8), "\n")
    cat("R^2.......................:", round(R, 6), "\n")
    cat("Adjusted r^2..............:", round(R.adj, 6), "\n")
    cat("Sample size of AE DB......:", AE_N, "\n")
    cat("Sample size of model......:", sample_size, "\n")
    cat("Missing data %............:", round(Perc_Miss, 6), "\n")
  }
  return(output)
  print(output)
}

GLM.BIN <- function(fit, DATASET, x_name, y){
  cat("Analyzing in dataset '", DATASET ,"' the association of '", x_name ,"' with '", y ,"' ...\n")
  if (nrow(summary(fit)$coefficients) == 1) {
    output = c(DATASET, x_name, y, rep(NA,9))
    cat("Model not fitted; probably singular.\n")
  }else {
    cat("Collecting data...\n")
    effectsize = summary(fit)$coefficients[2,1]
    SE = summary(fit)$coefficients[2,2]
    OReffect = exp(summary(fit)$coefficients[2,1])
    CI_low = exp(effectsize - 1.96 * SE)
    CI_up = exp(effectsize + 1.96 * SE)
    zvalue = summary(fit)$coefficients[2,3]
    pvalue = summary(fit)$coefficients[2,4]
    dev <- fit$deviance
    nullDev <- fit$null.deviance
    modelN <- length(fit$fitted.values)
    R.l <- 1 - dev / nullDev
    R.cs <- 1 - exp(-(nullDev - dev) / modelN)
    R.n <- R.cs / (1 - (exp(-nullDev/modelN)))
    sample_size = nrow(model.frame(fit))
    AE_N = AEData_2016.update.DF.CEA.SampleSize
    Perc_Miss = 100 - ((sample_size * 100)/AE_N)
    
    output = c(DATASET, x_name, y, effectsize, SE, OReffect, CI_low, CI_up, zvalue, pvalue, R.l, R.cs, R.n, AE_N, sample_size, Perc_Miss)
    cat("We have collected the following:\n")
    cat("Dataset...................:", DATASET, "\n")
    cat("Score.....................:", x_name, "\n")
    cat("Trait.....................:", y, "\n")
    cat("Effect size...............:", round(effectsize, 6), "\n")
    cat("Standard error............:", round(SE, 6), "\n")
    cat("Odds ratio (effect size)..:", round(OReffect, 3), "\n")
    cat("Lower 95% CI..............:", round(CI_low, 3), "\n")
    cat("Upper 95% CI..............:", round(CI_up, 3), "\n")
    cat("Z-value...................:", round(zvalue, 6), "\n")
    cat("P-value...................:", signif(pvalue, 8), "\n")
    cat("Hosmer and Lemeshow r^2...:", round(R.l, 6), "\n")
    cat("Cox and Snell r^2.........:", round(R.cs, 6), "\n")
    cat("Nagelkerke's pseudo r^2...:", round(R.n, 6), "\n")
    cat("Sample size of AE DB......:", AE_N, "\n")
    cat("Sample size of model......:", sample_size, "\n")
    cat("Missing data %............:", round(Perc_Miss, 6), "\n")
  }
  return(output)
  print(output)
}

TRAITS.CON = c("MacrophagesPercBC", "SMCPercBC", "VesselDensityBC")

TRAITS.BIN = c("Calc.bin", "Collagen.bin", "Fat.bin_40", "IPH.bin")

cat("* Analysis of continuous/quantitative plaque traits as a function of current smoking...")
GLM.results <- data.frame(matrix(NA, ncol = 15, nrow = 0))
for (trait in 1:length(TRAITS.CON)) {
  TRAIT = TRAITS.CON[trait]
  print(TRAIT)
  currentDF = AEData_2016.update.DF.CEA.TFvar
  ### univariate
  # fit <- lm(currentDF[,TRAITS.CON[trait]] ~ SmokerCurrent,
  #           data  =  currentDF)
  ### EWAS-like
  # fit <- lm(currentDF[,TRAITS.CON[trait]] ~ SmokerCurrent + Age + Gender,
  #            data  =  currentDF)
  ### Multivariate
  fit <- lm(currentDF[,TRAITS.CON[trait]] ~ SmokerCurrent + Age + Gender + Hospital + BMI + DM.composite + Hypertension1 + CAD_history + PAOD + TC_final + LDL_final + GFR_MDRD + Med.anticoagulants + Med.diuretic + RAAS_med + Med.bblocker + ORyear,
            data  =  currentDF)
  
  GLM.results.TEMP <- data.frame(matrix(NA, ncol = 15, nrow = 0))
  GLM.results.TEMP[1,] = GLM.CON(fit, "AEData_2016.update.DF.CEA.TFvar", "SmokerCurrent", TRAIT)
  GLM.results = rbind(GLM.results, GLM.results.TEMP)
}

cat("- Edit the column names...")
colnames(GLM.results) = c("Dataset", "Predictor", "Trait", 
                          "Beta", "s.e.m.", 
                          "OR", "low95CI", "up95CI", 
                          "T-value", "P-value", "r^2", "r^2_adj", "AE_N", "Model_N", "Perc_Miss")

cat("- Correct the variable types...")
GLM.results$Beta <- as.numeric(GLM.results$Beta)
GLM.results$s.e.m. <- as.numeric(GLM.results$s.e.m.)
GLM.results$OR <- as.numeric(GLM.results$OR)
GLM.results$low95CI <- as.numeric(GLM.results$low95CI)
GLM.results$up95CI <- as.numeric(GLM.results$up95CI)
GLM.results$`T-value` <- as.numeric(GLM.results$`T-value`)
GLM.results$`P-value` <- as.numeric(GLM.results$`P-value`)
GLM.results$`r^2` <- as.numeric(GLM.results$`r^2`)
GLM.results$`r^2_adj` <- as.numeric(GLM.results$`r^2_adj`)
GLM.results$`AE_N` <- as.numeric(GLM.results$`AE_N`)
GLM.results$`Model_N` <- as.numeric(GLM.results$`Model_N`)
GLM.results$`Perc_Miss` <- as.numeric(GLM.results$`Perc_Miss`)

# Save the data
cat("- Writing results to Excel-file...")
### Univariate
# write.xlsx(GLM.results,
#            file = paste0(OUT_loc, "/",Today,".aems450k.meta.Con.Univariate.PlaquePhenotypes.xlsx"),
#            row.names = FALSE, col.names = TRUE)
### EWAS-like
# write.xlsx(GLM.results,
#            file = paste0(OUT_loc, "/",Today,".aems450k.meta.Con.EWASlike.PlaquePhenotypes.xlsx"),
#            row.names = FALSE, col.names = TRUE)
### Multivariate
write.xlsx(GLM.results,
           file = paste0(OUT_loc, "/",Today,".aems450k.meta.Con.Multivariate.PlaquePhenotypes.xlsx"),
           row.names = FALSE, col.names = TRUE)
# Removing intermediates
cat("- Removing intermediate files...")
rm(TRAIT, trait, currentDF, GLM.results, GLM.results.TEMP, fit)

cat("* Analysis of binary/qualitative plaque traits as a function of current smoking...")
GLM.results <- data.frame(matrix(NA, ncol = 16, nrow = 0))
for (trait in 1:length(TRAITS.BIN)) {
  TRAIT = TRAITS.BIN[trait]
  print(TRAIT)
  currentDF = AEData_2016.update.DF.CEA.TFvar
  ### univariate
  # fit <- glm(as.factor(currentDF[,TRAITS.BIN[trait]]) ~ SmokerCurrent,
  #           data  =  currentDF, family = binomial())
  ### EWAS-like
  # fit <- glm(as.factor(currentDF[,TRAITS.BIN[trait]]) ~ SmokerCurrent + Age + Gender + Hospital,
  # data  =  currentDF, family = binomial())
  ### Multivariate
  fit <- glm(as.factor(currentDF[,TRAITS.BIN[trait]]) ~ SmokerCurrent + Age + Gender + Hospital + BMI + DM.composite + Hypertension1 + CAD_history + PAOD + TC_final + LDL_final + GFR_MDRD + Med.anticoagulants + Med.diuretic + RAAS_med + Med.bblocker + ORyear,
             data  =  currentDF, family = binomial())
  GLM.results.TEMP <- data.frame(matrix(NA, ncol = 16, nrow = 0))
  GLM.results.TEMP[1,] = GLM.BIN(fit, "AEData_2016.update.DF.CEA.TFvar", "SmokerCurrent", TRAIT)
  GLM.results = rbind(GLM.results, GLM.results.TEMP)
}

cat("- Edit the column names...")
colnames(GLM.results) = c("Dataset", "Predictor", "Trait", 
                          "Beta", "s.e.m.", 
                          "OR", "low95CI", "up95CI", 
                          "Z-value", "P-value", "r^2_l", "r^2_cs", "r^2_nagelkerke", "AE_N", "Model_N", "Perc_Miss")

cat("- Correct the variable types...")
GLM.results$Beta <- as.numeric(GLM.results$Beta)
GLM.results$s.e.m. <- as.numeric(GLM.results$s.e.m.)
GLM.results$OR <- as.numeric(GLM.results$OR)
GLM.results$low95CI <- as.numeric(GLM.results$low95CI)
GLM.results$up95CI <- as.numeric(GLM.results$up95CI)
GLM.results$`Z-value` <- as.numeric(GLM.results$`Z-value`)
GLM.results$`P-value` <- as.numeric(GLM.results$`P-value`)
GLM.results$`r^2_l` <- as.numeric(GLM.results$`r^2_l`)
GLM.results$`r^2_cs` <- as.numeric(GLM.results$`r^2_cs`)
GLM.results$`r^2_nagelkerke` <- as.numeric(GLM.results$`r^2_nagelkerke`)
GLM.results$`AE_N` <- as.numeric(GLM.results$`AE_N`)
GLM.results$`Model_N` <- as.numeric(GLM.results$`Model_N`)
GLM.results$`Perc_Miss` <- as.numeric(GLM.results$`Perc_Miss`)

# Save the data
cat("- Writing results to Excel-file...")

### Univariate
# write.xlsx(GLM.results,
#            file = paste0(OUT_loc, "/",Today,".aems450k.meta.Bin.Univariate.PlaquePhenotypes.xlsx"),
#            row.names = FALSE, col.names = TRUE)
### EWAS-like
# write.xlsx(GLM.results,
#            file = paste0(OUT_loc, "/",Today,".aems450k.meta.Bin.EWASlike.PlaquePhenotypes.xlsx"),
#            row.names = FALSE, col.names = TRUE)
### Multivariate
write.xlsx(GLM.results,
           file = paste0(OUT_loc, "/",Today,".aems450k.meta.Bin.Multivariate.PlaquePhenotypes.xlsx"),
           row.names = FALSE, col.names = TRUE)

# Removing intermediates
cat("- Removing intermediate files...")
rm(TRAIT, trait, currentDF, GLM.results, GLM.results.TEMP, fit)

cat("\n===========================================================================================")
cat("SAVE THE DATA")

save.image(paste0(ANALYSIS_loc,"/",Today,".aems450k.meta.analysis.ewas.",EWAS_trait,"_vs_PlaquePheno.RData"))

