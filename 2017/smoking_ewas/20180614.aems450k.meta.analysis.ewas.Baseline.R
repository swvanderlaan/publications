cat("===========================================================================================
                META-ANALYSIS ATHERO-EXPRESS METHYLATION STUDIES 450K 1 & 2
                                   --- BASELINE TABLES ---
    
    Version:      v2.9
    
    Last update:  2018-06-27
    Written by:   Sander W. van der Laan (s.w.vanderlaan-2@umcutrecht.nl);
    Marten A. Siemelink
    
    Description:  Script to make baseline tables.
    
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
# install.packages.auto("RMySQL") # install via RStudio's CRAN
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
install.packages.auto("VariantAnnotation")
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
ifelse(!dir.exists(file.path(ANALYSIS_loc, "/BASELINE")), 
       dir.create(file.path(ANALYSIS_loc, "/BASELINE")), 
       FALSE)
BASELINE_loc = paste0(ANALYSIS_loc, "/BASELINE")

cat("===========================================================================================")
cat("\nLOAD ATHERO-EXPRESS METHYLATION STUDY DATASETS")
setwd(INP_loc)
list.files()

cat("\n  - loading B/Mvalues of plaque samples...")
load(paste0(INP_AEMS450K1_loc,"/20171229.aems450k1.MvaluesQCIMP.blood.RData"))
load(paste0(INP_AEMS450K1_loc,"/20171229.aems450k1.MvaluesQCIMP.plaque.RData"))
load(paste0(INP_AEMS450K2_loc,"/20171229.aems450k2.MvaluesQCIMP.plaque.RData"))

cat("\n* Setup the analysis.")
require(FDb.InfiniumMethylation.hg19)
feats <- features(FDb.InfiniumMethylation.hg19)
chr.list <- levels(seqnames(feats))
regions <- feats[seqnames(feats) %in% chr.list]

cat("\n  - Setup the model, first covariate is the variable/phenotype of interest.")
cat("  > AEMS450K1...")
aems450k1.MvaluesQCplaqueClean <- aems450k1.MvaluesQCplaque
dim(aems450k1.MvaluesQCplaque)
metadata(aems450k1.MvaluesQCplaqueClean)$formula <- ~SmokerCurrent + Sample_Sex + Age + Hospital
aems450k1.MvaluesQCbloodClean <- aems450k1.MvaluesQCblood
dim(aems450k1.MvaluesQCblood)
metadata(aems450k1.MvaluesQCbloodClean)$formula <- ~SmokerCurrent + Sample_Sex + Age + Hospital
cat("  > AEMS450K2...")
aems450k2.MvaluesQCplaqueClean <- aems450k2.MvaluesQCplaque
dim(aems450k2.MvaluesQCplaque)
metadata(aems450k2.MvaluesQCplaqueClean)$formula <- ~SmokerCurrent + Sample_Sex + Age + Hospital

cat("===========================================================================================")
cat("CREATE BASELINE TABLE")

# Baseline table variables
basetable_vars = c("Age", "Gender", 
                   "TC_final", "LDL_final", "HDL_final", "TG_final", 
                   "systolic", "diastoli", "GFR_MDRD", "BMI", 
                   "KDOQI", "BMI_WHO",
                   "eCigarettes", "ePackYearsSmoking",
                   "DM.composite", "Hypertension.composite", 
                   "Hypertension.drugs", "Med.anticoagulants", "Med.all.antiplatelet", "Med.Statin.LLD", 
                   "Stroke_Dx", "Symptoms.4g", "restenos")
basetable_bin = c("Gender", 
                  "KDOQI", "BMI_WHO",
                  "DM.composite", "Hypertension.composite", 
                  "Hypertension.drugs", "Med.anticoagulants", "Med.all.antiplatelet", "Med.Statin.LLD", 
                  "Stroke_Dx", "Symptoms.4g", "restenos")
basetable_con = basetable_vars[!basetable_vars %in% basetable_bin]

cat("Check for Normality of continues variables.")
pdf(paste0(BASELINE_loc,"/",Today,".aems450k1.",EWAS_trait,".blood.BaselineNormality.pdf"), 
    paper = "a4", onefile = TRUE)
  par(mfrow = c(2,2))
  for (x in basetable_con) {
    qqnorm(colData(aems450k1.MvaluesQCbloodClean)[,x], 
           main = x, 
           cex.main = 1.0, cex.lab = 0.90, cex.axis = 0.80, cex = 1.25,
           col = "#1290D9", 
           pch = 20, bty = "n")
    qqline(colData(aems450k1.MvaluesQCbloodClean)[,x], col = "#E55738", lty = 2, lwd = 2)
  }
dev.off()

pdf(paste0(BASELINE_loc,"/",Today,".aems450k1.",EWAS_trait,".plaque.BaselineNormality.pdf"), 
    paper = "a4", onefile = TRUE)
  par(mfrow = c(2,2))
  for (x in basetable_con) {
    qqnorm(colData(aems450k1.MvaluesQCplaqueClean)[,x], 
           main = x, 
           cex.main = 1.0, cex.lab = 0.90, cex.axis = 0.80, cex = 1.25,
           col = "#1290D9", 
           pch = 20, bty = "n")
    qqline(colData(aems450k1.MvaluesQCplaqueClean)[,x], col = "#E55738", lty = 2, lwd = 2)
  }
dev.off()
pdf(paste0(BASELINE_loc,"/",Today,".aems450k2.",EWAS_trait,".plaque.BaselineNormality.pdf"), 
    paper = "a4", onefile = TRUE)
  par(mfrow = c(2,2))
  for (x in basetable_con) {
    qqnorm(colData(aems450k2.MvaluesQCplaqueClean)[,x], 
           main = x, 
           cex.main = 1.0, cex.lab = 0.90, cex.axis = 0.80, cex = 1.25,
           col = "#1290D9", 
           pch = 20, bty = "n")
    qqline(colData(aems450k2.MvaluesQCplaqueClean)[,x], col = "#E55738", lty = 2, lwd = 2)
  }
dev.off()

cat("* Plotting histograms for all data...")
pdf(paste0(BASELINE_loc,"/",Today,".aems450k1.",EWAS_trait,".blood.Transformation.Histogram.pdf"), 
    paper = "a4r", onefile = TRUE)
par(mfrow = c(2,2), mar = c(3,3,3,1))
for (x in basetable_con) {
  hist(colData(aems450k1.MvaluesQCbloodClean)[,x], main = x, 
       xlab = "",
       col = "#1290D9", 
       cex.main = 1.0, cex.lab = 0.90, cex.axis = 0.90, cex = 1.25)
  mtext(table(is.na(colData(aems450k1.MvaluesQCbloodClean)[,x]))["FALSE"], cex = 0.90)
  mtext(round(table(colData(aems450k1.MvaluesQCbloodClean)[,x] == 0)["TRUE"]/table(is.na(colData(aems450k1.MvaluesQCbloodClean)[,x]))["FALSE"]*100,2), 
        cex = 0.90, side = 4, line = -1)
}
dev.off()

cat("* Plotting histograms for all data...")
pdf(paste0(BASELINE_loc,"/",Today,".aems450k1.",EWAS_trait,".plaque.Transformation.Histogram.pdf"), 
    paper = "a4r", onefile = TRUE)
par(mfrow = c(2,2), mar = c(3,3,3,1))
for (x in basetable_con) {
  hist(colData(aems450k1.MvaluesQCplaqueClean)[,x], main = x, 
       xlab = "",
       col = "#1290D9", 
       cex.main = 1.0, cex.lab = 0.90, cex.axis = 0.90, cex = 1.25)
  mtext(table(is.na(colData(aems450k1.MvaluesQCplaqueClean)[,x]))["FALSE"], cex = 0.90)
  mtext(round(table(colData(aems450k1.MvaluesQCplaqueClean)[,x] == 0)["TRUE"]/table(is.na(colData(aems450k1.MvaluesQCplaqueClean)[,x]))["FALSE"]*100,2), 
        cex = 0.90, side = 4, line = -1)
}
dev.off()

cat("* Plotting histograms for all data...")
pdf(paste0(BASELINE_loc,"/",Today,".aems450k2.",EWAS_trait,".plaque.Transformation.Histogram.pdf"), 
    paper = "a4r", onefile = TRUE)
par(mfrow = c(2,2), mar = c(3,3,3,1))
for (x in basetable_con) {
  hist(colData(aems450k2.MvaluesQCplaqueClean)[,x], main = x, 
       xlab = "",
       col = "#1290D9", 
       cex.main = 1.0, cex.lab = 0.90, cex.axis = 0.90, cex = 1.25)
  mtext(table(is.na(colData(aems450k2.MvaluesQCplaqueClean)[,x]))["FALSE"], cex = 0.80)
  mtext(round(table(colData(aems450k2.MvaluesQCplaqueClean)[,x] == 0)["TRUE"]/table(is.na(colData(aems450k2.MvaluesQCplaqueClean)[,x]))["FALSE"]*100,2), 
        cex = 0.80, side = 4, line = -1)
}
dev.off()

# Create baseline tables
AMES450K1.tableOne_plaque = print(CreateTableOne(vars = basetable_vars, strata = EWAS_trait, 
                                                 data = as.data.frame(colData(aems450k1.MvaluesQCplaqueClean)), 
                                                 factorVars = basetable_bin), 
                                  nonnormal = c(), quote = FALSE, showAllLevels = FALSE, cramVars = "Gender",
                                  format = "p", contDigits = 1)[,1:3]
AMES450K2.tableOne_plaque = print(CreateTableOne(vars = basetable_vars, strata = EWAS_trait, 
                                                 data = as.data.frame(colData(aems450k2.MvaluesQCplaqueClean)), 
                                                 factorVars = basetable_bin), 
                                  nonnormal = c(), quote = FALSE, showAllLevels = FALSE, cramVars = "Gender",
                                  format = "p", contDigits = 1)[,1:3]
AMES450K1.tableOne_blood = print(CreateTableOne(vars = basetable_vars, strata = EWAS_trait, 
                                                data = as.data.frame(colData(aems450k1.MvaluesQCbloodClean)), 
                                                factorVars = basetable_bin), 
                                 nonnormal = c(), quote = FALSE, showAllLevels = FALSE, cramVars = "Gender",
                                 format = "p", contDigits = 1)[,1:3]

# Write basetable
write.xlsx(file = paste0(BASELINE_loc, "/",Today,".aems450k1.meta.",EWAS_trait,".plaque.Basetable.xlsx"), 
           format(AMES450K1.tableOne_plaque, digits = 5, scientific = FALSE) , row.names = TRUE, col.names = TRUE)
write.xlsx(file = paste0(BASELINE_loc, "/",Today,".aems450k2.meta.",EWAS_trait,".plaque.Basetable.xlsx"), 
           format(AMES450K2.tableOne_plaque, digits = 5, scientific = FALSE) , row.names = TRUE, col.names = TRUE)
write.xlsx(file = paste0(BASELINE_loc, "/",Today,".aems450k1.meta.",EWAS_trait,".blood.Basetable.xlsx"), 
           format(AMES450K1.tableOne_blood, digits = 5, scientific = FALSE) , row.names = TRUE, col.names = TRUE)

# Check for NA, add manually to table
CLIN = colData(aems450k1.MvaluesQCplaqueClean)[, as.integer(labels(colnames(colData(aems450k1.MvaluesQCplaqueClean))))]
CLIN = CLIN[!is.na(CLIN$SmokerCurrent),]
for (x in basetable_vars) {
  print(paste(x, "smoker", table(is.na(CLIN[CLIN$SmokerCurrent == "yes", x]))["FALSE"]))
  print(paste(x, "non-smoker", table(is.na(CLIN[CLIN$SmokerCurrent == "no", x]))["FALSE"]))
}
table(CLIN[CLIN$SmokerCurrent == "yes", "Symptoms.4g"], useNA = "always")
table(CLIN[CLIN$SmokerCurrent == "no", "Symptoms.4g"], useNA = "always")

CLIN = colData(aems450k2.MvaluesQCplaqueClean)[, as.integer(labels(colnames(colData(aems450k2.MvaluesQCplaqueClean))))]
CLIN = CLIN[!is.na(CLIN$SmokerCurrent),]
for (x in basetable_vars) {
  print(paste(x, "smoker", table(is.na(CLIN[CLIN$SmokerCurrent == "yes", x]))["FALSE"]))
  print(paste(x, "non-smoker", table(is.na(CLIN[CLIN$SmokerCurrent == "no", x]))["FALSE"]))
}
table(CLIN[CLIN$SmokerCurrent == "yes", "Symptoms.4g"], useNA = "always")
table(CLIN[CLIN$SmokerCurrent == "no", "Symptoms.4g"], useNA = "always")

CLIN = colData(aems450k1.MvaluesQCbloodClean)[, as.integer(labels(colnames(colData(aems450k1.MvaluesQCbloodClean))))]
CLIN = CLIN[!is.na(CLIN$SmokerCurrent),]
for (x in basetable_vars) {
  print(paste(x, "smoker", table(is.na(CLIN[CLIN$SmokerCurrent == "yes", x]))["FALSE"]))
  print(paste(x, "non-smoker", table(is.na(CLIN[CLIN$SmokerCurrent == "no", x]))["FALSE"]))
}
table(CLIN[CLIN$SmokerCurrent == "yes", "Symptoms.4g"], useNA = "always")
table(CLIN[CLIN$SmokerCurrent == "no", "Symptoms.4g"], useNA = "always")

cat("===========================================================================================")
cat("BOX-COX TRANSFORMATION OF PLAQUE PHENOTYPES")
cat("Add in Box-Cox-transformed phenotypes: Transform Phenotype Measurements with Zero-values due 
    to Detection-limits.")

install.packages.auto("Hmisc") #some helpfull stuff     
install.packages.auto("foreign") #for importing SPSS   
install.packages.auto("stats") #some statistics
install.packages.auto("geoR") #for Box-cox -- install manually
# library("geoR")
install.packages.auto("plyr") #for Box-Cox
install.packages.auto("MASS") #for Box-Cox
install.packages.auto("AID") #for Box-Cox

cat("* Select and Transform data.")

### Options for data-transformation containing zero's 
### Drop the zero's to missing data (NA)
### Log(x+1)
### Log(x+c) with c = estimated or very small
### Box-Cox transformation (Box and Cox 1964): boxcoxfit in package geoR
### Inverse Hyperbolic Sine transformation (Burbidge, Magee and Robb - 1988)
### Mixed model. Binary variable = zero/not-zero + continues variable for non-zero measurements (only for independent vars in regression)
### Poisson Regression seems INVALID, it concerns TRUE zero's, not values below detection limit
### Quantile-regression model. Uses non-parametric model centered on slope at 50-percentile (or other set percentile)

# Function to compute Box-Cox transformation, from the internet
BC.trans = function(x){
  t <- try(lambda.pair <- boxcoxfit(x, lambda2 = TRUE)$lambda)
  
  # Estimating both lambdas sometimes fails; if so, estimate lambda1 only
  if (inherits(t, "try-error")) {
    lambda1 <- boxcoxfit(x)$lambda
    lambda2 <- 0
  } else {
    lambda1 <- lambda.pair[1]
    lambda2 <- lambda.pair[2]
  }
  
  boxcox.f <- function(x, lambda1, lambda2) {
    if (lambda1 != 0) {
      return(((x + lambda2) ^ lambda1 - 1) / lambda1)
    } else {
      return(log(x + lambda2))
    }
  }
  bc.data <- boxcox.f(x, lambda1, lambda2)
  return(bc.data)
}

# Function for setting 3SD outliers to NA
remove.outlier = function(vector){
  vector[abs(mean(vector, na.rm = TRUE) - vector) > (3 * sd(vector, na.rm = TRUE))] = NA
  return(vector)
}

cat("* Load All Athero-Express Biobank Study phenotype data -- v2016-3; 20160519...")
# row number = study number
AEData_2016 <- read.xlsx(paste0(INP_AE_loc,"/AEMS450KCombo/SampleInformation/20170330_AEDB_v3_20160519.xlsx"),
                         sheet = 1, skipEmptyRows = TRUE)

# Make dataframe of raw data
# str(AEData_2016, list.len = ncol(AEData_2016))
AEData_2016.data.raw <- as.data.frame(AEData_2016)[,c(1:2, 105:106,876)]

# NB. data in column 35 fails due to infinite numbers in BC calculation. 
# Division by small factor helps (1.01). INVESTIGATE !!!
# This turns out to be the case for all variables when you want to 
# remove outliers AFTER transformation. So we decided to divide all 
# data by 1000
AEData_2016.data.raw[,3:ncol(AEData_2016.data.raw)] = AEData_2016.data.raw[,3:ncol(AEData_2016.data.raw)]/1000

# Perform transformations BEFORE removing outliers
AEData_2016.data.raw.log = apply(AEData_2016.data.raw[,3:ncol(AEData_2016.data.raw)], MARGIN = 2, FUN = function(vec){log(vec + 1)} )
AEData_2016.data.raw.log2 = apply(AEData_2016.data.raw[,3:ncol(AEData_2016.data.raw)], MARGIN = 2, FUN = function(vec){log2(vec + 1)} )
AEData_2016.data.log = cbind(AEData_2016[,1:2], apply(AEData_2016.data.raw.log, MARGIN = 2, FUN = remove.outlier ))
AEData_2016.data.log2 = cbind(AEData_2016[,1:2], apply(AEData_2016.data.raw.log2, MARGIN = 2, FUN = remove.outlier ))
# Perform BOX-COX transformations BEFORE removing outliers 
AEData_2016.data.raw.BC = apply(AEData_2016.data.raw[,3:ncol(AEData_2016.data.raw)], MARGIN = 2, FUN = BC.trans)
AEData_2016.data.BC = cbind(AEData_2016[,1:2],apply(AEData_2016.data.raw.BC, MARGIN = 2, FUN = remove.outlier))

# Get raw data without outliers
AEData_2016.data.rm = cbind(AEData_2016[,1:2], apply(AEData_2016.data.raw[,3:ncol(AEData_2016.data.raw)], MARGIN = 2, FUN = remove.outlier))
rm(AEData_2016.data.raw.BC, AEData_2016.data.raw.log, AEData_2016.data.raw.log2)

# Below is included a Hyperbolic-Sine Transform, this is only better in a very few cases, Excluded for now
# data.ihs = cbind(data[,1:2],apply(data.rm, MARGIN=2, FUN=function(vec){log(vec+ sqrt(vec^vec+1))} ))
# data.ihs[data.ihs == Inf] = NA

cat("* Plotting qq-plots for all data...")
pdf(paste0(BASELINE_loc,"/",Today,".AE.CEA.",EWAS_trait,".Transformation.QQ.pdf"), 
    paper = "a4r", onefile = TRUE)
  par(mfrow = c(3,5), mar = c(3,3,3,1))
  for (x in 3:5) {
    qqnorm(AEData_2016.data.raw[,x], main = paste0("Raw ",colnames(AEData_2016.data.raw)[x]), 
           cex.main = 0.60)
    qqline(AEData_2016.data.raw[,x], col = "#595A5C")
    qqnorm(AEData_2016.data.rm[,x], main = paste0("Raw ",colnames(AEData_2016.data.rm)[x],"\noutliers removed"), 
           cex.main = 0.60)
    qqline(AEData_2016.data.rm[,x], col = "#F59D10")
    qqnorm(AEData_2016.data.log[,x], main = paste0("LN-transformed ",colnames(AEData_2016.data.log)[x]), 
           cex.main = 0.60)
    qqline(AEData_2016.data.log[,x], col = "#2F8BC9")
    qqnorm(AEData_2016.data.log2[,x], main = paste0("Log2-transformed ",colnames(AEData_2016.data.log2)[x]), 
           cex.main = 0.60)
    qqline(AEData_2016.data.log2[,x], col = "#E55738")
    qqnorm(AEData_2016.data.BC[,x], main = paste0("BC-transformed ",colnames(AEData_2016.data.BC)[x]), 
           cex.main = 0.60)
    qqline(AEData_2016.data.BC[,x], col = "#9FC228")
  }
dev.off()

cat("* Plotting histograms for all data...")
pdf(paste0(BASELINE_loc,"/",Today,".AE.CEA.",EWAS_trait,".Transformation.Histogram.pdf"), 
    paper = "a4r", onefile = TRUE)
  par(mfrow = c(3,5), mar = c(3,3,3,1))
  for (x in 3:5) {
    hist(AEData_2016.data.raw[,x], main = paste0("Raw ", colnames(AEData_2016.data.raw)[x]), 
         col = "#595A5C", cex.main = 0.60)
    mtext(table(is.na(AEData_2016.data.raw[,x]))["FALSE"], cex = 0.60)
    mtext(round(table(AEData_2016.data.raw[,x] == 0)["TRUE"]/table(is.na(AEData_2016.data.raw[,x]))["FALSE"]*100,2), cex = 0.60, side = 4, line = -1)
    
    hist(AEData_2016.data.rm[,x], main = paste0("Raw ", colnames(AEData_2016.data.rm)[x], "\noutliers removed"), 
         col = "#F59D10", cex.main = 0.60)
    mtext(table(is.na(AEData_2016.data.rm[,x]))["FALSE"], cex = 0.60)
    mtext(round(table(AEData_2016.data.rm[,x] == 0)["TRUE"]/table(is.na(AEData_2016.data.rm[,x]))["FALSE"]*100,2), cex = 0.60, side = 4, line = -1)
    
    hist(AEData_2016.data.log[,x], main = paste0("LN-transformed ", colnames(AEData_2016.data.log)[x]), 
         col = "#2F8BC9", cex.main = 0.60)
    mtext(table(is.na(AEData_2016.data.log[,x]))["FALSE"], cex = 0.60)
    mtext(round(table(AEData_2016.data.log[,x] == 0)["TRUE"]/table(is.na(AEData_2016.data.log[,x]))["FALSE"]*100,2), cex = 0.60, side = 4, line = -1)
    
    hist(AEData_2016.data.log2[,x], main = paste0("Log2-transformed ", colnames(AEData_2016.data.log2)[x]), 
         col = "#E55738", cex.main = 0.60)
    mtext(table(is.na(AEData_2016.data.log2[,x]))["FALSE"], cex = 0.60)
    mtext(round(table(AEData_2016.data.log2[,x] == 0)["TRUE"]/table(is.na(AEData_2016.data.log2[,x]))["FALSE"]*100,2), cex = 0.60, side = 4, line = -1)
    
    hist(AEData_2016.data.BC[,x], main = paste0("BC-transformed ", colnames(AEData_2016.data.BC)[x]), 
         col = "#9FC228", cex.main = 0.60)
    mtext(table(is.na(AEData_2016.data.BC[,x]))["FALSE"], cex = 0.60)
    mtext(round(table(AEData_2016.data.BC[,x] == 0)["TRUE"]/table(is.na(AEData_2016.data.BC[,x]))["FALSE"]*100,2), cex = 0.60, side = 4, line = -1)
  }
dev.off()

cat("* Parsing new data...")
names(AEData_2016.data.log)[names(AEData_2016.data.log) == "macmean0"] <- "MacrophagesPercLN"
names(AEData_2016.data.log)[names(AEData_2016.data.log) == "smcmean0"] <- "SMCPercLN"
names(AEData_2016.data.log)[names(AEData_2016.data.log) == "vessel_density_averaged"] <- "VesselDensityLN"
names(AEData_2016.data.log2)[names(AEData_2016.data.log2) == "macmean0"] <- "MacrophagesPercLog2"
names(AEData_2016.data.log2)[names(AEData_2016.data.log2) == "smcmean0"] <- "SMCPercLog2"
names(AEData_2016.data.log2)[names(AEData_2016.data.log2) == "vessel_density_averaged"] <- "VesselDensityLog2"
names(AEData_2016.data.BC)[names(AEData_2016.data.BC) == "macmean0"] <- "MacrophagesPercBC"
names(AEData_2016.data.BC)[names(AEData_2016.data.BC) == "smcmean0"] <- "SMCPercBC"
names(AEData_2016.data.BC)[names(AEData_2016.data.BC) == "vessel_density_averaged"] <- "VesselDensityBC"

AEData_2016.data.new = cbind(AEData_2016[,1:2],AEData_2016.data.log[,3:5], AEData_2016.data.log2[,3:5], AEData_2016.data.BC[,3:5])

rm(AEData_2016.data.log, AEData_2016.data.log2, AEData_2016.data.BC, AEData_2016.data.raw, AEData_2016.data.rm)

cat("* Add these data to our SummarizedExperiment files...")
aems450k1.MvaluesQCbloodCleanBACKUP = aems450k1.MvaluesQCbloodClean
aems450k1.MvaluesQCplaqueCleanBACKUP = aems450k1.MvaluesQCplaqueClean
aems450k2.MvaluesQCplaqueCleanBACKUP = aems450k2.MvaluesQCplaqueClean
colData(aems450k1.MvaluesQCbloodClean) <- merge(colData(aems450k1.MvaluesQCbloodClean), AEData_2016.data.new, 
                                                 by.x = "STUDY_NUMBER", by.y = "STUDY_NUMBER", all.x = FALSE)
colData(aems450k1.MvaluesQCplaqueClean) <- merge(colData(aems450k1.MvaluesQCplaqueClean), AEData_2016.data.new, 
                                                 by.x = "STUDY_NUMBER", by.y = "STUDY_NUMBER", all.x = FALSE)

colData(aems450k2.MvaluesQCplaqueClean) <- merge(colData(aems450k2.MvaluesQCplaqueClean), AEData_2016.data.new, 
                                                 by.x = "STUDY_NUMBER", by.y = "STUDY_NUMBER", all.x = FALSE)

# We create these for ease of analyses downstream
aems450k1.DF = as.data.frame(colData(aems450k1.MvaluesQCplaqueClean))
aems450k2.DF = as.data.frame(colData(aems450k2.MvaluesQCplaqueClean))
aems450k.meta.DF <- rbind(aems450k1.DF, aems450k2.DF)

cat("* Add these data to our original AE database...")
AEData_2016.update.DF <- cbind(AEData_2016, AEData_2016.data.new)
cat("* Subset original AE database to only include CEA...")
AEData_2016.update.DF.CEA <- AEData_2016.update.DF[ which(AEData_2016.update.DF$Artery_summary == "carotid (left & right)"), ]
AEData_2016.update.DF.CEA.SampleSize <- nrow(AEData_2016.update.DF.CEA)

# By the way -- we need a baseline table of the whole cohort as well.

AEData_2016.update.DF.CEA.tableOne_blood = print(CreateTableOne(vars = basetable_vars, strata = EWAS_trait, 
                                                                data = as.data.frame(AEData_2016.update.DF.CEA), 
                                                                factorVars = basetable_bin), 
                                                 nonnormal = c(), quote = FALSE, showAllLevels = FALSE, 
                                                 format = "p", contDigits = 3)[,1:3]

# Write basetable
write.xlsx(file = paste0(BASELINE_loc, "/",Today,".AE.CEA.Basetable.xlsx"), 
           format(AEData_2016.update.DF.CEA.tableOne_blood, digits = 5, scientific = FALSE) , row.names = TRUE, col.names = TRUE)

# Get the missing numbers
CLIN = AEData_2016.update.DF.CEA[, as.integer(labels(colnames(AEData_2016.update.DF.CEA)))]
CLIN = CLIN[!is.na(CLIN$SmokerCurrent),]
for (x in basetable_vars) {
  print(paste(x, "smoker", table(is.na(CLIN[CLIN$SmokerCurrent == "yes", x]))["FALSE"]))
  print(paste(x, "non-smoker", table(is.na(CLIN[CLIN$SmokerCurrent == "no", x]))["FALSE"]))
}
table(CLIN[CLIN$SmokerCurrent == "yes", "Symptoms.4g"], useNA = "always")
table(CLIN[CLIN$SmokerCurrent == "no", "Symptoms.4g"], useNA = "always")

cat("\n*** saving UPDATE data in between steps ***")
aems450k1.MvaluesQCbloodCleanTFvar = aems450k1.MvaluesQCbloodClean
aems450k1.MvaluesQCplaqueCleanTFvar = aems450k1.MvaluesQCplaqueClean
aems450k2.MvaluesQCplaqueCleanTFvar = aems450k2.MvaluesQCplaqueClean
AEData_2016.update.DF.CEA.TFvar = AEData_2016.update.DF.CEA
save(AEData_2016.update.DF.CEA.TFvar, file = paste0(OUT_loc,"/",Today,".AEData_2016.update.DF.CEA.TFvar.RData"))
save(aems450k1.MvaluesQCbloodCleanTFvar, file = paste0(OUT_loc,"/",Today,".aems450k1.MvaluesQCbloodClean.TFvar.RData"))
save(aems450k1.MvaluesQCplaqueCleanTFvar, file = paste0(OUT_loc,"/",Today,".aems450k1.MvaluesQCplaqueClean.TFvar.RData"))
save(aems450k2.MvaluesQCplaqueCleanTFvar, file = paste0(OUT_loc,"/",Today,".aems450k2.MvaluesQCplaqueClean.TFvar.RData"))

cat("\nRemove intermediate files.")
rm(aems450k1.DF, aems450k2.DF, aems450k.meta.DF,
   aems450k1.MvaluesQCbloodCleanBACKUP, aems450k1.MvaluesQCplaqueCleanBACKUP, aems450k2.MvaluesQCplaqueCleanBACKUP,
   aems450k1.MvaluesQCblood, aems450k1.MvaluesQCplaque, aems450k2.MvaluesQCplaque,
   aems450k1.BvaluesQCblood, aems450k1.BvaluesQCplaque, aems450k2.BvaluesQCplaque,
   aems450k1.MvaluesQCbloodClean, aems450k1.MvaluesQCplaqueClean, aems450k2.MvaluesQCplaqueClean,
   CLIN, feats, regions,
   AEData_2016.data.new, AEData_2016.update.DF.CEA)

cat("\n===========================================================================================")
cat("SAVE THE DATA")

save.image(paste0(ANALYSIS_loc,"/",Today,".aems450k.meta.analysis.ewas.",EWAS_trait,".BASELINE.RData"))

