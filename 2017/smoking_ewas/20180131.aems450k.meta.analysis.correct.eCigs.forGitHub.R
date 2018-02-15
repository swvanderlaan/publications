cat("===========================================================================================
                META-ANALYSIS ATHERO-EXPRESS METHYLATION STUDIES 450K 1 & 2
                        --- correction of eCigarettes variable ---
    
    Version:      v1.0
    
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
install.packages.auto("RMySQL") # "RMySQL" will be phased out: https://github.com/r-dbi/RMySQL,
# will be replaced by RMariaDB (https://github.com/r-dbi/RMariaDB)
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

cat(paste0("\n\n* Load ",PROJECTDATASET," data..."))
cat("\n  - loading B/Mvalues of blood samples...") # not available for AEMS450K2
load(paste0(INP_AEMS450K1_loc,"/20170721.aems450k1.MvaluesQCIMP.blood.RData"))
load(paste0(INP_AEMS450K1_loc,"/20170721.aems450k1.BvaluesQCIMP.blood.RData"))
aems450k1.MvaluesQCblood = MvaluesQCblood
aems450k1.BvaluesQCblood = BvaluesQCblood

cat("\n  - loading B/Mvalues of plaque samples...")
load(paste0(INP_AEMS450K1_loc,"/20170721.aems450k1.MvaluesQCIMP.plaque.RData"))
load(paste0(INP_AEMS450K1_loc,"/20170721.aems450k1.BvaluesQCIMP.plaque.RData"))
cat("\n  - renaming these...")
aems450k1.MvaluesQCplaque = MvaluesQCplaque
aems450k1.BvaluesQCplaque = BvaluesQCplaque

load(paste0(INP_AEMS450K2_loc,"/20170822.aems450k2.MvaluesQCIMP.plaque.RData"))
load(paste0(INP_AEMS450K2_loc,"/20170822.aems450k2.BvaluesQCIMP.plaque.RData"))
cat("\n  - renaming these...")
aems450k2.MvaluesQCplaque = MvaluesQCplaque
aems450k2.BvaluesQCplaque = BvaluesQCplaque
rm(MvaluesQCplaque, BvaluesQCplaque, MvaluesQCblood, BvaluesQCblood) 


cat("===========================================================================================")
cat("\n[ CORRECTING ESTIMATED NUMBER OF PACK YEARS -- as a response of reviewer question ]")

GENOMIC_loc = "location to SPSS data on our UMC server"
AEDATA_loc = paste0(GENOMIC_loc,"/AE-AAA_GS_DBs")
AEdata = as.data.table(read_spss(paste0(AEDATA_loc,"/2017-4_AtheroExpressDatabase_ScientificAE_20171207.sav")))

AEdataForUpdate <- subset(AEdata, select = c(STUDY_NUMBER, ePackYearsSmoking))
setDF(AEdataForUpdate)
str(AEdataForUpdate)
AEdataForUpdate$STUDY_NUMBER <- as.numeric(AEdataForUpdate$STUDY_NUMBER)
AEdataForUpdate$ePackYearsSmoking <- as.numeric(AEdataForUpdate$ePackYearsSmoking)
str(AEdataForUpdate)

keep.aems450k1.MvaluesQCblood <- aems450k1.MvaluesQCblood$STUDY_NUMBER
keep.aems450k1.BvaluesQCblood <- aems450k1.BvaluesQCblood$STUDY_NUMBER
keep.aems450k1.MvaluesQCplaque <- aems450k1.MvaluesQCplaque$STUDY_NUMBER
keep.aems450k1.BvaluesQCplaque <- aems450k1.BvaluesQCplaque$STUDY_NUMBER
keep.aems450k2.MvaluesQCplaque <- aems450k2.MvaluesQCplaque$STUDY_NUMBER
keep.aems450k2.BvaluesQCplaque <- aems450k2.BvaluesQCplaque$STUDY_NUMBER

AEdataForUpdate.aems450k1.MvaluesQCblood <- AEdataForUpdate[keep.aems450k1.MvaluesQCblood,]
AEdataForUpdate.aems450k1.BvaluesQCblood <- AEdataForUpdate[keep.aems450k1.BvaluesQCblood,]
AEdataForUpdate.aems450k1.MvaluesQCplaque <- AEdataForUpdate[keep.aems450k1.MvaluesQCplaque,]
AEdataForUpdate.aems450k1.BvaluesQCplaque <- AEdataForUpdate[keep.aems450k1.BvaluesQCplaque,]
AEdataForUpdate.aems450k2.MvaluesQCplaque <- AEdataForUpdate[keep.aems450k2.MvaluesQCplaque,]
AEdataForUpdate.aems450k2.BvaluesQCplaque <- AEdataForUpdate[keep.aems450k2.BvaluesQCplaque,]

aems450k1.MvaluesQCblood$ePackYearsSmoking <- NULL
aems450k1.BvaluesQCblood$ePackYearsSmoking <- NULL
aems450k1.MvaluesQCplaque$ePackYearsSmoking <- NULL
aems450k1.BvaluesQCplaque$ePackYearsSmoking <- NULL
aems450k2.MvaluesQCplaque$ePackYearsSmoking <- NULL
aems450k2.BvaluesQCplaque$ePackYearsSmoking <- NULL

colData(aems450k1.MvaluesQCblood) <- cbind(colData(aems450k1.MvaluesQCblood), AEdataForUpdate.aems450k1.MvaluesQCblood)
colData(aems450k1.BvaluesQCblood) <- cbind(colData(aems450k1.BvaluesQCblood), AEdataForUpdate.aems450k1.BvaluesQCblood)
colData(aems450k1.MvaluesQCplaque) <- cbind(colData(aems450k1.MvaluesQCplaque), AEdataForUpdate.aems450k1.MvaluesQCplaque)
colData(aems450k1.BvaluesQCplaque) <- cbind(colData(aems450k1.BvaluesQCplaque), AEdataForUpdate.aems450k1.BvaluesQCplaque)
colData(aems450k2.MvaluesQCplaque) <- cbind(colData(aems450k2.MvaluesQCplaque), AEdataForUpdate.aems450k2.MvaluesQCplaque)
colData(aems450k2.BvaluesQCplaque) <- cbind(colData(aems450k2.BvaluesQCplaque), AEdataForUpdate.aems450k2.BvaluesQCplaque)

### SANITY CHECK
# summary(aems450k1.MvaluesQCblood$ePackYearsSmoking)
# summary(aems450k1.BvaluesQCblood$ePackYearsSmoking)
# summary(aems450k1.MvaluesQCplaque$ePackYearsSmoking)
# summary(aems450k1.BvaluesQCplaque$ePackYearsSmoking)
# summary(aems450k2.MvaluesQCplaque$ePackYearsSmoking)
# summary(aems450k2.BvaluesQCplaque$ePackYearsSmoking)

cat("\n*** saving UPDATE data in between steps ***")
save(aems450k1.MvaluesQCblood, file = paste0(INP_AEMS450K1_loc,"/",Today,".aems450k1.MvaluesQCIMP.blood.RData"))
save(aems450k1.BvaluesQCblood, file = paste0(INP_AEMS450K1_loc,"/",Today,".aems450k1.BvaluesQCIMP.blood.RData"))
save(aems450k1.MvaluesQCplaque, file = paste0(INP_AEMS450K1_loc,"/",Today,".aems450k1.MvaluesQCIMP.plaque.RData"))
save(aems450k1.BvaluesQCplaque, file = paste0(INP_AEMS450K1_loc,"/",Today,".aems450k1.BvaluesQCIMP.plaque.RData"))
save(aems450k2.MvaluesQCplaque, file = paste0(INP_AEMS450K2_loc,"/",Today,".aems450k2.MvaluesQCIMP.plaque.RData"))
save(aems450k2.BvaluesQCplaque, file = paste0(INP_AEMS450K2_loc,"/",Today,".aems450k2.BvaluesQCIMP.plaque.RData"))

cat("\n*** removing the B-value data - as we don't need that ***")

rm(aems450k1.BvaluesQCblood, aems450k1.BvaluesQCplaque, aems450k2.BvaluesQCplaque)

