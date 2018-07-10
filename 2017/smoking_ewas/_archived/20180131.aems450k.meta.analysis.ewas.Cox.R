cat("===========================================================================================
                META-ANALYSIS ATHERO-EXPRESS METHYLATION STUDIES 450K 1 & 2
                                    -- Cox regressions --
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

cat(paste0("\n\n* Load ",PROJECTDATASET," data...\n"))
cat("\n  - loading Mvalues of plaque samples...\n")
load(paste0(OUT_loc,"/20180207.aems450k1.MvaluesQCplaqueClean.TFvar.RData"))
load(paste0(OUT_loc,"/20180207.aems450k2.MvaluesQCplaqueClean.TFvar.RData"))

aems450k1.MvaluesQCplaqueClean = aems450k1.MvaluesQCplaqueCleanTFvar
aems450k2.MvaluesQCplaqueClean = aems450k2.MvaluesQCplaqueCleanTFvar

cat("\n  - loading AE database...\n")
load(paste0(OUT_loc,"/20180207.AEData_2016.update.DF.CEA.TFvar.RData"))

AEData_2016.update.DF.CEA.SampleSize = nrow(AEData_2016.update.DF.CEA.TFvar)

AEData_2016.update.DF.CEA = AEData_2016.update.DF.CEA.TFvar

cat("\n  - loading EWAS analyses *REPLICATED* results...\n")
load(paste0(OUT_loc,"/20180207.SmokerCurrent.aems450k.meta.resultspfCGIQC.RData"))

aems450k.meta.resultspfCGIQC = aems450k.meta.resultspfCGIQC.SmokerCurrent

cat("===========================================================================================")
cat("CpGs vs OUTCOME")

# Reference: https://bioconductor.org/packages/devel/bioc/vignettes/MultiAssayExperiment/inst/doc/QuickStartMultiAssay.html
# If you want to suppress warnings and messages when installing/loading packages
# suppressPackageStartupMessages({})
install.packages.auto("survival")
install.packages.auto("survminer")
install.packages.auto("Hmisc")
# listing all significant CpGs after meta-analysis and FDR-adjustment
hilight.top.cpg <- head(aems450k.meta.resultspfCGIQC[order(aems450k.meta.resultspfCGIQC$PvalCorAdj.meta),], 4)[,1]

cat("* Creating function to summarize Cox regression and prepare container for results.")
# Function to get summary statistics from Cox regression model
COX.STAT <- function(coxfit, DATASET, OUTCOME, cpg){
  cat("Summarizing Cox regression results for '", cpg ,"' and its association to '",OUTCOME,"' in '",DATASET,"'.\n")
  if (nrow(summary(coxfit)$coefficients) == 1) {
    output = c(cpg, rep(NA,8))
    cat("Model not fitted; probably singular.\n")
  }else {
    cat("Collecting data.\n\n")
    cox.sum <- summary(coxfit)
    cox.effectsize = cox.sum$coefficients[1,1]
    cox.SE = cox.sum$coefficients[1,3]
    cox.HReffect = cox.sum$coefficients[1,2]
    cox.CI_low = exp(cox.effectsize - 1.96 * cox.SE)
    cox.CI_up = exp(cox.effectsize + 1.96 * cox.SE)
    cox.zvalue = cox.sum$coefficients[1,4]
    cox.pvalue = cox.sum$coefficients[1,5]
    cox.sample_size = cox.sum$n
    cox.nevents = cox.sum$nevent
    
    output = c(DATASET, OUTCOME, cpg, cox.effectsize, cox.SE, cox.HReffect, cox.CI_low, cox.CI_up, cox.zvalue, cox.pvalue, cox.sample_size, cox.nevents)
    cat("We have collected the following:\n")
    cat("Dataset used..............:", DATASET, "\n")
    cat("Outcome analyzed..........:", OUTCOME, "\n")
    cat("CpG.......................:", cpg, "\n")
    cat("Effect size...............:", round(cox.effectsize, 6), "\n")
    cat("Standard error............:", round(cox.SE, 6), "\n")
    cat("Odds ratio (effect size)..:", round(cox.HReffect, 3), "\n")
    cat("Lower 95% CI..............:", round(cox.CI_low, 3), "\n")
    cat("Upper 95% CI..............:", round(cox.CI_up, 3), "\n")
    cat("T-value...................:", round(cox.zvalue, 6), "\n")
    cat("P-value...................:", signif(cox.pvalue, 8), "\n")
    cat("Sample size in model......:", cox.sample_size, "\n")
    cat("Number of events..........:", cox.nevents, "\n")
  }
  return(output)
  print(output)
} 

# only MACE 
times = c("ep_major_t_3years")

endpoints = c("epmajor.3years")

summary(AEData_2016.update.DF.CEA$EP_major_time)
summary(AEData_2016.update.DF.CEA$epmajor.3years)

cat("* Check distribution of events over time.")
# Check distribution of events over time
pdf(file = paste0(PLOT_loc, "/",Today,".aems450k1.",EWAS_trait,".EventDistributionPerYear.pdf"), onefile = TRUE)
for (x in times) {
  hist(colData(aems450k1.MvaluesQCplaqueClean)[,x], main = x, breaks = 15,
       xlab = "year", col = uithof_color[16])
}
dev.off()
pdf(file = paste0(PLOT_loc, "/",Today,".aems450k2.",EWAS_trait,".EventDistributionPerYear.pdf"), onefile = TRUE)
for (x in times) {
  hist(colData(aems450k2.MvaluesQCplaqueClean)[,x], main = x, breaks = 15,
       xlab = "year", col = uithof_color[16])
}
dev.off()

# Set up a dataframe to receive results
COX.results <- data.frame(matrix(NA, ncol = 12, nrow = 0))

# Looping over each endpoint/time combination
for (i in 1:length(times)) {
  eptime = times[i]
  ep = endpoints[i]
  cat(paste0("* Analyzing the effect of significant CpGs on [",ep,"].\n"))
  cat(" - creating temporary SE for this work.\n")
  # TEMP.SE = aems450k1.MvaluesQCplaqueClean
  TEMP.SE = aems450k2.MvaluesQCplaqueClean
  cat(" - making a 'Surv' object and adding this to colData().\n")
  TEMP.SE$event <- as.integer(colData(TEMP.SE)[,ep] == "Excluded")
  TEMP.SE$y <- Surv(colData(TEMP.SE)[,eptime], TEMP.SE$event)
  cat(" - making strata of each of the top CpG and start survival analysis.\n")
  
  for (cpg in 1:length(hilight.top.cpg)) {
    cat(paste0("   > processing [",hilight.top.cpg[cpg],"]; ",cpg," out of ",length(hilight.top.cpg)," probes.\n"))
    # splitting into two groups
    TEMP.SE[[ hilight.top.cpg[cpg] ]] <- cut2(assay(TEMP.SE)[hilight.top.cpg[cpg],], g = 2)
    cat(paste0("   > cross tabulation of ",hilight.top.cpg[cpg],"-stratum.\n"))
    show(table(TEMP.SE[[ hilight.top.cpg[cpg] ]]))
    cat(paste0("\n   > fitting the model for ",hilight.top.cpg[cpg],"-stratum.\n"))
    fit <- survfit(as.formula(paste0("y ~ ", hilight.top.cpg[cpg])), data = colData(TEMP.SE))
    cat(paste0("\n   > make a Kaplan-Meier-shizzle...\n"))
    # make Kaplan-Meier curve and save it
    show(ggsurvplot(fit, data = colData(TEMP.SE),
                    palette = c("#DB003F", "#1290D9"),
                    # palete = c("F59D10", "#DB003F", "#49A01D", "#1290D9"), 
                    linetype = c(1,2),
                    # linetype = c(1,2,3,4),
                    # conf.int = FALSE, conf.int.fill = "#595A5C", conf.int.alpha = 0.1,
                    pval = FALSE, pval.method = FALSE, pval.size = 4,
                    risk.table = TRUE, risk.table.y.text = FALSE, tables.y.text.col = TRUE, fontsize = 4,
                    censor = FALSE, 
                    legend = "right",
                    legend.title = paste0("",hilight.top.cpg[cpg],""),
                    legend.labs = c("demthylated", "methylated"), 
                    title = paste0("Risk of ",ep,""), xlab = "Time [years]", font.main = c(16, "bold", "black")))
    # dev.copy2pdf(file = paste0(COX_loc,"/",
    #                   Today,".aems450k1.meta.plaque.survival.",ep,".2G.",
    #                   hilight.top.cpg[cpg],".pdf"), width = 12, height = 10, onefile = FALSE)
    dev.copy2pdf(file = paste0(COX_loc,"/",
                               Today,".aems450k2.meta.plaque.survival.",ep,".2G.",
                               hilight.top.cpg[cpg],".pdf"), width = 12, height = 10, onefile = FALSE)

    cat(paste0("\n   > perform the Cox-regression fashizzle and plot it...\n"))
    ### Do Cox-regression and plot it
    # correct for smoking
    cox = coxph(Surv(colData(TEMP.SE)[,eptime], event) ~ TEMP.SE[[ hilight.top.cpg[cpg] ]]+Age+Gender+SmokerCurrent, data = colData(TEMP.SE))
    coxplot = coxph(Surv(colData(TEMP.SE)[,eptime], event) ~ strata(TEMP.SE[[ hilight.top.cpg[cpg] ]])+Age+Gender+SmokerCurrent, data = colData(TEMP.SE))
    plot(survfit(coxplot), main = paste0("Cox proportional hazard of [",ep,"] per [",eptime,"]."),
         # ylim = c(0.2, 1), xlim = c(0,3), col = c("#595A5C", "#DB003F", "#1290D9"),
         ylim = c(0, 1), xlim = c(0,3), col = c("#DB003F", "#1290D9"),
         lty = c(1,2), lwd = 2,
         ylab = "Suvival probability", xlab = "FU time [years]",
         mark.time = FALSE, axes = FALSE, bty = "n")
    legend("topright", 
           c("demethylated", "methylated"),
           title = paste0("",hilight.top.cpg[cpg],""),
           # levels(TEMP.SE[[ hilight.top.cpg[cpg] ]]), 
           # label(c("demethylated", "methylated")), 
           col = c("#DB003F", "#1290D9"),
           lty = c(1,2), lwd = 2, 
           bty = "n")
    axis(side = 1, at = seq(0, 3, by = 1))
    axis(side = 2, at = seq(0, 1, by = 0.2))
    # dev.copy2pdf(file = paste0(COX_loc,"/",
    #                            Today,".aems450k1.meta.plaque.Cox.",ep,".2G.",
    #                            # Today,".aems450k1.meta.plaque.Cox.",ep,".4G.",
    #                            hilight.top.cpg[cpg],".withSMOKING.pdf"), height = 12, width = 10, onefile = TRUE)
    dev.copy2pdf(file = paste0(COX_loc,"/",
                               Today,".aems450k2.meta.plaque.Cox.",ep,".2G.",
                               # Today,".aems450k2.meta.plaque.Cox.",ep,".4G.",
                               hilight.top.cpg[cpg],".withSMOKING.pdf"), height = 12, width = 10, onefile = TRUE)
    show(summary(cox))
    
    cat(paste0("\n   > writing the Cox-regression fashizzle to Excel...\n"))
    
    COX.results.TEMP <- data.frame(matrix(NA, ncol = 12, nrow = 0))
    # COX.results.TEMP[1,] = COX.STAT(cox, "aems450k1.MvaluesQCplaqueClean", ep, hilight.top.cpg[cpg])
    COX.results.TEMP[1,] = COX.STAT(cox, "aems450k2.MvaluesQCplaqueClean", ep, hilight.top.cpg[cpg])
    COX.results = rbind(COX.results, COX.results.TEMP)
    
  }
}

cat("- Edit the column names...\n")
colnames(COX.results) = c("Dataset", "Outcome", "CpG", 
                          "Beta", "s.e.m.", 
                          "HR", "low95CI", "up95CI", 
                          "Z-value", "P-value", "SampleSize", "N_events")

cat("- Correct the variable types...\n")
COX.results$Beta <- as.numeric(COX.results$Beta)
COX.results$s.e.m. <- as.numeric(COX.results$s.e.m.)
COX.results$HR <- as.numeric(COX.results$HR)
COX.results$low95CI <- as.numeric(COX.results$low95CI)
COX.results$up95CI <- as.numeric(COX.results$up95CI)
COX.results$`Z-value` <- as.numeric(COX.results$`Z-value`)
COX.results$`P-value` <- as.numeric(COX.results$`P-value`)
COX.results$SampleSize <- as.numeric(COX.results$SampleSize)
COX.results$N_events <- as.numeric(COX.results$N_events)

# aems450k1.COX.results.withSMOKING <- COX.results
aems450k2.COX.results.withSMOKING <- COX.results

# Save the data
cat("- Writing results to Excel-file...\n")
head.style <- createStyle(textDecoration = "BOLD")
# write.xlsx(aems450k1.COX.results.withSMOKING,
#            file = paste0(OUT_loc, "/",Today,".aems450k1.meta.plaque.Cox.2G.withSMOKING.xlsx"),
#            creator = "Sander W. van der Laan",
#            sheetName = "Results", headerStyle = head.style,
#            row.names = FALSE, col.names = TRUE, overwrite = TRUE)
write.xlsx(aems450k2.COX.results.withSMOKING,
           file = paste0(OUT_loc, "/",Today,".aems450k2.meta.plaque.Cox.2G.withSMOKING.xlsx"),
           creator = "Sander W. van der Laan",
           sheetName = "Results", headerStyle = head.style,
           row.names = FALSE, col.names = TRUE, overwrite = TRUE)

# Removing intermediates
cat("- Removing intermediate files...\n")
rm(TEMP.SE, cpg, fit, cox, coxplot, COX.results, COX.results.TEMP, head.style)

cat("\n * Correct for bias and inflation using 'BACON'. 
    (ref: M. van Iterson, E.W. van Zwet, and B.T. Heijmans. 2017. Genome Biol. 18 (1): 19.).")
library("bacon")

cat(" * create objects for 'bacon'...")
cat("   > with SMOKING...")
aems450k.meta.outcome.t <- cbind(limmaZ.aems450k1 = aems450k1.COX.results.withSMOKING[,9],
                                 limmaZ.aems450k2 = aems450k2.COX.results.withSMOKING[,9])
aems450k.meta.outcome.es <- cbind(limmaB.aems450k1 = aems450k1.COX.results.withSMOKING[,4],
                                  limmaB.aems450k2 = aems450k2.COX.results.withSMOKING[,4])
aems450k.meta.outcome.se <- cbind(limmaSE.aems450k1 = aems450k1.COX.results.withSMOKING[,5],
                                  limmaSE.aems450k2 = aems450k2.COX.results.withSMOKING[,5])

cat(" * run 'bacon'...")
library(BiocParallel)
register(MulticoreParam(2))
aems450k.meta.outcome.bcp.t <- bacon(teststatistics = aems450k.meta.outcome.t, 
                                     effectsizes = NULL, standarderrors = NULL,
                                     verbose = TRUE) 
aems450k.meta.outcome.bcp.b <- bacon(teststatistics = NULL,
                                     effectsizes = aems450k.meta.outcome.es, standarderrors = aems450k.meta.outcome.se,
                                     verbose = TRUE)
cat("   - 'bacon' model result.")
aems450k.meta.outcome.bcp.t
aems450k.meta.outcome.bcp.b
cat("\n - estimates of the mixture fit")
estimates(aems450k.meta.outcome.bcp.t)
estimates(aems450k.meta.outcome.bcp.b)
cat("   - the inflation")
inflation(aems450k.meta.outcome.bcp.t)
inflation(aems450k.meta.outcome.bcp.b)
cat("   - the bias")
bias(aems450k.meta.outcome.bcp.t) 
bias(aems450k.meta.outcome.bcp.b) 

cat("\n * Perform meta-analysis.")
aems450k.meta.outcome.bcpm.b <- meta(aems450k.meta.outcome.bcp.b)
aems450k.meta.outcome.pvalsbcpm.b <- pval(aems450k.meta.outcome.bcpm.b)
cat("\n * Uncorrected results.")
aems450k.meta.outcome.bcpm.pvalsp <- pval(aems450k.meta.outcome.bcpm.b, corrected = FALSE)
aems450k.meta.outcome.bcpm.zp = qnorm(aems450k.meta.outcome.bcpm.pvalsp[,1]/2)
aems450k.meta.outcome.bcpm.lambdap = round(median(aems450k.meta.outcome.bcpm.zp^2)/0.4549364,3)
cat(paste0("\n  - lambda is: [",aems450k.meta.outcome.bcpm.lambdap,"]."))
cat("\n * Corrected results.")
aems450k.meta.outcome.bcpm.pvalsbcp <- pval(aems450k.meta.outcome.bcpm.b, corrected = TRUE)
aems450k.meta.outcome.bcpm.zbcpp = qnorm(aems450k.meta.outcome.bcpm.pvalsbcp[,1]/2)
aems450k.meta.outcome.bcpm.lambdabcp = round(median(aems450k.meta.outcome.bcpm.zbcpp^2)/0.4549364,3)
cat(paste0("\n  - lambda is: [",aems450k.meta.outcome.bcpm.lambdabcp,"]."))

cat("\n * Merging with results...")
# add in effect sizes
temp.results <- cbind(aems450k1.COX.results.withSMOKING, aems450k2.COX.results.withSMOKING, es(aems450k.meta.outcome.bcpm.b))
names(temp.results)[names(temp.results) == "limmaB.aems450k1"] <- "EffectSize.aems450k1"
names(temp.results)[names(temp.results) == "limmaB.aems450k2"] <- "EffectSize.aems450k2"
names(temp.results)[names(temp.results) == "meta"] <- "EffectSize.meta"
# add in standard errors
temp.results2 <- cbind(temp.results, se(aems450k.meta.outcome.bcpm.b))
names(temp.results2)[names(temp.results2) == "limmaSE.aems450k1"] <- "SE.aems450k1"
names(temp.results2)[names(temp.results2) == "limmaSE.aems450k2"] <- "SE.aems450k2"
names(temp.results2)[names(temp.results2) == "meta"] <- "SE.meta"
# add in uncorrected p-values
temp.results3 <- cbind(temp.results2, pval(aems450k.meta.outcome.bcpm.b, corrected = FALSE))
names(temp.results3)[names(temp.results3) == "limmaB.aems450k1"] <- "Pval.aems450k1"
names(temp.results3)[names(temp.results3) == "limmaB.aems450k2"] <- "Pval.aems450k2"
temp.results4 <- temp.results3[, -which(names(temp.results3) %in% c("meta"))]
# add in corrected p-values
temp.results5 <- cbind(temp.results4, pval(aems450k.meta.outcome.bcpm.b, corrected = TRUE))
names(temp.results5)[names(temp.results5) == "limmaB.aems450k1"] <- "PvalCor.aems450k1"
names(temp.results5)[names(temp.results5) == "limmaB.aems450k2"] <- "PvalCor.aems450k2"
names(temp.results5)[names(temp.results5) == "meta"] <- "PvalCor.meta"
# add in genomic data
temp.results5.infoCpG <- cpgInfo(temp.results5$CpG, TxDb = "TxDb.Hsapiens.UCSC.hg19.knownGene")
temp.results.final.cox <- cbind(temp.results5, temp.results5.infoCpG, corrected = TRUE)

aems450k.meta.outcome.COX.withSMOKING <- temp.results.final.cox

head.style <- createStyle(textDecoration = "BOLD")
write.xlsx(aems450k.meta.outcome.COX.withSMOKING,
           file = paste0(OUT_loc, "/",Today,".aems450k.meta.plaque.Cox.2G.withSMOKING.xlsx"),
           creator = "Sander W. van der Laan",
           sheetName = "Results", headerStyle = head.style,
           row.names = FALSE, col.names = TRUE, overwrite = TRUE)

# Removing intermediates
cat("- Removing intermediate files...\n")
rm(temp.results, temp.results2, temp.results3, temp.results4, temp.results5, head.style,
   aems450k.meta.outcome.t, aems450k.meta.outcome.se, aems450k.meta.outcome.es,
   aems450k.meta.outcome.pvalsbcpm.b, aems450k.meta.outcome.bcpm.pvalsbcp,
   aems450k.meta.outcome.bcpm.pvalsp, aems450k.meta.outcome.bcpm.zbcpp, aems450k.meta.outcome.bcpm.zp,
   aems450k.meta.outcome.bcpm.lambdap, aems450k.meta.outcome.bcpm.lambdabcp, 
   temp.results5.infoCpG, temp.results.final.cox,
   hm450.controls, seqinfo.hg19,
   aems450k.meta.outcome.bcp.b, aems450k.meta.outcome.bcp.t, aems450k.meta.outcome.bcpm.b,
   AEData_2016.update.DF.CEA.TFvar, aems450k.meta.resultspfCGIQC.SmokerCurrent,
   aems450k1.MvaluesQCplaqueCleanTFvar, aems450k2.MvaluesQCplaqueCleanTFvar)

cat("===========================================================================================")
cat("SMOKING vs OUTCOME")

# Set up a dataframe to receive results
COX.results <- data.frame(matrix(NA, ncol = 12, nrow = 0))

for (i in 1:length(times)) {
  eptime = times[i]
  ep = endpoints[i]
  cat(paste0("* Analyzing the effect of current smoking on [",ep,"].\n"))
  cat(" - creating temporary DF for this work.\n")
  TEMP.DF = AEData_2016.update.DF.CEA
  cat(" - making a 'Surv' object and adding this to dataframe.\n")
  TEMP.DF$event <- as.integer(TEMP.DF[,ep] == "Excluded")
  TEMP.DF$y <- Surv(TEMP.DF[,eptime], TEMP.DF$event)
  fit <- survfit(y ~ SmokerCurrent, data = TEMP.DF)
  
  show(ggsurvplot(fit, data = TEMP.DF,
                  palette = c("#DB003F", "#1290D9"),
                  # palete = c("F59D10", "#DB003F", "#49A01D", "#1290D9"), 
                  linetype = c(1,2),
                  # linetype = c(1,2,3,4),
                  # conf.int = FALSE, conf.int.fill = "#595A5C", conf.int.alpha = 0.1,
                  pval = TRUE, pval.method = TRUE, pval.size = 4,
                  risk.table = TRUE, risk.table.y.text = FALSE, tables.y.text.col = TRUE, fontsize = 4,
                  censor = FALSE, 
                  cumevents = TRUE,
                  cumcensor = TRUE,
                  # ncensor.plot = TRUE, ncensor.plot.title = paste0("Censors of ",ep,""), ncensor.plot.height = 0.5,
                  # surv.plot.height = 0.2, 
                  legend = "right",
                  legend.title = "Current smoker",
                  legend.labs = c("no", "yes"),
                  title = paste0("Risk of ",ep,""), xlab = "Time [years]", font.main = c(16, "bold", "black")))
  dev.copy2pdf(file = paste0(COX_loc,"/",
                             Today,".ae.smoking.survival.",ep,".2G.pdf"), width = 12, height = 10, onefile = FALSE)
  dev.copy2eps(file = paste0(COX_loc,"/",
                             Today,".ae.smoking.survival.",ep,".2G.ps"), width = 12, height = 10, onefile = FALSE)
  
  cat(paste0("\n   > perform the Cox-regression fashizzle and plot it...\n"))
  ### Do Cox-regression and plot it
  # EWAS model
  # cox = coxph(Surv(TEMP.DF[,eptime], TEMP.DF$event) ~ SmokerCurrent+Age+Gender+Hospital, data = TEMP.DF)
  # coxplot = coxph(Surv(TEMP.DF[,eptime], TEMP.DF$event) ~ strata(SmokerCurrent)+Age+Gender+Hospital, data = TEMP.DF)
  # Full model
  cox = coxph(Surv(TEMP.DF[,eptime], TEMP.DF$event) ~ SmokerCurrent+Age+Gender+Hospital+DM.composite+Hypertension1+CAD_history+PAOD+TC_final+LDL_final+BMI+GFR_MDRD+Med.anticoagulants,
              data = TEMP.DF)
  coxplot = coxph(Surv(TEMP.DF[,eptime], TEMP.DF$event) ~
                    strata(SmokerCurrent)+Age+Gender+Hospital+DM.composite+Hypertension1+CAD_history+PAOD+TC_final+LDL_final+BMI+GFR_MDRD+Med.anticoagulants,
                  data = TEMP.DF)
  plot(survfit(coxplot), main = paste0("Cox proportional hazard of [",ep,"] per [",eptime,"]."),
       # ylim = c(0.2, 1), xlim = c(0,3), col = c("#595A5C", "#DB003F", "#1290D9"),
       ylim = c(0, 1), xlim = c(0,3), col = c("#DB003F", "#1290D9"),# Normal
       # ylim = c(0.7, 1), xlim = c(0,3), col = c("#DB003F", "#1290D9"), # Zoom
       lty = c(1,2), lwd = 2,
       ylab = "Suvival probability", xlab = "FU time [years]",
       mark.time = FALSE, axes = FALSE, bty = "n")
  legend("topright", 
         c("no", "yes"),
         title = "Current smoker",
         col = c("#DB003F", "#1290D9"),
         lty = c(1,2), lwd = 2, 
         bty = "n")
  axis(side = 1, at = seq(0, 3, by = 1))
  axis(side = 2, at = seq(0, 1, by = 0.2))
  # EWAS model
  # dev.copy2pdf(file = paste0(COX_loc,"/",
  #                            Today,".ae.smoking.cox.",ep,".2G.pdf"), width = 12, height = 10, onefile = FALSE)
  # dev.copy2eps(file = paste0(COX_loc,"/",
  #                            Today,".ae.smoking.cox.",ep,".2G.ps"), width = 12, height = 10, onefile = FALSE)
  # dev.copy2pdf(file = paste0(COX_loc,"/",
  #                            Today,".ae.smoking.cox.",ep,".2G.zoom.pdf"), width = 12, height = 10, onefile = FALSE)
  # dev.copy2eps(file = paste0(COX_loc,"/",
  #                            Today,".ae.smoking.cox.",ep,".2G.zoom.ps"), width = 12, height = 10, onefile = FALSE)
  
  
  # Full model
  dev.copy2pdf(file = paste0(COX_loc,"/",
                             Today,".ae.smoking.cox.",ep,".2G.fullmodel.pdf"), width = 12, height = 10, onefile = FALSE)
  dev.copy2eps(file = paste0(COX_loc,"/",
                             Today,".ae.smoking.cox.",ep,".2G.fullmodel.ps"), width = 12, height = 10, onefile = FALSE)
  # dev.copy2pdf(file = paste0(COX_loc,"/",
  #                            Today,".ae.smoking.cox.",ep,".2G.fullmodel.zoom.pdf"), width = 12, height = 10, onefile = FALSE)
  # dev.copy2eps(file = paste0(COX_loc,"/",
  #                            Today,".ae.smoking.cox.",ep,".2G.fullmodel.zoom.ps"), width = 12, height = 10, onefile = FALSE)
  
  show(summary(cox))
  
  # temporary matrix and dataframes to receive results
  COX.results.TEMP <- data.frame(matrix(NA, ncol = 12, nrow = 0))
  COX.results.TEMP[1,] = COX.STAT(cox, "AEDB2016", ep, "CurrentSmoking")
  COX.results = rbind(COX.results, COX.results.TEMP)
}

cat("- Edit the column names...\n")
colnames(COX.results) = c("Dataset", "Outcome", "Trait", 
                          "Beta", "s.e.m.", 
                          "HR", "low95CI", "up95CI", 
                          "Z-value", "P-value", "SampleSize", "N_events")

cat("- Correct the variable types...\n")
COX.results$Beta <- as.numeric(COX.results$Beta)
COX.results$s.e.m. <- as.numeric(COX.results$s.e.m.)
COX.results$HR <- as.numeric(COX.results$HR)
COX.results$low95CI <- as.numeric(COX.results$low95CI)
COX.results$up95CI <- as.numeric(COX.results$up95CI)
COX.results$`Z-value` <- as.numeric(COX.results$`Z-value`)
COX.results$`P-value` <- as.numeric(COX.results$`P-value`)
COX.results$SampleSize <- as.numeric(COX.results$SampleSize)
COX.results$N_events <- as.numeric(COX.results$N_events)

# EWAS model
# AEData_2016.update.DF.CEA.COX.resultsSMOKING <- COX.results
# Full model
AEData_2016.update.DF.CEA.COX.resultsSMOKING.fullmodel <- COX.results

# Removing intermediates
cat("- Removing intermediate files...\n")
rm(TEMP.DF, fit, cox, coxplot, COX.results, COX.results.TEMP, ep, eptime)

# Save the data
cat("- Writing Cox-regression fashizzle results to Excel-file...\n")
head.style <- createStyle(textDecoration = "BOLD")
# EWAS model
# write.xlsx(AEData_2016.update.DF.CEA.COX.resultsSMOKING,
#            file = paste0(OUT_loc, "/",Today,".AE.Smoking.Cox.2G.EWASmodel.xlsx"),
#            creator = "Sander W. van der Laan",
#            sheetName = "Results", headerStyle = head.style,
#            row.names = FALSE, col.names = TRUE, overwrite = TRUE)
# Full model
write.xlsx(AEData_2016.update.DF.CEA.COX.resultsSMOKING.fullmodel,
           file = paste0(OUT_loc, "/",Today,".AE.Smoking.Cox.2G.Fullmodel.xlsx"),
           creator = "Sander W. van der Laan",
           sheetName = "Results", headerStyle = head.style,
           row.names = FALSE, col.names = TRUE, overwrite = TRUE)
rm(head.style)

cat("\n===========================================================================================")
cat("SAVE THE DATA")

save.image(paste0(ANALYSIS_loc,"/",Today,".aems450k.meta.analysis.ewas.",EWAS_trait,".COX_regressions.RData"))

