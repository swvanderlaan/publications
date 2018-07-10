cat("===========================================================================================
                META-ANALYSIS ATHERO-EXPRESS METHYLATION STUDIES 450K 1 & 2
            -- Linear & Logistic regressions of CpGs vs. Plaque Phenotypes --
    
    Version:      v2.8
    
    Last update:  2018-06-27
    Written by:   Sander W. van der Laan (s.w.vanderlaan-2@umcutrecht.nl);
    Marten A. Siemelink
    
    Description:  Script to perform regression analyses of DNA methylation vs. plaque
    phenotypes.
    
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

cat("\n  - loading EWAS analyses *REPLICATED* results...\n")
load(paste0(OUT_loc,"/20180207.SmokerCurrent.aems450k.meta.resultspfCGIQC.RData"))

aems450k.meta.resultspfCGIQC = aems450k.meta.resultspfCGIQC.SmokerCurrent

cat("===========================================================================================")
cat("CpGs vs PLAQUE PHENOTYPES")

cat("* Subset the regions of significant CpGs *after* REPLICATION...")
cat("- listing all top regions; we're being greedy and will analyze each CpG mapped to gene X...")
head(aems450k.meta.resultspfCGIQC[order(aems450k.meta.resultspfCGIQC$PvalCorAdj.meta),], 10)
ahrr <- aems450k.meta.resultspfCGIQC$CpG[grep("AHRR", aems450k.meta.resultspfCGIQC$SYMBOL)]
itpk1 <- aems450k.meta.resultspfCGIQC$CpG[grep("ITPK1", aems450k.meta.resultspfCGIQC$SYMBOL)]

hilight.top <- c(ahrr, itpk1)

cat("- listing all EWAS significant probes, with meta-p < 0.05...")
head(aems450k.meta.resultspfCGIQC[order(aems450k.meta.resultspfCGIQC$PvalCorAdj.meta),], 11)
n_ewas_sign = nrow(aems450k.meta.resultspfCGI) # these are the total number of tested CpGs
n_ewas_sign_rep = nrow(aems450k.meta.resultspfCGIQC) # these are the total number of tested CpGs AFTER replication
p_ewas_sign = 0.05/n_ewas_sign
p_ewas_sug = (0.05/n_ewas_sign) * 10
cat(paste0("  > EWAS significance threshold for ", format(n_ewas_sign, big.mark = ","),
" CpG probes equals to p-value = ",format(p_ewas_sign, digit = 3),". After replication there are ", 
format(n_ewas_sign_rep, big.mark = ","),
" CpG probes with p-value < 0.05. We had set the suggestive significance threshold to p-value = ",
format(p_ewas_sug, digit = 3),"."))

cat("- subsetting the data...")
aems450k1.MvaluesQCplaqueCleanTop <- aems450k1.MvaluesQCplaqueCleanTFvar[hilight.top, ]
aems450k2.MvaluesQCplaqueCleanTop <- aems450k2.MvaluesQCplaqueCleanTFvar[hilight.top, ]

cat("\n * Run the analysis.")
# first combine all the phenotypes
TRAITS.CON = c("MacrophagesPercBC", "SMCPercBC", "VesselDensityBC")

TRAITS.BIN = c("Calc.bin", "Collagen.bin", "Fat.bin_40", "IPH.bin")

TRAITS.PLAQUE <- c(TRAITS.BIN, TRAITS.CON)

for (trait in 1:length(TRAITS.PLAQUE)) {
  TRAIT = TRAITS.PLAQUE[trait]
  cat(paste0("\n - analyzing: [ ",TRAIT," ] for AEMS450K1."))
  currentSE1 = aems450k1.MvaluesQCplaqueCleanTop
  cat(paste0("\n - the number of CpG vs samples to analyse initially..."))
  print(dim(currentSE1))
  metadata(currentSE1)$formula <- as.formula(paste0("~",TRAIT," + Sample_Sex + Age + Hospital"))
  currentSE1.covariates <- get_all_vars(metadata(currentSE1)$formula, data = colData(currentSE1))
  currentSE1.nas <- apply(currentSE1.covariates, 1, anyNA)
  currentSE1 <- currentSE1[, !currentSE1.nas]
  cat(paste0("\n - the number of are CpG vs samples to analyse after removing any NAs..."))
  print(dim(currentSE1))
  currentSE1.designp <- model.matrix(metadata(currentSE1)$formula, data = colData(currentSE1))
  currentSE1.datap <- assays(currentSE1)$data
  currentSE1.fitp <- limma::lmFit(currentSE1.datap, currentSE1.designp)
  cat("\n  > calculate T-statistic = beta/standard error...")
  currentSE1.effectsizep <- currentSE1.fitp$coefficients
  currentSE1.SEp <- currentSE1.fitp$stdev.unscaled * currentSE1.fitp$sigma
  currentSE1.tstatp <- currentSE1.effectsizep/currentSE1.SEp
  cat("\n  > calculate p-value from T-statistic, and adjust using FDR (Benjamin-Hochberg)...")
  currentSE1.pvalp <- 2*pnorm(-abs(currentSE1.tstatp[,2]))
  currentSE1.padjp <- p.adjust(sort(currentSE1.pvalp, decreasing = FALSE), method = "BH")
  print(head(currentSE1.padjp[currentSE1.padjp < 0.05]))
  
  cat(paste0("\n - analyzing: [ ",TRAIT," ] for AEMS450K2."))
  currentSE2 = aems450k2.MvaluesQCplaqueCleanTop
  cat(paste0("\n - there are CpG vs samples to analyse initially..."))
  print(dim(currentSE2))
  metadata(currentSE2)$formula <- as.formula(paste0("~",TRAIT," + Sample_Sex + Age + Hospital"))
  currentSE2.covariates <- get_all_vars(metadata(currentSE2)$formula, data = colData(currentSE2))
  currentSE2.nas <- apply(currentSE2.covariates, 1, anyNA)
  currentSE2 <- currentSE2[, !currentSE2.nas]
  cat(paste0("\n - there are CpG vs samples to analyse after removing any NAs..."))
  print(dim(currentSE2))
  currentSE2.designp <- model.matrix(metadata(currentSE2)$formula, data = colData(currentSE2))
  currentSE2.datap <- assays(currentSE2)$data
  currentSE2.fitp <- limma::lmFit(currentSE2.datap, currentSE2.designp)
  cat("\n  > calculate T-statistic = beta/standard error...")
  currentSE2.effectsizep <- currentSE2.fitp$coefficients
  currentSE2.SEp <- currentSE2.fitp$stdev.unscaled * currentSE2.fitp$sigma
  currentSE2.tstatp <- currentSE2.effectsizep/currentSE2.SEp
  cat("\n  > calculate p-value from T-statistic, and adjust using FDR (Benjamin-Hochberg)...")
  currentSE2.pvalp <- 2*pnorm(-abs(currentSE2.tstatp[,2]))
  currentSE2.padjp <- p.adjust(sort(currentSE2.pvalp, decreasing = FALSE), method = "BH")
  print(head(currentSE2.padjp[currentSE2.padjp < 0.05]))
  
  cat("\n * Correct for bias and inflation using 'BACON'. 
      (ref: M. van Iterson, E.W. van Zwet, and B.T. Heijmans. 2017. Genome Biol. 18 (1): 19.).")
  library("bacon")
  # plaque
  currentSE.meta.t <- cbind(limmaT.aems450k1 = currentSE1.tstatp[,2], limmaT.aems450k2 = currentSE2.tstatp[,2])
  currentSE.meta.es <- cbind(limmaB.aems450k1 = currentSE1.effectsizep[,2], limmaB.aems450k2 = currentSE2.effectsizep[,2])
  currentSE.meta.se <- cbind(limmaSE.aems450k1 = currentSE1.SEp[,2], limmaSE.aems450k2 = currentSE2.SEp[,2])
  
  library(BiocParallel)
  register(MulticoreParam(8))
  currentSE.meta.bcp.t <- bacon(teststatistics = currentSE.meta.t, 
                                effectsizes = NULL, standarderrors = NULL,
                                verbose = TRUE) 
  currentSE.meta.bcp.b <- bacon(teststatistics = NULL,
                                effectsizes = currentSE.meta.es, standarderrors = currentSE.meta.se,
                                verbose = TRUE)
  cat("\n   - Bacon model result.")
  print(currentSE.meta.bcp.t)
  print(currentSE.meta.bcp.b)
  cat("\n   - estimates of the mixture fit")
  print(estimates(currentSE.meta.bcp.t))
  print(estimates(currentSE.meta.bcp.b))
  cat("\n   - the inflation")
  print(inflation(currentSE.meta.bcp.t))
  print(inflation(currentSE.meta.bcp.b))
  cat("\n   - the bias")
  print(bias(currentSE.meta.bcp.t))
  print(bias(currentSE.meta.bcp.b))
  
  cat("\n * Perform meta-analysis.")
  currentSE.meta.bcpm.b <- meta(currentSE.meta.bcp.b)
  currentSE.meta.pvalsbcpm.b <- pval(currentSE.meta.bcpm.b)
  head(currentSE.meta.pvalsbcpm.b[order(currentSE.meta.pvalsbcpm.b[,1]),])
  print(bacon::topTable(currentSE.meta.bcpm.b))
  
  cat("\n * Uncorrected results.")
  currentSE.meta.pvalsp <- pval(currentSE.meta.bcpm.b, corrected = FALSE)
  head(currentSE.meta.pvalsp[order(currentSE.meta.pvalsp[,3]),])
  currentSE.meta.zp = qnorm(currentSE.meta.pvalsp[,1]/2)
  currentSE.meta.lambdap = round(median(currentSE.meta.zp^2)/0.4549364,3)
  cat(paste0("\n  - lambda is: [",currentSE.meta.lambdap,"]."))
  
  cat("\n * Corrected results.")
  currentSE.meta.pvalsbcp <- pval(currentSE.meta.bcpm.b, corrected = TRUE)
  head(currentSE.meta.pvalsbcp[order(currentSE.meta.pvalsbcp[,3]),])
  currentSE.meta.zbcp = qnorm(currentSE.meta.pvalsbcp[,1]/2)
  currentSE.meta.lambdabcp = round(median(currentSE.meta.zbcp^2)/0.4549364,3)
  cat(paste0("\n  - lambda is: [",currentSE.meta.lambdabcp,"]."))
  
  cat("\n * Annotating corrected results for plotting and other purposes.")
  cat("\n   - annotating...")
  # plaque
  currentSE.infobcp <- cpgInfo(rownames(currentSE.meta.pvalsbcp), TxDb = "TxDb.Hsapiens.UCSC.hg19.knownGene")
  print(head(currentSE.infobcp))
  cat("\n   - merging with results...")
  # add in effect sizes
  temp.results <- merge(currentSE.infobcp, es(currentSE.meta.bcpm.b), by = "row.names")
  rownames(temp.results) <- temp.results$Row.names
  temp.results2 <- temp.results[, -which(names(temp.results) %in% c("Row.names"))]
  names(temp.results2)[names(temp.results2) == "limmaB.aems450k1"] <- "EffectSize.aems450k1"
  names(temp.results2)[names(temp.results2) == "limmaB.aems450k2"] <- "EffectSize.aems450k2"
  names(temp.results2)[names(temp.results2) == "meta"] <- "EffectSize.meta"
  # add in standard errors
  temp.results3 <- merge(temp.results2, se(currentSE.meta.bcpm.b), by = "row.names")
  rownames(temp.results3) <- temp.results3$Row.names
  temp.results4 <- temp.results3[, -which(names(temp.results3) %in% c("Row.names"))]
  names(temp.results4)[names(temp.results4) == "limmaSE.aems450k1"] <- "SE.aems450k1"
  names(temp.results4)[names(temp.results4) == "limmaSE.aems450k2"] <- "SE.aems450k2"
  names(temp.results4)[names(temp.results4) == "meta"] <- "SE.meta"
  # add in uncorrected p-values
  temp.results5 <- merge(temp.results4, currentSE.meta.pvalsp, by = "row.names")
  rownames(temp.results5) <- temp.results5$Row.names
  temp.results6 <- temp.results5[, -which(names(temp.results5) %in% c("Row.names"))]
  names(temp.results6)[names(temp.results6) == "limmaB.aems450k1"] <- "Pval.aems450k1"
  names(temp.results6)[names(temp.results6) == "limmaB.aems450k2"] <- "Pval.aems450k2"
  temp.results7 <- temp.results6[, -which(names(temp.results6) %in% c("meta"))]
  # add in corrected p-values
  currentSE.meta.resultsp <- merge(temp.results7, currentSE.meta.pvalsbcp, by = "row.names")
  rownames(currentSE.meta.resultsp) <- currentSE.meta.resultsp$Row.names
  names(currentSE.meta.resultsp)[names(currentSE.meta.resultsp) == "Row.names"] <- "CpG"
  names(currentSE.meta.resultsp)[names(currentSE.meta.resultsp) == "limmaB.aems450k1"] <- "PvalCor.aems450k1"
  names(currentSE.meta.resultsp)[names(currentSE.meta.resultsp) == "limmaB.aems450k2"] <- "PvalCor.aems450k2"
  names(currentSE.meta.resultsp)[names(currentSE.meta.resultsp) == "meta"] <- "PvalCor.meta"
  currentSE.meta.resultsp$PvalCorAdj.meta <- p.adjust(currentSE.meta.resultsp$PvalCor.meta, method = "BH")
  print(head(currentSE.meta.resultsp))
  print(str(currentSE.meta.resultsp))
  print(dim(currentSE.meta.resultsp))
  
  cat("\n   - removing intermediate file...")
  rm(temp.results, temp.results2, temp.results3, temp.results4, temp.results5, temp.results6, temp.results7)
  
  cat("\n   - recreating the chromosomes...")
  currentSE.list.chr <- levels(currentSE.meta.resultsp$seqnames)[1:25]
  currentSE.list.chr
  currentSE.meta.resultspf <- currentSE.meta.resultsp[(currentSE.meta.resultsp$seqnames %in% currentSE.list.chr), ]
  currentSE.meta.resultspf$chr <- gsub("chr", "",currentSE.meta.resultspf$seqnames)
  currentSE.meta.resultspf$chr <- gsub("X", "23",currentSE.meta.resultspf$chr)
  currentSE.meta.resultspf$chr <- gsub("Y", "24",currentSE.meta.resultspf$chr)
  currentSE.meta.resultspf$chr <- gsub("XY", "25",currentSE.meta.resultspf$chr)
  currentSE.meta.resultspf$chr <- gsub("M", "26",currentSE.meta.resultspf$chr)
  currentSE.meta.resultspf$chr <- as.numeric(currentSE.meta.resultspf$chr)
  
  cat("\n* Saving results...")
  cat("\n  - writing plaque results...\n")
  head.style <- createStyle(textDecoration = "BOLD")
  write.xlsx(currentSE.meta.resultspf, 
             file = paste0(OUT_loc, "/", Today,".aems450k.meta.analysis.",EWAS_trait,".CpG_vs_",TRAIT,".xlsx"),
             creator = "Sander W. van der Laan",
             sheetName = paste0(TRAIT), tabColour = uithof_color[trait], headerStyle = head.style,
             row.names = FALSE, col.names = TRUE, overwrite = TRUE)
  cat("\n*******************************************************************************************\n\n")
}
rm(currentSE1, currentSE1.covariates, currentSE1.nas,
   currentSE2, currentSE2.covariates, currentSE2.nas,
   currentSE1.effectsizep, currentSE1.SEp, currentSE1.tstatp, 
   currentSE1.fitp, currentSE1.padjp, currentSE1.pvalp,
   currentSE2.effectsizep, currentSE2.SEp, currentSE2.tstatp, 
   currentSE2.fitp, currentSE2.padjp, currentSE2.pvalp,
   currentSE.meta.t, currentSE.meta.se, currentSE.meta.es,
   currentSE.meta.bcp.t, currentSE.meta.bcp.b,
   currentSE.meta.resultspf, currentSE.meta.resultsp,
   currentSE.meta.bcpm.b, currentSE.meta.pvalsbcpm.b, 
   currentSE.meta.pvalsp, currentSE.meta.zp, currentSE.meta.lambdap, 
   currentSE.meta.pvalsbcp, currentSE.meta.zbcp,currentSE.meta.lambdabcp,
   currentSE.list.chr, currentSE.infobcp, head.style,
   currentSE1.datap, currentSE1.designp, currentSE2.datap, currentSE2.designp,
   hm450.controls, seqinfo.hg19,
   aems450k1.MvaluesQCplaqueCleanTFvar, aems450k2.MvaluesQCplaqueCleanTFvar)

cat("\n===========================================================================================")
cat("SAVE THE DATA")

save.image(paste0(ANALYSIS_loc,"/",Today,".aems450k.meta.analysis.ewas.",EWAS_trait,".CpG_vs_PlaquePheno.RData"))

