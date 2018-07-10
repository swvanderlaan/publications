cat("===========================================================================================
                META-ANALYSIS ATHERO-EXPRESS METHYLATION STUDIES 450K 1 & 2
    
    Version:      v2.8
    
    Last update:  2018-01-31
    Written by:   Sander W. van der Laan (s.w.vanderlaan-2@umcutrecht.nl);
    Marten A. Siemelink
    
    Description:  Script to do meta-analysis of Athero-Express Methylation 
    Study 450K 1 (2013) & 2 (2016) EWAS on smoking.
    Based on DNAmArray: 
    - https://molepi.github.io/DNAmArray_workflow/index.html
    - https://github.com/molepi/DNAmArray
    
    Minimum requirements: R version 3.4.1 (2017-06-30) -- 'Single Candle', Mac OS X El Capitan
    
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
    eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE, repos = \"http://cran-mirror.cs.uu.nl/\")", x)))
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
# install.packages.auto("RMySQL") # install this one from CRAN!
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

cat("\n* Manhattan plotting function, based on library(\"qqman\")...")
manhattan.uithof <- function(x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", 
                              col = c("gray10", "gray60"), chrlabs = NULL, 
                              suggestiveline = -log10(1e-05),  genomewideline = -log10(5e-08), 
                              highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = TRUE, ...) {
  CHR = BP = P = index = NULL
  if (!(chr %in% names(x))) 
    stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) 
    stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) 
    stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x))) 
    warning(paste("No SNP/CpG column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]])) 
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) 
    stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) 
    stop(paste(p, "column should be numeric."))
  d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
  if (!is.null(x[[snp]])) 
    d = transform(d, SNP = x[[snp]])
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  }
  else {
    d$logp <- d$P
  }
  d$pos = NA
  d$index = NA
  ind = 0
  for (i in unique(d$CHR)) {
    ind = ind + 1
    d[d$CHR == i, ]$index = ind
  }
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    d$pos = d$BP
    ticks = floor(length(d$pos))/2 + 1
    xlabel = paste("Chromosome", unique(d$CHR), "position")
    labs = ticks
  }
  else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + tail(subset(d, index == 
                                            i - 1)$BP, 1)
        d[d$index == i, ]$pos = d[d$index == i, ]$BP + 
          lastbase
      }
      ticks = c(ticks, (min(d[d$index == i, ]$pos) + max(d[d$index == 
                                                             i, ]$pos))/2 + 1)
    }
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
                   las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0, 
                                                                     ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
                                            names(dotargs)]))
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  if (nchr == 1) {
    axis(1, ...)
  }
  else {
    axis(1, at = ticks, labels = labs, ...)
  }
  col = rep(col, max(d$CHR))
  if (nchr == 1) {
    with(d, points(pos, logp, pch = 20, col = col[1], ...))
  }
  else {
    icol = 1
    for (i in unique(d$index)) {
      with(d[d$index == unique(d$index)[i], ], points(pos, 
                                                      logp, col = col[icol], pch = 20, ...))
      icol = icol + 1
    }
  }
  if (suggestiveline) 
    abline(h = suggestiveline, col = "#1290D9", lty = 2, lwd = 2)
  if (genomewideline) 
    abline(h = genomewideline, col = "#E55738", lty = 2, lwd = 2)
  if (!is.null(highlight)) {
    if (any(!(highlight %in% d$SNP))) 
      warning("You're trying to highlight SNPs/CpGs that don't exist in your results.")
    d.highlight = d[which(d$SNP %in% highlight), ]
    with(d.highlight, points(pos, logp, col = "#595A5C", pch = 20, 
                             ...))
  }
  if (!is.null(annotatePval)) {
    topHits = subset(d, P <= annotatePval)
    par(xpd = TRUE)
    if (annotateTop == FALSE) {
      with(subset(d, P <= annotatePval), textxy(pos, -log10(P), 
                                                offset = 0.625, labs = topHits$SNP, cex = 1.25), 
           ...)
    }
    else {
      topHits <- topHits[order(topHits$P), ]
      topSNPs <- NULL
      for (i in unique(topHits$CHR)) {
        chrSNPs <- topHits[topHits$CHR == i, ]
        topSNPs <- rbind(topSNPs, chrSNPs[1, ])
      }
      library("calibrate")
      textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, 
             labs = topSNPs$SNP, cex = 1.25, ...)
    }
  }
  par(xpd = FALSE)
}

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

cat("\n  - loading B/Mvalues of plaque samples...")
load(paste0(INP_AEMS450K1_loc,"/20171229.aems450k1.BvaluesQCIMP.plaque.RData"))
load(paste0(INP_AEMS450K1_loc,"/20171229.aems450k1.MvaluesQCIMP.plaque.RData"))

load(paste0(INP_AEMS450K2_loc,"/20171229.aems450k2.BvaluesQCIMP.plaque.RData"))
load(paste0(INP_AEMS450K2_loc,"/20171229.aems450k2.MvaluesQCIMP.plaque.RData"))

cat("===========================================================================================")
cat("\n[ META-ANALYSIS of EPIGENOME-WIDE ASSOCIATION STUDY on ",EWAS_trait," in PLAQUE in AEMS450K1 & 2 ]")
# Reference: https://molepi.github.io/DNAmArray_workflow/06_EWAS.html

cat("===========================================================================================")
cat("\n[ SANNITY CHECKING THE DATA ]")

cat("\n * AEMS450K1, in plaque...")
pdf(paste0(QC_loc,"/",Today,".aems450k1.",EWAS_trait,".plaque.MethylationDensity.pdf"),
    width = 12, height = 8, onefile = TRUE)
  par(mfrow = c(1,2), oma = c(0, 0, 2, 0))
  densityPlot(assays(aems450k1.BvaluesQCplaque)$data, sampGroups = aems450k1.BvaluesQCplaque$SmokerCurrent, main = "Beta-values", 
              legend = FALSE, 
              xlab = "Beta-values", 
              pal = c("#9FC228", "#E55738"), 
              bty = "n")
  legend("topright", legend = levels(factor(aems450k1.BvaluesQCplaque$SmokerCurrent)), 
         text.col = c("#9FC228", "#E55738"), 
         lty = 1, lwd = 1, col = c("#9FC228", "#E55738"), 
         bty = "n")
  densityPlot(assays(aems450k1.MvaluesQCplaque)$data, sampGroups = aems450k1.MvaluesQCplaque$SmokerCurrent, main = "M-values", 
              legend = FALSE, 
              xlab = "M-values", 
              pal = c("#9FC228", "#E55738"), 
              bty = "n")
  legend("topright", legend = levels(factor(aems450k1.MvaluesQCplaque$SmokerCurrent)), 
         text.col = c("#9FC228", "#E55738"), 
         lty = 1, lwd = 1, col = c("#9FC228", "#E55738"), 
         bty = "n")
  mtext(paste0("Overall methylation density (plaque, ",SUBPROJECTNAME1,")"), outer = TRUE, cex = 1.5)
  par(mfrow = c(1,1), oma = c(0,0,0,0))
dev.off()

cat("\n * AEMS450K2, in plaque...")
pdf(paste0(QC_loc,"/",Today,".aems450k2.",EWAS_trait,".plaque.MethylationDensity.pdf"),
    width = 12, height = 8, onefile = TRUE)
  par(mfrow = c(1,2), oma = c(0, 0, 2, 0))
  densityPlot(assays(aems450k2.BvaluesQCplaque)$data, sampGroups = aems450k2.BvaluesQCplaque$SmokerCurrent, main = "Beta-values", 
              legend = FALSE, 
              xlab = "Beta-values", 
              pal = c("#9FC228", "#E55738"), 
              bty = "n")
  legend("topright", legend = levels(factor(aems450k2.BvaluesQCplaque$SmokerCurrent)), 
         text.col = c("#9FC228", "#E55738"), 
         lty = 1, lwd = 1, col = c("#9FC228", "#E55738"), 
         bty = "n")
  densityPlot(assays(aems450k2.MvaluesQCplaque)$data, sampGroups = aems450k2.MvaluesQCplaque$SmokerCurrent, main = "M-values", 
              legend = FALSE, 
              xlab = "M-values", 
              pal = c("#9FC228", "#E55738"), 
              bty = "n")
  legend("topright", legend = levels(factor(aems450k2.MvaluesQCplaque$SmokerCurrent)), 
         text.col = c("#9FC228", "#E55738"), 
         lty = 1, lwd = 1, col = c("#9FC228", "#E55738"), 
         bty = "n")
  mtext(paste0("Overall methylation density (plaque, ",SUBPROJECTNAME2,")"), outer = TRUE, cex = 1.5)
  par(mfrow = c(1,1), oma = c(0,0,0,0))
dev.off()

cat("\nRemoving BvaluesQCIMP objects - as we don't use these anymore.")
rm(aems450k1.BvaluesQCplaque, aems450k2.BvaluesQCplaque)

cat("===========================================================================================")
cat("\n[ CONTINUE (META-)ANALYSIS AEMS450K STUDIES ]")

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
cat("  > AEMS450K2...")
aems450k2.MvaluesQCplaqueClean <- aems450k2.MvaluesQCplaque
dim(aems450k2.MvaluesQCplaque)
metadata(aems450k2.MvaluesQCplaqueClean)$formula <- ~SmokerCurrent + Sample_Sex + Age + Hospital

cat("\n  - Next, we extract those samples having a complete set of covariates. 
    Notice that we subset the SummarizedExperiment-object!")
cat("  > AEMS450K1...")
aems450k1.covariates <- get_all_vars(metadata(aems450k1.MvaluesQCplaqueClean)$formula, data = colData(aems450k1.MvaluesQCplaqueClean))
aems450k1.nas <- apply(aems450k1.covariates, 1, anyNA)
aems450k1.MvaluesQCplaqueClean <- aems450k1.MvaluesQCplaqueClean[, !aems450k1.nas]
dim(aems450k1.MvaluesQCplaqueClean)
cat("  > AEMS450K2...")
aems450k2.covariates <- get_all_vars(metadata(aems450k2.MvaluesQCplaqueClean)$formula, data = colData(aems450k2.MvaluesQCplaqueClean))
aems450k2.nas <- apply(aems450k2.covariates, 1, anyNA)
aems450k2.MvaluesQCplaqueClean <- aems450k2.MvaluesQCplaqueClean[, !aems450k2.nas]
dim(aems450k2.MvaluesQCplaqueClean)

cat("\n  - Remove probes containing 
    SNPs or does mapping to multiple locations (ref: W. Zhou, P.W. Laird, and H. Shen. 2016. Nucleic Acids Res.).")
cat("  > next probes with SNPs...")
data(hm450.manifest.pop.GoNL) ##From DNAmArray
#hm450.manifest.pop.GoNL
hm450.manifest.pop.GoNL <- hm450.manifest.pop.GoNL[!is.na(hm450.manifest.pop.GoNL$MASK.general.GoNL) &
                                                     hm450.manifest.pop.GoNL$MASK.general.GoNL == TRUE, ]
cat("  > AEMS450K1...")
aems450k1.MvaluesQCplaqueClean <- aems450k1.MvaluesQCplaqueClean[!(names(aems450k1.MvaluesQCplaqueClean) %in% names(hm450.manifest.pop.GoNL)),]  
cat("  > AEMS450K2...")
aems450k2.MvaluesQCplaqueClean <- aems450k2.MvaluesQCplaqueClean[!(names(aems450k2.MvaluesQCplaqueClean) %in% names(hm450.manifest.pop.GoNL)),]  

cat("\n  - Subsetting AEMS450K1/AEMS450K2.")
cat("  > get a list of probes that difference between the two datasets")
aems450k1.ranges <- unlist(names(rowRanges(aems450k1.MvaluesQCplaqueClean)))
aems450k2.ranges <- unlist(names(rowRanges(aems450k2.MvaluesQCplaqueClean)))

cat("  > what is in both?")
aems450k.meta.intersect <- intersect(aems450k1.ranges,aems450k2.ranges)
length(aems450k.meta.intersect)

cat("  > what is different in AEMS450K1...?")
aems450k1.meta.difs <- setdiff(aems450k1.ranges,aems450k2.ranges)
length(aems450k1.meta.difs)

cat("  > what is different in AEMS450K2...?")
aems450k2.meta.difs <- setdiff(aems450k2.ranges,aems450k1.ranges)
length(aems450k2.meta.difs)

cat("  > subsetting...")
aems450k1.MvaluesQCplaqueClean <- aems450k1.MvaluesQCplaqueClean[(names(aems450k1.MvaluesQCplaqueClean) %in% aems450k.meta.intersect),]  
aems450k2.MvaluesQCplaqueClean <- aems450k2.MvaluesQCplaqueClean[(names(aems450k2.MvaluesQCplaqueClean) %in% aems450k.meta.intersect),]  

cat("  > we started with:")
dim(aems450k1.MvaluesQCplaque)  # 483731 CpGs, 485 samples
dim(aems450k2.MvaluesQCplaque)  # 484249 CpGs, 190 samples
# rowRanges(aems450k1.MvaluesQCplaqueClean)
cat("  > we end up with:")
dim(aems450k1.MvaluesQCplaqueClean)  # 443084 CpGs, 477 samples
dim(aems450k2.MvaluesQCplaqueClean)  # 443084 CpGs, 187 samples
# rowRanges(aems450k2.MvaluesQCplaqueClean)

cat("\n * Run the analysis.")
install.packages.auto("limma")
cat(" - AEMS450K1 *with* hospital.")
aems450k1.designph <- model.matrix(metadata(aems450k1.MvaluesQCplaqueClean)$formula, data = colData(aems450k1.MvaluesQCplaqueClean))
aems450k1.dataph <- assays(aems450k1.MvaluesQCplaqueClean)$data
aems450k1.fitph <- limma::lmFit(aems450k1.dataph, aems450k1.designph)
cat("  > calculate T-statistic = beta/standard error...")
aems450k1.effectsizeph <- aems450k1.fitph$coefficients
aems450k1.SEph <- aems450k1.fitph$stdev.unscaled * aems450k1.fitph$sigma
aems450k1.tstatph <- aems450k1.effectsizeph/aems450k1.SEph
cat("  > calculate p-value from T-statistic, and adjust using FDR (Benjamin-Hochberg)...")
aems450k1.pvalph <- 2*pnorm(-abs(aems450k1.tstatph[,2]))
aems450k1.padjph <- p.adjust(sort(aems450k1.pvalph, decreasing = FALSE), method = "BH")
head(aems450k1.padjph[aems450k1.padjph < 0.05])
# summary
summary(aems450k1.fitph)

cat(" - AEMS450K2 *with* hospital.")
aems450k2.designph <- model.matrix(metadata(aems450k2.MvaluesQCplaqueClean)$formula, data = colData(aems450k2.MvaluesQCplaqueClean))
aems450k2.dataph <- assays(aems450k2.MvaluesQCplaqueClean)$data
aems450k2.fitph <- limma::lmFit(aems450k2.dataph, aems450k2.designph)
cat("  > calculate T-statistic = beta/standard error...")
aems450k2.effectsizeph <- aems450k2.fitph$coefficients
aems450k2.SEph <- aems450k2.fitph$stdev.unscaled * aems450k2.fitph$sigma
aems450k2.tstatph <- aems450k2.effectsizeph/aems450k2.SEph
cat("  > calculate p-value from T-statistic, and adjust using FDR (Benjamin-Hochberg)...")
aems450k2.pvalph <- 2*pnorm(-abs(aems450k2.tstatph[,2]))
aems450k2.padjph <- p.adjust(sort(aems450k2.pvalph, decreasing = FALSE), method = "BH")
head(aems450k2.padjph[aems450k2.padjph < 0.05])
# summary
summary(aems450k2.fitph)

cat("\n * Correct for bias and inflation using 'BACON'. 
    (ref: M. van Iterson, E.W. van Zwet, and B.T. Heijmans. 2017. Genome Biol. 18 (1): 19.).")
library("bacon")
# plaque *with* hospital
aems450k.meta.th <- cbind(limmaT.aems450k1 = aems450k1.tstatph[,2], limmaT.aems450k2 = aems450k2.tstatph[,2])
aems450k.meta.esh <- cbind(limmaB.aems450k1 = aems450k1.effectsizeph[,2], limmaB.aems450k2 = aems450k2.effectsizeph[,2])
aems450k.meta.seh <- cbind(limmaSE.aems450k1 = aems450k1.SEph[,2], limmaSE.aems450k2 = aems450k2.SEph[,2])

cat("\n * Getting a handle on the test statistics of both datasets.")
cat("  - of the data 'as-is'...")
summary(aems450k.meta.th)
summary(aems450k.meta.esh)
summary(aems450k.meta.seh)
cat("  - of the data 'absolutely'...")
summary(abs(aems450k.meta.th))
summary(abs(aems450k.meta.esh))
summary(abs(aems450k.meta.seh))

# plaque *with* hospital
mean.betah = mean(aems450k.meta.esh)
sd.betah = sd(aems450k.meta.esh)

neg.cutoffh = mean.betah - 4 * sd.betah
pos.cutoffh = mean.betah + 4 * sd.betah

aems450k1.plothisth <- hist(aems450k.meta.esh[,1], breaks = 100, border = FALSE)
aems450k2.plothisth <- hist(aems450k.meta.esh[,2], breaks = 100, border = FALSE)

aems450k1.col <- adjustcolor(229,87,56, alpha.f = 0.25)
aems450k2.col <- adjustcolor(18,144,217, alpha.f = 0.25)

# PLOT FOR SANITY CHECK -- are distributions of the two cohorts similar?
pdf(paste0(QC_loc,"/",Today,".aems450k.cohorts.",EWAS_trait,".plaque.EffectSizeDistributions.pdf"),
    width = 10, height = 10, onefile = TRUE)
  par(mfrow = c(1,1), oma = c(0, 0, 2, 0))
  plot(aems450k1.plothisth, col = aems450k1.col, 
       ylim = c(0,80000),
       main = "Model corrected for age + sex + hospital", 
       xlab = bquote("effect size ["~beta~"]"))  # first histogram
  plot(aems450k2.plothisth, col = aems450k2.col, add = TRUE)  # second
  abline(v = mean.betah, col = "#595A5C", lty = 2, lwd = 1)
  abline(v = neg.cutoffh, col = "#595A5C", lty = 2, lwd = 1)
  abline(v = pos.cutoffh, col = "#595A5C", lty = 2, lwd = 1)
  
  legend("topright",c("AEMS450K1", "AEMS450K2"), 
         col = c(aems450k1.col, aems450k2.col), 
         pch = 16, bty = "n")
  mtext("Distribution of effect estimates", outer = TRUE, cex = 1.5)
  dev.off()
par(mfrow = c(1,1), oma = c(0, 0, 0, 0))

# actual bacon correction - run in parallel
library(BiocParallel)
register(MulticoreParam(8))
# plaque *with* hospital
aems450k.meta.bcp.th <- bacon(teststatistics = aems450k.meta.th, 
                              effectsizes = NULL, standarderrors = NULL,
                              verbose = TRUE) 
aems450k.meta.bcp.bh <- bacon(teststatistics = NULL,
                              effectsizes = aems450k.meta.esh, standarderrors = aems450k.meta.seh,
                              verbose = TRUE)

cat("   - Bacon model result.")
aems450k.meta.bcp.th
aems450k.meta.bcp.bh

cat("\n - estimates of the mixture fit")
estimates(aems450k.meta.bcp.th)
estimates(aems450k.meta.bcp.bh)

cat("\n - the inflation")
inflation(aems450k.meta.bcp.th)
inflation(aems450k.meta.bcp.bh)

cat("\n - the bias")
bias(aems450k.meta.bcp.th) 
bias(aems450k.meta.bcp.bh) 

# PLOT FOR SANITY CHECK -- is the meta-analysis behaving as expected?
# plaque *with* hospital
pdf(paste0(QC_loc,"/",Today,".aems450k.cohorts.traces.",EWAS_trait,".plaque.pdf"),
    width = 10, height = 10, onefile = TRUE)
  traces(aems450k.meta.bcp.th)
dev.off()
pdf(paste0(QC_loc,"/",Today,".aems450k.cohorts.posteriors.",EWAS_trait,".plaque.pdf"),
    width = 10, height = 10, onefile = TRUE)
  posteriors(aems450k.meta.bcp.th)
dev.off()
pdf(paste0(QC_loc,"/",Today,".aems450k.cohorts.fittedbacon.",EWAS_trait,".plaque.pdf"),
    width = 10, height = 10, onefile = TRUE)
  fit(aems450k.meta.bcp.th, n = 100) # visualization of the fitting using the Gibbs Sampling algorithm
dev.off()
pdf(paste0(QC_loc,"/",Today,".aems450k.cohorts.histogram.",EWAS_trait,".plaque.pdf"),
    width = 10, height = 10, onefile = TRUE)
  print(plot(aems450k.meta.bcp.th, type = "hist"))
dev.off()
pdf(paste0(QC_loc,"/",Today,".aems450k.cohorts.qq.",EWAS_trait,".plaque.pdf"),
    width = 10, height = 10, onefile = TRUE)
  print(plot(aems450k.meta.bcp.th, type = "qq"))
dev.off()

cat("\n * Perform meta-analysis.")
# actual meta-analysis in plaque *with* hospital
aems450k.meta.bcpm.bh <- meta(aems450k.meta.bcp.bh)
aems450k.meta.pvalsbcpm.bh <- pval(aems450k.meta.bcpm.bh)
head(aems450k.meta.pvalsbcpm.bh[order(aems450k.meta.pvalsbcpm.bh[,1]),])
print(bacon::topTable(aems450k.meta.bcpm.bh))

pdf(paste0(QC_loc,"/",Today,".aems450k.meta.qq.",EWAS_trait,".plaque.pdf"),
    width = 12, height = 12, onefile = TRUE)
  print(plot(aems450k.meta.bcpm.bh, type = "qq"))
dev.off()
png(paste0(QC_loc,"/",Today,".aems450k.meta.qq.",EWAS_trait,".plaque.png"),
    width = 800, height = 800)
  print(plot(aems450k.meta.bcpm.bh, type = "qq"))
dev.off()

cat("\n * Uncorrected results.")
# plaque *with* hospital
aems450k.meta.pvalsph <- pval(aems450k.meta.bcpm.bh, corrected = FALSE)
head(aems450k.meta.pvalsph[order(aems450k.meta.pvalsph[,3]),], 10)
aems450k.meta.zph = qnorm(aems450k.meta.pvalsph[,3]/2)
aems450k.meta.lambdaph = round(median(aems450k.meta.zph^2)/0.4549364,3)
cat(paste0("\n  - lambda is: [",aems450k.meta.lambdaph,"].")) # 1.2713076331

cat("\n * Corrected results.")
# plaque *with* hospital
aems450k.meta.pvalsbcph <- pval(aems450k.meta.bcpm.bh, corrected = TRUE)
head(aems450k.meta.pvalsbcph[order(aems450k.meta.pvalsbcph[,3]),], 10)
aems450k.meta.zbcph = qnorm(aems450k.meta.pvalsbcph[,3]/2)
aems450k.meta.lambdabcph = round(median(aems450k.meta.zbcph^2)/0.4549364,3)
cat(paste0("\n  - lambda is: [",aems450k.meta.lambdabcph,"].")) # 1.2713076331

cat("\n * Annotating corrected results for plotting and other purposes.")
cat("\n   - annotating...")
# plaque
# for some reason does this function take forever since the last update to R 3.4.3.
# infobcp <- cpgInfo(rownames(aems450k.meta.pvalsbcph), TxDb = "TxDb.Hsapiens.UCSC.hg19.knownGene")
infobcp <- as.data.frame(fread(paste0(OUT_loc, "/20171229.aems450k.meta.infobcp.hm450kannot.txt"),
                 header = TRUE, na.strings = "NA", dec = ".", verbose = TRUE, showProgress = TRUE), keep.rownames = TRUE)
rownames(infobcp) <- infobcp$V1
infobcp$V1 <- NULL
# head(infobcp)

cat("\n   - merging with results...")
# add in effect sizes
temp.results <- merge(infobcp, es(aems450k.meta.bcpm.bh), by = "row.names")
rownames(temp.results) <- temp.results$Row.names
temp.results2 <- temp.results[, -which(names(temp.results) %in% c("Row.names"))]
names(temp.results2)[names(temp.results2) == "limmaB.aems450k1"] <- "EffectSize.aems450k1"
names(temp.results2)[names(temp.results2) == "limmaB.aems450k2"] <- "EffectSize.aems450k2"
names(temp.results2)[names(temp.results2) == "meta"] <- "EffectSize.meta"
# add in standard errors
temp.results3 <- merge(temp.results2, se(aems450k.meta.bcpm.bh), by = "row.names")
rownames(temp.results3) <- temp.results3$Row.names
temp.results4 <- temp.results3[, -which(names(temp.results3) %in% c("Row.names"))]
names(temp.results4)[names(temp.results4) == "limmaSE.aems450k1"] <- "SE.aems450k1"
names(temp.results4)[names(temp.results4) == "limmaSE.aems450k2"] <- "SE.aems450k2"
names(temp.results4)[names(temp.results4) == "meta"] <- "SE.meta"
# add in uncorrected p-values
temp.results5 <- merge(temp.results4, aems450k.meta.pvalsph, by = "row.names")
rownames(temp.results5) <- temp.results5$Row.names
temp.results6 <- temp.results5[, -which(names(temp.results5) %in% c("Row.names"))]
names(temp.results6)[names(temp.results6) == "limmaB.aems450k1"] <- "Pval.aems450k1"
names(temp.results6)[names(temp.results6) == "limmaB.aems450k2"] <- "Pval.aems450k2"
temp.results7 <- temp.results6[, -which(names(temp.results6) %in% c("meta"))]
# add in corrected p-values
aems450k.meta.resultsp <- merge(temp.results7, aems450k.meta.pvalsbcph, by = "row.names")
rownames(aems450k.meta.resultsp) <- aems450k.meta.resultsp$Row.names
names(aems450k.meta.resultsp)[names(aems450k.meta.resultsp) == "Row.names"] <- "CpG"
names(aems450k.meta.resultsp)[names(aems450k.meta.resultsp) == "limmaB.aems450k1"] <- "PvalCor.aems450k1"
names(aems450k.meta.resultsp)[names(aems450k.meta.resultsp) == "limmaB.aems450k2"] <- "PvalCor.aems450k2"
names(aems450k.meta.resultsp)[names(aems450k.meta.resultsp) == "meta"] <- "PvalCor.meta"
aems450k.meta.resultsp$Pval.meta <- pval(aems450k.meta.bcpm.bh, corrected = FALSE)[,3]
aems450k.meta.resultsp$PvalCorAdj.meta <- p.adjust(aems450k.meta.resultsp$PvalCor.meta, method = "BH")
head(aems450k.meta.resultsp)
str(aems450k.meta.resultsp)
dim(aems450k.meta.resultsp)

cat("\n   - removing intermediate file...")
rm(temp.results, temp.results2, temp.results3, temp.results4, temp.results5, temp.results6, temp.results7)

cat("\n   - recreating the chromosomes...")
list.chr <- levels(as.factor(aems450k.meta.resultsp$seqnames))[1:25]
# list.chr
aems450k.meta.resultspf <- aems450k.meta.resultsp[(aems450k.meta.resultsp$seqnames %in% list.chr), ]
aems450k.meta.resultspf$chr <- gsub("chr", "",aems450k.meta.resultspf$seqnames)
aems450k.meta.resultspf$chr <- gsub("X", "23",aems450k.meta.resultspf$chr)
aems450k.meta.resultspf$chr <- gsub("Y", "24",aems450k.meta.resultspf$chr)
aems450k.meta.resultspf$chr <- gsub("XY", "25",aems450k.meta.resultspf$chr)
aems450k.meta.resultspf$chr <- gsub("M", "26",aems450k.meta.resultspf$chr)
aems450k.meta.resultspf$chr <- as.numeric(aems450k.meta.resultspf$chr)

cat("\n * Annotating results with additional CpG information.")
cat("\n   - annotating...")
# plaque
# install.packages.auto("IlluminaHumanMethylation450kmanifest")
# install.packages.auto("IlluminaHumanMethylation450kanno.ilmn12.hg19")
# ?IlluminaHumanMethylation450kmanifest
# data(IlluminaHumanMethylation450kmanifest)
# IlluminaHumanMethylation450kmanifest
# ?IlluminaHumanMethylation450kanno.ilmn12.hg19
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
IlluminaHumanMethylation450kanno.ilmn12.hg19
# data(Locations)
# data(Manifest)
# data(SNPs.137CommonSingle)
data(Islands.UCSC)
# data(Other)
# Locations
# Manifest
# SNPs.137CommonSingle
# Islands.UCSC
# Other

Islands.UCSC.forannotate <- as.data.frame(Islands.UCSC)
Islands.UCSC.forannotate$CpG <- rownames(Islands.UCSC.forannotate)
Other.forannotate <- as.data.frame(Other)

aems450k.meta.resultspfCGI.temp <- merge(aems450k.meta.resultspf, Islands.UCSC.forannotate, by = "row.names")
aems450k.meta.resultspfCGI.temp$Row.names <- NULL
aems450k.meta.resultspfCGI.temp$CpG.y <- NULL
rownames(aems450k.meta.resultspfCGI.temp) <- aems450k.meta.resultspfCGI.temp$CpG.x
names(aems450k.meta.resultspfCGI.temp)[names(aems450k.meta.resultspfCGI.temp) == "CpG.x"] <- "CpG"
aems450k.meta.resultspfCGI <- merge(aems450k.meta.resultspfCGI.temp, Other.forannotate, by = "row.names")
aems450k.meta.resultspfCGI$Row.names <- NULL
rownames(aems450k.meta.resultspfCGI) <- aems450k.meta.resultspfCGI$CpG
head(aems450k.meta.resultspfCGI)
str(aems450k.meta.resultspfCGI)
dim(aems450k.meta.resultspfCGI)
rm(aems450k.meta.resultspf, Islands.UCSC.forannotate, Other.forannotate, Islands.UCSC, aems450k.meta.resultspfCGI.temp)

cat("\n* Plotting...")
cat("\n  - QQ-plot in plaque...")
png(paste0(PLOT_loc,"/",Today,".aems450k.meta.",EWAS_trait,".plaque.QQPlot.png"),
    width = 1024, height = 800)
pdf(paste0(PLOT_loc,"/",Today,".aems450k.meta.",EWAS_trait,".plaque.QQPlot.pdf"),
    width = 12, height = 10, onefile = TRUE)
postscript(paste0(PLOT_loc,"/",Today,".aems450k.meta.",EWAS_trait,".plaque.QQPlot.ps"),
           width = 12, height = 10, onefile = TRUE, bg = "transparent", family = "Helvetica")
tiff(paste0(PLOT_loc,"/",Today,".aems450k.meta.",EWAS_trait,".plaque.QQPlot.tiff"),
     width = 1024, height = 800, units = "px", pointsize = 12,
     compression = "none", bg = "transparent",
     type = "quartz")

  par(mfrow = c(1,2), oma = c(0, 0, 2, 0), mar = c(5, 6, 4, 2) + 0.1)
  qq(aems450k.meta.resultspfCGI$Pval.meta, 
     main = bquote("Uncorrected " ~ lambda == .(aems450k.meta.lambdaph)),
     col = "#1290D9", pch = 16, xlim = c(0,8), ylim = c(0,25),
     cex = 1.75, cex.lab = 1.75, cex.axis = 1.75,  cex.main = 1.50,
     bty = "n")
  abline(0, 1, col = "#E55738")
  qq(aems450k.meta.resultspfCGI$PvalCor.meta, 
     main = bquote("Corrected" ~ lambda == .(aems450k.meta.lambdabcph)),
     col = "#1290D9", pch = 16, xlim = c(0,8), ylim = c(0,25),
     cex = 1.75, cex.lab = 1.75, cex.axis = 1.75, cex.main = 1.50,
     bty = "n")
  abline(0, 1, col = "#E55738")
  mtext("QQ-plots", outer = TRUE, cex = 2.0)
dev.off()
par(mfrow = c(1,1), oma = c(0, 0, 0, 0), mar = c(5, 4, 4, 2) + 0.1)

cat("\n  - Manhattan-plot in plaque...")
# plaque
# meta-analyse
# head(aems450k.meta.resultspfCGI[order(aems450k.meta.resultspfCGI[,31]),], 10)
# example of known smoking-related CpG
# (other) top hits: AHRR (cg25648203, cg05575921, cg03991871, cg12806681, cg02385153
#     resp. p = 5.823467e-19, p = 8.715543e-15, p = 1.101234e-10, p = 1.162238e-07, p = 1.705505e-07, chr5), 
# NTHL1 (cg16650073, cg01815912, p = 7.853257e-10, p = 1.103098e-06, chr16), 
# ALPI (cg05951221, p = 1.172635e-07, chr2)
# ITPK1 (cg05284742, p = 3.044496e-07, chr14), 
# CRLF1 (cg22702618, p = 1.103098e-06, chr1)
# EEFSEC (cg19505196, p = 5.787613e-07, chr3)
png(paste0(PLOT_loc,"/",Today,".aems450k.meta.",EWAS_trait,".plaque.ManhattanPlot.png"),
    width = 1920, height = 1080)
pdf(paste0(PLOT_loc,"/",Today,".aems450k.meta.",EWAS_trait,".plaque.ManhattanPlot.pdf"),
    width = 28, height = 8, onefile = TRUE)
postscript(paste0(PLOT_loc,"/",Today,".aems450k.meta.",EWAS_trait,".plaque.ManhattanPlot.ps"),
           width = 28, height = 10, onefile = TRUE, bg = "transparent", family = "Helvetica")
tiff(paste0(PLOT_loc,"/",Today,".aems450k.meta.",EWAS_trait,".plaque.ManhattanPlot.tiff"),
     width = 1920, height = 1080, units = "px", pointsize = 12,
     compression = "none", bg = "transparent",
     type = "quartz")

  par(mfrow = c(1,1), oma = c(0, 0, 0, 0), mai = c(1, 1, 1, 0))

  ahrr <- aems450k.meta.resultspfCGI$CpG[grep("AHRR", aems450k.meta.resultspfCGI$SYMBOL)]
  itpk1 <- aems450k.meta.resultspfCGI$CpG[grep("ITPK1", aems450k.meta.resultspfCGI$SYMBOL)]
  hilight.top <- c(ahrr, itpk1)
  manhattan.uithof(aems450k.meta.resultspfCGI, chr = "chr", bp = "start", p = "PvalCor.meta", snp = "CpG", 
                   highlight = hilight.top, suggestiveline = FALSE,
                   genomewideline = FALSE, 
                   cex = 2.5, cex.axis = 1.75, cex.lab = 1.75,
                   col = uithof_color, ylim = c(0,20),
                   chrlabs = c(1:22, "X", "Y"),
                   annotatePval = NULL, annotateTop = NULL)
dev.off()
rm(hilight.top, ahrr, itpk1)
par(mfrow = c(1,1), oma = c(0, 0, 0, 0), mai = c(1, 1, 1, 1))

cat("\n* Filter on *REPLICATED* signal only...")
aems450k.meta.resultspfCGIQC <- subset(aems450k.meta.resultspfCGI, PvalCor.aems450k2 < 0.05)
dim(aems450k.meta.resultspfCGIQC)

aems450k.meta.zph.QC = qnorm(aems450k.meta.resultspfCGIQC$Pval.meta/2)
aems450k.meta.lambdaph.QC = round(median(aems450k.meta.zph.QC^2)/0.4549364,3) # 1.293
aems450k.meta.zph.cor.QC = qnorm(aems450k.meta.resultspfCGIQC$PvalCor.meta/2)
aems450k.meta.lambdaph.cor.QC = round(median(aems450k.meta.zph.cor.QC^2)/0.4549364,3) # 3.273

cat("\n* Plotting...")
cat("\n  - QQ-plot in plaque...")
png(paste0(PLOT_loc,"/",Today,".aems450k.meta.",EWAS_trait,".plaque.QQPlot.replicated.png"),
    width = 1024, height = 800)
pdf(paste0(PLOT_loc,"/",Today,".aems450k.meta.",EWAS_trait,".plaque.QQPlot.replicated.pdf"),
    width = 12, height = 10, onefile = TRUE)
postscript(paste0(PLOT_loc,"/",Today,".aems450k.meta.",EWAS_trait,".plaque.QQPlot.replicated.ps"),
           width = 12, height = 10, onefile = TRUE, bg = "transparent", family = "Helvetica")
tiff(paste0(PLOT_loc,"/",Today,".aems450k.meta.",EWAS_trait,".plaque.QQPlot.replicated.tiff"),
     width = 1024, height = 800, units = "px", pointsize = 12,
     compression = "none", bg = "transparent",
     type = "quartz")

  par(mfrow = c(1,2), oma = c(0, 0, 2, 0), mar = c(5, 6, 4, 2) + 0.1)
  qq(aems450k.meta.resultspfCGIQC$Pval.meta, 
     main = bquote("Uncorrected " ~ lambda == .(aems450k.meta.lambdaph.QC)),
     # main = bquote("Uncorrected "),
     col = "#1290D9", pch = 16, xlim = c(0,8), ylim = c(0,25),
     cex = 1.75, cex.lab = 1.75, cex.axis = 1.75, cex.main = 1.50,
     bty = "n")
  abline(0, 1, col = "#E55738")
  qq(aems450k.meta.resultspfCGIQC$PvalCor.meta, 
     main = bquote("Corrected" ~ lambda == .(aems450k.meta.lambdaph.cor.QC)),
     # main = bquote("Corrected"),
     col = "#1290D9", pch = 16, xlim = c(0,8), ylim = c(0,25),
     cex = 1.75, cex.lab = 1.75, cex.axis = 1.75, cex.main = 1.50,
     bty = "n")
  abline(0, 1, col = "#E55738")
  mtext("QQ-plots", outer = TRUE, cex = 2.0)
dev.off()
par(mfrow = c(1,1), oma = c(0, 0, 0, 0), mar = c(5, 4, 4, 2) + 0.1)

cat("\n  - Manhattan-plot in plaque...")
# plaque
# head(aems450k.meta.resultspfCGIQC[order(aems450k.meta.resultspfCGIQC[,31]),], 10)
# example of known smoking-related CpG
# (other) top hits: AHRR (cg25648203, cg05575921, cg03991871, cg12806681, cg02385153
#     resp. p = 5.823467e-19, p = 8.715543e-15, p = 1.101234e-10, p = 1.162238e-07, p = 1.705505e-07, chr5), 
# NTHL1 (cg16650073, cg01815912, p = 7.853257e-10, p = 1.103098e-06, chr16), 
# ALPI (cg05951221, p = 1.172635e-07, chr2)
# ITPK1 (cg05284742, p = 3.044496e-07, chr14), 
# CRLF1 (cg22702618, p = 1.103098e-06, chr1)
# EEFSEC (cg19505196, p = 5.787613e-07, chr3)

png(paste0(PLOT_loc,"/",Today,".aems450k.meta.",EWAS_trait,".plaque.ManhattanPlot.replicated.png"),
    width = 1920, height = 1080)
pdf(paste0(PLOT_loc,"/",Today,".aems450k.meta.",EWAS_trait,".plaque.ManhattanPlot.replicated.pdf"),
    width = 28, height = 8, onefile = TRUE)
postscript(paste0(PLOT_loc,"/",Today,".aems450k.meta.",EWAS_trait,".plaque.ManhattanPlot.replicated.ps"),
           width = 28, height = 10, onefile = TRUE, bg = "transparent", family = "Helvetica")
tiff(paste0(PLOT_loc,"/",Today,".aems450k.meta.",EWAS_trait,".plaque.ManhattanPlot.replicated.tiff"),
     width = 1920, height = 1080, units = "px", pointsize = 12,
     compression = "none", bg = "transparent",
     type = "quartz")

  par(mfrow = c(1,1), oma = c(0, 0, 0, 0), mai = c(1, 1, 1, 0))
  
  ahrr <- aems450k.meta.resultspfCGIQC$CpG[grep("AHRR", aems450k.meta.resultspfCGIQC$SYMBOL)]
  itpk1 <- aems450k.meta.resultspfCGIQC$CpG[grep("ITPK1", aems450k.meta.resultspfCGIQC$SYMBOL)]
  hilight.top <- c(ahrr, itpk1)
  manhattan.uithof(aems450k.meta.resultspfCGIQC, chr = "chr", bp = "start", p = "PvalCor.meta", snp = "CpG", 
                   highlight = hilight.top, suggestiveline = FALSE,
                   genomewideline = FALSE, 
                   cex = 2.50, cex.axis = 1.75, cex.lab = 1.75,
                   col = uithof_color, ylim = c(0,20),
                   chrlabs = c(1:22, "X", "Y"),
                   annotatePval = NULL, annotateTop = NULL)
dev.off()
rm(hilight.top, ahrr, itpk1)
par(mfrow = c(1,1), oma = c(0, 0, 0, 0), mai = c(1, 1, 1, 1))

cat("\n  > in plaque...")
pdf(paste0(PLOT_loc,"/",Today,".aems450k.meta.",EWAS_trait,".plaque.TopCor.pdf"))
  top.symbol <- aems450k.meta.resultspfCGIQC[order(aems450k.meta.resultspfCGIQC$PvalCor.meta), "SYMBOL"][1]
  top <- aems450k.meta.resultspfCGIQC[order(aems450k.meta.resultspfCGIQC$PvalCor.meta), "CpG"][1]
  aems450k1.x <- aems450k1.dataph[rownames(aems450k1.dataph) == top, ]
  aems450k1.y <- aems450k1.designph[, 2]
  aems450k1.pearsonR <- signif(cor(aems450k1.x, aems450k1.y), 3)
  
  aems450k2.x <- aems450k2.dataph[rownames(aems450k2.dataph) == top, ]
  aems450k2.y <- aems450k2.designph[, 2]
  aems450k2.pearsonR <- signif(cor(aems450k2.x, aems450k2.y), 3)
  
  par(mfrow = c(1,2), oma = c(0, 0, 2, 0), mar = c(5, 5, 4, 0) + 0.1)
  # AEMS450K1
  boxplot(aems450k1.x~aems450k1.y, main = bquote(Pearson~r^2 == .(aems450k1.pearsonR)),
          sub = "(AEMS450K1)",
          cex.sub = 0.8,
          xlab = "Current smoker", 
          ylab = bquote(.(top)~" -"~italic(.(top.symbol))), 
          plot = TRUE, notch = FALSE, outline = TRUE, 
          names = c("no", "yes"), 
          whisklty = 1, staplelty = 0, # no whisker-ends, but with lines
          frame = FALSE)
  stripchart(aems450k1.x~aems450k1.y, vertical = TRUE,
             method = "jitter", add = TRUE, pch = 20, col = "#1290D9")
  # AEMS450K2
  boxplot(aems450k2.x~aems450k2.y, main = bquote(Pearson~r^2 == .(aems450k2.pearsonR)),
          sub = "(AEMS450K2)",
          cex.sub = 0.8,
          xlab = "Current smoker", 
          # ylab = bquote(.(top)~" methylation near "~italic(.(top.symbol))), 
          ylab = "",
          plot = TRUE, notch = FALSE, outline = TRUE, 
          names = c("no", "yes"), 
          whisklty = 1, staplelty = 0, # no whisker-ends, but with lines
          frame = FALSE)
  stripchart(aems450k2.x~aems450k2.y, vertical = TRUE,
             method = "jitter", add = TRUE, pch = 20, col = "#1290D9")
  mtext("Top correlated CpG", outer = TRUE, cex = 1.25)

dev.off()
rm(top.data, top.symbol, top,
   aems450k1.x, aems450k1.y, aems450k1.pearsonR,
   aems450k2.x, aems450k2.y, aems450k2.pearsonR)
par(mfrow = c(1,1), oma = c(0, 0, 0, 0), mar = c(5, 4, 4, 2) + 0.1)

pdf(paste0(PLOT_loc,"/",Today,".aems450k.meta.",EWAS_trait,".plaque.Top4.pdf"), paper = "a4r",
    width = 12, height = 8)
  # head(aems450k.meta.resultspfCGIQC[order(aems450k.meta.resultspfCGIQC[,31]),], 4)
  # example of known smoking-related CpG
  # No. CpG       seqnames SYMBOL PvalCor.meta    PvalCorAdj.meta
  # 1   cg05575921   chr5   AHRR  1.713638e-13    3.796429e-08
  # 2   cg03991871   chr5   AHRR  2.896744e-11    4.278337e-06
  # 3   cg12806681   chr5   AHRR  5.949100e-09    5.271902e-04
  # 4   cg05284742  chr14  ITPK1  2.045460e-07    1.510517e-02

  top.ahrr1 <- aems450k.meta.resultspfCGIQC[order(aems450k.meta.resultspfCGIQC$PvalCor.meta), "CpG"][1]
  top.ahrr2 <- aems450k.meta.resultspfCGIQC[order(aems450k.meta.resultspfCGIQC$PvalCor.meta), "CpG"][2]
  top.ahrr3 <- aems450k.meta.resultspfCGIQC[order(aems450k.meta.resultspfCGIQC$PvalCor.meta), "CpG"][3]
  top.itpk1 <- aems450k.meta.resultspfCGIQC[order(aems450k.meta.resultspfCGIQC$PvalCor.meta), "CpG"][4]
  
  par(mfrow = c(2,2), mar = c(5,6,4,0), oma = c(0, 0, 2, 0))
  # AHRR 1 - 3
  x.ahrr1 <- aems450k1.dataph[rownames(aems450k1.dataph) == top.ahrr1, ]
  y.ahrr1 <- aems450k1.designph[, 2]
  pearsonR.ahrr1 <- signif(cor(x.ahrr1, y.ahrr1), 3)
  x.ahrr2 <- aems450k1.dataph[rownames(aems450k1.dataph) == top.ahrr2, ]
  y.ahrr2 <- aems450k1.designph[, 2]
  pearsonR.ahrr2 <- signif(cor(x.ahrr2, y.ahrr2), 3)
  x.ahrr3 <- aems450k1.dataph[rownames(aems450k1.dataph) == top.ahrr3, ]
  y.ahrr3 <- aems450k1.designph[, 2]
  pearsonR.ahrr3 <- signif(cor(x.ahrr3, y.ahrr3), 3)
  
  boxplot(x.ahrr1~y.ahrr1, main = bquote(Pearson~r^2 == .(pearsonR.ahrr1)),
          xlab = "Current Smoking", 
          ylab = bquote(.(top.ahrr1)~" -"~italic("AHRR")), 
          # ylim = c(0, 4.5),
          plot = TRUE, notch = FALSE, outline = TRUE, 
          cex = 1.25, cex.lab = 1.50, cex.axis = 1.50, cex.main = 2.0,
          names = c("no", "yes"), 
          whisklty = 1, staplelty = 0, # no whisker-ends, but with lines
          frame = FALSE)
  stripchart(x.ahrr1~y.ahrr1, vertical = TRUE,
             method = "jitter", add = TRUE, pch = 20, col = "#1290D9")
  boxplot(x.ahrr2~y.ahrr2, main = bquote(Pearson~r^2 == .(pearsonR.ahrr2)),
          xlab = "Current Smoking", 
          ylab = bquote(.(top.ahrr2)~" -"~italic("AHRR")),
          # ylim = c(0, 4.5),
          plot = TRUE, notch = FALSE, outline = TRUE, 
          cex = 1.25, cex.lab = 1.50, cex.axis = 1.50, cex.main = 2.0,
          names = c("no", "yes"), 
          whisklty = 1, staplelty = 0, # no whisker-ends, but with lines
          frame = FALSE)
  stripchart(x.ahrr2~y.ahrr2, vertical = TRUE,
             method = "jitter", add = TRUE, pch = 20, col = "#1290D9")
  boxplot(x.ahrr3~y.ahrr3, main = bquote(Pearson~r^2 == .(pearsonR.ahrr3)),
          xlab = "Current Smoking", 
          ylab = bquote(.(top.ahrr3)~" -"~italic("AHRR")),
          # ylim = c(0, 4.5),
          plot = TRUE, notch = FALSE, outline = TRUE, 
          cex = 1.25, cex.lab = 1.50, cex.axis = 1.50, cex.main = 2.0,
          names = c("no", "yes"), 
          whisklty = 1, staplelty = 0, # no whisker-ends, but with lines
          frame = FALSE)
  stripchart(x.ahrr3~y.ahrr3, vertical = TRUE,
             method = "jitter", add = TRUE, pch = 20, col = "#1290D9")
  
  # ITPK1
  x.itpk1 <- aems450k1.dataph[rownames(aems450k1.dataph) == top.itpk1, ]
  y.itpk1 <- aems450k1.designph[, 2]
  pearsonR.itpk1 <- signif(cor(x.itpk1, y.itpk1), 3)
  boxplot(x.itpk1~y.itpk1, main = bquote(Pearson~r^2 == .(pearsonR.itpk1)),
          xlab = "Current Smoking", 
          ylab = bquote(.(top.itpk1)~" -"~italic("ITPK1")),
          # ylim = c(0, 4.5),
          plot = TRUE, notch = FALSE, outline = TRUE, 
          cex = 1.25, cex.lab = 1.50, cex.axis = 1.50, cex.main = 2.0,
          names = c("no", "yes"), 
          whisklty = 1, staplelty = 0, # no whisker-ends, but with lines
          frame = FALSE)
  stripchart(x.itpk1~y.itpk1, vertical = TRUE,
             method = "jitter", add = TRUE, pch = 20, col = "#1290D9")
  
  mtext("Top 4 CpGs", outer = TRUE, cex = 2.25)

dev.off()
rm(top.ahrr1,top.ahrr2,top.ahrr3,top.itpk1,
   x.ahrr1, y.ahrr1, pearsonR.ahrr1,
   x.ahrr2, y.ahrr2, pearsonR.ahrr2,
   x.ahrr3, y.ahrr3, pearsonR.ahrr3,
   x.itpk1, y.itpk1, pearsonR.itpk1)
par(mfrow = c(1,1), mar = c(5, 4, 4, 2), oma = c(0, 0, 0, 0))

cat("\n* Saving results...")
cat("\n  - writing ALL results...")
fwrite(aems450k.meta.resultspfCGI, 
       file = paste0(OUT_loc, "/", Today,".aems450k.meta.",EWAS_trait,".ResultsPlaqueCleaned.txt"),
       quote = FALSE, sep = ";", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE,
       showProgress = TRUE, verbose = TRUE)
cat("\n  - writing top 10 results...")
fwrite(utils::head(aems450k.meta.resultspfCGI[order(aems450k.meta.resultspfCGI[,31]),], 10), 
       file = paste0(OUT_loc, "/", Today,".aems450k.meta.",EWAS_trait,".resultspf.top10meta.txt"),
       quote = FALSE, sep = ";", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE,
       showProgress = TRUE, verbose = TRUE)
cat("\n  - writing top 20 discovery results...")
fwrite(utils::head(aems450k.meta.resultspfCGI[order(aems450k.meta.resultspfCGI[,27]),], 20), 
       file = paste0(OUT_loc, "/", Today,".aems450k.meta.",EWAS_trait,".resultspf.top20discovery.txt"),
       quote = FALSE, sep = ";", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE,
       showProgress = TRUE, verbose = TRUE)
cat("\n  - writing top 20 replication results...")
fwrite(utils::head(aems450k.meta.resultspfCGIQC[order(aems450k.meta.resultspfCGIQC[,31]),], 20), 
       file = paste0(OUT_loc, "/", Today,".aems450k.meta.",EWAS_trait,".resultspf.top20replication.txt"),
       quote = FALSE, sep = ";", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE,
       showProgress = TRUE, verbose = TRUE)

cat("\n*** saving final datasets data in between steps ***")
aems450k.meta.resultspfCGI.SmokerCurrent = aems450k.meta.resultspfCGI
aems450k.meta.resultspfCGIQC.SmokerCurrent = aems450k.meta.resultspfCGIQC
save(aems450k.meta.resultspfCGI.SmokerCurrent, file = paste0(OUT_loc,"/",Today,".",EWAS_trait,".aems450k.meta.resultspfCGI.RData"))
save(aems450k.meta.resultspfCGIQC.SmokerCurrent, file = paste0(OUT_loc,"/",Today,".",EWAS_trait,".aems450k.meta.resultspfCGIQC.RData"))

cat("\n* Let's clean up some old objects we do not need anymore...")
rm(list.chr, aems450k.meta.resultsp, aems450k.meta.pvalsph, aems450k.meta.zph, aems450k.meta.pvalsbcph, aems450k.meta.zbcph,
   aems450k.meta.bcp.bh, aems450k.meta.bcpm.bh, aems450k.meta.bcp.th, aems450k.meta.esh, aems450k.meta.seh,
   aems450k1.col, aems450k2.col, aems450k1.plothisth, aems450k2.plothisth,
   mean.betah, sd.betah, neg.cutoffh, pos.cutoffh,
   # aems450k1.dataph, aems450k2.dataph,
   # aems450k1.designph, aems450k2.designph,
   aems450k1.fitph, aems450k2.fitph,
   aems450k1.effectsizeph, aems450k2.effectsizeph,
   aems450k1.SEph, aems450k2.SEph,
   aems450k1.tstatph, aems450k2.tstatph,
   aems450k1.pvalph, aems450k2.pvalph,
   aems450k1.padjph, aems450k2.padjph, 
   aems450k1.meta.difs, aems450k2.meta.difs, aems450k.meta.intersect, 
   aems450k1.ranges, aems450k2.ranges, 
   hm450.manifest.pop.GoNL, 
   aems450k1.covariates, aems450k2.covariates,
   aems450k1.nas, aems450k2.nas,
   feats, chr.list, regions)

cat("\n===========================================================================================")
cat("SAVE THE DATA")

save.image(paste0(ANALYSIS_loc,"/",Today,".aems450k.meta.analysis.ewas.",EWAS_trait,".plaque.RData"))


