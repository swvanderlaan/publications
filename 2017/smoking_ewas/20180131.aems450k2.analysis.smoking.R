cat("===========================================================================================
                ANALYSIS ATHERO-EXPRESS METHYLATION STUDIES 450K 2
    
    Version:      v2.9
    
    Last update:  2018-01-31
    Written by:   Sander W. van der Laan (s.w.vanderlaan-2@umcutrecht.nl);
    Marten A. Siemelink
    
    Description:  Script to do meta-analysis of Athero-Express Methylation 
    Study 450K 2 (2016) EWAS on smoking.
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
install.packages.auto("limma")
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
install_github("molepi/DNAmArray", force = F)
library(DNAmArray)
install_github("molepi/omicsPrint", ref = "R3.4", force = F)
library(omicsPrint)
install_github("bbmri-nl/BBMRIomics", subdir = "BBMRIomics", force = F)
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
PROJECTDATASET = "AEMS450K2"
PROJECTNAME = "smoking"
EWAS_trait = "SmokerCurrent" # Phenotype

INP_AE_loc = paste0(ROOT_loc, "/PLINK/_AE_Originals")
INP_AEMS450K1_loc = paste0(INP_AE_loc, "/AEMS450K1")
INP_AEMS450K2_loc = paste0(INP_AE_loc, "/AEMS450K2")

### Mac
EPIGENETICS_loc = paste0(ROOT_loc, "/PLINK/analyses/epigenetics")

# ### HPC
# EPIGENETICS_loc = paste0(ROOT_loc, "/svanderlaan/projects/epigenetics")

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

cat("----------------------------------------------------------------------------")
cat("\nLOAD ATHERO-EXPRESS METHYLATION STUDY 2 DATASETS")
setwd(INP_loc)
list.files()

cat(paste0("\n\n* Load ",PROJECTDATASET," data..."))
cat("\n  - loading B/Mvalues of plaque samples...")
load(paste0(INP_AEMS450K2_loc,"/20171229.aems450k2.BvaluesQCIMP.plaque.RData"))
load(paste0(INP_AEMS450K2_loc,"/20171229.aems450k2.MvaluesQCIMP.plaque.RData"))

cat("----------------------------------------------------------------------------")
cat(paste0("\n[ EPIGENOME-WIDE ASSOCIATION STUDIES on SMOKING in PLAQUE in ",PROJECTDATASET," ]"))
# Reference: https://molepi.github.io/DNAmArray_workflow/06_EWAS.html

cat("\n* Sannity checking the data.")
# png(paste0(PLOT_loc,"/",Today,".aems450k2.smoking.plaque.MethylationDensity.png"),
#     width = 800, height = 600)
pdf(paste0(QC_loc,"/",Today,".aems450k2.smoking.plaque.MethylationDensity.pdf"),
    width = 12, height = 8, onefile = TRUE)
  par(mfrow = c(1,2), oma = c(0, 0, 2, 0))
  densityPlot(assays(aems450k2.BvaluesQCplaque)$data, sampGroups = aems450k2.BvaluesQCplaque$SmokerCurrent, main = "Beta-values", 
              legend = FALSE, 
              xlab = "Beta-values", 
              col = c("#9FC228", "#E55738"), 
              bty = "n")
  legend("topright", legend = levels(factor(aems450k2.BvaluesQCplaque$SmokerCurrent)), 
         text.col = c("#9FC228", "#E55738"), 
         lty = 1, lwd = 1, col = c("#9FC228", "#E55738"), 
         bty = "n")
  densityPlot(assays(aems450k2.MvaluesQCplaque)$data, sampGroups = aems450k2.MvaluesQCplaque$SmokerCurrent, main = "M-values", 
              legend = FALSE, 
              xlab = "M-values", 
              col = c("#9FC228", "#E55738"), 
              bty = "n")
  legend("topright", legend = levels(factor(aems450k2.MvaluesQCplaque$SmokerCurrent)), 
         text.col = c("#9FC228", "#E55738"), 
         lty = 1, lwd = 1, col = c("#9FC228", "#E55738"), 
         bty = "n")
  mtext(paste0("Overall methylation density (plaque, ",PROJECTDATASET,")"), outer = TRUE, cex = 1.5)
  par(mfrow = c(1,1), oma = c(0,0,0,0))
dev.off()

cat("\nRemoving BvaluesQCIMP objects - as we don't use these anymore.")
rm(aems450k2.BvaluesQCplaque)

cat("===========================================================================================")
cat("\n[ CONTINUE EPIGENOME-WIDE ASSOCIATION STUDIES on SMOKING in PLAQUE in ",PROJECTDATASET," ]")

cat("\n* Setup the analysis.")
require(FDb.InfiniumMethylation.hg19)
feats <- features(FDb.InfiniumMethylation.hg19)
chr.list <- levels(seqnames(feats))
regions <- feats[seqnames(feats) %in% chr.list]

cat("\n  - Setup the model, first covariate is the variable/phenotype of interest.")
#plaque
MvaluesQCplaqueClean <- aems450k2.MvaluesQCplaque
metadata(MvaluesQCplaqueClean)$formula <- ~SmokerCurrent + Sample_Sex + Age + Hospital

cat("\n  - Next, we extract those samples having a complete set of covariates. 
    Notice that we subset the SummarizedExperiment-object!")
# plaque
covariatesp <- get_all_vars(metadata(MvaluesQCplaqueClean)$formula, data = colData(MvaluesQCplaqueClean))
nasp <- apply(covariatesp, 1, anyNA)
MvaluesQCplaqueClean <- MvaluesQCplaqueClean[, !nasp]

cat("\n  - Remove probes on the X and Y chromosomes and probes containing 
    SNPs or does mapping to multiple locations (ref: W. Zhou, P.W. Laird, and H. Shen. 2016. Nucleic Acids Res.).")
cat("\n > first chromosome X and Y... [NOTE: we skip this as we are interested in X, Y, and M(t)]")
# plaque
# MvaluesQCplaqueClean <- MvaluesQCplaqueClean[!(seqnames(MvaluesQCplaqueClean) %in% c("chrX","chrY")),]
cat("\n > next probes with SNPs...")
data(hm450.manifest.pop.GoNL) ##From DNAmArray
#hm450.manifest.pop.GoNL
hm450.manifest.pop.GoNL <- hm450.manifest.pop.GoNL[!is.na(hm450.manifest.pop.GoNL$MASK.general.GoNL) &
                                                     hm450.manifest.pop.GoNL$MASK.general.GoNL == TRUE, ]
# plaque
MvaluesQCplaqueClean <- MvaluesQCplaqueClean[!(names(MvaluesQCplaqueClean) %in% names(hm450.manifest.pop.GoNL)),]  
cat("\n > we started with:")
dim(MvaluesQCplaque) # 484249 CpGs, 190 samples
cat("\n > we end up with:")
dim(MvaluesQCplaqueClean) # 443614 CpGs, 187 samples

cat("\n* Run the analysis.")
install.packages.auto("limma")
# plaque
designp <- model.matrix(metadata(MvaluesQCplaqueClean)$formula, data = colData(MvaluesQCplaqueClean))
datap <- assays(MvaluesQCplaqueClean)$data
fitp <- lmFit(datap, designp)
cat("  - calculate T-statistic = beta/standard error...")
tstatp <- fitp$coef/fitp$stdev.unscaled/fitp$sigma
effectsizep <- fitp$coef[,2]
cat("  - standard error is stdev.unscaled/sigma in limma prior to eBayes()...")
SEp <- fitp$stdev.unscaled[,2]/fitp$sigma
cat("  - calculate p-value from T-statistic, and adjust using FDR (Benjamin-Hochberg)...")
pvalp <- 2*pnorm(-abs(tstatp[,2]))
padjp <- p.adjust(sort(pvalp, decreasing = FALSE), method = "BH")
head(padjp[padjp < 0.05])

cat("\n* Correct for bias and inflation using 'BACON'. 
    (ref: M. van Iterson, E.W. van Zwet, and B.T. Heijmans. 2017. Genome Biol. 18 (1): 19.).")
library("bacon")
# plaque
# tstatsp <- cbind(limma = tstatp[,2], cate = as.vector(fap$beta.t))
tstatsp <- cbind(limmaT = tstatp[,2], limmaB = effectsizep, limmaSE = SEp)

library(BiocParallel)
register(MulticoreParam(8))
# plaque
bcp <- bacon(teststatistics = tstatsp[,1], effectsizes = tstatsp[,2], standarderrors = tstatsp[,3],
             verbose = TRUE) 
cat("\n* Bacon model result.")
bcp
cat("\n - estimates of the mixture fit")
estimates(bcp) 
cat("\n - the inflation")
inflation(bcp)  
cat("\n - the bias")
bias(bcp) 

# PLOT FOR SANITY CHECK -- is the meta-analysis behaving as expected?
# plaque *with* hospital
pdf(paste0(QC_loc,"/",Today,".aems450k2.traces.smoking.plaque.pdf"),
    width = 10, height = 10, onefile = TRUE)
  traces(bcp)
dev.off()
pdf(paste0(QC_loc,"/",Today,".aems450k2.posteriors.smoking.plaque.pdf"),
    width = 10, height = 10, onefile = TRUE)
  posteriors(bcp)
dev.off()
pdf(paste0(QC_loc,"/",Today,".aems450k2.fittedbacon.smoking.plaque.pdf"),
    width = 10, height = 10, onefile = TRUE)
  fit(bcp, n = 100) # visualization of the fitting using the Gibbs Sampling algorithm
dev.off()
pdf(paste0(QC_loc,"/",Today,".aems450k2.histogram.smoking.plaque.pdf"),
    width = 10, height = 10, onefile = TRUE)
  print(plot(bcp, type = "hist"))
dev.off()
pdf(paste0(QC_loc,"/",Today,".aems450k2.qq.smoking.plaque.pdf"),
    width = 10, height = 10, onefile = TRUE)
  print(plot(bcp, type = "qq"))
dev.off()

cat("\n * Uncorrected results.")
cat("\n - in plaque...")
pvalsp <- pval(bcp, corrected = FALSE)
head(pvalsp[order(pvalsp[,1]),])
zp = qnorm(pvalsp[,1]/2)
lambdap = round(median(zp^2)/0.4549364,3)
cat(paste0("\n  - lambda is: [",lambdap,"].")) # lambda: 1.433

cat("\n * Corrected results.")
cat("\n - in plaque...")
pvalsbcp <- pval(bcp, corrected = TRUE)
head(pvalsbcp[order(pvalsbcp[,1]),])
zbcp = qnorm(pvalsbcp[,1]/2)
lambdabcp = round(median(zbcp^2)/0.4549364,3)
cat(paste0("\n  - lambda is: [",lambdabcp,"].")) # lambda: 1.307

cat("\n * Annotating corrected results for plotting and other purposes.")
cat("\n   - annotating...")
# plaque
# for some reason does this function take forever since the last update to R 3.4.3.
# infobcp <- cpgInfo(rownames(pvalsbcp), TxDb = "TxDb.Hsapiens.UCSC.hg19.knownGene")
infobcp <- as.data.frame(fread(paste0(OUT_loc, "/20170828.aems450k2.infobcp.hm450kannot.txt"),
                               header = TRUE, na.strings = "NA", dec = ".", verbose = TRUE, showProgress = TRUE), keep.rownames = TRUE)
rownames(infobcp) <- infobcp$V1
infobcp$V1 <- NULL
# head(infobcp)

# plaque
cat("\n   - merging with results...")
cat("\n   > in plaque...")
temp.results <- merge(infobcp, tstatsp, by = "row.names")
rownames(temp.results) <- temp.results$Row.names
temp.results2 <- temp.results[ , -which(names(temp.results) %in% c("Row.names"))]
temp.results3 <- merge(temp.results2, pvalsp, by = "row.names")
names(temp.results3)[names(temp.results3) == "V1"] <- "Pval_limma"
resultsp <- merge(temp.results3, pvalsbcp, by.x = "Row.names", by.y = "row.names")
names(resultsp)[names(resultsp) == "V1"] <- "Pval_limma_bacon"
names(resultsp)[names(resultsp) == "Row.names"] <- "CpG"
rownames(resultsp) <- resultsp$CpG
head(resultsp)
str(resultsp)
dim(resultsp)

cat("\n   - removing intermediate file...")
rm(temp.results, temp.results2, temp.results3)

cat("\n   - recreating the chromosomes...")
cat("\n   > in plaque...")
list.chr <- levels(as.factor(resultsp$seqnames))[1:25]
list.chr
resultspf <- resultsp[(resultsp$seqnames %in% list.chr), ]
resultspf$chr <- gsub("chr", "",resultspf$seqnames)
resultspf$chr <- gsub("X", "23",resultspf$chr)
resultspf$chr <- gsub("Y", "24",resultspf$chr)
resultspf$chr <- gsub("XY", "25",resultspf$chr)
resultspf$chr <- gsub("M", "26",resultspf$chr)
resultspf$chr <- as.numeric(resultspf$chr)

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
data(Locations)
# data(Manifest)
data(SNPs.137CommonSingle)
data(Islands.UCSC)
data(Other)
Locations
# Manifest
SNPs.137CommonSingle
Islands.UCSC
Other

Islands.UCSC.forannotate <- as.data.frame(Islands.UCSC)
Locations.forannotate <- as.data.frame(Locations)
SNPs.137CommonSingle.forannotate <- as.data.frame(SNPs.137CommonSingle)
Other.forannotate <- as.data.frame(Other)
Islands.UCSC.forannotate$CpG <- rownames(Islands.UCSC.forannotate)
Locations.forannotate$CpG <- rownames(Locations.forannotate)
SNPs.137CommonSingle.forannotate$CpG <- rownames(SNPs.137CommonSingle.forannotate)
Other.forannotate$CpG <- rownames(Other.forannotate)

# PLAQUE
aems450k2.resultspfCGI.temp <- merge(resultspf, Islands.UCSC.forannotate, by = "row.names")
aems450k2.resultspfCGI.temp$Row.names <- NULL
aems450k2.resultspfCGI.temp$CpG.y <- NULL
rownames(aems450k2.resultspfCGI.temp) <- aems450k2.resultspfCGI.temp$CpG.x
names(aems450k2.resultspfCGI.temp)[names(aems450k2.resultspfCGI.temp) == "CpG.x"] <- "CpG"
aems450k2.resultspfCGI <- merge(aems450k2.resultspfCGI.temp, Other.forannotate, by = "row.names")
aems450k2.resultspfCGI$Row.names <- NULL
rownames(aems450k2.resultspfCGI) <- aems450k2.resultspfCGI$CpG.x
aems450k2.resultspfCGI$CpG.y <- NULL
names(aems450k2.resultspfCGI)[names(aems450k2.resultspfCGI) == "CpG.x"] <- "CpG"
utils::head(aems450k2.resultspfCGI)
str(aems450k2.resultspfCGI)
dim(aems450k2.resultspfCGI)
rm(resultspf, aems450k2.resultspfCGI.temp, Other, SNPs.137CommonSingle, Locations, Islands.UCSC)

cat("\n* Plotting...")
cat("\n  - QQ-plot in plaque...")
png(paste0(PLOT_loc,"/",Today,".aems450k2.smoking.plaque.QQPlot.png"),
    width = 1024, height = 800)
# pdf(paste0(PLOT_loc,"/",Today,".aems450k2.smoking.plaque.QQPlot.pdf"),
#     width = 12, height = 10, onefile = TRUE)
# postscript(paste0(PLOT_loc,"/",Today,".aems450k2.smoking.plaque.QQPlot.ps"),
#            width = 12, height = 10, onefile = TRUE, bg = "transparent", family = "Helvetica")
  par(mfrow = c(1,2), oma = c(0, 0, 2, 0), mar = c(5, 6, 4, 2) + 0.1)
  qq(aems450k2.resultspfCGI$Pval_limma, 
     main = bquote("Uncorrected " ~ lambda == .(lambdap)),
     col = "#1290D9", pch = 16, xlim = c(0,8), ylim = c(0,25),
     cex = 1.75, cex.lab = 1.75, cex.axis = 1.75,  cex.main = 1.75,
     bty = "n")
  abline(0, 1, col = "#E55738")
  qq(aems450k2.resultspfCGI$Pval_limma_bacon, 
     main = bquote("Corrected" ~ lambda == .(lambdabcp)),
     col = "#1290D9", pch = 16, xlim = c(0,8), ylim = c(0,25),
     cex = 1.75, cex.lab = 1.75, cex.axis = 1.75, cex.main = 1.75,
     bty = "n")
  abline(0, 1, col = "#E55738")
  mtext("QQ-plots", outer = TRUE, cex = 2.0)
dev.off()
par(mfrow = c(1,1), oma = c(0, 0, 0, 0), mar = c(5, 4, 4, 2) + 0.1)

cat("\n  - Manhattan-plot in plaque...")
  # 0.05/nrow(resultspf)*10
  # head(resultspf[order(resultspf[,23]),], 10)
  # example of known smoking-related CpG
  # (other) top hits: AHRR (cg05575921, p = 1.492649e-06, chr5), 
  # CHSY1 (cg24312730, p = 2.522552e-06, chr15), 
  # MTUS2 (cg24366680, p = 3.663666e-06, chr13),
  # MRPL54 (cg16260977, p = 5.718725e-06, chr19),
  # SF3B1 (cg08357651, p = 5.942498e-06, chr2), 
  # SNX1 (cg03431705, p = 6.228995e-06, chr15), 
  # ahrr <- resultspf$CpG[grep("AHRR", resultspf$SYMBOL)]
  # chsy1 <- resultspf$CpG[grep("CHSY1", resultspf$SYMBOL)]
  # mtus2 <- resultspf$CpG[grep("MTUS2", resultspf$SYMBOL)]
  # mrpl54 <- resultspf$CpG[grep("MRPL54", resultspf$SYMBOL)]
  # sf3b1 <- resultspf$CpG[grep("SF3B1", resultspf$SYMBOL)]
  # snx1 <- resultspf$CpG[grep("SNX1", resultspf$SYMBOL)]

    # hilight.top <- c(ahrr, chsy1, mtus2, mrpl54, sf3b1, snx1)
  # hilight.top <- c(ahrr, chsy1, mtus2)
png(paste0(PLOT_loc,"/",Today,".aems450k2.smoking.plaque.ManhattanPlot.png"),
    width = 1920, height = 1080)
# pdf(paste0(PLOT_loc,"/",Today,".aems450k2.smoking.plaque.ManhattanPlot.pdf"),
#     width = 28, height = 8, onefile = TRUE)
# postscript(paste0(PLOT_loc,"/",Today,".aems450k2.smoking.plaque.ManhattanPlot.ps"),
#            width = 28, height = 10, onefile = TRUE, bg = "transparent", family = "Helvetica")

  par(mfrow = c(1,1), oma = c(0, 0, 0, 0), mai = c(1, 1, 1, 0))
  
  # Loci identified in AEMS450K1
  ahrr <- aems450k2.resultspfCGI$CpG[grep("AHRR", aems450k2.resultspfCGI$SYMBOL)]
  nthl1 <- aems450k2.resultspfCGI$CpG[grep("NTHL1", aems450k2.resultspfCGI$SYMBOL)]
  crlf1 <- aems450k2.resultspfCGI$CpG[grep("CRLF1", aems450k2.resultspfCGI$SYMBOL)]
  glis1 <- aems450k2.resultspfCGI$CpG[grep("GLIS1", aems450k2.resultspfCGI$SYMBOL)]
  eefsec <- aems450k2.resultspfCGI$CpG[grep("EEFSEC", aems450k2.resultspfCGI$SYMBOL)]
  prdm16 <- aems450k2.resultspfCGI$CpG[grep("PRDM16", aems450k2.resultspfCGI$SYMBOL)]
  olfm2 <- aems450k2.resultspfCGI$CpG[grep("OLFM2", aems450k2.resultspfCGI$SYMBOL)]

  hilight.top <- c(ahrr, nthl1, crlf1, glis1, 
                   eefsec, prdm16, olfm2)
  manhattan.uithof(aems450k2.resultspfCGI, chr = "chr", bp = "start", p = "Pval_limma_bacon", snp = "CpG", 
                   highlight = hilight.top, suggestiveline = FALSE,
                   genomewideline = FALSE, 
                   cex = 2.5, cex.axis = 1.75, cex.lab = 1.75,
                   col = uithof_color, ylim = c(0,20),
                   chrlabs = c(1:22, "X", "Y"),
                   annotatePval = NULL, annotateTop = NULL)
dev.off()
rm(hilight.top, ahrr, nthl1, crlf1, glis1, eefsec, prdm16, olfm2)
par(mfrow = c(1,1), oma = c(0, 0, 0, 0), mai = c(1, 1, 1, 1))


cat("\n  - Top result correlation plot...")
pdf(paste0(PLOT_loc,"/",Today,".aems450k2.smoking.plaque.TopCor.pdf"))
  top.data <- aems450k2.resultspfCGI[order(aems450k2.resultspfCGI$Pval_limma_bacon), ][1,]
  top.symbol <- aems450k2.resultspfCGI[order(aems450k2.resultspfCGI$Pval_limma_bacon), "SYMBOL"][1]
  top <- aems450k2.resultspfCGI[order(aems450k2.resultspfCGI$Pval_limma_bacon), "CpG"][1]
  x <- datap[rownames(datap) == top, ]
  y <- designp[, 2]
  pearsonR <- signif(cor(x, y), 3)
  
  par(mfrow = c(1,1), mar = c(5,4,4,2))
  boxplot(x~y, main = bquote(Pearson~r^2 == .(pearsonR)),
          xlab = "Current Smoking", 
          ylab = bquote(.(top)~" methylation near "~italic(.(top.symbol))),
          plot = TRUE, notch = FALSE, outline = TRUE, 
          names = c("no smoker", "current smoker"), 
          whisklty = 1, staplelty = 0, # no whisker-ends, but with lines
          frame = FALSE)
  stripchart(x~y, vertical = TRUE,
             method = "jitter", add = TRUE, pch = 20, cex = 1.50, col = "#1290D9")
dev.off()
rm(top.data, top.symbol, top,
   x, y, pearsonR)

cat("\n* Saving results...")
cat("\n  - writing ALL results...")
fwrite(aems450k2.resultspfCGI, 
       file = paste0(OUT_loc, "/", Today,".aems450k2.ResultsPlaqueCleaned.txt"),
       quote = FALSE, sep = ";", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE,
       showProgress = TRUE, verbose = TRUE)

cat("\n* Let's clean up some old objects we do not need anymore...")
rm(resultsp, bcp, covariatesp, 
   nasp, fitp, tstatp, tstatsp, 
   feats, hm450.manifest.pop.GoNL, regions, chr.list, list.chr,
   SEp,
   pvalsbcp, pvalsp, zp, zbcp, 
   # datap, designp, # we keep this for future use
   infobcp, 
   padjp, pvalp, effectsizep)

cat("\n===========================================================================================")
cat("SAVE THE DATA")

save.image(paste0(ANALYSIS_loc,"/",Today,".aems450k2.analysis.ewas.smoking.plaque.RData"))

