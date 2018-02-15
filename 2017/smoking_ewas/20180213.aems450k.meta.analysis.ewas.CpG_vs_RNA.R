cat("===========================================================================================
                META-ANALYSIS ATHERO-EXPRESS METHYLATION STUDIES 450K 1 & 2
                        -- CpG vs. RNA expression: A Pilot Study --
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
ifelse(!dir.exists(file.path(ANALYSIS_loc, "/RNAseq")), 
       dir.create(file.path(ANALYSIS_loc, "/RNAseq")), 
       FALSE)
RNAOUT_loc = paste0(ANALYSIS_loc, "/RNAseq")
ifelse(!dir.exists(file.path(ANALYSIS_loc, "/CPGRNA")), 
       dir.create(file.path(ANALYSIS_loc, "/CPGRNA")), 
       FALSE)
CPGRNA_loc = paste0(ANALYSIS_loc, "/CPGRNA")

cat("===========================================================================================")
cat("\nLOAD ATHERO-EXPRESS METHYLATION STUDY DATASETS")
setwd(INP_loc)
list.files()

cat(paste0("\n\n* Load ",PROJECTDATASET," data...\n"))
cat("\n  - loading Mvalues of plaque samples...\n")
load(paste0(OUT_loc,"/20180207.aems450k1.MvaluesQCplaqueClean.TFvar.RData"))
aems450k1.MvaluesQCplaqueClean = aems450k1.MvaluesQCplaqueCleanTFvar

cat("\n  - loading EWAS analyses *REPLICATED* results...\n")
load(paste0(OUT_loc,"/20180207.SmokerCurrent.aems450k.meta.resultspfCGIQC.RData"))

aems450k.meta.resultspfCGIQC = aems450k.meta.resultspfCGIQC.SmokerCurrent

cat("\n* loading AE database...\n")
load(paste0(OUT_loc,"/20180207.AEData_2016.update.DF.CEA.TFvar.RData"))

AEData_2016.update.DF.CEA.SampleSize = nrow(AEData_2016.update.DF.CEA.TFvar)
AEData_2016.update.DF.CEA = AEData_2016.update.DF.CEA.TFvar

cat("\n* loading analyzed RNA pilot data...\n")
# load(paste0(RNAOUT_loc,"/20180208.aerna.pilot.",EWAS_trait,".dds_qcf.RData"))
# 
# load(paste0(RNAOUT_loc,"/20180208.aerna.pilot.",EWAS_trait,".se_fqc.RData"))
# 
# load(paste0(RNAOUT_loc,"/20180208.aerna.pilot.",EWAS_trait,".res_qcf.RData"))
# load(paste0(RNAOUT_loc,"/20180208.aerna.pilot.",EWAS_trait,".res_qcfdata.RData"))

load(paste0(RNAOUT_loc,"/20180208.aerna.pilot.",EWAS_trait,".sampleTable.RData"))
load(paste0(RNAOUT_loc,"/20180208.aerna.pilot.",EWAS_trait,".sampleTableFiltered.RData"))

cat("\n* raw, mapped data normalized for the read depth...")
cat("\nLoad datasets...")
INPRNA_loc = paste0(ROOT_loc, "/PLINK/_AE_Originals/AERNASMOKE")
RNAPROJECTNAME = "160223_NS500813_0095_AHT2FKBGXX"
DATADIR = paste0(INPRNA_loc, "/",RNAPROJECTNAME)
NORMDATA = read.table(paste0(DATADIR, "/read_counts/",RNAPROJECTNAME,"_readCounts_normalized.txt"), header = TRUE, 
                      sep = "\t", row.names = 1)

cat("===========================================================================================")
cat("CpGs vs RNAseq")

cat("\nCorrelation of significant CpGs to genes.")
### dds_qcf - normalized data, exclusion of PC outliers, filtered genes with zero-counts
### AHRR = ENSG00000063438; CpG = cg05575921
### Also: all of the 21 RNAseq samples overlap with AEMS450K1.

cat("\n* PART 1: compare all significant *and* replicated CpGs with the expression of all genes (annotated to those CpGs)...")
# - get all the result + RNAseq counts in one table
# - get all the result + methylation M-values in one table
# - filter away CpGs not mapped to a gene/transcript
# - get a list of all CpGs associated with DEmethylation, PvalCor.meta <= 0.05, beta < 0
# - get a list of all CpGs associated with methylation, PvalCor.meta <= 0.05, beta > 0
# - list all genes per list of CpGs
# - in RNAseq data filter out all genes not in above not
# - correlate DEmethylation and methylation with RNAseq expression
# - plot this correlation
# - nice-to-have group CpG lists per CpG location

# Remember: our methylation data and RNA data are in SummarizedExperiments 
# Ref: http://www.bioconductor.org/packages/3.7/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html#subsetting
cat("  > extract only relevant samples.\n")
rna.sample.list <- sampleTableFiltered$STUDYNUMBER
rna.sample.name.list <- sampleTableFiltered$RunID

aems450k1.MvaluesQCplaqueClean.21samples <- aems450k1.MvaluesQCplaqueClean[,aems450k1.MvaluesQCplaqueClean$STUDY_NUMBER %in% rna.sample.list]
# dim(aems450k1.MvaluesQCplaqueClean.21samples)
# aems450k1.MvaluesQCplaqueClean.21samples
rm(aems450k1.MvaluesQCplaqueClean, aems450k1.MvaluesQCplaqueCleanTFvar)

cat("  > extract only relevant CpGs, i.e. replicated.\n")
cpg.repl.list <- aems450k.meta.resultspfCGIQC$CpG

aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl <- aems450k1.MvaluesQCplaqueClean.21samples[cpg.repl.list,]
# dim(aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl)
# aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl

cat("  > we set the colnames to be the same in the AEMS450K1 data...")
colnames(aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl) <- aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl$STUDY_NUMBER
# aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl
# rowRanges(aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl)
genome(aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl) <- "hg19"
aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl <- sortSeqlevels(aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl)
aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl <- sort(aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl)
# rowRanges(aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl)
# metadata(aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl)
# aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl

cat("  > concatenate all the annotations and parse them (externally)...")
# Map location information to CpGs -- very nice to have
# source("http://bioconductor.org/biocLite.R")
# biocLite("IlluminaHumanMethylation450kanno.ilmn12.hg19")
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data("Locations")
data("Other")
data("Manifest")
data("SNPs.147CommonSingle")
data("SNPs.Illumina")
data("Islands.UCSC")
anno.450k <- IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other
# names(anno.450k)

anno.450k.Locations <- as.data.frame(Locations)
anno.450k.Locations$CpG <- rownames(anno.450k.Locations)
anno.450k.Manifest <- as.data.frame(Manifest)
anno.450k.Manifest$CpG <- rownames(anno.450k.Manifest)
anno.450k.Other <- as.data.frame(Other)
anno.450k.Other$CpG <- rownames(anno.450k.Other)
anno.450k.SNPs.147CommonSingle <- as.data.frame(SNPs.147CommonSingle)
anno.450k.SNPs.147CommonSingle$CpG <- rownames(anno.450k.SNPs.147CommonSingle)
anno.450k.SNPs.Illumina <- as.data.frame(SNPs.Illumina)
anno.450k.SNPs.Illumina$CpG <- rownames(anno.450k.SNPs.Illumina)
anno.450k.Islands.UCSC <- as.data.frame(Islands.UCSC)
anno.450k.Islands.UCSC$CpG <- rownames(anno.450k.Islands.UCSC)

anno.450k.temp1 <- merge.data.frame(anno.450k.Locations, anno.450k.Manifest, by = "CpG")
anno.450k.temp2 <- merge.data.frame(anno.450k.temp1, anno.450k.Other, by = "CpG")
anno.450k.temp3 <- merge.data.frame(anno.450k.temp2, anno.450k.SNPs.147CommonSingle, by = "CpG")
anno.450k.temp4 <- merge.data.frame(anno.450k.temp3, anno.450k.SNPs.Illumina, by = "CpG")
anno.450k.combined <- merge.data.frame(anno.450k.temp4, anno.450k.Islands.UCSC, by = "CpG")

rm(anno.450k.temp1, anno.450k.temp2, anno.450k.temp3, anno.450k.temp4)

fwrite(anno.450k.combined, file = paste0(CPGRNA_loc, "/anno.450k.combined.txt"), 
       sep = "\t", col.names = TRUE, row.names = FALSE, na = "NA",
       showProgress = TRUE, verbose = TRUE)

# Use a Perl script (ann2expann.pl) to split each line on GeneName/GeneID/GeneGroup by ";",
# and write out each line to a new file.
anno.450k.combined.EXPANDED <- fread(paste0(CPGRNA_loc, "/anno.450k.combined.EXPANDED.txt"), 
                                     na.strings = "NA", showProgress = TRUE, verbose = TRUE)
# dim(anno.450k.combined.EXPANDED)
# anno.450k.combined.EXPANDED[1:5,1:5]

cat("\nNow we will get lists of CpGs based on each Gene Group...\n")
# Now there are a couple of gene groups; we will split the annotation file based on these groups
# and plot the association of CpG with RNAseq based on these.
levels(as.factor(anno.450k.combined.EXPANDED$UCSC_RefGene_Group))
# Group 1: ""
# Group 2: "1stExon"
# Group 3: "3'UTR"
# Group 4: "5'UTR"
# Group 5: "Body"
# Group 6: "TSS1500"
# Group 7: "TSS200" 

temp <- subset(anno.450k.combined.EXPANDED, UCSC_RefGene_Group == "", select = c(CpG, UCSC_RefGene_Name, UCSC_RefGene_Group))
temp$CpG_Group <- "Intergenic"
anno.450k.combined.EXPANDED.rest.DF = temp
dim(anno.450k.combined.EXPANDED.rest.DF)
anno.450k.combined.EXPANDED.rest.DF[1:5, 1:4]
anno.450k.combined.EXPANDED.rest <- temp$CpG
temp <- subset(anno.450k.combined.EXPANDED, UCSC_RefGene_Group == "1stExon", select = c(CpG, UCSC_RefGene_Name, UCSC_RefGene_Group))
temp$CpG_Group <- "1st exon"
anno.450k.combined.EXPANDED.1stexon.DF = temp
anno.450k.combined.EXPANDED.1stexon <- temp$CpG
temp <- subset(anno.450k.combined.EXPANDED, UCSC_RefGene_Group == "3'UTR", select = c(CpG, UCSC_RefGene_Name, UCSC_RefGene_Group))
temp$CpG_Group <- "3'UTR"
anno.450k.combined.EXPANDED.3utr.DF = temp
anno.450k.combined.EXPANDED.3utr <- temp$CpG
temp <- subset(anno.450k.combined.EXPANDED, UCSC_RefGene_Group == "5'UTR", select = c(CpG, UCSC_RefGene_Name, UCSC_RefGene_Group))
temp$CpG_Group <- "5'UTR"
anno.450k.combined.EXPANDED.5utr.DF = temp
anno.450k.combined.EXPANDED.5utr <- temp$CpG
temp <- subset(anno.450k.combined.EXPANDED, UCSC_RefGene_Group == "Body", select = c(CpG, UCSC_RefGene_Name, UCSC_RefGene_Group))
temp$CpG_Group <- "Body"
anno.450k.combined.EXPANDED.body.DF = temp
anno.450k.combined.EXPANDED.body <- temp$CpG
temp <- subset(anno.450k.combined.EXPANDED, UCSC_RefGene_Group == "TSS200", select = c(CpG, UCSC_RefGene_Name, UCSC_RefGene_Group))
temp$CpG_Group <- "TSS200"
anno.450k.combined.EXPANDED.tss200.DF = temp
anno.450k.combined.EXPANDED.tss200 <- temp$CpG
temp <- subset(anno.450k.combined.EXPANDED, UCSC_RefGene_Group == "TSS1500", select = c(CpG, UCSC_RefGene_Name, UCSC_RefGene_Group))
temp$CpG_Group <- "TSS1500"
anno.450k.combined.EXPANDED.tss1500.DF = temp
anno.450k.combined.EXPANDED.tss1500 <- temp$CpG
rm(temp)

cat("\n...and we will subset the methylation data accordingly...\n")
cat("\nGet a list of probes that difference between the two datasets.\n")
ranges <- unlist(names(rowRanges(aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl)))

cat("\nWhat is in both?\n")
ranges.intersect.rest <- intersect(ranges, anno.450k.combined.EXPANDED.rest)
# length(ranges.intersect.rest) # 3688
ranges.intersect.1stexon <- intersect(ranges, anno.450k.combined.EXPANDED.1stexon)
# length(ranges.intersect.1stexon) # 1120
ranges.intersect.3utr <- intersect(ranges, anno.450k.combined.EXPANDED.3utr)
# length(ranges.intersect.3utr) # 646
ranges.intersect.5utr <- intersect(ranges, anno.450k.combined.EXPANDED.5utr)
# length(ranges.intersect.5utr) # 2064
ranges.intersect.body <- intersect(ranges, anno.450k.combined.EXPANDED.body)
# length(ranges.intersect.body) # 5911
ranges.intersect.tss200 <- intersect(ranges, anno.450k.combined.EXPANDED.tss200)
# length(ranges.intersect.tss200) # 1872
ranges.intersect.tss1500 <- intersect(ranges, anno.450k.combined.EXPANDED.tss1500)
# length(ranges.intersect.tss1500) # 2677

CpG.Type <- c("TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "3'UTR", "Intergenic")
CpG.N <- c(length(ranges.intersect.tss1500), length(ranges.intersect.tss200), length(ranges.intersect.5utr), 
           length(ranges.intersect.1stexon), length(ranges.intersect.body), length(ranges.intersect.3utr),
           length(ranges.intersect.rest))
CpG.TypeN.DF <- data.frame(CpG.Type, CpG.N)
rm(CpG.Type, CpG.N)

pdf(paste0(CPGRNA_loc,"/",Today,".aerna.pilot.signCpG_per_Type.pdf"))
# postscript(paste0(CPGRNA_loc,"/",Today,".aerna.pilot.signCpG_per_Type.ps"))
  # Simple Horizontal Bar Plot with Added Labels 
  par(las = 2) # make label text perpendicular to axis
  par(mar = c(5,8,4,2)) # increase y-axis margin.
  bp <- barplot(CpG.TypeN.DF$CpG.N, main = "Type of CpGs", horiz = TRUE,
                names.arg = c("TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "3'UTR", "Intergenic"),
                col = c("#FBB820", "#F59D10", "#C5D220", "#CC0071", "#1290D9", "#9FC228", "#595A5C"), border = FALSE,
                cex.names = 0.8)
  text(0, bp, format(CpG.TypeN.DF$CpG.N, big.mark = ","), cex = 1, pos = 4)
dev.off()
  par(las = 0) # make label text perpendicular to axis
  par(mar = c(5, 4, 4, 2) + 0.1) # increase y-axis margin.

cat("\n * subsetting the methylation data...\n")
aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl.rest <- aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl[(names(aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl) %in% ranges.intersect.rest),]  
aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl.1stexon <- aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl[(names(aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl) %in% ranges.intersect.1stexon),]  
aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl.3utr <- aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl[(names(aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl) %in% ranges.intersect.3utr),]  
aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl.5utr <- aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl[(names(aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl) %in% ranges.intersect.5utr),]  
aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl.body <- aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl[(names(aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl) %in% ranges.intersect.body),]  
aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl.tss200 <- aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl[(names(aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl) %in% ranges.intersect.tss200),]  
aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl.tss1500 <- aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl[(names(aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl) %in% ranges.intersect.tss1500),]  

cat("\n * subsetting...\n")
aems450k.meta.resultspfCGIQC.cleaned = aems450k.meta.resultspfCGIQC
aems450k.meta.resultspfCGIQC.cleaned$addressA <- NULL
aems450k.meta.resultspfCGIQC.cleaned$addressB <- NULL
aems450k.meta.resultspfCGIQC.cleaned$channel <- NULL
aems450k.meta.resultspfCGIQC.cleaned$platform <- NULL
aems450k.meta.resultspfCGIQC.cleaned$sourceSeq <- NULL
aems450k.meta.resultspfCGIQC.cleaned$probeType <- NULL
aems450k.meta.resultspfCGIQC.cleaned$probeStart <- NULL
aems450k.meta.resultspfCGIQC.cleaned$probeEnd <- NULL
aems450k.meta.resultspfCGIQC.cleaned$probeTarget <- NULL
aems450k.meta.resultspfCGIQC.cleaned$chr <- NULL
aems450k.meta.resultspfCGIQC.cleaned$SourceSeq <- NULL
aems450k.meta.resultspfCGIQC.cleaned$Forward_Sequence <- NULL
aems450k.meta.resultspfCGIQC.cleaned$Pval.aems450k1 <- NULL
aems450k.meta.resultspfCGIQC.cleaned$Pval.aems450k2 <- NULL
aems450k.meta.resultspfCGIQC.cleaned$Pval.meta <- NULL
aems450k.meta.resultspfCGIQC.cleaned$PvalCorAdj.meta <- NULL
aems450k.meta.resultspfCGIQC.cleaned$UCSC_RefGene_Name <- NULL
aems450k.meta.resultspfCGIQC.cleaned$UCSC_RefGene_Accession <- NULL
aems450k.meta.resultspfCGIQC.cleaned$UCSC_RefGene_Group <- NULL
# dim(aems450k.meta.resultspfCGIQC.cleaned)
# aems450k.meta.resultspfCGIQC.cleaned[1:5,1:29]

temp <- merge(subset(as.data.frame(aems450k.meta.resultspfCGIQC.cleaned), PvalCor.meta <= 0.05), 
              as.data.frame(assay(aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl.rest)), 
              by = "row.names", sort = FALSE)
temp$Row.names <- NULL
res_aems450k1.rest <- merge(temp, anno.450k.combined.EXPANDED.rest.DF, 
                            by = "CpG", sort = FALSE)
res_aems450k1.rest$Row.names <- NULL
temp <- merge(subset(as.data.frame(aems450k.meta.resultspfCGIQC.cleaned), PvalCor.meta <= 0.05), 
              as.data.frame(assay(aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl.1stexon)), 
              by = "row.names", sort = FALSE)
temp$Row.names <- NULL
res_aems450k1.1stexon <- merge(temp, anno.450k.combined.EXPANDED.1stexon.DF, 
                            by = "CpG", sort = FALSE)
res_aems450k1.1stexon$Row.names <- NULL
temp <- merge(subset(as.data.frame(aems450k.meta.resultspfCGIQC.cleaned), PvalCor.meta <= 0.05), 
              as.data.frame(assay(aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl.3utr)), 
              by = "row.names", sort = FALSE)
temp$Row.names <- NULL
res_aems450k1.3utr <- merge(temp, anno.450k.combined.EXPANDED.3utr.DF, 
                            by = "CpG", sort = FALSE)
res_aems450k1.3utr$Row.names <- NULL
temp <- merge(subset(as.data.frame(aems450k.meta.resultspfCGIQC.cleaned), PvalCor.meta <= 0.05), 
              as.data.frame(assay(aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl.5utr)), 
              by = "row.names", sort = FALSE)
temp$Row.names <- NULL
res_aems450k1.5utr <- merge(temp, anno.450k.combined.EXPANDED.5utr.DF, 
                               by = "CpG", sort = FALSE)
res_aems450k1.5utr$Row.names <- NULL
temp <- merge(subset(as.data.frame(aems450k.meta.resultspfCGIQC.cleaned), PvalCor.meta <= 0.05), 
              as.data.frame(assay(aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl.body)), 
              by = "row.names", sort = FALSE)
temp$Row.names <- NULL
res_aems450k1.body <- merge(temp, anno.450k.combined.EXPANDED.body.DF, 
                               by = "CpG", sort = FALSE)
res_aems450k1.body$Row.names <- NULL
temp <- merge(subset(as.data.frame(aems450k.meta.resultspfCGIQC.cleaned), PvalCor.meta <= 0.05), 
              as.data.frame(assay(aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl.tss200)), 
              by = "row.names", sort = FALSE)
temp$Row.names <- NULL
res_aems450k1.tss200 <- merge(temp, anno.450k.combined.EXPANDED.tss200.DF, 
                               by = "CpG", sort = FALSE)
res_aems450k1.tss200$Row.names <- NULL
temp <- merge(subset(as.data.frame(aems450k.meta.resultspfCGIQC.cleaned), PvalCor.meta <= 0.05), 
              as.data.frame(assay(aems450k1.MvaluesQCplaqueClean.21samplesCpGrepl.tss1500)), 
              by = "row.names", sort = FALSE)
temp$Row.names <- NULL
res_aems450k1.tss1500 <- merge(temp, anno.450k.combined.EXPANDED.tss1500.DF, 
                              by = "CpG", sort = FALSE)
res_aems450k1.tss1500$Row.names <- NULL
rm(temp)

# we are left with several dataframes like this:
# dim(res_aems450k1.rest)
# res_aems450k1.rest[1:5, 1:53]

# calculate the average Mvalue across all samples
res_aems450k1.rest$MedianCpGMvalue <- apply(subset(res_aems450k1.rest, 
                                                   select = c("696", "904", "936", "947", "1029",
                                                              "1140", "1267", "1343", "1373", "1376",
                                                              "1390", "1413", "1858", "1932", "1938", 
                                                              "2117", "2195", "2207", "2419", "2555", "2611")), 1, median, na.rm = TRUE)
res_aems450k1.1stexon$MedianCpGMvalue <- apply(subset(res_aems450k1.1stexon, 
                                                      select = c("696", "904", "936", "947", "1029",
                                                                 "1140", "1267", "1343", "1373", "1376",
                                                                 "1390", "1413", "1858", "1932", "1938", 
                                                                 "2117", "2195", "2207", "2419", "2555", "2611")), 1, median, na.rm = TRUE)
res_aems450k1.3utr$MedianCpGMvalue <- apply(subset(res_aems450k1.3utr, 
                                                   select = c("696", "904", "936", "947", "1029",
                                                              "1140", "1267", "1343", "1373", "1376",
                                                              "1390", "1413", "1858", "1932", "1938", 
                                                              "2117", "2195", "2207", "2419", "2555", "2611")), 1, median, na.rm = TRUE)
res_aems450k1.5utr$MedianCpGMvalue <- apply(subset(res_aems450k1.5utr, 
                                                   select = c("696", "904", "936", "947", "1029",
                                                              "1140", "1267", "1343", "1373", "1376",
                                                              "1390", "1413", "1858", "1932", "1938", 
                                                              "2117", "2195", "2207", "2419", "2555", "2611")), 1, median, na.rm = TRUE)
res_aems450k1.body$MedianCpGMvalue <- apply(subset(res_aems450k1.body, 
                                                   select = c("696", "904", "936", "947", "1029",
                                                              "1140", "1267", "1343", "1373", "1376",
                                                              "1390", "1413", "1858", "1932", "1938", 
                                                              "2117", "2195", "2207", "2419", "2555", "2611")), 1, median, na.rm = TRUE)
res_aems450k1.tss200$MedianCpGMvalue <- apply(subset(res_aems450k1.tss200, 
                                                     select = c("696", "904", "936", "947", "1029",
                                                                "1140", "1267", "1343", "1373", "1376",
                                                                "1390", "1413", "1858", "1932", "1938", 
                                                                "2117", "2195", "2207", "2419", "2555", "2611")), 1, median, na.rm = TRUE)
res_aems450k1.tss1500$MedianCpGMvalue <- apply(subset(res_aems450k1.tss1500, 
                                                      select = c("696", "904", "936", "947", "1029",
                                                                 "1140", "1267", "1343", "1373", "1376",
                                                                 "1390", "1413", "1858", "1932", "1938", 
                                                                 "2117", "2195", "2207", "2419", "2555", "2611")), 1, median, na.rm = TRUE)

cat("calculate the average counts across all samples.\n")
NORMDATA$AvgCount <- rowMeans(subset(NORMDATA, 
                                     select = c("nsm1267", "nsm1343", "nsm1373", "nsm1376", "nsm1858",
                                                "nsm1932", "nsm2117", "nsm2207", "nsm2419", "nsm2555",
                                                "nsm2611", "nsm696", "nsm904", "nsm947", "sm1029",
                                                "sm1140", "sm1390", "sm1413", "sm1938", "sm2195", "sm936")), 
                              na.rm = TRUE)
cat("\n  - annotating results...")
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
NORMDATA$ENSEMBLID <- row.names(NORMDATA)
NORMDATA$UCSC_RefGene_Name <- mapIds(org.Hs.eg.db,
                                     keys = NORMDATA$ENSEMBLID,
                                     column = "SYMBOL",
                                     keytype = "ENSEMBL",
                                     multiVals = "first")
# dim(NORMDATA)
# NORMDATA[1:5,1:33]

# get a dataset with only the genes and the average count
NORMDATA.CLEANED <- subset(NORMDATA, NORMDATA$AvgCount > 0, select = c("UCSC_RefGene_Name", "AvgCount"))
# dim(NORMDATA.CLEANED)
# NORMDATA.CLEANED[1:5,1:2]

# calculate the per-gene average of all the averaged M-values across samples
# Ref: https://stackoverflow.com/questions/21982987/mean-per-group-in-a-data-frame
library(plyr)

res_aems450k1.rest.perGene <- ddply(res_aems450k1.rest, .(SYMBOL), summarize,  MedianCpGMvaluePerGene = median(MedianCpGMvalue))
res_aems450k1.1stexon.perGene <- ddply(res_aems450k1.1stexon, .(UCSC_RefGene_Name), summarize,  MedianCpGMvaluePerGene = median(MedianCpGMvalue))
res_aems450k1.3utr.perGene <- ddply(res_aems450k1.3utr, .(UCSC_RefGene_Name), summarize,  MedianCpGMvaluePerGene = median(MedianCpGMvalue))
res_aems450k1.5utr.perGene <- ddply(res_aems450k1.5utr, .(UCSC_RefGene_Name), summarize,  MedianCpGMvaluePerGene = median(MedianCpGMvalue))
res_aems450k1.body.perGene <- ddply(res_aems450k1.body, .(UCSC_RefGene_Name), summarize,  MedianCpGMvaluePerGene = median(MedianCpGMvalue))
res_aems450k1.tss200.perGene <- ddply(res_aems450k1.tss200, .(UCSC_RefGene_Name), summarize,  MedianCpGMvaluePerGene = median(MedianCpGMvalue))
res_aems450k1.tss1500.perGene <- ddply(res_aems450k1.tss1500, .(UCSC_RefGene_Name), summarize,  MedianCpGMvaluePerGene = median(MedianCpGMvalue))

# concatenate the counts and the demethylation/methylations
res_rna_450k.rest <- merge(as.data.table(res_aems450k1.rest.perGene), 
                              as.data.table(NORMDATA.CLEANED),
                              by.x = "SYMBOL", by.y = "UCSC_RefGene_Name", sort = FALSE)
res_rna_450k.rest$CpG.Effect <- ifelse(res_rna_450k.rest$MedianCpGMvaluePerGene < 0,
                                          c("demethylated"), c("methylated"))
res_rna_450k.1stexon <- merge(as.data.table(res_aems450k1.1stexon.perGene), 
                              as.data.table(NORMDATA.CLEANED),
                              by.x = "UCSC_RefGene_Name", by.y = "UCSC_RefGene_Name", sort = FALSE)
res_rna_450k.1stexon$CpG.Effect <- ifelse(res_rna_450k.1stexon$MedianCpGMvaluePerGene < 0,
                                          c("demethylated"), c("methylated"))
res_rna_450k.3utr <- merge(as.data.table(res_aems450k1.3utr.perGene), 
                              as.data.table(NORMDATA.CLEANED),
                              by.x = "UCSC_RefGene_Name", by.y = "UCSC_RefGene_Name", sort = FALSE)
res_rna_450k.3utr$CpG.Effect <- ifelse(res_rna_450k.3utr$MedianCpGMvaluePerGene < 0,
                                          c("demethylated"), c("methylated"))
res_rna_450k.5utr <- merge(as.data.table(res_aems450k1.5utr.perGene), 
                              as.data.table(NORMDATA.CLEANED),
                              by.x = "UCSC_RefGene_Name", by.y = "UCSC_RefGene_Name", sort = FALSE)
res_rna_450k.5utr$CpG.Effect <- ifelse(res_rna_450k.5utr$MedianCpGMvaluePerGene < 0,
                                          c("demethylated"), c("methylated"))
res_rna_450k.body <- merge(as.data.table(res_aems450k1.body.perGene), 
                              as.data.table(NORMDATA.CLEANED),
                              by.x = "UCSC_RefGene_Name", by.y = "UCSC_RefGene_Name", sort = FALSE)
res_rna_450k.body$CpG.Effect <- ifelse(res_rna_450k.body$MedianCpGMvaluePerGene < 0,
                                          c("demethylated"), c("methylated"))
res_rna_450k.tss200 <- merge(as.data.table(res_aems450k1.tss200.perGene), 
                              as.data.table(NORMDATA.CLEANED),
                              by.x = "UCSC_RefGene_Name", by.y = "UCSC_RefGene_Name", sort = FALSE)
res_rna_450k.tss200$CpG.Effect <- ifelse(res_rna_450k.tss200$MedianCpGMvaluePerGene < 0,
                                          c("demethylated"), c("methylated"))
res_rna_450k.tss1500 <- merge(as.data.table(res_aems450k1.tss1500.perGene), 
                              as.data.table(NORMDATA.CLEANED),
                              by.x = "UCSC_RefGene_Name", by.y = "UCSC_RefGene_Name", sort = FALSE)
res_rna_450k.tss1500$CpG.Effect <- ifelse(res_rna_450k.tss1500$MedianCpGMvaluePerGene < 0,
                                          c("demethylated"), c("methylated"))

### START PLOTTING
pdf(paste0(CPGRNA_loc,"/",Today,".aerna.pilot.signCpG_vs_theirGenes.tss1500.pdf"))
# postscript(paste0(CPGRNA_loc,"/",Today,".aerna.pilot.signCpG_vs_theirGenes.tss1500.ps"))
# nice-to-have: stratify by type of CpG
par(mfrow = c(1,1), mar = c(6,6,4,0), oma = c(0, 0, 2, 0))
  
  ### TSS1500
  gene.counts.tss1500 <- res_rna_450k.tss1500$AvgCount
  gene.methylation.tss1500 <- as.factor(res_rna_450k.tss1500$MedianCpGMvaluePerGene)
  # levels(as.factor(res_rna_450k.tss1500$CpG.Effect))
  gene.methylation.tss1500 <- as.numeric(gene.methylation.tss1500)
  gene.cor.tss1500 <- signif(cor(gene.counts.tss1500, gene.methylation.tss1500, 
                                 use = "pairwise.complete.obs", method = "spearman"), 3)
  gene.cor.tss1500
  wilcox.test(log10(AvgCount) ~ as.factor(CpG.Effect), data = res_rna_450k.tss1500, 
                      alternative = "two.sided")
  pval.tss1500 = format(wilcox.test(log10(AvgCount) ~ as.factor(CpG.Effect), data = res_rna_450k.tss1500, 
                                    alternative = "two.sided")$`p.value`, digits = 3, scientific = FALSE)
  
  boxplot(log10(AvgCount) ~ CpG.Effect, res_rna_450k.tss1500,
          main = bquote(rho == .(gene.cor.tss1500)~"-"~italic("P") == .(pval.tss1500)),
          # sub = "1,500bp from TSS",
          cex.sub = 1.75, cex.axis = 1.75, cex.main = 2.0, cex.lab = 2,
          ylim = c(-2.5, 3.5),
          xlab = "CpG effect\n1,500bp from TSS",
          ylab = expression(log[10]~(normalized~counts)), 
          plot = TRUE, notch = FALSE, outline = TRUE, 
          whisklty = 1, staplelty = 0, # no whisker-ends, but with lines
          boxwex = 0.3,
          frame = FALSE)
  stripchart(log10(AvgCount) ~ CpG.Effect, data = res_rna_450k.tss1500, vertical = TRUE,
             method = "jitter", add = TRUE, pch = 20, col = "#1290D9", cex = 2.25)
dev.off()
  
pdf(paste0(CPGRNA_loc,"/",Today,".aerna.pilot.signCpG_vs_theirGenes.tss200.pdf"))
# postscript(paste0(CPGRNA_loc,"/",Today,".aerna.pilot.signCpG_vs_theirGenes.tss200.ps"))
# nice-to-have: stratify by type of CpG
par(mfrow = c(1,1), mar = c(6,6,4,0), oma = c(0, 0, 2, 0))

  ### TSS200
  gene.counts.tss200 <- res_rna_450k.tss200$AvgCount
  gene.methylation.tss200 <- as.factor(res_rna_450k.tss200$CpG.Effect)
  # levels(as.factor(res_rna_450k.tss200$CpG.Effect))
  gene.methylation.tss200 <- as.numeric(gene.methylation.tss200)
  gene.cor.tss200 <- signif(cor(gene.counts.tss200, gene.methylation.tss200, 
                                 use = "pairwise.complete.obs", method = "spearman"), 3)
  gene.cor.tss200
  wilcox.test(log10(AvgCount) ~ as.factor(CpG.Effect), data = res_rna_450k.tss200, 
              alternative = "two.sided")
  pval.tss200 = format(wilcox.test(log10(AvgCount) ~ as.factor(CpG.Effect), data = res_rna_450k.tss200, 
                                  alternative = "two.sided")$`p.value`, digits = 3, scientific = TRUE)
  
  boxplot(log10(AvgCount) ~ CpG.Effect, res_rna_450k.tss200,
          main = bquote(rho == .(gene.cor.tss200)~"-"~italic("P") == .(pval.tss200)),
          sub = "200bp from TSS",
          cex.sub = 1.75, cex.axis = 1.75, cex.main = 2.0, cex.lab = 2,
          ylim = c(-2.5, 3.5),
          xlab = "CpG effect",
          ylab = expression(log[10]~(normalized~counts)), 
          plot = TRUE, notch = FALSE, outline = TRUE, 
          whisklty = 1, staplelty = 0, # no whisker-ends, but with lines
          boxwex = 0.3,
          frame = FALSE)
  stripchart(log10(AvgCount) ~ CpG.Effect, data = res_rna_450k.tss200, vertical = TRUE,
             method = "jitter", add = TRUE, pch = 20, col = "#1290D9", cex = 2.25)
dev.off()
  
pdf(paste0(CPGRNA_loc,"/",Today,".aerna.pilot.signCpG_vs_theirGenes.5utr.pdf"))
# postscript(paste0(CPGRNA_loc,"/",Today,".aerna.pilot.signCpG_vs_theirGenes.5utr.ps"))
# nice-to-have: stratify by type of CpG
par(mfrow = c(1,1), mar = c(6,6,4,0), oma = c(0, 0, 2, 0))

  ### 5'UTR
  gene.counts.5utr <- res_rna_450k.5utr$AvgCount
  gene.methylation.5utr <- as.factor(res_rna_450k.5utr$CpG.Effect)
  # levels(as.factor(res_rna_450k.5utr$CpG.Effect))
  gene.methylation.5utr <- as.numeric(gene.methylation.5utr)
  gene.cor.5utr <- signif(cor(gene.counts.5utr, gene.methylation.5utr, 
                                 use = "pairwise.complete.obs", method = "spearman"), 3)
  gene.cor.5utr
  wilcox.test(log10(AvgCount) ~ as.factor(CpG.Effect), data = res_rna_450k.5utr, 
              alternative = "two.sided")
  pval.5utr = format(wilcox.test(log10(AvgCount) ~ as.factor(CpG.Effect), data = res_rna_450k.5utr, 
                                 alternative = "two.sided")$`p.value`, digits = 3, scientific = FALSE)
  
  boxplot(log10(AvgCount) ~ CpG.Effect, res_rna_450k.5utr,
          main = bquote(rho == .(gene.cor.5utr)~"-"~italic("P") == .(pval.5utr)),
          sub = "5'UTR",
          cex.sub = 1.75, cex.axis = 1.75, cex.main = 2.0, cex.lab = 2,
          ylim = c(-2.5, 3.5),
          xlab = "CpG effect",
          ylab = expression(log[10]~(normalized~counts)), 
          plot = TRUE, notch = FALSE, outline = TRUE, 
          whisklty = 1, staplelty = 0, # no whisker-ends, but with lines
          boxwex = 0.3,
          frame = FALSE)
  stripchart(log10(AvgCount) ~ CpG.Effect, data = res_rna_450k.5utr, vertical = TRUE,
             method = "jitter", add = TRUE, pch = 20, col = "#1290D9", cex = 2.25)
dev.off()
  
pdf(paste0(CPGRNA_loc,"/",Today,".aerna.pilot.signCpG_vs_theirGenes.1stexon.pdf"))
# postscript(paste0(CPGRNA_loc,"/",Today,".aerna.pilot.signCpG_vs_theirGenes.1stexon.ps"))
# nice-to-have: stratify by type of CpG
par(mfrow = c(1,1), mar = c(6,6,4,0), oma = c(0, 0, 2, 0))

  ### 1st Exon
  gene.counts.1stexon <- res_rna_450k.1stexon$AvgCount
  gene.methylation.1stexon <- as.factor(res_rna_450k.1stexon$CpG.Effect)
  # levels(as.factor(res_rna_450k.1stexon$CpG.Effect))
  gene.methylation.1stexon <- as.numeric(gene.methylation.1stexon)
  gene.cor.1stexon <- signif(cor(gene.counts.1stexon, gene.methylation.1stexon, 
                                 use = "pairwise.complete.obs", method = "spearman"), 3)
  gene.cor.1stexon
  wilcox.test(log10(AvgCount) ~ as.factor(CpG.Effect), data = res_rna_450k.1stexon, 
              alternative = "two.sided")
  pval.1stexon = format(wilcox.test(log10(AvgCount) ~ as.factor(CpG.Effect), data = res_rna_450k.1stexon, 
                                    alternative = "two.sided")$`p.value`, digits = 3, scientific = TRUE)
  
  boxplot(log10(AvgCount) ~ CpG.Effect, res_rna_450k.1stexon,
          main = bquote(rho == .(gene.cor.1stexon)~"-"~italic("P") == .(pval.1stexon)),
          sub = expression(1^st~exon),
          cex.sub = 1.75, cex.axis = 1.75, cex.main = 2.0, cex.lab = 2,
          ylim = c(-2.5, 3.5),
          xlab = "CpG effect",
          ylab = expression(log[10]~(normalized~counts)), 
          plot = TRUE, notch = FALSE, outline = TRUE, 
          whisklty = 1, staplelty = 0, # no whisker-ends, but with lines
          boxwex = 0.3,
          frame = FALSE)
  stripchart(log10(AvgCount) ~ CpG.Effect, data = res_rna_450k.1stexon, vertical = TRUE,
             method = "jitter", add = TRUE, pch = 20, col = "#1290D9", cex = 2.25)
dev.off()

pdf(paste0(CPGRNA_loc,"/",Today,".aerna.pilot.signCpG_vs_theirGenes.body.pdf"))
# postscript(paste0(CPGRNA_loc,"/",Today,".aerna.pilot.signCpG_vs_theirGenes.body.ps"))
# nice-to-have: stratify by type of CpG
par(mfrow = c(1,1), mar = c(6,6,4,0), oma = c(0, 0, 2, 0))

  ### BODY
  gene.counts.body <- res_rna_450k.body$AvgCount
  gene.methylation.body <- as.factor(res_rna_450k.body$CpG.Effect)
  # levels(as.factor(res_rna_450k.body$CpG.Effect))
  gene.methylation.body <- as.numeric(gene.methylation.body)
  gene.cor.body <- signif(cor(gene.counts.body, gene.methylation.body, 
                                 use = "pairwise.complete.obs", method = "spearman"), 3)
  gene.cor.body
  wilcox.test(log10(AvgCount) ~ as.factor(CpG.Effect), data = res_rna_450k.body, 
              alternative = "two.sided")
  pval.body = format(wilcox.test(log10(AvgCount) ~ as.factor(CpG.Effect), data = res_rna_450k.body, 
                                 alternative = "two.sided")$`p.value`, digits = 3, scientific = FALSE)
  
  boxplot(log10(AvgCount) ~ CpG.Effect, res_rna_450k.body,
          main = bquote(rho == .(gene.cor.body)~"-"~italic("P") == .(pval.body)),
          sub = "Body",
          cex.sub = 1.75, cex.axis = 1.75, cex.main = 2.0, cex.lab = 2,
          ylim = c(-2.5, 3.5),
          xlab = "CpG effect",
          ylab = expression(log[10]~(normalized~counts)),  
          plot = TRUE, notch = FALSE, outline = TRUE, 
          whisklty = 1, staplelty = 0, # no whisker-ends, but with lines
          boxwex = 0.3,
          frame = FALSE)
  stripchart(log10(AvgCount) ~ CpG.Effect, data = res_rna_450k.body, vertical = TRUE,
             method = "jitter", add = TRUE, pch = 20, col = "#1290D9", cex = 2.25)
dev.off()
  
pdf(paste0(CPGRNA_loc,"/",Today,".aerna.pilot.signCpG_vs_theirGenes.3utr.pdf"))
# postscript(paste0(CPGRNA_loc,"/",Today,".aerna.pilot.signCpG_vs_theirGenes.3utr.ps"))
# nice-to-have: stratify by type of CpG
par(mfrow = c(1,1), mar = c(6,6,4,0), oma = c(0, 0, 2, 0))

  ### 3'UTR
  gene.counts.3utr <- res_rna_450k.3utr$AvgCount
  gene.methylation.3utr <- as.factor(res_rna_450k.3utr$CpG.Effect)
  # levels(as.factor(res_rna_450k.3utr$CpG.Effect))
  gene.methylation.3utr <- as.numeric(gene.methylation.3utr)
  gene.cor.3utr <- signif(cor(gene.counts.3utr, gene.methylation.3utr, 
                                use = "pairwise.complete.obs", method = "spearman"), 3)
  gene.cor.3utr
  wilcox.test(log10(AvgCount) ~ as.factor(CpG.Effect), data = res_rna_450k.3utr, 
              alternative = "two.sided")
  pval.3utr = format(wilcox.test(log10(AvgCount) ~ as.factor(CpG.Effect), data = res_rna_450k.3utr, 
                                 alternative = "two.sided")$`p.value`, digits = 3, scientific = FALSE)
  
  boxplot(log10(AvgCount) ~ CpG.Effect, res_rna_450k.3utr,
          main = bquote(rho == .(gene.cor.3utr)~"-"~italic("P") == .(pval.3utr)),
          sub = "3'UTR",
          cex.sub = 1.75, cex.axis = 1.75, cex.main = 2.0, cex.lab = 2,
          ylim = c(-2.5, 3.5),
          xlab = "CpG effect",
          ylab = expression(log[10]~(normalized~counts)), 
          plot = TRUE, notch = FALSE, outline = TRUE, 
          whisklty = 1, staplelty = 0, # no whisker-ends, but with lines
          boxwex = 0.3,
          frame = FALSE)
  stripchart(log10(AvgCount) ~ CpG.Effect, data = res_rna_450k.3utr, vertical = TRUE,
             method = "jitter", add = TRUE, pch = 20, col = "#1290D9", cex = 2.25)
dev.off()
  
pdf(paste0(CPGRNA_loc,"/",Today,".aerna.pilot.signCpG_vs_theirGenes.intergenic.pdf"))
# postscript(paste0(CPGRNA_loc,"/",Today,".aerna.pilot.signCpG_vs_theirGenes.intergenic.ps"))
# nice-to-have: stratify by type of CpG
par(mfrow = c(1,1), mar = c(6,6,4,0), oma = c(0, 0, 2, 0))

  ### INTERGENIC
  gene.counts.rest <- res_rna_450k.rest$AvgCount
  gene.methylation.rest <- as.factor(res_rna_450k.rest$CpG.Effect)
  # levels(as.factor(res_rna_450k.rest$CpG.Effect))
  gene.methylation.rest <- as.numeric(gene.methylation.rest)
  gene.cor.rest <- signif(cor(gene.counts.rest, gene.methylation.rest, 
                              use = "pairwise.complete.obs", method = "spearman"), 3)
  gene.cor.rest
  wilcox.test(log10(AvgCount) ~ as.factor(CpG.Effect), data = res_rna_450k.rest, 
              alternative = "two.sided")
  pval.rest = format(wilcox.test(log10(AvgCount) ~ as.factor(CpG.Effect), data = res_rna_450k.rest, 
                                 alternative = "two.sided")$`p.value`, digits = 3, scientific = FALSE)
  
  boxplot(log10(AvgCount) ~ CpG.Effect, res_rna_450k.rest,
          main = bquote(rho == .(gene.cor.rest)~"-"~italic("P") == .(pval.rest)),
          sub = "Intergenic",
          cex.sub = 1.75, cex.axis = 1.75, cex.main = 2.0, cex.lab = 2,
          ylim = c(-2.5, 3.5),
          xlab = "CpG effect",
          ylab = expression(log[10]~(normalized~counts)), 
          plot = TRUE, notch = FALSE, outline = TRUE, 
          whisklty = 1, staplelty = 0, # no whisker-ends, but with lines
          boxwex = 0.3,
          frame = FALSE)
  stripchart(log10(AvgCount) ~ CpG.Effect, data = res_rna_450k.rest, vertical = TRUE,
             method = "jitter", add = TRUE, pch = 20, col = "#1290D9", cex = 2.25)
dev.off()
  
par(mfrow = c(1,1), mar = c(5, 4, 4, 2), oma = c(0, 0, 0, 0))

cat("\n===========================================================================================")
cat("SAVE THE DATA")

save.image(paste0(ANALYSIS_loc,"/",Today,".aems450k.meta.analysis.ewas.",EWAS_trait,".CpG_vs_RNA.RData"))

