cat("===========================================================================================
                META-ANALYSIS ATHERO-EXPRESS METHYLATION STUDIES 450K 1 & 2
                  -- RNA Expression in Carotid Plaques: A Pilot Study --
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

cat("\n===========================================================================================")
cat("LOAD RNAseq pilot data")
cat("Experimental setup: 30 samples were selected based on the extremes of the most significant
    CpG; so extremely demethylated 15 samples and methylated 15 samples.")  

# References: 
# - https://f1000research.com/articles/4-1070/v2
# - http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
cat("\nLoad datasets...")
INPRNA_loc = paste0(ROOT_loc, "/PLINK/_AE_Originals/AERNASMOKE")
PROJECTNAME = "160223_NS500813_0095_AHT2FKBGXX"
DATADIR = paste0(INPRNA_loc, "/",PROJECTNAME)

cat("\n* raw, mapped data...")
RAWDATA = read.table(paste0(DATADIR, "/read_counts/",PROJECTNAME,"_readCounts_raw.txt"), header = TRUE, 
                     sep = "\t", row.names = 1)
cat("\n* raw, mapped data normalized for the read depth and gene length...")
RPKMDATA = read.table(paste0(DATADIR, "/read_counts/",PROJECTNAME,"_readCounts_RPKM.txt"), header = TRUE, 
                      sep = "\t", row.names = 1)
cat("\n* raw, mapped data normalized for the read depth...")
NORMDATA = read.table(paste0(DATADIR, "/read_counts/",PROJECTNAME,"_readCounts_normalized.txt"), header = TRUE, 
                      sep = "\t", row.names = 1)

cat("\n* sample information...")
coldata = read.table(paste0(INPRNA_loc, "/samples.phenotypes.txt"), header = TRUE, sep = "\t",
                     row.names = 1)

cat("\n* loading AE database...\n")
load(paste0(OUT_loc,"/20180207.AEData_2016.update.DF.CEA.TFvar.RData"))

AEData_2016.update.DF.CEA.SampleSize = nrow(AEData_2016.update.DF.CEA.TFvar)

AEData_2016.update.DF.CEA = AEData_2016.update.DF.CEA.TFvar

cat("\n* read in BAM-files and associated sample information...")
cat("\n  - sample information...")
sampleTable = read.table(paste0(INPRNA_loc, "/samples.phenotypes.txt"), header = TRUE, sep = "\t",
                         row.names = 1)

cat("\n  - list BAM-files...")
# Allignment was done using GATK/STAR and "Homo_sapiens.GRCh37.74.gtf"
filenames <- file.path(paste0(DATADIR, "/", # root folder
                              sampleTable$RunID, # sample folder
                              "/mapping/", # sample subfolder
                              sampleTable$RunID, "_sorted.bam")) # sample BAM-file
cat("\n  - check existence of BAM-files...")
file.exists(filenames)
library("Rsamtools")
bamfiles <- BamFileList(filenames, yieldSize = 2000000)
# check the chromosome-nomenclature: "1" or "chr1" under 'seqnames'
seqinfo(bamfiles[1])

cat("\n===========================================================================================")
cat("DEFINE GENE MODELS")
cat("\nDefining the gene models for the BAM-files...")
# Here we will demonstrate loading from a GTF file:
library("GenomicFeatures")
# We indicate that none of our sequences (chromosomes) are circular using a 0-length character vector.
# gzipped files can be processed - saves space.
cat("\n* get the human genome data, GRch37.74")
gtffile <- file.path(paste0(DATADIR, "/Homo_sapiens.GRCh37.74.gtf"))
txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character(),
                        organism = "Homo sapiens", taxonomyId = 9)
txdb
columns(txdb)
keytypes(txdb)

cat("\n* get the exons")
(ebg <- exonsBy(txdb, by = "gene"))

cat("\n===========================================================================================")
cat("READ COUNTING")
cat("\nRead counting the shizzle...")
library("GenomicAlignments")
install.packages.auto("Rsubread")

# ### Version 1 - very fast counting, approx. 1 minute per BAM file
# countdata <- featureCounts(files = filenames, 
#                            annot.ext = gtffile,
#                            isGTFAnnotationFile = TRUE,
#                            GTF.featureType = "exon",
#                            GTF.attrType = "gene_id", 
#                            nthreads = 8)
# 
# ### Changing the column headers of the targets, for matching with 'coldata' downstream.
# countdataORIGINAL = countdata
# countdata$targets <- sub("X.Users.swvanderlaan..PLINK._AE_Originals.AERNASMOKE.160223_NS500813_0095_AHT2FKBGXX.", "", countdata$targets)
# countdata$targets <- sub(pattern = ".mapping.*", replacement = "", x = countdata$targets, perl = TRUE)
# colnames(countdata$counts) <- sub("X.Users.swvanderlaan..PLINK._AE_Originals.AERNASMOKE.160223_NS500813_0095_AHT2FKBGXX.", "", colnames(countdata$counts))
# colnames(countdata$counts) <- sub(pattern = ".mapping.*", replacement = "", x = colnames(countdata$counts), perl = TRUE)
# 
### Version 2 - substantially slower, also BiocParallel is needed - prefered as it creates a SE
library("BiocParallel")
register(MulticoreParam(workers = 2))
register(MulticoreParam())
se <- summarizeOverlaps(features = ebg, reads = bamfiles,
                        mode = "Union",
                        singleEnd = TRUE,
                        ignore.strand = FALSE,
                        inter.feature = TRUE,
                        fragments = FALSE )
(colData(se) <- DataFrame(sampleTable))

se.ignore.strand <- summarizeOverlaps(features = ebg, reads = bamfiles,
                                      mode = "Union",
                                      singleEnd = TRUE,
                                      ignore.strand = TRUE,
                                      inter.feature = TRUE,
                                      fragments = FALSE )
(colData(se.ignore.strand) <- DataFrame(sampleTable))

cat("\n===========================================================================================")
cat("RNA SEQUENCING ANALYSIS")
### References:
### - http://www.bioconductor.org/help/workflows/rnaseqGene/#exploratory-analysis-and-visualization

cat("\nQuality Control.\n")
# FUTURE VERSIONS:
# - check % mapped mRNAs as a function of 3" bias
# - check sample sex: RNA-seq vs. clinical information using 'omniPrint' 
# - visualize sex: plot X- and Y-chromosome specific genes per sex, e.g. median/mean
#                  which one would not expect for females/males for certain genes
# - check sample mix-ups: use available genotype data
cat("\n* % mapped mRNAs...\n")
  plot(sampleTable$PCT_CORRECT_STRAND_READS, sampleTable$PCT_MRNA_BASES, main = "",
       ylab = "% mapped mRNAs",
       xlab = "% correct strand reads",
       xlim = c(0.5,1),
       axes = FALSE, col = "#1290D9", pch = 19,
       bty = "n")
  # draw an axis on the left 
  axis(2, at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), 
       labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), 
       las = 2)
  # draw an axis on the bottom 
  axis(1, at = c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), 
       labels = c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), 
       las = 2)
  # add dashed for minimal mRNA % threshold
  abline(h = 0.05, lty = 2, col = "#595A5C")
  abline(v = 0.90, lty = 2, col = "#E55738")
dev.copy2pdf(file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.Percentage_mRNA_mapped_vs_corr_reads.pdf"), onefile = TRUE,
             width = 10, height = 10,
             paper = "a4r")

cat("\n* density plot of raw read counts (log10)...\n")
samples <- as.character(sampleTable$STUDYNUMBER)
logcounts <- log(assay(se)[,1] + 1, 10)
d <- density(logcounts)
plot(d, main = "", 
     xlim = c(1, 6), ylim = c(0, 1),
     xlab = "raw read counts per gene [log10]", 
     ylab = "density", 
     col = "#FBB820",
     lty = 1, lwd = 2,
     bty = "n")
for (s in 2:length(samples)) {
  logcounts <- log(assay(se)[,s] + 1, 10) 
  d <- density(logcounts)
  lines(d, col = uithof_color_legend[s], 
        lty = c(1,1,1,1,1,1,1,1,1,
                1,1,1,1,1,1,1,1,1,1,
                1,1,1,1,1,1,1,1,2,3), # gives the legend appropriate symbols (lines)
        lwd = c(2,2,2,2,2,2,2,2,2,2,
                2,2,2,2,2,2,2,2,2,2,
                2,2,2,2,2,2,2,2,2))
}
legend("topright", # places a legend at the appropriate place  
       legend = samples, # puts text in the legend
       lty = c(1,1,1,1,1,1,1,1,1,1,
               1,1,1,1,1,1,1,1,1,1,
               1,1,1,1,1,1,1,1,2,3), # gives the legend appropriate symbols (lines)
       lwd = c(2,2,2,2,2,2,2,2,2,2,
               2,2,2,2,2,2,2,2,2,2,
               2,2,2,2,2,2,2,2,2,2), # gives the legend appropriate symbols (lines)
       col = uithof_color_legend,
       title = "Samples", bty = "n") # gives the legend lines the correct color and width
dev.copy2pdf(file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.Raw_read_counts_per_gene.density.pdf"), onefile = TRUE,
             width = 10, height = 10,
             paper = "a4")
# remove temporary files
rm(d, s, samples, logcounts)

cat("\n* box plot of raw read counts (log10)...\n")
logcounts <- log(assay(se) + 1,10)
samples <- as.character(sampleTable$STUDYNUMBER)
boxplot(logcounts, main = "", 
        xlab = "", ylab = "Raw read counts per gene [log10]", 
        col = uithof_color,
        axes = FALSE)
axis(2)
axis(1, at = c(1:length(samples)),
     labels = colnames(logcounts), 
     las = 2, cex.axis = 0.8)
dev.copy2pdf(file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.Raw_read_counts_per_gene.boxplot.pdf"), onefile = TRUE,
             width = 10, height = 10,
             paper = "a4")
# remove temporary files
rm(logcounts, samples)

cat("\n* filter samples with very low % mapped mRNAs...\n")
sampleTableFiltered <- subset(sampleTable, PCT_MRNA_BASES > 0.05 & PCT_CORRECT_STRAND_READS > 0.9 )
sampleTableFiltered
keeps <- rownames(sampleTableFiltered)
library(SummarizedExperiment)
metadata(se)
# assays(se)$counts
rowRanges(se)
# colData(se)
se_f <- subset(se)[,keeps]
se_f
dim(se_f)

cat("\n* filter out genes with low counts...\n")
# Reference: https://cdn.rawgit.com/DarwinAwardWinner/resume/master/examples/Salomon/Teaching/RNA-Seq%20Lab.html
library(edgeR)
total_counts <- colSums(assays(se)$counts)
print(total_counts)
summary(total_counts)
nf <- calcNormFactors(assays(se)$counts, lib.size = total_counts, method = "TMM") # Trimmed Mean of M-values (TMM) method
print(nf)
summary(nf)
normalized_total_counts <- total_counts * nf
print(normalized_total_counts)
summary(normalized_total_counts)
mean_log_cpm <- aveLogCPM(assays(se)$counts, normalized_total_counts)
print(mean_log_cpm)
summary(mean_log_cpm)
filter_threshold <- 1
hist(mean_log_cpm, main = "Mean counts per million",
     xlim = c(0, 20),
     xlab = "mean log(CpM)",
     ylab = "counts",
     col = "#1290D9", bty = "n", breaks = 50)
abline(v = filter_threshold, col = "#E55738", lwd = 1, lty = 2)
dev.copy2pdf(file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.MeanCounts.QC.Histogram.pdf"), onefile = TRUE,
             width = 10, height = 10,
             paper = "a4")
qqnorm(mean_log_cpm, bty = "n")
abline(h = filter_threshold, col = "#E55738", lwd = 1, lty = 2)
dev.copy2pdf(file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.MeanCounts.QC.QQplot.pdf"), onefile = TRUE,
             width = 10, height = 10,
             paper = "a4")

keep_genes <- mean_log_cpm >= 2.5
keep_genes
se_fqc <- subset(se_f)[keep_genes,]
dim(se) # 63,677 genes; 30 samples
dim(se_f) # 63,677 genes; 21 samples
dim(se_fqc) # 30,488 genes; 21 samples
allgenescounts = nrow(se)
cat(paste0("\n  - all genes after sample QC: ", format(allgenescounts, big.mark = ","),""))
genecountsfiltered = nrow(se_fqc)
cat(paste0("\n  - genes after sample and expression QC: ", format(genecountsfiltered, big.mark = ","),""))

rm(total_counts, nf, normalized_total_counts, filter_threshold)

cat("\nGet some baseline characteristics of the RNAseq pilot data...\n")
# colnames(coldata)

# Baseline table variables
basetable_vars = c("Age", "sex", 
                   "Creat", "Homocyst", "Glucose", "FolicAcid", "VitB12", "CRP", 
                   "SBP", "DBP", 
                   "TC_final", "LDL_final", "HDL_final", "TG_final", 
                   "MAP", "PulsePressure", "eGFR", "BMI", 
                   "CKD", 
                   "T2D", "Hypertension", 
                   "Antihypertensives", "Anticoagulants", "Antiplatelets", "OralGlucInh", "LLDs", 
                   "Stroke_history", "Symptoms.4g", "EP_major")
basetable_bin = c("sex", 
                  "CKD", 
                  "T2D", "Hypertension", 
                  "Antihypertensives", "Anticoagulants", "Antiplatelets", "OralGlucInh", "LLDs", 
                  "Stroke_history", "Symptoms.4g", "EP_major")
basetable_con = basetable_vars[!basetable_vars %in% basetable_bin]

cat("Check for Normality of continues variables.")
pdf(paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,".BaselineNormality.RAW.pdf"), 
    paper = "a4", onefile = TRUE)
par(mfrow = c(2,2))
for (x in basetable_con) {
  cat(paste0("* plotting data for: [ ",x," ].\n"))
  qqnorm(colData(se)[,x], 
         main = x, 
         cex.main = 1.0, cex.lab = 0.90, cex.axis = 0.80, cex = 1.25,
         col = "#1290D9", 
         pch = 20, bty = "n")
  qqline(colData(se)[,x], col = "#E55738", lty = 2, lwd = 2)
}
dev.off()

cat("* Plotting histograms for all data...")
pdf(paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,".Transformation.Histogram.pdf"), 
    paper = "a4r", onefile = TRUE)
par(mfrow = c(2,2), mar = c(3,3,3,1))
for (x in basetable_con) {
  hist(colData(se)[,x], main = x, 
       xlab = "",
       col = "#1290D9", 
       cex.main = 1.0, cex.lab = 0.90, cex.axis = 0.90, cex = 1.25)
  mtext(table(is.na(colData(se)[,x]))["FALSE"], cex = 0.90)
  mtext(round(table(colData(se)[,x] == 0)["TRUE"]/table(is.na(colData(se)[,x]))["FALSE"]*100,2), 
        cex = 0.90, side = 4, line = -1)
}
dev.off()

# Create baseline tables
se.tableOne = print(CreateTableOne(vars = basetable_vars, strata = "CurSmoking", 
                                   data = as.data.frame(colData(se)), 
                                   factorVars = basetable_bin), 
                    nonnormal = c(), quote = FALSE, showAllLevels = FALSE, cramVars = "sex",
                    format = "p", contDigits = 1)[,1:3]

# Write basetable
write.xlsx(file = paste0(RNAOUT_loc, "/",Today,".aerna.pilot.",EWAS_trait,".Basetable.xlsx"), 
           format(se.tableOne, digits = 5, scientific = FALSE) , row.names = TRUE, col.names = TRUE)

# Check for NA, add manually to table
CLIN = colData(se)[, as.integer(labels(colnames(colData(se))))]
CLIN = CLIN[!is.na(CLIN$CurSmoking),]
for (x in basetable_vars) {
  print(paste(x, "smoker", table(is.na(CLIN[CLIN$CurSmoking == "yes", x]))["FALSE"]))
  print(paste(x, "non-smoker", table(is.na(CLIN[CLIN$CurSmoking == "no", x]))["FALSE"]))
}
table(CLIN[CLIN$CurSmoking == "yes", "Symptoms.4g"], useNA = "always")
table(CLIN[CLIN$CurSmoking == "no", "Symptoms.4g"], useNA = "always")

cat("\nExploratory data and differential expression analysis...\n")
library("DESeq2")
library("dplyr")
library("ggplot2")
cat("\n  - all raw data")
dds <- DESeqDataSet(se, design = ~ CurSmoking)
dds <- DESeq(dds)
dds_alt <- dds[ rowSums(counts(dds)) > 1, ] # remove genes with low counts (<1)
dds_alt <- DESeq(dds_alt)
cat("\n  - after sample QC")
dds_qc <- DESeqDataSet(se_f, design = ~ CurSmoking)
dds_qc <- DESeq(dds_qc)
dds_qc_alt <- dds_qc[ rowSums(counts(dds_qc)) > 1, ] # remove genes with low counts (<1)
dds_qc_alt <- DESeq(dds_qc_alt)
cat("\n  - after sample and expression QC")
dds_qcf <- DESeqDataSet(se_fqc, design = ~ CurSmoking)
dds_qcf <- DESeq(dds_qcf)
dds_qcf_alt <- dds_qcf[ rowSums(counts(dds_qcf)) > 1, ] # remove genes with low counts (<1)
dds_qcf_alt <- DESeq(dds_qcf_alt)

cat("\n* Plotting dispersions...")
# More on dispersions can be found here: https://support.bioconductor.org/p/75260/

png(file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.Dispersions.QC.png"),
    width = 1920, height = 1080)
pdf(file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.Dispersions.QC.pdf"), onefile = TRUE,
        width = 10, height = 10,
        paper = "a4")

  par(mfrow = c(3,2), mar = c(5,8,4,2), oma = c(0, 0, 2, 0))

  plotDispEsts(dds, main = "Raw data")
  plotDispEsts(dds_alt, main = "Raw data, filtered gene counts < 1")
  plotDispEsts(dds_qc, main = "QC'd data")
  plotDispEsts(dds_qc_alt, main = "QC'd data, filtered gene counts < 1")
  plotDispEsts(dds_qcf, main = "QC'd & filtered data")
  plotDispEsts(dds_qcf_alt, main = "QC'd & filtered data, filtered gene counts < 1")
  
  mtext("Expression dispersions", outer = TRUE, cex = 1.1)

dev.off()
par(mfrow = c(1,1), mar = c(5, 4, 4, 2), oma = c(0, 0, 0, 0))

cat("\n* Normalise data based on 'regularized-logarithm transformation' (only required for PC plotting!)...")
rld <- rlog(dds, blind = FALSE)
rld_alt <- rlog(dds_alt, blind = FALSE)
rld_qc <- rlog(dds_qc, blind = FALSE)
rld_qc_alt <- rlog(dds_qc_alt, blind = FALSE)
rld_qcf <- rlog(dds_qcf, blind = FALSE)
rld_qcf_alt <- rlog(dds_qcf_alt, blind = FALSE)

vsd <- vst(dds, blind = FALSE)
vsd_alt <- vst(dds_alt, blind = FALSE)
vsd_qc <- vst(dds_qc, blind = FALSE)
vsd_qc_alt <- vst(dds_qc_alt, blind = FALSE)
vsd_qcf <- vst(dds_qcf, blind = FALSE)
vsd_qcf_alt <- vst(dds_qcf_alt, blind = FALSE)

cat("\n* Sample distances...")
cat("\n  - pheatmap and 'Euclidean' distances")
library("pheatmap")
library("RColorBrewer")
pdf(file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.heatmap.QC.pdf"), onefile = FALSE,
    width = 10, height = 10,
    paper = "a4")
  sampleDists <- dist(t(assay(vsd_qc)))
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- paste(dds_qc$STUDYNUMBER, dds_qc$IPH, sep = "|")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors)
dev.off()
rm(sampleDists, sampleDistMatrix, colors)

library("PoiClaClu")
library("gplots")
cat("\n  - pheatmap and 'Poisson' distances based on counts")
pdf(file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.heatmap.poisson.QC.pdf"), onefile = FALSE,
    width = 10, height = 10,
    paper = "a4")
  poisd <- PoissonDistance(t(counts(dds)))
  samplePoisDistMatrix <- as.matrix( poisd$dd )
  rownames(samplePoisDistMatrix) <- paste(dds$STUDYNUMBER, dds$IPH, sep = "|")
  colnames(samplePoisDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  hc <- hclust(poisd$dd)
  heatmap.2( samplePoisDistMatrix, Rowv = as.dendrogram(hc),
             symm = TRUE, trace = "none", col = colors,
             margins = c(2,10), labCol = FALSE )
dev.off()
rm(poisd, hc, samplePoisDistMatrix, colors)

cat("\n* Principal component analysis...")
# MAYBE ADD IN: PCA with "age", and with "% mRNA mapped"
cat("\n  - using build in function")
pdf(file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.pca.Sex.pdf"), onefile = TRUE,
    width = 10, height = 10,
    paper = "a4")
  plotPCA(rld_qc, intgroup = c("CurSmoking", "sex"))
dev.off()
pdf(file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.pca.Age.pdf"), onefile = TRUE,
    width = 10, height = 10,
    paper = "a4")
  plotPCA(rld_qc, intgroup = c("CurSmoking", "Age"))
dev.off()
pdf(file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.pca.PCT_MRNA_BASES.pdf"), onefile = TRUE,
    width = 10, height = 10,
    paper = "a4")
  plotPCA(rld_qc, intgroup = c("CurSmoking", "PCT_MRNA_BASES"))
dev.off()
pdf(file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.pca.CurSmoking.pdf"), onefile = TRUE,
    width = 10, height = 10,
    paper = "a4")
  plotPCA(rld_qc, intgroup = c("CurSmoking"))
dev.off()

cat("\n* Get differential expression results...")
cat("\n  - collecting results...")
res <- results(dds)
res
# table(res$padj<0.05)
res_alt <- results(dds_alt)
res_alt
# table(res_alt$padj<0.05)
res_qc <- results(dds_qc)
res_qc
# table(res_qc$padj<0.05)
res_qc_alt <- results(dds_qc_alt)
res_qc_alt
# table(res_qc_alt$padj<0.05)
res_qcf <- results(dds_qcf)
res_qcf
# table(res_qcf$padj<0.05)
res_qcf_alt <- results(dds_qcf_alt)
res_qcf_alt
# table(res_qcf_alt$padj<0.05)

cat("\n  - ordering results on adjusted p-value...")
res <- res[order(res$padj), ]
res_alt <- res[order(res_alt$padj), ]
res_qc <- res[order(res_qc$padj), ]
res_qc_alt <- res[order(res_qc_alt$padj), ]
res_qcf <- res[order(res_qcf$padj), ]
res_qcf_alt <- res[order(res_qcf_alt$padj), ]

cat("\n  - merge with normalized count data...")
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized = TRUE)), 
                 by = "row.names", sort = FALSE)
names(resdata)[1] <- "ENSEMBLID"
res_altdata <- merge(as.data.frame(res_alt), as.data.frame(counts(dds_alt, normalized = TRUE)), 
                     by = "row.names", sort = FALSE)
names(res_altdata)[1] <- "ENSEMBLID"
res_qcdata <- merge(as.data.frame(res_qc), as.data.frame(counts(dds_qc, normalized = TRUE)), 
                    by = "row.names", sort = FALSE)
names(res_qcdata)[1] <- "ENSEMBLID"
res_qc_altdata <- merge(as.data.frame(res_qc_alt), as.data.frame(counts(dds_qc_alt, normalized = TRUE)), 
                        by = "row.names", sort = FALSE)
names(res_qc_altdata)[1] <- "ENSEMBLID"
res_qcfdata <- merge(as.data.frame(res_qcf), as.data.frame(counts(dds_qcf, normalized = TRUE)), 
                     by = "row.names", sort = FALSE)
names(res_qcfdata)[1] <- "ENSEMBLID"
res_qcf_altdata <- merge(as.data.frame(res_qcf_alt), as.data.frame(counts(dds_qcf_alt, normalized = TRUE)), 
                         by = "row.names", sort = FALSE)
names(res_qcf_altdata)[1] <- "ENSEMBLID"

cat("\n  - annotating results...")
### DESeqRESQC - normalized data, exclusion of PC outlier, age+sex+mRNA%+3Prime%
### DESeqRESQCmin - normalized data, exclusion of PC outlier
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)

resdata$symbol <- mapIds(org.Hs.eg.db,
                         keys = resdata$ENSEMBLID,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
resdata$entrezid <- mapIds(org.Hs.eg.db,
                           keys = resdata$ENSEMBLID,
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")
resdata$genename <- mapIds(org.Hs.eg.db,
                           keys = resdata$ENSEMBLID,
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")

res_altdata$symbol <- mapIds(org.Hs.eg.db,
                             keys = res_altdata$ENSEMBLID,
                             column = "SYMBOL",
                             keytype = "ENSEMBL",
                             multiVals = "first")
res_altdata$entrezid <- mapIds(org.Hs.eg.db,
                               keys = res_altdata$ENSEMBLID,
                               column = "SYMBOL",
                               keytype = "ENSEMBL",
                               multiVals = "first")
res_altdata$genename <- mapIds(org.Hs.eg.db,
                               keys = res_altdata$ENSEMBLID,
                               column = "SYMBOL",
                               keytype = "ENSEMBL",
                               multiVals = "first")

res_qcdata$symbol <- mapIds(org.Hs.eg.db,
                            keys = res_qcdata$ENSEMBLID,
                            column = "SYMBOL",
                            keytype = "ENSEMBL",
                            multiVals = "first")
res_qcdata$entrezid <- mapIds(org.Hs.eg.db,
                              keys = res_qcdata$ENSEMBLID,
                              column = "SYMBOL",
                              keytype = "ENSEMBL",
                              multiVals = "first")
res_qcdata$genename <- mapIds(org.Hs.eg.db,
                              keys = res_qcdata$ENSEMBLID,
                              column = "SYMBOL",
                              keytype = "ENSEMBL",
                              multiVals = "first")

res_qc_altdata$symbol <- mapIds(org.Hs.eg.db,
                                keys = res_qc_altdata$ENSEMBLID,
                                column = "SYMBOL",
                                keytype = "ENSEMBL",
                                multiVals = "first")
res_qc_altdata$entrezid <- mapIds(org.Hs.eg.db,
                                  keys = res_qc_altdata$ENSEMBLID,
                                  column = "SYMBOL",
                                  keytype = "ENSEMBL",
                                  multiVals = "first")
res_qc_altdata$genename <- mapIds(org.Hs.eg.db,
                                  keys = res_qc_altdata$ENSEMBLID,
                                  column = "SYMBOL",
                                  keytype = "ENSEMBL",
                                  multiVals = "first")

res_qcfdata$symbol <- mapIds(org.Hs.eg.db,
                             keys = res_qcfdata$ENSEMBLID,
                             column = "SYMBOL",
                             keytype = "ENSEMBL",
                             multiVals = "first")
res_qcfdata$entrezid <- mapIds(org.Hs.eg.db,
                               keys = res_qcfdata$ENSEMBLID,
                               column = "SYMBOL",
                               keytype = "ENSEMBL",
                               multiVals = "first")
res_qcfdata$genename <- mapIds(org.Hs.eg.db,
                               keys = res_qcfdata$ENSEMBLID,
                               column = "SYMBOL",
                               keytype = "ENSEMBL",
                               multiVals = "first")

res_qcf_altdata$symbol <- mapIds(org.Hs.eg.db,
                                 keys = res_qcf_altdata$ENSEMBLID,
                                 column = "SYMBOL",
                                 keytype = "ENSEMBL",
                                 multiVals = "first")
res_qcf_altdata$entrezid <- mapIds(org.Hs.eg.db,
                                   keys = res_qcf_altdata$ENSEMBLID,
                                   column = "SYMBOL",
                                   keytype = "ENSEMBL",
                                   multiVals = "first")
res_qcf_altdata$genename <- mapIds(org.Hs.eg.db,
                                   keys = res_qcf_altdata$ENSEMBLID,
                                   column = "SYMBOL",
                                   keytype = "ENSEMBL",
                                   multiVals = "first")

cat("\n  - writing this shizzle...")
library(data.table)
library(ReportingTools)
fwrite(resdata, file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.raw.",EWAS_trait,".txt"),
       sep = ";", na = "NA", dec = ".",
       row.names = FALSE, col.names = TRUE,
       showProgress = TRUE, verbose = TRUE)
htmlRep <- HTMLReport(shortName = "RawDataReport", 
                      basePath = RNAOUT_loc,
                      title = paste0("Athero-Express RNAseq Pilot: Current smoking (normalized, raw data, ",Today,")"),
                      reportDirectory = "./ae.rna.pilot.report")
publish(resdata, htmlRep)
url <- finish(htmlRep)

fwrite(res_qcdata, file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.qc.",EWAS_trait,".txt"),
       sep = ";", na = "NA", dec = ".",
       row.names = FALSE, col.names = TRUE,
       showProgress = TRUE, verbose = TRUE)
htmlRepQC <- HTMLReport(shortName = "QCDataReport", 
                        basePath = RNAOUT_loc,
                        title = paste0("Athero-Express RNAseq Pilot: Current smoking (normalized, sample QC data, ",Today,")"),
                        reportDirectory = "./ae.rna.pilot.report.QC")
publish(res_qcdata, htmlRepQC)
url <- finish(htmlRepQC)

fwrite(res_qcfdata, file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.qc.filtered.",EWAS_trait,".txt"),
       sep = ";", na = "NA", dec = ".",
       row.names = FALSE, col.names = TRUE,
       showProgress = TRUE, verbose = TRUE)
htmlRepQCF <- HTMLReport(shortName = "QCFDataReport", 
                         basePath = RNAOUT_loc,
                         title = paste0("Athero-Express RNAseq Pilot: Current smoking (normalized, sample QC, filtered gene counts data, ",Today,")"),
                         reportDirectory = "./ae.rna.pilot.report.QCF")
publish(res_qcfdata, htmlRepQCF)
url <- finish(htmlRepQCF)

## Examine plot of p-values
png(file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,"PvalHistograms.png"),
    width = 1920, height = 1080)
pdf(file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,"PvalHistograms.pdf"), onefile = TRUE,
        paper = "a4r")
  par(mfrow = c(2,2), mar = c(5,8,4,2), oma = c(0, 0, 2, 0))
  # example code
  # main = bquote(Pearson~r^2 == .(pearsonR.ahrr4))
  # ylab = bquote(.(top.ahrr4)~" -"~italic("AHRR")) 
  hist(-log10(resdata$pvalue), main = "Raw data\n(gene ~ CurSmoke)", 
       xlab = expression(log[10]~italic((P))), breaks = 50, col = uithof_color[1])
  hist(-log10(res_qcdata$pvalue), main = "After sample QC\n(gene ~ CurSmoke)", 
       xlab = expression(log[10]~italic((P))), breaks = 50, col = uithof_color[4])
  hist(-log10(res_qcfdata$pvalue), main = "After sample & gene-count QC\n(gene ~ CurSmoke)", 
       xlab = expression(log[10]~italic((P))), breaks = 50, col = uithof_color[16])

  mtext("Examination of p-values", outer = TRUE, cex = 1.25)

dev.off()
par(mfrow = c(1,1), mar = c(5, 4, 4, 2), oma = c(0, 0, 0, 0))

# ## Examine independent filtering
# attr(res, "filterThreshold")
# plot(attr(res,"filterNumRej"), type = "b", xlab = "quantiles of baseMean", ylab = "number of rejections")

cat("\n  - MA plots...")
# An MA plot is an application of a Bland–Altman plot for visual representation of genomic data. 
# The plot visualises the differences between measurements taken in two samples, 
# by transforming the data onto M (log ratio) and A (mean average) scales, then plotting these values.
## Could do with built-in DESeq2 function:
# DESeq2::plotMA(dds_qcf, cex=1)
# function to make MA plot
library(calibrate)
maplot <- function(res, thresh = 0.05, labelsig = TRUE, textcx = s1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch = 20, cex = 0.5, log = "x", ...))
  with(subset(res, padj < thresh), points(baseMean, log2FoldChange, col = uithof_color[3], 
                                        pch = 20, cex = 1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj < thresh), textxy(baseMean, log2FoldChange, labs = symbol, 
                                          cex = textcx, col = 2))
  }
}
png(file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,".MAplots.png"),
    width = 1920, height = 1080)
pdf(file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,".MAplots.pdf"), onefile = TRUE,
        paper = "a4r")
  par(mfrow = c(2,2), mar = c(5,8,4,2), oma = c(0, 0, 2, 0))
  maplot(resdata, threshold = 0.1, labelsign = TRUE, main = "Raw data", bty = "n",
         xlab = "mean counts over all samples", ylab = expression(log[2]~(Fold~change)),
         col = uithof_color[1])
  abline(h = 0, lty = 1, lwd = 2, col = uithof_color[30]) 
  
  maplot(res_qcdata, threshold = 0.1, labelsign = TRUE, main = "After sample QC", bty = "n",
         xlab = "mean counts over all samples", ylab = expression(log[2]~(Fold~change)),
         col = uithof_color[4])
  abline(h = 0, lty = 1, lwd = 2, col = uithof_color[30]) 
  
  maplot(res_qcfdata, threshold = 0.1, labelsign = TRUE, main = "After sample & gene-count QC", bty = "n",
         xlab = "mean counts over all samples", ylab = expression(log[2]~(Fold~change)),
         col = uithof_color[16])
  abline(h = 0, lty = 1, lwd = 2, col = uithof_color[30]) 
  
  mtext("MA-plots", outer = TRUE, cex = 1.25)

dev.off()
par(mfrow = c(1,1), mar = c(5, 4, 4, 2), oma = c(0, 0, 0, 0))

cat("\n  - Volcano plot with 'significant' genes labeled...")
# function to make volcanoplots
volcanoplot <- function(res, lfcthresh = 2, sigthresh = 0.05, 
                         main = "Volcano Plot", legendpos = "topleft", labelsig = TRUE, textcx = 1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch = 20, main = main, ...))
  with(subset(res, padj < sigthresh ), 
       points(log2FoldChange, -log10(pvalue), pch = 20, col = uithof_color[4], ...))
  with(subset(res, abs(log2FoldChange) > lfcthresh), 
       points(log2FoldChange, -log10(pvalue), pch = 20, col = uithof_color[1], ...))
  with(subset(res, padj < sigthresh & abs(log2FoldChange) > lfcthresh), 
       points(log2FoldChange, -log10(pvalue), pch = 20, col = uithof_color[20], ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj < sigthresh & abs(log2FoldChange) > lfcthresh), 
         textxy(log2FoldChange, -log10(pvalue), labs = symbol, cex = textcx, ...))
  }
  legend(legendpos, xjust = 1, yjust = 1, 
         legend = c(paste("FDR<", sigthresh, sep = ""), 
                    paste("|LogFC|>", lfcthresh, sep = ""), "both"), 
         pch = 20, col = c(uithof_color[4], uithof_color[1], uithof_color[20]), bty = "n")
}

png(file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,".VolcanoPlots.png"),
    width = 1920, height = 1080)
pdf(file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,".VolcanoPlots.pdf"), onefile = TRUE,
    paper = "a4r")
    
  par(mfrow = c(1,3), mar = c(5,8,4,2), oma = c(0, 0, 2, 0))
  
  volcanoplot(resdata, lfcthresh = 2, sigthresh = 0.05, textcx = 0.8,  main = "Raw data", bty = "n",
              xlab = expression(log[2]~(Fold~change)), ylab = expression(log[10]~italic((P)))
  )
  volcanoplot(res_qcdata, lfcthresh = 2, sigthresh = 0.05, textcx = 0.8,  main = "After sample QC", bty = "n",
              xlab = expression(log[2]~(Fold~change)), ylab = expression(log[10]~italic((P)))
  )
  volcanoplot(res_qcfdata, lfcthresh = 2, sigthresh = 0.05, textcx = 0.8,  main = "After sample & gene-count QC", bty = "n",
              xlab = expression(log[2]~Fold~change), ylab = expression(log[10]~italic((P)))
  )

    mtext("Volcano-plots", outer = TRUE, cex = 1.25)
dev.off()
par(mfrow = c(1,1), mar = c(5, 4, 4, 2), oma = c(0, 0, 0, 0))

cat("\n  - Top results...")
pdf(paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,".Top.pdf"))
  top.data <- resdata[order(resdata$pvalue), ][1,]
  top.symbol <- resdata[order(resdata$pvalue), "ENSEMBLID"][1]
  top <- ifelse(is.na(resdata[order(resdata$pvalue), "symbol"][1]) == TRUE, 
                "no_gene_name",
                resdata[order(resdata$pvalue), "symbol"][1])
  
  dds.cnts <- counts(estimateSizeFactors(dds))[top.symbol,]
  dds.status <- as.factor(colData(dds)[["CurSmoking"]])
  levels(dds.status) <- 1:length(levels(dds.status))
  dds.status <- as.numeric(dds.status)
  dds.pearson <- signif(cor(dds.cnts, dds.status), 3)
  
  par(mfrow = c(1,1), mar = c(5,8,4,2), oma = c(0, 0, 2, 0))
  boxplot(dds.cnts~dds.status, main = bquote(Pearson~r^2 == .(dds.pearson)),
          # sub = "(AEMS450K1)",
          cex.sub = 0.8,
          xlab = "Current smoker", 
          ylab = bquote(.(top.symbol)~" -"~italic(.(top))), 
          plot = TRUE, notch = FALSE, outline = TRUE, 
          names = c("no", "yes"), 
          whisklty = 1, staplelty = 0, # no whisker-ends, but with lines
          frame = FALSE)
  stripchart(dds.cnts~dds.status, vertical = TRUE,
             method = "jitter", add = TRUE, pch = 20, col = "#1290D9")
  
  mtext("Top correlated transcript (gene)", outer = TRUE, cex = 1.1)

dev.off()
rm(top.data, top.symbol, top,
   dds.cnts, dds.status, dds.pearson)
par(mfrow = c(1,1), mar = c(5, 4, 4, 2), oma = c(0, 0, 0, 0))

pdf(paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,".Top.QC.pdf"))
  top.data <- res_qcdata[order(res_qcdata$pvalue), ][1,]
  top.symbol <- res_qcdata[order(res_qcdata$pvalue), "ENSEMBLID"][1]
  top <- ifelse(is.na(res_qcdata[order(res_qcdata$pvalue), "symbol"][1]) == TRUE, 
                "no_gene_name",
                res_qcdata[order(res_qcdata$pvalue), "symbol"][1])
  
  dds.cnts <- counts(estimateSizeFactors(dds))[top.symbol,]
  dds.status <- as.factor(colData(dds)[["CurSmoking"]])
  levels(dds.status) <- 1:length(levels(dds.status))
  dds.status <- as.numeric(dds.status)
  dds.pearson <- signif(cor(dds.cnts, dds.status), 3)
  
  par(mfrow = c(1,1), mar = c(5,8,4,2), oma = c(0, 0, 2, 0))
  boxplot(dds.cnts~dds.status, main = bquote(Pearson~r^2 == .(dds.pearson)),
          # sub = "(AEMS450K1)",
          cex.sub = 0.8,
          xlab = "Current smoker", 
          ylab = bquote(.(top.symbol)~" -"~italic(.(top))), 
          plot = TRUE, notch = FALSE, outline = TRUE, 
          names = c("no", "yes"), 
          whisklty = 1, staplelty = 0, # no whisker-ends, but with lines
          frame = FALSE)
  stripchart(dds.cnts~dds.status, vertical = TRUE,
             method = "jitter", add = TRUE, pch = 20, col = "#1290D9")
  
  mtext("Top correlated transcript (gene)", outer = TRUE, cex = 1.1)

dev.off()
rm(top.data, top.symbol, top,
   dds.cnts, dds.status, dds.pearson)
par(mfrow = c(1,1), mar = c(5, 4, 4, 2), oma = c(0, 0, 0, 0))

pdf(paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,".Top.QCF.pdf"))
  top.data <- res_qcfdata[order(res_qcfdata$pvalue), ][1,]
  top.symbol <- res_qcfdata[order(res_qcfdata$pvalue), "ENSEMBLID"][1]
  top <- ifelse(is.na(res_qcfdata[order(res_qcfdata$pvalue), "symbol"][1]) == TRUE, 
                "no_gene_name",
                res_qcfdata[order(res_qcfdata$pvalue), "symbol"][1])
  
  dds.cnts <- counts(estimateSizeFactors(dds))[top.symbol,]
  dds.status <- as.factor(colData(dds)[["CurSmoking"]])
  levels(dds.status) <- 1:length(levels(dds.status))
  dds.status <- as.numeric(dds.status)
  dds.pearson <- signif(cor(dds.cnts, dds.status), 3)
  
  par(mfrow = c(1,1), mar = c(5,8,4,2), oma = c(0, 0, 2, 0))
  boxplot(dds.cnts~dds.status, main = bquote(Pearson~r^2 == .(dds.pearson)),
          # sub = "(AEMS450K1)",
          cex.sub = 0.8,
          xlab = "Current smoker", 
          ylab = bquote(.(top.symbol)~" -"~italic(.(top))), 
          plot = TRUE, notch = FALSE, outline = TRUE, 
          names = c("no", "yes"), 
          whisklty = 1, staplelty = 0, # no whisker-ends, but with lines
          frame = FALSE)
  stripchart(dds.cnts~dds.status, vertical = TRUE,
             method = "jitter", add = TRUE, pch = 20, col = "#1290D9")
  
  mtext("Top correlated transcript (gene)", outer = TRUE, cex = 1.1)

dev.off()
rm(top.data, top.symbol, top,
   dds.cnts, dds.status, dds.pearson)
par(mfrow = c(1,1), mar = c(5, 4, 4, 2), oma = c(0, 0, 0, 0))

cat("\n===========================================================================================")
cat("\n*** saving UPDATE data in between steps ***")
save(dds, file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,".dds.RData"))
save(dds_qc, file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,".dds_qc.RData"))
save(dds_qcf, file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,".dds_qcf.RData"))

save(se, file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,".se.RData"))
save(se_f, file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,".se_f.RData"))
save(se_fqc, file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,".se_fqc.RData"))

save(res, file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,".res.RData"))
save(res_qc, file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,".res_qc.RData"))
save(res_qcf, file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,".res_qcf.RData"))

save(resdata, file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,".resdata.RData"))
save(res_qcdata, file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,".res_qcdata.RData"))
save(res_qcfdata, file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,".res_qcfdata.RData"))

save(sampleTable, file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,".sampleTable.RData"))
save(sampleTableFiltered, file = paste0(RNAOUT_loc,"/",Today,".aerna.pilot.",EWAS_trait,".sampleTableFiltered.RData"))

cat("\n===========================================================================================")
cat("SAVE THE DATA")

save.image(paste0(ANALYSIS_loc,"/",Today,".aems450k.meta.analysis.ewas.",EWAS_trait,"_vs_RNA.RData"))

