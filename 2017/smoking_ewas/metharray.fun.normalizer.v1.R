#!/hpc/local/CentOS7/dhl_ec/software/R-3.3.3/bin/Rscript --vanilla

# Alternative shebang for local Mac OS X: "#!/usr/local/bin/Rscript --vanilla"
# Linux version for HPC: #!/hpc/local/CentOS7/dhl_ec/software/R-3.3.1/bin/Rscript --vanilla
cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    METHYLATION ARRAY FUNCTIONAL NORMALIZATION -- FOR ATHERO-EXPRESS METHYLATION STUDY 1 & 2\n
    * Version: v1.2
    * Last edit: 2017-04-27
    * Created by: Sander W. van der Laan | s.w.vanderlaan-2@umcutrecht.nl;
                  Annelot M. Dekker | a.m.dekker-7@umcutrecht.nl\n
    * Description: Load iDatfiles in parallel, remove MethylAid outliers and save RGset/snpprobes. Based on 
                   DNAmArray:
                   - https://molepi.github.io/DNAmArray_workflow/index.html
                   - https://github.com/molepi/DNAmArray\n
    The script should be usuable on both any Linux distribution with R 3+ installed, Mac OS X and Windows.
    
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
cat("--------------------------------------------------------------------------------------------------------------\n")
cat("CLEAR THE BOARD")
rm(list = ls())

cat("\n--------------------------------------------------------------------------------------------------------------\n")
cat("GENERAL R SETUP")
### FUNCTION TO INSTALL PACKAGES
### This function will automatically check in both CRAN and Bioconductor. This is 
### a function found by Sander W. van der Laan online from @Samir: 
### http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
### 
cat("\n* Creating functions to... ")
cat("\n  - install and load packages...")
install.packages.auto <- function(x) { 
  x <- as.character(substitute(x)) 
  if(isTRUE(x %in% .packages(all.available = TRUE))) { 
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else { 
    # Update installed packages - this may mean a full upgrade of R, which in turn
    # may not be warrented. 
    #update.packages(ask = FALSE) 
    eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE, repos = \"http://cran-mirror.cs.uu.nl/\")", x)))
  }
  if(isTRUE(x %in% .packages(all.available = TRUE))) { 
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

cat("\n  - load RData file to a specific named variable...")
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

### In this case I'm keeping track of the various packages, as versions and 
### actual loading of the libraries gave issues before.

cat("\n* Loading packages (or installing if necessary)...\n")
cat("- General packages...\n")
suppressMessages(install.packages.auto("optparse"))
suppressMessages(install.packages.auto("tools"))
suppressMessages(install.packages.auto("openxlsx"))
suppressMessages(install.packages.auto("devtools"))
suppressMessages(install.packages.auto("pheatmap"))
suppressMessages(install.packages.auto("irlba"))
suppressMessages(devtools::install_github("rstudio/shiny"))
suppressMessages(install.packages.auto('rsconnect'))
suppressMessages(install.packages.auto("lattice"))
suppressMessages(install.packages.auto("pheatmap"))
suppressMessages(install.packages.auto("fastcluster"))

cat("\n- DNAmArray package...\n")
# Also refer to: 
# - https://molepi.github.io/DNAmArray_workflow/index.html
# - https://github.com/molepi/DNAmArray
suppressMessages(install_github("molepi/DNAmArray"))

cat("\n- MethylAid/minfi packages...\n")
suppressMessages(install.packages.auto("shinyMethyl"))
suppressMessages(install.packages.auto("minfi"))
suppressMessages(install.packages.auto("minfiData"))
suppressMessages(install.packages.auto("MethylAid"))
suppressMessages(install.packages.auto("MethylAidData"))
suppressMessages(install.packages.auto("illuminaio"))

cat("\n- Parallelisation packages...\n")
suppressMessages(install.packages.auto("BiocParallel"))

cat("\n- Genomics packages...\n")
suppressMessages(install.packages.auto("Biobase"))
suppressMessages(install.packages.auto("BiocGenerics"))

cat("\n- Load libraries...\n")
suppressMessages(library("BiocParallel"))
suppressMessages(library("Biobase"))
suppressMessages(library("BiocGenerics"))
suppressMessages(library("minfi"))
suppressMessages(library("illuminaio"))
suppressMessages(library("DNAmArray"))

cat("\n- Create datestamp...\n") 
Today = format(as.Date(as.POSIXlt(Sys.time())), "%Y%m%d")
Today.Report = format(as.Date(as.POSIXlt(Sys.time())), "%A, %B %d, %Y")

### -----------------------------------------------------------------------------------
###                        UTRECHT SCIENCE PARK COLOURS SCHEME
### -----------------------------------------------------------------------------------
### Website to convert HEX to RGB: http://hex.colorrrs.com.
### For some functions you should divide these numbers by 255.
###
### -----------------------------------------------------------------------------------
### Color.............HEX.....RGB...............CHR...MAF/INFO
### -----------------------------------------------------------------------------------
### yellow............#FBB820 (251,184,32)	=>	1.....or 1.0 > INFO
###	gold..............#F59D10 (245,157,16)	=>	2		
###	salmon............#E55738 (229,87,56) 	=>	3.....or 0.05 < MAF < 0.2 or 0.4 < INFO < 0.6
###	darkpink..........#DB003F ((219,0,63)		=>	4		
###	lightpink.........#E35493 (227,84,147)	=>	5.....or 0.8 < INFO < 1.0
###	pink..............#D5267B (213,38,123)	=>	6		
###	hardpink..........#CC0071 (204,0,113)		=>	7		
###	lightpurple.......#A8448A (168,68,138)	=>	8		
###	purple............#9A3480 (154,52,128)	=>	9		
###	lavendel..........#8D5B9A (141,91,154)	=>	10		
###	bluepurple........#705296 (112,82,150)	=>	11		
###	purpleblue........#686AA9 (104,106,169)	=>	12		
###	lightpurpleblue...#6173AD (97,115,173)	=>	13		
###	seablue...........#4C81BF (76,129,191)	=>	14		
###	skyblue...........#2F8BC9 (47,139,201)	=>	15		
###	azurblue..........#1290D9 (18,144,217)	=>	16.....or 0.01 < MAF < 0.05 or 0.2 < INFO < 0.4
###	lightazurblue.....#1396D8 (19,150,216)	=>	17		
###	greenblue.........#15A6C1 (21,166,193)	=>	18		
###	seaweedgreen......#5EB17F (94,177,127)	=>	19		
###	yellowgreen.......#86B833 (134,184,51)	=>	20		
###	lightmossgreen....#C5D220 (197,210,32)	=>	21		
###	mossgreen.........#9FC228 (159,194,40)	=>	22.....or MAF > 0.20 or 0.6 < INFO < 0.8
###	lightgreen........#78B113 (120,177,19)	=>	23/X
###	green.............#49A01D (73,160,29)		=>	24/Y
###	grey..............#595A5C (89,90,92)		=>	25/XY..or MAF < 0.01 or 0.0 < INFO < 0.2
###	lightgrey.........#A2A3A4	(162,163,164)	=> 	26/MT
uithof_color=c("#FBB820","#F59D10","#E55738","#DB003F","#E35493","#D5267B",
               "#CC0071","#A8448A","#9A3480","#8D5B9A","#705296","#686AA9",
               "#6173AD","#4C81BF","#2F8BC9","#1290D9","#1396D8","#15A6C1",
               "#5EB17F","#86B833","#C5D220","#9FC228","#78B113","#49A01D",
               "#595A5C","#A2A3A4")
###	
### COLOR DEFINITION SAMPLE GROUPS
### 
### AE              azurblue  #1290D9
### TWIECE          hardpink  #9A3480
### HUVECS          gold      #F59D10
### CMPC            darkpink  #E55738
### SXS_Controls    green     #49A01D
###	
uithof_color_sample_groups=c("#1290D9", "#9A3480", "#F59D10", "#E55738", "#49A01D")
###	
### -----------------------------------------------------------------------------------

#--------------------------------------------------------------------------
### OPTION LISTING
option_list = list(
  make_option(c("-p", "--projectdir"), action="store", default=NA, type='character',
              help="Path to the project directory."),
  make_option(c("-i", "--inputdir"), action="store", default=NA, type='character',
              help="Path to the input directory."),
  make_option(c("-o", "--outputdir"), action="store", default=NA, type='character',
              help="Path to the output directory."),
  make_option(c("-r", "--rgsetfile"), action="store", default=NA, type='character',
              help="The RGset Rdata-file ('RGsetQC.Rdata'); contains the list of iDAT-file locations. Relative to the project directory."),
  make_option(c("-f", "--filename"), action="store", default=NA, type='character',
              help="The generic filename used to name the genotype Rdata-file ('snpprobes.Rdata', which contains the 65 SNP probes for genotype checks), and the extended RGChannelSet Rdata-file after removal of MethylAid outliers ('RGsetQC.Rdata'). Relative to the project directory."),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Should the program print extra stuff out? [default %default]"),
  make_option(c("-s", "--silent"), action="store_false", dest="verbose",
              help="Make the program not be verbose.")
)
opt = parse_args(OptionParser(option_list=option_list))
if (opt$verbose) {
  # You can use either the long or short name; so opt$a and opt$avar are the same.
  # Show the user what the variables are.
  cat("--------------------------------------------------------------------------------------------------------------")
  cat("\n* Verbose reporting. Checking the settings as given through the flags.")
  cat(paste0("\nThe project directory.............: ", opt$projectdir,"."))
  cat(paste0("\nThe input directory...............: ", opt$inputdir,"."))
  cat(paste0("\nThe output directory..............: ", opt$outputdir,"."))
  cat(paste0("\nThe RGset Rdata-file..............: ", opt$rgsetfile,"."))
  cat(paste0("\nThe generic filename..............: ", opt$filename,"."))
  cat("\n--------------------------------------------------------------------------------------------------------------")
  cat("\n\n")
}
cat("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
cat("Starting \"Functional Normalization\".")

if(!is.na(opt$projectdir) & !is.na(opt$inputdir) & !is.na(opt$outputdir) & !is.na(opt$rgsetfile) & !is.na(opt$filename)) {
  cat("\n* Setting directories.")
  
  cat("\n  - setting project directories and variables...")
  ROOT_loc = paste0(opt$projectdir)
  INP_loc = paste0(ROOT_loc, "/", opt$inputdir)
  OUT_loc = paste0(INP_loc, "/", opt$outputdir)
  
  QC_loc = paste0(INP_loc,"/QC") # already made by `metharray.targets.creator.v1.R`
  
  PLOT_loc = paste0(QC_loc,"/Plots") # already made by `metharray.targets.creator.v1.R`
  
  load_rgset <- paste0(QC_loc, "/", opt$rgsetfile)
  outfile_GRset <- paste0(QC_loc, "/",opt$filename,".GRChannelSetQC")

  cat(paste("\nLoad RGChannelSet for functional normalization. 
  Analysing these results..............................: [ ", load_rgset, " ].
  File names will be ..................................: [ ", opt$filename,".* ].
  Functional normalized results will be saved here.....: [ ", OUT_loc," ].\n"))
  
  cat(paste0("Today's: ",Today.Report,".\n"))
  
  cat("\n* Reading in RGChannelSet...\n")

  load(load_rgset)
  cat("\nDone! RGChannelSet loaded.\n")
  
  cat("\n* Functional normalization...\n")
  GRsetQC <- preprocessFunnorm.DNAmArray(RGsetQC, nPCs = 4, keepCN = FALSE)
  
  cat(paste("  > There are", dim(GRsetQC)[[2]], "samples in the GRsetQC after functional normalization.\n", sep = " "))
  cat("\nDone! GRChannelSet created.\n")
  
  cat("\n* Saving extended RGChannelSet...\n")
  # Save the (extended) RGset
  save(GRsetQC, file = paste0(outfile_GRset,".Rdata"))
  
  cat(paste("\nDone! Saved output to:",outfile_GRset,".Rdata\n", sep = ""))

} else {

cat("*** ERROR *** 
You didn't specify all variables:
- --p/projectdir      : path to the project directory
- --i/inputdir        : path to the input directory
- --o/outputdir       : path to the output directory
- --r/rgsetfile      : the target Rdata-file
- --f/filename        : the generic filename\n\n", 
file=stderr()) # print error messages to stderr
}

#--------------------------------------------------------------------------
### CLOSING MESSAGE
cat("\n Wow! I'm all done. That was a lot of work! \n")
cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

#--------------------------------------------------------------------------
### SAVE ENVIRONMENT | FOR DEBUGGING
#save.image(paste0(OUT_loc, "/",Today,".metharray.fun.normalizer.RData"))
