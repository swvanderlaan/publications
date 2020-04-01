#!/bin/bash
#

#$ -S /bin/bash 																			# the type of BASH you'd like to use
#$ -N create.aegs.dataset  																	# the name of this script
# -hold_jid some_other_basic_bash_script  													# the current script (basic_bash_script) will hold until some_other_basic_bash_script has finished
#$ -o /hpc/dhl_ec/svanderlaan/projects/polygenicscores/create.aegs.dataset.log  								# the log file of this job
#$ -e /hpc/dhl_ec/svanderlaan/projects/polygenicscores/create.aegs.dataset.errors 							# the error file of this job
#$ -l h_rt=00:15:00  																		# h_rt=[max time, e.g. 02:02:01] - this is the time you think the script will take
#$ -l h_vmem=8G  																			#  h_vmem=[max. mem, e.g. 45G] - this is the amount of memory you think your script will use
# -l tmpspace=64G  																		# this is the amount of temporary space you think your script will use
#$ -M s.w.vanderlaan-2@umcutrecht.nl  														# you can send yourself emails when the job is done; "-M" and "-m" go hand in hand
#$ -m ea  																					# you can choose: b=begin of job; e=end of job; a=abort of job; s=suspended job; n=no mail is send
#$ -cwd  																					# set the job start to the current directory - so all the things in this script are relative to the current directory!!!
#
### INTERACTIVE SHELLS
# You can also schedule an interactive shell, e.g.:
#
# qlogin -N "basic_bash_script" -l h_rt=02:00:00 -l h_vmem=24G -M s.w.vanderlaan-2@umcutrecht.nl -m ea
#
# You can use the variables above (indicated by "#$") to set some things for the submission system.
# Another useful tip: you can set a job to run after another has finished. Name the job 
# with "-N SOMENAME" and hold the other job with -hold_jid SOMENAME". 
# Further instructions: https://wiki.bioinformatics.umcutrecht.nl/bin/view/HPC/HowToS#Run_a_job_after_your_other_jobs
#
# It is good practice to properly name and annotate your script for future reference for
# yourself and others. Trust me, you'll forget why and how you made this!!!

### Creating display functions
### Setting colouring
NONE='\033[00m'
BOLD='\033[1m'
OPAQUE='\033[2m'
FLASHING='\033[5m'
UNDERLINE='\033[4m'

RED='\033[01;31m'
GREEN='\033[01;32m'
YELLOW='\033[01;33m'
PURPLE='\033[01;35m'
CYAN='\033[01;36m'
WHITE='\033[01;37m'
### Regarding changing the 'type' of the things printed with 'echo'
### Refer to: 
### - http://askubuntu.com/questions/528928/how-to-do-underline-bold-italic-strikethrough-color-background-and-size-i
### - http://misc.flogisoft.com/bash/tip_colors_and_formatting
### - http://unix.stackexchange.com/questions/37260/change-font-in-echo-command

### echo -e "\033[1mbold\033[0m"
### echo -e "\033[3mitalic\033[0m" ### THIS DOESN'T WORK ON MAC!
### echo -e "\033[4munderline\033[0m"
### echo -e "\033[9mstrikethrough\033[0m"
### echo -e "\033[31mHello World\033[0m"
### echo -e "\x1B[31mHello World\033[0m"

# for i in $(seq 0 5) 7 8 $(seq 30 37) $(seq 41 47) $(seq 90 97) $(seq 100 107) ; do 
# 	echo -e "\033["$i"mYou can change the font...\033[0m"; 
# done
### Creating some function
# function echobold { #'echobold' is the function name
#     echo -e "\033[1m${1}\033[0m" # this is whatever the function needs to execute.
# }
function echocyan { #'echobold' is the function name
    echo -e "${CYAN}${1}${NONE}" # this is whatever the function needs to execute.
}
function echobold { #'echobold' is the function name
    echo -e "${BOLD}${1}${NONE}" # this is whatever the function needs to execute, note ${1} is the text for echo
}
function echoitalic { #'echobold' is the function name
    echo -e "\033[3m${1}\033[0m" # this is whatever the function needs to execute.
}

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echobold "                            CREATE AEGS DATASET"
echo ""
echoitalic "* Written by  : Sander W. van der Laan"
echoitalic "* E-mail      : s.w.vanderlaan-2@umcutrecht.nl"
echoitalic "* Last update : 2020-03-31"
echoitalic "* Version     : v1.0.1"
echo ""
echoitalic "* Description : This script will set some directories, and will then submit a "
echoitalic "                job to filter AEGS data."
echo ""
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Today's: "$(date)
TODAY=$(date +"%Y%m%d")
echo ""
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echobold "The following directories are set."
AEDATADIR="/hpc/dhl_ec/data/_ae_originals"
ORIGINALDIR="${AEDATADIR}/AEGS_COMBINED_IMPUTE2_1000Gp3_GoNL5"
PROJECTDIR="/hpc/dhl_ec/svanderlaan/projects/polygenicscores"
SOFTWARE="/hpc/local/CentOS7/dhl_ec/software"
QCTOOL="${SOFTWARE}/qctool_v204"
SNPTEST="${SOFTWARE}/snptest_v2.5.4beta3"

echo "Original data directory____ ${ORIGINALDIR}"
echo "Project directory__________ ${PROJECTDIR}"
echo "Dosage directory___________ ${DOSAGES}"
echo "Software directory_________ ${SOFTWARE}"
echo "Where \"qctool\" resides_____ ${QCTOOL}"
echo ""

### Make directories for script if they do not exist yet (!!!PREREQUISITE!!!)
if [ ! -d ${PROJECTDIR}/AEGSCOMBO_IMPUTE2_1000Gp3_GoNL5_1439/ ]; then
	mkdir -v ${PROJECTDIR}/AEGSCOMBO_IMPUTE2_1000Gp3_GoNL5_1439/
fi
PROJECTDATADIR="${PROJECTDIR}/AEGSCOMBO_IMPUTE2_1000Gp3_GoNL5_1439"

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echobold "Extracting some data in a for loop example and submitting this as a job."
# for CHR in $(seq 1 22) X; do 
# 	echo "* processing chromosome [ ${CHR} ]..."
# 	echo "$QCTOOL -g $ORIGINALDIR/aegs_combo_1kGp3GoNL5_RAW_chr${CHR}.vcf.gz -s $AEDATADIR/pheno_cov_exclusions/aegscombo_phenocov.sample -excl-samples $PROJECTDIR/aegscombo.exclusion_nonCEA.list -og $PROJECTDATADIR/aegs_combo_1kGp3GoNL5_1439pts_InouyeVariants_chr${CHR}.vcf.gz -os $PROJECTDATADIR/aegs_combo_1kGp3GoNL5_1439pts_InouyeVariants_chr${CHR}.sample -incl-rsids $PROJECTDIR/Inouye_bioRxiv_2018/metaGRS_hg19_20180205.4PRSICE.listvariants.txt " > $PROJECTDATADIR/create.aegs.dataset.chr${CHR}.sh
# 	echo "$QCTOOL -g $PROJECTDATADIR/aegs_combo_1kGp3GoNL5_1439pts_InouyeVariants_chr${CHR}.vcf.gz -s $PROJECTDATADIR/aegs_combo_1kGp3GoNL5_1439pts_InouyeVariants_chr${CHR}.sample -snp-stats -osnp $PROJECTDATADIR/aegs_combo_1kGp3GoNL5_1439pts_InouyeVariants_chr${CHR}.sumstats.txt " >> $PROJECTDATADIR/create.aegs.dataset.chr${CHR}.sh
# 	qsub -S /bin/bash -e $PROJECTDATADIR/create.aegs.dataset.chr${CHR}.errors -o $PROJECTDATADIR/create.aegs.dataset.chr${CHR}.output -l h_rt=02:00:00 -l h_vmem=8G -M s.w.vanderlaan-2@umcutrecht.nl -m ea -cwd $PROJECTDATADIR/create.aegs.dataset.chr${CHR}.sh
# done

for CHR in $(seq 1 22) X; do 
	echo "* processing chromosome [ ${CHR} ]..."
	echo "$SNPTEST -data $PROJECTDATADIR/aegs_combo_1kGp3GoNL5_1439pts_InouyeVariants_chr${CHR}.vcf.gz $PROJECTDATADIR/aegs_combo_1kGp3GoNL5_1439pts_InouyeVariants_chr${CHR}.sample -summary_stats_only -o $PROJECTDATADIR/aegs_combo_1kGp3GoNL5_1439pts_InouyeVariants_chr${CHR}.allstats.txt -log $PROJECTDATADIR/aegs_combo_1kGp3GoNL5_1439pts_InouyeVariants_chr${CHR}.allstats.log -hwe" >> $PROJECTDATADIR/stats.aegs.dataset.chr${CHR}.sh
	qsub -S /bin/bash -N stats.aegs.dataset.chr${CHR} -e $PROJECTDATADIR/stats.aegs.dataset.chr${CHR}.errors -o $PROJECTDATADIR/stats.aegs.dataset.chr${CHR}.output -l h_rt=00:30:00 -l h_vmem=8G -M s.w.vanderlaan-2@umcutrecht.nl -m ea -cwd $PROJECTDATADIR/stats.aegs.dataset.chr${CHR}.sh
done

echo ""
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echobold "Wow. I'm all done buddy. What a job! let's have a beer!"
date

