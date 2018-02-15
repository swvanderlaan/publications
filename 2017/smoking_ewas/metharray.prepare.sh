#!/bin/bash
#
#$ -S /bin/bash 																			# the type of BASH you'd like to use
#$ -N metharray.prepare 																	# the name of this script
# -hold_jid some_other_basic_bash_script  													# the current script (basic_bash_script) will hold until some_other_basic_bash_script has finished
#$ -o /hpc/dhl_ec/data/_ae_originals/AEMS450KCombo/metharray.getsnps.log  								# the log file of this job
#$ -e /hpc/dhl_ec/data/_ae_originals/AEMS450KCombo/metharray.getsnps.errors 							# the error file of this job
#$ -l h_rt=00:15:00  																		# h_rt=[max time, e.g. 02:02:01] - this is the time you think the script will take
#$ -l h_vmem=4G  																			#  h_vmem=[max. mem, e.g. 45G] - this is the amount of memory you think your script will use
# -l tmpspace=64G  																		# this is the amount of temporary space you think your script will use
#$ -M s.w.vanderlaan-2@umcutrecht.nl  														# you can send yourself emails when the job is done; "-M" and "-m" go hand in hand
#$ -m ea  																					# you can choose: b=begin of job; e=end of job; a=abort of job; s=suspended job; n=no mail is send
#$ -cwd  																					# set the job start to the current directory - so all the things in this script are relative to the current directory!!!
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
### NONE='\033[00m'
### BOLD='\033[1m'
### OPAQUE='\033[2m'
### FLASHING='\033[5m'
### UNDERLINE='\033[4m'
### 
### RED='\033[01;31m'
### GREEN='\033[01;32m'
### YELLOW='\033[01;33m'
### PURPLE='\033[01;35m'
### CYAN='\033[01;36m'
### WHITE='\033[01;37m'
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

#for i in $(seq 0 5) 7 8 $(seq 30 37) $(seq 41 47) $(seq 90 97) $(seq 100 107) ; do 
#	echo -e "\033["$i"mYou can change the font...\033[0m"; 
#done
### Creating some function
function echobold { #'echobold' is the function name
    echo -e "\033[1m${1}\033[0m" # this is whatever the function needs to execute.
}
function echobold { #'echobold' is the function name
    echo -e "${BOLD}${1}${NONE}" # this is whatever the function needs to execute, note ${1} is the text for echo
}
function echoitalic { #'echobold' is the function name
    echo -e "\033[3m${1}\033[0m" # this is whatever the function needs to execute.
}

echobold "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echobold "                           METHYLATION ARRAY PREPARATION"
echo "                                      version 1.0"
echo ""
echoitalic "* Written by  : Sander W. van der Laan"
echoitalic "* E-mail      : s.w.vanderlaan-2@umcutrecht.nl"
echoitalic "* Last update : 2017-07-27"
echoitalic "* Version     : prepare.methylation.array.sh"
echo ""
echoitalic "* Description : This script will prepare methylation array data for downstream "
echoitalic "                analyses."
echo ""
echobold "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Today's: "$(date)
TODAY=$(date +"%Y%m%d")
echo ""
echobold "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echobold "The following directories are set."
SOFTWARE=/hpc/local/CentOS7/dhl_ec/software
ORIGINALS=/hpc/dhl_ec/data/_ae_originals
AEMS450KCOMBO="AEMS450KCombo"
PROJECTNAME="aems450k"

### QSUB Settings R-scripts
TARGETTIME="00:15:00"
TARGETMEM="8G"
TARGETCORES="8"
SUMTIME="00:30:00"
SUMMEM="16G"
SUMCORES="8"
RGTIME="02:00:00"
RGMEM="256G"
RCORES="8"
GRTIME="02:00:00"
GRMEM="256G"
GRCORES="8"

echo "Original data directory..... ${ORIGINALS}"
echo "Project directory........... ${AEMS450KCOMBO}"
echo ""

### Make directories for script if they do not exist yet (!!!PREREQUISITE!!!)
#if [ ! -d ${ORIGINALS}/somedir/ ]; then
#  mkdir -v ${ORIGINALS}/somedir/
#fi
#RAWDIR=${ORIGINALS}/somedir

echobold "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

###COMBO VERSION
echobold "Getting target file."
echo "Rscript metharray.target.creator.v1.R -p ${ORIGINALS} -i AEMS450KCombo -r RAWDATA -o OUTPUT -f ${PROJECTNAME} -n 8" > ${ORIGINALS}/${AEMS450KCOMBO}/metharray.target.creator.v1.sh
echo "" >> ${ORIGINALS}/${AEMS450KCOMBO}/metharray.target.creator.v1.sh
qsub -S /bin/bash -N TARGET.MA.${PROJECTNAME} -e ${ORIGINALS}/${AEMS450KCOMBO}/metharray.target.creator.v1.errors -o ${ORIGINALS}/${AEMS450KCOMBO}/metharray.target.creator.v1.log -l h_rt=${TARGETTIME} -l h_vmem=${TARGETMEM} -pe threaded ${TARGETCORES} -M s.w.vanderlaan-2@umcutrecht.nl -m ea -cwd ${ORIGINALS}/${AEMS450KCOMBO}/metharray.target.creator.v1.sh

echobold "Summarizing methylation array data."
echo "Rscript metharray.summarizer.v1.R -p ${ORIGINALS} -i AEMS450KCombo -r RAWDATA -o OUTPUT -f ${PROJECTNAME} -n 8" > ${ORIGINALS}/${AEMS450KCOMBO}/metharray.summarizer.v1.sh
echo "" >> ${ORIGINALS}/${AEMS450KCOMBO}/metharray.summarizer.v1.sh
qsub -S /bin/bash -N SUM.MA.${PROJECTNAME} -hold_jid TARGET.MA.${PROJECTNAME} -e ${ORIGINALS}/${AEMS450KCOMBO}/metharray.summarizer.v1.errors -o ${ORIGINALS}/${AEMS450KCOMBO}/metharray.summarizer.v1.log -l h_rt=${SUMTIME} -l h_vmem=${SUMMEM} -pe threaded ${SUMCORES} -M s.w.vanderlaan-2@umcutrecht.nl -m ea -cwd ${ORIGINALS}/${AEMS450KCOMBO}/metharray.summarizer.v1.sh


### PER DATASET VERSION
#echobold "Getting target file."
#echo "Rscript metharray.target.creator.v1.R -p ${ORIGINALS} -i AEMS450KCombo/AEMS450K1 -r RAWDATA -o OUTPUT -f ${PROJECTNAME} -n 8" > ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k1.target.creator.v1.sh
#qsub -S /bin/bash -N TARGET.AEMS450K1.MA.${PROJECTNAME} -e ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k1.target.creator.v1.errors -o ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k1.target.creator.v1.log -l h_rt=${TARGETTIME} -l h_vmem=${TARGETMEM} -pe threaded ${TARGETCORES} -M s.w.vanderlaan-2@umcutrecht.nl -m ea -cwd ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k1.target.creator.v1.sh
#echo "Rscript metharray.target.creator.v1.R -p ${ORIGINALS} -i AEMS450KCombo/AEMS450K2 -r RAWDATA -o OUTPUT -f ${PROJECTNAME} -n 8" > ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k2.target.creator.v1.sh
#qsub -S /bin/bash -N TARGET.AEMS450K2.MA.${PROJECTNAME} -e ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k2.target.creator.v1.errors -o ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k2.target.creator.v1.log -l h_rt=${TARGETTIME} -l h_vmem=${TARGETMEM} -pe threaded ${TARGETCORES} -M s.w.vanderlaan-2@umcutrecht.nl -m ea -cwd ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k2.target.creator.v1.sh
#
#echobold "Summarizing methylation array data."
#echo "Rscript metharray.summarizer.v1.R -p ${ORIGINALS} -i AEMS450KCombo/AEMS450K1 -r RAWDATA -o OUTPUT -f ${PROJECTNAME} -n 8" > ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k1.summarizer.v1.sh
#qsub -S /bin/bash -N SUM.AEMS450K1.MA.${PROJECTNAME} -hold_jid TARGET.AEMS450K1.MA.${PROJECTNAME} -e ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k1.summarizer.v1.errors -o ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k1.summarizer.v1.log -l h_rt=${SUMTIME} -l h_vmem=${SUMMEM} -pe threaded ${SUMCORES} -M s.w.vanderlaan-2@umcutrecht.nl -m ea -cwd ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k1.summarizer.v1.sh
#echo "Rscript metharray.summarizer.v1.R -p ${ORIGINALS} -i AEMS450KCombo/AEMS450K2 -r RAWDATA -o OUTPUT -f ${PROJECTNAME} -n 8" > ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k2.summarizer.v1.sh
#qsub -S /bin/bash -N SUM.AEMS450K2.MA.${PROJECTNAME} -hold_jid TARGET.AEMS450K2.MA.${PROJECTNAME} -e ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k2.summarizer.v1.errors -o ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k2.summarizer.v1.log -l h_rt=${SUMTIME} -l h_vmem=${SUMMEM} -pe threaded ${SUMCORES} -M s.w.vanderlaan-2@umcutrecht.nl -m ea -cwd ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k2.summarizer.v1.sh

echobold "Interactively visualize and print outliers. For convenience: DO THIS LOCALLY!!!"
Rscript metharray.visualizer.v1.R -p ${ORIGINALS} -i AEMS450KCombo -r RAWDATA -o OUTPUT -f ${PROJECTNAME} -n 8


###COMBO VERSION
echobold "Creating extended RGChannelSet."
###Rscript metharray.rgset.creator.parallel.v1.R -p ${ORIGINALS} -i AEMS450KCombo -o OUTPUT -t ${PROJECTNAME}.targets.Rdata -r ${PROJECTNAME}.outliers.csv -f ${PROJECTNAME} -n ${RCORES}
echo "Rscript metharray.rgset.creator.parallel.v1.R -p ${ORIGINALS} -i AEMS450KCombo -o OUTPUT -t ${PROJECTNAME}.targets.Rdata -r ${PROJECTNAME}.outliers.csv -f ${PROJECTNAME} -n ${RCORES}" > ${ORIGINALS}/${AEMS450KCOMBO}/metharray.rgset.creator.parallel.v1.sh
echo "" >> ${ORIGINALS}/${AEMS450KCOMBO}/metharray.rgset.creator.parallel.v1.sh
qsub -S /bin/bash -N RGSET.MA.${PROJECTNAME} -hold_jid SUM.MA.${PROJECTNAME} -e ${ORIGINALS}/${AEMS450KCOMBO}/metharray.rgset.creator.parallel.v1.errors -o ${ORIGINALS}/${AEMS450KCOMBO}/metharray.rgset.creator.parallel.v1.log -l h_rt=${RGTIME} -l h_vmem=${RGMEM} -M s.w.vanderlaan-2@umcutrecht.nl -m ea -cwd ${ORIGINALS}/${AEMS450KCOMBO}/metharray.rgset.creator.parallel.v1.sh

echobold "Functional normalization of the extended RGChannelSet."
###Rscript metharray.fun.normalizer.v1.R -p ${ORIGINALS} -i AEMS450KCombo -o OUTPUT -r aems450k.RGChannelSetQC.Rdata -f ${PROJECTNAME}
echo "Rscript metharray.fun.normalizer.v1.R -p ${ORIGINALS} -i AEMS450KCombo -o OUTPUT -r aems450k.RGChannelSetQC.Rdata -f ${PROJECTNAME}" > ${ORIGINALS}/${AEMS450KCOMBO}/metharray.fun.normalizer.v1.sh
echo "" >> ${ORIGINALS}/${AEMS450KCOMBO}/metharray.fun.normalizer.v1.sh
qsub -S /bin/bash -N GRSET.MA.${PROJECTNAME} -hold_jid RGSET.MA.${PROJECTNAME} -e ${ORIGINALS}/${AEMS450KCOMBO}/metharray.fun.normalizer.v1.errors -o ${ORIGINALS}/${AEMS450KCOMBO}/metharray.fun.normalizer.v1.log -l h_rt=${GRTIME} -l h_vmem=${GRMEM} -M s.w.vanderlaan-2@umcutrecht.nl -m ea -cwd ${ORIGINALS}/${AEMS450KCOMBO}/metharray.fun.normalizer.v1.sh


### PER DATASET VERSION
#echobold "Creating extended RGChannelSet."
#echo "Rscript metharray.rgset.creator.parallel.v1.R -p ${ORIGINALS} -i AEMS450KCombo/AEMS450K1 -o OUTPUT -t ${PROJECTNAME}.targets.Rdata -r ${PROJECTNAME}.outliers.csv -f ${PROJECTNAME} -n ${RCORES}" > ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k1.rgset.creator.parallel.v1.sh
#qsub -S /bin/bash -N RGSET.AEMS450K1.MA.${PROJECTNAME} -hold_jid SUM.AEMS450K1.MA.${PROJECTNAME} -e ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k1.rgset.creator.parallel.v1.errors -o ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k1.rgset.creator.parallel.v1.log -l h_rt=${RGTIME} -l h_vmem=${RGMEM} -M s.w.vanderlaan-2@umcutrecht.nl -m ea -cwd ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k1.rgset.creator.parallel.v1.sh
#echo "Rscript metharray.rgset.creator.parallel.v1.R -p ${ORIGINALS} -i AEMS450KCombo/AEMS450K2 -o OUTPUT -t ${PROJECTNAME}.targets.Rdata -r ${PROJECTNAME}.outliers.csv -f ${PROJECTNAME} -n ${RCORES}" > ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k2.rgset.creator.parallel.v1.sh
#qsub -S /bin/bash -N RGSET.AEMS450K2.MA.${PROJECTNAME} -hold_jid SUM.AEMS450K2.MA.${PROJECTNAME} -e ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k2.rgset.creator.parallel.v1.errors -o ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k2.rgset.creator.parallel.v1.log -l h_rt=${RGTIME} -l h_vmem=${RGMEM} -M s.w.vanderlaan-2@umcutrecht.nl -m ea -cwd ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k2.rgset.creator.parallel.v1.sh


#echobold "Functional normalization of the extended RGChannelSet."
#echo "Rscript metharray.fun.normalizer.v1.R -p ${ORIGINALS} -i AEMS450KCombo/AEMS450K1 -o OUTPUT -r aems450k.RGChannelSetQC.Rdata -f ${PROJECTNAME}" > ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k1.fun.normalizer.v1.sh
#qsub -S /bin/bash -N GRSET.AEMS450K1.MA.${PROJECTNAME} -hold_jid RGSET.AEMS450K1.MA.${PROJECTNAME} -e ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k1.fun.normalizer.v1.errors -o ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k1.fun.normalizer.v1.log -l h_rt=${GRTIME} -l h_vmem=${GRMEM} -M s.w.vanderlaan-2@umcutrecht.nl -m ea -cwd ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k1.fun.normalizer.v1.sh
#echo "Rscript metharray.fun.normalizer.v1.R -p ${ORIGINALS} -i AEMS450KCombo/AEMS450K2 -o OUTPUT -r aems450k.RGChannelSetQC.Rdata -f ${PROJECTNAME}" > ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k2.fun.normalizer.v1.sh
#qsub -S /bin/bash -N GRSET.AEMS450K2.MA.${PROJECTNAME} -hold_jid RGSET.AEMS450K2.MA.${PROJECTNAME} -e ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k2.fun.normalizer.v1.errors -o ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k2.fun.normalizer.v1.log -l h_rt=${GRTIME} -l h_vmem=${GRMEM} -M s.w.vanderlaan-2@umcutrecht.nl -m ea -cwd ${ORIGINALS}/${AEMS450KCOMBO}/metharray.aems450k2.fun.normalizer.v1.sh

echobold "Getting relevant SNP data from Athero-Express."
# for CHR in $(seq 1 22) X; do 
# 	echo "* processing chromosome [ ${CHR} ] ..."
# 	echo "${SOFTWARE}/qctool_v1.5 -g ${ORIGINALS}/AEGS_COMBINED_IMPUTE2_1000Gp3_GoNL5/aegs_combo_1kGp3GoNL5_RAW_chr${CHR}.gen.gz -og ${ORIGINALS}/${AEMS450KCOMBO}/OUTPUT/aegs_combo_1kGp3GoNL5_RAW_chr${CHR}.gen -incl-rsids ${ORIGINALS}/${AEMS450KCOMBO}/aems450k.65SNPlist.txt " > ${ORIGINALS}/${AEMS450KCOMBO}/metharray.getgeneticdata.chr${CHR}.v1.sh
# 	qsub -S /bin/bash -N GETSNP.chr${CHR}.${PROJECTNAME} -e ${ORIGINALS}/${AEMS450KCOMBO}/metharray.getgeneticdata.chr${CHR}.v1.errors -o ${ORIGINALS}/${AEMS450KCOMBO}/metharray.getgeneticdata.chr${CHR}.v1.log -l h_rt=${TARGETTIME} -l h_vmem=${TARGETMEM} -M s.w.vanderlaan-2@umcutrecht.nl -m ea -cwd ${ORIGINALS}/${AEMS450KCOMBO}/metharray.getgeneticdata.chr${CHR}.v1.sh
# done
cat aegs_combo_1kGp3GoNL5_RAW_chr1.gen > aegs_combo_1kGp3GoNL5_RAW.gen
for CHR in $(seq 2 22) X; do 
	echo "processing chromosome [ ${CHR} ]"
	cat aegs_combo_1kGp3GoNL5_RAW_chr${CHR}.gen >> aegs_combo_1kGp3GoNL5_RAW.gen
done

echo ""
echobold "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echobold "Wow. I'm all done buddy. What a job! let's have a beer!"
date

###	UtrechtSciencePark Colours Scheme
###
### Website to convert HEX to RGB: http://hex.colorrrs.com.
### For some functions you should divide these numbers by 255.
###
###	No.	Color				HEX		RGB							CMYK					CHR		MAF/INFO
### --------------------------------------------------------------------------------------------------------------------
###	1	yellow				#FBB820 (251,184,32)				(0,26.69,87.25,1.57) 	=>	1 		or 1.0 > INFO
###	2	gold				#F59D10 (245,157,16)				(0,35.92,93.47,3.92) 	=>	2		
###	3	salmon				#E55738 (229,87,56) 				(0,62.01,75.55,10.2) 	=>	3 		or 0.05 < MAF < 0.2 or 0.4 < INFO < 0.6
###	4	darkpink			#DB003F ((219,0,63)					(0,100,71.23,14.12) 	=>	4		
###	5	lightpink			#E35493 (227,84,147)				(0,63,35.24,10.98) 	=>	5 		or 0.8 < INFO < 1.0
###	6	pink				#D5267B (213,38,123)				(0,82.16,42.25,16.47) 	=>	6		
###	7	hardpink			#CC0071 (204,0,113)					(0,0,0,0) 	=>	7		
###	8	lightpurple			#A8448A (168,68,138)				(0,0,0,0) 	=>	8		
###	9	purple				#9A3480 (154,52,128)				(0,0,0,0) 	=>	9		
###	10	lavendel			#8D5B9A (141,91,154)				(0,0,0,0) 	=>	10		
###	11	bluepurple			#705296 (112,82,150)				(0,0,0,0) 	=>	11		
###	12	purpleblue			#686AA9 (104,106,169)				(0,0,0,0) 	=>	12		
###	13	lightpurpleblue		#6173AD (97,115,173/101,120,180)	(0,0,0,0) 	=>	13		
###	14	seablue				#4C81BF (76,129,191)				(0,0,0,0) 	=>	14		
###	15	skyblue				#2F8BC9 (47,139,201)				(0,0,0,0) 	=>	15		
###	16	azurblue			#1290D9 (18,144,217)				(0,0,0,0) 	=>	16		 or 0.01 < MAF < 0.05 or 0.2 < INFO < 0.4
###	17	lightazurblue		#1396D8 (19,150,216)				(0,0,0,0) 	=>	17		
###	18	greenblue			#15A6C1 (21,166,193)				(0,0,0,0) 	=>	18		
###	19	seaweedgreen		#5EB17F (94,177,127)				(0,0,0,0) 	=>	19		
###	20	yellowgreen			#86B833 (134,184,51)				(0,0,0,0) 	=>	20		
###	21	lightmossgreen		#C5D220 (197,210,32)				(0,0,0,0) 	=>	21		
###	22	mossgreen			#9FC228 (159,194,40)				(0,0,0,0) 	=>	22		or MAF > 0.20 or 0.6 < INFO < 0.8
###	23	lightgreen			#78B113 (120,177,19)				(0,0,0,0) 	=>	23/X
###	24	green				#49A01D (73,160,29)					(0,0,0,0) 	=>	24/Y
###	25	grey				#595A5C (89,90,92)					(0,0,0,0) 	=>	25/XY	or MAF < 0.01 or 0.0 < INFO < 0.2
###	26	lightgrey			#A2A3A4	(162,163,164)				(0,0,0,0) 	=> 	26/MT
### --------------------------------------------------------------------------------------------------------------------



