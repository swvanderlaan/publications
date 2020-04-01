#!/bin/bash

# Created by		Sander W. van der Laan | UMC Utrecht | s.w.vanderlaan[at]gmail[dot]com
# Last edit			2018-10-09
# Version			1.3.0

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

function echocyan { #'echobold' is the function name
    echo -e "${CYAN}${1}${NONE}" # this is whatever the function needs to execute.
}
function echobold { #'echobold' is the function name
    echo -e "${BOLD}${1}${NONE}" # this is whatever the function needs to execute, note ${1} is the text for echo
}
function echoitalic { #'echobold' is the function name
    echo -e "\033[3m${1}\033[0m" # this is whatever the function needs to execute.
}

script_copyright_message() {
	echo ""
	THISYEAR=$(date +'%Y')
	echoitalic "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echoitalic "+ The MIT License (MIT)                                                                                 +"
	echoitalic "+ Copyright (c) 1979-${THISYEAR} Sander W. van der Laan                                                        +"
	echoitalic "+                                                                                                       +"
	echoitalic "+ Permission is hereby granted, free of charge, to any person obtaining a copy of this software and     +"
	echoitalic "+ associated documentation files (the \"Software\"), to deal in the Software without restriction,         +"
	echoitalic "+ including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, +"
	echoitalic "+ and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, +"
	echoitalic "+ subject to the following conditions:                                                                  +"
	echoitalic "+                                                                                                       +"
	echoitalic "+ The above copyright notice and this permission notice shall be included in all copies or substantial  +"
	echoitalic "+ portions of the Software.                                                                             +"
	echoitalic "+                                                                                                       +"
	echoitalic "+ THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT     +"
	echoitalic "+ NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND                +"
	echoitalic "+ NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES  +"
	echoitalic "+ OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN   +"
	echoitalic "+ CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                            +"
	echoitalic "+                                                                                                       +"
	echoitalic "+ Reference: http://opensource.org.                                                                     +"
	echoitalic "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
}


echocyan "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echocyan "                           POLYGENIC SCORE CALCULATIONS"
echocyan ""
echocyan ""
echocyan "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echobold "*Setting the environment."

PROJECTDIR="/hpc/dhl_ec/svanderlaan/projects/polygenicscores"
HERCULESTOOLKIT="/hpc/local/CentOS7/dhl_ec/software/HerculesToolKit"
AEGSDATA="/hpc/dhl_ec/data/_ae_originals"

echobold "* Making some directories."
### Make directories for script if they do not exist yet (!!!PREREQUISITE!!!)
if [ ! -d ${PROJECTDIR}/PRSICEINOUYE/ ]; then
	echo "The PRSICE directory does not exist, Mr. Bourne will make it for you!"
	mkdir -v ${PROJECTDIR}/PRSICEINOUYE/
fi

chmod -R a+rwx ${PROJECTDIR}

echobold "* Setting some variables specific for PRSice calculations."
echoitalic " > data files ..."

echo "  - creating AEGS to GWAS data variantID update list ..."
# SNPID RSID Chr BP A_allele B_allele MinorAllele MajorAllele AA AB BB AA_calls AB_calls BB_calls MAF HWE missing missing_calls Info CAF
# --- 1:10177:A:AC 01 10177 A AC AC A 554.34 731.78 239.87 46 56 8 0.39696 0.018899 6.3195e-06 0.92792 0.35024 0.396962
# --- 1:10235:T:TA 01 10235 T TA TA T 1524.2 1.8055 0 1523 0 0 0.00059159 4.8216e-17 0 0.0019659 0.26078 0.000591577
# --- rs145072688:10352:T:TA 01 10352 T TA TA T 490.88 755.85 279.26 46 55 15 0.43066 0.14578 1.0239e-05 0.92398 0.34431 0.430661

# echo "variantid rsid_aegs" > ${PROJECTDIR}/aegs_combo_1kGp3GoNL5_RAW.matching_variantID.list
# zcat ${AEGSDATA}/AEGS_COMBINED_IMPUTE2_1000Gp3_GoNL5/aegs_combo_1kGp3GoNL5_RAW.stats.gz | \
# awk '{ print $3, $4, $2 }' | \
# awk '{ printf "%2i %s %s\n" , $1, $2 , $3 }' | \
# awk '{ print $1":"$2, $3 }' | \
# tail -n +2 >> ${PROJECTDIR}/aegs_combo_1kGp3GoNL5_RAW.matching_variantID.list
# gzip -fv ${PROJECTDIR}/aegs_combo_1kGp3GoNL5_RAW.matching_variantID.list

echo "  - creating 1000G phase 3 update list (for ISGC data)..."
### Header
### #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
### 1	10177	rs367896724	A	AC	100	PASS	AC=2130;AF=0.425319;AN=5008;NS=2504;DP=103152;EAS_AF=0.3363;AMR_AF=0.3602;AFR_AF=0.4909;EUR_AF=0.4056;SAS_AF=0.4949;AA=|||unknown(NO_COVERAGE);VT=INDEL
### 1	10235	rs540431307	T	TA	100	PASS	AC=6;AF=0.00119808;AN=5008;NS=2504;DP=78015;EAS_AF=0;AMR_AF=0.0014;AFR_AF=0;EUR_AF=0;SAS_AF=0.0051;AA=|||unknown(NO_COVERAGE);VT=INDEL
### 1	10352	rs555500075	T	TA	100	PASS	AC=2191;AF=0.4375;AN=5008;NS=2504;DP=88915;EAS_AF=0.4306;AMR_AF=0.4107;AFR_AF=0.4788;EUR_AF=0.4264;SAS_AF=0.4192;AA=|||unknown(NO_COVERAGE);VT=INDEL

RAW1000Gdata="/hpc/dhl_ec/data/references/1000G/Phase3/VCF_format/"
# echo "rsid variantid chr position" > ${PROJECTDIR}/1000G.phase3version5b.update_VARIANTIDCHRBP.txt
# zcat ${RAW1000Gdata}/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz | \
# grep -v "##" | awk '$3 ~ '/^rs/' || $3 == "ID" ' | \
# awk '{ print $3, $1":"$2, $1, $2 }' | \
# tail -n +2 >> ${PROJECTDIR}/1000G.phase3version5b.update_VARIANTIDCHRBP.txt
# gzip -fv ${PROJECTDIR}/1000G.phase3version5b.update_VARIANTIDCHRBP.txt

echoitalic "* Inouye M et al. metaGRS ..."
if [ ! -d ${PROJECTDIR}/Inouye_bioRxiv_2018/ ]; then
	echo "The INOUYE directory does not exist, Mr. Bourne will make it for you!"
	mkdir -v ${PROJECTDIR}/Inouye_bioRxiv_2018/
fi
INOUYEDIR=${PROJECTDIR}/Inouye_bioRxiv_2018

# Head INOUYE
# chr	position	rsid	allele1	allele2	effect_allele	beta
#  1	  2245570	rs2843152	C	G	G	-2.76009e-02
#  1	 22132518	rs35465346	A	G	G	 2.39340e-02
#  1	 38386727	rs28470722	A	G	G	-1.74935e-02

# echo "variantid rsid rsid_aegs chr position effect_allele other_allele beta P_ukbb" > ${INOUYEDIR}/metaGRS_hg19_20180205.foo
# zcat ${INOUYEDIR}/metaGRS_hg19_20180205.txt.gz | \
# parseTable --col chr,position,rsid,allele1,allele2,effect_allele,beta | \
# awk '{ if($6 == $5) { print $1":"$2, $3, $3, $1, $2, $6, $4, $7, "NA" } else { print $1":"$2, $3, $3, $1, $2, $6, $5, $7, "NA" } }' >> \
# ${INOUYEDIR}/metaGRS_hg19_20180205.foo

# mergeTablesv2 \
# --file1 ${UKBBDIR}/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.4pvalupdate.txt \
# --file2 ${INOUYEDIR}/metaGRS_hg19_20180205.foo \
# --index variantid --format NORM --replace > ${INOUYEDIR}/metaGRS_hg19_20180205.4PRSICE.foo
# 
# mergeTablesv2 \
# --file1 ${PROJECTDIR}/aegs_combo_1kGp3GoNL5_RAW.matching_variantID.list.gz \
# --file2 ${INOUYEDIR}/metaGRS_hg19_20180205.4PRSICE.foo \
# --index variantid --format GZIP1 --replace > ${INOUYEDIR}/metaGRS_hg19_20180205.4PRSICE.txt

INOUYEBASEDATA="${INOUYEDIR}/metaGRS_hg19_20180205.4PRSICE.txt"

echoitalic " > setting the TARGETDATA (AEGS) ..."
TARGETDATA="/hpc/dhl_ec/data/_ae_originals/AEGS_COMBINED_IMPUTE2_1000Gp3_GoNL5/aegs_combo_1kGp3GoNL5_RAW_chr# "

echoitalic " > getting a phenotype file and exclusion list ..."
echo "FID IID AEGS_type COHORT STUDY_TYPE sex Age AgeSQR OR_year OR_year_C PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 Calcification_bin Collagen_bin Fat10_bin Fat40_bin Macrophages_bin SMC_bin IPH Macrophages_BC Mastcells_BC Neutrophils_BC SMC_BC VesselDensityAvg_BC" > ${PROJECTDIR}/aegscombo_phenocov.pheno
cat ${AEGSDATA}/pheno_cov_exclusions/aegscombo_phenocov.sample | \
parseTable --col ID_1,ID_2,AEGS_type,COHORT,STUDY_TYPE,sex,Age,AgeSQR,OR_year,OR_year_C,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,Calcification_bin,Collagen_bin,Fat10_bin,Fat40_bin,Macrophages_bin,SMC_bin,IPH,Macrophages_BC,Mastcells_BC,Neutrophils_BC,SMC_BC,VesselDensityAvg_BC | tail -n +3 | \
sed 's/FEMALE/2/g' | sed 's/MALE/1/g' >> ${PROJECTDIR}/aegscombo_phenocov.pheno

echoitalic " > PRSice general and plotting settings ..."
# PRSICESETTINGS="--no-clump --print-snp --extract PRSice.valid --score sum --missing center --all-score --perm 10000"
# PRSICESETTINGS="--no-clump --print-snp --extract PRSice.valid --score sum --all-score "

PRSICEVALID="${INOUYEDIR}/PRSice.valid"

#--extract PRSice.valid 
PRSICESETTINGS="--print-snp --extract ${PRSICEVALID} --score sum --all-score "
PRSICEPLOTTING="--fastscore --bar-col-high #E55738 --bar-col-low #1290D9 --quantile 100 --quant-break 2.5,5,10,20,40,60,80,90,95,97.5,100 --quant-ref 60"
PRSICETHREADS="8"
PRSICESEED="91149214"
PRSICEBARLEVELS="0.00000005,0.0000005,0.000005,0.00005,0.0005,0.005,0.05,0.1,0.2,0.3,0.4,0.5,1"
# PRSICEBARLEVELS="0.00000005,0.000005,0.0005,0.001,0.01,0.05,0.1,0.2,0.5,1"
SCORETYPE="PRSsum"
DATATYPE="BED"
PERMUTATION="NOPERM_NOCENTER"

PHENOTYPEFILE="${PROJECTDIR}/aegscombo_phenocov.pheno"
COVARIATES="sex,Age,OR_year_C,@PC[1-4],COHORT"
COVARIATESFACTOR="sex,COHORT"
COVARIATESFILE="${PROJECTDIR}/aegscombo_phenocov.pheno"
EXCLUSION="${PROJECTDIR}/aegscombo.exclusion_nonCEA.list"

STATTYPE="--beta" 
SNPID="rsid_aegs"
CHRID="chr"
BPID="position"
A1ID="effect_allele"
A2ID="other_allele"
STATID="beta"

# PVALUEID="P_ukbb"
PVALUEID="P"

if [ ! -d ${PROJECTDIR}/PRSICEINOUYE ]; then
	echo "The metaGRS directory does not exist, Mr. Bourne will make it for you!"
	mkdir -v "${PROJECTDIR}/PRSICEINOUYE"
fi
BASEDATA="${INOUYEBASEDATA}"
PRSICEDIR="${PROJECTDIR}/PRSICEINOUYE"

for PHENO in IPH; do 
# for PHENO in Calcification_bin Collagen_bin Fat10_bin Fat40_bin Macrophages_bin SMC_bin IPH; do 
	echoitalic " > statistics ..."
	TARGETTYPE="T"

	echoitalic " > phenotypes and covariates ..."
	PHENOTYPE="${PHENO}"
	PRSICEOUTPUTNAME="--out PRSice.${PHENOTYPE}.${SCORETYPE}.${DATATYPE}.${PERMUTATION}"

	cd ${PRSICEDIR}

	prsice.R --prsice $(command -v prsice) \
	--dir ${PRSICEDIR} \
	--seed ${PRSICESEED} \
	--bar-levels ${PRSICEBARLEVELS} \
	--base ${BASEDATA} \
	--target ${TARGETDATA} \
	--thread ${PRSICETHREADS} \
	${STATTYPE} \
	--binary-target ${TARGETTYPE} \
	--snp ${SNPID} \
	--chr ${CHRID} \
	--bp ${BPID} \
	--A1 ${A1ID} \
	--A2 ${A2ID} \
	--stat ${STATID} \
	--pvalue ${PVALUEID} \
	--cov-file ${COVARIATESFILE} \
	--cov-col ${COVARIATES} \
	--cov-factor ${COVARIATESFACTOR} \
	--pheno-file ${PHENOTYPEFILE} \
	--pheno-col ${PHENOTYPE} \
	--remove ${EXCLUSION} \
	${PRSICEPLOTTING} \
	${PRSICESETTINGS} \
	${PRSICEOUTPUTNAME}

done

# for PHENO in Macrophages_BC Mastcells_BC Neutrophils_BC SMC_BC VesselDensityAvg_BC; do 
# 	echoitalic " > statistics ..."
# 	TARGETTYPE="F"
# 
# 	echoitalic " > phenotypes and covariates ..."
# 	PHENOTYPE="${PHENO}"
# 	PRSICEOUTPUTNAME="--out PRSice.${PHENOTYPE}.${SCORETYPE}.${DATATYPE}.${PERMUTATION}"
# 
# 	cd ${PRSICEDIR}
# 
# 	prsice.R --prsice $(command -v prsice) \
# 	--dir ${PRSICEDIR} \
# 	--seed ${PRSICESEED} \
# 	--bar-levels ${PRSICEBARLEVELS} \
# 	--base ${BASEDATA} \
# 	--target ${TARGETDATA} \
# 	--thread ${PRSICETHREADS} \
# 	${STATTYPE} \
# 	--binary-target ${TARGETTYPE} \
# 	--snp ${SNPID} \
# 	--chr ${CHRID} \
# 	--bp ${BPID} \
# 	--A1 ${A1ID} \
# 	--A2 ${A2ID} \
# 	--stat ${STATID} \
# 	--pvalue ${PVALUEID} \
# 	--cov-file ${COVARIATESFILE} \
# 	--cov-col ${COVARIATES} \
# 	--cov-factor ${COVARIATESFACTOR} \
# 	--pheno-file ${PHENOTYPEFILE} \
# 	--pheno-col ${PHENOTYPE} \
# 	--remove ${EXCLUSION} \
# 	${PRSICEPLOTTING} \
# 	${PRSICESETTINGS} \
# 	${PRSICEOUTPUTNAME}
# 
# done

echo ""
echocyan "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echobold "Wow. I'm all done buddy. What a job 😱 ! let's have a 🍻🍻 ... 🖖 "

script_copyright_message
