Publications
============
[![DOI](https://zenodo.org/badge/111704113.svg)](https://zenodo.org/badge/latestdoi/111704113)

#### Disclaimer
It is my sincere intend to be as open and forthcoming as I can be, and therefore to do my utmost best to share my scripts whenever is possible. 

#### This repository
This repositories contains scripts organized in folders for published projects (as of August 30, 2016). All scripts are annotated for debugging purposes - and future reference. Therefore, most scripts will be self-explanatory and this page only serves as a table of contents.

#### EGA
All (raw) genotype and methylation data will be made available through EGA. When the EGA-identifiers are known, I will put them here.

#### Notes
Scripts will work within the context of a certain Linux environment, for example a CentOS7 system on a SUN Grid Engine background or macOS X Lion+ (version 10.7.[x]+). 


--------------

## 2018

### CTMM: eQTL analysis in monocytes

[more information to follow]

For the eQTL analyses in monocytes from CTMM we developed the [QTLToolKit](https://github.com/swvanderlaan/QTLToolKit) which is using [QTLTools](https://qtltools.github.io/qtltools/).

Link to the publication:.


### Meta-analysis of GWAS for serum FABP4 levels

[more information to follow]

For this meta-analysis we redesigned and updated MANTEL (see [here](http://debakker.med.harvard.edu/resources.html) and [here](https://www.ncbi.nlm.nih.gov/pubmed/18852200)) and created [MetaGWASToolKit](https://github.com/swvanderlaan/MetaGWASToolKit) to this end. **MetaGWASToolKit** works with 1000G data (phase 1 and phase 3) as well as the Dutch 'optimized' reference 1000G+GoNL5.

Link to the publication:.


--------------

## 2017

### Smoking is associated to DNA methylation in atherosclerotic carotid lesions.

For this study we performed an **E**pigenome-**W**ide **A**ssociation **S**tudy (EWAS) of DNA derived from whole-blood and carotid plaques.

#### Description of scripts
##### Analyses

| **Script**														| **Description** |
| :---																| :--- |
| `20180614.aems450k.meta.analysis.ewas.smoking.plaque.R` 			| Meta-analysis of EWAS on current smoking-status in Athero-Express Methylation Study 1 (AEMS450K1) and Athero-Express Methylation Study 2 (AEMS450K2). |
| `20180614.aems450k.meta.analysis.ewas.ePackYears.plaque.R` 		| Meta-analysis of EWAS of number of Estimated Pack Years in AEMS450K1 and AEMS450K2. This was a sensitivity analysis. |
| `20180614.aems450k.meta.analysis.correct.eCigs.forGitHub.R` 		| Re-calculation of the number of Estimated Pack Years in AEMS450K1 and AEMS450K2. |
| `20180614.aems450k1.analysis.smoking.R` 							| EWAS on current smoking-status in the AEMS450K1 (blood and plaque). |
| `20180614.aems450k2.analysis.smoking.R` 							| EWAS on current smoking-status in the AEMS450K2 (plaque only). |
| `20180614.aems450k.meta.analysis.ewas.Baseline.R` 				| Calculate baseline characteristics of AEMS450K1 and AEMS450K2 |
| `20180614.aems450k.meta.analysis.ewas.Smoking_vs_RNA.R` 			| RNA-sequencing pilot study |
| `20180614.aems450k.meta.analysis.ewas.CpG_vs_PlaquePheno.R`		| Regression of CpGs with plaque phenotypes. |
| `20180614.aems450k.meta.analysis.ewas.CpG_vs_RNA.R`				| Compare CpG methylation and RNA expression. |
| `20180614.aems450k.meta.analysis.ewas.Smoking_vs_PlaquePheno.R`	| Linear and logistic regression of current smoking status with plaque phenotypes. _Not used in the article, but maybe useful for some._ |
| `ann2expann.pl`													| To edit the annotation file of Illumina Methylation 450K for annotation of results. |

##### 'Calling' and normalisation of methylation array data.
For the quality control we used [DNAmArray](https://github.com/molepi/DNAmArray) and followed the instructions in the [DNAmArray Workflow](https://molepi.github.io/DNAmArray_workflow/). Below some of the scripts used to 'call' and 'normalise' the methylation data.

| **Script**														| **Description** |
| :---																| :--- |
| `metharray.prepare.sh`					| Script to perform initial 'calling' of the methylation data. A CentOS7-style high-performance computer cluster is required.|
| `metharray.target.creator.v1.R`			| Creates 'targets' files.|
| `metharray.summarizer.v1.R`				| Summarise the methylation array data.|
| `metharray.visualizer.v1.R`				| Interactively visualize the methylation data to detect outliers. Best to do this locally with the required files.|
| `metharray.rgset.creator.parallel.v1.R`	| Create the RGChannelSet file.|
| `metharray.fun.normalizer.v1.R`			| Functional normalisation of the RGChannelSet file.|

For the mQTL analyses in the [Athero-Express](http://www.atheroexpress.nl) I developed the [fastQTLToolKit](https://github.com/swvanderlaan/fastQTLToolKit) which is using [fastQTL](https://www.ncbi.nlm.nih.gov/pubmed/26708335).
All (raw) genotype and methylation data will be made available soon via EGA; a link is supplied above. For access to BiKE and STAGE data contact with the respective 'owners' should be made.

Link to the publication: revision is submitted; link follows [asap].


### CVD loci

For this I used mainly [GWASToolKit](https://github.com/swvanderlaan/GWASToolKit). Some additional useful scripts I used are found in the [HerculesToolKit](https://github.com/swvanderlaan/HerculesToolKit). The polygenic (risk) scores were made with a script we (Jessica van Setten and I) had gotten via Eli Stahl. All (raw) genotype data will be made available soon via EGA; a link is supplied above.

Link to the publication: revision is submitted; link follows [asap].


--------------

## 2016

### Cystatin C and Cardiovascular Disease: A Mendelian Randomization Study.

The aim of this study was to use Mendelian randomization to investigate whether cystatin C is causally related to CVD in the general population. We published this in [JACC](https://www.ncbi.nlm.nih.gov/pubmed/?term=27561768).

Here you can find the R-script used to calculate power.


--------------

#### The MIT License (MIT)
##### Copyright (c) 2016-2018 Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:   

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Reference: http://opensource.org.
