Epigenome-Wide Association Study of Current Tobacco Smoking
===========================================================

#### This readme
This readme accompanies the paper "Smoking is associated to DNA methylation in atherosclerotic carotid lesions." by Siemelink M. and Van der Laan S.W. *et al*.

#### Abstract
*Background:* Tobacco smoking is a major risk factor for atherosclerotic disease and has been associated with DNA methylation (DNAm) changes in blood cells. However, whether smoking influences DNAm in the diseased vascular wall is unknown but may prove crucial in understanding the pathophysiology of atherosclerosis. In the this study we associated current tobacco smoking to epigenome-wide DNAm in atherosclerotic plaques from patients undergoing carotid endarterectomy (CEA).

*Methods:* DNAm at commonly methylated sites (CpGs) was assessed in atherosclerotic plaque samples and peripheral blood samples from 485 CEA patients. We tested the association of current tobacco smoking with DNAm corrected for age, and sex. To control for bias and inflation due to cellular heterogeneity we applied a Bayesian method to estimate an empirical null distribution. Replication of the smoking associated methylated CpGs in atherosclerotic plaques was executed in a second sample of 190 CEA patients, and results were meta-analyzed using a fixed-effects model.

*Results:* Tobacco smoking was significantly associated to differential DNAm in atherosclerotic lesions of 4 CpGs (FDR < 0.05) mapped to 2 different genes (*AHRR*, *ITPK1*). The strongest associations were found for CpGs mapped to the gene *AHRR*, a repressor of the aryl hydrocarbon receptor transcription factor involved in xenobiotic detoxification. One of these methylated CpGs were found to be regulated by local genetic variation.

*Conclusion:* Tobacco smoking associates with DNA methylation at multiple loci in carotid atherosclerotic lesions. These observations support further investigation of the relationship between risk factors and epigenetic regulation in atherosclerotic disease.


#### Quality Control
We followed the [Leiden DNAmArray Workflow](https://github.com/molepi/DNAmArray) for Quality Control (QC) which is developed by the team of Bas T. Heijmans at the LUMC (Leiden, the Netherlands). To this end we created several, fully annotated scripts.

- *metharray.prepare.sh*</br>
BASH script to perform the QC on a High Performance Computer Cluster (HPC). This script also includes some lines to extract relevant SNPs from imputed genotyping data - the same SNPs present on the Illumina Methylation 450K BeadArray.
- *metharray.target.creator.v1.R*</br>
R-script to create the target list (i.e. the list of IDAT-files) for downstream QC.
- *metharray.summarizer.v1.R*</br>
R-script to summarise the Illumina Methylation 450K data for downstream QC.
- *metharray.visualizer.v1.R*</br>
R-script to interactively visualize the QC using a browser. QC images and a table with excluded samples can be downloaded here as well. Depending on the capabilities and settings of the local HPC, the user could do this locally as it does not require a massive amount of computing power or memory. I did this on my MacBook Pro mid-2012 i5 2.4GHz with 16Gb RAM and 1Tb SSD.
- *metharray.rgset.creator.parallel.v1.R*</br>
R-script to create an extended RGChannelSet.
- *metharray.fun.normalizer.v1.R*</br>
R-script to normalize the RGChannelSet to a GR-object.


#### Analysis Scripts
Surely these scripts will not work immediately on your systems, but they may be used and edited for local use.
 
- *20180131.aems450k1.analysis.ewas.smoking.R*</br>
EWAS in the discovery sample (AEMS450K1). Includes analyses on DNA methylation (DNAm) in carotid plaques and whole blood.
- *20180131.aems450k2.analysis.ewas.smoking.R*</br>
EWAS in the replication sample (AEMS450K2). Only includes analyses on DNA methylation (DNAm) in carotid plaques.
- *20180131.aems450k.meta.analysis.ewas.smoking.plaque.R*</br>
Meta-analysis of EWAS using data from the discovery sample (AEMS450K1) and replication sample (AEMS450K2). Only includes analyses on DNA methylation (DNAm) in carotid plaques.
- *20180131.aems450k.meta.analysis.correct.eCigs.forGitHub.R*</br>
Calculate estimated number of pack years smoking in the data.
- *20180131.aems450k.meta.analysis.ewas.ePackYears.plaque.R*</br>
EWAS on estimated number of pack years smoking. Sensitivity analysis.
- *20180131.aems450k.meta.analysis.ewas.Smoking_vs_PlaquePheno.R*</br>
R-script to perform analysis of 'current smoking' on histological plaque characteristics.
- *20180131.aems450k.meta.analysis.ewas.Baseline.R*</br>
Calculate the baseline summary statistics of the sample sets.
- *20180131.aems450k.meta.analysis.ewas.Cox.R*</br>
Perform Cox-regression using the significant CpGs.
- *20180131.aems450k.meta.analysis.ewas.CpG_vs_PlaquePheno.R*</br>
Association of the significant CpGs with histological plaque characteristics.
- *20180131.aems450k.meta.analysis.ewas.Smoking_vs_RNA.R*</br>
Association of 'current smoking' with gene counts in carotid plaques. A pilot RNA-sequencing experiment.
- *20180213.aems450k.meta.analysis.ewas.CpG_vs_RNA.R*</br>
Association of significant CpGs with gene counts from the RNA-sequencing experiment. 

#### Notes
Scripts will work within the context of a certain Linux environment, for example a CentOS7 system on a SUN Grid Engine background or macOS X Lion+ (version 10.7.[x]+). 


--------------

#### The MIT License (MIT)
##### Copyright (c) 1979-present Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:   

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Reference: http://opensource.org.
