Publications
============
[![DOI](https://zenodo.org/badge/111704113.svg)](https://zenodo.org/badge/latestdoi/111704113)

#### Disclaimer
It is my sincere intend to be as open and forthcoming as I can be, and therefore to do my utmost best to share my scripts whenever is possible. 

#### This repository
This repositories contains scripts organized in folders for ongoing projects, once published a stand-alone repository will be made and linked here (as of August 30, 2016). All scripts are annotated (to some extend) for debugging purposes - and future reference. Therefore, most scripts will be self-explanatory and this page only serves as a table of contents.

#### EGA
All (raw) genotype and methylation data will be made available through EGA at some point in time. We are actively pursuiting this, but it takes time; watch this space. When the EGA-identifiers are known, I will put them here.

#### Notes
Scripts will work within the context of a certain Linux environment, for example a CentOS7 system on a SUN Grid Engine background or macOS X Lion+ (version 10.7.[x]+). 


--------------

## 2019

### [Family history and polygenic risk of cardiovascular disease: independent factors associated with secondary cardiovascular manifestations in patients undergoing carotid endarterectomy](https://github.com/swvanderlaan/2020_pid_timmerman_n_fhx_prs_cad)

Here we calculated a polygenic risk score for coronary artery disease based on the MetaGRS described in [Inouye M. *et al.*](http://www.onlinejacc.org/content/72/16/1883){target="_blank"}, and correlated this with family history (FHx) of cardiovascular disease, secondary major adverse clinical events after carotid endarterectomy, and atherosclerotic plaque characteristics.

Link to the publication: [Timmerman N. *et al.*](https://www.medrxiv.org/content/10.1101/19006718v2).


### [Common variants in OSMR contribute to carotid plaque vulnerability](2019/vankeulen_d_oncostatin)

In this study we studied the effects of *OSMR* on plaque composition in the Athero-Express Biobank Study.

Link to the publication: [Van Keulen D. *et al.*](https://www.biorxiv.org/content/10.1101/576793v1).


--------------

## 2018

### [CTMM: eQTL analysis in monocytes](2018/ctmm)

[more information to follow]

For the eQTL analyses in monocytes from CTMM we developed the [QTLToolKit](https://github.com/swvanderlaan/QTLToolKit) which is using [QTLTools](https://qtltools.github.io/qtltools/).

Link to the publication:.


### [Meta-analysis of GWAS for serum FABP4 levels](2018/fabp4)

[more information to follow]

For this meta-analysis we redesigned and updated MANTEL (see [here](http://debakker.med.harvard.edu/resources.html) and [here](https://www.ncbi.nlm.nih.gov/pubmed/18852200)) and created [MetaGWASToolKit](https://github.com/swvanderlaan/MetaGWASToolKit) to this end. **MetaGWASToolKit** works with 1000G data (phase 1 and phase 3) as well as the Dutch 'optimized' reference 1000G+GoNL5.

Link to the publication:.


--------------

## 2017

### [Smoking is associated to DNA methylation in atherosclerotic carotid lesions.](https://github.com/swvanderlaan/2018_circgenomprecismed_30354327_siemelink_m_smoking_ewas)

For this study we performed an **E**pigenome-**W**ide **A**ssociation **S**tudy (EWAS) of DNA derived from whole-blood and carotid plaques.

Link to the publication: [_Circ Genom Precis Med._ 2018;11:e002030](https://www.ahajournals.org/doi/full/10.1161/CIRCGEN.117.002030).


### [CVD loci](https://github.com/swvanderlaan/2018_circgenomprecmed_30354329_van_der_laan_sw_cvd_loci)

For this I used mainly [GWASToolKit](https://github.com/swvanderlaan/GWASToolKit). Some additional useful scripts I used are found in the [HerculesToolKit](https://github.com/swvanderlaan/HerculesToolKit). The polygenic (risk) scores were made with a script we (Jessica van Setten and I) had gotten via Eli Stahl. All (raw) genotype data will be made available soon via EGA; a link is supplied above.

Link to the publication: [_Circ Genom Precis Med._ 2018;11:e002115](https://www.ahajournals.org/doi/full/10.1161/CIRCGEN.118.002115).


--------------

## 2016

### [Cystatin C and Cardiovascular Disease: A Mendelian Randomization Study.](https://github.com/swvanderlaan/2016_jacc_27561768_vanderlaan_sw_cystatinc)

The aim of this study was to use Mendelian randomization to investigate whether cystatin C is causally related to CVD in the general population. We published this in [JACC](https://www.ncbi.nlm.nih.gov/pubmed/?term=27561768).


--------------

#### The MIT License (MIT)
##### Copyright (c) 1979-present Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:   

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Reference: http://opensource.org.
