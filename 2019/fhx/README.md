Family history and polygenic risk of cardiovascular disease: independent factors associated with secondary cardiovascular manifestations in patients undergoing carotid endarterectomy
===========================================================

#### This readme
This readme accompanies the paper "Family history and polygenic risk of cardiovascular disease: independent factors associated with secondary cardiovascular manifestations in patients undergoing carotid endarterectomy." by Timmerman N. and [Timmerman N. *et al*. **bioRxiv**](https://doi.org/10.1101/19006718).

--------------

#### Abstract

**Background** Family history (FHx) of cardiovascular disease (CVD) is a risk factor for CVD and a proxy for cardiovascular heritability. Polygenic risk scores (PRS) summarizing >1 million variants for coronary artery disease (CAD) are associated with incident and recurrent CAD events. However, little is known about the influence of FHx or PRS on secondary cardiovascular events (sCVE) in patients undergoing carotid endarterectomy (CEA). 

**Methods** We included 1,788 CEA patients from the Athero-Express Biobank ([www.atheroexpress.nl](http://www.atheroexpress.nl)). A weighted PRS for CAD including 1.7 million variants was calculated (MetaGRS). The composite endpoint of sCVE during three years follow-up included coronary, cerebrovascular and peripheral events and cardiovascular death. We assessed the impact of FHx and MetaGRS on sCVE and carotid plaque composition. 

**Results** Positive FHx was associated with a higher 3-year risk of sCVE independent of cardiovascular risk factors and MetaGRS (adjusted HR 1.40, 95%CI 1.07-1.82, p=0.013). Patients in the highest MetaGRS quintile had a higher 3-year risk of sCVE compared to the rest of the cohort independent of cardiovascular risk factors including FHx (adjusted HR 1.35, 95%CI 1.01-1.79, p=0.043), and their atherosclerotic plaques contained more fat (adjusted OR 1.59, 95%CI, 1.11-2.29, p=0.013) and more macrophages (OR 1.49, 95%CI 1.12-1.99, p=0.006). 

**Conclusion** In CEA patients, both positive FHx and higher MetaGRS were independently associated with increased risk of sCVE. Moreover, higher MetaGRS was associated with vulnerable plaque characteristics. Future studies should unravel underlying mechanisms and focus on the added value of PRS and FHx in individual risk prediction for sCVE.

--------------

*Competing Interest Statement*

The authors have declared no competing interest.

*Funding Statement*

Dr. Van der Laan was funded through grants from the Netherlands CardioVascular Research Initiative of the Netherlands Heart Foundation (CVON 2011/B019 and CVON 2017-20: Generating the best evidence-based pharmaceutical targets for atherosclerosis [GENIUS I&II]). Dr. Timmerman, Dr. De Borst, dr. De Kleijn, and dr. Pasterkamp are funded through the EU 755320 Taxinomisis grant. We acknowledge the European Research Area Network on Cardiovascular diseases (ERA-CVD, grant number 01KL1802). There are no competing interests declared.


*Author Declarations*

All relevant ethical guidelines have been followed and any necessary IRB and/or ethics committee approvals have been obtained.

  > Yes

All necessary patient/participant consent has been obtained and the appropriate institutional forms have been archived.

  > Yes

Any clinical trials involved have been registered with an ICMJE-approved registry such as ClinicalTrials.gov and the trial ID is included in the manuscript.

  > Not Applicable

I have followed all appropriate research reporting guidelines and uploaded the relevant Equator, ICMJE or other checklist(s) as supplementary files, if applicable.

  > Yes


#### Analysis Scripts

Surely these scripts will not work immediately on your systems, but they may be used and edited for local use.

- *[create.aegs.dataset.sh](create.aegs.dataset.sh)*</br>
This script was made during the revision to extract only the 1.7 million variants from the MetaGRS and calculate summary statistics of all variants in the AEGS dataset.

- *[prsice.fhx.metagrs.sh](prsice.fhx.metagrs.sh)*</br>
This script was used to calculate the MetaGRS in the Athero-Express Biobank Study using [PRSice-2](http://www.prsice.info)[^1].

- *[syntax.sas]()*</br>
This SPSS syntax was used for overal analyses.


#### Notes
Scripts will work within the context of a certain Linux environment, for example a CentOS7 system on a SUN Grid Engine background or macOS X Lion+ (version 10.7.[x]+). 


#### References
[^1] Choi SW, and O’Reilly PF. "PRSice-2: Polygenic Risk Score Software for Biobank-Scale Data." **GigaScience** 8, no. 7 (July 1, 2019). [https://doi.org/10.1093/gigascience/giz082](https://doi.org/10.1093/gigascience/giz082).

--------------

#### The MIT License (MIT)
##### Copyright (c) 1979-present Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:   

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Reference: http://opensource.org.
