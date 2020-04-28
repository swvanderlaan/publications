* Encoding: UTF-8.

** definitieve syntax 
- FH composite vs secundary events and plaque fenotypes
- Polygenic risk score Inouye vs secundary events and plaque fenotypes.


******************************************************************************************************************************************************

*Inouye exploratives
******************************************************************************************************************************************************

*explore Inouye. 
EXAMINE VARIABLES=CAD_PRSInouye
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

*is er een effect per chip, geslacht of ziekenhuis?.

EXAMINE VARIABLES=CAD_PRSInouye BY Gender AEGS_type Hospital 
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.


*boxplot per jaar.
EXAMINE VARIABLES=CAD_PRSInouye BY ORyear
  /PLOT=BOXPLOT
  /STATISTICS=NONE
  /NOTOTAL.

*boxplot year7 (before 2007 vs after 2007).

EXAMINE VARIABLES=CAD_PRSInouye BY Year7
  /PLOT=BOXPLOT
  /STATISTICS=NONE
  /NOTOTAL.


******************************************************************************************************************************************************
*Inouye vs CAD (=validatie)
******************************************************************************************************************************************************.

**********************************Inouye vs CAD_history = Validation*********************************.

  LOGISTIC REGRESSION VARIABLES CAD_history
  /METHOD=ENTER CAD_PRSInouye
  /METHOD=ENTER Age
   /METHOD=ENTER Gender
   /METHOD=ENTER AEGS_type
   /METHOD=ENTER PC1 PC2 PC2 PC3
   /CONTRAST (Gender)=Indicator(1)
   /CONTRAST (AEGS_type)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


******************************************************************************************************************************************************
*correleert Inouye score met FH composite? Dit zouden we wel verwachten namelijk = validatie.
******************************************************************************************************************************************************

**FH_composite2**.
*Relatie PRS en FH_composite (eerste graads familie leden)- dit zou wel moeten want je verwacht dat deze twee gecorreleerd zijn.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT CAD_PRSInouye
  /METHOD=ENTER FH_CVD_comp2
  /METHOD=ENTER Age
  /METHOD=ENTER Gender 
  /METHOD=ENTER AEGS_type
  /METHOD=ENTER  PC1 PC2 PC3 PC4.


LOGISTIC REGRESSION VARIABLES FH_CVD_comp2
  /METHOD=ENTER CAD_PRSInouye
  /METHOD=ENTER Age
   /METHOD=ENTER Gender
      /METHOD=ENTER  AEGS_type
      /METHOD=ENTER PC1 PC2 PC3 PC4 
  /CONTRAST (Gender)=Indicator(1)
   /CONTRAST (AEGS_type)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).



*is er een chip-effect in relatie CAD_PRSInouye vs FH_composite ?.

SORT CASES  BY AEGS_type.
SPLIT FILE LAYERED BY AEGS_type.

EXAMINE VARIABLES=CAD_PRSInouye BY FH_CVD_comp2
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

T-TEST GROUPS=FH_CVD_comp2(0 1)
  /MISSING=ANALYSIS
  /VARIABLES=CAD_PRSInouye
  /CRITERIA=CI(.95).


REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT CAD_PRSInouye
  /METHOD=ENTER FH_CVD_comp2
  /METHOD=ENTER Age
  /METHOD=ENTER Gender 
  /METHOD=ENTER  PC1 PC2 PC3 PC4.


LOGISTIC REGRESSION VARIABLES FH_CVD_comp2
  /METHOD=ENTER CAD_PRSInouye
  /METHOD=ENTER Age
   /METHOD=ENTER Gender
      /METHOD=ENTER PC1 PC2 PC3 PC4 
  /CONTRAST (Gender)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

SPLIT FILE off.

*Is er nog per geslacht een verschil tussen de relatieCAD_PRSInouye en FH_composite 2?.

SORT CASES  BY Gender.
SPLIT FILE LAYERED BY Gender.

EXAMINE VARIABLES= CAD_PRSInouye BY FH_CVD_comp2
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

T-TEST GROUPS=FH_CVD_comp2(0 1)
  /MISSING=ANALYSIS
  /VARIABLES=CAD_PRSInouye
  /CRITERIA=CI(.95).


*Correctie voor Age, PC1-PC4 + Chiptype.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT CAD_PRSInouye
  /METHOD=ENTER FH_CVD_comp2
  /METHOD=ENTER Age
  /METHOD=ENTER AEGS_type
  /METHOD=ENTER  PC1 PC2 PC3 PC4.


LOGISTIC REGRESSION VARIABLES FH_CVD_comp2
  /METHOD=ENTER CAD_PRSInouye
  /METHOD=ENTER Age
      /METHOD=ENTER  AEGS_type
      /METHOD=ENTER PC1 PC2 PC3 PC4 
   /CONTRAST (AEGS_type)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

SPLIT file off.

******************************************************************************************************************************************************
*Recode variables- Inouye Percentile scores - per chiptye berekenen aangezien de chip effect heeft op PRS.
******************************************************************************************************************************************************.
*chiptype 0.

USE ALL.
COMPUTE filter_$=(AEGS_type = 0).
VARIABLE LABELS filter_$ 'AEGS_type = 0 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

*568 cases met Inouye chip 0. 

*Bereken percentiele Inouye van alleen chip 0.

RANK VARIABLES=CAD_PRSInouye (A)
  /RANK
  /PRINT=YES
  /TIES=MEAN.

*verander naam in RCAD_PRSChip0.


*total cases 568.
COMPUTE Percentiles_Inouyechip0=((RCAD_PRSChip0-0.5)/568) * 100.
EXECUTE.
*decimals on 0, 0 veranderen in 1.

*in Quintiles.
RECODE Percentiles_Inouyechip0 (SYSMIS=SYSMIS) (MISSING=SYSMIS) (0.0000 thru 20.0000=1) (20.0001 thru 
    40.0000=2) (40.0001 thru 60.0000=3) (60.0001 thru 80.0000=4) (80.0001 thru 100.0000=5) INTO 
    Inouyechip0_quintiles.
VARIABLE LABELS  Inouyechip0_quintiles 'quintiles Inouye chip 0'.
EXECUTE.



****zelfde voor chip 1.

USE ALL.
COMPUTE filter_$=(AEGS_type = 1).
VARIABLE LABELS filter_$ 'AEGS_type = 1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

*cases with chip 1 = 868. 

*Bereken percentiele Inouye van alleen chip 1.

RANK VARIABLES=CAD_PRSInouye (A)
  /RANK
  /PRINT=YES
  /TIES=MEAN.

*verander naam in RCAD_PRSChip1.

*total cases 868.
COMPUTE Percentiles_Inouyechip1=((RCAD_PRSChip1-0.5)/868) * 100.
EXECUTE.
*decimals on 0, 0 veranderen in 1.



*in Quintiles.
RECODE Percentiles_Inouyechip1 (SYSMIS=SYSMIS) (MISSING=SYSMIS) (0.0000 thru 20.0000=1) (20.0001 thru 
    40.0000=2) (40.0001 thru 60.0000=3) (60.0001 thru 80.0000=4) (80.0001 thru 100.0000=5) INTO 
    Inouyechip1_quintiles.
VARIABLE LABELS  Inouyechip1_quintiles 'quintiles Inouye chip 1'.
EXECUTE.

FILTER OFF.
USE ALL.
EXECUTE.


*NU variabelen maken met top 20% en lowest 20% of scores door zowel dat te doen van Chip 1 als Chip 0.

COMPUTE Inouye20th=-999.
IF (Inouyechip1_quintiles=1) Inouye 20th=0.
IF (Inouyechip1_quintiles=5) Inouye 20th=1.
IF (Inouyechip0_quintiles=1) Inouye 20th=0.
IF (Inouyechip0_quintiles=5) Inouye 20th=1.
EXECUTE.



**addendum 19-3-2019**.
*1 Quintile variabele maken.

COMPUTE Inouye_Quintile_total=-999.
IF ((Inouyechip1_quintiles=1) OR (Inouyechip0_quintiles=1)) Inouye_Quintile_total=1.
IF ((Inouyechip1_quintiles=2) OR (Inouyechip0_quintiles=2)) Inouye_Quintile_total=2.
IF ((Inouyechip1_quintiles=3) OR (Inouyechip0_quintiles=3)) Inouye_Quintile_total=3.
IF ((Inouyechip1_quintiles=4) OR (Inouyechip0_quintiles=4)) Inouye_Quintile_total=4.
IF ((Inouyechip1_quintiles=5) OR (Inouyechip0_quintiles=5)) Inouye_Quintile_total=5.
EXECUTE.


*inouye80vs20 betekent bovenste quintile vergeleken met onderste 4 quintielen, oftewel top 20% vs onderste 80%.

COMPUTE Inouye80vs20=-999.
IF (Inouye_Quintile_total=1)  Inouye80vs20=0.
IF (Inouye_Quintile_total=2)  Inouye80vs20=0.
IF (Inouye_Quintile_total=3)  Inouye80vs20=0.
IF (Inouye_Quintile_total=4)  Inouye80vs20=0.
IF (Inouye_Quintile_total=5)  Inouye80vs20=1.
EXECUTE.
VARIABLE LABELS   Inouye80vs20 'top 20% (1) versus bottom 80% (0) '.
EXECUTE.
ADD VALUE LABELS Inouye80vs20 1 'top20' 0 'bottom80' -999 'No data available/missing'.
EXECUTE.


******************************************************************************************************************************************************
*ECHTE ANALYSES Starten hier: 

*is there an association between FH composite and secondary events and plaque fenotype? 
*and is this association also present for PRS Inouye (continuous) and percentiles (lowest 20th and upper20th?).

******************************************************************************************************************************************************.

******************************************************************************************************************************************************
*Baseline tabel- FH Composite.
******************************************************************************************************************************************************
*baseline tabel maken.

*baseline tables with p-values-  categorical variables.
*Gender SmokerCurrent DM.composite risk614 Hypertension1 CAD_history Stroke_history PAOD 
    Symptoms.4g stenosis_ipsilateral stenosis_con_bin Med.statin.derived Med.all.antiplatelet YearORnew.

CTABLES
  /VLABELS VARIABLES=Gender SmokerCurrent DM.composite risk614 Hypertension1 CAD_history Stroke_history PAOD 
    Symptoms.4g stenosis_ipsilateral stenosis_con_bin Med.statin.derived Med.all.antiplatelet YearORnew FH_CVD_comp2 DISPLAY=DEFAULT 
  /TABLE Gender [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1]
+ SmokerCurrent [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] + DM.composite [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] 
+ risk614 [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] + Hypertension1 [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] 
+ CAD_history [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] + Stroke_history [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1]
+ PAOD [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1]
+ symptoms.4g [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] + Stenosis_ipsilateral [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] + stenosis_con_bin [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1]
+ Med.statin.derived [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] + Med.all.antiplatelet [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] 
+ YearORnew [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1]  BY FH_CVD_comp2 [C]
  /SLABELS POSITION=ROW
  /CATEGORIES VARIABLES= Gender SmokerCurrent DM.composite risk614 Hypertension1 CAD_history Stroke_history PAOD Symptoms.4g stenosis_ipsilateral stenosis_con_bin Med.statin.derived 
  Med.all.antiplatelet YearORnew FH_CVD_comp2 ORDER=D KEY=VALUE EMPTY=INCLUDE
  /TITLES
    TITLE='Table 1: Categorical characteristics at baseline, by FH CVD comp2'
    /SIGTEST TYPE=CHISQUARE ALPHA=0.05 INCLUDEMRSETS=YES CATEGORIES=ALLVISIBLE.

*baseline tables with p-values-  continuous variables.

*continue variabelen: Age BMI GFR_MDRD totalchol  triglyceriden LDL HDL.
CTABLES
  /VLABELS VARIABLES=Age BMI GFR_MDRD totalchol triglyceriden LDL HDL FH_CVD_comp2 DISPLAY=DEFAULT 
  /TABLE Age [S][MEAN F40.2, STDDEV 'SD' F40.2, MEDIAN F40.2] + BMI [S][MEAN F40.2, STDDEV 'SD' F40.2, MEDIAN F40.2] + GFR_MDRD [S][MEAN F40.2, STDDEV 'SD' F40.2, MEDIAN F40.2] 
  +  totalchol [S][MEAN F40.2, STDDEV 'SD' F40.2, MEDIAN F40.2] + triglyceriden [S][MEAN F40.2, STDDEV 'SD' F40.2, MEDIAN F40.2] 
  + LDL [S][MEAN F40.2, STDDEV 'SD' F40.2, MEDIAN F40.2] + HDL [S][MEAN F40.2, STDDEV 'SD' F40.2, MEDIAN F40.2] BY FH_CVD_comp2 [C]
  /SLABELS POSITION=ROW
  /CATEGORIES VARIABLES= FH_CVD_comp2 ORDER=D KEY=VALUE EMPTY=INCLUDE
  /TITLES
    TITLE='Table 1: Continuous baseline characteristics , by FH CVD comp2'.


*for median and interquartile range.

EXAMINE VARIABLES=totalchol BY FH_CVD_comp2
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

EXAMINE VARIABLES=triglyceriden BY FH_CVD_comp2
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

EXAMINE VARIABLES=LDL BY FH_CVD_comp2
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

EXAMINE VARIABLES=HDL BY FH_CVD_comp2
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

*voor lipiden eerst LN transformatie.
T-TEST GROUPS=FH_CVD_comp2(0 1)
  /MISSING=ANALYSIS
  /VARIABLES=Age BMI GFR_MDRD LN_HDL LN_totalchol LN_triglyceriden LN_LDL
  /CRITERIA=CI(.95).

*Uitgebreider met p-waardes en missing values.

  CROSSTABS
  /TABLES=Gender SmokerCurrent DM.composite risk614 Hypertension1 CAD_history Stroke_history PAOD 
    Symptoms.4g Stenosis_ipsilateral stenosis_con_bin Med.statin.derived Med.all.antiplatelet YearORnew BY FH_CVD_comp2
  /FORMAT=AVALUE TABLES
  /STATISTICS=CHISQ 
  /CELLS=COUNT COLUMN 
  /COUNT ROUND CELL.

** baseline values for total cohort for overview****.

FREQUENCIES VARIABLES=Gender Gender SmokerCurrent DM.composite risk614 Hypertension1 CAD_history Stroke_history PAOD 
    Symptoms.4g Stenosis_ipsilateral stenosis_con_bin Med.statin.derived Med.all.antiplatelet YearORnew FH_CVD_comp2
  /ORDER=ANALYSIS.

EXAMINE VARIABLES= Age BMI GFR_MDRD HDL LDL triglyceriden totalchol
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING PAIRWISE
  /NOTOTAL.


******************************************************************************************************************************************************
*FH composite versus plaque characteristics-univariable.

******************************************************************************************************************************************************.

EXAMINE VARIABLES=macmean0 smcmean0 vessel_density_averaged Mast_cells_plaque neutrophils BY FH_CVD_comp2
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING PAIRWISE
  /NOTOTAL.

*check distributions. 
EXAMINE VARIABLES=LN_macmean LN_smcmean LN_vd LN_ms LN_neutro BY FH_CVD_comp2
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING PAIRWISE
  /NOTOTAL.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_macmean
  /METHOD=ENTER FH_CVD_comp2.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_smcmean
  /METHOD=ENTER FH_CVD_comp2.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_vd
  /METHOD=ENTER FH_CVD_comp2.

*REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_ms
  /METHOD=ENTER FH_CVD_comp2.

*REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_neutro
  /METHOD=ENTER FH_CVD_comp2.

**binary.

CROSSTABS
  /TABLES=Calc.bin Collagen.bin Fat.bin_40 Fat.bin_10 IPH.bin Macrophages.bin SMC.bin  BY FH_CVD_comp2
  /FORMAT=AVALUE TABLES
  /STATISTICS=CHISQ 
  /CELLS=COUNT COLUMN 
  /COUNT ROUND CELL.

LOGISTIC REGRESSION VARIABLES Calc.bin
  /METHOD=ENTER FH_CVD_comp2 
  /CONTRAST (FH_CVD_comp2)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Collagen.bin
  /METHOD=ENTER FH_CVD_comp2 
  /CONTRAST (FH_CVD_comp2)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Fat.bin_40
  /METHOD=ENTER FH_CVD_comp2 
  /CONTRAST (FH_CVD_comp2)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Fat.bin_10
  /METHOD=ENTER FH_CVD_comp2 
  /CONTRAST (FH_CVD_comp2)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES IPH.bin
  /METHOD=ENTER FH_CVD_comp2 
  /CONTRAST (FH_CVD_comp2)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Macrophages.bin 
  /METHOD=ENTER FH_CVD_comp2 
  /CONTRAST (FH_CVD_comp2)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES SMC.bin
  /METHOD=ENTER FH_CVD_comp2 
  /CONTRAST (FH_CVD_comp2)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


******************************************************************************************************************************************************
*FH composite versus plaque characteristics-Multivariable.

  
REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_macmean
  /METHOD=ENTER FH_CVD_comp2 Age Gender OKyear Symptoms.4g
    /METHOD=ENTER Hypertension1 med.statin.derived PAOD Stenosis_ipsilateral.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_smcmean
 /METHOD=ENTER FH_CVD_comp2 Age Gender OKyear Symptoms.4g
  /METHOD=ENTER risk614 Stroke_history GFR_MDRD.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_vd
 /METHOD=ENTER FH_CVD_comp2 Age Gender OKyear Symptoms.4g
 /METHOD=ENTER BMI stenosis_ipsilateral hypertension1 CAD_history LN_triglyceriden.
 
 
**binary.
LOGISTIC REGRESSION VARIABLES Collagen.bin
  /METHOD=ENTER FH_CVD_comp2 Age Gender OKyear Symptoms.4g
   /METHOD=ENTER Hypertension1 risk614 PAOD
  /CONTRAST (FH_CVD_comp2)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
      /CONTRAST (Symptoms.4g)=Indicator(1)
         /CONTRAST (Hypertension1)=Indicator(1)
              /CONTRAST (PAOD)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


LOGISTIC REGRESSION VARIABLES  Calc.bin
  /METHOD=ENTER FH_CVD_comp2 Age Gender OKyear Symptoms.4g
   /METHOD=ENTER Hypertension1 CAD_history Stroke_history PAOD LN_totalchol 
  /CONTRAST (FH_CVD_comp2)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
      /CONTRAST (Symptoms.4g)=Indicator(1)
         /CONTRAST (Hypertension1)=Indicator(1)
          /CONTRAST (CAD_history)=Indicator(1)
           /CONTRAST (stroke_history)=Indicator(1)
              /CONTRAST (PAOD)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Fat.bin_40 
/METHOD=ENTER FH_CVD_comp2 Age Gender OKyear Symptoms.4g
   /METHOD=ENTER BMI Stenosis_ipsilateral Risk614 Stroke_history PAOD 
  /CONTRAST (FH_CVD_comp2)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
      /CONTRAST (Symptoms.4g)=Indicator(1)
         /CONTRAST (Stenosis_ipsilateral)=Indicator(1)
          /CONTRAST (Risk614)=Indicator(1)
           /CONTRAST (stroke_history)=Indicator(1)
              /CONTRAST (PAOD)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES  Fat.bin_10
   /METHOD=ENTER FH_CVD_comp2 Age Gender OKyear Symptoms.4g
   /METHOD=ENTER Stroke_history PAOD med.statin.derived GFR_MDRD
  /CONTRAST (FH_CVD_comp2)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
      /CONTRAST (Symptoms.4g)=Indicator(1)
           /CONTRAST (stroke_history)=Indicator(1)
              /CONTRAST (PAOD)=Indicator(1)
               /CONTRAST (med.statin.derived)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES IPH.bin
 /METHOD=ENTER FH_CVD_comp2 Age Gender OKyear Symptoms.4g
   /METHOD=ENTER BMI Hypertension1 CAD_history Stroke_history LN_LDL LN_totalchol  
  /CONTRAST (FH_CVD_comp2)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
      /CONTRAST (Symptoms.4g)=Indicator(1)
         /CONTRAST (hypertension1)=Indicator(1)
          /CONTRAST (CAD_history)=Indicator(1)
           /CONTRAST (stroke_history)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*zelfde maar dan statine ipv LDL en totalchol (lipiden maakt wel meer sense).
LOGISTIC REGRESSION VARIABLES IPH.bin
 /METHOD=ENTER FH_CVD_comp2 Age Gender OKyear Symptoms.4g
   /METHOD=ENTER BMI Hypertension1 CAD_history Stroke_history med.statin.derived  
  /CONTRAST (FH_CVD_comp2)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
      /CONTRAST (Symptoms.4g)=Indicator(1)
         /CONTRAST (hypertension1)=Indicator(1)
          /CONTRAST (CAD_history)=Indicator(1)
           /CONTRAST (stroke_history)=Indicator(1)
                /CONTRAST (med.statin.derived)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES SMC.bin
  /METHOD=ENTER FH_CVD_comp2 Age Gender OKyear Symptoms.4g
   /METHOD=ENTER stroke_history risk614 GFR_MDRD
  /CONTRAST (FH_CVD_comp2)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
      /CONTRAST (Symptoms.4g)=Indicator(1)
         /CONTRAST (stroke_history)=Indicator(1)
          /CONTRAST (risk614)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES macrophages.bin
 /METHOD=ENTER FH_CVD_comp2 Age Gender OKyear Symptoms.4g
   /METHOD=ENTER Hypertension1 Stenosis_ipsilateral Stroke_history PAOD med.statin.derived
  /CONTRAST (FH_CVD_comp2)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
      /CONTRAST (Symptoms.4g)=Indicator(1)
         /CONTRAST (Stenosis_ipsilateral)=Indicator(1)
          /CONTRAST (Hypertension1)=Indicator(1)
           /CONTRAST (stroke_history)=Indicator(1)
              /CONTRAST (PAOD)=Indicator(1)
               /CONTRAST (med.statin.derived)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*zelfde maar dan TG ipv statine=used in paper.
LOGISTIC REGRESSION VARIABLES macrophages.bin
 /METHOD=ENTER FH_CVD_comp2 Age Gender OKyear Symptoms.4g
   /METHOD=ENTER Hypertension1 Stenosis_ipsilateral Stroke_history PAOD LN_triglyceriden
  /CONTRAST (FH_CVD_comp2)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
      /CONTRAST (Symptoms.4g)=Indicator(1)
         /CONTRAST (Stenosis_ipsilateral)=Indicator(1)
          /CONTRAST (Hypertension1)=Indicator(1)
           /CONTRAST (stroke_history)=Indicator(1)
              /CONTRAST (PAOD)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

******************************************************************************************************************************************************
*FH composite vs EP composite (secundary events) - univariate analysis

******************************************************************************************************************************************************.

KM ep_com_t_3years BY FH_CVD_comp2
  /STATUS=epcom.3years(1)
  /PRINT TABLE MEAN
  /PLOT SURVIVAL
  /TEST LOGRANK
  /COMPARE OVERALL POOLED.

**plaatje.
COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /STRATA=FH_CVD_comp2
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

**p-val from cox regression univariate.
COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /CONTRAST (FH_CVD_comp2)=Indicator(1)
  /METHOD=ENTER FH_CVD_comp2 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

*EP major and EP_cvdeath analysen- te weinig power.

******************************************************************************************************************************************************
*FH composite vs EP composite (secundary events) - Multivariate analysis.
******************************************************************************************************************************************************.

*traditional risk factors: Age, Gender, smoking, DM, BMI, hypertension1, hypercholesterolemie.

COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /CONTRAST (FH_CVD_comp2)=Indicator(1)
  /CONTRAST (risk614)=Indicator(1)
    /CONTRAST(Gender)=Indicator(1)
        /CONTRAST(SmokerCurrent)=Indicator(1)
            /CONTRAST(DM.composite)=Indicator(1)
                /CONTRAST(Hypertension1)=Indicator(1)
  /METHOD=ENTER FH_CVD_comp2 Age Gender risk614 SmokerCurrent DM.composite BMI Hypertension1
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

*Plus additional covariates: CAD, PAOD, Symptoms 4g en GFR= used model in paper. 

COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /CONTRAST (FH_CVD_comp2)=Indicator(1)
  /CONTRAST (risk614)=Indicator(1)
    /CONTRAST(Gender)=Indicator(1)
    /CONTRAST(SmokerCurrent)=Indicator(1)
     /CONTRAST(DM.composite)=Indicator(1)
     /CONTRAST(Hypertension1)=Indicator(1)
    /CONTRAST(CAD_history)=Indicator(1)
  /CONTRAST(PAOD)=Indicator(1)
   /CONTRAST(Symptoms.4g)=Indicator(1)  
  /METHOD=ENTER FH_CVD_comp2 Age Gender risk614 SmokerCurrent DM.composite BMI Hypertension1 CAD_history PAOD Symptoms.4g GFR_MDRD
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).



*Plus additional covariates: CAD, PAOD, Symptoms 4g en GFR= used model in paper + MetaGRS.

COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /CONTRAST (FH_CVD_comp2)=Indicator(1)
  /CONTRAST (AEGS_type)=Indicator(1)
  /CONTRAST (risk614)=Indicator(1)
    /CONTRAST(Gender)=Indicator(1)
    /CONTRAST(SmokerCurrent)=Indicator(1)
     /CONTRAST(DM.composite)=Indicator(1)
     /CONTRAST(Hypertension1)=Indicator(1)
    /CONTRAST(CAD_history)=Indicator(1)
  /CONTRAST(PAOD)=Indicator(1)
   /CONTRAST(Symptoms.4g)=Indicator(1)  
  /METHOD=ENTER FH_CVD_comp2 Age Gender risk614 SmokerCurrent DM.composite BMI Hypertension1 CAD_history PAOD Symptoms.4g GFR_MDRD CAD_PRSInouye AEGS_type PC1 PC2 PC3 PC4 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

 


*corrigeren voor LDL, TG  en totalchol ipv hypercholesterolemie.

COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /CONTRAST (FH_CVD_comp2)=Indicator(1)
    /CONTRAST(Gender)=Indicator(1)
    /CONTRAST(SmokerCurrent)=Indicator(1)
     /CONTRAST(DM.composite)=Indicator(1)
     /CONTRAST(Hypertension1)=Indicator(1)
    /CONTRAST(CAD_history)=Indicator(1)
  /CONTRAST(PAOD)=Indicator(1)
   /CONTRAST(Symptoms.4g)=Indicator(1)  
  /METHOD=ENTER FH_CVD_comp2 Age Gender SmokerCurrent DM.composite BMI Hypertension1 CAD_history PAOD Symptoms.4g LN_LDL LN_triglyceriden LN_totalchol
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).


*****************************************************************************************************************************************************
*INOUYE METAGRS ANALYSES 
*****************************************************************************************************************************************************

*****************************************************************************************************************************************************
**Baseline Inouye polygenic score.
*****************************************************************************************************************************************************.
* 1 rij met total cohort (hele gegenotypeerde cohort).
* dan lowest 20th vs highest 20th. 

**********************************************hele genotyped cohort ******************************************************************************************


CTABLES
  /VLABELS VARIABLES=Gender SmokerCurrent DM.composite risk614 Hypertension1 CAD_history Stroke_history PAOD 
    Symptoms.4g stenosis_ipsilateral stenosis_con_bin Med.statin.derived Med.all.antiplatelet YearORnew DISPLAY=DEFAULT 
  /TABLE Gender [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1]
+ SmokerCurrent [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] + DM.composite [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] 
+ risk614 [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] + Hypertension1 [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] 
+ CAD_history [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] + Stroke_history [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1]
+ PAOD [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1]
+ symptoms.4g [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] + Stenosis_ipsilateral [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] + stenosis_con_bin [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1]
+ Med.statin.derived [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] + Med.all.antiplatelet [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] 
+ YearORnew [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] 
  /SLABELS POSITION=ROW
  /CATEGORIES VARIABLES= Gender SmokerCurrent DM.composite risk614 Hypertension1 CAD_history Stroke_history PAOD Symptoms.4g stenosis_ipsilateral stenosis_con_bin Med.statin.derived 
  Med.all.antiplatelet YearORnew ORDER=D KEY=VALUE EMPTY=INCLUDE
  /TITLES
    TITLE='Table: Categorical characteristics at baseline, total cohort'.


*baseline tables with p-values-  continuous variables.

*continue variabelen: Age BMI GFR_MDRD totalchol  triglyceriden LDL HDL.
CTABLES
  /VLABELS VARIABLES=Age BMI GFR_MDRD totalchol triglyceriden LDL HDL DISPLAY=DEFAULT 
  /TABLE Age [S][MEAN F40.2, STDDEV 'SD' F40.2, MEDIAN F40.2] + BMI [S][MEAN F40.2, STDDEV 'SD' F40.2, MEDIAN F40.2] + GFR_MDRD [S][MEAN F40.2, STDDEV 'SD' F40.2, MEDIAN F40.2] 
  +  totalchol [S][MEAN F40.2, STDDEV 'SD' F40.2, MEDIAN F40.2] + triglyceriden [S][MEAN F40.2, STDDEV 'SD' F40.2, MEDIAN F40.2] 
  + LDL [S][MEAN F40.2, STDDEV 'SD' F40.2, MEDIAN F40.2] + HDL [S][MEAN F40.2, STDDEV 'SD' F40.2, MEDIAN F40.2] 
  /SLABELS POSITION=ROW
  /TITLES
    TITLE='Table 1: Continuous baseline characteristics , total cohort'.


*median and interquartile ranges.
EXAMINE VARIABLES=LDL
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

EXAMINE VARIABLES=HDL
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

EXAMINE VARIABLES=totalchol
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

EXAMINE VARIABLES=triglyceriden
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.


*voor totalen (gehele genotypeerd cohort).
FREQUENCIES VARIABLES=Gender SmokerCurrent DM.composite risk614 Hypertension1 CAD_history Stroke_history PAOD 
    Symptoms.4g Stenosis_ipsilateral stenosis_con_bin Med.statin.derived Med.all.antiplatelet YearORnew
  /ORDER=ANALYSIS.


******************************************Baselinecharacteristics compared with PRS Inouye Continuous******************************************************************************************
*needed to assess additional confounders for multivariate plaque and outcome analysis.

**Inouye vs CAD_history = Validation want dit moet wel*****************.

  LOGISTIC REGRESSION VARIABLES CAD_history
  /METHOD=ENTER CAD_PRSInouye
  /METHOD=ENTER Age
   /METHOD=ENTER Gender
   /METHOD=ENTER AEGS_type
   /METHOD=ENTER PC1 PC2 PC2 PC3
   /CONTRAST (Gender)=Indicator(1)
   /CONTRAST (AEGS_type)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

********************Correleert PRS_Inouye met klassieke risicofactoren?*****************************.

* dan correctie age, gender, chip, PC104. 


REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT CAD_PRSInouye
  /METHOD=ENTER Hypertension1
  /METHOD=ENTER Age Gender AEGS_type PC1 PC2 PC2 PC3 PC4.
   

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT CAD_PRSInouye
  /METHOD=ENTER  Hypertension.drugs
  /METHOD=ENTER Age Gender AEGS_type PC1 PC2 PC2 PC3 PC4. 


REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT CAD_PRSInouye
  /METHOD=ENTER  risk614
  /METHOD=ENTER Age Gender AEGS_type PC1 PC2 PC2 PC3 PC4. 

 
REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT CAD_PRSInouye
  /METHOD=ENTER  DM.composite
  /METHOD=ENTER Age Gender AEGS_type PC1 PC2 PC2 PC3 PC4.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT CAD_PRSInouye
  /METHOD=ENTER  SmokerCurrent
  /METHOD=ENTER Age Gender AEGS_type PC1 PC2 PC2 PC3 PC4.

 REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT CAD_PRSInouye
  /METHOD=ENTER  PAOD
   /METHOD=ENTER Age Gender AEGS_type PC1 PC2 PC2 PC3 PC4.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT CAD_PRSInouye
  /METHOD=ENTER  Stroke_history
   /METHOD=ENTER Age Gender AEGS_type PC1 PC2 PC2 PC3 PC4.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT CAD_PRSInouye
  /METHOD=ENTER  med.statin.derived
  /METHOD=ENTER Age Gender AEGS_type PC1 PC2 PC2 PC3 PC4.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT CAD_PRSInouye
  /METHOD=ENTER  med.all.antiplatelet
  /METHOD=ENTER Age Gender AEGS_type PC1 PC2 PC2 PC3 PC4.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT CAD_PRSInouye
  /METHOD=ENTER  BMI30
  /METHOD=ENTER Age Gender AEGS_type PC1 PC2 PC2 PC3 PC4.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT CAD_PRSInouye
  /METHOD=ENTER  BMI
   /METHOD=ENTER Age Gender AEGS_type PC1 PC2 PC2 PC3 PC4.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT CAD_PRSInouye
  /METHOD=ENTER  stenosis_con_bin
  /METHOD=ENTER Age Gender AEGS_type PC1 PC2 PC2 PC3 PC4.


REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT CAD_PRSInouye
  /METHOD=ENTER  ln_triglyceriden
  /METHOD=ENTER Age Gender AEGS_type PC1 PC2 PC2 PC3 PC4.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT CAD_PRSInouye
  /METHOD=ENTER  ln_HDL
  /METHOD=ENTER Age Gender AEGS_type PC1 PC2 PC2 PC3 PC4.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT CAD_PRSInouye
  /METHOD=ENTER  ln_LDL
  /METHOD=ENTER Age Gender AEGS_type PC1 PC2 PC2 PC3 PC4.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT CAD_PRSInouye
  /METHOD=ENTER  LN_totalchol
  /METHOD=ENTER Age Gender AEGS_type PC1 PC2 PC2 PC3 PC4.


REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT CAD_PRSInouye
  /METHOD=ENTER Age Gender AEGS_type PC1 PC2 PC2 PC3 PC4 Stenosis_ipsilateral.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT CAD_PRSInouye
  /METHOD=ENTER Age Gender AEGS_type PC1 PC2 PC2 PC3 PC4 GFR_MDRD.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT CAD_PRSInouye
  /METHOD=ENTER Age Gender AEGS_type PC1 PC2 PC2 PC3 PC4 OKyear.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT CAD_PRSInouye
  /METHOD=ENTER Age Gender AEGS_type PC1 PC2 PC2 PC3 PC4 Symptoms.4g.


*****************************************************************************************************************************************************.
** STANDARDIZATION OF METAGRS per SD 

*association with sCVE? Independent of FHx and IPH?
*association with plaque characteristics?
*gender stratified
*****************************************************************************************************************************************************.

* reviewer EHJ:

**************************************************************************************************************************************************************************.
*********************************************************ADDENDUM*********************************************************.

**Reviewer EHJ:
 
 *Genotyping was performed in two batches (which would also presumably have difference imputation efficiencies)
 did the metaGRS differ between the batches, adjusting for phenotype characteristics which also may differ? 
 If so, it would be advisable to standardise the metaGRS within each batch and merged. **.

****To make the metaGRS analyses more comparable with other studies, the metaGRS should be standardised to mean zero and unit variance***.

** dus per chip eerst standardizeren en dan mergen**.

*select cases: no missing for CAD_PRS en chiptype .
USE ALL.
COMPUTE filter_$=( ~ (SYSMIS(CAD_PRSInouye))).
VARIABLE LABELS filter_$ ' ~ (SYSMIS(CAD_PRSInouye)) (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.


*AEGS_type=0.
USE ALL.
COMPUTE filter_$=(AEGS_type = 0).
VARIABLE LABELS filter_$ 'AEGS_type = 0 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.


DESCRIPTIVES VARIABLES=CAD_PRSInouye /SAVE.
*rename: ZCAD_PRSInouye_chip0.


FILTER OFF.
USE ALL.
EXECUTE.


*AEGS_type= 1.
USE ALL.
COMPUTE filter_$=(AEGS_type = 1).
VARIABLE LABELS filter_$ 'AEGS_type = 1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

DESCRIPTIVES VARIABLES=CAD_PRSInouye /SAVE.
*rename: ZCAD_PRSInouye_chip1.

FILTER OFF.
USE ALL.
EXECUTE.

*merge standardized scores for chip 1 and chip 0.
COMPUTE MetaGRSstandard = ZCAD_PRSInouye_chip0.
IF (AEGS_type = 1) MetaGRSstandard=ZCAD_PRSInouye_chip1.
EXECUTE.

** MetaGRSstandard = goede variabelen= standaridezed MetaGRS variabele.


*********************************EP composite analysen Total Cohort************************************************************************************************************************************.

*chiptype is not in the model because this is included in MetaGRS standardization, so already corrected per chip.

*univariate.
COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /METHOD=ENTER MetaGRSstandard Age Gender PC1 PC2 PC3 PC4 
    /CONTRAST (Gender)=Indicator(1)
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).


*traditonal risk factors added.
COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /METHOD=ENTER MetaGRSstandard  Age Gender PC1 PC2 PC3 PC4 
      /METHOD=ENTER DM.composite risk614 Hypertension1 BMI SmokerCurrent
    /CONTRAST (Gender)=Indicator(1)
        /CONTRAST (DM.composite)=Indicator(1)
            /CONTRAST (risk614)=Indicator(1)
                /CONTRAST (SmokerCurrent)=Indicator(1)
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

*association independent from IPH?.
*traditonal risk factors added+ IPH.

COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /METHOD=ENTER MetaGRSstandard  Age Gender PC1 PC2 PC3 PC4 
      /METHOD=ENTER DM.composite risk614 Hypertension1 BMI SmokerCurrent IPH.bin
    /CONTRAST (Gender)=Indicator(1)
        /CONTRAST (DM.composite)=Indicator(1)
            /CONTRAST (risk614)=Indicator(1)
                /CONTRAST (SmokerCurrent)=Indicator(1)
                 /CONTRAST (IPH.bin)=Indicator(1)
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

*traditional risk factors including FHx.
COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /METHOD=ENTER MetaGRSstandard Age Gender PC1 PC2 PC3 PC4 
      /METHOD=ENTER DM.composite risk614 Hypertension1 BMI SmokerCurrent FH_CVD_comp2
    /CONTRAST (Gender)=Indicator(1)
        /CONTRAST (DM.composite)=Indicator(1)
            /CONTRAST (risk614)=Indicator(1)
                /CONTRAST (SmokerCurrent)=Indicator(1)
                    /CONTRAST (FH_CVD_comp2)=Indicator(1)
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).


*association independent from IPH?.
*traditional risk factors including FHx + Iph.bin.
COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /METHOD=ENTER MetaGRSstandard Age Gender PC1 PC2 PC3 PC4 
      /METHOD=ENTER DM.composite risk614 Hypertension1 BMI SmokerCurrent FH_CVD_comp2 IPH.bin
    /CONTRAST (Gender)=Indicator(1)
        /CONTRAST (DM.composite)=Indicator(1)
            /CONTRAST (risk614)=Indicator(1)
                /CONTRAST (SmokerCurrent)=Indicator(1)
                    /CONTRAST (FH_CVD_comp2)=Indicator(1)
                     /CONTRAST (IPH.bin)=Indicator(1)
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).


*********************************EP composite analysen- only Females************************************************************************************************************************************.

USE ALL.
COMPUTE filter_$=(Gender=0).
VARIABLE LABELS filter_$ 'Gender=0 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

*do not have correct for chiptype because this is included in MetaGRS standardization.
*also now no correction for gender as we only look at females.

*univariate.
COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /METHOD=ENTER MetaGRSstandard Age PC1 PC2 PC3 PC4 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).


*traditonal risk factors added.
COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /METHOD=ENTER MetaGRSstandard  Age PC1 PC2 PC3 PC4 
      /METHOD=ENTER DM.composite risk614 Hypertension1 BMI SmokerCurrent
        /CONTRAST (DM.composite)=Indicator(1)
            /CONTRAST (risk614)=Indicator(1)
                /CONTRAST (SmokerCurrent)=Indicator(1)
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).


*traditonal risk factors added+ IPH.
*association independent from IPH?.

COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /METHOD=ENTER MetaGRSstandard  Age  PC1 PC2 PC3 PC4 
      /METHOD=ENTER DM.composite risk614 Hypertension1 BMI SmokerCurrent IPH.bin
        /CONTRAST (DM.composite)=Indicator(1)
            /CONTRAST (risk614)=Indicator(1)
                /CONTRAST (SmokerCurrent)=Indicator(1)
                 /CONTRAST (IPH.bin)=Indicator(1)
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

*traditional risk factors including FHx.
COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /METHOD=ENTER MetaGRSstandard Age PC1 PC2 PC3 PC4 
      /METHOD=ENTER DM.composite risk614 Hypertension1 BMI SmokerCurrent FH_CVD_comp2
        /CONTRAST (DM.composite)=Indicator(1)
            /CONTRAST (risk614)=Indicator(1)
                /CONTRAST (SmokerCurrent)=Indicator(1)
                    /CONTRAST (FH_CVD_comp2)=Indicator(1)
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).


*association independent from IPH?.
*traditional risk factors including FHx + Iph.bin.
COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /METHOD=ENTER MetaGRSstandard Age PC1 PC2 PC3 PC4 
      /METHOD=ENTER DM.composite risk614 Hypertension1 BMI SmokerCurrent FH_CVD_comp2 IPH.bin
        /CONTRAST (DM.composite)=Indicator(1)
            /CONTRAST (risk614)=Indicator(1)
                /CONTRAST (SmokerCurrent)=Indicator(1)
                    /CONTRAST (FH_CVD_comp2)=Indicator(1)
                     /CONTRAST (IPH.bin)=Indicator(1)
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

FILTER OFF.
USE ALL.
EXECUTE.


*********************************EP composite analysen- Men***********************************************************************************************************************************.


USE ALL.
COMPUTE filter_$=(Gender=1).
VARIABLE LABELS filter_$ 'Gender=1 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.


*do not have correct for chiptype because this is included in MetaGRS standardization.

*univariate.
COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /METHOD=ENTER MetaGRSstandard Age PC1 PC2 PC3 PC4
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).


*traditonal risk factors added.
COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /METHOD=ENTER MetaGRSstandard  Age PC1 PC2 PC3 PC4 
      /METHOD=ENTER DM.composite risk614 Hypertension1 BMI SmokerCurrent
        /CONTRAST (DM.composite)=Indicator(1)
            /CONTRAST (risk614)=Indicator(1)
                /CONTRAST (SmokerCurrent)=Indicator(1)
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

*association independent from IPH?.
*traditonal risk factors added+ IPH.

COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /METHOD=ENTER MetaGRSstandard  Age PC1 PC2 PC3 PC4 
      /METHOD=ENTER DM.composite risk614 Hypertension1 BMI SmokerCurrent IPH.bin
        /CONTRAST (DM.composite)=Indicator(1)
            /CONTRAST (risk614)=Indicator(1)
                /CONTRAST (SmokerCurrent)=Indicator(1)
                 /CONTRAST (IPH.bin)=Indicator(1)
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

*traditional risk factors including FHx.
COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /METHOD=ENTER MetaGRSstandard Age PC1 PC2 PC3 PC4 
      /METHOD=ENTER DM.composite risk614 Hypertension1 BMI SmokerCurrent FH_CVD_comp2
        /CONTRAST (DM.composite)=Indicator(1)
            /CONTRAST (risk614)=Indicator(1)
                /CONTRAST (SmokerCurrent)=Indicator(1)
                    /CONTRAST (FH_CVD_comp2)=Indicator(1)
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).


*association independent from IPH?.
*traditional risk factors including FHx + Iph.bin.
COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /METHOD=ENTER MetaGRSstandard Age PC1 PC2 PC3 PC4 
      /METHOD=ENTER DM.composite risk614 Hypertension1 BMI SmokerCurrent FH_CVD_comp2 IPH.bin
        /CONTRAST (DM.composite)=Indicator(1)
            /CONTRAST (risk614)=Indicator(1)
                /CONTRAST (SmokerCurrent)=Indicator(1)
                    /CONTRAST (FH_CVD_comp2)=Indicator(1)
                     /CONTRAST (IPH.bin)=Indicator(1)
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

FILTER OFF.
USE ALL.
EXECUTE.


*******************************Plaque Characteristics-  Total Cohort************************************************************************************************************************************.

*plaque analyses in correlatie met Inouye PRS: met Age, Gender, PC1-4 en Chiptype als confounders.
*chiptype zit al in gestandardisserde score.

**plaque vulnerability score**.
REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT Instability 
  /METHOD=ENTER MetaGRSstandard Age Gender  PC1 PC2 PC3 PC4  
  /METHOD=ENTER MetaGRSstandard Age Gender  PC1 PC2 PC3 PC4 Symptoms.4g OKyear.


*Inouye80vs20 = upper quintile versus rest.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT Instability 
  /METHOD=ENTER Inouye80vs20 Age Gender  PC1 PC2 PC3 PC4  
  /METHOD=ENTER Inouye80vs20 Age Gender  PC1 PC2 PC3 PC4 Symptoms.4g OKyear.


*separate plaque characteristics.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_macmean
  /METHOD=ENTER MetaGRSstandard
  /METHOD=ENTER Age
  /METHOD=ENTER Gender
  /METHOD=ENTER PC1 PC2 PC3 PC4.


REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_smcmean
  /METHOD=ENTER MetaGRSstandard
  /METHOD=ENTER Age
  /METHOD=ENTER Gender
  /METHOD=ENTER PC1 PC2 PC3 PC4.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_vd
  /METHOD=ENTER MetaGRSstandard
  /METHOD=ENTER Age
  /METHOD=ENTER Gender
  /METHOD=ENTER PC1 PC2 PC3 PC4.


**binary.
LOGISTIC REGRESSION VARIABLES Collagen.bin
  /METHOD=ENTER MetaGRSstandard
  /METHOD=ENTER Age Gender 
   /METHOD=ENTER PC1 PC2 PC3 PC4
  /CONTRAST (Gender)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES calc.bin
  /METHOD=ENTER MetaGRSstandard
  /METHOD=ENTER Age Gender 
   /METHOD=ENTER PC1 PC2 PC3 PC4
  /CONTRAST (Gender)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Fat.bin_40
  /METHOD=ENTER MetaGRSstandard
  /METHOD=ENTER Age Gender 
   /METHOD=ENTER PC1 PC2 PC3 PC4
  /CONTRAST (Gender)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Fat.bin_10
  /METHOD=ENTER MetaGRSstandard
  /METHOD=ENTER Age Gender 
   /METHOD=ENTER PC1 PC2 PC3 PC4
  /CONTRAST (Gender)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


LOGISTIC REGRESSION VARIABLES iph.bin
  /METHOD=ENTER MetaGRSstandard
  /METHOD=ENTER Age Gender 
   /METHOD=ENTER PC1 PC2 PC3 PC4
  /CONTRAST (Gender)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


LOGISTIC REGRESSION VARIABLES smc.bin
  /METHOD=ENTER MetaGRSstandard
  /METHOD=ENTER Age Gender 
   /METHOD=ENTER PC1 PC2 PC3 PC4
  /CONTRAST (Gender)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES macrophages.bin
  /METHOD=ENTER MetaGRSstandard
  /METHOD=ENTER Age Gender 
   /METHOD=ENTER PC1 PC2 PC3 PC4
  /CONTRAST (Gender)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


*MVA*****.

*pragmatisch univariaat model Age, Gender, PC1-4 + OK year + Symptoms.4g.
* daaraan toevoegen confounders p<0.2 met Inouye_CAD score: dit zijn HDL, LLD, LDL, hyperchol (zie confoundertable).


*macmean: alleen correlatie met  LLD.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_macmean
  /METHOD=ENTER MetaGRSstandard Age Gender  PC1 PC2 PC3 PC4
  /METHOD=ENTER OKyear Symptoms.4g Med.Statin.LLD.

*smc mean: hyperchol.
REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_smcmean
  /METHOD=ENTER MetaGRSstandard Age Gender  AEGS_type PC1 PC2 PC3 PC4
  /METHOD=ENTER OKyear Symptoms.4g risk614.


*vd: Geen extra confounders dus alleen OK year en symptoms toevoegen=pragmatisch. 
REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_vd
  /METHOD=ENTER MetaGRSstandard Age Gender  AEGS_type PC1 PC2 PC3 PC4
  /METHOD=ENTER OKyear Symptoms.4g.


**binary************.

*collagen: + hyperchol.
LOGISTIC REGRESSION VARIABLES Collagen.bin
  /METHOD=ENTER MetaGRSstandard Age Gender PC1 PC2 PC3 PC4
   /METHOD=ENTER OKyear Symptoms.4g risk614
  /CONTRAST (Gender)=Indicator(1)
          /CONTRAST (Symptoms.4g)=Indicator(1)
               /CONTRAST (risk614)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*calc:  + HDL.
LOGISTIC REGRESSION VARIABLES calc.bin
  /METHOD=ENTER MetaGRSstandard Age Gender PC1 PC2 PC3 PC4
   /METHOD=ENTER OKyear Symptoms.4g LN_HDL
  /CONTRAST (Gender)=Indicator(1)
          /CONTRAST (Symptoms.4g)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*Fat 40: +hyperchol, +statin, meerdere varianten mogelijk.

*met hyperchol.
LOGISTIC REGRESSION VARIABLES Fat.bin_40
   /METHOD=ENTER MetaGRSstandard Age Gender PC1 PC2 PC3 PC4
   /METHOD=ENTER OKyear Symptoms.4g risk614
  /CONTRAST (Gender)=Indicator(1)
        /CONTRAST (risk614)=Indicator(1)
          /CONTRAST (Symptoms.4g)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*met  LLD.
LOGISTIC REGRESSION VARIABLES Fat.bin_40
  /METHOD=ENTER MetaGRSstandard Age Gender PC1 PC2 PC3 PC4
   /METHOD=ENTER OKyear Symptoms.4g Med.Statin.LLD
  /CONTRAST (Gender)=Indicator(1)
          /CONTRAST (Symptoms.4g)=Indicator(1)
                /CONTRAST (Med.Statin.LLD)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


*fat 10, +statin.
LOGISTIC REGRESSION VARIABLES Fat.bin_10
  /METHOD=ENTER MetaGRSstandard Age Gender PC1 PC2 PC3 PC4
   /METHOD=ENTER OKyear Symptoms.4g Med.Statin.LLD
  /CONTRAST (Gender)=Indicator(1)
          /CONTRAST (Symptoms.4g)=Indicator(1)
                /CONTRAST (Med.Statin.LLD)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*iph.bin *met LLD.
LOGISTIC REGRESSION VARIABLES iph.bin
 /METHOD=ENTER MetaGRSstandard Age Gender PC1 PC2 PC3 PC4
   /METHOD=ENTER OKyear Symptoms.4g Med.Statin.LLD
  /CONTRAST (Gender)=Indicator(1)
          /CONTRAST (Symptoms.4g)=Indicator(1)
                /CONTRAST (Med.Statin.LLD)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES iph.bin
 /METHOD=ENTER MetaGRSstandard Age Gender PC1 PC2 PC3 PC4
   /METHOD=ENTER OKyear Symptoms.4g Med.Statin.derived
  /CONTRAST (Gender)=Indicator(1)
          /CONTRAST (Symptoms.4g)=Indicator(1)
                /CONTRAST (Med.Statin.derived)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


*smc bin,  + hyperchol of LDL.

*met hyperchol.
LOGISTIC REGRESSION VARIABLES smc.bin
  /METHOD=ENTER MetaGRSstandard Age Gender PC1 PC2 PC3 PC4
   /METHOD=ENTER OKyear Symptoms.4g risk614
  /CONTRAST (Gender)=Indicator(1)
          /CONTRAST (Symptoms.4g)=Indicator(1)
                /CONTRAST (risk614)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


*macrophages: met LLD.

LOGISTIC REGRESSION VARIABLES macrophages.bin
  /METHOD=ENTER MetaGRSstandard Age Gender PC1 PC2 PC3 PC4
   /METHOD=ENTER OKyear Symptoms.4g Med.Statin.LLD
  /CONTRAST (Gender)=Indicator(1)
          /CONTRAST (Symptoms.4g)=Indicator(1)
                /CONTRAST (Med.Statin.LLD)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


*******************************Plaque Characteristics-  gender stratified**********.

SORT CASES  BY Gender.
SPLIT FILE LAYERED BY Gender.


REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_macmean
  /METHOD=ENTER MetaGRSstandard
  /METHOD=ENTER Age
  /METHOD=ENTER PC1 PC2 PC3 PC4.


REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_smcmean
  /METHOD=ENTER MetaGRSstandard
  /METHOD=ENTER Age
  /METHOD=ENTER PC1 PC2 PC3 PC4.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_vd
  /METHOD=ENTER MetaGRSstandard
  /METHOD=ENTER Age
  /METHOD=ENTER PC1 PC2 PC3 PC4.


**binary.
LOGISTIC REGRESSION VARIABLES Collagen.bin
  /METHOD=ENTER MetaGRSstandard
  /METHOD=ENTER Age 
   /METHOD=ENTER PC1 PC2 PC3 PC4
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES calc.bin
  /METHOD=ENTER MetaGRSstandard
  /METHOD=ENTER Age 
   /METHOD=ENTER PC1 PC2 PC3 PC4
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Fat.bin_40
  /METHOD=ENTER MetaGRSstandard
  /METHOD=ENTER Age 
   /METHOD=ENTER PC1 PC2 PC3 PC4
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Fat.bin_10
  /METHOD=ENTER MetaGRSstandard
  /METHOD=ENTER Age 
   /METHOD=ENTER PC1 PC2 PC3 PC4
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


LOGISTIC REGRESSION VARIABLES iph.bin
  /METHOD=ENTER MetaGRSstandard
  /METHOD=ENTER Age 
   /METHOD=ENTER PC1 PC2 PC3 PC4
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


LOGISTIC REGRESSION VARIABLES smc.bin
  /METHOD=ENTER MetaGRSstandard
  /METHOD=ENTER Age
   /METHOD=ENTER PC1 PC2 PC3 PC4
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES macrophages.bin
  /METHOD=ENTER MetaGRSstandard
  /METHOD=ENTER Age 
   /METHOD=ENTER PC1 PC2 PC3 PC4
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


*MVA*.

*macmean: alleen correlatie met  LLD.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_macmean
  /METHOD=ENTER MetaGRSstandard Age  PC1 PC2 PC3 PC4
  /METHOD=ENTER OKyear Symptoms.4g Med.Statin.LLD.

*smc mean: hyperchol.
REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_smcmean
  /METHOD=ENTER MetaGRSstandard Age AEGS_type PC1 PC2 PC3 PC4
  /METHOD=ENTER OKyear Symptoms.4g risk614.

*vd: Geen extra confounders dus alleen OK year en symptoms toevoegen. 
REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_vd
  /METHOD=ENTER MetaGRSstandard Age AEGS_type PC1 PC2 PC3 PC4
  /METHOD=ENTER OKyear Symptoms.4g.


**binary************.

*collagen: + hyperchol.
LOGISTIC REGRESSION VARIABLES Collagen.bin
  /METHOD=ENTER MetaGRSstandard Age PC1 PC2 PC3 PC4
   /METHOD=ENTER OKyear Symptoms.4g risk614
          /CONTRAST (Symptoms.4g)=Indicator(1)
               /CONTRAST (risk614)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*calc:  + HDL.
LOGISTIC REGRESSION VARIABLES calc.bin
  /METHOD=ENTER MetaGRSstandard Age PC1 PC2 PC3 PC4
   /METHOD=ENTER OKyear Symptoms.4g LN_HDL
          /CONTRAST (Symptoms.4g)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*Fat 40: +hyperchol, +statin, meerdere varianten mogelijk.

*met hyperchol.
*LOGISTIC REGRESSION VARIABLES Fat.bin_40
   /METHOD=ENTER MetaGRSstandard Age PC1 PC2 PC3 PC4
   /METHOD=ENTER OKyear Symptoms.4g risk614
        /CONTRAST (risk614)=Indicator(1)
          /CONTRAST (Symptoms.4g)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*met  LLD.
*LOGISTIC REGRESSION VARIABLES Fat.bin_40
  /METHOD=ENTER MetaGRSstandard Age PC1 PC2 PC3 PC4
   /METHOD=ENTER OKyear Symptoms.4g Med.Statin.LLD
          /CONTRAST (Symptoms.4g)=Indicator(1)
                /CONTRAST (Med.Statin.LLD)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


*fat 10, +statin.
LOGISTIC REGRESSION VARIABLES Fat.bin_10
  /METHOD=ENTER MetaGRSstandard Age PC1 PC2 PC3 PC4
   /METHOD=ENTER OKyear Symptoms.4g Med.Statin.LLD
          /CONTRAST (Symptoms.4g)=Indicator(1)
                /CONTRAST (Med.Statin.LLD)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*iph.bin *met LLD.
LOGISTIC REGRESSION VARIABLES iph.bin
 /METHOD=ENTER MetaGRSstandard Age PC1 PC2 PC3 PC4
   /METHOD=ENTER OKyear Symptoms.4g Med.Statin.LLD
          /CONTRAST (Symptoms.4g)=Indicator(1)
                /CONTRAST (Med.Statin.LLD)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


*smc bin,  + hyperchol of LDL.

*met hyperchol.
LOGISTIC REGRESSION VARIABLES smc.bin
  /METHOD=ENTER MetaGRSstandard Age PC1 PC2 PC3 PC4
   /METHOD=ENTER OKyear Symptoms.4g risk614
          /CONTRAST (Symptoms.4g)=Indicator(1)
                /CONTRAST (risk614)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


*macrophages: met LLD.

LOGISTIC REGRESSION VARIABLES macrophages.bin
  /METHOD=ENTER MetaGRSstandard Age PC1 PC2 PC3 PC4
   /METHOD=ENTER OKyear Symptoms.4g Med.Statin.LLD
          /CONTRAST (Symptoms.4g)=Indicator(1)
                /CONTRAST (Med.Statin.LLD)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).



*****************************************************************************************************************************************************
*****************************************************************************************************************************************************
*Inouye- percentiles- lowest quintile compared to highest quintile.

*baseline table
*association with sCVE
*association with plaque characteristics.
*****************************************************************************************************************************************************.

*****************************************************************************************************************************************************
***********************************************baseline characteristics divided of 1st quintile PRS and 5th quintile******************************************************************************************.

* 1st quintile vs 5th quintile.

*baseline tables with p-values-  categorical variables.
*Gender SmokerCurrent DM.composite risk614 Hypertension1 CAD_history Stroke_history PAOD 
    Symptoms.4g stenosis_ipsilateral stenosis_con_bin Med.statin.derived Med.all.antiplatelet YearORnew.

CTABLES
  /VLABELS VARIABLES=Gender SmokerCurrent DM.composite risk614 Hypertension1 CAD_history Stroke_history PAOD 
    Symptoms.4g stenosis_ipsilateral stenosis_con_bin Med.statin.derived Med.all.antiplatelet YearORnew Inouye20th DISPLAY=DEFAULT 
  /TABLE Gender [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1]
+ SmokerCurrent [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] + DM.composite [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] 
+ risk614 [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] + Hypertension1 [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] 
+ CAD_history [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] + Stroke_history [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1]
+ PAOD [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1]
+ symptoms.4g [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] + Stenosis_ipsilateral [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] + stenosis_con_bin [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1]
+ Med.statin.derived [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] + Med.all.antiplatelet [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1] 
+ YearORnew [C][COUNT 'N' F40.0, LAYERCOLPCT.VALIDN ' (%)' PCTPAREN40.1]  BY Inouye20th [C]
  /SLABELS POSITION=ROW
  /CATEGORIES VARIABLES= Gender SmokerCurrent DM.composite risk614 Hypertension1 CAD_history Stroke_history PAOD Symptoms.4g stenosis_ipsilateral stenosis_con_bin Med.statin.derived 
  Med.all.antiplatelet YearORnew Inouye20th ORDER=D KEY=VALUE EMPTY=INCLUDE
  /TITLES
    TITLE='Table: Categorical characteristics at baseline, by Inouye20th'
    /SIGTEST TYPE=CHISQUARE ALPHA=0.05 INCLUDEMRSETS=YES CATEGORIES=ALLVISIBLE.

*baseline tables with p-values-  continuous variables.

*continue variabelen: Age BMI GFR_MDRD totalchol  triglyceriden LDL HDL.
CTABLES
  /VLABELS VARIABLES=Age BMI GFR_MDRD totalchol triglyceriden LDL HDL Inouye20th DISPLAY=DEFAULT 
  /TABLE Age [S][MEAN F40.2, STDDEV 'SD' F40.2, MEDIAN F40.2] + BMI [S][MEAN F40.2, STDDEV 'SD' F40.2, MEDIAN F40.2] + GFR_MDRD [S][MEAN F40.2, STDDEV 'SD' F40.2, MEDIAN F40.2] 
  +  totalchol [S][MEAN F40.2, STDDEV 'SD' F40.2, MEDIAN F40.2] + triglyceriden [S][MEAN F40.2, STDDEV 'SD' F40.2, MEDIAN F40.2] 
  + LDL [S][MEAN F40.2, STDDEV 'SD' F40.2, MEDIAN F40.2] + HDL [S][MEAN F40.2, STDDEV 'SD' F40.2, MEDIAN F40.2] BY Inouye20th [C]
  /SLABELS POSITION=ROW
  /CATEGORIES VARIABLES= Inouye20th ORDER=D KEY=VALUE EMPTY=INCLUDE
  /TITLES
    TITLE='Table 1: Continuous baseline characteristics , by Inouye20th'.


*for median and interquartile range.

EXAMINE VARIABLES=totalchol BY Inouye20th
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

EXAMINE VARIABLES=triglyceriden BY Inouye20th
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

EXAMINE VARIABLES=LDL BY Inouye20th
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

EXAMINE VARIABLES=HDL BY Inouye20th
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

*voor lipiden eerst LN transformatie.
T-TEST GROUPS=Inouye20th(0 1)
  /MISSING=ANALYSIS
  /VARIABLES=Age BMI GFR_MDRD LN_HDL LN_totalchol LN_triglyceriden LN_LDL
  /CRITERIA=CI(.95).

*Uitgebreider met p-waardes en missing values.

  CROSSTABS
  /TABLES=Gender SmokerCurrent DM.composite risk614 Hypertension1 CAD_history Stroke_history PAOD 
    Symptoms.4g Stenosis_ipsilateral stenosis_con_bin Med.statin.derived Med.all.antiplatelet YearORnew BY Inouye20th
  /FORMAT=AVALUE TABLES
  /STATISTICS=CHISQ 
  /CELLS=COUNT COLUMN 
  /COUNT ROUND CELL.

*select alleen cases Inouye20th.
USE ALL.
COMPUTE filter_$=( ~ (SYSMIS(Inouye20th))).
VARIABLE LABELS filter_$ ' ~ (SYSMIS(Inouye20th)) (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

CROSSTABS
  /TABLES= YearORnew BY AEGS_type
  /FORMAT=AVALUE TABLES
  /STATISTICS=CHISQ 
  /CELLS=COUNT COLUMN 
  /COUNT ROUND CELL.


SORT CASES  BY AEGS_type.
SPLIT FILE SEPARATE BY AEGS_type.

CROSSTABS
  /TABLES= YearORnew BY Inouye20th
  /FORMAT=AVALUE TABLES
  /STATISTICS=CHISQ 
  /CELLS=COUNT COLUMN 
  /COUNT ROUND CELL.


SPLIT FILE OFF.
FILTER OFF.
USE ALL.
EXECUTE.

*chiptype over beide groepen anders?.
CROSSTABS
  /TABLES= AEGS_type BY Inouye20th
  /FORMAT=AVALUE TABLES
  /STATISTICS=CHISQ 
  /CELLS=COUNT COLUMN 
  /COUNT ROUND CELL.


*****************************************************************************************************************************************************
*Inouye 20th  vs secondary outcome measures (Composite events). 
******************************************************************************************************************************************************

********************************************************Inouye 1 vs 5 quintile**********************************************************************************************.

KM ep_com_t_3years BY  Inouye20th
  /STATUS=epcom.3years(1)
  /PRINT MEAN
  /PLOT SURVIVAL
  /TEST LOGRANK
  /COMPARE OVERALL POOLED.

COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /METHOD=ENTER Inouye20th
    /CONTRAST (Inouye20th)=Indicator(1)
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

*gecorrigeerd voor Age Gender.
COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /METHOD=ENTER Inouye20th Age Gender
   /CONTRAST (Inouye20th)=Indicator(1)
    /CONTRAST (Gender)=Indicator(1)
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).


*multivariate- traditional risk factors added.
COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /METHOD=ENTER Inouye20th Age Gender 
  /METHOD=ENTER SmokerCurrent Hypertension1 risk614 DM.composite BMI 
   /CONTRAST (Inouye20th)=Indicator(1)
    /CONTRAST (Gender)=Indicator(1)
        /CONTRAST (SmokerCurrent)=Indicator(1)
      /CONTRAST (Hypertension1)=Indicator(1)
       /CONTRAST (risk614)=Indicator(1)
        /CONTRAST (DM.composite)=Indicator(1) 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).



*traditional risk factors + FH composite.
COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /METHOD=ENTER Inouye20th Age Gender 
  /METHOD=ENTER SmokerCurrent Hypertension1 risk614 DM.composite BMI FH_CVD_comp2
   /CONTRAST (Inouye20th)=Indicator(1)
    /CONTRAST (Gender)=Indicator(1)
        /CONTRAST (SmokerCurrent)=Indicator(1)
      /CONTRAST (Hypertension1)=Indicator(1)
       /CONTRAST (risk614)=Indicator(1)
        /CONTRAST (DM.composite)=Indicator(1) 
        /CONTRAST (FH_CVD_comp2)=Indicator(1) 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).


******************************************************************************************************************************************************
*inouye 1vs5 quintile vs plaque characteristics
******************************************************************************************************************************************************.

************************************Univariate***************************.


EXAMINE VARIABLES=macmean0 smcmean0 vessel_density_averaged Mast_cells_plaque neutrophils BY Inouye20th
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING PAIRWISE
  /NOTOTAL.

*check distributions. 
EXAMINE VARIABLES=LN_macmean LN_smcmean LN_vd LN_ms LN_neutro BY Inouye20th
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING PAIRWISE
  /NOTOTAL.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_macmean
  /METHOD=ENTER Inouye20th.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_smcmean
  /METHOD=ENTER Inouye20th.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_vd
  /METHOD=ENTER Inouye20th.


**binary.

CROSSTABS
  /TABLES=Calc.bin Collagen.bin Fat.bin_40 Fat.bin_10 IPH.bin Macrophages.bin SMC.bin  BY Inouye20th
  /FORMAT=AVALUE TABLES
  /STATISTICS=CHISQ 
  /CELLS=COUNT COLUMN 
  /COUNT ROUND CELL.

LOGISTIC REGRESSION VARIABLES Calc.bin
  /METHOD=ENTER Inouye20th 
  /CONTRAST (Inouye20th)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Collagen.bin
  /METHOD=ENTER Inouye20th 
  /CONTRAST (Inouye20th)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Fat.bin_40
  /METHOD=ENTER Inouye20th 
  /CONTRAST (Inouye20th)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Fat.bin_10
  /METHOD=ENTER Inouye20th 
  /CONTRAST (Inouye20th)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES IPH.bin
  /METHOD=ENTER Inouye20th 
  /CONTRAST (Inouye20th)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Macrophages.bin 
  /METHOD=ENTER Inouye20th 
  /CONTRAST (Inouye20th)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES SMC.bin
  /METHOD=ENTER Inouye20th 
  /CONTRAST (Inouye20th)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


*univariate, but corrected for age and gender.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_macmean 
  /METHOD=ENTER Inouye20th
  /METHOD=ENTER Age Gender.


REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_smcmean 
  /METHOD=ENTER Inouye20th
  /METHOD=ENTER Age Gender.


REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_VD
  /METHOD=ENTER Inouye20th
  /METHOD=ENTER Age Gender.


LOGISTIC REGRESSION VARIABLES Calc.bin
  /METHOD=ENTER Inouye20th Age Gender
  /CONTRAST (Inouye20th)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Collagen.bin
    /METHOD=ENTER Inouye20th Age Gender
  /CONTRAST (Inouye20th)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Fat.bin_40
   /METHOD=ENTER Inouye20th Age Gender 
  /CONTRAST (Inouye20th)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Fat.bin_10
    /METHOD=ENTER Inouye20th Age Gender
  /CONTRAST (Inouye20th)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES IPH.bin
   /METHOD=ENTER Inouye20th Age Gender
  /CONTRAST (Inouye20th)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Macrophages.bin 
    /METHOD=ENTER Inouye20th Age Gender
  /CONTRAST (Inouye20th)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES SMC.bin
    /METHOD=ENTER Inouye20th Age Gender 
  /CONTRAST (Inouye20th)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).



*******************************Multivariate only age, gender,  OR year, symptoms.4g + baseline differences p<0.2*************************.

*only pragmatic.
REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_macmean 
  /METHOD=ENTER Inouye20th
  /METHOD=ENTER Age Gender ORyear Symptoms.4g.

*+ p<0.2 baseline characteristics.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_macmean 
  /METHOD=ENTER Inouye20th
  /METHOD=ENTER Age Gender ORyear Symptoms.4g med.statin.derived DM.composite Stenosis_con_bin.

*only pragmatic.
REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_smcmean 
  /METHOD=ENTER Inouye20th
  /METHOD=ENTER Age Gender ORyear Symptoms.4g.

* p<0.2 added.
REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_smcmean 
  /METHOD=ENTER Inouye20th
  /METHOD=ENTER Age Gender ORyear Symptoms.4g stenosis_con_bin risk614.

*VD geen andere confounders.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_VD
  /METHOD=ENTER Inouye20th
  /METHOD=ENTER Age Gender ORyear Symptoms.4g.


*pragmatic.
LOGISTIC REGRESSION VARIABLES Calc.bin
  /METHOD=ENTER Inouye20th Age Gender ORyear Symptoms.4g
  /CONTRAST (Inouye20th)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*with P<0.2.
LOGISTIC REGRESSION VARIABLES Calc.bin
  /METHOD=ENTER Inouye20th Age Gender ORyear Symptoms.4g LN_totalchol
  /CONTRAST (Inouye20th)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


*pragmatic.
LOGISTIC REGRESSION VARIABLES Collagen.bin
    /METHOD=ENTER Inouye20th Age Gender ORyear Symptoms.4g
  /CONTRAST (Inouye20th)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*with covariates p<0.2.

LOGISTIC REGRESSION VARIABLES Collagen.bin
    /METHOD=ENTER Inouye20th Age Gender ORyear Symptoms.4g risk614 stenosis_con_bin
  /CONTRAST (Inouye20th)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
        /CONTRAST (risk614)=Indicator(1)
            /CONTRAST (stenosis_con_bin)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*pragmatic.
LOGISTIC REGRESSION VARIABLES Fat.bin_40
   /METHOD=ENTER Inouye20th Age Gender ORyear Symptoms.4g
  /CONTRAST (Inouye20th)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*with p<0.2 added (let op max 7 covariates).

*with statine.
LOGISTIC REGRESSION VARIABLES Fat.bin_40
   /METHOD=ENTER Inouye20th Age Gender ORyear Symptoms.4g DM.Composite Stenosis_con_bin med.statin.derived
  /CONTRAST (Inouye20th)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
     /CONTRAST (DM.composite)=Indicator(1)
          /CONTRAST (stenosis_con_bin)=Indicator(1)
             /CONTRAST (med.statin.derived)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*with risk614 in case of statine.
LOGISTIC REGRESSION VARIABLES Fat.bin_40
   /METHOD=ENTER Inouye20th Age Gender ORyear Symptoms.4g DM.Composite Stenosis_con_bin risk614
  /CONTRAST (Inouye20th)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
     /CONTRAST (DM.composite)=Indicator(1)
          /CONTRAST (stenosis_con_bin)=Indicator(1)
             /CONTRAST (risk614)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*pragmatic.
LOGISTIC REGRESSION VARIABLES Fat.bin_10
    /METHOD=ENTER Inouye20th Age Gender ORyear Symptoms.4g
  /CONTRAST (Inouye20th)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*p<0.2 added.
LOGISTIC REGRESSION VARIABLES Fat.bin_10
    /METHOD=ENTER Inouye20th Age Gender ORyear Symptoms.4g DM.Composite Stenosis_con_bin med.statin.derived
  /CONTRAST (Inouye20th)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
     /CONTRAST (DM.composite)=Indicator(1)
          /CONTRAST (stenosis_con_bin)=Indicator(1)
             /CONTRAST (med.statin.derived)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*pragmatic.
LOGISTIC REGRESSION VARIABLES IPH.bin
   /METHOD=ENTER Inouye20th Age Gender ORyear Symptoms.4g
  /CONTRAST (Inouye20th)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


*with p<0.2 with LDL and total chol.

LOGISTIC REGRESSION VARIABLES IPH.bin
   /METHOD=ENTER Inouye20th Age Gender ORyear Symptoms.4g LN_LDL LN_Totalchol
  /CONTRAST (Inouye20th)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


*empirically.
LOGISTIC REGRESSION VARIABLES Macrophages.bin 
    /METHOD=ENTER Inouye20th Age Gender ORyear Symptoms.4g
  /CONTRAST (Inouye20th)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*with p<0.2- statine.
LOGISTIC REGRESSION VARIABLES Macrophages.bin 
    /METHOD=ENTER Inouye20th Age Gender ORyear Symptoms.4g med.statin.derived
  /CONTRAST (Inouye20th)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
     /CONTRAST (med.statin.derived)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


*empircally.
LOGISTIC REGRESSION VARIABLES SMC.bin
    /METHOD=ENTER Inouye20th Age Gender ORyear Symptoms.4g
  /CONTRAST (Inouye20th)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*added p<0.2 with risk614.
LOGISTIC REGRESSION VARIABLES SMC.bin
    /METHOD=ENTER Inouye20th Age Gender ORyear Symptoms.4g stenosis_con_bin risk614
  /CONTRAST (Inouye20th)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
     /CONTRAST (stenosis_con_bin)=Indicator(1)
      /CONTRAST (risk614)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*****************************************************************************************************************************************************.
*****************************************************************************************************************************************************.
*****************************************************************************************************************************************************
*Inouye- percentiles- lowest quintile compared to rest (lower 80%).

*baseline table
*sCVE
*plaque characteristics.

*variable = Inouye80vs20.
*****************************************************************************************************************************************************.
*****************************************************************************************************************************************************.
**************************************BASELINE TABLE************************************************.


*for median and interquartile range.

EXAMINE VARIABLES=totalchol BY Inouye80vs20
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

EXAMINE VARIABLES=triglyceriden BY Inouye80vs20
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

EXAMINE VARIABLES=LDL BY Inouye80vs20
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

EXAMINE VARIABLES=HDL BY Inouye80vs20
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

*voor lipiden eerst LN transformatie.
T-TEST GROUPS=Inouye80vs20(0 1)
  /MISSING=ANALYSIS
  /VARIABLES=Age BMI GFR_MDRD LN_HDL LN_totalchol LN_triglyceriden LN_LDL
  /CRITERIA=CI(.95).


  CROSSTABS
  /TABLES=Gender SmokerCurrent DM.composite risk614 Hypertension1 CAD_history Stroke_history PAOD 
    Symptoms.4g Stenosis_ipsilateral stenosis_con_bin Med.statin.derived Med.all.antiplatelet YearORnew BY Inouye80vs20
  /FORMAT=AVALUE TABLES
  /STATISTICS=CHISQ 
  /CELLS=COUNT COLUMN 
  /COUNT ROUND CELL.


*******************************************Inouye  upper 20th vs lowest 80th and EP Composite*********************************************************.

KM ep_com_t_3years BY  Inouye80vs20
  /STATUS=epcom.3years(1)
  /PRINT MEAN
  /PLOT SURVIVAL
  /TEST LOGRANK
  /COMPARE OVERALL POOLED.

COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /METHOD=ENTER  Inouye80vs20
  /METHOD=ENTER Age Gender
    /CONTRAST (Gender)=Indicator(1)
   /CONTRAST (Inouye80vs20)=Indicator(1)
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).

COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /METHOD=ENTER  Inouye80vs20 Age Gender SmokerCurrent Hypertension1 risk614 DM.composite BMI 
   /CONTRAST (Inouye80vs20)=Indicator(1)
    /CONTRAST (Gender)=Indicator(1)
     /CONTRAST (SmokerCurrent)=Indicator(1)
      /CONTRAST (Hypertension1)=Indicator(1)
       /CONTRAST (risk614)=Indicator(1)
        /CONTRAST (DM.composite)=Indicator(1) 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).


*alle risicofactoren en nu met FH composite.
COXREG ep_com_t_3years
  /STATUS=epcom.3years(1)
  /METHOD=ENTER  Inouye80vs20 Age Gender SmokerCurrent Hypertension1 risk614 DM.composite BMI FH_CVD_comp2
   /CONTRAST (Inouye80vs20)=Indicator(1)
    /CONTRAST (Gender)=Indicator(1)
     /CONTRAST (SmokerCurrent)=Indicator(1)
      /CONTRAST (Hypertension1)=Indicator(1)
       /CONTRAST (risk614)=Indicator(1)
        /CONTRAST (DM.composite)=Indicator(1) 
        /CONTRAST (FH_CVD_comp2)=Indicator(1) 
  /PLOT SURVIVAL
  /PRINT=CI(95)
  /CRITERIA=PIN(.05) POUT(.10) ITERATE(20).



******************************************************************************************************************************************************
*inouye top 20 vs bottom 80 vs plaque characteristics
******************************************************************************************************************************************************.

************************************Univariate***************************.


EXAMINE VARIABLES=macmean0 smcmean0 vessel_density_averaged Mast_cells_plaque neutrophils BY Inouye80vs20
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING PAIRWISE
  /NOTOTAL.

*check distributions. 
EXAMINE VARIABLES=LN_macmean LN_smcmean LN_vd LN_ms LN_neutro BY Inouye80vs20
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING PAIRWISE
  /NOTOTAL.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_macmean
  /METHOD=ENTER Inouye80vs20.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_smcmean
  /METHOD=ENTER Inouye80vs20.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_vd
  /METHOD=ENTER Inouye80vs20.


**binary.

CROSSTABS
  /TABLES=Calc.bin Collagen.bin Fat.bin_40 Fat.bin_10 IPH.bin Macrophages.bin SMC.bin  BY Inouye80vs20
  /FORMAT=AVALUE TABLES
  /STATISTICS=CHISQ 
  /CELLS=COUNT COLUMN 
  /COUNT ROUND CELL.

LOGISTIC REGRESSION VARIABLES Calc.bin
  /METHOD=ENTER Inouye80vs20 
  /CONTRAST (Inouye80vs20)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Collagen.bin
  /METHOD=ENTER Inouye80vs20 
  /CONTRAST (Inouye80vs20)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Fat.bin_40
  /METHOD=ENTER Inouye80vs20 
  /CONTRAST (Inouye80vs20)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Fat.bin_10
  /METHOD=ENTER Inouye80vs20 
  /CONTRAST (Inouye80vs20)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES IPH.bin
  /METHOD=ENTER Inouye80vs20 
  /CONTRAST (Inouye80vs20)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Macrophages.bin 
  /METHOD=ENTER Inouye80vs20 
  /CONTRAST (Inouye80vs20)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES SMC.bin
  /METHOD=ENTER Inouye80vs20 
  /CONTRAST (Inouye80vs20)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


*univariate, but corrected for age and gender.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_macmean 
  /METHOD=ENTER Inouye80vs20
  /METHOD=ENTER Age Gender.


REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_smcmean 
  /METHOD=ENTER Inouye80vs20
  /METHOD=ENTER Age Gender.


REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_VD
  /METHOD=ENTER Inouye80vs20
  /METHOD=ENTER Age Gender.


LOGISTIC REGRESSION VARIABLES Calc.bin
  /METHOD=ENTER Inouye80vs20 Age Gender
  /CONTRAST (Inouye80vs20)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Collagen.bin
    /METHOD=ENTER Inouye80vs20 Age Gender
  /CONTRAST (Inouye80vs20)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Fat.bin_40
   /METHOD=ENTER Inouye80vs20 Age Gender 
  /CONTRAST (Inouye80vs20)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Fat.bin_10
    /METHOD=ENTER Inouye80vs20 Age Gender
  /CONTRAST (Inouye80vs20)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES IPH.bin
   /METHOD=ENTER Inouye80vs20 Age Gender
  /CONTRAST (Inouye80vs20)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES Macrophages.bin 
    /METHOD=ENTER Inouye80vs20 Age Gender
  /CONTRAST (Inouye80vs20)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

LOGISTIC REGRESSION VARIABLES SMC.bin
    /METHOD=ENTER Inouye80vs20 Age Gender 
  /CONTRAST (Inouye80vs20)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).



*******************************Multivariate only age, gender,  OR year, symptoms.4g + baseline differences p<0.2*************************.

*only pragmatic.
REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_macmean 
  /METHOD=ENTER Inouye80vs20
  /METHOD=ENTER Age Gender ORyear Symptoms.4g.

*+ p<0.2 baseline characteristics.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_macmean 
  /METHOD=ENTER Inouye80vs20
  /METHOD=ENTER Age Gender ORyear Symptoms.4g DM.composite Stenosis_con_bin.

*only pragmatic.
REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_smcmean 
  /METHOD=ENTER Inouye80vs20
  /METHOD=ENTER Age Gender ORyear Symptoms.4g.

* p<0.2 added.
REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_smcmean 
  /METHOD=ENTER Inouye80vs20
  /METHOD=ENTER Age Gender ORyear Symptoms.4g stenosis_con_bin.


*VD.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LN_VD
  /METHOD=ENTER Inouye80vs20
  /METHOD=ENTER Age Gender ORyear Symptoms.4g BMI.


*pragmatic, no additonal p<0.2.
LOGISTIC REGRESSION VARIABLES Calc.bin
  /METHOD=ENTER Inouye80vs20 Age Gender ORyear Symptoms.4g
  /CONTRAST (Inouye80vs20)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


*pragmatic.
LOGISTIC REGRESSION VARIABLES Collagen.bin
    /METHOD=ENTER Inouye80vs20 Age Gender ORyear Symptoms.4g
  /CONTRAST (Inouye80vs20)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*with covariates p<0.2.

LOGISTIC REGRESSION VARIABLES Collagen.bin
    /METHOD=ENTER Inouye80vs20 Age Gender ORyear Symptoms.4g stenosis_con_bin
  /CONTRAST (Inouye80vs20)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
            /CONTRAST (stenosis_con_bin)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*pragmatic.
LOGISTIC REGRESSION VARIABLES Fat.bin_40
   /METHOD=ENTER Inouye80vs20 Age Gender ORyear Symptoms.4g
  /CONTRAST (Inouye80vs20)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*with p<0.2 added (let op max 7 covariates).

LOGISTIC REGRESSION VARIABLES Fat.bin_40
   /METHOD=ENTER Inouye80vs20 Age Gender ORyear Symptoms.4g BMI DM.Composite Stenosis_con_bin
  /CONTRAST (Inouye80vs20)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
     /CONTRAST (DM.composite)=Indicator(1)
          /CONTRAST (stenosis_con_bin)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*pragmatic.
LOGISTIC REGRESSION VARIABLES Fat.bin_10
    /METHOD=ENTER Inouye80vs20 Age Gender ORyear Symptoms.4g
  /CONTRAST (Inouye80vs20)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*p<0.2 added.
LOGISTIC REGRESSION VARIABLES Fat.bin_10
    /METHOD=ENTER Inouye80vs20 Age Gender ORyear Symptoms.4g DM.Composite Stenosis_con_bin 
  /CONTRAST (Inouye80vs20)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
     /CONTRAST (DM.composite)=Indicator(1)
          /CONTRAST (stenosis_con_bin)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*pragmatic.
LOGISTIC REGRESSION VARIABLES IPH.bin
   /METHOD=ENTER Inouye80vs20 Age Gender ORyear Symptoms.4g
  /CONTRAST (Inouye80vs20)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*with p<0.2 BMI.

LOGISTIC REGRESSION VARIABLES IPH.bin
   /METHOD=ENTER Inouye80vs20 Age Gender ORyear Symptoms.4g BMI
  /CONTRAST (Inouye80vs20)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*empirically, no additional p<0.2 confounders.
LOGISTIC REGRESSION VARIABLES Macrophages.bin 
    /METHOD=ENTER Inouye80vs20 Age Gender ORyear Symptoms.4g
  /CONTRAST (Inouye80vs20)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*empircally.
LOGISTIC REGRESSION VARIABLES SMC.bin
    /METHOD=ENTER Inouye80vs20 Age Gender ORyear Symptoms.4g
  /CONTRAST (Inouye80vs20)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

*added p<0.2 with risk614.
LOGISTIC REGRESSION VARIABLES SMC.bin
    /METHOD=ENTER Inouye80vs20 Age Gender ORyear Symptoms.4g stenosis_con_bin 
  /CONTRAST (Inouye80vs20)=Indicator(1)
   /CONTRAST (Gender)=Indicator(1)
    /CONTRAST (Symptoms.4g)=Indicator(1)
     /CONTRAST (stenosis_con_bin)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).




*********************************************************************************************************************************************************************.

****************************************Inouye vs plaque vulnerability score****************************************************************************************.

* plaque score gebruiken om multiple testing theorie tegen te gaan. 

*** syntax- Plaque vulnerability**.
COMPUTE Macro_instab = -999.
IF macrophages.bin=2 Macro_instab=1.
IF macrophages.bin=1 Macro_instab=0.
EXECUTE.

COMPUTE Fat10_instab = -999.
IF Fat.bin_10=2 Fat10_instab=1.
IF Fat.bin_10=1 Fat10_instab=0.
EXECUTE.

COMPUTE coll_instab=-999.
IF Collagen.bin=2 coll_instab=0.
IF Collagen.bin=1 coll_instab=1.
EXECUTE.


COMPUTE SMC_instab=-999.
IF SMC.bin=2 SMC_instab=0.
IF SMC.bin=1 SMC_instab=1.
EXECUTE.

COMPUTE IPH_instab=-999.
IF IPH.bin=0 IPH_instab=0.
IF IPH.bin=1 IPH_instab=1.
EXECUTE.

COMPUTE Instability=Macro_instab + Fat10_instab +  coll_instab + SMC_instab + IPH_instab.
EXECUTE.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT Instability
  /METHOD=ENTER MetaGRSstandard Age Gender PC1 PC2 PC3 PC4
  /METHOD=ENTER Symptoms.4g ORyear.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT Instability
  /METHOD=ENTER Inouye20th
  /METHOD=ENTER Age Gender
    /METHOD=ENTER Symptoms.4g ORyear.


*idem maar dan log regressie. 


LOGISTIC REGRESSION VARIABLES Inouye20th 
   /METHOD=ENTER  Instability
  /METHOD=ENTER Age Gender OKyear Symptoms.4g
  /CONTRAST (Gender)=Indicator(1)
   /CONTRAST (Symptoms.4g)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


*plaque score met top20 vs 80.


LOGISTIC REGRESSION VARIABLES Inouye80vs20
   /METHOD=ENTER  Instability
  /METHOD=ENTER Age Gender OKyear Symptoms.4g
  /CONTRAST (Gender)=Indicator(1)
   /CONTRAST (Symptoms.4g)=Indicator(1)
  /PRINT=CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).


REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT Instability
  /METHOD=ENTER MetaGRSstandard Age Gender PC1 PC2 PC3 PC4
  /METHOD=ENTER Symptoms.4g ORyear.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT Instability
  /METHOD=ENTER Inouye20th
  /METHOD=ENTER Age Gender
    /METHOD=ENTER Symptoms.4g ORyear.


*********************************************************************************************************************************************************************.
*******************************************************KM curves & LIFETABLES**************************************************************************************************************.


SURVIVAL TABLE=ep_com_t_3years BY FH_CVD_comp2(0 1)
  /INTERVAL=THRU 3.2 BY 0.2
  /STATUS=epcom.3years(1)
  /PRINT=TABLE
  /PLOTS (SURVIVAL)=ep_com_t_3years BY  FH_CVD_comp2.

KM ep_com_t_3years BY FH_CVD_comp2
  /STATUS=epcom.3years(1)
  /PRINT TABLE MEAN
  /PLOT SURVIVAL
  /TEST LOGRANK
  /COMPARE OVERALL POOLED.



USE ALL.
COMPUTE filter_$=( ~ (MISSING(Inouye80vs20))).
VARIABLE LABELS filter_$ ' ~ (MISSING(Inouye80vs20)) (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

KM ep_com_t_3years BY  Inouye80vs20
  /STATUS=epcom.3years(1)
  /PRINT  TABLE MEAN
  /PLOT SURVIVAL
  /TEST LOGRANK
  /COMPARE OVERALL POOLED.



SURVIVAL TABLE=ep_com_t_3years BY Inouye80vs20(0 1)
  /INTERVAL=THRU 3.2 BY 0.2
  /STATUS=epcom.3years(1)
  /PRINT=TABLE
  /PLOTS (SURVIVAL)=ep_com_t_3years BY  Inouye80vs20.



KM ep_com_t_3years BY  Inouye20th
  /STATUS=epcom.3years(1)
  /PRINT  TABLE MEAN
  /PLOT SURVIVAL
  /TEST LOGRANK
  /COMPARE OVERALL POOLED.

SURVIVAL TABLE=ep_com_t_3years BY Inouye20th(0 1)
  /INTERVAL=THRU 3.2 BY 0.2
  /STATUS=epcom.3years(1)
  /PRINT=TABLE
  /PLOTS (SURVIVAL)=ep_com_t_3years BY  Inouye20th.





*********************************************************************************************************************************************************************.
* extra: nu Asymp meegenomen, hoeveel van deze asymp hebben wel een CVD event gehad?.
*********************************************************************************************************************************************************************.


*select: gegenotypeerd cohort en asymp (symptoms.4g =0).

USE ALL.
COMPUTE filter_$=(((~ (SYSMIS(CAD_PRSInouye)))  AND (symptoms.4g=0))).
VARIABLE LABELS filter_$ '((~ (SYSMIS(CAD_PRSInouye)))  AND (symptoms.4g=0)) (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.


*voor totalen (gehele genotypeerd cohort).
FREQUENCIES VARIABLES=Gender SmokerCurrent DM.composite risk614 Hypertension1 CAD_history Stroke_history PAOD 
    Symptoms.4g Stenosis_ipsilateral stenosis_con_bin Med.statin.derived Med.all.antiplatelet YearORnew
  /ORDER=ANALYSIS.


COMPUTE asymp_comorb=0.
IF ((CAD_history=1) OR (Stroke_history=1) OR (PAOD=1) OR (Peripheral.interv=1) OR (StrokeTIA_history=1)) asymp_comorb=1.
EXECUTE.


COMPUTE asymp_comorb_DM=0.
IF ((CAD_history=1) OR (Stroke_history=1) OR (PAOD=1) OR (Peripheral.interv=1) OR (StrokeTIA_history=1) OR (DM.composite=1)) asymp_comorb_DM=1.
EXECUTE.

COMPUTE asymp_comorb_risk=0.
IF ((CAD_history=1) OR (Stroke_history=1) OR (PAOD=1) OR (Peripheral.interv=1) OR (StrokeTIA_history=1) OR (DM.composite=1) OR (hypertension1=1) OR (risk614=1) OR (smokercurrent=1)) asymp_comorb_risk=1.
EXECUTE.

COMPUTE asymp_comorb_risk1=0.
IF ((CAD_history=1) OR (Stroke_history=1) OR (PAOD=1) OR (Peripheral.interv=1) OR (StrokeTIA_history=1) OR (DM.composite=1) OR (hypertension1=1) OR (risk614=1)) asymp_comorb_risk1=1.
EXECUTE.
