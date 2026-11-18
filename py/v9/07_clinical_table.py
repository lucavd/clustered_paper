# 07_clinical_table.py
# v9 rewrite (2026-11-08): Table 5 of v9 (clinical applications), entered from the v8 supplement and checked against it (cell 45).
# Cells extracted verbatim, in execution order, from the session that produced them (see py/v9/provenance/);
# only file paths were changed. Run from the repository root.
# Inputs: manuscript/v9/sources/MMM_v8_SMMR_manuscript.md and _supplement.md (copied as ms.md, sup.md in a temporary folder).
# Outputs: tables/v9/clinical_table.csv and .md
import os, shutil, tempfile, pandas as pd, numpy as np
ROOT = os.getcwd(); TMP = tempfile.mkdtemp()
shutil.copy('manuscript/v9/sources/MMM_v8_SMMR_manuscript.md', os.path.join(TMP, 'ms.md'))
shutil.copy('manuscript/v9/sources/MMM_v8_SMMR_supplement.md', os.path.join(TMP, 'sup.md'))
os.chdir(TMP)

# ---- cell 44: Building clinical applications table
NA='n.a.'
C='Neonatal lung-ultrasound cohort (207 patients; posterior vs anterior zones)'
R='geepack respiratory trial (111 patients; treatment effect, full mean model)'
rows=[
 # dataset, method, scale, est(SE), usable 10..50, ratio 10..50, ratio as reported, source
 [C,'Marginalized model (MMM)','Marginal','0.35 (0.10)', '77.2','85.4','94.6','96.6','100.0', NA,NA,NA,NA,NA,'0.56 to 0.62 (robust SE: 0.85 to 0.99)','Ms 3.11; S9; Table S15'],
 [C,'GEE exchangeable, uncorrected','Marginal','0.38 (0.17)', '70.0','80.8','91.0','95.4','100.0', NA,NA,NA,NA,NA,NA,'Ms 3.11; Table S15'],
 [C,'GEE exchangeable, Fay-Graubard','Marginal','n.a. (same fit as uncorrected; corrected SE not reported)', '34.2','61.2','80.0','88.4','99.8', NA,NA,NA,NA,NA,NA,'Table S15'],
 [C,'GEE exchangeable, Mancl-DeRouen','Marginal','n.a. (same fit as uncorrected; corrected SE not reported)', '15.0','37.2','62.2','80.0','98.8', NA,NA,NA,NA,NA,NA,'Table S15'],
 [C,'GLMM (random intercept)','Conditional','0.48 (0.12)', '78.2','85.4','94.6','96.6','100.0', NA,NA,NA,NA,NA,'0.51 to 0.58','Ms 3.11; Table S15'],
 [R,'Marginalized model (MMM)','Marginal','1.3 (0.34)', '50','84','94','99.2','100.0', '0.7','0.7','0.7','0.7','0.9','0.70 at 10, 0.85 at 50 (robust SE: 0.55 to 0.87)','Table S29; Table S30; S9'],
 [R,'GEE exchangeable, uncorrected','Marginal','1.3 (0.35)', '47','83','94','99.2','100.0', '0.5','0.6','0.7','0.7','0.9',NA,'Table S29; Table S30'],
 [R,'GEE exchangeable, Fay-Graubard','Marginal','1.3 (0.35)', '23','68','88','98','100.0', '2.2','1.3','1.0','0.9','0.9',NA,'Table S29; Table S30'],
 [R,'GEE exchangeable, Mancl-DeRouen','Marginal','1.3 (0.37)', '23','68','88','97','100.0', '2.2','1.4','1.2','1.0','1.0',NA,'Table S29; Table S30'],
 [R,'GLMM (random intercept)','Conditional','2.1 (0.54)', '53','88','97','100.0','100.0', '0.7','0.7','0.6','0.7','0.7',NA,'Table S29; Table S30'],
]
Ns=[10,15,20,25,50]
cols=['Dataset','Method','Scale','Full-data estimate (SE)']+[f'Usable % N={n}' for n in Ns]+[f'SE/SD ratio (FPC) N={n}' for n in Ns]+['SE/SD ratio as reported in text','Source (v8)']
ct=pd.DataFrame(rows,columns=cols); ct.to_csv('clinical_table.csv',index=False)
print(ct.shape)

# ---- cell 45: Verifying clinical numbers against supplement
import re
sup=open('sup.md').read().split('\n')
def parse(lo,hi,ncol):
    out={}; cur=None
    for l in sup[lo:hi]:
        t=l.split()
        if not t or t[0].startswith('-') or t[0]=='Method': continue
        nums=[x for x in t if re.fullmatch(r'-?\d+(\.\d+)?',x)]
        name=' '.join(x for x in t if not re.fullmatch(r'-?\d+(\.\d+)?',x))
        if name: cur=name
        out[(cur,nums[0])]=nums[1:]
    return out
s15=parse(428,491,4); s30=parse(900,963,6)
mapC={'Marginalized model (MMM)':'Marginalized model','GEE exchangeable, uncorrected':'GEE, uncorrected','GEE exchangeable, Fay-Graubard':'GEE, Fay-Graubard','GEE exchangeable, Mancl-DeRouen':'GEE, Mancl-DeRouen','GLMM (random intercept)':'GLMM'}
bad=[]
for _,r in ct.iterrows():
    src=s15 if r.Dataset==C else s30
    for n in Ns:
        v=src[(mapC[r.Method],str(n))]
        if r[f'Usable % N={n}']!=v[3]: bad.append((r.Method,n,'usable',r[f'Usable % N={n}'],v[3]))
        if r.Dataset==R and r[f'SE/SD ratio (FPC) N={n}']!=v[5]: bad.append((r.Method,n,'ratio',r[f'SE/SD ratio (FPC) N={n}'],v[5]))
print(bad, len(s15), len(s30))

# ---- cell 46: Writing clinical markdown and paragraph
para=("In both datasets the two marginal methods agree on the full data (cohort: 0.35 by the marginalized model, 0.38 by GEE; trial: 1.3 and 1.3). "
"The mixed model returns a larger coefficient on its conditional scale (0.48; 2.1). This is the estimand gap of Section 1, not a disagreement about the data: "
"with a patient variance of 3.5 the attenuation factor is about 1.5, the factor used to build scenario M1c. "
"The subsamples reproduce three features of the simulation. Availability falls as patients become fewer, for every method. "
"The uncorrected sandwich understates the spread of the estimates at 10 to 15 patients (trial ratio 0.5 to 0.6), the anticonservatism seen on the grid. "
"The Fay-Graubard and Mancl-DeRouen corrections overstate it at 10 patients (ratio 2.2), the direction of their conservatism in the simulation. "
"The clinical data also expose what the simulation underplays. With richer mean models (seven coefficients in the trial; disease, ventilation and gestational age in the cohort), "
"the corrections often returned an estimate without a standard error. They were usable in 15% to 34% of 10-patient draws, against about 50% to 78% for the other methods; "
"in the simulation, with four or five coefficients, they lost little availability (Table 4). "
"The model-based standard errors of the marginalized and mixed models fell short of the spread in both datasets. In the cohort part of this comes from between-patient heterogeneity of the zone contrast "
"(random-slope variance 10), which none of the four mechanisms contains. The robust standard error of the marginalized model was close to the spread there (0.85 to 0.99). "
"The applications cannot show type I error or coverage. The reference is the full-data estimate, not a known truth, and the generating mechanism is unknown. "
"Nor can they test the separation failure of GEE: the trial had 365 of 5,000 subsamples with an infinite contrast and the cohort none, but without a known null the false positives of Section 3.4 cannot be counted. "
"Each rests on one dataset with one outcome.")
print(len(para.split()))

# ---- cell 47: Trimming paragraph and saving clinical table
para=("In both datasets the two marginal methods agree on the full data (cohort: 0.35 by the marginalized model, 0.38 by GEE; trial: 1.3 by both), "
"while the mixed model returns a larger coefficient on its conditional scale (0.48; 2.1). This is the estimand gap of Section 1, not a disagreement about the data; "
"with a patient variance of 3.5 the attenuation factor is about 1.5, as built into scenario M1c. "
"The subsamples reproduce three simulation features: availability falls with fewer patients for every method; "
"the uncorrected sandwich understates the spread at 10 to 15 patients (trial ratio 0.5 to 0.6); "
"and the Fay-Graubard and Mancl-DeRouen corrections overstate it at 10 patients (ratio 2.2), the direction of their conservatism on the grid. "
"The applications also expose what the simulation underplays. With richer mean models, the corrections often returned an estimate without a standard error "
"and were usable in only 15% to 34% of 10-patient draws, against 47% to 78% for the other methods; with four or five coefficients in the simulation they lost little availability (Table 4). "
"Model-based standard errors of the marginalized and mixed models fell short of the spread in both datasets; in the cohort, part of this reflects between-patient heterogeneity of the zone contrast "
"(random-slope variance 10), which no simulated mechanism contains, and the robust standard error of the marginalized model was close (0.85 to 0.99). "
"The applications cannot show type I error or coverage: the reference is the full-data estimate, not a known truth, and the mechanism is unknown. "
"Nor can they count the separated-data false positives of GEE (Section 3.4), although 365 of 5,000 trial subsamples had an infinite contrast.")
print(len(para.split()))
def row_md(r): return '| '+' | '.join(str(x) for x in r)+' |'
short=ct.copy()
short['Dataset']=short.Dataset.map({C:'Neonatal cohort (207)',R:'Respiratory trial (111)'})
mdcols=['Dataset','Method','Scale','Full-data estimate (SE)']+[f'Usable % {n}' for n in Ns]+[f'SE/SD {n}' for n in Ns]+['SE/SD as reported','Source (v8)']
L=['# Clinical applications: full-data estimates, usable share and SE calibration in subsamples','',
'Subsamples: 500 draws of whole patients at each size N = 10, 15, 20, 25, 50 patients. Usable % = draws returning a finite estimate and SE with |estimate| and SE ≤ 10 (all draws in the denominator). SE/SD = mean standard error (model-based for MMM and GLMM, sandwich for GEE) divided by the empirical SD of the estimates over subsamples, multiplied by √(1 − N/P) (finite-population correction, P = patients in the full dataset).','',
row_md(mdcols),'|'+'---|'*len(mdcols)]
for _,r in short.iterrows(): L.append(row_md(r.values))
L+=['','Notes.',
'1. Cohort contrast: posterior vs anterior zones, adjusted for disease, invasive ventilation and gestational age; trial contrast: treatment effect adjusted for baseline status, visit, age, sex and centre (full mean model, Table S29; reduced-model estimates in Table S29).',
'2. GLMM estimates are on the conditional (subject-specific) scale and are not comparable in size with the marginal estimates.',
'3. Cohort, N = 10: 22% of draws contained no invasively ventilated patient, so no method could fit the mean model; these draws are in the denominator (S17). Among returned fits, usable: GLMM 100%, MMM 100%, uncorrected GEE 92%, FG 45%, MD 20%.',
'4. Cohort SE/SD ratios are reported in v8 only as ranges over subsample sizes (§3.11, S9); per-N values are n.a. Cohort SE ratios for the GEEs and full-data SEs of the corrected GEEs are not reported in v8 (n.a.). With a patient-timepoint random intercept the GLMM ratio was 0.50 to 0.57; with a random slope for posterior zones 0.66 to 0.73 at ≥25 patients (S17).',
'5. Most failures of the corrections were an estimate without a standard error (Tables S15, S30). In the cohort FG and MD were computed with geesmv and agreed with fewclust within 1% in 99% of draws (S9). Kauermann-Carroll (not shown): cohort usable 11.2/34.8/61.2/79.6/98.8%; trial 1.3 (0.36), usable 22/69/89/97/100.0%.',
'6. Monte Carlo SE: at most 2.2 points for usable shares; at most 0.15 for the trial SE ratio (Table S30).',
'7. Trial: 365 of 5,000 subsamples had an infinite treatment contrast; cohort 0 (S19).',
'','## Linking the clinical findings to the simulation','',para,'']
open('clinical_table.md','w').write('\n'.join(L))

for f in ['clinical_table.csv', 'clinical_table.md']: shutil.copy(f, os.path.join(ROOT, 'tables/v9', f))
