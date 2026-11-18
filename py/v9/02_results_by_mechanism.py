# 02_results_by_mechanism.py
# v9 rewrite (2026-11-08): type I error, usable share and power by mechanism and region; best and robust strategy per region; sensitivity sets.
# Cells extracted verbatim, in execution order, from the session that produced them (see py/v9/provenance/);
# only file paths were changed (artifact store -> tables/v9/). Run from the repository root.
# Input: tables/v9/v8_dp_tidy.csv.gz (py/v9/01_build_dp_tidy.py). Cell 13 is omitted: it failed and was replaced by cell 14.

# ---- cell 2: Loading design-point table and inspecting
import pandas as pd, numpy as np, os
p='tables/v9/v8_dp_tidy.csv.gz'
d=pd.read_csv(p)
print(d.shape); print(d.dtypes.to_dict())
print(d.groupby(['set','dgp']).size()/19)
print(d.Nband.unique(), d.Tband.unique(), d.strategy.unique())
print(d.head(3).T)

# ---- cell 3: Computing region and mechanism summaries
os.makedirs('tables/v9',exist_ok=True)
print(d.groupby(['dgp','Nband','Tband']).size().loc['GenBinaryY']/19 if False else d[d.set=='main'].groupby(['dgp','Nband','Tband']).size().unstack()/19)
print('rate NaN:',d.rate.isna().sum(),' n<20:',(d.n<20).sum(), ' hi&n<20:',(d.hi&(d.n<20)).sum())
print(d[['n','nd']].describe().T)
# check hi/lo consistent with Clopper-Pearson
from scipy.stats import beta
a=d[d.n>=20]
lo_ci=np.where(a.k>0,beta.ppf(.025,a.k,a.n-a.k+1),0); up_ci=np.where(a.k<a.n,beta.ppf(.975,a.k+1,a.n-a.k),1)
print('hi check mismatches',((lo_ci>0.075)!=a.hi).sum(),' lo check',((up_ci<0.025)!=a.lo).sum())

# ---- cell 4: Building summary CSVs for main set
labels={'gee_exch_unc':'GEE, uncorrected','gee_exch_md':'GEE, Mancl-DeRouen','gee_exch_fg':'GEE, Fay-Graubard',
 'gee_ar1_fg':'GEE AR(1), Fay-Graubard','glmer':'GLMM, Wald','mmm_cc':'MMM (with check)','flip':'Sign-flip test'}
main7=list(labels)
mechs=['GenBinaryY','GLMM-native','GLMM-AR1','Copula-AR1']
Nb=['N=10','N=15-25','N=35-65','N=75-145']; Tb=['T=2','T=3-4','T=6-10']
all19=main7+[s for s in sorted(d.strategy.unique()) if s not in main7]
def summ(g):
    a=g[g.n>=20]; na=len(a)
    return pd.Series({'n_points':len(g),'n_assessed':na,
      'mean_rate':a.rate.mean() if na else np.nan,
      'mean_abs_dep':(a.rate-5).abs().mean() if na else np.nan,
      'mean_mcse':a.mcse.mean() if na else np.nan,
      'max_rate':a.rate.max() if na else np.nan,
      'n_above':int(a.hi.sum()),'pct_above':100*a.hi.mean() if na else np.nan,
      'n_below':int(a.lo.sum()),'pct_below':100*a.lo.mean() if na else np.nan,
      'mean_pu':g.pu.mean(),'min_pu':g.pu.min(),'mean_upow':g.upow.mean()})
M=d[d.set=='main'].copy()
M['dgp']=pd.Categorical(M.dgp,mechs,ordered=True); M['Nband']=pd.Categorical(M.Nband,Nb,ordered=True)
M['Tband']=pd.Categorical(M.Tband,Tb,ordered=True); M['strategy']=pd.Categorical(M.strategy,all19,ordered=True)
reg=M.groupby(['dgp','Tband','Nband','strategy'],observed=True).apply(summ).reset_index()
reg.insert(4,'label',reg.strategy.map(labels).fillna(reg.strategy.astype(str)))
mech=M.groupby(['dgp','strategy'],observed=True).apply(summ).reset_index()
mech.insert(2,'label',mech.strategy.map(labels).fillna(mech.strategy.astype(str)))
r7=reg[reg.strategy.isin(main7)]; m7=mech[mech.strategy.isin(main7)]
for df,f in [(r7,'region_summary_main7.csv'),(m7,'mechanism_summary_main7.csv'),(mech,'mechanism_summary_all19.csv'),(reg,'region_summary_all19.csv')]:
    df.round(3).to_csv('tables/v9/'+f,index=False)
print(len(r7),len(m7),len(mech))
pd.set_option('display.width',250)
print(m7[['dgp','strategy','n_assessed','mean_rate','n_above','pct_above','n_below','mean_pu','mean_upow']].round(1).to_string())

# ---- cell 5: Identifying best calibrated strategy per cell
interval_all=[s for s in all19 if s not in ('flip','glmer_lrt')]
interval_main=[s for s in main7 if s!='flip']
reg['calibrated']=(reg.n_assessed>0)&(reg.n_above==0)&reg.mean_rate.between(2.5,7.5)
def lab(s): return labels.get(s,s)
rows=[]
for cset,cands in [('main6',interval_main),('all17',interval_all)]:
    for (m,t,n),g in reg[reg.strategy.isin(cands)].groupby(['dgp','Tband','Nband'],observed=True):
        c=g[g.calibrated].sort_values('mean_pu',ascending=False)
        if len(c)==0:
            rows.append(dict(candidate_set=cset,dgp=m,Tband=t,Nband=n,n_calibrated=0,best='none',best_label='none',best_pu=np.nan,best_rate=np.nan,best_upow=np.nan,equivalent='',equivalent_labels=''))
            continue
        b=c.iloc[0]; eq=c[(c.mean_pu>=b.mean_pu-2.5)&(c.strategy!=b.strategy)]
        rows.append(dict(candidate_set=cset,dgp=m,Tband=t,Nband=n,n_calibrated=len(c),best=b.strategy,best_label=lab(b.strategy),
            best_pu=b.mean_pu,best_rate=b.mean_rate,best_upow=b.mean_upow,
            equivalent=';'.join(eq.strategy.astype(str)),equivalent_labels='; '.join(eq.strategy.astype(str).map(lab))))
best=pd.DataFrame(rows)
best.round(2).to_csv('tables/v9/best_by_mechanism_region.csv',index=False)
for cset in ['main6','all17']:
    b=best[best.candidate_set==cset]
    w=b.pivot_table(index=['Tband','Nband'],columns='dgp',values='best',aggfunc='first',observed=True)[mechs]
    print(cset); print(w.to_string())


# ---- cell 6: Inspecting best-strategy detail and equivalents
b=best[best.candidate_set=='main6']
print(b[['Tband','Nband','dgp','n_calibrated','best','best_pu','best_rate','equivalent']].round(1).to_string())
r=reg[reg.strategy.isin(interval_main)]
print(r[(r.Nband.isin(['N=10','N=15-25']))][['dgp','Tband','Nband','strategy','n_assessed','mean_rate','n_above','mean_pu','calibrated']].round(1).to_string())

# ---- cell 8: Quantifying mechanism dependence of best strategy
best['n_assessed_best']=[reg.set_index(['dgp','Tband','Nband','strategy']).n_assessed.get((r.dgp,r.Tband,r.Nband,r.best),np.nan) if r.best!='none' else np.nan for r in best.itertuples()]
best['low_usable_flag']=best.best_pu<50
best.round(2).to_csv('tables/v9/best_by_mechanism_region.csv',index=False)
hl=[]
for cset in ['main6','all17']:
    for (t,n),g in best[best.candidate_set==cset].groupby(['Tband','Nband'],sort=False):
        g=g.set_index('dgp').loc[mechs]
        sets=[set([r.best]+[e for e in r.equivalent.split(';') if e]) for r in g.itertuples()]
        common=set.intersection(*sets)
        hl.append(dict(candidate_set=cset,Tband=t,Nband=n,best_GenBinaryY=g.best.iloc[0],best_GLMM_native=g.best.iloc[1],best_GLMM_AR1=g.best.iloc[2],best_Copula_AR1=g.best.iloc[3],
          n_distinct_best=g.best.nunique(),common_best_or_equiv=';'.join(sorted(common)) or 'none',min_best_pu=g.best_pu.min()))
hl=pd.DataFrame(hl)
hl.to_csv('tables/v9/best_agreement_by_region.csv',index=False)
print(hl.round(1).to_string())
for cset in ['main6','all17']:
    h=hl[hl.candidate_set==cset]
    print(cset,'regions with >1 distinct best:',(h.n_distinct_best>1).sum(),'/12; with no common best-or-equivalent:',(h.common_best_or_equiv=='none').sum())

# ---- cell 9: Deriving robust choice across mechanisms
rows=[]
R=reg.copy(); R['strategy']=R.strategy.astype(str)
for cset,cands in [('main6',interval_main),('all17',interval_all)]:
    for t in Tb:
        for n in Nb:
            g=R[(R.Tband==t)&(R.Nband==n)&R.strategy.isin(cands)]
            agg=g.groupby('strategy').agg(n_mech_cal=('calibrated','sum'),min_pu=('mean_pu','min'),max_pu=('mean_pu','max'),
                    min_pow=('mean_upow','min'),max_pow=('mean_upow','max'),max_rate=('mean_rate','max'),min_rate=('mean_rate','min'),tot_above=('n_above','sum'))
            rob=agg[agg.n_mech_cal==4].sort_values('min_pu',ascending=False)
            fl=R[(R.Tband==t)&(R.Nband==n)&(R.strategy=='flip')]
            base=dict(candidate_set=cset,Tband=t,Nband=n,n_robust=len(rob),
                      flip_tot_above=int(fl.n_above.sum()),flip_min_pu=fl.mean_pu.min(),flip_pow_range=f"{fl.mean_upow.min():.1f}-{fl.mean_upow.max():.1f}")
            if len(rob)==0:
                # least-bad: most mechanisms calibrated, then min_pu
                lb=agg.sort_values(['n_mech_cal','min_pu'],ascending=False).iloc[0]
                rows.append({**base,'robust_choice':'none','label':'none','min_pu':np.nan,'max_pu':np.nan,'min_upow':np.nan,'max_upow':np.nan,
                             'rate_range':'','other_robust':'','closest':f"{lb.name} (calibrated under {int(lb.n_mech_cal)}/4; min pu {lb.min_pu:.1f})"})
                continue
            b=rob.iloc[0]
            rows.append({**base,'robust_choice':b.name,'label':lab(b.name),'min_pu':b.min_pu,'max_pu':b.max_pu,'min_upow':b.min_pow,'max_upow':b.max_pow,
                'rate_range':f"{b.min_rate:.1f}-{b.max_rate:.1f}",
                'other_robust':'; '.join(f"{s} ({r.min_pu:.1f})" for s,r in rob.iloc[1:].iterrows()),'closest':''})
rc=pd.DataFrame(rows); rc.round(2).to_csv('tables/v9/robust_choice.csv',index=False)
print(rc.drop(columns=['other_robust']).round(1).to_string())
print(rc[['candidate_set','Tband','Nband','other_robust']].to_string())

# ---- cell 10: Adding tie-break and robustness counts
R['absdep']=R.mean_abs_dep
rows=[];cnt=[]
for cset,cands in [('main6',interval_main),('all17',interval_all)]:
    for t in Tb:
        for n in Nb:
            g=R[(R.Tband==t)&(R.Nband==n)&R.strategy.isin(cands)]
            agg=g.groupby('strategy').agg(n_mech_cal=('calibrated','sum'),min_pu=('mean_pu','min'),max_pu=('mean_pu','max'),
                    min_pow=('mean_upow','min'),max_pow=('mean_upow','max'),max_rate=('mean_rate','max'),min_rate=('mean_rate','min'),worst_absdep=('absdep','max'))
            for s,r in agg.iterrows(): cnt.append(dict(candidate_set=cset,Tband=t,Nband=n,strategy=s,robust=r.n_mech_cal==4))
            rob=agg[agg.n_mech_cal==4].copy(); rob['pu_r']=rob.min_pu.round(0)
            rob=rob.sort_values(['pu_r','worst_absdep'],ascending=[False,True])
            fl=R[(R.Tband==t)&(R.Nband==n)&(R.strategy=='flip')]
            base=dict(candidate_set=cset,Tband=t,Nband=n,n_robust=len(rob),
                      flip_points_above=int(fl.n_above.sum()),flip_min_pu=fl.mean_pu.min(),flip_max_pu=fl.mean_pu.max(),flip_min_upow=fl.mean_upow.min(),flip_max_upow=fl.mean_upow.max())
            if len(rob)==0:
                lb=agg.sort_values(['n_mech_cal','min_pu'],ascending=False).iloc[0]
                rows.append({**base,'robust_choice':'none','label':'none','min_pu':np.nan,'max_pu':np.nan,'min_upow':np.nan,'max_upow':np.nan,
                             'min_rate':np.nan,'max_rate':np.nan,'other_robust':'','closest_nonrobust':f"{lb.name} (calibrated under {int(lb.n_mech_cal)}/4 mechanisms; min pu {lb.min_pu:.1f})"}); continue
            b=rob.iloc[0]
            rows.append({**base,'robust_choice':b.name,'label':lab(b.name),'min_pu':b.min_pu,'max_pu':b.max_pu,'min_upow':b.min_pow,'max_upow':b.max_pow,
                'min_rate':b.min_rate,'max_rate':b.max_rate,'other_robust':'; '.join(f"{s} (min pu {r.min_pu:.1f})" for s,r in rob.iloc[1:].iterrows()),'closest_nonrobust':''})
rc=pd.DataFrame(rows); rc.round(2).to_csv('tables/v9/robust_choice.csv',index=False)
cnt=pd.DataFrame(cnt)
print(rc[['candidate_set','Tband','Nband','n_robust','robust_choice','min_pu','max_pu','min_upow','max_upow','min_rate','max_rate']].round(1).to_string())
print(cnt[cnt.candidate_set=='main6'].groupby('strategy').robust.sum())
print(cnt[cnt.candidate_set=='all17'].groupby('strategy').robust.sum().sort_values())
c6=cnt[cnt.candidate_set=='all17']; print(c6[(c6.strategy=='gee_ar1_fg')&~c6.robust])

# ---- cell 11: Summarising M1c and M6-M7 sensitivity sets
S=d[d.set.isin(['M1c','M6-M7'])].copy()
S['strategy']=pd.Categorical(S.strategy,all19,ordered=True)
sens=S.groupby(['set','dgp','strategy'],observed=True).apply(summ,include_groups=False).reset_index()
sens.insert(3,'label',sens.strategy.astype(str).map(lab))
sens.round(3).to_csv('tables/v9/sensitivity_M1c_M67.csv',index=False)
print(sens[sens.strategy.isin(main7)][['set','dgp','strategy','n_points','n_assessed','mean_rate','max_rate','n_above','n_below','mean_pu','mean_upow']].round(1).to_string())
# region-level check of robust choices in sensitivity sets
Sr=S.groupby(['set','dgp','Tband','Nband','strategy'],observed=True).apply(summ,include_groups=False).reset_index()
Sr['calibrated']=(Sr.n_assessed>0)&(Sr.n_above==0)&Sr.mean_rate.between(2.5,7.5)
chk=rc.merge(Sr.assign(strategy=Sr.strategy.astype(str)),left_on=['Tband','Nband','robust_choice'],right_on=['Tband','Nband','strategy'])
print(chk[['candidate_set','set','dgp','Tband','Nband','robust_choice','n_points','n_assessed','mean_rate','n_above','mean_pu','calibrated']].round(1).to_string())

# ---- cell 14: Saving sensitivity check of robust choices
print([c for c in chk.columns if c.endswith(('_x','_y'))])
chk2=chk.rename(columns={'max_rate_y':'max_rate'})
chk2[['candidate_set','set','dgp','Tband','Nband','robust_choice','n_points','n_assessed','mean_rate','max_rate','n_above','n_below','mean_pu','mean_upow','calibrated']].round(2).to_csv('tables/v9/sensitivity_robust_choice_check.csv',index=False)
print(chk2.groupby(['candidate_set','set']).agg(cells=('calibrated','size'),calibrated=('calibrated','sum'),above=('n_above','sum')))
print(sens[sens.strategy.isin(['gee_ar1_fg','gee_exch_fg','gee_ind_fg','gee_ind_cr2','glmer','gee_exch_md'])][['set','dgp','strategy','n_assessed','n_above','pct_above']].round(1).to_string())

