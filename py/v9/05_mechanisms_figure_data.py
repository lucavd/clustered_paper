# 05_mechanisms_figure_data.py
# v9 rewrite (2026-11-08): data of Figure 1 (lag correlation and marginal prevalence by mechanism).
# Cells extracted verbatim, in execution order, from the session that produced them (see py/v9/provenance/);
# only file paths were changed. Run from the repository root.
# Inputs: tables/v9/mech_*.csv, written by R/v9/mechanisms_generate.R (run on the server, ~2 min on 4 cores; rerun on 2026-11-18 gave
# mech_pair_corr.csv identical byte for byte to the original). Output: tables/v9/mechanisms_lag_correlation.csv
# The plotting cells (30 and the first lines of 31) are in py/v9/provenance/mechanisms_and_clinical_cells.json; they need the figure-style helpers.

# ---- cell 29: Aggregating lag correlations and prevalence
import pandas as pd, numpy as np
d='tables/v9/'
pairs=pd.read_csv(d+'mech_pair_corr.csv'); prev=pd.read_csv(d+'mech_prev_counts.csv'); pool=pd.read_csv(d+'mech_pooled_corr.csv'); info=pd.read_csv(d+'mech_dataset_info.csv')
# per replicate: n-weighted mean of pair/arm-specific correlations at each lag
rep=pairs.groupby(['mechanism','rep','lag']).apply(lambda g: pd.Series({'r':np.average(g.r,weights=g.n),'n_pairs':g.n.sum()})).reset_index()
lagt=rep.groupby(['mechanism','lag']).agg(r_mean=('r','mean'),r_sd=('r','std'),n_rep=('r','size'),n_pairs_total=('n_pairs','sum')).reset_index()
lagt['mcse']=lagt.r_sd/np.sqrt(lagt.n_rep); lagt['lo']=lagt.r_mean-2*lagt.mcse; lagt['hi']=lagt.r_mean+2*lagt.mcse
pl=pool.groupby(['mechanism','lag']).r_pooled.mean().reset_index()
lagt=lagt.merge(pl,on=['mechanism','lag'])
pr=prev.groupby(['mechanism','rep','occasion'])[['n','y1']].sum().reset_index(); pr['p']=pr.y1/pr.n
prt=pr.groupby(['mechanism','occasion']).agg(prev_mean=('p','mean'),prev_sd=('p','std'),n_rep=('p','size'),n_obs_total=('n','sum')).reset_index()
prt['mcse']=prt.prev_sd/np.sqrt(prt.n_rep); prt['lo']=prt.prev_mean-2*prt.mcse; prt['hi']=prt.prev_mean+2*prt.mcse
# by arm too
pa=prev.groupby(['mechanism','occasion','arm'])[['n','y1']].sum().reset_index(); pa['p']=pa.y1/pa.n
print(lagt.pivot(index='lag',columns='mechanism',values='r_mean').round(3))
print(lagt.mcse.max().round(4), (lagt.r_mean-lagt.r_pooled).abs().max().round(3))
print(prt.pivot(index='occasion',columns='mechanism',values='prev_mean').round(3))
print(info.groupby('mechanism')[['n_subj','n_obs','mean_y']].mean())

# ---- cell 31 (data part): Retitling panel b and writing CSV
a=lagt.rename(columns={'lag':'x','r_mean':'value','lo':'lo_2mcse','hi':'hi_2mcse','n_pairs_total':'n_total'})[['mechanism','x','value','mcse','lo_2mcse','hi_2mcse','n_rep','n_total','r_pooled']]
a.insert(0,'panel','a_lag_correlation')
b=prt.rename(columns={'occasion':'x','prev_mean':'value','lo':'lo_2mcse','hi':'hi_2mcse','n_obs_total':'n_total'})[['mechanism','x','value','mcse','lo_2mcse','hi_2mcse','n_rep','n_total']]
b.insert(0,'panel','b_prevalence'); b['r_pooled']=np.nan
out=pd.concat([a,b]); out['scenario']='M1'; out['hypothesis']='H0'; out['N_per_dataset']=2000; out['T']=10
out.to_csv('tables/v9/mechanisms_lag_correlation.csv',index=False,float_format='%.5g')
print(out.shape, out.groupby('panel').size().to_dict())
