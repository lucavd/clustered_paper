# 03_guidance_by_region.py
# v9 rewrite (2026-11-08/10): region guidance table (Table 4 of v9).
# Cells 46 and 50 of the rewrite session, extracted verbatim (see provenance/); paths made relative. Run from the repository root.
# Input: tables/v9/region_summary_main7.csv (02_results_by_mechanism.py). Output: tables/v9/guidance_by_region.csv
import os, numpy as np, pandas as pd
D = 'tables/v9'
# cell 46: region x strategy aggregates (main-text strategies)
rs = pd.read_csv(f'{D}/region_summary_main7.csv')
main7 = ['gee_exch_unc', 'gee_exch_md', 'gee_exch_fg', 'gee_ar1_fg', 'glmer', 'mmm_cc', 'flip']
lab = dict(zip(main7, ['GEE unc', 'GEE MD', 'GEE FG', 'GEE-AR1 FG', 'GLMM', 'MMM', 'Flip']))
r = rs[rs.strategy.isin(main7)].copy()
r['reg'] = r.Tband + ' ' + r.Nband
g = r.groupby(['strategy', 'reg']).agg(above=('n_above', 'sum'), maxmeanrate=('mean_rate', 'max'), minmeanrate=('mean_rate', 'min'),
                                       minpu=('min_pu', 'min'), meanpu=('mean_pu', 'mean'), pow=('mean_upow', 'mean')).reset_index()
order = [f'{t} {n}' for t in ['T=2', 'T=3-4', 'T=6-10'] for n in ['N=10', 'N=15-25', 'N=35-65', 'N=75-145']]
# cell 50: strategies with no design point clearly above 7.5% under any mechanism, by region
iv = ['gee_exch_unc', 'gee_exch_md', 'gee_exch_fg', 'gee_ar1_fg', 'glmer', 'mmm_cc']
rows = []
for reg in order:
    x = g[g.reg == reg].set_index('strategy')
    valid = [s for s in iv if x.loc[s, 'above'] == 0]
    rows.append(dict(region=reg, valid=', '.join(lab[s] for s in valid) or 'none',
        best_valid=max(valid, key=lambda s: x.loc[s, 'meanpu']) if valid else None,
        **{f'pu_{lab[s]}': x.loc[s, 'meanpu'] for s in iv + ['flip']},
        **{f'pow_{lab[s]}': x.loc[s, 'pow'] for s in iv + ['flip']},
        **{f'above_{lab[s]}': int(x.loc[s, 'above']) for s in iv},
        **{f'minrate_{lab[s]}': x.loc[s, 'minmeanrate'] for s in iv},
        **{f'maxrate_{lab[s]}': x.loc[s, 'maxmeanrate'] for s in iv}))
gt = pd.DataFrame(rows)
gt.to_csv(f'{D}/guidance_by_region.csv', index=False); print(gt[['region', 'valid']].to_string())
