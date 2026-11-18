# 01_build_dp_tidy.py
# v9 rewrite (2026-11-08/10): per-design-point tidy table .
# Cells extracted verbatim from the rewrite session (cells 34, 36, 37, 46, 50; see provenance/), with file paths
# made relative to the repository root. Run from the repository root:  python py/v9/01_build_dp_tidy.py
# Inputs:  tables/v9/S11_degen_table5_points_versions.csv.gz  (copy of results_v4/S11/degen/table5_points_versions.csv.gz)
#          tables/v9/S10_tabelle_s4_all_dp.csv                (copy of results_v4/S10/tabelle/s4_all_dp.csv; consistency check only)
# Output:  tables/v9/v8_dp_tidy.csv.gz
import os, numpy as np, pandas as pd
D = 'tables/v9'
v = pd.read_csv(f'{D}/S11_degen_table5_points_versions.csv.gz')
s4 = pd.read_csv(f'{D}/S10_tabelle_s4_all_dp.csv')
# cell 34: base version; consistency check against the S10 normal/t-reference counts
b = v[v.version == 'base'].drop(columns='version').copy()
s4h0 = s4[s4.label == 'manual'][['dgp', 'scen', 'N', 'ntime', 'strategy', 'n', 'kn', 'kt']].rename(columns={'n': 'n_s4'})
m = b.merge(s4h0, on=['dgp', 'scen', 'N', 'ntime', 'strategy'], how='left')
chk = m.dropna(subset=['kt'])
print('n match', (chk.n == chk.n_s4).mean().round(3), 'k match', (chk.k == chk.kt).mean().round(3))
# cell 36-37: region bands and tidy table
d = b.copy()
d['mcse'] = np.sqrt(d.rate * (100 - d.rate) / d.n.where(d.n > 0))
d['Nband'] = pd.cut(d.N, [0, 10, 25, 65, 200], labels=['N=10', 'N=15-25', 'N=35-65', 'N=75-145'])
d['Tband'] = pd.cut(d.ntime, [0, 2, 4, 10], labels=['T=2', 'T=3-4', 'T=6-10'])
d['set'] = d.grp.map({'main': 'main', 'copula': 'main', 'm67': 'M6-M7'})
d.loc[d.scen == 'M1c', 'set'] = 'M1c'
d.to_csv(f'{D}/v8_dp_tidy.csv.gz', index=False); print('tidy', d.shape)
