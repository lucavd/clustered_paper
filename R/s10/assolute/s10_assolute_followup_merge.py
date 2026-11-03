#!/usr/bin/env python3
"""S10 'assolute' follow-up: replace the rows of tables/V7_S10_assolute_check.csv (written by s10_assolute_checks.py)
with the follow-up verdicts of tables/V7_S10_assolute_followup.csv (same columns plus 'followup').
Run from the repository root after s10_assolute_checks.py."""
import pandas as pd
ck = pd.read_csv("tables/V7_S10_assolute_check.csv", keep_default_na=False, dtype=str)
fu = pd.read_csv("tables/V7_S10_assolute_followup.csv", keep_default_na=False, dtype=str)
for _, r in fu.iterrows():
    i = ck.index[ck.cand_id == r.cand_id]; assert len(i) == 1
    for c in ck.columns: ck.at[i[0], c] = r[c]
ck.to_csv("tables/V7_S10_assolute_check.csv", index=False)
print(ck.verdict.value_counts().to_dict())
