# 06_method_differences_table.py
# v9 rewrite (2026-11-08): Table 1 of v9 (key differences between the strategies), entered from the v8 manuscript and supplement.
# Cells extracted verbatim, in execution order, from the session that produced them (see py/v9/provenance/);
# only file paths were changed. Run from the repository root.
# Outputs: tables/v9/method_differences.csv and .md
import os, pandas as pd, numpy as np
os.chdir('tables/v9')

# ---- cell 37: Building method-differences table
cols=['Method family','Target estimand','Within-subject dependence','Small-sample device','Reference distribution','Returns','Typical failure with few clusters / short series','Source (v8)']
rows=[
['GEE, uncorrected sandwich',
 'Marginal (population-averaged) log odds ratio',
 'Working correlation (exchangeable, independence or AR(1)) as nuisance; robust sandwich variance',
 'None',
 't with N - p df (normal in Table 3)',
 'Estimate + Wald CI',
 'Anticonservative: clearly above 7.5% at 242 main/M1c, 69 M6-M7, 100 Copula-AR1 design points',
 'Ms 2.4; Table 5; S14'],
['GEE + Mancl-DeRouen (MD) or Fay-Graubard (FG); also Kauermann-Carroll (KC), CR2',
 'Marginal (population-averaged) log odds ratio',
 'As uncorrected GEE; level of FG changed <0.5 points across working correlations',
 'Leverage-adjusted sandwich: MD (I-H)^-1, KC/CR2 (I-H)^-1/2, FG bounded inflation (b = 0.75)',
 't with N - p df (CR2 not with Satterthwaite df)',
 'Estimate + Wald CI',
 'Separated data: plausible estimate, no warning (FG-AR(1) at N = 10, T = 2); KC, MD stay anticonservative; SE missing in sparse subsamples',
 'Ms 2.4, 3.4, 3.5, 3.11, 4.1, 4.3; S4'],
['GLMM, Wald test, Laplace approximation (lme4)',
 'Conditional (subject-specific); larger than marginal by about (1 + k^2 sigma^2)^1/2',
 'Normal random intercept, random slopes in M3-M5; conditional independence given random effects',
 'None in main run; adaptive quadrature (10 nodes) examined post hoc, random intercept only',
 't with N - p df; LRT against chi-square(1)',
 'Estimate + Wald CI on conditional scale',
 'Laplace anticonservative (515 points); random slope at T = 2 not identifiable (usable 0.002%); misses serial dependence (Copula-AR1)',
 'Ms 1, 2.4, 3.3, 3.4, 3.6, 3.10, 4.1, 4.4; S12'],
['Marginalized multilevel model (MMLB::mm) with convergence check',
 'Marginal (population-averaged) log odds ratio',
 'Likelihood: latent variable (SD may trend over time) plus first-order transition term',
 'Check: boundary fit (nlm code not 1/2, or latent-SD slope < 0.001) replaced by converged fit without trend',
 't with N - p df',
 'Estimate + Wald CI (model-based SE)',
 'Boundary non-convergence: no result in 16% of M3-M5 fits (32% at T <= 3); conservative at N = 10',
 'Ms 2.4, 3.2, 3.4, 4.2; S16, S20'],
['Sign-flip score test',
 'No estimate; tests zero marginal interaction via effective score of logistic model without interaction',
 'Subjects treated as independent units; dependence absorbed in per-subject score sums',
 'Randomization: all sign flips up to 12 subjects, else 1,999 random plus identity',
 'Sign-flip (randomization) distribution',
 'p-value only; inverted interval often unbounded at T <= 3',
 'Not computable if null fit diverges or scores vanish (e.g. one arm all zero); power 4.2% at T <= 3, N <= 25',
 'Ms 2.4, 3.7, 4.1; S7'],
]
md=pd.DataFrame(rows,columns=cols)
wc=md.iloc[:,1:7].map(lambda s: len(s.replace('/',' / ').split()))
print(wc.max().to_dict(), wc.values.max())

# ---- cell 38: Tightening cells and saving method table
F='Typical failure with few clusters / short series'
md.loc[0,F]='Anticonservative: clearly above 7.5% at 242 main/M1c, 69 M6-M7 and 100 Copula-AR1 design points'
md.loc[1,F]='Separated data: plausible estimate, no warning (FG-AR(1) fails at N = 10, T = 2); KC, MD exceed 7.5% at 138-207, 20-63 points'
md.loc[2,F]='Laplace anticonservative (above 7.5% at 515 points); random slope unidentifiable at T = 2 (usable 0.002%); ignores serial dependence'
md.loc[3,F]='Boundary non-convergence: no result in 16% of M3-M5 fits (32% at T <= 3); conservative at N = 10; rarely answers when separated'
md.loc[4,F]='Not computable if null fit diverges or scores vanish (one arm all zero); power 4.2% at T <= 3, N <= 25'
wc=md.iloc[:,1:7].map(lambda s: len(s.split()))
print(wc.values.max(), wc.max().to_dict())
md.to_csv('method_differences.csv',index=False)
def esc(s): return s.replace('|','\\|')
lines=['# Key differences between the method families (v8 manuscript and supplement)','',
 '| '+' | '.join(cols)+' |','|'+'---|'*len(cols)]
for _,r in md.iterrows(): lines.append('| '+' | '.join(esc(x) for x in r)+' |')
lines+=['','Notes. N = number of clusters (subjects); p = number of mean-model coefficients; T = maximum number of occasions; df = degrees of freedom; k = 16√3/(15π), σ² = random-intercept variance. "Clearly above 7.5%": lower limit of the exact 95% interval of the null rejection rate above 7.5%, counted over the 2,376 design points with a t reference (Table 5); for GEE corrections, the range spans the three working correlations. Kauermann-Carroll counts use *fewclust* (published formula); *geesmv* computes a different quantity (S4). The shrinkage GLM (rescaled GLMM) and the GLMM likelihood-ratio test (no interval; above 7.5% at 412 design points) are not separate rows. Adaptive quadrature lowered the GLMM null rejection rate by 3.1-4.0 points but left excess at T = 10, mainly on Copula-AR1 (§3.10). "Rarely answers when separated": 4 usable MMM results on datasets with an infinite interaction estimate in the whole study (§3.4).']
open('method_differences.md','w').write('\n'.join(lines)+'\n')
print('\n'.join(lines[:4]))

# ---- cell 39: Shortening over-length failure cells
print(wc[F].to_dict())
md.loc[1,F]='Separated data: plausible estimate, no warning (FG-AR(1) failure, N = 10, T = 2); KC, MD exceed 7.5% often'
md.loc[2,F]='Laplace anticonservative (515 points above 7.5%); random slope unidentifiable at T = 2; ignores serial dependence'
md.loc[3,F]='Boundary non-convergence: no result in 16% of M3-M5 fits (32% at T <= 3); conservative at N = 10'
md.loc[4,F]='Not computable if null fit diverges or scores vanish (one arm all zero); power 4.2% at T <= 3, N <= 25'
wc=md.iloc[:,1:7].map(lambda s: len(s.split())); print(wc.values.max(), wc[F].to_dict())
md.to_csv('method_differences.csv',index=False)
lines=lines[:4]
for _,r in md.iterrows(): lines.append('| '+' | '.join(esc(x) for x in r)+' |')
note=('Notes. N = number of clusters (subjects); p = number of mean-model coefficients; T = maximum number of occasions; df = degrees of freedom; k = 16√3/(15π), σ² = random-intercept variance. '
 '"Clearly above 7.5%": lower limit of the exact 95% interval of the null rejection rate above 7.5%, over the 2,376 design points with a t reference (Table 5). '
 'With a t reference KC was clearly above 7.5% at 138-207 and MD at 20-63 design points across working correlations; FG with exchangeable correlation at 15 (§3.4, §4.3). KC counts use *fewclust* (published formula); *geesmv* computes a different quantity (S4). '
 'GLMM with random slopes at T = 2: usable Wald test in 0.002% of attempts (Table 4). Adaptive quadrature lowered the GLMM null rejection rate by 3.1-4.0 points but left excess at T = 10, mainly on Copula-AR1 (§3.10). '
 'MMM returned 4 usable results on separated datasets (infinite interaction estimate) in the whole study (§3.4). '
 'The shrinkage GLM (rescaled GLMM) and the GLMM likelihood-ratio test (no interval; above 7.5% at 412 design points) are not separate rows.')
open('method_differences.md','w').write('\n'.join(lines)+'\n\n'+note+'\n')

# ---- cell 40: Final trim and save of method table
md.loc[4,F]='Not computable when null fit diverges or scores vanish (one arm all zero); power 4.2% at T <= 3, N <= 25'
md.loc[4,F]='Not computable when null fit diverges or scores vanish; power 4.2% at T <= 3, N <= 25'
wc=md.iloc[:,1:7].map(lambda s: len(s.split())); print(wc.values.max())
md.to_csv('method_differences.csv',index=False)
lines=lines[:4]
for _,r in md.iterrows(): lines.append('| '+' | '.join(esc(x) for x in r)+' |')
note2=note+' Sign-flip test not computable e.g. when one arm has only zero outcomes (S7).'
open('method_differences.md','w').write('\n'.join(lines)+'\n\n'+note2+'\n')

