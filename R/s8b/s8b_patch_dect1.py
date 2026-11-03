
import re,io
def sub1(path, old, new):
    s=open(path).read(); assert s.count(old)==1, (path, old, s.count(old)); open(path,'w').write(s.replace(old,new))
p='R/s5/s5_build.R'
sub1(p, 'PV <- V[, .(dgp, scen, N, ntime, arm, method, sim_id, p, u0)]',
 '# DEC-T1 (S8b, D-73): size-adjusted rule of Section 2.5 = 95th percentile (type 7) of the usable null |estimate/se|.\n# sa_one is unchanged: it receives -|z|; type 7 is reflection-symmetric, so its 5th percentile is -q95(|z0|) and\n# p1 < q0 <=> |z1| > q95. The likelihood-ratio test has no z and stays on p.\nPV <- V[, .(dgp, scen, N, ntime, arm, method, sim_id, p = fifelse(method == "glmer_lrt", p, -abs(estimate / se)), u0)]')
sub1(p, '}, mc.cores = 40, mc.preschedule = FALSE)', '}, mc.cores = as.integer(Sys.getenv("S5_CORES", "40")), mc.preschedule = FALSE)')
for p,old in [('R/s5/s5_run.R','mc.cores = 45, mc.preschedule = FALSE)'),('R/s6/s6_run.R','mc.cores = 45, mc.preschedule = FALSE)'),('R/s6/s6_run_boby.R','mc.cores = 45, mc.preschedule = FALSE)')]:
    sub1(p, old, 'mc.cores = as.integer(Sys.getenv("S5_CORES", "45")), mc.preschedule = FALSE)')
p='R/s5/s5_run.R'
sub1(p, 'binom.test(k, n, 0.05)$p.value, 0), mc.cores = 40)', 'binom.test(k, n, 0.05)$p.value, 0), mc.cores = as.integer(Sys.getenv("S5_CORES", "40")))')
p='R/s7f/s7f_table5_recompute.R'
s=open(p).read()
old='Q0 <- V[label == "manual" & u0 %in% TRUE & is.finite(p), .(q0 = if (.N >= 2) quantile(p, 0.05, type = 7, names = FALSE) else NA_real_), by = k]'
assert s.count(old)==1
s=s.replace(old, '# DEC-T1 (S8b, D-73): Wald tests use the declared rule, 95th percentile (type 7) of the usable null |estimate/se|,\n# written as the 5th percentile of -|z| (reflection-symmetric); the likelihood-ratio test stays on p.\nV[, pv := fifelse(method == "glmer_lrt", p, -abs(estimate / se))]\nQ0 <- V[label == "manual" & u0 %in% TRUE & is.finite(pv), .(q0 = if (.N >= 2) quantile(pv, 0.05, type = 7, names = FALSE) else NA_real_), by = k]')
old2='H1 <- V[label == "mmlb", .(n = sum(u0 %in% TRUE), kk = sum(u0 %in% TRUE & is.finite(p) & p < 0.05), p = list(p[u0 %in% TRUE & is.finite(p)])), by = k]'
assert s.count(old2)==1
s=s.replace(old2,'H1 <- V[label == "mmlb", .(n = sum(u0 %in% TRUE), kk = sum(u0 %in% TRUE & is.finite(p) & p < 0.05), p = list(pv[u0 %in% TRUE & is.finite(pv)])), by = k]')
s=s.replace('# 5th percentile, type 7, of the usable null p-values of the same cell)','# 95th percentile, type 7, of the usable null |estimate/se| of the same cell, D-73)')
open(p,'w').write(s)
p='R/postrun/agg_v4.R'
s=open(p).read(); old='          q <- quantile(g0$p[s0], 0.05, type = 7, names = FALSE)'
assert s.count(old)==1
s=s.replace(old,'          # Superseded for the paper by D-73 (S8b): the size-adjusted power printed in the manuscript comes from\n          # R/s5/s5_build.R with the 95th percentile of |estimate/se|; crit_p_sa/power_sa here keep the old p rule.\n'+old)
open(p,'w').write(s)
print("ok")
