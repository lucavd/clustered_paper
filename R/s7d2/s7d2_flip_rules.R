# s7d2_flip_rules.R -- the sign-flip test at b0 = 0 (null of the paper) as saved in S7d / S7e (p_flip) and with the
# D-57 rules (p0_guard: not computable when the reduced fit diverges or the score is numerically zero), on the same
# 4,752,000 datasets: per point and arm, share answering (denominator: estimable designs), rejection rate, MCSE and
# the flags "clearly above 7.5%" / "clearly below 2.5%" (>= 20 answers). Written after the run (open doubt of S7d2).
# Outputs: tables/V4_S7d2_flip_rules_dp.csv, tables/V4_S7d2_flip_rules.csv
source(path.expand("~/MMM/R/s7d2/s7d2_common.R")); setDTthreads(8L)
TA <- s7d2_tasks()
X <- rbindlist(parallel::mclapply(TA$file, function(f) readRDS(f)[, .(dgp, scen, N, ntime, label, sim_id, design_degenerate, old_p_flip, p0_guard)], mc.cores = 16L))
X[, nd := !design_degenerate]
L <- rbind(X[, .(dgp, scen, N, ntime, label, version = "as in S7d", nd, use = nd & is.finite(old_p_flip), rej = nd & is.finite(old_p_flip) & old_p_flip < CI_ALPHA)],
           X[, .(dgp, scen, N, ntime, label, version = "with D-57 rules", nd, use = nd & is.finite(p0_guard), rej = nd & is.finite(p0_guard) & p0_guard < CI_ALPHA)])
D <- L[, .(nd = sum(nd), n = sum(use), pu = 100 * sum(use) / sum(nd), rate = 100 * sum(rej) / sum(use)), by = .(dgp, scen, N, ntime, label, version)]
D[, mcse := sqrt(rate * (100 - rate) / n)]; D[, `:=`(hi = n >= 20 & rate - 2 * mcse > 7.5, lo = n >= 20 & rate + 2 * mcse < 2.5)]
D[, small := N <= 25 & ntime <= 3]
fwrite(D, file.path(TABDIR, "V4_S7d2_flip_rules_dp.csv"))
S <- rbind(D[, .(set = "all 2,376 points", points = .N, mean_pu = mean(pu), mean_rate = mean(rate, na.rm = TRUE), n_hi = sum(hi), n_lo = sum(lo), points_lt20 = sum(n < 20)), by = .(label, version)],
           D[small == TRUE, .(set = "small designs (T 2-3, N <= 25)", points = .N, mean_pu = mean(pu), mean_rate = mean(rate, na.rm = TRUE), n_hi = sum(hi), n_lo = sum(lo), points_lt20 = sum(n < 20)), by = .(label, version)],
           D[dgp == DGP_COP, .(set = "copula (480 points)", points = .N, mean_pu = mean(pu), mean_rate = mean(rate, na.rm = TRUE), n_hi = sum(hi), n_lo = sum(lo), points_lt20 = sum(n < 20)), by = .(label, version)],
           D[dgp != DGP_COP, .(set = "study mechanisms (1,896 points)", points = .N, mean_pu = mean(pu), mean_rate = mean(rate, na.rm = TRUE), n_hi = sum(hi), n_lo = sum(lo), points_lt20 = sum(n < 20)), by = .(label, version)])
S[, label := fifelse(label == "manual", "H0", "H1")]; setorder(S, set, label, version)
print(S, digits = 4); fwrite(S, file.path(TABDIR, "V4_S7d2_flip_rules.csv"))
print(D[label == "manual" & version == "with D-57 rules" & small == TRUE, .(pu = mean(pu), rate = mean(rate, na.rm = TRUE)), by = .(N, ntime)][order(ntime, N)], digits = 3)
