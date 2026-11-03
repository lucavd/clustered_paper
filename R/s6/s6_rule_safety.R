# =============================================================================
# s6_rule_safety.R -- S6: evidence for the decision rule, judged design point by design point.
# A strategy is a method with a reference distribution (normal "n" or t with N - p degrees of
# freedom "t"; the LRT has only its chi-square). At each design point (main grid and M1c) it is
#   clearly anticonservative  if rate - 2 MCSE > 7.5   (Bradley's upper limit for a 5% test)
#   clearly conservative      if rate + 2 MCSE < 2.5
#   not assessable            if fewer than 20 usable null replicates
# with the binomial MCSE of the design point. Usable share and power per attempt (rejections on
# the H1 arm over non-degenerate replicates, failures counted as no answer) by the same groups.
# Input : tables/V4_S6_dp_rates.csv.  Outputs: tables/V4_S6_rule_safety.csv (by strategy x N x
#         series band, and by strategy x N over all T), tables/V4_S6_rule_safety_dp.csv (flags)
# =============================================================================
source(path.expand("~/MMM/R/s6/s6_common.R"))
X <- fread(file.path(TABDIR, "V4_S6_dp_rates.csv"))
L <- rbind(X[, .(dgp, scen, N, ntime, arm, method, n, nd, kk = k, ref = "n")],
           X[method != "glmer_lrt", .(dgp, scen, N, ntime, arm, method, n, nd, kk = kt, ref = "t")])
L[, strategy := paste(method, ref, sep = "|")]
H0 <- L[arm == "H0"][, rate := fifelse(n > 0, 100 * kk / n, NA_real_)]
H0[, mc := fifelse(n > 0, 100 * sqrt((rate / 100) * (1 - rate / 100) / n), NA_real_)]
H0[, `:=`(assess = n >= 20, pu = 100 * n / nd)]
H0[, `:=`(hi = assess & rate - 2 * mc > 7.5, lo = assess & rate + 2 * mc < 2.5)]
H1 <- L[arm == "H1", .(dgp, scen, N, ntime, strategy, upow = 100 * kk / nd)]
D <- merge(H0, H1, by = c("dgp", "scen", "N", "ntime", "strategy"))
D[, band := fifelse(ntime <= 3, "T2-3", fifelse(ntime <= 6, "T4-6", "T8-10"))]
fwrite(D[, .(dgp, scen, N, ntime, band, strategy, n, nd, rate, mc, assess, hi, lo, pu, upow)], file.path(TABDIR, "V4_S6_rule_safety_dp.csv"))
sm <- function(by) D[, .(n_dp = .N, n_hi = sum(hi), n_lo = sum(lo), n_na = sum(!assess), max_rate = max(rate[assess]),
                         n_point_above = sum(assess & rate > 7.5), pu = mean(pu), upow = mean(upow)), by = by]
S <- rbind(sm(c("strategy", "N", "band")), sm(c("strategy", "N"))[, band := "all"], sm(c("strategy", "band"))[, N := NA_integer_],
           sm("strategy")[, `:=`(N = NA_integer_, band = "all")], use.names = TRUE)
setorder(S, strategy, band, N); fwrite(S, file.path(TABDIR, "V4_S6_rule_safety.csv"))
print(S[is.na(N) & band == "all"][order(n_hi)], digits = 3)
