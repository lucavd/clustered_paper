# s6_mmm_tail.R -- is the gap between MMM empirical SE and its SEs due to a few extreme estimates?
# Stored results only (results_v4/S4/rows_v4.rds), all replicates, H0, T = 10, N = 145, time:trt.
# Compares, on the same replicates, the MMM and uncorrected-GEE estimates: SD, robust SD (MAD),
# quantiles, and the MMM - GEE difference. Output: tables/V4_S6_mmm_tail.csv
source(path.expand("~/MMM/R/s6/s6_common.R"))
V <- readRDS(file.path(OUT_S4, "rows_v4.rds"))
V <- V[method %in% c("mmm", "gee") & N == 145 & ntime == 10 & label == "manual" & term == "time:trt" & is.finite(estimate) & is.finite(se)]
V[, model := as.integer(sub("^M([1-5]).*$", "\\1", scen))]
W <- dcast(V, dgp + model + sim_id ~ method, value.var = c("estimate", "se"))
W <- W[is.finite(estimate_mmm) & is.finite(estimate_gee)]
T <- W[, .(n = .N, sd_mmm = sd(estimate_mmm), mad_mmm = mad(estimate_mmm), sd_gee = sd(estimate_gee), mad_gee = mad(estimate_gee),
           med_se_mmm = median(se_mmm), sd_diff = sd(estimate_mmm - estimate_gee), cor = cor(estimate_mmm, estimate_gee),
           q001 = quantile(estimate_mmm, 0.001), q999 = quantile(estimate_mmm, 0.999), n_out4 = sum(abs(estimate_mmm - median(estimate_mmm)) > 4 * mad(estimate_mmm)),
           sd_mmm_trim = sd(estimate_mmm[abs(estimate_mmm - median(estimate_mmm)) <= 4 * mad(estimate_mmm)])), by = .(dgp, model)]
fwrite(T, file.path(TABDIR, "V4_S6_mmm_tail.csv")); print(T, digits = 3)
