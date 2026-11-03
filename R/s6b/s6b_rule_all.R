# s6b_rule_all.R -- S6b: evidence for the Figure 1 rule over all simulated design points: the 1,512 of the study
# (tables/variant_cc/V4_S6_rule_safety_dp.csv, MMM with the convergence check) and the 384 of the new scenarios
# (tables/V4_S6b_ns_rule_safety_dp.csv). Same flags (Bradley limits with the MCSE of the design point).
# Also: GB3 with the corrected completeness check (an analysis returning a fit without the time:trt term is an outcome).
# Outputs: tables/V4_S6b_rule_all.csv, tables/V4_S6b_ns_gates_post.csv (row GB3b added)
source(path.expand("~/MMM/R/s4/s4_common.R"))
A <- fread(file.path(REPO, "tables/variant_cc/V4_S6_rule_safety_dp.csv"))[, src := "studio"]
B <- fread(file.path(TABDIR, "V4_S6b_ns_rule_safety_dp.csv"))[!startsWith(strategy, "mmm_nocc")][, src := "nuovi"]
D <- rbind(A, B, fill = TRUE)
D[, Ng := fifelse(N <= 25, "N<=25", "N>=35")]
sm <- function(X, by) X[, .(n_dp = .N, n_hi = sum(hi), n_lo = sum(lo), max_rate = max(rate[assess], na.rm = TRUE), pu = mean(pu), upow = mean(upow)), by = by]
S <- rbind(sm(D, c("strategy", "src"))[, sub := "all"], sm(D, "strategy")[, `:=`(src = "tutti", sub = "all")],
           sm(D[band == "T2-3"], c("strategy", "Ng"))[, `:=`(src = "tutti", sub = "T2-3")][, sub := paste(sub, Ng)][, Ng := NULL], fill = TRUE)
S <- S[strategy %in% c("mmm|t", "mmm|n", "gee_fg|t", "gee_fg|n", "gee_md|t", "glmer_lrt|n", "glmer|t", "gee_kc|t")]
setorder(S, strategy, src, sub); fwrite(S, file.path(TABDIR, "V4_S6b_rule_all.csv")); print(S, digits = 3)
G <- fread(file.path(TABDIR, "V4_S6b_ns_gates_post.csv"))[gate == "GB3"]
R3 <- fread(file.path(TABDIR, "V4_S6b_ns_diag_gb3.csv"))
ok <- nrow(R3) > 0 && R3[, all(has_timetrt == FALSE & !grepl("ERROR", terms))]
G <- rbind(G, data.table(gate = "GB3b", pass = ok && G$pass == FALSE,
  detail = sprintf("corrected check: the %d replicates without a row for some analysis were re-run (%d analyses): every one returns a fit without the time:trt term (collinear design), 0 errors; as in the study these count as not usable", uniqueN(R3[, .(dgp, scen, N, ntime, label, sim_id)]), nrow(R3))))
fwrite(G, file.path(TABDIR, "V4_S6b_ns_gates_post.csv")); print(G)
