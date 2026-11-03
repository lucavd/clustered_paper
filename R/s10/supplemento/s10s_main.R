
# s10s_main.R -- S10 "supplemento": independent per-cell aggregates for the main grid and M1c (Tables S3-S13,
# S16-S25, S36 and related text). Definitions from the rendered text (main 2.2-2.7, supplement); no project
# aggregation function is used. Inputs: compact copies made by s10s_build.R plus S6cc/cc_all.rds and truth files.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(3)
RP <- "/home/user/MMM/results_v4"; D <- file.path(RP, "S10/supplemento"); TB <- "/home/user/MMM/tables"
K <- c("dgp","scen","N","ntime","label","sim_id"); CELL <- c("dgp","scen","N","ntime","label")
SC <- c(paste0("M", 1:5), "M1c"); diag <- list()
t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%5.1f min] %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), s))
# ---- replicate universe: every generated dataset, its design and whether the outcome is constant
R <- readRDS(file.path(D, "recs_s7d.rds"))[scen %in% SC]
R[, gen := generated %in% TRUE][, deg := design_degenerate %in% TRUE]
diag$recs_n <- nrow(R); diag$not_gen <- R[gen == FALSE, .N]; diag$not_gen_not_deg <- R[!gen & !deg, .N]
diag$deg_per_dgp_label <- R[deg | !gen, .N, by = .(dgp, label, grid = scen == "M1c")]
R[, estd := gen & !deg]; R[, const := estd & (sum_y == 0 | sum_y == n_obs)]
diag$const_n <- R[const == TRUE, .N]
KEYS <- R[estd == TRUE, c(K, "const"), with = FALSE]
# ---- study rows (variant_cc): mmm (check), gee, gee_md/gee_fg (geesmv), glmer, shrink, glmer_lrt, glmer_bobyqa
S <- readRDS(file.path(D, "study_cc.rds"))[method != "generation"]; tm("study rows loaded")
S[method == "gee_kc", method := "gee_kc_smv"]
diag$dup_rows <- S[, .(cnt = .N), by = c(K, "method")][cnt > 1L, .N]
# Kauermann-Carroll from fewclust (S7d, exchangeable): estimate and SE from the stored per-replicate GEE rows
G <- readRDS(file.path(D, "gee_s7d.rds"))[cor == "exch" & scen %in% SC]
gst <- S[method == "gee", c(K, "estimate", "status"), with = FALSE]
G <- merge(G, gst, by = K, all.x = TRUE)
diag$kc_est_vs_study_gee_maxdiff <- G[is.finite(est) & is.finite(estimate), max(abs(est - estimate))]
diag$kc_gee_ok_but_no_study_est <- G[gee_ok %in% TRUE & !is.finite(estimate), .N]
KC <- G[, .(dgp, scen, N, ntime, label, sim_id, method = "gee_kc", estimate = est, se = se_kc,
            p = 2 * pnorm(-abs(est / se_kc)), converged = gee_ok, status, U0 = NA, U1 = NA)]
rm(G, gst); gc()
# MMM without the check (study S4 rows)
MN <- readRDS(file.path(D, "study_mmm_nocc.rds"))[scen %in% SC][, method := "mmm_nocc"]
# MMM with the check rebuilt from the text definition (main 2.4) on the S6cc fits of M3-M5
CC <- readRDS(file.path(RP, "S6cc/cc_all.rds"))
CC[, bnd_text := !(A_code %in% c(1L, 2L)) | (is.finite(A_slope_min) & A_slope_min < 0.001)]
CC[, useB := bnd_text & (B_ok %in% TRUE)]
diag$cc_tab <- CC[, .N, by = .(A_ok, bnd_text, bnd, B_ok, used, useB)]
diag$cc_Bcode_when_A_kept <- CC[bnd_text & B_ok %in% TRUE & used == "A", .N, by = .(B_code, B_se_fin = is.finite(B_se), A_code)]
diag$cc_Bcode_when_B_used <- CC[used == "B", .N, by = .(B_code)]
CC[, `:=`(e = fifelse(useB, B_est, A_est), s = fifelse(useB, B_se, A_se), pp = fifelse(useB, B_p, A_p))]
diag$cc_vs_stored <- CC[, .(n = .N, est_diff = sum(!(e == cc_est | (is.na(e) & is.na(cc_est)))), se_diff = sum(!(s == cc_se | (is.na(s) & is.na(cc_se)))))]
mstudy <- S[method == "mmm" & scen %in% c("M3","M4","M5"), c(K, "estimate", "se", "converged", "status"), with = FALSE]
cmp <- merge(CC[, c(K, "e", "s", "A_ok", "useB", "B_status", "A_status"), with = FALSE], mstudy, by = K, all = TRUE)
diag$cc_vs_variant_cc_mmm <- cmp[, .(n = .N, only_cc = sum(is.na(converged)), only_study = sum(is.na(A_ok)),
   est_diff = sum(is.finite(e) & is.finite(estimate) & abs(e - estimate) > 1e-12),
   fin_mismatch = sum(is.finite(e) != is.finite(estimate), na.rm = TRUE))]
MC <- cmp[, .(dgp, scen, N, ntime, label, sim_id, method = "mmm_cc_text", estimate = e, se = s,
              p = 2 * pnorm(-abs(e / s)), converged = A_ok %in% TRUE | useB %in% TRUE,
              status = fifelse(useB %in% TRUE, B_status, A_status), U0 = NA, U1 = NA)]
MC <- rbind(MC, S[method == "mmm" & scen %in% c("M1","M2","M1c")][, method := "mmm_cc_text"], use.names = TRUE, fill = TRUE)
rm(CC, cmp, mstudy); gc()
X <- rbindlist(list(S, KC, MN, MC), use.names = TRUE, fill = TRUE); rm(S, KC, MN, MC); gc(); tm("rows assembled")
X <- merge(X, KEYS, by = K)              # keep estimable designs only
diag$rows_by_method <- X[, .N, by = method]
# ---- per-replicate definitions (main 2.4-2.5)
Z <- qnorm(0.975); WALDM <- setdiff(unique(X$method), "glmer_lrt")
X[, z := estimate / se]
X[, fin := converged %in% TRUE & is.finite(estimate) & is.finite(se) & is.finite(p) & abs(estimate) < 10 & se > 0 & se < 10]
X[method == "glmer_lrt", fin := converged %in% TRUE & is.finite(p)]
X[, use := fin & !const]; X[, use_u1 := use & is.na(status)]
X[, np := fifelse(scen %in% c("M2","M5"), 5L, 4L)]
X[, `:=`(rn = abs(z) > Z, rt = abs(z) > qt(0.975, N - np), rp = p < 0.05)]
X[method == "glmer_lrt", `:=`(rn = p < 0.05, rt = p < 0.05)]
diag$p_vs_z_disagree <- X[use & method != "glmer_lrt", .(n = .N, disagree = sum(rn != rp)), by = method]
diag$const_fin_by_method <- X[const == TRUE, .(n = .N, fin = sum(fin)), by = method]
# truths: marginal per cell (generator / GH table), conditional per cell (NAT) or per replicate (AR1)
TT <- fread(file.path(TB, "V4_S4_truth_by_cell.csv"))
TR <- readRDS(file.path(RP, "S4/truth_ar1_rep.rds")); TR[, `:=`(N = as.integer(N), ntime = as.integer(ntime))]
X <- merge(X, TT[, .(dgp, scen, ntime, label, tm = truth_marginal, tc = truth_conditional)], by = c("dgp","scen","ntime","label"), all.x = TRUE)
X <- merge(X, TR[, .(dgp, scen, N, ntime, label, sim_id, tcr = truth_cond)], by = K, all.x = TRUE)
X[dgp == "GLMM-AR1", tc := tcr]; X[, tcr := NULL]
X[, truth := fifelse(method %in% c("glmer","glmer_bobyqa"), tc, tm)]
X[dgp == "GenBinaryY" & method %in% c("glmer","glmer_bobyqa"), truth := NA_real_]
X[, `:=`(cvn = abs(estimate - truth) <= Z * se, cvt = abs(estimate - truth) <= qt(0.975, N - np) * se)]
tm("per-replicate flags done")
# ---- common set: replicates usable for all eight primary methods (Table S6 / S19)
PRIM <- c("mmm","gee","gee_md","gee_kc","gee_fg","glmer","shrink","glmer_lrt")
CM <- X[method %in% PRIM, .(nu = sum(use), nm = .N), by = K][nu == 8L, c(K), with = FALSE][, com := TRUE]
X <- merge(X, CM, by = K, all.x = TRUE); X[is.na(com), com := FALSE]
# ---- size-adjusted power: per design point and method, critical value from the usable null replicates
Q <- X[label == "manual" & use == TRUE, .(crit = if (method[1] == "glmer_lrt") quantile(p, 0.05, names = FALSE) else quantile(abs(z), 0.95, names = FALSE),
                                         n0 = .N), by = .(dgp, scen, N, ntime, method)]
X <- merge(X, Q, by = c("dgp","scen","N","ntime","method"), all.x = TRUE)
X[, `:=`(sa = fifelse(method == "glmer_lrt", p < crit, abs(z) > crit), sa_ge = fifelse(method == "glmer_lrt", p <= crit, abs(z) >= crit))]
# ---- cell aggregates
A <- X[, .(n_row = .N, n_fin = sum(fin), n_use = sum(use), k_n = sum(use & rn), k_t = sum(use & rt), k_p = sum(use & rp),
           n_u1 = sum(use_u1), k_n_u1 = sum(use_u1 & rn), k_t_u1 = sum(use_u1 & rt),
           n_com = sum(use & com), k_n_com = sum(use & com & rn),
           n_tr = sum(use & is.finite(truth)), cov_n = sum(use & cvn, na.rm = TRUE), cov_t = sum(use & cvt, na.rm = TRUE),
           s_est = sum(estimate[use]), s_est2 = sum(estimate[use]^2), s_se = sum(se[use]), s_se2 = sum(se[use]^2), s_err = sum((estimate - truth)[use & is.finite(truth)]),
           s_err2 = sum(((estimate - truth)^2)[use & is.finite(truth)]),
           k_sa = sum(use & sa, na.rm = TRUE), k_sa_ge = sum(use & sa_ge, na.rm = TRUE), k_sa37 = sum(use & sa & abs(z) > 37, na.rm = TRUE),
           k_n37 = sum(use & rn & abs(z) > 37), crit = crit[1], n0 = n0[1],
           n_use_inclconst = sum(fin), k_n_inclconst = sum(fin & rn), k_t_inclconst = sum(fin & rt)), by = c(CELL, "method")]
NE <- R[, .(n_all = .N, n_gen = sum(gen), n_est = sum(estd), n_const = sum(const, na.rm = TRUE)), by = CELL]
MM <- unique(A$method); NE2 <- NE[rep(seq_len(.N), each = length(MM))][, method := rep(MM, nrow(NE))]
A <- merge(NE2, A, by = c(CELL, "method"), all.x = TRUE)
for (v in setdiff(names(A), c(CELL, "method", "crit", "n0", "s_est", "s_est2", "s_se", "s_err", "s_err2"))) set(A, which(is.na(A[[v]])), v, 0)
fwrite(A, file.path(D, "s10s_main_cells.csv"))
# per design point sd of the estimates for SE ratio: need sum over usable of (est - mean)^2 -> from s_est, s_est2
sink(file.path(D, "s10s_main_diag.txt")); print(diag); sink()
tm("MAIN OK")
