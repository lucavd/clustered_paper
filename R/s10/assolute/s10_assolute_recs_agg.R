## S10 track 'assolute': aggregate the per-replicate facts written by s10_assolute_recs.R.
## Usage: Rscript --vanilla s10_assolute_recs_agg.R <dir>   (reads recs_all.rds, fgpt_rows.rds; writes CSVs in <dir>)
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(1L)
d <- commandArgs(TRUE)[1]
R <- readRDS(file.path(d, "recs_all.rds"))
R[, ok_design := generated %in% TRUE & !(design_degenerate %in% TRUE)]
R[, const := ok_design & (sum_y == 0 | sum_y == n_obs)]
R[, flip_seen := !is.na(p_flip) | !is.na(flip_why)]
sc <- intersect(c("n_gee", "n_glmer", "n_glmer_lrt", "n_mmm", "n_shrink"), names(R))
A <- R[, c(list(n_recs = .N, n_gen = sum(generated %in% TRUE), n_deg = sum(design_degenerate %in% TRUE),
          n_notgen = sum(!(generated %in% TRUE)), n_notgen_nondeg = sum(!(generated %in% TRUE) & !(design_degenerate %in% TRUE)),
          n_ok_design = sum(ok_design), n_const = sum(const), n_const_all0 = sum(const & sum_y == 0),
          gee3_on_ok = sum(ok_design & n_gee_cor == 3), gee_missing_on_ok = sum(ok_design & n_gee_cor < 3),
          gee_on_notok = sum(!ok_design & n_gee_cor > 0), flip_seen_on_ok = sum(ok_design & flip_seen),
          flip_on_notok = sum(!ok_design & flip_seen)),
       lapply(setNames(sc, paste0(sc, "_missing_on_ok")), function(v) sum(ok_design & (is.na(get(v)) | get(v) == 0))),
       lapply(setNames(sc, paste0(sc, "_on_notok")), function(v) sum(!ok_design & !is.na(get(v)) & get(v) > 0))),
       by = .(run, dgp, scen, N, ntime, label)]
fwrite(A, file.path(d, "recs_by_dp.csv"))
K <- readRDS("/home/user/MMM/results_v4/S9/const_keys.rds")
C <- R[const == TRUE, .(dgp, scen, N, ntime, label, sim_id)]
both <- merge(C, K, by = c("dgp", "scen", "N", "ntime", "label", "sim_id"))
cmp <- data.table(recomputed = nrow(C), keys = nrow(K), both = nrow(both), only_recomputed = nrow(C) - nrow(both), only_keys = nrow(K) - nrow(both),
                  small_recomputed = C[ntime <= 3 & N <= 25, .N], all_datasets = nrow(R),
                  const_any_design = R[generated %in% TRUE & (sum_y == 0 | sum_y == n_obs), .N])
fwrite(cmp, file.path(d, "const_compare.csv"))
G <- readRDS(file.path(d, "fgpt_rows.rds"))
G <- merge(G, R[, .(dgp, scen, N, ntime, label, sim_id, const, ok_design)], by = c("dgp", "scen", "N", "ntime", "label", "sim_id"), all.x = TRUE)
G[, df := N - 5L]
G[, p_fg := 2 * pt(-abs(est / se_fg), df)]
G[, usable := gee_ok %in% TRUE & is.finite(est) & is.finite(se_fg) & is.finite(p_fg) & abs(est) < 10 & se_fg > 0 & se_fg < 10]
G[, rej := usable & p_fg < 0.05]
F <- G[, .(n_pre = sum(usable), k_pre = sum(rej), n_post = sum(usable & !const), k_post = sum(rej & !const),
           const_usable = sum(usable & const), const_usable_rej = sum(rej & const)), by = .(dgp, scen, N, ntime, label)]
fwrite(F, file.path(d, "fgpt_check.csv"))
print(cmp); print(F)
