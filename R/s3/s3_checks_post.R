# =============================================================================
# s3_checks_post.R -- S3 post-run expectations Q4-Q9 on results_v4/S3/.
# Reads every S3 chunk once (and, for blocks F and L, the matching results_full
# chunk), then evaluates the pre-registered criteria of docs/sessions/2026-10-16_S3.md.
# =============================================================================
source(path.expand("~/MMM/R/s3/s3_common.R"))
library(data.table)
TAB <- file.path(REPO, "tables")
tasks <- as.data.table(s3_tasks())
WALD <- c("mmm", "gee", "gee_fg", "gee_md", "gee_kc", "glmer", "shrink"); LRT <- "glmer_lrt"
deg_s2 <- fread(file.path(TAB, "V4_S2_degenerate_designs.csv"))[deg == TRUE, .(model, N, ntime, sim_id)]
cmp_cols <- c("estimate", "se", "p", "status")
per_chunk <- function(k) {
  t <- tasks[k]; if (!file.exists(t$file)) return(list(chunk = data.table(t[, .(key, block)], exists = FALSE)))
  x <- readRDS(t$file); s <- as.data.table(attr(x, "sims")); x <- as.data.table(x)
  ids <- (t$chunk - 1L) * 100L + 1:100
  ch <- data.table(t[, .(key, block, dgp, scen, model, N, ntime, label)], exists = TRUE,
    sims_ok = nrow(s) == 100L && identical(sort(as.integer(s$sim_id)), ids),
    every_sim_has_rows = all(ids %in% x$sim_id),
    gen_rows_only_on_failures = setequal(x[method == "generation", sim_id], s[generated == FALSE, sim_id]),
    n_gen_fail = sum(!s$generated), n_gen_fail_nondeg = sum(!s$generated & !s$design_degenerate),
    n_trace_analysis = x[method != "generation" & grepl("^\\[v4 trace\\]", status), .N],
    design_ref_checked = sum(!is.na(s$design_identical_ref)), design_ref_identical = sum(s$design_identical_ref %in% TRUE),
    re_ok = all(x[method == LRT, re_structure] == c("(1|id)", "(1|id)", "(1 + time|id)", "(1 + time|id)", "(1 + time + dose|id)")[t$model]),
    n_lrt_rows = x[method == LRT, .N])
  # Q5: gee identity in block L (all terms)
  w <- NULL; q9 <- NULL
  if (t$block == "L" || t$block %in% c("F1", "F2")) o <- as.data.table(readRDS(t$orig_file))
  if (t$block == "L") {
    a <- x[method == "gee" & !grepl("^\\[v4 trace\\]", status %||% "")]; b <- o[method == "gee"]
    setkey(a, sim_id, term); setkey(b, sim_id, term)
    ch[, `:=`(gee_sims_equal = setequal(unique(a$sim_id), unique(b$sim_id)), gee_rows_v4 = nrow(a), gee_rows_orig = nrow(b),
              gee_identical = nrow(a) == nrow(b) && identical(a$term, b$term) && identical(as.integer(a$sim_id), as.integer(b$sim_id)) &&  # sim_id is integer in v4, numeric in results_full
                all(vapply(cmp_cols, function(v) identical(a[[v]], b[[v]]), TRUE)))]
  }
  # Q8 controprova: LRT vs Wald glmer, time:trt, both without warning
  wald <- if (t$block == "L") o[method == "glmer" & term == "time:trt"] else x[method == "glmer" & term == "time:trt"]
  lr <- x[method == LRT & term == "time:trt"]
  m <- merge(lr[, .(sim_id, e_lrt = estimate, s_lrt = status, c_lrt = converged)], wald[, .(sim_id, e_w = estimate, s_w = status, c_w = converged)], by = "sim_id")
  m <- m[c_lrt & c_w & is.na(s_lrt) & is.na(s_w)]
  ch[, `:=`(q8_pairs = nrow(m), q8_close = sum(abs(m$e_lrt - m$e_w) < 1e-3))]
  # Q9 inputs: F1 time:trt rows, paired with the original H0
  if (t$block == "F1") {
    u <- function(D) { D <- D[!is.na(term) & term == "time:trt" & method %in% c(WALD, LRT)]
      D[, U := fifelse(method %in% WALD, converged & is.finite(estimate) & is.finite(se) & abs(estimate) < 10 & se > 0 & se < 10 & is.finite(p), converged & is.finite(p))]
      D[is.na(U), U := FALSE]; D[, R := U & p < 0.05]; D[, .(method, sim_id = as.integer(sim_id), U, R)] }
    q9 <- merge(u(x), u(o), by = c("method", "sim_id"), all = TRUE, suffixes = c("_new", "_old"))
    q9[, `:=`(model = t$model, N = t$N, ntime = t$ntime)]
  }
  degk <- s[design_degenerate == TRUE, .(sim_id = as.integer(sim_id))][, `:=`(dgp = t$dgp, scen = t$scen, model = t$model, N = t$N, ntime = t$ntime, label = t$label)]
  list(chunk = ch, deg = degk, q9 = q9)
}
`%||%` <- function(a, b) ifelse(is.na(a), "", a)
R <- parallel::mclapply(seq_len(nrow(tasks)), per_chunk, mc.cores = 40, mc.preschedule = FALSE)
bad <- !vapply(R, is.list, TRUE); if (any(bad)) { print(head(R[bad])); stop(sprintf("%d chunk readers failed", sum(bad))) }
CH <- rbindlist(lapply(R, `[[`, "chunk"), fill = TRUE); DEG <- rbindlist(lapply(R, `[[`, "deg")); Q9 <- rbindlist(lapply(R, `[[`, "q9"))
fwrite(CH, file.path(TAB, "V4_S3_chunk_checks.csv"))
out <- list()
out$Q6 <- CH[, .(expected = .N, exist = sum(exists), sims_ok = sum(sims_ok, na.rm = TRUE), every_sim_rows = sum(every_sim_has_rows, na.rm = TRUE),
                 gen_rows_match = sum(gen_rows_only_on_failures, na.rm = TRUE), gen_fail = sum(n_gen_fail, na.rm = TRUE),
                 trace_analysis_rows = sum(n_trace_analysis, na.rm = TRUE)), by = block]
out$Q4 <- CH[block %in% c("F1", "F2", "C"), .(checked = sum(design_ref_checked, na.rm = TRUE), identical = sum(design_ref_identical, na.rm = TRUE)), by = block]
out$Q5 <- CH[block == "L", .(chunks = .N, gee_identical = sum(gee_identical, na.rm = TRUE), gee_sims_equal = sum(gee_sims_equal, na.rm = TRUE), rows_v4 = sum(gee_rows_v4, na.rm = TRUE), rows_orig = sum(gee_rows_orig, na.rm = TRUE))]
# Q7: degenerate flags vs S2 keys
exp_main <- deg_s2[model %in% 3:5]; exp_c <- deg_s2[model == 1 & N <= 45 & ntime <= 6]
cellkeys <- unique(tasks[, .(dgp, scen, model, N, ntime, label)])
q7 <- rbindlist(lapply(seq_len(nrow(cellkeys)), function(k) { c0 <- cellkeys[k]
  got <- DEG[dgp == c0$dgp & scen == c0$scen & N == c0$N & ntime == c0$ntime & label == c0$label, sort(sim_id)]
  e <- if (c0$scen == "M1c") exp_c[N == c0$N & ntime == c0$ntime, sort(sim_id)] else exp_main[model == c0$model & N == c0$N & ntime == c0$ntime, sort(sim_id)]
  data.table(c0, n_flag = length(got), n_expected = length(e), bijection = identical(as.integer(got), as.integer(e))) }))
out$Q7 <- q7[, .(cells = .N, bijection = sum(bijection), flagged = sum(n_flag), expected = sum(n_expected)), by = .(block = ifelse(scen == "M1c", "C", "main"))]
out$Q7_genfail_nondeg <- CH[, .(gen_fail_on_nondegenerate = sum(n_gen_fail_nondeg, na.rm = TRUE)), by = block]
fwrite(q7, file.path(TAB, "V4_S3_q7_degenerate.csv"))
# Q8
out$Q8_re <- CH[, .(chunks = .N, re_ok = sum(re_ok, na.rm = TRUE), lrt_rows = sum(n_lrt_rows, na.rm = TRUE)), by = block]
q8 <- CH[exists == TRUE, .(pairs = sum(q8_pairs), close = sum(q8_close)), by = .(dgp, scen)][, share := close / pairs][, pass := share >= 0.95]
out$Q8 <- q8; fwrite(q8, file.path(TAB, "V4_S3_q8_lrt_vs_wald.csv"))
# Q9
for (v in c("U_old", "U_new", "R_old", "R_new")) Q9[is.na(get(v)), (v) := FALSE]
S9 <- Q9[, { du <- as.numeric(U_new) - as.numeric(U_old); b <- U_new & U_old; dr <- as.numeric(R_new[b]) - as.numeric(R_old[b])
  .(n = .N, usable_old = mean(U_old), usable_new = mean(U_new), d_usable = mean(du), mcse_d_usable = sd(du) / sqrt(.N),
    t1_old = sum(R_old) / max(sum(U_old), 1), t1_new = sum(R_new) / max(sum(U_new), 1),
    n_both = sum(b), d_t1 = if (sum(b)) mean(dr) else NA_real_, mcse_d_t1 = if (sum(b) > 1) sd(dr) / sqrt(sum(b)) else NA_real_) },
  by = .(model, N, ntime, method)]
S9[, `:=`(z_usable = d_usable / mcse_d_usable, z_t1 = d_t1 / mcse_d_t1)]
fwrite(S9, file.path(TAB, "V4_S3_q9_h0_paired.csv"))
k <- S9[model %in% 3:4]
out$Q9 <- data.table(measure = c("usable", "type I"),
  evaluable = c(k[mcse_d_usable > 0, .N], k[mcse_d_t1 > 0, .N]),
  within_3 = c(k[mcse_d_usable > 0 & abs(z_usable) <= 3, .N], k[mcse_d_t1 > 0 & abs(z_t1) <= 3, .N]))[, share := within_3 / evaluable][, pass := share >= 0.99]
out$Q9_M5_explore <- S9[model == 5, .(evaluable_t1 = sum(mcse_d_t1 > 0, na.rm = TRUE), over3_t1 = sum(abs(z_t1) > 3, na.rm = TRUE),
                                       evaluable_u = sum(mcse_d_usable > 0, na.rm = TRUE), over3_u = sum(abs(z_usable) > 3, na.rm = TRUE))]
options(width = 220)
for (n in names(out)) { cat("\n====", n, "\n"); print(out[[n]]) }
saveRDS(out, file.path(TAB, "V4_S3_post_checks.rds"))
