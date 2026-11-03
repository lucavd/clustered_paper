# =============================================================================
# s7d_gate_reprod.R [file_list] -- gates GE2 (same datasets as the study), GE3 (same sign-flip test as the S7c pilot)
# and GE5 (geesmv / fewclust coincidence with independence) on the chunk files of results_v4/S7d/rows/
# (all files, or those listed one per line in file_list, for the dry run). Writes tables/V4_S7d_gates_reprod.csv
# (or _dry.csv) and results_v4/S7d/gates/reprod_detail*.rds.
# =============================================================================
source(path.expand("~/MMM/R/s4/s4_common.R"))
suppressMessages(library(data.table)); setDTthreads(4L)
a <- commandArgs(TRUE); OUT <- file.path(REPO, "results_v4/S7d"); dir.create(file.path(OUT, "gates"), showWarnings = FALSE)
files <- if (length(a)) readLines(a[1]) else list.files(file.path(OUT, "rows"), pattern = "\\.rds$", full.names = TRUE)
tag <- if (length(a)) "_dry" else ""
L <- lapply(files, readRDS)
ROWS <- rbindlist(lapply(L, `[[`, "rows"), fill = TRUE); RECS <- rbindlist(lapply(L, `[[`, "recs"), fill = TRUE); rm(L)
cat("chunks", length(files), "rows", nrow(ROWS), "recs", nrow(RECS), "\n")
keys <- unique(RECS[, .(dgp, scen, N, ntime, label)])
KEEP <- c("dgp", "scen", "N", "ntime", "label", "method", "sim_id", "estimate", "se")
SA <- readRDS(file.path(REPO, "results_v4/variant_cc/S4/rows_v4.rds"))[method %in% c("gee", "gee_md", "gee_kc", "gee_fg") & term %in% "time:trt", ..KEEP]
SB <- readRDS(file.path(REPO, "results_v4/S6b/ns/rows_ns.rds"))[method %in% c("gee", "gee_md", "gee_kc", "gee_fg") & term %in% "time:trt", ..KEEP]
S <- rbind(SA, SB); rm(SA, SB); invisible(gc())
S <- S[unique(RECS[, .(dgp, scen, N, ntime, label, sim_id)]), on = c("dgp", "scen", "N", "ntime", "label", "sim_id"), nomatch = NULL]
rel <- function(x, y) abs(x - y) / pmax(abs(y), 1e-12)
# GE2a: exchangeable geeglm estimate and SE vs saved "gee"
E <- merge(ROWS[cor == "exch", .(dgp, scen, N, ntime, label, sim_id, est, se_unc)], S[method == "gee"], by = c("dgp", "scen", "N", "ntime", "label", "sim_id"), all.y = TRUE)
fin <- E[is.finite(estimate) & is.finite(se)]
ge2a <- fin[, .(n = .N, missing_new = sum(!is.finite(est) | !is.finite(se_unc)), max_d_est = max(abs(est - estimate), na.rm = TRUE), max_d_se = max(abs(se_unc - se), na.rm = TRUE))]
# GE2b: geesmv exchangeable (replicates 1-50) vs saved gee_md / gee_kc / gee_fg SE
SM <- melt(ROWS[cor == "exch" & sim_id <= 50, .(dgp, scen, N, ntime, label, sim_id, gee_md = smv_md, gee_kc = smv_kc, gee_fg = smv_fg)],
           id.vars = c("dgp", "scen", "N", "ntime", "label", "sim_id"), variable.name = "method", value.name = "se_new")
SM <- merge(SM, S[method != "gee" & sim_id <= 50, .(dgp, scen, N, ntime, label, sim_id, method, se)], by = c("dgp", "scen", "N", "ntime", "label", "sim_id", "method"))
ge2b <- SM[, .(n_both = sum(is.finite(se_new) & is.finite(se)), only_saved = sum(!is.finite(se_new) & is.finite(se)), only_new = sum(is.finite(se_new) & !is.finite(se)),
               max_rel = suppressWarnings(max(rel(se_new, se)[is.finite(se_new) & is.finite(se)]))), by = method]
# GE3: sign-flip vs S7c pilot (M1, N 10/25, T 2/3/10)
pil <- RECS[scen == "M1" & N %in% c(10L, 25L) & ntime %in% c(2L, 3L, 10L)]
ge3 <- data.table(n = 0L, mism_primary = 0L, mism_strict = 0L)
if (nrow(pil)) {
  P <- rbindlist(lapply(list.files(file.path(REPO, "results_v4/S7c/posthoc"), full.names = TRUE), readRDS), fill = TRUE)[, .(dgp, N, ntime, label, sim_id, p_flip_relaxed)]
  A <- rbindlist(lapply(list.files(file.path(REPO, "results_v4/S7c/pilot"), full.names = TRUE), function(f) readRDS(f)$tests), fill = TRUE)[test == "flip", .(dgp, N, ntime, label, sim_id, pA = p)]
  M <- merge(merge(pil[, .(dgp, N, ntime, label, sim_id, p_flip, p_flip_strict)], P, by = c("dgp", "N", "ntime", "label", "sim_id"), all.x = TRUE), A, by = c("dgp", "N", "ntime", "label", "sim_id"), all.x = TRUE)
  same <- function(x, y) (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)
  ge3 <- M[, .(n = .N, mism_primary = sum(!same(p_flip, p_flip_relaxed)), mism_strict = sum(!same(p_flip_strict, pA)))]
}
# GE5: independence, geesmv vs fewclust, MD and FG
I <- ROWS[cor == "ind" & sim_id <= 50]
ge5 <- rbindlist(lapply(c("md", "kc", "fg"), function(ty) { x <- I[[paste0("smv_", ty)]]; y <- I[[paste0("se_", ty)]]; b <- is.finite(x) & is.finite(y)
  data.table(corr = ty, n_both = sum(b), only_geesmv = sum(is.finite(x) & !is.finite(y)), only_fewclust = sum(!is.finite(x) & is.finite(y)),
             max_rel = if (any(b)) max(rel(x[b], y[b])) else NA_real_, share_rel_lt_1e8 = if (any(b)) mean(rel(x[b], y[b]) < 1e-8) else NA_real_) }))
# exchangeable coincidence distribution (reported, not gate)
X <- ROWS[cor == "exch" & sim_id <= 50]
exq <- rbindlist(lapply(c("md", "kc", "fg"), function(ty) { r <- X[[paste0("smv_", ty)]] / X[[paste0("se_", ty)]]; r <- r[is.finite(r)]
  data.table(corr = ty, n = length(r), q05 = quantile(r, 0.05), q50 = quantile(r, 0.5), q95 = quantile(r, 0.95), share_within_1pct = mean(abs(r - 1) < 0.01)) }))
print(ge2a); print(ge2b); print(ge3); print(ge5); print(exq)
pass2 <- ge2a$n > 0 && ge2a$missing_new == 0 && ge2a$max_d_est < 1e-8 && ge2a$max_d_se < 1e-8 && all(ge2b$max_rel < 1e-8 | ge2b$n_both == 0)
pass3 <- ge3$mism_primary == 0 && ge3$mism_strict == 0
pass5 <- all(ge5[corr %in% c("md", "fg")]$max_rel < 1e-8)
G <- data.table(gate = c("GE2", "GE3", "GE5"), pass = c(pass2, pass3, pass5),
  detail = c(sprintf("gee exch: %d saved finite rows, %d without new value, max |d est| %.2e, max |d se| %.2e; geesmv (i<=50): %s", ge2a$n, ge2a$missing_new, ge2a$max_d_est, ge2a$max_d_se,
                     paste(sprintf("%s both %d, only saved %d, only new %d, max rel %.2e", ge2b$method, ge2b$n_both, ge2b$only_saved, ge2b$only_new, ge2b$max_rel), collapse = "; ")),
             sprintf("%d pilot replicates, mismatches primary %d, strict %d", ge3$n, ge3$mism_primary, ge3$mism_strict),
             paste(sprintf("%s: both %d, only geesmv %d, only fewclust %d, max rel %.2e", ge5$corr, ge5$n_both, ge5$only_geesmv, ge5$only_fewclust, ge5$max_rel), collapse = "; ")))
print(G)
fwrite(G, file.path(REPO, sprintf("tables/V4_S7d_gates_reprod%s.csv", tag))); fwrite(exq, file.path(REPO, sprintf("tables/V4_S7d_exch_coincidence%s.csv", tag)))
saveRDS(list(ge2a = ge2a, ge2b = ge2b, ge3 = ge3, ge5 = ge5, exq = exq), file.path(OUT, sprintf("gates/reprod_detail%s.rds", tag)))
