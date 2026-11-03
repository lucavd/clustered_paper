# =============================================================================
# s4_load.R -- S4: the v4 row set (D-15). results_full/ rows come from
# results_v4/S2/rows_timetrt.rds (V-15, md5 checked); S3 chunks are read here.
# Output (not tracked; md5 in the S4 report):
#   results_v4/S4/rows_v4.rds      time:trt rows of the v4 set, with U0/U1/U2
#   results_v4/S4/lrt_v3_m3m5.rds  glmer_lrt rows of results_full/ in M3-M5 (D-08 effect only)
#   results_v4/S4/sims_s3.rds      per-replicate records of the S3 chunks
# =============================================================================
source(path.expand("~/MMM/R/s4/s4_common.R"))
t0 <- Sys.time()
f2 <- file.path(REPO, "results_v4/S2/rows_timetrt.rds")
md5 <- unname(tools::md5sum(f2)); cat("rows_timetrt.rds md5", md5, "\n")
stopifnot(md5 == "092a5ad60b60e0ea28cb6b78ba2adbeb")
A <- readRDS(f2)
A <- A[method != "gee_ss"]                       # 'fit returned NULL' rows: absent = unusable (as S2)
A[, scen := paste0("M", model)]; A[, src := "full"]
A[, `:=`(chunk = NULL, D = NULL)]
cat(sprintf("results_full rows (time:trt or failure, gee_ss dropped): %d\n", nrow(A)))

# ---- S3 chunks ---------------------------------------------------------------
fs <- list.files(file.path(REPO, "results_v4/S3"), pattern = "[.]rds$", recursive = TRUE, full.names = TRUE)
stopifnot(length(fs) == 18720L)
rd <- function(f) {
  x <- as.data.table(readRDS(f)); s <- as.data.table(attr(x, "sims"))
  x <- x[is.na(term) | term == "time:trt",
         .(method, term, estimate, se, stat, p, converged, status, sim_id = as.integer(sim_id),
           model = as.integer(model), N = as.integer(N), ntime = as.integer(ntime), label, dgp,
           scen = scenario, block, design_degenerate, re_structure)]
  s[, `:=`(file = basename(f), block = basename(dirname(f)))]
  list(x = x, s = s)
}
L <- parallel::mclapply(fs, rd, mc.cores = 40, mc.preschedule = TRUE)
bad <- which(vapply(L, function(z) !is.list(z) || !is.data.table(z$x), logical(1)))
if (length(bad)) stop("unreadable S3 chunks: ", length(bad))
S <- rbindlist(lapply(L, `[[`, "x")); SIMS <- rbindlist(lapply(L, `[[`, "s")); rm(L); gc(verbose = FALSE)
S <- S[method != "gee_ss"]
cat(sprintf("S3 rows read: %d (F1 %d, F2 %d, L %d, C %d)\n", nrow(S), S[block == "F1", .N], S[block == "F2", .N], S[block == "L", .N], S[block == "C", .N]))

# ---- substitutions (D-15) ----------------------------------------------------
key_f12 <- unique(S[block %in% c("F1", "F2"), .(dgp, scen, N, ntime, label)])
LV3 <- A[method == LRT & model %in% 3:5]          # set aside BEFORE the F1/F2 removal (all M3-M5 cells)
n_before <- nrow(A)
A <- A[!key_f12, on = c("dgp", "scen", "N", "ntime", "label")]
cat(sprintf("rows of results_full removed for F1/F2 cells: %d (cells %d)\n", n_before - nrow(A), nrow(key_f12)))
A <- A[!(method == LRT & model %in% 3:5)]
cat(sprintf("glmer_lrt rows of results_full in M3-M5 set aside: %d\n", nrow(LV3)))
S[, src := block]
V <- rbindlist(list(A, S[block %in% c("F1", "F2", "C")],
                    S[block == "L" & method == LRT]), use.names = TRUE, fill = TRUE)
rm(A); gc(verbose = FALSE)
# D-16: glmer_bobyqa = estimate/se of the full model of the LRT, Wald p
B <- V[method == LRT & !is.na(term)]
B[, `:=`(method = SENS, p = 2 * pnorm(-abs(estimate / se)), stat = estimate / se)]
V <- rbindlist(list(V, B), use.names = TRUE)
add_usability(V)
LV3[, re_structure := "(1|id)"]; add_usability(LV3)
setkeyv(V, c(CELL, "method", "sim_id"))
saveRDS(V, file.path(OUT_S4, "rows_v4.rds"), compress = FALSE)
saveRDS(LV3, file.path(OUT_S4, "lrt_v3_m3m5.rds"), compress = FALSE)
saveRDS(SIMS, file.path(OUT_S4, "sims_s3.rds"), compress = FALSE)
cat(sprintf("v4 rows %d | by src: %s\n", nrow(V), paste(V[, .N, by = src][, sprintf("%s %d", src, N)], collapse = ", ")))
cat(sprintf("elapsed %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
