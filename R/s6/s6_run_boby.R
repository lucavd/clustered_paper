# =============================================================================
# s6_run_boby.R -- S6 follow-up: the claim engine on M1c with glmer replaced by glmer_bobyqa
# (sensitivity of D-16), on every bootstrap replicate. Needed because in M1c 2-5.5% of the
# usable null-arm glmer fits (default optimizer) are degenerate (p ~ 0, tiny SE, max|grad|
# warnings), which moves the size-adjusted critical value. Output: results_v4/S6/s6_boot_m1c_boby.rds
# =============================================================================
source(path.expand("~/MMM/R/s6/s6_common.R"))
t0 <- Sys.time()
AG <- readRDS(file.path(OUT_S5, "agg.rds")); META <- AG$META; B <- AG$B
ksa_col <- match(META$j, as.integer(colnames(AG$KSA))); PBIN <- readRDS(file.path(OUT_S5, "pbin.rds"))
mk <- function(b) { cl <- META[, .(d, dgp, scen, N, ntime, arm, method, set, tok)]
  cl[, `:=`(n = AG$A$n[b, ], k = AG$A$k[b, ], kt = AG$A$kt[b, ], kc = AG$A$kc[b, ], e2 = AG$A$e2[b, ], ksa = AG$KSA[b, ksa_col], nd = AG$ND[b, d])]; cl }
boby <- function(cl) { x <- cl[method != "glmer"]; x[method == "glmer_bobyqa", method := "glmer"]; x }   # as R/s5/s5_run.R
R <- parallel::mclapply(seq_len(B + 1L), function(b) tryCatch(engine(boby(relabel_m1c(mk(b))), PBIN = PBIN), error = function(e) conditionMessage(e)),
                        mc.cores = as.integer(Sys.getenv("S5_CORES", "45")), mc.preschedule = FALSE)
bad <- which(vapply(R, is.character, logical(1))); if (length(bad)) stop("engine failures: ", length(bad))
nm <- unique(unlist(lapply(R, names)))
M <- t(vapply(R, function(z) { v <- z[nm]; names(v) <- nm; v }, numeric(length(nm)))); colnames(M) <- nm
saveRDS(M, file.path(OUT_S6, "s6_boot_m1c_boby.rds"))
cat(sprintf("M1c bobyqa: %d columns x %d replicates | %.1f min\n", ncol(M), nrow(M), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
