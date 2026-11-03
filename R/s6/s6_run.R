# =============================================================================
# s6_run.R -- S6: the S5 claim engine and the S6 extras on every bootstrap replicate
# (b = 0 is the estimate), for the main grid (1,440 design points) and for M1c (72),
# with the S5 aggregates and bootstrap indices. Output: results_v4/S6/s6_boot.rds
# =============================================================================
source(path.expand("~/MMM/R/s6/s6_common.R")); source(file.path(REPO, "R/s6/s6_extras.R"))
t0 <- Sys.time()
AG <- readRDS(file.path(OUT_S5, "agg.rds")); META <- AG$META; B <- AG$B
cat(sprintf("agg loaded (B = %d) %.1f min\n", B, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
ksa_col <- match(META$j, as.integer(colnames(AG$KSA)))
PBIN <- readRDS(file.path(OUT_S5, "pbin.rds"))
mk <- function(b) {                                   # identical to R/s5/s5_run.R
  cl <- META[, .(d, dgp, scen, N, ntime, arm, method, set, tok)]
  cl[, `:=`(n = AG$A$n[b, ], k = AG$A$k[b, ], kt = AG$A$kt[b, ], kc = AG$A$kc[b, ], e2 = AG$A$e2[b, ],
            ksa = AG$KSA[b, ksa_col], nd = AG$ND[b, d])]
  cl
}
one <- function(b) {
  cl <- mk(b); c1 <- relabel_m1c(cl)
  list(main = engine(cl, PBIN = PBIN), m1c = engine(c1, PBIN = PBIN),
       xmain = extras(cl[scen != "M1c"]), xm1c = extras(c1))
}
R <- parallel::mclapply(seq_len(B + 1L), function(b) tryCatch(one(b), error = function(e) conditionMessage(e)),
                        mc.cores = as.integer(Sys.getenv("S5_CORES", "45")), mc.preschedule = FALSE)
bad <- which(vapply(R, is.character, logical(1)))
if (length(bad)) { print(head(unlist(R[bad]))); stop("engine failures: ", length(bad)) }
tomat <- function(part) { nm <- unique(unlist(lapply(R, function(z) names(z[[part]]))))
  M <- t(vapply(R, function(z) { v <- z[[part]][nm]; names(v) <- nm; v }, numeric(length(nm)))); colnames(M) <- nm; M }
OUT <- list(main = tomat("main"), m1c = tomat("m1c"), xmain = tomat("xmain"), xm1c = tomat("xm1c"), B = B)
saveRDS(OUT, file.path(OUT_S6, "s6_boot.rds"))
cat(sprintf("S6 boot: main %d, M1c %d, extras %d / %d columns; %d replicates | %.1f min\n", ncol(OUT$main), ncol(OUT$m1c),
            ncol(OUT$xmain), ncol(OUT$xm1c), nrow(OUT$main), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
