## S10 'assolute' follow-up: read-only timing of the S5 claim engine on the current row set (variant nc), 4 bootstrap
## replicates on 1 core, to size an exact rerun of the threshold check. Writes nothing in the repository.
Sys.setenv(MMM_VARIANT = "nc")
for (p in c("/home/user/MMM/results_v4/variant_nc/S5/agg.rds", "/home/user/MMM/results_v4/variant_nc/S5/pbin.rds",
            "/home/user/MMM/reports/morris/variant_nc/S5_logs", "/home/user/MMM/reports/morris/variant_nc/S4_logs",
            "/home/user/MMM/results_v4/S5", "/home/user/MMM/reports/morris/S5_logs", "/home/user/MMM/results_v4/S4", "/home/user/MMM/reports/morris/S4_logs"))
  if (!file.exists(p)) stop("missing (would be created): ", p)
source(path.expand("~/MMM/R/s5/s5_common.R")); source(file.path(REPO, "R/s5/s5_engine.R"))
t0 <- Sys.time()
AG <- readRDS(file.path(OUT_S5, "agg.rds")); META <- AG$META; B <- AG$B
PBIN <- readRDS(file.path(OUT_S5, "pbin.rds")); ksa_col <- match(META$j, as.integer(colnames(AG$KSA)))
cat("loaded", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
mk <- function(b) { cl <- META[, .(d, dgp, scen, N, ntime, arm, method, set, tok)]
  cl[, `:=`(n = AG$A$n[b, ], k = AG$A$k[b, ], kt = AG$A$kt[b, ], kc = AG$A$kc[b, ], e2 = AG$A$e2[b, ], ksa = AG$KSA[b, ksa_col], nd = AG$ND[b, d])]; cl }
boby <- function(cl) { x <- cl[method != "glmer"]; x[method == "glmer_bobyqa", method := "glmer"]; x }
cl0 <- mk(1L)
D0 <- cl0[scen != "M1c" & arm == "H0" & set == "own" & method %in% PRIM8 & n > 0][, dev := abs(100 * k / n - 5)]
D0[, mo := match(method, PRIM8)]; setorder(D0, dgp, scen, N, ntime, dev, mo)
PAIRS <- D0[, .(method_w = method[1], method_r = method[2]), by = .(dgp, scen, N, ntime)][!is.na(method_r)]; PAIRS[, pid := .I]
tt <- sapply(2:5, function(b) { s <- Sys.time(); cl <- mk(b); m <- engine(cl, PBIN = PBIN, pairs_c05 = PAIRS); bb <- engine(boby(cl), PBIN = PBIN); m1 <- engine(cl, grid = "M1c")
  as.numeric(difftime(Sys.time(), s, units = "secs")) })
cat("seconds per bootstrap replicate (main + bobyqa + M1c engines):", round(tt, 1), "\n")
cat("projected s5_run on 1 core for 2,001 replicates:", round(mean(tt) * 2001 / 3600, 2), "h\n")
