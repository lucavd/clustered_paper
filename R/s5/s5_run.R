# =============================================================================
# s5_run.R -- S5: the claim engine on every bootstrap replicate (b = 0 is the estimate).
# Variants: main (1,440 design points), bobyqa (glmer replaced by glmer_bobyqa, D-16),
# M1c (descriptive, C01-C03 and C20). Output: results_v4/S5/claims_boot.rds
# =============================================================================
source(path.expand("~/MMM/R/s5/s5_common.R")); source(file.path(REPO, "R/s5/s5_engine.R"))
t0 <- Sys.time()
AG <- readRDS(file.path(OUT_S5, tag_out("agg.rds"))); META <- AG$META; B <- AG$B
cat(sprintf("agg loaded (B = %d) %.1f min\n", B, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
ksa_col <- match(META$j, as.integer(colnames(AG$KSA)))
pbf <- file.path(OUT_S5, "pbin.rds")
if (file.exists(pbf)) PBIN <- readRDS(pbf) else {
  pl <- parallel::mclapply(1:1000, function(n) vapply(0:n, function(k) binom.test(k, n, 0.05)$p.value, 0), mc.cores = as.integer(Sys.getenv("S5_CORES", "40")))
  PBIN <- list(p = unlist(pl), off = c(0, cumsum(lengths(pl)))[1:1000]); saveRDS(PBIN, pbf) }
mk <- function(b) {
  cl <- META[, .(d, dgp, scen, N, ntime, arm, method, set, tok)]
  cl[, `:=`(n = AG$A$n[b, ], k = AG$A$k[b, ], kt = AG$A$kt[b, ], kc = AG$A$kc[b, ], e2 = AG$A$e2[b, ],
            ksa = AG$KSA[b, ksa_col], nd = AG$ND[b, d])]
  cl
}
boby <- function(cl) { x <- cl[method != "glmer"]; x[method == "glmer_bobyqa", method := "glmer"]; x }
# C05 pairs: winner and runner-up per design point at b = 0 (8 primary methods, order of PRIM8 on ties)
cl0 <- mk(1L)
D0 <- cl0[scen != "M1c" & arm == "H0" & set == "own" & method %in% PRIM8 & n > 0][, dev := abs(100 * k / n - 5)]
D0[, mo := match(method, PRIM8)]; setorder(D0, dgp, scen, N, ntime, dev, mo)
PAIRS <- D0[, .(method_w = method[1], method_r = method[2]), by = .(dgp, scen, N, ntime)][!is.na(method_r)]
PAIRS[, pid := .I]
one <- function(b) {
  cl <- mk(b)
  m <- engine(cl, PBIN = PBIN, pairs_c05 = PAIRS)
  list(main = m, c05 = attr(m, "C05"), boby = engine(boby(cl), PBIN = PBIN), m1c = engine(cl, grid = "M1c"))
}
R <- parallel::mclapply(seq_len(B + 1L), function(b) tryCatch(one(b), error = function(e) conditionMessage(e)),
                        mc.cores = as.integer(Sys.getenv("S5_CORES", "45")), mc.preschedule = FALSE)
bad <- which(vapply(R, is.character, logical(1)))
if (length(bad)) { print(head(unlist(R[bad]))); stop("engine failures: ", length(bad)) }
tomat <- function(part) { nm <- unique(unlist(lapply(R, function(z) names(z[[part]]))))
  M <- t(vapply(R, function(z) { v <- z[[part]][nm]; names(v) <- nm; v }, numeric(length(nm)))); colnames(M) <- nm; M }
OUT <- list(main = tomat("main"), boby = tomat("boby"), m1c = tomat("m1c"),
            c05 = do.call(rbind, lapply(R, `[[`, "c05")), pairs = PAIRS, B = B)
saveRDS(OUT, file.path(OUT_S5, tag_out("claims_boot.rds")))
cat(sprintf("claims: main %d, bobyqa %d, M1c %d columns; %d replicates | %.1f min\n", ncol(OUT$main), ncol(OUT$boby), ncol(OUT$m1c),
            nrow(OUT$main), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(round(OUT$main[1, grep("^C0[1-4]", colnames(OUT$main))], 3))
