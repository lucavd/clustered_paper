# =============================================================================
# s5_build.R -- S5: per-replicate indicators of every cell and their block-bootstrap
# aggregates, with the bootstrap index matrix of S4 (D-19, D-23).
# Input : results_v4/S4/rows_v4.rds, results_v4/S4/boot_idx.rds, results_v4/S4/truth_ar1_rep.rds,
#         tables/V4_S4_truth_by_cell.csv, tables/V4_S2_degenerate_designs.csv
# Output: results_v4/S5/agg.rds (not tracked; md5 in the report)
#   META  one row per cell (dp x arm x method x set), column j of every aggregate matrix
#   A     list of (B+1) x K matrices: n (usable), k (rejections at 0.05), kt (t reference),
#         kc (covered), e2 (sum of squared errors); row 1 = all replicates once (b = 0)
#   ND    (B+1) x n_dp matrix of non-degenerate replicates per design point
#   KSA   (B+1) x K_H1 size-adjusted rejection counts (H1 own cells, rule of D-19)
#   LIN   per-replicate sums for G3; MCN McNemar counts (b, c, n) for P3
# =============================================================================
source(path.expand("~/MMM/R/s5/s5_common.R"))
t0 <- Sys.time()
V <- readRDS(file.path(OUT_S4, "rows_v4.rds"))
V <- V[!is.na(term)]
V[, model := as.integer(sub("^M([1-5]).*$", "\\1", scen))]
V <- V[, .(dgp, scen, model, N, ntime, label, method, sim_id, estimate, se, p, U0)]
DEG <- fread(file.path(TABDIR, "V4_S2_degenerate_designs.csv"))[, .(model, N, ntime, sim_id)]
V[, deg := FALSE]; V[DEG, on = c("model", "N", "ntime", "sim_id"), deg := TRUE]
V[, u0 := U0 & !deg]; V[, U0 := NULL]
V[, arm := fifelse(label == "mmlb", "H1", "H0")]
cat(sprintf("rows %d | load %.1f s\n", nrow(V), as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# ---- truths (agg_v4.R rules: D-12, D-18) ----------------------------------------
TT <- fread(file.path(TABDIR, "V4_S4_truth_by_cell.csv"))[, .(dgp, scen, ntime, label, tm = truth_marginal, tc = truth_conditional)]
TR <- readRDS(file.path(OUT_S4, "truth_ar1_rep.rds"))[, .(dgp, scen, N, ntime, label, sim_id, tr = truth_cond)]
V[TT, on = c("dgp", "scen", "ntime", "label"), `:=`(tm = i.tm, tc = i.tc)]
V[TR, on = c("dgp", "scen", "N", "ntime", "label", "sim_id"), tr := i.tr]
V[, tok := method %in% MARGINAL_METHODS | (method %in% COND_METHODS & dgp != "GenBinaryY")]
V[, th := fifelse(method %in% MARGINAL_METHODS, tm, fifelse(dgp == "GLMM-native", tc, tr))]
V[tok == FALSE, th := NA_real_]
V[, c("tm", "tc", "tr") := NULL]
stopifnot(V[u0 & tok & !is.finite(th), .N] == 0L)

# ---- indicators ------------------------------------------------------------------
V[, r  := u0 & fifelse(is.na(p), FALSE, p < 0.05)]
V[, rt := u0 & method %in% WALD_ALL & fifelse(is.finite(estimate / se), 2 * pt(-abs(estimate / se), df = N - NP[model]) < 0.05, FALSE)]
V[, cv := u0 & tok & fifelse(is.finite(th), abs(estimate - th) <= Z975 * se, FALSE)]
V[, e2 := fifelse(u0 & tok, (estimate - th)^2, 0)]
V[, za := abs(estimate / se)]   # DEC-T1 (S8b, D-73): |estimate/se| kept for the size-adjusted rule below
V[, c("estimate", "se", "th") := NULL]
CK <- c(DPK, "label", "sim_id")
INT <- V[method %in% PRIM8, .(inter = sum(u0) == 8L), by = CK]
V[INT, on = CK, inter := i.inter]; V[is.na(inter), inter := FALSE]
PR <- V[method %in% c("mmm", "gee_fg"), .(pair = sum(u0) == 2L), by = CK]
V[PR, on = CK, pair := i.pair]; V[is.na(pair), pair := FALSE]
rm(INT, PR); gc(verbose = FALSE)

# ---- McNemar counts for P3 (main grid, H0, pairs usable for both) --------------------
W <- dcast(V[arm == "H0" & scen != "M1c" & pair == TRUE & method %in% c("mmm", "gee_fg")],
           dgp + scen + N + ntime + sim_id ~ method, value.var = "r")
MCN <- W[, .(b = sum(mmm & !gee_fg), c = sum(!mmm & gee_fg), n = .N), by = dgp]
rm(W); print(MCN)

# ---- cell layout -------------------------------------------------------------------
DP <- unique(V[, ..DPK]); setorderv(DP, DPK); DP[, d := .I]
META <- rbind(
  CJ(d = DP$d, arm = c("H0", "H1"), method = METHODS, set = c("own", "intersection")),
  CJ(d = DP$d, arm = c("H0", "H1"), method = c("mmm", "gee_fg"), set = "pair"))
META <- merge(META, DP, by = "d"); setorder(META, model, d, arm, method, set); META[, j := .I]
META[, tok := method %in% MARGINAL_METHODS | (method %in% COND_METHODS & dgp != "GenBinaryY")]
cat(sprintf("design points %d | cells %d\n", nrow(DP), nrow(META)))

# ---- bootstrap counts, exactly as agg_v4.R ------------------------------------------------
BIDX <- readRDS(file.path(OUT_S4, "boot_idx.rds"))
CNT <- lapply(BIDX, function(M) { C <- matrix(0L, B_FULL, 1000L); for (b in seq_len(B_FULL)) C[b, ] <- tabulate(M[b, ], 1000L); C })
CNT <- lapply(CNT, function(C) rbind(1L, C[seq_len(B_USE), , drop = FALSE]))     # row 1: b = 0
storage <- function(C) { storage.mode(C) <- "double"; C }
CNTd <- lapply(CNT, storage)

TYPES <- c(n = "u0", k = "r", kt = "rt", kc = "cv", e2 = "e2")
A <- setNames(vector("list", length(TYPES)), names(TYPES)); for (t in names(TYPES)) A[[t]] <- vector("list", 5)
ND <- vector("list", 5); LIN <- vector("list", 5)
for (m in 1:5) {
  Mm <- META[model == m]; jl <- setNames(seq_len(nrow(Mm)), Mm$j)
  Vm <- V[model == m]
  X <- lapply(TYPES, function(z) matrix(0, 1000L, nrow(Mm)))
  for (s in c("own", "intersection", "pair")) {
    Vs <- if (s == "own") Vm else if (s == "intersection") Vm[inter == TRUE] else Vm[pair == TRUE & method %in% c("mmm", "gee_fg")]
    Vs <- Vs[u0 == TRUE]
    key <- Mm[set == s, .(dgp, scen, N, ntime, arm, method, j)]
    Vs <- key[Vs, on = c("dgp", "scen", "N", "ntime", "arm", "method"), nomatch = NULL]
    stopifnot(!anyNA(Vs$j))
    idx <- cbind(Vs$sim_id, jl[as.character(Vs$j)])
    for (t in names(TYPES)) X[[t]][idx] <- as.numeric(Vs[[TYPES[[t]]]])
  }
  for (t in names(TYPES)) A[[t]][[m]] <- CNTd[[m]] %*% X[[t]]
  DPm <- DP[model == m]
  NDX <- matrix(1, 1000L, nrow(DPm))
  Dg <- DEG[model == m]
  if (nrow(Dg)) {
    hit <- merge(DPm[, .(d, N, ntime)], Dg, by = c("N", "ntime"), allow.cartesian = TRUE)
    NDX[cbind(hit$sim_id, match(hit$d, DPm$d))] <- 0
  }
  ND[[m]] <- CNTd[[m]] %*% NDX
  # G3: per-replicate sums of the main-grid H0 own pooled rates and of the pairs
  g <- Mm[scen != "M1c" & arm == "H0" & method %in% PRIM8 & set %in% c("own", "pair")]
  LIN[[m]] <- g[, .(a = list(rowSums(X$k[, jl[as.character(j)], drop = FALSE])),
                    u = list(rowSums(X$n[, jl[as.character(j)], drop = FALSE]))), by = .(dgp, method, set)]
  rm(X, NDX); gc(verbose = FALSE)
  cat(sprintf("model %d done (%d cells) %.1f min\n", m, nrow(Mm), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
A <- lapply(A, function(L) do.call(cbind, L))
NDf <- matrix(NA_real_, B_USE + 1L, nrow(DP))
for (m in 1:5) NDf[, DP[model == m, d]] <- ND[[m]]           # column d = design point d
ND <- NDf; rm(NDf); stopifnot(!anyNA(ND), identical(DP$d, seq_len(nrow(DP))))
stopifnot(identical(META$j, seq_len(nrow(META))), all(diff(META$model) >= 0))   # A columns = META$j

# ---- size-adjusted power, bootstrap rule of agg_v4.R (D-19) -------------------------
sa_one <- function(C, p0, s0, p1, s1) {
  # C: B x 1000 counts (without row b = 0). Returns k_sa for b = 0..B and n1, n0 for b = 0..B
  o0 <- s0[order(p0[s0])]; v0 <- p0[o0]; o1 <- s1[order(p1[s1])]; v1 <- p1[o1]
  nb <- nrow(C) + 1L
  if (length(o0) < 2 || length(o1) < 1) return(rep(NA_real_, nb))
  q0 <- quantile(p0[s0], 0.05, type = 7, names = FALSE); k0 <- sum(p1[s1] < q0)
  rc <- function(M) if (ncol(M) == 1L) M else t(apply(M, 1, cumsum))
  cs0 <- rc(C[, o0, drop = FALSE]); n0 <- cs0[, ncol(cs0)]
  h <- (n0 - 1) * 0.05 + 1; lo <- floor(h); fr <- h - lo
  j_lo <- rowSums(cs0 < lo) + 1L; j_hi <- rowSums(cs0 < pmin(lo + 1, n0)) + 1L
  q <- v0[j_lo] + fr * (v0[j_hi] - v0[j_lo]); q[n0 < 2] <- NA
  cs1 <- cbind(0L, rc(C[, o1, drop = FALSE]))
  k_sa <- cs1[cbind(seq_len(nrow(C)), findInterval(q, v1, left.open = TRUE) + 1L)]
  k_sa[is.na(q)] <- NA
  c(k0, k_sa)
}
H1c <- META[arm == "H1" & set == "own"]
# DEC-T1 (S8b, D-73): size-adjusted rule of Section 2.5 = 95th percentile (type 7) of the usable null |estimate/se|.
# sa_one is unchanged: it receives -|z|; type 7 is reflection-symmetric, so its 5th percentile is -q95(|z0|) and
# p1 < q0 <=> |z1| > q95. The likelihood-ratio test has no z and stays on p.
PV <- V[, .(dgp, scen, N, ntime, arm, method, sim_id, p = fifelse(method == "glmer_lrt", p, -za), u0)]
setkey(PV, dgp, scen, N, ntime, method, arm)
jobs <- split(H1c[, .(dgp, scen, model, N, ntime, method, j)], by = c("dgp", "scen", "N", "ntime"))
KS <- parallel::mclapply(jobs, function(J) {
  k <- J[1]; C <- CNT[[k$model]][-1, , drop = FALSE]
  out <- lapply(seq_len(nrow(J)), function(r) {
    mth <- J$method[r]
    g0 <- PV[.(k$dgp, k$scen, k$N, k$ntime, mth, "H0")]; g1 <- PV[.(k$dgp, k$scen, k$N, k$ntime, mth, "H1")]
    p0 <- rep(NA_real_, 1000); p0[g0$sim_id] <- g0$p; s0 <- g0$sim_id[g0$u0 %in% TRUE]
    p1 <- rep(NA_real_, 1000); p1[g1$sim_id] <- g1$p; s1 <- g1$sim_id[g1$u0 %in% TRUE]
    sa_one(C, p0, s0, p1, s1)
  })
  M <- do.call(cbind, out); colnames(M) <- J$j; M
}, mc.cores = as.integer(Sys.getenv("S5_CORES", "40")), mc.preschedule = FALSE)
bad <- which(vapply(KS, function(z) !is.matrix(z), logical(1)))
if (length(bad)) { print(head(KS[bad])); stop("size-adjusted failures: ", length(bad)) }
KSA <- do.call(cbind, KS); KSA <- KSA[, order(as.integer(colnames(KSA))), drop = FALSE]
cat(sprintf("KSA %d x %d | %.1f min\n", nrow(KSA), ncol(KSA), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
saveRDS(list(META = META, DP = DP, A = A, ND = ND, KSA = KSA, LIN = LIN, MCN = MCN, B = B_USE),
        file.path(OUT_S5, tag_out("agg.rds")), compress = FALSE)
cat(sprintf("saved | elapsed %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
