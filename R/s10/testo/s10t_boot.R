# s10t_boot.R -- S10 "Testo": block bootstrap written from Section 2.5 (B = 2000, seed 20261015): replicate indices are
# resampled jointly across all cells that share a seed (seed 20260805 + 1000 m + i: M1c shares model 1; the four
# mechanisms share the seeds of their scenario). Input: results_v4/S10/testo/ind_boot.qs2 (from s10t_dp.R).
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages({library(data.table); library(qs2)}); setDTthreads(2L)
setwd("/home/user/MMM"); OUT <- "results_v4/S10/testo"
t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%.1f min] %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), s))
B <- 2000L; DPK <- c("dgp", "scen", "N", "ntime")
BI <- qs_read(file.path(OUT, "ind_boot.qs2")); tm(paste("loaded", paste(names(BI), collapse = ",")))
R <- qs_read(file.path(OUT, "x_recs.qs2")); R[, nd := generated & !design_degenerate]
fam <- function(scen) fifelse(scen == "M1c", 1L, as.integer(sub("^M", "", scen)))
CELLS <- unique(R[, .(dgp, scen, N, ntime, label)])[order(dgp, scen, N, ntime, label)]; CELLS[, fam := fam(scen)]; CELLS[, cid := .I]
set.seed(20261015L)
W <- lapply(1:7, function(m) { Cm <- matrix(0, 1000L, B + 1L); Cm[, 1] <- 1; for (b in 1:B) Cm[, b + 1] <- tabulate(sample.int(1000L, 1000L, TRUE), 1000L); Cm })
tot <- function(X, col) {   # cells x (B+1) bootstrap totals of a per-replicate 0/1 column
  out <- matrix(NA_real_, nrow(CELLS), B + 1L)
  X <- merge(X[, c("dgp", "scen", "N", "ntime", "label", "sim_id", col), with = FALSE], CELLS, by = c("dgp", "scen", "N", "ntime", "label"))
  for (m in 1:7) { cc <- CELLS[fam == m]$cid; x <- X[fam == m]; M <- matrix(0, length(cc), 1000L)
    M[cbind(match(x$cid, cc), x$sim_id)] <- as.numeric(x[[col]]); out[cc, ] <- M %*% W[[m]] }
  out }
ND <- tot(R, "nd"); tm("nd")
TT <- list(); for (s in names(BI)) { TT[[s]] <- list(n = tot(BI[[s]], "use"), kn = tot(BI[[s]], "rn"), kt = tot(BI[[s]], "rt")); tm(s) }
rm(BI); invisible(gc())
h0 <- CELLS$label == "manual"; h1 <- CELLS$label == "mmlb"
K0 <- CELLS[h0]; K1 <- CELLS[h1]; stopifnot(identical(K0[, ..DPK], K1[, ..DPK]))
rate <- function(s, ref = "kt") 100 * TT[[s]][[ref]][h0, ] / TT[[s]]$n[h0, ]
pu <- function(s) 100 * TT[[s]]$n[h0, ] / ND[h0, ]
pu_both <- function(s) 100 * (TT[[s]]$n[h0, ] + TT[[s]]$n[h1, ]) / (ND[h0, ] + ND[h1, ])
upow <- function(s) 100 * TT[[s]]$kt[h1, ] / ND[h1, ]
cm <- function(M, sel) colMeans(M[sel, , drop = FALSE], na.rm = TRUE)
summ <- function(v) c(est = v[1], mcse = sd(v[-1]), lo = unname(quantile(v[-1], 0.025)), hi = unname(quantile(v[-1], 0.975)))
OUTL <- list(); put <- function(group, region, q, v, n_dp = NA) OUTL[[length(OUTL) + 1]] <<- data.table(group = group, region = region, q = q, n_dp = n_dp, t(summ(v)))
# ---- direct comparison (Section 3.4, Figure 1): 2,376 design points, t reference ----------------------------------
G <- "gee_ar1_fg"; M <- "st_mmm"
RG <- list(all = rep(TRUE, nrow(K0)), N10_T3p = K0$N == 10 & K0$ntime >= 3, N15p_T2 = K0$N >= 15 & K0$ntime == 2, N15p_T3p = K0$N >= 15 & K0$ntime >= 3,
           "map|10|2" = K0$N == 10 & K0$ntime == 2, small = K0$N <= 25 & K0$ntime <= 3)
rg <- rate(G); rm_ <- rate(M); pg <- pu(G); pm <- pu(M); ug <- upow(G); um <- upow(M)
for (r in names(RG)) { sel <- RG[[r]]; nd <- sum(sel)
  put("h2h", r, "mad_gee", cm(abs(rg - 5), sel), nd); put("h2h", r, "mad_mmm", cm(abs(rm_ - 5), sel), nd)
  put("h2h", r, "mad_gee_minus_mmm", cm(abs(rg - 5), sel) - cm(abs(rm_ - 5), sel), nd)
  put("h2h", r, "rate_gee", cm(rg, sel), nd); put("h2h", r, "rate_mmm", cm(rm_, sel), nd)
  put("h2h", r, "pu_gee", cm(pg, sel), nd); put("h2h", r, "pu_mmm", cm(pm, sel), nd); put("h2h", r, "pu_gee_minus_mmm", cm(pg, sel) - cm(pm, sel), nd)
  put("h2h", r, "upow_gee", cm(ug, sel), nd); put("h2h", r, "upow_mmm", cm(um, sel), nd); put("h2h", r, "upow_mmm_minus_gee", cm(um, sel) - cm(ug, sel), nd) }
tm("h2h")
# ---- replication: resamples with at least one design point clearly above 7.5% (Section 4.5) ------------------------
hi_count <- function(s) { n <- TT[[s]]$n[h0, ]; k <- TT[[s]]$kt[h0, ]
  lo <- matrix(qbeta(0.025, pmax(k, 1e-300), pmax(n - k + 1, 1e-300)), nrow(n)); lo[k == 0] <- 0
  colSums(n >= 20 & lo > 0.075, na.rm = TRUE) }
for (s in c(M, G)) { h <- hi_count(s); put("s37", "all", paste0("hi_count|", s), h)
  OUTL[[length(OUTL) + 1]] <- data.table(group = "s37", region = "all", q = paste0("p_any_hi|", s), n_dp = nrow(K0), est = 100 * mean(h[-1] >= 1), mcse = NA, lo = NA, hi = NA) }
tm("s37")
# ---- Table 3 MCSE: normal reference, main grid M1-M5, by mechanism -------------------------------------------------
T3M <- c(st_gee_fg = "FG", st_gee_md = "MD", st_mmm = "MMM", st_glmer = "GLMM", st_shrink = "Shrink", st_gee = "GEE", st_glmer_lrt = "LRT", gee_exch_kc = "KC")
for (s in names(T3M)) for (g in c("GenBinaryY", "GLMM-native", "GLMM-AR1")) {
  sel <- K0$dgp == g & K0$scen %in% paste0("M", 1:5); n <- TT[[s]]$n[h0, ][sel, ]; k <- TT[[s]]$kn[h0, ][sel, ]
  put("t3", g, paste0("rate|", T3M[[s]]), 100 * colSums(k) / colSums(n), sum(sel))
  Ns <- K0$N[sel]; rN <- 100 * rowsum(k, Ns) / rowsum(n, Ns)
  put("t3", g, paste0("worst|", T3M[[s]]), apply(abs(rN - 5), 2, max), sum(sel)); put("t3", g, paste0("mad|", T3M[[s]]), colMeans(abs(rN - 5)), sum(sel)) }
tm("t3")
# ---- Table 4 MCSE: usable share pooled over the main grid by random-effects structure and T --------------------------
T4M <- c(st_glmer = "GLMM", st_glmer_lrt = "LRT", st_mmm = "MMM", st_gee_fg = "FG", st_gee = "GEE")
for (s in names(T4M)) for (re in c("RI", "RS")) for (tt in c(2, 3, 4, 6, 8, 10)) {
  sc <- if (re == "RI") c("M1", "M2") else c("M3", "M4", "M5")
  sel <- CELLS$scen %in% sc & CELLS$ntime == tt & CELLS$dgp != "Copula-AR1"
  put("t4", paste(re, tt), paste0("pu_both3|", T4M[[s]]), 100 * colSums(TT[[s]]$n[sel, ]) / colSums(ND[sel, ]), sum(sel) / 2)
  sel0 <- sel & h0; put("t4", paste(re, tt), paste0("pu_h0_3|", T4M[[s]]), 100 * colSums(TT[[s]]$n[sel0, ]) / colSums(ND[sel0, ]), sum(sel0))
  sel4 <- CELLS$scen %in% sc & CELLS$ntime == tt; put("t4", paste(re, tt), paste0("pu_both4|", T4M[[s]]), 100 * colSums(TT[[s]]$n[sel4, ]) / colSums(ND[sel4, ]), sum(sel4) / 2) }
tm("t4")
# ---- working correlation of FG (fewclust), 1,896 points of the three study mechanisms, t reference; copula vs study ----
st3 <- K0$dgp != "Copula-AR1"
put("wc", "1896", "FG ind - exch", cm(rate("gee_ind_fg") - rate("gee_exch_fg"), st3), sum(st3))
put("wc", "1896", "FG ar1 - exch", cm(rate("gee_ar1_fg") - rate("gee_exch_fg"), st3), sum(st3))
cop <- which(K0$dgp == "Copula-AR1"); key <- paste(K0$scen, K0$N, K0$ntime)
cmp <- function(s) { r <- rate(s); a <- r[cop, , drop = FALSE]
  b <- t(sapply(cop, function(i) colMeans(r[which(key == key[i] & K0$dgp %in% c("GenBinaryY", "GLMM-native", "GLMM-AR1")), , drop = FALSE], na.rm = TRUE)))
  colMeans(a, na.rm = TRUE) - colMeans(b, na.rm = TRUE) }
for (s in c("st_glmer", "st_mmm", "gee_ar1_fg", "flip")) put("cop", "480", paste0("copula - study|", s), cmp(s), length(cop))
tm("misc")
fwrite(rbindlist(OUTL), file.path(OUT, "boot_summary.csv")); tm("done")
