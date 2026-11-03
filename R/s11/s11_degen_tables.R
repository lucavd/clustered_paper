# S11 copy (D-90): reads the indicators of results_v4/S11/degen/ind (R/s11/s11_degen_ind_cv.R) and the classes of S10;
# writes tables/V8_S11_degen_*.csv and results_v4/S11/degen/. Otherwise identical to R/s10/degen/s10_degen_tables.R.
# =============================================================================
# s10_degen_tables.R -- S10 "Degeneri", steps 2b-4: (1) Tables 3, 4 and 5 recomputed from the per-replicate indicators
# (results_v4/S10/degen/ind, s10_degen_ind.R) and compared with the current files (tolerance 1e-9) BEFORE any removal;
# (2) usable results and rejections by class x strategy x mechanism x region x arm; (3) sensitivity versions
#   a: Si removed from the usable results of the Wald strategies (GEE, MMM, GLMM Wald, shrinkage); LRT, sign-flip unchanged
#   b: Si removed from every strategy      c: as a, with S instead of Si      d: as a, plus V1, E1, B0, B1
# (same denominators: nd = replicates with an estimable design), Tables 3-5 and the rule recomputed with the code of
# R/s7f/s7f_table4.R (exact Clopper-Pearson flags of R/s7f/s7f_cpflags.R, >= 20 usable), R/s5/s5_engine.R (rate, wc, mad)
# and R/s6/s6_extras.R (puT); (4) predictions PD1-PD6 judged as pre-registered (value of the change against the margin).
# Without results_v4/S10/degen/classes.rds only (1) runs. Outputs: tables/V7_S10_degen_*.csv
# =============================================================================
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(4L); setwd("/home/user/MMM")
D <- "results_v4/S11/degen"; IND <- file.path(D, "ind"); TB <- "tables"; D10 <- "results_v4/S10/degen"   # S11: D-90 copy of R/s10/degen/s10_degen_tables.R
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id"); K5 <- K6[1:5]; DPK <- K6[1:4]
src <- readLines("R/s7f/s7f_cpflags.R"); eval(parse(text = grep("^cp_(lo|up) <- function", src, value = TRUE)))   # verbatim
SHORT <- c(GenBinaryY = "GB", `GLMM-native` = "NAT", `GLMM-AR1` = "AR1")
EST14 <- c("gee_ar1_fg", "gee_ind_fg", "gee_exch_fg", "gee_ar1_md", "gee_ind_md", "gee_exch_md", "gee_ar1_kc", "gee_ind_kc",
           "gee_exch_kc", "gee_ind_cr2", "gee_exch_cr2", "gee_ar1_unc", "gee_ind_unc", "gee_exch_unc")
T5S <- c(EST14, "mmm_cc", "glmer", "shrink", "glmer_lrt", "flip")
T34S <- c("mmm", "gee", "gee_md", "gee_kc", "gee_fg", "glmer", "shrink", "glmer_lrt")
NOWALD <- c("glmer_lrt", "flip")
BITS <- c("const", "S", "Si", "V1", "E1", "B0", "B1")
HAVE_CL <- file.exists(file.path(D10, "classes.rds"))
G <- list(); gate <- function(name, ok, detail = "") { G[[length(G) + 1L]] <<- data.table(check = name, pass = isTRUE(ok), detail = detail); cat(sprintf("[%s] %s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name, detail)) }
## ---- denominators: estimable designs per point (from the regenerated datasets, else from the S7d2 rows) -------------
if (HAVE_CL) {
  CL <- readRDS(file.path(D10, "classes.rds"))
  ND <- CL[, .(nd = sum(!design_degenerate)), by = K5]
  CLB <- CL[elig == TRUE, c(K6, BITS[-1]), with = FALSE]
  CLC <- CL[const == TRUE & design_degenerate == FALSE, K6, with = FALSE][, const := TRUE]
  CLB <- merge(CLB, CLC, by = K6, all = TRUE); for (b in BITS) set(CLB, which(is.na(CLB[[b]])), b, FALSE)
  rm(CL); invisible(gc())
} else {
  f2 <- list.files("results_v4/S7d2/rows", pattern = "\\.rds$", full.names = TRUE)
  ND <- rbindlist(parallel::mclapply(f2, function(f) readRDS(f)[, .(nd = sum(!design_degenerate)), by = .(dgp, scen, N, ntime, label)], mc.cores = 4L))[, .(nd = sum(nd)), by = K5]
}
stopifnot(nrow(ND) == 4752L)
## ---- aggregate the indicators by point x arm x class profile -----------------------------------------------------------
agg_one <- function(family, s) {
  Y <- readRDS(file.path(IND, sprintf("%s_%s.rds", family, s)))
  wald <- !(s %in% NOWALD)
  Y[, L := if (wald) (abs(est) >= 5 | se >= 5) else FALSE]
  if (HAVE_CL) { Y <- merge(Y, CLB, by = K6, all.x = TRUE); for (b in BITS) set(Y, which(is.na(Y[[b]])), b, FALSE) } else for (b in BITS) set(Y, j = b, value = FALSE)
  A <- Y[, .(n = .N, k = sum(rej)), by = c(K5, BITS, "L")]
  A[, `:=`(family = family, strategy = s, wald = wald)]; A
}
AG <- rbind(rbindlist(lapply(T5S, function(s) agg_one("T5", s))), rbindlist(lapply(T34S, function(s) agg_one("T34", s))))
cat("aggregated profiles:", nrow(AG), "\n")
## ---- versions ------------------------------------------------------------------------------------------------------------
VERS <- c(base = "no removal", a = "Si removed, Wald strategies", b = "Si removed, all strategies", c = "S removed, Wald strategies",
          d = "Si, V1, E1, B0, B1 removed, Wald strategies")
drop_of <- function(v, A) switch(v, base = rep(FALSE, nrow(A)), a = A$wald & A$Si, b = A$Si, c = A$wald & A$S,
                                 d = A$wald & (A$Si | A$V1 | A$E1 | A$B0 | A$B1))
counts <- function(v, family) {
  fam <- family; A <- AG[AG$family == fam]; keep <- !drop_of(v, A)
  X <- A[keep, .(n = sum(n), k = sum(k)), by = c(K5, "strategy")]
  pts <- unique(ND[, ..K5]); if (family == "T34") pts <- pts[dgp %in% names(SHORT) & scen %in% sprintf("M%d", 1:5)]
  GR <- CJ(i = seq_len(nrow(pts)), strategy = unique(A$strategy)); GR <- cbind(pts[GR$i], GR[, .(strategy)])
  X <- merge(GR, X, by = c(K5, "strategy"), all.x = TRUE); X[is.na(n), `:=`(n = 0L, k = 0L)]
  merge(X, ND, by = K5)
}
## Table 5 (R/s7f/s7f_table4.R)
t5 <- function(v) {
  X <- counts(v, "T5")
  H0 <- X[label == "manual"]; H0[, `:=`(assess = n >= 20, pu = 100 * n / nd, rate = fifelse(n > 0, 100 * k / n, NA_real_))]
  H0[, `:=`(hi = assess & 100 * cp_lo(k, n) > 7.5, lo = assess & 100 * cp_up(k, n) < 2.5)]; H0[is.na(hi), hi := FALSE]; H0[is.na(lo), lo := FALSE]
  H1 <- X[label == "mmlb", .(dgp, scen, N, ntime, strategy, upow = 100 * k / nd)]
  P <- merge(H0[, .(dgp, scen, N, ntime, strategy, n, k, nd, rate, hi, lo, pu)], H1, by = c(DPK, "strategy"))
  P[, grp := fifelse(dgp == "Copula-AR1", "copula", fifelse(scen %in% c("M6", "M7"), "m67", "main"))]; P[, small := ntime <= 3 & N <= 25]
  T <- P[, .(n_dp = .N, hi_main = sum(hi & grp == "main"), hi_m67 = sum(hi & grp == "m67"), hi_copula = sum(hi & grp == "copula"),
             hi_all = sum(hi), lo_all = sum(lo), lo_N10 = sum(lo & N == 10), n_N10 = sum(N == 10),
             pu = mean(pu), upow = mean(upow, na.rm = TRUE), pu_small = mean(pu[small]), n_small = sum(small),
             n_main = sum(grp == "main"), n_m67 = sum(grp == "m67"), n_copula = sum(grp == "copula")), by = strategy]
  T[, has_interval := !(strategy %in% NOWALD)]; T[, eligible := has_interval & hi_all == 0L]
  T[, chosen := eligible & pu == max(pu[eligible])]; setorder(T, -has_interval, -eligible, -pu)
  list(P = P, T = T)
}
## Tables 3 and 4 (R/s5/s5_engine.R, R/s6/s6_extras.R): H0, own set, main grid
t34 <- function(v) {
  X <- counts(v, "T34")[label == "manual"]
  bN <- X[, .(n = sum(n), k = sum(k)), by = .(dgp, strategy, N)][, rate := 100 * k / n]
  S <- bN[, .(wc = max(abs(rate - 5), na.rm = TRUE), mad = mean(abs(rate - 5), na.rm = TRUE)), by = .(dgp, strategy)]
  P <- X[, .(rate = 100 * sum(k) / sum(n)), by = .(dgp, strategy)]
  T3 <- merge(P, S, by = c("dgp", "strategy"))
  X[, `:=`(pu = 100 * n / nd, grp = fifelse(scen %in% c("M1", "M2"), "RI", "RS"))]
  T4 <- X[, .(pu = mean(pu)), by = .(grp, strategy, ntime)]
  list(T3 = T3, T4 = T4)
}
R <- lapply(names(VERS), function(v) c(t5(v), t34(v))); names(R) <- names(VERS)
## ---- (1) baseline against the current files --------------------------------------------------------------------------
T5f <- fread(file.path(TB, "V5_S7f_table4.csv")); B5 <- R$base$T
m <- merge(B5, T5f, by = "strategy", suffixes = c("", ".f"))
numc <- c("n_dp", "hi_main", "hi_m67", "hi_copula", "hi_all", "lo_all", "lo_N10", "n_N10", "pu", "upow", "pu_small", "n_small", "n_main", "n_m67", "n_copula")
d5 <- max(sapply(numc, function(cc) max(abs(m[[cc]] - m[[paste0(cc, ".f")]]))))
l5 <- all(sapply(c("has_interval", "eligible", "chosen"), function(cc) identical(as.logical(m[[cc]]), as.logical(m[[paste0(cc, ".f")]]))))
gate("baseline Table 5: 19 strategies, every numeric column within 1e-9, flags identical", nrow(m) == 19L && d5 < 1e-9 && l5, sprintf("rows %d, max |diff| %.3g", nrow(m), d5))
RD <- fread(file.path(TB, "V5_S7f_rule_dp.csv")); mp <- merge(R$base$P, RD, by = c(DPK, "strategy"), suffixes = c("", ".f"))
gate("baseline Table 5 per point: 45,144 rows, pu/upow within 1e-9, hi/lo identical", nrow(mp) == 45144L && mp[, max(abs(pu - pu.f))] < 1e-9 && mp[, max(abs(upow - upow.f))] < 1e-9 &&
     mp[, all(hi == hi.f & lo == lo.f)], sprintf("rows %d, max |pu| %.3g, max |upow| %.3g, hi diff %d, lo diff %d", nrow(mp), mp[, max(abs(pu - pu.f))], mp[, max(abs(upow - upow.f))], mp[, sum(hi != hi.f)], mp[, sum(lo != lo.f)]))
gate("baseline rule: same chosen strategy", identical(B5[chosen == TRUE, strategy], T5f[chosen == TRUE, strategy]), paste(B5[chosen == TRUE, strategy]))
VAL <- fread(file.path(TB, "variant_cv/V4_S6_values.csv"))[variant == "main"]
b3 <- melt(R$base$T3, id.vars = c("dgp", "strategy"), variable.name = "what")[, name := sprintf("%s|%s|%s", what, SHORT[dgp], strategy)]
m3 <- merge(b3, VAL[, .(name, est)], by = "name")
gate("baseline Table 3: 72 values (rate, wc, mad x 3 mechanisms x 8 methods) within 1e-9", nrow(m3) == 72L && m3[, max(abs(value - est))] < 1e-9, sprintf("rows %d, max |diff| %.3g", nrow(m3), m3[, max(abs(value - est))]))
b4 <- R$base$T4[, name := sprintf("puT|%s|%s|%d", grp, strategy, ntime)]; m4 <- merge(b4, VAL[, .(name, est)], by = "name")
gate("baseline Table 4 (puT, all methods and T): values within 1e-9", nrow(m4) == 96L && m4[, max(abs(pu - est))] < 1e-9, sprintf("rows %d, max |diff| %.3g", nrow(m4), m4[, max(abs(pu - est))]))
fwrite(rbindlist(G), file.path(TB, "V8_S11_degen_baseline.csv"))
if (!HAVE_CL) { cat("classes not available: baseline only\n"); quit(save = "no") }
saveRDS(list(R = R, AG = AG, ND = ND), file.path(D, "tables_state.rds"))

## ---- (2) usable results and rejections by class -------------------------------------------------------------------------
AG[, `:=`(region = fifelse(ntime <= 3 & N <= 25, "T2-3 & N<=25", "other"), arm = fifelse(label == "mmlb", "H1", "H0"), anyc = S | V1 | E1 | B0 | B1)]
CLS <- c("S", "Si", "V1", "E1", "B0", "B1", "L", "const", "anyc")
ub <- function(by) rbindlist(lapply(CLS, function(cc) AG[, .(usable_total = sum(n), rej_total = sum(k), usable_class = sum(n[get(cc)]), rej_class = sum(k[get(cc)])), by = by][, class := cc]))
U <- rbind(ub(c("family", "strategy", "dgp", "region", "arm")), ub(c("family", "strategy", "arm"))[, `:=`(dgp = "all", region = "all")],
           ub(c("family", "strategy", "region", "arm"))[, dgp := "all"], ub(c("family", "strategy"))[, `:=`(dgp = "all", region = "all", arm = "both")])
U[class == "anyc", class := "S|V1|E1|B0|B1"]; U[, pct_usable_in_class := 100 * usable_class / usable_total]
setnames(U, "dgp", "mechanism"); setcolorder(U, c("family", "strategy", "mechanism", "region", "arm", "class")); setorder(U, family, strategy, mechanism, region, arm, class)
fwrite(U, file.path(TB, "V8_S11_degen_usable.csv"))
## ---- (3) sensitivity: changes against the baseline, value and printed value -------------------------------------------
p1 <- function(x) sprintf("%.1f", x); p0 <- function(x) sprintf("%d", as.integer(round(x)))
SEN <- list()
B5 <- R$base$T; B3 <- R$base$T3; B4 <- R$base$T4
T4PRINT <- c("glmer", "glmer_lrt", "mmm", "gee_fg", "gee"); T3ORDER <- c("gee_fg", "gee_md", "mmm", "glmer", "shrink", "gee", "glmer_lrt", "gee_kc")
for (v in setdiff(names(VERS), "base")) {
  N5 <- R[[v]]$T
  for (cc in c("hi_main", "hi_m67", "hi_copula", "hi_all", "lo_all", "lo_N10", "pu", "upow", "pu_small")) {
    m <- merge(B5[, .(strategy, b = get(cc))], N5[, .(strategy, x = get(cc))], by = "strategy")
    fm <- if (cc %in% c("pu", "upow", "pu_small")) p1 else p0
    SEN[[length(SEN) + 1L]] <- m[, .(version = v, table = "Table 5", row = strategy, column = cc, base = b, new = x, diff = x - b,
                                       printed_base = fm(b), printed_new = fm(x), in_paper = TRUE)]
  }
  m <- merge(B5[, .(strategy, b = as.numeric(chosen))], N5[, .(strategy, x = as.numeric(chosen))], by = "strategy")
  SEN[[length(SEN) + 1L]] <- m[, .(version = v, table = "Table 5", row = strategy, column = "chosen", base = b, new = x, diff = x - b, printed_base = p0(b), printed_new = p0(x), in_paper = TRUE)]
  m <- merge(B5[, .(strategy, b = as.numeric(eligible))], N5[, .(strategy, x = as.numeric(eligible))], by = "strategy")
  SEN[[length(SEN) + 1L]] <- m[, .(version = v, table = "Table 5", row = strategy, column = "eligible", base = b, new = x, diff = x - b, printed_base = p0(b), printed_new = p0(x), in_paper = FALSE)]
  for (cc in c("rate", "wc", "mad")) {
    m <- merge(B3[, .(dgp, strategy, b = get(cc))], R[[v]]$T3[, .(dgp, strategy, x = get(cc))], by = c("dgp", "strategy"))
    SEN[[length(SEN) + 1L]] <- m[, .(version = v, table = "Table 3", row = strategy, column = paste(cc, SHORT[dgp]), base = b, new = x, diff = x - b,
                                       printed_base = p1(b), printed_new = p1(x), in_paper = TRUE)]
  }
  m <- merge(B4[, .(grp, strategy, ntime, b = pu)], R[[v]]$T4[, .(grp, strategy, ntime, x = pu)], by = c("grp", "strategy", "ntime"))
  SEN[[length(SEN) + 1L]] <- m[, .(version = v, table = "Table 4", row = strategy, column = sprintf("%s T=%d", grp, ntime), base = b, new = x, diff = x - b,
                                     printed_base = p1(b), printed_new = p1(x), in_paper = strategy %in% T4PRINT & ntime != 8L)]
}
SEN <- rbindlist(SEN); SEN[, `:=`(version_label = VERS[version], printed_changes = printed_base != printed_new)]
fwrite(SEN, file.path(TB, "V8_S11_degen_sensitivity.csv"))
T5V <- rbindlist(lapply(names(VERS), function(v) cbind(version = v, R[[v]]$T))); fwrite(T5V, file.path(TB, "V8_S11_degen_table5_versions.csv"))
fwrite(rbindlist(lapply(names(VERS), function(v) cbind(version = v, R[[v]]$P))), file.path(D, "table5_points_versions.csv.gz"))
## ---- (4) predictions -----------------------------------------------------------------------------------------------------
PDL <- list(); pd <- function(code, quantity, expected, relevance, observed, verdict, explanation) PDL[[length(PDL) + 1L]] <<- data.table(code, quantity, expected, relevance, observed, verdict, explanation)
CLr <- readRDS(file.path(D10, "classes.rds"))[elig == TRUE, .(N, ntime, S)]
s_all <- 100 * CLr[, mean(S)]; s_sm <- 100 * CLr[ntime <= 3 & N <= 25, mean(S)]; s_45 <- 100 * CLr[N >= 45, mean(S)]; rm(CLr)
pd("PD1", "share of separated datasets (S) among eligible datasets", "2-10% overall; 20-60% at T 2-3 & N<=25; <1% at N>=45", "descriptive; > 10% overall -> Section 2.5 must say how separated datasets are counted",
   sprintf("overall %.2f%%; T 2-3 & N<=25 %.2f%%; N>=45 %.3f%%", s_all, s_sm, s_45),
   if (s_all <= 10) "confermata" else "smentita",
   sprintf("expected ranges met: overall %s, small %s, N>=45 %s", s_all >= 2 & s_all <= 10, s_sm >= 20 & s_sm <= 60, s_45 < 1))
sh <- function(s, arms = c("manual", "mmlb")) AG[family == "T5" & strategy == s & label %in% arms, 100 * sum(n[Si]) / sum(n)]
hp <- R$base$P[strategy == "gee_ar1_fg" & hi == TRUE, .(dgp, scen, N, ntime)]
hpn <- AG[family == "T5" & strategy == "gee_ar1_fg" & label == "manual" & Si][hp, on = DPK, nomatch = NULL][, sum(n)]
o2 <- c(mmm_cc = sh("mmm_cc"), gee_ar1_fg = sh("gee_ar1_fg"))
pd("PD2", "usable results on Si datasets, share of usable results (both arms, 2,376 points)", "below 1% for both; at most 10 at the 2 points where GEE AR(1) FG is clearly above 7.5%", "0-1% overall",
   sprintf("MMM with check %.3f%% (H0 %.3f%%); GEE AR(1) FG %.3f%% (H0 %.3f%%); Si usable at its %d flagged points (H0): %d", o2[1], sh("mmm_cc", "manual"), o2[2], sh("gee_ar1_fg", "manual"), nrow(hp), hpn),
   if (all(o2 <= 1)) "confermata" else "smentita", sprintf("flagged points: %s", paste(hp[, sprintf("%s %s N%d T%d", dgp, scen, N, ntime)], collapse = "; ")))
for (v in c("a", "b")) { T <- R[[v]]$T
  c1 <- identical(T[chosen == TRUE, strategy], "mmm_cc"); c2 <- T[strategy == "mmm_cc", hi_all] == 0L; g <- T[strategy == "gee_ar1_fg", hi_all]; c3 <- g >= 1 & g <= 3
  pd(paste0("PD3", v), sprintf("rule of Table 5 after version %s (%s)", v, VERS[v]), "chooses MMM with check; MMM 0 points clearly above 7.5%; GEE AR(1) FG 1-3", "categorical: confirmed if all three true",
     sprintf("chosen %s; MMM hi %d; GEE AR(1) FG hi %d", paste(T[chosen == TRUE, strategy], collapse = ","), T[strategy == "mmm_cc", hi_all], g),
     if (c1 && c2 && c3) "confermata" else "smentita", "") }
d3 <- SEN[version == "a" & table == "Table 3"]; d3[, what := sub(" .*$", "", column)]
r_ag <- d3[what %in% c("rate", "mad"), range(diff)]; r_wc <- d3[what == "wc", range(diff)]
pd("PD4", "Table 3 after version a: change of rate and MAD (each method x mechanism); change of worst case", "-0.5..+0.5 points; worst case -1..+1.5", "-1..+1 for rate and MAD; -2.5..+2.5 for the worst case",
   sprintf("rate/MAD change %.3f to %.3f; worst-case change %.3f to %.3f", r_ag[1], r_ag[2], r_wc[1], r_wc[2]),
   if (r_ag[1] >= -1 && r_ag[2] <= 1 && r_wc[1] >= -2.5 && r_wc[2] <= 2.5) "confermata" else "smentita",
   sprintf("within the expected ranges: %s", r_ag[1] >= -0.5 && r_ag[2] <= 0.5 && r_wc[1] >= -1 && r_wc[2] <= 1.5))
d4 <- SEN[version == "a" & table == "Table 4"]; d4[, Tn := as.integer(sub("^.*T=", "", column))]
r4 <- d4[, range(diff)]; r4p <- d4[in_paper == TRUE, range(diff)]
pd("PD5", "Table 4 after version a: change of the usable share per cell", "-3..0 at T = 2; -0.5..0 at T >= 4", "-2.5..+2.5 per cell",
   sprintf("all 96 cells %.3f to %.3f (printed cells %.3f to %.3f); T=2 %.3f to %.3f; T>=4 %.3f to %.3f", r4[1], r4[2], r4p[1], r4p[2], d4[Tn == 2, min(diff)], d4[Tn == 2, max(diff)], d4[Tn >= 4, min(diff)], d4[Tn >= 4, max(diff)]),
   if (r4[1] >= -2.5 && r4[2] <= 2.5) "confermata" else "smentita",
   sprintf("within the expected ranges: %s", d4[Tn == 2, min(diff) >= -3 & max(diff) <= 0] && d4[Tn >= 4, min(diff) >= -0.5 & max(diff) <= 0]))
Td <- R$d$T; d6 <- SEN[version == "d" & table == "Table 3"]; d6[, what := sub(" .*$", "", column)]; r6 <- d6[what %in% c("rate", "mad"), range(diff)]
c6 <- identical(Td[chosen == TRUE, strategy], B5[chosen == TRUE, strategy])
pd("PD6", "version d (as a, plus V1, E1, B0, B1): rule choice; Table 3 rate and MAD change", "choice unchanged; Table 3 within +/-0.5", "categorical for the rule; -1..+1 for Table 3",
   sprintf("chosen %s; rate/MAD change %.3f to %.3f; worst-case change %.3f to %.3f", paste(Td[chosen == TRUE, strategy], collapse = ","), r6[1], r6[2], d6[what == "wc", min(diff)], d6[what == "wc", max(diff)]),
   if (c6 && r6[1] >= -1 && r6[2] <= 1) "confermata" else "smentita", sprintf("Table 3 within +/-0.5: %s", r6[1] >= -0.5 && r6[2] <= 0.5))
PD <- rbindlist(PDL); fwrite(PD, file.path(TB, "V8_S11_degen_predictions.csv")); print(PD)
print(SEN[printed_changes == TRUE & in_paper == TRUE, .N, by = .(version, table)])
cat("done\n")
