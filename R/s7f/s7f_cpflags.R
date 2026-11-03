#!/usr/bin/env Rscript
# =============================================================================
# s7f_cpflags.R -- S7f (PI decision of 2026-10-27, from the second adversarial review): the flags "clearly above 7.5%"
# and "clearly below 2.5%" use the exact two-sided 95% Clopper-Pearson interval of the rejection rate instead of
# rate -/+ 2 Wald MCSE (with 0 rejections the Wald MCSE is 0, so a cell with 23 usable replicates was "clearly below").
# For every per-point table that carries flags: (1) k = n * rate / 100 is recovered and checked integer; (2) the Wald
# rule is reapplied and must reproduce the stored flags exactly (gate: same rows, same assess rule); (3) the flags are
# replaced by the exact ones; every other column is copied unchanged. Aggregates are rebuilt with the code of the step
# that first produced them. Inputs read-only; outputs tables/V5_S7f_cp_*.csv. Run from ~/MMM before s7f_table4.R.
# =============================================================================
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(1L); setwd("/home/user/MMM")
cp_lo <- function(k, n) fifelse(k <= 0, 0, qbeta(0.025, pmax(k, 1e-12), n - k + 1))
cp_up <- function(k, n) fifelse(k >= n, 1, qbeta(0.975, k + 1, pmax(n - k, 1e-12)))
G <- list()
gate <- function(name, ok, detail = "") { G[[length(G) + 1]] <<- data.table(gate = name, ok = isTRUE(ok), detail = detail); if (!isTRUE(ok)) stop("gate failed: ", name, " ", detail) }
reflag <- function(X, name, mc_col = "mc") {
  X <- copy(X)
  X[, k_ := n * rate / 100]
  gate(paste0(name, ": k integer"), X[assess == TRUE & !is.na(rate), max(abs(k_ - round(k_)))] < 1e-6)
  X[, k_ := round(k_)]
  w_hi <- X[, assess & rate - 2 * get(mc_col) > 7.5]; w_lo <- X[, assess & rate + 2 * get(mc_col) < 2.5]
  w_hi[is.na(w_hi)] <- FALSE; w_lo[is.na(w_lo)] <- FALSE
  gate(paste0(name, ": Wald reproduces stored hi"), identical(as.logical(w_hi), as.logical(X$hi)), sprintf("%d rows", nrow(X)))
  gate(paste0(name, ": Wald reproduces stored lo"), identical(as.logical(w_lo), as.logical(X$lo)))
  h <- X[, assess & 100 * cp_lo(k_, n) > 7.5]; l <- X[, assess & 100 * cp_up(k_, n) < 2.5]
  h[is.na(h)] <- FALSE; l[is.na(l)] <- FALSE
  CH <- data.table(table = name, hi_wald = sum(X$hi), hi_cp = sum(h), lo_wald = sum(X$lo), lo_cp = sum(l), hi_changed = sum(h != X$hi), lo_changed = sum(l != X$lo))
  X[, `:=`(hi = h, lo = l)][, k_ := NULL]
  attr(X, "change") <- CH; X
}
CHG <- list(); keep <- function(X) { CHG[[length(CHG) + 1]] <<- attr(X, "change"); X }
## ---- per-point tables used by Table 4, the N - 2 rule, Figure 1, Table S28 and Table S31 -------------------------
A  <- keep(reflag(fread("tables/V4_S7d_dp_rates.csv"), "S7d dp_rates"));  fwrite(A, "tables/V5_S7f_cp_S7d_dp_rates.csv")
B  <- keep(reflag(fread("tables/V4_S7e_dp_rates.csv"), "S7e dp_rates"));  fwrite(B, "tables/V5_S7f_cp_S7e_dp_rates.csv")
C  <- keep(reflag(fread("tables/V4_S7b_rule_dp.csv"), "S7b rule_dp"));    fwrite(C, "tables/V5_S7f_cp_S7b_rule_dp.csv")
FL <- keep(reflag(fread("tables/V4_S7d2_flip_rules_dp.csv")[, assess := n >= 20], "S7d2 flip_rules_dp", mc_col = "mcse")); FL[, assess := NULL]
fwrite(FL, "tables/V5_S7f_cp_S7d2_flip_rules_dp.csv")
## ---- Table S28 (fourth mechanism): rebuilt as in R/s7e/s7e_summary.R ----------------------------------------------
NEW <- B[ref == "t"]
RS <- A[ref == "t" & scen %in% sprintf("M%d", 1:5)]; RS[strategy == "study_mmm", strategy := "mmm_cc"]
RB <- C[ref == "t" & scen %in% sprintf("M%d", 1:5) & strategy %in% c("glmer|t", "shrink|t")]; RB[, strategy := sub("\\|t$", "", strategy)]
S28o <- fread("tables/V4_S7e_summary.csv")
SMc <- NEW[, .(clearly_above = sum(hi), clearly_below = sum(lo)), by = .(strategy, dgp)]
S28 <- merge(S28o[, !c("clearly_above", "clearly_below")], SMc, by = c("strategy", "dgp"), all.x = TRUE)
gate("S28: every row matched", S28[dgp == "Copula-AR1", !anyNA(clearly_above)])
setcolorder(S28, names(S28o)); fwrite(S28, "tables/V5_S7f_cp_S7e_summary.csv")

## ---- Table S28, sign-flip row: the version with the non-computability rules of Table 4 (second review) -------------
f0 <- FL[version == "with D-57 rules" & label == "manual" & dgp == "Copula-AR1"]
f1 <- FL[version == "with D-57 rules" & label == "mmlb" & dgp == "Copula-AR1"]
gate("S28 flip: 480 points per arm", nrow(f0) == 480L && nrow(f1) == 480L)
fr <- merge(f0[, .(scen, N, ntime, rate, hi, lo, pu)], f1[, .(scen, N, ntime, upow = pu * rate / 100)], by = c("scen", "N", "ntime"))
S28[strategy == "flip" & dgp == "Copula-AR1", `:=`(mean_type1 = fr[, mean(rate, na.rm = TRUE)], clearly_above = fr[, sum(hi)], clearly_below = fr[, sum(lo)],
                                                    mean_usable = fr[, mean(pu)], mean_power = fr[, mean(upow)])]
fwrite(S28, "tables/V5_S7f_cp_S7e_summary.csv"); print(S28[strategy == "flip"])
## ---- scenarios M6-M7 (nested-loop tables): rebuilt as in R/s6b/s6b_ns_agg.R ---------------------------------------
NS <- keep(reflag(fread("tables/V4_S6b_ns_rule_safety_dp.csv"), "S6b ns rule_safety_dp"))
sm <- function(D, by) D[, .(n_dp = .N, n_hi = sum(hi), n_lo = sum(lo), n_na = sum(!assess), max_rate = max(rate[assess]),
                            n_point_above = sum(assess & rate > 7.5), pu = mean(pu), upow = mean(upow)), by = by]
NSS <- rbind(sm(NS, c("strategy", "dgp", "scen")), sm(NS, c("strategy", "dgp", "scen", "band")), sm(NS, "strategy")[, `:=`(dgp = "all", scen = "all")], fill = TRUE)
NSo <- fread("tables/V4_S6b_ns_rule_safety.csv")
NSS[is.na(band), band := ""]   # fread reads the empty band of the stored table as ""
gate("S6b ns: non-flag columns identical", isTRUE(all.equal(NSo[, !c("n_hi", "n_lo")], NSS[, !c("n_hi", "n_lo")], check.attributes = FALSE)))
fwrite(NSS, "tables/V5_S7f_cp_S6b_ns_rule_safety.csv")
## ---- main grid and M1c, variant with KC from fewclust: rebuilt as in R/s6/s6_rule_safety.R -------------------------
RK <- keep(reflag(fread("tables/variant_cv/V4_S6_rule_safety_dp.csv"), "kcfc rule_safety_dp"))
smk <- function(by) RK[, .(n_dp = .N, n_hi = sum(hi), n_lo = sum(lo), n_na = sum(!assess), max_rate = max(rate[assess]),
                           n_point_above = sum(assess & rate > 7.5), pu = mean(pu), upow = mean(upow)), by = by]
SK <- rbind(smk(c("strategy", "N", "band")), smk(c("strategy", "N"))[, band := "all"], smk(c("strategy", "band"))[, N := NA_integer_],
            smk("strategy")[, `:=`(N = NA_integer_, band = "all")], use.names = TRUE)
setorder(SK, strategy, band, N)
SKo <- fread("tables/variant_cv/V4_S6_rule_safety.csv")
gate("kcfc rule_safety: non-flag columns identical", isTRUE(all.equal(SKo[, !c("n_hi", "n_lo")], SK[, !c("n_hi", "n_lo")], check.attributes = FALSE)))
fwrite(SK, "tables/V5_S7f_cp_kcfc_rule_safety.csv")
## ---- convergence-check comparison cells (Section 3.3): k and n stored per cell -----------------------------------
CC <- fread("tables/V4_S6cc_cells.csv")
h0 <- CC$arm == "H0"
for (s in c("c", "s")) {   # R/s6/s6cc_compare.R: H0 rows only, no minimum of usable replicates
  n_ <- CC[[paste0("n_", s)]]; k_ <- CC[[paste0("k_", s)]]; r_ <- 100 * k_ / n_; m_ <- 100 * sqrt((r_ / 100) * (1 - r_ / 100) / n_)
  st <- CC[[paste0("hi_", s)]]; w <- r_ - 2 * m_ > 7.5
  gate(paste0("S6cc ", s, ": Wald reproduces hi"), identical(fifelse(is.na(w[h0]), FALSE, w[h0]), fifelse(is.na(st[h0]), FALSE, as.logical(st[h0]))))
  hv <- 100 * cp_lo(k_, n_) > 7.5; lv <- 100 * cp_up(k_, n_) < 2.5
  hv[!h0] <- NA; lv[!h0] <- NA
  CC[[paste0("hi_", s)]] <- hv; CC[[paste0("lo_", s)]] <- lv
}
fwrite(CC, "tables/V5_S7f_cp_S6cc_cells.csv")
## ---- S7d rule (Kauermann-Carroll, fewclust against geesmv; Section S4) ----------------------------------------------
R7o <- fread("tables/V4_S7d_rule.csv")
ST <- C[ref == "t", .(strategy = paste0("paper_", sub("\\|t$", "", strategy)), hi, lo)]
R7n <- rbind(ST[, .(src = "study (saved rows)", n_hi = sum(hi), n_lo = sum(lo)), by = strategy],
             A[ref == "t", .(src = "S7d run", n_hi = sum(hi), n_lo = sum(lo)), by = strategy])
R7 <- merge(R7o[, !c("n_hi", "n_lo")], R7n, by = c("strategy", "src"), all.x = TRUE, sort = FALSE)
gate("S7d rule: rows matched", R7[, !anyNA(n_hi)] && nrow(R7) == nrow(R7o))
setcolorder(R7, names(R7o)); fwrite(R7, "tables/V5_S7f_cp_S7d_rule.csv")
## ---- null shift (Section 4.4: flags that change after subtracting the excess) --------------------------------------
NSH <- fread("tables/V4_S7b_null_shift.csv")
NSH[, k_ := n * rate / 100]; NSH[, ka_ := n * rate_adj / 100]
NSH[, `:=`(hi = assess & 100 * cp_lo(k_, n) > 7.5, lo = assess & 100 * cp_up(k_, n) < 2.5,
           hi_adj = assess & 100 * cp_lo(pmax(ka_, 0), n) > 7.5, lo_adj = assess & 100 * cp_up(pmax(ka_, 0), n) < 2.5)]
for (v in c("hi", "lo", "hi_adj", "lo_adj")) NSH[is.na(get(v)), (v) := FALSE]
NF <- data.table(k = "null_flips_cp", value = NSH[, sum(hi != hi_adj | lo != lo_adj)]); setnames(NF, "k", "key")
fwrite(NF, "tables/V5_S7f_cp_null_flips.csv")
## ---- report --------------------------------------------------------------------------------------------------------
CH <- rbindlist(CHG); fwrite(CH, "tables/V5_S7f_cp_changes.csv"); print(CH)
fwrite(rbindlist(G), "tables/V5_S7f_cp_gates.csv"); cat("gates:", length(G), "passed\n"); print(NF)
