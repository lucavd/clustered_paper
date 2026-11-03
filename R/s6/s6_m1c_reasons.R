# =============================================================================
# s6_m1c_reasons.R -- S6, follow-up of predictions P1-P2 (smentite): why the Wald routes stay
# far less usable than the LRT at N <= 20 in M1c. Reasons for a non-usable replicate, H0 arm,
# usability U0 of R/s4/s4_common.R, non-degenerate designs only.
# Output: tables/V4_S6_m1c_unusable_reasons.csv
# =============================================================================
source(path.expand("~/MMM/R/s6/s6_common.R"))
V <- readRDS(file.path(OUT_S4, "rows_v4.rds"))
print(names(V))
V <- V[!is.na(term) & N <= 20 & label == "manual" & method %in% PRIM8]
V[, grid := fifelse(scen == "M1c", "M1c", fifelse(scen %in% c("M1", "M2"), "M1-M2", "M3-M5"))]
V[, model := as.integer(sub("^M([1-5]).*$", "\\1", scen))]
DEG <- fread(file.path(TABDIR, "V4_S2_degenerate_designs.csv"))[, .(model, N, ntime, sim_id)]
V[, deg := FALSE]; V[DEG, on = c("model", "N", "ntime", "sim_id"), deg := TRUE]; V <- V[deg == FALSE]
V[, reason := fifelse(U0, "usable",
  fifelse(is.na(converged) | !converged, "error / not converged",
  fifelse(!is.finite(p), "p not finite",
  fifelse(method != LRT & (!is.finite(estimate) | !is.finite(se)), "estimate or SE not finite",
  fifelse(method != LRT & abs(estimate) >= 10, "|estimate| >= 10",
  fifelse(method != LRT & (se >= 10 | se <= 0), "SE >= 10 or <= 0", "other"))))))]
R <- V[, .N, by = .(grid, method, reason)][, pct := 100 * N / sum(N), by = .(grid, method)]
setorder(R, grid, method, -pct); fwrite(R, file.path(TABDIR, "V4_S6_m1c_unusable_reasons.csv"))
print(dcast(R, grid + method ~ reason, value.var = "pct", fill = 0), digits = 3)
