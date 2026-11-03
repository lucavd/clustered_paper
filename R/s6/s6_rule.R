# =============================================================================
# s6_rule.R -- S6: evidence for the branches of the decision rule (Figure 1a).
# Per branch (cluster count N x series band: T = 2-3 or T >= 4), pooled over mechanisms and
# scenarios (main grid) or over mechanisms (M1c), from tables/V4_S6_nestedloop.csv:
#   leader      method with the highest unconditional validity (criterion of v3 sec. 4.2)
#   lead        its margin over the second (points)
#   rate        null rejection rate of the leader in the branch (%)
#   floor_best  method with the highest unconditional validity among those whose branch
#               rejection rate is inside Bradley's 2.5-7.5% (option with a calibration floor)
#   cost        unconditional validity given up by the floor (points)
#   n_mech_lead mechanisms (of 3) in which the pooled leader is also the leader
# Outputs: tables/V4_S6_rule_branches.csv, tables/V4_S6_rule_rates.csv (rate and validity of every method per branch)
# =============================================================================
source(path.expand("~/MMM/R/s6/s6_common.R"))
NL <- fread(file.path(TABDIR, "V4_S6_nestedloop.csv"))
stopifnot(nrow(NL) == 1512L * 8L, !anyNA(NL$n_den), all(NL[n_obs > 0, is.finite(rej)]))   # every cell present; empty cells count as failures
NL[, `:=`(k = fifelse(n_obs > 0, rej / 100 * n_obs, 0), Tb = TCODE(ntime), grid = fifelse(scen == "M1c", "M1c", "main"))]
agg <- function(by) NL[, .(k = sum(k), n = sum(n_obs), nd = sum(n_den)), by = by][, `:=`(rate = 100 * k / n, va = 100 * (n - k) / nd)]
P <- agg(c("grid", "N", "Tb", "method")); PM <- agg(c("grid", "dgp", "N", "Tb", "method"))
lm <- PM[, .(lead_m = method[which.max(va)]), by = .(grid, dgp, N, Tb)]
R <- P[, { o <- order(-va); ld <- method[o[1]]; ok <- which(rate >= 2.5 & rate <= 7.5); fb <- if (length(ok)) method[ok[which.max(va[ok])]] else NA_character_
          .(leader = ld, lead = va[o[1]] - va[o[2]], rate = rate[o[1]], va_leader = va[o[1]],
            floor_best = fb, rate_floor = if (is.na(fb)) NA_real_ else rate[method == fb], cost = if (is.na(fb)) NA_real_ else va[o[1]] - va[method == fb]) },
       by = .(grid, N, Tb)]
R[lm[, .(n_mech_lead = .N), by = .(grid, N, Tb, lead_m)], on = .(grid, N, Tb, leader = lead_m), n_mech_lead := i.n_mech_lead]
R[is.na(n_mech_lead), n_mech_lead := 0L]
setorder(R, grid, Tb, N); fwrite(R, file.path(TABDIR, "V4_S6_rule_branches.csv"))
setorder(P, grid, Tb, N, method); fwrite(P[, .(grid, N, Tb, method, n, nd, rate, va)], file.path(TABDIR, "V4_S6_rule_rates.csv"))   # every method, every branch
print(R, digits = 3)
