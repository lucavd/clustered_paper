# =============================================================================
# s5_engine.R -- S5: every comparative claim of the v3 manuscript as a function of the
# cell-level counts of ONE bootstrap replicate (inventory reports/morris/S5_claims_inventory.csv).
# cl: one row per cell with dgp, scen, N, ntime, arm, method, set, tok, n, k, kt, kc, e2, ksa, nd
#     (n usable, k rejections, kt t-reference rejections, kc covered, e2 sum of squared errors,
#      ksa size-adjusted rejections (H1 own), nd non-degenerate replicates of the design point).
# Returns a named numeric vector. Names: <claim>|<mechanism or scope>[|<extra>].
# tie = "share" (D-23) or "first" (v3, order of the rows of cl inside the design point).
# =============================================================================
`%||%` <- function(a, b) if (is.null(a)) b else a
gapmin <- function(v, target) { o <- v[names(v) != target]; o <- o[is.finite(o)]; if (!length(o) || !is.finite(v[target])) NA_real_ else min(o) - v[[target]] }
gapmax_over <- function(v, target) { o <- v[names(v) != target]; o <- o[is.finite(o)]; if (!length(o) || !is.finite(v[target])) NA_real_ else v[[target]] - max(o) }
maxother <- function(v, target) { o <- v[names(v) != target]; o <- o[is.finite(o)]; if (!length(o) || !is.finite(v[target])) NA_real_ else v[[target]] - max(o) }
vec <- function(dt, val, nm = "method") setNames(dt[[val]], dt[[nm]])
sp <- function(a, b) { ok <- is.finite(a) & is.finite(b); if (sum(ok) < 3) NA_real_ else suppressWarnings(cor(a[ok], b[ok], method = "spearman")) }
SHORT <- c("GenBinaryY" = "GB", "GLMM-native" = "NAT", "GLMM-AR1" = "AR1")

engine <- function(cl, tie = "share", grid = "main", PBIN = NULL, pairs_c05 = NULL, v3 = FALSE) {
  out <- list()
  put <- function(nm, v) { if (length(v) > 1) stop("non-scalar value for ", nm); out[[nm]] <<- if (length(v) == 0) NA_real_ else as.numeric(v) }
  cl <- if (grid == "main") cl[scen != "M1c"] else cl[scen == "M1c"]
  has <- function(col) col %in% names(cl) && any(is.finite(cl[[col]]))
  H0 <- cl[arm == "H0" & set == "own" & method %in% PRIM8]
  H0[, dev := fifelse(n > 0, abs(100 * k / n - 5), NA_real_)]
  bN <- H0[, .(n = sum(n), k = sum(k)), by = .(dgp, method, N)][, rate := 100 * k / n]
  S  <- bN[, .(wc = max(abs(rate - 5), na.rm = TRUE), mad = mean(abs(rate - 5), na.rm = TRUE),
               argN = N[which.max(abs(rate - 5))]), by = .(dgp, method)]
  P  <- H0[, .(rate = 100 * sum(k) / sum(n)), by = .(dgp, method)]
  for (g in DGPS) {
    s <- S[dgp == g]; p <- P[dgp == g]; G <- SHORT[[g]]
    if (!nrow(s)) next
    for (mm in s$method) { put(sprintf("rate|%s|%s", G, mm), p[method == mm, rate]); put(sprintf("wc|%s|%s", G, mm), s[method == mm, wc]); put(sprintf("mad|%s|%s", G, mm), s[method == mm, mad]) }
    put(paste0("C01|", G), gapmin(vec(s, "wc"), "mmm"))
    put(paste0("C02|", G), gapmin(vec(s, "mad"), "gee_fg"))
    put(paste0("C22|", G), s[method == "mmm", argN])
  }
  # ---- C03 paired, C03b intersection --------------------------------------------
  pr <- cl[arm == "H0" & set == "pair", .(rate = 100 * sum(k) / sum(n)), by = .(dgp, method)]
  it <- cl[arm == "H0" & set == "intersection" & method %in% PRIM8, .(rate = 100 * sum(k) / sum(n)), by = .(dgp, method)]
  for (g in DGPS) { G <- SHORT[[g]]
    if (nrow(pr[dgp == g]) == 2) put(paste0("C03|", G), pr[dgp == g & method == "mmm", rate] - pr[dgp == g & method == "gee_fg", rate])
    if (nrow(it[dgp == g])) put(paste0("C03b|", G), it[dgp == g & method == "mmm", rate] - it[dgp == g & method == "gee_fg", rate])
  }
  if (grid == "M1c") {
    H1 <- cl[arm == "H1" & set == "own" & method %in% WALD7 & tok == TRUE & n >= 2]
    cv <- H1[, .(dev = abs(100 * sum(kc) / sum(n) - 95)), by = .(dgp, method)]
    for (g in DGPS) { G <- SHORT[[g]]; v <- vec(cv[dgp == g], "dev"); put(paste0("C20a|", G), gapmin(v, "gee_fg")); put(paste0("C20b|", G), maxother(v, "gee_kc")) }
    return(unlist(out))
  }
  # ---- C04 wins, C05 winner vs runner-up -------------------------------------------
  D <- H0[is.finite(dev), .(dgp, scen, N, ntime, method, dev)]
  if (tie == "share") {
    W <- D[, { mn <- min(dev); w <- method[dev == mn]; .(method = w, w = 1 / length(w)) }, by = .(dgp, scen, N, ntime)]
  } else {
    W <- D[, .(method = method[which.min(dev)], w = 1), by = .(dgp, scen, N, ntime)]
  }
  wins <- W[, .(w = sum(w)), by = method]
  for (mm in PRIM8) put(paste0("wins|", mm), wins[method == mm, sum(w)])
  put("C04|all", wins[method == "mmm", sum(w)] - wins[method == "gee_fg", sum(w)])
  C05 <- NULL
  if (!is.null(pairs_c05)) {
    x <- merge(pairs_c05, D[, .(dgp, scen, N, ntime, method_w = method, dev_w = dev)], by = c("dgp", "scen", "N", "ntime", "method_w"), all.x = TRUE)
    x <- merge(x, D[, .(dgp, scen, N, ntime, method_r = method, dev_r = dev)], by = c("dgp", "scen", "N", "ntime", "method_r"), all.x = TRUE)
    setorder(x, pid); C05 <- x$dev_r - x$dev_w
  }
  # ---- C06, C07, C13d by series length -----------------------------------------------
  D[, Ts := fifelse(ntime == 2, "T2", fifelse(ntime == 3, "T3", "T4plus"))]
  mT <- D[, .(mad = mean(dev)), by = .(method, Ts)]
  for (ts in c("T2", "T3", "T4plus")) for (mm in PRIM8) put(sprintf("madT|%s|%s", ts, mm), mT[method == mm & Ts == ts, mad])
  put("C06a|all", gapmin(vec(mT[Ts == "T2"], "mad"), "mmm"))
  put("C06b|all", gapmin(vec(mT[Ts == "T3"], "mad"), "gee_fg"))
  put("C06c|all", gapmin(vec(mT[Ts == "T4plus"], "mad"), "gee_fg"))
  put("C06d|all", mT[method == "mmm" & Ts == "T2", mad] - mT[method == "mmm" & Ts == "T4plus", mad])
  put("C13d|all", mT[method == "glmer" & Ts == "T2", mad] - mT[method == "mmm" & Ts == "T2", mad])
  m2 <- D[ntime == 2, .(mad = mean(dev)), by = method]; m10 <- D[ntime == 10, .(mad = mean(dev)), by = method]
  mm2 <- merge(m2, m10, by = "method", suffixes = c("2", "10"))
  put("C07|all", sp(mm2$mad2, mm2$mad10))
  # ---- C08, C09 across mechanisms ------------------------------------------------------
  pairs <- list(c("GenBinaryY", "GLMM-native"), c("GenBinaryY", "GLMM-AR1"), c("GLMM-native", "GLMM-AR1"))
  SP <- merge(S, P, by = c("dgp", "method"))
  SP <- merge(SP, D[, .(madcell = mean(dev)), by = .(dgp, method)], by = c("dgp", "method"), all.x = TRUE)   # weighting of V-92
  for (crit in c("wc", "mad", "rate", "madcell")) {
    r <- sapply(pairs, function(pp) { a <- SP[dgp == pp[1]][order(method)]; b <- SP[dgp == pp[2]][order(method)]; mmx <- intersect(a$method, b$method); sp(a[method %in% mmx][[crit]], b[method %in% mmx][[crit]]) })
    put(paste0("C08|", crit), min(r))
  }
  a <- SP[dgp == "GLMM-AR1"][order(method)]; b <- SP[dgp == "GLMM-native"][order(method)]
  put("C09|all", mean(abs(a$rate - b$rate)))
  put("C09|rank_wc_identical", as.numeric(identical(rank(a$wc), rank(b$wc))))
  put("C09|rank_mad_identical", as.numeric(identical(rank(a$mad), rank(b$mad))))
  # ---- C10 multiplicity -------------------------------------------------------------------
  if (!is.null(PBIN)) {
    mc <- H0[n > 0, .(method, n, k)]
    mc[, pb := PBIN$p[PBIN$off[n] + k + 1L]]
    put("C10a|raw", 100 * mean(mc$pb < 0.05))
    put("C10a|BH", 100 * mean(p.adjust(mc$pb, "BH") < 0.05)); put("C10a|BY", 100 * mean(p.adjust(mc$pb, "BY") < 0.05))
    put("C10a|bonferroni", 100 * mean(p.adjust(mc$pb, "bonferroni") < 0.05)); put("C10a|ncells", nrow(mc))
    for (adj in c("BH", "BY")) {
      f <- mc[, .(pct = 100 * mean(p.adjust(pb, adj) < 0.05)), by = method]
      for (mm in f$method) put(sprintf("flag|%s|%s", adj, mm), f[method == mm, pct])
      put(paste0("C10b|", adj), gapmin(vec(f, "pct"), "gee_fg")); put(paste0("C10c|", adj), maxother(vec(f, "pct"), "gee_kc"))
    }
  }
  # ---- C11 intersection shift -------------------------------------------------------------
  sh <- merge(P, it, by = c("dgp", "method"), suffixes = c("_own", "_int"))[, shift := rate_int - rate_own]
  for (g in DGPS) { G <- SHORT[[g]]; s <- sh[dgp == g]
    for (mm in s$method) put(sprintf("shift|%s|%s", G, mm), s[method == mm, shift])
    put(paste0("C11a|", G), -s[method == "glmer_lrt", shift]); put(paste0("C11b|", G), s[method == "gee_kc", shift])
    put(paste0("C11c|", G), s[!method %in% c("glmer_lrt", "gee_kc"), max(abs(shift))]) }
  put("C11d|GB", it[dgp == "GenBinaryY" & method == "glmer_lrt", rate] - it[dgp == "GenBinaryY" & method == "glmer", rate])
  # ---- C12 uncorrected GEE -------------------------------------------------------------------
  put("C12a|GB", bN[dgp == "GenBinaryY" & method == "gee" & N == 10, rate] - 7.5)
  put("C12b|GB", bN[dgp == "GenBinaryY" & method == "gee" & N == 145, rate] - 5)
  for (g in DGPS) put(paste0("C12c|", SHORT[[g]]), bN[dgp == g & method == "gee", min(rate - 5)])
  # ---- usability (H0 arm, cells never fitted at zero) -----------------------------------------
  U <- cl[arm == "H0" & set == "own" & method %in% PRIM8, .(dgp, scen, N, ntime, method, n, k, nd)]
  U[, pu := 100 * n / nd]; U[, va := 100 * (n - k) / nd]
  put("C13a|all", U[method == "glmer" & scen %in% c("M3", "M4", "M5") & ntime == 2, mean(pu)])
  put("C13c|all", U[method == "glmer_lrt" & scen %in% c("M1", "M2") & ntime == 2, mean(pu)] - U[method == "glmer" & scen %in% c("M1", "M2") & ntime == 2, mean(pu)])
  u20 <- U[N <= 20, .(pu = mean(pu), va = mean(va)), by = method]
  for (mm in u20$method) { put(paste0("use20|", mm), u20[method == mm, pu]); put(paste0("val20|", mm), u20[method == mm, va]) }
  put("C14a|all", maxother(vec(u20, "pu"), "glmer_lrt")); put("C14b|all", maxother(vec(u20, "va"), "glmer_lrt"))
  put("C14c|all", gapmin(vec(u20, "pu"), "glmer"))
  put("C14d|all", max(abs(u20[method == "mmm", pu] - u20[method %in% c("gee_md", "gee_kc", "gee_fg"), pu])))
  # ---- C15 strata by estimability ------------------------------------------------------------
  PM <- U[method %in% MARG5, .(pu_marg = mean(n / nd)), by = .(dgp, scen, N, ntime)]
  PM[, stratum := cut(pu_marg, breaks = c(-Inf, 0.5, 0.8, 0.95, 0.999, Inf), labels = c("s1", "s2", "s3", "s4", "s5"), right = TRUE)]
  for (s in c("s1", "s2", "s3", "s4", "s5")) put(paste0("ndp|", s), PM[stratum == s, .N])
  DS <- merge(H0[, .(dgp, scen, N, ntime, method, dev)], PM, by = c("dgp", "scen", "N", "ntime"))
  ms <- DS[, .(mad = mean(dev, na.rm = TRUE)), by = .(method, stratum)]
  for (s in c("s1", "s2", "s3", "s4", "s5")) for (mm in PRIM8) put(sprintf("madS|%s|%s", s, mm), ms[method == mm & stratum == s, mad])
  for (s in c("s1", "s2", "s3")) put(paste0("C15a|", s), gapmin(vec(ms[stratum == s], "mad"), "mmm"))
  put("C15b|s5", gapmin(vec(ms[stratum == "s5"], "mad"), "gee_md")); put("C15c|s4", gapmin(vec(ms[stratum == "s4"], "mad"), "gee_fg"))
  GP <- dcast(DS[method %in% c("mmm", "gee_fg")], dgp + scen + N + ntime + stratum + pu_marg ~ method, value.var = "dev")
  GP <- GP[is.finite(mmm) & is.finite(gee_fg)][, gap := mmm - gee_fg]
  sg <- c(s1 = -1, s2 = -1, s3 = -1, s4 = 1, s5 = 1)
  for (s in names(sg)) { x <- GP[stratum == s, mean(gap)]; put(paste0("gapS|", s), x); put(paste0("C15d|", s), sg[[s]] * x) }
  put("C15e|all", sp(GP$pu_marg, GP$gap))
  # ---- C16 power (H1 own) -------------------------------------------------------------------------
  if (has("ksa")) {
    H1 <- cl[arm == "H1" & set == "own" & method %in% PRIM8]
    pw <- H1[, .(raw = 100 * sum(k) / sum(n), sa = 100 * sum(ksa[is.finite(ksa)]) / sum(n[is.finite(ksa)])), by = .(dgp, method)]
    for (g in DGPS) { G <- SHORT[[g]]; x <- pw[dgp == g]
      for (mm in x$method) { put(sprintf("raw|%s|%s", G, mm), x[method == mm, raw]); put(sprintf("sa|%s|%s", G, mm), x[method == mm, sa]) }
      put(paste0("C16a|", G), x[method == "glmer", sa] - x[method == "gee_fg", sa]); put(paste0("C16e|", G), max(x$sa) - min(x$sa)) }
    put("C16b|all", out[["C16a|GB"]] - out[["C16a|NAT"]])
    x <- pw[dgp == "GenBinaryY"]; put("C16c|GB", maxother(vec(x, "raw"), "gee_kc")); put("C16d|GB", -maxother(vec(x, "sa"), "gee_kc"))
  }
  # ---- C17, C24 t reference ---------------------------------------------------------------------
  if (has("kt")) {
    Ht <- H0[method %in% WALD7]
    tN <- Ht[, .(n = sum(n), k = sum(k), kt = sum(kt)), by = .(dgp, method, N)][, `:=`(rn = 100 * k / n, rt = 100 * kt / n)]
    tP <- Ht[, .(rn = 100 * sum(k) / sum(n), rt = 100 * sum(kt) / sum(n)), by = .(dgp, method)]
    tW <- tN[, .(wct = max(abs(rt - 5))), by = .(dgp, method)]
    for (g in DGPS) { G <- SHORT[[g]]; x <- tP[dgp == g]; w <- tW[dgp == g]
      for (mm in x$method) { put(sprintf("ratet|%s|%s", G, mm), x[method == mm, rt]); put(sprintf("wct|%s|%s", G, mm), w[method == mm, wct]) }
      put(paste0("C17a|", G), min(x$rn - x$rt))
      put(paste0("C17c|", G), gapmin(vec(w, "wct"), "gee_fg"))
      put(paste0("C17d|", G), w[method == "mmm", wct] - w[method == "gee_fg", wct])
      put(paste0("C17e|", G), if (g == "GenBinaryY") w[method == "gee_md", wct] - w[method == "mmm", wct] else w[method == "mmm", wct] - w[method == "gee_md", wct])
      put(paste0("C17f|", G), abs(x[method == "gee_fg", rt - rn]) - abs(x[method == "gee_fg", rn] - x[method == "gee_md", rn]))
      put(paste0("C24a|", G), tN[dgp == g & method == "gee_kc" & N == 10, rt])
      put(paste0("C24c|", G), x[method == "gee_md", rt]) }
    s25 <- Ht[N <= 25, .(d = 100 * sum(k) / sum(n) - 100 * sum(kt) / sum(n)), by = method]
    put("C17b|all", min(s25$d))
    put("C24b|all", Ht[method == "gee_kc" & N == 145, 100 * sum(kt) / sum(n)] - 7.5)
  }
  # ---- C19 RMSE, C20 coverage (H1 own, estimand-specific truth, n >= 2) --------------------------
  if (has("kc")) {
    H1 <- if (v3) cl[arm == "H1" & set == "own" & method %in% PRIM8 & nc > 0] else cl[arm == "H1" & set == "own" & method %in% WALD7 & tok == TRUE & n >= 2]
    nc <- if ("nc" %in% names(H1)) H1$nc else H1$n
    H1[, ncv := nc]
    cv <- H1[, .(cov = 100 * sum(kc) / sum(ncv)), by = .(dgp, method)][, dev := abs(cov - 95)]
    for (g in DGPS) { G <- SHORT[[g]]; x <- cv[dgp == g]
      for (mm in x$method) put(sprintf("cov|%s|%s", G, mm), x[method == mm, cov])
      put(paste0("C20a|", G), gapmin(vec(x, "dev"), "gee_fg")); put(paste0("C20b|", G), maxother(vec(x, "dev"), "gee_kc")) }
    if (has("e2")) {
      H1[, rmse := sqrt(e2 / n)]
      for (g in DGPS) { G <- SHORT[[g]]; x <- H1[dgp == g]
        ms <- unique(x$method)
        cm <- function(meths) { y <- x[method %in% meths]; ok <- y[, .(m = uniqueN(method)), by = .(scen, N, ntime)][m == length(meths)]
                                y[ok, on = c("scen", "N", "ntime")][, .(r = mean(rmse)), by = method] }
        r <- cm(ms); for (mm in r$method) put(sprintf("rmse|%s|%s", G, mm), r[method == mm, r])
        rr <- vec(r, "r"); put(paste0("C19a|", G), min(rr[names(rr) != "shrink"] / rr[["shrink"]]) - 1)
        if (g != "GenBinaryY") { r2 <- vec(cm(c("glmer", "mmm")), "r"); put(paste0("C19b|", G), r2[["glmer"]] / r2[["mmm"]] - 1) }
        r3 <- vec(cm(c("mmm", "gee_md", "gee_kc", "gee_fg")), "r"); put(paste0("C19c|", G), max(r3) / min(r3) - 1) }
    }
    # ---- C21 combined rank -------------------------------------------------------------------------
    lev <- P[, .(dgp, method, lev = abs(rate - 5))]
    usb <- U[, .(usable = mean(pu)), by = .(dgp, method)]
    crank <- function(lev, cvv, usb) {
      x <- merge(merge(lev, cvv[, .(dgp, method, cov = dev)], by = c("dgp", "method")), usb, by = c("dgp", "method"))
      x[, `:=`(r1 = frank(lev), r2 = frank(cov), r3 = frank(-usable)), by = dgp][, mr := (r1 + r2 + r3) / 3]
      x[, .(mr = mean(mr)), by = method] }
    cr <- crank(lev, cv, usb)
    for (mm in cr$method) put(paste0("crank|", mm), cr[method == mm, mr])
    put("C21a|all", gapmin(vec(cr, "mr"), "gee_fg"))
    H0b <- H0[N <= 25]; lev2 <- H0b[, .(lev = abs(100 * sum(k) / sum(n) - 5)), by = .(dgp, method)]
    cv2 <- H1[N <= 25, .(dev = abs(100 * sum(kc) / sum(ncv) - 95)), by = .(dgp, method)]
    cr2 <- crank(lev2, cv2, U[N <= 25, .(usable = mean(pu)), by = .(dgp, method)])
    for (mm in cr2$method) put(paste0("crank25|", mm), cr2[method == mm, mr])
  }
  # ---- C23 thresholds (unconditional validity and usability per N) -----------------------------
  VN <- U[, .(va = mean(va), pu = mean(pu)), by = .(dgp, method, N)]
  code <- setNames(seq_along(PRIM8), PRIM8)
  for (g in DGPS) { G <- SHORT[[g]]; x <- VN[dgp == g]; Ns <- sort(unique(x$N))
    am <- sapply(Ns, function(nn) { y <- x[N == nn]; y$method[which.max(y$va)] })
    sw <- if (am[1] != "glmer_lrt") 0 else Ns[max(which(cumprod(am == "glmer_lrt") == 1))]
    put(paste0("C23a|", G), sw)
    nxt <- if (sw == 0 || sw == max(Ns)) NA else am[which(Ns == sw) + 1]
    put(paste0("C23b|", G), if (is.na(nxt)) NA else code[[nxt]])
    put(paste0("C23c|", G), maxother(vec(x[N == 10], "va"), "glmer_lrt"))
    for (nn in c(25, 35)) put(sprintf("C23d|%s|N%d", G, nn), maxother(vec(x[N == nn], "pu"), "glmer_lrt"))
    w20 <- bN[dgp == g & N <= 20 & method %in% MARG5, .(wc = max(abs(rate - 5))), by = method]
    put(paste0("C23e|", G), gapmin(vec(w20, "wc"), "mmm")) }
  # ---- C25 large-N behaviour ------------------------------------------------------------------------
  for (g in DGPS) { G <- SHORT[[g]]
    if (g != "GenBinaryY") put(paste0("C25a|", G), bN[dgp == g & method == "gee_fg" & N >= 55, max(abs(rate - 5))])
    put(paste0("C25b|", G), bN[dgp == g & method == "mmm" & N >= 55, mean(rate - 5)])
    for (mm in c("gee_fg", "gee_md")) { x <- bN[dgp == g & method == mm & N >= 55]
      put(sprintf("C25c|%s|%s", G, mm), mean(abs(x$rate - 5) < 0.1))
      for (nn in x$N) put(sprintf("rateN|%s|%s|%d", G, mm, nn), x[N == nn, rate]) } }
  v <- unlist(out)
  if (!is.null(C05)) attr(v, "C05") <- C05
  v
}
