# =============================================================================
# s6_extras.R -- quantities for the v4 tables and Figure 1 that the S5 engine does not
# return, as a function of the cell-level counts of ONE bootstrap replicate (same input
# as engine()). Returns a named numeric vector.
#   vaN|G|m|N, puN|G|m|N       unconditional validity / usable share by cluster count (Fig 1c)
#   vaNT|G|m|N|Tb              the same by series-length band (T23: T = 2-3, T4p: T >= 4)
#   argN|G|N, argNT|G|N|Tb     index in PRIM8 of the method maximising unconditional validity
#   argNT|all|N|Tb             the same pooled over mechanisms; gapNT|.. lead over the second
#   madT1|m|T                  mean absolute departure by exact T over all mechanisms (Fig 1b)
#   wc96|G|m, mad96|G|m        worst case and mean departure over the N x T combinations (Fig 1d)
#   puT|grp|m|T                usable share by T, random intercept only (RI) or random slopes (RS) (Table 3.8)
#   ncell20|m                  design points at N <= 20 with at least one usable replicate (Table 3.9)
#   tsh25|m                    normal minus t rejection rate at N <= 25, pooled (Table 3.12)
#   irate|G|m                  rejection rate on the intersection set (Table 3.3)
#   tshift|G|m                 t minus normal rejection rate, own set (Table 3.12)
# =============================================================================
extras <- function(cl) {
  nv <- function(nm, v) setNames(as.numeric(v), nm)
  S <- SHORT
  H0 <- cl[arm == "H0" & set == "own" & method %in% PRIM8]
  H0[, `:=`(pu = 100 * n / nd, va = 100 * (n - k) / nd, dev = fifelse(n > 0, abs(100 * k / n - 5), NA_real_), Tb = TCODE(ntime))]
  out <- list()
  x <- H0[, .(va = mean(va), pu = mean(pu)), by = .(dgp, method, N)]
  out$vaN <- nv(sprintf("vaN|%s|%s|%d", S[x$dgp], x$method, x$N), x$va)
  out$puN <- nv(sprintf("puN|%s|%s|%d", S[x$dgp], x$method, x$N), x$pu)
  am <- x[, .(code = match(method[which.max(va)], PRIM8)), by = .(dgp, N)]
  out$argN <- nv(sprintf("argN|%s|%d", S[am$dgp], am$N), am$code)
  y <- H0[, .(va = mean(va)), by = .(dgp, method, N, Tb)]
  out$vaNT <- nv(sprintf("vaNT|%s|%s|%d|%s", S[y$dgp], y$method, y$N, y$Tb), y$va)
  lead <- function(d) d[, { o <- order(-va); .(code = match(method[o[1]], PRIM8), gap = va[o[1]] - va[o[2]]) }, by = .(g, N, Tb)]
  a1 <- lead(y[, .(g = S[dgp], method, N, Tb, va)])
  a2 <- lead(H0[, .(va = mean(va)), by = .(method, N, Tb)][, g := "all"])
  a <- rbind(a1, a2)
  out$argNT <- nv(sprintf("argNT|%s|%d|%s", a$g, a$N, a$Tb), a$code)
  out$gapNT <- nv(sprintf("gapNT|%s|%d|%s", a$g, a$N, a$Tb), a$gap)
  m <- H0[is.finite(dev), .(mad = mean(dev)), by = .(method, ntime)]
  out$madT1 <- nv(sprintf("madT1|%s|%d", m$method, m$ntime), m$mad)
  r <- H0[, .(rate = 100 * sum(k) / sum(n)), by = .(dgp, method, N, ntime)]
  s <- r[, .(wc = max(abs(rate - 5), na.rm = TRUE), mad = mean(abs(rate - 5), na.rm = TRUE)), by = .(dgp, method)]
  out$wc96 <- nv(sprintf("wc96|%s|%s", S[s$dgp], s$method), s$wc)
  out$mad96 <- nv(sprintf("mad96|%s|%s", S[s$dgp], s$method), s$mad)
  H0[, grp := fifelse(scen %in% c("M1", "M2"), "RI", "RS")]
  p <- H0[, .(pu = mean(pu)), by = .(grp, method, ntime)]
  out$puT <- nv(sprintf("puT|%s|%s|%d", p$grp, p$method, p$ntime), p$pu)
  c20 <- H0[N <= 20, .(nc = sum(n > 0)), by = method]
  out$ncell20 <- nv(sprintf("ncell20|%s", c20$method), c20$nc)
  t25 <- H0[N <= 25 & method %in% WALD7, .(d = 100 * (sum(k) - sum(kt)) / sum(n)), by = method]
  out$tsh25 <- nv(sprintf("tsh25|%s", t25$method), t25$d)
  it <- cl[arm == "H0" & set == "intersection" & method %in% PRIM8, .(rate = 100 * sum(k) / sum(n)), by = .(dgp, method)]
  out$irate <- nv(sprintf("irate|%s|%s", S[it$dgp], it$method), it$rate)
  ts <- H0[method %in% WALD7, .(d = 100 * (sum(kt) - sum(k)) / sum(n)), by = .(dgp, method)]
  out$tshift <- nv(sprintf("tshift|%s|%s", S[ts$dgp], ts$method), ts$d)
  unlist(unname(out))
}
