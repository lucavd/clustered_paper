# =============================================================================
# s6_tables.R -- S6: the v4 tables of the manuscript (main grid and M1c), every number read
# from tables/V4_S6_values.csv with its block-bootstrap MCSE, rounded to the MCSE.
# Rounding rule (pre-registered in docs/steps/S6.md): last digit = first significant digit of
# the MCSE, the largest MCSE of the column fixing the column; MCSE in the caption when the
# largest MCSE of the column is at most twice the smallest, otherwise "estimate (MCSE)".
# Outputs: tables/V4_S6_table_cells.csv, reports/morris/S6_tables.md
# =============================================================================
source(path.expand("~/MMM/R/s6/s6_common.R"))
VAL <- fread(file.path(TABDIR, "V4_S6_values.csv"))
LAB <- c(gee_fg = "GEE, Fay-Graubard", gee_md = "GEE, Mancl-DeRouen", mmm = "Marginalized model", glmer = "GLMM, Wald",
         shrink = "Shrinkage GLM", gee = "GEE, uncorrected", glmer_lrt = "GLMM, likelihood ratio", gee_kc = "GEE, Kauermann-Carroll")
ORD <- names(LAB); ORDW <- setdiff(ORD, "glmer_lrt")
MECH <- c(GB = "GenBinaryY", NAT = "GLMM-native", AR1 = "GLMM-AR1")
VARLAB <- c(main = "main grid, 1,440 design points", m1c = "cohort-anchored scenario M1c, 72 design points")
digits_of <- function(se) { se <- se[is.finite(se) & se > 0]; if (!length(se)) return(1L); as.integer(min(3, max(0, -floor(log10(max(se)))))) }
ceil1 <- function(x) { e <- floor(log10(x)); v <- ceiling(x / 10^e); if (v < 10) v * 10^e else ceiling(x / 10^(e - 1)) * 10^(e - 1) }
sfmt <- function(x) trimws(formatC(x, format = "fg", digits = 2))
se1 <- function(x) ifelse(x == 0, "0", trimws(formatC(signif(x, 1), format = "fg", digits = 1)))
fmt <- function(x, d) { d <- rep_len(d, length(x)); unname(mapply(function(a, b) if (is.finite(a)) formatC(a, format = "f", digits = b) else "n.a.", x, d)) }
CELLS <- list(); MD <- character()
# spec: rows = data.table(rk, label); cols = list of list(label, fun(key) -> value name, count = FALSE)
emit <- function(id, title, caption, variant_, rows, cols, rowhead = "Method") {
  M <- rbindlist(lapply(seq_along(cols), function(j) { cc <- cols[[j]]
    nm <- vapply(rows$rk, cc$fun, character(1))
    v <- VAL[variant == variant_, .(name, est, se)][data.table(name = nm), on = "name"]
    data.table(table_id = id, variant = variant_, row = rows$label, col = cc$label, j = j, i = seq_len(nrow(rows)), name = nm, est = v$est, se = v$se,
               count = isTRUE(cc$count)) }))
  M[grepl("^(ndp|ncell20)\\|", name), `:=`(count = TRUE, se = NA_real_)]
  M[, digits := if (all(count)) 0L else digits_of(se[!count]), by = j]
  M[count == TRUE, digits := 0L]
  M[, mode := { s <- se[is.finite(se) & se > 0]; if (all(count) || !length(s)) "none" else if (max(s) <= 2 * min(s)) "caption" else "inline" }, by = j]
  M[, text := fmt(est, digits)]
  M[count == TRUE, mode := "none"]
  M[mode == "inline" & is.finite(se), text := sprintf("%s (%s)", text, se1(se))]
  CELLS[[id]] <<- M
  capm <- M[mode == "caption" & is.finite(se), .(mx = max(se, na.rm = TRUE)), by = col]
  mc <- if (nrow(capm)) paste0(" Monte Carlo SE (block bootstrap, B = 2000): at most ", paste(sprintf("%s for %s", vapply(capm$mx, function(z) sfmt(ceil1(z)), ""), capm$col), collapse = "; "), ".") else ""
  inl <- if (any(M$mode == "inline")) " Where an entry carries a value in parentheses, it is its Monte Carlo SE." else ""
  W <- dcast(M, i + row ~ j, value.var = "text"); setorder(W, i)
  hdr <- c(rowhead, vapply(cols, `[[`, character(1), "label"))
  MD <<- c(MD, sprintf("### %s. %s (%s)", id, title, VARLAB[[variant_]]), "",
           paste0("| ", paste(hdr, collapse = " | "), " |"), paste0("|", paste(rep("---", length(hdr)), collapse = "|"), "|"),
           apply(W[, -1], 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), "",
           paste0("*", caption, mc, inl, "*"), "")
}
mrows <- function(keys) data.table(rk = keys, label = LAB[keys])
for (vr in c("main", "m1c")) {
  sfx <- if (vr == "m1c") "c" else ""
  Tlev <- if (vr == "main") c(2, 3, 4, 6, 8, 10) else c(2, 3, 4, 6)
  # Table 3: type I error
  cols <- c(lapply(names(MECH), function(G) list(label = paste("rate", G), fun = function(m) sprintf("rate|%s|%s", G, m))),
            lapply(names(MECH), function(G) list(label = paste("worst", G), fun = function(m) sprintf("wc|%s|%s", G, m))),
            lapply(names(MECH), function(G) list(label = paste("MAD", G), fun = function(m) sprintf("mad|%s|%s", G, m))))
  emit(paste0("T3", sfx), "Type I error at the nominal 5% (%)",
       "Rejection rate pooled over the grid, own usable replicates (usability U0). The rate is replicate-weighted within each cluster count; the worst case is the largest absolute departure from 5 over the cluster counts and MAD their mean. GB = GenBinaryY, NAT = GLMM-native, AR1 = GLMM-AR1.",
       vr, mrows(ORD), cols)
  emit(paste0("T4", sfx), "Mean absolute departure from 5% by number of occasions (percentage points)",
       sprintf("Cell-weighted mean over design points and mechanisms. The last column pools T = %s.", paste(Tlev[Tlev >= 4], collapse = ", ")),
       vr, mrows(ORD), list(list(label = "T = 2", fun = function(m) sprintf("madT|T2|%s", m)), list(label = "T = 3", fun = function(m) sprintf("madT|T3|%s", m)),
                            list(label = "T >= 4", fun = function(m) sprintf("madT|T4plus|%s", m))))
  emit(paste0("T5", sfx), "Method-cells detectably off nominal (%)",
       "Exact binomial test of each method-cell with at least one usable replicate against 5%, false discovery rate controlled within method at 5% by Benjamini-Hochberg (BH) or Benjamini-Yekutieli (BY).",
       vr, mrows(ORD), list(list(label = "BH", fun = function(m) sprintf("flag|BH|%s", m)), list(label = "BY", fun = function(m) sprintf("flag|BY|%s", m))))
  cols <- c(lapply(names(MECH), function(G) list(label = paste("own", G), fun = function(m) sprintf("rate|%s|%s", G, m))),
            lapply(names(MECH), function(G) list(label = paste("common", G), fun = function(m) sprintf("irate|%s|%s", G, m))),
            lapply(names(MECH), function(G) list(label = paste("shift", G), fun = function(m) sprintf("shift|%s|%s", G, m))))
  emit(paste0("T6", sfx), "Type I error on each method's own replicates and on the common set (%)",
       "Common set: the replicates on which all eight methods returned a usable result. Shift = common minus own, in percentage points.", vr, mrows(ORD), cols)
  grps <- if (vr == "main") c("RI", "RS") else "RI"
  rr <- CJ(grp = grps, T = Tlev, sorted = FALSE)[, `:=`(rk = paste(grp, T, sep = ":"), label = sprintf("%s, T = %d", fifelse(grp == "RI", "random intercept", "random slopes"), T))]
  cols <- lapply(ORD, function(m) list(label = LAB[[m]], fun = function(k) { p <- strsplit(k, ":")[[1]]; sprintf("puT|%s|%s|%s", p[1], m, p[2]) }))
  emit(paste0("T7", sfx), "Replicates returning a usable result, by number of occasions (%)",
       if (vr == "main") "Pooled over mechanisms and cluster counts, every design point in the denominator. Random intercept only: scenarios M1-M2; random slopes: M3-M5." else "Pooled over mechanisms and cluster counts, every design point in the denominator. M1c has a random intercept only.",
       vr, rr[, .(rk, label)], cols, rowhead = "Structure, occasions")
  emit(paste0("T8", sfx), "Usable result and unconditional validity at 20 clusters or fewer (%)",
       "Means over design points at N <= 20, pooled over mechanisms and scenarios, cells a method never fits counted at zero. Unconditional validity = probability of a usable result that does not falsely reject.",
       vr, mrows(ORD), list(list(label = "usable", fun = function(m) sprintf("use20|%s", m)), list(label = "unconditional validity", fun = function(m) sprintf("val20|%s", m)),
                            list(label = "design points with a usable replicate", fun = function(m) sprintf("ncell20|%s", m), count = TRUE)))
  slab <- c(s1 = "<50%", s2 = "50-80%", s3 = "80-95%", s4 = "95-99.9%", s5 = "100%")
  emit(paste0("T9", sfx), "Mean absolute departure from 5% by estimability stratum (percentage points)",
       "Design points stratified by the mean usable proportion of the five marginal methods (right-closed strata). The last row gives the design points per stratum.",
       vr, rbind(mrows(ORD), data.table(rk = "ndp", label = "design points")),
       lapply(names(slab), function(s) list(label = slab[[s]], fun = function(m) if (m == "ndp") sprintf("ndp|%s", s) else sprintf("madS|%s|%s", s, m))))
  rr <- CJ(G = names(MECH), m = ORD, sorted = FALSE)[, `:=`(rk = paste(G, m, sep = ":"), label = paste(MECH[G], LAB[m], sep = ", "))]
  emit(paste0("T10", sfx), "Raw and size-adjusted power at the interaction contrast (%)",
       "Pooled by mechanism. Size-adjusted: each cell's own empirical 5% critical value from the null arm (for the likelihood-ratio test on the p-value scale).",
       vr, rr[, .(rk, label)], list(list(label = "raw", fun = function(k) { p <- strsplit(k, ":")[[1]]; sprintf("raw|%s|%s", p[1], p[2]) }),
                                     list(label = "size-adjusted", fun = function(k) { p <- strsplit(k, ":")[[1]]; sprintf("sa|%s|%s", p[1], p[2]) })),
       rowhead = "Mechanism, method")
  cols <- c(lapply(names(MECH), function(G) list(label = paste("normal", G), fun = function(m) sprintf("rate|%s|%s", G, m))),
            lapply(names(MECH), function(G) list(label = paste("t", G), fun = function(m) sprintf("ratet|%s|%s", G, m))),
            lapply(names(MECH), function(G) list(label = paste("shift", G), fun = function(m) sprintf("tshift|%s|%s", G, m))),
            list(list(label = "shift at N <= 25", fun = function(m) sprintf("tsh25|%s", m))))
  emit(paste0("T11", sfx), "Type I error under a normal and a t reference (%)",
       "t reference with N - p degrees of freedom. Shift = t minus normal, percentage points; the last column is normal minus t pooled over mechanisms at N <= 25 (a positive value: the t reference rejects less).",
       vr, mrows(ORDW), cols)
  emit(paste0("T12", sfx), "Coverage of the 95% confidence interval (%)",
       "Scored against each method's own estimand (marginal for the marginal methods, conditional for the GLMM, which has no conditional true value on GenBinaryY: not available); H1 arm, cells with at least two usable replicates.",
       vr, mrows(ORDW), lapply(names(MECH), function(G) list(label = MECH[[G]], fun = function(m) sprintf("cov|%s|%s", G, m))))
}
CELLS <- rbindlist(CELLS)
fwrite(CELLS, file.path(TABDIR, "V4_S6_table_cells.csv"))
writeLines(c("Generated by `R/s6/s6_tables.R` from `tables/V4_S6_values.csv`; every entry is listed with its source in `tables/V4_S6_table_cells.csv`.", "", MD),
           file.path(REPDIR, "S6_tables.md"))
cat(sprintf("%d cells in %d tables; missing values %d; inline-MCSE columns %d of %d\n", nrow(CELLS), uniqueN(CELLS$table_id), CELLS[!is.finite(est), .N],
            uniqueN(CELLS[mode == "inline", .(table_id, j)]), uniqueN(CELLS[, .(table_id, j)])))
print(CELLS[!is.finite(est), .N, by = .(table_id, col)])
