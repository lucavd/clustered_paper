# =============================================================================
# s6_check_g3.R -- S6 gate G3: every number of the S6 tables and of the S6 figure data is
# re-read from its source file, independently of the code that wrote it.
#  tables : parse reports/morris/S6_tables.md; each cell is matched by position to
#           tables/V4_S6_table_cells.csv (source name) and compared with tables/V4_S6_values.csv
#           (estimate within half a unit of the last printed digit; MCSE within half a unit of
#           its single printed significant digit)
#  figures: reports/morris/S6_figures/S6_figures_data.csv against V4_S6_values.csv
#           (estimate, estimate -/+ 2 MCSE), V4_S6_nestedloop.csv, V4_S6_rule_safety.csv and V4_S6_zip.csv
# Output: tables/V4_S6_g3.csv
# =============================================================================
source(path.expand("~/MMM/R/s6/s6_common.R"))
VAL <- fread(file.path(TABDIR, "V4_S6_values.csv")); CEL <- fread(file.path(TABDIR, "V4_S6_table_cells.csv"))
L <- readLines(file.path(REPDIR, "S6_tables.md"))
hd <- grep("^### ", L); res <- list()
for (h in hd) {
  id <- sub("^### ([^.]+)\\..*$", "\\1", L[h]); k <- h + 4; rows <- character()
  while (k <= length(L) && startsWith(L[k], "|")) { rows <- c(rows, L[k]); k <- k + 1 }
  C <- CEL[table_id == id]
  for (ri in seq_along(rows)) {
    f <- trimws(strsplit(sub("^\\|", "", sub("\\|$", "", rows[ri])), "|", fixed = TRUE)[[1]])[-1]
    for (cj in seq_along(f)) {
      cc <- C[which(C$i == ri & C$j == cj)]
      v <- VAL[variant == cc$variant & name == cc$name]
      if (!nrow(v)) v <- data.table(est = NA_real_, se = NA_real_)   # value not defined (printed "n.a.")
      num <- suppressWarnings(as.numeric(sub("^(-?[0-9.]+).*$", "\\1", f[cj])))
      se_p <- suppressWarnings(as.numeric(sub("^.*\\(([0-9.]+)\\)$", "\\1", f[cj])))
      if (!grepl("\\(", f[cj])) se_p <- NA_real_
      d <- if (grepl("\\.", sub(" .*$", "", f[cj]))) nchar(sub("^-?[0-9]+\\.([0-9]+).*$", "\\1", f[cj])) else 0
      ok_est <- if (f[cj] == "n.a.") !is.finite(v$est) else is.finite(num) && abs(num - v$est) <= 0.5 * 10^(-d) + 1e-9
      ok_se <- if (is.na(se_p)) TRUE else if (!is.finite(v$se)) FALSE else if (v$se == 0) se_p == 0 else abs(se_p - v$se) <= 0.5 * 10^floor(log10(v$se)) + 1e-12
      res[[length(res) + 1]] <- data.table(kind = "table", id = id, i = ri, j = cj, name = cc$name, printed = f[cj], est = v$est, se = v$se, ok = ok_est && ok_se)
    }
  }
}
TB <- rbindlist(res)
FD <- fread(file.path(FIGDIR, "S6_figures_data.csv"))
NL <- fread(file.path(TABDIR, "V4_S6_nestedloop.csv")); SF <- fread(file.path(TABDIR, "V4_S6_rule_safety.csv")); Z <- fread(file.path(TABDIR, "V4_S6_zip.csv"))[usable == TRUE]
chk <- function(r) {
  if (startsWith(r$source, "V4_S6_nestedloop:")) { p <- strsplit(sub("^V4_S6_nestedloop:", "", r$source), "|", fixed = TRUE)[[1]]
    x <- NL[dgp == p[1] & scen == p[2] & N == as.integer(p[3]) & ntime == as.integer(p[4]) & method == p[5]]
    if (x$n_obs < 20) return(!is.finite(r$est)); return(abs(x$rej - r$est) < 1e-9 && abs(x$rej - 2 * x$mcse_rej - r$lo) < 1e-9) }
  if (startsWith(r$source, "V4_S6_rule_safety:")) { p <- strsplit(sub("^V4_S6_rule_safety:", "", r$source), "|", fixed = TRUE)[[1]]
    x <- SF[strategy == paste(p[1], p[2], sep = "|") & N == as.integer(p[3]) & band == "all"]; return(nrow(x) == 1 && abs(100 * x$n_hi / x$n_dp - r$est) < 1e-9) }
  if (startsWith(r$source, "V4_S6_zip:")) { p <- strsplit(sub("^V4_S6_zip:", "", r$source), "|", fixed = TRUE)[[1]]
    x <- Z[N == as.integer(p[1]) & method == p[2]]; cv <- 100 * mean(abs(x$estimate - x$truth) <= qnorm(0.975) * x$se); return(abs(cv - r$est) < 1e-6) }
  v <- VAL[variant == "main" & name == r$source]
  abs(v$est - r$est) < 1e-9 && abs(v$est - 2 * v$se - r$lo) < 1e-9 && abs(v$est + 2 * v$se - r$hi) < 1e-9 }
FD[, ok := vapply(seq_len(.N), function(k) isTRUE(chk(FD[k])), logical(1))]
G3 <- rbind(TB[, .(kind, id, name, printed = as.character(printed), ok)], FD[, .(kind = "figure", id = paste(figure, panel), name = source, printed = as.character(est), ok)])
fwrite(G3, file.path(TABDIR, "V4_S6_g3.csv"))
cat(sprintf("G3: table cells %d checked, %d failing; figure values %d checked, %d failing\n", nrow(TB), TB[ok == FALSE, .N], nrow(FD), FD[ok == FALSE, .N]))
print(head(G3[ok == FALSE], 20))
