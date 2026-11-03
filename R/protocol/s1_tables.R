# =============================================================================
# s1_tables.R -- builds the markdown tables of the S1 report and of
# docs/PROTOCOL_ADEMP.md from the CSV files written by the extraction and the
# tests, so that no number in either document is typed by hand.
#   in : tables/V4_dgm_parameters.csv, tables/V4_S1_tests_P1_P4.csv,
#        tables/V4_S1_tests_P5.csv, tables/V4_S1_tests_P5_mersenne.csv
#   out: reports/morris/S1_params.md, reports/morris/S1_results.md,
#        and the block between <!-- PARAMS:BEGIN --> / <!-- PARAMS:END --> in
#        docs/PROTOCOL_ADEMP.md
# =============================================================================
REPO <- path.expand("~/MMM")
rd <- function(f) read.csv(file.path(REPO, f), stringsAsFactors = FALSE)
tab <- rd("tables/V4_dgm_parameters.csv")
t14 <- rd("tables/V4_S1_tests_P1_P4.csv")
t5  <- rd("tables/V4_S1_tests_P5.csv")
t5m <- rd("tables/V4_S1_tests_P5_mersenne.csv")
md5c <- rd("tables/V4_S1_md5check.csv")
v <- function(g, m, a, p) {
  x <- tab$value[tab$dgp == g & tab$model == m & tab$arm == a & tab$parameter == p]
  if (length(x) == 1L) x else ""
}
tm <- function(g, m, p) {
  x <- tab$terms[tab$dgp == g & tab$model == m & tab$arm == "H1" & tab$parameter == p]
  if (length(x) == 1L && !is.na(x)) x else ""
}
cell <- function(s) gsub("|", "\\|", s, fixed = TRUE)

# ---- parameter summary --------------------------------------------------------------
L <- c("| DGM | M | mean model (coefficient order) | coefficients H1 | dependence (H1) | H0 differs from H1 in |",
       "|---|---|---|---|---|---|")
for (g in c("GenBinaryY", "GLMM-native", "GLMM-AR1")) for (m in 1:5) {
  bp <- switch(g, "GenBinaryY" = "beta", "GLMM-native" = "beta_c", "GLMM-AR1" = "target_beta_m")
  dep <- switch(g,
    "GenBinaryY"  = sprintf("lv %s, sigma (%s); t %s, gamma (%s); q = %s",
                            v(g, m, "H1", "lv.formula"), v(g, m, "H1", "sigma"),
                            v(g, m, "H1", "t.formula"), v(g, m, "H1", "gamma"), v(g, m, "H1", "q")),
    "GLMM-native" = sprintf("b_i on %s, D = [%s]", v(g, m, "H1", "re.formula"), v(g, m, "H1", "D")),
    "GLMM-AR1"    = sprintf("b_i on %s, D = [%s]; AR(1) rho = %s, sigma2_ar = %s",
                            v(g, m, "H1", "re.formula"), v(g, m, "H1", "D"),
                            v(g, m, "H1", "rho"), v(g, m, "H1", "sigma2_ar")))
  par <- t14$detail[t14$test == "P2" & t14$item == sprintf("%s M%d parity", g, m)]
  par <- sub("^also differ: ", "time:trt and ", par)
  if (g == "GenBinaryY" && m >= 3) {
    par <- sprintf("%s (H0: sigma (%s), gamma (%s))", par, v(g, m, "H0", "sigma"), v(g, m, "H0", "gamma"))
  }
  L <- c(L, sprintf("| %s | %d | %s (%s) | %s %s | %s | %s |", g, m,
                    v(g, m, "H1", "mean.formula"), tm(g, m, bp),
                    switch(g, "GenBinaryY" = "marginal", "GLMM-native" = "conditional",
                           "GLMM-AR1" = "marginal target"),
                    v(g, m, "H1", bp), cell(dep), par))
}
g <- "GenBinaryY"
dr <- tab[tab$dgp == "GenBinaryY" & tab$model == 5 & tab$arm == "H1" & grepl("^draw_", tab$parameter), ]
# draws per scenario must be the same in the three mechanisms and both arms
draws_of <- function(gg, m, a) {
  x <- tab[tab$dgp == gg & tab$model == m & tab$arm == a & grepl("^draw_", tab$parameter), ]
  # same variable, same function and arguments, with the sample size written as its level
  # (GenBinaryY/native write rnorm(n = length(id)), AR1 rnorm(n = nrow(d)): same draw)
  lvl <- sub(".*level=([a-z]+);.*", "\\1", x$note)
  val <- mapply(function(v, l) sub("(n|size) = [^,)]+", paste0("\\1 = <", l, ">"), v), x$value, lvl)
  paste(x$terms, val, sub("; assigned_as=.*", "", x$note), sep = "=", collapse = " ; ")
}
for (m in 1:5) stopifnot(length(unique(c(sapply(c("GenBinaryY", "GLMM-native", "GLMM-AR1"),
  function(gg) c(draws_of(gg, m, "H1"), draws_of(gg, m, "H0")))))) == 1L)
per_m <- vapply(1:5, function(m) {
  x <- tab[tab$dgp == "GenBinaryY" & tab$model == m & tab$arm == "H1" & grepl("^draw_", tab$parameter), ]
  sprintf("M%d: %s", m, paste(unique(x$terms), collapse = ", "))
}, "")
D <- c("", "Design draws, in execution order (table rows `draw_*`). They are the same in the three mechanisms and in both arms (checked); the scenarios differ only in how many covariates they draw:", "",
       paste0("- ", per_m), "", "Full list (M5):", "",
       sprintf("- follow-up: `nclust = %s` (p_full_followup = %s)",
               v(g, 5, "H1", "followup_rule"), v(g, 5, "H1", "p_full_followup")))
for (k in seq_len(nrow(dr)))
  D <- c(D, sprintf("- `%s` <- `%s` (%s)", dr$terms[k], dr$value[k], dr$note[k]))
U <- c("", "Covariates drawn but entering no part of the DGM (they consume random numbers):", "")
for (gg in c("GenBinaryY", "GLMM-native", "GLMM-AR1")) for (m in 1:5) {
  u <- v(gg, m, "H1", "unused_covariates")
  if (u != "none") U <- c(U, sprintf("- %s M%d: %s", gg, m, u))
}
R <- c("", sprintf("Seed: `set.seed(%s)` per replicate; RNGkind in the workers: `%s` (%s).",
                   v(g, 1, "H1", "seed_rule"), v(g, 1, "H1", "RNGkind"),
                   tab$note[tab$parameter == "RNGkind"][1]))
params_md <- c(L, D, U, R)
writeLines(params_md, file.path(REPO, "reports/morris/S1_params.md"))

pf <- file.path(REPO, "docs/PROTOCOL_ADEMP.md")
if (file.exists(pf)) {
  p <- readLines(pf)
  b <- grep("<!-- PARAMS:BEGIN -->", p, fixed = TRUE); e <- grep("<!-- PARAMS:END -->", p, fixed = TRUE)
  stopifnot(length(b) == 1L, length(e) == 1L, e > b)
  writeLines(c(p[seq_len(b)], params_md, p[e:length(p)]), pf)
}

# ---- expectation results ------------------------------------------------------------
cnt <- function(t) sprintf("%d/%d", sum(t$pass), nrow(t))
p1f <- t14[t14$test == "P1" & !t14$pass, ]
p2f <- t14[t14$test == "P2" & !t14$pass, ]
p3 <- t14[t14$test == "P3", ]; p4 <- t14[t14$test == "P4", ]
p3err <- sum(grepl("generation error", p3$detail))
res <- c("| # | attesa (pre-registrata, commit b0826c5) | osservato | esito |", "|---|---|---|---|",
  sprintf("| P1 | 30 combinazioni; valori = audit B1, B3, B4 (esatti) | %d combinazioni; %s controlli; falliti: %s | %s |",
          nrow(unique(tab[, c("dgp", "model", "arm")])), cnt(t14[t14$test == "P1", ]),
          paste(sprintf("%s (%s)", p1f$item, p1f$detail), collapse = "; "),
          if (nrow(p1f)) "non superato alla lettera (vedi sez. 3)" else "superato"),
  sprintf("| P2 | 12/15 coppie con parita; eccezioni GenBinaryY M3, M4, M5 | %s; eccezioni: %s | %s |",
          cnt(t14[t14$test == "P2", ]), paste(sprintf("%s (%s)", sub(" parity", "", p2f$item), p2f$detail), collapse = "; "),
          if (identical(sort(p2f$item), sprintf("GenBinaryY M%d parity", 3:5))) "superato" else "non superato"),
  sprintf("| P3 | 240/240 dataset bit-identici dal solo contenuto della tabella | %s (errori di generazione in entrambi i percorsi: %d) | %s |",
          cnt(p3), p3err, if (all(p3$pass) && nrow(p3) == 240) "superato" else "non superato"),
  sprintf("| P4 | 96/96 dataset M1-M2 identici fra codice del job 23035d84 e `.rundir/` | %s | %s |",
          cnt(p4), if (all(p4$pass) && nrow(p4) == 96) "superato" else "non superato"),
  sprintf("| P5 | 30 rifit i = 1, N = 25, T = 4: max abs diff <= 1e-8 | %d/%d; %d righe confrontate; max abs diff stima %.3g, SE %.3g; status identici %d/%d | %s |",
          sum(t5$pass), nrow(t5), sum(t5$rows_common), max(t5$max_abs_diff_estimate),
          max(t5$max_abs_diff_se), sum(t5$status_same), nrow(t5),
          if (all(t5$pass) && nrow(t5) == 30) "superato (secondo run, vedi sez. 3)" else "non superato"),
  sprintf("| P5' | (non pre-registrato) stesso rifit sotto Mersenne-Twister | %d/%d; max abs diff stima %.3g | controprova del RNGkind |",
          sum(t5m$pass), nrow(t5m), max(t5m$max_abs_diff_estimate)),
  sprintf("| P6 | `md5sum -c` su `results_full/` prima e dopo: 28.800/28.800 | %s | %s |",
          paste(sprintf("%s (%s): %d OK, %d falliti", md5c$stage, md5c$time, md5c$ok, md5c$failed), collapse = "; "),
          if (all(md5c$ok == 28800 & md5c$failed == 0)) "superato" else "non superato"))
writeLines(res, file.path(REPO, "reports/morris/S1_results.md"))
cat("written S1_params.md (", length(params_md), "lines), S1_results.md, protocol block:", file.exists(pf), "\n")
cat(res, sep = "\n")

# ---- expected MCSE for the n_sim justification (protocol sez. 6) --------------------
mc <- function(p, n) 100 * sqrt(p * (1 - p) / n)
mcse <- c("| quantita' | p | n | MCSE (punti percentuali) |", "|---|---|---|---|",
  sprintf("| errore di I tipo / copertura | %.2f | 1000 | %.3f |", 0.05, mc(0.05, 1000)),
  sprintf("| errore di I tipo, 50%% usabili | %.2f | 500 | %.3f |", 0.05, mc(0.05, 500)),
  sprintf("| potenza (massimo) | %.2f | 1000 | %.3f |", 0.5, mc(0.5, 1000)),
  sprintf("| potenza | %.2f | 1000 | %.3f |", 0.8, mc(0.8, 1000)),
  sprintf("| n per MCSE = 0,5 pp al 5%% | 0.05 | %d | 0.500 |", as.integer(round(0.05 * 0.95 / 0.005^2))))
writeLines(mcse, file.path(REPO, "reports/morris/S1_mcse.md"))
cat(mcse, sep = "\n")
