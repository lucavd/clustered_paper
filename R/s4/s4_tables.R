# =============================================================================
# s4_tables.R -- builds reports/morris/S4_results.md and tables/V4_S4_expectations.csv
# from the files written by s4_q0.sh, s4_truth.R, agg_v4.R and s4_checks.R.
# Every number is read from a file.
# =============================================================================
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")); suppressMessages(library(data.table))
REPO <- path.expand("~/MMM"); TAB <- file.path(REPO, "tables"); LOGD <- file.path(REPO, "reports/morris/S4_logs")
rd <- function(f) fread(file.path(TAB, f)); lg <- function(f) readLines(file.path(LOGD, f), warn = FALSE)
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d, decimal.mark = ",")
int <- function(x) formatC(x, format = "d", big.mark = ".")
sci <- function(x) sprintf("%.1e", x)
md <- function(dt) { dt <- as.data.frame(dt); dt[] <- lapply(dt, function(v) gsub("\\|", "\\\\|", as.character(v)))
  c(paste0("| ", paste(names(dt), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(dt)), collapse = "|"), "|"),
    apply(dt, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))) }
yn <- function(b) ifelse(b, "si", "no")
L <- character(); add <- function(...) L <<- c(L, ...)

q0pre <- grep("^Q0 (results|rows)", lg("S4_q0_pre.log"), value = TRUE)
q0post <- if (file.exists(file.path(LOGD, "S4_q0_post.log"))) grep("^Q0 (results|rows)", lg("S4_q0_post.log"), value = TRUE) else "non eseguito"
q0ok <- function(x) all(grepl(" 0 FAIL lines|092a5ad60b60e0ea28cb6b78ba2adbeb \\(expected", x)) && length(x) == 4
Q12 <- rd("V4_S4_q1_q2.csv"); Q3 <- rd("V4_S4_q3_reproduction.csv"); Q4 <- rd("V4_S4_q4_usable_vs_S2.csv")
Q5 <- rd("V4_S4_q5_formulas.csv"); Q6 <- rd("V4_S4_q6_truth_checks.csv"); Q7 <- rd("V4_S4_q7_bootstrap.csv")
Q9 <- rd("V4_S4_q9_bobyqa.csv"); D8 <- rd("V4_S4_d08_lrt_effect.csv"); D8S <- rd("V4_S4_d08_summary.csv")
INV <- fread(file.path(LOGD, "S4_q1c_investigation.csv"))
P <- fread(file.path(TAB, "V4_perf_all.csv.gz"))
g <- function(lbl) Q12[check == lbl]
q5r <- Q5[grepl("^rsimsum", check)]; q5i <- Q5[grepl("^identity", check)]
q8 <- D8[group == "D-08 only", .(cells = .N, ok = sum(n_both >= 1))]
E <- rbindlist(list(
  list("Q0", "md5 di results_full/ (28.800), results_v4 S2 (126) e S3 (18.720) e di rows_timetrt.rds invariati prima e dopo", "zero falliti",
       paste(c(paste("prima:", paste(sub("^Q0 ", "", q0pre), collapse = "; ")), paste("dopo:", paste(sub("^Q0 ", "", q0post), collapse = "; "))), collapse = " || "),
       yn(q0ok(q0pre) && q0ok(q0post))),
  list("Q1a", "nessuna riga glmer_lrt v3 in M3-M5; re_structure come D-08 in M3-M5 e (1|id) in M1c", "esatta",
       sprintf("righe v3: %s; re_structure inattese: %s", int(g("Q1a glmer_lrt rows from results_full in M3-M5")$observed),
               int(g("Q1a glmer_lrt rows with unexpected re_structure (M3-M5, M1c)")$observed)),
       yn(g("Q1a glmer_lrt rows from results_full in M3-M5")$pass && g("Q1a glmer_lrt rows with unexpected re_structure (M3-M5, M1c)")$pass)),
  list("Q1b", "celle F1 e F2 interamente da S3", "esatta", sprintf("righe non da S3: %s", int(g("Q1b F1/F2 cell rows not from S3")$observed)), yn(g("Q1b F1/F2 cell rows not from S3")$pass)),
  list("Q1c", "3.024 celle, 1.512 design point, ogni cella con righe per tutti gli 8 metodi primari", "esatta",
       sprintf("celle %s, design point %s; celle con meno di 8 metodi: %s (qualunque riga), %s (righe time:trt); duplicati %s",
               int(g("Q1c cells")$observed), int(g("Q1c design points")$observed),
               int(g("Q1c cells with fewer than 8 primary methods (any row)")$observed),
               int(g("Q1c cells with fewer than 8 primary methods (time:trt rows; descriptive)")$observed),
               int(g("Q1 duplicate (cell, method, sim_id) rows")$observed)),
       yn(g("Q1c cells")$pass && g("Q1c design points")$pass && g("Q1c cells with fewer than 8 primary methods (any row)")$pass && g("Q1 duplicate (cell, method, sim_id) rows")$pass)),
  list("Q1d", "tabella lunga 3.024 x 9 x 3 x 2 = 163.296 righe", "esatta", int(nrow(P)), yn(nrow(P) == 163296L)),
  list("Q2", "n_den: 2.878.518 (griglia) e 143.742 (M1c); zero righe U0 su chiavi degeneri; ogni chiave non degenere con una riga di analisi", "esatta",
       sprintf("%s; %s; U0 su degeneri %s; chiavi senza righe %s", int(g("Q2 sum n_den main grid")$observed), int(g("Q2 sum n_den M1c")$observed),
               int(g("Q2 U0 rows on degenerate keys")$observed), int(g("Q2 non-degenerate keys without any analysis row")$observed)),
       yn(all(Q12[grepl("^Q2", check), pass]))),
  list("Q3", "righe non toccate da S3 riproducono O_level_all (n, k_rej) e O_coverage_all (n, cov_gen)", "100%, 1e-12",
       paste(Q3[, sprintf("%s: %s/%s n uguali, diff. max %s, mancanti %s%s", table, int(n_equal), int(rows), sci(max_diff), int(missing),
                          ifelse(is.na(reverse_missing), "", sprintf(", inverse mancanti %s", int(reverse_missing))))], collapse = "; "), yn(all(Q3$pass))),
  list("Q4", "n_U0, n_U1, n_U2 = V4_S2_usable_by_cell nelle combinazioni non toccate", "esatta",
       Q4[, sprintf("%s righe; U0 %s, U1 %s, U2 %s uguali; mancanti %s", int(rows), int(U0_equal), int(U1_equal), int(U2_equal), int(missing))], yn(Q4$pass)),
  list("Q5a", "rsimsum 0.13.1 su 300 combinazioni: 17 misure e MCSE entro 1e-10", "1e-10",
       sprintf("%s confronti, oltre tolleranza %s, diff. max %s", int(sum(q5r$compared)), int(sum(q5r$over_tol)), sci(max(q5r$max_abs_diff))), yn(sum(q5r$over_tol) == 0)),
  list("Q5b", "MSE = bias^2 + SE empirico^2 (n - 1)/n in tutte le combinazioni con n_obs >= 2", "1e-12",
       paste(q5i[, sprintf("%s: %s/%s oltre, diff. max %s", sub("identity MSE, ", "", check), int(over_tol), int(combos), sci(max_abs_diff))], collapse = "; "),
       yn(sum(q5i$over_tol) == 0)),
  list("Q6", "verita' AR1 per replicato (riproduzione S2, biiezione con i degeneri, replay = generatore); quadratura M1", "1e-12; < 1e-8",
       paste(Q6[, sprintf("%s: %s", sub("^Q6[abc] ", "", check), ifelse(pass, "ok", "NO"))], collapse = "; "), yn(all(Q6$pass))),
  list("Q7", "SE bootstrap della potenza grezza entro il 10% dalla MCSE binomiale", ">= 95%",
       Q7[, sprintf("%s/%s (%s); rapporto mediano %s", int(within_10pct), int(combos), fmt(share, 4), fmt(median_ratio, 3))], yn(Q7$pass)),
  list("Q8", "D-08: ogni cella del blocco L con almeno un replicato usabile in entrambe le versioni", "100%",
       sprintf("%s/%s", int(q8$ok), int(q8$cells)), yn(q8$ok == q8$cells)),
  list("Q9", "glmer_bobyqa: p = 2 Phi(-abs(z)); n_obs <= n_obs(glmer_lrt) per cella", "esatta",
       Q9[, sprintf("p diversi %s su %s righe; celle %s/%s (eccesso massimo %s)", int(p_mismatch), int(rows), int(cells_bobyqa_le_lrt), int(cells), int(max_excess))], yn(Q9$pass))
))
setnames(E, c("id", "attesa", "tolleranza", "osservato", "superato"))
fwrite(E, file.path(TAB, "V4_S4_expectations.csv"))
add("### Tabella E. Attese pre-registrate e osservato", "", md(E), "")

# E1: Q1c investigation
iv <- INV[, .(`coppie cella x metodo` = .N, `con righe` = sum(n_rows > 0), `solo righe di fallimento` = sum(n_rows > 0 & n_rows == n_failure_rows),
              `glmer mai convergente nella cella` = sum(!glmer_converged_any), `motivo piu' frequente` = top_status[which.max(n_rows)]),
          by = .(metodo = method, scenario = scen, T = ntime)][order(metodo, scenario, T)]
add("### E1. Q1c: celle senza righe time:trt per un metodo primario", "", md(iv), "")

# E2: D-04 usability, pooled main grid
U <- P[set == "own" & scen != "M1c", .(n_obs = sum(n_obs), n_den = sum(n_den)), by = .(method, usab)]
U <- dcast(U[, .(method, usab, share = fmt(n_obs / n_den, 3))], method ~ usab, value.var = "share")
add("### E2. Quota usabile per definizione (griglia principale, entrambi gli arm, denominatore D-11)", "", md(U), "")

# E3: MCSE overview, primary analysis, cells with n_obs >= 100
pr <- P[usab == "U0" & set == "own" & n_obs >= 100]
mc <- rbindlist(lapply(list(c("rej", "mcse_rej", "H0"), c("rej", "mcse_rej", "H1"), c("cover", "mcse_cover", "H1"), c("becover", "mcse_becover", "H1"),
                            c("bias", "mcse_bias", "H1"), c("empse", "mcse_empse", "H1"), c("relerr_modse", "mcse_relerr", "H1"),
                            c("power_sa", "mcse_power_sa", "H1"), c("p_usable", "mcse_p_usable", "H0")), function(v) {
  x <- pr[arm == v[3] & is.finite(get(v[2]))][[v[2]]]
  data.table(misura = v[1], arm = v[3], `combinazioni` = int(length(x)), `MCSE mediana` = sprintf("%.4f", median(x)), `MCSE massima` = sprintf("%.4f", max(x)))
}))
add("### E3. MCSE nell'analisi primaria (U0, own; celle con n_obs >= 100)", "",
    "Scala: proporzioni (I tipo, copertura, potenza, quota usabile), unita' del coefficiente (bias, SE empirico), punti percentuali (errore relativo dello SE).", "", md(mc), "")

# E4: truths
TT <- fread(file.path(TAB, "V4_S4_truth_m1c_native.csv"))
add("### E4. Valori veri di M1c in GLMM-native (quadratura, Q6c descrittivo)", "",
    md(TT[, .(arm = ifelse(label == "mmlb", "H1", "H0"), T = ntime, condizionale = fmt(truth_conditional, 6), marginale = fmt(truth_marginal, 6), `errore 40 vs 80 nodi` = sci(quad_err_40v80))]), "")
add(md(Q6[, .(controllo = check, chiavi = int(keys), valore = sci(value), `discordanze NA` = int(na_mismatch), superato = yn(pass))]), "")

# E5: rsimsum detail
add("### E5. Controprova rsimsum (Q5a) e identita' (Q5b)", "", md(Q5[, .(controllo = check, combinazioni = int(combos), confronti = int(compared), `diff. max` = sci(max_abs_diff), `oltre tolleranza` = int(over_tol), tolleranza = sci(tol))]), "")

# E6: D-08
add("### E6. Effetto di D-08 sul LRT: differenza appaiata v4 - v3 (U0), media sulle celle", "",
    "d I tipo/potenza: differenza del tasso di rifiuto sui replicati usabili in entrambe le versioni; d usabile: differenza della quota usabile sui replicati con disegno stimabile. Celle oltre 2 MCSE appaiate.", "",
    md(D8S[, .(gruppo = group, meccanismo = dgp, scenario = scen, arm, celle = int(cells), `d rifiuto medio` = fmt(mean_d_rej, 4),
               `celle oltre 2 MCSE (rifiuto)` = int(cells_rej_over2), `d usabile medio` = fmt(mean_d_usable, 4),
               `celle oltre 2 MCSE (usabile)` = int(cells_usable_over2), `n appaiati minimo` = int(min_n_both))]), "")

# E7: Wald default vs bobyqa (descriptive)
wb <- P[usab == "U0" & set == "own" & method %in% c("glmer", "glmer_bobyqa") & n_obs >= 100 & scen != "M1c"]
wb <- dcast(wb[, .(dgp, scen, N, ntime, arm, method, rej, mcse_rej, p_usable)], dgp + scen + N + ntime + arm ~ method, value.var = c("rej", "mcse_rej", "p_usable"))
wb[, z := (rej_glmer_bobyqa - rej_glmer) / sqrt(mcse_rej_glmer^2 + mcse_rej_glmer_bobyqa^2)]
wbs <- wb[is.finite(z), .(celle = .N, `rifiuto glmer` = fmt(mean(rej_glmer), 4), `rifiuto bobyqa` = fmt(mean(rej_glmer_bobyqa), 4),
                          `celle con abs(z) > 2 (z non appaiato, conservativo)` = sum(abs(z) > 2), `quota usabile glmer` = fmt(mean(p_usable_glmer), 3),
                          `quota usabile bobyqa` = fmt(mean(p_usable_glmer_bobyqa), 3)), by = .(meccanismo = dgp, arm)][order(meccanismo, arm)]
add("### E7. Wald glmer (default) e Wald bobyqa (D-16), media sulle celle con n_obs >= 100 per entrambi", "", md(wbs), "")
writeLines(L, file.path(REPO, "reports/morris/S4_results.md"))
cat(sprintf("S4_results.md written: %d lines; expectations passed %d/%d\n", length(L), sum(E$superato == "si"), nrow(E)))
print(E[, .(id, superato)])
