# =============================================================================
# s2_tables.R -- builds reports/morris/S2_results.md from the S2 tables and
# logs. Every number in the S2 report tables comes from here.
# =============================================================================
source(path.expand("~/MMM/R/s2/s2_common.R"))
LOGS <- file.path(REPO, "reports/morris/S2_logs")
rd <- function(f) fread(file.path(TABDIR, f))
f3 <- function(x, d = 3) formatC(x, format = "f", digits = d, big.mark = ".", decimal.mark = ",")
fi <- function(x) formatC(x, format = "d", big.mark = ".")
fg <- function(x, d = 2) gsub("e", " x 10^", formatC(x, format = "e", digits = d))
esc <- function(x) gsub("|", "\\|", as.character(x), fixed = TRUE)
md <- function(D) { D <- as.data.frame(D); D[] <- lapply(D, esc); names(D) <- esc(names(D)); h <- paste0("| ", paste(names(D), collapse = " | "), " |")
  s <- paste0("|", paste(rep("---", ncol(D)), collapse = "|"), "|")
  b <- apply(D, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")); c(h, s, b, "") }
loggrep <- function(f, pat) { L <- readLines(file.path(LOGS, f)); L[grepl(pat, L)] }
num_after <- function(s, pat) as.numeric(sub(paste0(".*", pat, " *([0-9.e+-]+).*"), "\\1", s))
OUT <- c()

# ---- Q0 ------------------------------------------------------------------------
b <- readLines(file.path(LOGS, "S2_md5_before.txt")); a <- readLines(file.path(LOGS, "S2_md5_after.txt"))
q0_ok <- grepl("^ *28800 OK$", b[1]) && grepl("^ *28800 OK$", a[1]) && length(grep("FAILED", c(a, b))) == 0

# ---- Q1 ------------------------------------------------------------------------
q11 <- loggrep("S2_a_b.log", "^Q1.1")
m11 <- regmatches(q11, regexpr("matched [0-9]+, mismatched [0-9]+", q11))
q11_ok <- all(grepl("mismatched 0; absent in O with n_U0 > 0: 0", q11))
q12 <- loggrep("S2_a_b.log", "^Q1.2 usable rows")
q12_ok <- grepl("'other': 0$", q12)
P <- rd("V4_S2_usable_pooled.csv")

# ---- Q2 ------------------------------------------------------------------------
G <- rd("V4_S2_replicates_by_cell.csv")
DP <- G[, .(m = min(n_gen)), by = .(dgp, model, N, ntime)]
n_dp_gen <- sum(DP$m < 1000); min_gen <- min(DP$m)
q21_man <- loggrep("S2_b_check.log", "^manuscript definition")
n_dp_man <- num_after(q21_man, "at least one arm:"); min_man <- num_after(q21_man, "min")
q22 <- loggrep("S2_a_b.log", "^Q2.2 generated sim_ids with fewer")
n_part <- num_after(q22, "techniques:")
q22b <- loggrep("S2_b_check.log", "^Q2.2 partial ids that are also missing")
q22c <- loggrep("S2_b_rerun.log", "^Q2.2 partial ids re-run")
M <- rd("V4_S2_missing_rerun.csv")
NB <- fread(file.path(OUTDIR, "missing_neighbours_rerun.csv"))
LD <- rd("V4_S2_lost_replicates_design.csv")

# ---- Q3 ------------------------------------------------------------------------
TV <- rd("V4_S2_truth_values.csv"); NT <- rd("V4_S2_truth_glmmnative.csv"); AS <- rd("V4_S2_truth_ar1_conditional.csv")
q31 <- loggrep("S2_c_truth.log", "^Q3.1"); q33 <- loggrep("S2_c_truth.log", "^Q3.3"); q34 <- loggrep("S2_c_truth.log", "^Q3.4")
FA <- fread(file.path(OUTDIR, "ar1_F_truth_reproduction.csv"))
AP <- fread(file.path(OUTDIR, "ar1_population_gee_limit_exploratory.csv"))

# ---- Q4 ------------------------------------------------------------------------
S <- rd("V4_S2_series_length.csv"); CK <- rd("V4_S2_series_length_check.csv"); CK[gen_error == "", gen_error := NA_character_]
pf <- S[, sum(p_full * n_subj) / sum(n_subj)]; mc_pf <- sqrt(0.21 / S[, sum(n_subj)])
T2 <- S[ntime == 2]; S[, ratio := mean_obs / ntime]
ck_ev <- CK[dgp != "across" & is.na(gen_error)]; ck_ac <- CK[dgp == "across"]

# ---- Q5 ------------------------------------------------------------------------
C1 <- rd("V4_S2_corr_arms.csv"); C2 <- rd("V4_S2_corr_mechanisms.csv")
k <- C1[method %in% c("mmm", "gee", "glmer") & n_pairs >= 100]
med <- function(pr, mt) C2[pair == pr & model <= 4 & method == mt & n_pairs >= 100, median(rho)]
medall <- function(pr, mt) C2[pair == pr & method == mt & n_pairs >= 100, median(rho)]

# ---- Q6 ------------------------------------------------------------------------
PAR <- rd("V4_S2_pilot_parity.csv"); REP <- rd("V4_S2_pilot_reproduction.csv"); PH <- rd("V4_S2_pilot_h0.csv")
run1 <- loggrep("S2_f_pilot_run1.log", "^Q6.1 parity")
q63 <- loggrep("S2_f_pilot_analysis.log", "^Q6.3 design identical")
PH34 <- PH[model %in% 3:4]; PH5 <- PH[model == 5]
zero_mcse <- PH34[!is.na(pass_t1) & mcse_d_t1 == 0, .N]
PHz <- PH[!is.na(mcse_d_t1) & mcse_d_t1 > 0][, z := d_t1_paired / mcse_d_t1]

# ---- Q7 ------------------------------------------------------------------------
LU <- rd("V4_S2_lus_vs_dgm.csv"); g <- function(it, col = "value") LU[item == it][[col]]

# ---- expectation table ------------------------------------------------------------
E <- rbindlist(list(
  list("Q0", "md5 di results_full/ prima e dopo", "28.800/28.800 OK, zero falliti",
       sprintf("prima: %s; dopo: %s", trimws(b[1]), trimws(a[1])), if (q0_ok) "si" else "no"),
  list("Q1.1", "conteggi usabili = colonna n di O_level_all (H0) e O_coverage_all (H1)", "zero discordanze",
       paste(m11, collapse = "; "), if (q11_ok) "si" else "no"),
  list("Q1.2", "status non NA classificati", "100%, classe altro vuota", sub("^Q1.2 ", "", q12), if (q12_ok) "si" else "no"),
  list("Q2.1", "design point con < 1.000 replicati generati in almeno un arm", "132 su 1.440, minimo 971",
       sprintf("%d su 1.440 (solo GLMM-AR1), minimo %d; con la definizione del manoscritto (almeno una riga time:trt) %d, minimo %d",
               n_dp_gen, min_gen, n_dp_man, min_man), if (n_dp_gen == 132 && min_gen == 971) "si" else "no"),
  list("Q2.2", "replicati generati con meno di 5 tecniche", "0", sprintf("%d", n_part), if (n_part == 0) "si" else "no"),
  list("Q2.3", "mancanti rieseguiti riproducono un errore di generazione con messaggio", "100%; vicini senza errore",
       sprintf("%d/%d (%s); vicini %d/%d senza errore", sum(!M$ok), nrow(M), paste(unique(M$msg), collapse = "; "), sum(NB$ok), nrow(NB)),
       if (all(!M$ok) && all(NB$ok)) "si" else "no"),
  list("Q3.1", "verita' marginale GenBinaryY = beta del generatore", "esatta", sub("^Q3.1 ", "", q31), if (grepl("10/10", q31)) "si" else "no"),
  list("Q3.2", "verita' condizionale GLMM-native = beta_c", "esatta",
       paste(unique(NT[, sprintf("M%d %s %s", model, label, f3(truth_conditional))]), collapse = "; "),
       if (all(NT[label == "mmlb", truth_conditional == fifelse(model == 4, 0.058, 0.232)]) && all(NT[label == "manual", truth_conditional == 0])) "si" else "no"),
  list("Q3.3", "quadratura riproduce F_truth_glmmnative.csv", "<= 1e-10", sub(".*F_truth_glmmnative = ([0-9.e+-]+).*", "max |diff| = \\1", q33),
       if (grepl("TRUE\\)$", q33)) "si" else "no"),
  list("Q3.4", "errore di quadratura 40 contro 80 nodi", "<= 1e-10", sub(".*nodes\\| = ([0-9.e+-]+).*", "\\1", q34), if (grepl("TRUE\\)$", q34)) "si" else "no"),
  list("Q3.5a", "verita' condizionale AR1 costante fra replicati", "sd < 1e-8 in ogni (M, arm, T, N)",
       sprintf("%d su %d; sd massima %s", sum(AS$sd < 1e-8, na.rm = TRUE), nrow(AS), f3(max(AS$sd, na.rm = TRUE), 4)),
       if (all(AS$sd < 1e-8, na.rm = TRUE)) "si" else "no"),
  list("Q3.5b", "riproduzione di F_truth_ar1.csv", "<= 1e-8",
       paste(FA[, .(d = max(d)), by = model][, sprintf("M%d %s", model, fg(d))], collapse = "; "), if (max(FA$d) <= 1e-8) "si" else "no"),
  list("Q4.1", "quota di serie piene aggregata", "0,70 +/- 0,01", sprintf("%s (MCSE %s)", f3(pf, 5), f3(mc_pf, 5)), if (abs(pf - 0.7) <= 0.01) "si" else "no"),
  list("Q4.2", "quota per cella", "0,70 +/- 0,01; nessun superamento a N >= 25",
       sprintf("%d superamenti su %d; |z| massimo %s", sum(!S$pass_q42), nrow(S), f3(max(abs(S$z_full)), 2)), if (all(S$pass_q42)) "si" else "no"),
  list("Q4.3", "T = 2: quota con una osservazione; media/T", "0,30 +/- 0,01; 0,85 +/- 0,02",
       sprintf("%s a %s; %s a %s", f3(min(T2$p_single), 4), f3(max(T2$p_single), 4), f3(min(S$ratio), 4), f3(max(S$ratio), 4)),
       if (all(abs(T2$p_single - 0.3) <= 0.01) && all(abs(S$ratio - 0.85) <= 0.02)) "si" else "no"),
  list("Q4.4", "nclust rigenerato = make_all_data; nclust e trt uguali fra meccanismi", "150/150",
       sprintf("%d/%d valutabili identici (%d non valutabili: generazione AR1 fallita); fra meccanismi %d/%d",
               sum(ck_ev$nclust_identical), nrow(ck_ev), CK[dgp != "across" & !is.na(gen_error), .N], sum(ck_ac$nclust_identical, na.rm = TRUE), nrow(ck_ac)),
       if (nrow(ck_ev) == 150 && all(ck_ev$nclust_identical)) "si" else sprintf("no alla lettera: %d/150 valutabili", nrow(ck_ev))),
  list("Q5.1", "rho fra arm > 0,5 (mmm, gee, glmer, >= 100 coppie); IC > 0", ">= 95%; tutte",
       sprintf("%d/%d (%s%%); IC > 0 in %d/%d", sum(k$rho > 0.5), nrow(k), f3(100 * mean(k$rho > 0.5), 1), sum(k$lo > 0), nrow(k)),
       if (mean(k$rho > 0.5) >= 0.95 && all(k$lo > 0)) "si" else "no"),
  list("Q5.2", "mediana rho native-AR1 > GenBinaryY-native (M1-M4); |mediana GB-native| < 0,3", "gee e glmer",
       sprintf("gee %s contro %s; glmer %s contro %s; |GB-native| %s e %s", f3(med("native_vs_AR1", "gee")), f3(med("GenBinaryY_vs_native", "gee")),
               f3(med("native_vs_AR1", "glmer")), f3(med("GenBinaryY_vs_native", "glmer")), f3(abs(medall("GenBinaryY_vs_native", "gee"))), f3(abs(medall("GenBinaryY_vs_native", "glmer")))),
       if (med("native_vs_AR1", "gee") > med("GenBinaryY_vs_native", "gee") && med("native_vs_AR1", "glmer") > med("GenBinaryY_vs_native", "glmer")) "si" else "no"),
  list("Q6.1", "parita' H0 corretto / H1: differisce solo l'ultimo elemento di beta", "12/12",
       sprintf("primo controllo %s (strumento difettoso, vedi sez. 2); controllo corretto %d/12, controprova dello strumento %d/12",
               sub(".*: ", "", run1), sum(PAR$pass), sum(PAR$instr_ok)), if (all(PAR$pass) && all(PAR$instr_ok)) "si (al secondo controllo)" else "no"),
  list("Q6.2", "il banco riproduce results_full (i = 1-5)", "60/60", sprintf("%d/%d", sum(REP$identical), nrow(REP)), if (all(REP$identical)) "si" else "no"),
  list("Q6.3", "disegno identico fra H0 corretto e H0 originale", "100%", sub(".*: ", "", q63), if (grepl("12000/12000", q63)) "si" else "no"),
  list("Q6.4", "M3-M4: differenza appaiata di usabilita' e di errore di I tipo entro 3 MCSE", "circa 0,2 superamenti su circa 70",
       sprintf("usabilita' %d/%d; I tipo %d/%d (%d NA; %d confronti con MCSE appaiata 0)", sum(PH34$pass_usable, na.rm = TRUE), sum(!is.na(PH34$pass_usable)),
               sum(PH34$pass_t1, na.rm = TRUE), sum(!is.na(PH34$pass_t1)), sum(is.na(PH34$pass_t1)), zero_mcse),
       if (all(PH34$pass_usable, na.rm = TRUE) && all(PH34$pass_t1, na.rm = TRUE)) "si" else "no"),
  list("Q7.1", "riproduzione della coorte: 2.296 esiti; varianza paziente 3,545; contrasto 0,378 (0,168)", "+/- 0,001",
       sprintf("%d esiti; %s; MMM %s (SE %s)", g("n_zone_outcomes"), f3(g("var_id_modelA"), 4), f3(g("contrast_MMM_regionposterior"), 4), f3(g("contrast_MMM_regionposterior", "se"), 4)),
       if (g("n_zone_outcomes") == 2296 && abs(g("var_id_modelA") - 3.545) <= 0.001 && abs(g("contrast_MMM_regionposterior") - 0.378) <= 0.001 && abs(g("contrast_MMM_regionposterior", "se") - 0.168) <= 0.001) "si" else "no"),
  list("Q7.2", "DGM dentro l'IC 95% della coorte (prevalenza, varianze)", "coerente / scarto dichiarato",
       sprintf("prevalenza %s [%s, %s] contro %s: %s; varianza paziente (solo paziente) %s [%s, %s] contro 1: %s",
               f3(g("prevalence_lus_severe")), f3(g("prevalence_lus_severe", "lo")), f3(g("prevalence_lus_severe", "hi")), f3(plogis(-2)), g("prevalence_lus_severe", "note"),
               f3(g("var_id_modelA")), f3(g("var_id_modelA", "lo")), f3(g("var_id_modelA", "hi")), sub(";.*", "", g("var_id_modelA", "note"))), "criterio applicato")
))
setnames(E, c("id", "attesa", "tolleranza", "osservato", "superato"))
OUT <- c(OUT, "<!-- generated by R/s2/s2_tables.R; do not edit by hand -->", "", "### Tabella E. Attese pre-registrate e osservato", "", md(E))

# ---- detail tables ------------------------------------------------------------------
OUT <- c(OUT, "### Tabella A1. Quota usabile per metodo, sull'intera griglia (2.880 celle, 2.879.506 replicati generati)", "",
  "U0 = filtro di `agg_all5.R` (warning ammessi); U1 = U0 senza alcun warning; U2 = U0 senza warning di convergenza (lme4 o glm). Le combinazioni cella x metodo senza righe contano 0.", "",
  md(P[order(match(method, c("mmm", "gee", "gee_md", "gee_kc", "gee_fg", "glmer", "shrink", "glmer_lrt"))),
       .(method, U0 = f3(U0), U1 = f3(U1), U2 = f3(U2), `warning fra gli usabili` = f3(warn_among_U0), `mediana per cella` = f3(cell_med_warn),
         `celle > 10%` = cells_warn_gt10, `celle > 50%` = cells_warn_gt50, `celle con U0 = 0` = cells_zero_U0)]))
X <- rd("V4_S2_typeI_U0_vs_U1_summary.csv")
OUT <- c(OUT, "### Tabella A2. Esplorativa, non pre-registrata (per D-04): errore di I tipo con U0 e con U1, celle H0 con almeno 100 usabili in entrambe", "",
  md(X[, .(method, celle = cells, `I tipo U0` = f3(mean_t1_U0, 4), `I tipo U1` = f3(mean_t1_U1, 4), `differenza media` = f3(mean_d, 4),
           `5° pct` = f3(q05_d, 4), `95° pct` = f3(q95_d, 4), `celle con |diff| > 2 MCSE` = cells_absd_gt_2mcse)]))
R <- rd("V4_S2_unusable_reasons.csv")
OUT <- c(OUT, "### Tabella A3. Motivi di non usabilita' (righe time:trt o di fallimento)", "", md(dcast(R, method ~ reason, value.var = "N", fill = 0)))
miss <- rd("V4_S2_missing_ids.csv")
OUT <- c(OUT, "### Tabella B1. Replicati mancanti in GLMM-AR1 per scenario e T (entrambi gli arm, stessi sim_id)", "",
  md(dcast(miss[, .N, by = .(model, ntime)], model ~ ntime, value.var = "N", fill = 0)),
  "### Tabella B2. Disegno dei replicati mancanti (arm H1; l'arm H0 ha lo stesso disegno)", "",
  md(LD[, .N, by = .(tipo = type)]), md(dcast(LD[, .(cnt = .N), by = .(N, ntime)], N ~ ntime, value.var = "cnt", fill = 0)))
OUT <- c(OUT, "### Tabella C1. Verita' marginale di time:trt in GLMM-native (limite della GEE di indipendenza, disegno realizzato)", "",
  md(dcast(NT[, .(model, label, ntime, v = f3(truth_marginal_weighted, 4))], model + label ~ ntime, value.var = "v")),
  sprintf("Scarto massimo fra versione pesata e griglia bilanciata: %s.", f3(max(abs(NT$truth_marginal_weighted - NT$truth_marginal_balanced)), 5)), "",
  "### Tabella C2. Verita' condizionale di time:trt in GLMM-AR1 (media su 50 replicati a N = 145; sd massima su N in {10, 25, 75, 145})", "",
  md(merge(dcast(AS[N == 145, .(model, label, ntime, v = f3(mean, 4))], model + label ~ ntime, value.var = "v"),
           AS[, .(`sd massima` = f3(max(sd, na.rm = TRUE), 4)), by = .(model, label)], by = c("model", "label"))),
  sprintf("Esplorativo: il limite di popolazione della GEE di indipendenza implicato dal beta_c calibrato dista dal target al massimo %s.", f3(max(abs(AP$gee_limit_pop - AP$target)), 4)), "",
  "### Tabella C3. Riproduzione di F_truth_ar1.csv con la sua stessa procedura", "",
  md(FA[, .(`max |diff|` = fg(max(d)), righe = .N), by = model]))
OUT <- c(OUT, "### Tabella D1. Lunghezza delle serie per T (media su scenari e N)", "",
  md(S[, .(`quota piene` = f3(sum(p_full * n_subj) / sum(n_subj), 4), `quota con 1 osservazione` = f3(sum(p_single * n_subj) / sum(n_subj), 4),
           `osservazioni medie` = f3(sum(mean_obs * n_subj) / sum(n_subj), 3)), by = ntime]))
OUT <- c(OUT, "### Tabella E1. Correlazione delle stime fra arm (mediana [min, max] sulle celle con almeno 100 coppie)", "",
  md(dcast(C1[n_pairs >= 100, .(v = sprintf("%s [%s, %s]", f3(median(rho)), f3(min(rho)), f3(max(rho)))), by = .(dgp, method)], method ~ dgp, value.var = "v")),
  "### Tabella E2. Correlazione delle stime fra meccanismi (mediana sulle celle con almeno 100 coppie)", "",
  md(dcast(C2[n_pairs >= 100, .(v = f3(median(rho))), by = .(pair, method)], method ~ pair, value.var = "v")))
OUT <- c(OUT, "### Tabella F1. Pilota H0 corretto: errore di I tipo e quota usabile, originale e corretto", "",
  md(PH[, .(M = model, N, T = ntime, method, `usabile orig.` = f3(usable_old), `usabile corr.` = f3(usable_new),
            `I tipo orig.` = f3(t1_old), `I tipo corr.` = f3(t1_new), `diff. appaiata` = f3(d_t1_paired), MCSE = f3(mcse_d_t1),
            entro = fifelse(is.na(pass_t1), "NA", fifelse(pass_t1, "si", "no")))]),
  sprintf("|z| massimo della differenza appaiata: %s (M%d, N = %d, T = %d, %s).", f3(max(abs(PHz$z)), 2),
          PHz[which.max(abs(z))]$model, PHz[which.max(abs(z))]$N, PHz[which.max(abs(z))]$ntime, PHz[which.max(abs(z))]$method), "")
LV <- LU[grepl("^var_|^icc_|^prevalence|^dgm_marginal|^LRT", item)]
OUT <- c(OUT, "### Tabella G1. Coorte LUS contro DGM: prevalenza e struttura di varianza", "",
  "Modello A: `lus_severe ~ region + (1|id)`; B: + `(1|pt)`; C: + `(1|hemi)`. IC 95% di profilo sulla SD, al quadrato; prevalenza con bootstrap per paziente (B = 2000, seed 20261015).", "",
  md(LV[, .(voce = item, valore = f3(value, 4), lo = fifelse(is.na(lo), "", f3(lo, 4)), hi = fifelse(is.na(hi), "", f3(hi, 4)), DGM = fifelse(is.na(dgm), "", f3(dgm, 4)), nota = note)]))
LC <- LU[grepl("^contrast_", item)]
LC[, c("metodo", "termine") := tstrsplit(sub("^contrast_", "", item), "_", keep = 1:2)]
LC[, termine := sub("^[A-Z]+_", "", sub("^contrast_", "", item))]
OUT <- c(OUT, "### Tabella G2. Contrasti della coorte sulla scala logit (riferimento di grandezza, D-10b); DGM: 0,2 per occasione", "",
  md(dcast(LC[, .(termine, metodo, v = sprintf("%s (%s)", f3(value), f3(se)))], termine ~ metodo, value.var = "v")))
writeLines(OUT, file.path(REPO, "reports/morris/S2_results.md"), useBytes = TRUE)
fwrite(E, file.path(TABDIR, "V4_S2_expectations.csv"))
cat("wrote S2_results.md,", length(OUT), "lines\n"); print(E[, .(id, superato)])
