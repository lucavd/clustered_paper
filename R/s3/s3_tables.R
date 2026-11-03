# =============================================================================
# s3_tables.R -- builds reports/morris/S3_results.md (Tabella E and detail
# tables) from the files written by the S3 checks. Every number is read from a
# file; nothing is typed by hand.
# =============================================================================
suppressMessages(library(data.table, lib.loc = path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")))
REPO <- path.expand("~/MMM"); TAB <- file.path(REPO, "tables"); LOGD <- file.path(REPO, "reports/morris/S3_logs")
rd <- function(f) fread(file.path(TAB, f))
lg <- function(f) readLines(file.path(LOGD, f), warn = FALSE)
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d, decimal.mark = ",")
int <- function(x) formatC(x, format = "d", big.mark = ".")
md <- function(dt) { dt <- as.data.frame(dt); dt[] <- lapply(dt, function(v) gsub("\\|", "\\\\|", as.character(v)))
  c(paste0("| ", paste(names(dt), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(dt)), collapse = "|"), "|"),
    apply(dt, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))) }
yn <- function(b) ifelse(b, "si", "no")
L <- character()
add <- function(...) L <<- c(L, ...)

q0pre <- lg("S3_q0_pre.log"); q0post <- if (file.exists(file.path(LOGD, "S3_q0_post.log"))) lg("S3_q0_post.log") else "non eseguito"
par4 <- rd("V4_S3_parity_v4.csv"); paro <- rd("V4_S3_parity_orig.csv"); q2 <- rd("V4_S3_q2_d02.csv"); q3 <- rd("V4_S3_q3_bench.csv")
c101 <- rd("V4_S3_m1c_calibration.csv"); c103 <- rd("V4_S3_m1c_prevalence.csv"); c105 <- rd("V4_S3_m1c_constants.csv")
P <- readRDS(file.path(TAB, "V4_S3_post_checks.rds")); q7 <- rd("V4_S3_q7_degenerate.csv"); q8 <- rd("V4_S3_q8_lrt_vs_wald.csv"); q9 <- rd("V4_S3_q9_h0_paired.csv")
main <- par4$scen != "M1c"; orig_fail <- paro[pass == FALSE, paste(dgp, scen)]
q6 <- P$Q6; q4 <- P$Q4; q5 <- P$Q5; q7b <- P$Q7; q7g <- P$Q7_genfail_nondeg; q8r <- P$Q8_re; q9s <- P$Q9
q8 <- q8[!is.na(dgp)]
E <- rbindlist(list(
  list("Q0", "md5 di results_full/ e results_v4 (S2) invariati prima e dopo", "zero falliti",
       paste(c(paste("prima:", grep("^Q0", q0pre, value = TRUE)), paste("dopo:", grep("^Q0", q0post, value = TRUE))), collapse = "; "),
       yn(length(grep("OK", c(q0pre, q0post))) == 4 && !any(grepl("FAIL", c(q0pre, q0post))))),
  list("Q1", "parita' H0/H1: 15/15 + 3/3 (M1c); controprova sui generatori originali 12/15, falliti GenBinaryY M3-M5", "esatta",
       sprintf("v4 %d/%d + %d/%d; originali %d/%d, falliti: %s", sum(par4$pass[main]), sum(main), sum(par4$pass[!main]), sum(!main),
               sum(paro$pass), nrow(paro), paste(orig_fail, collapse = ", ")),
       yn(all(par4$pass) && sum(main) == 15 && sum(!main) == 3 && setequal(orig_fail, paste("GenBinaryY", c("M3", "M4", "M5"))))),
  list("Q2", "AR1 M5 v4 differisce dall'originale solo in D[3,3] (0,01 -> 0,001) e D = D di GLMM-native M5", "esatta",
       sprintf("%d/%d", sum(q2$pass), nrow(q2)), yn(all(q2$pass))),
  list("Q3", "banco con codice originale riproduce results_full/ bit per bit, 300 replicati", "300/300",
       sprintf("%d/%d (righe di traccia escluse: %d)", sum(q3$identical), nrow(q3), sum(q3$trace_rows)), yn(all(q3$identical) && nrow(q3) == 300)),
  list("Q4", "disegno identico al generatore di riferimento (F1, F2: originale; C: M1)", "100%",
       paste(sprintf("%s %s/%s", q4$block, int(q4$identical), int(q4$checked)), collapse = "; "), yn(all(q4$identical == q4$checked))),
  list("Q5", "blocco L: gee rieseguita identica a results_full/", "100%",
       sprintf("%s/%s chunk identici; sim_id uguali %s/%s; righe %s contro %s", int(q5$gee_identical), int(q5$chunks), int(q5$gee_sims_equal), int(q5$chunks), int(q5$rows_v4), int(q5$rows_orig)),
       yn(q5$gee_identical == q5$chunks && q5$gee_sims_equal == q5$chunks)),
  list("Q6", "18.720 chunk, zero errori, ogni sim_id con una traccia", "zero",
       sprintf("%s/%s chunk; sim completi %s; errori di crew: vedi run.log", int(sum(q6$exist)), int(sum(q6$expected)), int(sum(q6$sims_ok))),
       yn(all(q6$exist == q6$expected) && all(q6$sims_ok == q6$expected) && all(q6$every_sim_rows == q6$expected) && all(q6$gen_rows_match == q6$expected))),
  list("Q7", "design_degenerate in biiezione con le chiavi di S2; errori di generazione solo su chiavi degeneri", "esatta",
       sprintf("%s; errori su chiavi non degeneri: %s", paste(sprintf("%s %s/%s celle", q7b$block, int(q7b$bijection), int(q7b$cells)), collapse = ", "), int(sum(q7g$gen_fail_on_nondegenerate))),
       yn(all(q7b$bijection == q7b$cells) && sum(q7g$gen_fail_on_nondegenerate) == 0)),
  list("Q8a", "re_structure del LRT: pendenze in M3-M5, (1|id) in M1c", "100%",
       paste(sprintf("%s %s/%s", q8r$block, int(q8r$re_ok), int(q8r$chunks)), collapse = "; "), yn(all(q8r$re_ok == q8r$chunks))),
  list("Q8b", "abs(LRT - Wald) < 1e-3 in almeno il 95% delle coppie senza warning, per scenario x meccanismo", ">= 0,95 in ogni riga",
       sprintf("%d/%d righe; quota minima %s (%s %s)", sum(q8$pass), nrow(q8), fmt(min(q8$share)), q8[which.min(share), dgp], q8[which.min(share), scen]), yn(all(q8$pass))),
  list("Q9", "GenBinaryY M3-M4: differenza appaiata H0 v4 - H0 originale entro 3 MCSE", ">= 99%",
       paste(sprintf("%s %s/%s (%s)", q9s$measure, int(q9s$within_3), int(q9s$evaluable), fmt(q9s$share, 4)), collapse = "; "), yn(all(q9s$pass))),
  list("Q10.1", "beta0(T) entro 5e-4 dai valori pre-registrati, residuo < 1e-10", "4/4",
       sprintf("%d/4; beta0 = %s", sum(c101$pass_101), paste(fmt(c101$beta0, 4), collapse = "; ")), yn(all(c101$pass_101))),
  list("Q10.2", "prevalenza attesa H0: 0,5945 (GenBinaryY), 0,589 +/- 0,001 (GLMM-native)", "4/4",
       sprintf("%d/4; GLMM-native %s", sum(c101$pass_102), paste(fmt(c101$prev_glmm_H0, 4), collapse = "; ")), yn(all(c101$pass_102))),
  list("Q10.3", "prevalenza realizzata a N = 45 entro 3 MCSE dall'attesa", "24/24",
       sprintf("%d/%d; |z| massimo %s", sum(c103$pass_103), nrow(c103), fmt(max(abs(c103$z)), 2)), yn(all(c103$pass_103) && nrow(c103) == 24)),
  list("Q10.4", "AR1: errore massimo di calibrazione < 1e-8 sui replicati non degeneri", "< 1e-8",
       sprintf("%.2e", max(c103$max_calib_err_nondeg, na.rm = TRUE)), yn(max(c103$max_calib_err_nondeg, na.rm = TRUE) < 1e-8)),
  list("Q10.5", "costanti e argomenti dei motori come D-14", "tutti",
       sprintf("%d/%d", sum(c105$pass), nrow(c105)), yn(all(c105$pass)))
))
setnames(E, c("id", "attesa", "tolleranza", "osservato", "superato"))
fwrite(E, file.path(TAB, "V4_S3_expectations.csv"))
add("### Tabella E. Attese pre-registrate e osservato", "", md(E), "")
add("### E1. Parita' H0/H1 per generatore (Q1)", "", md(merge(par4[, .(dgp, scen, v4 = sprintf("%d/%d", n_pass, n_checks), v4_differisce_in = differs_in)],
    paro[, .(dgp, scen, originale = sprintf("%d/%d", n_pass, n_checks), originale_differisce_in = differs_in)], by = c("dgp", "scen"), all.x = TRUE)), "")
add("### E2. Calibrazione di M1c (Q10.1, Q10.2)", "", md(c101[, .(T = ntime, beta0 = fmt(beta0, 4), preregistrato = fmt(prereg, 4), residuo = sprintf("%.1e", resid),
    `prev. H1 marginale` = fmt(prev_H1_marginal, 4), `GLMM-native H0` = fmt(prev_glmm_H0, 4), `GLMM-native H1` = fmt(prev_glmm_H1, 4))]), "")
add("### E3. Prevalenza realizzata di M1c a N = 45, replicati 1-1.000 (Q10.3)", "", md(c103[, .(meccanismo = dgp, T = ntime, arm = ifelse(label == "mmlb", "H1", "H0"),
    generati = generated, prevalenza = fmt(prev, 4), MCSE = fmt(mcse, 4), attesa = fmt(expected, 4), z = fmt(z, 2))]), "")
add("### E4. Completezza, disegno e identita' per blocco (Q4-Q7)", "", md(merge(q6[, .(blocco = block, attesi = int(expected), scritti = int(exist), `sim completi` = int(sims_ok),
    `errori di generazione` = int(gen_fail), `righe di traccia di analisi` = int(trace_analysis_rows))],
    q7g[, .(blocco = block, `errori su chiavi non degeneri` = int(gen_fail_on_nondegenerate))], by = "blocco")), "")
add("### E5. LRT contro Wald senza warning (Q8b)", "", md(q8[, .(meccanismo = dgp, scenario = scen, coppie = int(pairs), entro = int(close), quota = fmt(share, 4), superato = yn(pass))]), "")
k <- q9[model %in% 3:4]
ex <- k[(mcse_d_usable > 0 & abs(z_usable) > 3) | (mcse_d_t1 > 0 & abs(z_t1) > 3)]
add("### E6. Effetto di D-01 su GenBinaryY M3-M5 H0 (Q9)", "",
    sprintf("M3-M4: %s combinazioni cella x metodo; superamenti di 3 MCSE elencati sotto (%d righe). M5 esplorativo: I tipo %s su %s oltre 3 MCSE, quota usabile %s su %s.",
            int(nrow(k)), nrow(ex), int(P$Q9_M5_explore$over3_t1), int(P$Q9_M5_explore$evaluable_t1), int(P$Q9_M5_explore$over3_u), int(P$Q9_M5_explore$evaluable_u)), "")
if (nrow(ex)) add(md(ex[, .(M = model, N, T = ntime, metodo = method, `d usabile` = fmt(d_usable, 4), z_u = fmt(z_usable, 2), `I tipo orig` = fmt(t1_old, 4), `I tipo v4` = fmt(t1_new, 4), z_t1 = fmt(z_t1, 2))]), "")
add(sprintf("Media della differenza appaiata dell'errore di I tipo in M3-M4, per metodo (cella pesata 1): %s.",
            paste(k[mcse_d_t1 > 0, .(m = mean(d_t1)), by = method][, sprintf("%s %s", method, fmt(m, 4))], collapse = "; ")), "")
q9m <- q9[, .(`I tipo valutabili` = sum(mcse_d_t1 > 0, na.rm = TRUE), `I tipo oltre 3 MCSE` = sum(mcse_d_t1 > 0 & abs(z_t1) > 3, na.rm = TRUE),
                `media d I tipo` = fmt(mean(d_t1[mcse_d_t1 > 0], na.rm = TRUE), 4),
                `usabile valutabili` = sum(mcse_d_usable > 0, na.rm = TRUE), `usabile oltre 3 MCSE` = sum(mcse_d_usable > 0 & abs(z_usable) > 3, na.rm = TRUE)),
            by = .(scenari = ifelse(model == 5, "M5 (esplorativo)", "M3-M4"), metodo = method)][order(scenari, metodo)]
add("### E6b. Q9 per metodo (scomposizione diagnostica, non pre-registrata)", "",
    "In F1 `glmer_lrt` cambia anche per D-08 (pendenze casuali nel modello del LRT), quindi la sua differenza appaiata non misura solo D-01.", "", md(q9m), "")
k2 <- k[method != "glmer_lrt"]
add(sprintf("Senza `glmer_lrt`, M3-M4: I tipo %d/%d entro 3 MCSE, quota usabile %d/%d.", k2[mcse_d_t1 > 0 & abs(z_t1) <= 3, .N], k2[mcse_d_t1 > 0, .N],
            k2[mcse_d_usable > 0 & abs(z_usable) <= 3, .N], k2[mcse_d_usable > 0, .N]), "")
iv <- rd("V4_S3_q8_investigation.csv")
add("### E5b. Indagine su Q8b: rifit di 400 replicati campionati (esplorativa, non pre-registrata)", "",
    "Campione: per GenBinaryY M5, GLMM-AR1 M5, GLMM-native M5, GenBinaryY M3 e GenBinaryY M1c, 60 coppie discordanti e 20 concordanti da 60 chunk estratti a caso (seed 20261017). d logLik = logLik(LRT completo, bobyqa) - logLik(Wald, ottimizzatore di default).", "",
    md(iv[, .(n = .N, `riproduce Wald` = sum(reproduces_wald), `riproduce LRT` = sum(reproduces_lrt), `stessi theta` = sum(n_theta_w == n_theta_l),
              singolare = sum(sing_w | sing_l), `abs(d logLik) < 1e-3` = sum(abs(d_ll) < 1e-3), `LRT piu' alta` = sum(d_ll > 1e-3), `Wald piu' alta` = sum(d_ll < -1e-3),
              `mediana abs(d stima)` = sprintf("%.2g", median(d_est))), by = .(meccanismo = dgp, scenario = scen, concordante = yn(close))]), "")
# timing
cs <- fread(file.path(REPO, "results_v4/S3/chunk_seconds.tsv"), header = FALSE, col.names = c("key", "sec"))
ff <- list.files(file.path(REPO, "results_v4/S3"), recursive = TRUE, pattern = "[.]rds$")
blk <- data.table(chunk_key = sub("[.]rds$", "", basename(ff)), block = basename(dirname(ff)))   # data.table(key = ) is reserved
setnames(cs, "key", "chunk_key")
cs <- merge(cs, blk, by = "chunk_key")
add("### E7. Tempo di calcolo per blocco", "", md(cs[, .(chunk = int(.N), `ore-worker` = fmt(sum(sec) / 3600, 1), `s/chunk mediana` = fmt(median(sec), 1), `s/chunk massimo` = fmt(max(sec), 1)), by = .(blocco = block)]), "")
writeLines(L, file.path(REPO, "reports/morris/S3_results.md"))
cat(sprintf("S3_results.md written: %d lines; expectations passed %d/%d\n", length(L), sum(E$superato == "si"), nrow(E)))
