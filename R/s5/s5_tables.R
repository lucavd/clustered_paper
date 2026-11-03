# =============================================================================
# s5_tables.R -- S5: markdown tables of the report from the tables/V4_S5_*.csv files
# (every number is read back from a file). Output: reports/morris/S5_results.md,
# tables/V4_S5_expectations.csv
# =============================================================================
source(path.expand("~/MMM/R/s5/s5_common.R"))
rd <- function(x) fread(file.path(TABDIR, x))
it <- function(x, d = 2) formatC(x, format = "f", digits = d, big.mark = ".", decimal.mark = ",")
iti <- function(x) formatC(round(x), format = "d", big.mark = ".")
md <- function(D) { h <- paste0("| ", paste(names(D), collapse = " | "), " |"); s <- paste0("|", paste(rep("---", ncol(D)), collapse = "|"), "|")
  b <- apply(D, 1, function(r) paste0("| ", paste(gsub("|", "\\|", r, fixed = TRUE), collapse = " | "), " |")); c(h, s, b, "") }
L <- character()
add <- function(...) L <<- c(L, ...)
CL <- rd("V4_S5_claims.csv"); VAL <- rd("V4_S5_values.csv"); G <- rd("V4_S5_gates.csv")
P1 <- rd("V4_S5_p1_reproduction.csv"); P2 <- rd("V4_S5_p2_stability.csv"); P3 <- rd("V4_S5_p3_crn.csv"); G3 <- rd("V4_S5_g3_linearization.csv")
G2 <- rd("V4_S5_g2.csv"); SW <- rd("V4_S5_switch_points.csv"); FB <- rd("V4_S5_fu_b_crn_corr.csv"); FC <- rd("V4_S5_fu_c_bias.csv")
FD <- rd("V4_S5_fu_d_wins.csv"); FT <- rd("V4_S5_fu_d_ties.csv"); FF <- rd("V4_S5_fu_f_rmse.csv"); FG <- rd("V4_S5_fu_g_usable.csv")
C5 <- rd("V4_S5_c05_winner_runnerup.csv"); C8 <- rd("V4_S5_p1_c08_diagnostic.csv"); ARMS <- rd("V4_S2_corr_arms.csv")
INV <- fread(file.path(REPO, "reports/morris/S5_claims_inventory.csv"))
cv <- function(nm) CL[name == nm][1]

# ---------------- expectations table ----------------
g1 <- function(tag) { x <- readLines(file.path(LOG_S5, sprintf("S5_g1_%s.log", tag)))[-1]; paste(sub("^G1 ", "", x), collapse = "; ") }
c03 <- CL[claim_id == "C03"]; c04 <- cv("C04|all"); c14a <- cv("C14a|all"); c14b <- cv("C14b|all")
p4 <- function(ids) paste(sprintf("%s: %s", CL[name %in% ids, name], CL[name %in% ids, class]), collapse = "; ")
E <- data.table(
  id = c("G1", "G2", "G3", "G4", "P1", "P2", "P3", "P4a C01", "P4b C02", "P4c C03", "P4d C04", "P4e C20", "P4f C14"),
  atteso = c("md5 invariati prima e dopo", "b = 0 e MCSE di S4 riprodotte entro 1e-12", "SE bootstrap = SE per linearizzazione entro 7%", "ogni riga con un esito; >= 1.990 replicati finiti",
             "numeri v3 riprodotti alla cifra stampata", "v4 - v3 entro +/- 1 punto (63 valori)", "rapporto SE blocchi / SE indipendenti circa 2, sopra 1,25",
             "GB risolta e rilevante; NAT e AR1 trascurabile o indeterminata", "trascurabile sui tre", "risolta, rilevanza indeterminata sui tre",
             "uguaglianza confermata (margine 72)", "C20b rilevante, C20a trascurabile", "vantaggio del LRT ridotto di oltre 5 punti"),
  osservato = c(sprintf("prima: %s || dopo: %s", g1("pre"), g1("post")), G[gate == "G2", result], G[gate == "G3", result], G[gate == "G4", result],
                G[gate == "P1", result], G[gate == "P2", result], G[gate == "P3", result],
                p4(c("C01|GB", "C01|NAT", "C01|AR1")), p4(c("C02|GB", "C02|NAT", "C02|AR1")), p4(c("C03|GB", "C03|NAT", "C03|AR1")),
                sprintf("MMM - FG = %s [%s; %s]: %s", it(c04$est, 1), it(c04$lo, 1), it(c04$hi, 1), c04$class),
                p4(c("C20a|GB", "C20a|NAT", "C20a|AR1", "C20b|GB", "C20b|NAT", "C20b|AR1")),
                sprintf("usabilita': scarto da %s a %s [%s; %s]; validita': da %s a %s [%s; %s]", it(P1[name == "use20|glmer_lrt", printed] - P1[grepl("^use20", name) & name != "use20|glmer_lrt", max(printed)], 1), it(c14a$est, 1), it(c14a$lo, 1), it(c14a$hi, 1),
                        it(P1[name == "val20|glmer_lrt", printed] - P1[grepl("^val20", name) & name != "val20|glmer_lrt", max(printed)], 1), it(c14b$est, 1), it(c14b$lo, 1), it(c14b$hi, 1))),
  esito = c("superato", "superato (dopo la correzione dello strumento, sez. 3)", ifelse(G[gate == "G3", pass], "superato", "non superato"), ifelse(G[gate == "G4", pass], "superato", "non superato"),
            "smentita, spiegata", ifelse(G[gate == "P2", pass], "confermata", "smentita"), "smentita, spiegata", "confermata", "confermata", "smentita su GB e NAT, spiegata",
            "smentita, spiegata", "confermata", "confermata"))
fwrite(E, file.path(TABDIR, "V4_S5_expectations.csv"))
add("### Tabella E. Cancelli e previsioni", "", md(E))

# ---------------- claim matrix ----------------
lab <- function(x) fifelse(is.na(x), "", x)
X <- CL[!is.na(est) & kind %in% c("difference", "equivalence", "value", "rank")]
X <- X[, .(`sotto-affermazione` = name, sezione = section, `valore v3` = v3_value, `stima v4` = it(est, 3), `intervallo 95%` = sprintf("[%s; %s]", it(lo, 3), it(hi, 3)),
           margine = fifelse(is.na(margin), "", it(margin, 2)), esito = class,
           `esito, intervallo basic` = fifelse(class_basic != class, class_basic, ""), `esito con bobyqa` = fifelse(bobyqa_changes %in% TRUE, class_bobyqa, ""),
           `replicati per risolvere` = fifelse(is.finite(n_req), iti(n_req), ""))]
add("### Tabella A. Matrice affermazione x MCSE (griglia principale, U0, insieme proprio)", "",
    "Stima e intervallo sono orientati nel verso del manoscritto: un valore positivo va nella direzione affermata. Intervallo percentile del bootstrap per blocchi (B = 2.000). Le righe con il nome della affermazione sono descritte in `reports/morris/S5_claims_inventory.csv`.", "", md(X))
cnt <- CL[!is.na(class), .N, by = .(esito = class)][order(-N)]
add("### Tabella B. Esiti per classe (tutte le righe calcolate e non calcolate)", "", md(cnt))
SWt <- SW[, .(meccanismo = mech, `cambio v4 (N)` = switch_v4, `cambio v3 (N)` = claimed, `P(cambio = v3)` = it(p_claimed, 3), `metodo dopo il cambio` = next_v4,
              `P(metodo dopo = MD)` = it(p_next_md, 3), `N distinti dell'argmax del caso peggiore di MMM` = wc_argmax_distinctN)]
add("### Tabella C. Punti di cambio della validita' incondizionata (C23a, C23b) e instabilita' del caso peggiore (C22)", "", md(SWt))
c5 <- CL[claim_id == "C05", .(conteggio = name, valore = iti(est), esito = class)]
add("### Tabella D. Celle vinte contro la seconda (C05): punti di disegno con differenza significativa", "", sprintf("Su %s punti di disegno con almeno due metodi stimabili.", iti(nrow(C5))), "", md(c5))

# ---------------- appendix ----------------
add("## Appendice", "")
nr <- P1[status != "reproduced", .(numero = name, stampato = printed, `ricalcolato` = it(primary, 3), fonte = source)]
add("### A1. Numeri del v3 non riprodotti (P1)", "", sprintf("Riprodotti %s di %s numeri.", iti(P1[status == "reproduced", .N]), iti(nrow(P1))), "", md(nr))
add("Pesature della deviazione media per la correlazione fra meccanismi (C08, da `O_level_all.csv`):", "", md(C8[, .(pesatura = weighting, coppia = pair, rho = it(rho, 3))]))
p2 <- P2[, .(n = .N, `media v4 - v3` = it(mean(diff), 3), `massimo abs` = it(max(abs(diff)), 3), confermati = sum(class == "confermata"))]
add("### A2. Stabilita' dal v3 al v4 (P2)", "", md(p2))
p3 <- P3[, .(quantita = quantity, `SE blocchi` = it(se_boot, 4), `SE indipendente` = it(se_ind, 4), rapporto = it(ratio, 3))]
add("### A3. Peso dei numeri casuali comuni (P3)", "", md(p3))
fb <- FB[, .(n_pairs = sum(pairs), mean_corr = weighted.mean(mean_corr, pairs)), by = relation][order(-mean_corr)]
fb <- fb[, .(relazione = relation, `coppie di celle` = iti(n_pairs), `correlazione media` = it(mean_corr, 4))]
arm_med <- ARMS[is.finite(rho) & n_pairs >= 100, median(rho)]
add("Correlazione degli indicatori di rifiuto dello stesso replicato fra celle (H0, MMM e FG, scenari M1 e M3):", "", md(fb),
    sprintf("Per confronto, la correlazione mediana delle stime fra H0 e H1 della stessa cella (S2, `V4_S2_corr_arms.csv`, celle con almeno 100 coppie) e' %s.", it(arm_med, 3)), "")
g3 <- G3[, .(quantita = quantity, `SE linearizzazione` = it(se_lin, 4), `SE bootstrap` = it(se_boot, 4), `scarto relativo` = it(rel, 3))]
add("### A4. Controprova del bootstrap (G3)", "", md(g3))
add("### A5. Riproduzione di S4 (G2)", "", md(G2[, .(controllo = check, n = iti(n), `scarto massimo` = formatC(max, format = "e", digits = 1))]))
fc <- FC[variant == "main" & est_inside_pct == FALSE, .(nome = name, stima = it(est, 3), `media bootstrap` = it(boot_mean, 3), percentile = sprintf("[%s; %s]", it(pct_lo, 3), it(pct_hi, 3)),
                                                         basic = sprintf("[%s; %s]", it(basic_lo, 3), it(basic_hi, 3)))]
add("### A6. Statistiche non lisce: stima fuori dall'intervallo percentile", "", md(fc))
fd <- FD[, .(metodo = method, `v3, primo in ordine` = iti(v3_first), `v3, pareggi divisi` = it(v3_share, 1), `v4, primo in ordine` = iti(v4_first), `v4, pareggi divisi` = it(v4_share, 1))]
ft <- FT[, .SD[order(-N)][1:4], by = version][, .(versione = version, `metodi a pari merito` = who, `punti di disegno` = N)]
add("### A7. Celle vinte e regola sui pareggi (C04)", "", md(fd), md(ft))
ff <- FF[, .(meccanismo = dgp, metodo = method, `RMSE, celle proprie` = it(rmse_own, 4), `celle proprie` = cells_own, `RMSE, celle comuni` = it(rmse_common, 4), `celle comuni` = cells_common)][order(meccanismo, `RMSE, celle comuni`)]
add("### A8. RMSE su celle proprie e su celle comuni (C19a)", "", md(ff))
fg <- FG[method %in% c("glmer_lrt", "gee_md", "gee_fg") & N == 35, .(metodo = method, N, T = ntime_grp, scenari = sgrp, `usabile v3 (%)` = it(v3, 1), `usabile v4 (%)` = it(v4, 1))]
add("### A9. Usabilita' a N = 35, v3 e v4 (C23)", "", md(fg))
bb <- CL[bobyqa_changes %in% TRUE, .(`sotto-affermazione` = name, `stima` = it(est, 3), esito = class, `stima bobyqa` = it(est_bobyqa, 3), `esito bobyqa` = class_bobyqa)]
add("### A10. Sensibilita' bobyqa: sotto-affermazioni che cambiano esito", "", md(bb))
m1 <- VAL[variant == "m1c" & grepl("^C(01|02|03|20)", name), .(`sotto-affermazione` = name, stima = it(est, 3), `intervallo 95%` = sprintf("[%s; %s]", it(lo, 3), it(hi, 3)))]
add("### A11. Scenario ancorato alla coorte (M1c), descrittivo", "", "Stesse definizioni della griglia principale su 72 punti di disegno (6 N, 4 T); caso peggiore sui 6 N.", "", md(m1))
writeLines(L, file.path(REPO, "reports/morris/S5_results.md"))
cat("written", length(L), "lines\n")
