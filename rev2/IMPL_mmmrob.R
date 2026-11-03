.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
# =============================================================================
# IMPL_mmmrob.R -- ricostruzione dello script perduto che ha prodotto
# rev2/MMMROB_raw.csv (D2-37, M-30).
#
# DOMANDA. summary.MMLongit costruisce mean.table da sqrt(diag(mod.cov)) e passa
# a rob.cov solo se control[7] (analisi a verosimiglianza pesata, samp.probi)
# e' vero. La pipeline legge summary(mod)$mean.table (analysis_model1.R riga 15
# e omologhe), quindi l'MMM del manoscritto usa SEMPRE la varianza model-based;
# rob.cov e' presente nell'oggetto MMLongit e non viene mai letta. Qui si
# rifitta la stessa griglia estraendo ENTRAMBE le varianze DALLO STESSO fit, e
# si misura il livello dell'MMM sotto le due.
#
# DISEGNO. Braccio H0 (label="manual"), scenari M1 e M3, N in {10,25,45,145},
# T in {4,10}, tre meccanismi (GenBinaryY, GLMM-native, GLMM-AR1) = 48 celle,
# 250 repliche per cella (12.000 tentativi).
#
# SEME. RNGkind("L'Ecuyer-CMRG") e set.seed(20260805L + 1000L*model + i) prima
# della generazione, i indice di replica.
#
# IL GENERATORE PSEUDOCASUALE VA DICHIARATO, NON EREDITATO. run_full.R gira
# dentro worker crew/mirai, che impostano RNGkind("L'Ecuyer-CMRG"); set.seed
# NON cambia il tipo di generatore, quindi lo stesso intero produce flussi
# diversi sotto Mersenne-Twister (default di Rscript --vanilla) e sotto
# L'Ecuyer. Verificato contro results_full/GenBinaryY_m1_N10_T4_manual_c01.rds,
# repliche 1-4, coefficiente time:trt dell'MMM:
#   chunk del run principale : 0.824334  0.058120  16.253979  -0.874875
#   L'Ecuyer-CMRG            : 0.824334  0.058120  16.253979  -0.874875  <- coincide
#   Mersenne-Twister        : -19.007885 -0.000000 18.836008   0.023643
# Con L'Ecuyer le repliche sono LE STESSE del run principale, quindi il
# confronto con l'MMM pubblicato (REV_cells.csv) e' appaiato replica per
# replica. MMMROB_raw.csv, il file non riproducibile, non coincide con i chunk:
# era stato prodotto sotto Mersenne-Twister.
#
# CONVENZIONI. Fit utilizzabile: est ed SE finiti, |est| < 10, 0 < SE < 10
# (la stessa regola di REV_agg.R). Il livello e' 100*k_rej/n_usable, il rifiuto
# e' |est/SE| > qnorm(0.975).
# Due convenzioni, entrambe riportate, perche' i due estimatori di varianza
# FALLISCONO SU REPLICHE DIVERSE (a N=10 la model-based e' utilizzabile in ~66%
# delle repliche e la robusta in ~76%):
#   (A) FILTRI SEPARATI -- ciascun livello sul proprio denominatore. E' la
#       convenzione del confronto pubblicato: e' quella che riproduce i valori
#       di MMMROB_raw.csv, ed e' anche la lettura corretta della domanda, perche'
#       un analista che scegliesse la varianza robusta avrebbe a disposizione
#       quel campione di repliche e non un altro.
#   (B) FILTRO CONGIUNTO -- solo le repliche in cui ENTRAMBE le SE sono
#       utilizzabili, cioe' un confronto appaiato. Sposta i livelli di alcuni
#       decimi di punto; riportato accanto ad (A).
# Il worst case e' il massimo su N di |tasso - 5|, con il tasso aggregato per
# repliche entro N (pooling su scenario e T), come in Sezione 3.1.
# =============================================================================
suppressMessages({library(data.table); library(MMLB); library(MASS); library(parallel)
                  library(dplyr); library(fastGHQuad)})
SRC <- path.expand("~/MMM/R"); RD <- path.expand("~/MMM/.rundir")
OUT <- path.expand("~/MMM/rev2")

## --- generatori: uno per meccanismo, in ambienti separati (come run_full.R) ---
eG <- new.env(parent = globalenv())
for (f in c("generation_model1.R","generation_model3.R","generate_data.R","make_all_data.R"))
  sys.source(file.path(SRC, f), envir = eG)
mk_gby <- get("make_all_data", envir = eG)

eN <- new.env(parent = globalenv())
sys.source(file.path(RD, "gen_glmm_data.R"), envir = eN)
for (f in c("manual/generation_model1.R","manual/generation_model3.R",
            "manual/generate_data.R","manual/make_all_data.R"))
  sys.source(file.path(RD, f), envir = eN)
mk_nat <- get("make_all_data", envir = eN)

eA <- new.env(parent = globalenv())
sys.source(file.path(RD, "gen_ar1_data.R"), envir = eA)
for (f in c("ar1/generation_models_ar1.R","manual/generate_data.R","ar1/make_all_data.R"))
  sys.source(file.path(RD, f), envir = eA)
mk_ar1 <- get("make_all_data", envir = eA)

gen <- function(dgp, m, N, T, i, label) {
  suppressWarnings(RNGkind("L'Ecuyer-CMRG"))   # come i worker crew di run_full.R
  set.seed(20260805L + 1000L*m + i)
  f <- switch(dgp, "GenBinaryY" = mk_gby, "GLMM-native" = mk_nat, "GLMM-AR1" = mk_ar1)
  as.data.frame(f(N, T, i, as.character(m), label))
}

## --- fit MMM con la specificazione della pipeline per lo scenario ---
fit_mmm <- function(d, m) {
  lvf <- if (m == 1) ~1 else ~time
  suppressWarnings(mm(binY ~ time*trt, lv.formula = lvf, t.formula = ~1,
                      data = d, id = id, q = 20))
}

one <- function(dgp, m, N, T, i, label) {
  d <- tryCatch(gen(dgp, m, N, T, i, label), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  mo <- tryCatch(fit_mmm(d, m), error = function(e) NULL)
  if (is.null(mo)) return(NULL)
  tb <- tryCatch(summary(mo)$mean.table, error = function(e) NULL)
  if (is.null(tb) || !"time:trt" %in% rownames(tb)) return(NULL)
  k <- which(rownames(tb) == "time:trt")
  data.table(dgp = dgp, model = m, N = N, ntime = T, label = label, rep = i,
             est = tb[k,1], se_tab = tb[k,2],
             se_mod = tryCatch(sqrt(diag(mo$mod.cov))[k], error = function(e) NA_real_),
             se_rob = tryCatch(sqrt(diag(mo$rob.cov))[k], error = function(e) NA_real_))
}

CORES <- 25L
NREP  <- 250L
grid  <- CJ(dgp = c("GenBinaryY","GLMM-native","GLMM-AR1"), model = c(1L,3L),
            N = c(10L,25L,45L,145L), ntime = c(4L,10L), sorted = FALSE)
cat("celle:", nrow(grid), " tentativi previsti:", nrow(grid)*NREP, "\n")

RAW <- file.path(OUT, "IMPL_MMMROB_raw_full.csv")
t0 <- Sys.time()
if (file.exists(RAW) && !nzchar(Sys.getenv("MMMROB_REFIT"))) {
  # Riesecuzione della sola aggregazione: i fit sono deterministici dato il seme,
  # quindi rileggerli e' equivalente a rifarli. MMMROB_REFIT=1 forza il rifit.
  D <- fread(RAW); cat("rilettura di", RAW, ":", nrow(D), "fit\n")
} else {
acc <- vector("list", nrow(grid))
for (k in seq_len(nrow(grid))) {
  g <- grid[k]
  rr <- mclapply(seq_len(NREP),
                 function(i) tryCatch(one(g$dgp, g$model, g$N, g$ntime, i, "manual"),
                                      error = function(e) NULL),
                 mc.cores = CORES)
  acc[[k]] <- rbindlist(Filter(is.data.table, rr))
  cat(sprintf("[%s] %-11s m%d N=%3d T=%2d : %3d fit (%.1f min)\n", format(Sys.time(), "%H:%M:%S"),
              g$dgp, g$model, g$N, g$ntime, nrow(acc[[k]]),
              as.numeric(difftime(Sys.time(), t0, units = "mins")))); flush.console()
}
D <- rbindlist(acc)
D[, `:=`(se_mod = as.numeric(se_mod), se_rob = as.numeric(se_rob))]
fwrite(D[, .(dgp, model, N, ntime, rep, est, se_mod, se_rob)],
       file.path(OUT, "IMPL_MMMROB_raw.csv"))
fwrite(D, RAW)
}

REP <- file.path(OUT, "IMPL_MMMROB_report.txt")
sink(REP, split = TRUE)
zq <- qnorm(0.975)
ok_est <- function(x) is.finite(x$est) & abs(x$est) < 10
ok_se  <- function(s) is.finite(s) & s > 0 & s < 10
cat("\n===== 0. controllo del seme contro i chunk del run principale =====\n")
# I chunk di results_full sono la fonte primaria: se il seme e il generatore
# pseudocasuale sono quelli giusti, le stime dell'MMM coincidono replica per
# replica. Nessun valore di questo script si consegna se questo controllo fallisce.
ck <- path.expand("~/MMM/results_full/GenBinaryY_m1_N10_T4_manual_c01.rds")
if (file.exists(ck)) {
  cq <- as.data.table(readRDS(ck))[term=="time:trt" & method=="mmm", .(rep=sim_id, est_chunk=estimate)]
  cn <- D[dgp=="GenBinaryY" & model==1 & N==10 & ntime==4, .(rep, est_new=est)]
  MC <- merge(cq, cn, by="rep")
  cat(sprintf("repliche appaiate col chunk c01: %d | identiche a 1e-8: %d | max |diff| %.3e\n",
              nrow(MC), sum(abs(MC$est_chunk-MC$est_new) < 1e-8, na.rm=TRUE),
              max(abs(MC$est_chunk-MC$est_new), na.rm=TRUE)))
  print(head(MC, 4))
} else cat("chunk di riferimento assente: controllo non eseguito\n")

cat("\n===== 1. mean.table = model-based, non robusta =====\n")
cat(sprintf("fit con mean.table: %d ; max |se_tab - sqrt(diag(mod.cov))| = %.3e\n",
            nrow(D), max(abs(D$se_tab - D$se_mod), na.rm = TRUE)))
cat(sprintf("max |se_tab - sqrt(diag(rob.cov))| = %.3e ; mediana = %.3e\n",
            max(abs(D$se_tab - D$se_rob), na.rm = TRUE),
            median(abs(D$se_tab - D$se_rob), na.rm = TRUE)))
cat(sprintf("rob.cov disponibile in %d fit su %d (%.2f%%)\n",
            sum(is.finite(D$se_rob)), nrow(D), 100*mean(is.finite(D$se_rob))))

## ---- convenzione (A): FILTRI SEPARATI -------------------------------------
E  <- D[ok_est(D)]
Um <- E[ok_se(se_mod)]; Ur <- E[ok_se(se_rob)]
cat(sprintf("\ntentativi: %d ; fit restituiti: %d ; est utilizzabile: %d\n",
            nrow(grid)*NREP, nrow(D), nrow(E)))
cat(sprintf("utilizzabili con SE model-based: %d (%.2f%% dei tentativi)\n",
            nrow(Um), 100*nrow(Um)/(nrow(grid)*NREP)))
cat(sprintf("utilizzabili con SE robusta    : %d (%.2f%% dei tentativi)\n",
            nrow(Ur), 100*nrow(Ur)/(nrow(grid)*NREP)))

cat("\n===== 2A. LIVELLO (%) per N, FILTRI SEPARATI (convenzione pubblicata) =====\n")
am <- Um[, .(n_mod = .N, lev_mod = 100*mean(abs(est/se_mod) > zq)), by = .(dgp, N)]
ar <- Ur[, .(n_rob = .N, lev_rob = 100*mean(abs(est/se_rob) > zq)), by = .(dgp, N)]
byN <- merge(am, ar, by = c("dgp","N"))
byN[, `:=`(shift = lev_rob - lev_mod,
           mcse_mod = 100*sqrt(0.05*0.95/n_mod), mcse_rob = 100*sqrt(0.05*0.95/n_rob))]
print(as.data.frame(byN[order(dgp, N), .(dgp, N, n_mod, lev_mod = round(lev_mod,2),
      mcse_mod = round(mcse_mod,2), n_rob, lev_rob = round(lev_rob,2),
      mcse_rob = round(mcse_rob,2), shift = round(shift,2))]), row.names = FALSE)
cat("\n-- pooled per meccanismo, filtri separati --\n")
PA <- merge(Um[, .(n_mod = .N, lev_mod = round(100*mean(abs(est/se_mod) > zq),2)), by = dgp],
            Ur[, .(n_rob = .N, lev_rob = round(100*mean(abs(est/se_rob) > zq),2)), by = dgp], by = "dgp")
PA[, shift := round(lev_rob - lev_mod, 2)]; print(as.data.frame(PA), row.names = FALSE)

cat("\n===== 2B. LIVELLO (%) per N, FILTRO CONGIUNTO (confronto appaiato) =====\n")
U <- E[ok_se(se_mod) & ok_se(se_rob)]
cat(sprintf("fit con ENTRAMBE le SE utilizzabili: %d su %d tentativi\n", nrow(U), nrow(grid)*NREP))
U[, `:=`(rej_mod = abs(est/se_mod) > zq, rej_rob = abs(est/se_rob) > zq)]
byNj <- U[, .(n = .N, lev_mod = 100*mean(rej_mod), lev_rob = 100*mean(rej_rob)), by = .(dgp, N)]
byNj[, `:=`(shift = lev_rob - lev_mod, mcse = 100*sqrt(0.05*0.95/n))]
print(as.data.frame(byNj[order(dgp, N), .(dgp, N, n, lev_mod = round(lev_mod,2),
      lev_rob = round(lev_rob,2), shift = round(shift,2), mcse = round(mcse,2))]), row.names = FALSE)
cat("\n-- pooled per meccanismo, filtro congiunto --\n")
P <- U[, .(n = .N, lev_mod = round(100*mean(rej_mod),2), lev_rob = round(100*mean(rej_rob),2)), by = dgp]
P[, shift := round(lev_rob - lev_mod, 2)]; print(as.data.frame(P), row.names = FALSE)

cat("\n-- rapporto SE_rob/SE_mod (sui fit con entrambe le SE) --\n")
RR <- U[, .(n = .N, ratio_mean = mean(se_rob/se_mod), ratio_med = median(se_rob/se_mod)), by = .(dgp, N)]
print(as.data.frame(RR[order(dgp,N), .(dgp, N, n, ratio_mean = round(ratio_mean,4),
      ratio_med = round(ratio_med,4))]), row.names = FALSE)
cat(sprintf("a N=10, tutti i meccanismi: rapporto medio %.4f (robusta %+.1f%% rispetto alla model-based, n=%d)\n",
            mean(U[N==10]$se_rob/U[N==10]$se_mod),
            100*(mean(U[N==10]$se_rob/U[N==10]$se_mod)-1), nrow(U[N==10])))

cat("\n===== 3. WORST CASE (max su N di |tasso - 5|), filtri separati =====\n")
W <- byN[, .(worst_mod = max(abs(lev_mod-5)), N_mod = N[which.max(abs(lev_mod-5))],
             worst_rob = max(abs(lev_rob-5)), N_rob = N[which.max(abs(lev_rob-5))]), by = dgp]
print(as.data.frame(W[, .(dgp, worst_mod = round(worst_mod,2), N_mod,
                          worst_rob = round(worst_rob,2), N_rob)]), row.names = FALSE)
cat(sprintf("intervallo fra meccanismi: model-based %.2f-%.2f pp ; robusta %.2f-%.2f pp\n",
            min(W$worst_mod), max(W$worst_mod), min(W$worst_rob), max(W$worst_rob)))
cat("-- lo stesso, filtro congiunto --\n")
Wj <- byNj[, .(worst_mod = round(max(abs(lev_mod-5)),2), worst_rob = round(max(abs(lev_rob-5)),2)), by = dgp]
print(as.data.frame(Wj), row.names = FALSE)
cat("-- variante per cella (dgp,model,N,T), non aggregata, filtri separati --\n")
bcm <- Um[, .(n = .N, lev = 100*mean(abs(est/se_mod) > zq)), by = .(dgp, model, N, ntime)]
bcr <- Ur[, .(n = .N, lev = 100*mean(abs(est/se_rob) > zq)), by = .(dgp, model, N, ntime)]
print(as.data.frame(merge(bcm[, .(worst_mod = round(max(abs(lev-5)),2)), by = dgp],
                          bcr[, .(worst_rob = round(max(abs(lev-5)),2)), by = dgp], by = "dgp")),
      row.names = FALSE)
fwrite(byN,  file.path(OUT, "IMPL_MMMROB_byN.csv"))
fwrite(byNj, file.path(OUT, "IMPL_MMMROB_byN_joint.csv"))
fwrite(merge(bcm, bcr, by = c("dgp","model","N","ntime"), suffixes = c("_mod","_rob")),
       file.path(OUT, "IMPL_MMMROB_bycell.csv"))

## --- M-30: calibrazione della SE model-based a N=60, T=6 -------------------
cat("\n===== 4. M-30: calibrazione a GenBinaryY M3, N=60, T=6, 400 repliche =====\n")
for (lb in c("manual","mmlb")) {
  rr <- rbindlist(Filter(is.data.table, mclapply(seq_len(400L),
        function(i) tryCatch(one("GenBinaryY", 3L, 60L, 6L, i, lb), error = function(e) NULL),
        mc.cores = CORES)))
  V <- rr[is.finite(est) & abs(est) < 10 & is.finite(se_mod) & se_mod > 0 & se_mod < 10]
  cat(sprintf("label=%-6s : fit %3d/400, utilizzabili %3d | SE media (model-based) %.5f | SD empirica %.5f | rapporto %.4f\n",
              lb, nrow(rr), nrow(V), mean(V$se_mod), sd(V$est), mean(V$se_mod)/sd(V$est)))
  Vr <- V[is.finite(se_rob) & se_rob > 0 & se_rob < 10]
  cat(sprintf("               robusta: n %3d | SE media %.5f | SD empirica %.5f | rapporto %.4f | SE_rob/SE_mod medio %.4f\n",
              nrow(Vr), mean(Vr$se_rob), sd(Vr$est), mean(Vr$se_rob)/sd(Vr$est), mean(Vr$se_rob/Vr$se_mod)))
  fwrite(rr, file.path(OUT, sprintf("IMPL_MMMROB_n60_%s.csv", lb)))
}
## --- 5. confronto col worst case pubblicato in Sezione 3.1 ------------------
# §3.1 calcola il worst case sui 16 livelli di N dell'intera griglia (5 scenari,
# 6 valori di T). Questa griglia ne copre un sottoinsieme (M1, M3; N in 4 valori;
# T in {4,10}), quindi il confronto onesto e' col worst case di REV_cells.csv
# RISTRETTO allo stesso sottoinsieme, riportato accanto al valore pubblicato.
cat("\n===== 5. worst case dell'MMM: pubblicato, sottogriglia, robusta =====\n")
RC <- fread(file.path(OUT, "REV_cells.csv"))
H  <- RC[label=="manual" & set=="own" & method=="mmm"]
pub <- H[, .(n = sum(n), k = sum(k_rej)), by = .(dgp, N)][, rate := 100*k/n][
        , .(worst_pubblicato = round(max(abs(rate-5)),2)), by = dgp]
sub <- H[model %in% c(1,3) & N %in% c(10,25,45,145) & ntime %in% c(4,10)][
        , .(n = sum(n), k = sum(k_rej)), by = .(dgp, N)][, rate := 100*k/n][
        , .(worst_sottogriglia = round(max(abs(rate-5)),2)), by = dgp]
cmp <- merge(merge(pub, sub, by = "dgp"),
             byN[, .(worst_ricalcolato_mod = round(max(abs(lev_mod-5)),2),
                     worst_ricalcolato_rob = round(max(abs(lev_rob-5)),2)), by = dgp], by = "dgp")
print(as.data.frame(cmp), row.names = FALSE)
cat(sprintf("intervallo fra meccanismi -- pubblicato %.2f-%.2f | sottogriglia %.2f-%.2f | mod ricalcolato %.2f-%.2f | ROBUSTA %.2f-%.2f pp\n",
            min(cmp$worst_pubblicato), max(cmp$worst_pubblicato),
            min(cmp$worst_sottogriglia), max(cmp$worst_sottogriglia),
            min(cmp$worst_ricalcolato_mod), max(cmp$worst_ricalcolato_mod),
            min(cmp$worst_ricalcolato_rob), max(cmp$worst_ricalcolato_rob)))
fwrite(cmp, file.path(OUT, "IMPL_MMMROB_worstcase.csv"))
cat(sprintf("\ntempo totale %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("DONE_IMPL_MMMROB\n")
sink()
