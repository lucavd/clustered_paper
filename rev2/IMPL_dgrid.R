.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
# =============================================================================
# IMPL_dgrid.R -- ricostruzione dello script perduto che ha prodotto
# rev2/DGRID_raw.csv e rev2/DGRID_publishable.csv (D2-16r, D2-36).
#
# DOMANDA. Il disegno principale (run_full.R) non ha alcun fattore che misuri la
# forza della dipendenza: i suoi fattori sono meccanismo, scenario, N, T, braccio.
# La varianza degli effetti casuali e' fissata dentro ciascun generatore. Qui la
# si fa variare come fattore esplicito, tenendo fisso tutto il resto, per vedere
# se l'ordinamento dei metodi dipende dalla forza della dipendenza piu' che dal
# meccanismo.
#
# PROTOCOLLO. Meccanismo GLMM-native, scenario M1, braccio H0 (label "manual",
# time:trt marginale = 0); N in {10,25,45}; T in {4,10}; D in {0.25,1,4},
# corrispondenti a ICC latente D/(D+pi^2/3) = 0.0705 / 0.2330 / 0.5485.
# 18 celle x 300 repliche tentate. Tutti e otto i metodi con la pipeline
# invariata (R/analysis_model1.R, R/analysis_smallsample.R, R/do_analysis.R).
#
# beta_c RISCALATO. La pipeline a D=1 usa beta_c = (-2.32, 0.116, 0.232, 0),
# cioe' il bersaglio marginale (-2, 0.1, 0.2, 0) moltiplicato per il fattore di
# Zeger sqrt(1 + c^2 D) con c^2 = 0.346. Qui lo stesso riscalamento e' applicato
# a ogni D, cosi' il modello medio marginale resta (-2, 0.1, 0.2, 0) e l'unica
# cosa che cambia fra le tre condizioni e' la forza della dipendenza. (A D=1 il
# valore esatto -2.32047 differisce dalla costante arrotondata -2.32 del
# generatore di 5e-4: differenza dichiarata, non corretta.)
#
# SEME. RNGkind("L'Ecuyer-CMRG") e set.seed(500000 + 7919*i + j), i indice di
# cella (1..18, ordine della griglia stampata sotto), j indice di replica.
# Il tipo di generatore va DICHIARATO: set.seed non lo cambia, e run_full.R gira
# in worker crew/mirai che usano L'Ecuyer-CMRG mentre Rscript --vanilla parte in
# Mersenne-Twister. Verificato sui chunk di results_full (vedi IMPL_mmmrob.R,
# sezione 0): con L'Ecuyer le stime dell'MMM coincidono replica per replica col
# run principale, con Mersenne no. Cosi' questa griglia condivide con il disegno
# principale il generatore pseudocasuale, non soltanto il codice.
#
# CONVENZIONI. Utilizzabilita' e livello esattamente come REV_agg.R: per i sette
# metodi di Wald converged & est,se,p finiti & |est| < 10 & 0 < se < 10; per il
# rapporto di verosimiglianza converged & p finito. Livello = 100*k/n_usable con
# k = #(p < 0.05); MCSE = 100*sqrt(p(1-p)/n_usable); denominatore dichiarato in
# ogni tabella. n_att = repliche in cui il metodo ha restituito una riga.
# =============================================================================
suppressMessages({library(data.table); library(MMLB); library(lme4); library(geepack)
                  library(geesmv); library(dplyr); library(tibble); library(MASS)
                  library(fastGHQuad); library(parallel)})
SRC <- path.expand("~/MMM/R"); RD <- path.expand("~/MMM/.rundir")
OUT <- path.expand("~/MMM/rev2")

## --- pipeline di analisi, invariata ---
for (f in c("analysis_model1.R","analysis_model2.R","analysis_model3.R",
            "analysis_model4.R","analysis_model5.R","analysis_smallsample.R",
            "do_analysis.R")) sys.source(file.path(RD, f), envir = globalenv())
sys.source(file.path(RD, "gen_glmm_data.R"), envir = globalenv())

## --- generatore GLMM-native M1 con D come argomento -------------------------
C2      <- 0.346                      # c^2 di Zeger, come nella pipeline a D=1
TARGET0 <- c(-2, 0.1, 0.2, 0)         # bersaglio MARGINALE, braccio H0
gen_dgrid <- function(N, ntime, Dvar) {
  beta_c <- TARGET0 * sqrt(1 + C2*Dvar)
  # struttura di follow-up identica a manual/generation_model1.R
  nclust <- ifelse(runif(N) < 0.7, ntime, sample(seq_len(ntime - 1), N, replace = TRUE))
  id   <- rep(seq(N), nclust)
  trt  <- rep(rbinom(N, 1, 0.5), nclust)
  time <- unlist(lapply(nclust, function(x) seq(x) - 1))
  data <- data.frame(id, time, trt) |> arrange(id, time)
  gen_glmm_data(mean.formula = ~ time*trt, re.formula = ~1, beta_c = beta_c,
                D = Dvar, id = id, data = data, Yname = "binY") |>
    dplyr::mutate(sim_id = NA_integer_)
}

TECHS <- c("mmm","gee","glmer","gee_ss","glmer_lrt")   # glmer restituisce anche shrink
one_rep <- function(i, j, N, ntime, Dvar) {
  suppressWarnings(RNGkind("L'Ecuyer-CMRG"))   # come i worker crew di run_full.R
  set.seed(500000L + 7919L*i + j)
  d <- tryCatch(gen_dgrid(N, ntime, Dvar), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  d$sim_id <- j
  rows <- lapply(TECHS, function(tc)
    tryCatch(do_analysis(d, N, ntime, "manual", 1, tc, j), error = function(e) NULL))
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  r <- as.data.table(dplyr::bind_rows(rows))
  r <- r[term == "time:trt"]
  if (!nrow(r)) return(NULL)
  r[, .(method, estimate, se, p, converged, sim_id = j, Ncl = N, Tt = ntime, Dvar = Dvar)]
}

grid <- CJ(Dvar = c(0.25,1,4), Ncl = c(10L,25L,45L), Tt = c(4L,10L), sorted = FALSE)
setorder(grid, Dvar, Ncl, Tt)
grid[, i := .I]
cat("=== griglia (i = indice di cella nel seme) ===\n"); print(as.data.frame(grid), row.names = FALSE)
NREP  <- 300L; CORES <- 20L
# modalita' di collaudo: DGRID_SMOKE=n riduce le repliche e tiene solo la prima
# cella di ogni D, per esercitare il percorso completo (otto metodi, forking del
# timeout di geesmv) prima di impegnare l'ora di calcolo.
if (nzchar(Sys.getenv("DGRID_SMOKE"))) {
  NREP <- as.integer(Sys.getenv("DGRID_SMOKE"))
  grid <- grid[Ncl == 10L & Tt == 4L]
  cat("*** SMOKE:", NREP, "repliche su", nrow(grid), "celle ***\n")
}

t0 <- Sys.time(); acc <- vector("list", nrow(grid))
for (k in seq_len(nrow(grid))) {
  g <- grid[k]
  rr <- mclapply(seq_len(NREP),
                 function(j) tryCatch(one_rep(g$i, j, g$Ncl, g$Tt, g$Dvar),
                                     error = function(e) NULL),
                 mc.cores = CORES)
  acc[[k]] <- rbindlist(Filter(is.data.table, rr), fill = TRUE)
  cat(sprintf("[%s] i=%2d D=%.2f N=%2d T=%2d : %5d righe, %3d repliche generate (%.1f min)\n",
              format(Sys.time(), "%H:%M:%S"), g$i, g$Dvar, g$Ncl, g$Tt, nrow(acc[[k]]),
              uniqueN(acc[[k]]$sim_id), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  flush.console()
}
D <- rbindlist(acc, fill = TRUE)
D[, converged := as.logical(converged)][is.na(converged), converged := FALSE]
fwrite(D, file.path(OUT, "IMPL_DGRID_raw.csv"))
att <- D[, .(n_att = uniqueN(sim_id)), by = .(Dvar,Ncl,Tt,method)]

## --- aggregazione -----------------------------------------------------------
# geesmv e geepack scrivono a livello C sullo stdout del processo: il rapporto
# viene deviato su un file con sink(), che raccoglie solo l'output R-level.
REP <- file.path(OUT, "IMPL_DGRID_report.txt")
sink(REP, split = TRUE)
cat(sprintf("righe grezze: %d ; celle: %d ; metodi: %s\n", nrow(D),
            uniqueN(D[, .(Dvar,Ncl,Tt)]), paste(sort(unique(D$method)), collapse = ",")))
cat("repliche generate per cella (min/mediana/max): ",
    paste(quantile(D[, uniqueN(sim_id), by = .(Dvar,Ncl,Tt)]$V1, c(0,0.5,1)), collapse = " / "), "\n")
cat(sprintf("tentativi per cella-metodo: mediana %.1f, intervallo %d-%d\n",
            median(att$n_att), min(att$n_att), max(att$n_att)))
LRT <- "glmer_lrt"
D[, usable := ifelse(method == LRT, converged & is.finite(p),
                     converged & is.finite(estimate) & is.finite(se) & is.finite(p) &
                       abs(estimate) < 10 & se > 0 & se < 10)]
A <- D[, .(n_att = .N, n_use = sum(usable), k = sum(usable & p < 0.05)), by = .(method, Dvar)]
A[, `:=`(usable = 100*n_use/n_att, level = 100*k/n_use)]
A[, `:=`(mcse = 100*sqrt((k/n_use)*(1-k/n_use)/n_use), dev = abs(level - 5))]
A[, icc := Dvar/(Dvar + pi^2/3)]
setorder(A, Dvar, dev)
A[, rk := seq_len(.N), by = Dvar]
fwrite(A[, .(method, Dvar, n_att, n_use, k, usable, level, mcse, dev, icc = round(icc,3), rk)],
       file.path(OUT, "IMPL_DGRID_publishable.csv"))
cat("\n===== LIVELLO (%) e MCSE (pp) per metodo e forza di dipendenza =====\n")
print(as.data.frame(A[, .(Dvar, icc = round(icc,4), method, n_att, n_use,
      usable = round(usable,2), level = round(level,2), mcse = round(mcse,3),
      dev = round(dev,2), rk)]), row.names = FALSE)

cat("\n===== ORDINAMENTI e correlazione di rango fra livelli di dipendenza =====\n")
W <- dcast(A, method ~ Dvar, value.var = "rk")
setnames(W, c("method","rk_D025","rk_D1","rk_D4")); print(as.data.frame(W), row.names = FALSE)
sp <- function(a,b) cor(a, b, method = "spearman")
cat(sprintf("Spearman fra ordinamenti: ICC 0.071 vs 0.233 = %.4f ; 0.071 vs 0.549 = %.4f ; 0.233 vs 0.549 = %.4f\n",
            sp(W$rk_D025,W$rk_D1), sp(W$rk_D025,W$rk_D4), sp(W$rk_D1,W$rk_D4)))
Wd <- dcast(A, method ~ Dvar, value.var = "dev")
setnames(Wd, c("method","d025","d1","d4"))
cat(sprintf("Spearman sugli scostamenti |livello-5|: 0.071 vs 0.233 = %.4f ; 0.071 vs 0.549 = %.4f ; 0.233 vs 0.549 = %.4f\n",
            sp(Wd$d025,Wd$d1), sp(Wd$d025,Wd$d4), sp(Wd$d1,Wd$d4)))

cat("\n===== TEST: il livello cambia con la forza della dipendenza? =====\n")
cat("prop.test a 2 campioni su k/n_use, D=0.25 contro D=4, per metodo\n")
TT <- rbindlist(lapply(sort(unique(A$method)), function(mm) {
  a <- A[method==mm & Dvar==0.25]; b <- A[method==mm & Dvar==4]
  ht <- prop.test(c(a$k,b$k), c(a$n_use,b$n_use), correct = FALSE)
  z  <- A[method==mm][order(Dvar)]
  tr <- tryCatch(prop.trend.test(z$k, z$n_use, score = z$icc), error = function(e) NULL)
  data.table(method = mm, lev_low = a$level, lev_high = b$level,
             delta = b$level - a$level, p_2camp = ht$p.value,
             p_trend = if (is.null(tr)) NA_real_ else tr$p.value)
}))
setorder(TT, p_2camp)
print(as.data.frame(TT[, .(method, lev_low = round(lev_low,2), lev_high = round(lev_high,2),
      delta = round(delta,2), p_2camp = signif(p_2camp,3), p_trend = signif(p_trend,3))]),
      row.names = FALSE)
fwrite(TT, file.path(OUT, "IMPL_DGRID_tests.csv"))

## --- confronto con l'ordinamento fra meccanismi del disegno principale ------
cat("\n===== confronto: correlazione di rango FRA MECCANISMI (REV_cells.csv) =====\n")
R <- fread(file.path(OUT, "REV_cells.csv"))
H <- R[label=="manual" & set=="own"]
sub <- H[model==1 & N %in% c(10,25,45) & ntime %in% c(4,10)]
mk_rank <- function(x) {
  z <- x[, .(n = sum(n), k = sum(k_rej)), by = .(dgp, method)][, dev := abs(100*k/n - 5)]
  z[, rk := frank(dev, ties.method = "first"), by = dgp]
  dcast(z, method ~ dgp, value.var = "rk")
}
for (tag in c("sottoinsieme M1,N<=45,T in {4,10}", "griglia intera")) {
  Q <- mk_rank(if (tag == "griglia intera") H else sub)
  cat("--", tag, "\n"); print(as.data.frame(Q), row.names = FALSE)
  nmv <- setdiff(names(Q), "method")
  for (a in seq_along(nmv)) for (b in seq_along(nmv)) if (a < b)
    cat(sprintf("   Spearman %s vs %s = %.4f\n", nmv[a], nmv[b],
                cor(Q[[nmv[a]]], Q[[nmv[b]]], method = "spearman")))
}

## --- il disegno principale non ha un fattore di dipendenza -----------------
cat("\n===== fattori del disegno in run_full.R =====\n")
rf <- readLines(file.path(SRC, "run_full.R"))
cat(grep("^(N_LEVELS|T_LEVELS|MODELS|LABELS|DGPS|N_SIM|CHUNK|TECHS)", rf, value = TRUE), sep = "\n")
cat("\nrighe di expand.grid che definiscono le celle:\n")
cat(grep("expand.grid", rf, value = TRUE), sep = "\n")
cat("\nmenzioni di rho / ICC / varianza fra i fattori:",
    length(grep("rho|ICC|icc|Dvar|variance component", rf, value = TRUE)), "\n")
## --- parametri fissi dei tre generatori ------------------------------------
cat("\n===== parametri fissi dei generatori (grep sul sorgente) =====\n")
for (m in 1:5) {
  g <- readLines(file.path(SRC, sprintf("generation_model%d.R", m)))
  cat(sprintf("-- GenBinaryY M%d: %s | %s\n", m,
      paste(trimws(grep("mmlb *= *(c\\()?[0-9]", g, value = TRUE)), collapse = " ;; "),
      paste(trimws(grep("^ *(lv|t)\\.formula", g, value = TRUE)), collapse = " ; ")))
}
for (m in 1:5) {
  g <- readLines(file.path(RD, sprintf("manual/generation_model%d.R", m)))
  cat(sprintf("-- GLMM-native M%d: D = %s\n", m,
      paste(trimws(grep("^ *(D <-|1\\.0|0\\.01|0\\.001)", g, value = TRUE)), collapse = " ")))
}
ga <- readLines(file.path(RD, "ar1/generation_models_ar1.R"))
cat("-- GLMM-AR1: ", paste(trimws(grep("^AR1_(RHO|SIGMA2)", ga, value = TRUE)), collapse = " ; "), "\n")
cat("-- GLMM-AR1 D per scenario:\n")
cat(paste(trimws(grep("D = 1.0|matrix\\(c\\(|^ *0\\.01|^ *1\\.0", ga, value = TRUE)), collapse = "\n"), "\n")
cat(sprintf("\ntempo totale %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("DONE_IMPL_DGRID\n")
sink()
