.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(dplyr); library(lme4); library(geepack); library(MMLB); library(data.table)})
setwd(path.expand("~/MMM"))
d0 <- readRDS("ARTICOLO 2026/LURID_bis/lurid_long.rds")

## ---------- 1. ricostruzione timepoint a blocchi di dieci zone ----------
d0$row_in_id <- ave(seq_len(nrow(d0)), d0$id, FUN = seq_along)
d0$tp  <- (d0$row_in_id - 1L) %/% 10L + 1L
d0$pt  <- factor(paste(d0$id, d0$tp, sep = "_"))
chk <- tapply(as.character(d0$section), d0$pt, function(s) length(unique(s)))
cat("=== T1: blocchi di dieci zone ===\n")
cat("n patient-timepoint:", nlevels(d0$pt), " blocchi con 10 sezioni distinte:", sum(chk==10), "/", length(chk), "\n")
cat("colonna timepoint esplicita nei dati:", any(grepl("^(tp|timepoint|visit|time)$", names(readRDS('ARTICOLO 2026/LURID_bis/lurid_long.rds')))), "\n")
cat("righe per pt:", paste(names(table(table(d0$pt))), table(table(d0$pt)), sep="x", collapse=" "), "\n")

## ---------- 2. region / lus_severe (R3-V11) ----------
d <- d0 %>% filter(!is.na(lus_severe)) %>%
  mutate(region = factor(case_when(grepl("A$", as.character(section)) ~ "anterior",
                                   grepl("P$", as.character(section)) ~ "posterior",
                                   TRUE ~ "lateral"),
                         levels = c("anterior","lateral","posterior")),
         id = droplevels(factor(id)), pt = droplevels(factor(pt)),
         hemi = factor(paste(pt, side, sep="_")),
         post_bin = as.integer(region == "posterior"),
         region_lin = as.numeric(region))
cat("\n=== R3-V11 ===\n")
cat("righe con outcome:", nrow(d), " pazienti:", n_distinct(d$id), "\n")
print(table(d$region))
cat("prevalenza lus_severe per regione:\n"); print(round(tapply(d$lus_severe, d$region, mean), 4))
cat("prevalenza complessiva:", round(mean(d$lus_severe),4), "\n")
cat("mappa lus_score->lus_severe:\n"); print(table(d$lus_score, d$lus_severe))

## ---------- 3. invasive: quanti pazienti (R3-S08 / R3-S07 / R3-E03) ----------
cat("\n=== invasive per paziente ===\n")
inv <- d %>% group_by(id) %>% summarise(any_inv = any(invasive == 1), all_inv = all(invasive == 1),
                                        nlev = n_distinct(invasive), .groups="drop")
cat("pazienti totali:", nrow(inv), " con almeno una riga invasive=1:", sum(inv$any_inv),
    " con tutte invasive=1:", sum(inv$all_inv), " con entrambi i livelli entro paziente:", sum(inv$nlev==2), "\n")
cat("righe invasive=1:", sum(d$invasive==1), " patient-timepoint invasive=1:", n_distinct(d$pt[d$invasive==1]), "\n")
K <- sum(inv$any_inv); P <- nrow(inv)
cat("\n=== ipergeometrica: P(nessun paziente invasivo) con", K, "invasivi su", P, "===\n")
for (n in c(10,15,20,25,50)) {
  p0 <- exp(lchoose(P-K, n) - lchoose(P, n))
  cat(sprintf("N=%2d  P0=%.6f  tetto fit stimabili=%.4f%%  (dhyper=%.6f)\n", n, p0, 100*(1-p0), dhyper(0, K, P-K, n)))
}
cat("\n(controllo con K=25 come dichiara il piano)\n")
for (n in c(10,15,20,25,50)) cat(sprintf("N=%2d  P0(K=25)=%.6f\n", n, dhyper(0, 25, 207-25, n)))

## ---------- 4. estrazioni con il seme dichiarato: degenerazione ----------
cat("\n=== estrazioni seed 20260807, ordine dell'autore ===\n")
set.seed(20260807)
ids <- levels(d$id); NLEV <- c(10,15,20,25,50); REPS <- 500
draws <- list(); k <- 0L
for (n in NLEV) for (r in seq_len(REPS)) { k <- k+1L; draws[[k]] <- list(N=n, rep=r, sel=sample(ids, n)) }
saveRDS(draws, "rev2/T3_draws.rds")
dg <- rbindlist(lapply(draws, function(z) {
  s <- d[d$id %in% z$sel, ]
  data.table(N=z$N, rep=z$rep,
             n_inv_pat = length(unique(s$id[s$invasive==1])),
             lev_inv = length(unique(as.character(s$invasive))),
             lev_dis = length(unique(as.character(s$disease))),
             lev_age = length(unique(s$age_weeks)),
             nrows = nrow(s), nev = sum(s$lus_severe))
}))
fwrite(dg, "rev2/T3_degeneracy.csv")
cat("percentuale di estrazioni con invasive a un solo livello:\n")
print(dg[, .(n=.N, degen=sum(lev_inv==1), pct_degen=round(100*mean(lev_inv==1),2),
             pct_dis_lt4=round(100*mean(lev_dis<4),2), pct_age_lt2=round(100*mean(lev_age<2),2),
             pct_ident_full=round(100*mean(lev_inv==2 & lev_dis>=2),2)), by=N][order(N)])

## ---------- 5. R3-V08: chi2 del livello patient-timepoint ----------
cat("\n=== R3-V08 ===\n")
multi <- inv$id[0]; nt <- d %>% group_by(id) %>% summarise(k=n_distinct(pt), .groups="drop")
multi <- as.character(nt$id[nt$k > 1]); single <- as.character(nt$id[nt$k == 1])
cat("pazienti con >1 timepoint:", length(multi), " loro righe:", sum(d$id %in% multi),
    " su", nrow(d), "\n")
cat("distribuzione timepoint per paziente:\n"); print(table(nt$k))
lrt2 <- function(dat) {
  dat <- droplevels(dat)
  m1 <- try(glmer(lus_severe ~ region + (1|id), family=binomial, data=dat,
                  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))), silent=TRUE)
  m2 <- try(glmer(lus_severe ~ region + (1|id) + (1|pt), family=binomial, data=dat,
                  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))), silent=TRUE)
  if (inherits(m1,"try-error")||inherits(m2,"try-error")) return(c(NA,NA,NA))
  ch <- as.numeric(2*(logLik(m2)-logLik(m1)))
  c(ch, unname(unlist(VarCorr(m1))[1]), unname(sum(unlist(VarCorr(m2)))))
}
a_all <- lrt2(d); a_mul <- lrt2(d[d$id %in% multi,]); a_sin <- lrt2(d[d$id %in% single,])
cat(sprintf("chi2 totale        = %.6f\n", a_all[1]))
cat(sprintf("chi2 sui 22 multi  = %.6f  (quota = %.4f%%)\n", a_mul[1], 100*a_mul[1]/a_all[1]))
cat(sprintf("chi2 sui 185 single= %.8f\n", a_sin[1]))
cat(sprintf("p boundary 50:50 (tot) = %.3e   p naive chi2_1 = %.3e\n",
            0.5*pchisq(a_all[1],1,lower.tail=FALSE), pchisq(a_all[1],1,lower.tail=FALSE)))

## ---------- 6. decomposizione della varianza (R3-V12 / R3-D01) ----------
cat("\n=== R3-V12: ICC per specificazione ===\n")
specs <- list(
  "region 3 livelli"      = lus_severe ~ region,
  "posterior binario"     = lus_severe ~ post_bin,
  "region lineare"        = lus_severe ~ region_lin,
  "section 10 livelli"    = lus_severe ~ section,
  "solo intercetta"       = lus_severe ~ 1,
  "region+side"           = lus_severe ~ region + side,
  "region+dis+inv+age"    = lus_severe ~ region + disease + invasive + age_weeks,
  "region+age_weeks"      = lus_severe ~ region + age_weeks)
S <- 3.289868
out <- list()
for (nm in names(specs)) {
  f <- specs[[nm]]
  fid <- update(f, . ~ . + (1|id)); fpt <- update(f, . ~ . + (1|pt))
  fbo <- update(f, . ~ . + (1|id) + (1|pt))
  g <- function(ff) try(glmer(ff, family=binomial, data=d,
        control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))), silent=TRUE)
  m_id <- g(fid); m_pt <- g(fpt); m_bo <- g(fbo)
  v <- function(m, lab) if (inherits(m,"try-error")) NA_real_ else {
    x <- unlist(VarCorr(m)); vv <- x[grep(lab, names(x))]; if (!length(vv)) NA_real_ else unname(vv[1]) }
  vid1 <- v(m_id,"id"); vpt1 <- v(m_pt,"pt"); vidB <- v(m_bo,"id"); vptB <- v(m_bo,"pt")
  ch <- if (inherits(m_bo,"try-error")||inherits(m_id,"try-error")) NA_real_ else
        as.numeric(2*(logLik(m_bo)-logLik(m_id)))
  out[[nm]] <- data.table(spec=nm, var_id_only=vid1, var_pt_only=vpt1,
    var_id_both=vidB, var_pt_both=vptB,
    icc_pt_A = vptB/(vptB+S),                 # var_pt / (var_pt + pi^2/3)
    icc_pt_B = vptB/(vidB+vptB+S),            # var_pt / totale
    icc_pt_C = (vidB+vptB)/(vidB+vptB+S),     # correlazione entro pt
    icc_pt_only = vpt1/(vpt1+S),
    icc_id_only = vid1/(vid1+S),
    icc_id_both = vidB/(vidB+vptB+S),
    chi2_pt = ch, p_naive = pchisq(ch,1,lower.tail=FALSE))
}
O <- rbindlist(out)
fwrite(O, "rev2/T2_icc.csv")
print(O[, .(spec, var_id_only, var_pt_both, var_id_both, icc_pt_A, icc_pt_B, icc_pt_only, icc_id_only, chi2_pt)], digits=6)
cat("\nchi2/p per la specificazione con covariate:\n")
print(O[spec=="region+dis+inv+age", .(chi2_pt, p_naive, p_bound=0.5*p_naive)], digits=6)
cat("\nintervallo ICC livello paziente (icc_id_only):", sprintf("%.6f - %.6f", min(O$icc_id_only,na.rm=TRUE), max(O$icc_id_only,na.rm=TRUE)), "\n")
cat("0.403 dentro l'intervallo:", 0.403 >= min(O$icc_id_only,na.rm=TRUE) & 0.403 <= max(O$icc_id_only,na.rm=TRUE), "\n")

## ---------- 7. fit a coorte piena: specificazione piena vs ridotta ----------
cat("\n=== fit coorte piena: contrasto posterior vs anterior ===\n")
TERM <- "regionposterior"
fitone <- function(dat, f, tag) {
  dat <- dat %>% arrange(id) %>% mutate(id=droplevels(factor(id)),
      across(where(is.factor), droplevels))
  r <- list()
  m <- try(mm(f, t.formula = ~1, data=as.data.frame(dat), id=id, verbose=FALSE), silent=TRUE)
  if (!inherits(m,"try-error") && TERM %in% names(m$beta)) {
    V <- if (!is.null(m$rob.cov)) m$rob.cov else m$mod.cov
    r[[1]] <- data.table(tag=tag, method="MMM", est=unname(m$beta[TERM]), se=sqrt(V[TERM,TERM]))
  } else r[[1]] <- data.table(tag=tag, method="MMM", est=NA_real_, se=NA_real_)
  g <- try(geeglm(f, family=binomial, id=id, corstr="exchangeable", data=as.data.frame(dat)), silent=TRUE)
  if (!inherits(g,"try-error")) { s <- summary(g)$coefficients
    r[[2]] <- data.table(tag=tag, method="GEE", est=s[TERM,1], se=s[TERM,2]) }
  q <- try(suppressWarnings(glmer(update(f, . ~ . + (1|id)), family=binomial, data=dat,
        control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))), silent=TRUE)
  if (!inherits(q,"try-error")) { s <- summary(q)$coefficients
    r[[3]] <- data.table(tag=tag, method="GLMM", est=s[TERM,1], se=s[TERM,2]) }
  rbindlist(r)
}
dB <- d %>% filter(region != "lateral") %>% mutate(region=droplevels(region))
cat("subset A (tutte le righe):", nrow(d), "obs,", n_distinct(d$id), "pazienti\n")
cat("subset B (solo A$/P$):", nrow(dB), "obs,", n_distinct(dB$id), "pazienti\n")
FS <- list(pieno = lus_severe ~ region + disease + invasive + age_weeks,
           ridotto_region = lus_severe ~ region,
           ridotto_region_age = lus_severe ~ region + age_weeks)
allf <- list(); k <- 0L
for (sub in c("A","B")) for (fn in names(FS)) {
  k <- k+1L
  allf[[k]] <- fitone(if (sub=="A") d else dB, FS[[fn]], paste0(sub,":",fn))
}
FF <- rbindlist(allf)
fwrite(FF, "rev2/T2_fullcohort.csv")
print(FF[order(tag,method)], digits=6)
cat("\nDONE_T2\n")
