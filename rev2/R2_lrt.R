# R2_lrt.R -- GIRO 2, contestazione 1: riconciliazione della controprova LRT.
# Registra le CAUSE di fallimento e nobs vs numero di effetti casuali, non solo i tassi.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({library(dplyr); library(tidyr); library(lme4); library(parallel); library(MASS)})
RD <- path.expand("~/MMM/.rundir")
sys.source(file.path(RD,"gen_glmm_data.R"), envir=globalenv())
for (f in c("manual/generation_model3.R","manual/generate_data.R","manual/make_all_data.R"))
  sys.source(file.path(RD,f), envir=globalenv())
CTRL <- lme4::glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e5))

# filtro dichiarato per il LRT: converged & is.finite(p). NESSUNA condizione su SE o stima.
lrt <- function(d, re) {
  ff <- as.formula(paste("binY ~ time*trt +", re))
  fr <- as.formula(paste("binY ~ time + trt +", re))
  full <- tryCatch(suppressMessages(glmer(ff, data=d, family=binomial, control=CTRL)),
                   error=function(e) structure(conditionMessage(e), class="ee"))
  if (inherits(full,"ee")) return(list(ok=0, why="full_error", msg=substr(full,1,70), p=NA))
  red <- tryCatch(suppressMessages(glmer(fr, data=d, family=binomial, control=CTRL)),
                  error=function(e) structure(conditionMessage(e), class="ee"))
  if (inherits(red,"ee")) return(list(ok=0, why="reduced_error", msg=substr(red,1,70), p=NA))
  a <- tryCatch(anova(red, full), error=function(e) NULL)
  if (is.null(a)) return(list(ok=0, why="anova_error", msg=NA, p=NA))
  p <- a$`Pr(>Chisq)`[2]
  if (!is.finite(p)) return(list(ok=0, why="p_nonfinite", msg=NA, p=NA))
  list(ok=1, why="ok", msg=NA, p=p)
}

NREP <- 500L; MODEL <- 3L
grid <- expand.grid(N=c(10,20,45), ntime=c(2,3,4), label=c("manual","mmlb"), stringsAsFactors=FALSE)
cell <- function(k) {
  g <- grid[k,]
  rows <- lapply(1:NREP, function(i) {
    set.seed(20260805L + 1000L*MODEL + i)
    d <- tryCatch(make_all_data(g$N, g$ntime, i, as.character(MODEL), g$label), error=function(e) NULL)
    if (is.null(d)) return(NULL)
    d <- as.data.frame(d); d$id <- factor(d$id)
    nobs <- nrow(d); nre_rs <- 2L*nlevels(d$id); nre_ri <- nlevels(d$id)
    a <- lrt(d, "(1|id)"); b <- lrt(d, "(1 + time|id)")
    data.frame(rep=i, nobs=nobs, nre_ri=nre_ri, nre_rs=nre_rs,
               all_full = as.integer(nobs == g$N*g$ntime),
               ok_ri=a$ok, p_ri=a$p, why_ri=a$why,
               ok_rs=b$ok, p_rs=b$p, why_rs=b$why, msg_rs=ifelse(is.na(b$msg),"",b$msg))
  })
  r <- bind_rows(rows)
  r$N <- g$N; r$ntime <- g$ntime; r$label <- g$label
  r
}
res <- mclapply(seq_len(nrow(grid)), function(k) try(cell(k), silent=TRUE), mc.cores=18)
bad <- !sapply(res, is.data.frame); cat("celle fallite:", sum(bad), "\n")
R <- bind_rows(res[!bad])
saveRDS(R, path.expand("~/MMM/rev2/R2_lrt_raw.rds"))

cat("\n===== CONTEGGI GREZZI per (N,T,arm): tentativi, utilizzabili (1|id), utilizzabili (1+time|id) =====\n")
tb <- R %>% group_by(N, ntime, label) %>%
  summarise(tentativi=n(), us_ri=sum(ok_ri), us_rs=sum(ok_rs),
            rej_ri=sum(ok_ri==1 & p_ri<0.05), rej_rs=sum(ok_rs==1 & p_rs<0.05),
            nobs_eq_nre = sum(nobs == nre_rs), all_full=sum(all_full), .groups="drop")
print(as.data.frame(tb), row.names=FALSE)

cat("\n===== CAUSE di non-utilizzabilita con (1+time|id) =====\n")
print(as.data.frame(R %>% filter(ok_rs==0) %>% count(ntime, why_rs, name="repliche")), row.names=FALSE)
cat("\nmessaggi di errore distinti (troncati), con frequenza:\n")
print(as.data.frame(R %>% filter(ok_rs==0, nzchar(msg_rs)) %>% count(msg_rs, name="n") %>% arrange(-n) %>% head(6)),
      row.names=FALSE)

cat("\n===== il meccanismo del T=2: quando nobs raggiunge 2N il fit passa =====\n")
m2 <- R %>% filter(ntime==2) %>% group_by(N) %>%
  summarise(tentativi=n(), nobs_eq_2N=sum(nobs==nre_rs), us_rs=sum(ok_rs),
            us_rs_tra_nobs_eq=sum(ok_rs[nobs==nre_rs]), .groups="drop")
m2$p_teorica_tutti_completi <- round(0.7^m2$N, 5)
m2$pct_us_rs <- round(100*m2$us_rs/m2$tentativi, 3)
print(as.data.frame(m2), row.names=FALSE)
cat(sprintf("\nT=2 pooled sui tre N e i due bracci: utilizzabili (1+time|id) %d su %d = %.3f%% ; con (1|id) %d su %d = %.2f%%\n",
  sum(R$ok_rs[R$ntime==2]), sum(R$ntime==2), 100*mean(R$ok_rs[R$ntime==2]),
  sum(R$ok_ri[R$ntime==2]), sum(R$ntime==2), 100*mean(R$ok_ri[R$ntime==2])))

cat("\n===== LIVELLO (braccio H0) con i conteggi grezzi =====\n")
h0 <- R %>% filter(label=="manual")
lv <- h0 %>% group_by(ntime) %>%
  summarise(n_ri=sum(ok_ri), k_ri=sum(ok_ri==1 & p_ri<0.05),
            n_rs=sum(ok_rs), k_rs=sum(ok_rs==1 & p_rs<0.05), .groups="drop") %>%
  mutate(level_ri=round(100*k_ri/n_ri,2), mcse_ri=round(100*sqrt((k_ri/n_ri)*(1-k_ri/n_ri)/n_ri),2),
         level_rs=ifelse(n_rs>0, round(100*k_rs/n_rs,2), NA),
         mcse_rs=ifelse(n_rs>0, round(100*sqrt((k_rs/n_rs)*(1-k_rs/n_rs)/n_rs),2), NA))
print(as.data.frame(lv), row.names=FALSE)
cat(sprintf("\npooled T>=3, H0: (1|id) %.2f%% (k=%d, n=%d, MCSE %.2f) | (1+time|id) %.2f%% (k=%d, n=%d, MCSE %.2f)\n",
  100*sum(lv$k_ri[lv$ntime>=3])/sum(lv$n_ri[lv$ntime>=3]), sum(lv$k_ri[lv$ntime>=3]), sum(lv$n_ri[lv$ntime>=3]),
  100*sqrt((sum(lv$k_ri[lv$ntime>=3])/sum(lv$n_ri[lv$ntime>=3]))*(1-sum(lv$k_ri[lv$ntime>=3])/sum(lv$n_ri[lv$ntime>=3]))/sum(lv$n_ri[lv$ntime>=3])),
  100*sum(lv$k_rs[lv$ntime>=3])/sum(lv$n_rs[lv$ntime>=3]), sum(lv$k_rs[lv$ntime>=3]), sum(lv$n_rs[lv$ntime>=3]),
  100*sqrt((sum(lv$k_rs[lv$ntime>=3])/sum(lv$n_rs[lv$ntime>=3]))*(1-sum(lv$k_rs[lv$ntime>=3])/sum(lv$n_rs[lv$ntime>=3]))/sum(lv$n_rs[lv$ntime>=3]))))
cat(sprintf("pooled TUTTI i T, H0: (1|id) %.2f%% (n=%d) | (1+time|id) %.2f%% (n=%d)\n",
  100*sum(h0$ok_ri==1 & h0$p_ri<0.05, na.rm=TRUE)/sum(h0$ok_ri), sum(h0$ok_ri),
  100*sum(h0$ok_rs==1 & h0$p_rs<0.05, na.rm=TRUE)/sum(h0$ok_rs), sum(h0$ok_rs)))

cat("\n===== la cifra dell'Abstract: 98,8% a T=2 diventa cosa? =====\n")
cat(sprintf("T=2, struttura dichiarata (1+time|id): %d utilizzabili su %d tentativi = %.3f%% [IC binomiale esatto %.3f-%.3f%%]\n",
  sum(R$ok_rs[R$ntime==2]), sum(R$ntime==2), 100*mean(R$ok_rs[R$ntime==2]),
  100*binom.test(sum(R$ok_rs[R$ntime==2]), sum(R$ntime==2))$conf.int[1],
  100*binom.test(sum(R$ok_rs[R$ntime==2]), sum(R$ntime==2))$conf.int[2]))
for (tt in c(3,4)) cat(sprintf("T=%d, struttura dichiarata: %.2f%% (%d su %d) ; come eseguito %.2f%%\n", tt,
  100*mean(R$ok_rs[R$ntime==tt]), sum(R$ok_rs[R$ntime==tt]), sum(R$ntime==tt),
  100*mean(R$ok_ri[R$ntime==tt])))
cat("\nDONE_LRT\n")
