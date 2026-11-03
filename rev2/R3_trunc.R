# R3_trunc.R -- GIRO 3 punto 5: troncamento e identificabilita' a T=2.
# La domanda: il GLMM a pendenza casuale non fitta per il TRONCAMENTO o per la SATURAZIONE?
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(data.table); library(lme4); library(parallel); library(MASS)})
RD <- path.expand("~/MMM/.rundir")
sys.source(file.path(RD,"gen_glmm_data.R"), envir=globalenv())
for (f in c("manual/generation_model3.R","manual/generate_data.R","manual/make_all_data.R"))
  sys.source(file.path(RD,f), envir=globalenv())
CTRL <- lme4::glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e5))

cat("########## A. composizione dei cluster a T=2, misurata\n")
cnt <- rbindlist(lapply(1:2000, function(i) {
  set.seed(20260805L+3000L+i); N <- 200
  nc <- ifelse(runif(N) < 0.7, 2L, sample(seq_len(1), N, replace=TRUE))
  data.table(rep=i, n1=sum(nc==1), n2=sum(nc==2), N=N)
}))
cat(sprintf("su 2000 repliche x 200 cluster: quota a 1 osservazione %.4f, a 2 osservazioni %.4f\n",
  sum(cnt$n1)/sum(cnt$N), sum(cnt$n2)/sum(cnt$N)))
cat("(a T=2 il ramo troncato campiona da seq_len(T-1) = {1}, quindi 30% esattamente a 1 osservazione)\n")

cat("\n########## B. identificabilita': quanti parametri a livello di soggetto sono stimabili?\n")
for (N in c(10,20,45)) {
  # composizione attesa
  n2 <- 0.7*N; n1 <- 0.3*N; nobs <- 2*n2 + n1
  cat(sprintf("N=%d: cluster a 2 oss %.1f, a 1 oss %.1f, osservazioni %.1f\n", N, n2, n1, nobs))
  cat(sprintf("   effetti casuali richiesti da (1+time|id): 2N = %d ; osservazioni disponibili %.1f -> deficit %.1f\n",
    N, 2*N, nobs, 2*N - nobs))
  cat(sprintf("   contributi ENTRO cluster (cluster con >=2 oss): %.1f su %d cluster (%.0f%%)\n",
    n2, N, 100*n2/N))
  cat(sprintf("   gradi di liberta' entro-cluster: nobs - N - 2 = %.1f\n", nobs - N - 2))
}
cat("\nsenza troncamento (tutti a 2 oss): nobs = 2N, effetti casuali 2N, deficit 0 -> saturazione ESATTA\n")
cat("con troncamento (70/30): nobs = 1.7N, effetti casuali 2N, deficit 0.3N -> nobs < numero di effetti casuali\n")

cat("\n########## C. il GLMM a pendenza casuale fitta nelle repliche SENZA troncamento?\n")
# genero forzando tutti i cluster a T osservazioni, cioe' rimuovo il troncamento
gen_full <- function(N, T, i, label) {
  set.seed(20260805L+1000L*3+i)
  bc <- c(-2.32, 0.116, 0.232, if (label=="manual") 0 else 0.232)
  id <- rep(seq_len(N), each=T); tt <- rep(seq_len(T)-1, N)
  trt <- rep(rbinom(N,1,0.5), each=T)
  dat <- data.frame(id=id, time=tt, trt=trt)
  D2 <- matrix(c(1,.01,.01,.001),2,byrow=TRUE)
  gen_glmm_data(~time*trt, ~1+time, beta_c=bc, D=D2, id=id, data=dat, Yname="binY")
}
try_fit <- function(d, re) {
  d <- as.data.frame(d); d$id <- factor(d$id)
  f <- as.formula(paste("binY ~ time*trt +", re))
  m <- tryCatch(suppressMessages(glmer(f, data=d, family=binomial, control=CTRL)),
                error=function(e) structure(conditionMessage(e), class="ee"))
  if (inherits(m,"ee")) return(list(ok=0, why=substr(m,1,64), est=NA, se=NA))
  s <- tryCatch(summary(m)$coefficients, error=function(e) NULL)
  if (is.null(s) || !"time:trt" %in% rownames(s)) return(list(ok=0, why="no_coef", est=NA, se=NA))
  e <- s["time:trt",1]; v <- s["time:trt",2]
  ok <- as.integer(is.finite(e) & is.finite(v) & v>0 & v<10 & abs(e)<10)
  list(ok=ok, why=ifelse(ok==1,"ok","filtered"), est=e, se=v)
}
NREP <- 400L
res <- rbindlist(lapply(c(10,20,45), function(N) {
  rr <- mclapply(1:NREP, function(i) {
    dtr <- tryCatch(as.data.frame(make_all_data(N,2,i,"3","manual")), error=function(e) NULL)  # con troncamento
    dfu <- tryCatch(as.data.frame(gen_full(N,2,i,"manual")), error=function(e) NULL)           # senza
    if (is.null(dtr) || is.null(dfu)) return(NULL)
    a <- try_fit(dtr, "(1 + time|id)"); b <- try_fit(dfu, "(1 + time|id)")
    a2 <- try_fit(dtr, "(1|id)");       b2 <- try_fit(dfu, "(1|id)")
    data.table(N=N, rep=i, nobs_tr=nrow(dtr), nobs_full=nrow(dfu),
      rs_tr=a$ok, rs_tr_why=a$why, rs_full=b$ok, rs_full_why=b$why,
      ri_tr=a2$ok, ri_full=b2$ok)
  }, mc.cores=20)
  rbindlist(Filter(Negate(is.null), rr))
}))
saveRDS(res, path.expand("~/MMM/rev2/R3_trunc_raw.rds"))
cat("\n===== GLMM a pendenza casuale (1+time|id) a T=2: con e senza troncamento =====\n")
S <- res[, .(tentativi=.N, nobs_tr=round(mean(nobs_tr),1), nobs_full=round(mean(nobs_full),1),
             usable_tr=sum(rs_tr), pct_tr=round(100*mean(rs_tr),2),
             usable_full=sum(rs_full), pct_full=round(100*mean(rs_full),2)), by=N]
print(as.data.frame(S), row.names=FALSE)
cat("\n===== la stessa cosa con la sola intercetta casuale (controllo) =====\n")
S2 <- res[, .(tentativi=.N, ri_tr=round(100*mean(ri_tr),2), ri_full=round(100*mean(ri_full),2)), by=N]
print(as.data.frame(S2), row.names=FALSE)
cat("\n===== cause di fallimento senza troncamento =====\n")
print(as.data.frame(res[rs_full==0, .N, by=.(N, rs_full_why)][order(N,-N)]), row.names=FALSE)
cat("\n===== cause di fallimento con troncamento =====\n")
print(as.data.frame(res[rs_tr==0, .N, by=.(N, rs_tr_why)][order(N,-N)]), row.names=FALSE)
cat(sprintf("\nVERDETTO: pendenza casuale a T=2, con troncamento %.2f%% utilizzabili (%d su %d) ; senza troncamento %.2f%% (%d su %d)\n",
  100*mean(res$rs_tr), sum(res$rs_tr), nrow(res), 100*mean(res$rs_full), sum(res$rs_full), nrow(res)))
cat(sprintf("intercetta casuale: con troncamento %.2f%% ; senza %.2f%%\n",
  100*mean(res$ri_tr), 100*mean(res$ri_full)))
fwrite(S, path.expand("~/MMM/rev2/R3_trunc_glmm.csv"))
cat("\nDONE_TRUNC\n")
