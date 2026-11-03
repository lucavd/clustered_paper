# R3_df.R -- GIRO 3 punto 2: §3.12 ricalcolata con i gdl ENTRO-cluster per time:trt.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(data.table); library(parallel)})
RES <- path.expand("~/MMM/results_full")
WALD <- c("mmm","gee","gee_fg","gee_md","gee_kc","glmer","shrink")
NP <- c(4,5,4,4,5)                      # parametri del modello medio per scenario
# covariate ENTRO cluster nel modello medio: time e time:trt (trt, age, dose sono fra cluster)
NP_WITHIN <- c(2,2,2,2,2)
one <- function(f) {
  x <- try(as.data.table(readRDS(f)), silent=TRUE); if (inherits(x,"try-error")) return(NULL)
  x <- x[term=="time:trt" & method %in% WALD]; if (!nrow(x)) return(NULL)
  x[, model := as.integer(model)]
  x[, converged := as.logical(converged)][is.na(converged), converged := FALSE]
  y <- x[converged & is.finite(estimate) & is.finite(se) & abs(estimate)<10 & se>0 & se<10]
  if (!nrow(y)) return(NULL)
  # nobs per replica: ricavabile solo se la colonna esiste; altrimenti la stimo dal disegno
  cols <- intersect(c("nobs","n_obs","N_obs"), names(y))
  y[, nobs_col := if (length(cols)) get(cols[1]) else NA_real_]
  y[, .(dgp,method,model,N,ntime,label,sim_id,estimate,se,nobs_col)]
}
D <- rbindlist(Filter(Negate(is.null), mclapply(list.files(RES, full.names=TRUE), one, mc.cores=45)), fill=TRUE)
cat("righe:", nrow(D), " colonna nobs nei chunk:", all(is.na(D$nobs_col)), "\n")
# nobs atteso per cella: N * (0.7*T + 0.3*mean(1..T-1)) -- il troncamento e' uniforme in 1..T-1
D[, nobs_exp := N*(0.7*ntime + 0.3*(ntime)/2)]
cat("nobs atteso per T: "); print(unique(D[, .(ntime, per_cluster=round(nobs_exp/N,3))])[order(ntime)])
D[, `:=`(z = estimate/se,
         df_bc  = N - NP[model],                          # cio' che e' stato usato (between-cluster)
         df_wc  = round(nobs_exp - N - NP_WITHIN[model]))] # residuo entro-cluster
D <- D[df_bc >= 1 & df_wc >= 1]
zq <- qnorm(0.975)
D[, `:=`(rej_norm = abs(z) > zq,
         rej_t_bc = abs(z) > qt(0.975, df_bc),
         rej_t_wc = abs(z) > qt(0.975, df_wc))]
cat("\n===== §3.12 con i due sistemi di gdl: tassi per metodo x meccanismo (braccio H0) =====\n")
H <- D[label=="manual"]
T1 <- H[, .(n=.N, norm=100*mean(rej_norm), t_bc=100*mean(rej_t_bc), t_wc=100*mean(rej_t_wc)), by=.(method,dgp)]
T1[, `:=`(shift_bc=round(t_bc-norm,2), shift_wc=round(t_wc-norm,2))]
W <- dcast(T1, method~dgp, value.var=c("norm","t_bc","t_wc","shift_bc","shift_wc"))
ord <- T1[, .(s=mean(shift_bc)), by=method][order(s)]$method
cat("\nformato del manoscritto: GenBinaryY / GLMM-native / GLMM-AR1\n")
for (mth in ord) {
  r <- T1[method==mth]
  g <- function(d,c) round(r[dgp==d][[c]],2)
  cat(sprintf("%-8s normale %5.2f / %5.2f / %5.2f | t (N-p) %5.2f / %5.2f / %5.2f | t entro-cluster %5.2f / %5.2f / %5.2f | shift(N-p) %+.2f / %+.2f / %+.2f | shift(entro) %+.2f / %+.2f / %+.2f\n",
    mth, g("GenBinaryY","norm"), g("GLMM-native","norm"), g("GLMM-AR1","norm"),
    g("GenBinaryY","t_bc"), g("GLMM-native","t_bc"), g("GLMM-AR1","t_bc"),
    g("GenBinaryY","t_wc"), g("GLMM-native","t_wc"), g("GLMM-AR1","t_wc"),
    g("GenBinaryY","shift_bc"), g("GLMM-native","shift_bc"), g("GLMM-AR1","shift_bc"),
    g("GenBinaryY","shift_wc"), g("GLMM-native","shift_wc"), g("GLMM-AR1","shift_wc")))
}
cat("\n===== shift a N<=25 (la colonna che il manoscritto evidenzia) =====\n")
S25 <- H[N<=25, .(n=.N, norm=100*mean(rej_norm), t_bc=100*mean(rej_t_bc), t_wc=100*mean(rej_t_wc)), by=method]
S25[, `:=`(shift_bc=round(t_bc-norm,2), shift_wc=round(t_wc-norm,2))]
print(as.data.frame(S25[order(shift_bc), .(method,n,norm=round(norm,2),t_bc=round(t_bc,2),
  t_wc=round(t_wc,2),shift_bc,shift_wc)]), row.names=FALSE)
cat(sprintf("\nintervallo dello shift a N<=25: (N-p) da %.2f a %.2f pp | entro-cluster da %.2f a %.2f pp\n",
  min(S25$shift_bc), max(S25$shift_bc), min(S25$shift_wc), max(S25$shift_wc)))
cat(sprintf("pooled sulla griglia: (N-p) da %.2f a %.2f pp | entro-cluster da %.2f a %.2f pp\n",
  min(T1$shift_bc), max(T1$shift_bc), min(T1$shift_wc), max(T1$shift_wc)))

cat("\n===== WORST CASE per cella: chi ha il primato sotto i tre riferimenti =====\n")
CELL <- H[, .(n=.N, norm=100*mean(rej_norm), t_bc=100*mean(rej_t_bc), t_wc=100*mean(rej_t_wc)),
          by=.(method,dgp,model,N,ntime)]
WC <- CELL[, .(worst_norm=max(abs(norm-5)), worst_t_bc=max(abs(t_bc-5)), worst_t_wc=max(abs(t_wc-5))),
           by=.(method,dgp)]
for (dg in c("GenBinaryY","GLMM-native","GLMM-AR1")) {
  cat(sprintf("\n-- %s --\n", dg))
  x <- WC[dgp==dg][order(worst_norm)]
  print(as.data.frame(x[, .(method, worst_norm=round(worst_norm,2), worst_t_bc=round(worst_t_bc,2),
    worst_t_wc=round(worst_t_wc,2))]), row.names=FALSE)
  cat(sprintf("   primato: normale -> %s (%.2f) | t (N-p) -> %s (%.2f) | t entro-cluster -> %s (%.2f)\n",
    x[which.min(worst_norm)]$method, min(x$worst_norm),
    x[which.min(worst_t_bc)]$method, min(x$worst_t_bc),
    x[which.min(worst_t_wc)]$method, min(x$worst_t_wc)))
}
cat("\n===== la conclusione di §3.12 sopravvive? mmm contro gee_fg sui tre riferimenti =====\n")
for (dg in c("GenBinaryY","GLMM-native","GLMM-AR1")) {
  a <- WC[dgp==dg & method=="mmm"]; b <- WC[dgp==dg & method=="gee_fg"]
  cat(sprintf("%-12s normale: mmm %.2f vs fg %.2f (%s) | t(N-p): mmm %.2f vs fg %.2f (%s) | t(entro): mmm %.2f vs fg %.2f (%s)\n",
    dg, a$worst_norm, b$worst_norm, ifelse(a$worst_norm<b$worst_norm,"mmm","fg"),
    a$worst_t_bc, b$worst_t_bc, ifelse(a$worst_t_bc<b$worst_t_bc,"mmm","fg"),
    a$worst_t_wc, b$worst_t_wc, ifelse(a$worst_t_wc<b$worst_t_wc,"mmm","fg")))
}
cat("\n===== valori dei gdl a confronto =====\n")
G <- unique(D[, .(model,N,ntime,df_bc,df_wc)])[order(model,N,ntime)]
print(as.data.frame(G[model %in% c(1,5) & N %in% c(10,25,45,145) & ntime %in% c(2,4,10)]), row.names=FALSE)
cat(sprintf("\nquantile t al 97.5%%: a 5 gdl %.4f | a 6 gdl %.4f | a 22 gdl %.4f | a 336 gdl %.4f | normale %.4f\n",
  qt(.975,5), qt(.975,6), qt(.975,22), qt(.975,336), zq))
fwrite(T1, path.expand("~/MMM/rev2/R3_df_bymethod.csv"))
fwrite(WC, path.expand("~/MMM/rev2/R3_df_worstcase.csv"))
fwrite(S25, path.expand("~/MMM/rev2/R3_df_N25.csv"))
cat("\nDONE_DF\n")
