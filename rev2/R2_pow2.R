# R2_pow2.R -- riparazione del pooling (n mascherato dentro summarise) e verdetto su §3.11.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(dplyr); library(tidyr); library(data.table)})
WALD <- c("mmm","gee","gee_fg","gee_md","gee_kc","glmer","shrink"); LRT <- "glmer_lrt"
keep_rows <- function(d) bind_rows(
  d %>% filter(method %in% WALD, converged, is.finite(estimate), is.finite(se),
               abs(estimate) < 10, se > 0, se < 10, is.finite(p)),
  d %>% filter(method %in% LRT, converged, is.finite(p)))
DGPS <- c(GenBinaryY="GenBinaryY", GLMMnative="GLMM-native", GLMMAR1="GLMM-AR1")
run_pow <- function(RES, models) {
  pw <- list()
  for (pref in names(DGPS)) {
    lab <- DGPS[[pref]]
    fs <- list.files(RES, pattern=paste0("^",pref,"_"), full.names=TRUE)
    if (!length(fs)) next
    k <- sub("_c[0-9]+[.]rds$","",basename(fs)); f <- names(which(table(k)==10))
    fs <- fs[k %in% f]; if (!length(fs)) next
    d <- bind_rows(lapply(fs, readRDS)); d$model <- as.integer(d$model)
    d <- d %>% filter(term=="time:trt", model %in% models)
    if (!nrow(d)) next
    H0 <- keep_rows(d %>% filter(label=="manual")); H1 <- keep_rows(d %>% filter(label=="mmlb"))
    cut <- H0 %>% group_by(model, N, ntime, method) %>%
      summarise(q05 = quantile(p, 0.05, na.rm=TRUE), .groups="drop")
    pw[[lab]] <- H1 %>% inner_join(cut, by=c("model","N","ntime","method")) %>%
      group_by(method, N) %>%
      summarise(nrep=n(), power_raw=100*mean(p<0.05), power_adj=100*mean(p<q05), .groups="drop") %>%
      mutate(dgp=lab)
    rm(d,H0,H1); gc(verbose=FALSE)
  }
  as.data.table(bind_rows(pw))
}
pool <- function(x) x[, .(n=sum(nrep), raw=round(sum(power_raw*nrep)/sum(nrep),2),
                          adj=round(sum(power_adj*nrep)/sum(nrep),2)), by=.(dgp,method)]
pub <- fread(path.expand("~/MMM/tables/P_power_sizeadj.csv"))
pubp <- pub[, .(n_pub=sum(n), raw_pub=round(sum(power_raw*n)/sum(n),2),
                adj_pub=round(sum(power_adj*n)/sum(n),2)), by=.(dgp,method)]
cmp <- function(P, lab) {
  m <- merge(P, pubp, by=c("dgp","method"))
  m[, `:=`(d_adj=round(adj-adj_pub,2), d_raw=round(raw-raw_pub,2), d_n=n-n_pub)]
  cat(sprintf("\n--- %s ---\n", lab))
  print(as.data.frame(m[order(dgp,-adj), .(dgp,method,n,adj,adj_pub,d_adj,raw,raw_pub,d_raw,d_n)]),
        row.names=FALSE)
  cat(sprintf("max |d_adj| %.2f | max |d_raw| %.2f | max |d_n| %d | somma |d_n| %d\n",
    max(abs(m$d_adj)), max(abs(m$d_raw)), max(abs(m$d_n)), sum(abs(m$d_n))))
  m
}
A <- pool(run_pow(path.expand("~/MMM/results_full"), 1:5))
mA <- cmp(A, "A. chunk ATTUALI, scenari 1-5")
QD <- path.expand("~/MMM/results_superseded_m3m5")
D <- pool(rbind(run_pow(path.expand("~/MMM/results_full"), 1:2), run_pow(QD, 3:5)))
mD <- cmp(D, "D. RICOSTRUZIONE: 1-2 attuali + 3-5 quarantenati")
cat("\n############ verdetto\n")
cat(sprintf("conteggi di replica: attuale scarto totale %d ; ricostruzione scarto totale %d\n",
  sum(abs(mA$d_n)), sum(abs(mD$d_n))))
cat(sprintf("potenza aggiustata: attuale max |scarto| %.2f pp ; ricostruzione max |scarto| %.2f pp\n",
  max(abs(mA$d_adj)), max(abs(mD$d_adj))))
cat("\n############ la cifra di §3.11 su GenBinaryY nelle tre versioni\n")
for (nm in c("attuale","ricostruzione","pubblicato")) {
  x <- switch(nm, attuale=A, ricostruzione=D,
              pubblicato=pubp[, .(dgp,method,adj=adj_pub)])
  y <- x[dgp=="GenBinaryY"]
  cat(sprintf("  %-14s glmer %.2f | gee_fg %.2f | gap %.2f pp | spread %.2f-%.2f\n", nm,
    y[method=="glmer"]$adj, y[method=="gee_fg"]$adj,
    y[method=="glmer"]$adj-y[method=="gee_fg"]$adj, min(y$adj), max(y$adj)))
}
cat("\nspread sugli altri due meccanismi:\n")
for (dg in c("GLMM-native","GLMM-AR1")) {
  cat(sprintf("  %-12s attuale %.2f-%.2f | ricostruzione %.2f-%.2f | pubblicato %.2f-%.2f\n", dg,
    min(A[dgp==dg]$adj), max(A[dgp==dg]$adj), min(D[dgp==dg]$adj), max(D[dgp==dg]$adj),
    min(pubp[dgp==dg]$adj_pub), max(pubp[dgp==dg]$adj_pub)))
}
fwrite(mA, path.expand("~/MMM/rev2/R2_power_current.csv"))
fwrite(mD, path.expand("~/MMM/rev2/R2_power_reconstruction.csv"))
cat("\nDONE_POW2\n")
