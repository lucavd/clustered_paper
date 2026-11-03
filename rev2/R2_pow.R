# R2_pow.R -- GIRO 2, contestazione 5: post_gan.R E' il codice della potenza. Riproduzione
# esatta con la sua convenzione (quantile type 7, cut per model+N+ntime+method, pooling per N),
# sui chunk ATTUALI e sui chunk QUARANTENATI (results_superseded_m3m5).
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(dplyr); library(tidyr); library(data.table)})
OUTD <- path.expand("~/MMM/rev2")
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
      summarise(q05 = quantile(p, 0.05, na.rm=TRUE), .groups="drop")   # type 7, come post_gan.R
    pw[[lab]] <- H1 %>% inner_join(cut, by=c("model","N","ntime","method")) %>%
      group_by(method, N) %>%
      summarise(n=n(), power_raw=100*mean(p<0.05), power_adj=100*mean(p<q05), .groups="drop") %>%
      mutate(dgp=lab)
    rm(d,H0,H1); gc(verbose=FALSE)
  }
  bind_rows(pw)
}
pool <- function(x) x %>% group_by(dgp, method) %>%
  summarise(n=sum(n), raw=round(sum(power_raw*n)/sum(n),2), adj=round(sum(power_adj*n)/sum(n),2),
            .groups="drop")
pub <- fread(path.expand("~/MMM/tables/P_power_sizeadj.csv"))
pubp <- pub[, .(n_pub=sum(n), raw_pub=round(sum(power_raw*n)/sum(n),2),
                adj_pub=round(sum(power_adj*n)/sum(n),2)), by=.(dgp,method)]

cat("############ A. chunk ATTUALI, tutti i cinque scenari, convenzione di post_gan.R\n")
A <- pool(run_pow(path.expand("~/MMM/results_full"), 1:5))
m <- merge(as.data.table(A), pubp, by=c("dgp","method"))
m[, `:=`(d_adj=round(adj-adj_pub,2), d_raw=round(raw-raw_pub,2), d_n=n-n_pub)]
print(as.data.frame(m[order(dgp,-adj)]), row.names=FALSE)
cat(sprintf("max |d_adj| %.2f | max |d_raw| %.2f | max |d_n| %d\n",
  max(abs(m$d_adj)), max(abs(m$d_raw)), max(abs(m$d_n))))

cat("\n############ B. SOLO scenari 1-2 (non rigenerati): la tabella combacia qui?\n")
B <- pool(run_pow(path.expand("~/MMM/results_full"), 1:2))
print(as.data.frame(B[order(B$dgp,-B$adj),]), row.names=FALSE)

cat("\n############ C. scenari 3-5 dai chunk QUARANTENATI (disegno pre-correzione)\n")
QD <- path.expand("~/MMM/results_superseded_m3m5")
cat("file quarantenati:", length(list.files(QD, pattern="[.]rds$")), "\n")
C <- pool(run_pow(QD, 3:5))
print(as.data.frame(C[order(C$dgp,-C$adj),]), row.names=FALSE)

cat("\n############ D. RICOSTRUZIONE: scenari 1-2 attuali + scenari 3-5 quarantenati\n")
D <- bind_rows(run_pow(path.expand("~/MMM/results_full"), 1:2), run_pow(QD, 3:5))
Dp <- pool(D)
md <- merge(as.data.table(Dp), pubp, by=c("dgp","method"))
md[, `:=`(d_adj=round(adj-adj_pub,2), d_raw=round(raw-raw_pub,2), d_n=n-n_pub)]
print(as.data.frame(md[order(dgp,-adj)]), row.names=FALSE)
cat(sprintf("max |d_adj| %.2f | max |d_raw| %.2f | max |d_n| %d\n",
  max(abs(md$d_adj)), max(abs(md$d_raw)), max(abs(md$d_n))))
cat("\nse D combacia e A no, §3.11 e' calcolata sul disegno PRE-correzione per M3-M5.\n")
fwrite(md, file.path(OUTD,"R2_power_reconstruction.csv"))
fwrite(m,  file.path(OUTD,"R2_power_current.csv"))

cat("\n############ E. la cifra di §3.11: GLMM vs FG su GenBinaryY nelle tre versioni\n")
for (nm in c("attuale (1-5)","ricostruzione (1-2 attuali + 3-5 quarantenati)")) {
  x <- if (nm=="attuale (1-5)") as.data.table(A) else as.data.table(Dp)
  g <- x[dgp=="GenBinaryY" & method=="glmer"]$adj; fg <- x[dgp=="GenBinaryY" & method=="gee_fg"]$adj
  cat(sprintf("  %-46s glmer %.2f | gee_fg %.2f | gap %.2f pp | spread %.2f-%.2f\n", nm, g, fg, g-fg,
    min(x[dgp=="GenBinaryY"]$adj), max(x[dgp=="GenBinaryY"]$adj)))
}
p <- pubp[dgp=="GenBinaryY"]
cat(sprintf("  %-46s glmer %.2f | gee_fg %.2f | gap %.2f pp | spread %.2f-%.2f\n", "pubblicato",
  p[method=="glmer"]$adj_pub, p[method=="gee_fg"]$adj_pub,
  p[method=="glmer"]$adj_pub-p[method=="gee_fg"]$adj_pub, min(p$adj_pub), max(p$adj_pub)))
cat("\nDONE_POW\n")
