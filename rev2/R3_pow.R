# R3_pow.R -- GIRO 3 punto 3: post_gan.R RIESEGUITO sui chunk attuali. Tabella §3.11 pronta.
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(dplyr); library(tidyr); library(data.table)})
RES <- path.expand("~/MMM/results_full")
DGPS <- c(GenBinaryY="GenBinaryY", GLMMnative="GLMM-native", GLMMAR1="GLMM-AR1")
WALD <- c("mmm","gee","gee_fg","gee_md","gee_kc","glmer","shrink"); LRT <- "glmer_lrt"
keep_rows <- function(d) bind_rows(
  d %>% filter(method %in% WALD, converged, is.finite(estimate), is.finite(se),
               abs(estimate) < 10, se > 0, se < 10, is.finite(p)),
  d %>% filter(method %in% LRT, converged, is.finite(p)))
pw <- list(); pwcell <- list()
for (pref in names(DGPS)) {
  lab <- DGPS[[pref]]
  fs <- list.files(RES, pattern=paste0("^",pref,"_"), full.names=TRUE)
  k <- sub("_c[0-9]+[.]rds$","",basename(fs)); f <- names(which(table(k)==10))
  d <- bind_rows(lapply(fs[k %in% f], readRDS)); d$model <- as.integer(d$model)
  d <- d %>% filter(term == "time:trt")
  H0 <- keep_rows(d %>% filter(label == "manual")); H1 <- keep_rows(d %>% filter(label == "mmlb"))
  cut <- H0 %>% group_by(model, N, ntime, method) %>%
    summarise(q05 = quantile(p, 0.05, na.rm = TRUE), .groups = "drop")
  J <- H1 %>% inner_join(cut, by = c("model","N","ntime","method"))
  pw[[lab]] <- J %>% group_by(method, N) %>%
    summarise(n = n(), power_raw = 100*mean(p < 0.05), power_adj = 100*mean(p < q05),
              mcse = 100*sqrt(mean(p<q05)*(1-mean(p<q05))/n()), .groups="drop") %>% mutate(dgp = lab)
  pwcell[[lab]] <- J %>% group_by(method, model, N, ntime) %>%
    summarise(n=n(), power_raw=100*mean(p<0.05), power_adj=100*mean(p<q05), .groups="drop") %>%
    mutate(dgp=lab)
  rm(d,H0,H1,J); gc(verbose=FALSE)
}
PW <- as.data.table(bind_rows(pw)); PC <- as.data.table(bind_rows(pwcell))
fwrite(PW, path.expand("~/MMM/rev2/P_power_sizeadj_CORRECTED.csv"))
fwrite(PC, path.expand("~/MMM/rev2/P_power_sizeadj_CORRECTED_bycell.csv"))

cat("########## TABELLA §3.11 CORRETTA (chunk attuali) ##########\n")
cat("\n=== potenza aggiustata per metodo x N, per meccanismo (formato della tabella pubblicata) ===\n")
for (dg in c("GenBinaryY","GLMM-native","GLMM-AR1")) {
  cat(sprintf("\n-- %s --\n", dg))
  w <- dcast(PW[dgp==dg], method~N, value.var="power_adj")
  nm <- names(w)[-1]
  for (j in nm) w[[j]] <- round(w[[j]],1)
  ord <- PW[dgp==dg, .(a=sum(power_adj*n)/sum(n)), by=method][order(-a)]$method
  print(as.data.frame(w[match(ord, w$method)]), row.names=FALSE)
}
cat("\n=== pooled per meccanismo: potenza grezza e aggiustata, pesate per replica ===\n")
P <- PW[, .(n=sum(n), raw=round(sum(power_raw*n)/sum(n),2), adj=round(sum(power_adj*n)/sum(n),2)),
        by=.(dgp,method)]
print(as.data.frame(P[order(dgp,-adj)]), row.names=FALSE)
cat("\n=== spread size-adjusted per meccanismo ===\n")
S <- P[, .(min_adj=min(adj), max_adj=max(adj), best=method[which.max(adj)], worst=method[which.min(adj)],
           min_raw=min(raw), max_raw=max(raw)), by=dgp]
print(as.data.frame(S), row.names=FALSE)
cat("\n=== divario GLMM - Fay-Graubard, per meccanismo ===\n")
G <- P[method %in% c("glmer","gee_fg"), .(glmer=adj[method=="glmer"], gee_fg=adj[method=="gee_fg"]), by=dgp]
G[, gap := round(glmer-gee_fg,2)]
print(as.data.frame(G), row.names=FALSE)
cat("\n=== confronto con i valori PUBBLICATI ===\n")
pub <- fread(path.expand("~/MMM/tables/P_power_sizeadj.csv"))
pubp <- pub[, .(raw_pub=round(sum(power_raw*n)/sum(n),2), adj_pub=round(sum(power_adj*n)/sum(n),2),
                n_pub=sum(n)), by=.(dgp,method)]
m <- merge(P, pubp, by=c("dgp","method"))
m[, `:=`(d_adj=round(adj-adj_pub,2), d_raw=round(raw-raw_pub,2), d_n=n-n_pub)]
print(as.data.frame(m[order(dgp,-adj), .(dgp,method,n,adj,adj_pub,d_adj,raw,raw_pub,d_raw,d_n)]), row.names=FALSE)
cat(sprintf("\nGenBinaryY: spread %.2f-%.2f (pubblicato 21,66-26,27) ; gap GLMM-FG %.2f pp (pubblicato 3,98)\n",
  S[dgp=="GenBinaryY"]$min_adj, S[dgp=="GenBinaryY"]$max_adj, G[dgp=="GenBinaryY"]$gap))
cat(sprintf("GLMM-native: spread %.2f-%.2f (pubblicato 23,14-24,55)\n",
  S[dgp=="GLMM-native"]$min_adj, S[dgp=="GLMM-native"]$max_adj))
cat(sprintf("GLMM-AR1: spread %.2f-%.2f (pubblicato 23,00-25,22)\n",
  S[dgp=="GLMM-AR1"]$min_adj, S[dgp=="GLMM-AR1"]$max_adj))
cat("\n=== l'ordinamento dei metodi cambia? ===\n")
for (dg in c("GenBinaryY","GLMM-native","GLMM-AR1")) {
  a <- P[dgp==dg][order(-adj)]$method; b <- pubp[dgp==dg][order(-adj_pub)]$method
  cat(sprintf("%-12s attuale: %s\n%-12s pubbl.: %s\n", dg, paste(a,collapse=" > "), "", paste(b,collapse=" > ")))
  cat(sprintf("             Spearman fra i due ordinamenti: %.3f\n",
    cor(match(a,a), match(a,b), method="spearman")))
}
cat("\nDONE_POW3\n")
