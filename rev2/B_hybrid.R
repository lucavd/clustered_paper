
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(dplyr); library(tidyr); library(parallel)})
CUR <- path.expand("~/MMM/results_full")
SUP <- path.expand("~/MMM/results_superseded_m3m5")
OUT <- path.expand("~/MMM/rev2")
DGPS <- c(GenBinaryY="GenBinaryY", GLMMnative="GLMM-native", GLMMAR1="GLMM-AR1")
WALD <- c("mmm","gee","gee_fg","gee_md","gee_kc","glmer","shrink"); LRT <- "glmer_lrt"
keep_rows <- function(d) bind_rows(
  d %>% filter(method %in% WALD, converged, is.finite(estimate), is.finite(se),
               abs(estimate) < 10, se > 0, se < 10, is.finite(p)),
  d %>% filter(method %in% LRT, converged, is.finite(p)))
rd <- function(x) { z <- readRDS(x); z[z$term=="time:trt",
   c("estimate","se","p","method","N","ntime","label","model","converged")] }

pw <- list()
for (pref in names(DGPS)) {
  lab <- DGPS[[pref]]
  ## models 1-2 from the CURRENT tree, models 3-5 from the QUARANTINED tree
  f1 <- list.files(CUR, pattern=paste0("^",pref,"_m[12]_"), full.names=TRUE)
  f2 <- list.files(SUP, pattern=paste0("^",pref,"_m[345]_"), full.names=TRUE)
  comp <- function(fs) { k <- sub("_c[0-9]+[.]rds$","",basename(fs))
                         fs[k %in% names(which(table(k)==10))] }
  fs <- c(comp(f1), comp(f2))
  cat(lab, "hybrid files:", length(comp(f1)), "(m1-2 current) +", length(comp(f2)),
      "(m3-5 superseded) =", length(fs), "\n")
  d <- bind_rows(mclapply(fs, rd, mc.cores=20)); d$model <- as.integer(d$model)
  H0 <- keep_rows(d %>% filter(label=="manual")); H1 <- keep_rows(d %>% filter(label=="mmlb"))
  cut <- H0 %>% group_by(model,N,ntime,method) %>%
    summarise(q05 = quantile(p, 0.05, na.rm=TRUE), .groups="drop")
  pw[[lab]] <- H1 %>% inner_join(cut, by=c("model","N","ntime","method")) %>%
    group_by(method, N) %>%
    summarise(n=n(), power_raw=100*mean(p<0.05), power_adj=100*mean(p<q05),
              mcse=100*sqrt(mean(p<q05)*(1-mean(p<q05))/n()), .groups="drop") %>%
    mutate(dgp=lab)
  rm(d,H0,H1); gc(verbose=FALSE)
}
HY <- bind_rows(pw)
write.csv(HY, file.path(OUT,"R2_power_hybrid_repro.csv"), row.names=FALSE)

PUB <- read.csv(path.expand("~/MMM/tables/P_power_sizeadj.csv"), stringsAsFactors=FALSE)
cmp <- PUB %>% rename(n_pub=n, raw_pub=power_raw, adj_pub=power_adj) %>%
  inner_join(HY %>% rename(n_hy=n, raw_hy=power_raw, adj_hy=power_adj),
             by=c("method","N","dgp"))
cat("\n=== ROW-LEVEL match, published P_power_sizeadj.csv vs hybrid recomputation ===\n")
cat("published rows:", nrow(PUB), " matched rows:", nrow(cmp), "\n")
cat("max |adj diff| =", format(max(abs(cmp$adj_pub-cmp$adj_hy)), digits=6),
    "  max |raw diff| =", format(max(abs(cmp$raw_pub-cmp$raw_hy)), digits=6),
    "  max |n diff| =", max(abs(cmp$n_pub-cmp$n_hy)), "\n")
cat("rows identical to 1e-9 on power_adj:", sum(abs(cmp$adj_pub-cmp$adj_hy)<1e-9),
    "/", nrow(cmp), "\n")

poolf <- function(x) x %>% group_by(dgp, method) %>%
  summarise(nH1=sum(n), raw=weighted.mean(power_raw,n), adj=weighted.mean(power_adj,n),
            .groups="drop")
PP <- poolf(PUB) %>% rename(raw_pub=raw, adj_pub=adj, n_pub=nH1)
PH <- poolf(HY)  %>% rename(raw_hy=raw,  adj_hy=adj,  n_hy=nH1)
J <- PP %>% inner_join(PH, by=c("dgp","method"))
cat("\n=== 24 POOLED rows: published (from tables/) vs hybrid ===\n")
print(as.data.frame(J %>% mutate(adj_pub=round(adj_pub,4), adj_hy=round(adj_hy,4),
      d_adj=round(adj_pub-adj_hy,10), raw_pub=round(raw_pub,3), raw_hy=round(raw_hy,3)) %>%
      select(dgp,method,n_pub,n_hy,raw_pub,raw_hy,adj_pub,adj_hy,d_adj) %>%
      arrange(dgp,method)), row.names=FALSE)
cat("\nall 24 pooled adj identical to 1e-9:", all(abs(J$adj_pub-J$adj_hy)<1e-9), "\n")
cat("hybrid spread per dgp:\n")
print(as.data.frame(PH %>% group_by(dgp) %>% summarise(min=round(min(adj_hy),3),
      max=round(max(adj_hy),3))), row.names=FALSE)
cat("hybrid glmer/gee_fg:\n")
print(as.data.frame(PH %>% filter(method %in% c("glmer","gee_fg")) %>%
      mutate(adj_hy=round(adj_hy,4), raw_hy=round(raw_hy,4))), row.names=FALSE)
cat("\nDONE_B\n")
