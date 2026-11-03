
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(dplyr); library(tidyr); library(data.table)})
RES <- path.expand("~/MMM/results_full")
DGPS <- c(GenBinaryY="GenBinaryY", GLMMnative="GLMM-native", GLMMAR1="GLMM-AR1")
WALD <- c("mmm","gee","gee_fg","gee_md","gee_kc","glmer","shrink"); LRT <- "glmer_lrt"
keep <- function(d) bind_rows(
  d %>% filter(method %in% WALD, converged, is.finite(estimate), is.finite(se),
               abs(estimate) < 10, se > 0, se < 10, is.finite(p)),
  d %>% filter(method == LRT, converged, is.finite(p)))

cat("=== mtime dei chunk contro la data di F_paired (2026-08-08 10:18) ===\n")
fi <- file.info(list.files(RES, full.names=TRUE, pattern="[.]rds$"))
cat("mtime min:", format(min(fi$mtime)), " max:", format(max(fi$mtime)), "\n")
cat("chunk modificati DOPO 2026-08-08 10:18:", sum(fi$mtime > as.POSIXct("2026-08-08 10:18:00")), "su", nrow(fi), "\n")
cat("F_paired mtime:", format(file.info(path.expand("~/MMM/tables/F_paired_mmm_fg.csv"))$mtime), "\n")

cat("\n=== RIESECUZIONE VERBATIM del blocco pr[[lab]] di fix_tables.R ===\n")
out <- list()
for (pref in names(DGPS)) {
  lab <- DGPS[[pref]]
  fs <- list.files(RES, pattern=paste0("^",pref,"_"), full.names=TRUE)
  k <- sub("_c[0-9]+[.]rds$","",basename(fs)); f <- names(which(table(k)==10))
  d <- bind_rows(lapply(fs[k %in% f], readRDS)); d$model <- as.integer(d$model)
  d <- d %>% filter(term == "time:trt")
  H0 <- keep(d %>% filter(label == "manual"))
  pr <- H0 %>% filter(method %in% c("mmm","gee_fg")) %>%
    select(model, N, ntime, sim_id, method, p) %>%
    mutate(rej = p < 0.05) %>% select(-p) %>%
    pivot_wider(names_from=method, values_from=rej) %>%
    filter(!is.na(mmm), !is.na(gee_fg)) %>%
    summarise(n_matched=n(), b=sum(mmm & !gee_fg), c=sum(!mmm & gee_fg)) %>%
    mutate(dgp=lab, diff_pp=100*(b-c)/n_matched,
           chisq=(abs(b-c)-1)^2/(b+c), p=pchisq(chisq,1,lower.tail=FALSE))
  print(as.data.frame(pr))
  out[[lab]] <- pr
  # diagnostica: chiave senza chunk -> collisioni?
  H0m <- H0 %>% filter(method %in% c("mmm","gee_fg"))
  dup <- H0m %>% count(model, N, ntime, sim_id, method) %>% filter(n > 1)
  cat("  collisioni di chiave (model,N,ntime,sim_id,method) con n>1:", nrow(dup),
      " righe coinvolte:", sum(dup$n), "\n")
  if (nrow(dup)) print(head(as.data.frame(dup),5))
  cat("  righe mmm:", sum(H0m$method=="mmm"), " righe gee_fg:", sum(H0m$method=="gee_fg"), "\n")
  # e con la chiave che include il chunk?
  rm(d,H0,H0m); gc(verbose=FALSE)
}
cat("\n=== confronto con il file pubblicato ===\n")
P <- fread(path.expand("~/MMM/tables/F_paired_mmm_fg.csv"))
print(P)
R <- rbindlist(out)
print(R)
M <- merge(P[, .(dgp, npair_file=n_matched, b_file=b, c_file=c)],
           R[, .(dgp, npair_new=n_matched, b_new=b, c_new=c)], by="dgp")
M[, `:=`(d_npair=npair_new-npair_file, d_b=b_new-b_file, d_c=c_new-c_file,
         d_bc=(b_new+c_new)-(b_file+c_file))]
print(M)
cat("\nDONE_D\n")
