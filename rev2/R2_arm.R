# R2_arm.R -- GIRO 2, contestazione 3: livello su GenBinaryY M3 con il braccio H0
# (a) come eseguito, sigma=c(1,0)  e  (b) con la dipendenza di H1, sigma=c(1,0.01).
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages({library(dplyr); library(MMLB); library(lme4); library(geepack); library(geesmv)
                  library(parallel); library(MASS); library(tibble)})
SRC <- path.expand("~/MMM/R")
for (f in c("analysis_model3.R","analysis_smallsample.R","do_analysis.R",
            "generation_model3.R","generate_data.R","make_all_data.R"))
  sys.source(file.path(SRC,f), envir=globalenv())

# generatore ibrido: media di H0 (time:trt = 0) con la dipendenza di H1 (sigma = c(1,0.01)).
txt <- readLines(file.path(SRC,"generation_model3.R"))
i <- grep("mmlb   = c\\(-2, 0.1, 0.2, 0.2\\)", txt); stopifnot(length(i)==1)
txt[i] <- sub("c(-2, 0.1, 0.2, 0.2)", "c(-2, 0.1, 0.2, 0)", txt[i], fixed=TRUE)
e <- new.env(parent=globalenv()); eval(parse(text=paste(txt, collapse="\n")), envir=e)
gen_hyb <- get("generate_data_model3", envir=e)
cat("ibrido: riga patchata ->", trimws(txt[i]), "\n")

TECHS <- c("mmm","gee","glmer","gee_ss","glmer_lrt")
WALD  <- c("mmm","gee","gee_fg","gee_md","gee_kc","glmer","shrink")
one_rep <- function(N, T, i, mode) {
  set.seed(20260805L + 1000L*3L + i)
  d <- tryCatch(if (mode=="hybrid") gen_hyb(N,T,i,"mmlb") else make_all_data(N,T,i,"3","manual"),
                error=function(e) NULL)
  if (is.null(d)) return(NULL)
  rows <- lapply(TECHS, function(tc)
    tryCatch(do_analysis(d, N, T, "manual", 3, tc, i), error=function(e) NULL))
  r <- bind_rows(Filter(Negate(is.null), rows))
  if (!nrow(r)) return(NULL)
  r <- r %>% filter(term=="time:trt")
  if (!nrow(r)) return(NULL)
  r$converged <- ifelse(is.na(r$converged), FALSE, as.logical(r$converged))
  r$usable <- ifelse(r$method=="glmer_lrt", r$converged & is.finite(r$p),
    r$converged & is.finite(r$estimate) & is.finite(r$se) & is.finite(r$p) &
    abs(r$estimate)<10 & r$se>0 & r$se<10)
  r %>% transmute(method, usable, rej = is.finite(p) & p<0.05, N=N, ntime=T, rep=i, mode=mode)
}
NREP <- 300L
grid <- expand.grid(N=c(10,20,45), ntime=c(2,4,10), mode=c("asrun","hybrid"), stringsAsFactors=FALSE)
t0 <- Sys.time()
out <- list()
for (k in seq_len(nrow(grid))) {
  g <- grid[k,]
  if (as.numeric(difftime(Sys.time(), t0, units="mins")) > 150) { cat("CAP DI TEMPO raggiunto a k=",k,"\n"); break }
  rr <- mclapply(1:NREP, function(i) tryCatch(one_rep(g$N, g$ntime, i, g$mode), error=function(e) NULL),
                 mc.cores=20)
  out[[k]] <- bind_rows(Filter(Negate(is.null), rr))
  cat(sprintf("[%s] N=%d T=%d %s: righe %d (%.1f min)\n", format(Sys.time(),"%H:%M:%S"),
              g$N, g$ntime, g$mode, nrow(out[[k]]),
              as.numeric(difftime(Sys.time(), t0, units="mins")))); flush.console()
}
D <- bind_rows(out)
saveRDS(D, path.expand("~/MMM/rev2/R2_arm_raw.rds"))
S <- D %>% filter(usable) %>% group_by(method, N, ntime, mode) %>%
  summarise(n=n(), k=sum(rej), .groups="drop") %>%
  mutate(rate=100*k/n, mcse=100*sqrt((k/n)*(1-k/n)/n))
W <- S %>% select(method,N,ntime,mode,n,k,rate,mcse) %>%
  tidyr::pivot_wider(names_from=mode, values_from=c(n,k,rate,mcse)) %>%
  mutate(diff_pp = rate_hybrid - rate_asrun,
         mcse_diff = sqrt(mcse_asrun^2 + mcse_hybrid^2),
         z = diff_pp/mcse_diff)
cat("\n===== LIVELLO GenBinaryY M3, H0 come eseguito vs H0 con la dipendenza di H1 =====\n")
print(as.data.frame(W %>% mutate(across(where(is.numeric), ~round(.x,2))) %>% arrange(method,ntime,N)),
      row.names=FALSE)
cat("\n===== pooled per metodo sulle celle disponibili =====\n")
P <- S %>% group_by(method, mode) %>% summarise(n=sum(n), k=sum(k), .groups="drop") %>%
  mutate(rate=100*k/n, mcse=100*sqrt((k/n)*(1-k/n)/n)) %>%
  select(method,mode,n,k,rate,mcse) %>%
  tidyr::pivot_wider(names_from=mode, values_from=c(n,k,rate,mcse)) %>%
  mutate(diff_pp=round(rate_hybrid-rate_asrun,2),
         mcse_diff=round(sqrt(mcse_asrun^2+mcse_hybrid^2),2),
         z=round((rate_hybrid-rate_asrun)/sqrt(mcse_asrun^2+mcse_hybrid^2),2))
print(as.data.frame(P %>% mutate(across(where(is.numeric), ~round(.x,2)))), row.names=FALSE)
cat(sprintf("\n|diff| massima %.2f pp ; MCSE tipica della differenza %.2f pp ; celle con |z|>1.96: %d su %d\n",
  max(abs(W$diff_pp), na.rm=TRUE), median(W$mcse_diff, na.rm=TRUE),
  sum(abs(W$z)>1.96, na.rm=TRUE), sum(is.finite(W$z))))
write.csv(W, path.expand("~/MMM/rev2/R2_arm_level_bycell.csv"), row.names=FALSE)
write.csv(P, path.expand("~/MMM/rev2/R2_arm_level_pooled.csv"), row.names=FALSE)
cat("\nDONE_ARM\n")
