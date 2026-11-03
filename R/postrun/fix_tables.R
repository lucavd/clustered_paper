# =============================================================================
# fix_tables.R — regenerate the three tables the final check found stale or absent.
#
# The bootstrap and by-scenario tables were computed BEFORE the per-method filter
# fix, so `gee` and `glmer_lrt` carry values from the wrong replicate set
# (verified: gee worst 21.527 vs 20.900 corrected; glmer_lrt 5.256 vs 6.761).
# Coverage is promised in Methods and never reported in Results at all.
#
# The filter, restated: Wald-based methods need a usable SE; the likelihood-ratio
# test does not use one, so filtering its rows on `se` discards valid results.
# =============================================================================
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(dplyr); library(tidyr)})
RES <- path.expand("~/MMM/results_full"); OUT <- path.expand("~/MMM/tables")
set.seed(20260808)
DGPS <- c(GenBinaryY="GenBinaryY", GLMMnative="GLMM-native", GLMMAR1="GLMM-AR1")
WALD <- c("mmm","gee","gee_fg","gee_md","gee_kc","glmer","shrink")
LRT  <- "glmer_lrt"
METH <- c(WALD, LRT)
MARG <- c("1"=0.200,"2"=0.200,"3"=0.200,"4"=0.050,"5"=0.200)   # per-scenario truth
COND <- c("1"=0.232,"2"=0.232,"3"=0.232,"4"=0.058,"5"=0.232)

keep <- function(d) bind_rows(
  d %>% filter(method %in% WALD, converged, is.finite(estimate), is.finite(se),
               abs(estimate) < 10, se > 0, se < 10, is.finite(p)),
  d %>% filter(method == LRT, converged, is.finite(p)))

cw <- list(); sc <- list(); cv <- list(); pr <- list()
for (pref in names(DGPS)) {
  lab <- DGPS[[pref]]
  fs <- list.files(RES, pattern=paste0("^",pref,"_"), full.names=TRUE)
  k <- sub("_c[0-9]+[.]rds$","",basename(fs)); f <- names(which(table(k)==10))
  d <- bind_rows(lapply(fs[k %in% f], readRDS)); d$model <- as.integer(d$model)
  d <- d %>% filter(term == "time:trt")
  H0 <- keep(d %>% filter(label == "manual"))
  H1 <- keep(d %>% filter(label == "mmlb"))

  cw[[lab]] <- H0 %>% group_by(method, N) %>%
    summarise(n=n(), k=sum(p<0.05), .groups="drop") %>% mutate(dgp=lab)
  sc[[lab]] <- H0 %>% group_by(model, method, N) %>%
    summarise(rate=100*mean(p<0.05), .groups="drop") %>%
    group_by(model, method) %>% summarise(worst=max(abs(rate-5)), .groups="drop") %>%
    mutate(dgp=lab)

  # coverage: Wald methods only (an LRT yields no interval), per-scenario truth
  cv[[lab]] <- H1 %>% filter(method %in% WALD) %>%
    mutate(truth = ifelse(method=="glmer", COND[as.character(model)],
                                           MARG[as.character(model)]),
           covers = abs(estimate - truth) <= 1.959964*se) %>%
    group_by(method, N) %>%
    summarise(n=n(), cov=100*mean(covers),
              mcse=100*sqrt(mean(covers)*(1-mean(covers))/n()), .groups="drop") %>%
    mutate(dgp=lab)

  # paired rejections, MMM vs Fay-Graubard, matched replicate by replicate
  pr[[lab]] <- H0 %>% filter(method %in% c("mmm","gee_fg")) %>%
    select(model, N, ntime, sim_id, method, p) %>%
    mutate(rej = p < 0.05) %>% select(-p) %>%
    pivot_wider(names_from=method, values_from=rej) %>%
    filter(!is.na(mmm), !is.na(gee_fg)) %>%
    summarise(n_matched=n(), b=sum(mmm & !gee_fg), c=sum(!mmm & gee_fg)) %>%
    mutate(dgp=lab, diff_pp=100*(b-c)/n_matched,
           chisq=(abs(b-c)-1)^2/(b+c), p=pchisq(chisq,1,lower.tail=FALSE))
  rm(d,H0,H1); gc(verbose=FALSE)
}
CW <- bind_rows(cw); SC <- bind_rows(sc); CV <- bind_rows(cv); PR <- bind_rows(pr)

# ---- bootstrap of the worst-case maximum, on the CORRECTED counts ----
B <- 10000L
bw <- function(s) {
  ns <- s$n; ph <- s$k/s$n
  m <- matrix(rbinom(B*length(ns), rep(ns, each=B), rep(ph, each=B)), nrow=B) /
       matrix(rep(ns, each=B), nrow=B)
  apply(abs(100*m - 5), 1, max)
}
rows <- list(); i <- 0L
for (lab in unique(CW$dgp)) for (m in METH) {
  s <- CW %>% filter(dgp==lab, method==m) %>% arrange(N)
  if (!nrow(s)) next
  bs <- bw(s); obs <- max(abs(100*s$k/s$n - 5))
  i <- i+1L
  rows[[i]] <- data.frame(dgp=lab, method=m, worst_obs=obs, boot_mean=mean(bs),
                          lo=quantile(bs,.025), hi=quantile(bs,.975),
                          upward_bias=mean(bs)-obs)
}
W <- bind_rows(rows); rownames(W) <- NULL
write.csv(W,  file.path(OUT,"F_worstcase_bootstrap.csv"), row.names=FALSE)
write.csv(SC, file.path(OUT,"F_worstcase_by_scenario.csv"), row.names=FALSE)
write.csv(CV, file.path(OUT,"F_coverage.csv"), row.names=FALSE)
write.csv(PR, file.path(OUT,"F_paired_mmm_fg.csv"), row.names=FALSE)

cat("===== BOOTSTRAP on corrected counts (normal reference) =====\n")
print(as.data.frame(W %>% mutate(across(where(is.numeric), ~round(.x,3))) %>%
  arrange(dgp, worst_obs)), row.names=FALSE)
cat("\n-- P(MMM worst < FG worst) --\n")
for (lab in unique(CW$dgp)) {
  a <- bw(CW %>% filter(dgp==lab, method=="mmm") %>% arrange(N))
  b2 <- bw(CW %>% filter(dgp==lab, method=="gee_fg") %>% arrange(N))
  cat(sprintf("  %-12s P = %.3f  (%s)\n", lab, mean(a<b2),
      ifelse(mean(a<b2)>0.975|mean(a<b2)<0.025,"resolvable","NOT resolvable")))
}
cat("\n-- noise floor: median upward bias of the maximum --\n")
cat(sprintf("  %.3f pp (range %.3f-%.3f)\n", median(W$upward_bias),
            min(W$upward_bias), max(W$upward_bias)))

cat("\n===== PAIRED MMM vs Fay-Graubard (direction matters) =====\n")
print(as.data.frame(PR %>% mutate(across(c(diff_pp,chisq), ~round(.x,3)), p=signif(p,3))),
      row.names=FALSE)
cat("(diff_pp > 0 means MMM rejects MORE often than FG under the null, i.e. WORSE)\n")

cat("\n===== COVERAGE, pooled (nominal 95%) =====\n")
cvp <- CV %>% group_by(dgp, method) %>%
  summarise(cov=sum(cov*n)/sum(n), n=sum(n),
            mcse=100*sqrt((sum(cov*n/100)/sum(n))*(1-sum(cov*n/100)/sum(n))/sum(n)),
            .groups="drop")
print(as.data.frame(cvp %>% mutate(across(c(cov,mcse), ~round(.x,3))) %>%
  arrange(dgp, desc(cov))), row.names=FALSE)
cat("\n-- coverage at the small-N corner --\n")
print(as.data.frame(CV %>% filter(N %in% c(10,15,25)) %>%
  select(dgp, method, N, cov) %>% mutate(cov=round(cov,2)) %>%
  pivot_wider(names_from=N, values_from=cov, names_prefix="N=") %>%
  arrange(dgp, method)), row.names=FALSE)

cat("\n===== per-scenario worst case, corrected =====\n")
cat("times each method is best across the 15 (mechanism x scenario) cells:\n")
print(as.data.frame(SC %>% group_by(dgp, model) %>% slice_min(worst, n=1, with_ties=FALSE) %>%
  ungroup() %>% count(method, name="times_best") %>% arrange(desc(times_best))), row.names=FALSE)
cat("\nwrote F_worstcase_bootstrap.csv, F_worstcase_by_scenario.csv, F_coverage.csv, F_paired_mmm_fg.csv\n")
