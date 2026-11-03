# =============================================================================
# post_gan.R — the analyses the adversarial verdicts made necessary.
#
# P1 REBUILD every headline number with the CORRECT per-method filter. The LRT
#    p-value does not use the Wald SE, so filtering LRT rows on `se` discarded
#    ~44,000 valid results and understated its type I error (9.32 -> 10.92%).
#    Wald-based methods still need a usable SE; LRT-based ones do not.
# P2 SIZE-ADJUSTED POWER. A method with nominal level but no power is useless at
#    few clusters, and we never reported power. Because raw power is confounded
#    with size, compare at a COMMON empirical size: take each method's own
#    empirical 5% critical value from the H0 replicates in that cell, then apply
#    it to the H1 replicates of the same cell.
# P3 UNCONDITIONAL VALIDITY: P(usable result) x P(correct decision | usable).
#    This is the endpoint that separates the methods most sharply, and it treats
#    "returns nothing" as the failure it is rather than dropping it.
# P4 MEAN ABSOLUTE DEVIATION across N, alongside the worst case — the reviewer
#    showed every non-maximum summary favours Fay-Graubard, so both must appear.
# =============================================================================
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB); suppressMessages({library(dplyr); library(tidyr)})
RES <- path.expand("~/MMM/results_full"); OUT <- path.expand("~/MMM/tables")
DGPS <- c(GenBinaryY="GenBinaryY", GLMMnative="GLMM-native", GLMMAR1="GLMM-AR1")
WALD <- c("mmm","gee","gee_fg","gee_md","gee_kc","glmer","shrink")   # need a usable SE
LRT  <- c("glmer_lrt")                                               # does not
MARG <- c("1"=0.200,"2"=0.200,"3"=0.200,"4"=0.050,"5"=0.200)
COND <- c("1"=0.232,"2"=0.232,"3"=0.232,"4"=0.058,"5"=0.232)

keep_rows <- function(d) {
  w <- d %>% filter(method %in% WALD, converged, is.finite(estimate), is.finite(se),
                    abs(estimate) < 10, se > 0, se < 10, is.finite(p))
  l <- d %>% filter(method %in% LRT,  converged, is.finite(p))
  bind_rows(w, l)
}

t1 <- list(); pw <- list(); uv <- list()
for (pref in names(DGPS)) {
  lab <- DGPS[[pref]]
  fs <- list.files(RES, pattern=paste0("^",pref,"_"), full.names=TRUE)
  k <- sub("_c[0-9]+[.]rds$","",basename(fs)); f <- names(which(table(k)==10))
  d <- bind_rows(lapply(fs[k %in% f], readRDS)); d$model <- as.integer(d$model)
  d <- d %>% filter(term == "time:trt")
  attempts <- d %>% group_by(model, N, ntime, method) %>%
    summarise(n_attempt = 1000L, .groups="drop")     # 10 chunks x 100 replicates

  H0 <- keep_rows(d %>% filter(label == "manual"))
  H1 <- keep_rows(d %>% filter(label == "mmlb"))

  # ---- P1 / P4: level, per N, with the correct filter ----
  t1[[lab]] <- H0 %>% group_by(method, N) %>%
    summarise(n = n(), rate = 100*mean(p < 0.05),
              mcse = 100*sqrt(mean(p<0.05)*(1-mean(p<0.05))/n()), .groups="drop") %>%
    mutate(dgp = lab)

  # ---- P2: size-adjusted power, per cell ----
  # critical |z| (or LRT p) achieving exactly 5% rejection on THIS cell's H0 draws
  cut <- H0 %>% group_by(model, N, ntime, method) %>%
    summarise(q05 = quantile(p, 0.05, na.rm = TRUE), .groups = "drop")
  pw[[lab]] <- H1 %>% inner_join(cut, by = c("model","N","ntime","method")) %>%
    group_by(method, N) %>%
    summarise(n = n(),
              power_raw = 100*mean(p < 0.05),
              power_adj = 100*mean(p < q05),
              mcse = 100*sqrt(mean(p<q05)*(1-mean(p<q05))/n()), .groups="drop") %>%
    mutate(dgp = lab)

  # ---- P3: unconditional validity under H0 ----
  # P(usable) from attempts; P(correct decision) = P(not rejecting | usable) at 5%
  usable <- H0 %>% group_by(model, N, ntime, method) %>%
    summarise(n_usable = n(), rej = mean(p < 0.05), .groups="drop") %>%
    left_join(attempts, by = c("model","N","ntime","method")) %>%
    mutate(p_usable = n_usable / n_attempt,
           uncond_valid = p_usable * (1 - rej))
  uv[[lab]] <- usable %>% group_by(method, N) %>%
    summarise(p_usable = 100*mean(p_usable),
              uncond = 100*mean(uncond_valid), .groups="drop") %>% mutate(dgp = lab)
  rm(d, H0, H1); gc(verbose=FALSE)
}
T1 <- bind_rows(t1); PW <- bind_rows(pw); UV <- bind_rows(uv)
write.csv(T1, file.path(OUT,"P_level_corrected.csv"), row.names=FALSE)
write.csv(PW, file.path(OUT,"P_power_sizeadj.csv"),   row.names=FALSE)
write.csv(UV, file.path(OUT,"P_unconditional.csv"),   row.names=FALSE)

cat("=====================================================================\n")
cat("P1 — TYPE I ERROR with the correct per-method filter (nominal 5%)\n")
cat("=====================================================================\n")
lev <- T1 %>% group_by(dgp, method) %>%
  summarise(pooled = 100*sum(rate*n/100)/sum(n),
            worst  = max(abs(rate - 5)),
            mad    = mean(abs(rate - 5)), .groups="drop")
print(as.data.frame(lev %>% mutate(across(where(is.numeric), ~round(.x,2))) %>%
  arrange(dgp, worst)), row.names=FALSE)
cat("\n(mad = mean absolute deviation across the 16 cluster counts; the reviewer\n")
cat(" showed every non-maximum summary favours Fay-Graubard, so both appear.)\n")

cat("\n=====================================================================\n")
cat("P2 — SIZE-ADJUSTED POWER (each method judged at its own empirical 5%)\n")
cat("=====================================================================\n")
p2 <- PW %>% filter(N %in% c(10,15,25,45,85,145)) %>%
  select(dgp, method, N, power_adj) %>% mutate(power_adj = round(power_adj,1)) %>%
  pivot_wider(names_from = N, values_from = power_adj, names_prefix = "N=")
print(as.data.frame(p2 %>% arrange(dgp, method)), row.names=FALSE)
cat("\nraw vs size-adjusted, pooled (raw power is confounded with size):\n")
print(as.data.frame(PW %>% group_by(dgp, method) %>%
  summarise(raw = round(weighted.mean(power_raw, n),1),
            adj = round(weighted.mean(power_adj, n),1), .groups="drop") %>%
  arrange(dgp, desc(adj))), row.names=FALSE)

cat("\n=====================================================================\n")
cat("P3 — UNCONDITIONAL VALIDITY: P(usable) x P(no false rejection | usable)\n")
cat("=====================================================================\n")
p3 <- UV %>% filter(N %in% c(10,15,25,45,145)) %>%
  select(dgp, method, N, uncond) %>% mutate(uncond = round(uncond,1)) %>%
  pivot_wider(names_from = N, values_from = uncond, names_prefix = "N=")
print(as.data.frame(p3 %>% arrange(dgp, method)), row.names=FALSE)
cat("\nP(usable result) alone, at the small-N corner:\n")
print(as.data.frame(UV %>% filter(N %in% c(10,15,25)) %>%
  select(dgp, method, N, p_usable) %>% mutate(p_usable = round(p_usable,1)) %>%
  pivot_wider(names_from = N, values_from = p_usable, names_prefix = "N=") %>%
  arrange(dgp, method)), row.names=FALSE)
cat("\nwrote P_level_corrected.csv, P_power_sizeadj.csv, P_unconditional.csv\n")
