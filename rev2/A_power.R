
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(dplyr); library(tidyr); library(parallel)})
RES <- path.expand("~/MMM/results_full"); OUT <- path.expand("~/MMM/rev2")
DGPS <- c(GenBinaryY="GenBinaryY", GLMMnative="GLMM-native", GLMMAR1="GLMM-AR1")
WALD <- c("mmm","gee","gee_fg","gee_md","gee_kc","glmer","shrink")
LRT  <- c("glmer_lrt")

keep_rows <- function(d) {
  w <- d %>% filter(method %in% WALD, converged, is.finite(estimate), is.finite(se),
                    abs(estimate) < 10, se > 0, se < 10, is.finite(p))
  l <- d %>% filter(method %in% LRT,  converged, is.finite(p))
  bind_rows(w, l)
}

cells <- list(); pw <- list(); sizechk <- list(); var_cell <- list()
for (pref in names(DGPS)) {
  lab <- DGPS[[pref]]
  fs <- list.files(RES, pattern=paste0("^",pref,"_"), full.names=TRUE)
  k <- sub("_c[0-9]+[.]rds$","",basename(fs)); f <- names(which(table(k)==10))
  fs <- fs[k %in% f]
  cat(lab, ": files used =", length(fs), " cells complete =", length(f), "\n")
  d <- bind_rows(mclapply(fs, function(x) {
         z <- readRDS(x); z <- z[z$term=="time:trt", c("estimate","se","p","method","N",
              "ntime","label","model","converged")]; z }, mc.cores=40))
  d$model <- as.integer(d$model)
  H0 <- keep_rows(d %>% filter(label == "manual"))
  H1 <- keep_rows(d %>% filter(label == "mmlb"))
  cat("  H0 usable rows =", nrow(H0), " H1 usable rows =", nrow(H1), "\n")

  ## ---------- VARIANT C : DEC-CRIT, quantile on p  (post_gan.R P2 verbatim) ----
  cutC <- H0 %>% group_by(model, N, ntime, method) %>%
    summarise(q05 = quantile(p, 0.05, na.rm = TRUE), nH0 = n(), .groups = "drop")
  PWc <- H1 %>% inner_join(cutC, by = c("model","N","ntime","method")) %>%
    group_by(method, N) %>%
    summarise(n = n(), power_raw = 100*mean(p < 0.05), power_adj = 100*mean(p < q05),
              mcse = 100*sqrt(mean(p<q05)*(1-mean(p<q05))/n()), .groups="drop") %>%
    mutate(dgp = lab)
  pw[[lab]] <- PWc

  ## ---------- SIZE VALIDATION : apply q05 back to the SAME H0 replicates -------
  sizechk[[lab]] <- H0 %>% inner_join(cutC, by=c("model","N","ntime","method")) %>%
    group_by(method) %>% summarise(n = n(), size_emp = 100*mean(p < q05),
      size_cellmean = 100*mean(tapply(p<q05, interaction(model,N,ntime), mean)[
                                 as.character(interaction(model,N,ntime))]),
      .groups="drop") %>% mutate(dgp = lab)

  ## ---------- VARIANTS A,B (|z| quantiles), D, E : per-cell power --------------
  H0z <- H0 %>% mutate(z = ifelse(method %in% LRT, NA_real_, abs(estimate/se)))
  H1z <- H1 %>% mutate(z = ifelse(method %in% LRT, NA_real_, abs(estimate/se)))
  cutAB <- H0z %>% group_by(model, N, ntime, method) %>%
    summarise(q95_t7 = if (all(is.na(z))) NA_real_ else quantile(z, 0.95, na.rm=TRUE, type=7),
              q95_t1 = if (all(is.na(z))) NA_real_ else quantile(z, 0.95, na.rm=TRUE, type=1),
              q05_p7 = quantile(p, 0.05, na.rm=TRUE, type=7),
              q05_p1 = quantile(p, 0.05, na.rm=TRUE, type=1),
              nH0 = n(), .groups="drop")
  vc <- H1z %>% inner_join(cutAB, by=c("model","N","ntime","method")) %>%
    group_by(dgp=lab, model, N, ntime, method) %>%
    summarise(nH1 = n(), nH0 = nH0[1],
              pA = 100*mean(if (all(is.na(z))) p < q05_p7 else z > q95_t7),
              pB = 100*mean(if (all(is.na(z))) p < q05_p1 else z > q95_t1),
              pC = 100*mean(p < q05_p7),
              praw = 100*mean(p < 0.05), .groups="drop")
  var_cell[[lab]] <- vc
  rm(d, H0, H1, H0z, H1z); gc(verbose=FALSE)
}
PW <- bind_rows(pw); SZ <- bind_rows(sizechk); VC <- bind_rows(var_cell)
write.csv(PW, file.path(OUT,"R2_power_sizeadj_current.csv"), row.names=FALSE)
write.csv(SZ, file.path(OUT,"R2_power_sizecheck.csv"), row.names=FALSE)
write.csv(VC, file.path(OUT,"R2_power_variants_cell.csv"), row.names=FALSE)

cat("\n================ TABLE 3.11 : pooled raw and size-adjusted power ==========\n")
tab <- PW %>% group_by(dgp, method) %>%
  summarise(nH1 = sum(n), raw = weighted.mean(power_raw, n),
            adj = weighted.mean(power_adj, n), .groups="drop")
print(as.data.frame(tab %>% mutate(raw=round(raw,3), adj=round(adj,3)) %>%
      arrange(dgp, desc(adj))), row.names=FALSE)
cat("\n---- spread of size-adjusted power per mechanism ----\n")
print(as.data.frame(tab %>% group_by(dgp) %>% summarise(min=round(min(adj),3),
      max=round(max(adj),3), argmin=method[which.min(adj)], argmax=method[which.max(adj)])),
      row.names=FALSE)
cat("\n---- GLMM vs Fay-Graubard gap (size-adjusted) ----\n")
print(as.data.frame(tab %>% filter(method %in% c("glmer","gee_fg")) %>%
      select(dgp, method, adj) %>% pivot_wider(names_from=method, values_from=adj) %>%
      mutate(gap = round(glmer-gee_fg,4), glmer=round(glmer,4), gee_fg=round(gee_fg,4))),
      row.names=FALSE)
cat("\n---- RAW power ranking, GenBinaryY (R2-08) ----\n")
print(as.data.frame(tab %>% filter(dgp=="GenBinaryY") %>% arrange(desc(raw)) %>%
      select(method, raw, nH1) %>% mutate(raw=round(raw,3))), row.names=FALSE)

cat("\n================ SIZE VALIDATION of the protocol (R2-09) ==================\n")
cat("empirical size when q05 is applied to the H0 replicates it came from:\n")
print(as.data.frame(SZ %>% mutate(size_emp=round(size_emp,4)) %>%
      select(dgp, method, n, size_emp) %>% arrange(dgp, method)), row.names=FALSE)
cat("range over the 24 (dgp,method):", round(min(SZ$size_emp),4), "-",
    round(max(SZ$size_emp),4), "\n")

cat("\n================ VARIANTS (R2-09) ========================================\n")
pool <- function(v, w) weighted.mean(v, w)
vv <- VC %>% group_by(dgp, method) %>% summarise(
  A = pool(pA, nH1), B = pool(pB, nH1), C = pool(pC, nH1),
  D = { i <- nH0 >= 20 & nH1 >= 20; if (any(i)) pool(pC[i], nH1[i]) else NA_real_ },
  E = mean(pC), ncell = n(),
  ncell_drop = sum(!(nH0 >= 20 & nH1 >= 20)), .groups="drop")
print(as.data.frame(vv %>% mutate(across(c(A,B,C,D,E), ~round(.x,3)))), row.names=FALSE)
cat("\n---- GLMM-vs-FG gap under each variant ----\n")
gg <- vv %>% filter(method %in% c("glmer","gee_fg")) %>%
  select(dgp, method, A,B,C,D,E) %>%
  pivot_longer(A:E, names_to="variant", values_to="pw") %>%
  pivot_wider(names_from=method, values_from=pw) %>%
  mutate(gap = round(glmer - gee_fg, 4))
print(as.data.frame(gg %>% mutate(glmer=round(glmer,3), gee_fg=round(gee_fg,3))), row.names=FALSE)
cat("\nA minus C gap difference (GenBinaryY):",
    round(gg$gap[gg$dgp=="GenBinaryY" & gg$variant=="A"] -
          gg$gap[gg$dgp=="GenBinaryY" & gg$variant=="C"], 4), "\n")
write.csv(vv, file.path(OUT,"R2_power_variants_pooled.csv"), row.names=FALSE)
write.csv(tab, file.path(OUT,"R2_power_311_pooled.csv"), row.names=FALSE)
cat("\nDONE_A\n")
