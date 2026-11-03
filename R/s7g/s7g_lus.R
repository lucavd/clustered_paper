# =============================================================================
# s7g_lus.R -- S7g, doubt 8 (neonatal cohort), prediction PG5.
# The mixed model of the cohort subsamples (rev2/T3_rerun.R: glmer, full mean model, (1 | id), bobyqa, maxfun 2e5) is refitted
# on the same 2,500 draws (rev2/T3_draws.rds, seed 20260807) with (1 | id) and with (1 | id) + (1 | pt), pt = patient-timepoint.
# Check: the (1 | id) refit reproduces n, mean SE and SD of tables/V5_S7f_seratio_fpc.csv (GLMM, full model) to 1e-6.
# Ratio = mean SE / SD of usable estimates x sqrt(1 - N/P); usable = finite estimate and SE, |estimate| <= 10, SE < 10.
# Output: tables/V5_S7g_lus_glmm.csv. Read-only on inputs.
# Explanation run (after PG5 was refuted; S7G_SLOPE=1): RE adds (1 + post | id), post = posterior zone, on draws with N >= 25 -> tables/V5_S7g_lus_glmm_slope.csv.
# =============================================================================
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages({library(dplyr); library(lme4); library(data.table); library(parallel)}); setDTthreads(1L)
setwd("/home/user/MMM"); TERM <- "regionposterior"; NCORE <- as.integer(Sys.getenv("S7G_CORES", "12"))
d0 <- readRDS("ARTICOLO 2026/LURID_bis/lurid_long.rds")
d0$row_in_id <- ave(seq_len(nrow(d0)), d0$id, FUN = seq_along); d0$tp <- (d0$row_in_id - 1L) %/% 10L + 1L
d <- d0 %>% filter(!is.na(lus_severe)) %>%
  mutate(region = factor(case_when(grepl("A$", as.character(section)) ~ "anterior", grepl("P$", as.character(section)) ~ "posterior", TRUE ~ "lateral"),
                         levels = c("anterior", "lateral", "posterior")),
         pt = factor(paste(id, tp, sep = "_")), id = droplevels(factor(id)))
f <- lus_severe ~ region + disease + invasive + age_weeks
fitg <- function(dat, re) {
  dat <- dat %>% arrange(id) %>% mutate(id = droplevels(factor(id)), across(where(is.factor), droplevels))
  q <- try(suppressWarnings(glmer(update(f, as.formula(paste(". ~ . +", re))), family = binomial, data = dat,
                                  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))), silent = TRUE)
  if (inherits(q, "try-error")) return(c(NA_real_, NA_real_, NA_real_))
  s <- try(summary(q)$coefficients, silent = TRUE)
  if (inherits(s, "try-error") || !TERM %in% rownames(s)) return(c(NA_real_, NA_real_, NA_real_))
  vc <- as.data.frame(VarCorr(q)); c(s[TERM, 1], s[TERM, 2], if ("pt" %in% vc$grp) vc$vcov[vc$grp == "pt"] else if (any(vc$var1 %in% "post" & is.na(vc$var2))) vc$vcov[vc$var1 %in% "post" & is.na(vc$var2)] else NA_real_)
}
SLOPE <- Sys.getenv("S7G_SLOPE") == "1"
d$post <- as.numeric(d$region == "posterior")
RE <- if (SLOPE) c(id = "(1 | id)", id_slope = "(1 + post | id)") else c(id = "(1 | id)", id_pt = "(1 | id) + (1 | pt)")
full <- rbindlist(lapply(names(RE), function(k) { r <- fitg(d, RE[[k]]); data.table(re = k, est = r[1], se = r[2], var_pt = r[3]) }))
print(full)
draws <- readRDS("rev2/T3_draws.rds"); P <- nlevels(d$id)
if (SLOPE) draws <- draws[vapply(draws, function(z) z$N >= 25, TRUE)]
t0 <- Sys.time()
res <- mclapply(seq_along(draws), function(i) { z <- draws[[i]]; s <- d[d$id %in% z$sel, ]
  rbindlist(lapply(names(RE), function(k) { r <- fitg(s, RE[[k]]); data.table(N = z$N, rep = z$rep, re = k, est = r[1], se = r[2], var_pt = r[3]) })) },
  mc.cores = NCORE, mc.preschedule = TRUE)
stopifnot(all(vapply(res, is.data.table, TRUE)))
R <- rbindlist(res); cat(sprintf("fits: %d in %.1f min\n", nrow(R), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
R[, use := is.finite(est) & is.finite(se) & abs(est) <= 10 & se < 10]
S <- R[use == TRUE, .(n = .N, mean_se = mean(se), sd_est = sd(est), median_var_pt = median(var_pt, na.rm = TRUE)), by = .(re, N)]
S[, `:=`(P = P, ratio = mean_se / sd_est, ratio_fpc = mean_se / sd_est * sqrt(1 - N / P))][, mcse_ratio := ratio / sqrt(2 * (n - 1))]
OLD <- fread("tables/V5_S7f_seratio_fpc.csv")[application == "neonatal cohort" & spec == "pieno" & method == "GLMM", .(N, n_old = n, mean_se_old = mean_se, sd_old = sd_est)]
CHK <- merge(S[re == "id"], OLD, by = "N"); CHK[, ok := n == n_old & abs(mean_se - mean_se_old) < 1e-6 & abs(sd_est - sd_old) < 1e-6]
print(CHK[, .(N, n, n_old, mean_se, mean_se_old, sd_est, sd_old, ok)])
S[, reproduces_stored := re == "id" & N %in% CHK[ok == TRUE, N]]
S <- merge(S, full[, .(re, full_est = est, full_se = se, full_var_pt = var_pt)], by = "re")
fwrite(S, if (SLOPE) "tables/V5_S7g_lus_glmm_slope.csv" else "tables/V5_S7g_lus_glmm.csv"); print(S, digits = 3)
cat(sprintf("all reproduce: %s\n", all(CHK$ok)))
