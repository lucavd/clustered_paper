
# s10s_lus_to.R -- S10 "supplemento": (1) independent refit of the cohort mixed model with a random slope for posterior
# zones (Section S17) on the stored subsample draws with N >= 25; (2) share of geesmv time-outs per cell in the study rows.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages({library(data.table); library(lme4); library(parallel)}); setDTthreads(1)
setwd("/home/user/MMM"); D <- "results_v4/S10/supplemento"
d <- as.data.frame(readRDS("ARTICOLO 2026/LURID_bis/lurid_long.rds")); cat("cols:", names(d), "\n")
d <- d[!is.na(d$lus_severe), ]
sec <- as.character(d$section)
d$region <- factor(ifelse(grepl("A$", sec), "anterior", ifelse(grepl("P$", sec), "posterior", "lateral")), levels = c("anterior","lateral","posterior"))
d$post <- as.numeric(d$region == "posterior"); d$id <- droplevels(factor(d$id)); P <- nlevels(d$id)
cat("patients", P, "rows", nrow(d), "\n")
fit <- function(s) { s$id <- droplevels(factor(s$id))
  q <- try(suppressWarnings(glmer(lus_severe ~ region + disease + invasive + age_weeks + (1 + post | id), family = binomial, data = s,
           control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))), silent = TRUE)
  if (inherits(q, "try-error")) return(c(NA, NA)); cf <- summary(q)$coefficients
  if (!"regionposterior" %in% rownames(cf)) return(c(NA, NA)); cf["regionposterior", 1:2] }
full <- fit(d); cat("full:", full, "\n")
dr <- readRDS("rev2/T3_draws.rds"); dr <- dr[vapply(dr, function(z) z$N >= 25, TRUE)]
R <- rbindlist(mclapply(dr, function(z) { r <- fit(d[d$id %in% z$sel, ]); data.table(N = z$N, rep = z$rep, est = r[1], se = r[2]) }, mc.cores = 2))
R[, use := is.finite(est) & is.finite(se) & abs(est) < 10 & se > 0 & se < 10]
S <- R[use == TRUE, .(n = .N, mean_se = mean(se), sd_est = sd(est)), by = N][, ratio_fpc := mean_se / sd_est * sqrt(1 - N / P)]
print(S); fwrite(S, file.path(D, "s10s_lus_slope.csv"))
# geesmv time-outs in the study rows
S2 <- readRDS(file.path(D, "study_cc.rds"))[method %in% c("gee_md", "gee_fg", "gee_kc")]
S2[, to := grepl("timeout|time-out|timed out", status, ignore.case = TRUE)]
print(S2[, .N, by = .(method, to)])
TO <- S2[, .(n = .N, pct_to = 100 * mean(to)), by = .(dgp, scen, N, ntime, label, method)]
print(TO[order(-pct_to)][1:8]); fwrite(TO, file.path(D, "s10s_timeouts.csv"))
cat("LUS_TO OK\n")
