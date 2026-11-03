# s10t_clin.R -- S10 "Testo": independent refit of the full-data numbers of Sections 2.7 and 3.11 from the raw data
# (neonatal cohort, ARTICOLO 2026/LURID_bis/lurid_long.rds; respiratory trial, geepack::respiratory). Prints aggregates only.
.libPaths(c("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu", "/home/user/R/s7d_lib"))
suppressMessages({library(data.table); library(lme4); library(geepack); library(MMLB); library(geesmv)}); setDTthreads(1L)
setwd("/home/user/MMM"); OUT <- "results_v4/S10/testo"
res <- list(); put <- function(k, v, how) { res[[length(res) + 1]] <<- data.table(item = k, value = as.numeric(v), how = how); cat(k, v, "\n") }
d0 <- readRDS("ARTICOLO 2026/LURID_bis/lurid_long.rds"); d0 <- as.data.frame(d0)
d0$row_in_id <- ave(seq_len(nrow(d0)), d0$id, FUN = seq_along); d0$tp <- (d0$row_in_id - 1L) %/% 10L + 1L
d <- d0[!is.na(d0$lus_severe), ]
sec <- as.character(d$section)
d$region <- factor(ifelse(grepl("A$", sec), "anterior", ifelse(grepl("P$", sec), "posterior", "lateral")), levels = c("anterior", "lateral", "posterior"))
d$id <- droplevels(factor(d$id)); d$pt <- droplevels(factor(paste(d$id, d$tp, sep = "_")))
put("lus_npat", nlevels(d$id), "patients with an outcome"); put("lus_npt", nlevels(d$pt), "patient-timepoints"); put("lus_nobs", nrow(d), "zone-level outcomes")
k <- tapply(as.character(d$pt), d$id, function(x) length(unique(x))); put("lus_nmulti", sum(k > 1), "patients with more than one timepoint")
put("lus_prev", mean(d$lus_severe), "prevalence of severe loss of aeration"); put("lus_rows_per_pt_10", mean(table(d$pt) == 10), "share of patient-timepoints with 10 zones")
ctl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
m1 <- glmer(lus_severe ~ region + (1 | id), family = binomial, data = d, control = ctl)
m2 <- glmer(lus_severe ~ region + (1 | id) + (1 | pt), family = binomial, data = d, control = ctl)
vc <- as.data.frame(VarCorr(m2))
put("lus_var_id", as.data.frame(VarCorr(m1))$vcov[1], "patient variance, region + (1|id)")
put("lus_var_pt", vc$vcov[vc$grp == "pt"], "patient-timepoint variance, + (1|pt)"); put("lus_var_id_in_m2", vc$vcov[vc$grp == "id"], "patient variance in the two-level model")
put("lus_chi2", 2 * (as.numeric(logLik(m2)) - as.numeric(logLik(m1))), "likelihood-ratio chi2 for the patient-timepoint level")
D <- d[order(d$id), ]; D$id <- droplevels(D$id)
f <- lus_severe ~ region + disease + invasive + age_weeks
g <- geeglm(f, family = binomial, id = id, corstr = "exchangeable", data = D)
put("lus_est_gee", coef(g)["regionposterior"], "GEE exch"); put("lus_se_gee", sqrt(diag(vcov(g)))["regionposterior"], "GEE exch uncorrected sandwich SE")
q <- glmer(update(f, . ~ . + (1 | id)), family = binomial, data = D, control = ctl)
put("lus_est_glmm", fixef(q)["regionposterior"], "GLMM (1|id)"); put("lus_se_glmm", sqrt(diag(as.matrix(vcov(q))))["regionposterior"], "GLMM SE")
m <- mm(f, lv.formula = ~1, t.formula = ~1, data = D, id = id, q = 20, verbose = FALSE)
put("lus_est_mmm", m$beta["regionposterior"], "MMM lv ~1, t ~1, q = 20"); put("lus_se_mmm", sqrt(m$mod.cov["regionposterior", "regionposterior"]), "MMM model-based SE")
put("lus_se_mmm_rob", sqrt(m$rob.cov["regionposterior", "regionposterior"]), "MMM robust SE")
put("lus_lateral_mmm", m$beta["regionlateral"], "MMM lateral vs anterior")
m0 <- mm(f, t.formula = ~1, data = D, id = id, verbose = FALSE)
put("lus_est_mmm_nolv", m0$beta["regionposterior"], "MMM without latent variable (rev2 specification)")
inv <- tapply(D$invasive == 1 | D$invasive == "1", D$id, any); put("lus_n_inv_pat", sum(inv), "patients with an invasive row")
put("lus_p_noinv_N10", dhyper(0, sum(inv), length(inv) - sum(inv), 10), "hypergeometric P(no invasive patient) in 10 patients")
# ---- respiratory trial ------------------------------------------------------------------------------------------
data(respiratory, package = "geepack"); R <- as.data.frame(respiratory)
R$pid <- factor(paste(R$center, R$id, sep = "_")); R$id <- as.integer(R$pid)
R$treat <- factor(as.character(R$treat), levels = c("P", "A")); R$sex <- factor(as.character(R$sex)); R$center <- factor(R$center)
R <- R[order(R$pid, R$visit), ]
put("resp_P", nlevels(R$pid), "patients"); put("resp_visits", length(unique(R$visit)), "visits")
fr <- outcome ~ treat + baseline + visit + age + sex + center
gr <- geeglm(fr, family = binomial, id = id, corstr = "exchangeable", data = R)
put("resp_est_gee", coef(gr)["treatA"], "GEE exch"); put("resp_se_gee", sqrt(diag(vcov(gr)))["treatA"], "uncorrected SE")
kk <- match("treatA", names(coef(gr)))
for (ty in c("md", "kc", "fg")) { v <- do.call(getExportedValue("geesmv", paste0("GEE.var.", ty)), list(fr, id = "id", family = binomial, data = R, corstr = "exchangeable"))
  put(paste0("resp_se_", ty, "_geesmv"), sqrt(as.numeric(v$cov.beta)[kk]), paste("geesmv", ty)) }
if (requireNamespace("fewclust", quietly = TRUE)) for (ty in c("md", "kc", "fg")) put(paste0("resp_se_", ty, "_fewclust"), sqrt(fewclust::gee_vcov(gr, type = ty)[kk, kk]), paste("fewclust", ty))
qr <- glmer(update(fr, . ~ . + (1 | id)), family = binomial, data = R, control = ctl)
put("resp_est_glmm", fixef(qr)["treatA"], "GLMM (1|id)")
mr <- mm(fr, lv.formula = ~1, t.formula = ~1, data = R, id = id, q = 20, verbose = FALSE)
put("resp_est_mmm", mr$beta["treatA"], "MMM"); put("resp_se_mmm", sqrt(mr$mod.cov["treatA", "treatA"]), "MMM model-based SE"); put("resp_se_mmm_rob", sqrt(mr$rob.cov["treatA", "treatA"]), "MMM robust SE")
put("resp_ncoef", length(coef(gr)), "mean-model coefficients")
fwrite(rbindlist(res), file.path(OUT, "clin.csv"))
