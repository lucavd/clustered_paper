
# s10s_resp.R -- S10 "supplemento": independent refit of the respiratory trial (Table S29, Section S9) from
# geepack::respiratory, with the models described in the rendered text. No project script is sourced.
.libPaths(c("/home/user/R/s7d_lib", "/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(geepack); library(lme4); library(MMLB); library(geesmv)})
data(respiratory, package = "geepack"); d <- respiratory
d$id <- as.integer(factor(paste(d$center, d$id)))          # patients are numbered within centre
d$trtA <- as.integer(d$treat == "A"); d$center <- factor(d$center)
d <- d[order(d$id, d$visit), ]
cat("patients", length(unique(d$id)), "rows", nrow(d), "visits", paste(sort(unique(d$visit)), collapse = ","), "\n")
cat("treat levels", levels(respiratory$treat), "\n")
ff <- outcome ~ trtA + baseline + visit + age + sex + center; fr <- outcome ~ trtA + baseline
out <- list()
for (sp in c("full", "reduced")) {
  f <- if (sp == "full") ff else fr
  g <- geeglm(f, id = id, data = d, family = binomial, corstr = "exchangeable")
  est <- coef(g)[["trtA"]]; se_unc <- sqrt(vcov(g)["trtA", "trtA"])
  dd <- d; dd$id <- d$id
  md <- GEE.var.md(f, id = "id", family = binomial, data = dd, corstr = "exchangeable")
  fg <- GEE.var.fg(f, id = "id", family = binomial, data = dd, corstr = "exchangeable")
  j <- which(names(coef(g)) == "trtA")
  kc <- tryCatch({ if (requireNamespace("fewclust", quietly = TRUE)) { fc <- fewclust::vcov_gee(g, type = "KC"); sqrt(fc[j, j]) } else NA }, error = function(e) { cat("fewclust:", conditionMessage(e), "\n"); NA })
  out[[paste(sp, "gee")]] <- c(est = est, se_unc = se_unc, se_md = sqrt(md$cov.beta[j]), se_fg = sqrt(fg$cov.beta[j]), se_kc = kc)
  gl <- glmer(update(f, . ~ . + (1 | id)), data = d, family = binomial)
  out[[paste(sp, "glmm")]] <- c(est = fixef(gl)[["trtA"]], se = sqrt(vcov(gl)["trtA", "trtA"]))
  for (tf in c("with_transition", "no_transition")) {
    m <- tryCatch(if (tf == "with_transition") mm(f, lv.formula = ~1, t.formula = ~1, id = id, data = d, q = 20)
                  else mm(f, lv.formula = ~1, id = id, data = d, q = 20), error = function(e) e)
    if (inherits(m, "error")) { cat(sp, tf, "error", conditionMessage(m), "\n"); next }
    b <- m$beta; k <- which(names(b) == "trtA")
    vm <- if (!is.null(m$mod.cov)) m$mod.cov else NA; vr <- if (!is.null(m$rob.cov)) m$rob.cov else NA
    out[[paste(sp, "mmm", tf)]] <- c(est = b[[k]], se_mod = if (is.matrix(vm)) sqrt(vm[k, k]) else NA, se_rob = if (is.matrix(vr)) sqrt(vr[k, k]) else NA)
    if (sp == "full" && tf == "with_transition") { cat("mm object names:", names(m), "\n") }
  }
}
print(out)
cat("fewclust exports:", if (requireNamespace("fewclust", quietly = TRUE)) paste(getNamespaceExports("fewclust"), collapse = " ") else "not installed", "\n")
