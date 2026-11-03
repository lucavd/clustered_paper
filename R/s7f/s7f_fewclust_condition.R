# s7f_fewclust_condition.R -- S7f, condition for citing fewclust as implementation of the paper's rule (pre-registered 2026-10-25).
# Applies auto_cluster_fit() of the registered fewclust (0.2.0, commit in ~/R/s7d_lib/fewclust_commit.txt) to the package data
# fc_trial and to one simulated dataset with 25 subjects and 4 occasions; records the chosen strategy.
.libPaths(c(path.expand("~/R/s7d_lib"), path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"), .libPaths()))
suppressMessages({ library(fewclust); library(data.table) })
setwd(path.expand("~/MMM"))
data(fc_trial, package = "fewclust")
cat("fc_trial columns:", paste(names(fc_trial), collapse = ", "), "\n")
sim <- simulate_clustered_binary(n_groups = 25, n_occasions = 4, seed = 20261027)
cat("sim columns:", paste(names(sim), collapse = ", "), "\n")
pick <- function(nm, cand) cand[cand %in% nm][1]
run <- function(d, label) {
  nm <- names(d); y <- pick(nm, c("y", "outcome", "resp")); g <- pick(nm, c("id", "cluster", "group", "subject"))
  tt <- pick(nm, c("time", "occasion", "visit")); tr <- pick(nm, c("trt", "treatment", "arm", "x"))
  f <- as.formula(sprintf("%s ~ %s * %s", y, tr, tt))
  fit <- tryCatch(auto_cluster_fit(f, data = d, group = g, time = tt), error = function(e) e)
  chosen <- if (inherits(fit, "error")) paste("error:", conditionMessage(fit)) else {
    s <- fit$strategy; if (is.null(s)) s <- attr(fit, "strategy"); if (is.null(s)) s <- class(fit)[1]; paste(unlist(s), collapse = "/") }
  data.table(dataset = label, formula = deparse(f), group = g, time = tt, chosen = chosen)
}
R <- rbind(run(as.data.frame(fc_trial), "fc_trial"), run(as.data.frame(sim), "simulated N=25 T=4 seed 20261027"))
R[, paper_rule := "gee_ar1_fg_t"]; R[, condition_met := grepl("ar1", chosen, ignore.case = TRUE) & grepl("fg|fay", chosen, ignore.case = TRUE)]
R[, fewclust_commit := readLines(path.expand("~/R/s7d_lib/fewclust_commit.txt"))[1]]
fwrite(R, "tables/V5_S7f_fewclust_condition.csv"); print(R)
