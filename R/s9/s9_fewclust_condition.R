# s9_fewclust_condition.R -- S9, condition for citing fewclust as implementation of the paper's rule (D-43),
# re-run for the rule of the corrected study (D-82: marginalized model with the convergence check + t).
# Same data and same call as R/s7f/s7f_fewclust_condition.R (fc_trial; simulated N = 25, T = 4, seed 20261027);
# fewclust 0.3.1 installed from commit c688ecc of branch package into ~/R/s9_fewclust_lib (by s9_fewclust_condition.sh).
lib <- path.expand("~/R/s9_fewclust_lib")
.libPaths(c(lib, path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"), .libPaths()))
suppressMessages({ library(fewclust, lib.loc = lib); library(data.table) })
setwd(path.expand("~/MMM"))
cat("fewclust", as.character(packageVersion("fewclust", lib.loc = lib)), "\n")
data(fc_trial, package = "fewclust")
sim <- simulate_clustered_binary(n_groups = 25, n_occasions = 4, seed = 20261027)
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
R[, paper_rule := "mmm_cc_t"]; R[, condition_met := tolower(chosen) == "mmm"]
R[, fewclust_version := as.character(packageVersion("fewclust", lib.loc = lib))]
R[, fewclust_commit := readLines(file.path(lib, "fewclust_commit.txt"))[1]]
fwrite(R, "tables/V6_S9_fewclust_condition.csv"); print(R)
