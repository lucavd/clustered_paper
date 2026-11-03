# s11_thr_compare.R -- S11 (finding F13, prediction PC7): the sub-claims of the claim table (V4_S5_claims.csv) and the
# rule safety of the marginalized model under the thresholds 1e-4 (variant t4cv) and 1e-2 (variant t2cv) of the convergence
# check, against the primary threshold 1e-3 (variant cv), all with the rules D-82 and D-90. Also: which sub-claims are
# cited by a claim tag (cl:) in the manuscript and the supplement. Outputs: tables/V8_S11_thr_claims.csv (one row per
# sub-claim and variant), tables/V8_S11_thr_summary.csv (key, value) read by R/s11/s11_keys.R.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setwd("/home/user/MMM")
C0 <- fread("tables/variant_cv/V4_S5_claims.csv")
src <- paste(readLines("manuscript/v5/MMM_v5.qmd"), readLines("manuscript/v5/MMM_v5_supplement.qmd"), collapse = "\n")
cited <- unique(sub(":bobyqa$", "", regmatches(src, gregexpr("(?<=cl:)[^ ]+?(?= as=)", src, perl = TRUE))[[1]]))
OUT <- list(); SUM <- list()
for (tg in c("t4cv", "t2cv")) {
  C1 <- fread(sprintf("tables/variant_%s/V4_S5_claims.csv", tg))
  M <- merge(C0[, .(name, class_cv = class, est_cv = est, lo_cv = lo, hi_cv = hi)], C1[, .(name, class_v = class, est_v = est, lo_v = lo, hi_v = hi)], by = "name")
  stopifnot(nrow(M) == nrow(C0), nrow(M) == nrow(C1))
  M[, `:=`(variant = tg, changed = class_cv != class_v, cited = name %in% cited)]
  OUT[[tg]] <- M
  RS <- fread(sprintf("tables/variant_%s/V4_S6_rule_safety.csv", tg))[is.na(N) & band == "all" & strategy == "mmm|t"]
  SUM[[tg]] <- setnames(data.table(k = paste0("thr|", tg, c("|n", "|changed", "|cited", "|cited_changed", "|mmm_t_hi")),
                          value = c(nrow(M), M[, sum(changed)], M[, sum(cited)], M[, sum(changed & cited)], sum(RS$n_hi))), "k", "key")   # data.table(key = ...) would set a key
}
X <- rbindlist(OUT); fwrite(X, "tables/V8_S11_thr_claims.csv"); S <- rbindlist(SUM); fwrite(S, "tables/V8_S11_thr_summary.csv")
print(S); print(X[changed == TRUE, .(variant, name, class_cv, class_v, cited)])
