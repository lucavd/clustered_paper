
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({library(data.table)})
OUT <- path.expand("~/MMM/rev2")
CUT <- fread(file.path(OUT,"R2_311_variantA_cut.csv"))
Z <- CUT[scale == "z"]
cat("cells with q95 exactly 0, by method:\n")
print(Z[, list(ncell = .N, q95_eq0 = sum(q95 == 0), q95_tiny = sum(q95 < 1e-8)), by = method][order(method)])
cat("\nthe cells with a degenerate critical value:\n")
print(Z[q95 < 1e-8])
cat("\nsmallest ten q95 for shrink and glmer:\n")
print(head(Z[method %in% c("shrink","glmer")][order(q95)], 10))
