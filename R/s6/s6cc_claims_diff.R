# s6cc_claims_diff.R -- D-32: comparative sub-claims of the manuscript whose classification changes when MMM in
# M3-M5 is analysed with the convergence check. Inputs: tables/V4_S5_claims.csv (study),
# tables/variant_cc/V4_S5_claims.csv (variant), reports/morris/S5_claims_inventory.csv.
# Output: tables/V4_S6cc_claims_diff.csv
source(path.expand("~/MMM/R/s4/s4_common.R"))
A <- fread(file.path(REPO, "tables/V4_S5_claims.csv")); B <- fread(file.path(REPO, "tables/variant_cc/V4_S5_claims.csv"))
M <- merge(A[, .(name, class_study = class, est_study = est)], B[, .(name, class_cc = class, est_cc = est)], by = "name")
M[, changed := class_study != class_cc]
fwrite(M[changed == TRUE][order(name)], file.path(REPO, "tables/V4_S6cc_claims_diff.csv"))
cat("sub-claims", nrow(M), "changed", M[, sum(changed)], "\n"); print(M[changed == TRUE, .N, by = .(class_study, class_cc)][order(-N)])
