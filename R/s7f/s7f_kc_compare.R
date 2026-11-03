# s7f_kc_compare.R -- S7f, predictions PF3 and PF4: what changes when Kauermann-Carroll comes from fewclust (variant kcfc)
# instead of geesmv (variant_cc). Read-only on both variant trees.
# Outputs: tables/V5_S7f_claims_diff.csv (every sub-claim of V4_S5_claims.csv and V4_S6_m1c_claims.csv: values and class in
#          both variants), tables/V5_S7f_values_diff.csv (every value of V4_S6_values.csv that changed), tables/V5_S7f_PF.csv.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(1L); setwd("/home/user/MMM")
cmp <- function(a, b, what) {
  A <- fread(a); B <- fread(b); stopifnot(identical(sort(A$name), sort(B$name)))
  M <- merge(A, B, by = "name", suffixes = c("_cc", "_kcfc"))
  num <- function(x) suppressWarnings(as.numeric(x))
  M[, changed := !((is.na(num(est_cc)) & is.na(num(est_kcfc))) | abs(num(est_cc) - num(est_kcfc)) < 1e-12) |
                 !((is.na(num(lo_cc)) & is.na(num(lo_kcfc))) | abs(num(lo_cc) - num(lo_kcfc)) < 1e-12) |
                 !((is.na(num(hi_cc)) & is.na(num(hi_kcfc))) | abs(num(hi_cc) - num(hi_kcfc)) < 1e-12)]
  M[, class_change := class_cc != class_kcfc]
  if ("class_bobyqa_cc" %in% names(M)) M[, class_bobyqa_change := class_bobyqa_cc != class_bobyqa_kcfc] else M[, class_bobyqa_change := FALSE]
  M[, file := what]
  M[, .(file, name, est_cc, lo_cc, hi_cc, class_cc, est_kcfc, lo_kcfc, hi_kcfc, class_kcfc, changed, class_change, class_bobyqa_change)]
}
CD <- rbind(cmp("tables/variant_cc/V4_S5_claims.csv", "tables/variant_kcfc/V4_S5_claims.csv", "main grid"),
            cmp("tables/V4_S6_m1c_claims.csv", "tables/variant_kcfc/V4_S6_m1c_claims.csv", "M1c"))
fwrite(CD, "tables/V5_S7f_claims_diff.csv")
cat("sub-claims:", nrow(CD), "| value changed:", CD[, sum(changed)], "| class changed:", CD[, sum(class_change)],
    "| class changed with unchanged value (PF3):", CD[, sum(class_change & !changed)], "\n")
print(CD[changed == TRUE, .(file, name, est_cc, est_kcfc, class_cc, class_kcfc)], nrows = 200)
VA <- fread("tables/variant_cc/V4_S6_values.csv"); VB <- fread("tables/variant_kcfc/V4_S6_values.csv")
V <- merge(VA, VB, by = c("variant", "name"), suffixes = c("_cc", "_kcfc"))
V[, changed := abs(est_cc - est_kcfc) > 1e-12 | xor(is.na(est_cc), is.na(est_kcfc))]
fwrite(V[changed == TRUE], "tables/V5_S7f_values_diff.csv")
cat("values:", nrow(V), "| changed:", V[, sum(changed, na.rm = TRUE)], "| changed names without gee_kc:",
    V[changed == TRUE & !grepl("gee_kc", name), .N], "\n")
print(V[changed == TRUE & !grepl("gee_kc", name), .(variant, name, est_cc, est_kcfc)], nrows = 60)
# PF4: pooled KC rate on GenBinaryY, normal reference; conservative interval of the difference: +/- 2 (se_cc + se_kcfc)
k <- V[variant == "main" & name == "rate|GB|gee_kc"]
pf4 <- data.table(q = "rate|GB|gee_kc fewclust - geesmv", est = k$est_kcfc - k$est_cc,
                  lo = k$est_kcfc - k$est_cc - 2 * (k$se_cc + k$se_kcfc), hi = k$est_kcfc - k$est_cc + 2 * (k$se_cc + k$se_kcfc),
                  rate_cc = k$est_cc, rate_kcfc = k$est_kcfc)
pf4[, verdict := fifelse(hi < -1, "confermata", fifelse(lo > -1, "smentita", "indeterminata"))]
PF <- rbind(data.table(code = "PF3", est = CD[, sum(class_change & !changed)], lo = NA, hi = NA,
                       verdict = fifelse(CD[, sum(class_change & !changed)] == 0, "confermata", "smentita")),
            pf4[, .(code = "PF4", est, lo, hi, verdict)])
fwrite(PF, "tables/V5_S7f_PF.csv"); print(PF); print(pf4)
