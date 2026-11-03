# =============================================================================
# s6_checks.R -- S6: gates G2 and G5, predictions P1-P4, classification of the claims in
# M1c with the classification code of R/s5/s5_checks.R (evaluated verbatim from that file).
# Inputs : results_v4/S6/s6_boot.rds, results_v4/S5/claims_boot.rds, tables/V4_S5_values.csv,
#          tables/V4_S5_claims.csv, reports/morris/S5_claims_inventory.csv
# Outputs: tables/V4_S6_values.csv, V4_S6_m1c_claims.csv, V4_S6_m1c_switch_points.csv,
#          V4_S6_gates.csv, V4_S6_expectations.csv
# =============================================================================
source(path.expand("~/MMM/R/s6/s6_common.R"))
S6 <- readRDS(file.path(OUT_S6, "s6_boot.rds")); CB <- readRDS(file.path(OUT_S5, "claims_boot.rds"))
INV <- fread(file.path(REPO, "reports/morris/S5_claims_inventory.csv"))
src <- readLines(file.path(REPO, "R/s5/s5_checks.R"))
a <- grep("^# =+ classification =+$", src); b <- grep("^CLm <- classify", src)
eval(parse(text = src[(a + 1):(b - 1)]))                      # margin_of, summ, cls_diff, cls_equiv, classify
cat(sprintf("classification code: lines %d-%d of s5_checks.R\n", a + 1, b - 1))
G <- list()
# ---------------- G2a: main grid reproduces claims_boot (S5) ----------------
nm <- colnames(CB$main); A1 <- S6$main[, nm, drop = FALSE]; A0 <- CB$main
na_mis <- sum(xor(is.na(A1), is.na(A0))); dd <- abs(A1 - A0); mx <- max(dd[!is.na(dd)])
G$G2a <- data.table(gate = "G2a", pass = na_mis == 0 && mx <= 1e-12 && identical(dim(A1), dim(A0)),
                    result = sprintf("%d x %d values, max |diff| %.3g, NA mismatches %d", nrow(A1), ncol(A1), mx, na_mis))
# ---------------- G2b: M1c reproduces the 93 m1c values of S5 ----------------
V5 <- fread(file.path(TABDIR, "V4_S5_values.csv"))[variant == "m1c"]
g2b <- rbindlist(lapply(V5$name, function(z) { s <- summ(S6$m1c, z); v <- V5[name == z]
  data.table(name = z, d_est = abs(s[["est"]] - v$est), d_se = abs(s[["se"]] - v$se), d_lo = abs(s[["lo"]] - v$lo), d_hi = abs(s[["hi"]] - v$hi)) }))
mx2 <- max(as.matrix(g2b[, -1]))
G$G2b <- data.table(gate = "G2b", pass = nrow(g2b) == 93 && mx2 <= 1e-12, result = sprintf("%d values, max |diff| over est, se, lo, hi %.3g", nrow(g2b), mx2))
# ---------------- G2c: the copied classification reproduces S5 ----------------
CL5 <- fread(file.path(TABDIR, "V4_S5_claims.csv")); CLm <- classify(CB$main)
cmpc <- merge(CLm[, .(name, class6 = class)], CL5[, .(name, class5 = class)], by = "name")
G$G2c <- data.table(gate = "G2c", pass = nrow(cmpc) == nrow(CLm) && all(cmpc$class6 == cmpc$class5),
                    result = sprintf("%d of %d classes identical", sum(cmpc$class6 == cmpc$class5), nrow(CLm)))
print(rbindlist(G))
# ---------------- values: every engine value and extra, main and M1c ----------------
vals <- function(M, variant) { nms <- colnames(M)
  rbindlist(lapply(nms, function(z) { x <- unname(M[-1, z]); e <- unname(M[1, z])
    isarg <- grepl("^arg", z)
    data.table(variant = variant, name = z, est = e, se = if (isarg) NA_real_ else sd(x, na.rm = TRUE),
               lo = if (isarg) NA_real_ else unname(quantile(x, 0.025, na.rm = TRUE)), hi = if (isarg) NA_real_ else unname(quantile(x, 0.975, na.rm = TRUE)),
               n_finite = sum(is.finite(x)), p_same = if (isarg) mean(x == e, na.rm = TRUE) else NA_real_) })) }
VAL <- rbind(vals(S6$main, "main"), vals(S6$m1c, "m1c"), vals(S6$xmain, "main"), vals(S6$xm1c, "m1c"))
VAL <- VAL[!duplicated(VAL[, .(variant, name)])]
fwrite(VAL, file.path(TABDIR, "V4_S6_values.csv"))
# ---------------- M1c classification with the pre-registered applicability ----------------
NOT_APPL <- c(C07 = "usa T = 10, assente in M1c", C12b = "usa N = 145, assente in M1c", C24b = "usa N = 145, assente in M1c",
              C25a = "usa N >= 55, assente in M1c", C25b = "usa N >= 55, assente in M1c", C25c = "usa N >= 55, assente in M1c",
              C13a = "usa gli scenari M3-M5, assenti in M1c")
EXCL <- c(C04 = "conteggio delle celle vinte fuori dal manoscritto (D-25)", C05 = "conteggio delle celle vinte fuori dal manoscritto (D-25)")
cn <- grep("^C[0-9]", colnames(S6$m1c), value = TRUE)
nonfin <- cn[!is.finite(S6$m1c[1, cn])]                      # estimate not defined on the M1c grid
cat("M1c claim columns with a non-finite estimate:", paste(nonfin, collapse = ", "), "\n")
CM <- classify(S6$m1c[, setdiff(colnames(S6$m1c), nonfin), drop = FALSE])
CM <- rbind(CM, data.table(name = nonfin, claim_id = sub("\\|.*$", "", nonfin), est = unname(S6$m1c[1, nonfin]), class = "non calcolabile"), fill = TRUE)
CM <- merge(CM, CL5[, .(name, class_main = class, est_main = est, lo_main = lo, hi_main = hi)], by = "name", all.x = TRUE)
CM[, applicability := "applicabile"]; CM[, reason := ""]
CM[claim_id %in% names(NOT_APPL), `:=`(applicability = "non applicabile", reason = NOT_APPL[claim_id])]
CM[claim_id %in% names(EXCL), `:=`(applicability = "esclusa", reason = EXCL[claim_id])]
ndp <- S6$m1c[1, grep("^ndp\\|", colnames(S6$m1c))]; names(ndp) <- sub("^ndp\\|", "", names(ndp))
CM[, stratum := fifelse(grepl("^C15[a-d]\\|s[1-5]$", name), sub("^.*\\|", "", name), NA_character_)]
CM[!is.na(stratum) & ndp[stratum] == 0, `:=`(applicability = "non applicabile", reason = "strato vuoto in M1c")]
CM[!is.na(stratum), ndp_m1c := ndp[stratum]]
CM[applicability == "applicabile" & class == "non calcolabile", reason := "stima non definita sulla griglia di M1c"]
CM[applicability == "applicabile" & n_finite < 1990 & class != "non calcolabile", reason := sprintf("valori finiti in %d replicati su 2.000", n_finite)]
CM[applicability != "applicabile", class := applicability]
bf <- file.path(OUT_S6, "s6_boot_m1c_boby.rds")
if (file.exists(bf)) {                                       # D-16 sensitivity in M1c
  MB <- readRDS(bf); nb <- grep("^C[0-9]", colnames(MB), value = TRUE); nb <- nb[is.finite(MB[1, nb])]
  CB2 <- classify(MB[, c(setdiff(colnames(MB), grep("^C[0-9]", colnames(MB), value = TRUE)), nb), drop = FALSE])
  CM <- merge(CM, CB2[, .(name, est_bobyqa = est, lo_bobyqa = lo, hi_bobyqa = hi, class_bobyqa = class)], by = "name", all.x = TRUE)
  CM[, bobyqa_changes := applicability == "applicabile" & !is.na(class_bobyqa) & class_bobyqa != class]
  cat(sprintf("M1c bobyqa sensitivity: %d applicable sub-claims change class\n", CM[bobyqa_changes == TRUE, .N]))
}
CM <- merge(CM, INV[, .(claim_id, section, claim)], by = "claim_id", all.x = TRUE)
NC <- INV[!claim_id %in% CM$claim_id, .(claim_id, name = claim_id, kind, class = kind, applicability = kind,
         reason = fifelse(claim_id %in% names(EXCL), EXCL[claim_id], "riga non calcolata anche in S5"), section, claim)]
NC[claim_id %in% names(EXCL), `:=`(class = "esclusa", applicability = "esclusa")]
CM <- rbind(CM, NC, fill = TRUE); setorder(CM, claim_id, name)
fwrite(CM, file.path(TABDIR, "V4_S6_m1c_claims.csv"))
print(CM[, .N, by = .(applicability, class)][order(applicability, -N)])
# switch points in M1c (C23a/b), as in S5
SWm <- rbindlist(lapply(c("GB", "NAT", "AR1"), function(Gm) { a <- unname(S6$m1c[, paste0("C23a|", Gm)]); bb <- unname(S6$m1c[, paste0("C23b|", Gm)])
  tb <- table(a[-1]); data.table(mech = Gm, switch_m1c = a[1], modal_switch = as.numeric(names(which.max(tb))), p_modal = max(tb) / sum(is.finite(a[-1])),
                                 p_same = mean(a[-1] == a[1], na.rm = TRUE), next_m1c = if (is.finite(bb[1])) PRIM8[bb[1]] else NA_character_) }))
print(SWm); fwrite(SWm, file.path(TABDIR, "V4_S6_m1c_switch_points.csv"))
# ---------------- G5: completeness ----------------
cov_ids <- unique(CM$claim_id); miss <- setdiff(INV$claim_id, cov_ids)
nofin <- CM[applicability == "applicabile" & is.finite(n_finite) & n_finite < 1990 & reason == ""]
G$G5 <- data.table(gate = "G5", pass = length(miss) == 0 && nrow(nofin) == 0,
                   result = sprintf("%d of %d inventory rows with an outcome or a reason; applicable sub-claims %d, of which %d below 1,990 finite replicates (all with the reason recorded)",
                                    length(intersect(INV$claim_id, cov_ids)), nrow(INV), CM[applicability == "applicabile", .N], CM[applicability == "applicabile" & n_finite < 1990, .N]))
# ---------------- predictions ----------------
iv <- function(z) CM[name == z, .(est, lo, hi)]
below <- function(z, m) { x <- iv(z); if (x$hi < m) "confermata" else if (x$lo >= m) "smentita" else "indeterminata" }
p1 <- iv("C14a|all"); p3 <- iv("C17b|all")
P <- rbindlist(list(
  data.table(id = "P1", quantity = "C14a|all (M1c)", expected = "circa 0, sotto 5 punti", est = p1$est, lo = p1$lo, hi = p1$hi, outcome = below("C14a|all", 5)),
  rbindlist(lapply(c("GB", "NAT", "AR1"), function(Gm) { z <- paste0("C23c|", Gm); x <- iv(z)
    data.table(id = "P2", quantity = paste0(z, " (M1c)"), expected = "sotto 5 punti", est = x$est, lo = x$lo, hi = x$hi, outcome = below(z, 5)) })),
  data.table(id = "P3", quantity = "C17b|all (M1c)", expected = "1,39 +/- 0,5", est = p3$est, lo = p3$lo, hi = p3$hi,
             outcome = if (p3$lo > 0.889 && p3$hi < 1.889) "confermata" else if (p3$hi < 0.889 || p3$lo > 1.889) "smentita" else "indeterminata"),
  rbindlist(lapply(c("GB", "NAT", "AR1"), function(Gm) { z <- paste0("C19c|", Gm); x <- iv(z)
    data.table(id = "P4", quantity = paste0(z, " (M1c)"), expected = "sotto 0,10", est = x$est, lo = x$lo, hi = x$hi, outcome = below(z, 0.10)) }))))
print(P); fwrite(P, file.path(TABDIR, "V4_S6_expectations.csv"))
fwrite(rbindlist(G), file.path(TABDIR, "V4_S6_gates.csv"))
