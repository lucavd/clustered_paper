# =============================================================================
# s5_checks.R -- S5: gates G2-G4, predictions P1-P4, classification of every claim (D-23).
# Inputs : results_v4/S5/agg.rds, results_v4/S5/claims_boot.rds, tables/V4_perf_all.csv.gz,
#          v3 cell tables (tables/O_level_all.csv, O_coverage_all.csv, O_reference_dist.csv, rev2/REV_cells.csv),
#          reports/morris/S5_claims_inventory.csv
# Outputs: tables/V4_S5_*.csv
# =============================================================================
source(path.expand("~/MMM/R/s5/s5_common.R")); source(file.path(REPO, "R/s5/s5_engine.R"))
t0 <- Sys.time()
AG <- readRDS(file.path(OUT_S5, tag_out("agg.rds"))); META <- AG$META; B <- AG$B
CB <- readRDS(file.path(OUT_S5, tag_out("claims_boot.rds")))
PBIN <- readRDS(file.path(OUT_S5, "pbin.rds"))
TO <- function(x) if (DEV) file.path(OUT_S5, tag_out(x)) else file.path(TABDIR, x)
INV <- fread(file.path(REPO, "reports/morris/S5_claims_inventory.csv"))

# ============================== G2: same data, same indices ==============================
P <- fread(file.path(TABDIR, "V4_perf_all.csv.gz"))[usab == "U0"]
P[, arm := fifelse(label == "mmlb", "H1", "H0")]
M <- META[set %in% c("own", "intersection")]
M[, `:=`(n0 = AG$A$n[1, j], k0 = AG$A$k[1, j], kc0 = AG$A$kc[1, j], nd0 = AG$ND[1, d])]
ksa_col <- match(M$j, as.integer(colnames(AG$KSA)))
M[, ksa0 := AG$KSA[1, ksa_col]]
if (!DEV) {
  M[, se_raw := apply(AG$A$k[-1, j, drop = FALSE] / AG$A$n[-1, j, drop = FALSE], 2, sd, na.rm = TRUE)]
  ii <- which(!is.na(ksa_col)); se_sa <- rep(NA_real_, nrow(M))
  se_sa[ii] <- apply(AG$KSA[-1, ksa_col[ii], drop = FALSE] / AG$A$n[-1, M$j[ii], drop = FALSE], 2, sd, na.rm = TRUE)
  M[, se_sa := se_sa]   # (first pass indexed j[ii] inside M[ii, ...], i.e. twice: a tool error, fixed here)
}
G2 <- merge(M, P, by = c("dgp", "scen", "N", "ntime", "arm", "method", "set"))
cmp <- function(a, b) { ok <- is.finite(a) | is.finite(b); x <- abs(a[ok] - b[ok]); x[is.na(x)] <- Inf; c(n = sum(ok), max = if (any(ok)) max(x) else 0) }
g2 <- rbind(
  data.table(check = "n_obs", t(cmp(G2$n0, G2$n_obs))),
  data.table(check = "rej", t(cmp(fifelse(G2$n0 > 0, G2$k0 / G2$n0, NA_real_), G2$rej))),
  data.table(check = "p_usable", t(cmp(G2$n0 / G2$nd0, G2$p_usable))),
  data.table(check = "cover", t(cmp(G2[is.finite(cover), kc0 / n0], G2[is.finite(cover), cover]))),
  data.table(check = "power_sa (H1 own)", t(cmp(G2[arm == "H1" & set == "own", ksa0 / n0], G2[arm == "H1" & set == "own", power_sa]))))
if (!DEV) g2 <- rbind(g2,
  data.table(check = "boot_se_power_raw", t(cmp(G2[is.finite(boot_se_power_raw), se_raw], G2[is.finite(boot_se_power_raw), boot_se_power_raw]))),
  data.table(check = "mcse_power_sa (H1 own)", t(cmp(G2[arm == "H1" & set == "own", se_sa], G2[arm == "H1" & set == "own", mcse_power_sa]))))
g2[, pass := max <= 1e-12]
cat("cells matched:", nrow(G2), "of", nrow(M), "\n"); print(g2)

# ============================== G3: linearization vs bootstrap ==============================
L <- rbindlist(lapply(seq_along(AG$LIN), function(m) AG$LIN[[m]][, model := m]))
lin_ratio <- function(a_list, u_list) {        # a, u: lists of per-model 1000-vectors
  R <- sum(unlist(a_list)) / sum(unlist(u_list))
  v <- sum(mapply(function(a, u) { z <- a - R * u; sum((z - mean(z))^2) }, a_list, u_list))
  c(est = 100 * R, se = 100 * sqrt(v) / sum(unlist(u_list))) }
mainB <- CB$main
bse <- function(nm) sd(mainB[-1, nm], na.rm = TRUE)
g3 <- list()
for (g in DGPS) for (mm in PRIM8) {
  x <- L[dgp == g & method == mm & set == "own"]; r <- lin_ratio(x$a, x$u)
  nm <- sprintf("rate|%s|%s", SHORT[[g]], mm)
  g3[[length(g3) + 1]] <- data.table(quantity = nm, est_lin = r[["est"]], est_boot = mainB[1, nm], se_lin = r[["se"]], se_boot = bse(nm)) }
for (g in DGPS) {
  a <- L[dgp == g & method == "mmm" & set == "pair"]; f <- L[dgp == g & method == "gee_fg" & set == "pair"]
  setorder(a, model); setorder(f, model)
  dl <- mapply(function(x, y) x - y, a$a, f$a, SIMPLIFY = FALSE)
  r <- lin_ratio(dl, a$u); nm <- paste0("C03|", SHORT[[g]])
  g3[[length(g3) + 1]] <- data.table(quantity = nm, est_lin = r[["est"]], est_boot = mainB[1, nm], se_lin = r[["se"]], se_boot = bse(nm)) }
g3 <- rbindlist(g3)[, rel := se_boot / se_lin - 1][, pass := abs(rel) <= 0.07 & abs(est_lin - est_boot) < 1e-9]
print(g3)
fwrite(g3, TO("V4_S5_g3_linearization.csv"))

# ============================== P3: weight of the common random numbers ==============================
cl0 <- META[, .(d, dgp, scen, N, ntime, arm, method, set)][, `:=`(n = AG$A$n[1, ], k = AG$A$k[1, ])]
ind <- cl0[scen != "M1c" & arm == "H0" & set == "own" & method %in% PRIM8 & n > 0,
           .(se_ind = 100 * sqrt(sum(k * (n - k) / n)) / sum(n)), by = .(dgp, method)]
p3 <- merge(g3[grepl("^rate", quantity)][, .(quantity, se_boot)],
            ind[, .(quantity = sprintf("rate|%s|%s", SHORT[dgp], method), se_ind)], by = "quantity")
mc <- AG$MCN[, .(quantity = paste0("C03|", SHORT[dgp]), se_ind = 100 * sqrt((b + c) - (b - c)^2 / n) / n)]
p3 <- rbind(p3, merge(g3[grepl("^C03", quantity)][, .(quantity, se_boot)], mc, by = "quantity"))
p3[, ratio := se_boot / se_ind][, confirmed := ratio > 1.25]
print(p3); fwrite(p3, TO("V4_S5_p3_crn.csv"))

# ============================== P1: v3 definitions reproduce the manuscript ==============================
REV <- fread(file.path(REPO, "rev2/REV_cells.csv")); REV[, ord := as.numeric(.I)]
mkv3 <- function(src) {   # src: data.table with dgp, model, N, ntime, label/arm, method, set, n, k, (kt, kc, nc), ord
  att <- REV[, .(nd = max(n_attempt)), by = .(dgp, model, N, ntime)]
  grid <- CJ(dgp = DGPS, model = 1:5, N = sort(unique(REV$N)), ntime = sort(unique(REV$ntime)), arm = c("H0", "H1"),
             method = PRIM8, set = c("own", "intersection"))
  x <- merge(grid, src, by = c("dgp", "model", "N", "ntime", "arm", "method", "set"), all.x = TRUE)
  x <- merge(x, att, by = c("dgp", "model", "N", "ntime"), all.x = TRUE)
  for (cc in intersect(c("n", "k", "kt", "kc", "nc", "ord"), names(x))) set(x, NULL, cc, as.numeric(x[[cc]]))
  for (cc in c("n", "k", "kt", "kc", "nc")) if (cc %in% names(x)) set(x, which(is.na(x[[cc]])), cc, 0) else set(x, NULL, cc, NA_real_)
  x[is.na(ord), ord := Inf]; setorder(x, ord)
  x[, `:=`(scen = paste0("M", model), tok = TRUE, e2 = NA_real_, ksa = NA_real_)]
  x }
r3 <- REV[, .(dgp, model, N, ntime, arm = fifelse(label == "mmlb", "H1", "H0"), method, set, n = as.numeric(n), k = as.numeric(k_rej), kt = as.numeric(k_rej_t), kc = as.numeric(k_cov), nc = as.numeric(n_cov), ord)]
E_rev <- engine(mkv3(r3), tie = "first", PBIN = PBIN, v3 = TRUE)
OL <- fread(file.path(TABDIR, "O_level_all.csv"))[, .(dgp, model, N, ntime, arm = "H0", method, set, n = as.numeric(n), k = as.numeric(k_rej), ord = as.numeric(.I))]
E_ol <- tryCatch(engine(mkv3(OL), tie = "first", PBIN = PBIN, v3 = TRUE), error = function(e) { cat("O_level engine:", conditionMessage(e), "\n"); NULL })
OC <- fread(file.path(TABDIR, "O_coverage_all.csv"))[set == "own", .(dgp, model, N, ntime, arm = "H1", method, set, nc = n, kc = cov_gen * n / 100)]
r3c <- copy(r3); r3c[arm == "H1" & set == "own", `:=`(kc = NA_real_, nc = NA_real_)]
r3c[OC, on = c("dgp", "model", "N", "ntime", "arm", "method", "set"), `:=`(kc = i.kc, nc = i.nc)]
r3c[arm == "H1" & set == "own" & is.na(nc), `:=`(kc = 0, nc = 0)]
E_oc <- engine(mkv3(r3c), tie = "first", PBIN = PBIN, v3 = TRUE)
ORD <- fread(file.path(TABDIR, "O_reference_dist.csv"))
ord_wct <- ORD[, .(wct = max(abs(tref - 5))), by = .(dgp, method)]
ord_rt <- ORD[, .(rt = sum(tref * n) / sum(n)), by = .(dgp, method)]

M3 <- c(GB = "GB", NAT = "NAT", AR1 = "AR1")
tg <- list(); addt <- function(name, printed, digits, src) tg[[length(tg) + 1]] <<- data.table(name = name, printed = printed, digits = digits, source = src)
three <- function(kind, mm, vals, d, src) for (i in 1:3) addt(sprintf("%s|%s|%s", kind, M3[i], mm), vals[i], d, src)
T31 <- list(gee_fg = list(c(5.79,5.32,5.34), c(4.55,2.83,2.91), c(1.00,0.52,0.55)), gee_md = list(c(6.05,5.51,5.56), c(8.33,5.90,6.16), c(1.41,0.88,0.97)),
            mmm = list(c(6.59,6.19,6.38), c(2.06,1.50,1.55), c(1.61,1.20,1.37)), glmer = list(c(8.51,6.30,6.54), c(4.43,1.83,2.18), c(3.38,1.38,1.57)),
            shrink = list(c(9.17,6.64,6.90), c(4.86,2.30,2.49), c(4.09,1.61,1.85)), gee = list(c(8.31,7.61,7.69), c(20.90,17.40,18.36), c(3.80,3.00,3.11)),
            glmer_lrt = list(c(11.06,7.48,8.11), c(7.16,3.85,4.09), c(6.06,2.48,3.11)), gee_kc = list(c(13.10,10.43,10.71), c(9.23,6.48,6.72), c(7.93,5.32,5.57)))
for (mm in names(T31)) { three("rate", mm, T31[[mm]][[1]], 2, "3.1"); three("wc", mm, T31[[mm]][[2]], 2, "3.1"); three("mad", mm, T31[[mm]][[3]], 2, "3.1") }
for (x in list(c("mmm",321), c("gee_fg",308), c("glmer",216), c("gee_md",210), c("gee",172), c("glmer_lrt",117), c("gee_kc",57), c("shrink",39))) addt(paste0("wins|", x[1]), as.numeric(x[2]), 0, "3.1")
MT <- list(mmm = c(1.80,1.50,1.74), gee = c(2.17,3.55,3.78), gee_fg = c(2.42,1.02,1.18), gee_kc = c(2.47,3.45,7.83), gee_md = c(2.53,1.14,1.60),
           glmer = c(5.59,3.70,2.06), shrink = c(10.83,4.66,2.15), glmer_lrt = c(11.09,3.68,2.25))
for (mm in names(MT)) for (i in 1:3) addt(sprintf("madT|%s|%s", c("T2","T3","T4plus")[i], mm), MT[[mm]][i], 2, "3.1")
addt("C07|all", -0.357, 3, "V-90"); addt("C08|wc", 0.976, 3, "V-92"); addt("C08|mad", 0.833, 3, "V-92"); addt("C08|rate", 0.857, 3, "V-92"); addt("C08|madcell", 0.833, 3, "V-92"); addt("C09|all", 0.217, 3, "V-91")
addt("C10a|raw", 52.25, 2, "V-81"); addt("C10a|BH", 47.13, 2, "V-81"); addt("C10a|bonferroni", 27.82, 2, "V-81"); addt("C10a|ncells", 11163, 0, "3.2")
for (x in list(c("gee_fg",23.2), c("gee_md",28.1), c("mmm",29.1), c("gee",45.4), c("glmer",50.2), c("shrink",52.7), c("glmer_lrt",60.6), c("gee_kc",87.2))) addt(paste0("flag|BH|", x[1]), as.numeric(x[2]), 1, "3.2")
SH <- list(glmer_lrt = c(-2.13,-1.78,-1.84), gee_kc = c(0.45,0.51,0.48), shrink = c(-0.69,-0.35,-0.36), gee = c(-0.53,-0.19,-0.26), mmm = c(-0.32,-0.17,-0.17),
           gee_md = c(-0.03,0.15,0.12), gee_fg = c(-0.06,0.12,0.08), glmer = c(-0.02,-0.01,0.00))
for (mm in names(SH)) three("shift", mm, SH[[mm]], 2, "3.3")
U9 <- list(glmer_lrt = c(98.7,89.3), gee = c(78.7,65.2), gee_fg = c(69.0,63.9), gee_kc = c(68.8,62.3), gee_md = c(68.7,62.4), shrink = c(66.9,62.7), mmm = c(66.6,62.3), glmer = c(64.5,61.0))
for (mm in names(U9)) { addt(paste0("use20|", mm), U9[[mm]][1], 1, "3.9"); addt(paste0("val20|", mm), U9[[mm]][2], 1, "3.9") }
MS <- list(mmm = c(2.62,1.56,2.00,2.34,1.02), gee_kc = c(3.06,1.95,3.07,6.15,8.15), gee_fg = c(3.89,2.83,2.75,1.13,0.63), gee_md = c(4.24,4.25,4.12,1.32,0.61),
           glmer = c(5.46,2.22,2.47,2.66,2.38), glmer_lrt = c(5.35,8.79,8.36,4.34,1.59), gee = c(8.36,9.95,8.69,2.79,1.04), shrink = c(10.92,5.37,3.88,2.94,2.40))
for (mm in names(MS)) for (i in 1:5) addt(sprintf("madS|s%d|%s", i, mm), MS[[mm]][i], 2, "3.10")
for (i in 1:5) { addt(sprintf("ndp|s%d", i), c(83,100,147,502,608)[i], 0, "3.10"); addt(sprintf("gapS|s%d", i), c(-1.27,-1.27,-0.75,1.21,0.40)[i], 2, "3.10") }
addt("C15e|all", 0.245, 3, "3.10")
WT <- list(gee_fg = c(0.80,1.08,0.95), mmm = c(2.02,2.50,2.63), gee_md = c(3.05,1.16,1.39), glmer = c(4.23,2.67,2.66), shrink = c(4.43,1.90,2.19), gee_kc = c(8.85,6.14,6.34), gee = c(12.87,9.64,10.60))
for (mm in names(WT)) three("wct", mm, WT[[mm]], 2, "V-83")
RT <- list(gee = c(7.16,6.49,6.56), gee_kc = c(12.00,9.45,9.72), gee_md = c(5.24,4.73,4.82), mmm = c(5.79,5.39,5.58), gee_fg = c(5.03,4.58,4.63), shrink = c(8.48,6.10,6.32), glmer = c(7.86,5.77,5.97))
for (mm in names(RT)) three("ratet", mm, RT[[mm]], 2, "3.12")
CV <- list(gee_fg = c(94.55,94.83,94.91), gee_md = c(94.24,94.59,94.64), mmm = c(93.17,93.24,93.24), glmer = c(88.32,93.03,92.28), gee = c(91.90,92.44,92.47), shrink = c(90.06,92.94,92.59), gee_kc = c(86.36,88.94,88.84))
for (mm in names(CV)) three("cov", mm, CV[[mm]], 2, "3.15")
for (x in list(c("gee_fg",2.33), c("gee_md",3.67), c("gee",4.00), c("shrink",4.33), c("glmer",4.67), c("mmm",4.67), c("glmer_lrt",5.00), c("gee_kc",7.33))) addt(paste0("crank|", x[1]), as.numeric(x[2]), 2, "V-67")
for (x in list(c("glmer",2.67), c("glmer_lrt",2.67), c("shrink",3.33), c("gee_fg",3.67))) addt(paste0("crank25|", x[1]), as.numeric(x[2]), 2, "V-67")
TG <- rbindlist(tg)
getv <- function(E, nm) if (is.null(E) || !nm %in% names(E)) NA_real_ else E[[nm]]
alt_t <- function(nm) { p <- strsplit(nm, "|", fixed = TRUE)[[1]]; g <- names(SHORT)[SHORT == p[2]]
  if (p[1] == "wct") ord_wct[dgp == g & method == p[3], wct] else if (p[1] == "ratet") ord_rt[dgp == g & method == p[3], rt] else NA_real_ }
TG[, `:=`(rev = vapply(name, function(z) getv(E_rev, z), 0), olev = vapply(name, function(z) getv(E_ol, z), 0),
          ocov = vapply(name, function(z) getv(E_oc, z), 0), oref = vapply(name, function(z) { v <- alt_t(z); if (length(v)) v else NA_real_ }, 0))]
TG[, primary := fifelse(source %in% c("3.1", "3.3", "V-90", "V-92", "V-91", "V-81", "3.2"), olev,
                fifelse(source == "3.15", ocov, fifelse(source %in% c("V-83", "3.12"), oref, rev)))]
tol <- function(d) 0.5 * 10^(-d) + 1e-9
TG[, status := fifelse(is.na(primary), "not computed", fifelse(abs(primary - printed) <= tol(digits), "reproduced",
               fifelse(abs(primary - printed) <= 10^(-digits) + 1e-9, "one unit of the last digit", "not reproduced")))]
TG[, rev_status := fifelse(is.na(rev), "", fifelse(abs(rev - printed) <= tol(digits), "reproduced", "differs"))]
print(TG[, .N, by = .(source, status)][order(source)])
print(TG[status != "reproduced"], nrows = 200)
fwrite(TG, TO("V4_S5_p1_reproduction.csv"))

# ---- P1 diagnostics: C08 on the MAD under three weightings (O_level_all, own) ------------------------
OLd <- OL[set == "own"][, dev := abs(100 * k / n - 5)]
madv <- list(perN = OLd[, .(n = sum(n), k = sum(k)), by = .(dgp, method, N)][, .(v = mean(abs(100 * k / n - 5))), by = .(dgp, method)],
             perCell = OLd[, .(v = mean(dev)), by = .(dgp, method)],
             perNT = OLd[, .(n = sum(n), k = sum(k)), by = .(dgp, method, N, ntime)][, .(v = mean(abs(100 * k / n - 5))), by = .(dgp, method)])
c08d <- rbindlist(lapply(names(madv), function(w) { x <- madv[[w]]
  pr <- list(c("GenBinaryY", "GLMM-native"), c("GenBinaryY", "GLMM-AR1"), c("GLMM-native", "GLMM-AR1"))
  r <- sapply(pr, function(pp) { a <- x[dgp == pp[1]][order(method)]; b <- x[dgp == pp[2]][order(method)]; cor(a$v, b$v, method = "spearman") })
  data.table(weighting = w, pair = sapply(pr, paste, collapse = " vs "), rho = r) }))
print(c08d); fwrite(c08d, TO("V4_S5_p1_c08_diagnostic.csv"))


# ============================== classification ==============================
margin_of <- function(id) {
  mg <- INV[claim_id == id, margin]
  if (!length(mg) || mg == "") return(NA_real_)
  if (grepl("^1 pp", mg)) 1 else if (grepl("^2 pp", mg)) 2 else if (grepl("^5 points", mg)) 5 else if (grepl("^72", mg)) 72 else
  if (grepl("10% relative", mg)) 0.10 else if (grepl("0.3", mg, fixed = TRUE)) 0.3 else if (grepl("0.1", mg, fixed = TRUE)) 0.1 else
  if (grepl("^Bradley 2.5", mg)) 2.5 else if (grepl("^1 point", mg)) 1 else NA_real_ }
summ <- function(MAT, nm) { x <- unname(MAT[-1, nm]); c(est = unname(MAT[1, nm]), se = sd(x, na.rm = TRUE), lo = unname(quantile(x, 0.025, na.rm = TRUE)),
                                              hi = unname(quantile(x, 0.975, na.rm = TRUE)), n_fin = sum(is.finite(x)), p_pos = mean(x > 0, na.rm = TRUE)) }
cls_diff <- function(lo, hi, m) { if (!is.finite(lo) || !is.finite(hi)) return("non calcolabile"); if (is.na(m)) m <- 0
  if (hi < 0) return(if (hi < -m && m > 0) "smentita e rilevante" else "smentita")
  if (lo > 0) return(if (m == 0) "risolta" else if (lo > m) "risolta e rilevante" else if (hi < m) "trascurabile" else "risolta, rilevanza indeterminata")
  if (lo > -m && hi < m) "trascurabile" else "non risolta" }
cls_equiv <- function(lo, hi, m) { if (!is.finite(lo) || !is.finite(hi)) return("non calcolabile")
  if (lo > -m && hi < m) "confermata" else if (lo >= m || hi <= -m) "smentita" else "indeterminata" }
classify <- function(MAT) {
  nms <- grep("^C[0-9]", colnames(MAT), value = TRUE)
  rows <- lapply(nms, function(nm) {
    id <- strsplit(nm, "|", fixed = TRUE)[[1]][1]; inv <- INV[claim_id == id]; kind <- if (nrow(inv)) inv$kind else "descriptive"
    if (grepl("rank_", nm, fixed = TRUE)) kind <- "descriptive"
    s <- summ(MAT, nm); m <- margin_of(id); sg <- if (id == "C07") -1 else 1
    est <- sg * s[["est"]]; lo0 <- if (sg < 0) -s[["hi"]] else s[["lo"]]; hi0 <- if (sg < 0) -s[["lo"]] else s[["hi"]]
    one_class <- function(lo, hi) switch(kind,
      difference = cls_diff(lo, hi, m),
      equivalence = cls_equiv(lo, hi, m),
      value = if (id == "C08") { if (lo >= 0.7) "confermata" else if (hi < 0.7) "smentita" else "indeterminata" }
              else if (id %in% c("C12a", "C24b")) cls_diff(lo, hi, 0)
              else if (id %in% c("C12b", "C12c")) cls_diff(lo, hi, 2.5)
              else if (id == "C13a") { if (hi < 1) "confermata" else if (lo >= 1) "smentita" else "indeterminata" }
              else if (id == "C24a") { if (lo >= 2.5 && hi <= 7.5) "confermata" else if (hi < 2.5 || lo > 7.5) "smentita" else "indeterminata" }
              else if (id == "C24c") { if (lo >= 3.6 && hi <= 6.4) "confermata" else if (hi < 3.6 || lo > 6.4) "smentita" else "indeterminata" }
              else "descrittiva",
      rank = if (id == "C21a") cls_diff(lo, hi, 0) else "vedi tabella dei punti di cambio",
      "descrittiva")
    lo <- lo0; hi <- hi0; cl <- one_class(lo0, hi0)
    blo <- 2 * est - hi0; bhi <- 2 * est - lo0; clb <- one_class(blo, bhi)     # basic (reverse percentile) interval, sensitivity
    data.table(name = nm, claim_id = id, kind = kind, est = est, se = s[["se"]], lo = lo, hi = hi, n_finite = s[["n_fin"]], margin = m, class = cl,
               boot_mean = sg * mean(MAT[-1, nm], na.rm = TRUE), basic_lo = blo, basic_hi = bhi, class_basic = clb) })
  X <- rbindlist(rows)
  X[, p := fifelse(kind == "difference" & se > 0, 2 * pnorm(-abs(est / se)), NA_real_)]
  ii <- which(is.finite(X$p)); X[, `:=`(p_BH = NA_real_, p_BY = NA_real_)]
  X[ii, p_BH := p.adjust(p, "BH")]; X[ii, p_BY := p.adjust(p, "BY")]
  X[, flag := fifelse(grepl("^risolta", class) & is.finite(p_BY) & p_BY >= 0.05, "risolta solo senza correzione", "")]
  # replicates per cell that would resolve the claim (2 MCSE = distance to the boundary that the interval crosses)
  X[, boundary := NA_real_]
  X[kind == "difference" & class == "non risolta", boundary := 0]
  X[kind == "difference" & class == "risolta, rilevanza indeterminata", boundary := margin]
  X[kind == "equivalence" & class == "indeterminata", boundary := fifelse(hi >= margin, margin, -margin)]
  X[, n_req := fifelse(is.finite(boundary) & abs(est - boundary) > 0, 1000 * (2 * se / abs(est - boundary))^2, NA_real_)]
  X }
CLm <- classify(CB$main); CLb <- classify(CB$boby)
CL <- merge(CLm, CLb[, .(name, est_bobyqa = est, lo_bobyqa = lo, hi_bobyqa = hi, class_bobyqa = class)], by = "name", all.x = TRUE)
CL[, bobyqa_changes := !is.na(class_bobyqa) & class_bobyqa != class]
CL <- merge(CL, INV[, .(claim_id, section, claim, v3_value)], by = "claim_id", all.x = TRUE)
setorder(CL, claim_id, name)
# switch points (C23a/b) and argmax spread (C22)
claimed_sw <- c(GB = 35, NAT = 45, AR1 = 45); md_code <- match("gee_md", PRIM8)
SW <- rbindlist(lapply(names(claimed_sw), function(G) { a <- unname(CB$main[, paste0("C23a|", G)]); bb <- unname(CB$main[, paste0("C23b|", G)])
  data.table(mech = G, switch_v4 = a[1], claimed = claimed_sw[[G]], p_claimed = mean(a[-1] == claimed_sw[[G]], na.rm = TRUE),
             modal_switch = as.numeric(names(which.max(table(a[-1])))), p_modal = max(table(a[-1])) / sum(is.finite(a[-1])),
             next_v4 = PRIM8[bb[1]], p_next_md = mean(bb[-1] == md_code, na.rm = TRUE),
             wc_argmax_distinctN = length(unique(CB$main[-1, paste0("C22|", G)]))) }))
print(SW); fwrite(SW, TO("V4_S5_switch_points.csv"))
# C05: winner vs runner-up per design point
C5 <- copy(CB$pairs); C5[, `:=`(gap = CB$c05[1, ], se = apply(CB$c05[-1, , drop = FALSE], 2, sd, na.rm = TRUE))]
C5[, p := fifelse(se > 0, 2 * pnorm(-abs(gap / se)), 1)][, `:=`(p_BH = p.adjust(p, "BH"), p_BY = p.adjust(p, "BY"))]
c05 <- data.table(name = c("C05|uncorrected", "C05|BH", "C05|BY"), claim_id = "C05", kind = "count",
                  est = c(sum(C5$p < 0.05), sum(C5$p_BH < 0.05), sum(C5$p_BY < 0.05)))
c05[, class := fifelse(name == "C05|uncorrected", "descrittiva", fifelse(est == 0, "confermata", "smentita"))]
CL <- rbind(CL, merge(c05, INV[, .(claim_id, section, claim, v3_value)], by = "claim_id"), fill = TRUE)
fwrite(C5, TO("V4_S5_c05_winner_runnerup.csv"))
# non-computed rows of the inventory
NC <- INV[kind %in% c("superseded", "out_of_scope", "not_operationalized"), .(claim_id, name = claim_id, kind, class = kind, section, claim, v3_value, notes)]
CL <- rbind(CL, NC, fill = TRUE)
fwrite(CL, TO("V4_S5_claims.csv"))
print(CL[, .N, by = class][order(-N)])

# ============================== secondary values with MCSE ==============================
VAL <- rbindlist(lapply(c("main", "boby", "m1c"), function(v) { MAT <- CB[[v]]
  rbindlist(lapply(colnames(MAT), function(nm) { s <- summ(MAT, nm); data.table(variant = v, name = nm, est = s[["est"]], se = s[["se"]], lo = s[["lo"]], hi = s[["hi"]], n_finite = s[["n_fin"]]) })) }))
fwrite(VAL, TO("V4_S5_values.csv"))

# ============================== P2: stability v3 -> v4 (7 non-LRT methods) ==============================
P2 <- TG[source == "3.1" & grepl("^(rate|wc|mad)\\|", name) & !grepl("glmer_lrt", name), .(name, v3 = printed)]
P2 <- merge(P2, VAL[variant == "main", .(name, v4 = est, se)], by = "name")
P2[, `:=`(diff = v4 - v3, lo = v4 - v3 - 2 * se, hi = v4 - v3 + 2 * se)]
P2[, class := fifelse(lo > -1 & hi < 1, "confermata", fifelse(lo >= 1 | hi <= -1, "smentita", "indeterminata"))]
p2agg <- P2[, .(mean_diff = mean(diff), mean_abs = mean(abs(diff)), max_abs = max(abs(diff)), n_conf = sum(class == "confermata"), n = .N)]
print(p2agg); print(P2[class != "confermata"])
fwrite(P2, TO("V4_S5_p2_stability.csv"))
lrt <- TG[source == "3.1" & grepl("glmer_lrt", name), .(name, v3 = printed)]
lrt <- merge(lrt, VAL[variant == "main", .(name, v4 = est, se)], by = "name"); print(lrt)

# ============================== G4 and the gate table ==============================
inv_ids <- INV$claim_id
g4_missing <- setdiff(inv_ids, unique(CL$claim_id))
g4_nf <- CL[kind %in% c("difference", "equivalence", "value", "rank") & !is.na(n_finite) & n_finite < (B - 10), name]
GATES <- data.table(
  gate = c("G2", "G3", "G4", "P1", "P2", "P3"),
  what = c("b = 0 reproduces V4_perf_all and the S4 bootstrap SEs (<= 1e-12)", "bootstrap SE vs linearization SE within 7% (27 quantities)",
           "every inventory row classified; finite values in >= B - 10 replicates", "v3 numbers reproduced at the printed digit",
           "v4 - v3 within +/- 1 pp (63 values)", "block SE / independent SE > 1.25 (27 quantities)"),
  result = c(sprintf("%d of %d checks pass (max %.2g)", sum(g2$pass), nrow(g2), max(g2$max)),
             sprintf("%d of %d pass; rel. diff. %.3f to %.3f", sum(g3$pass), nrow(g3), min(g3$rel), max(g3$rel)),
             sprintf("missing %d; sub-claims with < B - 10 finite: %d", length(g4_missing), length(g4_nf)),
             sprintf("%d of %d reproduced, %d within one unit, %d not reproduced", TG[status == "reproduced", .N], nrow(TG), TG[status == "one unit of the last digit", .N], TG[status == "not reproduced", .N]),
             sprintf("mean diff %.3f; %d of %d confirmed", p2agg$mean_diff, p2agg$n_conf, p2agg$n),
             sprintf("ratio %.2f to %.2f; %d of %d above 1.25", min(p3$ratio), max(p3$ratio), sum(p3$confirmed), nrow(p3))),
  pass = c(all(g2$pass), all(g3$pass), length(g4_missing) == 0 && length(g4_nf) == 0, TG[status != "reproduced", .N] == 0,
           p2agg$n_conf == p2agg$n, all(p3$confirmed)))
print(GATES); fwrite(GATES, TO("V4_S5_gates.csv")); fwrite(g2, TO("V4_S5_g2.csv"))
cat(sprintf("elapsed %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
