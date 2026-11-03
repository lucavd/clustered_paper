# =============================================================================
# s2_b_check_partial.R -- S2(b), Q2.2 follow-up: for the generated sim_ids that
# lack a time:trt (or failure) row for some technique, read the raw chunk rows
# for ALL terms; check the overlap with the GLMM-AR1 missing ids; recount the
# manuscript's "attempted" definition (replicates with >= 1 time:trt row).
# =============================================================================
source(path.expand("~/MMM/R/s2/s2_common.R"))
A <- readRDS(file.path(OUTDIR, "rows_timetrt.rds"))
CELL <- c("dgp", "model", "N", "ntime", "label")
miss <- fread(file.path(TABDIR, "V4_S2_missing_ids.csv"))
A[, tech := fifelse(method %in% c("gee_md", "gee_kc", "gee_fg", "gee_ss"), "gee_ss", fifelse(method == "shrink", "glmer", method))]
TECHS <- c("mmm", "gee", "glmer", "gee_ss", "glmer_lrt")
have <- unique(A[, c(CELL, "sim_id", "tech"), with = FALSE])
part <- have[, .(missing_tech = paste(setdiff(TECHS, tech), collapse = ",")), by = c(CELL, "sim_id")][missing_tech != ""]
# (1) raw rows of those replicates, all terms
part[, chunk := (sim_id - 1L) %/% 100L + 1L]
pref <- c("GenBinaryY" = "GenBinaryY", "GLMM-native" = "GLMMnative", "GLMM-AR1" = "GLMMAR1")
part[, file := sprintf("%s/results_full/%s_m%d_N%d_T%d_%s_c%02d.rds", REPO, pref[dgp], model, N, ntime, label, chunk)]
raw <- rbindlist(lapply(unique(part$file), function(f) {
  x <- as.data.table(readRDS(f)); k <- part[file == f]
  x[sim_id %in% k$sim_id, .(file = f, sim_id = as.integer(sim_id), method, term, converged, status)]
}))
raw <- merge(raw, part[, .(file, sim_id, dgp, model, N, ntime, label, missing_tech)], by = c("file", "sim_id"))
raw[, tech := fifelse(method %in% c("gee_md", "gee_kc", "gee_fg", "gee_ss"), "gee_ss", fifelse(method == "shrink", "glmer", method))]
pres <- raw[, .(techs_any_term = paste(sort(unique(tech)), collapse = ",")), by = .(file, sim_id)]
cat(sprintf("Q2.2 partial ids %d; with all 5 techniques present on SOME term: %d\n", nrow(part),
            sum(vapply(strsplit(pres$techs_any_term, ","), function(z) all(TECHS %in% z), TRUE))))
cat("terms returned by the 'missing' technique:\n")
print(raw[mapply(grepl, tech, missing_tech), .N, by = .(dgp, method, terms = term)][order(dgp, method, -N)][, head(.SD, 8), by = .(dgp, method)])
# (2) overlap with the GLMM-AR1 missing ids (same model, N, T, arm, sim_id)
ov <- merge(unique(part[, .(model, N, ntime, label, sim_id, dgp)]), miss[, .(model, N, ntime, label, sim_id, ar1_missing = TRUE)],
            by = c("model", "N", "ntime", "label", "sim_id"), all.x = TRUE)
cat(sprintf("Q2.2 partial ids that are also missing in GLMM-AR1: %d of %d\n", sum(!is.na(ov$ar1_missing)), nrow(ov)))
# per AR1-missing id: which methods returned a time:trt row in GenBinaryY / GLMM-native?
tt <- A[term == "time:trt", .(n_tt = uniqueN(method)), by = c(CELL, "sim_id")]
z <- merge(CJ_dt <- miss[, .(model, N, ntime, label, sim_id)][, .(dgp = c("GenBinaryY", "GLMM-native")), by = .(model, N, ntime, label, sim_id)],
           tt, by = c(CELL, "sim_id"), all.x = TRUE)
z[is.na(n_tt), n_tt := 0L]
cat("AR1-missing ids, in the other mechanisms: number of methods with a time:trt row\n"); print(z[, .N, by = .(dgp, n_tt)][order(dgp, n_tt)])
# (3) manuscript definition: attempted = replicates with at least one time:trt row
att <- A[term == "time:trt", .(n_att = uniqueN(sim_id)), by = CELL]
dp <- att[, .(m = min(n_att)), by = .(dgp, model, N, ntime)]
cat(sprintf("manuscript definition: design points with < 1000 'attempted' in at least one arm: %d (by dgp: %s); min %d\n",
            sum(dp$m < 1000), paste(dp[m < 1000, .N, by = dgp][, paste(dgp, N)], collapse = "; "), min(dp$m)))
# are the same (model, N, T, sim_id) involved across mechanisms?
lost <- A[, .(any_tt = any(term == "time:trt", na.rm = TRUE)), by = c(CELL, "sim_id")][any_tt == FALSE]
lost <- rbind(lost[, .(dgp, model, N, ntime, label, sim_id)], miss[, .(dgp, model, N, ntime, label, sim_id)])
w <- dcast(lost[, .(cnt = .N), by = .(model, N, ntime, label, sim_id, dgp)], model + N + ntime + label + sim_id ~ dgp, value.var = "cnt", fill = 0)
cat(sprintf("replicates lost under the manuscript definition: %d (model,N,T,arm,sim_id) keys; lost in all 3 mechanisms: %d\n",
            nrow(w), sum(rowSums(w[, .(`GLMM-AR1`, `GLMM-native`, GenBinaryY)]) == 3)))
fwrite(w, file.path(TABDIR, "V4_S2_lost_replicates_by_mechanism.csv"))
fwrite(pres, file.path(OUTDIR, "partial_ids_techniques_any_term.csv"))
