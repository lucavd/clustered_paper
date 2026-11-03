# =============================================================================
# s2_b_rerun_missing.R -- S2(b), Q2.2/Q2.3: re-run, with their own seed and the
# production RNGkind (D-07), (i) every missing sim_id (generation only) and
# (ii) every generated sim_id that lacks one or more techniques (generation +
# the missing technique, without run_chunk's outer tryCatch), recording the
# error message. Controprova: the nearest present sim_id generates cleanly.
# =============================================================================
source(path.expand("~/MMM/R/s2/s2_common.R"))
A <- readRDS(file.path(OUTDIR, "rows_timetrt.rds"))
CELL <- c("dgp", "model", "N", "ntime", "label")
miss <- fread(file.path(TABDIR, "V4_S2_missing_ids.csv"))
A[, tech := fifelse(method %in% c("gee_md", "gee_kc", "gee_fg", "gee_ss"), "gee_ss", fifelse(method == "shrink", "glmer", method))]
TECHS <- c("mmm", "gee", "glmer", "gee_ss", "glmer_lrt")
have <- unique(A[, c(CELL, "sim_id", "tech"), with = FALSE])
part <- have[, .(missing_tech = paste(setdiff(TECHS, tech), collapse = ",")), by = c(CELL, "sim_id")][missing_tech != ""]
present <- unique(A[, c(CELL, "sim_id"), with = FALSE])
rm(A, have); invisible(gc())
cat(sprintf("missing ids %d | partial ids %d\n", nrow(miss), nrow(part)))
print(part[, .N, by = .(dgp, missing_tech)][order(-N)])

ANALYSIS_FILES <- c("analysis_model1.R", "analysis_model2.R", "analysis_model3.R", "analysis_model4.R",
                    "analysis_model5.R", "analysis_smallsample.R", "do_analysis.R")
worker_env <- function(dgp) {
  env <- load_family(dgp, RUNDIR)
  for (f in ANALYSIS_FILES) sys.source(file.path(RUNDIR, f), envir = env)
  env
}
ENVS <- setNames(lapply(DGPS, worker_env), DGPS)   # loaded before fork, inherited by workers
job <- function(dgp, model, N, ntime, label, sim_id, techs) {
  env <- ENVS[[dgp]]
  set_rep_seed(model, sim_id)
  d <- tryCatch(env$make_all_data(N, ntime, sim_id, model, label), error = function(e) e)
  if (inherits(d, "error"))
    return(data.table(stage = "generation", tech = NA_character_, ok = FALSE,
                      msg = conditionMessage(d), call = paste(deparse(conditionCall(d))[1], collapse = "")))
  if (!length(techs)) return(data.table(stage = "generation", tech = NA_character_, ok = TRUE, msg = NA_character_, call = NA_character_))
  rbindlist(lapply(techs, function(tc) {
    r <- tryCatch(env$do_analysis(d, N, ntime, label, model, tc, sim_id), error = function(e) e)
    if (inherits(r, "error")) data.table(stage = "analysis_outer", tech = tc, ok = FALSE, msg = conditionMessage(r),
                                         call = paste(deparse(conditionCall(r))[1], collapse = ""))
    else data.table(stage = "analysis_outer", tech = tc, ok = TRUE, msg = NA_character_, call = NA_character_)
  }))
}
run_set <- function(S, techs_col = NULL) {
  out <- parallel::mcMap(function(k) {
    r <- S[k]
    tc <- if (is.null(techs_col)) character(0) else strsplit(r[[techs_col]], ",")[[1]]
    cbind(r[, c(CELL, "sim_id"), with = FALSE], job(r$dgp, r$model, r$N, r$ntime, r$label, r$sim_id, tc))
  }, seq_len(nrow(S)), mc.cores = 40, mc.preschedule = TRUE)
  rbindlist(out, fill = TRUE)
}
t0 <- Sys.time()
M <- run_set(miss)
cat(sprintf("Q2.3 missing ids re-run: %d; generation errors reproduced: %d (%.1f%%)\n",
            nrow(M), sum(!M$ok), 100 * mean(!M$ok)))
print(M[, .N, by = .(ok, msg = substr(msg, 1, 80), call = substr(call, 1, 60))])
# controprova: nearest present neighbour of each missing id, same cell
nb <- miss[, .(sim_id = { s <- sim_id; cand <- setdiff(1:1000, s); vapply(s, function(x) cand[which.min(abs(cand - x))], 1L) }), by = CELL]
nb <- unique(nb)[present, on = c(CELL, "sim_id"), nomatch = 0]
NB <- run_set(nb)
cat(sprintf("Q2.3 controprova: neighbours %d, generated cleanly %d\n", nrow(NB), sum(NB$ok)))
P2 <- run_set(part, "missing_tech")
cat(sprintf("Q2.2 partial ids re-run: %d ids, %d technique runs; outer errors reproduced %d\n",
            nrow(part), sum(P2$stage == "analysis_outer"), sum(!P2$ok)))
print(P2[, .N, by = .(dgp, tech, ok, msg = substr(msg, 1, 90), call = substr(call, 1, 50))][order(-N)])
fwrite(M, file.path(TABDIR, "V4_S2_missing_rerun.csv"))
fwrite(NB, file.path(OUTDIR, "missing_neighbours_rerun.csv"))
fwrite(P2, file.path(TABDIR, "V4_S2_partial_rerun.csv"))
cat(sprintf("elapsed %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
