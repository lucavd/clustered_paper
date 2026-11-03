# S3 Q3: the S3 replicate function, loaded with the ORIGINAL code (.rundir),
# reproduces results_full/ bit for bit (row keys, estimate, se, p, status) for
# replicates 1-5 of 15 generators x 2 arms x (N, T) in {(10, 2), (145, 10)}.
# Trace rows (E3 addition, status "[v4 trace] ...") are not in results_full by
# construction and are excluded before the comparison (counted separately).
source(path.expand("~/MMM/R/s3/s3_common.R"))
G <- expand.grid(dgp = DGPS, model = 1:5, label = c("mmlb", "manual"), NT = 1:2, stringsAsFactors = FALSE)
G$N <- c(10, 145)[G$NT]; G$ntime <- c(2, 10)[G$NT]
jobs <- G[rep(seq_len(nrow(G)), each = 5), ]; jobs$i <- rep(1:5, nrow(G))
res <- parallel::mclapply(seq_len(nrow(jobs)), function(k) {
  j <- jobs[k, ]; env <- load_env(j$dgp, RUNDIR_ORIG)
  f <- file.path(REPO, "results_full", sprintf("%s_m%d_N%d_T%d_%s_c01.rds", dgp_tag(j$dgp), j$model, j$N, j$ntime, j$label))
  st <- as.data.frame(readRDS(f)); st <- st[st$sim_id == j$i, ]
  o <- one_rep(env, sprintf("M%d", j$model), j$model, j$N, j$ntime, j$label, j$i, TECHS_ALL)
  rf <- as.data.frame(o$rows); ntr <- sum(grepl("^\\[v4 trace\\]", rf$status)); rf <- rf[!grepl("^\\[v4 trace\\]", rf$status) | is.na(rf$status), ]
  key <- function(x) paste(x$method, x$term)
  ok <- nrow(st) == nrow(rf) && identical(key(st), key(rf)) && identical(st$estimate, rf$estimate) &&
    identical(st$se, rf$se) && identical(st$p, rf$p) && identical(st$status, rf$status)
  data.frame(j[c("dgp", "model", "label", "N", "ntime", "i")], rows_stored = nrow(st), rows_bench = nrow(rf),
             trace_rows = ntr, identical = ok)
}, mc.cores = 40, mc.preschedule = FALSE)
bad <- !vapply(res, is.data.frame, TRUE); if (any(bad)) { print(res[bad]); stop("bench task errors") }
res <- do.call(rbind, res); print(res[!res$identical, ])
cat(sprintf("Q3 bench reproduces results_full: %d/%d (trace rows excluded: %d)\n", sum(res$identical), nrow(res), sum(res$trace_rows)))
write.csv(res, file.path(REPO, "tables/V4_S3_q3_bench.csv"), row.names = FALSE)
