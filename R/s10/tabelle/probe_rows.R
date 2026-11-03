
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(2L)
R <- "/home/user/MMM/results_v4/"
for (f in c("variant_cc/S4/rows_v4.rds", "S6b/ns/rows_ns.rds", "S6cc/cc_all.rds")) {
  cat("=========", f, "\n"); t0 <- Sys.time(); A <- readRDS(file.path(R, f)); setDT(A)
  cat("read", format(Sys.time()-t0), "rows", nrow(A), "\n"); str(A, vec.len = 2)
  for (v in intersect(c("method","dgp","scen","label","term","ntime","N"), names(A))) { cat(v, ":"); print(table(A[[v]], useNA="ifany")) }
  if (all(c("method","converged") %in% names(A))) {
    print(A[, .(n=.N, conv=sum(converged %in% TRUE), convNA=sum(is.na(converged)), stNA=sum(is.na(status)),
      finite_ok=sum(is.finite(estimate)&is.finite(se)&abs(estimate)<10&se>0&se<10&is.finite(p)),
      notconv_but_ok=sum(!(converged %in% TRUE) & is.finite(estimate)&is.finite(se)&abs(estimate)<10&se>0&se<10&is.finite(p)),
      U0=if ("U0" %in% names(A)) sum(U0) else NA_integer_), by=.(method, term)][order(method)])
    if ("status" %in% names(A)) { cat("status of notconv_but_ok, top:\n"); print(head(A[!(converged %in% TRUE) & is.finite(estimate)&is.finite(se)&abs(estimate)<10&se>0&se<10&is.finite(p), .N, by=.(method, s=substr(status,1,70))][order(-N)], 15)) }
  }
  rm(A); invisible(gc())
}
x <- readRDS(file.path(R, "S7e/rows/copula_m3_N10_T10_manual_c01.rds")); cat("S7e cc example\n"); str(x$cc, vec.len=2); print(table(x$study$method))
