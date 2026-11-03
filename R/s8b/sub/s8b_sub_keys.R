#!/usr/bin/env Rscript
# =============================================================================
# s8b_sub_keys.R -- S8b sub-task: one keys file for the registry, tables/V5_S8b_sub_keys.csv (key, value, mcse, kind,
# source, how), the row-bind of the per-task key files written by
#   R/s8b/sub/s8b_table_s36.R (s36|...), R/s8b/sub/s8b_table_s37.R (s37|...), R/s8b/sub/s8b_mdkc_s26.R (mdkc|...),
#   py/s8b/s8b_mcse_captions.py (cap|...), py/s8b/s8b_figS1.py (figS1|...).
# Checks: every part present, unique keys, allowed kinds, numeric values for int/sig2/plain, int values integer.
# Run from ~/MMM last: Rscript --vanilla R/s8b/sub/s8b_sub_keys.R
# =============================================================================
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(1L); setwd("/home/user/MMM")
PARTS <- c(s36 = "tables/V5_S8b_sub_s36_keys.csv", s37 = "tables/V5_S8b_sub_s37_keys.csv", mdkc = "tables/V5_S8b_sub_mdkc_keys.csv",
           cap = "tables/V5_S8b_sub_cap_keys.csv", figS1 = "tables/V5_S8b_sub_figS1_keys.csv")
stopifnot(all(file.exists(PARTS)))
K <- rbindlist(lapply(names(PARTS), function(p) { x <- fread(PARTS[[p]], colClasses = "character", na.strings = NULL)
  stopifnot(identical(names(x), c("key", "value", "mcse", "kind", "source", "how")), all(startsWith(x$key, paste0(p, "|")))); x }))
stopifnot(!anyDuplicated(K$key), all(K$kind %in% c("int", "sig2", "plain", "text")))
v <- suppressWarnings(as.numeric(K$value)); isnum <- K$kind %in% c("int", "sig2", "plain")
stopifnot(all(is.finite(v[isnum])), all(abs(v[K$kind == "int"] - round(v[K$kind == "int"])) < 1e-9))
fwrite(K, "tables/V5_S8b_sub_keys.csv")
print(K[, .N, by = .(part = sub("\\|.*$", "", key), kind)])
cat("keys written:", nrow(K), "\n")
