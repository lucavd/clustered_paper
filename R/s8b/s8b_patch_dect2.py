import json
p='R/s7f/s7f_h2h_boot.R'
s=open(p).read()
old="ND0 <- 1 - mats(\"gee_ar1_fg\", \"H0\", \"design_degenerate\"); ND1 <- 1 - mats(\"gee_ar1_fg\", \"H1\", \"design_degenerate\")"
assert s.count(old)==1
s=s.replace(old,"# DEC-T2 (S8b, D-74): the D-11 denominator (replicates with an estimable design) comes from the design records of every\n# replicate (RECS), not from the GEE rows: replicates with a degenerate design and no GEE row were counted as estimable.\nRU <- unique(RECS[, c(K6, \"design_degenerate\"), with = FALSE], by = K6)\nRU[, arm := fifelse(label == \"mmlb\", \"H1\", \"H0\")][, key := paste(dgp, scen, N, ntime)]\nndm <- function(arm_) { Y <- RU[arm == arm_ & key %in% KEYS$key]\n  stopifnot(nrow(Y) == 1000L * nrow(KEYS), Y[, .(n = .N), by = key][, all(n == 1000L)])\n  Mx <- matrix(NA_real_, nrow(KEYS), 1000L); Mx[cbind(match(Y$key, KEYS$key), Y$sim_id)] <- 1 - as.numeric(Y$design_degenerate)\n  stopifnot(!anyNA(Mx)); Mx }\nND0 <- ndm(\"H0\"); ND1 <- ndm(\"H1\")\ncat(sprintf(\"D-11 denominators from the design records: H0 %d, H1 %d estimable of %d; GEE-row denominators were H0 %d, H1 %d\\n\",\n            sum(ND0), sum(ND1), 1000L * nrow(KEYS), sum(1 - mats(\"gee_ar1_fg\", \"H0\", \"design_degenerate\")), sum(1 - mats(\"gee_ar1_fg\", \"H1\", \"design_degenerate\"))))")
open(p,'w').write(s); print('ok')
