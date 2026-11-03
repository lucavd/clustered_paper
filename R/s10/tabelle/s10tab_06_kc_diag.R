
source("/home/user/MMM/R/s10/tabelle/s10tab_common.R"); setDTthreads(2L)
KC <- readRDS(file.path(OUT, "s1_kc_rep.rds"))[label == "manual" & dgp == "GLMM-AR1"]
V <- readRDS(file.path(REPO, "results_v4/variant_nc/S4/rows_v4.rds"))
V <- V[method %in% c("gee_kc", "gee") & dgp == "GLMM-AR1" & label == "manual" & scen %in% c("M1","M2","M3","M4","M5")]
V[, scen := norm_scen(scen)]
G <- V[method == "gee_kc" & term %in% "time:trt", .(dgp, scen, N, ntime, label, sim_id, estimate, se, p, U0)]
X <- merge(KC, G, by = K6, all = TRUE)
X[, rej_proj := U0 %in% TRUE & p < 0.05]
print(X[, .(n = .N, own_use = sum(use %in% TRUE), proj_use = sum(U0 %in% TRUE), use_diff = sum(xor(use %in% TRUE, U0 %in% TRUE)), rej_diff = sum(xor(rn %in% TRUE, rej_proj)), proj_row_missing = sum(is.na(U0)))])
d <- X[xor(use %in% TRUE, U0 %in% TRUE) | xor(rn %in% TRUE, rej_proj)]
print(d[, .(scen, N, ntime, sim_id, use, rn, estimate, se, p, U0)])
GE <- V[method == "gee" & term %in% "time:trt", .(dgp, scen, N, ntime, label, sim_id, gee_conv = converged, gee_se = se, gee_status = status)]
print(merge(d[, ..K6], GE, by = K6, all.x = TRUE))
