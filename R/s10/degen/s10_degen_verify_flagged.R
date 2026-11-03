# s10_degen_verify_flagged.R -- S10 "Degeneri": independent recount at the two design points where GEE AR(1) FG is clearly
# above 7.5% (M2, N = 10, T = 2, GLMM-AR1 and GLMM-native), from the indicator files and classes.rds; and a divergence test:
# each usable Si dataset of GEE AR(1) FG there is regenerated and refitted with geese maxit 25 (default, as the study) and 200.
source(path.expand("~/MMM/R/s10/degen/s10_degen_common.R")); setDTthreads(2L)
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
CL <- readRDS(file.path(OUT_DG, "classes.rds"))[scen == "M2" & N == 10L & ntime == 2L & dgp %in% c("GLMM-AR1", "GLMM-native")]
out <- list()
for (s in c("gee_ar1_fg", "gee_ind_fg", "gee_exch_fg", "mmm_cc", "glmer_lrt", "flip")) {
  y <- readRDS(file.path(OUT_DG, sprintf("ind/T5_%s.rds", s)))[scen == "M2" & N == 10L & ntime == 2L & dgp %in% c("GLMM-AR1", "GLMM-native")]
  y <- merge(y, CL[, c(K6, "S", "Si", "elig"), with = FALSE], by = K6, all.x = TRUE)
  out[[s]] <- y[, .(strategy = s, usable = .N, rej = sum(rej), usable_Si = sum(Si %in% TRUE), rej_Si = sum(rej & Si %in% TRUE), usable_notS = sum(!(S %in% TRUE)), rej_notS = sum(rej & !(S %in% TRUE))), by = .(dgp, label)]
  if (s == "gee_ar1_fg") SI <- y[Si %in% TRUE & label == "manual"]
}
R <- rbindlist(out); print(R); fwrite(R, file.path(OUT_DG, "verify_flagged_points.csv"))
cat("eligible at the two points:", CL[, sum(elig)], " S:", CL[, sum(S, na.rm = TRUE)], " Si:", CL[, sum(Si, na.rm = TRUE)], "\n")
TK <- dg_tasks()
fit <- function(d, cs, mx) { f <- tryCatch(geepack::geeglm(binY ~ time * trt + age, id = id, data = d, family = binomial, corstr = cs, scale.fix = TRUE,
                                                          control = geepack::geese.control(maxit = mx)), error = function(e) NULL)
  if (is.null(f)) return(c(NA, NA)); c(unname(coef(f)["time:trt"]), sqrt(diag(vcov(f)))[["time:trt"]]) }
DV <- rbindlist(lapply(seq_len(nrow(SI)), function(r) { k <- SI[r]
  t <- TK[src == "S7d" & dgp == k$dgp & scen == "M2" & N == 10L & ntime == 2L & label == "manual" & chunk == (k$sim_id - 1L) %/% 100L + 1L]
  gen <- dg_gen(t); set_rep_seed(t$model, k$sim_id); d <- as.data.frame(gen(t$N, t$ntime, k$sim_id, t$label)); d <- d[order(d$id, d$time), ]
  a <- fit(d, "ar1", 25L); b <- fit(d, "ar1", 200L); g <- fit(d, "independence", 25L); h <- fit(d, "independence", 200L)
  data.table(dgp = k$dgp, sim_id = k$sim_id, est_stored = k$est, se_stored = k$se, rej = k$rej, est25 = a[1], se25 = a[2], est200 = b[1], se200 = b[2],
             ind25 = g[1], ind200 = h[1]) }))
DV[, `:=`(same_as_stored = abs(est25 - est_stored) < 1e-8, moved_out = abs(est200) > abs(est25) + 0.5)]
print(DV[, .(n = .N, rej = sum(rej), same_as_stored = sum(same_as_stored, na.rm = TRUE), est200_na = sum(is.na(est200)), moved_out = sum(moved_out, na.rm = TRUE),
             med_abs25 = median(abs(est25), na.rm = TRUE), med_abs200 = median(abs(est200), na.rm = TRUE), med_se25 = median(se25, na.rm = TRUE), med_se200 = median(se200, na.rm = TRUE),
             med_abs_ind25 = median(abs(ind25), na.rm = TRUE), med_abs_ind200 = median(abs(ind200), na.rm = TRUE)), by = dgp])
fwrite(DV, file.path(OUT_DG, "verify_flagged_divergence.csv"))
