
suppressMessages(library(data.table))
N <- fread("tables/V4_S7e_dp_rates.csv")[ref == "t" & strategy == "glmer", .(dgp, scen, N, ntime, rate)]
O <- fread("tables/V4_S7b_rule_dp.csv")[ref == "t" & strategy == "glmer|t" & scen %in% sprintf("M%d", 1:5), .(dgp, scen, N, ntime, rate)]
A <- rbind(N, O)[, .(rate = mean(rate, na.rm = TRUE)), by = .(dgp, scen, ntime)]
W <- dcast(A, scen + ntime ~ dgp, value.var = "rate"); print(W, digits = 3); fwrite(W, "tables/V4_S7e_glmm_by_scen_T.csv")
