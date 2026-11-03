# s8b_lus_lateral.R -- S8b, D-70: lateral-vs-anterior contrast of the marginalized model on the full cohort, with the latent
# variable (lv.formula = ~1, q = 20), the specification of R/s8b/s8b_clin_mmm.R; replaces the value without latent variable
# of tables/V4_S2_lus_vs_dgm.csv (item contrast_MMM_regionlateral) in Section 2.2. Output: tables/V5_S8b_lus_lateral.csv.
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({ library(data.table); library(MMLB); library(dplyr) }); setwd(path.expand("~/MMM"))
d0 <- readRDS("ARTICOLO 2026/LURID_bis/lurid_long.rds")
d0$row_in_id <- ave(seq_len(nrow(d0)), d0$id, FUN = seq_along); d0$tp <- (d0$row_in_id - 1L) %/% 10L + 1L
dL <- d0 %>% filter(!is.na(lus_severe)) %>%
  mutate(region = factor(case_when(grepl("A$", as.character(section)) ~ "anterior", grepl("P$", as.character(section)) ~ "posterior", TRUE ~ "lateral"),
                         levels = c("anterior", "lateral", "posterior")), id = droplevels(factor(id)))
D <- as.data.frame(dL %>% arrange(id) %>% mutate(id = droplevels(factor(id)), across(where(is.factor), droplevels)))
m <- mm(lus_severe ~ region + disease + invasive + age_weeks, lv.formula = ~1, t.formula = ~1, data = D, id = id, q = 20, verbose = FALSE)
O <- data.table(term = c("regionlateral", "regionposterior"), estimate = unname(m$beta[c("regionlateral", "regionposterior")]),
                se_mod = sqrt(diag(m$mod.cov)[c("regionlateral", "regionposterior")]), code = unname(m$control["convergence_code"]))
O <- O[term == "regionlateral"]; print(O); fwrite(O, "tables/V5_S8b_lus_lateral.csv")
