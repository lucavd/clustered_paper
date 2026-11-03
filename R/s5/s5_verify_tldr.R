# =============================================================================
# s5_verify_tldr.R -- S5: deterministic read-back of every number quoted in the TL;DR and in the
# decisions of docs/steps/S5.md (attested computation: file, md5, value, quoted text).
# Output: tables/V4_S5_tldr_receipts.csv ; stops if any quoted value differs from the file.
# =============================================================================
source(path.expand("~/MMM/R/s5/s5_common.R"))
f <- function(x) file.path(TABDIR, x)
md5 <- function(p) unname(tools::md5sum(p))
CL <- fread(f("V4_S5_claims.csv")); VAL <- fread(f("V4_S5_values.csv")); SW <- fread(f("V4_S5_switch_points.csv"))
FD <- fread(f("V4_S5_fu_d_wins.csv")); FT <- fread(f("V4_S5_fu_d_ties.csv")); P3 <- fread(f("V4_S5_p3_crn.csv")); P1 <- fread(f("V4_S5_p1_reproduction.csv"))
INV <- fread(file.path(REPO, "reports/morris/S5_claims_inventory.csv"))
v3log <- readLines(file.path(LOG_S5, "S5_v3_ranking_script.log"))
it <- function(x, d) formatC(x, format = "f", digits = d, decimal.mark = ",", big.mark = ".")
rng <- function(x, d) paste0(it(min(x), d), "-", it(max(x), d))
est <- function(nm, variant = "main") if (variant == "main") CL[name == nm, est] else VAL[variant == variant & name == nm, est]
R <- list(); add <- function(id, quoted, computed, file) R[[length(R) + 1]] <<- data.table(V = id, quoted = quoted, computed = computed, file = file, md5 = md5(file))
fc <- f("V4_S5_claims.csv"); fv <- f("V4_S5_values.csv")
# V-53 coverage of the inventory
add("V-53", "71", as.character(nrow(INV)), file.path(REPO, "reports/morris/S5_claims_inventory.csv"))
add("V-53", "64", as.character(CL[!is.na(est), uniqueN(claim_id)]), fc)
add("V-53", "149", as.character(CL[!is.na(est), .N]), fc)
add("V-53", "7", as.character(CL[kind %in% c("superseded", "out_of_scope", "not_operationalized"), .N]), fc)
add("V-53", "0", as.character(CL[class == "non risolta", .N]), fc)
add("V-53", "17", as.character(CL[class == "risolta, rilevanza indeterminata", .N]), fc)
# V-54 negligible differences
add("V-54", "0,25-0,42", rng(CL[claim_id == "C02", est], 2), fc)
add("V-54", "0,24-0,32", rng(CL[claim_id == "C20a", est], 2), fc)
add("V-54", "0,35-0,47", rng(CL[name %in% c("C15a|s1", "C15a|s2", "C15a|s3"), est], 2), fc)
add("V-54", "0,24-0,26", rng(CL[name %in% c("C01|NAT", "C01|AR1"), est], 2), fc)
add("V-54", "0,73", it(est("C03|GB"), 2), fc); add("V-54", "0,85", it(est("C03|NAT"), 2), fc); add("V-54", "1,03", it(est("C03|AR1"), 2), fc)
# V-55 relevant differences
add("V-55", "1,49", it(est("C01|GB"), 2), fc)
add("V-55", "3,4-3,7", rng(CL[claim_id == "C20b", est], 1), fc)
add("V-55", "1,2-1,5", rng(CL[claim_id == "C25b", est], 1), fc)
# V-56 conclusions changed by v4
add("V-56", "25", paste(unique(SW$switch_v4), collapse = ","), f("V4_S5_switch_points.csv"))
add("V-56", "35 e 45", paste(unique(SW$claimed), collapse = " e "), f("V4_S5_switch_points.csv"))
add("V-56", "85,8", it(VAL[variant == "main" & name == "use20|glmer_lrt", est], 1), fv)
add("V-56", "98,7", it(P1[name == "use20|glmer_lrt", printed], 1), f("V4_S5_p1_reproduction.csv"))
add("V-56", "-1,76", it(est("C11d|GB"), 2), fc)
add("V-56", "smentita x3", paste(unique(CL[claim_id == "C19a", class]), length(CL[claim_id == "C19a", class])), fc)
add("V-56", "323", it(FD[method == "mmm", v4_share], 0), f("V4_S5_fu_d_wins.csv")); add("V-56", "188", it(FD[method == "gee_fg", v4_share], 0), f("V4_S5_fu_d_wins.csv"))
add("V-56", "321", it(FD[method == "mmm", v3_first], 0), f("V4_S5_fu_d_wins.csv")); add("V-56", "308", it(FD[method == "gee_fg", v3_first], 0), f("V4_S5_fu_d_wins.csv"))
add("V-56", "138", as.character(FT[version == "v4" & who == "gee_fg+gee_md", N]), f("V4_S5_fu_d_ties.csv"))
add("V-56", "158", as.character(FT[version == "v3" & who == "gee_fg+gee_md", N]), f("V4_S5_fu_d_ties.csv"))
# V-57 common random numbers
add("V-57", "0,97-1,04", rng(P3$ratio, 2), f("V4_S5_p3_crn.csv"))
# V-58 reproduction of v3
add("V-58", "275", as.character(P1[status == "reproduced", .N]), f("V4_S5_p1_reproduction.csv"))
add("V-58", "288", as.character(nrow(P1)), f("V4_S5_p1_reproduction.csv"))
add("V-58", "12", as.character(P1[status != "reproduced" & grepl("^crank", name), .N]), f("V4_S5_p1_reproduction.csv"))
l <- grep("^[0-9]+: +gee_fg +[0-9.]+$", v3log, value = TRUE)[1]
add("V-58", "1,67", it(as.numeric(sub(".* ", "", l)), 2), file.path(LOG_S5, "S5_v3_ranking_script.log"))
# V-59 bobyqa
add("V-59", "6", as.character(CL[bobyqa_changes %in% TRUE, .N]), fc)
add("V-59", "AR1 trascurabile, NAT smentita", paste0("AR1 ", CL[name == "C16a|AR1", class_bobyqa], ", NAT ", CL[name == "C16a|NAT", class_bobyqa]), fc)
# V-60 M1c
add("V-60", "-1,2 a -1,7", paste0(it(max(VAL[variant == "m1c" & grepl("^C01\\|", name), est]), 1), " a ", it(min(VAL[variant == "m1c" & grepl("^C01\\|", name), est]), 1)), fv)
add("V-60", "1,8-1,9", rng(VAL[variant == "m1c" & grepl("^C03\\|", name), est], 1), fv)
# V-61 numbers of the decisions
nr <- CL[is.finite(n_req), n_req]
add("V-61", "40.000", it(round(max(nr[nr < 1e5]), -3), 0), fc)
add("V-61", "500.000", it(floor(max(nr) / 1e5) * 1e5, 0), fc)
add("V-61", "1,33", it(VAL[variant == "main" & name == "crank|gee_fg", est], 2), fv)
add("V-61", "2,56", it(VAL[variant == "main" & name == "crank|gee_md", est], 2), fv)
add("V-61", "0,12-0,14", rng(VAL[variant == "main" & grepl("^rateN\\|.*\\|gee_(fg|md)\\|", name), se], 2), fv)
X <- rbindlist(R)
X[, match := mapply(function(q, cc) {
  if (q == "smentita x3") return(cc == "smentita 3")
  if (q == "35 e 45") return(cc %in% c("35 e 45", "45 e 35"))
  if (q == "40.000") return(cc == "40.000")
  if (q == "500.000") return(cc == "500.000")
  q == cc }, quoted, computed)]
fwrite(X, f("V4_S5_tldr_receipts.csv"))
print(X[, .(V, quoted, computed, match)], nrows = 100)
if (!all(X$match)) stop("quoted values not matching the files: ", sum(!X$match))
cat("all", nrow(X), "quoted values match\n")
