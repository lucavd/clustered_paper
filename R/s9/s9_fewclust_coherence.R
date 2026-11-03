# s9_fewclust_coherence.R -- S9 (PI decision, 2026-11-02): fewclust is cited as implementation of the recommendation,
# with one sentence on its extra coherence condition. Reads the check-3 diagnosis of fewclust 0.3.1 from the committed
# package branch (git show, no copy by hand) and writes tables/V6_S9_fewclust_coherence.csv.
suppressMessages(library(data.table)); setwd(path.expand("~/MMM"))
C <- "626393bd5703b8f459a5b13655b843dbee78855d"
x <- system2("git", c("-C", "/home/user/MMM-package", "show", paste0(C, ":fewclust/dev/checks_0.3.1/out/c9b_summary.txt")), stdout = TRUE)
g <- function(pat) { l <- grep(pat, x, value = TRUE)[1]; stopifnot(!is.na(l)); l }
num <- function(l, re) as.integer(sub(re, "\\1", l))
l1 <- g("^compared:"); l3 <- g("^mismatches among datasets usable"); l5 <- g("^description, default control")
X <- data.table(
  item = c("n_cmp", "n_mismatch", "n_usable_study", "n_mismatch_usable", "n_coh"),
  value = c(num(l1, "^compared: ([0-9]+).*"), num(l1, ".*mismatches: ([0-9]+).*"),
            num(l3, ".*: [0-9]+ of ([0-9]+).*"), num(l3, ".*: ([0-9]+) of [0-9]+.*"),
            num(l5, ".*reported fit different in ([0-9]+) of.*")),
  source = paste0("MMM-package@", substr(C, 1, 7), ":fewclust/dev/checks_0.3.1/out/c9b_summary.txt"))
stopifnot(num(l5, ".*different in [0-9]+ of ([0-9]+) .*") == X[item == "n_cmp", value])
fwrite(X, "tables/V6_S9_fewclust_coherence.csv"); print(X)
