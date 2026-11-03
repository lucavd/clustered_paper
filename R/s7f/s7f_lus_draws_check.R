
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setwd("/home/user/MMM")
d <- as.data.table(readRDS("ARTICOLO 2026/LURID_bis/lurid_long.rds"))
cat("patients", uniqueN(d$id), "\n")
inv <- d[, .(any_inv = any(as.character(invasive) == "1"), all_inv = all(as.character(invasive) == "1")), by = id]
K <- inv[, sum(any_inv)]; P <- nrow(inv); cat("patients with any invasive row:", K, " with all rows invasive:", inv[, sum(all_inv)], "\n")
for (n in c(10, 15, 20)) cat(sprintf("N=%d expected P(no invasive patient) = %.4f\n", n, dhyper(0, K, P - K, n)))
dr <- readRDS("rev2/T3_draws.rds")
S <- rbindlist(lapply(dr, function(z) data.table(N = z$N, rep = z$rep, n_inv = sum(z$sel %in% inv[any_inv == TRUE, id]))))
print(S[, .(draws = .N, share_none = mean(n_inv == 0), mean_inv = mean(n_inv), expected_mean = N[1] * K / P), by = N])
x <- S[N == 10, sum(n_inv == 0)]; cat("binomial test N=10:", x, "/", S[N == 10, .N], " p =", binom.test(x, S[N == 10, .N], dhyper(0, K, P - K, 10))$p.value, "\n")
cat("class of id:", class(d$id), " ids in draws found:", mean(unlist(lapply(dr, `[[`, "sel")) %in% d$id), "\n")
