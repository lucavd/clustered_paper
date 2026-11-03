
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(1L)
x <- readRDS("/home/user/MMM/results_v4/S6b/ns/rows/GenBinaryY_m6_N10_T2_manual_c01.rds")
cat("NS class", class(x), "\n"); if (is.list(x) && !is.data.frame(x)) { print(names(x)); for (n in names(x)) {cat("==",n,"\n"); str(x[[n]], vec.len=1, max.level=1)} } else str(x, vec.len=1)
f <- list.files("/home/user/MMM/results_v4/S7d2/rows")[1:3]; print(f); print(length(list.files("/home/user/MMM/results_v4/S7d2/rows")))
y <- readRDS(file.path("/home/user/MMM/results_v4/S7d2/rows", f[1])); cat("S7d2 class", class(y), "\n"); if (is.list(y) && !is.data.frame(y)) { print(names(y)); for (n in names(y)) {cat("==",n,"\n"); str(y[[n]], vec.len=1, max.level=1)} } else str(y, vec.len=1)
k <- readRDS("/home/user/MMM/results_v4/S9/const_keys.rds"); str(k, vec.len=1); print(nrow(k))
