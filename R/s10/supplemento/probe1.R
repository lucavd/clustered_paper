
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(1)
f <- list.files("~/MMM/results_v4/S7d/rows", full.names=TRUE); cat(length(f), basename(f[1:3]), "\n")
x <- readRDS(f[1]); str(x, max.level=2, list.len=40)
print(head(as.data.table(x$recs),3)); print(head(as.data.table(x$rows),3))
g <- list.files("~/MMM/results_v4/S7e/rows", full.names=TRUE); cat(length(g), basename(g[1:3]), "\n")
y <- readRDS(g[1]); str(y, max.level=1)
for (n in names(y)) {cat("==",n,"\n"); print(head(as.data.table(y[[n]]),2))}
