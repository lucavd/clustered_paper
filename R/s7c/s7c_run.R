# s7c_run.R [workers] -- pilot run, resumable (skips chunks already written)
source(path.expand("~/MMM/R/s7c/s7c_common.R"))
a <- commandArgs(TRUE); W <- if (length(a)) as.integer(a[1]) else 20L
tk <- pilot_tasks(); todo <- tk[!file.exists(file)]
cat(format(Sys.time()), "tasks", nrow(tk), "todo", nrow(todo), "workers", W, "\n")
t0 <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(todo)), function(k) tryCatch(run_chunk_s7c(todo[k]), error = function(e) paste("ERROR", conditionMessage(e))),
                          mc.cores = W, mc.preschedule = FALSE)
err <- Filter(is.character, res)
cat(format(Sys.time()), "done; errors", length(err), "; wall", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
if (length(err)) print(head(unlist(err)))
cat("files present", sum(file.exists(tk$file)), "of", nrow(tk), "\n")
