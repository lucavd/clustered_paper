
import re, os
os.makedirs("/tmp/gc1/tables", exist_ok=True); os.makedirs("/tmp/gc1/scripts", exist_ok=True)
spec = {"R/s7b/s7b_df.R": ('file.path(OUT_S7B, "rows_h0_marg.rds")', '"/tmp/gc1/rows_h0_marg.rds"'),
        "R/s7e/s7e_agg.R": ('file.path(OUT, "agg/boot_mats.rds")', '"/tmp/gc1/boot_mats.rds"'),
        "R/s7f/s7f_h2h_boot.R": ('file.path(OUT, "per_replicate_gee_ar1_fg_mmm_cc.rds")', '"/tmp/gc1/per_replicate_gee_ar1_fg_mmm_cc.rds"')}
for f, (a, b) in spec.items():
    s = open("/home/user/MMM/" + f).read()
    assert s.count(a) == 1, f
    s = s.replace(a, b)
    s = s.replace('source(path.expand("~/MMM/R/s11/s11_cv.R"))', 's11_drop_cv <- function(...) invisible(0L)   # GC1 check: D-90 off')
    s = re.sub(r'(fwrite\(.*?), file\.path\(TABDIR, ', r'\1, file.path("/tmp/gc1/tables", ', s)
    s = re.sub(r'(fwrite\([^\n]*?)"tables/', r'\1"/tmp/gc1/tables/', s)
    left = [l for l in s.splitlines() if ("fwrite(" in l and "/tmp/gc1" not in l) or ("saveRDS(" in l and "/tmp/gc1" not in l and "rows_h0" not in l)]
    print(f, "unredirected writes:", left)
    open("/tmp/gc1/scripts/" + os.path.basename(f), "w").write(s)
