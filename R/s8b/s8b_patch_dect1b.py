p='R/s5/s5_build.R'
s=open(p).read()
old1="V[, c(\"estimate\", \"se\", \"th\") := NULL]"
assert s.count(old1)==1
s=s.replace(old1,"V[, za := abs(estimate / se)]   # DEC-T1 (S8b, D-73): |estimate/se| kept for the size-adjusted rule below\nV[, c(\"estimate\", \"se\", \"th\") := NULL]")
old2="p = fifelse(method == \"glmer_lrt\", p, -abs(estimate / se)), u0)]"
assert s.count(old2)==1
s=s.replace(old2,"p = fifelse(method == \"glmer_lrt\", p, -za), u0)]")
open(p,'w').write(s); print('ok')
