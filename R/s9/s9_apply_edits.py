
import json,sys
E=json.load(open(sys.argv[1] if len(sys.argv)>1 else '/tmp/s9_edits.json')); log=[]
cache={}
for kind,f,a,b in E:
    s=cache.get(f) or open(f).read()
    if kind=='rep':
        n=s.count(a); assert n==1, (f,'rep count',n,a[:80]); s=s.replace(a,b)
    elif kind=='after':
        lines=s.split('\n'); idx=[i for i,l in enumerate(lines) if l.startswith(a)]
        assert len(idx)==1,(f,'after count',len(idx),a); lines.insert(idx[0]+1,b); s='\n'.join(lines)
    elif kind=='sub':
        n=s.count(a); assert n>=1,(f,'sub none'); s=s.replace(a,b); log.append(f'{f}: {n} x {a}->{b}')
    cache[f]=s
for f,s in cache.items(): open(f,'w').write(s)
print('\n'.join(log)); print('files edited:',len(cache))
