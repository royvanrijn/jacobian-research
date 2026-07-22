from pathlib import Path
import re,pickle,time,math
from flint import fmpz
R=Path('.');pat=re.compile(r'^C\|(\d+)\|(\d+)\|(\d+)\|(\d+)\|\((.*)\)$')
def pp(s,p):
 s=s.replace(' ','');o=[0]*5
 for t in re.findall(r'[+-]?[^+-]+',s):
  if '*w' in t:a,r=t.split('*w',1);e=int(r[1:]) if r.startswith('^') else 1
  elif 'w' in t:
   m=re.fullmatch(r'([+-]?)w(?:\^(\d+))?',t);a='-1' if m.group(1)=='-' else '1';e=int(m.group(2) or 1)
  else:a=t;e=0
  o[e]=int(a)%p
 return o
def rf(f,p):
 d=[]
 for l in f.read_text(errors='ignore').splitlines():
  m=pat.match(l)
  if m:d.append((tuple(map(int,m.group(1,2,3,4))),pp(m.group(5),p)))
 return d
fs=[]
for f in R.glob('hcert_syz_p*.out'):
 t=f.read_text(errors='ignore')
 if 'CHECK_PASS' not in t:continue
 p=int(re.search(r'p(\d+)\.out$',f.name).group(1));d=rf(f,p)
 if len(d)==385:fs.append((p,d))
fs.sort();base=[k for k,a in fs[0][1]];xs=[fmpz(0)]*1925;M=fmpz(1);t=time.time()
for z,(p,d) in enumerate(fs):
 if [k for k,a in d]!=base:raise SystemExit('support')
 inv=pow(int(M%p),-1,p);vals=[v for k,a in d for v in a]
 for j,r in enumerate(vals):xs[j]+=M*(((r-int(xs[j]%p))*inv)%p)
 M*=p
pickle.dump((base,M,xs),open('crt_cache.pkl','wb'),protocol=4)
print('files',len(fs),'digits',int(M.bit_length()/math.log2(10))+1,'sec',time.time()-t,'bytes',Path('crt_cache.pkl').stat().st_size)
