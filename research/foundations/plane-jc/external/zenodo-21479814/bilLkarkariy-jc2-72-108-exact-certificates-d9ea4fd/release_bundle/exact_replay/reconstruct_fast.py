from pathlib import Path
import re,sys,pickle,math,time
from flint import fmpz,fmpq,fmpq_poly
R=Path('.');HN=pickle.load(open('hne0_polred.pkl','rb'));F=fmpq_poly([26,0,3,3,-1,1]);pat=re.compile(r'^C\|(\d+)\|(\d+)\|(\d+)\|(\d+)\|\((.*)\)$')
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
fs.sort();print('files',len(fs),flush=True);base=[k for k,a in fs[0][1]];xs=[fmpz(0)]*1925;M=fmpz(1);t0=time.time()
for z,(p,d) in enumerate(fs):
 if [k for k,a in d]!=base:raise SystemExit('support')
 inv=pow(int(M%p),-1,p);vals=[v for k,a in d for v in a]
 for j,r in enumerate(vals):xs[j]+=M*(((r-int(xs[j]%p))*inv)%p)
 M*=p
 if (z+1)%100==0:print('crt',z+1,'sec',round(time.time()-t0,1),flush=True)
digits=int(M.bit_length()/math.log2(10))+1;print('digits',digits,'crtsec',round(time.time()-t0,1),flush=True)
def rr(x,m):
 x%=m;B=(m//2).isqrt();r0,r1=m,x;t0,t1=fmpz(0),fmpz(1)
 while abs(r1)>B:q=r0//r1;r0,r1=r1,r0-q*r1;t0,t1=t1,t0-q*t1
 a,b=r1,t1
 if b<0:a,b=-a,-b
 if not b or abs(a)>B or b>B or a.gcd(b)!=1 or (x*b-a)%m:return None
 return int(a),int(b)
res=[]
for j,x in enumerate(xs):
 res.append(rr(x,M))
 if (j+1)%300==0:print('rr',j+1,flush=True)
miss=sum(z is None for z in res);print('reconstructed',1925-miss,'missing',miss,flush=True)
if miss:sys.exit(2)
Ts=[{} for _ in range(4)];pos=0
for key in base:
 i,a,b,c=key;co=[]
 for e in range(5):n,d=res[pos];pos+=1;co.append(fmpq(n,d))
 Ts[i-1][(a,b,c)]=fmpq_poly(co)%F
E=[{m:fmpq_poly([fmpq(*c.get(i,(0,1))) for i in range(5)])%F for m,c in e.items()} for e in HN];S={}
def add(m,c):
 c%=F
 if not c:return
 z=(S.get(m,fmpq_poly([0]))+c)%F
 if z:S[m]=z
 elif m in S:del S[m]
for i in range(4):
 for mt,ct in Ts[i].items():
  for me,ce in E[i].items():add(tuple(mt[k]+me[k] for k in range(3)),ct*ce)
add((1,0,0),fmpq_poly([-1]));print('residual',len(S));
if S:sys.exit(3)
print('EXACT_PASS')
with open('h_certificate_exact.txt','w') as fh:
 fh.write('minpoly=w^5-w^4+3*w^3+3*w^2+26\nidentity h=sum(Ti*Ei), i=1..4\n')
 for i,T in enumerate(Ts,1):
  for m,c in sorted(T.items()):
   ts=[]
   for e in range(len(c)):
    q=c[e]
    if q:ts.append(str(q)+('' if e==0 else ('*w' if e==1 else f'*w^{e}')))
   fh.write(f'C|{i}|{m[0]}|{m[1]}|{m[2]}|'+('+'.join(ts).replace('+-','-') or '0')+'\n')
