"""Exact Hensel check at the retained discriminant-square local witnesses."""
from fractions import Fraction as F
from hashlib import sha256
from math import floor
from pathlib import Path
import json,resource,time

resource.setrlimit(resource.RLIMIT_CPU,(30,35))
ROOT=Path('/home/royvanrijn/src/jacobian-research/research')
SOURCE=ROOT/'artifacts/generated-results/elkies-k3-mestre-parent-branch-cancellation-v1/input.json'
OUT=ROOT/'artifacts/generated-results/elkies-k3-common-quartic-splitting-local-gate-v1'

def vp(n,p):
 if not n:return 1000000
 n=F(n);v=0;a,b=n.numerator,n.denominator
 while a%p==0:a//=p;v+=1
 while b%p==0:b//=p;v-=1
 return v

def ev(poly,t):return sum(F(c)*t**i for i,c in enumerate(poly))
def red(q,p,n):
 q=F(q);m=p**n
 return q.numerator*pow(q.denominator,-1,m)%m

def square(q,p):
 v=vp(q,p)
 if v==1000000 or v%2:return False
 unit=q/F(p)**v
 return red(unit,p,3)==1 if p==2 else pow(red(unit,p,1),(p-1)//2,p)==1

def root(a,b,p):
 scale=min(vp(a,p)//2,vp(b,p)//3)
 a,b=a/F(p)**(2*scale),b/F(p)**(3*scale)
 active=[(0,0)];visits=0
 while active and visits<20000:
  r,n=active.pop(0);visits+=1
  value=F(r)**3+a*r+b;der=3*F(r)**2+a
  if value==0 or vp(value,p)>2*vp(der,p):
   return {'status':'ROOT','root_scale':scale,'root_residue':r,'f_valuation':vp(value,p),'derivative_valuation':vp(der,p),'visits':visits}
  if n==12:return {'status':'DEFERRED','visits':visits}
  modulus=p**(n+1);aa,bb=red(a,p,n+1),red(b,p,n+1)
  for j in range(p):
   rr=r+j*p**n
   if (rr**3+aa*rr+bb)%modulus==0:active.append((rr,n+1))
 return {'status':'NO_ROOT' if not active else 'DEFERRED','visits':visits}

start=time.process_time();parents={p['name']:p for p in json.loads(SOURCE.read_text())['parents']}
preview=json.loads((OUT/'preview.json').read_text());rows=[]
for row in preview['rows']:
 parent=parents[row['parent']];p=row['prime'];coord=row['coordinate']
 t=F(coord) if row['chart']=='t' else F(1,coord) if coord else None
 a,b=(F(parent['A'][-1]),F(parent['B'][-1])) if t is None else (ev(parent['A'],t),ev(parent['B'],t))
 disc=-4*a**3-27*b**2
 assert square(disc,p),(row,disc)
 result=root(a,b,p)
 rows.append({'parent':row['parent'],'p':p,'chart':row['chart'],'coordinate':coord,**result})
 print(row['parent'],p,result['status'],result.get('root_residue'),flush=True)
result={'schema':'splitting-hensel-preflight-v1','preview_sha256':sha256((OUT/'preview.json').read_bytes()).hexdigest(),
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'rows':rows,'cpu_seconds':time.process_time()-start}
with (OUT/'hensel-preview.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
