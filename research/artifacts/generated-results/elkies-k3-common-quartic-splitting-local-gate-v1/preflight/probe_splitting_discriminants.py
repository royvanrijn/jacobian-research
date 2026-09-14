"""Bounded local discriminant-square preflight; no global point search."""
from fractions import Fraction
from hashlib import sha256
from math import comb, lcm
from pathlib import Path
import json
import resource
import time

resource.setrlimit(resource.RLIMIT_CPU,(30,35))
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
ROOT=Path('/home/royvanrijn/src/jacobian-research/research')
SOURCE=ROOT/'artifacts/generated-results/elkies-k3-mestre-parent-branch-cancellation-v1/input.json'
OUT=ROOT/'artifacts/generated-results/elkies-k3-common-quartic-splitting-local-gate-v1'
PRIMES=[2,3,5,7,11,13,17,19,23,29,31,37,41,43]

def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c

def vp(n,p):
 if not n:return 1000000
 v=0
 while n%p==0:n//=p;v+=1
 return v

def sub(a,r,m):
 return [sum(a[i]*comb(i,j)*r**(i-j)*m**j for i in range(j,len(a))) for j in range(len(a))]

def ev(a,t):return sum(c*t**i for i,c in enumerate(a))

def probe(H,p):
 # Every projective Qp base value: t integral or v=1/t in pZp.
 todo=[('t',H,0,1,0),('v',sub(H[::-1],0,p),0,p,0)]
 leaves=[];visits=0;deferred=[]
 while todo and visits<2000:
  chart,f,a,m,depth=todo.pop(0);visits+=1
  val=min(vp(c,p) for c in f); q=[c//p**val for c in f]
  modulus=8 if p==2 else p
  for r in range(modulus):
   b=ev(q,r)
   if b%p:
    square=(b%8==1) if p==2 else pow(b%p,(p-1)//2,p)==1
    if val%2==0 and square:
     z=a+m*r
     return {'status':'LOCAL_POINT','chart':chart,'coordinate':z,'valuation':val,
             'unit_mod':b%modulus,'modulus':modulus,'visits':visits}
    leaves.append({'chart':chart,'a':a+m*r,'m':m*modulus,'v':val,'unit_mod':b%modulus})
   elif depth==6:
    deferred.append({'chart':chart,'a':a+m*r,'m':m*modulus})
   else:
    todo.append((chart,sub(q,r,modulus),a+m*r,m*modulus,depth+1))
    # Removed val is only even when further zeros can matter. Track its parity.
    if val%2:todo[-1]=(chart,[p*c for c in todo[-1][1]],a+m*r,m*modulus,depth+1)
 if todo or deferred:
  return {'status':'DEFERRED','visits':visits,'queued':len(todo),'depth_leaves':len(deferred)}
 return {'status':'NO_LOCAL_POINT','visits':visits,'leaves':leaves}

start=time.process_time();pkt=json.loads(SOURCE.read_text())
rows=[]
for parent in pkt['parents']:
 A=list(map(Fraction,parent['A']));B=list(map(Fraction,parent['B']))
 a3=mul(A,mul(A,A));b2=mul(B,B)
 H=[-4*a3[i]-27*b2[i] for i in range(25)]
 denominator=lcm(*(c.denominator for c in H))
 H=[int(c*denominator**2) for c in H]
 common=0
 # A global square scale is irrelevant to the discriminant squareclass.
 for p in PRIMES:
  v=min(vp(c,p) for c in H);H=[c//p**(2*(v//2)) for c in H]
 for p in PRIMES:
  result=probe(H,p);row={'parent':parent['name'],'prime':p,**result};rows.append(row)
  print(parent['name'],p,result['status'],result.get('visits'),flush=True)
OUT.mkdir(parents=True,exist_ok=True)
result={'schema':'splitting-discriminant-preflight-v1',
        'source_sha256':sha256(SOURCE.read_bytes()).hexdigest(),
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope':'All four retained literal parents; fixed primes2 through43; projective Qp base cover with depth6 and2000 nodes per prime. Discriminant squareclasses only; no global point or low-degree-divisor conclusion from local solubility.',
        'cpu_seconds':time.process_time()-start,'rows':rows}
with (OUT/'preview.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
