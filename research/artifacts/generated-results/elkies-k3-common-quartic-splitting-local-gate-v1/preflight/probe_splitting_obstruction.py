"""Resolve projective local base cells for complete cubic splitting, p>=5."""
from fractions import Fraction as F
from hashlib import sha256
from math import comb,lcm
from pathlib import Path
import json,resource,time

resource.setrlimit(resource.RLIMIT_CPU,(30,35))
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
ROOT=Path('/home/royvanrijn/src/jacobian-research/research')
SOURCE=ROOT/'artifacts/generated-results/elkies-k3-mestre-parent-branch-cancellation-v1/input.json'
OUT=ROOT/'artifacts/generated-results/elkies-k3-common-quartic-splitting-local-gate-v1'

def vp(n,p):
 if not n:return 1000000
 n=F(n);v=0;a,b=n.numerator,n.denominator
 while a%p==0:a//=p;v+=1
 while b%p==0:b//=p;v-=1
 return v
def red(n,p):
 n=F(n);return n.numerator*pow(n.denominator,-1,p)%p
def ev(poly,t):return sum(c*t**i for i,c in enumerate(poly))
def sub(a,r,m):
 return [sum(a[i]*comb(i,j)*r**(i-j)*m**j for i in range(j,len(a))) for j in range(len(a))]
def mul(a,b):
 c=[F(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c
def discriminant(A,B):
 a3,b2=mul(A,mul(A,A)),mul(B,B)
 return [-4*a3[i]-27*b2[i] for i in range(25)]

def square_on_cell(H,p):
 den=lcm(*(c.denominator for c in H));H=[int(c*den**2) for c in H]
 todo=[(H,0,1,0)];leaves=[];visits=0;deferred=[]
 while todo and visits<2000:
  f,a,m,depth=todo.pop(0);visits+=1
  v=min(vp(c,p) for c in f);q=[c//p**v for c in f]
  for r in range(p):
   b=ev(q,r)%p
   if b:
    if v%2==0 and pow(b,(p-1)//2,p)==1:return {'status':'SQUARE','coordinate':a+m*r,'visits':visits}
    leaves.append({'a':a+m*r,'m':m*p,'v':v,'unit':b})
   elif depth==7:deferred.append((a+m*r,m*p))
   else:todo.append(([p**(v%2)*c for c in sub(q,r,p)],a+m*r,m*p,depth+1))
 return {'status':'NO_SQUARE' if not todo and not deferred else 'DEFERRED','visits':visits,'leaves':leaves}

def probe(A,B,p):
 todo=[('t',A,B,0,1,0),('v',sub(A[::-1],0,p),sub(B[::-1],0,p),0,p,0)]
 leaves=[];visits=0;deferred=[]
 while todo and visits<2000:
  chart,A,B,a,m,depth=todo.pop(0);visits+=1
  va,vb=min(vp(c,p) for c in A),min(vp(c,p) for c in B)
  scale=min(va//2,vb//3);A=[c/F(p)**(2*scale) for c in A];B=[c/F(p)**(3*scale) for c in B]
  va-=2*scale;vb-=3*scale
  for r in range(p):
   aa,bb=red(ev(A,r),p),red(ev(B,r),p)
   leaf={'chart':chart,'a':a+m*r,'m':m*p,'depth':depth,'scale':scale,'A_residue':aa,'B_residue':bb}
   if aa==bb==0:
    # Full splitting of z^3+a*z+b with triple residue root0 requires all
    # three roots in pZp, hence a in p^2 and b in p^3.
    fail=(va<2 and red(ev(A,r)/F(p)**va,p)!=0) or (vb<3 and red(ev(B,r)/F(p)**vb,p)!=0)
    if fail:leaves.append({**leaf,'reason':'triple-root coefficient valuations'});continue
    if depth==7:deferred.append(leaf)
    else:todo.append((chart,sub(A,r,p),sub(B,r,p),a+m*r,m*p,depth+1))
    continue
   delta=(-4*aa**3-27*bb**2)%p
   roots=[x for x in range(p) if (x**3+aa*x+bb)%p==0]
   if delta:
    if len(roots)==3:return {'status':'SPLITS','chart':chart,'coordinate':a+m*r,'roots_mod_p':roots,'visits':visits}
    leaves.append({**leaf,'reason':'residue cubic does not split','roots':roots});continue
   # Double plus simple residue root. The simple root lifts; the quadratic
   # complement splits precisely when the nonzero cubic discriminant is square.
   assert len(roots)==2
   H=discriminant(sub(A,r,p),sub(B,r,p));sq=square_on_cell(H,p)
   if sq['status']=='SQUARE':return {'status':'SPLITS','chart':chart,'coordinate':a+m*(r+p*sq['coordinate']),'via':'simple root and discriminant square','visits':visits}
   if sq['status']=='DEFERRED':deferred.append({**leaf,'reason':'discriminant square unresolved'})
   else:leaves.append({**leaf,'reason':'discriminant nonsquare','square_exclusion':sq})
 return {'status':'NO_SPLITTING' if not todo and not deferred else 'DEFERRED','visits':visits,'leaves':leaves,'deferred':deferred}

start=time.process_time();pkt=json.loads(SOURCE.read_text());rows=[]
for c in pkt['parents']:
 for p in [7,17,19,29]:
  result=probe(list(map(F,c['A'])),list(map(F,c['B'])),p)
  rows.append({'parent':c['name'],'p':p,**result})
  print(c['name'],p,result['status'],'visits',result['visits'],'coordinate',result.get('coordinate'),flush=True)
result={'schema':'full-splitting-cell-preflight-v1','source_sha256':sha256(SOURCE.read_bytes()).hexdigest(),
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'rows':rows,'cpu_seconds':time.process_time()-start,
        'scope':'Four fixed parents; primes7,17,19,29; complete projective Qp base cells, depth7 and2000 nodes for each stage. Any DEFERRED remains unresolved. No global rational-point or coefficient search.'}
with (OUT/'cells-preview.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
