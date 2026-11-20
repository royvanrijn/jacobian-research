#!/usr/bin/env python3
"""Small exact controls for the local coalescing-binomial bridge only."""
from __future__ import annotations
from fractions import Fraction as F
from math import factorial,comb
from pathlib import Path
import argparse,hashlib,json
if not __debug__:raise SystemExit('Refusing optimized execution: exact assertions required.')
HERE=Path(__file__).resolve().parent
class Quad:
 def __init__(self,a,d):self.a=F(a);self.d=d;self.D=F(a+d);self.t=self.D+self.a-1;self.p=self.a*self.D
 def add(self,x,y):return(x[0]+y[0],x[1]+y[1])
 def neg(self,x):return(-x[0],-x[1])
 def sub(self,x,y):return self.add(x,self.neg(y))
 def mul(self,x,y):return(x[0]*y[0]-x[1]*y[1]*self.p,x[0]*y[1]+x[1]*y[0]+x[1]*y[1]*self.t)
 def scale(self,x,c):return(x[0]*c,x[1]*c)
 def conj(self,x):return(x[0]+self.t*x[1],-x[1])
 def norm(self,x):return x[0]*x[0]+self.t*x[0]*x[1]+self.p*x[1]*x[1]
 def trace(self,x):return 2*x[0]+self.t*x[1]
 def inv(self,x):
  n=self.norm(x);assert n;return self.scale(self.conj(x),1/n)
 def div(self,x,y):return self.mul(x,self.inv(y))
 def power(self,x,n):
  if n<0:return self.power(self.inv(x),-n)
  out=(F(1),F(0))
  while n:
   if n&1:out=self.mul(out,x)
   x=self.mul(x,x);n//=2
  return out
ONE=(F(1),F(0));ZERO=(F(0),F(0));S=(F(0),F(1))
def falling(n,k):
 p=1
 for j in range(k):p*=n-j
 return p

def polymul(R,A,B):
 out=[ZERO]*(len(A)+len(B)-1)
 for i,u in enumerate(A):
  for j,v in enumerate(B):out[i+j]=R.add(out[i+j],R.mul(u,v))
 return out

def actual_shift(a,d,N):
 n=N-1;D=a+d
 J=[F(comb(n,k)*factorial(a*n+d*k),N**(a*n+d*k+1))for k in range(n+1)]
 # Coefficient index k corresponds to w^(n-k), before and after substitution.
 I=[F(comb(n,k)*factorial(a*n+d*k),n**(a*n+d*k+1))for k in range(n+1)]
 rhs=[F(n,N)**(D*n+1)*c*F(N,n)**(d*(n-k))for k,c in enumerate(I)]
 assert J==rhs
 return J


MANDATORY_SOURCES = {'airy_proof': ('packet', 'COALESCING_BINOMIAL_SADDLE_BRIDGE.md'), 'bulk_proof': ('packet', 'COALESCING_BINOMIAL_HERMITE_BULK.md'), 'airy_checker': ('packet', 'verify_coalescing_binomial_saddle_bridge.py'), 'bulk_checker': ('packet', 'verify_coalescing_binomial_hermite_bulk.py'), 'large_lower': ('repository', 'research/factorial/image-rigidity/GAMMA_PRODUCT_ROOT_MODULUS.tex'), 'lower_four': ('repository', 'research/factorial/image-rigidity/LOWER_FOUR_QUADRATIC_SECTOR_SFC.tex'), 'lower_five': ('repository', 'research/factorial/image-rigidity/LOWER_FIVE_QUADRATIC_SECTOR_SFC.tex'), 'gamma_obstruction': ('repository', 'research/factorial/image-rigidity/GAMMA_CONVOLUTION_APPELL_OBSTRUCTION.tex'), 'etna_context': ('repository', 'research/factorial/image-rigidity/sources/dunster_gil_ruiz_segura_etna65_2026_140_155.pdf'), 'wolfs_context': ('repository', 'research/factorial/image-rigidity/sources/wolfs-2024-multiple-orthogonal-2401.08312v1.pdf')}

def verify_sources(repo_arg):
    repo=repo_arg.resolve() if repo_arg else next((p for p in (HERE,*HERE.parents,Path.cwd(),*Path.cwd().parents) if (p/'research/MATH_STATUS.json').is_file()),None)
    if repo is None:raise SystemExit('Use --repo-root outside the checkout.')
    ledger_path=HERE/'coalescing_binomial_local_source_ledger.json'
    ledger=json.loads(ledger_path.read_text())
    assert ledger['schema_version']==1 and set(ledger['files'])==set(MANDATORY_SOURCES),'Missing or extra mandatory source record'
    for key,(scope,name) in MANDATORY_SOURCES.items():
        rec=ledger['files'][key]
        assert rec['scope']==scope and rec['path']==name
        p=(HERE if scope=='packet' else repo)/name
        assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==rec['sha256'],'Stale source: '+key
    return ledger_path

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path);ap.add_argument('--output',type=Path,default=HERE/'coalescing_binomial_saddle_bridge_verification.json');args=ap.parse_args();ledger_path=verify_sources(args.repo_root)
 cases=[];jet_controls=0;gap2=0
 for d in range(1,11):
  boundary=F((d-1)**2,4)
  parameters={int(boundary)+1,d*d+1,10*d*d+1}
  if d==2:parameters.add(4)
  for a in sorted(parameters):
   assert a>boundary;R=Quad(a,d);D=a+d
   sm_a=R.sub(S,(F(a),F(0)));Dm_s=R.sub((F(D),F(0)),S)
   w=R.div(R.mul(R.power(S,d),Dm_s),sm_a)
   w_reduced=R.scale(R.mul(R.power(S,d),R.sub((F(D*(d-1)),F(0)),R.scale(S,d))),F(1,a))
   assert w==w_reduced and all(v.denominator==1 for v in w)
   norm=F(a)**(d-1)*F(D)**(d+1)
   assert R.norm(w)==norm and R.trace(w).denominator==1
   L=[F(2),R.t]
   for k in range(2,d+2):L.append(R.t*L[-1]-R.p*L[-2])
   assert R.trace(w)==F(D*(d-1)*L[d]-d*L[d+1],a)
   assert R.trace(w)**2-4*norm<0
   if d==2:
    assert w==(F(3*a*D),F(1-4*a))
    assert R.trace(w)==-2*a*a+10*a+1
    assert R.trace(w)**2-4*norm==(1-4*a)**3;gap2+=1
   g=[]
   for k in range(5):
    lo=R.scale(R.mul(w,R.power(S,a-k)),F(falling(a,k)))if a>=k else ZERO
    hi=R.scale(R.power(S,D-k),F(falling(D,k)))if D>=k else ZERO
    g.append(R.add(lo,hi))
   ratios=[R.div(v,g[0])for v in g]
   assert ratios[1]==ratios[2]==ONE
   phi3=R.add(R.sub(ratios[3],R.scale(R.mul(ratios[2],ratios[1]),3)),R.scale(R.power(ratios[1],3),2))
   expected=R.neg(R.div(R.sub(R.scale(S,2),(R.t,F(0))),R.power(S,2)))
   assert phi3==expected and phi3!=ZERO
   assert R.norm(phi3)==F(4*a-(d-1)**2,(a*D)**2)
   phi_sw=R.neg(R.div(R.scale(R.power(S,d-1),d),R.power(R.add(w,R.power(S,d)),2)))
   assert phi_sw==R.neg(R.div(R.power(sm_a,2),R.scale(R.power(S,d+1),d)))
   assert R.neg(R.mul(phi_sw,w))==(F(1,d),F(0))
   gcrit=R.div(R.scale(R.power(S,D),d),sm_a);assert gcrit==g[0]
   invg=R.inv(g[0]);dinvg=R.neg(R.div(g[1],R.power(g[0],2)))
   assert dinvg==R.neg(invg)
   det=R.scale(R.inv(R.mul(g[0],phi3)),-2)
   assert R.norm(det)>0
   # Formal amplitude-jet determinant for arbitrary b,c: b(c/g-b²/g)-c*b/g=-b³/g.
   # Compare all coefficients in the b polynomial, with a nonzero arbitrary c.
   c=(F(7),F(3));q0=[ZERO,invg];q1=[R.mul(c,invg),ZERO,dinvg]
   left=polymul(R,[ZERO,ONE],q1)
   right=[R.mul(c,u)for u in q0]
   difference=left[:]
   for k,v in enumerate(right):difference[k]=R.sub(difference[k],v)
   assert difference==[ZERO,ZERO,ZERO,R.neg(invg)];jet_controls+=1
   h=F(d-1,2);m=F(a)+h;v2=F(a)-h*h
   assert m*m-(d*d-1)*v2==(F(a)-F(d*(d-1),2))**2
   if d>=2:assert m*m-v2>0
   cases.append({'a':a,'d':d,'trace':str(R.trace(w)),'norm':str(norm),'critical_value_integer_coefficients':[str(q)for q in w]})
 quartics=0
 for r in range(1,7):
  a=r*r;d=2*r+1;D=(r+1)**2;s=F(r*(r+1));w=F(r**(2*r)*(r+1)**(2*r+2))
  g=[w*falling(a,k)*s**(a-k)+falling(D,k)*s**(D-k)for k in range(5)]
  q=[z/g[0]for z in g]
  assert q[1]==q[2]==q[3]==1
  phi4=q[4]-4*q[3]*q[1]-3*q[2]**2+12*q[2]*q[1]**2-6*q[1]**4
  assert phi4==-2/(s*s);quartics+=1
 shifts=0
 for a in range(1,6):
  for d in range(1,6):
   for N in range(2,7):actual_shift(a,d,N);shifts+=1
 # Corruption guards distinguish actual shape shifts and wrong critical norms.
 true=actual_shift(4,2,5);bad=[F(comb(4,k)*factorial(4*5+2*k),5**(4*5+2*k+1))for k in range(5)]
 assert true!=bad
 R4=Quad(4,2);w4=R4.div(R4.mul(R4.power(S,2),R4.sub((F(6),F(0)),S)),R4.sub(S,(F(4),F(0))))
 assert R4.norm(w4)!=4**2*6**3 and R4.trace(w4)==9 and R4.norm(w4)==864
 # A real cubic-root direction cannot satisfy the nonreal simple-turning cube.
 R=Quad(4,2);p3=R.neg(R.div(R.sub(R.scale(S,2),(R.t,F(0))),R.power(S,2)))
 assert R.inv(p3)[1]!=0
 proof=HERE/'COALESCING_BINOMIAL_SADDLE_BRIDGE.md';assert proof.is_file()
 def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
 out={'status':'PASS','scope':'Local exact critical-circle/integrality/Lucas/turning/actual-amplitude and radial-shift controls only. No global cycle, uniform remainder or SFC theorem.','counts':{'quadratic_algebra_cases':len(cases),'gap_two_explicit_cases':gap2,'amplitude_jet_controls':jet_controls,'quartic_boundary_cases':quartics,'actual_phase_shift_identities':shifts,'negative_controls':3},'files':[{'path':p.name,'sha256':sha(p)}for p in[proof,Path(__file__).resolve(),ledger_path]],'cases':cases}
 args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','counts':out['counts']}))
if __name__=='__main__':main()
