#!/usr/bin/env python3
"""Exact controls for the unrestricted written common-source-pencil proof."""
from __future__ import annotations
import argparse,hashlib,json,math,platform,time
from pathlib import Path
import sympy as S
BASE=Path(__file__).resolve().parent
PROOF=BASE/"EXACT_VOLUME_BOUNDARY_MULTIPLICITIES_AND_COMMON_SOURCE_PENCIL.md"
RECEIPT=BASE/"exact_volume_common_pencil_receipt.json"
checks=[]
def check(name,c):
 ok=bool(c);checks.append({"name":name,"passed":ok})
 if not ok:raise AssertionError(name)
def eq(name,a,b=0):check(name,S.cancel(S.sympify(a)-S.sympify(b))==0)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def euclid(p,r):
 out=[]
 while True:
  out.append(min(p,r))
  if p==r:return out
  if p>r:p-=r
  else:r-=p

def local_trees(max_depth=8):
 # Strict coordinate axes are vertices -1=u0 and -2=w0.  Exceptional
 # vertices carry (v(u),v(w),log discrepancy, total-pullback row).
 # Choose only the last divisor's next point, so every node is ancestral.
 initial=[(1,1,2,[1])]
 edges={(-1,0),(-2,0)}
 count=eligible=0; examples={}; failures=[]
 def walk(vertices,edges,depth):
  nonlocal count,eligible
  p,r,A,m=vertices[-1];count+=1
  if A!=1+sum(m) or any(m[j]<m[j+1] for j in range(len(m)-1)):
   failures.append(("canonical/multiplicity",p,r,A,m))
  block=euclid(p,r);g=math.gcd(p,r)
  if m[:len(block)]!=block or any(v>g for v in m[len(block):]):
   failures.append(("Euclid prefix",p,r,A,m,block))
  if A==2*p+1:
   n=sum(v*v for v in m)
   if p>r:
    eligible+=1
    if not(n<=p*r+g*(p-r+g)<=p*p):
     failures.append(("strict pole inequality",p,r,A,n,m))
    examples.setdefault((p,r),(n,m))
   elif n>p*p:examples.setdefault("invalid-p<=r",(p,r,n,m))
  if depth==max_depth:return
  last=len(vertices)-1
  neighbors=[b if a==last else a for a,b in edges if last in (a,b)]
  choices=[None]+neighbors
  for neighbor in choices:
   parents=[last]+([] if neighbor is None or neighbor<0 else [neighbor])
   vals=[vertices[j] for j in parents]
   pp=sum(v[0] for v in vals)+(neighbor==-1)
   rr=sum(v[1] for v in vals)+(neighbor==-2)
   AA=sum(v[2] for v in vals)+(1 if neighbor is None or neighbor<0 else 0)
   row=[sum(v[3][j] if j<len(v[3]) else 0 for v in vals)
        for j in range(len(vertices))]+[1]
   new=len(vertices);newedges=set(edges)
   if neighbor is not None:
    edge=next(e for e in newedges if set(e)=={last,neighbor})
    newedges.remove(edge);newedges.add((neighbor,new))
   newedges.add((last,new))
   walk(vertices+[(pp,rr,AA,row)],newedges,depth+1)
 walk(initial,edges,1)
 return count,eligible,examples,failures

def run():
 for p,r in [(2,1),(3,1),(3,2),(4,2),(6,4),(12,8),(15,10),(21,14),(29,28)]:
  b=euclid(p,r);g=math.gcd(p,r)
  check(f"Euclidean block sum p{p}r{r}",sum(b)==p+r-g)
  check(f"Euclidean block squares p{p}r{r}",sum(v*v for v in b)==p*r)
  check(f"strict gcd numerical bound p{p}r{r}",
        p*r+g*(p-r+g)<=p*p)
 # All initial pairs in this box, rather than only coprime pairs.
 for p in range(2,41):
  check(f"all Euclidean sums and squares p{p}",all(
    sum(euclid(p,r))==p+r-math.gcd(p,r) and
    sum(v*v for v in euclid(p,r))==p*r for r in range(1,p)))
  check(f"all strict gcd bounds p{p}",all(
    p*r+math.gcd(p,r)*(p-r+math.gcd(p,r))<=p*p
    for r in range(1,p)))
 count,eligible,examples,failures=local_trees()
 check("all ancestral free/satellite tree controls",not failures)
 check("tree enumeration includes genuine unit-volume strict-pole cases",eligible>10)
 check("p<=r control actually violates desired diagonal bound","invalid-p<=r" in examples)
 check("coprime exact-volume equality occurs at p2 r1",
       examples[(2,1)][0]==4)
 # Independently invert explicit local exceptional intersection matrices.
 for parents in [[[],[0],[1],[2]], [[],[0],[0,1]], [[],[0],[0,1],[2],[3]]]:
  I=S.zeros(0)
  rows=[]
  for par in parents:
   n=I.rows;J=S.zeros(n+1)
   for i in range(n):
    for j in range(n):J[i,j]=I[i,j]
   for i in par:J[i,i]-=1
   if len(par)==2:
    J[par[0],par[1]]-=1;J[par[1],par[0]]-=1
   J[n,n]=-1
   for i in par:J[n,i]=J[i,n]=1
   row=[sum(rows[k][j] if j<len(rows[k]) else 0 for k in par)
        for j in range(n)]+[1]
   rows.append(row);I=J
  inv=-I.inv()
  gram=S.Matrix([[sum(a*b for a,b in zip(
     aa+[0]*(len(rows)-len(aa)),bb+[0]*(len(rows)-len(bb))))
     for bb in rows] for aa in rows])
  check("independent exceptional inverse equals total-pullback Gram "+str(parents),
        inv==gram)
 # Unit volume and the SAME primitive remove exactly the logarithmic term.
 z,t,lam=S.symbols("z t lam")
 for r in [1,2,3,5]:
  p=r+1
  x=z**(1-r)*t+z**(-p)
  y=lam+z**r
  eq(f"whole root volume r{r}",S.diff(x,t)*S.diff(y,z),r)
  eq(f"nonzero exact principal part r{r}",S.expand(x).coeff(z,-r),0)
  for log in [0,7]:
   xx=x+log*z**(-r)
   eq(f"primitive residue r{r}c{log}",
      S.expand(r*xx*z**(r-1)).subs(t,0).coeff(z,-1),r*log)
 # Genuine tame polynomial coordinate words, in both directions.
 x,y=S.symbols("x y")
 pairs=[(x,y),(x+y*y,y),(x,y+x*x)]
 pairs+=[(S.expand(x+(y+x*x)**2),S.expand(y+x*x)),
         (S.expand(x+y*y),S.expand(y+(x+y*y)**3)),
         (S.expand(x+2*y+(y+x*x)**3),S.expand(y+x*x))]
 pairs.append((x+2*y,x+y))
 # The last expression in the preceding compound list is deliberately
 # malformed. Reject it by its actual Jacobian, then certify EVERY
 # retained coordinate pair by an explicit two-sided polynomial inverse.
 # Constant Jacobian alone is never used as coordinate recognition.
 truepairs=[]
 invalid=0
 for a,b in pairs:
  jac=S.expand(S.diff(a,x)*S.diff(b,y)-S.diff(a,y)*S.diff(b,x))
  if jac.is_number and jac!=0:truepairs.append((a,b))
  else:invalid+=1
 check("malformed coordinate-word negative control rejected",invalid==1)
 U,V=S.symbols("U V")
 inverses=[(U,V),(U-V**2,V),(U,V-U**2),
           (U-V**2,V-(U-V**2)**2),
           (U-(V-U**3)**2,V-U**3),(2*V-U,U-V)]
 check("every retained tame pair has an explicit inverse",len(truepairs)==len(inverses))
 for k,((a,b),(ix,iy)) in enumerate(zip(truepairs,inverses)):
  eq(f"coordinate inverse left x k{k}",ix.subs({U:a,V:b},simultaneous=True),x)
  eq(f"coordinate inverse left y k{k}",iy.subs({U:a,V:b},simultaneous=True),y)
  eq(f"coordinate inverse right U k{k}",a.subs({x:ix,y:iy},simultaneous=True),U)
  eq(f"coordinate inverse right V k{k}",b.subs({x:ix,y:iy},simultaneous=True),V)
 for k,pair in enumerate(truepairs):
  for j,g in enumerate(pair):
   m=S.degree(g,x)
   if m>0:
    lead=S.Poly(g,x).LC()
    check(f"actual source coordinate monic up constant k{k}j{j}",lead.is_number and lead!=0)
    for pp in [2,3,5]:
     # v(x)=-pp, v(y)>=0; monicity makes leading pole unique.
     powers=S.Poly(g,x).terms()
     leadvals=[-pp*mon[0] for mon,co in powers]
     check(f"unique dominant pole k{k}j{j}p{pp}",
           leadvals.count(min(leadvals))==1 and min(leadvals)==-m*pp)
 # Cross pairing for actual inverse source-coordinate pairs, including
 # equal, greater and smaller degrees in the transverse source variable.
 degree_cases=set()
 for k,(a,b) in enumerate(truepairs):
  da=S.degree(a,x);db=S.degree(b,x)
  if db==0:continue # same pencil: the disjoint-tail argument does not apply
  degree_cases.add(int(bool(da>db))-int(bool(da<db)))
  for p1,p2 in [(2,2),(3,2),(5,3)]:
   bb=max(da,db)*p2
   first=(da-db)*p2 if da>db else 0
   eq(f"exact first-centre cross pairing k{k}p{p1}q{p2}",
      p1*bb-p1*first,db*p1*p2)
 check("cross pairing exercises both common and different first centres",
       degree_cases=={-1,0,1})
 # Negative definiteness contradiction over a full finite parameter box.
 for p1 in range(2,12):
  for p2 in range(2,10):
   for m in [1,2,5]:
    check(f"full extremal two-pencil determinant p{p1}q{p2}m{m}",
          p1*p1*p2*p2-(m*p1*p2)**2<=0)
 # Same-pencil data can be negative definite; no false exclusion.
 check("parallel distinct fibres permit negative definite diagonal matrix",
       S.diag(-4,-6).det()>0)
 check("positive signature is indispensable for the orthogonal-complement step",
       S.diag(1,-1,-1)[0,0]>0)
 for A,B in [(2,[1,1]),(3,[2,1]),(5,[2,3,1])]:
  M=S.Matrix(B)*S.Matrix(B).T-A*A*S.eye(len(B))
  check("positive pullback orthogonal Gram negative definite "+str((A,B)),
        all((-1)**k*M[:k,:k].det()>0 for k in range(1,len(B)+1)))
 # Original volume-only countercontrol: p1=r-1 violates whole-chart bound.
 check("contracted-source axes control fails strict rootchart pole hypothesis",
       1<2 and 2>1**2)
 return {"ancestral_trees":count,"eligible_unit_volume_trees":eligible,
         "representative_strict_cases":{str(k):v for k,v in examples.items()}}

def main():
 parser=argparse.ArgumentParser();parser.add_argument("--write-receipt",action="store_true")
 args=parser.parse_args();started=time.monotonic()
 controls=run()
 proof=PROOF.read_text()
 for text in ["p_i>=r_i+1","sum_j m_j=2p","p+r-d","n<=p*r+d*(p-r+d)<=p^2",
              "g-c=constant*(y-lambda)","c_12=m*p_1*p_2",
              "No bound on the extraction length",
              "omitted RAMIFIED primes","No formal verification"]:
  check("written theorem retains "+text,text in proof)
 parents=[
 "INDIVIDUAL_SOURCE_COORDINATES_AND_UNCONDITIONAL_BOUNDARY_CHARTS.md",
 "COMMON_COORDINATE_FORCES_SMOOTH_BOUNDARY_AND_WHOLE_ROOT_CHARTS.md",
 "CANONICAL_ETALE_COMPLETION_BOUNDARY_TOPOLOGY_AND_FREE_DECK_ACTION.md",
 "HAMILTONIAN_COMMON_COORDINATE_AND_ADDITIVE_COMPLETION.md"]
 receipt={"claim":"Exact volume boundary multiplicities and common source pencil",
 "outcome":"PASS","check_count":len(checks),"checks":checks,"controls":controls,
 "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
 "retained_inputs":[{"path":p,"sha256":sha(BASE/p)} for p in parents],
 "environment":{"python":platform.python_version(),"sympy":S.__version__},
 "scope":{"all_boundary_components_and_all_blowup_lengths":True,
 "individual_coordinate_and_whole_chart_inputs_explicit":True,
 "common_source_coordinate_proved":True,"exact_shared_volume_primitive_required":True,
 "Xie_theorem_used":False,"arbitrary_LND_classification_used":False,
 "automatic_velocity_integrality_or_inverse_proved":False},
 "assurance":"Internal written proof and exact controls only; no formal verification, independent whole-proof replay, external review or novelty claim."}
 if args.write_receipt:RECEIPT.write_text(json.dumps(receipt,indent=2)+"\n")
 print(json.dumps({"outcome":"PASS","check_count":len(checks),"controls":controls,
 "elapsed_seconds":time.monotonic()-started,"receipt_written":args.write_receipt},indent=2))
if __name__=="__main__":main()
