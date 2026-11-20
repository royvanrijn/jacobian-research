#!/usr/bin/env python3
"""Exact support for global ACTUAL ordinary gap-two zero localization."""
if not __debug__:raise SystemExit("Unoptimized Python is required for exact checks.")
from pathlib import Path
from fractions import Fraction as F
from math import comb
import argparse,hashlib,json,copy
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def padd(a,b):
 out=[F(0)]*max(len(a),len(b))
 for i,v in enumerate(a):out[i]+=v
 for i,v in enumerate(b):out[i]+=v
 return out
def scale(a,c):return [v*c for v in a]
def pmul(a,b):
 out=[F(0)]*(len(a)+len(b)-1)
 for i,u in enumerate(a):
  for j,v in enumerate(b):out[i+j]+=u*v
 return out
def peval(a,x):
 out=F(0)
 for c in a[::-1]:out=out*x+c
 return out
def eq(a,b):return not any(padd(a,scale(b,-1)))
def rising(x,k):
 v=F(1)
 for j in range(k):v*=x+j
 return v
def gamma_poly(n,g):return [F(comb(n,j))*rising(g,n-j)for j in range(n+1)]
def actual(n,a):return [F(comb(n,j))*rising(a*n+1,2*(n-j))for j in range(n+1)]
def diamond(p,q,n):return [p[k]*q[k]/comb(n,k)for k in range(n+1)]
def esym(xs):
 out=[F(1)]+[F(0)]*len(xs)
 for x in xs:
  for k in range(len(xs),0,-1):out[k]+=x*out[k-1]
 return out
def cadd(x,y):return (x[0]+y[0],x[1]+y[1])
def cneg(x):return(-x[0],-x[1])
def cmul(x,y):return(x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def norm(x):return x[0]*x[0]+x[1]*x[1]
def cinv(x):
 q=norm(x);assert q
 return(x[0]/q,-x[1]/q)
# Polynomial identities in independent formal a,t.
def poly(*terms):
 out={}
 for i,j,c in terms:out[(i,j)]=out.get((i,j),F(0))+F(c)
 return{k:v for k,v in out.items()if v}
def add(*args):
 out={}
 for p in args:
  for e,c in p.items():out[e]=out.get(e,F(0))+c
 return{k:v for k,v in out.items()if v}
def times(p,q):
 out={}
 for u,c in p.items():
  for v,d in q.items():
   e=(u[0]+v[0],u[1]+v[1]);out[e]=out.get(e,F(0))+c*d
 return{k:v for k,v in out.items()if v}
def ps(p,c):return{k:v*c for k,v in p.items()if v*c}
def validate_input(d):
 assert d=={"schema_version":1,"minimum_lower_exponent":40000,
 "minimum_terminal":2,"physical_gap":2,"gamma_parameter_samples":["1/2","1","5/3","10"],
 "matrix_maximum_dimension":16,"actual_fixture_lower_exponents":[1,2,4,6,40000],
 "actual_fixture_orders":list(range(1,9)),"polar_identity_maximum_degree":7,
 "proof_scope":"Every actual upper-moment zero lies in the global Hermite rectangle; no consecutive coprimality claim."}
def run_controls(inp):
 recurrences=0;duplications=0;polar=0;circle=0
 # Every successive principal minor: fixed alpha=-gamma-n throughout.
 for n in range(1,inp["matrix_maximum_dimension"]+1):
  for g in map(F,inp["gamma_parameter_samples"]):
   prev=[F(1)];before=None
   for k in range(1,n+1):
    b=2*k-1-g-n
    current=pmul([-b,F(1)],prev)
    if k>1:current=padd(current,scale(before,(k-1)*(g+n-k+1)))
    assert eq(current,gamma_poly(k,g+n-k));recurrences+=1
    before,prev=prev,current
   assert eq(prev,gamma_poly(n,g))
 for a in inp["actual_fixture_lower_exponents"]:
  for n in inp["actual_fixture_orders"]:
   p=gamma_poly(n,F(a*n+1,2));q=gamma_poly(n,F(a*n+2,2))
   h=diamond(p,q,n)
   assert actual(n,a)==[c*4**(n-j)for j,c in enumerate(h)];duplications+=1
 # Full polynomial polar identity in zeta, via degree+1 exact rational controls.
 for n in range(1,inp["polar_identity_maximum_degree"]+1):
  roots=[F(j+1)for j in range(n)];q=[F(1)]
  for b in roots:q=pmul(q,[-b,F(1)])
  p=[F((-1)**j*(j+2),j+1)for j in range(n+1)]
  h=diamond(p,q,n);assert q[0]
  for zeta in range(n+1):
   sym=esym([-F(zeta)/b for b in roots])
   value=sum((p[k]/comb(n,k)*sym[k]for k in range(n+1)),F(0))
   assert value==peval(h,F(zeta))/q[0];polar+=1
 # Exact inverse-disc geometry, including closed-boundary points.
 points=[(F(0),F(0)),(F(1),F(0)),(F(-1),F(0)),(F(0),F(1)),
         (F(0),F(-1)),(F(3,5),F(4,5)),(F(-3,5),F(4,5))]
 for z in [(F(2),F(0)),(F(3,2),F(1)),(F(-2),F(1,2))]:
  K=norm(z)-1;assert K>0
  center=(z[0]/K,-z[1]/K);radius2=1/(K*K);images=[]
  for r in points:
   t=cinv(cadd(z,cneg(r)));images.append(t)
   assert norm(cadd(t,cneg(center)))<=radius2;circle+=1
  avg=(sum(t[0]for t in images)/len(images),sum(t[1]for t in images)/len(images))
  assert norm(cadd(avg,cneg(center)))<=radius2 and norm(avg)>0
  assert norm(cadd(z,cneg(cinv(avg))))<=1;circle+=1
 # Exact complete coefficient identities; no sampling a or N establishes these.
 a=poly((1,0,1));t=poly((0,1,1));one=poly((0,0,1))
 h=add(one,ps(t,-1));g1=ps(add(a,t),F(1,2));g2=ps(add(a,ps(t,2)),F(1,2))
 M=times(g1,g2);E=add(ps(times(h,add(g1,g2)),3),ps(times(h,h),9))
 U=add(times(add(ps(g1,5),ps(g2,3)),h),ps(times(h,h),9))
 identities=[
 (U,poly((1,0,4),(0,0,9),(1,1,-4),(0,1,F(-25,2)),(0,2,F(7,2)))),
 (add(E,U),poly((1,0,7),(0,0,18),(1,1,-7),(0,1,-26),(0,2,8))),
 (add(ps(M,2),E),poly((2,0,F(1,2)),(1,0,3),(0,0,9),(1,1,F(-3,2)),(0,1,F(-27,2)),(0,2,F(11,2)))),
 (add(times(a,a),ps(M,-4),ps(U,4)),poly((1,0,16),(1,1,-19),(0,0,36),(0,1,-50),(0,2,12))),
 (add(times(a,a),ps(M,-4),ps(E,-4)),poly((1,0,-12),(1,1,9),(0,0,-36),(0,1,54),(0,2,-20))),
 (ps(times(poly((1,0,7),(0,0,18)),poly((2,0,F(1,2)),(1,0,3),(0,0,9))),4),
  poly((3,0,14),(2,0,120),(1,0,468),(0,0,648)))]
 assert all(x==y for x,y in identities)
 guards=[]
 def ck(label,l,r):
  l,r=F(l),F(r);assert l<r,(label,l,r)
  guards.append({"label":label,"left":str(l),"relation":"<","right":str(r)})
 A=F(inp["minimum_lower_exponent"])
 ck("M exceeds U in the all-a range",4*A+9,A*A/4)
 ck("M-U lower-bound derivative positive",4,A/2)
 ck("all-root real width",F(8,200)+F(18,40000*200),F(1,20))
 ck("all-root imaginary width squared",14+120/A+468/(A*A)+648/(A**3),16)
 # Sign of every removed t term, for 0<t<=1/2 and positive a.
 ck("U removed t coefficient",F(7,4),4*A+F(25,2))
 ck("E+U removed t coefficient",4,7*A+26)
 ck("2M+E removed t coefficient",F(11,4),F(3,2)*A+F(27,2))
 # Four actual corruptions/invalid extensions, rather than fresh sampled roots.
 negatives=0
 p=gamma_poly(1,F(3));q=gamma_poly(1,F(5))
 correct=diamond(p,q,1)
 assert correct==[F(15),F(1)] and correct!=[F(-15),F(1)];negatives+=1
 a,n=4,2
 wrong=gamma_poly(n,F(a*(n-1)+1,2));wrong2=gamma_poly(n,F(a*(n-1)+2,2))
 assert actual(n,a)!=[c*4**(n-j)for j,c in enumerate(diamond(wrong,wrong2,n))];negatives+=1
 a,n=1,2;h=F(n-1);g1=F(a*n+1,2);g2=F(a*n+2,2)
 assert (5*g1+3*g2)*h+9*h*h>g1*g2;negatives+=1
 # Wrong polar sign retains degree n rather than dropping it.
 p=[F(-1),F(0),F(1)];dp=[F(0),F(2)];alpha=F(3)
 correct_polar=padd(scale(p,2),pmul([alpha,F(-1)],dp))
 wrong_polar=padd(scale(p,2),pmul([alpha,F(1)],dp))
 assert eq(correct_polar,[F(-2),F(6)]) and eq(wrong_polar,[F(-2),F(6),F(4)])
 assert correct_polar[2]==0 and wrong_polar[2]!=0;negatives+=1
 return {"matrix_minor_coefficient_identities":recurrences,
 "actual_duplication_vectors":duplications,"polar_identity_values":polar,
 "closed_inverse_disc_controls":circle,"universal_cleared_identities":len(identities),
 "rational_domain_guards":len(guards),"negative_controls":negatives},guards
def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument("--input",type=Path,default=HERE/"gap_two_grace_global_localization_inputs.json")
 ap.add_argument("--ledger",type=Path,default=HERE/"gap_two_grace_global_localization_source_ledger.json")
 ap.add_argument("--output",type=Path,default=HERE/"gap_two_grace_global_localization_verification.json")
 args=ap.parse_args();ledger=json.loads(args.ledger.read_text());assert ledger["schema_version"]==1
 assert len(ledger["sources"])==3 and {s["role"]for s in ledger["sources"]}=={"proof","checker","input"}
 paths={}
 for row in ledger["sources"]:
  p=(args.ledger.parent/row["path"]).resolve();assert p.is_file()and sha(p)==row["sha256"],"Source pin mismatch: "+str(p)
  paths[row["role"]]=p
 assert paths["checker"]==Path(__file__).resolve()and paths["input"]==args.input.resolve()
 inp=json.loads(args.input.read_text());validate_input(inp)
 counts,guards=run_controls(inp)
 result={"schema_version":1,"status":"PASS","source_ledger_sha256":sha(args.ledger),
 "scope":"Exact support for the full written all-N global ACTUAL gap-two zero-localization proof. The finite controls do not by themselves prove the universal matrix/Grace argument. No SFC closure.",
 "counts":counts,"rational_guards":guards}
 args.output.parent.mkdir(parents=True,exist_ok=True)
 args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
 print(json.dumps({"status":"PASS","counts":counts}))
if __name__=="__main__":main()
