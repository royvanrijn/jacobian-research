#!/usr/bin/env python3
"""Exact covariant Airy controls on k[t] and its crossing pullback t=XQ."""
from fractions import Fraction as Q
from pathlib import Path
import argparse,hashlib,json,platform,tempfile
HERE=Path(__file__).resolve().parent
def add(*ps):
 out={}
 for p in ps:
  for e,c in p.items():out[e]=out.get(e,Q(0))+c
 return {e:c for e,c in out.items() if c}
def scale(p,c):return {e:c*v for e,v in p.items() if c*v}
def mul(p,q):
 out={}
 for e,c in p.items():
  for f,d in q.items():
   ef=tuple(x+y for x,y in zip(e,f));out[ef]=out.get(ef,Q(0))+c*d
 return {e:c for e,c in out.items() if c}
def deriv(p,i):
 out={}
 for e,c in p.items():
  if e[i]:
   ee=list(e);ee[i]-=1;out[tuple(ee)]=c*e[i]
 return out
def vadd(*vs):return tuple(add(*(v[i] for v in vs)) for i in range(3))
def vscale(v,c):return tuple(scale(p,c) for p in v)
def vmul(v,p):return tuple(mul(z,p) for z in v)
def vd(v,i):return tuple(deriv(z,i) for z in v)
def deg(v):return max((sum(e) for p in v for e in p),default=-1)
def airyD(v):
 a,b,c=v;t={(1,):Q(1)}
 return (add(deriv(a,0),c,scale(mul(t,b),-1)),
         add(deriv(b,0),scale(a,-2)),
         add(deriv(c,0),scale(mul(t,a),2)))
def airy_inverse(v):
 u,vv,w=v;t={(1,):Q(1)}
 R=add(w,scale(deriv(u,0),-1),scale(deriv(deriv(vv,0),0),Q(-1,2)),mul(t,vv))
 N=max((e[0] for e in R),default=-1);b={}
 for n in range(N,-1,-1):
  c=(R.get((n,),Q(0))+Q((n+3)*(n+2)*(n+1),2)*b.get((n+3,),Q(0)))/Q(2*n+1)
  if c:b[(n,)]=c
 a=scale(add(deriv(b,0),scale(vv,-1)),Q(1,2))
 c=add(u,scale(deriv(a,0),-1),mul(t,b))
 out=(a,b,c);assert airyD(out)==v
 assert deg(out)<=deg(v)+2 if any(v) else not any(out)
 return out
def pull(v):return tuple({(e[0],e[0]):c for e,c in p.items()} for p in v)
X={(1,0):Q(1)};Qpos={(0,1):Q(1)};t={(1,1):Q(1)}
def adM(v):
 a,b,c=v
 return (add(c,scale(mul(t,b),-1)),scale(a,-2),scale(mul(t,a),2))
def nabla(v,i):
 return vadd(vd(v,i),vmul(adM(v),Qpos if i==0 else X))
def d0(v):return (nabla(v,0),nabla(v,1))
def d1(alpha):return vadd(nabla(alpha[1],0),vscale(nabla(alpha[0],1),-1))
def weight_piece(v,weight,formweight=0):
 return tuple({e:c for e,c in p.items() if e[0]-e[1]+formweight==weight} for p in v)
def onevar(v,shift):
 # shift QdX=(0,1), or XdQ=(1,0), or a function/2-form=(0,0).
 out=[]
 for p in v:
  q={}
  for e,c in p.items():
   ee=(e[0]-shift[0],e[1]-shift[1]);assert ee[0]==ee[1] and ee[0]>=0
   q[(ee[0],)]=c
  out.append(q)
 return tuple(out)
def solve_one(alpha):
 assert not any(d1(alpha));F,G=alpha;out=({}, {}, {})
 weights=set(e[0]-e[1]+1 for p in F for e in p)|set(e[0]-e[1]-1 for p in G for e in p)
 for m in sorted(weights-{0}):
  Fm=weight_piece(F,m,1);Gm=weight_piece(G,m,-1)
  out=vadd(out,vscale(vadd(vmul(Fm,X),vscale(vmul(Gm,Qpos),-1)),Q(1,m)))
 A=onevar(weight_piece(F,0,1),(0,1));B=onevar(weight_piece(G,0,-1),(1,0))
 assert A==B
 out=vadd(out,pull(airy_inverse(A)))
 assert d0(out)==alpha
 assert deg(out)<=max(deg(F),deg(G))+3 if any(F) or any(G) else not any(out)
 return out
def solve_two(R):
 A=({}, {}, {});B=({}, {}, {})
 weights=set(e[0]-e[1] for p in R for e in p)
 for m in sorted(weights-{0}):
  piece=weight_piece(R,m)
  A=vadd(A,vscale(vmul(piece,Qpos),Q(1,m)));B=vadd(B,vscale(vmul(piece,X),Q(1,m)))
 R0=onevar(weight_piece(R,0),(0,0));U=airy_inverse(R0)
 obstruction=tuple(p.get((0,),Q(0)) for p in U)
 if any(obstruction):return obstruction,None
 V=tuple({(e[0]-1,):c for e,c in p.items()} for p in U)
 B=vadd(B,vmul(pull(V),X))
 assert d1((A,B))==R
 assert max(deg(A),deg(B))<=deg(R)+3 if any(R) else not any(A) and not any(B)
 return obstruction,(A,B)
def encode(v):return [[[list(e),str(c)] for e,c in sorted(p.items())] for p in v]
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--output",type=Path);args=parser.parse_args()
 E12=({}, {(0,):Q(1)}, {});E21=({}, {}, {(0,):Q(1)});H=({(0,):Q(1)}, {}, {})
 # Exact low generators and a high input exercise the finite descending solver.
 for v in (E12,E21,H,({(8,):Q(3)},{(7,):Q(-2)},{(9,):Q(5),(0,):Q(7)})):
  assert airyD(airy_inverse(v))==v
 assert airy_inverse(H)==E21
 basis=[E12,E21,H];classes=[]
 for v in basis:
  R=pull(airyD(v));obstruction,solution=solve_two(R)
  assert solution is None
  expect=tuple(p.get((0,),Q(0)) for p in v);assert obstruction==expect
  classes.append({"representative":encode(R),"obstruction":[str(z) for z in obstruction]})
 # A nontrivial closed 1-form has the unique polynomial primitive.
 C=({(3,1):Q(1)},{(2,0):Q(1),(0,3):Q(2)},{(1,2):Q(-3)})
 alpha=d0(C);assert not any(d1(alpha));assert solve_one(alpha)==C
 # Positive exact 2-form, including nonzero weights and the crossing weight.
 R=vadd((t,{}, {(0,0):Q(1)}),({(2,0):Q(1)},{(0,3):Q(1)},{}))
 obstruction,solution=solve_two(R);assert obstruction==(0,0,0) and solution is not None
 # In particular H*dX wedge dQ is closed but has no polynomial primitive.
 obstruction,badsolution=solve_two(pull(H));assert obstruction==(0,0,1) and badsolution is None
 def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
 receipt={"claim":"Adjoint Airy polynomial solver and crossing cohomology H0=H1=0,H2=sl2(k)",
  "status":"PASS","environment":{"python":platform.python_version(),"arithmetic":"Fraction polynomial differentiation and products"},
  "one_variable_operator":"D=partial_t+ad([[0,1],[t,0]])",
  "scalar_elimination":"L=-partial_t^3/2+2*t*partial_t+1",
  "inverse_recurrence":"b_n=(r_n+(n+3)(n+2)(n+1)*b_(n+3)/2)/(2n+1)",
  "one_variable_degree_bound":"inputdegree+2","crossing":"t=XQ",
  "cohomology_basis_controls":classes,
  "constant_H_curvature_obstruction":["0","0","1"],
  "polynomial_closed_one_form_solver_verified":True,
  "compatible_two_form_solver_verified":True,
  "position_primitive_degree_bound":"inputdegree+3, written proof",
  "Long_covariant_cochain_bound":"inputdegree+69 for canonical degree23 frame, conditional on closedness and the three vanishing crossing coordinates",
  "scope":"Written all-polynomial cohomology and degree theorem; controls verify explicit inverse/primitive identities. No actual new finite Long tuple or bounded formal gauge is claimed.",
  "hashes":{"ADJOINT_AIRY_CROSSING_COHOMOLOGY_AND_SOLVER.md":sha(HERE/"ADJOINT_AIRY_CROSSING_COHOMOLOGY_AND_SOLVER.md"),
            "check_adjoint_airy_crossing.py":sha(Path(__file__))},
  "independent_replay":False,"formal_verification":False,"external_review":False}
 out=args.output or Path(tempfile.mkdtemp(prefix="adjoint-airy-crossing-"))/"adjoint_airy_crossing_controls.json"
 out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(receipt,indent=2)+"\n")
 print("PASS: Airy inverse, crossing obstruction basis, closed-one and compatible-two polynomial solvers")
 print("Receipt:",out)
if __name__=="__main__":main()
