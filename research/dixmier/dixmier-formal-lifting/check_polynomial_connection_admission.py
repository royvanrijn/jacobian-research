#!/usr/bin/env python3
"""Exact finite controls for the global XQ-pullback Airy source embedding."""
from fractions import Fraction as Q
from math import comb,factorial
from pathlib import Path
import argparse,hashlib,itertools,json,platform,tempfile
HERE=Path(__file__).resolve().parent
ZERO=(0,0,0,0,0)
def clean(p):return {e:c for e,c in p.items() if c}
def add(*polys):
 out={}
 for poly in polys:
  for e,c in poly.items():out[e]=out.get(e,Q(0))+c
 return clean(out)
def scale(p,c):return clean({e:c*v for e,v in p.items()})
def mul(a,b):
 """Normal PBW product X^a Q^b P^c Z^d h^m, [P,X]=[Z,Q]=h."""
 out={}
 for e,c in a.items():
  for f,d in b.items():
   for i in range(min(e[2],f[0])+1):
    ci=comb(e[2],i)*factorial(f[0])//factorial(f[0]-i)
    for j in range(min(e[3],f[1])+1):
     cj=comb(e[3],j)*factorial(f[1])//factorial(f[1]-j)
     key=(e[0]+f[0]-i,e[1]+f[1]-j,e[2]+f[2]-i,e[3]+f[3]-j,e[4]+f[4]+i+j)
     out[key]=out.get(key,Q(0))+c*d*ci*cj
 return clean(out)
def derivative(poly,i):
 out={}
 for e,c in poly.items():
  if e[i]:
   ee=list(e);ee[i]-=1;out[tuple(ee)]=c*e[i]
 return out
def madd(*matrices):
 out={}
 for mat in matrices:
  for ij,p in mat.items():out[ij]=add(out.get(ij,{}),p)
 return {ij:p for ij,p in out.items() if p}
def mscale(mat,c):return {ij:q for ij,p in mat.items() if (q:=scale(p,c))}
def mmul(a,b):
 out={}
 for (i,j),p in a.items():
  for (k,l),q in b.items():
   if j==k:out[i,l]=add(out.get((i,l),{}),mul(p,q))
 return {ij:p for ij,p in out.items() if p}
def mcomm(a,b):return madd(mmul(a,b),mscale(mmul(b,a),-1))
def diag(p):return {(0,0):p,(1,1):p} if p else {}
def mpow(a,n):
 out=diag({ZERO:Q(1)})
 for _ in range(n):out=mmul(out,a)
 return out
def mat_derivative(a,i):return {ij:q for ij,p in a.items() if (q:=derivative(p,i))}
def degree(mat,weighted=False):
 weights=(1,1,3,3,0) if weighted else (1,1,1,1,0)
 return max((sum(x*w for x,w in zip(e,weights)) for poly in mat.values() for e in poly),default=-1)
def encode(mat):
 return {f"{i},{j}":[[list(e),str(c)] for e,c in sorted(poly.items())]
         for (i,j),poly in sorted(mat.items())}
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--output",type=Path);args=parser.parse_args()
 variables=[{tuple(int(i==j) for i in range(5)):Q(1)} for j in range(5)]
 X,Qpos,P,Z,h=variables;one={ZERO:Q(1)}
 M={(0,1):one,(1,0):mul(X,Qpos)}
 MX={ij:mul(Qpos,p) for ij,p in M.items()}
 MQ={ij:mul(X,p) for ij,p in M.items()}
 flat=madd(mat_derivative(MQ,0),mscale(mat_derivative(MX,1),-1),mcomm(MX,MQ))
 assert not flat
 images=[diag(X),diag(Qpos),madd(diag(P),{ij:mul(h,p) for ij,p in MX.items()}),
         madd(diag(Z),{ij:mul(h,p) for ij,p in MQ.items()})]
 for i,j in itertools.combinations(range(4),2):
  expected=mscale(diag(h),-1) if (i,j) in ((0,2),(1,3)) else {}
  assert mcomm(images[i],images[j])==expected,(i,j)
 # A wrong relative connection sign genuinely violates a source Weyl relation.
 badZ=madd(diag(Z),mscale({ij:mul(h,p) for ij,p in MQ.items()},-1))
 assert mcomm(images[2],badZ)
 assert [degree(z) for z in images]==[1,1,3,3]
 assert [degree(z,True) for z in images]==[1,1,3,3]
 # One nontrivial filtered polynomial control; the universal theorem
 # is the written parameter-order induction, not a finite sample.
 value=madd(mpow(images[2],2),images[3],diag(mul(mul(X,X),Qpos)))
 assert degree(value,True)==6
 residue={ij:{ee:c for ee,c in poly.items() if ee[4]==0} for ij,poly in value.items()}
 residue={ij:poly for ij,poly in residue.items() if poly}
 assert residue==diag(add(mul(P,P),Z,mul(mul(X,X),Qpos)))
 # X and Q stay scalar positions; upper/lower matrices are genuinely nonconstant.
 assert MX[(0,1)]==Qpos and MX[(1,0)]==mul(X,mul(Qpos,Qpos))
 assert MQ[(0,1)]==X and MQ[(1,0)]==mul(mul(X,X),Qpos)
 def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
 receipt={"claim":"Polynomial position-flat connections reflect uniform PBW boundedness; global XQ Airy source datum",
  "status":"PASS","environment":{"python":platform.python_version(),"arithmetic":"Fraction, independent normal PBW multiplication"},
  "source_order":["X","Q","P","Z"],"brackets":"[P,X]=[Z,Q]=h",
  "datum":"M(t)=[[0,1],[t,0]],t=XQ;M_X=Q*M(t),M_Q=X*M(t)",
  "actual_global_source_images":[encode(mat) for mat in images],
  "all_six_source_relations_verified_exactly":True,"position_flatness_verified_exactly":True,
  "image_PBW_degrees":[degree(z) for z in images],"L":3,
  "weighted_generator_bounds_verified":True,
  "nontrivial_filtered_polynomial_control":"P^2+Z+X^2Q has image weighted degree6 and the prescribed scalar residue",
  "wrong_relative_sign_rejected":True,
  "universal_theorem_proof":"Written filteredness and parameter-order induction, not finite controls",
  "irreducibility_scope":"Written rational Riccati infinity argument over position field k(X,Q) only",
  "composition_scope":"Only compositional images under polynomial position-flat source representations; arbitrary finite noncompositional graph actions remain open",
  "hashes":{"POLYNOMIAL_CONNECTION_BOUNDEDNESS_REFLECTION.md":sha(HERE/"POLYNOMIAL_CONNECTION_BOUNDEDNESS_REFLECTION.md"),
            "check_polynomial_connection_admission.py":sha(Path(__file__))},
  "independent_replay":False,"formal_verification":False,"external_review":False}
 out=args.output or Path(tempfile.mkdtemp(prefix="polynomial-connection-"))/"polynomial_connection_admission_controls.json"
 out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(receipt,indent=2)+"\n")
 print("PASS: global XQ Airy source embedding; all six Weyl relations and flatness exact")
 print("Weighted generator bounds and nontrivial polynomial control PASS; receipt:",out)
if __name__=="__main__":main()
