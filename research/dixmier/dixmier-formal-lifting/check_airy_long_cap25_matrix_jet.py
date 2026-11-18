#!/usr/bin/env python3
"""Exact admitted Airy-Long cap25 matrix jet and fixed h5 extension certificate."""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial,comb
import argparse,hashlib,importlib.util,json,platform,tempfile
HERE=Path(__file__).resolve().parent
def load(path,name):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--output-dir",type=Path);args=ap.parse_args()
 alg=load(HERE/"check_long_graph_curvature.py","independent_polynomial_engine")
 cert=load(HERE/"check_graded_long_cap23.py","complete_even_coefficient_basis")
 source=HERE.parent/"dixmier-audit/verify_long_pbw.py";prod=load(source,"canonical_Long_input")
 f=[{tuple(e):Q(c) for e,c in poly.items()} for poly in prod.long_images()]
 first_path=HERE/"graded_long_cap23_witness.json";first=json.loads(first_path.read_text())
 g0=[{} for _ in range(4)]
 for name,e,c in first["first_particular"]:g0[cert.LABELS.index(name)][tuple(e)]=Q(c)
 wp=HERE/"airy_long_cap25_matrix_jet_witness.json";w=json.loads(wp.read_text())
 def madd(*ms):
  out={}
  for m in ms:
   for ij,p in m.items():out[ij]=alg.add(out.get(ij,{}),p)
  return {ij:p for ij,p in out.items() if p}
 def scale(m,c):return {ij:q for ij,p in m.items() if (q:=alg.scale(p,c))}
 def diag(p):return {(0,0):p,(1,1):p} if p else {}
 def bidiff(a,b,n):
  out={}
  for (i,j),p in a.items():
   for (k,l),q in b.items():
    if j==k:out[i,l]=alg.add(out.get((i,l),{}),alg.bidifferential(p,q,n))
  return {ij:p for ij,p in out.items() if p}
 def md(m,i):
  out={}
  for ij,p in m.items():
   q=alg.derivative(p,tuple(int(j==i) for j in range(4)))
   if q:out[ij]=q
  return out
 def coefficient(jets,i,j,n):
  out={}
  for a in range(len(jets)):
   for b in range(len(jets)):
    ell=n-a-b
    if ell>=0:
     value=madd(bidiff(jets[a][i],jets[b][j],ell),scale(bidiff(jets[a][j],jets[b][i],ell),-1))
     out=madd(out,scale(value,Q(1,2**ell*factorial(ell))))
  return out
 def parse(rows):
  out={}
  for i,j,e,c in rows:
   e=tuple(e);assert e not in out.setdefault((i,j),{});out[i,j][e]=Q(c)
  return out
 def degree(m):return max((sum(e) for p in m.values() for e in p),default=-1)
 X={(1,0,0,0):Q(1)};Qp={(0,1,0,0):Q(1)};P={(0,0,1,0):Q(1)};Z={(0,0,0,1):Q(1)}
 one={(0,0,0,0):Q(1)};t=alg.multiply(X,Qp);M={(0,1):one,(1,0):t}
 V=lambda poly:alg.add(alg.multiply(Qp,alg.derivative(poly,(0,0,1,0))),
                       alg.multiply(X,alg.derivative(poly,(0,0,0,1))))
 power=lambda p,n:cert.power(alg,p,n);r=f[0]
 B=alg.add(alg.scale(alg.multiply(power(X,2),P),3),
   alg.scale(alg.multiply(alg.add(one,alg.scale(alg.multiply(X,Qp),-3)),Z),2))
 A0=alg.add(Qp,alg.scale(alg.multiply(X,power(Qp,2)),3),alg.scale(alg.multiply(X,B),Q(-1,3)))
 Y=alg.add(alg.multiply(r,A0),alg.scale(alg.multiply(X,Qp),-3))
 G4=alg.multiply(power(r,3),alg.add(Y,alg.scale(one,2)))
 K4=[alg.bidifferential(G4,fi,1) for fi in f]
 delta=Q(-2,27);g2=[alg.add(gi,alg.scale(ki,delta)) for gi,ki in zip(g0,K4)]
 assert not cert.product(cert.linear_columns(alg,f,cert.basis(25,1)),
                         cert.vector(cert.basis(25,1),K4))
 F0=[diag(p) for p in f]
 F1=[{ij:alg.multiply(p,V(fi)) for ij,p in M.items()} for fi in f]
 F2=[diag(alg.add(gi,alg.scale(alg.multiply(t,V(V(fi))),Q(1,2)))) for fi,gi in zip(f,g2)]
 # Finite ordered gauge-derivative recurrence reconstructs the composition
 # coefficients, without storing any analytic gauge function.
 Mi=[{ij:alg.multiply(Qp,p) for ij,p in M.items()},{ij:alg.multiply(X,p) for ij,p in M.items()}]
 U={(0,0):diag(one)};W={(0,0):diag(one)}
 for n in range(1,4):
  for a in range(n+1):
   e=(a,n-a);i=0 if a else 1;prev=list(e);prev[i]-=1;prev=tuple(prev)
   U[e]=madd(md(U[prev],i),scale(bidiff(U[prev],Mi[i],0),-1))
   W[e]=madd(bidiff(Mi[i],W[prev],0),md(W[prev],i))
 T={}
 for n in range(4):
  for a in range(n+1):
   e=(a,n-a);value={}
   for aa in range(e[0]+1):
    for bb in range(e[1]+1):
     alpha=(aa,bb);beta=(e[0]-aa,e[1]-bb)
     value=madd(value,scale(bidiff(U[alpha],W[beta],0),
                            (-1)**(aa+bb)*comb(e[0],aa)*comb(e[1],bb)))
   T[e]=value
 def J(poly,n):
  out={}
  for a in range(n+1):
   der=alg.derivative(poly,(0,0,a,n-a))
   term={ij:alg.scale(alg.multiply(p,der),Q(1,2**n*factorial(a)*factorial(n-a)))
         for ij,p in T[a,n-a].items()}
   out=madd(out,term)
  return out
 for i,fi in enumerate(f):
  assert not madd(J(fi,1),scale(F1[i],-1))
  assert not madd(madd(J(fi,2),diag(g2[i])),scale(F2[i],-1))
 F3=[madd(J(fi,3),J(gi,1)) for fi,gi in zip(f,g2)]
 jets=[F0,F1,F2,F3]
 for m in (1,2,3):
  for actual,rows in zip(jets[m],w["point"][f"F{m}"]):
   assert not madd(actual,scale(parse(rows),-1))
 for m,frame in enumerate(jets):
  for i,mat in enumerate(frame):
   assert degree(mat)<=25
   assert all(-e[0]+e[1]+4*e[2]+2*e[3]==cert.WEIGHTS[i]-3*m
              for p in mat.values() for e in p)
 for i,j in cert.PAIRS:
  for n in range(5):
   expected=diag({(0,0,0,0):Q(-1)}) if n==1 and (i,j) in ((0,2),(1,3)) else {}
   assert coefficient(jets,i,j,n)==expected,(i,j,n)
 # Known complete h5 defect for precisely these lower coefficients.
 rhs={}
 for pi,(i,j) in enumerate(cert.PAIRS):
  for (a,b),poly in scale(coefficient(jets,i,j,5),-1).items():
   for e,c in poly.items():rhs[(pi,a,b,e)]=c
 left={(pi,a,b,tuple(e)):Q(c) for pi,a,b,e,c in w["fixed_h5"]["functional"]}
 u4=cert.basis(25,2);assert len(u4)==178
 def covariant(i,C):
  return madd(bidiff(diag(f[i]),C,1),bidiff(F1[i],C,0),scale(bidiff(C,F1[i],0),-1))
 count=0;nnz=0;rowset=set(rhs)
 for component,e in u4:
  for a,b in ((0,0),(0,1),(1,0),(1,1)):
   C={(a,b):{e:Q(1)}};col={}
   for pi,(i,j) in enumerate(cert.PAIRS):
    value=covariant(i,C) if component==j else (scale(covariant(j,C),-1) if component==i else {})
    for (aa,bb),poly in value.items():
     for ee,c in poly.items():col[(pi,aa,bb,ee)]=c
   assert sum((c*col.get(key,Q(0)) for key,c in left.items()),Q(0))==0
   rowset.update(col);nnz+=len(col);count+=1
 assert count==712 and nnz==239744 and len(rowset)==7100
 assert sum((c*rhs.get(key,Q(0)) for key,c in left.items()),Q(0))==1
 assert len(rhs)==402 and len(left)==274
 # Deliberately corrupt the exact obstruction functional.
 bad=dict(left);key=next(iter(rhs));bad[key]=bad.get(key,Q(0))+1
 assert sum((c*rhs.get(key,Q(0)) for key,c in bad.items()),Q(0))!=1
 def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
 receipt={"claim":"Actual global canonical Airy-Long cap25 matrix jet through h4; no cap25 h5 extension for these fixed lower coefficients",
  "status":"PASS","environment":{"python":platform.python_version(),"arithmetic":"Fraction full matrix Moyal products and complete covariant coefficient enumeration"},
  "cap":25,"rank":2,"delta_R4a":"-2/27","point":{"all_six_relations_through_h4":True,
   "F1_degrees":[degree(m) for m in F1],"F2_degrees":[degree(m) for m in F2],
   "F3_degrees":[degree(m) for m in F3],"source_global_polynomial_entries":True,
   "graded_weights_verified":True,"compositional_origin_explicit":True},
  "fixed_h5":{"all_matrix_entries_and_images_allowed":True,"columns":count,"rows":len(rowset),
   "nonzeros":nnz,"RHS_terms":len(rhs),"exact_functional_support":len(left),
   "all_columns_annihilated":True,"RHS_evaluation":"1","scope":"Only these fixed F1,F2,F3; homogeneous F3 freedoms or other F2 choices not excluded"},
  "invalid_control":{"corrupted_functional_rejected":True},
  "hashes":{"AIRY_LONG_CAP25_MATRIX_JET_AND_FIXED_EXTENSION_OBSTRUCTION.md":sha(HERE/"AIRY_LONG_CAP25_MATRIX_JET_AND_FIXED_EXTENSION_OBSTRUCTION.md"),
   "check_airy_long_cap25_matrix_jet.py":sha(Path(__file__)),
   "airy_long_cap25_matrix_jet_witness.json":sha(wp),
   "graded_long_cap23_witness.json":sha(first_path),
   "check_graded_long_cap23.py":sha(HERE/"check_graded_long_cap23.py"),
   "check_long_graph_curvature.py":sha(HERE/"check_long_graph_curvature.py"),
   "../dixmier-audit/verify_long_pbw.py":sha(source)},
  "scope":"Finite admitted matrix jet, not an exact full Weyl tuple or a new noncompositional carrier. The h5 obstruction fixes every lower coefficient and does not classify all Airy matrix jets. No bounded formal-gauge or DC1/DC2 claim.",
  "independent_replay":False,"formal_verification":False,"external_review":False}
 out=args.output_dir or Path(tempfile.mkdtemp(prefix="airy-long-cap25-"));out.mkdir(parents=True,exist_ok=True)
 dest=out/"airy_long_cap25_matrix_jet_controls.json";dest.write_text(json.dumps(receipt,indent=2)+"\n")
 print("PASS: actual cap25 matrix jet through h4; exact fixed h5 functional on all712columns")
 print("Receipt:",dest)
if __name__=="__main__":main()
