#!/usr/bin/env python3
"""Exact replay for the complete even graded Long cap25 two-jet obstruction."""
from fractions import Fraction as Q
from pathlib import Path
from math import isqrt
import argparse,hashlib,importlib.util,json,platform,tempfile
HERE=Path(__file__).resolve().parent
def load(path,name):
 spec=importlib.util.spec_from_file_location(name,path)
 mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument("--output-dir",type=Path);args=parser.parse_args()
 cert=load(HERE/"check_graded_long_cap23.py","graded_polynomial_engine")
 param=load(HERE/"check_graded_long_cap24.py","parameter_rhs_engine")
 alg=load(HERE/"check_long_graph_curvature.py","independent_polynomial_arithmetic")
 source=HERE.parent/"dixmier-audit/verify_long_pbw.py"
 producer=load(source,"canonical_long_input")
 f=[{tuple(e):Q(c) for e,c in poly.items()} for poly in producer.long_images()]
 witness_path=HERE/"graded_long_cap25_witness.json"
 witness=json.loads(witness_path.read_text())
 first_path=HERE/witness["first_particular_source"]
 first=json.loads(first_path.read_text());p=witness["prime"]
 assert p==1000003 and all(p%d for d in range(2,isqrt(p)+1))
 assert witness["cap"]==25
 g0=[{} for _ in range(4)]
 for name,e,c in first["first_particular"]:
  i=cert.LABELS.index(name);e=tuple(e);assert e not in g0[i];g0[i][e]=Q(c)
 power=lambda poly,n:cert.power(alg,poly,n)
 x={(1,0,0,0):Q(1)};q={(0,1,0,0):Q(1)}
 mom={(0,0,1,0):Q(1)};z={(0,0,0,1):Q(1)}
 one={(0,0,0,0):Q(1)};r,s=f[0],f[3]
 beta=alg.add(alg.scale(alg.multiply(power(x,2),mom),3),
  alg.scale(alg.multiply(alg.add(one,alg.scale(alg.multiply(x,q),-3)),z),2),
  alg.scale(power(q,2),-9))
 y=alg.add(q,alg.scale(alg.multiply(x,beta),Q(-1,3)))
 r4=power(r,4);assert all(e[0]>=4 for e in r4)
 # G4=R^4*a=R^4*y+R^4/X is an ACTUAL polynomial.
 G4=alg.add(alg.multiply(r4,y),
             {(e[0]-1,e[1],e[2],e[3]):c for e,c in r4.items()})
 k1=[{},{},alg.scale(power(r,2),-3),{}]
 k2=[{},r4,alg.scale(alg.multiply(power(r,3),s),-4),{}]
 k3=[alg.bidifferential(power(x,3),fi,1) for fi in f]
 k4=[alg.bidifferential(G4,fi,1) for fi in f]
 kernels=[k1,k2,k3,k4]
 u3=cert.basis(25,1);a3=cert.linear_columns(alg,f,u3)
 b3=cert.lower_rhs(alg,[f],1)
 rows3=sorted(set(b3).union(*(set(c) for c in a3)))
 assert [sum(i==j for i,e in u3) for j in range(4)]==[71,102,127,91]
 assert cert.product(a3,cert.vector(u3,g0))==b3
 bad=[dict(poly) for poly in g0];bad[3][(5,0,0,0)]+=1
 assert cert.product(a3,cert.vector(u3,bad))!=b3
 for kernel in kernels:assert not cert.product(a3,cert.vector(u3,kernel))
 assert all(not kernel[0] for kernel in kernels) and not g0[0]
 assert k1[2] and not k1[1] and not k1[3]
 assert k2[1] and not k2[3]
 assert any(e[2]+e[3]>0 for e in k3[3])
 assert k4[3] and all(e[2]+e[3]==0 for e in k4[3])
 assert k4[3]==alg.scale(alg.multiply(x,r4),2)
 det3=cert.minor(a3,b3,rows3,witness["h3_minor"],p)
 assert len(witness["h3_minor"]["rows"])==387
 gauge0=alg.add(alg.multiply(power(r,3),g0[3]),
               alg.scale(alg.bidifferential(power(r,3),s,2),Q(1,8)))
 gaugeg=alg.multiply(power(r,3),k3[3])
 gauged=alg.multiply(power(r,3),k4[3])
 assert max(map(sum,gauge0),default=0)<=16
 assert max(map(sum,gaugeg),default=0)<=21
 assert max(map(sum,gauged),default=0)<=22
 # Exact flow controls. The all-n non-LND conclusion is a written induction,
 # not an inference from this prefix, and is only about this exponential.
 assert alg.bidifferential(G4,r,1)=={}
 assert alg.bidifferential(G4,x,1)==alg.scale(alg.multiply(r4,power(x,3)),-1)
 second_s=alg.bidifferential(G4,k4[3],1)
 assert second_s==alg.scale(alg.multiply(power(r,8),power(x,3)),-2)
 assert max(map(sum,second_s))==27

 b0,bg,bgg=param.parameter_rhs(cert,alg,f,g0,k3)
 _,bd,bdd=param.parameter_rhs(cert,alg,f,g0,k4)
 bgd={}
 for pi,(i,j) in enumerate(cert.PAIRS):
  poly=alg.add(alg.bidifferential(k3[i],k4[j],1),
               alg.bidifferential(k4[i],k3[j],1))
  for e,c in poly.items():bgd[(pi,e)]=-c
 rhses=[b0,bg,bd,bgg,bgd,bdd]
 u5=cert.basis(25,2);a5=cert.linear_columns(alg,f,u5)
 assert [sum(i==j for i,e in u5) for j in range(4)]==[29,47,62,40]
 rows5=sorted(set().union(*(set(c) for c in a5),*(set(b) for b in rhses)))
 left=witness["h5_left_functional"]
 assert left["coefficient_columns"]==len(a5)==178 and left["row_count"]==len(rows5)
 functional={}
 for index,key,c in left["functional"]:
  key=(key[0],tuple(key[1]))
  assert rows5[index]==key and key not in functional
  functional[key]=Q(c)
 def apply(poly):return sum((c*poly.get(key,Q(0)) for key,c in functional.items()),Q(0))
 assert all(not apply(column) for column in a5)
 values=[apply(b) for b in rhses]
 assert values==[Q(1),Q(0),Q(0),Q(0),Q(0),Q(0)]
 badleft=dict(functional);key=next(iter(b0))
 badleft[key]=badleft.get(key,Q(0))+1
 assert sum((c*b0.get(key,Q(0)) for key,c in badleft.items()),Q(0))!=1

 def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
 receipt={
  "claim":"No even equivariant canonical Long full Weyl jet of common degree25 through h5",
  "status":"PASS","environment":{"python":platform.python_version(),
    "arithmetic":"Exact Fraction polynomial identities, bare modular determinant and exact rational functional"},
  "cap":25,"image_order":list(cert.LABELS),"image_weights":list(cert.WEIGHTS),
  "source_weights":[-1,1,4,2],"h_weight":3,
  "h3":{"counts":[71,102,127,91],"unknowns":391,"rows":len(rows3),
    "exact_particular_from_cap23":True,
    "exact_independent_kernels":["R^3","R^4*S","X^3","R^4*a"],
    "minor_size":387,"minor_det_modp":det3,"rank_Q":387,
    "all_first_R_corrections_zero":True},
  "new_polynomial_hamiltonian":{"formula":"R^4*a=R^4*y+R^4/X",
    "terms":len(G4),"degree":max(map(sum,G4)),
    "kernel_degrees":{name:max(map(sum,poly),default=None)
                     for name,poly in zip(cert.LABELS,k4)},
    "second_S_variation_degree":27,"source_parameter_not_removed":True},
  "target_gauge_h4_D_degrees":[max(map(sum,poly),default=0)
                              for poly in (gauge0,gaugeg,gauged)],
  "h5":{"counts":[29,47,62,40],"unknowns":178,"rows":len(rows5),
    "rhs_nonzeros":[len(b) for b in rhses],
    "parameters_retained":["gamma*{X^3,f}","delta*{R^4*a,f}"],
    "functional_support":len(functional),"all178columns_annihilated":True,
    "rhs_order":["b0","bgamma","bdelta","bgamma2","bgamma_delta","bdelta2"],
    "exact_evaluations":[str(c) for c in values],
    "compatibility":"0=1 for every gamma,delta"},
  "invalid_controls":{"corrupted_particular_rejected":True,
                       "corrupted_functional_rejected":True},
  "prime":p,"hashes":{"check_graded_long_cap25.py":sha(Path(__file__)),
    "GRADED_LONG_CAP25_JET_OBSTRUCTION.md":sha(HERE/"GRADED_LONG_CAP25_JET_OBSTRUCTION.md"),
    "graded_long_cap25_witness.json":sha(witness_path),
    "graded_long_cap23_witness.json":sha(first_path),
    "check_graded_long_cap23.py":sha(HERE/"check_graded_long_cap23.py"),
    "check_graded_long_cap24.py":sha(HERE/"check_graded_long_cap24.py"),
    "check_long_graph_curvature.py":sha(HERE/"check_long_graph_curvature.py"),
    "../dixmier-audit/verify_long_pbw.py":sha(source)},
  "coverage_boundary":"Complete finite even equivariant cap25 two-jet family, "
   "four first-kernel directions exhausted, target gauges removed with verified "
   "two-jet cap control, both source parameters retained and contradicted by "
   "a single exact rational functional. No nongraded/odd-jet theorem, cap26 "
   "search, all-gauge no-go or DC1/DC2 conclusion.",
  "independent_replay":False,"formal_verification":False,"external_review":False}
 out=args.output_dir or Path(tempfile.mkdtemp(prefix="graded-long-cap25-"))
 out.mkdir(parents=True,exist_ok=True)
 dest=out/"graded_long_cap25_controls.json"
 dest.write_text(json.dumps(receipt,indent=2)+"\n")
 print("PASS: h3 rank387/391; h5 exact functional gives0=1 for all gamma,delta")
 print("Minor determinant",det3,"prime",p,"functional support",len(functional))
 print("Receipt:",dest)

if __name__=="__main__":main()
