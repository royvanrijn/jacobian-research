#!/usr/bin/env python3
"""Exact cap24 Long jet obstruction: one modular minor and a rational left functional."""
from fractions import Fraction as Q
from pathlib import Path
from math import isqrt
import argparse,hashlib,importlib.util,json,platform,tempfile

HERE=Path(__file__).resolve().parent
def load(path,name):
 spec=importlib.util.spec_from_file_location(name,path)
 mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod

def parameter_rhs(cert,alg,f,g0,v):
 b0=cert.lower_rhs(alg,[f,g0],2);b1={};b2={}
 for pi,(i,j) in enumerate(cert.PAIRS):
  linear=alg.add(alg.bidifferential(v[i],g0[j],1),
                 alg.bidifferential(g0[i],v[j],1),
    alg.scale(alg.add(alg.bidifferential(f[i],v[j],3),
                      alg.bidifferential(v[i],f[j],3)),Q(1,24)))
  quadratic=alg.bidifferential(v[i],v[j],1)
  for e,c in linear.items():b1[(pi,e)]=-c
  for e,c in quadratic.items():b2[(pi,e)]=-c
 return b0,b1,b2

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument("--output-dir",type=Path);args=parser.parse_args()
 cert=load(HERE/"check_graded_long_cap23.py","shared_graded_polynomial_engine")
 alg=load(HERE/"check_long_graph_curvature.py","independent_polynomial_arithmetic")
 source=HERE.parent/"dixmier-audit/verify_long_pbw.py"
 producer=load(source,"canonical_long_input")
 images=[{tuple(e):Q(c) for e,c in f.items()} for f in producer.long_images()]
 witness_path=HERE/"graded_long_cap24_witness.json"
 witness=json.loads(witness_path.read_text())
 first_path=HERE/witness["first_particular_source"]
 first=json.loads(first_path.read_text())
 p=witness["prime"]
 assert p==1000003 and all(p%d for d in range(2,isqrt(p)+1))
 assert witness["cap"]==24
 g0=[{} for _ in range(4)]
 for label,e,c in first["first_particular"]:
  i=cert.LABELS.index(label);e=tuple(e);assert e not in g0[i];g0[i][e]=Q(c)
 unknowns3=cert.basis(24,1)
 assert [sum(i==j for i,e in unknowns3) for j in range(4)]==[62,91,114,80]
 a3=cert.linear_columns(alg,images,unknowns3)
 b3=cert.lower_rhs(alg,[images],1)
 rows3=sorted(set(b3).union(*(set(c) for c in a3)))
 assert cert.product(a3,cert.vector(unknowns3,g0))==b3
 bad_particular=[dict(poly) for poly in g0]
 bad_particular[3][(5,0,0,0)]+=1
 assert cert.product(a3,cert.vector(unknowns3,bad_particular))!=b3
 r,s=images[0],images[3]
 power=lambda f,n:cert.power(alg,f,n)
 kernels=[
  [{},{},alg.scale(power(r,2),-3),{}],
  [{},power(r,4),alg.scale(alg.multiply(power(r,3),s),-4),{}],
  [alg.bidifferential({(3,0,0,0):Q(1)},f,1) for f in images]]
 for kernel in kernels:
  assert not cert.product(a3,cert.vector(unknowns3,kernel))
 assert kernels[0][2] and not kernels[0][1] and not kernels[0][3]
 assert kernels[1][1] and not kernels[1][3]
 assert kernels[2][3]
 assert all(not kernel[0] for kernel in kernels) and not g0[0]
 det3=cert.minor(a3,b3,rows3,witness["h3_minor"],p)
 assert len(witness["h3_minor"]["rows"])==344
 # 347 columns, three exact independent kernels, and a344 minor: rank344
 # and the complete first family g0+alpha*K1+beta*K2+gamma*K3.
 v=kernels[2]
 gauge0=alg.add(alg.multiply(power(r,3),g0[3]),
                alg.scale(alg.bidifferential(power(r,3),s,2),Q(1,8)))
 gauge1=alg.multiply(power(r,3),v[3])
 assert max(map(sum,gauge0),default=0)<=16
 assert max(map(sum,gauge1),default=0)<=21

 unknowns5=cert.basis(24,2)
 assert [sum(i==j for i,e in unknowns5) for j in range(4)]==[24,40,54,34]
 a5=cert.linear_columns(alg,images,unknowns5)
 rhses=parameter_rhs(cert,alg,images,g0,v)
 rows5=sorted(set().union(*(set(c) for c in a5),*(set(b) for b in rhses)))
 left=witness["h5_left_functional"]
 assert left["coefficient_columns"]==len(a5)==152 and left["row_count"]==len(rows5)
 functional={}
 for row_index,key,c in left["functional"]:
  key=(key[0],tuple(key[1]))
  assert rows5[row_index]==key and key not in functional
  functional[key]=Q(c)
 def apply(poly):return sum((c*poly.get(key,Q(0)) for key,c in functional.items()),Q(0))
 assert all(not apply(column) for column in a5)
 evaluations=[apply(b) for b in rhses]
 assert evaluations==[Q(1),Q(0),Q(0)]
 bad_functional=dict(functional)
 bad_key=next(iter(rhses[0]))
 bad_functional[bad_key]=bad_functional.get(bad_key,Q(0))+1
 assert sum((c*rhses[0].get(key,Q(0)) for key,c in bad_functional.items()),Q(0))!=1
 # Therefore for every gamma, applying this SAME rational functional to
 # A*k=b0+gamma*b1+gamma^2*b2 gives the impossible equation0=1.
 def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
 receipt={
  "claim":"No even equivariant canonical Long Weyl jet of common degree24 through h5",
  "status":"PASS","environment":{"python":platform.python_version(),
    "arithmetic":"Exact Fraction polynomial identities, bare modular determinant, exact rational left functional"},
  "cap":24,"image_order":list(cert.LABELS),"image_weights":list(cert.WEIGHTS),
  "source_weights":[-1,1,4,2],"h_weight":3,
  "h3":{"basis_counts":[62,91,114,80],"unknowns":347,"rows":len(rows3),
    "nonzeros":sum(map(len,a3)),"exact_particular_from_cap23":True,
    "exact_independent_kernels":["G=R^3","G=R^4*S","G=X^3"],
    "minor_size":344,"minor_det_modp":det3,"rank_Q":344,
    "first_family":"g0+alpha*K1+beta*K2+gamma*K3",
    "all_first_R_corrections_zero":True},
  "target_gauge_h4_D_degrees":{"constant":max(map(sum,gauge0),default=0),
                              "gamma_coefficient":max(map(sum,gauge1),default=0)},
  "h5":{"basis_counts":[24,40,54,34],"unknowns":152,"rows":len(rows5),
    "rhs_nonzeros":[len(b) for b in rhses],
    "parameter":"gamma retained without source-gauge normalization",
    "left_functional_support":len(functional),"annihilates_every_A_column":True,
    "evaluations_on_b0_b1_b2":[str(c) for c in evaluations],
    "universal_compatibility":"0=1 for every gamma"},
  "invalid_controls":{"corrupted_particular_rejected":True,"corrupted_functional_rejected":True},
  "prime":p,"hashes":{"check_graded_long_cap24.py":sha(Path(__file__)),
    "GRADED_LONG_CAP24_JET_OBSTRUCTION.md":sha(HERE/"GRADED_LONG_CAP24_JET_OBSTRUCTION.md"),
    "graded_long_cap24_witness.json":sha(witness_path),
    "graded_long_cap23_witness.json":sha(first_path),
    "check_graded_long_cap23.py":sha(HERE/"check_graded_long_cap23.py"),
    "check_long_graph_curvature.py":sha(HERE/"check_long_graph_curvature.py"),
    "../dixmier-audit/verify_long_pbw.py":sha(source)},
  "coverage_boundary":"Complete finite even equivariant cap24 first/second jet family; "
    "all first parameters classified, two target gauges removed cap-preservingly, "
    "extra source parameter retained and contradicted for every value by an exact "
    "rational functional. No nongraded or odd-jet exclusion, cap25 conclusion or "
    "bounded normalization/DC1/DC2 theorem is claimed.",
  "independent_replay":False,"formal_verification":False,"external_review":False}
 out=args.output_dir or Path(tempfile.mkdtemp(prefix="graded-long-cap24-"))
 out.mkdir(parents=True,exist_ok=True)
 destination=out/"graded_long_cap24_controls.json"
 destination.write_text(json.dumps(receipt,indent=2)+"\n")
 print("PASS: h3 rank344/347; h5 exact functional gives0=1 for every gamma")
 print("First minor determinant:",det3,"prime",p,"functional support",len(functional))
 print("Receipt:",destination)

if __name__=="__main__":main()
