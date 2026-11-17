#!/usr/bin/env python3
"""Replay the degree27 positive even Long h5 jet and fixed-g0 h7 obstruction."""
from pathlib import Path
from fractions import Fraction as Q
from math import isqrt
import argparse,hashlib,importlib.util,json,platform,tempfile
HERE=Path(__file__).resolve().parent
def load(path,name):
 spec=importlib.util.spec_from_file_location(name,path)
 module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument("--output-dir",type=Path);args=parser.parse_args()
 cert=load(HERE/"check_graded_long_cap23.py","linear_column_engine")
 alg=load(HERE/"check_long_graph_curvature.py","independent_polynomial_engine")
 source=HERE.parent/"dixmier-audit/verify_long_pbw.py"
 producer=load(source,"canonical_polynomial_input")
 f=[{tuple(e):Q(c) for e,c in poly.items()} for poly in producer.long_images()]
 firstpath=HERE/"graded_long_cap23_witness.json"
 first=json.loads(firstpath.read_text())
 g=[{} for _ in range(4)]
 for label,e,c in first["first_particular"]:
  g[cert.LABELS.index(label)][tuple(e)]=Q(c)
 assert sum(map(len,g))==59
 path=HERE/"minimal_even_long_jet_witness.json"
 witness=json.loads(path.read_text())
 p=witness["prime"];assert p==1000003 and all(p%d for d in range(2,isqrt(p)+1))
 def deserialize(rows):
  out=[{} for _ in range(4)]
  for label,e,c in rows:
   i=cert.LABELS.index(label);e=tuple(e)
   assert e not in out[i];out[i][e]=Q(c)
  return out
 K=deserialize(witness["positive_F4_particular"])
 kernel=deserialize(witness["positive_F4_kernel"])
 r,T,D,S=f
 add,mul,scale=alg.add,alg.multiply,alg.scale
 power=lambda poly,n:cert.power(alg,poly,n)
 variables=[{tuple(int(i==j) for i in range(4)):Q(1)} for j in range(4)]
 for i,j in cert.PAIRS:
  assert alg.bidifferential(f[i],f[j],1)==alg.bidifferential(variables[i],variables[j],1)
 for poly,weight in zip(f,cert.WEIGHTS):
  assert all(-e[0]+e[1]+4*e[2]+2*e[3]==weight for e in poly)
 u2=cert.basis(27,1);a2=cert.linear_columns(alg,f,u2)
 assert cert.product(a2,cert.vector(u2,g))==cert.lower_rhs(alg,[f],1)
 b5=cert.lower_rhs(alg,[f,g],2)
 u26=cert.basis(26,2);a26=cert.linear_columns(alg,f,u26)
 rows26=sorted(set(b5).union(*(set(c) for c in a26)))
 assert [sum(i==j for i,e in u26) for j in range(4)]==[34,54,71,47]
 assert len(u26)==witness["fixed_C26_unknowns"]==206
 assert len(rows26)==witness["fixed_C26_rows"]
 lowdet=cert.minor(a26,b5,rows26,witness["fixed_C26_augmented_minor"],p,True)
 assert len(witness["fixed_C26_augmented_minor"]["rows"])==206
 expected=[{},{},scale(power(r,8),-9),{}]
 assert kernel==expected==[alg.bidifferential(power(r,9),z,1) for z in f]
 assert kernel[2] and not cert.product(a26,cert.vector(u26,kernel))
 # Thus rank(A26)<=205, while augmented rank>=206: fixed-g0 inconsistency.
 u27=cert.basis(27,2);a27=cert.linear_columns(alg,f,u27)
 rows27=sorted(set(b5).union(*(set(c) for c in a27)))
 assert [sum(i==j for i,e in u27) for j in range(4)]==[40,62,80,54]
 assert len(u27)==witness["positive_F4_unknowns"]==236
 assert len(rows27)==witness["positive_F4_rows"]
 positive_det=cert.minor(a27,b5,rows27,witness["positive_rank_minor"],p)
 assert len(witness["positive_rank_minor"]["rows"])==235
 assert cert.product(a27,cert.vector(u27,K))==b5
 assert not cert.product(a27,cert.vector(u27,kernel))
 assert [len(poly) for poly in K]==[0,12,59,3]
 degrees=[max(map(sum,poly),default=-1) for poly in K]
 assert degrees==[-1,19,27,15]
 # Two compact global formulas independently identify positive components.
 X,Qpos,P,Z=variables;one={(0,0,0,0):Q(1)}
 t=mul(X,Qpos)
 B=add(scale(mul(power(X,2),P),3),scale(mul(add(one,scale(t,-3)),Z),2))
 A0=add(Qpos,scale(mul(X,power(Qpos,2)),3),scale(mul(X,B),Q(-1,3)))
 assert K[1]==scale(mul(mul(power(X,8),power(r,2)),add(one,mul(X,A0))),Q(49,64))
 assert K[3]==scale(mul(power(X,9),power(r,2)),Q(-49,16))
 # Full F4 family is K+w*kernel; h7 is affine in w, not sampled.
 b7=cert.lower_rhs(alg,[f,g,K],3)
 bw={}
 for pi,(i,j) in enumerate(cert.PAIRS):
  value=add(alg.bidifferential(g[i],kernel[j],1),alg.bidifferential(kernel[i],g[j],1),
   scale(add(alg.bidifferential(f[i],kernel[j],3),alg.bidifferential(kernel[i],f[j],3)),Q(1,24)))
  for e,c in value.items():bw[(pi,e)]=-c
 assert not bw
 u6=cert.basis(27,3);a6=cert.linear_columns(alg,f,u6)
 rows6=sorted(set(b7).union(*(set(c) for c in a6)))
 assert [sum(i==j for i,e in u6) for j in range(4)]==[13,24,34,20]
 assert len(u6)==witness["h7_F6_unknowns"]==91
 assert len(rows6)==witness["h7_rows"]==974
 assert [len(b7),len(bw)]==witness["h7_rhs_nonzeros"]==[114,0]
 h7det=cert.minor(a6,b7,rows6,witness["h7_rank_minor"],p)
 assert len(witness["h7_rank_minor"]["rows"])==91
 functional={}
 for index,key,c in witness["h7_functional"]:
  key=(key[0],tuple(key[1]));assert rows6[index]==key and key not in functional
  functional[key]=Q(c)
 def apply(poly):return sum((q*poly.get(key,Q(0)) for key,q in functional.items()),Q(0))
 assert all(not apply(column) for column in a6)
 assert apply(b7)==1 and len(functional)==72
 # This rules out all 91 F6 coefficients and every w for the fixed g0.
 def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
 files=["MINIMAL_EVEN_LONG_JET_AND_FIXED_LOWER_OBSTRUCTION.md",
        "check_minimal_even_long_jet.py","minimal_even_long_jet_witness.json",
        "graded_long_cap23_witness.json","check_graded_long_cap23.py",
        "check_long_graph_curvature.py","GRADED_LONG_CAP26_EVEN_JET_OBSTRUCTION.md",
        "check_graded_long_cap26_even.py","graded_long_cap26_even_witness.json"]
 hashes={name:sha(HERE/name) for name in files}
 hashes["../dixmier-audit/verify_long_pbw.py"]=sha(source)
 receipt={"claim":"Sharp even graded h5 spatial minimum27; actual positive jet; fixed-g0 entire F4 family fails h7 at cap27",
  "status":"PASS","environment":{"python":platform.python_version(),"arithmetic":"Python integers and Fraction; no algebraic/nullspace solver"},
  "prime":p,"fixed_F2_coefficients":59,
  "fixed_cap26_augmented_minor_size":206,"fixed_cap26_augmented_minor_det":lowdet,
  "positive_cap":27,"positive_F4_unknowns":236,"positive_F4_rank_Q":235,
  "positive_minor_size":235,"positive_minor_det":positive_det,
  "positive_F4_component_degrees":degrees,"positive_F4_component_terms":[0,12,59,3],
  "full_fixed_F4_family":"K+w*(0,0,-9R^8,0)",
  "universal_h5_minimum":27,
  "universal_lower_bound_source":"GRADED_LONG_CAP26_EVEN_JET_OBSTRUCTION.md; imported theorem, not replayed by this checker",
  "h7_fixed_lower_scope":{"F6_unknowns":91,"F6_rank_Q":91,"minor_size":91,"minor_det":h7det,
   "rhs_nonzeros":[114,0],"free_F4_parameter_retained":True,
   "functional_support":72,"functional_all91columns":0,"functional_constant_rhs":1},
  "hashes":hashes,
  "coverage_boundary":"Positive cap27 h5 point, complete F4 family for fixed retained g0, fixed-g0 cap26 inconsistency and cap27 h7 obstruction. Universal h5 minimum imports the separately pinned full cap26 theorem. No universal cap27 h7 or full-lift claim.",
  "independent_replay":False,"formal_verification":False,"external_review":False}
 output=args.output_dir or Path(tempfile.mkdtemp(prefix="minimal-even-long-jet-"))
 output.mkdir(parents=True,exist_ok=True)
 target=output/"minimal_even_long_jet_controls.json"
 target.write_text(json.dumps(receipt,indent=2)+"\n")
 print("PASS: actual minimum-cap27 h5 jet; complete fixed-g0 F4 family; h7 obstruction for every w")
 print("Minors:",lowdet,positive_det,h7det,"mod",p,"functional72")
 print("Receipt:",target)
if __name__=="__main__":main()
