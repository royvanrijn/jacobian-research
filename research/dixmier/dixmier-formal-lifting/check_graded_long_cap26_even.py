#!/usr/bin/env python3
"""Replay the complete even graded Long cap26 h3/h5 obstruction."""
from fractions import Fraction as Q
from pathlib import Path
from math import isqrt
import argparse,hashlib,importlib.util,itertools,json,platform,tempfile
HERE=Path(__file__).resolve().parent
def load(path,name):
 spec=importlib.util.spec_from_file_location(name,path)
 mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument("--output-dir",type=Path);args=parser.parse_args()
 cert=load(HERE/"check_graded_long_cap23.py","complete_coefficient_engine")
 alg=load(HERE/"check_long_graph_curvature.py","exact_polynomial_arithmetic")
 eng=load(HERE/"check_graded_long_cap23_full_jet.py","parameter_coefficient_engine")
 source=HERE.parent/"dixmier-audit/verify_long_pbw.py"
 producer=load(source,"canonical_Long_input")
 f=[{tuple(e):Q(c) for e,c in poly.items()} for poly in producer.long_images()]
 witness_path=HERE/"graded_long_cap26_even_witness.json";w=json.loads(witness_path.read_text())
 assert w["cap"]==26 and w["parity"]=="even";first=w["first_family"]
 p=1000003;assert all(p%d for d in range(2,isqrt(p)+1))
 for poly,weight in zip(f,cert.WEIGHTS):
  assert all(-e[0]+e[1]+4*e[2]+2*e[3]==weight for e in poly)
 def parse(rows):
  out=[{} for _ in range(4)]
  for name,e,c in rows:
   i=cert.LABELS.index(name);e=tuple(e);assert e not in out[i];out[i][e]=Q(c)
  return out
 g0=parse(first["particular"]);kernels=[parse(z) for z in first["kernels"]]
 prior_path=HERE/"graded_long_cap23_witness.json";prior=json.loads(prior_path.read_text())
 assert g0==parse(prior["first_particular"]) and len(first["particular"])==59
 u2=cert.basis(26,1);a2=cert.linear_columns(alg,f,u2);b3=cert.lower_rhs(alg,[f],1)
 assert [sum(i==j for i,e in u2) for j in range(4)]==[80,114,141,102]
 assert len(u2)==437 and cert.product(a2,cert.vector(u2,g0))==b3
 assert len(kernels)==6 and all(not g[0] for g in [g0]+kernels)
 vectors=[cert.vector(u2,k) for k in kernels]
 for k,v,Grows in zip(kernels,vectors,first["hamiltonians"]):
  assert not cert.product(a2,v)
  G={tuple(e):Q(c) for e,c in Grows}
  assert all(-e[0]+e[1]+4*e[2]+2*e[3]==-3 for e in G)
  assert [alg.bidifferential(G,fi,1) for fi in f]==k
 free=first["free_columns"];assert len(free)==6 and len(set(free))==6
 assert [[v[j] for j in free] for v in vectors]==[
        [Q(int(i==j)) for j in range(6)] for i in range(6)]
 linear_rows2=sorted(set().union(*(set(c) for c in a2)))
 det2=cert.minor(a2,{},linear_rows2,w["F2_rank_minor"],p)
 assert len(w["F2_rank_minor"]["rows"])==431
 # A compact alternative basis, verified against the retained normalized basis.
 add,mul,scale=alg.add,alg.multiply,alg.scale
 power=lambda poly,n:cert.power(alg,poly,n)
 X,Qpos,P,Z=[{tuple(int(i==j) for i in range(4)):Q(1)} for j in range(4)]
 one={(0,0,0,0):Q(1)};r,T,D,S=f;t=mul(X,Qpos)
 B=add(scale(mul(power(X,2),P),3),scale(mul(add(one,scale(t,-3)),Z),2))
 A0=add(Qpos,scale(mul(X,power(Qpos,2)),3),scale(mul(X,B),Q(-1,3)))
 Y=add(mul(r,A0),scale(t,-3))
 J4=add(scale(power(Y,4),9),scale(power(Y,3),66),
         scale(mul(t,power(Y,2)),36),scale(mul(t,Y),150),
         scale(power(t,2),171),scale(t,-380))
 compact_G=[power(r,3),mul(power(r,4),S),power(X,3),
            mul(power(r,3),add(Y,scale(one,2))),mul(power(X,4),Qpos),
            mul(power(r,3),J4)]
 compact_jets=[[alg.bidifferential(G,z,1) for z in f] for G in compact_G]
 compact_vectors=[cert.vector(u2,z) for z in compact_jets]
 transition=[[v[j] for j in free] for v in compact_vectors]
 for compact_v,coordinates in zip(compact_vectors,transition):
  assert [sum((c*v[j] for c,v in zip(coordinates,vectors)),Q(0))
          for j in range(len(u2))]==compact_v
 detmatrix=[row[:] for row in transition];detcompact=Q(1)
 for j in range(6):
  pivot=next(i for i in range(j,6) if detmatrix[i][j])
  if pivot!=j:detmatrix[j],detmatrix[pivot]=detmatrix[pivot],detmatrix[j];detcompact=-detcompact
  q=detmatrix[j][j];detcompact*=q
  for i in range(j+1,6):
   factor=detmatrix[i][j]/q
   for k in range(j+1,6):detmatrix[i][k]-=factor*detmatrix[j][k]
 assert detcompact
 # Nonzero rank431 minor and six independent exact kernels prove completeness.
 zero=(0,)*6;F0={zero:f};F2={zero:g0}
 for j,k in enumerate(kernels):F2[tuple(int(i==j) for i in range(6))]=k
 rhs=eng.rhs_coefficients(cert,alg,[F0,{},F2],5)
 monos=sorted(rhs,key=lambda m:(sum(m),m))
 assert [list(m) for m in monos]==w["h5"]["monomials"] and len(monos)==20
 u4=cert.basis(26,2);a4=cert.linear_columns(alg,f,u4)
 assert [sum(i==j for i,e in u4) for j in range(4)]==[34,54,71,47]
 assert len(u4)==206
 rows4=eng.rows_for(a4,rhs);assert len(rows4)==1590
 det4=cert.minor(a4,{},sorted(set().union(*(set(c) for c in a4))),
                 w["F4_rank_minor"],p)
 KR9=[{},{},alg.scale(cert.power(alg,f[0],8),-9),{}]
 assert not cert.product(a4,cert.vector(u4,KR9))
 assert len(w["F4_rank_minor"]["rows"])==205
 functionals=[eng.functional(w["h5"][name],rows4) for name in
              ("functional_p0","functional_p1")]
 alpha=Q(32,3969);k3k4=(0,0,0,1,1,0);k4sq=(0,0,0,0,2,0)
 expected=[{zero:Q(1),k3k4:-alpha},{k4sq:Q(1)}]
 for left,want in zip(functionals,expected):
  assert all(not eng.apply(left,column) for column in a4)
  assert {m:eng.apply(left,b) for m,b in rhs.items() if eng.apply(left,b)}==want
 def pmul(a,b):
  out={}
  for m,c in a.items():
   for n,d in b.items():
    key=eng.plus(m,n);out[key]=out.get(key,Q(0))+c*d
  return {m:c for m,c in out.items() if c}
 # Exact polynomial Bezout certificate: no radical or root-field assumption.
 part0=pmul({zero:Q(1),k3k4:alpha},expected[0])
 part1=pmul({(0,0,0,2,0,0):alpha**2},expected[1])
 for m,c in part1.items():part0[m]=part0.get(m,Q(0))+c
 assert {m:c for m,c in part0.items() if c}=={zero:Q(1)}
 # An ACTUAL polynomial solution of the same h5 equation beyond cap26.
 # Regenerate the weighted de Rham homotopy, rather than importing a
 # previously computed positive-control tuple.
 E=[[{} for j in range(4)] for i in range(4)]
 for (pi,e),c in rhs[zero].items():
  i,j=cert.PAIRS[pi];E[i][j][e]=c;E[j][i][e]=-c
 for i,j,k in itertools.combinations(range(4),3):
  assert not add(alg.bidifferential(f[i],E[j][k],1),
                 alg.bidifferential(f[j],E[k][i],1),
                 alg.bidifferential(f[k],E[i][j],1))
 contraction=[scale(D,-4),scale(S,-1),scale(r,-1),scale(T,2)]
 positive=[scale(add(*(mul(contraction[j],E[j][i]) for j in range(4))),Q(-1,9))
           for i in range(4)]
 for i,j in cert.PAIRS:
  assert add(alg.bidifferential(positive[i],f[j],1),
             alg.bidifferential(f[i],positive[j],1))==E[i][j]
 positive_degrees=[max(map(sum,z),default=-1) for z in positive]
 assert positive_degrees==[19,31,39,27]
 for poly,weight in zip(positive,cert.WEIGHTS):
  assert all(-e[0]+e[1]+4*e[2]+2*e[3]==weight-12 for e in poly)
 positive_h7=cert.lower_rhs(alg,[f,g0,positive],3)
 assert len(positive_h7)==434 and max(sum(e) for pi,e in positive_h7)==48
 # Deliberate corruptions must fail actual identities.
 bad=[dict(poly) for poly in g0];e=next(iter(bad[3]));bad[3][e]+=1
 assert cert.product(a2,cert.vector(u2,bad))!=b3
 badleft=dict(functionals[0]);key=next(iter(rhs[zero]))
 badleft[key]=badleft.get(key,Q(0))+1
 assert eng.apply(badleft,rhs[zero])!=1
 def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
 receipt={
  "claim":"No even graded canonical Long scalar full Weyl jet of common spatial degree26 through h5",
  "status":"PASS","environment":{"python":platform.python_version(),
    "arithmetic":"Fraction polynomial identities, bare modular determinants and polynomial Bezout identity"},
  "cap":26,"parity":"even","image_order":list(cert.LABELS),
  "h3":{"unknowns":437,"component_counts":[80,114,141,102],
    "rank_Q":431,"minor_det_modp":det2,"independent_exact_kernels":6,
    "polynomial_Hamiltonians_verified":True,
    "compact_basis":["R^3","R^4*S","X^3","R^3*(Y+2)","X^4*Q","R^3*J4"],
    "compact_basis_transition_det_Q":str(detcompact),
    "all_R2_corrections_zero_as_consequence":True,
    "first_particular_nonzero_coefficients":59},
  "h5":{"unknowns":206,"component_counts":[34,54,71,47],"rows":1590,
    "rank_Q":205,"minor_det_modp":det4,"exact_kernel":"K_R9",
    "all_first_parameters_retained":["k0","k1","k2","k3","k4","k5"],
    "RHS_parameter_coefficients":20,
    "functional_supports":[len(z) for z in functionals],
    "identities":["ell0*A=0; ell0*B=1-32*k3*k4/3969",
                  "ell1*A=0; ell1*B=k4^2"],
    "Bezout":"1=(1+32*k3*k4/3969)*p0+(32/3969)^2*k3^2*p1",
    "unit_ideal_verified_as_polynomial_identity":True,
    "all_206_columns_and_all_20_RHS_coefficients_verified":True},
  "positive_control_beyond_bound":{"construction":"Weighted de Rham homotopy of the closed h5 defect, with all first parameters zero",
    "F4_degrees":positive_degrees,"all_six_equations_verified":True,
    "truncated_h7_defect_nonzero":True,"truncated_h7_defect_terms":len(positive_h7),
    "truncated_h7_defect_max_degree":48,
    "scope":"Actual finite polynomial even jet through h5, common cap39; no all-order or bounded-lift conclusion"},
  "invalid_controls":{"corrupted_first_particular_rejected":True,
                      "corrupted_obstruction_functional_rejected":True},
  "prime":p,"hashes":{
   "check_graded_long_cap26_even.py":sha(Path(__file__)),
   "GRADED_LONG_CAP26_EVEN_JET_OBSTRUCTION.md":sha(HERE/"GRADED_LONG_CAP26_EVEN_JET_OBSTRUCTION.md"),
   "graded_long_cap26_even_witness.json":sha(witness_path),
   "graded_long_cap23_witness.json":sha(prior_path),
   "check_graded_long_cap23.py":sha(HERE/"check_graded_long_cap23.py"),
   "check_graded_long_cap23_full_jet.py":sha(HERE/"check_graded_long_cap23_full_jet.py"),
   "check_long_graph_curvature.py":sha(HERE/"check_long_graph_curvature.py"),
   "../dixmier-audit/verify_long_pbw.py":sha(source)},
  "coverage_boundary":"Complete EVEN graded cap26 family through h5; no kernel or parameter omitted, no gauge normalization used, all complex solutions excluded by the unit ideal. No odd-jet, nongraded, cap27, universal bounded-lift or DC1/DC2 conclusion.",
  "independent_replay":False,"formal_verification":False,"external_review":False}
 out=args.output_dir or Path(tempfile.mkdtemp(prefix="graded-long-cap26-"));out.mkdir(parents=True,exist_ok=True)
 dest=out/"graded_long_cap26_even_controls.json";dest.write_text(json.dumps(receipt,indent=2)+"\n")
 print("PASS: complete F2 rank431/437; all six parameters retained; h5 unit ideal from two exact functionals")
 print("Functional supports:",*[len(z) for z in functionals],"minors:",det2,det4,"receipt:",dest)
if __name__=="__main__":main()
