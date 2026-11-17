#!/usr/bin/env python3
"""Complete even graded Long cap27 h3/h5/h7 certificate replay."""
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
 cert=load(HERE/"check_graded_long_cap23.py","complete_even_columns")
 alg=load(HERE/"check_long_graph_curvature.py","independent_Fraction_polynomials")
 eng=load(HERE/"check_graded_long_cap23_full_jet.py","general_parameter_and_PBW_engine")
 source=HERE.parent/"dixmier-audit/verify_long_pbw.py"
 producer=load(source,"canonical_global_Long_input")
 f=[{tuple(e):Q(c) for e,c in poly.items()} for poly in producer.long_images()]
 witness_path=HERE/"graded_long_cap27_even_witness.json";w=json.loads(witness_path.read_text())
 assert w["cap"]==27 and w["parity"]=="even";first=w["first_family"]
 p=1000003;assert all(p%d for d in range(2,isqrt(p)+1))
 def parse(rows):
  out=[{} for _ in range(4)]
  for name,e,c in rows:
   i=cert.LABELS.index(name);e=tuple(e);assert e not in out[i];out[i][e]=Q(c)
  return out
 g0=parse(first["particular"]);kernels=[parse(z) for z in first["kernels"]]
 prior_path=HERE/"graded_long_cap23_witness.json";prior=json.loads(prior_path.read_text())
 assert g0==parse(prior["first_particular"]) and len(first["particular"])==59
 for poly,weight in zip(f,cert.WEIGHTS):
  assert all(-e[0]+e[1]+4*e[2]+2*e[3]==weight for e in poly)
 u2=cert.basis(27,1);a2=cert.linear_columns(alg,f,u2);b3=cert.lower_rhs(alg,[f],1)
 assert [sum(i==j for i,e in u2) for j in range(4)]==[91,127,156,114]
 assert len(u2)==488 and cert.product(a2,cert.vector(u2,g0))==b3
 vectors=[cert.vector(u2,k) for k in kernels];assert len(kernels)==8
 for k,v,Grows in zip(kernels,vectors,first["hamiltonians"]):
  assert not cert.product(a2,v)
  G={tuple(e):Q(c) for e,c in Grows}
  assert all(-e[0]+e[1]+4*e[2]+2*e[3]==-3 for e in G)
  assert [alg.bidifferential(G,fi,1) for fi in f]==k
 free=first["free_columns"]
 assert len(free)==8 and len(set(free))==8
 assert [[v[j] for j in free] for v in vectors]==[
        [Q(int(i==j)) for j in range(8)] for i in range(8)]
 det2=cert.minor(a2,{},sorted(set().union(*(set(c) for c in a2))),
                 w["F2_rank_minor"],p)
 assert len(w["F2_rank_minor"]["rows"])==480
 # Compact alternative Hamiltonians, including the new first-R direction.
 add,mul,scale=alg.add,alg.multiply,alg.scale;power=lambda poly,n:cert.power(alg,poly,n)
 X,Qpos,P,Z=[{tuple(int(i==j) for i in range(4)):Q(1)} for j in range(4)]
 one={(0,0,0,0):Q(1)};r,T,D,S=f;t=mul(X,Qpos)
 B=add(scale(mul(power(X,2),P),3),scale(mul(add(one,scale(t,-3)),Z),2))
 A0=add(Qpos,scale(mul(X,power(Qpos,2)),3),scale(mul(X,B),Q(-1,3)))
 Y=add(mul(r,A0),scale(t,-3))
 J4=add(scale(power(Y,4),9),scale(power(Y,3),66),scale(mul(t,power(Y,2)),36),
        scale(mul(t,Y),150),scale(power(t,2),171),scale(t,-380))
 compact_G=[power(r,3),mul(power(r,4),S),power(X,3),mul(power(r,3),add(Y,scale(one,2))),
            mul(power(X,4),Qpos),mul(power(r,3),J4),mul(power(r,5),T),mul(power(X,5),Z)]
 compact_jets=[[alg.bidifferential(G,z,1) for z in f] for G in compact_G]
 compact_vectors=[cert.vector(u2,z) for z in compact_jets]
 transition=[[v[j] for j in free] for v in compact_vectors]
 for compact_v,coordinates in zip(compact_vectors,transition):
  assert [sum((c*v[j] for c,v in zip(coordinates,vectors)),Q(0))
          for j in range(len(u2))]==compact_v
 dm=[row[:] for row in transition];detcompact=Q(1)
 for j in range(8):
  pivot=next(i for i in range(j,8) if dm[i][j])
  if pivot!=j:dm[j],dm[pivot]=dm[pivot],dm[j];detcompact=-detcompact
  q=dm[j][j];detcompact*=q
  for i in range(j+1,8):
   factor=dm[i][j]/q
   for k in range(j+1,8):dm[i][k]-=factor*dm[j][k]
 assert detcompact
 assert compact_jets[-1][0]==scale(power(X,7),-3)
 assert any(k[0] for k in kernels) # first image was NOT held fixed.
 # All eight first parameters are retained through the h5 equation.
 zero8=(0,)*8;F0={zero8:f};F2={zero8:g0}
 for j,k in enumerate(kernels):F2[tuple(int(i==j) for i in range(8))]=k
 rhs5=eng.rhs_coefficients(cert,alg,[F0,{},F2],5)
 mons5=sorted(rhs5,key=lambda m:(sum(m),m))
 assert [list(m) for m in mons5]==w["h5_all"]["monomials"] and len(mons5)==36
 u4=cert.basis(27,2);a4=cert.linear_columns(alg,f,u4)
 assert [sum(i==j for i,e in u4) for j in range(4)]==[40,62,80,54]
 assert len(u4)==236;rows5=eng.rows_for(a4,rhs5);assert len(rows5)==1722
 det4=cert.minor(a4,{},sorted(set().union(*(set(c) for c in a4))),
                 w["F4_rank_minor"],p)
 KR9=[{},{},scale(power(r,8),-9),{}]
 assert not cert.product(a4,cert.vector(u4,KR9))
 assert len(w["F4_rank_minor"]["rows"])==235
 expected_lower={}
 for index in (5,6,7):
  expected_lower[f"k{index}_square"]={tuple(2*int(j==index) for j in range(8)):Q(1)}
 expected_lower["linear"]={tuple(int(j==i) for j in range(8)):co for i,co in
      [(2,Q(1)),(4,Q(10,3)),(5,Q(3)),(6,Q(9)),(7,Q(-1340,27))]}
 lower_supports={}
 for name,want in expected_lower.items():
  left=eng.functional(w["h5_all"]["functionals"][name],rows5)
  assert all(not eng.apply(left,col) for col in a4)
  assert {m:eng.apply(left,b) for m,b in rhs5.items() if eng.apply(left,b)}==want
  lower_supports[name]=len(left)
 # Over fields the three individual squares force k5=k6=k7=0.
 # The linear identity then forces k2=-10*k4/3, without real/imag splitting.
 zero4=(0,)*4
 dirB=[add(x,scale(y,Q(-10,3))) for x,y in zip(kernels[4],kernels[2])]
 reduced_F0={zero4:f}
 reduced_F2={zero4:g0,(1,0,0,0):kernels[0],(0,1,0,0):kernels[1],
             (0,0,1,0):kernels[3],(0,0,0,1):dirB}
 reduced_rhs=eng.rhs_coefficients(cert,alg,[reduced_F0,{},reduced_F2],5)
 reduced_mons=sorted(reduced_rhs,key=lambda m:(sum(m),m))
 assert [list(m) for m in reduced_mons]==w["h5_reduced"]["monomials"]
 reduced_rows=eng.rows_for(a4,reduced_rhs)
 product_left=eng.functional(w["h5_reduced"]["product_functional"],reduced_rows)
 assert all(not eng.apply(product_left,col) for col in a4)
 assert {m:eng.apply(product_left,b) for m,b in reduced_rhs.items()
         if eng.apply(product_left,b)}=={(0,0,1,1):Q(1)}
 # Thus k3*k4=0 in a FIELD, covered by the two displayed components.
 u6=cert.basis(27,3);a6=cert.linear_columns(alg,f,u6)
 assert [sum(i==j for i,e in u6) for j in range(4)]==[13,24,34,20]
 assert len(u6)==91
 det6=cert.minor(a6,{},sorted(set().union(*(set(c) for c in a6))),
                 w["F6_rank_minor"],p)
 assert len(w["F6_rank_minor"]["rows"])==91
 branch_controls={};all_rhs={}
 for name,direction in (("A",kernels[3]),("B",dirB)):
  B0={zero4:f};B2={zero4:g0,(1,0,0,0):kernels[0],(0,1,0,0):kernels[1],(0,0,1,0):direction}
  branch=w["F4_branch_families"][name]
  B4={tuple(item["monomial"]):parse(item["coefficients"]) for item in branch["F4"]}
  assert B4[(0,0,0,1)]==KR9
  known5=eng.rhs_coefficients(cert,alg,[B0,{},B2],5)
  assert set(B4)==set(known5)|{(0,0,0,1)}
  for m,tup in B4.items():
   assert not tup[0]
   assert cert.product(a4,cert.vector(u4,tup))==known5.get(m,{})
  # rank235/236 + KR9 makes these particular families exhaustive.
  known7=eng.rhs_coefficients(cert,alg,[B0,{},B2,{},B4],7)
  ms=sorted(known7,key=lambda m:(sum(m),m))
  assert [list(m) for m in ms]==w["h7"][name]["monomials"]
  rows7=eng.rows_for(a6,known7);assert len(rows7)==974
  left=eng.functional(w["h7"][name]["functional"],rows7)
  assert all(not eng.apply(left,col) for col in a6)
  assert {m:eng.apply(left,b) for m,b in known7.items() if eng.apply(left,b)}=={zero4:Q(1)}
  assert all(m[0]==0 and m[3]==0 for m in known7) # u,w cancel, not removed.
  # Deliberate bad functional control on each component.
  badleft=dict(left);key=next(iter(known7[zero4]));badleft[key]=badleft.get(key,Q(0))+1
  assert eng.apply(badleft,known7[zero4])!=1
  values=(Q(2),Q(-1),Q(1),Q(3))
  sizes=eng.pbw_control(eng.evaluate_frames(alg,[B0,{},B2,{},B4],values),5)
  branch_controls[name]={"F4_monomials":[list(m) for m in B4],
    "F4_nonzero_coefficients":sum(sum(map(len,tup)) for tup in B4.values()),
    "all_R4_components_zero":True,"all_h5_coefficients_verified":True,
    "h7_functional_support":len(left),"h7_RHS_monomials":[list(m) for m in ms],
    "all_91_F6_columns_annihilated":True,"h7_constant_evaluation":"1",
    "all_nonconstant_evaluations":"0","u_w_cancel_without_normalization":True,
    "normal_PBW_positive_h5_point":{"parameters_u_v_t_w":[str(z) for z in values],
                                   "all_six_relations":True,"operator_terms":sizes}}
  all_rhs[name]=known7
 # Retained shared final functional, redundant with the two component proofs.
 shared_rhs={("constant",zero4):all_rhs["A"][zero4]};assert all_rhs["A"][zero4]==all_rhs["B"][zero4]
 for name,rr in all_rhs.items():
  for m,b in rr.items():
   if m!=zero4:shared_rhs[(name,m)]=b
 keys=list(shared_rhs)
 assert [[name,list(m)] for name,m in keys]==w["h7_common"]["RHS_key_order"]
 shared_rows=eng.rows_for(a6,{(i,):b for i,b in enumerate(shared_rhs.values())})
 shared_left=eng.functional(w["h7_common"]["functional"],shared_rows)
 assert all(not eng.apply(shared_left,col) for col in a6)
 assert [eng.apply(shared_left,b) for b in shared_rhs.values()]==[Q(1)]+[Q(0)]*(len(shared_rhs)-1)
 # Corrupt one retained particular and reject the exact identity.
 bad=parse(w["F4_branch_families"]["A"]["F4"][0]["coefficients"])
 e=next(iter(bad[2]));bad[2][e]+=1
 base5=eng.rhs_coefficients(cert,alg,[{zero4:f},{},{zero4:g0}],5)[zero4]
 assert cert.product(a4,cert.vector(u4,bad))!=base5
 def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
 receipt={"claim":"No even graded canonical Long scalar full Weyl jet of common degree27 through h7",
  "status":"PASS","environment":{"python":platform.python_version(),
    "arithmetic":"Fraction polynomial identities, bare modular determinants and separate normal PBW controls"},
  "cap":27,"parity":"even","h3":{"unknowns":488,"rank_Q":480,
    "exact_independent_kernels":8,"minor_det_modp":det2,
    "compact_basis_transition_det_Q":str(detcompact),"first_R2_deformations_included":True},
  "h5":{"unknowns":236,"rank_Q":235,"minor_det_modp":det4,"rows":1722,
    "all_eight_first_parameters_retained":True,"RHS_coefficients":36,
    "necessity_functional_supports":lower_supports,
    "reduced_product_support":len(product_left),
    "field_classification":"k5=k6=k7=0;k2=-10*k4/3;k3*k4=0",
    "scope":"Field-valued points, not an asserted nonreduced ideal equality",
    "both_entire_F4_families_verified":True,"full_kernel":"K_R9"},
  "h7":{"unknowns":91,"component_counts":[13,24,34,20],"rank_Q":91,"minor_det_modp":det6,
    "branches":branch_controls,"common_functional_support":len(shared_left),
    "conclusion":"0=1 on both complete components for all retained parameters"},
  "invalid_controls":{"corrupted_F4_particular_rejected":True,
                       "corrupted_branch_h7_functionals_rejected":True},
  "prime":p,"hashes":{
   "GRADED_LONG_CAP27_EVEN_JET_OBSTRUCTION.md":sha(HERE/"GRADED_LONG_CAP27_EVEN_JET_OBSTRUCTION.md"),
   "check_graded_long_cap27_even.py":sha(Path(__file__)),
   "graded_long_cap27_even_witness.json":sha(witness_path),
   "graded_long_cap23_witness.json":sha(prior_path),
   "check_graded_long_cap23.py":sha(HERE/"check_graded_long_cap23.py"),
   "check_graded_long_cap23_full_jet.py":sha(HERE/"check_graded_long_cap23_full_jet.py"),
   "check_long_graph_curvature.py":sha(HERE/"check_long_graph_curvature.py"),
   "../dixmier-audit/verify_long_pbw.py":sha(source)},
  "coverage_boundary":"Complete even graded scalar cap27 jet through h7; both field-valued h5 components and all later parameters retained; all first-image coefficients allowed. Positive h5 point included. No odd/nongraded/cap28/universal bounded-lift/DC1/DC2 theorem.",
  "independent_replay":False,"formal_verification":False,"external_review":False}
 out=args.output_dir or Path(tempfile.mkdtemp(prefix="graded-long-cap27-"));out.mkdir(parents=True,exist_ok=True)
 dest=out/"graded_long_cap27_even_controls.json";dest.write_text(json.dumps(receipt,indent=2)+"\n")
 print("PASS: complete cap27 even family; two h5 field components; both give h7 contradiction")
 print("h3/h5/h7 minors:",det2,det4,det6,"receipt:",dest)
if __name__=="__main__":main()
