#!/usr/bin/env python3
"""Full cap24 graded Long obstruction, retaining both quadratic-field branches."""
from fractions import Fraction as Q
from pathlib import Path
import argparse,hashlib,importlib.util,json,platform,tempfile
HERE=Path(__file__).resolve().parent
def load(path,name):
 spec=importlib.util.spec_from_file_location(name,path)
 mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
def basis(weights,m):
 out=[]
 for i,wi in enumerate(weights):
  for b in range(25):
   for c in range(25-b):
    for d in range(25-b-c):
     a=b+4*c+2*d+3*m-wi
     if a>=0 and a+b+c+d<=24:out.append((i,(a,b,c,d)))
 return sorted(out)
def parse_pairs(cert,rows):
 out=[({}, {}) for _ in range(4)]
 for label,e,co in rows:
  i=cert.LABELS.index(label);e=tuple(e)
  for component in (0,1):
   value=Q(co[component])
   if value:
    assert e not in out[i][component]
    out[i][component][e]=value
 return out
def pair_bidifferential(alg,a,b,n):
 real=alg.bidifferential(a[0],b[0],n) if a[0] and b[0] else {}
 imag={}
 if a[1] and b[1]:real=alg.add(real,alg.scale(alg.bidifferential(a[1],b[1],n),-2))
 if a[0] and b[1]:imag=alg.add(imag,alg.bidifferential(a[0],b[1],n))
 if a[1] and b[0]:imag=alg.add(imag,alg.bidifferential(a[1],b[0],n))
 return real,imag
def field_rhs(cert,alg,eng,frames,order):
 result={}
 for pi,(i,j) in enumerate(cert.PAIRS):
  for a,b,ell,factor in eng.rhs_indices(order,len(frames)):
   for ma,fa in frames[a].items():
    for mb,fb in frames[b].items():
     mon=eng.plus(ma,mb)
     value=pair_bidifferential(alg,fa[i],fb[j],ell)
     dest=result.setdefault(mon,({},{}))
     for component in (0,1):
      for e,c in value[component].items():
       key=(pi,e);dest[component][key]=dest[component].get(key,Q(0))-factor*c
 out={}
 for mon,(real,imag) in result.items():
  real={e:c for e,c in real.items() if c}
  imag={e:c for e,c in imag.items() if c}
  if real or imag:out[mon]=(real,imag)
 return out
def field_rows(columns,rhs):
 return sorted(set().union(*(set(c) for c in columns),
                           *(set(z) for pair in rhs.values() for z in pair)))
def conjugate(frames):
 return [{m:[(real,{e:-c for e,c in imag.items()}) for real,imag in tup]
          for m,tup in frame.items()} for frame in frames]

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument("--output-dir",type=Path);args=parser.parse_args()
 cert=load(HERE/"check_graded_long_cap23.py","complete_column_engine")
 eng=load(HERE/"check_graded_long_cap23_full_jet.py","general_coefficient_and_PBW_engine")
 alg=load(HERE/"check_long_graph_curvature.py","independent_polynomial_engine")
 source=HERE.parent/"dixmier-audit/verify_long_pbw.py"
 prod=load(source,"canonical_global_input")
 images=[{tuple(e):Q(c) for e,c in poly.items()} for poly in prod.long_images()]
 witness_path=HERE/"graded_long_cap24_full_jet_witness.json"
 witness=json.loads(witness_path.read_text());assert witness["cap"]==24
 lower_path=HERE/witness["complex_prior_witness"]
 lower=json.loads(lower_path.read_text())
 assert lower["cap"]==24 and len(lower["particulars"])==11
 add,mul,scale=alg.add,alg.multiply,alg.scale;power=lambda p,n:cert.power(alg,p,n)
 X,Qpos,P,Z=[{tuple(int(i==j) for i in range(4)):Q(1)} for j in range(4)]
 one={(0,0,0,0):Q(1)};r,T,D,S=images;t=mul(X,Qpos)
 B=add(scale(mul(power(X,2),P),3),scale(mul(add(one,scale(t,-3)),Z),2))
 A0=add(Qpos,scale(mul(X,power(Qpos,2)),3),scale(mul(X,B),Q(-1,3)))
 Y=add(mul(r,A0),scale(t,-3))
 H4=add(power(Y,4),scale(mul(t,power(Y,2)),4),
         scale(mul(t,Y),Q(-16,3)),scale(power(t,2),10))
 H5old=add(scale(power(Y,5),9),scale(mul(t,power(Y,3)),45),
           scale(mul(power(t,2),Y),45),scale(mul(t,power(Y,2)),-60),
           scale(mul(t,Y),110),scale(power(t,2),-720))
 H6=add(mul(power(r,4),D),scale(H5old,Q(1,25)))
 H5new=add(scale(power(Y,5),9),scale(mul(t,power(Y,3)),45),
           scale(mul(t,power(Y,2)),-60),scale(mul(power(t,2),Y),135),
           scale(mul(t,Y),-10),scale(power(t,2),-555))
 jets=lambda G:[alg.bidifferential(G,f,1) for f in images]
 jy,jb,jc,jd,jf=[jets(G) for G in (Y,mul(power(r,2),T),H4,H6,H5new)]
 zero7=(0,)*7
 F0={zero7:[(poly,{}) for poly in images]};F1={}
 for index,jet in enumerate((jy,jb,jc)):
  F1[tuple(int(index==j) for j in range(7))]=[(poly,{}) for poly in jet]
 # Parameter order a,b,c,z,u,v,gamma; d=(125+100omega)z.
 F1[(0,0,0,1,0,0,0)]=[
    (add(x,scale(y,125)),scale(y,100)) for x,y in zip(jf,jd)]
 monmap={"1":zero7,"a2":(2,0,0,0,0,0,0),"ab":(1,1,0,0,0,0,0),
  "ac":(1,0,1,0,0,0,0),"b2":(0,2,0,0,0,0,0),"bc":(0,1,1,0,0,0,0),
  "c2":(0,0,2,0,0,0,0),"az":(1,0,0,1,0,0,0),"bz":(0,1,0,1,0,0,0),
  "cz":(0,0,1,1,0,0,0),"z2":(0,0,0,2,0,0,0)}
 F2={}
 for item in lower["particulars"]:
  name=item["monomial"].replace("²","2").replace("^","").replace("*","")
  F2[monmap[name]]=parse_pairs(cert,item["coefficients"])
 for index,item in enumerate(lower["homogeneous_kernels"]):
  F2[tuple(int(j==index+4) for j in range(7))]=parse_pairs(cert,item["coefficients"])
 # Verify all rational and omega components of the input h3 identities.
 # This validates identities in Q(omega), not separate constraints on parameters.
 u2=basis(cert.WEIGHTS,2);a2=cert.linear_columns(alg,images,u2)
 prior_rhs=field_rhs(cert,alg,eng,[F0,F1],3)
 for mon in set(F2)|set(prior_rhs):
  tup=F2.get(mon,[({}, {}) for _ in range(4)])
  right=prior_rhs.get(mon,({},{}))
  for component in (0,1):
   assert cert.product(a2,cert.vector(u2,[z[component] for z in tup]))==right[component]

 rhs_field=field_rhs(cert,alg,eng,[F0,F1,F2],4)
 u3=basis(cert.WEIGHTS,3);a3=cert.linear_columns(alg,images,u3)
 rows3=field_rows(a3,rhs_field)
 assert len(u3)==236 and len(rows3)==1673
 left_field=eng.functional(witness["field_h4_z3"],rows3)
 assert all(not eng.apply(left_field,column) for column in a3)
 zcube=(0,0,0,3,0,0,0)
 real_values={m:eng.apply(left_field,z[0]) for m,z in rhs_field.items()
              if eng.apply(left_field,z[0])}
 imag_values={m:eng.apply(left_field,z[1]) for m,z in rhs_field.items()
              if eng.apply(left_field,z[1])}
 assert real_values=={zcube:Q(1)} and not imag_values
 # Explicitly check conjugation of the entire generated RHS. The SAME
 # rational functional gives z^3+omega*0 on BOTH Galois branches.
 conjugated_rhs=field_rhs(cert,alg,eng,conjugate([F0,F1,F2]),4)
 assert set(conjugated_rhs)==set(rhs_field)
 for mon,(real,imag) in rhs_field.items():
  assert conjugated_rhs[mon]==(real,{e:-c for e,c in imag.items()})
 assert {m:eng.apply(left_field,z[0]) for m,z in conjugated_rhs.items()
         if eng.apply(left_field,z[0])}=={zcube:Q(1)}
 assert all(not eng.apply(left_field,z[1]) for z in conjugated_rhs.values())

 # Now z=0 is FORCED over fields. Remaining frame coefficients are real
 # rational polynomials; a,b,c,u,v,gamma may still be arbitrary complex values.
 def set_z_zero(frame):
  out={}
  for m,tup in frame.items():
   if m[3]:continue
   assert all(not imag for real,imag in tup)
   key=m[:3]+m[4:]
   eng.add_tuple(alg,out,key,[real for real,imag in tup])
  return out
 Z0,Z1,Z2=[set_z_zero(frame) for frame in (F0,F1,F2)]
 rhs_rational=eng.rhs_coefficients(cert,alg,[Z0,Z1,Z2],4)
 rows_rational=eng.rows_for(a3,rhs_rational)
 left_L=eng.functional(witness["rational_z0_h4_L3"],rows_rational)
 assert all(not eng.apply(left_L,column) for column in a3)
 expectedL3={(0,0,3,0,0,0):Q(1),(1,0,2,0,0,0):Q(-81,544),
   (2,0,1,0,0,0):Q(2187,295936),(3,0,0,0,0,0):Q(-19683,160989184)}
 assert {m:eng.apply(left_L,b) for m,b in rhs_rational.items()
         if eng.apply(left_L,b)}==expectedL3
 # L^3=0 over a field implies a=544c/27, with no nilpotent-parameter assertion.
 det3=cert.minor(a3,{},sorted(set().union(*(set(c) for c in a3))),
                 witness["F3_rank_minor"],1000003)
 KR6=[{},{},scale(power(r,5),-6),{}]
 assert not cert.product(a3,cert.vector(u3,KR6))
 assert len(witness["F3_rank_minor"]["rows"])==235
 ratio=Q(544,27)
 def normalize(frame):
  out={}
  for m,tup in frame.items():
   key=(m[1],m[0]+m[2],m[3],m[4],m[5],0)
   eng.add_tuple(alg,out,key,eng.tuple_scale(alg,tup,ratio**m[0]))
  return out
 R0,R1,R2=[normalize(frame) for frame in (Z0,Z1,Z2)]
 reduced_rhs=eng.rhs_coefficients(cert,alg,[R0,R1,R2],4)
 R3={}
 for item in witness["F3_family"]["particulars"]:
  mon=tuple(item["monomial_b_c_u_v_gamma_w"]);tup=[{} for _ in range(4)]
  for name,e,c in item["coefficients"]:
   i=cert.LABELS.index(name);e=tuple(e);assert e not in tup[i];tup[i][e]=Q(c)
  assert not tup[0]
  assert cert.product(a3,cert.vector(u3,tup))==reduced_rhs[mon]
  R3[mon]=tup
 assert set(R3)==set(reduced_rhs) and len(R3)==7
 R3[(0,0,0,0,0,1)]=KR6
 assert sum(sum(map(len,tup)) for m,tup in R3.items() if m[-1]==0)==253

 rhs5=eng.rhs_coefficients(cert,alg,[R0,R1,R2,R3],5)
 u4=basis(cert.WEIGHTS,4);a4=cert.linear_columns(alg,images,u4)
 rows4=eng.rows_for(a4,rhs5)
 assert len(u4)==152 and len(rows4)==1289
 det4=cert.minor(a4,{},sorted(set().union(*(set(c) for c in a4))),
                 witness["F4_rank_minor"],1000003)
 KR9=[{},{},scale(power(r,8),-9),{}]
 assert not cert.product(a4,cert.vector(u4,KR9))
 assert len(witness["F4_rank_minor"]["rows"])==151
 left5=eng.functional(witness["h5"],rows4)
 assert all(not eng.apply(left5,column) for column in a4)
 zero6=(0,)*6
 assert {m:eng.apply(left5,b) for m,b in rhs5.items() if eng.apply(left5,b)}=={zero6:Q(1)}
 assert all(m[2]==0 and m[5]==0 for m in rhs5) # u,w cancelled, not removed.
 # Deliberately corrupted input controls.
 badleft=dict(left_field);key=next(iter(rhs_field[zcube][0]))
 badleft[key]=badleft.get(key,Q(0))+1
 assert eng.apply(badleft,rhs_field[zcube][0])!=1
 badleft=dict(left5);key=next(iter(rhs5[zero6]))
 badleft[key]=badleft.get(key,Q(0))+1
 assert eng.apply(badleft,rhs5[zero6])!=1
 mon=next(m for m,tup in R3.items() if tup[1]);bad=[dict(z) for z in R3[mon]]
 key=next(iter(bad[1]));bad[1][key]+=1
 assert cert.product(a3,cert.vector(u3,bad))!=reduced_rhs[mon]
 # Separate rational normal-PBW positive control, after the NECESSARY field reductions.
 values=tuple(Q(z) for z in (-2,1,3,-1,2,2))
 normal_sizes=eng.pbw_control(eng.evaluate_frames(alg,[R0,R1,R2,R3],values),4)
 def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
 receipt={
  "claim":"No graded canonical Long full scalar Weyl jet of degree24 through h5, arbitrary parity and both complex h3 branches",
  "status":"PASS","environment":{"python":platform.python_version(),
   "arithmetic":"Fraction coefficient-pair field identities, rational functionals, bare modular minors and normal PBW control"},
  "cap":24,"parity":"unrestricted","field":"omega^2=-2",
  "complex_h3_inputs_verified_coefficientwise":True,
  "field_h4":{"unknowns":236,"rows":1673,"functional_support":len(left_field),
   "identity":"ellA=0; ellB=z^3+omega*0",
   "both_Galois_signs_verified":True,
   "no_separate_real_imag_parameter_constraints":True},
  "z0_h4":{"functional_support":len(left_L),
   "identity":"ellB=(c-27a/544)^3","field_conclusion":"a=544c/27",
   "rank_F3_Q":235,"minor_det_modp":det3,
   "particular_tuples":7,"particular_nonzero_coefficients":253,
   "full_F3_family":"retained particular+w*K_R6"},
  "h5":{"unknowns":152,"component_counts":[24,40,54,34],"rows":1289,
   "rank_F4_Q":151,"minor_det_modp":det4,"functional_support":len(left5),
   "all152columns_annihilated":True,"parameters":["b","c","u","v","gamma","w"],
   "rhs_monomials":[list(m) for m in sorted(rhs5,key=lambda m:(sum(m),m))],
   "constant_evaluation":"1","all_nonconstant_evaluations":"0",
   "u_w_cancel_without_normalization":True,"contradiction":"0=1 for every parameter point",
   "Pi3_F1_F1_included":True},
  "normal_PBW_positive_control":{"parameters_b_c_u_v_gamma_w":[str(v) for v in values],
   "all_six_relations_through_h4":True,"operator_terms":normal_sizes,
   "scope":"Finite convention/family control only"},
  "invalid_controls":{"corrupted_field_functional_rejected":True,
   "corrupted_F3_particular_rejected":True,"corrupted_h5_functional_rejected":True},
  "hashes":{"check_graded_long_cap24_full_jet.py":sha(Path(__file__)),
   "GRADED_LONG_CAP24_FULL_JET_OBSTRUCTION.md":sha(HERE/"GRADED_LONG_CAP24_FULL_JET_OBSTRUCTION.md"),
   "graded_long_cap24_full_jet_witness.json":sha(witness_path),
   "GRADED_LONG_CAP24_COMPLEX_FIRST_JET_FAMILY.md":sha(HERE/"GRADED_LONG_CAP24_COMPLEX_FIRST_JET_FAMILY.md"),
   "check_graded_long_cap24_complex_family.py":sha(HERE/"check_graded_long_cap24_complex_family.py"),
   "graded_long_cap24_complex_witness.json":sha(lower_path),
   "check_graded_long_cap23_full_jet.py":sha(HERE/"check_graded_long_cap23_full_jet.py"),
   "check_graded_long_cap23.py":sha(HERE/"check_graded_long_cap23.py"),
   "check_long_graph_curvature.py":sha(HERE/"check_long_graph_curvature.py"),
   "../dixmier-audit/verify_long_pbw.py":sha(source)},
  "coverage_boundary":"Complete graded cap24 finite coefficient family, arbitrary parity, "
   "both actual quadratic-field h3 branches retained until exact h4 z^3 identity; "
   "then full rational F3 family and h5 contradiction with allsix parameters. "
   "No real/imag separation, nongraded normalization, cap25 full-parity conclusion "
   "or DC1/DC2 theorem.",
  "independent_replay":False,"formal_verification":False,"external_review":False}
 out=args.output_dir or Path(tempfile.mkdtemp(prefix="graded-long-cap24-full-"))
 out.mkdir(parents=True,exist_ok=True);dest=out/"graded_long_cap24_full_jet_controls.json"
 dest.write_text(json.dumps(receipt,indent=2)+"\n")
 print("PASS: both field branches forcez=0; rational h4 reduction; h5 gives0=1 for allsix parameters")
 print("Functional supports",len(left_field),len(left_L),len(left5),"rank minors",det3,det4)
 print("Independent normal-PBW surviving h4 point PASS; receipt:",dest)

if __name__=="__main__":main()
