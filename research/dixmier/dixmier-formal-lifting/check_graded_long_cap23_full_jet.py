#!/usr/bin/env python3
"""Full-parity graded Long cap23 obstruction: exact h4/h5 certificates."""
from fractions import Fraction as Q
from pathlib import Path
from math import factorial,comb
import argparse,hashlib,importlib.util,itertools,json,platform,tempfile
HERE=Path(__file__).resolve().parent
ZERO=(0,0,0,0,0)
def load(path,name):
 spec=importlib.util.spec_from_file_location(name,path)
 mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
def plus(a,b):return tuple(x+y for x,y in zip(a,b))
def tuple_add(alg,a,b):return [alg.add(x,y) for x,y in zip(a,b)]
def tuple_scale(alg,a,c):return [alg.scale(x,c) for x in a]
def add_tuple(alg,dic,key,tup):
 dic[key]=tuple_add(alg,dic.get(key,[{} for _ in range(4)]),tup)
 if not any(dic[key]):del dic[key]
def rhs_indices(order,available):
 out=[]
 for ell in range(1,order+1,2):
  for a in range(available):
   b=order-ell-a
   if 0<=b<available:out.append((a,b,ell,Q(1,2**(ell-1)*factorial(ell))))
 return out
def rhs_coefficients(cert,alg,frames,order):
 """Negative coefficient at h^order; frames omit unknown F_(order-1)."""
 out={}
 for pi,(i,j) in enumerate(cert.PAIRS):
  for a,b,ell,factor in rhs_indices(order,len(frames)):
   for ma,fa in frames[a].items():
    for mb,fb in frames[b].items():
     mon=plus(ma,mb)
     poly=alg.scale(alg.bidifferential(fa[i],fb[j],ell),-factor)
     dest=out.setdefault(mon,{})
     for e,c in poly.items():dest[(pi,e)]=dest.get((pi,e),Q(0))+c
 return {m:{e:c for e,c in p.items() if c} for m,p in out.items()
         if any(p.values())}
def rows_for(columns,rhs):
 return sorted(set().union(*(set(c) for c in columns),*(set(b) for b in rhs.values())))
def functional(data,rows):
 assert data["row_count"]==len(rows)
 out={}
 for index,key,c in data["functional"]:
  key=(key[0],tuple(key[1]))
  assert rows[index]==key and key not in out
  out[key]=Q(c)
 return out
def apply(left,poly):
 return sum((c*poly.get(key,Q(0)) for key,c in left.items()),Q(0))
def evaluate_frames(alg,frames,values):
 result=[]
 for frame in frames:
  tup=[{} for _ in range(4)]
  for mon,term in frame.items():
   factor=Q(1)
   for x,n in zip(values,mon):factor*=x**n
   tup=tuple_add(alg,tup,tuple_scale(alg,term,factor))
  result.append(tup)
 return [[result[k][i] for k in range(len(result))] for i in range(4)]

# A separate positive control: commutators below use ONLY normal PBW
# multiplication, not the Moyal bidifferential arithmetic.
def clean(poly):return {e:c for e,c in poly.items() if c}
def normal_order(jet,cutoff):
 out={}
 for hdegree,poly in enumerate(jet):
  for e,c in poly.items():
   for n in range(min(e[0],e[2],cutoff-hdegree)+1):
    for m in range(min(e[1],e[3],cutoff-hdegree-n)+1):
     factor=Q(factorial(e[0])*factorial(e[2])*factorial(e[1])*factorial(e[3]),
       factorial(e[0]-n)*factorial(e[2]-n)*factorial(e[1]-m)*factorial(e[3]-m)
       *2**(n+m)*factorial(n)*factorial(m))
     key=(e[0]-n,e[1]-m,e[2]-n,e[3]-m,hdegree+n+m)
     out[key]=out.get(key,Q(0))+c*factor
 return clean(out)
def pbw(left,right,cutoff):
 out={}
 for e,c in left.items():
  for f,d in right.items():
   budget=cutoff-e[4]-f[4]
   if budget<0:continue
   for n in range(min(e[2],f[0],budget)+1):
    cn=comb(e[2],n)*factorial(f[0])//factorial(f[0]-n)
    for m in range(min(e[3],f[1],budget-n)+1):
     cm=comb(e[3],m)*factorial(f[1])//factorial(f[1]-m)
     key=(e[0]+f[0]-n,e[1]+f[1]-m,e[2]+f[2]-n,e[3]+f[3]-m,e[4]+f[4]+n+m)
     out[key]=out.get(key,Q(0))+c*d*cn*cm
 return clean(out)
def difference(a,b):
 out=dict(a)
 for e,c in b.items():out[e]=out.get(e,Q(0))-c
 return clean(out)
def pbw_control(jets,cutoff):
 operators=[normal_order(jet,cutoff) for jet in jets]
 for i,j in itertools.combinations(range(4),2):
  actual=difference(pbw(operators[i],operators[j],cutoff),
                    pbw(operators[j],operators[i],cutoff))
  expected={(0,0,0,0,1):Q(-1)} if (i,j) in ((0,2),(1,3)) else {}
  assert actual==expected,(i,j,len(actual))
 return [len(op) for op in operators]

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument("--output-dir",type=Path);args=parser.parse_args()
 cert=load(HERE/"check_graded_long_cap23.py","complete_column_engine")
 tail=load(HERE/"check_graded_long_tail_rigidity.py","complete_tail_bases")
 alg=load(HERE/"check_long_graph_curvature.py","independent_polynomial_engine")
 source=HERE.parent/"dixmier-audit/verify_long_pbw.py"
 producer=load(source,"canonical_global_input")
 images=[{tuple(e):Q(c) for e,c in poly.items()} for poly in producer.long_images()]
 witness_path=HERE/"graded_long_cap23_full_jet_witness.json"
 witness=json.loads(witness_path.read_text());assert witness["cap"]==23
 odd_path=HERE/witness["prior_odd_h3_witness"]
 odd=json.loads(odd_path.read_text())
 tail_path=HERE/witness["tail_witness"];tail_data=json.loads(tail_path.read_text())
 add,mul,scale=alg.add,alg.multiply,alg.scale
 power=lambda p,n:cert.power(alg,p,n)
 X,Qpos,P,Z=[{tuple(int(i==j) for i in range(4)):Q(1)} for j in range(4)]
 one={(0,0,0,0):Q(1)};r,T,D,S=images;t=mul(X,Qpos)
 B=add(scale(mul(power(X,2),P),3),scale(mul(add(one,scale(t,-3)),Z),2))
 A0=add(Qpos,scale(mul(X,power(Qpos,2)),3),scale(mul(X,B),Q(-1,3)))
 Y=add(mul(r,A0),scale(t,-3))
 H4=add(power(Y,4),scale(mul(t,power(Y,2)),4),
        scale(mul(t,Y),Q(-16,3)),scale(power(t,2),10))
 g1={tuple(int(i==j) for i in range(5)):
      [alg.bidifferential(G,f,1) for f in images]
      for j,G in enumerate((Y,mul(power(r,2),T),H4))}
 mapping={"1":ZERO,"a^2":(2,0,0,0,0),"a*b":(1,1,0,0,0),
   "a*c":(1,0,1,0,0),"b^2":(0,2,0,0,0),"b*c":(0,1,1,0,0),"c^2":(0,0,2,0,0)}
 g2={}
 for item in odd["h3_compatibility"]["h2_particulars_on_d_zero"]:
  tup=[{} for _ in range(4)]
  for label,e,c in item["coefficients"]:tup[cert.LABELS.index(label)][tuple(e)]=Q(c)
  add_tuple(alg,g2,mapping[item["monomial"]],tup)
 g2[(0,0,0,1,0)]=[{},{},scale(power(r,2),-3),{}]
 g2[(0,0,0,0,1)]=[{},power(r,4),scale(mul(power(r,3),S),-4),{}]
 F0={ZERO:images}
 # Recheck every retained second coefficient, including both free kernels.
 rhs3=rhs_coefficients(cert,alg,[F0,g1],3)
 u2=tail.basis(cert.WEIGHTS,2);a2=cert.linear_columns(alg,images,u2)
 for mon in set(g2)|set(rhs3):
  assert cert.product(a2,cert.vector(u2,g2.get(mon,[{} for _ in range(4)])))==rhs3.get(mon,{})

 expected4=[(1,2,1,Q(1)),(2,1,1,Q(1)),
            (0,1,3,Q(1,24)),(1,0,3,Q(1,24))]
 expected5=[(1,3,1,Q(1)),(2,2,1,Q(1)),(3,1,1,Q(1)),
            (0,2,3,Q(1,24)),(1,1,3,Q(1,24)),(2,0,3,Q(1,24)),
            (0,0,5,Q(1,1920))]
 assert rhs_indices(4,3)==expected4 and rhs_indices(5,4)==expected5
 rhs4=rhs_coefficients(cert,alg,[F0,g1,g2],4)
 u3=tail.basis(cert.WEIGHTS,3);a3=cert.linear_columns(alg,images,u3)
 rows3=rows_for(a3,rhs4)
 assert len(u3)==206 and len(rows3)==1563
 det3=cert.minor(a3,{},rows3,tail_data[0]["rank_minor"],1000003)
 KR6=[{},{},scale(power(r,5),-6),{}]
 assert not cert.product(a3,cert.vector(u3,KR6))
 left4=functional(witness["h4"],rows3)
 assert witness["h4"]["coefficient_columns"]==len(a3)
 assert all(not apply(left4,column) for column in a3)
 expectedL={(0,0,1,0,0):Q(1),(1,0,0,0,0):Q(-27,544)}
 assert {mon:apply(left4,b) for mon,b in rhs4.items() if apply(left4,b)}==expectedL
 # So a=544c/27. Keep b,c,u,v and the homogeneous F3 parameter w.
 ratio=Q(544,27)
 def normalize(frame):
  out={}
  for m,tup in frame.items():
   key=(m[1],m[0]+m[2],m[3],m[4],0)
   add_tuple(alg,out,key,tuple_scale(alg,tup,ratio**m[0]))
  return out
 F1,F2=normalize(g1),normalize(g2)
 reduced_rhs4=rhs_coefficients(cert,alg,[F0,F1,F2],4)
 F3={}
 for item in witness["F3_family"]["particulars"]:
  mon=tuple(item["parameter_monomial_b_c_u_v"])+(0,)
  tup=[{} for _ in range(4)]
  for label,e,c in item["coefficients"]:
   i=cert.LABELS.index(label);e=tuple(e);assert e not in tup[i];tup[i][e]=Q(c)
  assert not tup[0] and not tup[3]
  assert cert.product(a3,cert.vector(u3,tup))==reduced_rhs4.get(mon,{})
  F3[mon]=tup
 assert set(F3)==set(reduced_rhs4)
 F3[(0,0,0,0,1)]=KR6
 assert sum(sum(map(len,tup)) for m,tup in F3.items() if m!=(0,0,0,0,1))==146

 rhs5=rhs_coefficients(cert,alg,[F0,F1,F2,F3],5)
 u4=tail.basis(cert.WEIGHTS,4);a4=cert.linear_columns(alg,images,u4)
 rows4=rows_for(a4,rhs5)
 assert len(u4)==130 and len(rows4)==1207
 det4=cert.minor(a4,{},rows4,tail_data[1]["rank_minor"],1000003)
 left5=functional(witness["h5"],rows4)
 assert witness["h5"]["coefficient_columns"]==130
 assert all(not apply(left5,column) for column in a4)
 assert {m:apply(left5,b) for m,b in rhs5.items() if apply(left5,b)}=={ZERO:Q(1)}
 assert all(m[2]==0 and m[4]==0 for m in rhs5)  # u,w cancel exactly.
 # The displayed first/third corrections have R/S3 zero; none of this
 # is a positional first-image assumption on the enumerated unknowns.
 assert all(not tup[0] for tup in F1.values())
 assert all(not tup[0] for tup in F2.values())
 assert all(not tup[0] and not tup[3] for tup in F3.values())
 # Deliberate invalid controls.
 badleft=dict(left5);key=next(iter(rhs5[ZERO]));badleft[key]=badleft.get(key,Q(0))+1
 assert apply(badleft,rhs5[ZERO])!=1
 firstmon=next(m for m,tup in F3.items() if tup[1])
 bad=[dict(poly) for poly in F3[firstmon]]
 key=next(iter(bad[1]));bad[1][key]+=1
 assert cert.product(a3,cert.vector(u3,bad))!=reduced_rhs4[firstmon]
 # Independent normal-PBW positive control at a nontrivial five-parameter point.
 assert normal_order([{(1,0,1,0):Q(1)}],1)=={
   (1,0,1,0,0):Q(1),(0,0,0,0,1):Q(1,2)}
 values=tuple(Q(c) for c in witness["PBW_positive_control"]["parameters_b_c_u_v_w"])
 normal_sizes=pbw_control(evaluate_frames(alg,[F0,F1,F2,F3],values),4)

 def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
 receipt={
  "claim":"No graded canonical Long full scalar Weyl jet of spatial degree23 through commutator h5, with arbitrary parity",
  "status":"PASS","environment":{"python":platform.python_version(),
   "arithmetic":"Exact Fraction identities, bare modular minors, rational functionals and separate normal-PBW positive control"},
  "cap":23,"parity":"unrestricted","source_weights":[-1,1,4,2],"h_weight":3,
  "h4":{"unknown_counts":[34,54,71,47],"unknowns":206,"rows":1563,
   "rank_Q":205,"minor_det_modp":det3,"functional_support":len(left4),
   "necessary_relation":"a=544*c/27","complete_F3_particulars":5,
   "particular_nonzero_coefficients":146,"homogeneous":"w*K_R6",
   "R3_S3_forced_zero":True},
  "h5":{"unknown_counts":[20,34,47,29],"unknowns":130,"rows":1207,
   "rank_Q":130,"minor_det_modp":det4,"functional_support":len(left5),
   "all130columns_annihilated":True,
   "parameter_order":["b","c","u","v","w"],
   "rhs_monomials":[list(m) for m in sorted(rhs5,key=lambda m:(sum(m),m))],
   "constant_evaluation":"1","all_nonconstant_evaluations":"0",
   "u_w_cancel_without_normalization":True,"contradiction":"0=1 for every parameter point",
   "Pi3_F1_F1_included":True},
  "coefficient_formula":"sum_{a+b+ell=N,ell odd} Pi^ell(Fa_i,Fb_j)/(2^(ell-1)*ell!)",
  "rhs_indices_h4":[[a,b,k,str(q)] for a,b,k,q in expected4],
  "rhs_indices_h5":[[a,b,k,str(q)] for a,b,k,q in expected5],
  "PBW_positive_control":{"parameters_b_c_u_v_w":[str(c) for c in values],
   "all_six_relations_through_h4":True,"normal_operator_terms":normal_sizes,
   "scope":"A finite convention/family control, not a universal certificate"},
  "invalid_controls":{"corrupted_F3_particular_rejected":True,
                       "corrupted_h5_functional_rejected":True},
  "hashes":{"check_graded_long_cap23_full_jet.py":sha(Path(__file__)),
   "GRADED_LONG_CAP23_FULL_JET_OBSTRUCTION.md":sha(HERE/"GRADED_LONG_CAP23_FULL_JET_OBSTRUCTION.md"),
   "graded_long_cap23_full_jet_witness.json":sha(witness_path),
   "ODD_GRADED_LONG_FIRST_JETS_AND_NORMALIZATION_GAP.md":sha(HERE/"ODD_GRADED_LONG_FIRST_JETS_AND_NORMALIZATION_GAP.md"),
   "check_odd_graded_long_first_jets.py":sha(HERE/"check_odd_graded_long_first_jets.py"),
   "odd_graded_long_first_jets_witness.json":sha(odd_path),
   "GRADED_LONG_TAIL_RIGIDITY.md":sha(HERE/"GRADED_LONG_TAIL_RIGIDITY.md"),
   "check_graded_long_tail_rigidity.py":sha(HERE/"check_graded_long_tail_rigidity.py"),
   "graded_long_tail_rigidity_witness.json":sha(tail_path),
   "check_graded_long_cap23.py":sha(HERE/"check_graded_long_cap23.py"),
   "check_long_graph_curvature.py":sha(HERE/"check_long_graph_curvature.py"),
   "../dixmier-audit/verify_long_pbw.py":sha(source)},
  "coverage_boundary":"Complete graded cap23 finite coefficient family with both parities, "
   "using the pinned complete odd h3 classification and homogeneous tail theorem. "
   "Exact h4 reduction and full F3 family, then universal h5 contradiction with all "
   "five parameters retained. No nongraded normalization, cap24 full-parity exclusion "
   "or DC1/DC2 theorem is asserted.",
  "independent_replay":False,"formal_verification":False,"external_review":False}
 out=args.output_dir or Path(tempfile.mkdtemp(prefix="graded-long-cap23-full-"))
 out.mkdir(parents=True,exist_ok=True);dest=out/"graded_long_cap23_full_jet_controls.json"
 dest.write_text(json.dumps(receipt,indent=2)+"\n")
 print("PASS: h4 exact relation and complete F3 family; h5 gives0=1 for all five parameters")
 print("Functional supports",len(left4),len(left5),"tail minor determinants",det3,det4)
 print("Independent normal-PBW h4 positive control PASS; receipt:",dest)

if __name__=="__main__":main()
