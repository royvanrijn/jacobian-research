#!/usr/bin/env python3
"""Complete arbitrary-parity graded Long cap25 h2/h3/h4/h5 certificate replay."""
from pathlib import Path
from fractions import Fraction as Q
from math import isqrt
import argparse,hashlib,importlib.util,json,platform,tempfile
HERE=Path(__file__).resolve().parent
def load(path,name):
 s=importlib.util.spec_from_file_location(name,path)
 m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def basis(weights,m):
 out=[]
 for i,w in enumerate(weights):
  for b in range(26):
   for c in range(26-b):
    for d in range(26-b-c):
     a=b+4*c+2*d+3*m-w
     if a>=0 and a+b+c+d<=25:out.append((i,(a,b,c,d)))
 return sorted(out)
def parse(cert,rows):
 out=[{} for _ in range(4)]
 for label,e,c in rows:
  i=cert.LABELS.index(label);e=tuple(e);assert e not in out[i];out[i][e]=Q(c)
 return out
def frame(cert,rows):
 out={}
 for item in rows:
  m=tuple(item["monomial"]);assert m not in out;out[m]=parse(cert,item["coefficients"])
 return out
def determinant_Q(matrix):
 a=[list(row) for row in matrix];n=len(a);result=Q(1)
 assert all(len(row)==n for row in a)
 for j in range(n):
  i=next((i for i in range(j,n) if a[i][j]),None)
  if i is None:return Q(0)
  if i!=j:a[j],a[i]=a[i],a[j];result=-result
  q=a[j][j];result*=q
  for i in range(j+1,n):
   r=a[i][j]/q
   for k in range(j+1,n):a[i][k]-=r*a[j][k]
 return result
# Exact, declared reduction modulo I=(d^2-250df+35625f^2,e*eta).
def reduce_monomial(m):
 if m[4] and m[6]:return {}
 if m[3]<2:return {m:Q(1)}
 x=list(m);x[3]-=1;x[5]+=1
 y=list(m);y[3]-=2;y[5]+=2
 out={}
 for z,q in [(tuple(x),Q(250)),(tuple(y),Q(-35625))]:
  for mon,c in reduce_monomial(z).items():out[mon]=out.get(mon,Q(0))+q*c
 return {mon:c for mon,c in out.items() if c}
def reduce_rhs(alg,rhs):
 out={}
 for m,p in rhs.items():
  for mon,c in reduce_monomial(m).items():out[mon]=alg.add(out.get(mon,{}),alg.scale(p,c))
 return {m:p for m,p in out.items() if p}
def unit(n,i=None,degree=1):
 e=[0]*n
 if i is not None:e[i]=degree
 return tuple(e)
def pscale(p,c):return {m:c*q for m,q in p.items() if c*q}
def padd(*polys):
 out={}
 for p in polys:
  for m,c in p.items():out[m]=out.get(m,Q(0))+c
 return {m:c for m,c in out.items() if c}
def pmul(a,b):
 out={}
 for m,c in a.items():
  for n,d in b.items():
   z=tuple(x+y for x,y in zip(m,n));out[z]=out.get(z,Q(0))+c*d
 return {m:c for m,c in out.items() if c}
def ppow(p,n,NP):
 out={unit(NP):Q(1)}
 for _ in range(n):out=pmul(out,p)
 return out
def substitute(alg,eng,source,images,NP):
 out={}
 for m,tup in source.items():
  polynomial={unit(NP):Q(1)}
  for n,img in zip(m,images):
   if n:polynomial=pmul(polynomial,ppow(img,n,NP))
  for mon,c in polynomial.items():eng.add_tuple(alg,out,mon,eng.tuple_scale(alg,tup,c))
 return out
def evaluations(eng,left,rhs):
 return {m:eng.apply(left,p) for m,p in rhs.items() if eng.apply(left,p)}
def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument("--output-dir",type=Path);args=parser.parse_args()
 cert=load(HERE/"check_graded_long_cap23.py","complete_columns")
 alg=load(HERE/"check_long_graph_curvature.py","independent_polynomials")
 eng=load(HERE/"check_graded_long_cap23_full_jet.py","general_parameter_engine")
 source=HERE.parent/"dixmier-audit/verify_long_pbw.py"
 prod=load(source,"canonical_input")
 f=[{tuple(e):Q(c) for e,c in p.items()} for p in prod.long_images()]
 path=HERE/"graded_long_cap25_full_jet_witness.json";w=json.loads(path.read_text())
 assert w["cap"]==25 and w["parity"]=="unrestricted"
 prime=w["prime"];assert prime==1000003 and all(prime%d for d in range(2,isqrt(prime)+1))
 assert w["source_weights"]==[-1,1,4,2] and w["image_weights"]==list(cert.WEIGHTS) and w["h_weight"]==3
 add,mul,scale=alg.add,alg.multiply,alg.scale;power=lambda p,n:cert.power(alg,p,n)
 variables=[{tuple(int(i==j) for i in range(4)):Q(1)} for j in range(4)]
 X,Qpos,P,Z=variables;one={(0,0,0,0):Q(1)};r,T,D,S=f;t=mul(X,Qpos)
 for i,j in cert.PAIRS:assert alg.bidifferential(f[i],f[j],1)==alg.bidifferential(variables[i],variables[j],1)
 B=add(scale(mul(power(X,2),P),3),scale(mul(add(one,scale(t,-3)),Z),2))
 A0=add(Qpos,scale(mul(X,power(Qpos,2)),3),scale(mul(X,B),Q(-1,3)))
 Y=add(mul(r,A0),scale(t,-3))
 H4=add(power(Y,4),scale(mul(t,power(Y,2)),4),scale(mul(t,Y),Q(-16,3)),scale(power(t,2),10))
 N5=add(scale(power(Y,5),9),scale(mul(t,power(Y,3)),45),scale(mul(power(t,2),Y),45),
         scale(mul(t,power(Y,2)),-60),scale(mul(t,Y),110),scale(power(t,2),-720))
 H6=add(mul(power(r,4),D),scale(N5,Q(1,25)))
 H5=add(scale(power(Y,5),9),scale(mul(t,power(Y,3)),45),scale(mul(t,power(Y,2)),-60),
         scale(mul(power(t,2),Y),135),scale(mul(t,Y),-10),scale(power(t,2),-555))
 M=mul(power(X,2),Z);J=power(mul(r,S),2)
 Gs=[t,Y,power(Y,2),mul(power(r,2),T),H4,H6,M,H5,J]
 jets=[[alg.bidifferential(G,z,1) for z in f] for G in Gs]
 u1=basis(cert.WEIGHTS,1);a1=cert.linear_columns(alg,f,u1)
 assert [sum(i==j for i,e in u1) for j in range(4)]==[102,141,171,127] and len(u1)==541
 vectors=[cert.vector(u1,p) for p in jets]
 assert all(not cert.product(a1,v) for v in vectors)
 assert jets[-1]==[{},scale(mul(power(r,2),S),2),scale(mul(r,power(S,2)),-2),{}]
 rows1=sorted(set().union(*(set(c) for c in a1)))
 det1=cert.minor(a1,{},rows1,w["h1_rank_minor"],prime)
 assert len(w["h1_rank_minor"]["rows"])==532
 indices=w["h1_independence_minor"]["indices"]
 assert len(indices)==len(set(indices))==9 and all(0<=i<len(u1) for i in indices)
 ind=determinant_Q([[v[i] for v in vectors] for i in indices])
 assert ind==Q(w["h1_independence_minor"]["det_Q"]) and ind
 # Only the uniform source XQ and target RS gauges have been removed.
 assert mul(r,S)==add(scale(t,-6),scale(power(Y,2),-3),scale(Y,-8))
 assert alg.bidifferential(t,P,1)==scale(Qpos,-1)
 assert alg.bidifferential(t,Z,1)==scale(X,-1)
 assert [alg.bidifferential(mul(r,S),z,1) for z in f]==[{},r,scale(S,-1),{}]
 NP=11;zero=unit(NP)
 F0={zero:f};F1={}
 for i,jet in enumerate([jets[j] for j in [1,3,4,5,6,7,8]]):F1[unit(NP,i)]=jet
 raw3=eng.rhs_coefficients(cert,alg,[F0,F1],3)
 u2=basis(cert.WEIGHTS,2);a2=cert.linear_columns(alg,f,u2)
 assert [sum(i==j for i,e in u2) for j in range(4)]==[71,102,127,91] and len(u2)==391
 rows3=eng.rows_for(a2,raw3)
 det2=cert.minor(a2,{},rows3,w["F2_rank_minor"],prime)
 assert len(w["F2_rank_minor"]["rows"])==387
 cone={unit(NP,3,2):Q(1),tuple(int(j==3)+int(j==5) for j in range(NP)):Q(-250),unit(NP,5,2):Q(35625)}
 eh={tuple(int(j==4)+int(j==6) for j in range(NP)):Q(1)}
 h3supports=[]
 for record,want in zip(w["h3_necessity_functionals"],[cone,eh]):
  left=eng.functional(record,rows3)
  assert all(not eng.apply(left,col) for col in a2)
  assert evaluations(eng,left,raw3)==want
  h3supports.append(len(left))
 F2=frame(cert,w["F2_mod_h3_ideal"])
 assert all(len(m)==NP and all(n>=0 for n in m) for m in F2)
 assert {m for m in F2 if any(m[7:])}=={unit(NP,i) for i in range(7,11)}
 reduced3=reduce_rhs(alg,raw3)
 for mon in set(F2)|set(reduced3):
  tup=F2.get(mon,[{} for _ in range(4)])
  assert cert.product(a2,cert.vector(u2,tup))==reduced3.get(mon,{})
 priorpath=HERE/"graded_long_cap23_witness.json"
 prior=json.loads(priorpath.read_text());assert F2[zero]==parse(cert,prior["first_particular"])
 kernelG=[power(r,3),mul(power(r,4),S),power(X,3),mul(power(r,3),add(Y,scale(one,2)))]
 kernel2=[[alg.bidifferential(G,z,1) for z in f] for G in kernelG]
 for i,k in enumerate(kernel2):
  assert F2[unit(NP,i+7)]==k
  assert not cert.product(a2,cert.vector(u2,k))
 assert kernel2[0][2] and not kernel2[0][1] and not kernel2[0][3]
 assert kernel2[1][1] and not kernel2[1][3]
 assert any(e[2]+e[3]==1 for e in kernel2[2][3])
 assert kernel2[3][3] and all(e[2]+e[3]==0 for e in kernel2[3][3])
 # Complete F2 family: rank387/391 and these four independent directions.
 full4=reduce_rhs(alg,eng.rhs_coefficients(cert,alg,[F0,F1,F2],4))
 u3=basis(cert.WEIGHTS,3);a3=cert.linear_columns(alg,f,u3)
 assert [sum(i==j for i,e in u3) for j in range(4)]==[47,71,91,62] and len(u3)==271
 rows4=eng.rows_for(a3,full4)
 det3=cert.minor(a3,{},rows4,w["h4_rank_minor"],prime)
 assert len(w["h4_rank_minor"]["rows"])==270
 KR6=[{},{},scale(power(r,5),-6),{}]
 assert not cert.product(a3,cert.vector(u3,KR6))
 cubes=[]
 for record,index in zip(w["h4_cube_functionals"],[5,4]):
  left=eng.functional(record,rows4)
  assert all(not eng.apply(left,col) for col in a3)
  assert evaluations(eng,left,full4)=={unit(NP,index,3):Q(1)}
  cubes.append(len(left))
 # Over fields f=e=0; the retained conic then gives d=0, including both
 # quadratic-field branches without a real/imaginary parameter split.
 forced=lambda F:{m:p for m,p in F.items() if not(m[3] or m[4] or m[5])}
 forced_frames=[forced(F) for F in [F0,F1,F2]]
 red4=eng.rhs_coefficients(cert,alg,forced_frames,4)
 redrows=eng.rows_for(a3,red4)
 a,b,c,eta,gamma,delta=[{unit(NP,j):Q(1)} for j in [0,1,2,6,9,10]]
 ell=padd(a,pscale(c,Q(-544,27)))
 KK=padd(pscale(ppow(c,2,NP),976),pscale(pmul(c,eta),5184),pscale(delta,27))
 targets=[pmul(ell,KK),pmul(padd(c,pscale(eta,36)),KK),
          padd(pscale(pmul(gamma,eta),864),pscale(ppow(ell,3,NP),-1))]
 factor_supports=[]
 for record,want in zip(w["h4_reduced_factor_functionals"],targets):
  left=eng.functional(record,redrows)
  assert all(not eng.apply(left,col) for col in a3)
  assert evaluations(eng,left,red4)==want
  factor_supports.append(len(left))
 # Four FIELD strata cover eta!=0/eta=0, including the denominator-zero
 # loci. Substitution is parameterization, not a source/target gauge.
 BP=7;z7=unit(BP);pars=[{unit(BP,i):Q(1)} for i in range(BP)]
 bb,cc,LL,hh,uu,vv,ww=pars;ratio=Q(544,27)
 configurations={
 "A_nonzero_h":[padd(pscale(cc,ratio),LL),bb,cc,{},{},{},hh,uu,vv,
   pscale(pmul(ppow(LL,3,BP),{unit(BP,3,-1):Q(1)}),Q(1,864)),
   pscale(padd(pscale(ppow(cc,2,BP),976),pscale(pmul(cc,hh),5184)),Q(-1,27))],
 "A_h_zero":[pscale(cc,ratio),bb,cc,{},{},{},{},uu,vv,LL,pscale(ppow(cc,2,BP),Q(-976,27))],
 "B_gamma_zero":[pscale(hh,Q(-2176,3)),bb,pscale(hh,-36),{},{},{},hh,uu,vv,{},LL],
 "B_h_zero":[{},bb,{},{},{},{},{},uu,vv,LL,hh]}
 u4=basis(cert.WEIGHTS,4);a4=cert.linear_columns(alg,f,u4)
 assert [sum(i==j for i,e in u4) for j in range(4)]==[29,47,62,40] and len(u4)==178
 branch_controls={}
 for name,images in configurations.items():
  record=w["field_strata"][name]
  B0={z7:f};B1=frame(cert,record["F1"]);B2=frame(cert,record["F2"]);B3=frame(cert,record["F3"])
  assert record["eta_nonzero"]==(name=="A_nonzero_h")
  assert B1==substitute(alg,eng,forced_frames[1],images,BP)
  assert B2==substitute(alg,eng,forced_frames[2],images,BP)
  rhs4=eng.rhs_coefficients(cert,alg,[B0,B1,B2],4)
  assert B3[unit(BP,6)]==KR6
  assert set(B3)==set(rhs4)|{unit(BP,6)}
  for m,tup in B3.items():assert cert.product(a3,cert.vector(u3,tup))==rhs4.get(m,{})
  rhs5=eng.rhs_coefficients(cert,alg,[B0,B1,B2,B3],5)
  assert [list(m) for m in sorted(rhs5)]==record["h5_monomials"]
  rows5=eng.rows_for(a4,rhs5)
  left=eng.functional(record["h5_constant_functional"],rows5)
  assert all(not eng.apply(left,col) for col in a4)
  assert evaluations(eng,left,rhs5)=={z7:Q(1)}
  bad=dict(left);key=next(iter(rhs5[z7]));bad[key]=bad.get(key,Q(0))+1
  assert eng.apply(bad,rhs5[z7])!=1
  if name=="A_nonzero_h":
   pbw_values=tuple(Q(x) for x in (0,0,1,1,0,0,0))
   pbw_sizes=eng.pbw_control(eng.evaluate_frames(alg,[B0,B1,B2,B3],pbw_values),4)
  branch_controls[name]={"eta_nonzero":record["eta_nonzero"],
    "coordinate_order":record["coordinate_order"],
    "F3_blocks":len(B3),"F3_coefficients":sum(sum(map(len,t)) for t in B3.values()),
    "h5_RHS_monomials":[list(m) for m in sorted(rhs5)],"h5_rows":len(rows5),
    "functional_support":len(left),"all178columns_annihilated":True,
    "constant_evaluation":"1","all_nonconstant_evaluations":"0"}
 assert eng.rhs_indices(4,3)==[(1,2,1,Q(1)),(2,1,1,Q(1)),(0,1,3,Q(1,24)),(1,0,3,Q(1,24))]
 assert eng.rhs_indices(5,4)==[(1,3,1,Q(1)),(2,2,1,Q(1)),(3,1,1,Q(1)),
  (0,2,3,Q(1,24)),(1,1,3,Q(1,24)),(2,0,3,Q(1,24)),(0,0,5,Q(1,1920))]
 def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
 hashes={name:sha(HERE/name) for name in [
 "GRADED_LONG_CAP25_FULL_JET_OBSTRUCTION.md","check_graded_long_cap25_full_jet.py",
 "graded_long_cap25_full_jet_witness.json","graded_long_cap23_witness.json",
 "check_graded_long_cap23.py","check_graded_long_cap23_full_jet.py","check_long_graph_curvature.py"]}
 hashes["../dixmier-audit/verify_long_pbw.py"]=sha(source)
 receipt={"claim":"No arbitrary-parity graded canonical Long scalar full jet of spatial degree25 through h5",
  "status":"PASS","environment":{"python":platform.python_version(),
    "arithmetic":"Fraction coefficient identities, explicit ideal reductions, bare modular minors, rational functionals and one separate normal-PBW positive control"},
  "cap":25,"parity":"unrestricted","prime":prime,
  "first_kernel":{"counts":[102,141,171,127],"unknowns":541,"rank_Q":532,
    "independent_kernels":9,"minor_det_modp":det1,"independence_det_Q":str(ind),
    "new_direction":"(R*S)^2; no cap-growing normalization used"},
  "h3":{"F2_unknowns":391,"rank_Q":387,"minor_det_modp":det2,
    "necessity_functional_supports":h3supports,"necessary_sufficient_over_fields":"d^2-250df+35625f^2=0; e*eta=0",
    "full_F2_ideal_reduction_verified":True,"all_four_kernel_parameters_kept":True},
  "h4":{"F3_unknowns":271,"rank_Q":270,"minor_det_modp":det3,
    "cube_functional_supports":cubes,"forced_over_fields":"f=e=d=0",
    "factor_functional_supports":factor_supports,
    "ell":"a-544c/27","K":"976c^2+5184c*eta+27delta",
    "remaining_conditions":["ell*K=0","(c+36eta)*K=0","864gamma*eta=ell^3"],
    "full_F3_kernel":"w*K_R6"},
  "h5":{"F4_counts":[29,47,62,40],"all_F4_unknowns":178,
    "field_strata":branch_controls,"all_strata_contradiction":"0=1",
    "Laurent_eta_nonzero_branch_and_eta_zero_strata_kept":True,
    "Pi3_F1_F1_included":True,"corrupted_four_functionals_rejected":True},
  "normal_PBW_positive_control":{
    "stratum":"A_nonzero_h","coordinate_order":["b","c","ell","eta","u","v","w"],
    "parameter_values":[str(x) for x in pbw_values],
    "source_parameters":{"a":"1","gamma":"1/864","delta":"0"},
    "all_six_relations_through_h4":True,"normal_operator_terms":pbw_sizes,
    "scope":"One finite convention/family control; universal proof remains the exact parameter identities."},
  "hashes":hashes,
  "coverage_boundary":"Complete FIELD classification for the stated cap25 grading, arbitrary parity, only uniform XQ/RS first gauges. Both complex conic branches retained until exact cube obstruction; all image coefficients allowed. No nonreduced ideal equality, evenization, cap26 conclusion, full nongraded-lift or DC1/DC2 theorem.",
  "independent_replay":False,"formal_verification":False,"external_review":False}
 output=args.output_dir or Path(tempfile.mkdtemp(prefix="graded-long-cap25-full-"));output.mkdir(parents=True,exist_ok=True)
 target=output/"graded_long_cap25_full_jet_controls.json";target.write_text(json.dumps(receipt,indent=2)+"\n")
 print("PASS: cap25 full-parity kernel, all field strata, four constant h5 contradictions")
 print("Minors:",det1,det2,det3,"supports:",[p["functional_support"] for p in branch_controls.values()])
 print("Receipt:",target)
if __name__=="__main__":main()
