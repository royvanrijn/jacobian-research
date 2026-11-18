#!/usr/bin/env python3
"""Complete arbitrary-parity graded Long cap26 h2/h3/h4/h5 certificate replay."""
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
  for b in range(27):
   for c in range(27-b):
    for d in range(27-b-c):
     a=b+4*c+2*d+3*m-w
     if a>=0 and a+b+c+d<=26:out.append((i,(a,b,c,d)))
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
def reduce_monomial(m):
 if m[7] and (m[3] or m[4] or m[5]):return {}
 if m[3]<2:return {m:Q(1)}
 x=list(m);x[3]-=1;x[5]+=1
 y=list(m);y[3]-=2;y[5]+=2
 out={}
 for z,q in [(tuple(x),Q(250)),(tuple(y),Q(-35625))]:
  for mon,c in reduce_monomial(z).items():out[mon]=out.get(mon,Q(0))+q*c
 return {m:c for m,c in out.items() if c}
def reduce_rhs(alg,rhs):
 out={}
 for m,p in rhs.items():
  for mon,c in reduce_monomial(m).items():out[mon]=alg.add(out.get(mon,{}),alg.scale(p,c))
 return {m:p for m,p in out.items() if p}
def poly_rows(rows):return {tuple(m):Q(q) for m,q in rows if Q(q)}
def check_functionals(eng,columns,rhs,records,targets):
 assert len(records)==len(targets);rows=eng.rows_for(columns,rhs);supports=[]
 for record,want in zip(records,targets):
  left=eng.functional(record,rows)
  assert all(not eng.apply(left,col) for col in columns)
  assert evaluations(eng,left,rhs)==want
  supports.append(len(left))
 return supports
def variables(n):return [{unit(n,i):Q(1)} for i in range(n)]
def configurations(name):
 names={
  "Efull":"a b c e u v gamma w","Ezero":"e b u v gamma delta w",
  "A0":"a b c h u v gamma w","HZres":"a b c h u v gamma chi w",
  "Bzero":"b u v gamma delta chi psi w","Cspec":"b c psi u v gamma w",
  "Bmain":"b z t K u v gamma w"}[name].split()
 NP=len(names);V=variables(NP);power=lambda p,n:ppow(p,n,NP)
 if name=="Efull":
  a,b,c,e,u,v,gamma=V[:7]
  imgs=[a,b,c,{},e,{},{},{},u,v,gamma,
   pscale(padd(pscale(power(c,2),976),pscale(pmul(c,e),216)),Q(-1,27)),
   pscale(pmul(c,e),-18),{}]
 elif name=="Ezero":
  e,b,u,v,gamma,delta=V[:6]
  imgs=[{},b,{},{},e,{},{},{},u,v,gamma,delta,{},{}]
 elif name=="A0":
  a,b,c,h,u,v,gamma=V[:7]
  imgs=[a,b,c,{},{},{},h,{},u,v,gamma,
   padd(pscale(power(c,2),Q(-976,27)),pscale(pmul(c,h),-192)),{},{}]
 elif name=="HZres":
  a,b,c,h,u,v,gamma,chi=V[:8]
  imgs=[a,b,c,{},{},{},h,pscale(h,6),u,v,gamma,
   padd(pscale(power(c,2),Q(-976,27)),pscale(pmul(c,h),Q(-1552,3)),pscale(power(h,2),-2052)),chi,{}]
 elif name=="Bzero":
  b,u,v,gamma,delta,chi,psi=V[:7]
  imgs=[{},b,{},{},{},{},{},{},u,v,gamma,delta,chi,psi]
 elif name=="Cspec":
  b,c,psi,u,v,gamma=V[:6]
  imgs=[pscale(c,Q(544,27)),b,c,{},{},{},pscale(c,Q(-22,45)),{},u,v,gamma,
   padd(pscale(power(c,2),Q(7792,135)),pscale(psi,Q(-908,3))),{},psi]
 else:
  b,z,t,K,u,v,gamma=V[:7]
  h=pscale(padd(z,pscale(t,-1)),Q(1,6));c=padd(pscale(t,6),pscale(z,Q(-3,2)))
  Lin=pscale(pmul(pmul(z,power(t,2)),{unit(NP,3,-1):Q(1)}),-576)
  a=padd(pscale(c,Q(544,27)),Lin)
  delta=padd(K,pscale(power(c,2),Q(-976,27)),pscale(pmul(c,h),-192),
   pscale(pmul(c,z),Q(-488,9)),pscale(pmul(h,z),-144),pscale(power(z,2),-33))
  chi=padd(pscale(pmul(pmul(K,z),{unit(NP,2,-1):Q(1)}),Q(-3,8)),
   pscale(pmul(c,z),Q(-100,3)),pscale(pmul(z,padd(z,pscale(t,48))),Q(9,2)))
  imgs=[a,b,c,{},{},{},h,z,u,v,gamma,delta,chi,{}]
 return names,imgs
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--output-dir",type=Path);args=ap.parse_args()
 cert=load(HERE/"check_graded_long_cap23.py","complete_columns")
 alg=load(HERE/"check_long_graph_curvature.py","independent_polynomials")
 eng=load(HERE/"check_graded_long_cap23_full_jet.py","general_parameter_engine")
 source=HERE.parent/"dixmier-audit/verify_long_pbw.py";prod=load(source,"canonical_input")
 f=[{tuple(e):Q(c) for e,c in p.items()} for p in prod.long_images()]
 wp=HERE/"graded_long_cap26_full_jet_witness.json";w=json.loads(wp.read_text())
 assert w["cap"]==26 and w["parity"]=="unrestricted"
 prime=w["prime"];assert prime==1000003 and all(prime%d for d in range(2,isqrt(prime)+1))
 assert w["source_weights"]==[-1,1,4,2] and w["image_weights"]==list(cert.WEIGHTS) and w["h_weight"]==3
 add,mul,scale=alg.add,alg.multiply,alg.scale;power=lambda p,n:cert.power(alg,p,n)
 Vs=[{tuple(int(i==j) for i in range(4)):Q(1)} for j in range(4)]
 X,Qpos,P,Z=Vs;one={(0,0,0,0):Q(1)};r,T,D,S=f;t=mul(X,Qpos)
 for i,j in cert.PAIRS:assert alg.bidifferential(f[i],f[j],1)==alg.bidifferential(Vs[i],Vs[j],1)
 B=add(scale(mul(power(X,2),P),3),scale(mul(add(one,scale(t,-3)),Z),2))
 A0=add(Qpos,scale(mul(X,power(Qpos,2)),3),scale(mul(X,B),Q(-1,3)))
 Y=add(mul(r,A0),scale(t,-3))
 H4=add(power(Y,4),scale(mul(t,power(Y,2)),4),scale(mul(t,Y),Q(-16,3)),scale(power(t,2),10))
 N5=add(scale(power(Y,5),9),scale(mul(t,power(Y,3)),45),scale(mul(power(t,2),Y),45),
  scale(mul(t,power(Y,2)),-60),scale(mul(t,Y),110),scale(power(t,2),-720))
 H6=add(mul(power(r,4),D),scale(N5,Q(1,25)))
 H5=add(scale(power(Y,5),9),scale(mul(t,power(Y,3)),45),scale(mul(t,power(Y,2)),-60),
  scale(mul(power(t,2),Y),135),scale(mul(t,Y),-10),scale(power(t,2),-555))
 E=mul(power(X,2),Z);J=power(mul(r,S),2);LL=mul(power(X,4),P);NN=mul(mul(power(X,3),Qpos),Z)
 Knew=mul(power(X,2),B)
 assert Knew==add(scale(LL,3),scale(E,2),scale(NN,-6))
 Gs=[t,Y,power(Y,2),mul(power(r,2),T),H4,H6,E,H5,J,LL,NN]
 jets=[[alg.bidifferential(G,z,1) for z in f] for G in Gs]
 u1=basis(cert.WEIGHTS,1);a1=cert.linear_columns(alg,f,u1)
 counts1=[sum(i==j for i,e in u1) for j in range(4)];assert counts1==[114,156,188,141]
 vectors=[cert.vector(u1,j) for j in jets];assert all(not cert.product(a1,v) for v in vectors)
 rows1=sorted(set().union(*(set(c) for c in a1)))
 det1=cert.minor(a1,{},rows1,w["h1_rank_minor"],prime);assert len(w["h1_rank_minor"]["rows"])==588
 indices=w["h1_independence_minor"]["indices"]
 assert len(indices)==len(set(indices))==11 and all(0<=i<len(u1) for i in indices)
 ind1=determinant_Q([[v[i] for v in vectors] for i in indices]);assert ind1==Q(w["h1_independence_minor"]["det_Q"]) and ind1
 # Only the two uniform cap-preserving gauges are removed.
 assert mul(r,S)==add(scale(t,-6),scale(power(Y,2),-3),scale(Y,-8))
 assert alg.bidifferential(t,P,1)==scale(Qpos,-1) and alg.bidifferential(t,Z,1)==scale(X,-1)
 assert [alg.bidifferential(mul(r,S),z,1) for z in f]==[{},r,scale(S,-1),{}]
 rawjets=[jets[j] for j in [1,3,4,5,6,7,8,9,10]]
 NR=15;zr=unit(NR);V=variables(NR)
 aa,bb,cc,dd,ee,ff,hh,jj,kk=V[:9]
 rawF1={unit(NR,i):j for i,j in enumerate(rawjets)}
 raw3=eng.rhs_coefficients(cert,alg,[{zr:f},rawF1],3)
 cone=padd(ppow(dd,2,NR),pscale(pmul(dd,ff),-250),pscale(ppow(ff,2,NR),35625))
 targets3=[cone,padd(pmul(ee,jj),pscale(ppow(kk,2,NR),Q(-1,6))),
  padd(pmul(ee,kk),pscale(pmul(jj,kk),Q(2,3)),pscale(ppow(kk,2,NR),Q(2,3))),
  pmul(ff,jj),pmul(ff,kk),
  padd(pmul(hh,jj),pscale(pmul(hh,kk),Q(1,2)),pscale(pmul(jj,kk),Q(1,36)),pscale(ppow(kk,2,NR),Q(1,72))),
  ppow(padd(jj,pscale(kk,Q(1,2))),2,NR)]
 u2=basis(cert.WEIGHTS,2);a2=cert.linear_columns(alg,f,u2);rows3=eng.rows_for(a2,raw3)
 assert [sum(i==j for i,e in u2) for j in range(4)]==[80,114,141,102]
 det2=cert.minor(a2,{},rows3,w["F2_rank_minor"],prime);assert len(w["F2_rank_minor"]["rows"])==431
 h3sup=check_functionals(eng,a2,raw3,w["h3_raw_necessity_functionals"],targets3)
 NP=14;zero=unit(NP);V=variables(NP)
 a,b,c,d,e,ff,h,z,u,v,gamma,delta,chi,psi=V
 Kjet=[alg.bidifferential(Knew,q,1) for q in f]
 reducedjets=rawjets[:7]+[Kjet]
 F0={zero:f};F1={unit(NP,i):j for i,j in enumerate(reducedjets)}
 F2=frame(cert,w["F2_field_reduction"]);red3=reduce_rhs(alg,eng.rhs_coefficients(cert,alg,[F0,F1],3))
 for m in set(F2)|set(red3):assert cert.product(a2,cert.vector(u2,F2.get(m,[{},{},{},{}])))==red3.get(m,{})
 J4=add(scale(power(Y,4),9),scale(power(Y,3),66),scale(mul(t,power(Y,2)),36),
  scale(mul(t,Y),150),scale(power(t,2),171),scale(t,-380))
 KG=[power(r,3),mul(power(r,4),S),power(X,3),mul(power(r,3),add(Y,scale(one,2))),mul(power(X,4),Qpos),mul(power(r,3),J4)]
 K2=[[alg.bidifferential(G,q,1) for q in f] for G in KG]
 assert {m for m in F2 if any(m[8:])}=={unit(NP,i) for i in range(8,14)}
 for i,j in enumerate(K2):assert F2[unit(NP,i+8)]==j and not cert.product(a2,cert.vector(u2,j))
 vec2=[cert.vector(u2,j) for j in K2];inds=w["F2_independence_minor"]["indices"]
 assert len(inds)==len(set(inds))==6 and all(0<=i<len(u2) for i in inds)
 ind2=determinant_Q([[v[i] for v in vec2] for i in inds]);assert ind2==Q(w["F2_independence_minor"]["det_Q"]) and ind2
 frames=[F0,F1,F2];full4=reduce_rhs(alg,eng.rhs_coefficients(cert,alg,frames,4))
 u3=basis(cert.WEIGHTS,3);a3=cert.linear_columns(alg,f,u3);rows4=eng.rows_for(a3,full4)
 assert [sum(i==j for i,e in u3) for j in range(4)]==[54,80,102,71]
 det3=cert.minor(a3,{},rows4,w["h4_rank_minor"],prime);assert len(w["h4_rank_minor"]["rows"])==306
 KR6=[{},{},scale(power(r,5),-6),{}];assert not cert.product(a3,cert.vector(u3,KR6))
 basic=[pmul(ppow(ff,2,NP),padd(d,pscale(ff,375))),pmul(ppow(e,2,NP),h),pmul(psi,z)]
 basic_sup=check_functionals(eng,a3,full4,w["h4_basic_functionals"],basic)
 # The next contexts are restrictions justified over FIELDS in the proof.
 L=padd(a,pscale(c,Q(-544,27)))
 K=padd(delta,pscale(ppow(c,2,NP),Q(976,27)),pscale(pmul(c,h),192),
  pscale(pmul(c,z),Q(488,9)),pscale(pmul(h,z),144),pscale(ppow(z,2,NP),33))
 tt=padd(z,pscale(h,-6));XX=padd(chi,pscale(pmul(c,z),Q(100,3)),
  pscale(pmul(z,padd(pscale(z,49),pscale(h,-288))),Q(-9,2)))
 De=padd(pscale(ppow(c,2,NP),976),pscale(pmul(c,e),216),pscale(delta,27))
 hpsi=pscale(c,Q(-22,45))
 KP=padd(delta,pscale(ppow(c,2,NP),Q(976,27)),pscale(pmul(c,hpsi),192),pscale(psi,Q(908,3)))
 contexts={
 "DF0":({3:{},5:{}},[pmul(psi,padd(pscale(c,22),pscale(h,45))),pmul(e,psi),pmul(e,padd(pscale(pmul(c,e),18),chi))]),
 "E":({3:{},5:{},6:{},7:{},13:{}},[pmul(L,De),pmul(c,De)]),
 "HZ":({3:{},5:{},4:{},13:{}},[padd(pmul(K,L),pscale(pmul(z,ppow(tt,2,NP)),576)),
  pmul(padd(pscale(c,2),pscale(h,72),pscale(z,-9)),K),padd(pscale(pmul(tt,XX),8),pscale(pmul(K,z),3))]),
 "Psi":({3:{},5:{},4:{},7:{},6:hpsi},[pmul(L,KP),padd(pmul(c,KP),pscale(pmul(L,psi),Q(-561600,414336)))]),
 "Psi_L0":({3:{},5:{},4:{},7:{},6:hpsi,0:pscale(c,Q(544,27))},[pmul(c,chi)])}
 context_controls={}
 for name,(sub,targets) in contexts.items():
  imgs=[sub.get(i,p) for i,p in enumerate(V)]
  record=w["h4_contexts"][name]
  assert [poly_rows(p) for p in record["substitution"]]==imgs
  Fs=[substitute(alg,eng,F,imgs,NP) for F in frames]
  rhs=eng.rhs_coefficients(cert,alg,Fs,4)
  assert [poly_rows(p) for p in record["targets"]]==targets
  context_controls[name]=check_functionals(eng,a3,rhs,record["functionals"],targets)
 u4=basis(cert.WEIGHTS,4);a4=cert.linear_columns(alg,f,u4)
 assert [sum(i==j for i,e in u4) for j in range(4)]==[34,54,71,47]
 branch_controls={}
 for name in ["Efull","Ezero","A0","HZres","Bzero","Cspec","Bmain"]:
  record=w["field_strata"][name];names,imgs=configurations(name);BP=len(names);zz=unit(BP)
  assert record["parameters"]==names and record["nonzero_parameters"]==(["t","K"] if name=="Bmain" else [])
  B0={zz:f};B1=frame(cert,record["F1"]);B2=frame(cert,record["F2"]);B3=frame(cert,record["F3"])
  nonzero_indices={names.index(p) for p in record["nonzero_parameters"]}
  for Bframe in [B1,B2,B3]:
   for mon in Bframe:
    assert len(mon)==BP
    assert all(exponent>=0 or index in nonzero_indices for index,exponent in enumerate(mon))
  assert all(not mon[BP-1] for Bframe in [B1,B2] for mon in Bframe)
  assert {mon for mon in B3 if mon[BP-1]}=={unit(BP,BP-1)}
  assert B3[unit(BP,BP-1)]==KR6
  assert B1==substitute(alg,eng,F1,imgs,BP) and B2==substitute(alg,eng,F2,imgs,BP)
  rhs4=eng.rhs_coefficients(cert,alg,[B0,B1,B2],4)
  assert B3[unit(BP,BP-1)]==KR6
  for m in set(B3)|set(rhs4):assert cert.product(a3,cert.vector(u3,B3.get(m,[{},{},{},{}])))==rhs4.get(m,{})
  rhs5=eng.rhs_coefficients(cert,alg,[B0,B1,B2,B3],5);VP=dict(zip(names,variables(BP)))
  if name in ["Efull","Ezero","A0","Cspec"]:targets=[{zz:Q(1)}]
  elif name=="Bzero":targets=[ppow(VP["chi"],2,BP),padd({zz:Q(1)},pscale(pmul(VP["chi"],VP["gamma"]),Q(-768,49)))]
  elif name=="Bmain":targets=[ppow(VP["K"],2,BP)]
  else:targets=[ppow(VP["h"],4,BP),ppow(padd(pscale(VP["a"],27),pscale(VP["c"],-544)),4,BP)]
  assert [poly_rows(p) for p in record["h5_targets"]]==targets
  supports=check_functionals(eng,a4,rhs5,record["h5_functionals"],targets)
  control={"parameters":names,"nonzero_parameters":record["nonzero_parameters"],"F3_blocks":len(B3),
   "F3_coefficients":sum(sum(map(len,t)) for t in B3.values()),"h5_rows":len(eng.rows_for(a4,rhs5)),
   "h5_monomials":len(rhs5),"functional_supports":supports,"all206columns_annihilated":True,
   "parameter_monomial_lengths_verified":True,"Laurent_domain_verified":True,
   "F1_F2_independent_of_w":True,"only_w_dependent_F3_block_is_w_KR6":True}
  if name=="HZres":
   subimgs=variables(BP);subimgs[0]=pscale(VP["c"],Q(544,27));subimgs[3]={}
   Fs=[substitute(alg,eng,F,subimgs,BP) for F in [B0,B1,B2,B3]]
   rr=eng.rhs_coefficients(cert,alg,Fs,5)
   targets=[ppow(VP["c"],4,BP),ppow(VP["chi"],2,BP),padd({zz:Q(1)},pscale(pmul(VP["chi"],VP["gamma"]),Q(-768,49)))]
   post=record["h5_after_h_ell_zero"];assert [poly_rows(p) for p in post["targets"]]==targets
   control["after_field_h_ell_zero_supports"]=check_functionals(eng,a4,rr,post["functionals"],targets)
   vals=tuple(Q(x) for x in [0,0,0,1,0,0,0,0,0])
   control["positive_h4_normal_PBW_control"]=eng.pbw_control(eng.evaluate_frames(alg,[B0,B1,B2,B3],vals),4)
  branch_controls[name]=control
 assert eng.rhs_indices(5,4)==[(1,3,1,Q(1)),(2,2,1,Q(1)),(3,1,1,Q(1)),
  (0,2,3,Q(1,24)),(1,1,3,Q(1,24)),(2,0,3,Q(1,24)),(0,0,5,Q(1,1920))]
 def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
 inputs=["GRADED_LONG_CAP26_FULL_JET_OBSTRUCTION.md","check_graded_long_cap26_full_jet.py",
  "graded_long_cap26_full_jet_witness.json","check_graded_long_cap23.py",
  "check_graded_long_cap23_full_jet.py","check_long_graph_curvature.py",
  "MINIMAL_EVEN_LONG_JET_AND_FIXED_LOWER_OBSTRUCTION.md"]
 hashes={name:sha(HERE/name) for name in inputs};hashes["../dixmier-audit/verify_long_pbw.py"]=sha(source)
 receipt={"claim":"No arbitrary-parity graded canonical Long scalar full jet of spatial degree26 through h5",
  "status":"PASS","environment":{"python":platform.python_version(),"arithmetic":"Fraction only; exact coefficient equations, explicit field cover, bare modular minors, rational row functionals and one independent normal-PBW positive control"},
  "cap":26,"parity":"unrestricted","prime":prime,
  "first_kernel":{"counts":counts1,"columns":599,"rank_Q":588,"kernels":11,"det_modp":det1,"independence_det_Q":str(ind1)},
  "h3":{"F2_columns":437,"rank_Q":431,"kernel_dimension":6,"det_modp":det2,"independence_det_Q":str(ind2),"necessity_supports":h3sup,"all9raw_first_parameters_kept":True,"j_plus_k_over2_square_is_field_elimination_only":True},
  "h4":{"F3_columns":307,"rank_Q":306,"kernel_dimension":1,"det_modp":det3,"basic_supports":basic_sup,"context_supports":context_controls,"field_strata":7,"both_complex_conic_signs_eliminated_by_equations":True},
  "h5":{"F4_columns":206,"Pi3_F1_F1_included":True,"all_seven_field_families":branch_controls},
  "scope":{"characteristic_zero_fields_only":True,"all_first_images_allowed":True,"only_uniform_t_RS_gauges":True,"no_evenization_assumed":True,"no_ungraded_bound_claim":True},
  "input_sha256":hashes}
 out=args.output_dir or Path(tempfile.mkdtemp(prefix="graded-long-cap26-full-"));out.mkdir(parents=True,exist_ok=True)
 dest=out/"graded_long_cap26_full_jet_controls.json";dest.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
 print(json.dumps({"status":"PASS","output":str(dest),"input_sha256":hashes}))
if __name__=="__main__":main()
