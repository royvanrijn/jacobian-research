#!/usr/bin/env python3
"""Complete cap26 matrix family with prescribed Airy first coefficient.
Requires python-flint; all parameter/functionals are verified over Fraction.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,importlib.util,json,platform,tempfile
import flint
HERE=Path(__file__).resolve().parent
def load(path,name):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def basis(weights,m):
 out=[]
 for i,w in enumerate(weights):
  for b in range(27):
   for c in range(27-b):
    for d in range(27-b-c):
     a=b+4*c+2*d+3*m-w
     if a>=0 and a+b+c+d<=26:out.append((i,(a,b,c,d)))
 return sorted(out)
def functional(data,rows):
 assert data["row_count"]==len(rows);out={}
 for index,key,c in data["functional"]:
  k=(key[0],tuple(key[1])) if len(key)==2 else (key[0],key[1],tuple(key[2]))
  assert rows[index]==k and k not in out;out[k]=Q(c)
 return out
def poly_add(a,b,c=Q(1)):
 out=dict(a)
 for key,d in b.items():out[key]=out.get(key,Q(0))+c*d
 return {key:d for key,d in out.items() if d}
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--output-dir",type=Path);args=ap.parse_args()
 eng=load(HERE/"check_airy_first_coefficient_cap25.py","matrix_covariant_arithmetic")
 alg=load(HERE/"check_long_graph_curvature.py","independent_Fraction_polynomials");eng.ALG=alg
 cert=load(HERE/"check_graded_long_cap23.py","classical_coefficient_engine")
 source=HERE.parent/"dixmier-audit/verify_long_pbw.py";prod=load(source,"canonical_Long_input")
 f=[{tuple(e):Q(c) for e,c in poly.items()} for poly in prod.long_images()]
 wp=HERE/"airy_first_coefficient_cap26_witness.json";w=json.loads(wp.read_text())
 pp=HERE/w["fixed_point_witness"];point=json.loads(pp.read_text())
 oldp=HERE/w["scalar_F2_witness"];old=json.loads(oldp.read_text())
 add,mul,scale=alg.add,alg.multiply,alg.scale
 madd,mscale,bidiff,diag=eng.madd,eng.mscale,eng.bidiff,eng.diag
 F0=[diag(p) for p in f];F1=[eng.parse(p) for p in point["point"]["F1"]]
 fixed2=[eng.parse(p) for p in point["point"]["F2"]];fixed3=[eng.parse(p) for p in point["point"]["F3"]]
 X={(1,0,0,0):Q(1)};Qp={(0,1,0,0):Q(1)};P={(0,0,1,0):Q(1)};Z={(0,0,0,1):Q(1)}
 one={(0,0,0,0):Q(1)};t=mul(X,Qp);r=f[0];power=lambda p,n:cert.power(alg,p,n)
 M={(0,1):one,(1,0):t};phi=power(X,3)
 V=lambda p:add(mul(Qp,alg.derivative(p,(0,0,1,0))),mul(X,alg.derivative(p,(0,0,0,1))))
 alpha=[V(fi) for fi in f]
 for a,actual in zip(alpha,F1):assert not madd(actual,mscale({ij:mul(p,a) for ij,p in M.items()},-1))
 B=add(scale(mul(power(X,2),P),3),scale(mul(add(one,scale(t,-3)),Z),2))
 A0=add(Qp,scale(mul(X,power(Qp,2)),3),scale(mul(X,B),Q(-1,3)))
 Y=add(mul(r,A0),scale(t,-3))
 J4=add(scale(power(Y,4),9),scale(power(Y,3),66),scale(mul(t,power(Y,2)),36),
         scale(mul(t,Y),150),scale(power(t,2),171),scale(t,-380))
 scalarG=[power(r,3),mul(power(r,4),f[3]),phi,mul(power(r,3),add(Y,scale(one,2))),
          mul(power(X,4),Qp),mul(power(r,3),J4)]
 SK=[[alg.bidifferential(G,fi,1) for fi in f] for G in scalarG]
 def cov(i,C):return madd(bidiff(F0[i],C,1),bidiff(F1[i],C,0),mscale(bidiff(C,F1[i],0),-1))
 Gs=[{(1,0):phi},{(0,1):phi,(1,0):mul(t,phi)},
     {(0,0):phi,(1,1):scale(phi,-1),(0,1):scale(mul(phi,t),2),(1,0):scale(mul(phi,power(t,2)),2)}]
 TK=[[cov(i,G) for i in range(4)] for G in Gs]
 p=1000003
 u2=basis(cert.WEIGHTS,2);A2sc=cert.linear_columns(alg,f,u2)
 assert len(u2)==437
 sv=[cert.vector(u2,k) for k in SK]
 for z in sv:assert not cert.product(A2sc,z)
 piv={}
 for row in zip(*sv):
  z={i:c for i,c in enumerate(row) if c}
  while z:
   j=min(z)
   if j not in piv:
    fac=z[j];piv[j]={k:c/fac for k,c in z.items()};break
   fac=z[j];z={k:z.get(k,Q(0))-fac*piv[j].get(k,Q(0)) for k in set(z)|set(piv[j])}
   z={k:c for k,c in z.items() if c}
 assert len(piv)==6
 old_vectors=[]
 for rows in old["first_family"]["kernels"]:
  tup=[{} for _ in range(4)]
  for label,e,c in rows:tup[cert.LABELS.index(label)][tuple(e)]=Q(c)
  old_vectors.append(cert.vector(u2,tup))
 scalar_free=old["first_family"]["free_columns"]
 assert [[z[j] for j in scalar_free] for z in old_vectors]==[
        [Q(int(i==j)) for j in range(6)] for i in range(6)]
 transition=[[z[j] for j in scalar_free] for z in sv]
 assert [[str(c) for c in row] for row in transition]==w["scalar_compact_to_RREF"]
 for actual,coords in zip(sv,transition):
  assert actual==[sum((c*z[j] for c,z in zip(coords,old_vectors)),Q(0)) for j in range(len(u2))]

 dsc2=eng.minor(A2sc,old["F2_rank_minor"],p)
 assert len(old["F2_rank_minor"]["rows"])==431
 A2sl=eng.sl2_columns(f,alpha,u2);assert len(A2sl)==1311
 dsl2=eng.minor(A2sl,w["traceless_homogeneous_minors"]["2"]["minor"],p)
 assert len(w["traceless_homogeneous_minors"]["2"]["minor"]["rows"])==1308
 lookup={key:i for i,key in enumerate(u2)};tv=[]
 for K in TK:
  z=[Q(0)]*1311
  for image,mat in enumerate(K):
   assert not add(mat.get((0,0),{}),mat.get((1,1),{}))
   for kind,ij in enumerate(((0,0),(0,1),(1,0))):
    for e,c in mat.get(ij,{}).items():z[3*lookup[image,e]+kind]=c
  assert not cert.product(A2sl,z);tv.append(z)
 free=w["traceless_kernel_independence_coordinates"];assert free==[1056,1057,1058]
 assert [[z[j] for j in free] for z in tv]==[[Q(-18),Q(0),Q(0)],[Q(0),Q(0),Q(36)],[Q(0),Q(108),Q(72)]]
 # Exact independent kernels give the opposite rank bound, hence nullity3.
 u3=basis(cert.WEIGHTS,3);A3sc=cert.linear_columns(alg,f,u3);A3sl=eng.sl2_columns(f,alpha,u3)
 assert len(u3)==307 and len(A3sl)==921
 dsc3=eng.minor(A3sc,w["scalar_F3_minor"],p);assert len(w["scalar_F3_minor"]["rows"])==306
 KR6=[{},{},scale(power(r,5),-6),{}]
 assert not cert.product(A3sc,cert.vector(u3,KR6))
 dsl3=eng.minor(A3sl,w["traceless_homogeneous_minors"]["3"]["minor"],p)
 assert len(w["traceless_homogeneous_minors"]["3"]["minor"]["rows"])==921
 # Complete h4 RHS: linearity in F2 plus ordinary commutators of the
 # traceless F2 parameters. Scalar F2 matrices commute pointwise.
 zero9=(0,)*9
 base=eng.parameter_rhs([{zero9:F0},{zero9:F1},{zero9:fixed2}],4).get(zero9,{})
 rhs4={zero9:base}
 def flatten(pi,mat):return {(pi,a,b,e):c for (a,b),poly in mat.items() for e,c in poly.items()}
 def linear_h4(K):
  out={}
  for pi,(i,j) in enumerate(cert.PAIRS):
   mat=mscale(madd(bidiff(F1[i],K[j],1),bidiff(K[i],F1[j],1),
                       mscale(bidiff(F1[j],K[i],1),-1),mscale(bidiff(K[j],F1[i],1),-1)),Q(-1,2))
   out.update(flatten(pi,mat))
  return out
 for j,K in enumerate([[diag(z) for z in k] for k in SK]+TK):
  value=linear_h4(K)
  if value:rhs4[tuple(int(i==j) for i in range(9))]=value
 for a in range(3):
  for b in range(a,3):
   value={}
   for pi,(i,j) in enumerate(cert.PAIRS):
    mat=madd(bidiff(TK[a][i],TK[b][j],0),mscale(bidiff(TK[b][j],TK[a][i],0),-1))
    if a!=b:mat=madd(mat,bidiff(TK[b][i],TK[a][j],0),mscale(bidiff(TK[a][j],TK[b][i],0),-1))
    value.update(flatten(pi,mscale(mat,-1)))
   if value:
    mon=[0]*9;mon[6+a]+=1;mon[6+b]+=1;rhs4[tuple(mon)]=value
 assert len(rhs4)==15 and sum(map(len,rhs4.values()))==7448
 trace={};sl={}
 for mon,poly in rhs4.items():
  tr={};ss={}
  for (pi,a,b,e),c in poly.items():
   if a==b:
    tr[pi,e]=tr.get((pi,e),Q(0))+c/2
    ss[pi,0,e]=ss.get((pi,0,e),Q(0))+(c/2 if a==0 else -c/2)
   else:ss[pi,1 if (a,b)==(0,1) else 2,e]=c
  tr={e:c for e,c in tr.items() if c};ss={e:c for e,c in ss.items() if c}
  if tr:trace[mon]=tr
  if ss:sl[mon]=ss
 rows_sc=sorted(set().union(*(set(c) for c in A3sc),*(set(b) for b in trace.values())))
 s3left=functional(w["h4_scalar_s3"]["functional"],rows_sc)
 assert all(not eng.apply(s3left,c) for c in A3sc)
 s3mon=tuple(int(i==8) for i in range(9))
 assert {m:eng.apply(s3left,b) for m,b in trace.items() if eng.apply(s3left,b)}=={s3mon:Q(1)}
 reduced={m:b for m,b in sl.items() if m[8]==0}
 rows_sl=sorted(set().union(*(set(c) for c in A3sl),*(set(b) for b in reduced.values())))
 wanted={name:tuple(int(i==j) for i in range(9)) for name,j in (("e",3),("eps",4),("zeta",5))}
 wanted.update({name:tuple(2*int(i==j) for i in range(9)) for name,j in (("s1_square",6),("s2_square",7))})
 supports={}
 for name,mon in wanted.items():
  left=functional(w["h4_after_s3_zero"]["functionals"][name],rows_sl)
  assert all(not eng.apply(left,c) for c in A3sl)
  assert {m:eng.apply(left,b) for m,b in reduced.items() if eng.apply(left,b)}=={mon:Q(1)}
  supports[name]=len(left)
 # Entire surviving family u,v,gamma,w; no gauge removes parameters.
 zero=(0,)*4
 times=lambda mat,poly:{ij:q for ij,p in mat.items() if (q:=mul(p,poly))}
 frames=[{zero:F0},{zero:F1},
  {zero:fixed2,(1,0,0,0):[diag(z) for z in SK[0]],(0,1,0,0):[diag(z) for z in SK[1]],(0,0,1,0):[diag(z) for z in SK[2]]},
  {zero:fixed3,(0,1,0,0):[times(M,V(z)) for z in SK[1]],(0,0,1,0):[times(M,V(z)) for z in SK[2]],(0,0,0,1):[diag(z) for z in KR6]}]
 for order in (3,4):assert not eng.parameter_rhs(frames[:3] if order==3 else frames,order)
 for mon,tup in frames[3].items():
  for i,mat in enumerate(tup):
   for poly in mat.values():
    for e in poly:assert (i,e) in set(u3)
 rhs5=eng.parameter_rhs(frames,5);mons=sorted(rhs5,key=lambda m:(sum(m),m))
 assert [list(m) for m in mons]==w["h5"]["monomials"]
 scalar5={}
 for m,poly in rhs5.items():
  z={}
  for (pi,a,b,e),c in poly.items():
   assert a==b
   if a==0:z[pi,e]=c
  for (pi,e),c in z.items():assert poly.get((pi,1,1,e),Q(0))==c
  scalar5[m]=z
 assert [len(scalar5[m]) for m in mons]==[201,273,96,301,138]
 u4=basis(cert.WEIGHTS,4);A4=cert.linear_columns(alg,f,u4);assert len(u4)==206
 rows5=sorted(set().union(*(set(c) for c in A4),*(set(b) for b in scalar5.values())))
 last=functional(w["h5"]["functional"],rows5)
 assert all(not eng.apply(last,c) for c in A4)
 assert [eng.apply(last,scalar5[m]) for m in mons]==[Q(1),Q(0),Q(0),Q(0),Q(0)]
 lift={}
 for (pi,e),c in last.items():lift[pi,0,0,e]=c/2;lift[pi,1,1,e]=c/2
 count=0
 for component,e in u4:
  for a,b in ((0,0),(0,1),(1,0),(1,1)):
   C={(a,b):{e:Q(1)}};col={}
   for pi,(i,j) in enumerate(cert.PAIRS):
    value=cov(i,C) if component==j else (mscale(cov(j,C),-1) if component==i else {})
    col.update(flatten(pi,value))
   assert eng.apply(lift,col)==0;count+=1
 assert count==824 and len(last)==158 and len(lift)==316
 assert [eng.apply(lift,rhs5[m]) for m in mons]==[Q(1),Q(0),Q(0),Q(0),Q(0)]
 bad=dict(lift);key=next(iter(rhs5[zero]));bad[key]=bad.get(key,Q(0))+1;assert eng.apply(bad,rhs5[zero])!=1
 def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
 receipt={"claim":"No graded cap26 2x2 Long matrix jet through h5 with prescribed F1=M(XQ)Vf",
  "status":"PASS","environment":{"python":platform.python_version(),"python_flint":flint.__version__,
   "arithmetic":"Fraction identities and exact FLINT modular determinants"},
  "cap":26,"F1_fixed":True,"all_later_entries_images_and_parameters_retained":True,
  "F2_ranks":{"scalar":"431/437","sl2":"1308/1311","new_exact_traceless_kernels":3},
  "scalar_compact_to_RREF":[[str(c) for c in row] for row in transition],
  "scalar_parameter_order":["u_R3","v_R4S","gamma_X3","delta_R4a=e-2/27","eps_X4Q","zeta_R3J4"],
  "F3_ranks":{"scalar":"306/307","sl2":"921/921"},
  "minor_determinants":{"scalar2":dsc2,"sl2_2":dsl2,"scalar3":dsc3,"sl2_3":dsl3},
  "h4":{"all9_first_parameters_kept":True,"s3_functional_rows":len(s3left),
    "after_s3_zero_functional_rows":supports,
    "exact_field_conclusions":"s1=s2=s3=0;e=eps=zeta=0",
    "scope":"Individual square identities over fields; no complex sum-of-squares inference"},
  "whole_h4_family":"F2=point2+uKR3I+vKR4SI+gammaKX3I;F3=point3+vMV(KR4S)+gammaMV(KX3)+wKR6I",
  "h5":{"columns":count,"scalar_functional_rows":len(last),"matrix_rows":len(lift),
    "parameter_order":["u","v","gamma","w"],"RHS_monomials":[list(m) for m in mons],
    "exact_values":["1","0","0","0","0"],"conclusion":"0=1 for every remaining parameter point"},
  "invalid_control":{"corrupted_final_functional_rejected":True},
  "hashes":{"AIRY_FIRST_COEFFICIENT_CAP26_OBSTRUCTION.md":sha(HERE/"AIRY_FIRST_COEFFICIENT_CAP26_OBSTRUCTION.md"),
   "check_airy_first_coefficient_cap26.py":sha(Path(__file__)),"airy_first_coefficient_cap26_witness.json":sha(wp),
   "check_airy_first_coefficient_cap25.py":sha(HERE/"check_airy_first_coefficient_cap25.py"),
   "airy_long_cap25_matrix_jet_witness.json":sha(pp),"graded_long_cap26_even_witness.json":sha(oldp),
   "check_graded_long_cap23.py":sha(HERE/"check_graded_long_cap23.py"),"check_long_graph_curvature.py":sha(HERE/"check_long_graph_curvature.py"),
   "../dixmier-audit/verify_long_pbw.py":sha(source)},
  "scope":"Complete graded cap26 family with fixed F1 only; no other connection, scaling, first coefficient, larger cap, nongraded, bounded gauge or DC1/DC2 conclusion",
  "independent_replay":False,"formal_verification":False,"external_review":False}
 out=args.output_dir or Path(tempfile.mkdtemp(prefix="airy-first-cap26-"));out.mkdir(parents=True,exist_ok=True)
 dest=out/"airy_first_coefficient_cap26_controls.json";dest.write_text(json.dumps(receipt,indent=2)+"\n")
 print("PASS: complete Airy-F1 cap26 family; all three new traceless freedoms obstructed; h5 gives0=1")
 print("Receipt:",dest)
if __name__=="__main__":main()
