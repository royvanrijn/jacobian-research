#!/usr/bin/env python3
"""Complete cap25 matrix family with prescribed Airy first coefficient.
Requires python-flint for the two large exact modular determinant controls.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial,isqrt
import argparse,hashlib,importlib.util,json,platform,tempfile
import flint
from flint import nmod_mat
HERE=Path(__file__).resolve().parent
ALG=None;PAIRS=((0,1),(0,2),(0,3),(1,2),(1,3),(2,3))
def load(path,name):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def madd(*ms):
 out={}
 for mat in ms:
  for ij,p in mat.items():out[ij]=ALG.add(out.get(ij,{}),p)
 return {ij:p for ij,p in out.items() if p}
def mscale(m,c):return {ij:q for ij,p in m.items() if (q:=ALG.scale(p,c))}
def diag(p):return {(0,0):p,(1,1):p} if p else {}
def bidiff(a,b,n):
 out={}
 for (i,j),p in a.items():
  for (k,l),q in b.items():
   if j==k:out[i,l]=ALG.add(out.get((i,l),{}),ALG.bidifferential(p,q,n))
 return {ij:p for ij,p in out.items() if p}
def parse(rows):
 out={}
 for a,b,e,c in rows:
  e=tuple(e);assert e not in out.setdefault((a,b),{});out[a,b][e]=Q(c)
 return out
def basis(weights,m):
 out=[]
 for i,w in enumerate(weights):
  for b in range(26):
   for c in range(26-b):
    for d in range(26-b-c):
     a=b+4*c+2*d+3*m-w
     if a>=0 and a+b+c+d<=25:out.append((i,(a,b,c,d)))
 return sorted(out)
def parameter_rhs(frames,n):
 out={}
 for pi,(i,j) in enumerate(PAIRS):
  for a in range(len(frames)):
   for b in range(len(frames)):
    ell=n-a-b
    if ell<0:continue
    for ma,fa in frames[a].items():
     for mb,fb in frames[b].items():
      mon=tuple(x+y for x,y in zip(ma,mb))
      value=mscale(madd(bidiff(fa[i],fb[j],ell),mscale(bidiff(fa[j],fb[i],ell),-1)),
                   Q(-1,2**ell*factorial(ell)))
      dest=out.setdefault(mon,{})
      for (aa,bb),p in value.items():
       for e,c in p.items():
        key=(pi,aa,bb,e);dest[key]=dest.get(key,Q(0))+c
 return {m:{key:c for key,c in p.items() if c} for m,p in out.items() if any(p.values())}
def traceless_projection(rhs):
 out={}
 for mon,p in rhs.items():
  val={}
  for (pi,a,b,e),c in p.items():
   if a==b:val[(pi,0,e)]=val.get((pi,0,e),Q(0))+(c/2 if a==0 else -c/2)
   else:val[(pi,1 if (a,b)==(0,1) else 2,e)]=c
  trace={}
  for (pi,a,b,e),c in p.items():
   if a==b:trace[pi,e]=trace.get((pi,e),Q(0))+c
  assert not any(trace.values())
  val={key:c for key,c in val.items() if c}
  if val:out[mon]=val
 return out
def sl2_columns(f,alpha,u):
 out=[]
 def shift(poly,e,fac=Q(1)):return {tuple(x+y for x,y in zip(k,e)):fac*c for k,c in poly.items()}
 for component,e in u:
  for kind in range(3):
   col={}
   for pi,(i,j) in enumerate(PAIRS):
    if component not in (i,j):continue
    other=i if component==j else j;sign=Q(1) if component==j else Q(-1)
    diff=ALG.scale(ALG.bidifferential(f[other],{e:Q(1)},1),sign)
    ac=shift(alpha[other],e,sign);tac=shift(ac,(1,1,0,0));terms={kind:diff}
    if kind==0:terms[1]=ALG.scale(ac,-2);terms[2]=ALG.scale(tac,2)
    elif kind==1:terms[0]=ALG.scale(tac,-1)
    else:terms[0]=ac
    for kk,p in terms.items():
     for ee,c in p.items():col[pi,kk,ee]=c
   out.append(col)
 return out
def minor(columns,data,p,rows=None):
 rows=rows or sorted(set().union(*(set(c) for c in columns)))
 assert data["shape"]==[len(rows),len(columns)]
 rr,cc=data["rows"],data["columns"]
 assert len(rr)==len(cc) and len(set(rr))==len(rr) and len(set(cc))==len(cc)
 def mod(c):
  assert c.denominator%p
  return c.numerator%p*pow(c.denominator%p,-1,p)%p
 selected=[[mod(columns[j].get(rows[i],Q(0))) for j in cc] for i in rr]
 det=int(nmod_mat(selected,p).det());assert det==data["det_modp"] and det
 return det
def functional(data,rows):
 assert data["row_count"]==len(rows)
 out={}
 for index,key,c in data["functional"]:
  key=(key[0],key[1],tuple(key[2]));assert rows[index]==key and key not in out;out[key]=Q(c)
 return out
def apply(left,p):return sum((c*p.get(key,Q(0)) for key,c in left.items()),Q(0))
def main():
 global ALG
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--output-dir",type=Path);args=ap.parse_args()
 ALG=load(HERE/"check_long_graph_curvature.py","independent_polynomial_arithmetic")
 cert=load(HERE/"check_graded_long_cap23.py","classical_coefficient_engine")
 source=HERE.parent/"dixmier-audit/verify_long_pbw.py";prod=load(source,"canonical_Long_input")
 f=[{tuple(e):Q(c) for e,c in poly.items()} for poly in prod.long_images()]
 wp=HERE/"airy_first_coefficient_cap25_witness.json";w=json.loads(wp.read_text())
 point_path=HERE/w["fixed_point_witness"];point=json.loads(point_path.read_text())
 old_path=HERE/w["scalar_F2_witness"];old=json.loads(old_path.read_text())
 F0=[diag(p) for p in f];F1=[parse(rows) for rows in point["point"]["F1"]]
 fixed2=[parse(rows) for rows in point["point"]["F2"]];fixed3=[parse(rows) for rows in point["point"]["F3"]]
 X={(1,0,0,0):Q(1)};Qp={(0,1,0,0):Q(1)};one={(0,0,0,0):Q(1)}
 t=ALG.multiply(X,Qp);M={(0,1):one,(1,0):t};r=f[0]
 power=lambda p,n:cert.power(ALG,p,n)
 V=lambda p:ALG.add(ALG.multiply(Qp,ALG.derivative(p,(0,0,1,0))),ALG.multiply(X,ALG.derivative(p,(0,0,0,1))))
 alpha=[V(fi) for fi in f]
 for actual,a in zip(F1,alpha):
  assert not madd(actual,mscale({ij:ALG.multiply(p,a) for ij,p in M.items()},-1))
 B=ALG.add(ALG.scale(ALG.multiply(power(X,2),{(0,0,1,0):Q(1)}),3),
           ALG.scale(ALG.multiply(ALG.add(one,ALG.scale(t,-3)),{(0,0,0,1):Q(1)}),2))
 A0=ALG.add(Qp,ALG.scale(ALG.multiply(X,power(Qp,2)),3),ALG.scale(ALG.multiply(X,B),Q(-1,3)))
 Y=ALG.add(ALG.multiply(r,A0),ALG.scale(t,-3))
 G4=ALG.multiply(power(r,3),ALG.add(Y,ALG.scale(one,2)))
 K=[
  [{},{},ALG.scale(power(r,2),-3),{}],
  [{},power(r,4),ALG.scale(ALG.multiply(power(r,3),f[3]),-4),{}],
  [ALG.bidifferential(power(X,3),fi,1) for fi in f],
  [ALG.bidifferential(G4,fi,1) for fi in f]]
 p=1000003;assert all(p%d for d in range(2,isqrt(p)+1))
 u2=basis(cert.WEIGHTS,2);a2=cert.linear_columns(ALG,f,u2)
 assert len(u2)==391
 # Exact independent scalar kernel vectors and the retained387minor.
 vectors=[cert.vector(u2,k) for k in K]
 for z in vectors:assert not cert.product(a2,z)
 piv={}
 for row in zip(*vectors):
  z={i:c for i,c in enumerate(row) if c}
  while z:
   j=min(z)
   if j not in piv:
    fac=z[j];piv[j]={k:c/fac for k,c in z.items()};break
   fac=z[j];z={k:z.get(k,Q(0))-fac*piv[j].get(k,Q(0)) for k in set(z)|set(piv[j])}
   z={k:c for k,c in z.items() if c}
 assert len(piv)==4
 b3=cert.lower_rhs(ALG,[f],1)
 detscalar2=minor(a2,old["h3_minor"],p,sorted(set(b3).union(*(set(c) for c in a2))))
 assert len(old["h3_minor"]["rows"])==387
 # Rebuild both large covariant operators; no stored rank flag is trusted.
 a2sl=sl2_columns(f,alpha,u2);assert len(a2sl)==1173
 detsl2=minor(a2sl,w["traceless_homogeneous_minors"]["2"]["minor"],p)
 assert len(w["traceless_homogeneous_minors"]["2"]["minor"]["rows"])==1173
 u3=basis(cert.WEIGHTS,3);a3sl=sl2_columns(f,alpha,u3);assert len(u3)==271 and len(a3sl)==813
 detsl3=minor(a3sl,w["traceless_homogeneous_minors"]["3"]["minor"],p)
 assert len(w["traceless_homogeneous_minors"]["3"]["minor"]["rows"])==813
 a3=cert.linear_columns(ALG,f,u3);KR6=[{},{},ALG.scale(power(r,5),-6),{}]
 assert not cert.product(a3,cert.vector(u3,KR6))
 detscalar3=minor(a3,w["scalar_F3_minor"],p);assert len(w["scalar_F3_minor"]["rows"])==270
 # Entire F2 family centered at the positive point, parameters u,v,gamma,e.
 zero4=(0,)*4;Full2={zero4:fixed2}
 for j,k in enumerate(K):Full2[tuple(int(i==j) for i in range(4))]=[diag(z) for z in k]
 frames=[{zero4:F0},{zero4:F1},Full2]
 assert not parameter_rhs(frames,3)
 h4=traceless_projection(parameter_rhs(frames,4))
 rows3=sorted(set().union(*(set(c) for c in a3sl),*(set(b) for b in h4.values())))
 assert len(rows3)==6677
 supports={}
 for name,index in (("gamma_X3",2),("e_R4a",3)):
  left=functional(w["F3_compatibility"][name],rows3)
  assert all(not apply(left,col) for col in a3sl)
  want={tuple(int(i==index) for i in range(4)):Q(1)}
  assert {m:apply(left,b) for m,b in h4.items() if apply(left,b)}==want
  supports[name]=len(left)
 # Entire surviving F2/F3 families, now parameters u,v,w.
 zero=(0,0,0);vpart=[{ij:ALG.multiply(p,V(ki)) for ij,p in M.items()} for ki in K[1]]
 frames=[{zero:F0},{zero:F1},
  {zero:fixed2,(1,0,0):[diag(z) for z in K[0]],(0,1,0):[diag(z) for z in K[1]]},
  {zero:fixed3,(0,1,0):vpart,(0,0,1):[diag(z) for z in KR6]}]
 assert not parameter_rhs(frames,4)
 rhs=parameter_rhs(frames,5)
 assert set(rhs)=={zero,(0,1,0)} and len(rhs[zero])==402 and len(rhs[0,1,0])==192
 left={(pi,a,b,tuple(e)):Q(c) for pi,a,b,e,c in point["fixed_h5"]["functional"]}
 assert len(left)==274 and apply(left,rhs[zero])==1 and apply(left,rhs[0,1,0])==0
 # Separate corollary: displayed scaled lower family for ALL lambda,u,v,w.
 # No classification of additional kernels at arbitrarylambda is assumed.
 first_path=HERE/old["first_particular_source"];first=json.loads(first_path.read_text())
 g0=[{} for _ in range(4)]
 for name,e,c in first["first_particular"]:g0[cert.LABELS.index(name)][tuple(e)]=Q(c)
 H={(0,0):one,(1,1):ALG.scale(one,-1)};Mp={(1,0):one}
 def times(mat,poly):
  return {ij:q for ij,p in mat.items() if (q:=ALG.multiply(p,poly))}
 scaled=[
  {(0,0,0,0):F0},
  {(1,0,0,0):F1},
  {(0,0,0,0):[diag(z) for z in g0],
   (2,0,0,0):[diag(ALG.add(ALG.scale(ki,Q(-2,27)),ALG.scale(ALG.multiply(t,V(V(fi))),Q(1,2))))
               for fi,ki in zip(f,K[3])],
   (0,1,0,0):[diag(z) for z in K[0]],(0,0,1,0):[diag(z) for z in K[1]]},
  {(1,0,0,0):[madd(times(M,V(gi)),mscale(times(Mp,V(ALG.derivative(fi,(0,0,1,1)))),Q(1,4)))
               for fi,gi in zip(f,g0)],
   (2,0,0,0):[mscale(times(H,V(V(V(fi)))),Q(1,12)) for fi in f],
   (3,0,0,0):[madd(mscale(times(M,ALG.multiply(t,V(V(V(fi))))),Q(1,6)),
                         mscale(times(M,V(ki)),Q(-2,27))) for fi,ki in zip(f,K[3])],
   (1,0,1,0):[times(M,V(ki)) for ki in K[1]],
   (0,0,0,1):[diag(z) for z in KR6]}]
 assert not parameter_rhs(scaled[:3],3) and not parameter_rhs(scaled,4)
 srhs=parameter_rhs(scaled,5);smons=sorted(srhs,key=lambda m:(sum(m),m))
 sd=w["scaled_displayed_family"]
 assert [list(m) for m in smons]==sd["RHS_monomials"]
 assert [len(srhs[m]) for m in smons]==sd["RHS_full_supports"]==[192,138,300,192,390]
 scalar_left={(pi,tuple(e)):Q(c) for pi,e,c in sd["exact_trace_functional"]}
 assert len(scalar_left)==139
 scaled_left={}
 for (pi,e),c in scalar_left.items():
  scaled_left[pi,0,0,e]=c/2;scaled_left[pi,1,1,e]=c/2
 assert len(scaled_left)==278
 assert [apply(scaled_left,srhs[m]) for m in smons]==[Q(1),Q(0),Q(0),Q(0),Q(0)]
 # The scalar trace functional annihilates the ordinary scalar derivative map.
 a4_scalar=cert.linear_columns(ALG,f,basis(cert.WEIGHTS,4))
 assert all(not apply(scalar_left,col) for col in a4_scalar)
 # Build all712 full matrix F4columns and test the exact final functional.
 u4=basis(cert.WEIGHTS,4);assert len(u4)==178
 def cov(i,C):return madd(bidiff(diag(f[i]),C,1),bidiff(F1[i],C,0),mscale(bidiff(C,F1[i],0),-1))
 count=0
 for component,e in u4:
  for a,b in ((0,0),(0,1),(1,0),(1,1)):
   C={(a,b):{e:Q(1)}};col={}
   for pi,(i,j) in enumerate(PAIRS):
    value=cov(i,C) if component==j else (mscale(cov(j,C),-1) if component==i else {})
    for (aa,bb),poly in value.items():
     for ee,c in poly.items():col[pi,aa,bb,ee]=c
   assert apply(left,col)==0 and apply(scaled_left,col)==0
   # Check the lambda^0 operator too. The difference is the lambda^1
   # commutator coefficient, so BOTH coefficients vanish for every lambda.
   col0={}
   for pi,(i,j) in enumerate(PAIRS):
    value=bidiff(diag(f[i]),C,1) if component==j else (
        mscale(bidiff(diag(f[j]),C,1),-1) if component==i else {})
    for (aa,bb),poly in value.items():
     for ee,c in poly.items():col0[pi,aa,bb,ee]=c
   assert apply(scaled_left,col0)==0
   count+=1
 assert count==712
 bad=dict(left);key=next(iter(rhs[zero]));bad[key]=bad.get(key,Q(0))+1;assert apply(bad,rhs[zero])!=1
 def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
 receipt={"claim":"No graded cap25 2x2 Long matrix jet through h5 with prescribed Airy F1=M(XQ)*(Q*dP+X*dZ)f",
  "status":"PASS","environment":{"python":platform.python_version(),"python_flint":flint.__version__,
    "arithmetic":"All coefficient and parameter identities Fraction; FLINT exact modular determinant replay for rank minors, none skipped"},
  "prescribed_first_coefficient":True,"all_later_matrix_entries_and_images_allowed":True,
  "scalar_F2_rank":"387/391","traceless_F2_rank":"1173/1173","scalar_F3_rank":"270/271","traceless_F3_rank":"813/813",
  "minor_determinants":{"scalar_F2":detscalar2,"traceless_F2":detsl2,"scalar_F3":detscalar3,"traceless_F3":detsl3},
  "h4_necessity_functional_supports":supports,"h4_conditions":"gamma=0;delta=-2/27",
  "entire_h4_family":"F2=point2+u*K_R3*I+v*K_R4S*I;F3=point3+v*M*V(K_R4S)+w*K_R6*I",
  "whole_family_h3_h4_verified_coefficientwise":True,
  "h5":{"all_columns":count,"functional_support":len(left),"RHS_coefficients":["constant402terms","v192terms"],
        "exact_functional_values":["1","0"],"u_w_cancel_without_removal":True,"conclusion":"0=1 for all complex or characteristic-zero u,v,w"},
  "scaled_displayed_family_corollary":{"parameters":["lambda","u_R3","v_R4S","w_R6"],
    "all_h3_h4_parameter_identities_verified":True,
    "h5_RHS_monomials":[list(m) for m in smons],"h5_RHS_supports":[len(srhs[m]) for m in smons],
    "scalar_functional_rows":139,"matrix_trace_lift_rows":278,
    "all_712_covariant_columns_annihilated_in_both_lambda_coefficients":True,
    "exact_values":["1","0","0","0","0"],
    "scope":"Only the displayed scaled F2/F3 family, not all lower kernels or all jets for arbitrarylambda"},
  "invalid_control":{"corrupted_final_functional_rejected":True},
  "hashes":{"AIRY_FIRST_COEFFICIENT_CAP25_OBSTRUCTION.md":sha(HERE/"AIRY_FIRST_COEFFICIENT_CAP25_OBSTRUCTION.md"),
   "check_airy_first_coefficient_cap25.py":sha(Path(__file__)),"airy_first_coefficient_cap25_witness.json":sha(wp),
   "AIRY_LONG_CAP25_MATRIX_JET_AND_FIXED_EXTENSION_OBSTRUCTION.md":sha(HERE/"AIRY_LONG_CAP25_MATRIX_JET_AND_FIXED_EXTENSION_OBSTRUCTION.md"),
   "airy_long_cap25_matrix_jet_witness.json":sha(point_path),"graded_long_cap25_witness.json":sha(old_path),
   "graded_long_cap23_witness.json":sha(first_path),
   "check_graded_long_cap23.py":sha(HERE/"check_graded_long_cap23.py"),"check_long_graph_curvature.py":sha(HERE/"check_long_graph_curvature.py"),
   "../dixmier-audit/verify_long_pbw.py":sha(source)},
  "scope":"Complete graded cap25 family with this fixed matrix F1, not all matrix first coefficients, connections, larger caps or noncompositional graph carriers; no bounded formal-gauge or DC1/DC2 claim",
  "independent_replay":False,"formal_verification":False,"external_review":False}
 out=args.output_dir or Path(tempfile.mkdtemp(prefix="airy-first-cap25-"));out.mkdir(parents=True,exist_ok=True)
 dest=out/"airy_first_coefficient_cap25_controls.json";dest.write_text(json.dumps(receipt,indent=2)+"\n")
 print("PASS: complete lambda1 Airy-F1 cap25 family; displayed scaled family fails for every lambda,u,v,w")
 print("Receipt:",dest)
if __name__=="__main__":main()
