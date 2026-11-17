#!/usr/bin/env python3
"""Replay exact minor certificates for the complete even graded Long cap23 jet."""
from fractions import Fraction as Q
from math import factorial,isqrt
from pathlib import Path
import argparse,hashlib,importlib.util,itertools,json,platform,tempfile

HERE=Path(__file__).resolve().parent
LABELS=("R","T","D","S")
WEIGHTS=(-1,2,4,1)
PAIRS=list(itertools.combinations(range(4),2))
CAP=23
def load(path,name):
 spec=importlib.util.spec_from_file_location(name,path)
 mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod

def basis(cap,k):
 """ALL source monomials of degree<=cap and weight w_i-6k."""
 out=[]
 for i,w in enumerate(WEIGHTS):
  for b in range(cap+1):
   for c in range(cap+1-b):
    for d in range(cap+1-b-c):
     a=b+4*c+2*d+6*k-w
     if a>=0 and a+b+c+d<=cap:out.append((i,(a,b,c,d)))
 return sorted(out)

def linear_columns(alg,images,unknowns):
 columns=[]
 for component,e in unknowns:
  monomial={e:Q(1)};column={}
  for pi,(i,j) in enumerate(PAIRS):
   if component==i:value=alg.bidifferential(monomial,images[j],1)
   elif component==j:value=alg.bidifferential(images[i],monomial,1)
   else:continue
   for exponent,c in value.items():column[(pi,exponent)]=c
  columns.append(column)
 return columns

def lower_rhs(alg,jets,k):
 """Negative known coefficient at h^(2k+1); unknown h^(2k) is excluded."""
 rhs={}
 for pi,(i,j) in enumerate(PAIRS):
  poly={}
  for a in range(min(k,len(jets)-1)+1):
   for b in range(min(k-a,len(jets)-1)+1):
    ell=k-a-b
    term=alg.scale(alg.bidifferential(jets[a][i],jets[b][j],2*ell+1),
                   Q(1,2**(2*ell)*factorial(2*ell+1)))
    poly=alg.add(poly,term)
  for e,c in poly.items():rhs[(pi,e)]=-c
 return rhs

def power(alg,f,n):
 out={(0,0,0,0):Q(1)}
 for _ in range(n):out=alg.multiply(out,f)
 return out

def vector(unknowns,polys):
 lookup={key:i for i,key in enumerate(unknowns)}
 out=[Q(0)]*len(unknowns)
 for i,poly in enumerate(polys):
  for e,c in poly.items():
   assert (i,e) in lookup,("support outside complete family",i,e)
   out[lookup[(i,e)]]=c
 return out

def product(columns,v):
 out={}
 for column,c in zip(columns,v):
  if c:
   for e,a in column.items():out[e]=out.get(e,Q(0))+c*a
 return {e:c for e,c in out.items() if c}

def determinant_mod(matrix,p):
 """Bare integer modular elimination; no FLINT or modular-rank library."""
 n=len(matrix);assert all(len(row)==n for row in matrix)
 a=[list(row) for row in matrix];det=1
 for j in range(n):
  pivot=next((i for i in range(j,n) if a[i][j]),None)
  if pivot is None:return 0
  if pivot!=j:a[j],a[pivot]=a[pivot],a[j];det=-det
  q=a[j][j];det=det*q%p;inverse=pow(q,-1,p)
  for i in range(j+1,n):
   if a[i][j]:
    factor=a[i][j]*inverse%p;a[i][j]=0
    for k in range(j+1,n):a[i][k]=(a[i][k]-factor*a[j][k])%p
 return det%p

def minor(columns,rhs,rows,certificate,p,augmented=False):
 def mod(c):
  assert c.denominator%p
  return c.numerator%p*pow(c.denominator%p,-1,p)%p
 n=len(columns)
 assert certificate["shape"]==[len(rows),n+int(augmented)]
 rr,cc=certificate["rows"],certificate["columns"]
 assert len(rr)==len(set(rr)) and len(cc)==len(set(cc))
 assert all(0<=i<len(rows) for i in rr)
 assert all(0<=j<n+int(augmented) for j in cc)
 matrix=[[mod(rhs.get(rows[i],Q(0)) if j==n
             else columns[j].get(rows[i],Q(0))) for j in cc] for i in rr]
 det=determinant_mod(matrix,p)
 assert det==certificate["det_modp"] and det
 return det

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument("--output-dir",type=Path);args=parser.parse_args()
 witness_path=HERE/"graded_long_cap23_witness.json"
 witness=json.loads(witness_path.read_text())
 p=witness["prime"]
 assert p>2 and all(p%d for d in range(2,isqrt(p)+1))
 assert p==1000003 and witness["cap"]==CAP
 alg=load(HERE/"check_long_graph_curvature.py","independent_graph_arithmetic")
 source=HERE.parent/"dixmier-audit/verify_long_pbw.py"
 producer=load(source,"canonical_long_input")
 images=[{tuple(e):Q(c) for e,c in f.items()} for f in producer.long_images()]
 for f,w in zip(images,WEIGHTS):
  assert all(-e[0]+e[1]+4*e[2]+2*e[3]==w for e in f)
 unknowns3=basis(CAP,1)
 assert [sum(i==j for i,e in unknowns3) for j in range(4)]==[54,80,102,71]
 columns3=linear_columns(alg,images,unknowns3)
 rhs3=lower_rhs(alg,[images],1)
 rows3=sorted(set(rhs3).union(*(set(c) for c in columns3)))
 particular=[{} for _ in range(4)]
 for name,e,c in witness["first_particular"]:
  i=LABELS.index(name);e=tuple(e);assert e not in particular[i]
  particular[i][e]=Q(c)
 assert product(columns3,vector(unknowns3,particular))==rhs3
 r,s=images[0],images[3]
 kernels=[
  [{},{},alg.scale(power(alg,r,2),-3),{}],
  [{},power(alg,r,4),alg.scale(alg.multiply(power(alg,r,3),s),-4),{}]]
 for kernel in kernels:assert not product(columns3,vector(unknowns3,kernel))
 assert not kernels[0][1] and kernels[1][1] and kernels[0][2]
 assert not particular[0]
 expected_s2={(5,0,0,0):Q(49,2),(6,1,0,0):Q(-165,4)}
 assert particular[3]==expected_s2
 det3=minor(columns3,rhs3,rows3,witness["h3_minor"],p)
 assert len(witness["h3_minor"]["rows"])==305
 # 307 unknowns, 305-minor and two exact independent kernels give EXACT
 # rational rank305 and exhaust the first-jet family over every char0 field.
 gauge_d4=alg.add(alg.multiply(power(alg,r,3),particular[3]),
                    alg.scale(alg.bidifferential(power(alg,r,3),s,2),Q(1,8)))
 assert max(map(sum,gauge_d4),default=0)<=16
 unknowns5=basis(CAP,2)
 assert [sum(i==j for i,e in unknowns5) for j in range(4)]==[20,34,47,29]
 columns5=linear_columns(alg,images,unknowns5)
 rhs5=lower_rhs(alg,[images,particular],2)
 rows5=sorted(set(rhs5).union(*(set(c) for c in columns5)))
 det5=minor(columns5,rhs5,rows5,witness["h5_augmented_minor"],p,True)
 assert len(unknowns5)==130 and len(witness["h5_augmented_minor"]["rows"])==131
 # Orientation of PAIRS has TS. Its negative RHS equals the actual ST defect.
 ts=PAIRS.index((1,3));st_defect={e:c for (pi,e),c in rhs5.items() if pi==ts}
 expected_st=alg.scale(alg.multiply({(10,0,0,0):Q(1)},power(alg,r,2)),Q(-49,4))
 assert st_defect==expected_st
 def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
 receipt={
  "claim":"No even equivariant full Long Weyl jet of common spatial degree23 through h5",
  "status":"PASS","environment":{"python":platform.python_version(),
   "arithmetic":"Fraction exact polynomial identities and bare integer modular determinants"},
  "cap":CAP,"source_weights":[-1,1,4,2],"h_weight":3,
  "image_order":list(LABELS),"image_weights":list(WEIGHTS),
  "h3":{"basis_counts":[54,80,102,71],"unknowns":307,"rows":len(rows3),
   "nonzeros":sum(map(len,columns3)),"rhs_nonzeros":len(rhs3),
   "particular_terms":sum(map(len,particular)),"exact_substitution":True,
   "exact_independent_kernels":["G=R^3","G=R^4*S"],
   "minor_size":305,"minor_det_modp":det3,"rank_Q":305,
   "all_solutions_first_R_correction_zero":True,
   "forced_S2":"49/2*X^5-165/4*X^6*Q"},
  "target_gauge_h4_D_term":"R^3*S2+Pi^2(R^3,S)/8",
  "target_gauge_h4_D_term_degree":max(map(sum,gauge_d4),default=0),
  "h5":{"basis_counts":[20,34,47,29],"unknowns":130,"rows":len(rows5),
   "rhs_nonzeros":len(rhs5),"augmented_minor_size":131,
   "augmented_minor_det_modp":det5,
   "rational_inconsistency":"augmented rank>=131 while coefficient rank<=130",
   "h2_truncation_ST_defect":"-49/4*X^10*R^2"},
  "prime":p,"hashes":{"check_graded_long_cap23.py":sha(Path(__file__)),
   "graded_long_cap23_witness.json":sha(witness_path),
   "GRADED_LONG_CAP23_JET_OBSTRUCTION.md":sha(HERE/"GRADED_LONG_CAP23_JET_OBSTRUCTION.md"),
   "check_long_graph_curvature.py":sha(HERE/"check_long_graph_curvature.py"),
   "../dixmier-audit/verify_long_pbw.py":sha(source)},
  "coverage_boundary":"Complete finite EVEN EQUIVARIANT degreecap23 first/second jet "
   "family, exact first-jet kernel classification, and nonzero augmented minor. "
   "No nongraded or odd-jet obstruction, cap24 exclusion, bounded normalization "
   "or DC1/DC2 theorem is claimed.",
  "independent_replay":False,"formal_verification":False,"external_review":False}
 output=args.output_dir or Path(tempfile.mkdtemp(prefix="graded-long-cap23-"))
 output.mkdir(parents=True,exist_ok=True)
 path=output/"graded_long_cap23_controls.json"
 path.write_text(json.dumps(receipt,indent=2)+"\n")
 print("PASS: h3 exact rank305/307; h5 augmented minor131 with130 unknowns")
 print("Modular determinants:",det3,det5,"prime",p)
 print("Receipt:",path)

if __name__=="__main__":main()
