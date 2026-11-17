#!/usr/bin/env python3
"""Complete cap24 odd first-jet field classification through commutator h3.

Uses Python integers/Fraction only; no algebraic or nullspace solver.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations_with_replacement
import json,hashlib,importlib.util,argparse,tempfile,platform
from math import isqrt
HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output-dir",type=Path)
args=parser.parse_args()
def load(path,name):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
cert=load(HERE/"check_graded_long_cap23.py","columns")
alg=load(HERE/"check_long_graph_curvature.py","algebra")
def determinant_Q(matrix):
 a=[list(row) for row in matrix];n=len(a);assert all(len(row)==n for row in a);result=Q(1)
 for j in range(n):
  pivot=next((i for i in range(j,n) if a[i][j]),None)
  if pivot is None:return Q(0)
  if pivot!=j:a[pivot],a[j]=a[j],a[pivot];result=-result
  q=a[j][j];result*=q
  for i in range(j+1,n):
   factor=a[i][j]/q
   for k in range(j+1,n):a[i][k]-=factor*a[j][k]
 return result
producer=load(HERE.parent/"dixmier-audit/verify_long_pbw.py","input")
F=[{tuple(e):Q(c) for e,c in p.items()} for p in producer.long_images()]
data_path=HERE/"graded_long_cap24_complex_witness.json"
data=json.loads(data_path.read_text())
assert data["cap"]==24 and data["source_weights"]==[-1,1,4,2]
assert data["image_order"]==list(cert.LABELS) and data["image_weights"]==list(cert.WEIGHTS)
assert data["h_weight"]==3
assert all(1000003%d for d in range(2,isqrt(1000003)+1))
add,mul,scale=alg.add,alg.multiply,alg.scale
power=lambda p,n:cert.power(alg,p,n)
r,T,D,S=F
variables=[{tuple(int(i==j) for i in range(4)):Q(1)} for j in range(4)]
for i,j in cert.PAIRS:
 assert alg.bidifferential(F[i],F[j],1)==alg.bidifferential(variables[i],variables[j],1)
for poly,weight in zip(F,cert.WEIGHTS):
 assert max(map(sum,poly))<=24
 assert all(-e[0]+e[1]+4*e[2]+2*e[3]==weight for e in poly)
X,Qpos,P,Z=variables;one={(0,0,0,0):Q(1)}
t=mul(X,Qpos)
B=add(scale(mul(power(X,2),P),3),scale(mul(add(one,scale(t,-3)),Z),2))
A0=add(Qpos,scale(mul(X,power(Qpos,2)),3),scale(mul(X,B),Q(-1,3)))
Y=add(mul(r,A0),scale(t,-3))
H4=add(power(Y,4),scale(mul(t,power(Y,2)),4),scale(mul(t,Y),Q(-16,3)),scale(power(t,2),10))
H6=add(mul(power(r,4),D),scale(add(scale(power(Y,5),9),scale(mul(t,power(Y,3)),45),
 scale(mul(power(t,2),Y),45),scale(mul(t,power(Y,2)),-60),scale(mul(t,Y),110),scale(power(t,2),-720)),Q(1,25)))
E=mul(power(X,2),Z)
H5=add(scale(power(Y,5),9),scale(mul(t,power(Y,3)),45),scale(mul(t,power(Y,2)),-60),
 scale(mul(power(t,2),Y),135),scale(mul(t,Y),-10),scale(power(t,2),-555))
def basis(m):
 out=[]
 for i,w in enumerate(cert.WEIGHTS):
  for b in range(25):
   for c in range(25-b):
    for d in range(25-b-c):
     a=b+4*c+2*d+3*m-w
     if a>=0 and a+b+c+d<=24:out.append((i,(a,b,c,d)))
 return sorted(out)
u1=basis(1);c1=cert.linear_columns(alg,F,u1)
rows1=sorted(set().union(*(set(c) for c in c1)))
assert len(u1)==487 and len(rows1)==2666
det1=cert.minor(c1,{},rows1,data["h1_minor"],1000003)
allgens=[t,Y,power(Y,2),mul(power(r,2),T),H4,H6,E,H5]
alljets=[[alg.bidifferential(G,f,1) for f in F] for G in allgens]
vec1=[cert.vector(u1,v) for v in alljets]
for v in vec1:assert not cert.product(c1,v)
ii=data["h1_independence"]["unknown_indices"]
assert len(ii)==len(set(ii))==8 and all(0<=i<len(u1) for i in ii)
inddet=determinant_Q([[v[i] for v in vec1] for i in ii])
assert inddet==Q(data["h1_independence"]["det_Q"])==3440640
assert len(data["h1_minor"]["rows"])==479 and len(vec1)==8
assert alljets[6][0]==scale(power(X,4),-3)
# Complete rational h3 compatibility inputs in normalizedsixparameterorder.
gens=[Y,mul(power(r,2),T),H4,H6,E,H5]
jets=[[alg.bidifferential(G,f,1) for f in F] for G in gens]
pairs=list(combinations_with_replacement(range(6),2))
assert data["quadratic_pairs"]==[list(p) for p in pairs]
rhs=[cert.lower_rhs(alg,[F],1)]
for a,b in pairs:
 value={}
 for pi,(i,j) in enumerate(cert.PAIRS):
  p=alg.bidifferential(jets[a][i],jets[b][j],1)
  if a!=b:p=add(p,alg.bidifferential(jets[b][i],jets[a][j],1))
  for e,q in p.items():value[(pi,e)]=-q
 rhs.append(value)
u2=basis(2);c2=cert.linear_columns(alg,F,u2)
rows2=sorted(set().union(*(set(c) for c in c2),*(set(b) for b in rhs)))
assert len(u2)==347 and len(rows2)==2137
det2=cert.minor(c2,{},rows2,data["h2_minor"],1000003)
index={pair:1+j for j,pair in enumerate(pairs)}
targets=[]
for spec in [{(3,3):Q(1),(3,5):Q(-250),(5,5):Q(35625)},{(4,4):Q(1)},{(4,5):Q(1)}]:
 v=[Q(0)]*len(rhs)
 for pair,q in spec.items():v[index[pair]]=q
 targets.append(v)
supports=[]
for item,target in zip(data["h3_functionals"],targets):
 functional={}
 for row,key,q in item["coefficients"]:
  key=(key[0],tuple(key[1]))
  assert rows2[row]==key and key not in functional
  functional[key]=Q(q)
 def apply(poly):return sum((q*poly.get(e,Q(0)) for e,q in functional.items()),Q(0))
 assert all(not apply(col) for col in c2)
 assert [apply(b) for b in rhs]==target
 supports.append(len(functional))
# Direct E=Q(q),q²=-2 bilinear calculation for four-parameter branch.
# Real/q components are expanded independently; no algebraic solver is used.
branchgens=[(Y,{}),(gens[1],{}),(H4,{}),(add(H5,scale(H6,125)),scale(H6,100))]
branchjets=[([alg.bidifferential(g[0],f,1) for f in F],
             [alg.bidifferential(g[1],f,1) for f in F]) for g in branchgens]
branchrhs=[(rhs[0],{})]
branchpairs=list(combinations_with_replacement(range(4),2))
for a,b in branchpairs:
 real={};imag={}
 for pi,(i,j) in enumerate(cert.PAIRS):
  ar,ai=branchjets[a];br,bi=branchjets[b]
  pr=add(alg.bidifferential(ar[i],br[j],1),scale(alg.bidifferential(ai[i],bi[j],1),-2))
  pq=add(alg.bidifferential(ar[i],bi[j],1),alg.bidifferential(ai[i],br[j],1))
  if a!=b:
   pr=add(pr,alg.bidifferential(br[i],ar[j],1),scale(alg.bidifferential(bi[i],ai[j],1),-2))
   pq=add(pq,alg.bidifferential(br[i],ai[j],1),alg.bidifferential(bi[i],ar[j],1))
  for e,q in pr.items():real[(pi,e)]=-q
  for e,q in pq.items():imag[(pi,e)]=-q
 branchrhs.append((real,imag))
# Retained block order places the seven old QQ blocks before the four newE blocks.
blockorder=[0,1,2,3,5,6,8,4,7,9,10]
def tuplepair(coefficients):
 out=([{} for _ in range(4)],[{} for _ in range(4)])
 for label,e,pair in coefficients:
  i=cert.LABELS.index(label);e=tuple(e)
  for j,q in enumerate(pair):
   q=Q(q)
   if q:out[j][i][e]=q
 return out
particulars=[]
for item,k in zip(data["particulars"],blockorder):
 candidate=tuplepair(item["coefficients"])
 for j in range(2):
  v=cert.vector(u2,candidate[j])
  assert cert.product(c2,v)==branchrhs[k][j]
 particulars.append(candidate)
assert len(particulars)==11
assert branchjets[-1][0][0]==scale(power(r,4),125)
assert branchjets[-1][1][0]==scale(power(r,4),100)
r7=power(r,7)
assert particulars[-1][0][0]==scale(r7,Q(-4800,7))
assert particulars[-1][1][0]==scale(r7,Q(-9600,7))
assert all(not candidate[j][0] for candidate in particulars[:-1] for j in range(2))
kernelgens=[power(r,3),mul(power(r,4),S),power(X,3)]
kerneljets=[]
for item,G in zip(data["homogeneous_kernels"],kernelgens):
 candidate=tuplepair(item["coefficients"])
 expected=[alg.bidifferential(G,f,1) for f in F]
 assert candidate[0]==expected and not any(candidate[1])
 assert not cert.product(c2,cert.vector(u2,expected))
 kerneljets.append(expected)
assert kerneljets[0][2] and not kerneljets[0][1] and not kerneljets[0][3]
assert kerneljets[1][1] and not kerneljets[1][3]
assert kerneljets[2][3]
assert len(data["h2_minor"]["rows"])==344 and len(u2)==347
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
out={"claim":"Complete cap24 odd first kernel and h3 field classification, retaining both complex conjugates",
 "status":"PASS","environment":{"python":platform.python_version(),"arithmetic":"Python integers/Fraction; rational coefficient pairs for q^2=-2"},
 "cap":24,"h1_basis_counts":[91,127,155,114],"h1_rank_Q":479,"h1_unknowns":487,"h1_rows":2666,"kernel_dimension":8,
 "h1_minor_size":479,"h1_minor_det_modp":det1,"h1_independence_size":8,"h1_independence_det_Q":str(inddet),
 "h2_basis_counts":[62,91,114,80],"h2_rank_Q":344,"h2_unknowns":347,"h3_rows":2137,
 "h2_minor_size":344,"h2_minor_det_modp":det2,"prime":1000003,
 "functional_supports":supports,"necessary_polynomials":["d^2-250*d*f+35625*f^2","e^2","e*f"],
 "particular_blocks_E":11,"particular_block_monomials":[p["monomial"] for p in data["particulars"]],
 "coefficient_pair_entries":sum(len(p["coefficients"]) for p in data["particulars"]),
 "homogeneous_kernels":["R^3","R^4*S","X^3"],
 "field":"E=Q(q),q^2=-2; both signs via q conjugation; z=0 included",
 "branch_hamiltonian":data["branch_hamiltonian"],
 "first_R":data["first_R"],"second_R":data["second_R"],
 "field_scope":"Characteristic-zero fields only. If -2 is not a square, d=f=e=0; otherwise both conjugate four-parameter branches are retained. No nonreduced-base classification.",
 "coverage_boundary":"Complete first kernel and commutator-h3 compatibility, with full F2 families. No h4 compatibility, full bounded lift, evenization or DC1/DC2 claim.",
 "hashes":{
  "GRADED_LONG_CAP24_COMPLEX_FIRST_JET_FAMILY.md":sha(HERE/"GRADED_LONG_CAP24_COMPLEX_FIRST_JET_FAMILY.md"),
  "check_graded_long_cap24_complex_family.py":sha(Path(__file__)),
  "graded_long_cap24_complex_witness.json":sha(data_path),
  "check_graded_long_cap23.py":sha(HERE/"check_graded_long_cap23.py"),
  "check_long_graph_curvature.py":sha(HERE/"check_long_graph_curvature.py"),
  "../dixmier-audit/verify_long_pbw.py":sha(HERE.parent/"dixmier-audit/verify_long_pbw.py")},
 "independent_replay":False,"formal_verification":False,"external_review":False}
output=args.output_dir or Path(tempfile.mkdtemp(prefix="graded-long-cap24-complex-"))
output.mkdir(parents=True,exist_ok=True)
receiptpath=output/"graded_long_cap24_complex_controls.json"
receiptpath.write_text(json.dumps(out,indent=2)+"\n")
print("PASS: cap24 odd first rank479/487; complete h3 field family, 11 E-points and 3 kernels",flush=True)
print("Minor determinants:",det1,det2,"mod1000003; independence:",inddet,flush=True)
print("Receipt:",receiptpath)
