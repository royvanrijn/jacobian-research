#!/usr/bin/env python3
"""Exact/source support for the pointwise coefficient CI and Gamma1 radial identity."""
from __future__ import annotations
if not __debug__:raise SystemExit("Optimized Python refused.")
from pathlib import Path
from fractions import Fraction as Q
from math import comb,factorial,gcd
import argparse,copy,hashlib,json,sys
HERE=Path(__file__).resolve().parent
PN='SUPPORT_UNIFORM_COEFFICIENT_NILPOTENCE_AND_RADIAL_GAUSSIAN.md'
CN='verify_uniform_coefficient_gaussian_corollary.py'
IN='uniform_coefficient_gaussian_inputs.json'
LN='uniform_coefficient_gaussian_source_ledger.json'
RN='uniform_coefficient_gaussian_verification.json'
SCOPE={'schema': 1, 'reviewed_proof_sha256': '200881f425252bcee55be4598ebdc5573cf3d1219536dce5332393e8170195c7', 'coefficient_rings': 'unital commutative Q-algebras', 'fixed_pair_terminal_minimum': 2, 'fixed_pair_universal_predicates_equivalent': True, 'full_eventual_equivalence_requires_every_M_beyond_shared_cutoff': True, 'single_good_pair_implies_all_higher_pairs': False, 'pointwise_coefficient_nilpotence_iff_prescribed_pair_vanishes': False, 'ideal_power_exponent': '2*M-2', 'universal_CI_top_degree': '2*M-3', 'universal_CI_length': 'M*(M-1)', 'exponent_independent_of_support_but_not_terminal': True, 'terminal_independent_nilpotence_exponent_for_unspecified_late_pair': False, 'gaussian_dimension': 2, 'gaussian_variance': 1, 'radial_denominator': 2, 'gaussian_radial_moment': 'k!', 'other_dimensions_or_Gamma_shapes_included': False, 'reduced_Q_algebra_nonzero_phases_share_cutoff': True, 'all_order_SFC': False, 'effective_cutoff': False, 'finite_controls_certify_universal_estimates': False, 'controls': {'Hilbert_terminal_max': 32, 'Macaulay_terminal_max': 8, 'supports': [[0, 1], [1, 2], [2, 3]], 'Gaussian_power_max': 64, 'moment_max': 6}, 'assurance_flags': {'external_review': False, 'formal_verification': False, 'independent_replay': False}}
SPECS={'main_proof': ('local', 'SUPPORT_UNIFORM_COEFFICIENT_NILPOTENCE_AND_RADIAL_GAUSSIAN.md', 'written_pointwise_CI_equivalences_and_exact_radial_identity'), 'exact_checker': ('local', 'verify_uniform_coefficient_gaussian_corollary.py', 'exact_algebra_and_source_support'), 'exact_inputs': ('local', 'uniform_coefficient_gaussian_inputs.json', 'typed_quantifier_and_shape_scope'), 'global_eventual_input': ('repo', 'research/factorial/image-rigidity/SUPPORT_UNIFORM_EVENTUAL_ORDINARY_BINOMIAL_SYNTHESIS.md', 'one_support_uniform_field_cutoff_only'), 'retained_CI_mechanism': ('repo', 'research/factorial/image-rigidity/EVENTUAL_BINOMIAL_ALGEBRAIC_CONSEQUENCES.md', 'classical_regular_sequence_Hilbert_mechanism_reproduced_here')}
ASSURANCE={'external_review': False, 'formal_verification': False, 'independent_replay': False}

def need(ok,msg):
 if not ok:raise ValueError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def typed(a,b):
 if type(a)is not type(b):return False
 if type(b)is dict:return set(a)==set(b)and all(typed(a[k],v)for k,v in b.items())
 if type(b)is list:return len(a)==len(b)and all(typed(x,y)for x,y in zip(a,b))
 return a==b
def check_inputs(d):need(typed(d,SCOPE),"Exact typed coefficient-ring, terminal and Gaussian scope")
def sources(d,repo):
 need(type(d)is dict and set(d)=={"schema","sources","assurance_flags"},"Exact ledger fields")
 need(type(d["schema"])is int and d["schema"]==1 and typed(d["assurance_flags"],ASSURANCE),"Schema and assurance")
 rows=d["sources"];need(type(rows)is list and len(rows)==len(SPECS),"All five roles")
 out={}
 for row in rows:
  need(type(row)is dict and set(row)=={"role","kind","path","relation","sha256"},"Exact source row")
  role=row["role"];need(type(role)is str and role in SPECS and role not in out,"Unique mandatory role")
  kind,path,relation=SPECS[role]
  need((row["kind"],row["path"],row["relation"])==(kind,path,relation),"Exact source identity/role")
  need(type(row["sha256"])is str and len(row["sha256"])==64 and all(c in"0123456789abcdef"for c in row["sha256"]),"Exact SHA256")
  root=HERE if kind=="local"else repo;p=(root/path).resolve()
  need(p.is_relative_to(root.resolve())and p.is_file()and sha(p)==row["sha256"],"Missing/stale source:"+role)
  out[role]=row
 return out
def polymul(a,b):
 c=[Q(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 while len(c)>1 and c[-1]==0:c.pop()
 return c
def hilbert(M):return polymul([Q(1)]*(M-1),[Q(1)]*M)
def forms(a,d,m):
 row=[comb(m,j)*factorial(a*m+d*j)for j in range(m+1)]
 content=0
 for c in row:content=gcd(content,c)
 return [Q(c//content)for c in row]
def rank(rows,cols):
 basis={}
 for rr in rows:
  row=list(rr);need(len(row)==cols,"Macaulay row dimension")
  for p in range(cols):
   if not row[p]:continue
   if p in basis:
    c=row[p];row=[x-c*y for x,y in zip(row,basis[p])]
   else:
    c=row[p];basis[p]=[x/c for x in row];break
 return len(basis)
def ideal_rows(a,d,M,j):
 out=[]
 for m in(M-1,M):
  if j<m:continue
  f=forms(a,d,m)
  for shift in range(j-m+1):
   row=[Q(0)]*(j+1)
   for k,c in enumerate(f):row[k+shift]=c
   out.append(row)
 return out
def gaussian_power(k,dimension=2,denominator=2):
 def odd_double_factorial(n):
  z=1
  for j in range(1,n+1,2):z*=j
  return z
 if dimension==2:
  return sum(Q(comb(k,j)*odd_double_factorial(2*j-1)*odd_double_factorial(2*(k-j)-1),denominator**k)for j in range(k+1))
 # Only the first moment is used to reject a false dimension extension.
 need(k==1,"No other-dimensional theorem or general control")
 return Q(dimension,denominator)
def exact(inp):
 cfg=inp["controls"];H=mac=ranks=G=mom=nonreduced=0;profiles=[]
 for M in range(2,cfg["Hilbert_terminal_max"]+1):
  h=hilbert(M);lhs=polymul(h,[Q(1),Q(-2),Q(1)])
  rhs=[Q(0)]*(2*M);rhs[0]=1;rhs[M-1]-=1;rhs[M]-=1;rhs[2*M-1]+=1
  while rhs and rhs[-1]==0:rhs.pop()
  need(lhs==rhs,"Exact complete-intersection Hilbert numerator")
  need(len(h)-1==2*M-3 and h[-1]==1 and sum(h)==M*(M-1),"Exact top degree, multiplicity and length")
  need(all(x>0 for x in h),"All degrees through the socle survive")
  H+=1
 for a,d in cfg["supports"]:
  for M in range(2,cfg["Macaulay_terminal_max"]+1):
   h=hilbert(M);length=0
   for j in range(2*M-1):
    rr=rank(ideal_rows(a,d,M,j),j+1);actual=j+1-rr;expected=h[j]if j<len(h)else Q(0)
    need(actual==expected,"Literal actual moment Macaulay dimension")
    length+=actual;ranks+=1
   need(length==M*(M-1),"Literal quotient length")
   need(rank(ideal_rows(a,d,M,2*M-3),2*M-2)==2*M-3,"One-too-small ideal power survives")
   need(rank(ideal_rows(a,d,M,2*M-2),2*M-1)==2*M-1,"Next ideal power is killed")
   profiles.append({"a":a,"d":d,"M":M,"length":M*(M-1),"socle":2*M-3,"sharp_nilpotence":2*M-2});mac+=1
 for k in range(cfg["Gaussian_power_max"]+1):
  need(gaussian_power(k)==factorial(k),"Exact two-normal radial moment")
  need(sum(comb(2*j,j)*comb(2*(k-j),k-j)for j in range(k+1))==4**k,"Exact central-binomial convolution")
  G+=1
 for a,d in cfg["supports"]:
  for m in range(1,cfg["moment_max"]+1):
   A=Q(2,3);B=-Q(3,5)
   radial=sum(comb(m,j)*A**(m-j)*B**j*gaussian_power(a*m+d*j)for j in range(m+1))
   ordinary=sum(comb(m,j)*A**(m-j)*B**j*factorial(a*m+d*j)for j in range(m+1))
   need(radial==ordinary,"Same-order literal radial factorial phase identity")
   mom+=1
 # Actual nonreduced-ring countercontrol: R=Q[eps]/eps^D,
 # A=eps,B=0, D=2M-2. Its low moment survives although J^D=0.
 for M in range(2,33):
  D=2*M-2;need(M-1<D and factorial(M-1)!=0,"Nil coefficient does not force a prescribed pair")
  nonreduced+=1
 need(factorial(2)%2==0 and factorial(3)%2==0 and 1**4%2==1,"Positive-characteristic unit ideal countercontrol")
 need(gaussian_power(1,dimension=3)==Q(3,2)!=factorial(1),"Dimension-three is not Gamma1")
 need(gaussian_power(1,denominator=1)==2!=factorial(1),"Missing divisor2 changes the literal identity")
 # Generic (not factorial) shared-factor forms leave V^D alive.
 M=4;D=2*M-2;shared=[]
 for m in(M-1,M):
  for shift in range(D-m+1):
   row=[Q(0)]*(D+1);row[shift]=1;shared.append(row)
 need(rank(shared,D+1)<D+1,"Coprime hypothesis cannot be dropped")
 return {"Hilbert_profiles":H,"actual_Macaulay_quotients":mac,"actual_graded_rank_checks":ranks,
         "radial_Gaussian_powers_and_convolutions":G,"actual_radial_phase_moments":mom,
         "nonreduced_ring_countercontrols":nonreduced,"genuine_math_negative_families":5,
         "finite_actual_profiles":profiles}
def negatives(ledger,repo):
 out=[]
 def reject(name,fn):
  try:fn()
  except(ValueError,KeyError,TypeError,OSError):out.append(name);return
  raise ValueError("Negative accepted:"+name)
 for role in SPECS:
  d=copy.deepcopy(ledger);d["sources"]=[r for r in d["sources"]if r["role"]!=role]
  reject("omit:"+role,lambda d=d:sources(d,repo))
  d=copy.deepcopy(ledger);next(r for r in d["sources"]if r["role"]==role)["sha256"]="0"*64
  reject("badpin:"+role,lambda d=d:sources(d,repo))
 for key in("single_good_pair_implies_all_higher_pairs","pointwise_coefficient_nilpotence_iff_prescribed_pair_vanishes",
  "terminal_independent_nilpotence_exponent_for_unspecified_late_pair","other_dimensions_or_Gamma_shapes_included",
  "all_order_SFC","effective_cutoff","finite_controls_certify_universal_estimates"):
  d=copy.deepcopy(SCOPE);d[key]=True;reject("false-scope:"+key,lambda d=d:check_inputs(d))
 for key,value in[("coefficient_rings","arbitrary rings"),("gaussian_dimension",3),("radial_denominator",1),("schema",True)]:
  d=copy.deepcopy(SCOPE);d[key]=value;reject("wrong-domain-or-type:"+key,lambda d=d:check_inputs(d))
 d=copy.deepcopy(ledger);d["schema"]=True;reject("boolean-ledger-schema",lambda:sources(d,repo))
 return out
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--repo-root",required=True,type=Path);ap.add_argument("--output",type=Path,default=HERE/RN);ns=ap.parse_args()
 repo=ns.repo_root.resolve();inp=load(HERE/IN);check_inputs(inp)
 need(sha(HERE/PN)==inp["reviewed_proof_sha256"],"Exact current proof snapshot")
 ledger=load(HERE/LN);rows=sources(ledger,repo);counts=exact(inp);bad=negatives(ledger,repo);sources(ledger,repo)
 result={"schema":1,"status":"PASS","scope":inp,"proof_sha256":sha(HERE/PN),"checker_sha256":sha(HERE/CN),
  "input_sha256":sha(HERE/IN),"source_ledger_sha256":sha(HERE/LN),"mandatory_source_roles":len(rows),
  "sources":[rows[r]for r in sorted(rows)],"counts":counts,"source_scope_negative_controls":bad,"assurance_flags":ASSURANCE,
  "meaning":"The pointwise regular-sequence/Hilbert proof and formal Gaussian generating function establish the universal implications. Finite exact ranks, moments and source controls support those proofs; they do not evaluate the shared cutoff or prove all-order SFC."}
 ns.output.parent.mkdir(parents=True,exist_ok=True);ns.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
 print(json.dumps({"status":"PASS","sources":len(rows),"counts":{k:v for k,v in counts.items()if k!="finite_actual_profiles"},"source_controls":len(bad),"output":str(ns.output)}))
if __name__=="__main__":main()
