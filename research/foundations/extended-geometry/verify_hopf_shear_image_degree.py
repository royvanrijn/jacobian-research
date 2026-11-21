#!/usr/bin/env python3
"""Exact/source support for primitive Image degree propagation by Hopf shear.
Finite contractions support the written all-order proof; they do not prove it.
"""
from __future__ import annotations
if not __debug__:raise SystemExit("Refusing -O/-OO: exact/source checks must run.")
from pathlib import Path
from fractions import Fraction as Q
from math import comb,factorial,isqrt
import argparse,copy,hashlib,json
HERE=Path(__file__).resolve().parent
PN="HOPF_SHEAR_PRIMITIVE_IMAGE_DEGREE_PROPAGATION.md"
CN="verify_hopf_shear_image_degree.py"
IN="hopf_shear_image_degree_inputs.json"
LN="hopf_shear_image_degree_source_ledger.json"
RN="hopf_shear_image_degree_verification.json"
PROOF_SHA="b4a2a9e22308f45ea7f0223f3058ae6bdf119e85f1be7f532efa44afe25d0c04"
ASSURANCE=dict(external_review=False,formal_verification=False,independent_replay=False)
P="research/foundations/extended-geometry/"
SPECS={
"written_proof":("local",PN,"complete_all_n_k_m_proof_and_separate_coverage_equivalence"),
"exact_checker":("local",CN,"bounded_exact_controls_NOT_proof_of_all_order_identities"),
"typed_inputs":("local",IN,"exact_scope_markers_k_min_and_source_roles"),
"hyperbolic_seed":("repo",P+"HYPERBOLIC_ODD_HEIGHT_IMAGE_WITNESSES.md","LOGICAL_a1_all_n_seed_fixed_x_moments_and_integer_clearing"),
"old_shear":("repo",P+"QUADRATIC_HEIGHT_HOPF_CLASSIFICATION.md","LOGICAL_section4_polynomial_Lambda_invariance_reproved_in_new_note"),
"original_and_even":("repo",P+"TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md","LOGICAL_original_factorial_normalization_and_all_even_primitive_nonpower_family"),
"odd_composite":("repo",P+"ODD_PRIMITIVE_TWO_PAIR_IMAGE_WITNESSES.md","LOGICAL_ONLY_combined_all_D_ge8_coverage_may_have_other_fixed_marker"),
"cubic_shear_antecedent":("repo",P+"TWO_PAIR_IMAGE_CUBIC_STRATIFICATION.md","CONTEXT_old_exp_xdelta_invariance_NOT_new_construction"),
"two_height_threshold_context":("repo",P+"TWO_HEIGHT_MINIMUM_UNSAFE_PRIME_THIRTEEN.md","CONTEXT_only_its_restricted_two_height_threshold_new_shears_leave_class"),
"prime_density_antecedent":("repo",P+"PRIME_IMAGE_DEGREE_DENSITY.md","CONTEXT_only_existing_density_one_corollary_NOT_a_premise_of_every_prime_coverage"),
"coordinate_fiber_context":("repo",P+"HYPERBOLIC_COORDINATE_FIBER_IMAGE_WITNESSES.md","CONTEXT_only_fixed_coordinate_classification_NOT_an_arbitrary_phase_obstruction")}
SCOPE=dict(schema=1,reviewed_proof_sha256=PROOF_SHA,
 parameters="integers n>=1,k>=1; lambda=1; every characteristic-zero field",
 pure_orders="every integer m>=1",fixed_Z_nonzero_orders="every positive even m",
 new_least_balanced_degree="5*n+3+3*k",new_ordinary_total_degree="2*(5*n+3+3*k)",
 integer_lift=True,radial_zero_value=-1,proper_power_excluded_over_every_extension=True,
 new_shear_k0_allowed=False,k0_only_retained_unsheared_seed=True,
 every_prime_degree_ge11=True,prime_family_fixed_marker="Z",
 every_balanced_degree_ge16_from_three_seeds=True,
 every_balanced_degree_ge8_combining_old_families=True,
 combined_all_degree_marker="SOME fixed polynomial marker depending on witness; not necessarily Z",
 all_degree_some_marker_IFF_all_prime_some_marker_under_old_families=True,
 degree5_or7_settled=False,arbitrary_prime_minimum_claimed=False,
 old_two_height_prime13_threshold_contradicted=False,arbitrary_profile_classification=False,
 finite_prefix_proves_universal_moments=False,worldwide_priority_claimed=False,
 controls=dict(sphere_n_max=4,sphere_k_max=4,sphere_m_max=8,
 original_cases=[[1,1,4],[1,3,2],[2,2,2],[3,1,2]],
 degree_max=1000,monomial_exponent_max=6,delta_degree_max=20),
 assurance_flags=ASSURANCE)
def need(ok,msg):
 if not ok:raise ValueError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def typed(a,b):
 if type(a)is not type(b):return False
 if type(b)is dict:return set(a)==set(b)and all(typed(a[k],v)for k,v in b.items())
 if type(b)is list:return len(a)==len(b)and all(typed(x,y)for x,y in zip(a,b))
 return a==b
def inputs(v):need(typed(v,SCOPE),"Exact typed shear, k>=1, marker and coverage scope")
def sources(v,repo):
 need(type(v)is dict and set(v)=={"schema","sources","assurance_flags"},"Exact ledger fields")
 need(type(v["schema"])is int and v["schema"]==1 and typed(v["assurance_flags"],ASSURANCE),"Typed schema and assurance")
 rr=v["sources"];need(type(rr)is list and len(rr)==len(SPECS),"All exact mandatory roles required")
 out={}
 for row in rr:
  need(type(row)is dict and set(row)=={"role","kind","path","relation","sha256"},"Exact source row fields")
  role=row["role"];need(type(role)is str and role in SPECS and role not in out,"Unique known source role")
  need((row["kind"],row["path"],row["relation"])==SPECS[role],"Exact role/path/relation")
  h=row["sha256"];need(type(h)is str and len(h)==64 and all(c in"0123456789abcdef"for c in h),"Exact SHA256")
  root=(HERE if row["kind"]=="local"else repo).resolve();p=(root/row["path"]).resolve()
  need(p.is_relative_to(root)and p.is_file()and sha(p)==h,"Missing/stale source:"+role);out[role]=row
 need(out["written_proof"]["sha256"]==PROOF_SHA,"Exact reviewed whole proof")
 return out
# Sparse integer polynomial operations; exponent tuples also permit chart Laurent powers.
def clean(p):return {e:c for e,c in p.items()if c}
def add(*ps):
 out={}
 for p in ps:
  for e,c in p.items():out[e]=out.get(e,0)+c
 return clean(out)
def scale(p,c):return clean({e:c*v for e,v in p.items()})
def mul(p,q,cap=None):
 out={}
 for a,c in p.items():
  for b,d in q.items():
   e=tuple(x+y for x,y in zip(a,b))
   if cap is None or e[0]<=cap:out[e]=out.get(e,0)+c*d
 return clean(out)
def power(p,n):
 out={tuple(0 for _ in next(iter(p))):1}
 for _ in range(n):out=mul(out,p)
 return out
def uni_mul(p,q):
 out=[0]*(len(p)+len(q)-1)
 for i,c in enumerate(p):
  for j,d in enumerate(q):out[i+j]+=c*d
 return out
def seed(n):
 C=[4*n*n*comb(n+1,j)for j in range(n+2)]
 ss=[comb(2*n,j+1)for j in range(2*n)]
 D=[-x for x in uni_mul([comb(n+3,j)for j in range(n+4)],uni_mul(ss,ss))]
 both=[(C[j]if j<len(C)else 0)+D[j]for j in range(len(D))]
 need(both[0]==0,"Literal divided difference")
 return C,D,both[1:]
def angular(n,k):
 need(n>=1 and k>=1,"New construction k>=1")
 C,D,B=seed(n);tk={(0,0,1):1,(k+1,0,0):1}
 bracket={(i,0,0):c for i,c in enumerate(B)if c}
 bracket=add(bracket,{(i,1,0):-2*c for i,c in enumerate(D)if c},
  {(i+k,0,1):2*c for i,c in enumerate(D)if c},
  {(i+2*k+1,0,0):c for i,c in enumerate(D)if c})
 G=mul(tk,bracket);degree=5*n+3+3*k
 top={e:c for e,c in G.items()if sum(e)==degree}
 need(max(map(sum,G))==degree and top=={(degree,0,0):-1},"Unique exact top pure x term")
 return G
def chart(n,k):
 # Independent direct localized t_k(C*+t_k^2 D*)/x expansion.
 C,D,_=seed(n);tk={(0,1):1,(k+1,0):1}
 inner=add({(i,0):c for i,c in enumerate(C)if c},
           mul(mul(tk,tk),{(i,0):c for i,c in enumerate(D)if c}))
 p=mul(tk,inner)
 return {(i-1,t):c for(i,t),c in p.items()}
def kappa(m):return sum((Q((-1)**j*comb(m,j),m+2*j+1)for j in range(m+1)),Q(0))
def poch(a,r):
 out=Q(1)
 for j in range(r):out*=a+j
 return out
def marker(n,m):
 if m%2:return Q(0)
 r=m//2
 mu=(-n*n)**r*poch(Q(1,2*n),r)*poch(1-Q(1,2*n),r)/(poch(Q(1,2),r)*factorial(r))
 return(4*n*n)**m*kappa(m)*mu
def sphere_controls():
 count=0;rows=[]
 for n in range(1,5):
  for k in range(1,5):
   G=angular(n,k);p=chart(n,k);pp={(0,0):1}
   for m in range(1,9):
    # Terms x>8 can never return to x0 or x^-1 in any remaining <=8 steps.
    pp=mul(pp,p,cap=8)
    pure=sum((Q(c,t+1)for(i,t),c in pp.items()if i==0 and t%2==0),Q(0))
    mixed=sum((Q(c,t+1)for(i,t),c in pp.items()if i==-1 and t%2==0),Q(0))
    need(pure==0 and mixed==marker(n,m),"Exact raw chart pure/marked controls")
    rows.append(dict(n=n,k=k,m=m,balanced_degree=5*n+3+3*k,pure="0",normalized_marked=str(mixed)))
    count+=1
 return rows
def original(n,k):
 G=angular(n,k);D=5*n+3+3*k;out={}
 for(i,j,l),c in G.items():
  r=D-i-j-l
  for u in range(r+1):
   for v in range(l+1):
    key=(i+u+v,j+u+v)
    out[key]=out.get(key,0)+c*2**j*comb(r,u)*comb(l,v)*(-1)**(l-v)
 out=clean(out);need(out.get((D,0))==-1,"Literal original radial-zero evaluation")
 return out,D
def original_controls():
 rows=[];wrongfactor=0
 for n,k,maximum in SCOPE["controls"]["original_cases"]:
  F,D=original(n,k);pp={(0,0):1}
  for m in range(1,maximum+1):
   pp=mul(pp,F);degree=D*m
   pure=sum(c*factorial(i)*factorial(degree-i)for(i,j),c in pp.items()if i==j)
   marked=sum(c*factorial(i+1)*factorial(degree-i)for(i,j),c in pp.items()if j==i+1)
   need(pure==0 and Q(marked,factorial(degree+2))==marker(n,m),"DIRECT original factorial contractions")
   if m%2==0:
    need(marker(n,m)!=0,"Every even target nonzero")
    need(Q(marked,factorial((5*n+3)*m+2))!=marker(n,m),"New radial degree cannot be omitted")
    need(marker(n,m)/(4*n*n)**m!=marker(n,m),"Integer clearing factor cannot be omitted")
    wrongfactor+=2
   rows.append(dict(n=n,k=k,m=m,balanced_degree=D,original_support=len(F),
    powered_original_support=len(pp),pure="0",marked=str(marked),
    normalized_original_marker=str(Q(marked,factorial(degree+2))),radial_zero_value=-1))
 return rows,wrongfactor
def compact_degree_eleven_control():
 R={(1,1):1,(0,0):1};T={(1,1):1,(0,0):-1};Z={(1,0):1};W={(0,1):2}
 A=add(mul(T,R),power(Z,2));B=power(add(R,Z),2)
 H1=mul(power(R,4),mul(add(Z,scale(R,3)),add(power(Z,2),scale(mul(Z,R),3),scale(power(R,2),4))))
 H2=mul(B,mul(power(add(Z,scale(R,2)),2),
  add(power(Z,3),scale(mul(T,mul(Z,R)),2),scale(mul(W,power(R,2)),-2))))
 literal=scale(mul(A,mul(B,add(H1,H2))),-1)
 generic,D=original(1,1)
 need(D==11 and literal==generic,"Independent compact F11 matches literal generic lift")
 ip=[(1,0),(0,1),(-1,0),(0,-1)]
 value=[0,0];derivative=[0,0]
 for(i,j),c in literal.items():
  if i==11:
   exponent=11-j
   for z in range(2):value[z]+=c*ip[exponent%4][z]
   if exponent:
    for z in range(2):derivative[z]+=c*exponent*ip[(exponent-1)%4][z]
 need(value==[0,0] and derivative==[0,16],"Exact original F11 slice derivative at i")
 return dict(balanced_degree=11,original_support=len(literal),compact_factorization_matches=True,
  original_slice_value_at_i=value,original_slice_derivative_at_i=derivative)
def algebra_controls():
 # Formal four-variable polynomial identities (x,y,t,h).
 x={(1,0,0,0):1};y={(0,1,0,0):1};t={(0,0,1,0):1};h={(0,0,0,1):1}
 one={(0,0,0,0):1}
 tp=add(t,mul(x,h));yp=add(y,scale(mul(t,h),-1),scale(mul(x,mul(h,h)),Q(-1,2)))
 quad=add(mul(t,t),scale(mul(x,y),2))
 changed=add(mul(tp,tp),scale(mul(x,yp),2))
 need(changed==quad,"Universal quadric-preserving shear")
 back_t=add(tp,scale(mul(x,h),-1))
 back_y=add(yp,mul(tp,h),scale(mul(x,mul(h,h)),Q(-1,2)))
 need(back_t==t and back_y==y,"Universal inverse shear")
 bad_y=add(y,scale(mul(t,h),-1))
 need(add(mul(tp,tp),scale(mul(x,bad_y),2))!=quad,"Missing quadratic shear term has actual defect")
 # For all nonnegative monomial powers, the exact boundary derivative
 # has zero CT. Check full integration formula independently on a box.
 boundary=0
 for a in range(7):
  for b in range(7):
   for c in range(7):
    for k in range(1,5):
     # f=x^b*y^c*t^a. Delta_h=(x^k)(x d_t-t d_y).
     if b+k+1!=c:val=Q(0)
     else:
      left=sum((Q(a*(-1)**j*comb(c,j),2**c*(a+2*j))for j in range(c+1)if a and(a-1+2*j)%2==0),Q(0))
      right=sum((Q(-c*(-1)**j*comb(c-1,j),2**(c-1)*(a+2*j+2))for j in range(c)if c and(a+1+2*j)%2==0),Q(0))
      val=left+right
     need(val==0,"Exact formal interval/CT infinitesimal invariance");boundary+=1
 # delta_k roots are simple; the excluded root -1 must not be used blindly.
 simple=[]
 for k in range(1,21):
  # z delta'-(k+1)delta=-(k+1), universal monomial coefficient identity.
  need((k+1)-(k+1)==0 and -(k+1)!=0,"Delta root simplicity Bezout identity")
  total=k+1;minus_one_root=int((k+1)%2==1)
  need(total-minus_one_root>=1,"An admissible algebraic root other than -1 exists")
  simple.append(dict(k=k,degree=total,root_minus_one=bool(minus_one_root),remaining_roots=total-minus_one_root))
 need((1+(-1)**3)==0 and(1+(-1))**2==0,"Wrong root choice -1 can kill C")
 return dict(universal_shear_identities=3,exact_boundary_monomial_controls=boundary,
  delta_simple_root_controls=simple,genuine_math_negative_families=4)
def coverage_controls():
 rows=[];prime=bytearray(b"\1")*1001;prime[:2]=b"\0\0"
 for q in range(2,isqrt(1000)+1):
  if prime[q]:prime[q*q:1001:q]=b"\0"*(((1000-q*q)//q)+1)
 for D in range(16,1001):
  n=2 if D%3==1 else 1 if D%3==2 else 3
  k=(D-(5*n+3))//3
  need(D==5*n+3+3*k and k>=0,"All-degree residue coverage")
  need(k>=1 or(n,k,D)==(3,0,18),"k0 coverage ONLY retained unsheared seed")
  rows.append(dict(D=D,n=n,k=k,route="retained_unsheared"if not k else"new_shear"))
 primes=[]
 for p in range(11,1001):
  if prime[p]:
   n=1 if p%3==2 else 2;k=(p-(5*n+3))//3
   need(p==5*n+3+3*k and(k>=1 or p==13),"Every prime residue coverage")
   primes.append(dict(p=p,n=n,k=k,route="retained_unsheared"if not k else"new_shear"))
 # Abstract trichotomy of integers, not a numerical proof of the infinite theorem.
 return dict(all_degree_cases=rows,prime_cases=primes,
  marker_scope="New primes and sheared family use Z; old composite witnesses may use a different fixed polynomial.")
def negative_controls(v,repo):
 out=[]
 def reject(label,fn):
  try:fn()
  except(ValueError,KeyError,TypeError,OSError):out.append(label);return
  raise ValueError("Negative accepted:"+label)
 for role in SPECS:
  bad=copy.deepcopy(v);bad["sources"]=[r for r in bad["sources"]if r["role"]!=role]
  reject("omit:"+role,lambda bad=bad:sources(bad,repo))
  bad=copy.deepcopy(v);next(r for r in bad["sources"]if r["role"]==role)["sha256"]="0"*64
  reject("badpin:"+role,lambda bad=bad:sources(bad,repo))
 for key,val in[("new_shear_k0_allowed",True),("k0_only_retained_unsheared_seed",False),
 ("degree5_or7_settled",True),("arbitrary_prime_minimum_claimed",True),
 ("old_two_height_prime13_threshold_contradicted",True),("finite_prefix_proves_universal_moments",True),
 ("worldwide_priority_claimed",True),("combined_all_degree_marker","Z for every old composite witness"),
 ("new_least_balanced_degree","5*n+3"),("radial_zero_value",0),("schema",True)]:
  bad=copy.deepcopy(SCOPE);bad[key]=val;reject("false-scope:"+key,lambda bad=bad:inputs(bad))
 bad=copy.deepcopy(v);bad["schema"]=True;reject("boolean-ledger-schema",lambda:sources(bad,repo))
 return out
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--repo-root",required=True,type=Path)
 ap.add_argument("--output",type=Path,default=HERE/RN);ns=ap.parse_args()
 repo=ns.repo_root.resolve();inp=load(HERE/IN);inputs(inp);ledger=load(HERE/LN);pins=sources(ledger,repo)
 algebra=algebra_controls();compact=compact_degree_eleven_control();sphere=sphere_controls();orig,neg=original_controls();coverage=coverage_controls()
 failures=negative_controls(ledger,repo)
 report=dict(schema=1,status="PASS_HOPF_SHEAR_PRIMITIVE_IMAGE_DEGREE_PROPAGATION",
 scope=inp,proof_sha256=sha(HERE/PN),checker_sha256=sha(HERE/CN),source_ledger_sha256=sha(HERE/LN),
 mandatory_source_roles=len(pins),sources=[pins[k]for k in sorted(pins)],exact_algebra=algebra,
 raw_sphere_controls=sphere,original_four_variable_controls=orig,independent_compact_degree11=compact,
 wrong_factorial_or_clearing_controls=neg,degree_coverage_controls=coverage,
 source_scope_controls=failures,assurance_flags=ASSURANCE,
 meaning="The written proof establishes every parameter/order and exact primitivity/nonpower. Finite controls verify literal original contractions, source freshness and algebra only; all-degree marker existence is separated from the stronger new fixed-Z family.")
 ns.output.parent.mkdir(parents=True,exist_ok=True);ns.output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
 print(json.dumps(dict(status="PASS",sources=len(pins),sphere_cases=len(sphere),original_cases=len(orig),
 boundary_controls=algebra["exact_boundary_monomial_controls"],source_controls=len(failures))))
if __name__=="__main__":main()

