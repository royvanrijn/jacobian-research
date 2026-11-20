#!/usr/bin/env python3
"""Source-guarded finite-regime synthesis for one support-uniform eventual cutoff."""
from __future__ import annotations
if not __debug__:raise SystemExit("Refusing -O/-OO: exact/source controls must execute.")
from pathlib import Path
from fractions import Fraction as Q
from math import comb,factorial
import argparse,copy,hashlib,importlib.util,json,re,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
PN="SUPPORT_UNIFORM_EVENTUAL_ORDINARY_BINOMIAL_SYNTHESIS.md"
CN="verify_support_uniform_ordinary_completion.py"
IN="support_uniform_ordinary_completion_inputs.json"
LN="support_uniform_ordinary_completion_source_ledger.json"
RN="support_uniform_ordinary_completion_verification.json"
PREFIX="research/factorial/image-rigidity/"
CHILDREN={
 "quadratic":("UNIFORM_QUADRATIC_STRONG_FACTORIAL_TAIL.md","verify_uniform_quadratic_strong_factorial_completion.py","uniform_quadratic_strong_factorial_completion_inputs.json","uniform_quadratic_strong_factorial_completion_source_ledger.json","uniform_quadratic_strong_factorial_completion_verification.json","files",13,"PASS_SOURCE_GUARDED_UNIFORM_QUADRATIC_SUPPORT",False,"source_ledger_sha256"),
 "small_ratio":("UNRESTRICTED_SMALL_RATIO_GAMMA_TWO_WELL_EXCLUSION.md","verify_unrestricted_small_ratio.py","unrestricted_small_ratio_inputs.json","unrestricted_small_ratio_source_ledger.json","unrestricted_small_ratio_verification.json","sources",4,"PASS",False,"source_ledger_sha256"),
 "subcritical":("COMPACT_SUBCRITICAL_RATIO_UNIFORM_EVENTUAL_SFC.md","verify_compact_subcritical_ratio_uniform.py","compact_subcritical_ratio_inputs.json","compact_subcritical_ratio_source_ledger.json","compact_subcritical_ratio_verification.json","sources",13,"PASS",True,"source_ledger_sha256"),
 "bounded_thermal":("UNIFORM_BOUNDED_THERMAL_TRANSITION_EVENTUAL_SFC.md","verify_bounded_thermal_transition.py","bounded_thermal_transition_inputs.json","bounded_thermal_transition_source_ledger.json","bounded_thermal_transition_verification.json","sources",27,"PASS",True,"source_ledger_sha256"),
 "large_thermal":("LARGE_THERMAL_QUARTIC_FOLD_MATCHING.md","verify_large_thermal_quartic_fold_matching.py","large_thermal_quartic_fold_inputs.json","large_thermal_quartic_fold_source_ledger.json","large_thermal_quartic_fold_verification.json","sources",23,"PASS",False,"source_ledger_sha256"),
 "supercritical":("COMPACT_SUPERCRITICAL_RATIO_UNIFORM_EVENTUAL_SFC.md","verify_compact_supercritical_ratio_uniform.py","compact_supercritical_ratio_universal_ledger.json","compact_supercritical_ratio_source_manifest.json","compact_supercritical_ratio_verification.json","sources",17,"PASS",True,"source_manifest_sha256")}
ASSURANCE={"formal_verification":False,"external_review":False,"independent_replay":False}
SPECS={"main_proof":("local",PN,"written_finite_cover_and_six_way_scope_argument"),
 "exact_checker":("local",CN,"exact_algebra_logical_and_source_support"),
 "exact_inputs":("local",IN,"typed_scope_and_universal_cover_ledger")}
for tag,record in CHILDREN.items():
 for role,name in zip(["proof","checker","inputs","ledger","receipt"],record[:5]):
  SPECS[tag+"_"+role]=("local",name,"logical_regional_theorem_"+tag if role=="proof" else "fresh_scoped_regional_evidence_"+tag)
for role,name in [("large_thermal_context","LARGE_THERMAL_PRIMARY_CONTEXT.md"),("large_thermal_primary_framework","sources/dlmf-36-12-large-thermal.html")]:
 SPECS[role]=("local",name,"declared_large_thermal_primary_context_only")
for role,name,relation in [
 ("raw_adjacent","UNIFORM_PADE_PARABOLA_SFC.tex","LOGICAL_thm_main_ordinary_d1_ONLY"),
 ("raw_primepowers","PRIME_POWER_TWO_TERM_SFC.tex","LOGICAL_thm_sfc_alpha1_d2_3_4_5_ONLY"),
 ("raw_constant_analytic","ALL_UNIVARIATE_TWO_TERM_SFC.tex","LOGICAL_thm_analytic_Nge10000_6ledlesqrtN_ONLY_NO_finite_prefix"),
 ("raw_constant_large","ALL_EXPONENTS_AND_LARGE_EXPONENT_CUTOFF.tex","LOGICAL_thm_large_2powerdge4Nsquared_ONLY"),
 ("raw_six_entropy","EFFECTIVE_EXPONENT_SIX_SIMPLICITY.tex","AUTHENTICATED_analytic_credit_source_closure_NOT_separate_theorem_NO_finite_prefix"),
 ("raw_terminal_highgap","lower_exponent_initial_and_boundary.tex","LOGICAL_thm_gapbound_fixed_terminal_support_finiteness_ONLY"),
 ("raw_Gamma_soft","GAMMA_PRODUCT_ROOT_MODULUS.tex","LOGICAL_old_soft_fixed_terminal_support_finiteness_ONLY"),
 ("exact_polynomial_kernel","verify_strict_subcritical_eventual_sfc.py","CODE_ONLY_sparse_rational_polynomial_kernel")]:
 SPECS[role]=("repo",PREFIX+name,relation)
for role,name in zip(["proof","checker","inputs","ledger","receipt"],[
 "ALL_ORDINARY_BINOMIAL_EVENTUAL_STRONG_FACTORIAL.md","verify_ordinary_binomial_eventual_completion.py",
 "ordinary_binomial_eventual_inputs.json","ordinary_binomial_eventual_source_ledger.json","ordinary_binomial_eventual_verification.json"]):
 SPECS["fixed_support_407_"+role]=("local",name,"LOGICAL_407_ONLY_finite_support_reverse_and_G7_bounded_gap_fallback"if role=="proof"else"source_407_already_freshly_replayed_inside_G7")
SCOPE={"schema":1,"parameters":"integers a>=0,d>=1,N>=2; c=a/(d-1)^2 ONLY when a>=1,d>=2",
 "main_reviewed_sha256":"2923c4c29f5e29dd02fb5f25d50e5167e6b9d71b1a66bf41d8e24453f7a1a722",
 "theorem":"one unevaluated integer Nstar for ALL ordinary supports and every N>=Nstar",
 "all_nonzero_coefficients_and_characteristic_zero_fields":True,
 "rational_Bezout_no_arbitrary_field_embedding":True,
 "small_ratio_G3_requires_RAW_old_sparse":True,
 "G4_compact_subcritical_interval":"fixed real compact subset of (0,1/4)",
 "G5_used_at_ONE_fixed_finite_U":True,
 "G6_sequential_to_fixed_relative_strip":True,
 "G7_Kabluchko_2_3_logical_and_407_bounded_gap_only":True,
 "quadratic_G2_domain":"a>=10^24*d^2,d>=2,N>=2",
 "K":4*10**24,
 "critical_halfwidth":"h=min(1/8,delta/16)",
 "small_cut":"epsilon=min(epsilon_small,1/16)",
 "constant_family_uniform_terminal":10000,
 "a0_boundary_uses_RAW_analytic_large_and_primepower_proofs_ONLY":True,
 "a0_92935_or9998_finite_certificate_import":False,
 "global_bad_triple_set_finite":True,
 "finite_projective_exceptional_phases":True,
 "exceptional_projective_ratios_algebraic_over_Q":True,
 "bad_support_projection_finite":True,
 "six_way_equivalence_under_declared_inputs":True,
 "QQ_tail_to_finite_B_requires_raw_fixed_terminal_input_NOT407":True,
 "finite_support_or_projective_phases_to_finite_B_REQUIRES407":True,
 "unscaled_coefficient_pairs_finite":False,
 "all_order_SFC":False,"global_bad_set_empty":False,"effective_Nstar_or_containing_box":False,
 "replication_of_one_bad_triple":False,"arbitrary_polynomial_supports":False,"general_Gamma_shapes":False,
 "global_analytic_estimates_newly_proved_by_finite_controls":False,
 "all_fast_regional_guards_replayed_fresh":True,
 "heavy_finite_boundary_runs_repeated":False,
 "controls":{"maximum_lower":4,"maximum_gap":8,"maximum_order":5,
  "terminal_count_maximum":64,"root_t_samples":[100,101,128,150,200],
  "critical_denominators":[1,2,3,7,19,64],
  "critical_deltas":["1/1000000000000","1/1000","1/10","1","10"],
  "small_epsilons":["1/1000000","1/100","1","100"]},
 "assurance_flags":ASSURANCE}
def need(ok,msg):
 if not ok:raise ValueError(msg)
def H(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def typed(a,b):
 if type(a)is not type(b):return False
 if type(b)is dict:return set(a)==set(b)and all(typed(a[k],v)for k,v in b.items())
 if type(b)is list:return len(a)==len(b)and all(typed(x,y)for x,y in zip(a,b))
 return a==b
def check_inputs(d):need(typed(d,SCOPE),"Exact global finite-regime and six-way input scope")
def check_sources(d,repo):
 need(type(d)is dict and set(d)=={"schema","assurance_flags","sources"},"Exact source-ledger schema")
 need(type(d["schema"])is int and d["schema"]==1 and typed(d["assurance_flags"],ASSURANCE),"Typed schema and assurance")
 entries=d["sources"];need(type(entries)is list and len(entries)==len(SPECS),"Every mandatory global role")
 roots={"local":HERE,"repo":repo};out={};seen=set()
 for r in entries:
  need(type(r)is dict and set(r)=={"role","kind","path","relation","sha256"},"Exact source record schema")
  role=r["role"];need(type(role)is str and role in SPECS and role not in out,"Unique exact role")
  need((r["kind"],r["path"],r["relation"])==SPECS[role],"Exact source kind/path/logical-role")
  need(type(r["sha256"])is str and re.fullmatch("[0-9a-f]{64}",r["sha256"]),"SHA256")
  root=roots[r["kind"]].resolve();p=(root/r["path"]).resolve()
  need(p.is_relative_to(root)and p.is_file()and H(p)==r["sha256"],"Missing/stale source:"+role)
  key=(r["kind"],r["path"]);need(key not in seen,"Unique source path");seen.add(key);out[role]=r
 need(out["main_proof"]["sha256"]==SCOPE["main_reviewed_sha256"],"Exact whole-reviewed global mathematical proof")
 forbidden={"all_exponent_prefix_certificates.json","exponent_six_prefix_certificates.json",
  "ordinary_lower_one_finite_bridge_verification.json","relative_gap_three_finite_box_verification.json"}
 need(not any(Path(n).name in forbidden for _,n in seen),"Heavy completed boundary data are NOT logical inputs")
 return out
def child_boundaries(repo,rows):
 out=[]
 direct={(r["path"]if r["kind"]=="repo"else PREFIX+r["path"]):r["sha256"]for r in rows.values()}
 for tag,rec in CHILDREN.items():
  proof,checker,inp,ledger,receipt,key,count,status,addsupport,hashfield=rec
  d=load(HERE/ledger);need(type(d.get("schema"))is int and d["schema"]==1,"Child schema1:"+tag)
  entries=d.get(key);need(type(entries)is list and len(entries)==count,"Exact child source count:"+tag)
  seen=set();core=set()
  for r in entries:
   need(type(r)is dict and type(r.get("path"))is str,"Child source record")
   kind=r.get("kind",r.get("root",r.get("location",r.get("scope"))))
   need(kind in{"local","support","packet","repo"},"Child source root")
   source=(repo if kind=="repo"else HERE)/r["path"]
   need(source.is_file()and H(source)==r["sha256"],"Complete child source freshness:"+tag+":"+r["path"])
   canonical=(kind,r["path"]);need(canonical not in seen,"Unique child path");seen.add(canonical)
   name=r["path"]if kind=="repo"else PREFIX+r["path"]
   if name in direct:need(direct[name]==r["sha256"],"All overlapping direct/child identities agree")
   for role,name in [("proof",proof),("checker",checker),("inputs",inp)]:
    if kind in{"local","support","packet"}and r["path"]==name:
     need(r["sha256"]==rows[tag+"_"+role]["sha256"],"Direct child-core crosslink");core.add(role)
  need(core=={"proof","checker","inputs"},"Every regional core proof/checker/input mandatory")
  out.append({"tag":tag,"direct_child_sources":count,"all_current":True})
 # The explicit407 import for the six-way reverse is authenticated independently
 # of any inference from the global conclusion. G7 freshly replays it.
 d=load(HERE/SPECS["fixed_support_407_ledger"][1])
 need(type(d.get("schema"))is int and d["schema"]==1 and len(d["sources"])==77,"Complete407 fallback source declaration")
 for r in d["sources"]:
  path=(repo if r["kind"]=="repo"else HERE)/r["path"]
  need(path.is_file()and H(path)==r["sha256"],"Entire declared407 fallback boundary current")
  name=r["path"]if r["kind"]=="repo"else PREFIX+r["path"]
  if name in direct:need(direct[name]==r["sha256"],"407 reverse shares exact current raw inputs")
 g7=load(HERE/CHILDREN["supercritical"][3]);parent={}
 for r in g7["sources"]:
  if r["role"].startswith("fallback_407_"):parent[r["role"]]=r
 for ours,theirs in [("proof","proof"),("checker","checker"),("inputs","inputs"),("ledger","manifest"),("receipt","receipt")]:
  r=parent["fallback_407_"+theirs];need(r["sha256"]==rows["fixed_support_407_"+ours]["sha256"],"Six-way407 identity matches freshly executed G7 fallback")
 out.append({"tag":"fixed_support407","direct_child_sources":77,"freshly_executed_by_G7":True})
 return out
def replay(repo,rows,python):
 out=[]
 with tempfile.TemporaryDirectory(prefix="support-uniform-regions-")as folder:
  for tag,rec in CHILDREN.items():
   proof,checker,inp,ledger,receipt,key,count,status,addsupport,hashfield=rec
   output=Path(folder)/(tag+".json")
   cmd=[str(python),"-B",str(HERE/checker),"--repo-root",str(repo)]
   if addsupport:cmd+=["--support-root",str(HERE)]
   cmd+=["--output",str(output)]
   r=subprocess.run(cmd,cwd=folder,capture_output=True,text=True,timeout=180)
   need(r.returncode==0,"Fresh regional guard failed:"+tag+";"+r.stdout+";"+r.stderr)
   raw=output.read_bytes();need(raw==(HERE/receipt).read_bytes(),"Fresh regional receipt byte identity:"+tag)
   got=json.loads(raw)
   need(got["status"]==status and got[hashfield]==rows[tag+"_ledger"]["sha256"],"Fresh child exact source/PASS scope")
   if tag=="quadratic":need(got["scope"]==load(HERE/inp)["scope"],"401 exact domain")
   else:need(got["scope"]==load(HERE/inp),"Exact regional input scope")
   out.append({"tag":tag,"receipt_byte_identical":True,"receipt_sha256":hashlib.sha256(raw).hexdigest()})
 g7=load(HERE/CHILDREN["supercritical"][4])
 need(g7["fresh_bounded_gap_407_fallback"]["receipt_sha256"]==rows["fixed_support_407_receipt"]["sha256"],"407 evidence is actually freshly executed, not merely declared")
 return out


def import_kernel(repo):
 sys.dont_write_bytecode=True
 p=repo/SPECS["exact_polynomial_kernel"][1]
 spec=importlib.util.spec_from_file_location("global_cover_QQ_kernel",p);M=importlib.util.module_from_spec(spec);spec.loader.exec_module(M);return M
def universal(M):
 C,X,add,sub,mul,power,scale=M.C,M.X,M.add,M.sub,M.mul,M.power,M.scale
 out=[]
 c,x=X(0),X(1);e=mul(sub(c,C(Q(1,4))),power(x,2));CC=add(mul(c,power(x,2)),scale(x,Q(1,2)))
 p=sub(mul(e,add(scale(mul(c,x),2),C(1))),scale(mul(x,mul(CC,sub(c,C(Q(1,4))))),2))
 need(not p,"Universal EXACT r=(c-1/4)/(c+1/(2(d-1)))");out.append("exact_thermal_ratio_identity")
 d=X(0);p=sub(sub(scale(power(sub(d,C(1)),2),4),power(d,2)),mul(sub(d,C(2)),sub(scale(d,3),C(2))))
 need(not p,"Universal high-c factor4 threshold");out.append("K4e24_threshold_factorization")
 return out
def logical_controls():
 rules=[
 ({"TailQ","RawFixedN"},"BadFinite"),({"BadFinite"},"TailQ"),
 ({"TailQ"},"FieldTail"),({"FieldTail"},"ComplexTail"),({"ComplexTail"},"TailQ"),
 ({"BadFinite"},"ProjectiveFinite"),({"ProjectiveFinite"},"SupportFinite"),
 ({"SupportFinite","FixedSupport407"},"BadFinite")]
 predicates={"TailQ","BadFinite","FieldTail","ComplexTail","ProjectiveFinite","SupportFinite"}
 def close(seed,fn,fs):
  reached=set(seed)
  if fn:reached.add("RawFixedN")
  if fs:reached.add("FixedSupport407")
  while True:
   after=reached|{c for p,c in rules if p<=reached}
   if after==reached:return reached
   reached=after
 out=[]
 for p in sorted(predicates):
  reached=close({p},True,True)
  need(predicates<=reached,"Six-way reverse has both explicit context premises")
  out.append("six-way seeded:"+p)
 need(not(close(set(),True,True)&predicates),"An unseeded implication component cannot prove itself");out.append("no-cycle-bootstrap")
 for p in ["ProjectiveFinite","SupportFinite"]:
  reached=close({p},True,False)
  need("BadFinite"not in reached and"TailQ"not in reached,"Projection finiteness reverse NEEDS407")
  out.append("407-required:"+p)
 reached=close({"TailQ"},False,True)
 need("BadFinite"not in reached,"Uniform terminal tail alone NEEDS fixed-terminal fiber finiteness")
 out.append("raw-fixedN-required-for-tail-to-B")
 reached=close({"TailQ"},True,False)
 need("BadFinite"in reached,"Tail-to-B reverse DOES NOT need407")
 out.append("raw-fixedN-suffices-without407")
 for p in sorted(predicates):
  reached=close({p},True,True)
  need(not({"BadEmpty","AllOrder","EffectiveCutoff"}&reached),"No forbidden backwards strengthening")
  out.append("no-allorder-empty-effective:"+p)
 # Abstract examples, deliberately NOT moment counterexamples.
 fixed_terminal={(j,j,3)for j in range(1,65)}
 need(len(fixed_terminal)==64 and len({n for a,d,n in fixed_terminal})==1,"Tail alone permits infinitely many supports at an early terminal in an abstract model")
 out.append("abstract-fixedN-infinite-support-pattern")
 fixed_support={(1,2,n)for n in range(2,66)}
 need(len(fixed_support)==64 and len({(a,d)for a,d,n in fixed_support})==1,"Finite support projection does not by itself bound terminals")
 out.append("abstract-fixed-support-infinite-terminal-pattern")
 pairs={(j,j)for j in range(1,65)};ratios={Q(a,b)for a,b in pairs}
 need(len(pairs)==64 and ratios=={Q(1)},"Projective finite is NOT unscaled coefficient finite")
 out.append("projective-versus-unscaled")
 return out
def controls(inp):
 cfg=inp["controls"];critical=cover=large=rawtail=terminal=moments=appell=scalings=0
 C0=10**24;K=inp["K"];need(K==4*C0,"Exact K constant")
 for delta0 in cfg["critical_deltas"]:
  delta=Q(delta0);h=min(Q(1,8),delta/16)
  need(h>0 and h<=Q(1,8) and 8*h<=delta/2<delta,"Critical strip strict margin")
  for denominator in cfg["critical_denominators"]:
   for cc in [Q(1,4)-h,Q(1,4),Q(1,4)+h]:
    rr=(cc-Q(1,4))/(cc+Q(1,2*denominator))
    need(abs(rr)<=8*h<delta,"Whole CLOSED critical ratio interval covered")
    critical+=1
  for eps0 in cfg["small_epsilons"]:
   eps=min(Q(eps0),Q(1,16));left=Q(1,4)-h;right=Q(1,4)+h
   need(0<eps<left<=right<K,"Ordered finite covering intervals")
   samples=[eps/2,eps,(eps+left)/2,left,Q(1,4),right,(right+K)/2,K-1,K,K+1,2*K]
   for cc in samples:
    regions=[0<cc<=eps,eps<=cc<=left,left<=cc<=right,right<=cc<=K,cc>=K]
    need(any(regions),"Every positive-ratio sample belongs to an explicitly proved region")
    cover+=1
 for d in range(2,129):
  a=K*(d-1)**2
  need(4*(d-1)**2>=d*d and a>=C0*d*d,"G2 applies with its ACTUAL d>=2 domain")
  large+=1
 need(2**100>4*100**4 and Q(1,2)-Q(4,100)>0,"UNIVERSAL a0 monotone inequality base and derivative")
 for t in cfg["root_t_samples"]:
  need(t>=100 and 2**t>=4*t**4 and Q(1,2)-Q(4,t)>0,"Exact a0 uniform analytic/large split bound")
  N=t*t
  for d in [1,2,3,4,5,6,t-1,t,t+1,2*t]:
   if d==1:case="adjacent"
   elif d<=5:case="primepower_shape1"
   elif d*d<=N:case="raw_analytic"
   else:case="raw_large"
   if case=="raw_analytic":need(N>=10000 and 6<=d and d*d<=N,"Actual thm:analytic hypotheses")
   if case=="raw_large":need(2**d>=4*N*N,"Actual thm:large hypotheses")
   rawtail+=1
 for N in range(2,cfg["terminal_count_maximum"]+1):
  target=4*N*(N-1);ee=(target-1).bit_length();DN=2*(N-1)*ee
  need(2**(ee-1)<target<=2**ee,"Old high-gap exact integer cutoff")
  positive=sum(d*(2*N-3)-1 for d in range(1,DN))
  need(positive==(2*N-3)*(DN-1)*DN//2-(DN-1),"Positive-support fixed-terminal containing count")
  e0=(4*N*N-1).bit_length()
  need(2**(e0-1)<4*N*N<=2**e0 and e0>=1,"a0 raw large-gap fixed-terminal complement")
  terminal+=1
 def rise(a,k):
  r=1
  for j in range(k):r*=a+j
  return r
 for a in range(cfg["maximum_lower"]+1):
  for d in range(1,cfg["maximum_gap"]+1):
   for m in range(1,cfg["maximum_order"]+1):
    A,B=Q(2,3),-Q(3,5);x=A/B
    direct=sum(Q(comb(m,j)*factorial(a*m+d*j))*A**(m-j)*B**j for j in range(m+1))
    pol=sum(Q(comb(m,j))*rise(a*m+1,d*j)*x**(m-j)for j in range(m+1))
    need(direct==B**m*factorial(a*m)*pol,"Actual ordinary coefficient and own-order Gamma shape")
    lam=Q(3,2)
    scaled=sum(Q(comb(m,j)*factorial(a*m+d*j))*(lam*A)**(m-j)*(lam*B)**j for j in range(m+1))
    need(scaled==lam**m*direct,"Projective phase scaling preserves zeros")
    moments+=1;scalings+=1
 for d in range(1,17):
  for N in range(2,9):
   upper=[Q(comb(N,k)*factorial(d*k))for k in reversed(range(N+1))]
   lower=[Q(comb(N-1,k)*factorial(d*k))for k in reversed(range(N))]
   need([j*upper[j]for j in range(1,N+1)]==[N*v for v in lower],"EXACT a0 Appell derivative")
   appell+=1
 # An algebraic-ratio toy, explicitly not an actual bad factorial pair:
 # q=x²+1, p=(x+2)q. The common projective ratios are finite and need not be rational.
 q=[Q(1),Q(0),Q(1)];p=[Q(2),Q(1),Q(2),Q(1)]
 product=[Q(0)]*4
 for i,a in enumerate(q):
  for j,b in enumerate([Q(2),Q(1)]):product[i+j]+=a*b
 need(product==p and len(q)-1==2,"Degree-three/two algebraic projective ratio toy")
 # Actual mathematical negative boundaries.
 need(2*rise(3,2)!=2*rise(2,2),"Appell derivative must NOT be moved to nonzero lower exponent")
 need(Q(4,(5-1)**2)==Q(1,4)and Q(4,5**2)<Q(1,4),"Wrong c normalization misclassifies actual cusp")
 def quadratic_domain(a,d,N):return d>=2 and N>=2 and a>=10**24*d*d
 need(not quadratic_domain(10**24,1,2) and quadratic_domain(4*10**24,2,2),"Actual G2 domain excludes d1; its coverage comes from G0")
 return {"critical_strip_margin_controls":critical,"finite_cover_controls":cover,
 "K4e24_quadratic_threshold_controls":large,"raw_constant_uniform_tail_controls":rawtail,
 "raw_fixed_terminal_support_counts":terminal,"actual_moment_normalizations":moments,
 "actual_projective_scaling_identities":scalings,"a0_Appell_vector_identities":appell,
 "generic_algebraic_projective_ratio_toy":1,"genuine_math_negative_controls":2,"explicit_G2_domain_controls":1}
def negatives(ledger,repo):
 out=[]
 def reject(name,fn):
  try:fn()
  except(ValueError,AssertionError,KeyError,TypeError,OSError):out.append(name);return
  raise ValueError("Negative accepted:"+name)
 for role in SPECS:
  d=copy.deepcopy(ledger);d["sources"]=[r for r in d["sources"]if r["role"]!=role]
  reject("omit:"+role,lambda d=d:check_sources(d,repo))
  d=copy.deepcopy(ledger);next(r for r in d["sources"]if r["role"]==role)["sha256"]="0"*64
  reject("badpin:"+role,lambda d=d:check_sources(d,repo))
 for key in ["a0_92935_or9998_finite_certificate_import","unscaled_coefficient_pairs_finite","all_order_SFC",
  "global_bad_set_empty","effective_Nstar_or_containing_box","replication_of_one_bad_triple",
  "arbitrary_polynomial_supports","general_Gamma_shapes","global_analytic_estimates_newly_proved_by_finite_controls",
  "heavy_finite_boundary_runs_repeated"]:
  d=copy.deepcopy(SCOPE);d[key]=True;reject("false-scope:"+key,lambda d=d:check_inputs(d))
 for key in ["all_nonzero_coefficients_and_characteristic_zero_fields","G5_used_at_ONE_fixed_finite_U",
  "G6_sequential_to_fixed_relative_strip","global_bad_triple_set_finite","finite_projective_exceptional_phases",
  "bad_support_projection_finite","six_way_equivalence_under_declared_inputs",
  "QQ_tail_to_finite_B_requires_raw_fixed_terminal_input_NOT407",
  "finite_support_or_projective_phases_to_finite_B_REQUIRES407"]:
  d=copy.deepcopy(SCOPE);d[key]=False;reject("missing-scope:"+key,lambda d=d:check_inputs(d))
 for key,value in [("K",10**24),("quadratic_G2_domain","d>=1"),("parameters","c=a/d^2 for all supports")]:
  d=copy.deepcopy(SCOPE);d[key]=value;reject("wrong-constant-or-domain:"+key,lambda d=d:check_inputs(d))
 d=copy.deepcopy(SCOPE);d["schema"]=True;reject("typed-schema",lambda:check_inputs(d))
 d=copy.deepcopy(ledger);d["sources"]=[];reject("empty-source-ledger",lambda:check_sources(d,repo))
 return out
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--repo-root",required=True,type=Path)
 ap.add_argument("--parent-python",type=Path);ap.add_argument("--output",type=Path,default=HERE/RN);ns=ap.parse_args()
 repo=ns.repo_root.resolve();need(repo.is_dir(),"Explicit repository root")
 python=ns.parent_python
 if python is None:python=next((p for p in [repo/".venv/bin/python",repo/".venv/Scripts/python.exe"]if p.is_file()),Path(sys.executable))
 python=python.absolute();need(python.is_file(),"Parent interpreter path; .absolute() preserves the venv symlink")
 inp=load(HERE/IN);check_inputs(inp);ledger=load(HERE/LN);rows=check_sources(ledger,repo)
 boundaries=child_boundaries(repo,rows);children=replay(repo,rows,python)
 identities=universal(import_kernel(repo));logic=logical_controls();counts=controls(inp);bad=negatives(ledger,repo)
 check_sources(ledger,repo)
 result={"schema":1,"status":"PASS","scope":inp,"proof_sha256":H(HERE/PN),"checker_sha256":H(HERE/CN),
 "input_sha256":H(HERE/IN),"source_ledger_sha256":H(HERE/LN),"mandatory_source_roles":len(rows),
 "sources":[rows[r]for r in sorted(rows)],"complete_child_boundaries":boundaries,"fresh_regional_replays":children,
 "fixed_support_407_freshly_replayed_inside_G7":True,"universal_cover_identities":identities,
 "logical_six_way_dependency_controls":logic,"counts":counts,"source_scope_negative_controls":bad,
 "assurance_flags":ASSURANCE,
 "meaning":"The whole written regional theorems and elementary finite cover prove uniform eventuality. Fresh source-closed replays authenticate the exact reused packets. Finite QQ/logic controls are supporting evidence, not a formal analytic proof, effective box, all-order theorem or proof that the finite bad set is empty."}
 ns.output.parent.mkdir(parents=True,exist_ok=True);ns.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
 print(json.dumps({"status":"PASS","mandatory_sources":len(rows),"fresh_regions":len(children),
 "universal_identities":len(identities),"logical_controls":len(logic),"counts":counts,
 "source_scope_controls":len(bad),"output":str(ns.output)}))
if __name__=="__main__":main()

