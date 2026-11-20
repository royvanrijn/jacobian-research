#!/usr/bin/env python3
"""Exact/source support for compact-subcritical ratio uniform eventuality."""
from __future__ import annotations
if not __debug__:raise SystemExit("Refusing -O/-OO: exact controls must execute.")
from pathlib import Path
from fractions import Fraction as Q
from math import comb,factorial
import argparse,copy,hashlib,importlib.util,json,re,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
PN="COMPACT_SUBCRITICAL_RATIO_UNIFORM_EVENTUAL_SFC.md"
CN="verify_compact_subcritical_ratio_uniform.py"
IN="compact_subcritical_ratio_inputs.json"
LN="compact_subcritical_ratio_source_ledger.json"
RN="compact_subcritical_ratio_verification.json"
PREFIX="research/factorial/image-rigidity/"
PARENT_SPECS={
 "proof":("local","SUBCRITICAL_ORDERED_PAIR_EVENTUAL_SFC.md"),
 "checker":("local","verify_strict_subcritical_eventual_sfc.py"),
 "inputs":("local","strict_subcritical_inputs.json"),
 "two_thimble_proof":("support","CUSP_TWO_THIMBLE_POSITIVE_CYCLE.md"),
 "global_complement_proof":("support","CUSP_GLOBAL_REGULAR_COMPLEMENT_AND_STOKES_DETOUR.md"),
 "retained_gamma_upper_bound":("repo",PREFIX+"GAMMA_PRODUCT_ROOT_MODULUS.tex")}
SPECS={
 "main_proof":("local",PN),"exact_checker":("local",CN),"exact_inputs":("local",IN),
 "parent_ledger":("support","strict_subcritical_source_ledger.json"),
 "parent_receipt":("support","strict_subcritical_verification.json"),
 "raw_old_high_gap":("repo",PREFIX+"lower_exponent_initial_and_boundary.tex"),
 "primary_kabluchko_context":("repo",PREFIX+"uniform_cusp_sources/kabluchko_2203_05533v3.pdf")}
for role,(kind,name) in PARENT_SPECS.items():SPECS["parent_"+role]=("repo" if kind=="repo" else "support",name)
ASSURANCE={"external_review":False,"formal_verification":False,"independent_replay":False}
SCOPE={"schema":1,
 "parameters":"integers a>=1,d>=2; fixed real 0<c_minus<=c_plus<1/4; c=a/(d-1)^2",
 "uniformity":"one non-effective terminal cutoff for every support in the fixed ratio interval",
 "same_cutoff_all_fields_and_coefficients":True,
 "entire_original_integer_integrands":True,
 "actual_leading_actions_no_m_over_d_restriction":True,
 "both_actual_endpoint_circle_capture_before_lambda_limit":True,
 "zero_disc_gaussian_replaces_lower_annulus":True,
 "model_all_circle_two_simple_global_maxima_including_q1":True,
 "exact_integer_partial_periods_not_model_circlelaw_transfer":True,
 "regional_finite_bad_iff_uniform_cutoff_and_fields":True,
 "old_high_gap_and_old_soft_only_for_regional_reverse":True,
 "kabluchko_2_4_is_off_circle_context_only":True,
 "fixed_d_fallback_405_strict_only":True,
 "uses_401_404_407_411_412":False,
 "uniform_as_c_tends_zero_or_quarter":False,
 "compact_supercritical_ratio_proved":False,
 "effective_cutoff":False,"all_order_SFC":False,"global_bad_triple_set_finite":False,
 "heavy_finite_prefix_replay":False,
 "reviewed_body_sha256":"fd1f3bd7680fc8d342b2fffe1676912d3f65d68aec8c89c38aacf2ba603cc014",
 "controls":{"intervals":[["1/100","1/8"],["1/16","1/5"],["1/8","3/16"]],
             "maximum_gap":16,"maximum_lower":48,"maximum_order":4,
             "tail_gaps":[129,257,1025],"fixed_terminal_maximum":64,"series_coefficients":16},
 "assurance_flags":ASSURANCE}
def need(ok,msg):
 if not ok:raise ValueError(msg)
def H(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def typed(a,b):
 if type(a) is not type(b):return False
 if type(b) is dict:return set(a)==set(b) and all(typed(a[k],v)for k,v in b.items())
 if type(b) is list:return len(a)==len(b) and all(typed(x,y)for x,y in zip(a,b))
 return a==b
def check_inputs(d):need(typed(d,SCOPE),"Exact typed compact-subcritical scope")
def check_sources(d,repo,support):
 need(type(d) is dict and set(d)=={"schema","assurance_flags","sources"},"Ledger schema")
 need(type(d["schema"]) is int and d["schema"]==1 and typed(d["assurance_flags"],ASSURANCE),"Typed schema and assurance")
 entries=d["sources"];need(type(entries) is list and len(entries)==len(SPECS)==13,"All thirteen roles")
 roots={"local":HERE,"support":support,"repo":repo};out={};paths=set()
 for row in entries:
  need(type(row) is dict and set(row)=={"role","kind","path","sha256"},"Row schema")
  role=row["role"];need(type(role) is str and role in SPECS and role not in out,"Unique exact role")
  need((row["kind"],row["path"])==SPECS[role],"Exact role selection")
  need(type(row["sha256"]) is str and re.fullmatch("[0-9a-f]{64}",row["sha256"]),"SHA256")
  root=roots[row["kind"]].resolve();p=(root/row["path"]).resolve()
  need(p.is_relative_to(root) and p.is_file() and H(p)==row["sha256"],"Missing/stale source:"+role)
  key=(row["kind"],row["path"]);need(key not in paths,"Unique path");paths.add(key);out[role]=row
 forbidden={"ALL_ORDINARY_BINOMIAL_EVENTUAL_STRONG_FACTORIAL.md","FIBERWISE_UNIFORM_EVENTUAL_BINOMIAL_SFC.md",
 "UNIFORM_QUADRATIC_STRONG_FACTORIAL_TAIL.md","UNIFORM_INTEGER_CUSP_EVENTUAL_STRONG_FACTORIAL.md",
 "GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md","all_exponent_prefix_certificates.json"}
 need(not any(Path(p).name in forbidden for _,p in paths),"Unneeded broad or sharp-lower route excluded")
 return out
def parent(repo,support,rows):
 sys.dont_write_bytecode=True
 path=support/PARENT_SPECS["checker"][1]
 spec=importlib.util.spec_from_file_location("compact_ratio_strict_parent",path)
 module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 need(module.SOURCES==PARENT_SPECS,"Exact parent source kernel")
 ledger=load(support/SPECS["parent_ledger"][1]);module.validate_ledger(ledger,repo,support)
 need(len(ledger["sources"])==6,"Exact six parent sources")
 for row in ledger["sources"]:
  need(row["sha256"]==rows["parent_"+row["role"]]["sha256"],"Direct and parent source pins agree")
 with tempfile.TemporaryDirectory(prefix="compact-subcritical-parent-")as folder:
  output=Path(folder)/"actual.json"
  r=subprocess.run([sys.executable,"-B",str(path),"--repo-root",str(repo),"--support-root",str(support),"--output",str(output)],
   cwd=folder,capture_output=True,text=True,timeout=60)
  need(r.returncode==0,"Fresh narrow parent failed:"+r.stdout+";"+r.stderr)
  raw=output.read_bytes();need(raw==(support/SPECS["parent_receipt"][1]).read_bytes(),"Parent receipt byte identity")
  d=json.loads(raw);need(d["status"]=="PASS" and d["source_ledger_sha256"]==rows["parent_ledger"]["sha256"],"Parent scope/ledger")
 return module,{"source_roles":6,"receipt_byte_identical":True,"receipt_sha256":hashlib.sha256(raw).hexdigest()}
def universal(M):
 C,X,add,sub,mul,power,scale=M.C,M.X,M.add,M.sub,M.mul,M.power,M.scale
 out=[]
 def ck(name,p):need(not p,"Universal polynomial identity:"+name);out.append(name)
 sigma,c,t,q,d,L=[X(j)for j in range(6)]
 one=C(1);two=C(2)
 # a*sigma^2=c(1-sigma)^2; the cancellation leaves -(1-sigma)*t.
 a_sigma2=mul(c,power(sub(one,sigma),2))
 L_sigma=add(scale(a_sigma2,2),mul(sigma,sub(one,sigma)))
 ck("exact central linear cancellation",sub(sub(scale(a_sigma2,2),L_sigma),scale(mul(sigma,sub(one,sigma)),-1)))
 # L/d-c=c(1-2sigma)^2+sigma[1-(1+2c)sigma].
 lhs=sub(L_sigma,c)
 rhs=add(mul(c,power(sub(one,scale(sigma,2)),2)),
         mul(sigma,sub(one,mul(add(one,scale(c,2)),sigma))))
 ck("uniform tail scale lower bound decomposition",sub(lhs,rhs))
 beta,E,t=[X(j)for j in range(3)]
 N=add(scale(mul(E,mul(sub(beta,t),add(beta,t))),2),scale(mul(E,add(beta,t)),-1),scale(mul(E,sub(beta,t)),-1))
 ck("model stationary derivative numerator",sub(N,scale(mul(E,sub(mul(beta,sub(beta,one)),power(t,2))),2)))
 A,B,b,z=[X(j)for j in range(4)]
 left=sub(mul(A,add(B,b)),scale(mul(z,add(one,mul(b,B))),2))
 right=mul(add(B,one),sub(A,scale(z,2)))
 ck("horizontal curvature strictness decomposition",sub(sub(left,right),mul(sub(b,one),sub(A,scale(mul(z,B),2)))))
 # Actual p=1+2t/d and C=L*d/2, a=C-(d-1)/2.
 d,L,p,HH=[X(j)for j in range(4)]
 a=scale(sub(mul(L,d),sub(d,one)),Q(1,2))
 C0=scale(mul(L,d),Q(1,2))
 raw=add(a,scale(mul(C0,p),-1),mul(d,HH))
 normalized=add(scale(sub(d,one),Q(-1,2)),scale(mul(C0,sub(p,one)),-1),mul(d,HH))
 ck("exact whole-tail numerator",sub(raw,normalized))
 S,a,d,w,Z=[X(j)for j in range(5)]
 crit=add(mul(sub(a,S),add(w,Z)),mul(d,Z))
 ck("stationary actual radial envelope",add(sub(mul(d,w),mul(sub(add(a,d),S),add(w,Z))),crit))
 a,d,m,j=[X(k)for k in range(4)]
 ck("old Gamma upper factor comparison",sub(add(mul(a,m),one,mul(d,sub(m,one)),j),
  sub(mul(add(a,d),m),sub(sub(d,one),j))))
 m,j=[X(k)for k in range(2)]
 ck("Gaussian unitary-Hermite coefficient exponent",sub(power(sub(m,scale(j,2)),2),sub(power(m,2),scale(mul(j,sub(m,j)),4))))
 return out
def controls(inp):
 cfg=inp["controls"];ints=[(Q(a),Q(b))for a,b in cfg["intervals"]]
 need(all(0<a<=b<Q(1,4)for a,b in ints),"Actual compact intervals")
 moments=own=tails=factors=terminal=series=models=0
 def rise(x,n):
  r=1
  for j in range(n):r*=x+j
  return r
 for d in range(2,cfg["maximum_gap"]+1):
  for a in range(1,cfg["maximum_lower"]+1):
   c=Q(a,(d-1)**2)
   if not any(lo<=c<=hi for lo,hi in ints):continue
   s=Q(d-1,2);C0=a+s;L=2*C0/d;kap=Q(d+1,d-1);w0=C0**d*kap
   need(4*a<(d-1)**2,"Exact strict405 fallback domain")
   for m in range(1,cfg["maximum_order"]+1):
    AA,BB=Q(2,3),-Q(3,5);x=AA/BB
    direct=sum(Q(comb(m,j)*factorial(a*m+d*j))*AA**(m-j)*BB**j for j in range(m+1))
    pol=sum(comb(m,j)*rise(a*m+1,d*j)*x**(m-j)for j in range(m+1))
    need(direct==BB**m*factorial(a*m)*pol,"Actual moment identity")
    moments+=1
    for j in range(d):
     factor=a*m+1+d*(m-1)+j
     need(factor==(a+d)*m-(d-1-j) and 0<factor<=(a+d)*m,"Actual old upper factor")
     factors+=1
   for N in [2,3,5,8]:
    n=N-1;qn=Q(2,3)/(n**d*w0);qN=Q(2,3)/(N**d*w0)
    need(qn==Q(N,n)**d*qN,"Both own-order arguments")
    own+=1
 for lo,hi in ints:
  R=max(Q(1),3/lo)
  for d in cfg["tail_gaps"]:
   sig=Q(1,d)
   for c in [lo,hi]:
    a=c*(d-1)**2;C0=a+Q(d-1,2);L=2*C0/d
    need(L/d>=lo,"Uniform tail scale")
    for t in [R,R+1,2*R,-R,-R-1,-2*R]:
     pp=1+2*t/d
     if pp<=0:continue
     for HH in [Q(0),Q(1,2),Q(1)]:
      actual=2*(-Q(d-1,2)-L*t+d*HH)/(d*pp)
      if t>=R:need(actual<=-lo/2,"Right whole original tail")
      else:need(actual>=1,"Left whole original tail")
      tails+=1
 for N in range(2,cfg["fixed_terminal_maximum"]+1):
  target=4*N*(N-1);ex=(target-1).bit_length();DN=2*(N-1)*ex
  need(2**(ex-1)<target<=2**ex,"Exact raw high-gap cutoff")
  count=sum(d*(2*N-3)-1 for d in range(2,DN))
  need(count==(2*N-3)*(DN*(DN-1)//2-1)-(DN-2),"Finite-terminal old support count")
  terminal+=1
 for j in range(1,cfg["series_coefficients"]+1):
  coeff=Q(2**(2*j+1),factorial(2*j+1))
  need(coeff>0,"sinh(2z)-2z coefficient")
  need(Q(2**(2*j+1),factorial(2*j))-coeff==2*j*coeff>0,
       "2z cosh(2z)-sinh(2z) coefficient")
  series+=2
 for n in range(12,81):
  c_left=Q(2*n*n,(n*n-1)**2)
  c_right=Q(n*n-1,4*n*n)
  need(0<c_left<Q(1,n) and 0<Q(1,4)-c_right<Q(1,n),"Abstract ratios can escape both compact boundaries")
  need(n*n<3*(2*n*n),"Left escaping model respects the earlier relative-gap complement")
  models+=2
 # Genuine scope/arithmetic negatives, not numerical claims about roots.
 need(Q(1,4) not in [hi for lo,hi in ints],"Quartic cusp excluded")
 need(Q(4,(5-1)**2)==Q(1,4) and Q(4,5**2)<Q(1,4),"Wrong d-squared ratio misclassifies an ACTUAL integer cusp")
 need(1+2*Q(-2)/2<=0,"Original left endpoint domain must not be crossed")
 return {"actual_moment_controls":moments,"actual_own_order_controls":own,
  "old_gamma_upper_factor_controls":factors,"whole_original_tail_inequalities":tails,
  "fixed_terminal_counts":terminal,"positive_hyperbolic_series_controls":series,
  "abstract_endpoint_escape_controls":models,"genuine_math_negative_controls":3}
def negatives(ledger,repo,support):
 out=[]
 def reject(label,fn):
  try:fn()
  except(ValueError,KeyError,TypeError,OSError):out.append(label);return
  raise ValueError("Negative accepted:"+label)
 for role in SPECS:
  doc=copy.deepcopy(ledger);doc["sources"]=[r for r in doc["sources"]if r["role"]!=role]
  reject("omit:"+role,lambda doc=doc:check_sources(doc,repo,support))
  doc=copy.deepcopy(ledger);next(r for r in doc["sources"]if r["role"]==role)["sha256"]="0"*64
  reject("badpin:"+role,lambda doc=doc:check_sources(doc,repo,support))
 for key in ["uses_401_404_407_411_412","uniform_as_c_tends_zero_or_quarter","compact_supercritical_ratio_proved",
             "effective_cutoff","all_order_SFC","global_bad_triple_set_finite","heavy_finite_prefix_replay"]:
  doc=copy.deepcopy(SCOPE);doc[key]=True;reject("false-scope:"+key,lambda doc=doc:check_inputs(doc))
 for key in ["zero_disc_gaussian_replaces_lower_annulus","both_actual_endpoint_circle_capture_before_lambda_limit",
             "regional_finite_bad_iff_uniform_cutoff_and_fields","actual_leading_actions_no_m_over_d_restriction"]:
  doc=copy.deepcopy(SCOPE);doc[key]=False;reject("missing-scope:"+key,lambda doc=doc:check_inputs(doc))
 doc=copy.deepcopy(SCOPE);doc["schema"]=True;reject("typed-schema",lambda:check_inputs(doc))
 doc=copy.deepcopy(ledger);doc["sources"]=[];reject("empty-ledger",lambda:check_sources(doc,repo,support))
 return out
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--repo-root",type=Path,required=True)
 ap.add_argument("--support-root",type=Path,default=HERE);ap.add_argument("--output",type=Path,default=HERE/RN);ns=ap.parse_args()
 repo,support=ns.repo_root.resolve(),ns.support_root.resolve()
 inp=load(HERE/IN);check_inputs(inp);ledger=load(HERE/LN);rows=check_sources(ledger,repo,support)
 M,upstream=parent(repo,support,rows);ids=universal(M);counts=controls(inp);bad=negatives(ledger,repo,support)
 check_sources(ledger,repo,support)
 result={"schema":1,"status":"PASS","scope":inp,"proof_sha256":H(HERE/PN),"checker_sha256":H(HERE/CN),
  "input_sha256":H(HERE/IN),"source_ledger_sha256":H(HERE/LN),"sources":[rows[r]for r in sorted(rows)],
  "mandatory_source_roles":len(rows),"fresh_narrow405_parent":upstream,"universal_identities":ids,
  "counts":counts,"source_scope_negative_controls":bad,"assurance_flags":ASSURANCE,
  "meaning":"Written actual whole-tail/packet estimates supply compact-parameter uniformity. Finite algebra, arithmetic and current-source controls support that proof, not formal verification or an effective cutoff."}
 ns.output.parent.mkdir(parents=True,exist_ok=True);ns.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
 print(json.dumps({"status":"PASS","mandatory_sources":len(rows),"universal_identities":len(ids),"counts":counts,
                   "source_scope_controls":len(bad),"output":str(ns.output)}))
if __name__=="__main__":main()

