#!/usr/bin/env python3
"""Exact and source support for compact-supercritical uniform eventual SFC."""
from __future__ import annotations
if not __debug__:raise SystemExit("Refusing -O/-OO: exact/source controls must execute.")
from pathlib import Path
from fractions import Fraction as Q
from math import factorial,comb
import argparse,copy,hashlib,importlib.util,json,re,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
PN="COMPACT_SUPERCRITICAL_RATIO_UNIFORM_EVENTUAL_SFC.md"
CN="verify_compact_supercritical_ratio_uniform.py"
IN="compact_supercritical_ratio_universal_ledger.json"
LN="compact_supercritical_ratio_source_manifest.json"
RN="compact_supercritical_ratio_verification.json"
PREFIX="research/factorial/image-rigidity/"
SPECS={
 "main_proof":("local",PN,"written_complete_actual_period_proof"),
 "exact_checker":("local",CN,"exact_and_source_support"),
 "universal_ledger":("local",IN,"typed_universal_identity_and_scope_ledger"),
 "fallback_407_proof":("support","ALL_ORDINARY_BINOMIAL_EVENTUAL_STRONG_FACTORIAL.md","logical_input_ONLY_bounded_gap_finite_support_fallback"),
 "fallback_407_checker":("support","verify_ordinary_binomial_eventual_completion.py","defining_existing_407_guard"),
 "fallback_407_inputs":("support","ordinary_binomial_eventual_inputs.json","defining_existing_407_scope"),
 "fallback_407_manifest":("support","ordinary_binomial_eventual_source_ledger.json","hierarchical_complete_77_source_boundary"),
 "fallback_407_receipt":("support","ordinary_binomial_eventual_verification.json","fresh_byte_identical_existing_407_replay"),
 "old_gamma_upper":("repo",PREFIX+"GAMMA_PRODUCT_ROOT_MODULUS.tex","logical_input_OLD_UPPER_bound_only_no_lower_annulus"),
 "simple_turn_mechanism":("repo",PREFIX+"SIMPLE_TURN_SHRINKING_WINDOW_CONSECUTIVE_EXCLUSION.md","logical_local_analytic_cubic_and_proper_arm_mechanisms_revalidated_here"),
 "large_thermal_cubic_mechanism":("repo",PREFIX+"LARGE_THERMAL_QUARTIC_FOLD_MATCHING.md","logical_ONLY_sections5_8_local_cubic_proper_arms_Cauchy_NOT_large_thermal_theorem"),
 "primary_kabluchko":("repo",PREFIX+"uniform_cusp_sources/kabluchko_2203_05533v3.pdf","logical_Theorem2_3_inverse_existence_ONLY_Theorem2_4_circlelaw_context"),
 "nist_airy_sector":("repo",PREFIX+"sources/dlmf-9-7-endpoint.html","logical_classical_Ai_sector_9_7_5"),
 "nist_airy_zeros":("repo",PREFIX+"sources/dlmf-9-9-endpoint.html","logical_classical_Ai_zero_location_9_9_i"),
 "nist_airy_ODE":("repo",PREFIX+"sources/dlmf-9-2-endpoint.html","logical_Airy_ODE_uniqueness_simplicity"),
 "nist_CFU_context":("repo",PREFIX+"sources/dlmf-2-4-endpoint.html","classical_CFU_context_does_NOT_select_actual_cycle"),
 "exact_polynomial_kernel":("repo",PREFIX+"verify_strict_subcritical_eventual_sfc.py","code_ONLY_sparse_Q_polynomial_kernel_NOT_subcritical_theorem")}
ASSURANCE={"formal_verification":False,"external_review":False,"independent_replay":False}
SCOPE={
 "schema":1,"parameters":"integers a>=1,d>=2; fixed real 1/4<c_minus<=c_plus<infinity; c=a/(d-1)^2",
 "large_gap_core":"there exist d0(I),N0(I); every d>=d0,N>=N0 on I is good",
 "full_corollary":"one unevaluated terminal cutoff for every support on the fixed ratio interval",
 "same_cutoff_all_nonzero_coefficients_and_characteristic_zero_fields":True,
 "large_gap_core_independent_407":True,"407_role":"ONLY finitely many bounded-d supports",
 "Kabluchko_2_3_logical_inverse_existence":True,
 "Kabluchko_2_4_and_circlelaw_context_only":True,
 "actual_zero_free_log_strip_and_complex_order_partials":True,
 "both_actual_endpoint_capture_before_step_or_label_assumptions":True,
 "actual_leading_actions_no_m_over_d_restriction":True,
 "right_single_Ai_class_NOT_Scorer":True,
 "proper_cubic_shared_left_sheet_low_midpoint":True,
 "bounded_fields_same_finite_Ai_label_and_derivative_free_comparison":True,
 "large_thermal_role":"ONLY local cubic/proper-arm/Cauchy construction with hypotheses reverified; NOT its sequential theorem",
 "raw_large_thermal_expected_sha256":"a9ca7dfb62b814d0637d1a56baac5084f951b479e8e64f9e5df5b36948d13c87",
 "primary_Kabluchko_v3_sha256":"fbf2a26f4424761bf0a75ac2846b27aade3cd4c65433ff11fdf6ca929073a02c",
 "reviewed_main_expected_sha256":"0a4945c8a491c20a628db6e1fd37af555b70fa5e4b850d4202087c6bc5c4219f",
 "sharp_lower_annulus_premise":False,"model_circlelaw_to_Gamma_transfer":False,
 "large_thermal_theorem_premise":False,"all_order_SFC":False,"global_bad_triple_set_finite":False,
 "effective_cutoff":False,"uniform_to_cusp_or_infinite_ratio":False,
 "heavy_finite_boundary_computation_rerun":False,
 "controls":{"intervals":[["1/3","1"],["1","4"],["4","16"]],
  "maximum_gap":12,"maximum_lower":48,"maximum_order":4,
  "pythagorean_parameter_max_denominator":12,"positive_series_terms":16},
 "assurance_flags":ASSURANCE}
def need(ok,msg):
 if not ok:raise ValueError(msg)
def H(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def typed(a,b):
 if type(a) is not type(b):return False
 if type(b) is dict:return set(a)==set(b)and all(typed(a[k],v)for k,v in b.items())
 if type(b) is list:return len(a)==len(b)and all(typed(x,y)for x,y in zip(a,b))
 return a==b
def check_inputs(d):need(typed(d,SCOPE),"Exact compact-supercritical scope and universal ledger")
def check_sources(d,repo,support):
 need(type(d)is dict and set(d)=={"schema","assurance_flags","sources"},"Source manifest schema")
 need(type(d["schema"])is int and d["schema"]==1,"Integer schema1")
 need(typed(d["assurance_flags"],ASSURANCE),"False typed assurance flags")
 need(type(d["sources"])is list and len(d["sources"])==len(SPECS)==17,"Every mandatory source role")
 roots={"local":HERE,"support":support,"repo":repo};out={};seen=set()
 for row in d["sources"]:
  need(type(row)is dict and set(row)=={"role","kind","path","relation","sha256"},"Typed source row")
  role=row["role"];need(type(role)is str and role in SPECS and role not in out,"Unique exact role")
  need((row["kind"],row["path"],row["relation"])==SPECS[role],"Exact kind/path/logical role")
  need(type(row["sha256"])is str and re.fullmatch("[0-9a-f]{64}",row["sha256"]),"SHA256")
  root=roots[row["kind"]].resolve();p=(root/row["path"]).resolve()
  need(p.is_relative_to(root)and p.is_file()and H(p)==row["sha256"],"Missing/stale source:"+role)
  key=(row["kind"],row["path"]);need(key not in seen,"Unique source path");seen.add(key);out[role]=row
 need(out["primary_kabluchko"]["sha256"]==SCOPE["primary_Kabluchko_v3_sha256"],"Exact inspected version3 primary PDF")
 need(out["main_proof"]["sha256"]==SCOPE["reviewed_main_expected_sha256"],"Frozen whole-reviewed main proof")
 need(out["large_thermal_cubic_mechanism"]["sha256"]==SCOPE["raw_large_thermal_expected_sha256"],"Exact independently audited local mechanism source")
 forbidden={"all_exponent_prefix_certificates.json","ordinary_lower_one_finite_bridge_verification.json",
 "relative_gap_three_finite_box_verification.json","GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md"}
 need(not any(Path(name).name in forbidden for _,name in seen),"Heavy boundary or new lower annulus route excluded")
 return out
def import_module(name,path):
 sys.dont_write_bytecode=True
 spec=importlib.util.spec_from_file_location(name,path);M=importlib.util.module_from_spec(spec);spec.loader.exec_module(M);return M
def parent(repo,support,rows):
 M=import_module("supercritical_bounded_gap_407",support/SPECS["fallback_407_checker"][1])
 ledger=load(support/SPECS["fallback_407_manifest"][1])
 need(type(ledger.get("schema"))is int and ledger["schema"]==1 and len(ledger["sources"])==77,"Exact complete retained407 source cardinality")
 actual=M.check_sources(ledger,repo,support);linked=M.links(repo,support,actual)
 need(len(actual)==77 and len(linked)==12,"All hierarchical407 boundaries")
 for outer,inner in [("fallback_407_proof","main_proof"),("fallback_407_checker","main_checker"),
  ("fallback_407_inputs","main_inputs"),("old_gamma_upper","gap405_retained_gamma_upper"),
  ("simple_turn_mechanism","turn_proof"),("exact_polynomial_kernel","gap405_strict_checker")]:
  need(rows[outer]["sha256"]==actual[inner]["sha256"],"Direct and parent pins agree:"+outer)
 with tempfile.TemporaryDirectory(prefix="supercritical-407-fallback-")as folder:
  output=Path(folder)/"actual.json"
  r=subprocess.run([sys.executable,"-B",str(support/SPECS["fallback_407_checker"][1]),
    "--repo-root",str(repo),"--support-root",str(support),"--output",str(output)],
    cwd=folder,capture_output=True,text=True,timeout=90)
  need(r.returncode==0,"Fresh bounded-gap407 parent failed:"+r.stdout+";"+r.stderr)
  raw=output.read_bytes();need(raw==(support/SPECS["fallback_407_receipt"][1]).read_bytes(),"Fresh407 receipt byte identity")
  got=json.loads(raw)
  need(got["status"]=="PASS"and got["source_ledger_sha256"]==rows["fallback_407_manifest"]["sha256"],"Retained407 PASS/source boundary")
  for key,role in [("proof_sha256","fallback_407_proof"),("checker_sha256","fallback_407_checker"),("input_sha256","fallback_407_inputs")]:
   need(got[key]==rows[role]["sha256"],"Fresh parent identity:"+key)
 return {"mandatory_parent_sources":77,"crosslinked_parent_ledgers":12,"fresh_children_inside_parent":4,
  "receipt_byte_identical":True,"receipt_sha256":hashlib.sha256(raw).hexdigest(),
  "logical_role":"ONLY bounded-d finite-support fallback in full corollary"}


def universal(M):
 C,X,add,sub,mul,power,scale=M.C,M.X,M.add,M.sub,M.mul,M.power,M.scale
 out=[]
 def ck(name,p):need(not p,"Universal identity:"+name);out.append(name)
 def g(p):return(p,{})
 def gi(p):return({},p)
 def ga(*ps):return(add(*(z[0]for z in ps)),add(*(z[1]for z in ps)))
 def gs(p,q):return(sub(p[0],q[0]),sub(p[1],q[1]))
 def gm(p,q):return(sub(mul(p[0],q[0]),mul(p[1],q[1])),add(mul(p[0],q[1]),mul(p[1],q[0])))
 def gc(p,c):return(scale(p[0],c),scale(p[1],c))
 def gp(p,n):
  q=g(C(1))
  for _ in range(n):q=gm(q,p)
  return q
 def ckG(name,p):need(not p[0]and not p[1],"Universal Gaussian identity:"+name);out.append(name)
 def diff(p,index):
  z={}
  for ex,a in p.items():
   if ex[index]:
    ey=list(ex);ey[index]-=1;z[tuple(ey)]=a*ex[index]
  return z
 def rb2(p):
  z={};beta=X(0);B=mul(beta,sub(C(1),beta))
  for ex,a in p.items():
   ey=list(ex);j=ey[1]//2;ey[1]%=2
   z=add(z,mul({tuple(ey):a},power(B,j)))
  return z
 def rbk(p):
  z={}
  for ex,a in p.items():
   ey=list(ex);j=min(ey[1],ey[2]//3);ey[1]-=j;ey[2]-=3*j;ey[0]+=2*j
   z[tuple(ey)]=z.get(tuple(ey),Q(0))+a
  return {ex:a for ex,a in z.items()if a}
 beta,t,E=[X(j)for j in range(3)];one=C(1)
 num=add(scale(mul(E,mul(sub(beta,t),add(beta,t))),2),
         scale(mul(E,add(beta,t)),-1),scale(mul(E,sub(beta,t)),-1))
 ck("Q_beta derivative with BOTH terms",sub(num,scale(mul(E,sub(mul(beta,sub(beta,one)),power(t,2))),2)))
 qnum=mul(E,sub(beta,t));qden=add(beta,t)
 ck("critical mixed field derivative cleared",add(
  scale(mul(power(beta,2),mul(E,mul(qnum,qden))),-4),
  mul(sub(power(beta,2),power(t,2)),power(add(mul(E,qden),qnum),2))))
 beta,b,k,x=[X(j)for j in range(4)];B2=power(beta,2);b2=power(b,2)
 ck("edge stationarity quadratic relation",rb2(add(B2,b2,scale(beta,-1))))
 ck("edge third derivative cleared",rb2(scale(mul(b,add(B2,b2,scale(beta,-1))),2)))
 ck("edge fourth derivative cleared",rb2(add(
  scale(mul(add(B2,b2),add(B2,scale(b2,3))),-2),
  scale(mul(B2,sub(C(3),scale(beta,2))),2))))
 ck("positive imaginary radial derivative",rb2(mul(k,add(B2,b2,scale(beta,-1)))))
 ck("radial Airy derivative cube",rbk(sub(B2,mul(b,power(k,3)))))
 # Work with Qz=12*b*k*z to avoid algebraic-root divisions.
 Qz=ga(gi(scale(mul(b,x),-12)),g(mul(sub(C(3),scale(beta,2)),power(x,2))))
 cube=gp(Qz,3)
 cube=tuple({ex:a for ex,a in part.items()if ex[3]<=4}for part in cube)
 lhs=gm(g(B2),cube)
 rhs=gm(g(power(scale(mul(b,k),12),3)),ga(
  gi(mul(b,power(x,3))),g(scale(mul(sub(C(3),scale(beta,2)),power(x,4)),Q(-1,4)))))
 jet=gs(lhs,rhs);ckG("right-decay normal-form cubic and quartic jets",tuple(rbk(p)for p in jet))
 Qwrong=ga(gi(scale(mul(b,x),-12)),g(scale(mul(sub(C(3),scale(beta,2)),power(x,2)),Q(3,2))))
 wrong=gp(Qwrong,3);wrong=tuple({ex:a for ex,a in part.items()if ex[3]<=4}for part in wrong)
 wrong=gs(gm(g(B2),wrong),rhs);need(any(rbk(p)for p in wrong),"Actual wrong1/8 right-bend control")
 x,Sh,Ch,beta=[X(j)for j in range(4)]
 E0=sub(mul(add(one,scale(power(x,2),2)),Sh),scale(mul(x,Ch),2))
 derivative=add(diff(E0,0),scale(mul(Ch,diff(E0,1)),2),scale(mul(Sh,diff(E0,2)),2))
 ck("E0 derivative has strict-positive integral",sub(derivative,scale(mul(power(x,2),Ch),4)))
 Tnum=sub(mul(beta,sub(scale(mul(x,Ch),2),Sh)),scale(mul(power(x,2),Sh),2))
 Trhs=add(scale(mul(beta,E0),-1),scale(mul(sub(one,beta),mul(power(x,2),Sh)),-2))
 ck("edge unique-global-max T derivative",sub(Tnum,Trhs))
 y=X(0);I=gi(one);P=g(add(scale(power(y,3),Q(1,3)),y))
 ckG("upper cubic critical fiber factorization",
  gs(gs(P,gi(C(Q(2,3)))),gc(gm(gp(gs(g(y),I),2),ga(g(y),gi(C(2)))),Q(1,3))))
 ckG("lower cubic critical fiber factorization",
  gs(ga(P,gi(C(Q(2,3)))),gc(gm(gp(ga(g(y),I),2),gs(g(y),gi(C(2)))),Q(1,3))))
 iv=gi(y);Piv=ga(gc(gp(iv,3),Q(1,3)),iv)
 ckG("shared-left-sheet boundary parameterization",gs(Piv,gi(sub(y,scale(power(y,3),Q(1,3))))))
 rho=X(0)
 ck("shared low midpoint phase",add(scale(power(scale(rho,-1),3),Q(1,3)),mul(power(rho,2),scale(rho,-1)),scale(power(rho,3),Q(4,3))))
 sigma,c=[X(j)for j in range(2)]
 a2=mul(c,power(sub(one,sigma),2))
 Lsigma=add(scale(a2,2),mul(sigma,sub(one,sigma)))
 ck("actual phase central linear cancellation",add(sub(scale(a2,2),Lsigma),mul(sigma,sub(one,sigma))))
 a,d,S,w,Z=[X(j)for j in range(5)]
 crit=add(mul(sub(a,S),add(w,Z)),mul(d,Z))
 ck("actual support-fixed stationary envelope",add(sub(mul(d,w),mul(sub(add(a,d),S),add(w,Z))),crit))
 a,d,m,j=[X(i)for i in range(4)]
 ck("old Gamma upper factor bound identity",sub(add(mul(a,m),one,mul(d,sub(m,one)),j),
  sub(mul(add(a,d),m),sub(sub(d,one),j))))
 m,j=X(0),X(1)
 ck("exact Gaussian-unitary-Hermite coefficient exponents",sub(power(sub(m,scale(j,2)),2),
  sub(power(m,2),scale(mul(j,sub(m,j)),4))))
 d,m,L,K,rho=[X(i)for i in range(5)]
 left=add(mul(d,mul(power(K,2),rho)),mul(d,K))
 right=add(mul(d,mul(m,mul(K,power(rho,4)))),mul(d,mul(m,power(rho,3))))
 ck("growing-fold amplitude-error scaling under K=m*rho^3",
  sub(sub(left,right),mul(d,mul(add(mul(K,rho),one),sub(K,mul(m,power(rho,3)))))))
 return out
def exact_controls(inp):
 cfg=inp["controls"];intervals=[(Q(a),Q(b))for a,b in cfg["intervals"]]
 need(all(Q(1,4)<a<=b for a,b in intervals),"Actual supercritical ratio interval fixtures")
 def rise(a,n):
  z=1
  for j in range(n):z*=a+j
  return z
 def qa(z,w):return(z[0]+w[0],z[1]+w[1])
 def qneg(z):return(-z[0],-z[1])
 def qm(z,w):return(z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0])
 def qs(z,c):return(z[0]*c,z[1]*c)
 def qi(z):
  norm=z[0]**2+z[1]**2;need(norm!=0,"Nonzero Gaussian denominator")
  return(z[0]/norm,-z[1]/norm)
 def qd(z,w):return qm(z,qi(w))
 def qp(z,n):
  w=(Q(1),Q(0))
  for _ in range(n):w=qm(w,z)
  return w
 moments=own=factors=edges=series=budgets=0
 for d in range(2,cfg["maximum_gap"]+1):
  for a in range(1,cfg["maximum_lower"]+1):
   c=Q(a,(d-1)**2)
   if not any(lo<=c<=hi for lo,hi in intervals):continue
   need(4*a>(d-1)**2,"ACTUAL strict abovecritical support")
   C=a+Q(d-1,2);kap=Q(d+1,d-1);w0=C**d*kap
   for m in range(1,cfg["maximum_order"]+1):
    A,B=Q(2,3),-Q(3,5);x=A/B
    direct=sum(Q(comb(m,j)*factorial(a*m+d*j))*A**(m-j)*B**j for j in range(m+1))
    poly=sum(Q(comb(m,j))*rise(a*m+1,d*j)*x**(m-j)for j in range(m+1))
    need(direct==B**m*factorial(a*m)*poly,"Actual factorial coefficient identity")
    moments+=1
    for j in range(d):
     f=a*m+1+d*(m-1)+j
     need(f==(a+d)*m-(d-1-j)and 0<f<=(a+d)*m,"Actual OLD upper factor")
     factors+=1
   for N in [2,3,5,8]:
    n=N-1;qN=Q(2,3)/(N**d*w0);qn=Q(2,3)/(n**d*w0)
    need(qn==Q(N,n)**d*qN,"TWO own-order arguments with support fixed")
    own+=1
 for den in range(2,cfg["pythagorean_parameter_max_denominator"]+1):
  for num in range(1,den):
   r=Q(num,den);cos=(1-r*r)/(1+r*r);sin=2*r/(1+r*r)
   beta=cos*cos;b=cos*sin
   need(0<beta<1 and b>0 and b*b==beta*(1-beta),"Exact positive compact-edge parameter")
   q=qd((beta,b),(beta,-b));norm=q[0]**2+q[1]**2
   need(norm==1,"Model edge target has unit modulus")
   one=(Q(1),Q(0));h=qd(qa(one,qneg(q)),qa(one,q))
   need(h==(0,-b/beta),"Exact stationary model logarithmic derivative")
   f2=qa((-1/beta+1,0),qneg(qp(h,2)));need(f2==(0,0),"Edge second derivative")
   f3=qs(qm(h,qa(one,qneg(qp(h,2)))),-2)
   need(f3==(0,2*b/(beta*beta)),"Edge third derivative")
   f4=qs(qm(qa(one,qneg(qp(h,2))),qa(one,qs(qp(h,2),-3))),-2)
   need(f4==(-2*(3-2*beta)/(beta*beta),0),"Edge fourth derivative")
   qFtq=qd(qs(q,-2),qp(qa(one,q),2))
   need(qFtq==(-1/(2*beta),0),"Positive imaginary radial parameter after Vz=i*k")
   need(1/(8*beta*b)>0,"Radial derivative cube has strict negative imaginary sign")
   need(qa(one,q)[0]**2+qa(one,q)[1]**2==4*beta,"Critical factor NONZERO despite nearby lower phase zero")
   edges+=1
 for j in range(1,cfg["positive_series_terms"]+1):
  a=Q(2**(2*j+1),factorial(2*j+1))
  need(a>0 and Q(2**(2*j+1),factorial(2*j))-a==2*j*a>0,"Positive hyperbolic strictness series")
  series+=2
 for t in range(2,41):
  N=t**3
  for d in [2,3,7,19]:
   need(Q(1,t**4)/Q(d,N)==Q(1,d*t),"Derivative-free bounded-field error divided by radial scale")
   need(Q(1,(t+1)*d)<Q(1,t*d),"The ratio tends monotonically to zero, not a claimed cutoff")
   budgets+=1
 # Genuinely false mathematical shortcuts.
 a,d,m=4,3,3
 co=[comb(m,j)*rise(a*m+1,d*j)for j in range(m+1)]
 need(co[2]**3!=co[1]**3*co[3],"ACTUAL Gamma cubic is not an exact model common-modulus polynomial")
 need(Q(4,(5-1)**2)==Q(1,4)and Q(4,5**2)<Q(1,4),"Wrong d-squared ratio misclassifies an integer cusp")
 need(2*rise(2*4+1,3)!=2*rise(4+1,3),"Actual adjacent Gamma shapes cannot be shared")
 return {"actual_moment_controls":moments,"actual_own_order_controls":own,
 "old_gamma_upper_factors":factors,"exact_rational_edge_jet_and_radial_controls":edges,
 "positive_hyperbolic_series_controls":series,"bounded_field_radial_error_controls":budgets,
 "genuine_math_negative_controls":4}
def negatives(ledger,repo,support):
 out=[]
 def reject(name,fn):
  try:fn()
  except(ValueError,AssertionError,KeyError,TypeError,OSError):out.append(name);return
  raise ValueError("Negative accepted:"+name)
 for role in SPECS:
  d=copy.deepcopy(ledger);d["sources"]=[r for r in d["sources"]if r["role"]!=role]
  reject("omit:"+role,lambda d=d:check_sources(d,repo,support))
  d=copy.deepcopy(ledger);next(r for r in d["sources"]if r["role"]==role)["sha256"]="0"*64
  reject("badpin:"+role,lambda d=d:check_sources(d,repo,support))
 for key in ["sharp_lower_annulus_premise","model_circlelaw_to_Gamma_transfer","large_thermal_theorem_premise",
  "all_order_SFC","global_bad_triple_set_finite","effective_cutoff","uniform_to_cusp_or_infinite_ratio","heavy_finite_boundary_computation_rerun"]:
  d=copy.deepcopy(SCOPE);d[key]=True;reject("false-scope:"+key,lambda d=d:check_inputs(d))
 for key in ["large_gap_core_independent_407","Kabluchko_2_3_logical_inverse_existence",
  "actual_zero_free_log_strip_and_complex_order_partials","both_actual_endpoint_capture_before_step_or_label_assumptions",
  "right_single_Ai_class_NOT_Scorer","proper_cubic_shared_left_sheet_low_midpoint"]:
  d=copy.deepcopy(SCOPE);d[key]=False;reject("missing-scope:"+key,lambda d=d:check_inputs(d))
 for key,value in [("407_role","logical input of large-gap core"),("parameters","all c>=1/4"),
  ("large_thermal_role","the entire large-thermal theorem is assumed")]:
  d=copy.deepcopy(SCOPE);d[key]=value;reject("wrong-role:"+key,lambda d=d:check_inputs(d))
 d=copy.deepcopy(SCOPE);d["schema"]=True;reject("typed-schema",lambda:check_inputs(d))
 d=copy.deepcopy(ledger);d["sources"]=[];reject("empty-manifest",lambda:check_sources(d,repo,support))
 return out
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--repo-root",type=Path,required=True)
 ap.add_argument("--support-root",type=Path,default=HERE);ap.add_argument("--output",type=Path,default=HERE/RN);ns=ap.parse_args()
 repo,support=ns.repo_root.resolve(),ns.support_root.resolve()
 inp=load(HERE/IN);check_inputs(inp);ledger=load(HERE/LN);rows=check_sources(ledger,repo,support)
 upstream=parent(repo,support,rows)
 M=import_module("supercritical_Q_poly_kernel",repo/SPECS["exact_polynomial_kernel"][1])
 identities=universal(M);counts=exact_controls(inp);bad=negatives(ledger,repo,support)
 check_sources(ledger,repo,support)
 result={"schema":1,"status":"PASS","scope":inp,"proof_sha256":H(HERE/PN),"checker_sha256":H(HERE/CN),
 "universal_ledger_sha256":H(HERE/IN),"source_manifest_sha256":H(HERE/LN),
 "mandatory_source_roles":len(rows),"sources":[rows[r]for r in sorted(rows)],
 "fresh_bounded_gap_407_fallback":upstream,"universal_polynomial_identities":identities,"counts":counts,
 "source_scope_negative_controls":bad,"assurance_flags":ASSURANCE,
 "meaning":"The complete written actual-cycle proof supplies uniformity. Exact algebra and current-source/portable checks support the jets, contours, scales and scoped imports; they are not formal verification or an effective or all-order certificate."}
 ns.output.parent.mkdir(parents=True,exist_ok=True);ns.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
 print(json.dumps({"status":"PASS","mandatory_sources":len(rows),"universal_identities":len(identities),
 "counts":counts,"source_scope_controls":len(bad),"output":str(ns.output)}))
if __name__=="__main__":main()

