#!/usr/bin/env python3
"""Exact/source support for actual large-thermal matching and its fixed-strip consequence."""
from __future__ import annotations
if not __debug__:
    raise SystemExit("Optimized Python refused: all exact/source controls must execute.")
from pathlib import Path
from fractions import Fraction as Q
from math import comb,factorial
import argparse,copy,hashlib,importlib.util,json,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
PN='LARGE_THERMAL_QUARTIC_FOLD_MATCHING.md'
CN='verify_large_thermal_quartic_fold_matching.py'
IN='large_thermal_quartic_fold_inputs.json'
LN='large_thermal_quartic_fold_source_ledger.json'
RN='large_thermal_quartic_fold_verification.json'
SPECS={'main_proof': ('local', 'LARGE_THERMAL_QUARTIC_FOLD_MATCHING.md', 'written_actual_complete_cycle_argument'), 'exact_checker': ('local', 'verify_large_thermal_quartic_fold_matching.py', 'exact_and_source_support'), 'exact_inputs': ('local', 'large_thermal_quartic_fold_inputs.json', 'exact_typed_fixture_and_scope'), 'primary_context': ('local', 'LARGE_THERMAL_PRIMARY_CONTEXT.md', 'primary_scope_audit'), 'dlmf_cubic_framework': ('local', 'sources/dlmf-36-12-large-thermal.html', 'classical_framework_context_only'), 'thermal_parent_proof': ('repo', 'research/factorial/image-rigidity/UNIFORM_BOUNDED_THERMAL_TRANSITION_EVENTUAL_SFC.md', 'small_T_and_global_capture_mechanisms_and_separate_strip_consequence'), 'thermal_parent_checker': ('repo', 'research/factorial/image-rigidity/verify_bounded_thermal_transition.py', 'fresh_parent_exact_and_source_support'), 'thermal_parent_inputs': ('repo', 'research/factorial/image-rigidity/bounded_thermal_transition_inputs.json', 'exact_parent_input_boundary'), 'thermal_parent_ledger': ('repo', 'research/factorial/image-rigidity/bounded_thermal_transition_source_ledger.json', 'complete_parent_source_boundary'), 'thermal_parent_receipt': ('repo', 'research/factorial/image-rigidity/bounded_thermal_transition_verification.json', 'byte_identical_fresh_parent_receipt'), 'exact_rational_kernel': ('repo', 'research/factorial/image-rigidity/verify_uniform_cusp_eventual.py', 'exact_QQ_rational_function_kernel_only'), 'actual_simple_turn_mechanism': ('repo', 'research/factorial/image-rigidity/SIMPLE_TURN_SHRINKING_WINDOW_CONSECUTIVE_EXCLUSION.md', 'analytic_local_normal_form_and_proper_arms_revalidated_here'), 'dlmf_airy_ode': ('repo', 'research/factorial/image-rigidity/sources/dlmf-9-2-endpoint.html', 'classical_model_ODE_and_uniqueness_only'), 'dlmf_airy_sector': ('repo', 'research/factorial/image-rigidity/sources/dlmf-9-7-endpoint.html', 'classical_model_Ai_sector_leading_coefficient_only'), 'dlmf_airy_zeros': ('repo', 'research/factorial/image-rigidity/sources/dlmf-9-9-endpoint.html', 'classical_model_Ai_negative_real_zero_location_only'), 'dlmf_cfu_framework': ('repo', 'research/factorial/image-rigidity/sources/dlmf-2-4-endpoint.html', 'classical_framework_context_allows_Ai_or_Scorer'), 'primary_kabluchko': ('repo', 'research/factorial/image-rigidity/uniform_cusp_sources/kabluchko_2203_05533v3.pdf', 'parent_classical_CW_context_only_not_actual_fold'), 'primary_newman': ('repo', 'research/factorial/image-rigidity/quartic_lee_yang_sources/newman_fourier_transforms_1976.pdf', 'model_location_for_separate_bounded_parent_only'), 'compact_402_tails': ('repo', 'research/factorial/image-rigidity/UNIFORM_BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md', 'actual_compact_phase_and_complete_tails_revalidated'), 'global_annulus_404': ('repo', 'research/factorial/image-rigidity/GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md', 'uniform_actual_root_annulus_for_global_capture'), 'uniform_regular_capture': ('repo', 'research/factorial/image-rigidity/UNIFORM_CUSP_CURIE_WEISS_REGULAR_COMPLEMENT.md', 'actual_regular_model_geometry_revalidated_for_r_tending_zero'), 'uniform_small_T_packets': ('repo', 'research/factorial/image-rigidity/UNIFORM_GROWING_QUARTIC_CUSP_NEIGHBORHOOD.md', 'small_temperature_quartic_packet_mechanism_revalidated'), 'old_Gamma_catch': ('repo', 'research/factorial/image-rigidity/GAMMA_PRODUCT_ROOT_MODULUS.tex', 'old_soft_lower_endpoint_local_catch_only')}
SCOPE={'schema': 1, 'integer_parameters': 'a>=1,d>=2,N>=2', 'sequential_theorem': 'no actual common-root sequence with N->infinity, 0<abs(r)->0 and sqrt(N)*abs(r)->infinity', 'quantified_large_thermal_region': 'exists unevaluated delta0>0,U0,N0 excluding 0<abs(r)<delta0 and sqrt(N)*abs(r)>U0 for N>=N0', 'nearcritical_strip_consequence': 'one unevaluated fixed relative strip and uniform eventual cutoff, using separate bounded-thermal parent at one U0', 'all_characteristic_zero_fields_and_coefficients': True, 'nonzero_detuning_lattice': 'abs(e)>=1/4', 'actual_critical_powers_retained': True, 'upper_high_saddle_automatically_active': False, 'actual_lower_fold_cycle': 'one nonzero right-exit Ai cycle', 'canonical_zeta_omega_at_limit': '-i', 'shared_midpoint': 'negative-real inverse sheet, canonical loss 4*rho^3/3', 'original_horizontal_midpoint_supplies_cubic_loss': False, 'endpoint_capture_precedes_small_adjacent_step': True, 'complex_m_uses_consistent_actual_logs': True, 'bounded_fold_uses_same_two_term_coefficients': True, 'remainder_differentiated_without_Cauchy': False, 'uses_407_or_411': False, 'support_global_uniform_cutoff': False, 'all_order_SFC': False, 'effective_width_or_cutoff': False, 'finite_checks_certify_uniform_estimates': False, 'controls': {'support_s_max': 9, 'detunings': [-1, 1], 'moment_max': 4, 'cone_grid': 12, 'comparison_minimum': 8, 'comparison_maximum': 40}, 'assurance_flags': {'external_review': False, 'formal_verification': False, 'independent_replay': False}}
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
def check_inputs(data):need(typed(data,SCOPE),"Exact typed input/scope snapshot")
def source_path(role,repo):
    kind,path,_=SPECS[role];base=HERE if kind=="local"else repo
    p=base/path;need(p.resolve().is_relative_to(base.resolve()),"Source escapes declared base")
    return p
def check_sources(ledger,repo):
    need(type(ledger)is dict and set(ledger)=={"schema","sources","assurance_flags"},"Exact ledger shape")
    need(type(ledger["schema"])is int and ledger["schema"]==1,"Exact integer ledger schema1")
    need(typed(ledger["assurance_flags"],ASSURANCE),"No assurance upgrade")
    rows=ledger["sources"]
    need(type(rows)is list and len(rows)==len(SPECS),"Every exact mandatory role")
    indexed={}
    for row in rows:
        need(type(row)is dict and set(row)=={"role","kind","path","relation","sha256"},"Exact source record shape")
        role=row["role"];need(type(role)is str and role in SPECS and role not in indexed,"Unique mandatory source role")
        kind,path,relation=SPECS[role]
        need((row["kind"],row["path"],row["relation"])==(kind,path,relation),"Exact kind/path/relation")
        need(type(row["sha256"])is str and len(row["sha256"])==64 and all(c in "0123456789abcdef"for c in row["sha256"]),"Exact SHA256 syntax")
        p=source_path(role,repo);need(p.is_file()and sha(p)==row["sha256"],"Missing/stale source:"+role)
        indexed[role]=row
    need(set(indexed)==set(SPECS),"All exact roles")
    return indexed
def import_kernel(repo):
    sys.dont_write_bytecode=True
    p=source_path("exact_rational_kernel",repo)
    spec=importlib.util.spec_from_file_location("large_thermal_QQ_kernel",p)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    need(module.VARS==("t","q","z","h","r","s","l","p","H","b","lam","y","x","d","m","rho","k","delta"),"Exact rational kernel variables")
    return module
def verify_parent_bytes(actual,expected,ledger_hash):
    need(actual==expected,"Fresh parent receipt must match byte for byte")
    d=json.loads(actual)
    need(type(d)is dict and d.get("status")=="PASS"and d.get("schema")==1,"Exact parent PASS")
    need(d.get("source_ledger_sha256")==ledger_hash,"Fresh parent ledger identity")
    s=d["scope"]
    need(s.get("all_order_SFC")is False and s.get("effective_cutoff")is False and s.get("unbounded_thermal_matching_solved")is False,"Parent remains bounded-thermal only")
    return d
def fresh_parent(repo,rows):
    checker=source_path("thermal_parent_checker",repo)
    with tempfile.TemporaryDirectory(prefix="large-thermal-parent-")as tmp:
        out=Path(tmp)/"fresh.json"
        result=subprocess.run([sys.executable,"-B",str(checker),"--repo-root",str(repo),"--support-root",str(checker.parent),"--output",str(out)],cwd=tmp,capture_output=True,text=True,timeout=90)
        need(result.returncode==0,"Fresh bounded parent failed:"+result.stdout+";"+result.stderr)
        actual=out.read_bytes();expected=source_path("thermal_parent_receipt",repo).read_bytes()
        d=verify_parent_bytes(actual,expected,rows["thermal_parent_ledger"]["sha256"])
    parent=load(source_path("thermal_parent_ledger",repo))
    need(type(parent["schema"])is int and parent["schema"]==1 and len(parent["sources"])==27,"Exact complete parent ledger")
    seen=set();crosslinks=0
    for row in parent["sources"]:
        role=row["role"];need(role not in seen,"Unique parent role");seen.add(role)
        path=(checker.parent if row["kind"]in("local","support")else repo)/row["path"]
        need(path.is_file()and sha(path)==row["sha256"],"Full parent source boundary:"+role)
        for ours in rows.values():
            if source_path(ours["role"],repo).resolve()==path.resolve():
                need(ours["sha256"]==row["sha256"],"Direct/parent source crosslink");crosslinks+=1
    need(d["mandatory_source_roles"]==27 and len(d["universal_identities"])==19,"Actual parent exact scope")
    return {"mandatory_roles":27,"source_boundary_verified":len(seen),"direct_crosslinks":crosslinks,
            "fresh_receipt_byte_identical":True,"receipt_sha256":hashlib.sha256(actual).hexdigest()}

def universal(module):
    R=module.RF;v=module.var;out=[]
    ii=module.VARS.index("delta")
    def reduce_i(poly):
        ans={}
        for e,c in poly.items():
            n=e[ii];a=list(e);a[ii]=n%2;a=tuple(a)
            ans[a]=ans.get(a,Q(0))+c*((-1)**(n//2))
        return {e:c for e,c in ans.items()if c}
    def c0(expr):
        need(bool(reduce_i(expr.d)),"Zero complex rational denominator")
        return not reduce_i(expr.n)
    def ck(name,expr):
        need(c0(expr),"Exact complex/rational identity:"+name);out.append(name)
    def parts(expr):
        den=reduce_i(expr.d);need(all(e[ii]==0 for e in den),"Real denominator needed")
        rr={};im={}
        for e,c in reduce_i(expr.n).items():
            ee=list(e);power=ee[ii];ee[ii]=0;ee=tuple(ee)
            dst=rr if power==0 else im;dst[ee]=dst.get(ee,Q(0))+c
        return R(rr,den),R(im,den)
    vv,x,y,k,z,w,sigma,wr,wi,I=(v(n)for n in("t","x","y","k","z","q","p","h","H","delta"))
    f=-vv**4/12+sigma*vv**2/2+w*vv
    real,imag=parts(f.subst(t=x+I*y,q=wr+I*wi))
    target=-x**4/12+(y*y+sigma)*x*x/2+wr*x
    ck("tilted horizontal x-dependent real phase",real-target-(-y**4/12-sigma*y*y/2-y*wi))
    ck("stationary cubic map",f.diff("t")+vv**3/3-sigma*vv-w)
    fold=f.subst(p=-1,q=-2*I/3)
    ck("actual lower-fold horizontal polynomial",fold.subst(t=x-I)+Q(1,4)-I*x**3/3+x**4/12)
    ck("lower-fold stationary factorization",fold.diff("t")+(vv+I)**2*(vv-2*I)/3)
    ck("lower stationary derivative",fold.diff("t").subst(t=-I))
    ck("lower second derivative vanishes",fold.diff("t").diff("t").subst(t=-I))
    ck("nonzero lower third derivative",fold.diff("t").diff("t").diff("t").subst(t=-I)-2*I)
    ck("inactive upper action is two",fold.subst(t=2*I)-2)
    ck("selected lower action is minus one quarter",fold.subst(t=-I)+Q(1,4))
    sr,si=parts(f.diff("t").subst(t=x-I*k,q=-I*z))
    ck("two-maxima real stationary equation",sr-x*(-x*x/3+k*k+sigma))
    ck("two-maxima imaginary stationary equation",si-(k*x*x-k**3/3-sigma*k-z))
    ck("two-maxima imaginary parameter relation",k*(3*(k*k+sigma))-k**3/3-sigma*k-z-(8*k**3/3+2*sigma*k-z))
    ck("two-maxima real horizontal curvature",-3*(k*k+sigma)+k*k+sigma+2*(k*k+sigma))
    canon=-I*x+x*x/12
    expr=canon**3/3-I*x**3/3+x**4/12
    ix=module.VARS.index("x")
    ck("right-exit cubic map through fourth degree",R({e:c for e,c in expr.n.items()if e[ix]<5},expr.d))
    linear=w*(-I+I*z)
    ck("canonical unfolding has zeta_omega minus i",linear.diff("z")-I*w)
    rho=v("rho");phi=z**3/3+rho*rho*z
    ck("common negative-real cubic midpoint gap",phi.subst(z=-rho)+4*rho**3/3)
    ck("upper cubic critical value",phi.subst(z=I*rho)-2*I*rho**3/3)
    ck("lower cubic critical value",phi.subst(z=-I*rho)+2*I*rho**3/3)
    tau,gamma,M,mm,L,thermalr,beta,gap,qq=(v(n)for n in("r","l","m","d","y","s","b","h","lam"))
    lhs=(beta*gap*qq/(mm*tau**3)).subst(b=-L*(1-thermalr)/gap)
    rhs=(-gamma*(1-thermalr)*qq/M).subst(l=L*tau,m=mm*tau**4)
    ck("exact own-order omega derivative",lhs-rhs)
    K=v("k")
    ck("first transported fold amplitude exponent",(gamma/(M*rho))/(gamma*rho)-rho/(M*rho**3))
    ck("second transported amplitude exponent",(1/(M*K*rho**3)-1/K**2).subst(m=K/rho**3))
    ck("complex-m Cauchy exponent",(1/(mm*K*gamma*rho)).subst(d=K/(tau**4*rho**3))-tau**4*rho**2/(gamma*K**2))
    n=v("s");eps=1/n**2
    weak=Q(3,4)*((1+eps)**2-1)**2*n**3
    ck("physical midpoint loss lacks cubic scale",weak-(3/n+3/n**3+Q(3,4)/n**5))
    ck("effective order differs by support-fixed tau fourth power",mm*tau**4-(mm-1)*tau**4-tau**4)
    # These are actual nonzero algebraic defects, not merely scope flags.
    need(not c0(linear.diff("z")+I*w),"Wrong unfolding orientation must fail")
    need(Q(2)>-Q(1,4),"The high inactive critical value is genuinely higher")
    need(not c0((fold.subst(t=2*I)-fold.subst(t=-I))),"All-saddle dominance is not a cycle rule")
    return out
def controls(inp):
    cfg=inp["controls"];jets=moments=maps=lattice=cones=labels=weak=0
    # Truncated raw-z=t/L Taylor coefficients use only rational arithmetic.
    def mul(a,b,n):
        c=[Q(0)]*(n+1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):
                if i+j<=n:c[i+j]+=x*y
        return c
    def logunit(a,n):
        need(a[0]==1,"Unit series")
        u=list(a);u[0]=0;p=[Q(1)]+[Q(0)]*n;out=[Q(0)]*(n+1)
        for j in range(1,n+1):
            p=mul(p,u,n)
            for i in range(n+1):out[i]+=Q((-1)**(j+1),j)*p[i]
        return out
    for ss in range(1,cfg["support_s_max"]+1):
        for detune in cfg["detunings"]:
            a=ss*ss+detune
            if a<1:continue
            d=2*ss+1;s=Q(ss);Sc=s*(s+1);C=Q(a)+s;e=Q(a)-s*s;r=e/C
            need(abs(e)>=Q(1,4)and C*(1-r)==Sc,"Ordinary nonzero detuning")
            need(4*e==4*a-(d-1)**2,"Exact quarter lattice")
            kap=(s+1)/s;n=4
            pp=[Q(comb(d,j))if j<=d else Q(0)for j in range(n+1)];pp[0]+=kap
            pp=[v/pp[0]for v in pp]
            lg=logunit(pp,n)
            phase=[Q(0)]+[Q(a*(-1)**(j+1),j)for j in range(1,n+1)];phase[1]-=C
            phase=[x+y for x,y in zip(phase,lg)]
            need(phase[1]==0 and phase[2]==-e/2 and phase[3]==e/3 and phase[4]==-Sc*Sc/12-e/4,"Actual detuned phase jets")
            jets+=1
            for m in range(1,cfg["moment_max"]+1):
                A=Q(2,3);B=-Q(3,5);x=A/B
                literal=sum(Q(comb(m,j)*factorial(a*m+d*j))*A**(m-j)*B**j for j in range(m+1))
                pol=sum(Q(comb(m,j)*factorial(a*m+d*j),factorial(a*m))*x**(m-j)for j in range(m+1))
                need(literal==B**m*factorial(a*m)*pol,"Literal factorial versus actual moment polynomial")
                moments+=1
            for N in(2,3,5,8):
                nn=N-1;qN=Q(2,3)/N**d;qn=Q(2,3)/nn**d;lam=Q(N,nn)**d
                need(qn==lam*qN,"Both own-order q maps")
                need(nn*r*r==Q(nn,N)*N*r*r,"Both effective order scalings")
                need(qn-1==lam*(qN-1)+(lam-1),"Exact affine field path")
                maps+=1
            lattice+=1
    eps=Q(1,32)
    for A in range(1,cfg["cone_grid"]+1):
        for j in range(-8,9):
            y=eps*A*Q(j,8)
            for h in range(0,9):
                T=eps*A*Q(h,4)
                # Both endpoints remain in their narrow negative-real cone.
                if abs(y+T)>eps*A:continue
                need(A*A+(y+T)**2>=A*A,"Increment cannot cancel the negative real component")
                need(T<=2*eps*A,"Both endpoint cones bound the increment")
                for k in range(5):
                    yi=y+Q(k,4)*T
                    need(abs(yi)<=eps*A,"Whole chord stays central")
                cones+=1
    for n in range(cfg["comparison_minimum"],cfg["comparison_maximum"]+1):
        # M=n^3, tau=1/n and gamma=1/2 make the exact remainder/main ratio rational.
        M=n**3;tau=Q(1,n);gamma=Q(1,2)
        ratio=(tau**4*Q(1,n**5)+Q(1,n**4))/(gamma/M)
        need(ratio==2*(Q(1,n)+Q(1,n**6)),"Derivative-free two-term remainder ratio")
        need(ratio<1,"Same-label remainder is below an actual radial step in the control family")
        labels+=1
        wk=3*Q(1,n)+3*Q(1,n**3)+Q(3,4*n**5)
        need(wk<4*Q(1,n),"Old physical midpoint has vanishing relative cubic gap")
        weak+=1
    return {"actual_detuned_jet_controls":jets,"literal_moment_controls":moments,
            "own_order_maps":maps,"nonzero_lattice_controls":lattice,
            "both_endpoint_cone_controls":cones,"bounded_two_term_comparison_controls":labels,
            "weak_midpoint_countercontrols":weak,"genuine_math_negative_controls":3}


def negatives(ledger,repo,parent_receipt):
    out=[]
    def reject(name,f):
        try:f()
        except (ValueError,KeyError,TypeError,FileNotFoundError):out.append(name);return
        raise ValueError("Negative unexpectedly accepted:"+name)
    for role in SPECS:
        d=copy.deepcopy(ledger);d["sources"]=[r for r in d["sources"]if r["role"]!=role]
        reject("omit:"+role,lambda d=d:check_sources(d,repo))
        d=copy.deepcopy(ledger);next(r for r in d["sources"]if r["role"]==role)["sha256"]="0"*64
        reject("badpin:"+role,lambda d=d:check_sources(d,repo))
    for key in("upper_high_saddle_automatically_active","original_horizontal_midpoint_supplies_cubic_loss",
               "remainder_differentiated_without_Cauchy","uses_407_or_411","support_global_uniform_cutoff",
               "all_order_SFC","effective_width_or_cutoff","finite_checks_certify_uniform_estimates"):
        d=copy.deepcopy(SCOPE);d[key]=True
        reject("false-scope:"+key,lambda d=d:check_inputs(d))
    for key in("actual_critical_powers_retained","endpoint_capture_precedes_small_adjacent_step",
               "complex_m_uses_consistent_actual_logs","bounded_fold_uses_same_two_term_coefficients",
               "all_characteristic_zero_fields_and_coefficients"):
        d=copy.deepcopy(SCOPE);d[key]=False
        reject("missing-scope:"+key,lambda d=d:check_inputs(d))
    d=copy.deepcopy(SCOPE);d["schema"]=True
    reject("boolean-input-schema",lambda:check_inputs(d))
    d=copy.deepcopy(SCOPE);d["controls"]["moment_max"]=True
    reject("boolean-integer-control",lambda:check_inputs(d))
    d=copy.deepcopy(ledger);d["schema"]=True
    reject("boolean-ledger-schema",lambda:check_sources(d,repo))
    d=copy.deepcopy(ledger);d["sources"]=[]
    reject("empty-ledger",lambda:check_sources(d,repo))
    d=copy.deepcopy(ledger);d["sources"].append(copy.deepcopy(d["sources"][0]))
    reject("duplicate-ledger",lambda:check_sources(d,repo))
    d=json.loads(parent_receipt);d["scope"]["unbounded_thermal_matching_solved"]=True
    bad=(json.dumps(d,indent=2,sort_keys=True)+"\n").encode()
    reject("resigned-false-parent-scope",lambda:verify_parent_bytes(bad,bad,d["source_ledger_sha256"]))
    return out
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--repo-root",required=True,type=Path)
    ap.add_argument("--output",type=Path,default=HERE/RN);ns=ap.parse_args()
    repo=ns.repo_root.resolve();inp=load(HERE/IN);check_inputs(inp)
    ledger=load(HERE/LN);rows=check_sources(ledger,repo)
    kernel=import_kernel(repo);parent=fresh_parent(repo,rows)
    identities=universal(kernel);counts=controls(inp)
    neg=negatives(ledger,repo,source_path("thermal_parent_receipt",repo).read_bytes())
    check_sources(ledger,repo)
    result={"schema":1,"status":"PASS","scope":inp,
            "proof_sha256":sha(HERE/PN),"checker_sha256":sha(HERE/CN),
            "input_sha256":sha(HERE/IN),"source_ledger_sha256":sha(HERE/LN),
            "mandatory_source_roles":len(rows),"sources":[rows[r]for r in sorted(rows)],
            "fresh_bounded_parent":parent,"universal_identities":identities,
            "counts":counts,"source_scope_negative_controls":neg,"assurance_flags":ASSURANCE,
            "meaning":"The written complete actual-cycle proof supplies the uniform estimates and exhaustion. Exact identities, finite controls and current-source replays support algebra and scope; they do not numerically certify analytic estimates, evaluate a cutoff or prove all-order SFC."}
    ns.output.parent.mkdir(parents=True,exist_ok=True)
    ns.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":"PASS","mandatory_sources":len(rows),"universal_identities":len(identities),
                     "fresh_parent_receipt_byte_identical":True,"counts":counts,
                     "source_scope_controls":len(neg),"output":str(ns.output)}))
if __name__=="__main__":main()

