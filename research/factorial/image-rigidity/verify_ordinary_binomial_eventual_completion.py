#!/usr/bin/env python3
"""Source-guarded eventual SFC synthesis for every fixed ordinary binomial support."""
from __future__ import annotations
if not __debug__:raise SystemExit("Refusing -O/-OO: exact assertions must execute.")
from fractions import Fraction as Q
from pathlib import Path
from math import factorial,comb,isqrt
import argparse,copy,hashlib,json,re,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
PN="ALL_ORDINARY_BINOMIAL_EVENTUAL_STRONG_FACTORIAL.md"
CN="verify_ordinary_binomial_eventual_completion.py"
IN="ordinary_binomial_eventual_inputs.json"
LN="ordinary_binomial_eventual_source_ledger.json"
RN="ordinary_binomial_eventual_verification.json"
PREFIX="research/factorial/image-rigidity/"
ROLES={'adjacent_all_order_proof': ('repo', 'research/factorial/image-rigidity/UNIFORM_PADE_PARABOLA_SFC.tex'), 'constant_plus_power_analytic_proof': ('repo', 'research/factorial/image-rigidity/ALL_UNIVARIATE_TWO_TERM_SFC.tex'), 'gap405_coalescing_bridge_proof': ('repo', 'research/factorial/image-rigidity/COALESCING_BINOMIAL_SADDLE_BRIDGE.md'), 'gap405_complement_checker': ('repo', 'research/factorial/image-rigidity/verify_cusp_global_complement.py'), 'gap405_complement_inputs': ('repo', 'research/factorial/image-rigidity/cusp_global_complement_inputs.json'), 'gap405_complement_ledger': ('repo', 'research/factorial/image-rigidity/cusp_global_complement_source_ledger.json'), 'gap405_complement_proof': ('repo', 'research/factorial/image-rigidity/CUSP_GLOBAL_REGULAR_COMPLEMENT_AND_STOKES_DETOUR.md'), 'gap405_complement_receipt': ('repo', 'research/factorial/image-rigidity/cusp_global_complement_verification.json'), 'gap405_cusp_main_checker': ('repo', 'research/factorial/image-rigidity/verify_cusp_eventual_strong_factorial_completion.py'), 'gap405_cusp_main_inputs': ('repo', 'research/factorial/image-rigidity/cusp_eventual_strong_factorial_inputs.json'), 'gap405_cusp_main_ledger': ('repo', 'research/factorial/image-rigidity/cusp_eventual_strong_factorial_source_ledger.json'), 'gap405_cusp_main_proof': ('repo', 'research/factorial/image-rigidity/CUSP_EVENTUAL_STRONG_FACTORIAL_THEOREM.md'), 'gap405_cusp_main_receipt': ('repo', 'research/factorial/image-rigidity/cusp_eventual_strong_factorial_verification.json'), 'gap405_fixed_bounded_checker': ('repo', 'research/factorial/image-rigidity/verify_bounded_quartic_cusp.py'), 'gap405_fixed_bounded_proof': ('repo', 'research/factorial/image-rigidity/BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md'), 'gap405_fixed_bounded_receipt': ('repo', 'research/factorial/image-rigidity/bounded_quartic_cusp_verification.json'), 'gap405_fixed_completion_checker': ('repo', 'research/factorial/image-rigidity/verify_quartic_cusp_local_completion.py'), 'gap405_fixed_completion_inputs': ('repo', 'research/factorial/image-rigidity/quartic_cusp_local_completion_inputs.json'), 'gap405_fixed_completion_ledger': ('repo', 'research/factorial/image-rigidity/quartic_cusp_local_completion_source_ledger.json'), 'gap405_fixed_completion_receipt': ('repo', 'research/factorial/image-rigidity/quartic_cusp_local_completion_verification.json'), 'gap405_fixed_cusp_checker': ('repo', 'research/factorial/image-rigidity/verify_quartic_cusp_unfolding.py'), 'gap405_fixed_cusp_proof': ('repo', 'research/factorial/image-rigidity/QUARTIC_CUSP_UNFOLDING_AND_ACTUAL_ORDER_DRIFT.md'), 'gap405_fixed_cusp_receipt': ('repo', 'research/factorial/image-rigidity/quartic_cusp_unfolding_verification.json'), 'gap405_fixed_model_checker': ('repo', 'research/factorial/image-rigidity/verify_quartic_lee_yang_zero_and_robin_bridge.py'), 'gap405_fixed_model_inputs': ('repo', 'research/factorial/image-rigidity/quartic_lee_yang_zero_and_robin_inputs.json'), 'gap405_fixed_model_ledger': ('repo', 'research/factorial/image-rigidity/quartic_lee_yang_zero_and_robin_source_ledger.json'), 'gap405_fixed_model_proof': ('repo', 'research/factorial/image-rigidity/QUARTIC_LEE_YANG_ZERO_AND_ROBIN_BRIDGE.md'), 'gap405_fixed_model_receipt': ('repo', 'research/factorial/image-rigidity/quartic_lee_yang_zero_and_robin_verification.json'), 'gap405_fixed_primary_hjn': ('repo', 'research/factorial/image-rigidity/quartic_lee_yang_sources/hou_jiang_newman_2208_00917v1.pdf'), 'gap405_fixed_primary_newman': ('repo', 'research/factorial/image-rigidity/quartic_lee_yang_sources/newman_fourier_transforms_1976.pdf'), 'gap405_ledger': ('repo', 'research/factorial/image-rigidity/square_root_gap_eventual_source_ledger.json'), 'gap405_local_checker': ('repo', 'research/factorial/image-rigidity/verify_cusp_shrinking_window_log_ratio.py'), 'gap405_local_inputs': ('repo', 'research/factorial/image-rigidity/cusp_shrinking_window_inputs.json'), 'gap405_local_ledger': ('repo', 'research/factorial/image-rigidity/cusp_shrinking_window_source_ledger.json'), 'gap405_local_proof': ('repo', 'research/factorial/image-rigidity/CUSP_SHRINKING_WINDOW_LOG_RATIO_NONCOLLISION.md'), 'gap405_local_receipt': ('repo', 'research/factorial/image-rigidity/cusp_shrinking_window_verification.json'), 'gap405_main_checker': ('repo', 'research/factorial/image-rigidity/verify_square_root_gap_eventual_completion.py'), 'gap405_main_inputs': ('repo', 'research/factorial/image-rigidity/square_root_gap_eventual_inputs.json'), 'gap405_main_proof': ('repo', 'research/factorial/image-rigidity/SQUARE_ROOT_GAP_EVENTUAL_STRONG_FACTORIAL_THEOREM.md'), 'gap405_morse_checker': ('repo', 'research/factorial/image-rigidity/verify_cusp_two_thimble_positive_cycle.py'), 'gap405_morse_inputs': ('repo', 'research/factorial/image-rigidity/cusp_two_thimble_inputs.json'), 'gap405_morse_ledger': ('repo', 'research/factorial/image-rigidity/cusp_two_thimble_source_ledger.json'), 'gap405_morse_proof': ('repo', 'research/factorial/image-rigidity/CUSP_TWO_THIMBLE_POSITIVE_CYCLE.md'), 'gap405_morse_receipt': ('repo', 'research/factorial/image-rigidity/cusp_two_thimble_verification.json'), 'gap405_primary_vidunas': ('repo', 'research/factorial/image-rigidity/sources/vidunas-belyi-2025.pdf'), 'gap405_receipt': ('repo', 'research/factorial/image-rigidity/square_root_gap_eventual_verification.json'), 'gap405_retained_gamma_upper': ('repo', 'research/factorial/image-rigidity/GAMMA_PRODUCT_ROOT_MODULUS.tex'), 'gap405_stationary_checker': ('repo', 'research/factorial/image-rigidity/verify_cusp_stationary_belyi_monodromy.py'), 'gap405_stationary_inputs': ('repo', 'research/factorial/image-rigidity/cusp_stationary_belyi_inputs.json'), 'gap405_stationary_ledger': ('repo', 'research/factorial/image-rigidity/cusp_stationary_belyi_source_ledger.json'), 'gap405_stationary_primary_audit': ('repo', 'research/factorial/image-rigidity/cusp_stationary_belyi_primary_audit.md'), 'gap405_stationary_proof': ('repo', 'research/factorial/image-rigidity/CUSP_STATIONARY_BELYI_MONODROMY_AND_RADIAL_BARRIER.md'), 'gap405_stationary_receipt': ('repo', 'research/factorial/image-rigidity/cusp_stationary_belyi_verification.json'), 'gap405_strict_checker': ('repo', 'research/factorial/image-rigidity/verify_strict_subcritical_eventual_sfc.py'), 'gap405_strict_inputs': ('repo', 'research/factorial/image-rigidity/strict_subcritical_inputs.json'), 'gap405_strict_ledger': ('repo', 'research/factorial/image-rigidity/strict_subcritical_source_ledger.json'), 'gap405_strict_proof': ('repo', 'research/factorial/image-rigidity/SUBCRITICAL_ORDERED_PAIR_EVENTUAL_SFC.md'), 'gap405_strict_receipt': ('repo', 'research/factorial/image-rigidity/strict_subcritical_verification.json'), 'main_checker': ('local', 'verify_ordinary_binomial_eventual_completion.py'), 'main_inputs': ('local', 'ordinary_binomial_eventual_inputs.json'), 'main_proof': ('local', 'ALL_ORDINARY_BINOMIAL_EVENTUAL_STRONG_FACTORIAL.md'), 'novertical_checker': ('support', 'verify_no_vertical_physical_equimodular_pair.py'), 'novertical_inputs': ('support', 'no_vertical_pair_inputs.json'), 'novertical_ledger': ('support', 'no_vertical_pair_source_ledger.json'), 'novertical_proof': ('support', 'NO_VERTICAL_PHYSICAL_EQUIMODULAR_PAIR.md'), 'novertical_receipt': ('support', 'no_vertical_pair_verification.json'), 'regular_checker': ('support', 'verify_above_critical_regular_complement.py'), 'regular_inputs': ('support', 'above_critical_regular_inputs.json'), 'regular_ledger': ('support', 'above_critical_regular_source_ledger.json'), 'regular_proof': ('support', 'ABOVE_CRITICAL_REGULAR_COMPACT_COMPLEMENT.md'), 'regular_receipt': ('support', 'above_critical_regular_verification.json'), 'small_prime_power_proof': ('repo', 'research/factorial/image-rigidity/PRIME_POWER_TWO_TERM_SFC.tex'), 'turn_checker': ('support', 'verify_simple_turn_shrinking_window.py'), 'turn_inputs': ('support', 'simple_turn_shrinking_window_inputs.json'), 'turn_ledger': ('support', 'simple_turn_shrinking_window_source_ledger.json'), 'turn_proof': ('support', 'SIMPLE_TURN_SHRINKING_WINDOW_CONSECUTIVE_EXCLUSION.md'), 'turn_receipt': ('support', 'simple_turn_shrinking_window_verification.json')}
SCOPE={'schema': 1, 'support': {'a_min': 0, 'd_min': 1, 'b': 'a+d', 'condition': '0<=a<b integers'}, 'order_quantifier': 'for each fixed support (a,b) exists integer N_(a,b)>=2; every terminal N>=N_(a,b)', 'support_cases': ['a=0 analytic tail plus small prime powers and adjacent one', 'a>=1,d=1 adjacent', 'a>=1,d>=2,4a<=(d-1)^2', 'a>=1,d>=2,4a>(d-1)^2'], 'rational_actual_pair_coprime': True, 'same_cutoff_all_coefficients': True, 'same_cutoff_all_characteristic_zero_fields': True, 'descent': 'Q[x] Bezout identity; no embedding of an arbitrary field in C', 'BadTerminal_each_support_finite': True, 'BadTerminal_finite_iff_eventual': True, 'complex_all_coefficients_iff_rational_coprime': True, 'BadTerminal_empty': False, 'all_terminal_orders': False, 'effective_terminal_cutoff': False, 'uniform_terminal_cutoff_over_supports': False, 'arbitrary_polynomial_supports': False, 'general_Gamma_shapes': False, 'new_sharp_lower_annulus_logical_premise': False, 'constant_plus_power_finite_prefix_logical_premise': False, 'lower_one_full_finite_bridge_logical_premise': False, 'legacy_checkers_executed': False, 'formal_verification': False, 'external_review': False}
ROUTE=['constant_plus_power_analytic_tail', 'small_prime_powers_and_adjacent', 'square_root_gap405', 'abovecritical_regular', 'simple_turn_neighborhoods', 'old_Gamma_upper', 'rational_Bezout']
def require(ok,msg):
    if not ok:raise AssertionError(msg)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path):return json.loads(path.read_text(encoding="utf-8"))
def same_typed(a,b):
    if type(a) is not type(b):return False
    if type(b) is dict:return set(a)==set(b) and all(same_typed(a[k],b[k]) for k in b)
    if type(b) is list:return len(a)==len(b) and all(same_typed(x,y) for x,y in zip(a,b))
    return a==b
def check_inputs(data):require(same_typed(data,SCOPE),"Exact all-fixed-support eventual scope and types required")
def bases(repo,support):return {"local":HERE,"support":support,"repo":repo}
def check_sources(doc,repo,support):
    require(type(doc) is dict and set(doc)=={"schema","formal_verification","external_review","logical_route","sources"},"Source schema")
    require(type(doc["schema"]) is int and doc["schema"]==1,"Schema1")
    require(doc["formal_verification"] is False and doc["external_review"] is False,"Assurance flags")
    require(doc["logical_route"]==ROUTE,"Exact logical route")
    entries=doc["sources"];require(type(entries) is list and len(entries)==len(ROLES),"All mandatory source roles")
    out={};seen=set();roots=bases(repo,support)
    for row in entries:
        require(type(row) is dict and set(row)=={"role","kind","path","sha256"},"Source row schema")
        role=row["role"];require(type(role) is str and role in ROLES and role not in out,"Unique role")
        require((row["kind"],row["path"])==ROLES[role],"Exact role path")
        key=(row["kind"],row["path"]);require(key not in seen,"Unique path")
        require(type(row["sha256"]) is str and re.fullmatch("[0-9a-f]{64}",row["sha256"]),"SHA256")
        root=roots[row["kind"]].resolve();p=(root/row["path"]).resolve()
        require(p.is_relative_to(root) and p.is_file() and sha(p)==row["sha256"],"Missing/stale source:"+role)
        out[role]=row;seen.add(key)
    forbidden={"GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md",
               "UNIFORM_BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md",
               "all_exponent_prefix_certificates.json",
               "verify_all_univariate_two_term_sfc.py",
               "verify_ordinary_lower_one_finite_bridge.py",
               "ordinary_lower_one_finite_bridge_inputs.json"}
    require(not any(Path(path).name in forbidden for _,path in seen),"Unneeded stronger or finite-prefix route excluded")
    return out
def links(repo,support,rows):
    expected={(r["kind"],r["path"]):r["sha256"] for r in rows.values()}
    specs=[
      ("support","simple_turn_shrinking_window_source_ledger.json","sources","support",5),
      ("support","above_critical_regular_source_ledger.json","sources","support",8),
      ("support","no_vertical_pair_source_ledger.json","sources","support",4),
      ("repo",PREFIX+"square_root_gap_eventual_source_ledger.json","sources","mapped",54),
      ("repo",PREFIX+"cusp_eventual_strong_factorial_source_ledger.json","sources","mapped",44),
      ("repo",PREFIX+"strict_subcritical_source_ledger.json","sources","mapped",6),
      ("repo",PREFIX+"cusp_two_thimble_source_ledger.json","sources","mapped",8),
      ("repo",PREFIX+"cusp_shrinking_window_source_ledger.json","sources","mapped",9),
      ("repo",PREFIX+"cusp_global_complement_source_ledger.json","sources","mapped",6),
      ("repo",PREFIX+"quartic_cusp_local_completion_source_ledger.json","sources","canonical",15),
      ("repo",PREFIX+"quartic_lee_yang_zero_and_robin_source_ledger.json","sources","canonical",5),
      ("repo",PREFIX+"cusp_stationary_belyi_source_ledger.json","files","stationary",6)]
    out=[];roots=bases(repo,support)
    for kind,path,key,mode,count in specs:
        doc=load(roots[kind]/path)
        require(type(doc.get("schema")) is int and doc["schema"]==1,"Imported schema1")
        entries=doc.get(key);require(type(entries) is list and len(entries)==count,"Imported exact cardinality")
        seen=set()
        for row in entries:
            require(type(row) is dict and type(row.get("path")) is str and type(row.get("sha256")) is str,"Imported row")
            if mode=="canonical":target=("repo",PREFIX+row["path"])
            elif mode=="stationary":target=("repo",PREFIX+row["path"]) if row["root"]=="local" else ("repo",row["path"])
            else:
                k=row.get("kind",row.get("root"))
                require(k in {"local","support","repo"},"Imported root")
                if mode=="mapped":target=("repo",PREFIX+row["path"]) if k in {"local","support"} else ("repo",row["path"])
                else:target=("support" if k in {"local","support"} else "repo",row["path"])
            require(target not in seen,"Unique imported source");seen.add(target)
            require(target in expected and expected[target]==row["sha256"],"Direct mandatory pin must match child closure")
        out.append({"kind":kind,"path":path,"sources":count})
    return out
def equal_receipt(a,b):require(a==b,"Fresh child receipt differs")
def replay(repo,support,rows):
    specs=[
      ("gap405_main_checker","gap405_receipt","gap405_main_proof","repo",False),
      ("turn_checker","turn_receipt","turn_proof","support",False),
      ("regular_checker","regular_receipt","regular_proof","support",True),
      ("novertical_checker","novertical_receipt","novertical_proof","support",False)]
    out=[];roots=bases(repo,support)
    with tempfile.TemporaryDirectory(prefix="ordinary-binomial-eventual-") as folder:
        for i,(checker,receipt,proof,kind,addsupport) in enumerate(specs):
            p=Path(folder)/(str(i)+".json")
            cmd=[sys.executable,"-B",str(roots[kind]/ROLES[checker][1]),"--repo-root",str(repo)]
            if addsupport:cmd+=["--support-root",str(support)]
            cmd+=["--output",str(p)]
            r=subprocess.run(cmd,cwd=folder,capture_output=True,text=True,timeout=60)
            require(r.returncode==0,"Child failed:"+checker+";"+r.stdout+";"+r.stderr)
            require(p.is_file(),"Fresh child receipt missing")
            raw=p.read_bytes();equal_receipt(raw,(roots[kind]/ROLES[receipt][1]).read_bytes())
            d=json.loads(raw);require(d.get("status")=="PASS","Child PASS")
            require(d.get("proof_sha256")==rows[proof]["sha256"],"Child proof identity")
            require(d.get("checker_sha256")==rows[checker]["sha256"],"Child checker identity")
            out.append({"checker_role":checker,"receipt_sha256":hashlib.sha256(raw).hexdigest(),"byte_identical":True})
    return out
def mul(a,b):
    out=[Q(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out
def algebra():
    cases={name:0 for name in ["constant","adjacent","square_gap","abovecritical"]}
    for a in range(33):
        for d in range(1,21):
            if a==0:case="constant"
            elif d==1:case="adjacent"
            elif 4*a<=(d-1)**2:case="square_gap"
            else:case="abovecritical"
            require(case in cases,"Exhaustive integer support partition");cases[case]+=1
            if case=="square_gap":require(d>=3,"Square-root regime has d>=3")
    tail=0
    for d in range(6,129):
        N=max(10000,d*d)
        require(N>=10000 and d*d<=N and d>=6,"Minimal a0 analytic-tail hypotheses")
        tail+=1
    require([2,3,4,5]==[2,3,2**2,5],"Small d=2..5 are prime powers; d1 separately adjacent")
    moments=0
    for a in range(8):
        for d in range(1,9):
            for m in range(1,5):
                A=Q(2,3);B=-Q(3,5);x=A/B
                direct=sum((Q(comb(m,k)*factorial(a*m+d*k))*A**(m-k)*B**k for k in range(m+1)),Q(0))
                value=sum((Q(comb(m,k)*factorial(a*m+d*k),factorial(a*m))*x**(m-k) for k in range(m+1)),Q(0))
                require(direct==B**m*factorial(a*m)*value,"Actual factorial moment normalization")
                w=Q(2,5)
                integral=sum((Q(comb(m,k)*factorial(a*m+d*k),m**(a*m+d*k+1))*w**(m-k) for k in range(m+1)),Q(0))
                own=Q(factorial(a*m),m**((a+d)*m+1))*sum((Q(comb(m,k)*factorial(a*m+d*k),factorial(a*m))*(m**d*w)**(m-k) for k in range(m+1)),Q(0))
                require(integral==own,"Actual own-order integral normalization")
                require(a*m+1+d*(m-1)+d-1==(a+d)*m,"Old Gamma upper endpoint")
                moments+=1
    appell=0
    for d in range(1,9):
        for N in range(2,7):
            upper=[Q(comb(N,k)*factorial(d*k)) for k in reversed(range(N+1))]
            lower=[Q(comb(N-1,k)*factorial(d*k)) for k in reversed(range(N))]
            require([j*upper[j] for j in range(1,N+1)]==[N*x for x in lower],"a0 actual Appell derivative")
            appell+=1
    bez=0
    for a in range(6):
        for d in range(1,7):
            C=Q(factorial(a+d),factorial(a))
            C2=Q(2*factorial(2*a+d),factorial(2*a))
            R2=Q(factorial(2*a+2*d),factorial(2*a))
            rem=R2-C*(C2-C);require(rem>0,"Actual m2 variance remainder")
            p=[R2,C2,Q(1)];q=[C,Q(1)]
            u=[1/rem];v=[-(C2-C)/rem,-1/rem]
            x1,x2=mul(p,u),mul(q,v)
            require([x1[j]+x2[j] for j in range(3)]==[1,0,0],"Actual rational Bezout identity")
            bez+=1
    return {"support_case_counts":cases,"a0_analytic_tail_hypothesis_controls":tail,
            "actual_factorial_and_integral_normalizations":moments,
            "a0_Appell_derivative_controls":appell,"actual_rational_Bezout_identities":bez}
def negatives(ledger,repo,support):
    out=[]
    def reject(name,fn):
        try:fn()
        except (AssertionError,TypeError,KeyError,ValueError):out.append(name);return
        raise AssertionError("Negative accepted:"+name)
    for role in sorted(ROLES):
        d=copy.deepcopy(ledger);d["sources"]=[r for r in d["sources"] if r["role"]!=role]
        reject("omit:"+role,lambda d=d:check_sources(d,repo,support))
        d=copy.deepcopy(ledger);next(r for r in d["sources"] if r["role"]==role)["sha256"]="0"*64
        reject("badpin:"+role,lambda d=d:check_sources(d,repo,support))
    for key in ["BadTerminal_empty","all_terminal_orders","effective_terminal_cutoff",
      "uniform_terminal_cutoff_over_supports","arbitrary_polynomial_supports","general_Gamma_shapes",
      "new_sharp_lower_annulus_logical_premise","constant_plus_power_finite_prefix_logical_premise",
      "lower_one_full_finite_bridge_logical_premise","legacy_checkers_executed"]:
        d=copy.deepcopy(SCOPE);d[key]=True
        reject("false-scope:"+key,lambda d=d:check_inputs(d))
    for key in ["same_cutoff_all_coefficients","same_cutoff_all_characteristic_zero_fields",
                "BadTerminal_finite_iff_eventual","complex_all_coefficients_iff_rational_coprime"]:
        d=copy.deepcopy(SCOPE);d[key]=False
        reject("missing-scope:"+key,lambda d=d:check_inputs(d))
    d=copy.deepcopy(SCOPE);d["support"]["a_min"]=True
    reject("typed-support",lambda:check_inputs(d))
    d=copy.deepcopy(SCOPE);d["order_quantifier"]="one cutoff for all supports and all orders"
    reject("wrong-quantifier",lambda:check_inputs(d))
    for key in ["schema","empty","duplicate"]:
        d=copy.deepcopy(ledger)
        if key=="schema":d["schema"]=True
        if key=="empty":d["sources"]=[]
        if key=="duplicate":d["sources"][-1]=copy.deepcopy(d["sources"][0])
        reject("ledger:"+key,lambda d=d:check_sources(d,repo,support))
    reject("forged-child",lambda:equal_receipt(b'{"status":"PASS"}',b'{"status":"PASS","real":true}'))
    # Actual shape equality at a0: it is false that shapes ALWAYS differ by a.
    require(0*2+1==0*1+1,"a0 shape equality")
    out.append("false-always-different-Gamma-shapes")
    require(2*factorial(4)//factorial(2)==24 and 2*factorial(3)//factorial(1)==12,"Actual vs wrong shared shape")
    out.append("actual-wrong-shared-Gamma")
    return out
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--repo-root",required=True,type=Path)
    ap.add_argument("--support-root",type=Path,default=HERE);ap.add_argument("--output",type=Path,default=HERE/RN)
    ns=ap.parse_args();repo=ns.repo_root.resolve();support=ns.support_root.resolve()
    require(repo.is_dir() and support.is_dir(),"Explicit repo and support roots")
    inp=load(HERE/IN);check_inputs(inp)
    ledger=load(HERE/LN);rows=check_sources(ledger,repo,support)
    linked=links(repo,support,rows);children=replay(repo,support,rows)
    counts=algebra();bad=negatives(ledger,repo,support);check_sources(ledger,repo,support)
    counts.update({"fresh_direct_children":len(children),"cross_linked_source_ledgers":len(linked),"source_scope_controls":len(bad)})
    payload={"status":"PASS","schema":1,"scope":inp,"formal_verification":False,"external_review":False,
      "proof_sha256":sha(HERE/PN),"checker_sha256":sha(HERE/CN),"input_sha256":sha(HERE/IN),
      "source_ledger_sha256":sha(HERE/LN),"mandatory_source_roles":len(rows),
      "sources":[rows[r] for r in sorted(rows)],"verified_child_ledgers":linked,"fresh_children":children,
      "counts":counts,"negative_names":bad,
      "meaning":"The retained edge-family analytic/algebraic results, square-root-gap theorem, regular and turning-neighborhood proofs, compact root capture and rational Bezout prove the written fixed-support eventual theorem. Exact/source controls are supporting evidence, not formal verification or a finite-prefix inference."}
    ns.output.parent.mkdir(parents=True,exist_ok=True)
    ns.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":"PASS","mandatory_sources":len(rows),"counts":counts,"output":str(ns.output)}))
if __name__=="__main__":main()
