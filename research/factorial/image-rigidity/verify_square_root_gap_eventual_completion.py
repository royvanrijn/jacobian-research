#!/usr/bin/env python3
"""Source-guarded eventual SFC synthesis for square-root separated supports."""
from __future__ import annotations
if not __debug__:raise SystemExit("Refusing -O/-OO: exact assertions must execute.")
from fractions import Fraction as Q
from pathlib import Path
from math import factorial,comb,isqrt
import argparse,copy,hashlib,json,re,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
PN="SQUARE_ROOT_GAP_EVENTUAL_STRONG_FACTORIAL_THEOREM.md"
CN="verify_square_root_gap_eventual_completion.py"
IN="square_root_gap_eventual_inputs.json"
LN="square_root_gap_eventual_source_ledger.json"
RN="square_root_gap_eventual_verification.json"
PREFIX="research/factorial/image-rigidity/"
ROLES={'coalescing_bridge_proof': ('repo', 'research/factorial/image-rigidity/COALESCING_BINOMIAL_SADDLE_BRIDGE.md'), 'complement_checker': ('support', 'verify_cusp_global_complement.py'), 'complement_inputs': ('support', 'cusp_global_complement_inputs.json'), 'complement_ledger': ('support', 'cusp_global_complement_source_ledger.json'), 'complement_proof': ('support', 'CUSP_GLOBAL_REGULAR_COMPLEMENT_AND_STOKES_DETOUR.md'), 'complement_receipt': ('support', 'cusp_global_complement_verification.json'), 'cusp_main_checker': ('support', 'verify_cusp_eventual_strong_factorial_completion.py'), 'cusp_main_inputs': ('support', 'cusp_eventual_strong_factorial_inputs.json'), 'cusp_main_ledger': ('support', 'cusp_eventual_strong_factorial_source_ledger.json'), 'cusp_main_proof': ('support', 'CUSP_EVENTUAL_STRONG_FACTORIAL_THEOREM.md'), 'cusp_main_receipt': ('support', 'cusp_eventual_strong_factorial_verification.json'), 'fixed_bounded_checker': ('repo', 'research/factorial/image-rigidity/verify_bounded_quartic_cusp.py'), 'fixed_bounded_proof': ('repo', 'research/factorial/image-rigidity/BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md'), 'fixed_bounded_receipt': ('repo', 'research/factorial/image-rigidity/bounded_quartic_cusp_verification.json'), 'fixed_completion_checker': ('repo', 'research/factorial/image-rigidity/verify_quartic_cusp_local_completion.py'), 'fixed_completion_inputs': ('repo', 'research/factorial/image-rigidity/quartic_cusp_local_completion_inputs.json'), 'fixed_completion_ledger': ('repo', 'research/factorial/image-rigidity/quartic_cusp_local_completion_source_ledger.json'), 'fixed_completion_receipt': ('repo', 'research/factorial/image-rigidity/quartic_cusp_local_completion_verification.json'), 'fixed_cusp_checker': ('repo', 'research/factorial/image-rigidity/verify_quartic_cusp_unfolding.py'), 'fixed_cusp_proof': ('repo', 'research/factorial/image-rigidity/QUARTIC_CUSP_UNFOLDING_AND_ACTUAL_ORDER_DRIFT.md'), 'fixed_cusp_receipt': ('repo', 'research/factorial/image-rigidity/quartic_cusp_unfolding_verification.json'), 'fixed_model_checker': ('repo', 'research/factorial/image-rigidity/verify_quartic_lee_yang_zero_and_robin_bridge.py'), 'fixed_model_inputs': ('repo', 'research/factorial/image-rigidity/quartic_lee_yang_zero_and_robin_inputs.json'), 'fixed_model_ledger': ('repo', 'research/factorial/image-rigidity/quartic_lee_yang_zero_and_robin_source_ledger.json'), 'fixed_model_proof': ('repo', 'research/factorial/image-rigidity/QUARTIC_LEE_YANG_ZERO_AND_ROBIN_BRIDGE.md'), 'fixed_model_receipt': ('repo', 'research/factorial/image-rigidity/quartic_lee_yang_zero_and_robin_verification.json'), 'fixed_primary_hjn': ('repo', 'research/factorial/image-rigidity/quartic_lee_yang_sources/hou_jiang_newman_2208_00917v1.pdf'), 'fixed_primary_newman': ('repo', 'research/factorial/image-rigidity/quartic_lee_yang_sources/newman_fourier_transforms_1976.pdf'), 'local_checker': ('support', 'verify_cusp_shrinking_window_log_ratio.py'), 'local_inputs': ('support', 'cusp_shrinking_window_inputs.json'), 'local_ledger': ('support', 'cusp_shrinking_window_source_ledger.json'), 'local_proof': ('support', 'CUSP_SHRINKING_WINDOW_LOG_RATIO_NONCOLLISION.md'), 'local_receipt': ('support', 'cusp_shrinking_window_verification.json'), 'main_checker': ('local', 'verify_square_root_gap_eventual_completion.py'), 'main_inputs': ('local', 'square_root_gap_eventual_inputs.json'), 'main_proof': ('local', 'SQUARE_ROOT_GAP_EVENTUAL_STRONG_FACTORIAL_THEOREM.md'), 'morse_checker': ('support', 'verify_cusp_two_thimble_positive_cycle.py'), 'morse_inputs': ('support', 'cusp_two_thimble_inputs.json'), 'morse_ledger': ('support', 'cusp_two_thimble_source_ledger.json'), 'morse_proof': ('support', 'CUSP_TWO_THIMBLE_POSITIVE_CYCLE.md'), 'morse_receipt': ('support', 'cusp_two_thimble_verification.json'), 'primary_vidunas': ('repo', 'research/factorial/image-rigidity/sources/vidunas-belyi-2025.pdf'), 'retained_gamma_upper': ('repo', 'research/factorial/image-rigidity/GAMMA_PRODUCT_ROOT_MODULUS.tex'), 'stationary_checker': ('repo', 'research/factorial/image-rigidity/verify_cusp_stationary_belyi_monodromy.py'), 'stationary_inputs': ('repo', 'research/factorial/image-rigidity/cusp_stationary_belyi_inputs.json'), 'stationary_ledger': ('repo', 'research/factorial/image-rigidity/cusp_stationary_belyi_source_ledger.json'), 'stationary_primary_audit': ('repo', 'research/factorial/image-rigidity/cusp_stationary_belyi_primary_audit.md'), 'stationary_proof': ('repo', 'research/factorial/image-rigidity/CUSP_STATIONARY_BELYI_MONODROMY_AND_RADIAL_BARRIER.md'), 'stationary_receipt': ('repo', 'research/factorial/image-rigidity/cusp_stationary_belyi_verification.json'), 'strict_checker': ('support', 'verify_strict_subcritical_eventual_sfc.py'), 'strict_inputs': ('support', 'strict_subcritical_inputs.json'), 'strict_ledger': ('support', 'strict_subcritical_source_ledger.json'), 'strict_proof': ('support', 'SUBCRITICAL_ORDERED_PAIR_EVENTUAL_SFC.md'), 'strict_receipt': ('support', 'strict_subcritical_verification.json')}
SCOPE={'schema': 1, 'support': {'a_min': 1, 'd_min': 3, 'b': 'a+d', 'condition': '4*a<=(d-1)^2', 'equivalent_condition': 'sqrt(b)-sqrt(a)>=1'}, 'order_quantifier': 'for each fixed allowed support (a,b) exists integer N_(a,b)>=2; every terminal N>=N_(a,b)', 'strict_and_integer_cusp_cases_exhaustive': True, 'rational_actual_pair_coprime': True, 'same_cutoff_all_coefficients': True, 'same_cutoff_all_characteristic_zero_fields': True, 'descent': 'Q[x] Bezout identity; no embedding of an arbitrary field in C', 'BadTerminal_each_allowed_support_finite': True, 'BadTerminal_finite_iff_eventual': True, 'BadTerminal_empty': False, 'all_terminal_orders': False, 'effective_terminal_cutoff': False, 'uniform_terminal_cutoff_over_supports': False, 'remaining_supports_proved': False, 'new_sharp_lower_annulus_logical_premise': False, 'uniform_H_bound_logical_premise': False, 'HJN_monotonicity_logical_premise': False, 'formal_verification': False, 'external_review': False}
ROUTE=["strict_subcritical","independent_integer_cusp","strict_or_equality_exhaustion","retained_upper","rational_Bezout"]
def require(ok,msg):
    if not ok:raise AssertionError(msg)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path):return json.loads(path.read_text(encoding="utf-8"))
def same_typed(a,b):
    if type(a) is not type(b):return False
    if type(b) is dict:return set(a)==set(b) and all(same_typed(a[k],b[k]) for k in b)
    if type(b) is list:return len(a)==len(b) and all(same_typed(x,y) for x,y in zip(a,b))
    return a==b
def check_inputs(data):require(same_typed(data,SCOPE),"Exact square-root-gap eventual scope and types required")
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
               "UNIFORM_BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md"}
    require(not any(Path(path).name in forbidden for _,path in seen),"Stronger unused theorem excluded")
    return out
def links(repo,support,rows):
    expected={(r["kind"],r["path"]):r["sha256"] for r in rows.values()}
    specs=[
      ("support","cusp_eventual_strong_factorial_source_ledger.json","sources","support",44),
      ("support","strict_subcritical_source_ledger.json","sources","support",6),
      ("support","cusp_two_thimble_source_ledger.json","sources","support",8),
      ("support","cusp_shrinking_window_source_ledger.json","sources","support",9),
      ("support","cusp_global_complement_source_ledger.json","sources","support",6),
      ("repo",PREFIX+"quartic_cusp_local_completion_source_ledger.json","sources","canonical",15),
      ("repo",PREFIX+"quartic_lee_yang_zero_and_robin_source_ledger.json","sources","canonical",5),
      ("repo",PREFIX+"cusp_stationary_belyi_source_ledger.json","files","stationary",6)]
    result=[];roots=bases(repo,support)
    for kind,path,key,mode,count in specs:
        data=load(roots[kind]/path)
        require(type(data.get("schema")) is int and data["schema"]==1,"Imported schema1")
        entries=data.get(key);require(type(entries) is list and len(entries)==count,"Exact imported cardinality")
        seen=set()
        for row in entries:
            require(type(row) is dict and type(row.get("path")) is str and type(row.get("sha256")) is str,"Imported row")
            if mode=="canonical":target=("repo",PREFIX+row["path"])
            elif mode=="stationary":target=("repo",PREFIX+row["path"]) if row["root"]=="local" else ("repo",row["path"])
            else:
                k=row.get("kind",row.get("root"))
                require(k in {"local","support","repo"},"Imported root")
                target=("support" if k in {"local","support"} else "repo",row["path"])
            require(target not in seen,"Unique imported source");seen.add(target)
            require(target in expected and expected[target]==row["sha256"],"Cross-linked source must match direct mandatory pin")
        result.append({"kind":kind,"path":path,"sources":count})
    return result
def equal_receipt(a,b):require(a==b,"Fresh child receipt differs")
def replay(repo,support,rows):
    result=[]
    specs=[("cusp_main_checker","cusp_main_receipt","cusp_main_proof"),
           ("strict_checker","strict_receipt","strict_proof")]
    with tempfile.TemporaryDirectory(prefix="sqrt-gap-eventual-") as folder:
        for index,(checker,receipt,proof) in enumerate(specs):
            out=Path(folder)/(str(index)+".json")
            command=[sys.executable,"-B",str(support/ROLES[checker][1]),
                     "--repo-root",str(repo),"--support-root",str(support),"--output",str(out)]
            process=subprocess.run(command,cwd=folder,capture_output=True,text=True,timeout=60)
            require(process.returncode==0,"Fresh child failed:"+checker+";"+process.stdout+";"+process.stderr)
            require(out.is_file(),"Fresh receipt missing")
            raw=out.read_bytes();equal_receipt(raw,(support/ROLES[receipt][1]).read_bytes())
            data=json.loads(raw)
            require(data.get("status")=="PASS","Child PASS")
            require(data.get("proof_sha256")==rows[proof]["sha256"],"Child proof identity")
            require(data.get("checker_sha256")==rows[checker]["sha256"],"Child checker identity")
            result.append({"checker_role":checker,"receipt_sha256":hashlib.sha256(raw).hexdigest(),"byte_identical":True})
    return result
def mul(a,b):
    out=[Q(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out
def algebra():
    support=0;equal=0;strict=0;outside=0
    # Two independent exact tests of sqrt(a+d)>=sqrt(a)+1:
    # classify the sign of d-1 first, then isolate sqrt(a) in rational intervals.
    scale=10**8
    for a in range(1,33):
        for d in range(1,33):
            discriminant=(d-1)**2-4*a
            first=(d>=3 and discriminant>=0)
            lower=Q(isqrt(a*scale*scale),scale)
            upper=lower if lower*lower==a else lower+Q(1,scale)
            if d-1>=2*upper:second=True
            elif d-1<2*lower:second=False
            else:
                require(isqrt(a)**2==a and d-1==2*isqrt(a),"Exact radical boundary")
                second=True
            require(first==second,"Exact square-root support equivalence")
            b=a+d
            require((a+b-1)**2-4*a*b==discriminant,"Universal support discriminant specialization")
            if first and discriminant==0:
                s=isqrt(a);require(s*s==a and d==2*s+1 and b==(s+1)**2,"Equality iff adjacent squares")
                equal+=1
            elif first:strict+=1
            else:outside+=1
            support+=1
    count=0
    for a in range(1,7):
        for d in range(3,9):
            for m in range(1,5):
                A=Q(2,3);B=-Q(3,5);x=A/B
                left=sum((Q(comb(m,k)*factorial(a*m+d*k))*A**(m-k)*B**k for k in range(m+1)),Q(0))
                value=sum((Q(comb(m,k)*factorial(a*m+d*k),factorial(a*m))*x**(m-k) for k in range(m+1)),Q(0))
                require(left==B**m*factorial(a*m)*value,"Actual factorial normalization")
                require(a*m+1+d*(m-1)+d-1==(a+d)*m,"Own-order upper-bound endpoint")
                count+=1
    # A rational polynomial Bezout identity transports under the unique Q->K map.
    x1=mul([Q(1),Q(0),Q(1)],[Q(1,2)])
    x2=mul([Q(1),Q(1)],[Q(1,2),-Q(1,2)])
    require([x1[j]+x2[j] for j in range(3)]==[1,0,0],"Rational Bezout identity")
    return {"exact_support_equivalence_controls":support,"adjacent_square_boundaries":equal,
            "strict_support_cases":strict,"excluded_support_cases":outside,
            "actual_factorial_normalizations":count,"rational_Bezout_identities":1}
def negatives(ledger,repo,support):
    names=[]
    def reject(name,fn):
        try:fn()
        except (AssertionError,TypeError,KeyError,ValueError):names.append(name);return
        raise AssertionError("Negative accepted:"+name)
    for role in sorted(ROLES):
        data=copy.deepcopy(ledger);data["sources"]=[r for r in data["sources"] if r["role"]!=role]
        reject("omit:"+role,lambda data=data:check_sources(data,repo,support))
        data=copy.deepcopy(ledger);next(r for r in data["sources"] if r["role"]==role)["sha256"]="0"*64
        reject("badpin:"+role,lambda data=data:check_sources(data,repo,support))
    for key in ["BadTerminal_empty","all_terminal_orders","effective_terminal_cutoff",
                "uniform_terminal_cutoff_over_supports","remaining_supports_proved",
                "new_sharp_lower_annulus_logical_premise","uniform_H_bound_logical_premise"]:
        data=copy.deepcopy(SCOPE);data[key]=True
        reject("false-scope:"+key,lambda data=data:check_inputs(data))
    for key in ["same_cutoff_all_coefficients","same_cutoff_all_characteristic_zero_fields",
                "BadTerminal_finite_iff_eventual","strict_and_integer_cusp_cases_exhaustive"]:
        data=copy.deepcopy(SCOPE);data[key]=False
        reject("missing-scope:"+key,lambda data=data:check_inputs(data))
    data=copy.deepcopy(SCOPE);data["support"]["condition"]="4*a>(d-1)^2"
    reject("wrong-support",lambda:check_inputs(data))
    data=copy.deepcopy(SCOPE);data["order_quantifier"]="one N0 for all supports and all terminal orders"
    reject("wrong-quantifier",lambda:check_inputs(data))
    data=copy.deepcopy(SCOPE);data["support"]["a_min"]=True
    reject("typed-support",lambda:check_inputs(data))
    for key in ["schema","empty","duplicate"]:
        data=copy.deepcopy(ledger)
        if key=="schema":data["schema"]=True
        if key=="empty":data["sources"]=[]
        if key=="duplicate":data["sources"][-1]=copy.deepcopy(data["sources"][0])
        reject("ledger:"+key,lambda data=data:check_sources(data,repo,support))
    reject("forged-child",lambda:equal_receipt(b'{"status":"PASS"}',b'{"status":"PASS","real":true}'))
    # Wrong support and fixed-shape substitutions are genuinely false controls.
    require((2-1)**2<4*1,"Gap two at a=1 is outside this theorem")
    names.append("actual-excluded-support")
    # The own-order m=2, a=1,d=3 first coefficient is 2*(3)_3=120,
    # while the false shared shape gamma=2 gives 2*(2)_3=48.
    require(2*factorial(5)//factorial(2)==120 and 2*factorial(4)//factorial(1)==48,"Actual/shared-shape coefficient")
    names.append("actual-wrong-shared-Gamma")
    return names
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
      "meaning":"The strict and equality analytic theorems plus exact support exhaustion and rational Bezout prove the written theorem. Finite/source controls authenticate supporting evidence; they are not formal verification or a finite-prefix inference."}
    ns.output.parent.mkdir(parents=True,exist_ok=True)
    ns.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":"PASS","mandatory_sources":len(rows),"counts":counts,"output":str(ns.output)}))
if __name__=="__main__":main()
