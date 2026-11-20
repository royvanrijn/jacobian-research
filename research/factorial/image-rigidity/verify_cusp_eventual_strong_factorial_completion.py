#!/usr/bin/env python3
"""Source-guarded synthesis of eventual SFC at each fixed integer cusp."""
from __future__ import annotations
if not __debug__:raise SystemExit("Refusing -O/-OO: exact assertions must execute.")
from fractions import Fraction as Q
from pathlib import Path
from math import factorial,comb
import argparse,copy,hashlib,json,re,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
PN="CUSP_EVENTUAL_STRONG_FACTORIAL_THEOREM.md"
CN="verify_cusp_eventual_strong_factorial_completion.py"
IN="cusp_eventual_strong_factorial_inputs.json"
LN="cusp_eventual_strong_factorial_source_ledger.json"
RN="cusp_eventual_strong_factorial_verification.json"
PREFIX="research/factorial/image-rigidity/"
ROLES={'coalescing_bridge_proof': ('repo', 'research/factorial/image-rigidity/COALESCING_BINOMIAL_SADDLE_BRIDGE.md'), 'complement_checker': ('support', 'verify_cusp_global_complement.py'), 'complement_inputs': ('support', 'cusp_global_complement_inputs.json'), 'complement_ledger': ('support', 'cusp_global_complement_source_ledger.json'), 'complement_proof': ('support', 'CUSP_GLOBAL_REGULAR_COMPLEMENT_AND_STOKES_DETOUR.md'), 'complement_receipt': ('support', 'cusp_global_complement_verification.json'), 'fixed_bounded_checker': ('repo', 'research/factorial/image-rigidity/verify_bounded_quartic_cusp.py'), 'fixed_bounded_proof': ('repo', 'research/factorial/image-rigidity/BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md'), 'fixed_bounded_receipt': ('repo', 'research/factorial/image-rigidity/bounded_quartic_cusp_verification.json'), 'fixed_completion_checker': ('repo', 'research/factorial/image-rigidity/verify_quartic_cusp_local_completion.py'), 'fixed_completion_inputs': ('repo', 'research/factorial/image-rigidity/quartic_cusp_local_completion_inputs.json'), 'fixed_completion_ledger': ('repo', 'research/factorial/image-rigidity/quartic_cusp_local_completion_source_ledger.json'), 'fixed_completion_receipt': ('repo', 'research/factorial/image-rigidity/quartic_cusp_local_completion_verification.json'), 'fixed_cusp_checker': ('repo', 'research/factorial/image-rigidity/verify_quartic_cusp_unfolding.py'), 'fixed_cusp_proof': ('repo', 'research/factorial/image-rigidity/QUARTIC_CUSP_UNFOLDING_AND_ACTUAL_ORDER_DRIFT.md'), 'fixed_cusp_receipt': ('repo', 'research/factorial/image-rigidity/quartic_cusp_unfolding_verification.json'), 'fixed_model_checker': ('repo', 'research/factorial/image-rigidity/verify_quartic_lee_yang_zero_and_robin_bridge.py'), 'fixed_model_inputs': ('repo', 'research/factorial/image-rigidity/quartic_lee_yang_zero_and_robin_inputs.json'), 'fixed_model_ledger': ('repo', 'research/factorial/image-rigidity/quartic_lee_yang_zero_and_robin_source_ledger.json'), 'fixed_model_proof': ('repo', 'research/factorial/image-rigidity/QUARTIC_LEE_YANG_ZERO_AND_ROBIN_BRIDGE.md'), 'fixed_model_receipt': ('repo', 'research/factorial/image-rigidity/quartic_lee_yang_zero_and_robin_verification.json'), 'fixed_primary_hjn': ('repo', 'research/factorial/image-rigidity/quartic_lee_yang_sources/hou_jiang_newman_2208_00917v1.pdf'), 'fixed_primary_newman': ('repo', 'research/factorial/image-rigidity/quartic_lee_yang_sources/newman_fourier_transforms_1976.pdf'), 'local_checker': ('support', 'verify_cusp_shrinking_window_log_ratio.py'), 'local_inputs': ('support', 'cusp_shrinking_window_inputs.json'), 'local_ledger': ('support', 'cusp_shrinking_window_source_ledger.json'), 'local_proof': ('support', 'CUSP_SHRINKING_WINDOW_LOG_RATIO_NONCOLLISION.md'), 'local_receipt': ('support', 'cusp_shrinking_window_verification.json'), 'main_checker': ('local', 'verify_cusp_eventual_strong_factorial_completion.py'), 'main_inputs': ('local', 'cusp_eventual_strong_factorial_inputs.json'), 'main_proof': ('local', 'CUSP_EVENTUAL_STRONG_FACTORIAL_THEOREM.md'), 'morse_checker': ('support', 'verify_cusp_two_thimble_positive_cycle.py'), 'morse_inputs': ('support', 'cusp_two_thimble_inputs.json'), 'morse_ledger': ('support', 'cusp_two_thimble_source_ledger.json'), 'morse_proof': ('support', 'CUSP_TWO_THIMBLE_POSITIVE_CYCLE.md'), 'morse_receipt': ('support', 'cusp_two_thimble_verification.json'), 'primary_vidunas': ('repo', 'research/factorial/image-rigidity/sources/vidunas-belyi-2025.pdf'), 'retained_gamma_upper': ('repo', 'research/factorial/image-rigidity/GAMMA_PRODUCT_ROOT_MODULUS.tex'), 'stationary_checker': ('repo', 'research/factorial/image-rigidity/verify_cusp_stationary_belyi_monodromy.py'), 'stationary_inputs': ('repo', 'research/factorial/image-rigidity/cusp_stationary_belyi_inputs.json'), 'stationary_ledger': ('repo', 'research/factorial/image-rigidity/cusp_stationary_belyi_source_ledger.json'), 'stationary_primary_audit': ('repo', 'research/factorial/image-rigidity/cusp_stationary_belyi_primary_audit.md'), 'stationary_proof': ('repo', 'research/factorial/image-rigidity/CUSP_STATIONARY_BELYI_MONODROMY_AND_RADIAL_BARRIER.md'), 'stationary_receipt': ('repo', 'research/factorial/image-rigidity/cusp_stationary_belyi_verification.json')}
SCOPE={
 "schema":1,"fixed_support":{"s_min":1,"a":"s^2","d":"2s+1"},
 "order_quantifier":"for each s exists integer N_s>=2; all terminal N>=N_s",
 "rational_actual_pair_coprime":True,
 "same_cutoff_all_coefficients":True,
 "same_cutoff_all_characteristic_zero_fields":True,
 "descent":"Q[x] Bezout identity; no embedding of an arbitrary field in C",
 "BadTerminal_s_finite":True,"BadTerminal_finite_iff_eventual":True,
 "BadTerminal_s_empty":False,"all_terminal_orders":False,
 "effective_Ns":False,"uniform_Ns_in_s":False,
 "new_sharp_lower_annulus_logical_premise":False,
 "uniform_H_bound_logical_premise":False,
 "HJN_monotonicity_logical_premise":False,
 "formal_verification":False,"external_review":False,
}
def require(ok,msg):
    if not ok:raise AssertionError(msg)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path):return json.loads(path.read_text(encoding="utf-8"))
def check_inputs(data):
    require(type(data) is dict and data==SCOPE,"Exact fixed-support eventual scope required")
    for k,v in SCOPE.items():require(type(data[k]) is type(v),"Scope type:"+k)
    require(type(data["fixed_support"]["s_min"]) is int,"Integer s_min")
def roots(repo,support):return {"local":HERE,"support":support,"repo":repo}
def check_sources(doc,repo,support):
    require(type(doc) is dict and set(doc)=={"schema","formal_verification","external_review","logical_route","sources"},"Source schema")
    require(type(doc["schema"]) is int and doc["schema"]==1,"Source schema1")
    require(doc["formal_verification"] is False and doc["external_review"] is False,"Assurance flags")
    require(doc["logical_route"]==["fixed_bounded_cusp","shrinking_cusp","two_thimble","global_complement","retained_upper","rational_Bezout"],"Minimal logical route")
    rows=doc["sources"];require(type(rows) is list and len(rows)==len(ROLES),"All mandatory source roles")
    out={};paths=set();base=roots(repo,support)
    for row in rows:
        require(type(row) is dict and set(row)=={"role","kind","path","sha256"},"Source row schema")
        r=row["role"];require(type(r) is str and r in ROLES and r not in out,"Unique source role")
        require((row["kind"],row["path"])==ROLES[r],"Exact source path")
        key=(row["kind"],row["path"]);require(key not in paths,"Unique source path")
        require(type(row["sha256"]) is str and re.fullmatch("[0-9a-f]{64}",row["sha256"]),"SHA256")
        root=base[row["kind"]].resolve();p=(root/row["path"]).resolve()
        require(p.is_relative_to(root) and p.is_file() and sha(p)==row["sha256"],"Missing/stale:"+r)
        paths.add(key);out[r]=row
    forbidden={"GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md",
               "UNIFORM_BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md"}
    require(not any(Path(path).name in forbidden for _,path in paths),"Unneeded stronger theorem is outside this route")
    return out
def child_links(repo,support,rows):
    expected={(r["kind"],r["path"]):r["sha256"] for r in rows.values()}
    specs=[
      ("support","cusp_two_thimble_source_ledger.json","sources","support",8),
      ("support","cusp_shrinking_window_source_ledger.json","sources","support",9),
      ("support","cusp_global_complement_source_ledger.json","sources","support",6),
      ("repo",PREFIX+"quartic_cusp_local_completion_source_ledger.json","sources","canonical",15),
      ("repo",PREFIX+"quartic_lee_yang_zero_and_robin_source_ledger.json","sources","canonical",5),
      ("repo",PREFIX+"cusp_stationary_belyi_source_ledger.json","files","stationary",6),
    ]
    summary=[];base=roots(repo,support)
    for kind,path,key,mode,count in specs:
        d=load(base[kind]/path)
        require(type(d.get("schema")) is int and d["schema"]==1,"Imported ledger schema")
        entries=d.get(key);require(type(entries) is list and len(entries)==count,"Imported exact cardinality")
        seen=set()
        for entry in entries:
            require(type(entry) is dict,"Imported row object")
            if mode=="canonical":target=("repo",PREFIX+entry["path"])
            elif mode=="stationary":
                target=("repo",PREFIX+entry["path"]) if entry["root"]=="local" else ("repo",entry["path"])
            else:
                e_kind=entry.get("kind",entry.get("root"))
                require(e_kind in {"local","support","repo"},"Imported source root")
                target=("support" if e_kind in {"local","support"} else "repo",entry["path"])
            require(target not in seen,"Imported unique source");seen.add(target)
            require(target in expected and expected[target]==entry["sha256"],"Imported hash agrees with direct mandatory pin")
        summary.append({"kind":kind,"path":path,"sources":count})
    return summary
def equal_receipt(a,b):require(a==b,"Fresh child receipt differs")
def replay(repo,support,rows):
    children=[
      ("morse_checker","morse_receipt","morse_proof",False),
      ("local_checker","local_receipt","local_proof",False),
      ("complement_checker","complement_receipt","complement_proof",True),
    ]
    result=[]
    with tempfile.TemporaryDirectory(prefix="cusp-eventual-main-") as folder:
        for index,(script_r,receipt_r,proof_r,add_support) in enumerate(children):
            out=Path(folder)/(str(index)+".json")
            cmd=[sys.executable,"-B",str(support/ROLES[script_r][1]),"--repo-root",str(repo)]
            if add_support:cmd+=["--support-root",str(support)]
            cmd+=["--output",str(out)]
            r=subprocess.run(cmd,cwd=folder,capture_output=True,text=True,timeout=30)
            require(r.returncode==0,"Child failed:"+script_r+";"+r.stdout+";"+r.stderr)
            require(out.is_file(),"Child receipt missing")
            raw=out.read_bytes();equal_receipt(raw,(support/ROLES[receipt_r][1]).read_bytes())
            d=json.loads(raw);require(d.get("status")=="PASS","Child PASS")
            require(d.get("proof_sha256")==rows[proof_r]["sha256"],"Child proof identity")
            require(d.get("checker_sha256")==rows[script_r]["sha256"],"Child checker identity")
            result.append({"checker_role":script_r,"receipt_sha256":hashlib.sha256(raw).hexdigest(),"byte_identical":True})
    return result
def algebra_controls():
    count=0
    for s in range(1,5):
        a=s*s;d=2*s+1
        for m in range(1,5):
            A=Q(2,3);B=-Q(3,5);x=A/B
            left=sum((Q(comb(m,k)*factorial(a*m+d*k))*A**(m-k)*B**k for k in range(m+1)),Q(0))
            poly=sum((Q(comb(m,k)*factorial(a*m+d*k),factorial(a*m))*x**(m-k) for k in range(m+1)),Q(0))
            require(left==B**m*factorial(a*m)*poly,"Exact actual moment normalization")
            require(a*m+1+d*(m-1)+d-1==(a+d)*m,"Upper bound endpoint normalization")
            count+=1
    # Rational Bezout identity valid after any characteristic-zero scalar extension.
    # (x²+1)/2+(1-x)(x+1)/2=1.
    p=[Q(1),Q(0),Q(1)];q=[Q(1),Q(1)];u=[Q(1,2)];v=[Q(1,2),-Q(1,2)]
    def mul(a,b):
        out=[Q(0)]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):out[i+j]+=x*y
        return out
    a1,a2=mul(p,u),mul(q,v)
    out=[(a1[i] if i<len(a1) else 0)+(a2[i] if i<len(a2) else 0) for i in range(max(len(a1),len(a2)))]
    require(out==[1,0,0],"Exact rational Bezout identity")
    return count
def negatives(ledger,repo,support):
    names=[]
    def reject(name,fn):
        try:fn()
        except (AssertionError,TypeError,KeyError,ValueError):names.append(name);return
        raise AssertionError("Negative accepted:"+name)
    for role in sorted(ROLES):
        d=copy.deepcopy(ledger);d["sources"]=[r for r in d["sources"] if r["role"]!=role]
        reject("omit:"+role,lambda d=d:check_sources(d,repo,support))
        d=copy.deepcopy(ledger);next(r for r in d["sources"] if r["role"]==role)["sha256"]="0"*64
        reject("badpin:"+role,lambda d=d:check_sources(d,repo,support))
    for key in ["BadTerminal_s_empty","all_terminal_orders","effective_Ns","uniform_Ns_in_s",
                "new_sharp_lower_annulus_logical_premise","uniform_H_bound_logical_premise"]:
        d=copy.deepcopy(SCOPE);d[key]=True
        reject("false-scope:"+key,lambda d=d:check_inputs(d))
    for key in ["same_cutoff_all_coefficients","same_cutoff_all_characteristic_zero_fields",
                "BadTerminal_finite_iff_eventual"]:
        d=copy.deepcopy(SCOPE);d[key]=False
        reject("missing-scope:"+key,lambda d=d:check_inputs(d))
    d=copy.deepcopy(SCOPE);d["order_quantifier"]="one uniform N0 for all s and all N>=2"
    reject("wrong-quantifier",lambda:check_inputs(d))
    for key in ["schema","empty","duplicate"]:
        d=copy.deepcopy(ledger)
        if key=="schema":d["schema"]=True
        if key=="empty":d["sources"]=[]
        if key=="duplicate":d["sources"][-1]=copy.deepcopy(d["sources"][0])
        reject("ledger:"+key,lambda d=d:check_sources(d,repo,support))
    reject("forged-child",lambda:equal_receipt(b'{"status":"PASS"}',b'{"status":"PASS","real":true}'))
    return names
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--repo-root",required=True,type=Path)
    ap.add_argument("--support-root",type=Path,default=HERE)
    ap.add_argument("--output",type=Path,default=HERE/RN);ns=ap.parse_args()
    repo=ns.repo_root.resolve();support=ns.support_root.resolve()
    require(repo.is_dir() and support.is_dir(),"Explicit repo/support roots must exist")
    inp=load(HERE/IN);check_inputs(inp)
    ledger=load(HERE/LN);rows=check_sources(ledger,repo,support)
    links=child_links(repo,support,rows)
    children=replay(repo,support,rows);actual=algebra_controls()
    bad=negatives(ledger,repo,support);check_sources(ledger,repo,support)
    payload={"status":"PASS","schema":1,"scope":inp,"formal_verification":False,"external_review":False,
      "proof_sha256":sha(HERE/PN),"checker_sha256":sha(HERE/CN),"input_sha256":sha(HERE/IN),
      "source_ledger_sha256":sha(HERE/LN),"mandatory_source_roles":len(rows),
      "sources":[rows[r] for r in sorted(rows)],"verified_child_ledgers":links,"fresh_children":children,
      "counts":{"fresh_direct_children":len(children),"recorded_source_ledgers":len(links),
                "exact_actual_moment_normalizations":actual,"rational_Bezout_identities":1,
                "source_scope_controls":len(bad)},
      "negative_names":bad,
      "meaning":"Written analytic inputs plus compact cover and rational Bezout prove the theorem; finite/source replay is supporting evidence, not formal verification or a finite-prefix inference."}
    ns.output.parent.mkdir(parents=True,exist_ok=True);ns.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":"PASS","mandatory_sources":len(rows),"counts":payload["counts"],"output":str(ns.output)}))
if __name__=="__main__":main()

