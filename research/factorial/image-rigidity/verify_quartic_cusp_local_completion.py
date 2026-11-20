#!/usr/bin/env python3
"""Source-guarded transport replay for the three scoped quartic-cusp proofs.

Finite controls support the written arguments; they do not formally verify
Newman's theorem, analytic implicit functions, or asymptotic quantifiers.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile

if not __debug__:
    raise SystemExit("Refusing optimized execution: exact assertions are required.")

HERE=Path(__file__).resolve().parent
SCRIPT="verify_quartic_cusp_local_completion.py"
INPUT="quartic_cusp_local_completion_inputs.json"
LEDGER="quartic_cusp_local_completion_source_ledger.json"
RECEIPT="quartic_cusp_local_completion_verification.json"
ROLES={
 "model_proof":"QUARTIC_LEE_YANG_ZERO_AND_ROBIN_BRIDGE.md",
 "model_checker":"verify_quartic_lee_yang_zero_and_robin_bridge.py",
 "model_inputs":"quartic_lee_yang_zero_and_robin_inputs.json",
 "model_ledger":"quartic_lee_yang_zero_and_robin_source_ledger.json",
 "model_receipt":"quartic_lee_yang_zero_and_robin_verification.json",
 "primary_newman":"quartic_lee_yang_sources/newman_fourier_transforms_1976.pdf",
 "primary_hjn":"quartic_lee_yang_sources/hou_jiang_newman_2208_00917v1.pdf",
 "cusp_proof":"QUARTIC_CUSP_UNFOLDING_AND_ACTUAL_ORDER_DRIFT.md",
 "cusp_checker":"verify_quartic_cusp_unfolding.py",
 "cusp_receipt":"quartic_cusp_unfolding_verification.json",
 "bounded_proof":"BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md",
 "bounded_checker":"verify_bounded_quartic_cusp.py",
 "bounded_receipt":"bounded_quartic_cusp_verification.json",
 "completion_checker":SCRIPT,
 "completion_inputs":INPUT,
}
SCOPE={
 "schema":1,
 "fixed_integer_cusp":True,
 "cusp_parameterization":{"s_min":1,"a":"s^2","d":"2s+1"},
 "every_fixed_finite_scaled_window":True,
 "actual_window":"abs(x/N^d-w_c) <= L*N^(-3/4)",
 "order_quantifier":"exists N0(s,L); all integers N>=N0(s,L)",
 "whole_original_cycle":True,
 "C2_one_sided_epsilon_extension":True,
 "model_Lee_Yang_requires_real_tau":True,
 "effective_N0":False,
 "uniform_in_s":False,
 "growing_scaled_window":False,
 "all_actual_roots_localized":False,
 "analytic_at_epsilon_zero":False,
 "global_cusp_SFC":False,
 "formal_verification":False,
 "external_review":False,
}
CHILDREN=[
 ("model_checker","model_receipt"),
 ("cusp_checker","cusp_receipt"),
 ("bounded_checker","bounded_receipt"),
]
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path):
    with path.open(encoding="utf-8") as stream:return json.load(stream)
def require(condition,message):
    if not condition:raise AssertionError(message)
def check_scope(data):
    require(type(data) is dict and data==SCOPE,"Exact fixed-cusp/fixed-window scope required")
    # Equality of bool and int is allowed by Python; reject that ambiguity.
    for key,value in SCOPE.items():
        require(type(data[key]) is type(value),"Scope type: "+key)
    require(type(data["cusp_parameterization"]["s_min"]) is int,"Integer s_min")
def check_ledger(data,root=HERE):
    require(type(data) is dict,"Source ledger object")
    require(set(data)=={"schema","package","formal_verification","external_review","sources"},
            "Source ledger top-level schema")
    require(type(data["schema"]) is int and data["schema"]==1,"Source schema 1 required")
    require(data["package"]=="quartic-cusp-local-completion","Source package identity")
    require(data["formal_verification"] is False and data["external_review"] is False,
            "No unsupported assurance flags")
    rows=data["sources"]
    require(type(rows) is list and len(rows)==len(ROLES),"All 15 mandatory roles required")
    by_role={}
    paths=set()
    for row in rows:
        require(type(row) is dict and set(row)=={"role","path","sha256"},"Exact source row schema")
        role=row["role"]
        require(type(role) is str and role in ROLES and role not in by_role,"Unique required role")
        require(row["path"]==ROLES[role] and row["path"] not in paths,"Exact unique relative path")
        require(type(row["sha256"]) is str and re.fullmatch(r"[0-9a-f]{64}",row["sha256"]),
                "Lowercase SHA-256")
        path=(root/row["path"]).resolve()
        require(path.is_relative_to(root.resolve()),"Source must remain inside package")
        require(path.is_file() and digest(path)==row["sha256"],"Missing or stale source: "+role)
        by_role[role]=row; paths.add(row["path"])
    require(set(by_role)==set(ROLES),"Complete role set")
    return by_role
def equal_receipt(actual,expected):
    require(actual==expected,"Fresh child receipt differs from pinned executed receipt")
def model_links(by_role):
    model=load(HERE/ROLES["model_ledger"])
    require(type(model) is dict and type(model.get("schema")) is int and model.get("schema")==1,"Model ledger schema")
    require(model.get("formal_verification") is False and model.get("external_review") is False,"Model assurance flags")
    rows=model.get("sources")
    expected={ROLES[r]:by_role[r]["sha256"] for r in
              ["model_proof","model_checker","model_inputs","primary_newman","primary_hjn"]}
    require(type(rows) is list and len(rows)==5,"All model source premises required")
    seen={}
    for row in rows:
        require(type(row) is dict and row.get("path") not in seen,"Unique model source")
        seen[row.get("path")]=row.get("sha256")
    require(seen==expected,"Direct model source pins must agree with completion pins")
def negative_controls(ledger):
    names=[]
    def rejected(name,fn):
        try:fn()
        except (AssertionError,KeyError,ValueError,TypeError):
            names.append(name);return
        raise AssertionError("Negative control accepted: "+name)
    for role in sorted(ROLES):
        trial=copy.deepcopy(ledger)
        trial["sources"]=[row for row in trial["sources"] if row["role"]!=role]
        rejected("omit-role:"+role,lambda trial=trial:check_ledger(trial))
        trial=copy.deepcopy(ledger)
        row=next(row for row in trial["sources"] if row["role"]==role)
        row["sha256"]="0"*64 if row["sha256"]!="0"*64 else "1"*64
        rejected("bad-hash:"+role,lambda trial=trial:check_ledger(trial))
    for name,mutate in [
      ("schema",lambda d:d.update(schema=2)),
      ("empty",lambda d:d.update(sources=[])),
      ("duplicate",lambda d:d["sources"].__setitem__(-1,copy.deepcopy(d["sources"][0]))),
    ]:
        trial=copy.deepcopy(ledger);mutate(trial)
        rejected(name,lambda trial=trial:check_ledger(trial))
    for key in ["effective_N0","uniform_in_s","growing_scaled_window",
                "all_actual_roots_localized","analytic_at_epsilon_zero","global_cusp_SFC"]:
        trial=copy.deepcopy(SCOPE);trial[key]=True
        rejected("false-scope:"+key,lambda trial=trial:check_scope(trial))
    for key in ["fixed_integer_cusp","every_fixed_finite_scaled_window","whole_original_cycle"]:
        trial=copy.deepcopy(SCOPE);trial[key]=False
        rejected("missing-scope:"+key,lambda trial=trial:check_scope(trial))
    trial=copy.deepcopy(SCOPE);trial["order_quantifier"]="all N>=2, all s, all windows"
    rejected("wrong-order-quantifier",lambda:check_scope(trial))
    rejected("forged-child-bytes",lambda:equal_receipt(b'{"status":"PASS"}\n',b'{"status":"PASS","actual":true}\n'))
    return names
def replay_children(repo_root,by_role):
    summaries=[]
    with tempfile.TemporaryDirectory(prefix="quartic-cusp-completion-") as directory:
        outdir=Path(directory)
        for index,(script_role,receipt_role) in enumerate(CHILDREN):
            output=outdir/(str(index)+".json")
            run=subprocess.run([sys.executable,"-B",str(HERE/ROLES[script_role]),
                                "--output",str(output)],
                               cwd=outdir,capture_output=True,text=True,timeout=30)
            require(run.returncode==0,
                    "Child replay failed: "+script_role+"; "+run.stdout+"; "+run.stderr)
            require(output.is_file(),"Child did not produce receipt")
            actual=output.read_bytes()
            equal_receipt(actual,(HERE/ROLES[receipt_role]).read_bytes())
            result=json.loads(actual)
            require(result.get("status")=="PASS","Child did not report PASS")
            require(result.get("checker_sha256")==by_role[script_role]["sha256"],
                    "Child checker identity")
            proof_role={"model_checker":"model_proof","cusp_checker":"cusp_proof",
                        "bounded_checker":"bounded_proof"}[script_role]
            require(result.get("proof_sha256")==by_role[proof_role]["sha256"],
                    "Child proof identity")
            summaries.append({"checker_role":script_role,
                              "receipt_sha256":hashlib.sha256(actual).hexdigest(),
                              "byte_identical":True})
    return summaries
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--repo-root",required=True,type=Path,
                    help="Explicit repository working context; no repository source is a hidden premise.")
    ap.add_argument("--output",type=Path,default=HERE/RECEIPT)
    args=ap.parse_args()
    require(args.repo_root.resolve().is_dir(),"Explicit repository root must exist")
    scope=load(HERE/INPUT);check_scope(scope)
    ledger=load(HERE/LEDGER);by_role=check_ledger(ledger)
    model_links(by_role)
    children=replay_children(args.repo_root.resolve(),by_role)
    negatives=negative_controls(ledger)
    # Ensure child execution did not change the pinned core files.
    check_ledger(ledger)
    payload={
      "status":"PASS","schema":1,"formal_verification":False,"external_review":False,
      "scope":scope,"mandatory_source_roles":len(ROLES),
      "source_ledger_sha256":digest(HERE/LEDGER),
      "checker_sha256":digest(HERE/SCRIPT),"input_sha256":digest(HERE/INPUT),
      "sources":[by_role[role] for role in sorted(by_role)],
      "fresh_children":children,"fresh_replays":len(children),
      "negative_controls":len(negatives),"negative_names":negatives,
      "interpretation":"Source and finite-algebra replay supports the complete written proofs; no effective N0 or global cusp SFC claim.",
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":"PASS","mandatory_sources":len(ROLES),
                      "fresh_replays":len(children),"negative_controls":len(negatives),
                      "output":str(args.output)}))
if __name__=="__main__":main()

