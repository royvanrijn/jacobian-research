"""Register the completed seven-term exclusion without reactivating archived queues."""
import copy
import hashlib
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research import claim_fingerprint
from research_partial_reviews import entry_fingerprint
HERE=Path(__file__).resolve().parent
BASE="jacobian/sparse-seven/"

def load(path):
    return json.loads((ROOT/path).read_text())
def save(path,value):
    (ROOT/path).write_text(json.dumps(value,indent=2,ensure_ascii=False)+"\n")

ledger=load("MATH_STATUS.json")
known={e["id"]:e for e in ledger["entries"]}
entry=copy.deepcopy(known["GLDC1"])
entry.update(
    id="GLDC2",title="Complete seven-term degree-seven Keller collision exclusion",
    scope="Over every characteristic-zero field, every three-variable Keller map in the collision frame F(0)=F(e_1)=0 and JF(0)=I, with raw coordinate degrees at most seven, has at least eight nonlinear monomial occurrences counted with coordinate multiplicity. The unimodular-axis necessity reduces all candidates to 90 minimal-axis-exponent shards. Complete determinant-bucket closure leaves 1876 seven-term support orbits and 450 six-term regression orbits. Every exact-support coefficient-torus ideal is (1) over Q in both Singular and SymPy. A separate standard-library implementation reconstructs the entire support cover using an unsharded traversal and an independent Leibniz determinant expansion, with exact set equality; an independent derivative-level compiler verifies every coefficient ideal, including all collision equations and rho*product(c)-1. Checkpoint resumption and agreement with the original degree-three enumerator are tested. Thus seven-term attainment in the normalized degree-at-most-seven range is impossible. Eight-term attainment, global degree minima, unrestricted-degree support minima, and special-form dimension minima remain open. This does not close any entire dense profile above maximum degree two and is neither proof-assistant formalization nor external review.",
    canonical_source=BASE+"PROOF.md",dependencies=["GLDC1","CAUF1"],
    checker=BASE+"verify.py",proof_type="exact_symbolic",independent_replay=True,
    formal_verification=False,external_review=False,
    artifact_hash="sha256:"+hashlib.sha256((HERE/"verify.py").read_bytes()).hexdigest(),
    software_lock=[BASE+"requirements.txt",BASE+"validation.json",BASE+"results/manifest.json",
                   BASE+"search.py",BASE+"run_all.py",BASE+"verify_coverage.py"],
    supersedes=[],supersedes_notes=["Support-seven attainment was previously unknown in the normalized raw-degree-at-most-seven census."],
    closes_problems=[],narrows_problems=["OP-GLDC"],consumers=[],
    invalidates_assumptions=[],replaced_by=[],priority="reference",programme_status="archived")
if "GLDC2" in known:
    known["GLDC2"].clear()
    known["GLDC2"].update(entry)
else:
    ledger["entries"].append(entry)
    known["GLDC2"]=entry
parent=known["OP-GLDC"]
parent["scope"]=(
    "Complete the support-first census of every three-dimensional Keller collision whose affine-invariant degree flag is lexicographically below (7,6,4). GLDC1 supplies the collision normalization, all 74 profiles, the determinant-bucket and stratified infinity-face compilers, and four dense maximum-degree-two profile exclusions. CAUF1 supplies the unimodular collision-axis frontend. GLDC2 completely excludes nonlinear support at most seven in raw coordinate degree at most seven: all 90 axis-exponent shards are covered and every one of 1876 seven-term and 450 six-term residual orbit ideals is (1) over Q in Singular and SymPy, with independent support enumeration and coefficient compilation. The current normalized support lower bound is eight; attainment at eight is unknown. The remaining census starts at eight, followed by filtered-output rank minors, leading-coefficient torus solvability and actual boundary compactification for genuine coefficient survivors. The 105-digit raw support space rules out flat unstructured enumeration. No profile beyond maximum degree two is eliminated in full without a support bound, and no absolute degree or special-form dimension minimum is established.")
parent["dependencies"]=list(dict.fromkeys(parent["dependencies"]+["GLDC2"]))
parent["forbidden_attack_classes"][0]=dict(
    attack="Repeating normalized raw-degree-at-most-seven searches with nonlinear support at most seven",
    reason="GLDC1 and GLDC2 exhaust this domain. All residual exact-support coefficient ideals are unit over Q, with two rational backends and independent enumeration and coefficient-input validation.",
    witnesses=["GLDC1","GLDC2"])
parent["forbidden_attack_classes"][1]["reason"]=(
    "The axis-only support minimum is three and is merely a frontend. The global normalized lower bound eight uses full determinant and coefficient compatibility. No dense profile above maximum degree two is fully excluded.")
save("MATH_STATUS.json",ledger)

lessons=load("knowledge/lessons.json")
identifier="METHOD-KELLER-AXIS-SEEDED-CLOSURE-20261114"
lesson=dict(
    id=identifier,area="geometry",title="Seed sparse Keller searches with necessary collision-axis terms",
    when="A normalized three-variable collision is searched by exact nonlinear support.",
    use="Start from the smallest F1 pure-axis exponent and the two smallest pure-axis exponents in one transverse row. Apply complete determinant-singleton repair, retain arbitrary balanced supersets, and solve every residual coefficient-torus ideal over Q. Keep a source-pinned DFS frontier and explicit terminal status. Verify set-level coverage with an independently implemented global traversal and independently compile the actual derivative coefficient systems.",
    avoid="Treating support-seven attainment as open in the completed raw-degree-at-most-seven domain; enumerating the entire 348-atom Boolean universe without the unimodular-axis frontend; equating matching counts or two solvers on unchecked common inputs with coverage.",
    boundary="The 90-shard cover and all 1876 seven-term residual ideals are closed. The lower bound eight concerns the normalized frame and raw degree at most seven. It gives neither eight-term attainment nor an unrestricted degree, support or dimension minimum.",
    reopen_when="A candidate has at least eight normalized nonlinear occurrences, exceeds the raw degree bound, or introduces a structural theorem outside the completed finite domain.",
    claims=["GLDC2","OP-GLDC"],sources=[BASE+"PROOF.md"],
    implementation=[BASE+"search.py",BASE+"verify_coverage.py",BASE+"verify.py"],
    keywords=["Keller","sparse","axis","seven","eight","checkpoint","coverage","saturation"],
    reviewed_claims={k:claim_fingerprint(known[k]) for k in ("GLDC2","OP-GLDC")},
    programme_status="archived")
old=next((e for e in lessons["lessons"] if e["id"]==identifier),None)
if old is None:lessons["lessons"].append(lesson)
else:old.clear();old.update(lesson)
for item in lessons["lessons"]:
    if "OP-GLDC" in item["reviewed_claims"]:
        item["reviewed_claims"]["OP-GLDC"]=claim_fingerprint(parent)
save("knowledge/lessons.json",lessons)
work=load("knowledge/work_items.json")
for item in work["actions"]:
    if "OP-GLDC" in item["reviewed_claims"]:
        item["reviewed_claims"]["OP-GLDC"]=claim_fingerprint(parent)
    if item["id"]=="WORK-GLDC":
        item.update(
            title="Find an eight-term Keller collision or a structural lower bound",
            next_step="Seven terms are completely excluded by GLDC2. Seek an eight-term normalized collision in raw degree at most seven, or prove a transverse-extension or other structural obstruction that closes the next necessary support shapes. Retain exact coefficient torus, whole-support coverage and checkpoint boundaries.",
            done_when="A certified eight-term witness, or a complete next exclusion with its full combinatorial and coefficient coverage.",
            prerequisites="The normalized lower bound is eight, with attainment unknown. Exactly four dense quadratic profiles are closed, while the other degree flags remain open. The requested seven-term follow-up is complete; repeated searches in support at most seven require a changed hypothesis or a concrete audit issue.",
        )
        item["claims"]=list(dict.fromkeys(item["claims"]+["GLDC2"]))
        item["sources"]=list(dict.fromkeys(item["sources"]+[BASE+"PROOF.md"]))
        item["reviewed_claims"]["GLDC2"]=claim_fingerprint(entry)
save("knowledge/work_items.json",work)
reviews=load("knowledge/partial_reviews.json")
for group in ("reviews","checker_absence_reviews"):
    for item in reviews.get(group,[]):
        if "OP-GLDC" in item.get("reviewed_claims",{}):
            item["reviewed_claims"]["OP-GLDC"]=entry_fingerprint(parent)
save("knowledge/partial_reviews.json",reviews)
print("Registered GLDC2; OP-GLDC remains open, lower bound eight.")

