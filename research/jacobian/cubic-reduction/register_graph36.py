"""Register exact all-degree target obstructions for the nine fixed36-term graph slices."""
from pathlib import Path
import copy,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
BASE="jacobian/cubic-reduction/"
sys.path.insert(0,str(ROOT/"scripts"))
from research import claim_fingerprint
def load(p):return json.loads((ROOT/p).read_text())
def save(p,d):(ROOT/p).write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n")
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
v=load(BASE+"graph36-validation.json")
assert v["status"]=="pass" and v["all_nine_arbitrary_degree_certificates"]
assert v["source_sha256"]==sha(HERE/"verify_graph36.py")
assert v["input_sha256"]==sha(HERE/"cubic12-integer-witness.json")
ledger=load("MATH_STATUS.json");known={e["id"]:e for e in ledger["entries"]}
entry=copy.deepcopy(known["KCUB3"]);identifier="KCUB4"
entry.update(id=identifier,title="All-degree elementary target obstruction for nine fixed cubic graph deletions",
    scope="For the fixed integral36-term twelve-variable cubic KCUB3 map, delete each carrier j=4,...,12 on its level K_j=0 by solving its triangular source graph. In every case, one raw retained component G_i cannot be made degree at most3 by subtracting P of the other ten raw retained outputs, for any polynomial P of any degree over any characteristic-zero field. Obstructed components are2,3,1,1,2,3,2,3,2 for pivots4 through12. Explicit rational coefficient functionals annihilate every source polynomial of degree<=3, evaluate to one on the bad raw output, and annihilate the entire polynomial algebra of the other raw outputs. Their maximum source support orders are6,7,4,4,6,7,6,7,5. A downward-closed monomial quotient preserves all sampled coefficients. Since every raw generator vanishes at zero, target monomials of degree above the support order vanish in that quotient. Independent standard-library Fraction multiplication actions check all lower-degree commutative target words, exact zero-branch coverage and every nonzero defect, with no reliance on primary RREF ranks, SymPy or FLINT. The checker also replays BCR2's successful quadratic square completion as a positive control. This closes increasing the raw target degree for these nine fixed graph routes, and also changes confined to the generated algebra of the other outputs. It does not exclude prior use of the bad output in other target coordinates, general sequential stages, different source coordinates or graph levels, temporary degree growth or other constructions. The separate bounded cubic-pullback kernel search is exploratory only. No minimum cubic dimension or global term count is proved; the positive bounds remain twelve variables/36 displayed terms and the separate normalized three-variable interval [10,48].",
    canonical_source=BASE+"GRAPH_DELETION_OBSTRUCTION.md",dependencies=["KCUB3"],
    checker=BASE+"verify_graph36.py",artifact_hash="sha256:"+sha(HERE/"verify_graph36.py"),
    software_lock=[BASE+"graph-requirements.txt",BASE+"verify.py",BASE+"graph36.py",BASE+"graph36_jet.py",
       BASE+"graph36_all_degree.py",BASE+"graph36_adaptive_wide.py",BASE+"cubic12-integer-witness.json",
       BASE+"graph36-validation.json",BASE+"graph36/checkpoint.json",BASE+"adaptive36-wide/checkpoint.json"]+
       [BASE+folder+"/pivot-%02d.json"%j for folder in ["graph36","adaptive36-wide"] for j in range(4,13)],
    supersedes_notes=["For the current36-term frame, finite dual certificates replace a target-degree cap by an arbitrary-polynomial-degree elementary-stage obstruction. This does not supersede the different-frame BCRDBGRAPH1 theorem or prove global dimension minimality."],
    supersedes=[],replaced_by=[],closes_problems=[],narrows_problems=[],consumers=[],invalidates_assumptions=[])
if identifier in known:known[identifier].clear();known[identifier].update(entry)
else:ledger["entries"].append(entry);known[identifier]=entry
save("MATH_STATUS.json",ledger)
lessons=load("knowledge/lessons.json")
record=dict(id="METHOD-KELLER-FINITE-JET-TARGET-ALGEBRA-20261114",area="geometry",
    title="Certify arbitrary target-degree obstructions by finite source coefficient functionals",
    when="A graph-restricted Keller map has a raw higher-degree component to be repaired by a polynomial in the other raw outputs.",
    use="Seek a rational coefficient functional supported above the allowed source degree, nonzero on the defect, and zero on the raw-output algebra. Close its monomial support downward to form a finite quotient. Origin-vanishing generators imply target words above the support order vanish, so exact finite coverage certifies every polynomial target degree. Check the dual equations, all commutative words, zero-prefix coverage and a known successful completion independently.",
    avoid="Calling finite-jet membership a full polynomial completion; confusing a capped target-degree miss with an arbitrary-degree theorem; accepting producer rank assertions instead of checking the explicit functional; redoing larger raw target-degree bounds after a full algebra obstruction is certified.",
    boundary="KCUB4 independently closes all nine fixed KCUB3 carrier graph slices for one raw-output elementary correction of arbitrary degree. Support orders<=7 suffice. Initial small quotients and the capped adaptive run were inconclusive; the final result uses explicit independently verified separators. Other target sequences, source graphs and global minima remain open.",
    reopen_when="Change the source graph or level, or mix the obstructed component into another output before correction. Transformations confined to the other raw-output algebra cannot escape the certificate. Higher intermediate degrees or a new carrier realization remain admissible.",
    claims=["KCUB3",identifier,"BCR2","BCRDBGRAPH1"],sources=[BASE+"GRAPH_DELETION_OBSTRUCTION.md",BASE+"EXPLORATION.md"],
    implementation=[BASE+"graph36.py",BASE+"graph36_jet.py",BASE+"graph36_all_degree.py",BASE+"graph36_adaptive.py",
                    BASE+"graph36_adaptive_wide.py",BASE+"verify_graph36.py",BASE+"stages36.py"],
    keywords=["Keller","graph deletion","finite jet","dual certificate","monomial quotient","target algebra","arbitrary degree"],
    programme_status="archived",reviewed_claims={})
old=next((r for r in lessons["lessons"] if r["id"]==record["id"]),None)
if old is None:lessons["lessons"].append(record)
else:old.clear();old.update(record)
for item in lessons["lessons"]:
    if item["id"]=="METHOD-KELLER-CUBIC-GROUPED-SHEAR-SPARSITY-20261114":
        addition=" KCUB4 now excludes arbitrary-degree raw elementary target correction for all nine fixed graph deletions in the36-term frame; changing that target-degree cap is not a next step."
        if addition not in item["boundary"]:item["boundary"]+=addition
        item["claims"]=list(dict.fromkeys(item["claims"]+[identifier]))
        item["sources"]=list(dict.fromkeys(item["sources"]+[BASE+"GRAPH_DELETION_OBSTRUCTION.md"]))
    for key in ["KCUB3",identifier,"BCR2","BCRDBGRAPH1"]:
        if key in item.get("claims",[]):item.setdefault("reviewed_claims",{})[key]=claim_fingerprint(known[key])
save("knowledge/lessons.json",lessons)
work=load("knowledge/work_items.json")
for item in work["actions"]:
    if item["id"]=="WORK-GLDC":
        addition=" KCUB4 closes arbitrary-degree elementary target repair for the nine fixed carrier graph deletions of the36-term cubic; any improvement through these slices must change the raw-output algebra, source graph, or permit genuine sequential/intermediate-degree changes. The bounded raw cubic-pullback and shared-cube searches are retained as exploratory evidence, not lower bounds."
        if addition not in item["next_step"]:item["next_step"]+=addition
        item["claims"]=list(dict.fromkeys(item["claims"]+[identifier]))
        item["sources"]=list(dict.fromkeys(item["sources"]+[BASE+"GRAPH_DELETION_OBSTRUCTION.md"]))
    if identifier in item.get("claims",[]):item.setdefault("reviewed_claims",{})[identifier]=claim_fingerprint(entry)
save("knowledge/work_items.json",work)
print("Registered KCUB4: all nine fixed raw graph target repairs excluded in every polynomial degree; global optimization remains open.")

