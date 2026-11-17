"""Register a checked cubic witness and retain all open optimization questions."""
from pathlib import Path
import copy,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
BASE="jacobian-followup/cubic-reduction/"
sys.path.insert(0,str(ROOT/"scripts"))
from research import claim_fingerprint
def load(p):return json.loads((ROOT/p).read_text())
def save(p,d):(ROOT/p).write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n")
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
receipt=load(BASE+"validation.json")
assert receipt["status"]=="pass" and receipt["dimension"]==17 and receipt["degree"]==3
assert receipt["source_sha256"]==sha(HERE/"verify.py")
assert receipt["witness_sha256"]==sha(HERE/"cubic17-witness.json")
assert receipt["primary_source_sha256"]==sha(HERE/"construct17.py")
ledger=load("MATH_STATUS.json");known={e["id"]:e for e in ledger["entries"]}
if known.get("KCUB1",{}).get("replaced_by"):
    raise SystemExit("Historical registration retired after the baseline correction; use register12.py.")
identifier="KCUB1";entry=copy.deepcopy(known["KFRAME3"])
entry.update(
    id=identifier,title="Seventeen-variable cubic Keller collision and the degree-three threshold",
    scope="An explicit rational polynomial map in17 variables has ordinary total degree3, determinant one, H(0)=0, JH(0)=I, and three distinct rational points mapping to (0,0,-1/4,0,...,0). Its displayed coordinate system has58 nonlinear monomial occurrences. The common collision fibre is not the origin fibre; affine collision normalization preserves degree but not that term count. The map removes two redundant carrier directions from the 19-variable cubic reduction posted by rhicksrad on20 July2026, attributed there to Alpoge's foundational example and the classical Bass-Connell-Wright stable reduction. The retained combinations are u=6*w11-w12 and v=3*w14-w15; complementary directions disappear from every retained output, leaving two removable identity coordinates after source shears. A standard-library Fraction polynomial replay independently verifies the entire substitution from the pinned external coefficient data, the identity derivative, three full rational collisions, nilpotence of the14-by14 carrier Jacobian correction with exponent3, both finite-series inverse identities, and the full symbolic Schur determinant equal to one. For any degree-at-most-two map Q in characteristic zero, Q(a)-Q(b)=JQ((a+b)/2)*(a-b); a Keller map therefore has no quadratic collision in any dimension. Hence the minimum ordinary degree for characteristic-zero Keller collisions across all finite dimensions is exactly3. This degree threshold is classical, not a novelty claim. Seventeen is only an upper bound on cubic dimension. No minimum dimension, cubic support optimum, three-variable degree improvement, cubic-homogeneous or symmetric-Jacobian normal form, global affine-orbit optimum, formal verification or external review is asserted. The three-variable normalized raw-degree-seven support interval remains [10,48]. Failed bounded monomial splitting and selected coordinate-graph slices are preserved as discovery evidence only.",
    canonical_source=BASE+"PROOF.md",dependencies=["KFRAME1"],
    checker=BASE+"verify.py",artifact_hash="sha256:"+sha(HERE/"verify.py"),
    software_lock=[BASE+"requirements.txt",BASE+"construct17.py",BASE+"cubic17-witness.json",
        BASE+"validation.json",BASE+"sources/degree3-map-external.json",BASE+"sources/provenance.json"],
    closes_problems=[],narrows_problems=[],supersedes=[],supersedes_notes=[],
    consumers=[],invalidates_assumptions=[],replaced_by=[],independent_replay=True)
if identifier in known:
    known[identifier].clear();known[identifier].update(entry)
else:ledger["entries"].append(entry)
save("MATH_STATUS.json",ledger)
lessons=load("knowledge/lessons.json")
record=dict(
    id="METHOD-KELLER-CUBIC-CARRIER-DESTABILIZATION-20261114",area="geometry",
    title="Remove redundant carrier directions before further cubic degree reduction",
    when="A stable cubic reduction introduces many auxiliary carrier variables.",
    use="Compute how top coordinates and all retained carrier outputs depend on carrier directions. Retain independent combinations, apply the same linear change to the corresponding targets, and check whether complementary outputs become identity coordinates plus polynomials of retained inputs. If so, source shears remove those directions exactly without raising retained degree. Verify the complete symbolic determinant using an explicitly checked nilpotent carrier-block inverse, and transport full collision points.",
    avoid="Calling every coordinate-graph slice cubic; conflating ordinary degree with generic fibre degree; comparing term counts across different collision normalizations; claiming a beam search or failed family of target corrections proves minimum dimension.",
    boundary="KCUB1 reduces the pinned external 19-variable construction to17 with58 displayed nonlinear occurrences and a fully independently verified determinant and three-point collision. The quadratic midpoint identity proves the global degree-three collision threshold, already classical. Dimension17 is an upper bound only. Direct slices and the retained bounded splitting beam do not improve it; selected quadratic carrier-target correction systems are discovery diagnostics, not global exclusions.",
    reopen_when="The carrier realization changes, simultaneous or general-direction slices become available, higher-degree target corrections are permitted, or a different cubic construction is supplied.",
    claims=[identifier],sources=[BASE+"PROOF.md",BASE+"EXPLORATION.md",BASE+"sources/provenance.json"],
    implementation=[BASE+"construct17.py",BASE+"verify.py",BASE+"search.py",BASE+"slice_external.py",BASE+"repair_slices.py"],
    keywords=["Keller","cubic","ordinary degree","dimension17","carrier","stable equivalence","Schur complement","midpoint"],
    reviewed_claims={identifier:claim_fingerprint(entry)},programme_status="archived")
old=next((x for x in lessons["lessons"] if x["id"]==record["id"]),None)
if old is None:lessons["lessons"].append(record)
else:old.clear();old.update(record)
save("knowledge/lessons.json",lessons)
work=load("knowledge/work_items.json")
for item in work["actions"]:
    if item["id"]=="WORK-GLDC":
        addition=" Across varying dimensions, KCUB1 verifies the degree-three collision threshold and an explicit cubic witness in17 variables with58 displayed nonlinear terms; its minimum cubic dimension and support are separate open targets. Seek a smaller carrier realization or transformations beyond the retained failed slices, rather than treating17 as optimal."
        if addition not in item["next_step"]:item["next_step"]+=addition
        item["claims"]=list(dict.fromkeys(item["claims"]+[identifier]))
        item["sources"]=list(dict.fromkeys(item["sources"]+[BASE+"PROOF.md",BASE+"EXPLORATION.md"]))
        item.setdefault("reviewed_claims",{})[identifier]=claim_fingerprint(entry)
save("knowledge/work_items.json",work)
print("Registered KCUB1: verified cubic dimension17 upper bound; full optimization goal remains open.")

