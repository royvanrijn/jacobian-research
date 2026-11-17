"""Register the verified thirteen-variable, 42-term cubic existence result."""
from pathlib import Path
import copy,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
BASE="jacobian/cubic-reduction/"
sys.path.insert(0,str(ROOT/"scripts"))
from research import claim_fingerprint
def load(p):return json.loads((ROOT/p).read_text())
def save(p,d):(ROOT/p).write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n")
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
receipt=load(BASE+"sparse-validation13.json")
assert receipt["status"]=="pass" and receipt["dimension"]==13 and receipt["nonlinear_total"]==42
assert receipt["source_sha256"]==sha(HERE/"verify_sparse13.py")
assert receipt["witness_sha256"]==sha(HERE/"cubic13-sparse-witness.json")
assert receipt["primary_source_sha256"]==sha(HERE/"construct_sparse13.py")
ledger=load("MATH_STATUS.json");known={e["id"]:e for e in ledger["entries"]}
if known.get("KCUB2",{}).get("replaced_by"):
    raise SystemExit("Historical registration retired after the baseline correction; use register12.py.")
entry=copy.deepcopy(known["KCUB1"]);identifier="KCUB2"
entry.update(
    id=identifier,title="Thirteen-variable cubic Keller collision with42 nonlinear occurrences",
    scope="An explicit rational polynomial map in13 variables has ordinary total degree3, determinant one, K(0)=0, JK(0)=I, and three distinct rational points mapping to (0,-2,0,0,...,0). Its displayed nonlinear monomial occurrence counts are (4,12,15,1,1,1,2,1,1,1,1,1,1), total42. The common collision fibre is not the origin fibre; affine normalization at a collision preserves degree but need not preserve42 terms. Starting from the tangent-normalized foundational degree-seven map,10 carrier extensions,13 elementary product target steps and4 linear target steps reuse existing output factors to give a47-term cubic in13 variables. Rational derivative normalization and the transverse conjugation S=I-3*E_(7,11), in one-based source indices, produce the42-term map. A separate SymPy construction and standard-library Fraction replay agree on every coefficient, every elementary transformation and all three transported rational points. The independent checker verifies the original full determinant, the constant derivative inverse, the final identity derivative, the entire final13-variable polynomial Schur determinant equal to one, and both inverse identities for the10-by10 carrier block I+N with N^3=0. This improves the constructed cubic dimension upper bound17 to13 and the displayed cubic occurrence count58 to42; no minimum dimension or support is proved. The degree-three collision threshold across varying dimensions remains the classical KCUB1 consequence of the quadratic midpoint identity. Beam widths60 and800,18 specified variable-order runs through dimension12, selected graph slices, and rational elementary-line descent are discovery evidence only, not affine-orbit or dimension coverage. No new counterexample mechanism, novelty priority, three-variable degree/profile improvement, cubic-homogeneous or symmetric-Jacobian normal form, formal verification or external review is claimed. The separate three-variable normalized raw-degree-seven support interval remains [10,48], with attainment of the lower bound ten and global optimality unresolved.",
    canonical_source=BASE+"PROOF13.md",dependencies=["KFRAME1","KCUB1"],
    checker=BASE+"verify_sparse13.py",artifact_hash="sha256:"+sha(HERE/"verify_sparse13.py"),
    software_lock=[BASE+"requirements.txt",BASE+"construct_sparse13.py",BASE+"verify13.py",BASE+"verify.py",
      BASE+"construct13.py",BASE+"search_reuse.py",BASE+"search.py",BASE+"optimize13.py",
      BASE+"reuse-witness.json",BASE+"cubic13-witness.json",BASE+"cubic13-sparse-candidate.json",
      BASE+"cubic13-sparse-witness.json",BASE+"sparse-validation13.json"],
    supersedes_notes=["The explicit constructed cubic upper bound dimension17 and displayed count58 from KCUB1 are improved to dimension13 and42 occurrences. KCUB1's theorem and evidence remain valid."],
    supersedes=[],closes_problems=[],narrows_problems=[],consumers=[],invalidates_assumptions=[],replaced_by=[])
if identifier in known:known[identifier].clear();known[identifier].update(entry)
else:ledger["entries"].append(entry)
save("MATH_STATUS.json",ledger)
lessons=load("knowledge/lessons.json")
record=dict(
    id="METHOD-KELLER-CUBIC-OUTPUT-FACTOR-REUSE-20261114",area="geometry",
    title="Reuse existing output factors in elementary cubic degree reduction",
    when="A polynomial Keller map is reduced in degree by adding auxiliary carrier coordinates.",
    use="Use existing output coordinates as factors in elementary target product steps, adding a carrier only for a missing factor. Verify the target row differs from both factor rows. Track all source extensions and collision coordinates, then normalize the constant derivative and recount support. A later linear conjugation can combine source carrier polynomials and cancel several top terms. Replay the entire path in separate arithmetic and compute the full final determinant.",
    avoid="Introducing fresh carriers for both factors by default; assuming a discovery score survives derivative normalization; identifying a greedy rational elementary-line minimum with the whole affine orbit; calling a bounded monomial beam a dimension lower bound; comparing term counts across different collision normalizations.",
    boundary="KCUB2 gives a13-variable cubic with47 nonlinear terms before S=I-3*E_(7,11), and42 after it. Primary symbolic arithmetic and independent rational replay verify the whole map and determinant. The18 variable-order runs through dimension12 and selected graph slices find no smaller witness only within their bounded scopes. The final cubic has no untested claim of optimal dimension, support or special form.",
    reopen_when="Use polynomial carriers beyond monomials, a different factor realization, simultaneous or general-direction slices, nonlinear transformations, or another Keller family. The three-variable strict lower-profile census remains a separate objective.",
    claims=[identifier],sources=[BASE+"PROOF13.md",BASE+"EXPLORATION.md",BASE+"orders-manifest.json"],
    implementation=[BASE+"search_reuse.py",BASE+"search_reuse_orders.py",BASE+"construct13.py",
       BASE+"verify13.py",BASE+"optimize13.py",BASE+"construct_sparse13.py",BASE+"verify_sparse13.py"],
    keywords=["Keller","cubic","dimension13","42 terms","output factor","carrier reuse","tame","Schur"],
    reviewed_claims={identifier:claim_fingerprint(entry)},programme_status="archived")
old=next((x for x in lessons["lessons"] if x["id"]==record["id"]),None)
if old is None:lessons["lessons"].append(record)
else:old.clear();old.update(record)
for item in lessons["lessons"]:
    if item["id"]=="METHOD-KELLER-CUBIC-CARRIER-DESTABILIZATION-20261114":
        addition=" KCUB2 subsequently improves the constructed upper bound to13 variables and42 displayed nonlinear occurrences by output-factor reuse and a transverse shear; neither is proved optimal."
        if addition not in item["boundary"]:item["boundary"]+=addition
        item["claims"]=list(dict.fromkeys(item["claims"]+[identifier]))
        item["sources"]=list(dict.fromkeys(item["sources"]+[BASE+"PROOF13.md"]))
        item.setdefault("reviewed_claims",{})[identifier]=claim_fingerprint(entry)
save("knowledge/lessons.json",lessons)
work=load("knowledge/work_items.json")
for item in work["actions"]:
    if item["id"]=="WORK-GLDC":
        old="Across varying dimensions, KCUB1 verifies the degree-three collision threshold and an explicit cubic witness in17 variables with58 displayed nonlinear terms; its minimum cubic dimension and support are separate open targets. Seek a smaller carrier realization or transformations beyond the retained failed slices, rather than treating17 as optimal."
        new="Across varying dimensions, KCUB1 supplies the classical degree-three collision threshold. KCUB2 improves the explicit cubic witness to13 variables with42 displayed nonlinear occurrences, outside the three-variable collision normalization. Seek a smaller or sparser cubic through polynomial carriers, another factor realization or transformations beyond the bounded searches; minimum cubic dimension and support remain open."
        assert old in item["next_step"] or new in item["next_step"]
        item["next_step"]=item["next_step"].replace(old,new)
        item["claims"]=list(dict.fromkeys(item["claims"]+[identifier]))
        item["sources"]=list(dict.fromkeys(item["sources"]+[BASE+"PROOF13.md"]))
        item.setdefault("reviewed_claims",{})[identifier]=claim_fingerprint(entry)
save("knowledge/work_items.json",work)
print("Registered KCUB2: cubic dimension13 upper bound and42 displayed nonlinear occurrences. Full simplicity goal remains open.")

