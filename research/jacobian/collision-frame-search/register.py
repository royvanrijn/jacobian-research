"""Register the verified 48-term frame and the two fully covered transverse families."""
from pathlib import Path
import copy,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
BASE="jacobian/collision-frame-search/"
sys.path.insert(0,str(ROOT/"scripts"))
from research import claim_fingerprint
from research_partial_reviews import entry_fingerprint
def load(p):return json.loads((ROOT/p).read_text())
def save(p,v):(ROOT/p).write_text(json.dumps(v,indent=2,ensure_ascii=False)+"\n")
receipt=load(BASE+"validation.json")
assert receipt["status"]=="pass" and receipt["normalized_support_interval"]==[10,48]
assert receipt["source_sha256"]==hashlib.sha256((HERE/"verify.py").read_bytes()).hexdigest()
assert receipt["singular_sha256"]==hashlib.sha256((HERE/"singular-validation.json").read_bytes()).hexdigest()
ledger=load("MATH_STATUS.json");known={e["id"]:e for e in ledger["entries"]}
entry=copy.deepcopy(known["KFRAME1"])
entry.update(
    id="KFRAME2",title="48-term normalized Keller frame and two sharp transverse minima",
    scope="For the foundational three-variable Keller counterexample, the transverse shear S=I+14*E_32 conjugates KFRAME1's 49-term normalized map to a rational map with determinant one, JF(0)=I, F(0)=F(e_1)=0, coordinate degrees (4,6,7), and nonlinear support counts (6,19,23), total48. This improves the broader normalized raw-degree-at-most-seven support interval to [10,48] using GLDC4's lower bound. The value48 is sharp in two specified families over every characteristic-zero field: all boundary collision pairs in KFRAME1 with arbitrary lower triangular transverse completions, and every GL2 transverse basis in the y,z-plane at the fixed pair (0,-2,-16),(1,-2,8). Nonzero column scalings and a column swap reduce the latter class to a product of two unipotent shears, without singular parameter strata. Complete coefficient-curve, leading-coefficient and resultant projections give469 and824 irreducible critical parameter factors. Exact quotient-field root bases cover every specialization. Singular independently verifies all1293 minima, the irreducible projection factors, every resultant and projection factor product, complete coefficient products and coprime root bases. Standard-library rational polynomial arithmetic independently reconstructs both input families and verifies the actual48-term determinant, identity derivative and collision. The discarded v1 field-specialization implementation is explicitly invalidated and preserved; the theorem uses v2 with distinct generators, exact modulus checks and independent Singular replay. Bounded beam exploration is only discovery evidence. Arbitrary transverse x-components, joint variation of the boundary pair with both GL2 parameters, interior pairs and global support optimality remain open. The profile is still exactly (7,6,4), so this is not a strictly lower-profile witness for OP-GLDC. No ordinary-degree improvement, new counterexample mechanism, formal verification, external review or novelty priority is claimed.",
    canonical_source=BASE+"PROOF.md",dependencies=["KFRAME1","GLDC4"],
    checker=BASE+"verify.py",proof_type="exact_symbolic",
    artifact_hash="sha256:"+hashlib.sha256((HERE/"verify.py").read_bytes()).hexdigest(),
    software_lock=[BASE+"requirements.txt",BASE+"validation.json",BASE+"singular-validation.json",
                   BASE+"verify_singular.py",BASE+"specialize.py",BASE+"two_parameter.py",
                   BASE+"full_transverse.py",BASE+"field-regression.json",
                   BASE+"normalized-48-witness.json",BASE+"specialization-manifest-v2.json",
                   BASE+"full-transverse/specialization-manifest-v2.json"],
    supersedes_notes=["The broader normalized raw-degree-seven upper bound49 is improved to48 outside KFRAME1's fixed transverse-axis frame class."],
    supersedes=[],closes_problems=[],narrows_problems=[],consumers=[],
    invalidates_assumptions=[],replaced_by=[],independent_replay=True)
if "KFRAME2" in known:
    assert known["KFRAME2"]["checker"]==entry["checker"]
    known["KFRAME2"].clear();known["KFRAME2"].update(entry)
else:ledger["entries"].append(entry);known["KFRAME2"]=entry
parent=known["OP-GLDC"]
old="The true minimum in the broader normalized raw-degree-at-most-seven domain is between ten and 49."
new=("KFRAME2 improves the broader normalized raw-degree-at-most-seven interval to [10,48] and proves48 sharp in two extended transverse frame families. Its profile remains exactly (7,6,4); affine frame changes preserve this flag and cannot themselves solve the strict lower-profile objective.")
assert old in parent["scope"] or new in parent["scope"]
parent["scope"]=parent["scope"].replace(old,new)
parent["dependencies"]=list(dict.fromkeys(parent["dependencies"]+["KFRAME2"]))
save("MATH_STATUS.json",ledger)

lessons=load("knowledge/lessons.json")
for item in lessons["lessons"]:
    if item["id"]=="METHOD-KELLER-BOUNDARY-FRAME-SPARSITY-20261114":
        item["use"]="Parameterize actual collision frames and use exact coefficient-zero arrangements. For two parameters, project every vertical coefficient factor, every moving-factor leading coefficient and every pair resultant; solve every critical parameter factor using coprime root bases in the remaining parameter. Verify all parameter inputs independently and replay every algebraic specialization in a separate backend."
        item["avoid"]="Calling a frame-family minimum global; sampling rational parameters instead of all algebraic support-drop strata; treating bounded beam exploration as affine-orbit coverage; using a profile-equality frame for OP-GLDC's strict lower-profile objective."
        item["boundary"]="KFRAME1 has sharp minimum49 in fixed-axis boundary frames. KFRAME2 has sharp minimum48 for boundary pairs with lower triangular transverse completions and, separately, for all y,z-plane transverse bases at the chosen fixed pair. The broader normalized degree-seven interval is [10,48]. Transverse x-components, joint boundary-position plus full-GL2 variation and interior pairs remain outside these complete results."
        item["reopen_when"]="The collision point stratum, transverse x-components, joint parameter family, Keller mechanism or permitted coordinate transformations change."
        item["claims"]=list(dict.fromkeys(item["claims"]+["KFRAME2"]))
        item["sources"]=list(dict.fromkeys(item["sources"]+[BASE+"PROOF.md"]))
        item["implementation"]=list(dict.fromkeys(item["implementation"]+
            [BASE+"two_parameter.py",BASE+"specialize.py",BASE+"verify_singular.py",BASE+"verify.py"]))
    if item["id"]=="METHOD-KELLER-AXIS-SEEDED-CLOSURE-20261114":
        item["boundary"]=item["boundary"].replace("upper bound 49","upper bound48 via KFRAME2")
        item["claims"]=list(dict.fromkeys(item["claims"]+["KFRAME2"]))
    for key in ("KFRAME2","OP-GLDC"):
        if key in item.get("claims",[]):item["reviewed_claims"][key]=claim_fingerprint(known[key])
identifier="METHOD-FORMAL-ALGEBRAIC-GENERATOR-CACHE-20261114"
record=dict(
    id=identifier,area="geometry",title="Keep algebraic field generators and coefficient moduli distinct",
    when="Many exact number fields are constructed from different minimal polynomials in one process.",
    use="Do not reuse one bare formal root expression across different supplied minimal polynomials. Use distinct formal generators for residue arithmetic, canonicalize the residue generator and assert the exact monic modulus both in the coefficient field and its polynomial ring. Independently replay the specializations from literal minimal polynomials in another backend and include known rational points as controls. Parenthesize polynomial powers in Singular expressions before rational division.",
    avoid="Accepting cached algebraic-field domains merely because their printed aliases look correct. The rejected v1 reuse identified QQ[w]/(w-2) with a request for QQ[w]/(w-14), producing wrong support minima. All v1 ledgers are invalidated rather than silently overwritten; old sparse proofs used no such field construction.",
    boundary="The repaired frame specialization uses distinct Dummy generators with modulus assertions; Singular verifies all1293 field cases. The explicit48-term witness is independently checked in Fraction polynomial arithmetic. This does not independently audit arbitrary SymPy number-field workflows.",
    reopen_when="A field constructor, algebraic alias scheme, polynomial coefficient domain or serialization layer changes.",
    claims=["KFRAME2"],sources=[BASE+"INVALIDATED_SPECIALIZATION.md",BASE+"PROOF.md"],
    implementation=[BASE+"check_field_generators.py",BASE+"specialize.py",BASE+"verify_singular.py"],
    keywords=["number field","cache","minimal polynomial","generator","Singular","rational exponent","invalidated"],
    reviewed_claims={"KFRAME2":claim_fingerprint(entry)},programme_status="archived")
existing=next((e for e in lessons["lessons"] if e["id"]==identifier),None)
if existing is None:lessons["lessons"].append(record)
else:existing.clear();existing.update(record)
save("knowledge/lessons.json",lessons)

work=load("knowledge/work_items.json")
for item in work["actions"]:
    if item["id"]=="WORK-GLDC":
        item.update(
            title="Separate the strict degree-profile census from broader support minimization",
            next_step="For OP-GLDC, seek size-ten candidates passing the coordinate gates in a flag strictly below (7,6,4), or close the next full support stratum. Separately, improve the broader degree-seven upper bound48 only by leaving KFRAME2's two completed frame classes: allow transverse x-components, jointly vary boundary position and full GL2, use interior collision pairs or another Keller family. Affine changes of the foundational map cannot lower its invariant degree flag.",
            done_when="A certified strictly lower-profile witness or complete next exclusion for OP-GLDC; separately a verified normalized degree-seven example below48 advances the broader sparsity goal. Absolute optimality requires matching bounds.",
            prerequisites="The broader normalized degree-seven interval is [10,48]. The upper witness has flag exactly (7,6,4), outside OP-GLDC's strict range. KFRAME1 and KFRAME2's specified frame minima are complete. Do not reuse invalid v1 number-field specializations or treat the bounded beam as complete affine-orbit coverage.")
        item["claims"]=list(dict.fromkeys(item["claims"]+["KFRAME2"]))
        item["sources"]=list(dict.fromkeys(item["sources"]+[BASE+"PROOF.md",BASE+"INVALIDATED_SPECIALIZATION.md"]))
    for key in ("KFRAME2","OP-GLDC"):
        if key in item.get("claims",[]):item["reviewed_claims"][key]=claim_fingerprint(known[key])
save("knowledge/work_items.json",work)
reviews=load("knowledge/partial_reviews.json")
for group in ("reviews","checker_absence_reviews"):
    for item in reviews.get(group,[]):
        if "OP-GLDC" in item.get("reviewed_claims",{}):
            item["reviewed_claims"]["OP-GLDC"]=entry_fingerprint(parent)
save("knowledge/partial_reviews.json",reviews)
print("Registered KFRAME2: normalized interval [10,48], two sharp frame families; strict OP-GLDC and the overall goal remain open.")

