"""Register the independent endpoint trace obstruction without changing support bounds."""
from pathlib import Path
import copy,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
BASE="jacobian/collision-frame-search/"
sys.path.insert(0,str(ROOT/"scripts"))
from research import claim_fingerprint
def load(p):return json.loads((ROOT/p).read_text())
def save(p,v):(ROOT/p).write_text(json.dumps(v,indent=2,ensure_ascii=False)+"\n")
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
receipt=load(BASE+"nonzero-target-validation.json")
assert receipt["status"]=="pass"
assert receipt["source_sha256"]==sha(HERE/"verify_nonzero_target.py")
assert receipt["probe_sha256"]==sha(HERE/"nonzero-target-probe.json")
assert receipt["primary_source_sha256"]==sha(HERE/"nonzero_target_probe.py")
ledger=load("MATH_STATUS.json");known={e["id"]:e for e in ledger["entries"]}
entry=copy.deepcopy(known["KFRAME2"])
entry.update(
    id="KFRAME3",title="Invariant trace poles exclude endpoint limits for a specified nonzero-target collision family",
    scope="For the foundational three-variable Keller map, let a(k)=((1-k)/2,-2,-4(1+3k)/k^2), b(k)=((1+k)/2,-2,-4(1-3k)/k^2), k!=0. Their common target is (0,-2,(1-k^2)/(2k^2)). In the normalized map H_k=(JF(a)L)^-1(F(a+LX)-F(a)), L=[b-a,e_y,e_z], no meromorphic matrix P(k) fixing e1 can give a finite coefficientwise limit of P^-1 H_k(PX) at k=0 or infinity, including after finite local ramification. Both trace JH(2e1) and trace JH(H(2e1)) are polynomial functions of the bounded-degree map coefficients and invariant under these conjugations. The first equals864*k^2-2736*k+723 and has a pole at infinity. The second has leading Laurent term -49460507836416/k^3 at zero. Exact primary symbolic arithmetic and an independent standard-library rational Laurent calculation agree on its entire Laurent polynomial. The pair family has105 generic nonlinear occurrences in its natural frame; no support minimum for its finite parameters is asserted. Other pair families, simultaneous variation of further parameters and finite-parameter sparsity remain open. The normalized raw-degree-seven support interval stays [10,48]; this is a closed degeneration route, not a degree improvement or global optimum.",
    canonical_source=BASE+"NONZERO_TARGET_LIMITS.md",
    dependencies=["KFRAME1"],checker=BASE+"verify_nonzero_target.py",
    proof_type="exact_symbolic",artifact_hash="sha256:"+sha(HERE/"verify_nonzero_target.py"),
    software_lock=[BASE+"requirements.txt",BASE+"nonzero_target_probe.py",
                   BASE+"nonzero-target-probe.json",BASE+"nonzero-target-validation.json"],
    closes_problems=[],narrows_problems=[],supersedes=[],supersedes_notes=[],
    consumers=[],invalidates_assumptions=[],replaced_by=[],independent_replay=True,
    formal_verification=False,external_review=False)
if "KFRAME3" in known:
    known["KFRAME3"].clear();known["KFRAME3"].update(entry)
else:ledger["entries"].append(entry)
save("MATH_STATUS.json",ledger)
lessons=load("knowledge/lessons.json")
record=dict(
    id="METHOD-KELLER-FRAME-LIMIT-TRACE-20261114",area="geometry",
    title="Test trace invariants on iterated collision-axis points before seeking a singular frame",
    when="A normalized bounded-degree collision family is conjugated by a parameter-dependent linear matrix fixing e1.",
    use="Evaluate trace JH(t*e1) and trace JH(H(t*e1)) at a fixed rational t. These are polynomial functions of map coefficients and invariant under every such conjugation. A pole at the proposed endpoint excludes any finite coefficientwise frame limit, even if individual characteristic polynomials along the straight collision axis stay finite. Check full Laurent polynomials independently from the original map and derivatives.",
    avoid="Treating failed diagonal scaling as an obstruction to all frames; assuming finite trace and characteristic coefficients on the straight axis suffice to make an arbitrary frame possible; revisiting a route after an invariant pole has closed it.",
    boundary="KFRAME3 excludes both endpoints of one specified nonzero-target pair family: trace JH(2e1) has leading864*k^2 at infinity and trace JH(H(2e1)) has leading -49460507836416/k^3 at zero. Natural and canonical-axis diagonal tests are retained as discovery evidence. This does not exclude other pair families or finite-parameter support drops and does not alter [10,48].",
    reopen_when="The collision pair family changes, extra parameters co-vary, or the permitted transformations go beyond linear conjugation preserving the normalized collision axis.",
    claims=["KFRAME3"],sources=[BASE+"NONZERO_TARGET_LIMITS.md",BASE+"NEXT_FAMILY.md"],
    implementation=[BASE+"nonzero_target_probe.py",BASE+"verify_nonzero_target.py"],
    keywords=["Keller","collision","singular limit","trace","Jacobian","conjugation","Laurent"],
    reviewed_claims={"KFRAME3":claim_fingerprint(entry)},programme_status="archived")
existing=next((r for r in lessons["lessons"] if r["id"]==record["id"]),None)
if existing is None:lessons["lessons"].append(record)
else:existing.clear();existing.update(record)
save("knowledge/lessons.json",lessons)
work=load("knowledge/work_items.json")
for item in work["actions"]:
    if item["id"]=="WORK-GLDC":
        item["next_step"]="For OP-GLDC, seek size-ten candidates passing coordinate gates in a flag strictly below (7,6,4), or close the next complete support stratum. Separately improve the broader upper bound48 outside KFRAME2's two completed frame classes: transverse x-components, joint boundary-position and full-GL2 variation, other interior pairs, nonlinear coordinates or another Keller family. Singular orbit limits may lower the degree flag, but KFRAME3 excludes both endpoints of the specified one-parameter nonzero-target pair family under every meromorphic frame. Test invariant traces at axis images before searching frames in a new degeneration."
        closure=" KFRAME3 closes only the specified nonzero-target endpoint family; do not revisit its diagonal or general meromorphic frames without changing the pair family."
        if closure not in item["prerequisites"]:item["prerequisites"]+=closure
        item["claims"]=list(dict.fromkeys(item["claims"]+["KFRAME3"]))
        item["sources"]=list(dict.fromkeys(item["sources"]+[BASE+"NONZERO_TARGET_LIMITS.md",BASE+"NEXT_FAMILY.md"]))
        item.setdefault("reviewed_claims",{})["KFRAME3"]=claim_fingerprint(entry)
save("knowledge/work_items.json",work)
print("Registered KFRAME3: both specified affine endpoint routes closed; support interval [10,48] unchanged.")

