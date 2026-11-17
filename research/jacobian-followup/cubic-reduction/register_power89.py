"""Register the independently checked89-variable integral cubic-linear witness."""
from pathlib import Path
import copy,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
BASE="jacobian-followup/cubic-reduction/"
sys.path.insert(0,str(ROOT/"scripts"))
from research import claim_fingerprint
def load(p):return json.loads((ROOT/p).read_text())
def save(p,d):(ROOT/p).write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n")
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
v=load(BASE+"power89/validation.json")
assert v["status"]=="pass" and v["dimension"]==89 and v["coefficient_ring"]=="Z"
assert v["source_sha256"]==sha(HERE/"verify_power89.py") and v["witness_sha256"]==sha(HERE/"power89/witness.json")
assert v["primary_source_sha256"]==sha(HERE/"build_power89.py")
ledger=load("MATH_STATUS.json");known={e["id"]:e for e in ledger["entries"]}
entry=copy.deepcopy(known["KCUB3"]);identifier="KCUB5"
entry.update(id=identifier,title="Eighty-nine-variable integral cubic-linear Keller collision",
    scope="An explicit integer89-by89 matrix A defines F(Y)=Y+(A*Y)^{circ3} with identity linear part, determinant identically one and two explicit distinct rational points sharing an image. This improves the audited BCRDB1 ninety-variable cubic-linear construction to89, without a minimum-dimension or priority claim. A common set of89 integer linear forms gives an exact18-output cube identity for the existing homogeneous companion plus the padding cube tau^3. The rational coefficient matrix has denominator LCM24 in its first17 rows, and its last row selects the pure tau form. Scaling the first17 rows by24, zeroing that form's row in L to obtain L0, and setting A=L0*B gives the integer matrix. The full projection identity is B F=S(G24(BY)), including the necessary elementary tau^3 target shear S; it is not the unsheared companion identity. The source12-variable determinant, the companion Schur identity and Sylvester's determinant identity prove det DF=1. A rational right-inverse certificate for B and explicit fibre correction produce and directly verify the complete89-coordinate collision. A separate standard-library Fraction implementation checks all18 cube identities, all18 projection identities, integer scaling and A factorization, right inverse, base determinant and both distinct points; it uses no SymPy or FLINT. The89-variable determinant is proved by exact factorization, not expanded. The fixed178-form candidate space spans dimension121. The initial support89 SMT query timed out UNKNOWN, a bounded plateau basis search stayed at90, and an exact basis-exchange search allowing temporary support growth to94 found89. These searches are not lower bounds or global optimality evidence. No new counterexample mechanism, formal proof assistant, external review, minimum cubic-linear dimension, general cubic dimension improvement or normalized three-variable support change is claimed.",
    canonical_source=BASE+"POWER89.md",dependencies=["BCRDB1"],
    checker=BASE+"verify_power89.py",artifact_hash="sha256:"+sha(HERE/"verify_power89.py"),
    software_lock=[BASE+"requirements.txt",BASE+"graph-requirements.txt",BASE+"verify.py",BASE+"verify_sparse12.py",
       BASE+"build_power89.py",BASE+"de-bondt12-witness.json",BASE+"power_basis_anneal.py",BASE+"power_basis.py",
       BASE+"power_search.py",BASE+"power-basis-anneal/candidate.json",BASE+"power-basis-anneal/checkpoint.json",
       BASE+"power89/witness.json",BASE+"power89/validation.json",
       "jacobian-followup/de-bondt-2610.24889/verify.py","jacobian-followup/de-bondt-2610.24889/polarize.py",
       "jacobian-followup/de-bondt-2610.24889/polarization.json"],
    supersedes_notes=["The audited cubic-linear construction bound90 is improved to89 by a different shared-cube basis. Independence and nondeletability of the original fixed90-form list remain valid."],
    supersedes=[],replaced_by=[],closes_problems=[],narrows_problems=[],consumers=[],invalidates_assumptions=[])
if identifier in known:known[identifier].clear();known[identifier].update(entry)
else:ledger["entries"].append(entry);known[identifier]=entry
save("MATH_STATUS.json",ledger)
lessons=load("knowledge/lessons.json")
record=dict(id="METHOD-KELLER-SHARED-CUBE-BASIS-EXCHANGE-20261114",area="geometry",
    title="Optimize the common cube basis jointly and allow temporary support growth",
    when="A cubic-linear lift is built from simultaneous power decompositions of a homogeneous companion.",
    use="Form a finite pool of alternate linear cubes, express all target rows in a full cube basis and pivot exact rational basis coordinates jointly. Score nonzero coefficient rows across every output. Allow bounded temporary support growth; verify the final cube identity independently rather than trusting the heuristic score. Retain the pure padding cube and its isolated last coefficient row, scale denominators, and preserve the projection target shear when building the integral lift.",
    avoid="Treating independence of one90-form list as a lower bound on all decompositions; restricting optimization to one polarization choice per term; reporting a timed-out SMT query as UNSAT; dropping the tau^3 target shear; replacing the determinant proof by linear-matrix nilpotence.",
    boundary="KCUB5 verifies an integer89-variable cubic-linear map and complete rational collision. The fixed178-form pool has rank121. A60-second support query is UNKNOWN, a1280-step plateau search stays at90, and a successful exact basis search crosses support94 before reaching89. None proves finite-pool or global optimality. The determinant uses checked cube identities, the12-variable base determinant, companion block elimination and Sylvester; no89-variable determinant expansion is claimed.",
    reopen_when="Start from the saved89-form decomposition and seek fewer forms, enlarge the linear-form pool, change the tensor or improve the underlying homogeneous companion. Preserve the exact projection and collision proof for every improvement.",
    claims=["BCRDB1",identifier],sources=[BASE+"POWER89.md",BASE+"EXPLORATION.md"],
    implementation=[BASE+"power_search.py",BASE+"power_basis.py",BASE+"power_basis_anneal.py",BASE+"build_power89.py",BASE+"verify_power89.py"],
    keywords=["Keller","Druzkowski","89 variables","shared cubes","basis exchange","temporary support","Sylvester","padding shear"],
    programme_status="archived",reviewed_claims={})
old=next((r for r in lessons["lessons"] if r["id"]==record["id"]),None)
if old is None:lessons["lessons"].append(record)
else:old.clear();old.update(record)
for item in lessons["lessons"]:
    if item["id"]=="METHOD-DB90-DIRECT-PAIRING-20261113":
        addition=" KCUB5 subsequently gives a verified89-variable integral cubic-linear lift using a different joint cube basis; the original90-form independence result is not a global lower bound."
        if addition not in item["boundary"]:item["boundary"]+=addition
        item["claims"]=list(dict.fromkeys(item["claims"]+[identifier]))
        item["sources"]=list(dict.fromkeys(item["sources"]+[BASE+"POWER89.md"]))
    for key in ["BCRDB1",identifier]:
        if key in item.get("claims",[]):item.setdefault("reviewed_claims",{})[key]=claim_fingerprint(known[key])
save("knowledge/lessons.json",lessons)
work=load("knowledge/work_items.json")
for item in work["actions"]:
    if item["id"]=="WORK-GLDC":
        addition=" For cubic-linear form, KCUB5 improves the audited90-variable map to an integer89-variable witness. Continue from the saved89-form decomposition to seek fewer forms or enlarge the finite pool; do not restart from the dominated90-form baseline or infer optimality from the SMT timeout. Ordinary cubic12/36 and three-variable normalized [10,48] remain separate targets."
        if addition not in item["next_step"]:item["next_step"]+=addition
        item["claims"]=list(dict.fromkeys(item["claims"]+[identifier]))
        item["sources"]=list(dict.fromkeys(item["sources"]+[BASE+"POWER89.md"]))
    if identifier in item.get("claims",[]):item.setdefault("reviewed_claims",{})[identifier]=claim_fingerprint(entry)
save("knowledge/work_items.json",work)
print("Registered KCUB5: integral89-variable cubic-linear Keller witness; minimum dimensions and full simplicity goal remain open.")

