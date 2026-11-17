"""Register the36-term twelve-variable map and correct the missed size baseline."""
from pathlib import Path
import copy,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
BASE="jacobian/cubic-reduction/"
sys.path.insert(0,str(ROOT/"scripts"))
from research import claim_fingerprint
def load(p):return json.loads((ROOT/p).read_text())
def save(p,d):(ROOT/p).write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n")
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
v=load(BASE+"integer-validation12.json")
assert v["status"]=="pass" and (v["dimension"],v["degree"],v["nonlinear_total"])==(12,3,36)
assert v["source_sha256"]==sha(HERE/"verify_integer12.py") and v["witness_sha256"]==sha(HERE/"cubic12-integer-witness.json")
ledger=load("MATH_STATUS.json");known={e["id"]:e for e in ledger["entries"]}
assert known["BCR2"]["state"]==known["BCRDB1"]["state"]=="proved"
entry=copy.deepcopy(known["KCUB2"]);identifier="KCUB3"
entry.update(
    id=identifier,title="Thirty-six-term integral conjugate of the established twelve-variable cubic Keller map",
    scope="The already independently verified BCRDB1 section4.6 twelve-variable degree-three map has40 nonlinear monomial occurrences in its displayed tangent-normalized coordinates. Conjugating first by I+(3/2)*E_(9,8), then by I+(3/2)*E_(11,10), with one-based indices, reduces the count40 to39 to36. A diagonal conjugation with entries1/2 in positions9 and11 makes all coefficients and both colliding points integral without changing support. The final map has K(0)=0, JK(0)=I, full constant Jacobian determinant one, degree3, and nonlinear row counts (4,10,10,2,1,1,1,1,2,1,1,2). Two distinct integer points map to (0,1,0,...,0). The collision fibre is not the origin fibre; this term count is separate from affine collision-axis normalization. A SymPy construction reconstructs the pinned BCRDB1 source and an independent standard-library rational-polynomial implementation verifies every coefficient and conjugation, the40-term input determinant, the36-term rational and integral determinants in all12 symbolic variables, the origin derivatives, integral coefficients and full integer collision. Each determinant uses an explicitly verified nilpotent9-by9 carrier correction N with N^3=0 and both inverse identities, followed by the complete polynomial Schur determinant. Dimension12 was already established by BCR2 and BCRDB1; this is a displayed sparsity improvement only. Earlier KCUB1 dimension17/count58 and KCUB2 dimension13/count42 examples are dominated by the pre-existing12/count40 baseline and retained as historical experiments, not improvements to known dimension bounds. No minimum dimension, global term optimum, new counterexample mechanism, cubic-homogeneous or cubic-linear improvement, novelty priority, formal verification or external review is claimed. The separate normalized three-variable raw-degree-seven support interval stays [10,48], and the full simplicity goal remains open.",
    canonical_source=BASE+"PROOF12.md",dependencies=["BCRDB1"],
    checker=BASE+"verify_integer12.py",artifact_hash="sha256:"+sha(HERE/"verify_integer12.py"),
    software_lock=[BASE+"requirements.txt",BASE+"construct_sparse12.py",BASE+"construct_integer12.py",
      BASE+"verify_sparse12.py",BASE+"verify.py",BASE+"optimize12.py",
      BASE+"de-bondt12-witness.json",BASE+"de-bondt12-sparse-candidate.json",
      BASE+"cubic12-sparse-witness.json",BASE+"cubic12-integer-witness.json",
      BASE+"sparse-validation12.json",BASE+"integer-validation12.json",
      "jacobian/de-bondt-2610.24889/verify.py"],
    supersedes_notes=["The established12-variable baseline was missed in the KCUB1/KCUB2 exploratory comparison. Those larger maps remain verified but do not improve known bounds. This result improves the displayed count40 to36 at unchanged dimension12."],
    supersedes=[],replaced_by=[],closes_problems=[],narrows_problems=[],consumers=[],invalidates_assumptions=[])
if identifier in known:known[identifier].clear();known[identifier].update(entry)
else:ledger["entries"].append(entry);known[identifier]=entry
known["KCUB1"]["scope"]="A verified rational cubic map in17 variables has determinant one, identity derivative at zero and a three-point rational collision, with58 displayed nonlinear monomial occurrences. It is obtained by removing two redundant carrier directions from the posted19-variable construction; the full coefficient substitution and symbolic Schur determinant are independently verified in rational sparse arithmetic. The quadratic midpoint identity also verifies the classical degree-three collision threshold across varying dimensions. This is a retained historical existence theorem, not a size improvement: BCR2 and BCRDB1 had already established twelve-variable cubic maps before this experiment, and BCRDB1's displayed count40 dominates17/count58. KCUB3 subsequently supplies an integral12-variable36-term conjugate. The missed baseline and superseded optimization role are explicitly recorded in BASELINE_CORRECTION.md. No optimality, new dimension record, priority, formalization or external review is claimed."
known["KCUB2"]["scope"]="A verified rational cubic map in13 variables has determinant one, K(0)=0, JK(0)=I and a three-point rational collision, with42 displayed nonlinear occurrences. Its output-factor-reuse construction and the final S=I-3*E_(7,11) conjugation are replayed independently, and the full13-variable polynomial Schur determinant is one. The construction is retained as a historical method experiment. It is dominated in both dimension and count by BCRDB1's pre-existing12-variable40-term map; the earlier comparison only to KCUB1's17-variable58-term example missed that known baseline and supplied no improvement to known size bounds. KCUB3 gives the current integral12-variable36-term conjugate. The search beams and order trials prove no lower bound or global optimum. The three-variable normalized interval remains [10,48]. No minimum dimension, support optimum, priority, formalization or external review is claimed."
for key in ("KCUB1","KCUB2"):
    known[key]["replaced_by"]=[identifier]
    known[key]["supersedes_notes"]=["Optimization role superseded: the existing BCRDB1 twelve-variable baseline dominates this experiment. KCUB3 gives a36-term integral conjugate. The original exact existence theorem and its certificates remain valid."]
save("MATH_STATUS.json",ledger)
lessons=load("knowledge/lessons.json")
for item in lessons["lessons"]:
    if item["id"] in ("METHOD-KELLER-CUBIC-CARRIER-DESTABILIZATION-20261114","METHOD-KELLER-CUBIC-OUTPUT-FACTOR-REUSE-20261114"):
        item["boundary"]="The exact17/58 and13/42 constructions remain verified, but BCR2 and BCRDB1 already supplied12-variable cubic maps. BCRDB1's count40 dominates both experiments; KCUB3 improves that established map to an integral36-term form. These older construction routes are historical method evidence, not current size improvements. All bounded beam, slice and rational-line limits remain local rather than optimality proofs."
        item["reopen_when"]="A further construction uses fewer than12 variables, improves the36-term target under comparable cubic/tangent normalization, or supplies a genuinely additional property absent from BCRDB1/KCUB3."
        item["claims"]=list(dict.fromkeys(item["claims"]+["BCR2","BCRDB1",identifier]))
        item["sources"]=list(dict.fromkeys(item["sources"]+[BASE+"BASELINE_CORRECTION.md",BASE+"PROOF12.md"]))
records=[
 dict(id="METHOD-KELLER-CUBIC-BASELINE-AUDIT-20261114",area="geometry",
      title="Audit degree-three and economic-stabilization baselines before claiming a smaller cubic",
      when="A new cubic Keller construction is compared with existing dimensions or support counts.",
      use="Read BCR2, BCRDB1 and BCRDBGRAPH1 and their canonical sources directly, including number-word variants such as twelve-variable and degree-three. A failed central search is not an absence result. Compare ordinary degree, cubic homogeneity, unipotence, cubic-linearity and collision/tangent normalizations separately.",
      avoid="Using the posted19-variable example as the baseline after12-variable examples were already verified in the same repository. The17/58 and13/42 constructions were exact but dominated, so their validation did not establish optimization progress.",
      boundary="The mistaken dimension framing is corrected in the live ledger and source notes. The retained experiments are not invalidated algebraically. KCUB3 advances the established12-variable map from40 to36 terms; dimension12 is unchanged.",
      reopen_when="A new literature or repository result changes the relevant baseline, or a different special form or normalization is imposed.",
      claims=["BCR2","BCRDB1","BCRDBGRAPH1","KCUB1","KCUB2",identifier],
      sources=[BASE+"BASELINE_CORRECTION.md","jacobian/de-bondt-2610.24889/README.md",BASE+"PROOF12.md"],
      implementation=[BASE+"register12.py"],
      keywords=["Keller","baseline","degree-three","economic stabilization","dominated","retrieval"]),
 dict(id="METHOD-KELLER-CUBIC-GROUPED-SHEAR-SPARSITY-20261114",area="geometry",
      title="Align repeated carrier combinations and clear denominators without changing cubic support",
      when="A known tangent-normalized cubic Keller map has repeated linear carrier combinations.",
      use="Conjugate by rational elementary shears that align the combinations, recount all coordinates including carrier outputs, and use diagonal conjugation to clear denominators. Independently verify complete coefficients, integrality, transported collisions and the full symbolic determinant.",
      avoid="Claiming a new dimension bound when only the coordinate support changed; treating a rational elementary-line minimum as full affine-orbit optimality; confusing the displayed collision fibre with an origin-axis collision normalization.",
      boundary="KCUB3 uses two3/2 shears to reduce BCRDB1's count40 to39 to36, then scales coordinates9 and11 by1/2 in the source and inversely in the target. The final map and collision points are integral. The current132-line rational descent finds no further single-shear decrease, but joint, algebraic and nonlinear changes remain open.",
      reopen_when="Joint coordinate changes, another basis, a nonliteral graph deletion, sequential target completion or another cubic family is available.",
      claims=["BCRDB1","BCRDBGRAPH1",identifier],
      sources=[BASE+"PROOF12.md",BASE+"BASELINE_CORRECTION.md"],
      implementation=[BASE+"optimize12.py",BASE+"construct_sparse12.py",BASE+"construct_integer12.py",
                      BASE+"verify_sparse12.py",BASE+"verify_integer12.py"],
      keywords=["Keller","36 terms","twelve variables","integral","conjugation","shear","sparsity"])]
for record in records:
    record["programme_status"]="archived";record["reviewed_claims"]={}
    old=next((r for r in lessons["lessons"] if r["id"]==record["id"]),None)
    if old is None:lessons["lessons"].append(record)
    else:old.clear();old.update(record)
for item in lessons["lessons"]:
    for key in ("KCUB1","KCUB2",identifier,"BCR2","BCRDB1","BCRDBGRAPH1"):
        if key in item.get("claims",[]):item.setdefault("reviewed_claims",{})[key]=claim_fingerprint(known[key])
save("knowledge/lessons.json",lessons)
work=load("knowledge/work_items.json")
for item in work["actions"]:
    if item["id"]=="WORK-GLDC":
        marker="Across varying dimensions,"
        assert marker in item["next_step"]
        item["next_step"]=item["next_step"][:item["next_step"].index(marker)]+"Across varying dimensions, BCR2 and BCRDB1 already establish cubic maps in12 variables; KCUB3 reduces BCRDB1's displayed count40 to an integral36-term form at unchanged dimension. Seek fewer than12 variables, fewer than36 cubic terms under comparable normalization, or a better special-form construction. Use nonliteral coordinates or sequential nonlinear stages beyond BCRDBGRAPH1's completed nine literal deletion obstructions. KCUB1/KCUB2 are dominated historical experiments and must not be treated as current size baselines."
        item["claims"]=list(dict.fromkeys(item["claims"]+["BCR2","BCRDB1","BCRDBGRAPH1",identifier]))
        item["sources"]=list(dict.fromkeys(item["sources"]+[BASE+"PROOF12.md",BASE+"BASELINE_CORRECTION.md",
             "jacobian/de-bondt-2610.24889/README.md"]))
    for key in ("KCUB1","KCUB2",identifier,"BCR2","BCRDB1","BCRDBGRAPH1"):
        if key in item.get("claims",[]):item.setdefault("reviewed_claims",{})[key]=claim_fingerprint(known[key])
save("knowledge/work_items.json",work)
print("Registered KCUB3: integral12-variable36-term cubic. Corrected dominated KCUB1/KCUB2 baselines. Full goal remains open.")

