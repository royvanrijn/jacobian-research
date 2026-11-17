#!/usr/bin/env python3
"""Register fixed scalar-cube-span optimality without claiming global minimality."""
from pathlib import Path
import copy,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent;BASE="jacobian-followup/cubic-reduction/"
sys.path.insert(0,str(ROOT/"scripts"))
from research import claim_fingerprint
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads((ROOT/p).read_text())
def save(p,d):(ROOT/p).write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n")
v=load(BASE+"power-span-split/validation.json");pair=load(BASE+"power-two-exchange/validation.json")
assert v["status"]=="pass" and v["source_sha256"]==sha(HERE/"verify_power_span_optimality.py")
assert v["classification_sha256"]==sha(HERE/"power-span-split/classification.json")
assert v["primary_checkpoint_sha256"]==sha(HERE/"power-span-split/checkpoint.json")
assert pair["status"]=="pass" and pair["source_sha256"]==sha(HERE/"verify_power_two_exchange.py")
assert pair["fixed_span_validation_sha256"]==sha(HERE/"power-span-split/validation.json")
ledger=load("MATH_STATUS.json");known={e["id"]:e for e in ledger["entries"]};identifier="KCUB8"
assert identifier not in known or known[identifier]["canonical_source"]==BASE+"SPAN_OPTIMALITY.md"
entry=copy.deepcopy(known["KCUB7"])
entry.update(id=identifier,title="Saved89-dimensional scalar cubic span forces89 forms",
 scope="Let P be KCUB5's fixed padded18-output homogeneous cubic tensor and S the scalar cubic span of its89 selected independent cubes. Every shared-cube representation of P using only cubes in S spans all of S and therefore uses at least89 forms, over every characteristic-zero extension. This permits every form and output coefficient to change. The projective cube locus in S consists exactly of65 isolated selected points and eight rational projective lines: span(z4,tau), span(z2,z4+tau), span(z2,z4-tau), and span(Wi,tau) for five Wi. Rational linear-factor splitting trees give408 nodes and236 checked factor identities across15 coordinate supports. A separate Fraction checker reconstructs membership row spaces, all2^18 source-support coverage, every factor branch and leaf, the complete73-component classification, and all rank premises. The65 isolated cubes are compulsory. An annihilator of a hypothetical proper representation span vanishes on the five full Wi binary cube spaces; the remaining annihilator conditions have rank87 and leave two explicit functionals. Their zero cubes project to a special space of dimension at most5, while the projected target has rank7, a contradiction. This is a fixed-span minimum, not global rank89 for P, a border-rank bound or a minimum cubic-linear Keller dimension. A shorter two-new-form exchange would require proportional new output columns and two nonzero residues outside S that cancel. In the930-form union pool, independent exact replay checks1774 such pairs and cancels at most two old forms, so no shortening occurs there; arbitrary outside forms remain open. The unrestricted fixed-tensor shared-cube interval remains[34,89], ordinary cubic12/36 and normalized three-variable[10,48] are unchanged, and all global simplicity minima remain open.",
 canonical_source=BASE+"SPAN_OPTIMALITY.md",dependencies=["KCUB5"],
 checker=BASE+"verify_power_span_optimality.py",artifact_hash="sha256:"+sha(HERE/"verify_power_span_optimality.py"),
 software_lock=[BASE+n for n in ["verify.py","verify_power_single_compress.py","verify_power89.py","power89/validation.json",
 "power-basis-anneal/candidate.json","power_span_split.py","power-span-split/checkpoint.json",
 "power-span-split/classification.json","power-span-split/validation.json","power_two_exchange.py","verify_power_two_exchange.py",
 "power-two-exchange/pool.json","power-two-exchange/checkpoint.json","power-two-exchange/validation.json",
 "power_search.py","power_line_compress.py","power-basis-pencils/pool.json","requirements.txt","graph-requirements.txt"]]
 +[BASE+"power-span-split/case-%02d.json"%i for i in range(15)]
 +[BASE+"power-single-compress/case-%02d.json"%i for i in range(15)],
 supersedes_notes=["The fixed-span theorem strengthens KCUB7's single-new-form obstruction and closes arbitrary simultaneous changes that keep every scalar cube in S. It does not close changes of scalar cubic span."],
 supersedes=[],replaced_by=[],closes_problems=[],narrows_problems=[],consumers=[],invalidates_assumptions=[])
if identifier in known:known[identifier].clear();known[identifier].update(entry)
else:ledger["entries"].append(entry);known[identifier]=entry
save("MATH_STATUS.json",ledger)
lessons=load("knowledge/lessons.json")
record=dict(id="METHOD-KELLER-CUBE-SPAN-CLASSIFICATION-20261114",area="geometry",title="Classify the cube locus before searching inside a fixed scalar span",
 when="A shared-cube Keller lift is to be compressed and one-new-form edits have already been excluded.",
 use="Build the coordinate co-occurrence graph and exact cube-span membership equations. Classify their zero sets by rational linear-factor splitting trees: verify each selected constraint combination, linear factor product, hyperplane kernel and terminal subspace independently. Merge leaves and identify compulsory isolated cubes. Use a linear annihilator and exact restriction polynomials on the remaining binary families to test whether a smaller representation can stay in the old scalar span. Then project a changed candidate pool to the quotient by that span; group proportional nonzero residues to enumerate exact two-new-cube exchanges.",
 avoid="Treating the saved span's minimum as the tensor's unrestricted rank or a Keller dimension minimum; proposing further within-span basis search after it is closed; trusting primary decomposition without a replayable cover; confusing pair-residue enumeration in one pool with arbitrary two-form coverage.",
 boundary="KCUB8 proves89 minimal among all representations whose cubes lie in the saved89-dimensional scalar span. The cube locus is65 isolated points plus eight projective lines, checked through408 tree nodes and236 factor identities. The annihilator rank87 leaves two functionals; their special zero-cube span has dimension at most5 against target rank7. The930-form union pool has rank154,151 inside-span cubes and445 nonzero residue directions. Its1774 admissible outside pairs cancel at most two old forms. Independent Fraction replay also finds209 rank3 proportional-column raw triples and five tangent binary monomials, so direct raw three-to-two replacement does not work. None proves global rank or dimension89.",
 reopen_when="Seek two new forms outside S with proportional output columns, their scalar residue cancellation lying in S, and at least three prescribed cancellation coordinates from one old coefficient direction; allow retained coefficient columns to change. Alternatively use at least three new forms or a different companion. Do not repeat the closed fixed-span or930-form two-exchange searches.",
 claims=[identifier,"KCUB5"],sources=[BASE+"SPAN_OPTIMALITY.md",BASE+"EXPLORATION.md"],
 implementation=[BASE+n for n in ["power_span_components.py","power_span_split.py","verify_power_span_optimality.py","power_two_exchange.py","verify_power_two_exchange.py"]],
 keywords=["Keller","cube span","Veronese","89 forms","factor splitting","annihilator","two new cubes","quotient residues"],programme_status="archived",reviewed_claims={})
for k in record["claims"]:record["reviewed_claims"][k]=claim_fingerprint(known[k])
old=next((r for r in lessons["lessons"] if r["id"]==record["id"]),None)
if old is None:lessons["lessons"].append(record)
else:old.clear();old.update(record)
for r in lessons["lessons"]:
 if r["id"] in ["METHOD-KELLER-ONE-CUBE-COMPRESSION-20261114","METHOD-KELLER-SHARED-CUBE-BASIS-EXCHANGE-20261114"]:
  addition=" KCUB8 completes the saved-span cube locus (65 isolated points and eight lines) and proves that every decomposition staying in that scalar span uses at least89 forms. A shorter two-new-form exchange must use proportional output columns and cancel two outside residues. The930-form pool has1774 admissible pairs, none canceling three old forms."
  if addition not in r["boundary"]:r["boundary"]+=addition
  r["reopen_when"]="Change the scalar cube span: seek two outside forms with proportional output columns whose residues cancel and whose combined coordinates delete at least three old forms in one coefficient direction, with retained columns re-solved; or use at least three new forms or a different companion. Fixed-span changes and the930-form two-exchange pool are closed by KCUB8."
  r["claims"]=list(dict.fromkeys(r["claims"]+[identifier]));r["sources"]=list(dict.fromkeys(r["sources"]+[BASE+"SPAN_OPTIMALITY.md"]))
  r.setdefault("reviewed_claims",{})[identifier]=claim_fingerprint(entry)
save("knowledge/lessons.json",lessons)
work=load("knowledge/work_items.json")
for r in work["actions"]:
 if r["id"]=="WORK-GLDC":
  addition=" KCUB8 classifies the saved-span cube locus and proves89 minimal within that89-dimensional scalar cubic span; arbitrary simultaneous inside-span changes are now closed. Any shorter two-new-form edit needs proportional output columns and two outside residues that cancel. The930-form union pool has1774 admissible pairs and no shortening under exact independent replay. Next solve the unrestricted rank-two scalar-cubic residue equations with at least three prescribed canceled coordinates from one old coefficient direction, or enlarge the form pool/change the companion; do not repeat within-span basis changes or the same930-form pair test. Global tensor rank, special-form dimension, ordinary cubic term and strict three-variable-profile minima remain open."
  if addition not in r["next_step"]:r["next_step"]+=addition
  r["claims"]=list(dict.fromkeys(r["claims"]+[identifier]));r["sources"]=list(dict.fromkeys(r["sources"]+[BASE+"SPAN_OPTIMALITY.md"]))
  r.setdefault("reviewed_claims",{})[identifier]=claim_fingerprint(entry)
save("knowledge/work_items.json",work)
print("Registered KCUB8: saved scalar-span minimum89; unrestricted tensor and Keller minima remain open.")
