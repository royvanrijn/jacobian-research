#!/usr/bin/env python3
"""Register the exact single-new-cube obstruction and bounded search boundaries."""
from pathlib import Path
import copy,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent;BASE="jacobian-followup/cubic-reduction/"
sys.path.insert(0,str(ROOT/"scripts"))
from research import claim_fingerprint
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads((ROOT/p).read_text())
def save(p,d):(ROOT/p).write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n")
v=load(BASE+"power-single-compress/validation.json")
assert v["status"]=="pass" and v["unit_identities"]==86 and v["coordinate_support_cases"]==15
assert v["source_sha256"]==sha(HERE/"verify_power_single_compress.py")
assert v["primary_checkpoint_sha256"]==sha(HERE/"power-single-compress/checkpoint.json")
assert v["certification_sha256"]==sha(HERE/"power-single-compress/certification.json")
ledger=load("MATH_STATUS.json");known={e["id"]:e for e in ledger["entries"]};identifier="KCUB7"
assert identifier not in known or known[identifier]["canonical_source"]==BASE+"SINGLE_CUBE_OBSTRUCTION.md"
entry=copy.deepcopy(known["KCUB6"])
entry.update(id=identifier,title="One arbitrary new cube cannot shorten the saved89-form identity",
 scope="For the exact KCUB5 padded18-output companion identity on its selected89 independent linear cubes, no shorter representation can retain a subset of those forms, add just one arbitrary new linear cube and re-solve every retained coefficient column. This holds over all characteristic-zero extension fields. The new cube must lie in the saved89-cube span; every deleted coefficient column must be proportional to its new coefficient column. There are134 eligible proportional pairs. A support graph has15 maximal coordinate cliques covering225 nonempty cliques. Cubic-span membership equations and coordinate polynomials on these supports reduce the remaining possible cancellations to86 polynomial systems. Each has an explicit rational unit-ideal identity, with multiplier degrees at most3. A separate standard-library Fraction checker reconstructs the89-cube rank, membership row spaces and coordinate maps, covers all2^18 source-coordinate subsets and all coefficient pairs, and multiplies out all86 identities. An explicit four-cube-to-one-cube positive control is also verified. The earlier3916-pair-line traversal and independent polynomial-gcd exclusions are retained as a weaker discovery route. A larger256-form quadratic-pencil pool has cube-span rank154 but a bounded10554-exchange run stays at89; a nonmonotone ordinary-cubic beam and coordinate-translation line tests stay at36 terms, without a global optimum claim. This theorem is an obstruction to one new cube relative to one saved identity, not a minimum89 for the tensor or for cubic-linear Keller maps. At least two new forms, intermediate support growth, other companions and unrestricted smaller constructions remain open; the fixed-tensor interval stays[34,89] and the other simplicity targets are unchanged.",
 canonical_source=BASE+"SINGLE_CUBE_OBSTRUCTION.md",dependencies=["KCUB5"],
 checker=BASE+"verify_power_single_compress.py",artifact_hash="sha256:"+sha(HERE/"verify_power_single_compress.py"),
 software_lock=[BASE+n for n in ["verify.py","verify_power89.py","power89/validation.json","power-basis-anneal/candidate.json",
 "power_single_compress.py","power_line_compress.py","certify_power_single_compress.py",
 "power-single-compress/checkpoint.json","power-single-compress/certification.json","power-single-compress/validation.json",
 "requirements.txt","graph-requirements.txt"]]+[BASE+"power-single-compress/"+name for i in range(15) for name in ["case-%02d.json"%i,"cert-%02d.json"%i,"cert-%02d.sing"%i]],
 supersedes_notes=["The unrestricted single-new-form obstruction includes and strengthens the retained3916-pair-line exclusion. Neither is a global decomposition minimum."],
 supersedes=[],replaced_by=[],closes_problems=[],narrows_problems=[],consumers=[],invalidates_assumptions=[])
if identifier in known:known[identifier].clear();known[identifier].update(entry)
else:ledger["entries"].append(entry);known[identifier]=entry
save("MATH_STATUS.json",ledger)
lessons=load("knowledge/lessons.json")
record=dict(id="METHOD-KELLER-ONE-CUBE-COMPRESSION-20261114",area="geometry",title="Use cube-span membership and unit identities to close arbitrary one-form compression",
 when="A simultaneous cubic identity has independent selected cubes and a smaller decomposition is sought by retaining a subset and adding one form.",
 use="If one new coefficient column is nonzero, its cube lies in the old cube span. Deleted coefficient columns must be proportional. Cover possible source supports by maximal cliques of the coordinate co-occurrence graph, build exact cubic membership and coordinate equations, normalize one nonzero deleted coordinate to1 over an algebraic closure, and export unit-ideal multipliers for every proportional pair. Independently reconstruct the membership row spaces and clique coverage, then multiply out each identity with Fraction arithmetic.",
 avoid="Equating independence with minimality; confusing no single new form with no smaller decomposition; taking Singular syntax errors for proper or unit ideals; sampling three generic values to exclude roots of a cubic denominator; treating primary Groebner status as independent proof.",
 boundary="KCUB7 excludes one arbitrary new cube relative to the saved89 identity over all characteristic-zero extensions, via15 coordinate supports,134 proportional pairs and86 exact unit identities of multiplier degree at most3. It does not exclude two new forms or all smaller tensors/maps. The independent pair-line replay covers3916 lines,178 fixed points and48 full-binary pairs (eight distinct lines); the unrestricted certificate is stronger. A first three-sample gap and an initial rational-expression Singular parse failure are preserved as inconclusive historical attempts, followed by corrected exact replays.",
 reopen_when="Use at least two new forms at once, classify the complete Veronese intersection of the saved cube span, permit intermediate support growth, or change the companion beyond its existing linear orbit. Preserve the current89 witness and require full matrix and collision verification for any improvement.",
 claims=[identifier,"KCUB5"],sources=[BASE+"SINGLE_CUBE_OBSTRUCTION.md",BASE+"EXPLORATION.md"],
 implementation=[BASE+n for n in ["power_line_compress.py","verify_power_line_compress.py","power_single_compress.py","certify_power_single_compress.py","verify_power_single_compress.py"]],
 keywords=["Keller","shared cubes","one new cube","Veronese","coordinate cliques","unit ideal","Nullstellensatz","89 forms"],programme_status="archived",reviewed_claims={})
for key in record["claims"]:record["reviewed_claims"][key]=claim_fingerprint(known[key])
old=next((r for r in lessons["lessons"] if r["id"]==record["id"]),None)
if old is None:lessons["lessons"].append(record)
else:old.clear();old.update(record)
for r in lessons["lessons"]:
 if r["id"]=="METHOD-KELLER-SHARED-CUBE-BASIS-EXCHANGE-20261114":
  addition=" Completing squares of quadratic-output pencils adds a256-form pool with rank154, genuinely larger than the old rank121. A60-second10554-exchange search remains at89. KCUB7 separately excludes any shortening with only one new form, even outside finite pools; simultaneous replacements remain open."
  if addition not in r["boundary"]:r["boundary"]+=addition
  r["reopen_when"]="Test at least two new forms simultaneously, classify the entire cube/Veronese intersection of the saved span, or change the companion beyond its linear orbit. Do not repeat terminal178,858 or256-form search caps, or propose a single-new-form compression closed by KCUB7."
  r["claims"]=list(dict.fromkeys(r["claims"]+[identifier]));r["sources"]=list(dict.fromkeys(r["sources"]+[BASE+"SINGLE_CUBE_OBSTRUCTION.md"]))
  r["implementation"]=list(dict.fromkeys(r["implementation"]+[BASE+"power_basis_pencils.py"]))
  r.setdefault("reviewed_claims",{})[identifier]=claim_fingerprint(entry)
affine=dict(id="METHOD-KELLER-CUBIC-NONMONOTONE-AFFINE-SEARCH-20261114",area="geometry",title="Distinguish nonmonotone affine searches from the established36-term construction",
 when="The36-term12-variable cubic is a local rational elementary-line minimum and an affine coordinate improvement is sought.",
 use="Retain the36-term integer input and full conjugation histories. Allow bounded temporary support increases, cancellation roots and generic rational values. For a translation centre p, normalize the target by JF(p)^-1 so the origin derivative staysI; compare the complete quadratic and cubic coefficient support.",
 avoid="Calling a bounded beam an affine-orbit optimum; inferring arbitrary-translation coverage from coordinate lines; counting a lower monomial presentation of a linearly conjugate companion as a new shared-cube mechanism.",
 boundary="The seeded width12 depth8 beam checks69747 elementary-line values with a support cap four above the current best and a coefficient-height bound, encounters1329 distinct coefficient maps, and stays at36. The primary SymPy translation computation examines all rational cancellation values on twelve coordinate lines after exact inverse-derivative normalization and stays at36. Neither is an independently replayed global support theorem.",
 reopen_when="Use simultaneous translation coordinates, another affine frame, or a nonlinear transformation with a verified degree-three final map. Do not rerun the same seeded beam and coordinate-line family as a new search.",
 claims=["KCUB3"],sources=[BASE+"SINGLE_CUBE_OBSTRUCTION.md",BASE+"EXPLORATION.md"],implementation=[BASE+"linear_beam12.py",BASE+"translation_lines12.py"],
 keywords=["Keller","36 terms","ordinary cubic","beam","translation","tangent normalization"],programme_status="archived",reviewed_claims={"KCUB3":claim_fingerprint(known["KCUB3"])})
old=next((r for r in lessons["lessons"] if r["id"]==affine["id"]),None)
if old is None:lessons["lessons"].append(affine)
else:old.clear();old.update(affine)
save("knowledge/lessons.json",lessons)
work=load("knowledge/work_items.json")
for r in work["actions"]:
 if r["id"]=="WORK-GLDC":
  addition=" KCUB7 now excludes shortening the saved89 identity with only one arbitrary new cube, with retained coefficients freely re-solved;86 unit identities cover15 source supports. Next use at least two new forms at once or classify all cubes in the saved span. The changed256-form quadratic-pencil pool increases rank to154 but stays at89 after10554 exchanges. The ordinary-cubic nonmonotone beam and coordinate-translation lines stay at36; simultaneous translations or nonlinear degree-preserving changes remain untested. These do not close global dimension, term or strict three-variable-profile minima."
  if addition not in r["next_step"]:r["next_step"]+=addition
  r["claims"]=list(dict.fromkeys(r["claims"]+[identifier]));r["sources"]=list(dict.fromkeys(r["sources"]+[BASE+"SINGLE_CUBE_OBSTRUCTION.md"]))
  r.setdefault("reviewed_claims",{})[identifier]=claim_fingerprint(entry)
save("knowledge/work_items.json",work)
print("Registered KCUB7: arbitrary single-new-cube compression excluded for the saved89 identity; full goal remains open.")
