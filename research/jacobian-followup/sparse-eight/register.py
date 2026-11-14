"""Register the eight-term exclusion and the separate bounded pairing theorem."""
import copy
import hashlib
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
BASE="jacobian-followup/sparse-eight/"
sys.path.insert(0,str(ROOT/"scripts"))
from research import claim_fingerprint
from research_partial_reviews import entry_fingerprint
def load(name):
    return json.loads((ROOT/name).read_text())
def save(name,data):
    (ROOT/name).write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n")
ledger=load("MATH_STATUS.json")
known={e["id"]:e for e in ledger["entries"]}
template=copy.deepcopy(known["GLDC2"])
entry=copy.deepcopy(template)
entry.update(
    id="GLDC3",title="Complete eight-term degree-seven Keller collision exclusion",
    scope="Over every characteristic-zero field, every three-variable Keller map in the collision frame F(0)=F(e_1)=0 and JF(0)=I, with raw coordinate degrees at most seven, has at least nine nonlinear monomial occurrences counted with coordinate multiplicity. Complete closure of all 90 minimal-axis-exponent shards leaves 20857 oriented size-eight supports, or 17482 residual-swap orbits. All 17482 size-eight coefficient-torus ideals and 2326 smaller regression ideals are (1) over Q in both Singular and SymPy. A separate global traversal with independently compiled monomials and Leibniz determinant contributions exactly reproduces all 23183 oriented supports through size eight after 1947529 visits. An independent derivative-level compiler verifies every full determinant, collision and exact-support saturation ideal; both rational backends are rerun. The coefficient proof uses no plane theorem. A separate structural classification finds 17435 size-eight univariate-coordinate orbits, 46 elementary-shear-coordinate orbits excluded by the degree-controlled plane reduction, and one hyperplane-normal-factor obstruction F2=y*(1+d*x). Thus eight-term attainment in the normalized raw-degree-at-most-seven range is impossible. Nine-term attainment, global degree minima, and unrestricted-degree support minima remain unknown. No entire dense profile above maximum degree two is excluded; no formal proof-assistant certificate, external review or novelty priority is claimed.",
    canonical_source=BASE+"PROOF.md",dependencies=["GLDC2","CAUF1"],
    checker=BASE+"verify.py",artifact_hash="sha256:"+hashlib.sha256((HERE/"verify.py").read_bytes()).hexdigest(),
    software_lock=[BASE+"requirements.txt",BASE+"validation.json",BASE+"results/manifest.json",
                   BASE+"verify_coverage.py",BASE+"verify_structure.py",
                   "jacobian-followup/sparse-seven/search.py",
                   "jacobian-followup/sparse-seven/run_all.py",
                   "jacobian-followup/sparse-seven/verify.py"],
    supersedes_notes=["GLDC2's normalized support lower bound eight is strengthened to nine."],
    supersedes=[],narrows_problems=["OP-GLDC"],closes_problems=[],
    consumers=[],invalidates_assumptions=[],replaced_by=[])
pair=copy.deepcopy(template)
pair.update(
    id="KPAIR1",kind="corollary",title="Degree-ten constant Jacobian pairing exclusion and sharp dimension four",
    scope="Let F be a three-variable Keller map over a characteristic-zero field with maximum ordinary coordinate degree D. If constant vectors v,w satisfy v^T*JF(X)*w=c nonzero identically, linear source and target changes give an affine-shear coordinate x1+H(x2,x3). A triangular source substitution makes it a parameter and bounds each plane fibre coordinate degree by D*max(1,deg H)<=D^2. Moh's established maximum-plane-component-degree-at-most-100 theorem therefore excludes every collision for D<=10, independently of support. The pairing property is affine-left-right invariant. For JF(0)=I it is exactly nonemptiness of the bilinear incidence v^T*M_alpha*w=0 for every nonconstant Jacobian coefficient matrix and v^T*w=1. A collision in raw degree at most ten must have empty incidence; emptiness is not sufficient for a collision. Univariate own-variable coordinates and nonconstant hyperplane normal derivatives are independently impossible in any degree. Among counterexamples with maximum coordinate degree at most ten admitting this pairing, the minimum dimension is four: dimensions one and two are directly injective, dimension three is excluded by the plane reduction, and the degree-seven foundational map normalized to determinant one and identity-padded has the pairing in dimension four. The checker exactly replays this four-variable collision and determinant, the nonlinear shear chain rule, and the normalized foundational map's 27-equation direction incidence as (1) in Singular and SymPy. The general reduction and dimension coverage are written proofs using Moh; that external theorem is not independently reproved. Unbounded-degree three-variable pairings, unrestricted dimension minima and cubic special-form minima remain outside this corollary. No novelty priority or external review is claimed.",
    canonical_source=BASE+"COORDINATE_GATES.md",
    dependencies=["F1","external: T.T. Moh plane degree-at-most-100 theorem, J.reine angew.Math.340(1983),140-212; primary author account https://www.math.purdue.edu/~ttm/jacobian"],
    checker=BASE+"verify_coordinate_gates.py",proof_type="hybrid",independent_replay=False,
    artifact_hash="sha256:"+hashlib.sha256((HERE/"verify_coordinate_gates.py").read_bytes()).hexdigest(),
    software_lock=[BASE+"requirements.txt",BASE+"coordinate-gates-validation.json"],
    supersedes_notes=[],supersedes=[],narrows_problems=["OP-GLDC"],
    closes_problems=[],consumers=[],invalidates_assumptions=[],replaced_by=[])
for new in (entry,pair):
    if new["id"] in known:
        known[new["id"]].clear();known[new["id"]].update(new)
    else:
        ledger["entries"].append(new);known[new["id"]]=new
old=known["GLDC2"]
old["scope"]=old["scope"].replace(
    "Eight-term attainment, global degree minima, unrestricted-degree support minima, and special-form dimension minima remain open.",
    "GLDC3 subsequently excludes eight-term attainment in the same raw-degree range; this predecessor alone determines no unrestricted degree, support or special-form dimension minimum.")
parent=known["OP-GLDC"]
parent["scope"]=(
    "Complete the support-first census of every three-dimensional Keller collision whose affine-invariant degree flag is lexicographically below (7,6,4). GLDC1 supplies the normalization, all 74 profiles, determinant-bucket and stratified face compilers, and four dense maximum-degree-two exclusions. CAUF1 supplies the axis frontend. GLDC2 closes seven terms, and GLDC3 closes eight terms in the normalized raw-degree-at-most-seven range with complete 90-shard coverage, 17482 size-eight orbit ideals (1) over Q in both Singular and SymPy, and independent enumeration and full coefficient-input compilation. The current normalized support lower bound is nine, with attainment unknown. KPAIR1 excludes every degree-at-most-ten three-variable candidate admitting a constant nonzero Jacobian pairing, by an affine-shear reduction to plane degree at most D^2<=100; univariate-coordinate and hyperplane-normal-factor gates work in arbitrary degree. The next search must start at nine and use these genuine coordinate obstructions before coefficient algebra, while preserving possible repairs of partial supports. The remaining dense census needs filtered-output rank minors, common leading-coefficient torus solvability and actual compactification of any surviving coefficient families. No dense profile above maximum degree two is fully excluded without a support bound; global and special-form complexity minima beyond KPAIR1's bounded four-dimensional statement remain open.")
parent["dependencies"]=list(dict.fromkeys(parent["dependencies"]+["GLDC3","KPAIR1"]))
parent["forbidden_attack_classes"][0]=dict(
    attack="Repeating normalized raw-degree-at-most-seven searches with nonlinear support at most eight",
    reason="GLDC1 through GLDC3 exhaust this domain. All residual exact-support ideals are unit over Q, with two rational backends, independent exact support coverage and derivative-level input validation.",
    witnesses=["GLDC1","GLDC2","GLDC3"])
parent["forbidden_attack_classes"][1]["reason"]=(
    "The axis-only support minimum three is merely a frontend. The normalized global support lower bound nine uses full support coverage and coefficient compatibility. No dense profile above maximum degree two is fully excluded.")
save("MATH_STATUS.json",ledger)
lessons=load("knowledge/lessons.json")
for item in lessons["lessons"]:
    if item["id"]=="METHOD-KELLER-AXIS-SEEDED-CLOSURE-20261114":
        item["avoid"]="Treating support-eight attainment as open in the completed raw-degree-at-most-seven range; traversing the entire Boolean support universe without the axis seed; using matching totals or two solvers on unchecked common inputs as coverage."
        item["boundary"]="The 90-shard cover through eight terms is closed, including all 17482 size-eight orbit ideals. The normalized lower bound nine is specific to raw degree at most seven; nine-term attainment and unrestricted complexity minima remain open."
        item["reopen_when"]="A candidate has at least nine normalized nonlinear occurrences, exceeds the raw degree bound, or introduces a structural theorem beyond the completed domain."
        item["claims"]=list(dict.fromkeys(item["claims"]+["GLDC3"]))
        item["sources"]=list(dict.fromkeys(item["sources"]+[BASE+"PROOF.md"]))
        item["implementation"]=list(dict.fromkeys(item["implementation"]+[BASE+"verify.py",BASE+"verify_coverage.py"]))
    if item["id"]=="METHOD-KELLER-COORDINATE-GATES-20261114":
        raise AssertionError("coordinate lesson already present; inspect before repeating registration")
    for key in ("OP-GLDC","GLDC2","GLDC3"):
        if key in item.get("claims",[]):
            item["reviewed_claims"][key]=claim_fingerprint(known[key])
lessons["lessons"].append(dict(
    id="METHOD-KELLER-COORDINATE-GATES-20261114",area="geometry",
    title="Exclude low-degree constant Jacobian pairings before sparse algebra",
    when="A normalized low-degree three-variable Keller collision candidate is considered.",
    use="Reject completed supports with a nonlinear univariate coordinate or a nonconstant normal derivative on a coordinate hyperplane. In degree at most ten, a constant nonzero v^T*JF*w is excluded by a triangular source reduction to a plane map of degree at most D^2<=100 and Moh. For normalized maps search the finite bilinear direction incidence. The complete size-eight cover has 17435 univariate, 46 elementary-shear, and one hyperplane-normal-factor orbit.",
    avoid="Discarding a repairable partial support because its current row is univariate or triangular; treating empty direction incidence as a counterexample; applying the D<=10 consequence at unbounded degree.",
    boundary="The univariate and hyperplane-normal-factor gates are arbitrary-degree necessities. The pairing gate imports Moh's degree-100 plane theorem. The special-form minimum four is for degree-at-most-ten counterexamples admitting a constant nonzero pairing, with a degree-seven stabilized witness; it is not an unrestricted dimension minimum.",
    reopen_when="The completed candidate passes all three gates, exceeds degree ten for the pairing reduction, or supplies a new plane theorem.",
    claims=["GLDC3","KPAIR1"],sources=[BASE+"PROOF.md",BASE+"COORDINATE_GATES.md"],
    implementation=[BASE+"verify_structure.py",BASE+"verify_coordinate_gates.py"],
    keywords=["Keller","support","Jacobian pairing","coordinate","Moh","degree 100","dimension four"],
    reviewed_claims={k:claim_fingerprint(known[k]) for k in ("GLDC3","KPAIR1")},
    programme_status="archived"))
save("knowledge/lessons.json",lessons)
work=load("knowledge/work_items.json")
for item in work["actions"]:
    if item["id"]=="WORK-GLDC":
        item.update(
            title="Find a nine-term Keller collision after the coordinate obstructions",
            next_step="GLDC3 closes eight terms. Seek a normalized size-nine collision in raw degree at most seven, first imposing the univariate, hyperplane-normal and constant-Jacobian-pairing gates of KPAIR1 on completed supports. Preserve every possible repair of partial supports, exact coefficient-torus coverage and actual boundary conditions.",
            done_when="A certified nine-term witness, or a complete next exclusion with independent set-level coverage and full coefficient input verification.",
            prerequisites="The normalized lower bound is nine, with attainment unknown. Exactly four dense quadratic profiles are closed. Constant nonzero Jacobian pairings are excluded in three variables through degree ten; dimension four is attained in that bounded special class by the stabilized degree-seven map. These results do not settle arbitrary dense low-degree or unrestricted special-form minima.")
        item["claims"]=list(dict.fromkeys(item["claims"]+["GLDC3","KPAIR1"]))
        item["sources"]=list(dict.fromkeys(item["sources"]+[BASE+"PROOF.md",BASE+"COORDINATE_GATES.md"]))
    for key in ("OP-GLDC","GLDC2","GLDC3","KPAIR1"):
        if key in item.get("claims",[]):
            item["reviewed_claims"][key]=claim_fingerprint(known[key])
save("knowledge/work_items.json",work)
reviews=load("knowledge/partial_reviews.json")
for group in ("reviews","checker_absence_reviews"):
    for item in reviews.get(group,[]):
        for key in ("OP-GLDC","GLDC2"):
            if key in item.get("reviewed_claims",{}):
                item["reviewed_claims"][key]=entry_fingerprint(known[key])
save("knowledge/partial_reviews.json",reviews)
print("Registered GLDC3 (lower bound nine) and KPAIR1 (bounded pairing and dimension-four corollary); OP-GLDC remains open.")

