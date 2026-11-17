"""Register verified nine-term exclusion and the boundary-frame support minimum."""
import copy,hashlib,json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
BASE="jacobian-followup/sparse-nine/"
sys.path.insert(0,str(ROOT/"scripts"))
from research import claim_fingerprint
from research_partial_reviews import entry_fingerprint
def load(p):return json.loads((ROOT/p).read_text())
def save(p,v):(ROOT/p).write_text(json.dumps(v,indent=2,ensure_ascii=False)+"\n")
validation=load(BASE+"validation.json")
assert validation["status"]=="pass" and validation["retained_candidates"]==0
assert validation["source_sha256"]==hashlib.sha256((HERE/"verify.py").read_bytes()).hexdigest()
assert validation["frame_certificate_independent"]["independent_coefficient_inputs"]
ledger=load("MATH_STATUS.json");known={e["id"]:e for e in ledger["entries"]}
assert not any(d.startswith("GLDC") and d[4:].isdigit() and int(d[4:])>4
               for d in known["OP-GLDC"]["dependencies"]),"Stale registration cannot lower a later frontier"
template=copy.deepcopy(known["GLDC3"])
entry=copy.deepcopy(template)
entry.update(
    id="GLDC4",title="Coordinate-gated nine-term degree-seven Keller collision exclusion",
    scope="Over every characteristic-zero field, every three-variable Keller map in the collision frame F(0)=F(e_1)=0 and JF(0)=I, with raw coordinate degrees at most seven, has at least ten nonlinear monomial occurrences counted with coordinate multiplicity. KPAIR1 forces every completed output row to contain a nonlinear own-variable term; the univariate-coordinate and hyperplane-normal-factor gates also apply. Block triangularity then forces the nonlinear variable-dependence graph to be strongly connected. A complete repair-preserving determinant-bucket search closes all 90 minimal-axis-exponent shards after 6040767 partial-support visits and retains no candidate through support nine. An independent global traversal reconstructs the universe and determinant triples, uses graph reachability instead of cut masks and repair-set intersections instead of the production budget case split, and agrees after 7228226 visits. A 49-term actual normalized degree-seven counterexample passes every gate and admits a complete branch path within its support. The original unseeded degree-three support-eight census supplies a 173-support regression, with actual checkpoint resumption tests. No residual nine-term ideal reaches algebra; no raw balanced-support census or new rational Groebner exclusion is claimed. The proof imports Moh's degree-100 plane theorem via KPAIR1 and is not a formal proof-assistant certificate or external review. Ten-term attainment, unrestricted support minima and every full dense profile above maximum degree two remain open.",
    canonical_source=BASE+"PROOF.md",checker=BASE+"verify.py",proof_type="hybrid",
    artifact_hash="sha256:"+hashlib.sha256((HERE/"verify.py").read_bytes()).hexdigest(),
    dependencies=["GLDC3","CAUF1","KPAIR1","KFRAME1"],
    software_lock=[BASE+"requirements.txt",BASE+"validation.json",BASE+"gated_search.py",
                   BASE+"run_all.py",BASE+"verify_coverage.py",BASE+"verify_controls.py",
                   BASE+"results/manifest.json",BASE+"coverage-validation.json",
                   BASE+"controls-validation.json","jacobian-followup/sparse-seven/search.py"],
    supersedes_notes=["The normalized raw-degree-at-most-seven lower bound rises from nine to ten using proved coordinate gates."],
    supersedes=[],closes_problems=[],narrows_problems=["OP-GLDC"],
    consumers=[],invalidates_assumptions=[],replaced_by=[],independent_replay=True)
frame=copy.deepcopy(template)
frame.update(
    id="KFRAME1",title="Sharp 49-term minimum in natural boundary collision frames",
    scope="For the foundational degree-seven three-variable Keller map F, consider every distinct collision pair a,b with a.x=0, and normalize by H(X)=(JF(a)*L)^(-1)*(F(a+L*X)-F(a)), L=[b-a,e_y,e_z]. The x=0 restriction is injective, so b.x is nonzero. The verified diagonal torus reduces b.x to one without changing support. All such frames are then the family a(w)=(0,4*w+6,-63*w^2-189*w-142), b(w)=(1,w,2-3*w). Its generic nonlinear support is 67. Complete exact coefficient factorization into 30 pairwise coprime monic factors gives the sharp minimum 49, attained only at w=-2 after torus normalization. At that parameter the map has determinant one, JH(0)=I, H(0)=H(e_1)=0, displayed coordinate degrees (4,6,7), and nonlinear counts (5,18,26). The previous w=-3/2 frame has 54 occurrences. Symbolic construction is independently replayed using standard-library Fraction polynomial arithmetic, including all family coefficient inputs, complete factor products, pairwise polynomial gcds, the witness determinant and exact collision; no irreducibility label is trusted. This is the known counterexample in a different affine collision frame, not a new mechanism. Sharpness is confined to source completions e_y,e_z and one collision point on x=0; arbitrary transverse completions, two interior points, and global support minima remain open. The profile is exactly (7,6,4), so this witness does not attain a flag strictly below the OP-GLDC target. No formal verification, external review or novelty priority is claimed.",
    canonical_source=BASE+"NORMALIZED_WITNESS.md",checker=BASE+"verify.py",proof_type="exact_symbolic",
    artifact_hash="sha256:"+hashlib.sha256((HERE/"verify.py").read_bytes()).hexdigest(),
    dependencies=["F1"],software_lock=[BASE+"requirements.txt",BASE+"validation.json",
                   BASE+"verify_controls.py",BASE+"frame-family-certificate.json",
                   BASE+"normalized-49-witness.json",BASE+"controls-validation.json"],
    supersedes_notes=[],supersedes=[],closes_problems=[],narrows_problems=[],
    consumers=[],invalidates_assumptions=[],replaced_by=[],independent_replay=True)
for new in (entry,frame):
    if new["id"] in known:
        assert known[new["id"]]["checker"]==new["checker"]
        known[new["id"]].clear();known[new["id"]].update(new)
    else:ledger["entries"].append(new);known[new["id"]]=new
old=known["GLDC3"]
old["scope"]=old["scope"].replace(
    "Nine-term attainment, global degree minima, and unrestricted-degree support minima remain unknown.",
    "GLDC4 subsequently excludes nine-term attainment using coordinate gates and Moh; global degree minima and unrestricted-degree support minima remain unknown.")
parent=known["OP-GLDC"]
parent["scope"]=(
    "Complete the support-first census of every three-dimensional Keller collision whose affine-invariant degree flag is lexicographically below (7,6,4). GLDC1 supplies the normalization, all 74 profiles and four dense maximum-degree-two exclusions. CAUF1 supplies the axis frontend. GLDC2 and GLDC3 close supports seven and eight by rational coefficient ideals. GLDC4 closes nine terms in raw coordinate degree at most seven using KPAIR1, own-variable and hyperplane-normal conditions, strong dependence connectivity and a complete 90-shard repair-preserving closure. Independent global traversal agrees; no residual nine-term ideal reaches algebra. The current normalized lower bound is ten, with attainment unknown. KFRAME1 supplies a verified 49-term normalized example and proves that 49 is sharp among the specified natural frames based on x=0; its profile is exactly (7,6,4), so it does not give an upper bound or witness within the strictly lower-profile census. The true minimum in the broader normalized raw-degree-at-most-seven domain is between ten and 49. Further work should combine size-ten candidates passing the coordinate gates with constructive frame changes outside the completed boundary class, preserving all repairs of partial supports. The dense census still requires filtered-output rank minors, leading-coefficient torus solvability and actual compactification. No entire dense profile above maximum degree two, unrestricted support minimum or new absolute degree minimum is settled.")
parent["dependencies"]=list(dict.fromkeys(parent["dependencies"]+["GLDC4","KFRAME1"]))
parent["forbidden_attack_classes"][0]=dict(
    attack="Repeating normalized raw-degree-at-most-seven searches with nonlinear support at most nine",
    reason="GLDC1 through GLDC4 close this domain. The nine-term proof uses necessary coordinate gates, complete repair-preserving coverage and independent traversal, rather than a raw coefficient-ideal census.",
    witnesses=["GLDC1","GLDC2","GLDC3","GLDC4"])
parent["forbidden_attack_classes"][1]["reason"]=(
    "The sharp axis-only minimum three is just a frontend. The current normalized lower bound ten uses full gated support coverage and the degree-controlled plane theorem; no dense profile above maximum degree two is fully excluded.")
save("MATH_STATUS.json",ledger)
lessons=load("knowledge/lessons.json")
for item in lessons["lessons"]:
    if item["id"]=="METHOD-KELLER-AXIS-SEEDED-CLOSURE-20261114":
        item["use"]="Start from the smallest F1 pure-axis exponent and two smallest pure-axis exponents of one transverse row. Use complete determinant repair and every necessary coordinate repair, retaining all possible additions inside a hypothetical final support. In degree at most ten enforce own-variable nonlinear terms, hyperplane-normal repairs and strongly connected dependence. Pin actual DFS frontiers and terminal status; verify coverage with independently implemented gates and budget bounds."
        item["avoid"]="Repeating support-at-most-nine searches in raw degree at most seven; discarding a repairable partial support; treating gated visit totals as a raw balanced-support census; claiming vacuous empty algebra ledgers are new Groebner certificates."
        item["boundary"]="The 90-shard gated cover through nine terms is empty. The normalized raw-degree-at-most-seven lower bound is ten. Moh is an explicit dependency of the own-variable gate. KFRAME1 supplies an upper bound 49 in this raw-degree domain, at the boundary profile (7,6,4); no lower-profile witness is constructed."
        item["reopen_when"]="A candidate has at least ten nonlinear occurrences, exceeds the raw degree bound, or changes the actual coordinate hypotheses. Constructive minimization must leave the fixed boundary frame class already optimized by KFRAME1."
        item["claims"]=list(dict.fromkeys(item["claims"]+["GLDC4","KFRAME1"]))
        item["sources"]=list(dict.fromkeys(item["sources"]+[BASE+"PROOF.md",BASE+"NORMALIZED_WITNESS.md"]))
        item["implementation"]=list(dict.fromkeys(item["implementation"]+[BASE+"gated_search.py",BASE+"verify_coverage.py",BASE+"verify.py"]))
    if item["id"]=="METHOD-KELLER-COORDINATE-GATES-20261114":
        item["use"]+=" GLDC4 also forces a nonlinear own-variable term in each row and hence a strongly connected three-vertex dependence graph by block triangularity. Missing row or cut conditions require explicit possible repairs; lower-bound repair costs may be summed only when their atom sets are disjoint."
        item["claims"]=list(dict.fromkeys(item["claims"]+["GLDC4"]))
        item["sources"]=list(dict.fromkeys(item["sources"]+[BASE+"PROOF.md"]))
        item["implementation"]=list(dict.fromkeys(item["implementation"]+[BASE+"gated_search.py",BASE+"verify_coverage.py"]))
    for key in ("GLDC3","GLDC4","KFRAME1","OP-GLDC"):
        if key in item.get("claims",[]):item["reviewed_claims"][key]=claim_fingerprint(known[key])
identifier="METHOD-KELLER-BOUNDARY-FRAME-SPARSITY-20261114"
record=dict(
    id=identifier,area="geometry",title="Optimize exact collision frames before adding search layers",
    when="A known low-degree Keller collision is normalized and support is being minimized.",
    use="Parameterize actual collision pairs, use verified monomial equivariance to reduce the family, and compute every coefficient polynomial of the normalized maps. A complete pairwise-coprime factor basis partitions all support-drop strata. Independently compile the coefficients and check factor products, gcds and the actual witness determinant using rational polynomial arithmetic.",
    avoid="Calling a frame-class minimum global; searching coefficient samples instead of all parameter-zero strata; counting generic support after a coefficient has specialized to zero; using a profile-equality witness for a strictly lower-profile objective.",
    boundary="KFRAME1 completely optimizes frames based on x=0 with transverse columns e_y,e_z for the foundational map: generic support 67 and unique normalized minimum 49 at w=-2. All coefficients and all 30 factor strata are retained. Arbitrary transverse completions and collision pairs with both x-coordinates nonzero remain outside this class.",
    reopen_when="The source completion, collision-point stratum, map family, or allowed coordinate transformations change.",
    claims=["KFRAME1"],sources=[BASE+"NORMALIZED_WITNESS.md"],
    implementation=[BASE+"verify_controls.py",BASE+"verify.py"],
    keywords=["Keller","frame","support","49","coefficient factorization","boundary","normalization"],
    reviewed_claims={"KFRAME1":claim_fingerprint(frame)},programme_status="archived")
existing=next((e for e in lessons["lessons"] if e["id"]==identifier),None)
if existing is None:lessons["lessons"].append(record)
else:existing.clear();existing.update(record)
save("knowledge/lessons.json",lessons)
work=load("knowledge/work_items.json")
for item in work["actions"]:
    if item["id"]=="WORK-GLDC":
        item.update(
            title="Close the normalized support gap ten to 49 beyond the optimized boundary frames",
            next_step="Use KPAIR1 and GLDC4's repair-safe row and connectivity conditions for size-ten candidates. Also seek a sparser actual normalized example by changing the transverse completion or using two interior collision points; the fixed natural boundary frame class has exact minimum 49 and must not be sampled again as an open route. Any witness for OP-GLDC must have profile strictly below (7,6,4).",
            done_when="A certified smaller normalized witness with its true profile, or a complete structural next exclusion with independent coverage; global sharpness requires matching bounds.",
            prerequisites="The raw-degree-at-most-seven normalized support interval is [10,49]. The upper witness has profile exactly (7,6,4), outside OP-GLDC's strict lower-profile range. Only four dense quadratic profiles are fully closed. The nine-term proof explicitly uses Moh and retains zero candidate ideals, not a raw support census.")
        item["claims"]=list(dict.fromkeys(item["claims"]+["GLDC4","KFRAME1"]))
        item["sources"]=list(dict.fromkeys(item["sources"]+[BASE+"PROOF.md",BASE+"NORMALIZED_WITNESS.md"]))
    for key in ("GLDC3","GLDC4","KFRAME1","OP-GLDC"):
        if key in item.get("claims",[]):item["reviewed_claims"][key]=claim_fingerprint(known[key])
save("knowledge/work_items.json",work)
reviews=load("knowledge/partial_reviews.json")
for group in ("reviews","checker_absence_reviews"):
    for item in reviews.get(group,[]):
        for key in ("GLDC3","OP-GLDC"):
            if key in item.get("reviewed_claims",{}):
                item["reviewed_claims"][key]=entry_fingerprint(known[key])
save("knowledge/partial_reviews.json",reviews)
print("Registered GLDC4 and KFRAME1; normalized raw-degree-seven interval [10,49]; OP-GLDC remains open.")

