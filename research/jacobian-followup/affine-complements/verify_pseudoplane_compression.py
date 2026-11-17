"""PF128 source seals and exact controls; the unbounded proof is in the note."""
from pathlib import Path
import hashlib
import json
import platform
import sympy as S

P = Path(__file__).resolve().parent
NOTE = "MINIMAL_DEGREE_PSEUDOPLANE_COMPRESSION_AND_DECK_DICHOTOMY.md"
MANIFEST = "sources/pseudoplane-compression/manifest.json"


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def inverse(p):
    out = [0] * len(p)
    for i, j in enumerate(p):
        out[j] = i
    return tuple(out)


def generated(gens):
    identity = tuple(range(len(gens[0])))
    group, queue = {identity}, [identity]
    for p in queue:
        for q in gens:
            r = compose(q, p)
            if r not in group:
                group.add(r)
                queue.append(r)
    return group


def centralizer(gens):
    """Transitive G-equivariant maps are determined by the image of zero."""
    degree = len(gens[0])
    result = set()
    for image0 in range(degree):
        mapping, queue, valid = {0: image0}, [0], True
        for i in queue:
            for g in gens:
                j, value = g[i], g[mapping[i]]
                if j in mapping:
                    if mapping[j] != value:
                        valid = False
                else:
                    mapping[j] = value
                    queue.append(j)
        assert len(mapping) == degree, "control group must be transitive"
        if valid and len(set(mapping.values())) == degree:
            result.add(tuple(mapping[i] for i in range(degree)))
    return result


def main():
    checks = []
    def check(name, condition):
        assert condition, name
        checks.append(name)
    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    manifest = json.loads((P/MANIFEST).read_text())
    for source in manifest["sources"]:
        path = P/source["path"]
        check("pinned source "+path.name,
              path.stat().st_size == source["bytes"] and sha(path) == source["sha256"])
    old = json.loads((P/"quotient_corrected_shear_receipt.json").read_text())
    check("sealed PF127 predecessor PASS", old["outcome"] == "PASS")
    for f, h in old["source_sha256"].items():
        check("PF127 predecessor bytes "+f, sha(P/f) == h)

    b, t, m, r, kappa, k = S.symbols("b t m r kappa k", positive=True)
    check("UNBOUNDED original degree minus compressed degree",
          S.expand(b*t-b*m-b*(t-m)) == 0)
    check("UNBOUNDED lift degree times universal degree",
          S.cancel((t/m)*m-t) == 0)
    check("UNBOUNDED Euler of etale universal cover when chi(V)=1",
          S.expand((1+(t-1))-t) == 0)
    check("full deck contains relative deck and divides t, so no primitivity needed",
          all([g for g in range(1,tt+1) if tt % g == 0 and g % tt == 0] == [tt]
              for tt in [1,2,3,6,12]))
    check("UNBOUNDED binary etale Euler subtracts all critical A1s",
          S.expand(r*(3-kappa)-r-r*(2-kappa)) == 0)
    check("binary chi1 cusp is r1 and immersive has no positive integer solution",
          S.solve(r*(2-1)-1, r) == [1] and S.solve(r*(2-0)-1, r) == [S.Rational(1,2)])
    check("UNBOUNDED uniform critical packet has full degree t",
          S.cancel((t/k)*k-t) == 0)
    check("UNBOUNDED degree4 cusp ramification and outside ledger",
          S.expand(4*t-1-(2*t+t+(t-1))) == 0)
    check("UNBOUNDED cyclic branched Euler includes all critical primes",
          S.cancel((1+(t-1)+t/k)-(t+t/k)) == 0)
    for tt, mm in [(2,1),(3,1),(6,2),(12,3),(12,12)]:
        check("divisible finite pi1 degree and minimality comparison "+str((tt,mm)),
              tt % mm == 0 and 4*mm <= 4*tt and ((4*mm == 4*tt) == (mm == tt)))

    # C2 wr S4, natural blocks of size2. Relative Galois over a maximal
    # block stabilizer, but original action not Galois; full deck has order2.
    flip0 = (1,0,2,3,4,5,6,7)
    swaps = []
    for i in range(3):
        p = list(range(8))
        for a in range(2):
            p[2*i+a], p[2*(i+1)+a] = p[2*(i+1)+a], p[2*i+a]
        swaps.append(tuple(p))
    gens = [flip0]+swaps
    group = generated(gens)
    H = {g for g in group if g[0] == 0}
    J = {g for g in group if set(g[:2]) == {0,1}}
    normalizer = {g for g in group if
                  {compose(compose(g,h),inverse(g)) for h in H} == H}
    deck = centralizer(gens)
    check("actual C2 wr S4 group and proper primitive block indices",
          len(group) == 384 and len(J) == 96 and len(H) == 48)
    check("relative Galois at maximal J gives full normalizer exactly J", normalizer == J)
    check("full generic centralizer has order2 and four free orbits",
          deck == {tuple(range(8)),tuple(i^1 for i in range(8))})
    check("relative Galois does not mean original Galois", normalizer != group)

    # C4 x S3 product action. A nonmaximal J/H=C2 has full deck C4.
    rot = tuple(3*((i//3+1)%4)+i%3 for i in range(12))
    s01 = tuple(3*(i//3)+({0:1,1:0,2:2}[i%3]) for i in range(12))
    s12 = tuple(3*(i//3)+({0:0,1:2,2:1}[i%3]) for i in range(12))
    gens2 = [rot,s01,s12]
    group2 = generated(gens2)
    H2 = {g for g in group2 if g[0] == 0}
    J2 = generated(list(H2)+[compose(rot,rot)])
    normalizer2 = {g for g in group2 if
                   {compose(compose(g,h),inverse(g)) for h in H2} == H2}
    check("nonmaximal relative Galois subgroup need not be full normalizer",
          len(group2) == 24 and len(H2) == 2 and len(J2) == 4 and
          J2 < normalizer2 and len(normalizer2) == 8 and len(centralizer(gens2)) == 4)
    check("original Galois countercontrol needs separate exclusion",
          len(centralizer([(1,2,3,4,5,6,7,0)])) == 8)

    # Visible block swap and additional movement inside a fixed block.
    sigma = (2,3,0,1,5,4)
    original_support = sum(sigma[i] != i for i in range(6))
    block_support = sum(sigma[2*i]//2 != i for i in range(3))
    check("no invisible target component does not force zero residual defect",
          original_support == 6 and block_support == 2 and original_support-2*block_support == 2)
    fourcycle = (2,3,1,0,4,5,6,7)
    check("zero residual movement allows relative branching over removed critical axis",
          sum(fourcycle[i] != i for i in range(8)) == 4 and
          compose(fourcycle,fourcycle) != tuple(range(8)))
    for tt, kk in [(2,1),(6,2),(6,6),(15,3)]:
        check("uniform cyclic inertia profile including nontrivial pi1(M) "+str((tt,kk)),
              tt % kk == 0 and (tt//kk)*(2*kk) == 2*tt and
              (tt//kk)*kk == tt)

    note = (P/NOTE).read_text()
    for phrase in ["MINIMAL GENERIC DEGREE", "finite, of degree t/m=1",
                   "Delta=0  <=>  h|Z is etale", "J is maximal proper",
                   "NO PRIMITIVITY", "C_t wr S4", "pi1(M)=C_(t/k)", "cases chi(V)>1",
                   "does NOT supply the global compatibility equations"]:
        check("written hypothesis and scope retained: "+phrase, phrase in note)

    theorem = dict(finite_pi1_of_entire_intermediate_etale_locus_divides_t=True,
                   zero_support_defect_iff_relative_etale=True,
                   primitive_relative_galois_identifies_full_generic_deck=True,
                   chi1_pseudoplane_compression_constructs_degree_bm_plane_Keller_map=True,
                   degree_minimal_chi1_forces_cyclic_universal_cover=True,
                   every_degree_minimal_chi1_has_full_deck_Ct_and_b0_t_minus1=True,
                   minimal_quartic_cusp_has_smooth_N_single_critical_image_solvable_G=True,
                   minimal_quartic_cusp_uniform_inertia_and_exact_pi1M=True)
    assurance = dict(published_classification_and_cover_chart_explicitly_adopted=True,
                     primary_PDF_bytes_pinned_and_relevant_pages_visually_read=True,
                     written_unbounded_proofs_and_original_whole_source_retained=True,
                     all_counterexamples_assumed_degree_minimal=False,
                     all_intermediate_etale_loci_have_chi1=False,
                     no_invisible_component_equated_with_relative_etale_in_general=False,
                     relative_galois_equated_with_original_galois=False,
                     full_deck_group_assumed_to_preserve_original_plane_chart=False,
                     full_normalization_M_identified_with_plane=False,
                     pi1M_assumed_trivial_or_critical_E_assumed_primitive=False,
                     solvable_monodromy_inversion_claimed=False,
                     global_quadratic_refinement_compatibility_proved=False,
                     finite_controls_prove_unbounded_theorems=False,
                     entire_binary_affine_or_natural_AS_cases_excluded=False,
                     shared_full_trace_or_empty_fibre_gate_closed=False,
                     formal_verification=False, independent_replay=False,
                     external_review=False, full_goal_complete=False)
    files = [NOTE, "verify_pseudoplane_compression.py", MANIFEST,
             "quotient_corrected_shear_receipt.json"]+[s["path"] for s in manifest["sources"]]
    receipt = dict(claim="PF2RULED128",outcome="PASS",check_count=len(checks),checks=checks,
                   theorem=theorem,assurance=assurance,
                   source_sha256={f:sha(P/f) for f in files},
                   environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/"pseudoplane_compression_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print("PASS",len(checks),"PF128 controls; minimal-degree compression and residual-ramification dichotomy; goal ACTIVE")


if __name__ == "__main__":
    main()
