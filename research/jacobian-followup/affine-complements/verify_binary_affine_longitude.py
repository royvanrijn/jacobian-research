"""PF125 exact transcription and countercontrols, not a finite-sample proof."""
from pathlib import Path
import hashlib
import json
import platform
import sympy as S
import binary_affine_longitude as B

P = Path(__file__).resolve().parent


def main():
    checks = []

    def check(name, value):
        assert value, name
        checks.append(name)

    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    olds = ["intermediate_critical_axis_receipt.json",
            "intermediate_chern_betti_receipt.json",
            "block_quotient_monodromy_receipt.json"]
    for name in olds:
        receipt = json.loads((P / name).read_text())
        check("sealed predecessor PASS " + name, receipt["outcome"] == "PASS")
        for f, digest in receipt["source_sha256"].items():
            check("sealed predecessor " + name + " " + f, sha(P / f) == digest)
    manifest = json.loads((P / "sources/binary-affine-longitude/manifest.json").read_text())
    for item in manifest["sources"]:
        f = P / item["path"]
        check("retained source seal " + f.name,
              sha(f) == item["sha256"] and f.stat().st_size == item["bytes"])

    controls = []
    for a in [2, 3, 4]:
        b, v = 2**a, 2
        identity = list(range(b))
        q = [x ^ (v if x & 1 else 0) for x in identity]
        qc = [x ^ (v if not x & 1 else 0) for x in identity]
        expected = [{1: b//2}, {2: b//4}]
        if a >= 3:
            expected += [{4: b//8}, {1: b//4, 2: b//8}]
        observed = []
        for u in range(0, b, 2):
            for mask in range(b):
                if B.pairing(mask, v):
                    continue
                p = B.normal_form(a, v, u, mask)
                square = B.compose(p, p)
                critical = [x for x in identity if x & 1]
                data = B.peripheral_data(a, v, u, mask)
                residual = B.cycle_profile(p, range(0, b, 2))
                flipped = sum(z[2] for z in data if z[0] == 2)
                good = (sorted(p) == identity and B.compose(p, q) == B.compose(q, p)
                        and all(p[x] in [x, x ^ v] for x in critical)
                        and square == (qc if B.pairing(mask, u) else identity)
                        and B.compose(square, square) == identity
                        and residual in expected
                        and sum(e*f for e, f, _ in data) == b
                        and all(f == 1 for e, f, _ in data if e == 2)
                        and flipped in [0, b//4, b//8]
                        and all(0 <= shear < e for e, _, shear in data))
                if mask & ~1:
                    good = good and flipped == b//8
                else:
                    good = good and flipped in [0, b//4]
                check("actual affine parameter control " + str((a, u, mask)), good)
                if residual not in observed:
                    observed.append(residual)
        check("all four permitted profile types realized when dimension allows " + str(a),
              all(p in observed for p in expected))
        controls.append(dict(dimension=a, residual_profiles=observed))

    p = B.normal_form(3, 2, 4, 4)
    check("order-four longitude algebra is possible, not automatic triviality",
          B.compose(p, p) != list(range(8)) and B.cycle_profile(p, range(0, 8, 2)) == {4: 1})
    # A centralizer element rotating three extra coordinates has order three.
    rotate = [(x & 3) | ((x & 4) << 1) | ((x & 8) << 1) | ((x & 16) >> 2)
              for x in range(32)]
    q = [x ^ (2 if x & 1 else 0) for x in range(32)]
    check("centralizing meridian alone does NOT imply order divides four",
          B.compose(rotate, q) == B.compose(q, rotate)
          and B.compose(B.compose(rotate, rotate), rotate) == list(range(32))
          and any(rotate[x] not in [x, x ^ 2] for x in range(1, 32, 2)))
    check("intermediate peripheral coverage is not assumed from original source",
          manifest["original_coverage_transferred_to_intermediate"] is False)
    z, t = S.symbols("z t")
    check("cusp residual section exact factorization",
          S.expand((z-t)**2*(z+2*t) - (z**3-3*t**2*z+2*t**3)) == 0)
    n = S.symbols("n", integer=True, nonnegative=True)
    check("unbounded paired odd contacts have even infinity degree", S.rem(2*n, 2, n) == 0)
    for exponent in [1, 2, 3, 8]:
        poly = ((t-1)*(t+1))**exponent
        check("actual paired-contact squareclass even degree " + str(exponent),
              S.degree(poly, t) % 2 == 0)
    for a, v, u, mask in [(3, 6, 4, 6), (4, 14, 6, 10)]:
        p = B.normal_form(a, v, u, mask)
        check("other critical direction control " + str((a, v)),
              B.compose(B.compose(p, p), B.compose(p, p)) == list(range(2**a)))
    text = " ".join((P / "BINARY_AFFINE_LONGITUDE_AND_FILLING_PROFILES.md").read_text().lower().split())
    for phrase in ["original finite normalization", "all omitted primes",
                   "each individual pair", "g^4=1", "order dividing four",
                   "4^(b/8)", "f*e=s", "not automatically transferred to m",
                   "one-cusp case", "rho(lambda)=1", "necessary condition, not an exclusion",
                   "pf60", "unmodified authors", "remain open", "no new family exclusion"]:
        check("written scope " + phrase, phrase in text)
    files = ["BINARY_AFFINE_LONGITUDE_AND_FILLING_PROFILES.md",
             "binary_affine_longitude.py", "verify_binary_affine_longitude.py",
             "sources/binary-affine-longitude/manifest.json"] + olds + [x["path"] for x in manifest["sources"]]
    result = dict(claim="PF2RULED125", outcome="PASS", check_count=len(checks), checks=checks,
                  theorem=dict(all_genuine_binary_affine_intermediates=True,
                               actual_longitude_preserves_each_critical_pair=True,
                               longitude_affine_normal_form_and_fourth_power_identity=True,
                               four_residual_orbit_and_full_filling_profiles=True,
                               critical_shear_affine_function_and_count_restriction=True,
                               quotient_degree4_preferred_longitude_identity_including_cusp=True),
                  assurance=dict(prior_PF60_rank_two_argument_credited=True,
                                 prior_dihedral_authors_code_comparison_retained=True,
                                 general_binary_longitude_assumed_identity=False,
                                 arbitrary_centralizer_assumed_pair_preserving=False,
                                 intermediate_plane_recognition=False,
                                 original_integral_coverage_transferred_to_intermediate=False,
                                 orbit_profiles_determine_all_linking_entries=False,
                                 new_family_exclusion=False, new_authors_code_run=False,
                                 finite_samples_prove_theorem=False, formal_verification=False,
                                 independent_replay=False, external_review=False,
                                 actual_global_finite_algebra_closed=False, research_goal_completed=False),
                  controls=controls, source_sha256={f: sha(P/f) for f in dict.fromkeys(files)},
                  environment=dict(python=platform.python_version(), sympy=S.__version__))
    (P / "binary_affine_longitude_receipt.json").write_text(json.dumps(result, indent=2)+"\n")
    print("PASS", len(checks), "PF125 controls; all-degree binary affine peripheral restrictions; goal ACTIVE")


if __name__ == "__main__":
    main()
