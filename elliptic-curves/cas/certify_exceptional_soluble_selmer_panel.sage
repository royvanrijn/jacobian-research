#!/usr/bin/env sage-python
"""Replay known soluble residual 2-Selmer subspaces on eleven fixed fibres.

No rational-point search, class group, or full descent is run. Small good-prime
square characters certify selected known Kummer classes; exact intersections
of quadrics and rational witnesses certify their global solubility. Checkpoint
after every fibre. Missing Selmer complements and their CT entries stay UNKNOWN.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

from sage.all import EllipticCurve, GF, Matrix, PolynomialRing, QQ, ZZ, lcm, prime_range
from sage.version import version as sage_version

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))
from build_bnf_free_two_covers import (  # noqa: E402
    cover_for, multiply_mod_cubic, verify_rational_cover_witness,
)
from icarm_curve398 import GENERAL_WEIERSTRASS_COEFFICIENTS, POINTS  # noqa: E402

OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/exceptional_soluble_selmer_panel_v1.json"
CHECKPOINT = ROOT / "artifacts/local/elliptic-curves/exceptional-soluble-selmer-panel-v1/checkpoint.json"
PUBLIC = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
LINEAGE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
PARENTS = ROOT / "elliptic-curves/data/icarm_mw16_parent_ladder_blind_inputs_v1.json"
MW16_PUBLIC = ROOT / "elliptic-curves/data/icarm_mw16_public_targets_v1.json"
MW16_RESULTS = ROOT / "artifacts/generated-results/elliptic-curves/icarm_mw16_parent_ladder_blind_v1.json"
REFRESH_INPUT = ROOT / "elliptic-curves/data/r17_refresh_jump_ladder_blind_inputs_v1.json"
REFRESH_RESULTS = ROOT / "artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_blind_v2.json"
CALIBRATION = ROOT / "artifacts/generated-results/elliptic-curves/icarm_mw16_blind_ladder_calibration_v1.json"
CURVE_IDS = (351, 356, 376, 377, 385, 398, 400, 401, 542, 543, 548)
RANKS = dict(zip(CURVE_IDS, (25, 29, 22, 23, 29, 30, 28, 27, 26, 29, 24)))


def read(path):
    return json.loads(path.read_text())


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def xy(row):
    return [QQ(row["x"]), QQ(row["y"])] if isinstance(row, dict) else list(map(QQ, row))


def f2_rank(rows):
    return int(Matrix(GF(2), rows).rank()) if rows else 0


def inputs():
    if read(LINEAGE)["status"] != "PROVED_EXACT_LINEAGE_REALIZATION_AND_DISPLAYED_QUOTIENTS":
        raise ArithmeticError("published-R17 generic-subgroup identification changed")
    if read(PARENTS)["status"] != "PASS_EXACT_COMPLEMENT_BLIND_NINE_PARENT_INPUTS":
        raise ArithmeticError("MW16 generic-subgroup input status changed")
    public = {int(r["id"]): r for r in read(PUBLIC)["records"]}
    parents = read(PARENTS)["parents"]
    other_public = {int(r["id"]): r for r in read(MW16_PUBLIC)["records"]}
    ladder = read(MW16_RESULTS)["parents"]
    result = []
    for curve_id in CURVE_IDS:
        if curve_id in (351, 356, 376, 377, 385):
            r = public[curve_id]
            model = list(map(QQ, r["ainvs"]))
            generic = list(map(xy, r["points"][:17]))
            candidates = list(map(xy, r["points"][17:]))
            labels = [f"public:P{i}" for i in range(18, len(r["points"]) + 1)]
            frame = "published-R17"
            source = str(PUBLIC.relative_to(ROOT))
        elif curve_id == 543:
            r = next(r for r in read(REFRESH_INPUT)["cases"] if int(r["curve_id"]) == curve_id)
            known = next(r for r in read(REFRESH_RESULTS)["results"] if int(r["curve_id"]) == curve_id)
            model = list(map(QQ, r["short_model"]))
            generic = list(map(xy, r["generic_points"]))
            candidates = list(map(xy, known["final_basis"]))
            labels = [f"blind-final:B{i + 1}" for i in range(len(candidates))]
            frame = r["native_chart"]
            source = str(REFRESH_RESULTS.relative_to(ROOT))
        else:
            p = min((p for p in parents if p["curve_id"] == curve_id), key=lambda p: p["priority_rank"])
            model = list(map(QQ, p["target_short_model"]))
            curve = EllipticCurve(QQ, model)
            generic = list(map(xy, p["specialized_generic_points"]))
            frame = p["parent_id"]
            if curve_id in (398, 400, 401):
                original_model = GENERAL_WEIERSTRASS_COEFFICIENTS if curve_id == 398 else other_public[curve_id]["ainvs"]
                original_points = POINTS if curve_id == 398 else other_public[curve_id]["points"]
                original_curve = EllipticCurve(QQ, list(map(QQ, original_model)))
                transport = original_curve.isomorphism_to(curve)
                candidates = [list(transport(original_curve(*xy(p))))[:2] for p in original_points]
                labels = [f"public:P{i + 1}" for i in range(len(candidates))]
                source = "elliptic-curves/cas/icarm_curve398.py" if curve_id == 398 else str(MW16_PUBLIC.relative_to(ROOT))
            else:
                known = next(r for r in ladder if r["parent_id"] == p["parent_id"])
                events = [e for e in known["discovered_group_saturation"]["events"] if e["type"] == "NEW_Q_INDEPENDENT_DIRECTION"]
                candidates = [xy(e["point"]) for e in events]
                labels = [f"blind-discovery:D{i + 1}" for i in range(len(candidates))]
                source = str(MW16_RESULTS.relative_to(ROOT))
        result.append(dict(curve_id=curve_id, model=model, generic=generic,
                           candidates=candidates, labels=labels, frame=frame, point_source=source))
    return result


def cubic_data(model, points):
    a1, a2, a3, a4, a6 = model
    b2, b4, b6 = a1*a1+4*a2, a1*a3+2*a4, a3*a3+4*a6
    ring = PolynomialRing(QQ, "theta")
    theta = ring.gen()
    f = theta**3 + b2*theta**2 + 8*b4*theta + 16*b6
    curve = EllipticCurve(QQ, model)
    if not f.is_irreducible():
        raise ArithmeticError("panel hypothesis E(Q)[2]=0 failed")
    transformed = []
    for x, y in points:
        curve(x, y)
        X, Y = 4*x, 4*(2*y+a1*x+a3)
        if Y**2 != f(X):
            raise ArithmeticError("point failed exact cubic transport")
        transformed.append((X, Y))
    return f, transformed


def signatures(f, points, generic_rank, target_rank, prime_bound):
    """Linear residue factors suffice: irreducible factors add no needed columns.

    Every column is a homomorphism on these unit squareclasses. We use it only
    for independence, never to certify dependence or unknown Selmer membership.
    """
    rows = [[] for _ in points]
    blocks = []
    disc = f.discriminant()
    denominators = [c.denominator() for c in f] + [X.denominator() for X, _ in points]
    for p in prime_range(3, prime_bound + 1):
        if disc.numerator() % p == 0 or any(d % p == 0 for d in denominators):
            continue
        field = GF(p)
        reduced = f.change_ring(field)
        roots = sorted(reduced.roots(multiplicities=False), key=int)
        if not roots:
            continue
        values = [[field(X)-root for root in roots] for X, _ in points]
        if any(v == 0 for row in values for v in row):
            continue
        extra = [[0 if v.is_square() else 1 for v in row] for row in values]
        trial = [row+new for row, new in zip(rows, extra)]
        if f2_rank(trial) == f2_rank(rows):
            continue
        rows = trial
        blocks.append(dict(prime=int(p), linear_roots=[int(r) for r in roots], rank_after=f2_rank(rows)))
        if f2_rank(rows) >= target_rank and f2_rank(rows[:generic_rank]) == generic_rank:
            break
    return rows, blocks


def primitive_quadric(poly):
    coefficients = poly.coefficients()
    den = lcm([c.denominator() for c in coefficients])
    ints = [ZZ(c*den) for c in coefficients]
    divisor = ZZ(0)
    for c in ints:
        divisor = divisor.gcd(c)
    normalized = poly * (den / divisor)
    return normalized, max(abs(ZZ(c)).nbits() for c in normalized.coefficients())


def build_cover(f, point, label, ring):
    X, Y = point
    alpha = [X, QQ(-1), QQ(0)]
    coeffs = f.list()
    witness = list(map(QQ, [1, 0, 0, 1]))
    if verify_rational_cover_witness(alpha, coeffs, witness, ring) != X:
        raise ArithmeticError("cover x-map failed")
    u, v, w, z = ring.gens()
    beta = [u, v, w]
    products = multiply_mod_cubic(alpha, multiply_mod_cubic(beta, beta, coeffs), coeffs)
    # det(multiplication by beta) is its cubic norm; verify the complete y-map.
    powers = ([1, 0, 0], [0, 1, 0], [0, 0, 1])
    norm_beta = Matrix(ring, [multiply_mod_cubic(beta, basis, coeffs) for basis in powers]).det()
    if Y * norm_beta(*witness) / witness[3]**3 != Y or f(X) != Y**2:
        raise ArithmeticError("cover y-map failed")
    q1, bits1 = primitive_quadric(products[1]+z*z)
    q2, bits2 = primitive_quadric(products[2])
    return dict(label=label, alpha_coefficients=list(map(str, alpha)),
                cubic_point=[str(X), str(Y)], norm_square_root=str(Y),
                quadrics=cover_for(alpha, coeffs, ring),
                primitive_quadrics=[str(q1), str(q2)],
                affine_y_map=f"({Y})*({norm_beta})/z^3",
                rational_witness=["1", "0", "0", "1"], witness_verified=True,
                local_solubility="PROVED_AT_EVERY_PLACE_BY_RATIONAL_WITNESS",
                sha_image="ZERO", ct_pairing_with_every_selmer_class="ZERO_BY_GLOBAL_KUMMER_EXACT_SEQUENCE",
                primitive_quadric_coefficient_max_bits=max(bits1, bits2))


def build_case(case, prime_bound):
    points = case["generic"] + case["candidates"]
    r = len(case["generic"])
    f, transformed = cubic_data(case["model"], points)
    rows, blocks = signatures(f, transformed, r, RANKS[case["curve_id"]], prime_bound)
    gr = f2_rank(rows[:r])
    if gr != r:
        raise ArithmeticError("generic mod-two independence not certified within prime bound")
    chosen = list(range(r))
    for i in range(r, len(rows)):
        if f2_rank([rows[j] for j in chosen]+[rows[i]]) > len(chosen):
            chosen.append(i)
    q = len(chosen)-r
    ring = PolynomialRing(QQ, names=("u", "v", "w", "z"))
    covers = [build_cover(f, transformed[i], case["labels"][i-r], ring) for i in chosen[r:]]
    bound = RANKS[case["curve_id"]]-r
    return dict(curve_id=case["curve_id"], frame=case["frame"], point_source=case["point_source"],
                model=list(map(str, case["model"])), generic_rank=r,
                cubic_model_transport=dict(
                    input_to_cubic_X="4*x",
                    input_to_cubic_Y=f"4*(2*y+({case['model'][0]})*x+({case['model'][2]}))",
                    cubic_to_input_x="X/4",
                    cubic_to_input_y=f"(Y-({case['model'][0]})*X-4*({case['model'][2]}))/8"),
                generic_points=[list(map(str, p)) for p in case["generic"]],
                candidate_points=[dict(label=label, point=list(map(str, p))) for label, p in zip(case["labels"], case["candidates"])],
                irreducible_cubic_coefficients_ascending=list(map(str, f.list())),
                exact_kummer_signature=dict(blocks=blocks, rows=rows, generic_rank=gr,
                                           selected_independent_point_indices=chosen, selected_rank=len(chosen)),
                known_soluble_residual_dimension_lower_bound=q,
                demonstrated_rank_gain_target=bound,
                target_recovered_by_mod2_certificate=(q >= bound),
                known_soluble_nonzero_class_count=2**q-1,
                basis_covers=covers,
                known_soluble_ct_matrix=[[0]*q for _ in range(q)],
                ct_method="EXACT_GLOBAL_POINT_WITNESSES_FORCE_ZERO_AGAINST_FULL_SELMER",
                full_residual_selmer_dimension=None, full_ct_matrix=None,
                full_ct_radical_dimension=None, sha_2_dimension=None,
                certified_insoluble_selmer_classes=[],
                unconstructed_complement_status="UNKNOWN",
                complexity=dict(model_dependent=True, minimized_over_equivalent_covers=False,
                                basis_max_coefficient_bits=[c["primitive_quadric_coefficient_max_bits"] for c in covers]),
                status="PASS_EXACT_SOLUBLE_SUBSPACE_LOWER_BOUND")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--prime-bound", type=int, default=2000)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--checkpoint", type=Path, default=CHECKPOINT)
    args = parser.parse_args()
    if not 3 <= args.prime_bound <= 5000:
        parser.error("fixed-panel auxiliary prime bound must be between 3 and 5000")
    paths = [Path(__file__).resolve(), CAS / "build_bnf_free_two_covers.py", CAS / "icarm_curve398.py",
             PUBLIC, LINEAGE, PARENTS, MW16_PUBLIC, MW16_RESULTS, REFRESH_INPUT, REFRESH_RESULTS, CALIBRATION]
    hashes = {str(p.relative_to(ROOT)): digest(p) for p in paths}
    runs = []
    for case in inputs():
        runs.append(build_case(case, args.prime_bound))
        args.checkpoint.parent.mkdir(parents=True, exist_ok=True)
        args.checkpoint.write_text(json.dumps(dict(input_hashes=hashes, prime_bound=args.prime_bound, curves=runs), indent=2, sort_keys=True)+"\n")
        print(f"SOLUBLESELMER|curve={case['curve_id']}|generic={runs[-1]['generic_rank']}|known_W_lower={runs[-1]['known_soluble_residual_dimension_lower_bound']}", flush=True)
    payload = dict(schema="elliptic-curves.exceptional-soluble-selmer-panel.v1",
                   status="PASS_ELEVEN_FIBRE_KNOWN_SOLUBLE_SUBSPACES", input_hashes=hashes,
                   software=dict(sage_version=sage_version), prime_bound=args.prime_bound, curves=runs,
                   budget=dict(curve_count=11, auxiliary_prime_bound=args.prime_bound, rational_point_searches=0,
                               full_descents=0, class_group_computations=0, checkpoint_after_each_fibre=True),
                   total_independent_soluble_basis_covers=sum(len(r["basis_covers"]) for r in runs),
                   mw16_existing_blind_calibration=read(CALIBRATION)["curve_results"],
                   claim_boundary=["These are lower bounds for W_G and the residual Selmer quotient, not their complete dimensions.",
                                   "The known-soluble CT block is rigorously zero; no arithmetic pairing on an unknown complement was computed.",
                                   "There are no certified insoluble Selmer controls in this panel. Missing witnesses and bounded misses are never insolubility certificates.",
                                   "Cover coefficients depend on the chosen model and known points. They give no prospective low-complexity predictor.",
                                   "Existing midpoint-search charts are birational to E; with cover map 2R-Q they have residual class zero when Q is generic."])
    text = json.dumps(payload, indent=2, sort_keys=True)+"\n"
    if args.check:
        if args.output.read_text() != text:
            raise ArithmeticError("stored soluble Selmer panel differs from exact replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(f"SOLUBLESELMER|status={payload['status']}|basis_covers={payload['total_independent_soluble_basis_covers']}", flush=True)


if __name__ == "__main__":
    main()
