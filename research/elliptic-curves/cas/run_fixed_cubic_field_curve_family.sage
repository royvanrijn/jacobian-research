#!/usr/bin/env sage
"""Run the class-group-free fixed-cubic-field local Kummer experiment.

The default anchor is the pinned Fermigier--Mestre rank-at-least-20 curve.
For u=-2,-1,0,1,2 this constructs

    alpha_u = theta + u*theta^2,
    E_u: y^2 = Norm(x-alpha_u),

computes every relevant local condition on the full known 20-dimensional
Kummer span, and takes one GF(2) kernel.  No class group is requested.

A local row is accepted only after exact squareclass tests find the known
dimension of E_u(Q_p)/2E_u(Q_p).  Any incomplete place makes the whole run
fail closed and suppresses W_u.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Any, Iterable, Sequence

from sage.all import (  # type: ignore
    AA,
    EllipticCurve,
    GF,
    PolynomialRing,
    QQ,
    ZZ,
    inverse_mod,
    kronecker,
    pari,
)
from sage.version import version as sage_version  # type: ignore


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from fixed_cubic_field_curve_family import (  # noqa: E402
    bounded_integer_parameters,
    cubic_discriminant,
    discriminant_multiplier,
    f2_kernel_masks,
    f2_rank,
    field_product,
    fixed_field_cubic_coefficients,
    inverse_theta_coefficients,
    mask_indices,
)
from run_fermigier_rank20_pari_descent import (  # noqa: E402
    KNOWN_RANK,
    load_descent_basis,
    sage_q,
)


PROTOCOL = "FIXEDCUBICLOCAL"
SCHEMA = "elliptic-curves.fixed-cubic-field-varying-curve-local-kummer.v1"
PASS_STATUS = "PASS_EXACT_FULL_SPAN_LOCAL_INTERSECTIONS_NO_CLASS_GROUP"
INCOMPLETE_STATUS = "INCOMPLETE_LOCAL_KUMMER_IMAGE_NO_WU_CLAIM"
DEFAULT_MANIFEST = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "fermigier_rank20_near_miss_v1.json"
)
DEFAULT_CANDIDATE = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "fixed_cubic_field_fermigier_rank20_local_kummer_u2_v1.json"
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def canonical_hash(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return sha256(payload.encode()).hexdigest()


def qtext(value: object) -> str:
    value = QQ(value)
    return str(value)


def qpari(value: object):
    value = QQ(value)
    return pari(int(value.numerator())) / pari(int(value.denominator()))


def factor_record(value: object, retained=None) -> list[dict[str, int]]:
    value = ZZ(abs(QQ(value).numerator()))
    if value in (0, 1):
        return []
    factors = ([{"prime": int(p), "exponent": int(e)} for p, e in value.factor()]
               if retained is None else retained)
    if (len({row["prime"] for row in factors}) != len(factors)
        or any(type(row["exponent"]) is not int or row["exponent"] <= 0
               or not ZZ(row["prime"]).is_prime(proof=True) for row in factors)):
        raise ArithmeticError("invalid retained factorization")
    product = ZZ(1)
    for row in factors:
        p = ZZ(row["prime"])
        if not p.is_prime():
            raise ArithmeticError(f"unproved factor in exact support: {p}")
        product *= p ** row["exponent"]
    if product != value:
        raise ArithmeticError("factorization did not reconstruct its input")
    return factors


def rational_square_in_qp(value: object, p: int) -> bool:
    value = QQ(value)
    if value == 0:
        return True
    numerator = ZZ(value.numerator())
    denominator = ZZ(value.denominator())
    vn = numerator.valuation(p)
    vd = denominator.valuation(p)
    if (vn - vd) & 1:
        return False
    numerator //= ZZ(p) ** vn
    denominator //= ZZ(p) ** vd
    if p == 2:
        unit = numerator * inverse_mod(denominator, 8)
        return int(unit % 8) == 1
    unit = numerator * inverse_mod(denominator, p)
    return kronecker(int(unit % p), p) == 1


# Compatibility export for consumers of the original fixed-field runner.
from research_runtime.local_kummer import LocalSquareclasses
from research_runtime.arithmetic import TwoTorsionContext
from research_runtime.sage_arithmetic import SageArithmetic
from research_runtime.subspace import local_intersection


def y_discriminant(curve: object, x_value: object):
    a1, a2, a3, a4, a6 = curve.a_invariants()
    return (a1 * x_value + a3) ** 2 + 4 * (
        x_value**3 + a2 * x_value**2 + a4 * x_value + a6
    )


def candidate_x_values(curve: object, p: int) -> Iterable[QQ]:
    """Fixed local search on a p-minimal (or raw) model.

    Completeness comes from reaching the independently known local Kummer
    dimension, not from exhausting this search box.
    """

    seen = set()

    def fresh(value: object):
        value = QQ(value)
        if value in seen:
            return None
        seen.add(value)
        return value

    integer_limit = 4096 if p == 2 or p <= 127 else 256
    for numerator in range(integer_limit + 1):
        values = (QQ(0),) if numerator == 0 else (QQ(numerator), QQ(-numerator))
        for value in values:
            value = fresh(value)
            if value is not None:
                yield value

    if p != 2:
        ring = PolynomialRing(GF(p), "z")
        z = ring.gen()
        _, a2, _, a4, a6 = curve.a_invariants()
        cubic = z**3 + GF(p)(a2) * z**2 + GF(p)(a4) * z + GF(p)(a6)
        repeated = cubic.gcd(cubic.derivative())
        for root, _ in repeated.roots():
            root = ZZ(root)
            for offset in range(-64, 65):
                value = fresh(QQ(root + offset * p))
                if value is not None:
                    yield value

        # Deterministic samples across the entire residue field prevent large
        # new primes from being represented only by tiny integer residues.
        for index in range(512):
            seed = sha256(f"{p}:{index}".encode()).digest()
            value = fresh(QQ(int.from_bytes(seed, "big") % p))
            if value is not None:
                yield value

    maximum_denominator_exponent = 8 if p == 2 else 3
    numerator_limit = 1024 if p == 2 else 128
    for exponent in range(1, maximum_denominator_exponent + 1):
        denominator = ZZ(p) ** exponent
        for numerator in range(1, numerator_limit + 1):
            if numerator % p == 0:
                continue
            for value in (QQ(numerator) / denominator, -QQ(numerator) / denominator):
                value = fresh(value)
                if value is not None:
                    yield value


# Preserve the public adapter and coordinate ordering used by retained witnesses.
from research_runtime.binary import quotient_rows as quotient_rows_binary


def build(args: argparse.Namespace, retained=None) -> dict[str, Any]:
    basis_data = load_descent_basis(args.manifest, args.candidate_record)
    if (
        basis_data.mod2_rank != KNOWN_RANK
        or not basis_data.mod2_certified
        or len(basis_data.points) != KNOWN_RANK
    ):
        raise ArithmeticError("the pinned rank-20 Kummer input is not certified")

    source_curve = EllipticCurve(QQ, [sage_q(value) for value in basis_data.model])
    short_curve = source_curve.integral_model().short_weierstrass_model()
    a1, a2, a3, A, B = short_curve.a_invariants()
    if (a1, a2, a3) != (0, 0, 0):
        raise ArithmeticError("the anchor did not produce y^2=x^3+A*x+B")
    A, B = ZZ(A), ZZ(B)
    to_short = source_curve.isomorphism_to(short_curve)
    points = [
        to_short(source_curve(sage_q(x), sage_q(y)))
        for x, y in basis_data.points
    ]
    if any(point.is_zero() or point not in short_curve for point in points):
        raise ArithmeticError("failed to transport the pinned point basis")

    polynomial_ring = PolynomialRing(QQ, "x")
    x = polynomial_ring.gen()
    base_polynomial = x**3 + A * x + B
    if not base_polynomial.is_irreducible():
        raise ArithmeticError("the anchor 2-division cubic is not irreducible")
    base_discriminant = ZZ(base_polynomial.discriminant())
    if base_discriminant == 0:
        raise ArithmeticError("singular anchor")

    pari_polynomial = pari(f"y^3+({A})*y+({B})")
    arithmetic = SageArithmetic()
    torsion_context = TwoTorsionContext(tuple(map(str, base_polynomial.list())))
    base_factorization = factor_record(base_discriminant, None if retained is None else
                                       retained["anchor"]["base_discriminant_factorization"])
    nf = arithmetic.nf(torsion_context,
        factor_primes=[row["prime"] for row in base_factorization], discover=True)
    theta = pari(f"Mod(y,{pari_polynomial})")
    betas = [qpari(point[0]) - theta for point in points]
    beta_rows = [
        [qtext(point[0]), "-1", "0"]
        for point in points
    ]
    for point, beta in zip(points, betas):
        if QQ(pari.nfeltnorm(nf, beta)) != QQ(point[1]) ** 2:
            raise ArithmeticError("a pinned Kummer representative has wrong norm")

    base_support = {2, *[row["prime"] for row in base_factorization]}
    local_cache: dict[int, LocalSquareclasses] = {}

    real_roots = [root for root, multiplicity in base_polynomial.roots(AA)]
    if len(real_roots) not in (1, 3):
        raise ArithmeticError("unexpected real-root count")
    beta_real_rows = [
        [1 if AA(point[0]) - root < 0 else 0 for root in real_roots]
        for point in points
    ]

    runs = []
    all_complete = True
    j_values = []
    for run_index, u in enumerate(bounded_integer_parameters(args.parameter_bound)):
        old_run = None if retained is None else retained["runs"][run_index]
        if old_run is not None and QQ(old_run["parameter_u"]) != QQ(u):
            raise ArithmeticError("retained parameter scope mismatch")
        coefficients = fixed_field_cubic_coefficients(A, B, u)
        family_polynomial = sum(QQ(value) * x**index for index, value in enumerate(coefficients))
        multiplier = discriminant_multiplier(A, B, u)
        expected_discriminant = QQ(base_discriminant) * QQ(multiplier) ** 2
        if multiplier == 0 or QQ(family_polynomial.discriminant()) != expected_discriminant:
            raise ArithmeticError(f"family discriminant identity failed at u={u}")

        alpha = theta + qpari(u) * theta**2
        alpha_minpoly = polynomial_ring(str(pari.lift(pari.minpoly(alpha))))
        if alpha_minpoly != family_polynomial:
            raise ArithmeticError(f"alpha_u did not have F_u as minpoly at u={u}")

        c0, c1, c2, _ = coefficients
        curve = EllipticCurve(QQ, [0, QQ(c2), 0, QQ(c1), QQ(c0)])
        j_invariant = QQ(curve.j_invariant())
        j_values.append(j_invariant)

        numerator_factors = factor_record(multiplier, None if old_run is None else
                                         old_run["multiplier_numerator_factorization"])
        denominator_factors = factor_record(QQ(multiplier).denominator(), None if old_run is None else
                                           old_run["multiplier_denominator_factorization"])
        support = set(base_support)
        support.update(row["prime"] for row in numerator_factors)
        support.update(row["prime"] for row in denominator_factors)

        constraints = [[] for _ in range(KNOWN_RANK)]
        finite_records = []
        run_complete = True
        for p in sorted(support):
            if p not in local_cache:
                local_cache[p] = LocalSquareclasses(
                    nf, p, arithmetic=arithmetic, context=torsion_context)
            local = local_cache[p]
            expected_dimension = local.point_kummer_dimension
            local_generators: list[object] = []
            witnesses = []

            minimal_curve = curve.local_data(p).minimal_model()
            back_to_raw = ~curve.isomorphism_to(minimal_curve)
            model_sources = ((curve, None, "raw"), (minimal_curve, back_to_raw, "p-minimal"))
            old_local = None if old_run is None else next(
                row for row in old_run["finite_local_conditions"] if row["prime"] == p)
            for search_curve, back_map, model_label in model_sources:
                candidates = (candidate_x_values(search_curve, p) if old_local is None else
                    (QQ(row["model_x"]) for row in old_local["basis_witnesses"] if row["model"] == model_label))
                for local_x in candidates:
                    if not rational_square_in_qp(y_discriminant(search_curve, local_x), p):
                        continue
                    raw_x = (
                        local_x
                        if back_map is None
                        else QQ(back_map.x_rational_map()(local_x))
                    )
                    local_class = qpari(raw_x) - alpha
                    trial_basis, _ = local.coordinates([*local_generators, local_class])
                    if len(trial_basis) > len(local_generators):
                        local_generators = trial_basis
                        witnesses.append(
                            {
                                "model": model_label,
                                "model_x": qtext(local_x),
                                "raw_x": qtext(raw_x),
                            }
                        )
                    if len(local_generators) == expected_dimension:
                        break
                if len(local_generators) == expected_dimension:
                    break

            complete = len(local_generators) == expected_dimension
            run_complete &= complete
            local_data = curve.local_data(p)
            record = {
                "prime": p,
                "cubic_completion_factor_count": len(local.primes),
                "expected_point_kummer_dimension": expected_dimension,
                "found_point_kummer_dimension": len(local_generators),
                "point_kummer_basis_complete": complete,
                "basis_witnesses": witnesses,
                "minimal_discriminant_valuation": int(local_data.discriminant_valuation()),
                "conductor_exponent": int(local_data.conductor_valuation()),
                "bad_reduction": bool(local_data.conductor_valuation()),
            }
            if complete:
                _, coordinate_rows = local.coordinates([*local_generators, *betas])
                quotient_rows = [
                    row[len(local_generators) :]
                    for row in coordinate_rows[len(local_generators) :]
                ]
                quotient_dimension = len(quotient_rows[0]) if quotient_rows else 0
                record["known_span_quotient_coordinate_dimension"] = quotient_dimension
                record["known_span_quotient_rows"] = quotient_rows
                for index, row in enumerate(quotient_rows):
                    constraints[index].extend(row)
            finite_records.append(record)

        real_record: dict[str, Any]
        if len(real_roots) == 1:
            real_record = {
                "real_root_count": 1,
                "point_kummer_dimension": 0,
                "known_span_quotient_coordinate_dimension": 0,
                "known_span_quotient_rows": [[] for _ in range(KNOWN_RANK)],
            }
        else:
            alpha_values = [root + QQ(u) * root**2 for root in real_roots]
            order = sorted(range(3), key=lambda index: alpha_values[index])
            local_generator = [1, 1, 1]
            local_generator[order[0]] = 0
            real_quotient_rows, real_quotient_dimension = quotient_rows_binary(
                [local_generator], beta_real_rows
            )
            for index, row in enumerate(real_quotient_rows):
                constraints[index].extend(row)
            real_record = {
                "real_root_count": 3,
                "alpha_root_order_by_theta_embedding": order,
                "point_kummer_dimension": 1,
                "point_kummer_sign_generator": local_generator,
                "known_span_quotient_coordinate_dimension": real_quotient_dimension,
                "known_span_quotient_rows": real_quotient_rows,
            }

        if run_complete:
            kernel_masks = list(local_intersection(KNOWN_RANK, [constraints]))
            constraint_rank = f2_rank(constraints)
            if len(kernel_masks) != KNOWN_RANK - constraint_rank:
                raise ArithmeticError("full-span kernel dimension mismatch")
            kernel_rows = [
                [(mask >> index) & 1 for index in range(KNOWN_RANK)]
                for mask in kernel_masks
            ]
            if f2_rank(kernel_rows) != len(kernel_masks):
                raise ArithmeticError("full-span kernel basis is dependent")
            for mask in kernel_masks:
                combined = [0] * (len(constraints[0]) if constraints else 0)
                for index in mask_indices(mask, KNOWN_RANK):
                    combined = [
                        left ^ right
                        for left, right in zip(combined, constraints[index])
                    ]
                if any(combined):
                    raise ArithmeticError("reported survivor violates a local condition")
            surviving_basis = []
            beta_power_rows = [
                (QQ(point[0]), QQ(-1), QQ(0)) for point in points
            ]
            for mask in kernel_masks:
                indices = mask_indices(mask, KNOWN_RANK)
                beta = field_product(
                    [beta_power_rows[index] for index in indices], QQ(A), QQ(B)
                )
                norm_square_root = QQ(1)
                for index in indices:
                    norm_square_root *= QQ(points[index][1])
                beta_pari = qpari(beta[0]) + qpari(beta[1]) * theta + qpari(beta[2]) * theta**2
                if QQ(pari.nfeltnorm(nf, beta_pari)) != norm_square_root**2:
                    raise ArithmeticError("surviving Kummer representative has wrong norm")
                surviving_basis.append(
                    {
                        "mask": mask,
                        "one_based_anchor_basis_indices": [index + 1 for index in indices],
                        "beta_power_basis_coefficients": [qtext(value) for value in beta],
                        "norm_square_root": qtext(norm_square_root),
                    }
                )
            wu_dimension: int | None = len(kernel_masks)
        else:
            constraint_rank = None
            wu_dimension = None
            surviving_basis = []

        newly_bad = [
            record["prime"]
            for record in finite_records
            if record["prime"] not in base_support and record["bad_reduction"]
        ]
        run = {
            "parameter_u": qtext(u),
            "discriminant_multiplier": qtext(multiplier),
            "family_polynomial_ascending": [qtext(value) for value in coefficients],
            "family_discriminant": qtext(family_polynomial.discriminant()),
            "raw_curve_ainvariants": [qtext(value) for value in curve.a_invariants()],
            "j_invariant": qtext(j_invariant),
            "theta_from_alpha_power_basis": [
                qtext(value) for value in inverse_theta_coefficients(A, B, u)
            ],
            "multiplier_numerator_factorization": numerator_factors,
            "multiplier_denominator_factorization": denominator_factors,
            "complete_finite_place_support": sorted(support),
            "newly_bad_primes_relative_to_anchor": newly_bad,
            "finite_local_conditions": finite_records,
            "real_local_condition": real_record,
            "all_local_kummer_images_complete": run_complete,
            "combined_condition_column_count": len(constraints[0]),
            "combined_condition_rank": constraint_rank,
            "W_u_dimension": wu_dimension,
            "W_u_basis": surviving_basis,
        }
        runs.append(run)
        all_complete &= run_complete
        print(
            f"{PROTOCOL}|u={u}|places={len(support)}|complete={str(run_complete).lower()}"
            f"|condition_rank={constraint_rank}|W_u_dimension={wu_dimension}",
            flush=True,
        )

    if len(set(j_values)) != len(j_values):
        raise ArithmeticError("the bounded parameter set contains repeated j-invariants")
    zero_run = next(run for run in runs if run["parameter_u"] == "0")
    if zero_run["W_u_dimension"] != KNOWN_RANK:
        raise ArithmeticError("the u=0 positive control did not recover the full known span")

    document = {
        "schema": SCHEMA,
        "status": PASS_STATUS if all_complete else INCOMPLETE_STATUS,
        "protocol": f"{PROTOCOL}_v1",
        "software": {
            "sage": str(sage_version),
            "pari": str(pari.version()),
        },
        "source_hashes": {
            str(args.manifest.relative_to(ROOT)): digest(args.manifest),
            str(args.candidate_record.relative_to(ROOT)): digest(args.candidate_record),
            str(Path(__file__).resolve().relative_to(ROOT)): digest(Path(__file__).resolve()),
            str((CAS / "fixed_cubic_field_curve_family.py").relative_to(ROOT)): digest(
                CAS / "fixed_cubic_field_curve_family.py"
            ),
        },
        "anchor": {
            "id": "fermigier-mestre-rank20-u28917/20",
            "source_model_ainvariants": [qtext(value) for value in basis_data.model],
            "short_model_ainvariants": [qtext(value) for value in short_curve.a_invariants()],
            "base_polynomial_ascending": [qtext(B), qtext(A), "0", "1"],
            "base_polynomial_irreducible": True,
            "base_polynomial_discriminant": str(base_discriminant),
            "base_discriminant_factorization": base_factorization,
            "known_kummer_dimension": KNOWN_RANK,
            "known_kummer_independence_source": (
                "pinned exact 20-dimensional E(Q)/2E(Q) finite-quotient certificate"
            ),
            "known_points_on_short_model": [
                [qtext(point[0]), qtext(point[1])] for point in points
            ],
            "known_kummer_basis_beta_power_coordinates": beta_rows,
            "known_kummer_norms_are_displayed_y_squares": True,
        },
        "family": {
            "definition": "alpha_u=theta+u*theta^2; F_u(x)=Norm_K/Q(x-alpha_u)",
            "coefficient_formula_ascending": [
                "B+A*B*u^2-B^2*u^3",
                "A+3*B*u+A^2*u^2",
                "2*A*u",
                "1",
            ],
            "discriminant_identity": "disc(F_u)=disc(f)*(1+A*u^2+B*u^3)^2",
            "excluded_parameter_equation": "1+A*u^2+B*u^3=0",
            "parameter_policy": f"all integers u with |u|<={args.parameter_bound}",
            "parameter_count": len(runs),
            "pairwise_distinct_j_invariants": True,
            "u_zero_is_anchor_short_model": runs[args.parameter_bound]["raw_curve_ainvariants"]
            == [qtext(value) for value in short_curve.a_invariants()],
        },
        "covering_template": {
            "scope": "apply to every surviving span class beta, not just the original basis beta_i",
            "variables": ["a", "b", "c", "d"],
            "gamma": "a+b*theta+c*theta^2",
            "equations": [
                "coeff_theta(beta*gamma^2)+d^2=0",
                "coeff_theta2(beta*gamma^2)+u*d^2=0",
            ],
            "affine_recovery_when_d_nonzero": {
                "x": "coeff_1(beta*gamma^2)/d^2",
                "y": "norm_square_root(beta)*Norm(gamma)/d^3",
            },
        },
        "runs": runs,
        "claim_boundary": [
            "Each W_u is the kernel on the whole certified 20-dimensional span; individual-basis survival is not substituted.",
            "Every prime in 2*disc(f) and every numerator or denominator prime of 1+A*u^2+B*u^3 is checked, together with the real place.",
            "Local Kummer-image completeness is certified by exact local squareclasses plus the theoretical local dimension; a search-box miss cannot remove a class.",
            "No class group, full global 2-Selmer group, Mordell--Weil rank, or rational point on a new curve is computed.",
            "Membership in W_u proves local admissibility of a fixed global cohomology class, not realization by a rational point; the displayed covering equations are the next gate.",
        ],
        "class_group_computation_performed": False,
        "point_realization_computation_performed": False,
    }
    document["result_sha256"] = canonical_hash(document)
    return document


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--candidate-record", type=Path, default=DEFAULT_CANDIDATE)
    parser.add_argument("--parameter-bound", type=int, default=2)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    args.manifest = args.manifest.resolve()
    args.candidate_record = args.candidate_record.resolve()
    args.output = args.output.resolve()
    if args.parameter_bound < 0:
        raise ValueError("parameter bound must be nonnegative")
    if args.check:
        from research_runtime.witnesses import compare_replay
        if not args.output.is_file():
            raise SystemExit(f"missing artifact: {args.output}")
        stored = json.loads(args.output.read_text())
        recorded_hash = stored.pop("result_sha256")
        if canonical_hash(stored) != recorded_hash:
            raise ArithmeticError("retained result hash mismatch")
        document = build(args, retained=stored)
        document.pop("result_sha256")
        compare_replay(stored, document, root=ROOT,
                       source_paths=[str(Path(__file__).resolve().relative_to(ROOT))])
        print(
            f"{PROTOCOL}|status=PASS|check=true|result_sha256={recorded_hash}",
            flush=True,
        )
        return
    document = build(args)
    rendered = json.dumps(document, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(
        f"{PROTOCOL}|status={document['status']}|output={args.output}"
        f"|result_sha256={document['result_sha256']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
