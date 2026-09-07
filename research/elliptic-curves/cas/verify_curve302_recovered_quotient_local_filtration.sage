#!/usr/bin/env sage-python
"""Exact post-search local filtration of the recovered curve-302 quotient.

This is deliberately a retrospective calculation.  It consumes the sealed
M17-to-M24 public-span identities only after the point-search waves ended.  It
does not select a chart, search for a point, run a class-group computation, or
claim an exact Mordell--Weil rank.

The output concerns the mod-two Kummer image of the displayed rank-31 group
D.  In particular, it does *not* choose a noncanonical integral complement to
M24 in D.
"""

from __future__ import annotations

import argparse
import json
import sys
from hashlib import sha256
from math import prod
from pathlib import Path

from sage.all import AA, GF, QQ, ZZ, EllipticCurve, PolynomialRing, VectorSpace, matrix, pari, prime_range, vector
from sage.version import version as SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
ART = ROOT / "artifacts" / "generated-results" / "elliptic-curves"
PARENT = ART / "curve302_recovered_mw17_parent_v1.json"
SPAN = ART / "curve302_recovered_public_span_v1"
STRICT = ART / "rank_jump_curve302_strict_constructor_arithmetic_v1.json"
STRICT_INPUT = ART / "rank_jump_curve302_strict_constructor_arithmetic_inputs_v1.json"
OUTPUT = ART / "curve302_recovered_quotient_local_filtration_v1.json"
LOCAL_KUMMER = CAS / "research_runtime" / "local_kummer.py"
PUBLIC_SOURCE = CAS / "icarm_curve302.py"

sys.path.insert(0, str(CAS))
from research_runtime.local_kummer import LocalSquareclasses  # noqa: E402
import icarm_curve302 as curve  # noqa: E402


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.relative_to(ROOT))


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def binary_row(values) -> list[int]:
    return [int(value) for value in values]


def binary_words(kernel, lift=None) -> list[list[int]]:
    """Return a deterministic row basis, optionally in public coordinates."""
    answer = []
    for row in kernel.basis_matrix().rows():
        value = vector(GF(2), row)
        if lift is not None:
            value = lift.change_ring(GF(2)) * value
        answer.append(binary_row(value))
    return answer


def point_kummer_polynomial(point, f, ring, nf):
    """The integral representative of [4x-theta] with square norm."""
    x, y = map(QQ, point.xy())
    X, Y = 4 * x, 8 * y + 4 * x + 4
    assert Y * Y == f(X)
    denominator = ZZ(X.denominator()).sqrt()
    assert denominator in ZZ and denominator * denominator == X.denominator()
    a, b = ZZ(X * denominator**2), ZZ(Y * denominator**3)
    beta = ring([a, -denominator**2])
    gamma = pari.Mod(pari(beta), pari(f))
    assert pari.nfeltnorm(nf, gamma) == b * b
    return beta, gamma


def local_signature_matrix(public, f, ring, nf, places):
    """Exact S-local squareclass signatures for every public point."""
    beta_polynomials = []
    betas = []
    for point in public:
        beta, gamma = point_kummer_polynomial(point, f, ring, nf)
        beta_polynomials.append(beta)
        betas.append(gamma)

    signatures = [[] for _ in public]
    reports = []
    for place in places:
        chars = LocalSquareclasses(nf, place)
        rows = [binary_row(chars.signature(beta)) for beta in betas]
        for aggregate, row in zip(signatures, rows):
            aggregate.extend(row)
        reports.append(
            {
                "place": int(place),
                "ambient_signature_width": len(rows[0]),
                "point_kummer_dimension": int(chars.point_kummer_dimension),
                "public_image_dimension": int(matrix(GF(2), rows).rank()),
            }
        )

    roots = f.roots(AA, multiplicities=False)
    assert len(roots) == 3
    real_rows = [[int(beta(root) < 0) for root in roots] for beta in beta_polynomials]
    for aggregate, row in zip(signatures, real_rows):
        aggregate.extend(row)
    reports.append(
        {
            "place": "infinity",
            "ambient_signature_width": 3,
            "point_kummer_dimension": 1,
            "public_image_dimension": int(matrix(GF(2), real_rows).rank()),
        }
    )
    return matrix(GF(2), signatures), beta_polynomials, reports


def complete_split_character_rank(beta_polynomials, f):
    """An independent finite good-prime Kummer-rank certificate for D/2D."""
    rows = [[] for _ in beta_polynomials]
    blocks = []
    current_rank = 0
    for p in prime_range(3, 752):
        if f.discriminant() % p == 0:
            continue
        ff = f.change_ring(GF(p))
        roots = [int(root) for root in ff.roots(multiplicities=False)]
        if len(roots) != 3:
            continue
        values = [[int(beta.change_ring(GF(p))(root)) for root in roots] for beta in beta_polynomials]
        # A zero needs a separate valuation normalization, so omit that prime
        # from this deliberately elementary finite character certificate.
        if any(value == 0 for row in values for value in row):
            continue
        for aggregate, row in zip(rows, values):
            aggregate.extend(int(pow(value, (p - 1) // 2, p) == p - 1) for value in row)
        rank = int(matrix(GF(2), rows).rank())
        blocks.append({"prime": int(p), "roots": roots, "rank_after_block": rank})
        assert rank >= current_rank
        current_rank = rank
        if rank == 31:
            break
    assert current_rank == 31
    return blocks


def explicit_strict_witness(beta_polynomials, f):
    """One short strict word and its good-prime nonsquare witness."""
    indices = [3, 11, 17, 22, 26]  # P4,P12,P18,P23,P27; public one-based labels.
    p, root = 47, 29
    assert [int(value) for value in f.change_ring(GF(p)).roots(multiplicities=False)] == [root]
    residue = 1
    for index in indices:
        residue = residue * int(beta_polynomials[index].change_ring(GF(p))(root)) % p
    assert residue == 23 and pow(residue, (p - 1) // 2, p) == p - 1
    word = [int(index in indices) for index in range(31)]
    return {
        "public_indices_one_based": [index + 1 for index in indices],
        "public_word_mod_2": word,
        "good_prime_nonsquare_witness": {"prime": p, "root": root, "residue": residue},
    }


def build() -> dict:
    parent = load(PARENT)
    span_input = load(SPAN / "input.json")
    span_result = load(SPAN / "result.json")
    strict = load(STRICT)
    strict_input = load(STRICT_INPUT)

    assert strict["status"] == "PASS" and span_result["status"] == "PASS"
    assert strict_input["model"] == [str(value) for value in curve.GENERAL_WEIERSTRASS_COEFFICIENTS]
    assert len(span_input["public_points"]) == 31 == len(curve.POINTS)
    assert len(span_input["recovered_points"]) == 24 == len(span_result["relations"])
    assert span_result["recovered_exceptional_directions"] == 7

    # The sealed public-span package uses the integral short model.  The
    # Kummer calculation below intentionally returns to the literal public
    # generalized model, where X=4x gives the retained cubic algebra.
    E_short = EllipticCurve(QQ, list(map(QQ, span_input["curve"])))
    public = [E_short(point) for point in span_input["public_points"]]
    recovered = [E_short(point) for point in span_input["recovered_points"]]
    assert span_input["public_points"] == [[str(value) for value in point] for point in curve.SHORT_POINTS]
    assert all(point in E_short for point in public + recovered)
    generic = matrix(ZZ, parent["basis_embedding_in_public_D"])
    recovered_words = matrix(ZZ, [relation["word"] for relation in span_result["relations"]]).transpose()
    assert generic.dimensions() == (31, 17) and recovered_words.dimensions() == (31, 24)
    assert all(relation["denominator"] == 1 for relation in span_result["relations"])
    assert recovered_words[:, :17] == generic
    for point, relation in zip(recovered, span_result["relations"]):
        image = sum((coefficient * source for coefficient, source in zip(relation["word"], public)), E_short(0))
        assert image == point

    generic_smith = [int(value) for value in generic.smith_form()[0].diagonal()]
    recovered_smith = [int(value) for value in recovered_words.smith_form()[0].diagonal()]
    assert generic_smith == [1] * 17 and recovered_smith == [1] * 24
    assert generic.rank() == 17 and recovered_words.rank() == 24

    pari.allocatemem(64_000_000, 268_435_456, silent=True)
    ring = PolynomialRing(QQ, "z")
    E_general = EllipticCurve(QQ, list(map(QQ, curve.GENERAL_WEIERSTRASS_COEFFICIENTS)))
    public_general = [E_general(point) for point in curve.POINTS]
    a1, a2, a3, a4, a6 = map(QQ, curve.GENERAL_WEIERSTRASS_COEFFICIENTS)
    assert (a1, a2, a3) == (1, 1, 1)
    f = ring([64 * a6 + 16, 16 * a4 + 8, 5, 1])
    assert list(map(str, f.list())) == strict["cubic_ascending"]
    places = [int(place) for place in strict["S_finite"]]
    factorization = [(ZZ(p), int(e)) for p, e in strict_input["discriminant_factors"]]
    assert all(p.is_prime(proof=True) for p, _ in factorization)
    assert prod(p**e for p, e in factorization) == abs(ZZ(E_general.discriminant()))
    assert places == [int(p) for p, _ in factorization]
    nf = pari.nfinit([pari(f), places])

    local, beta_polynomials, local_reports = local_signature_matrix(public_general, f, ring, nf, places)
    local_dimensions = [row["point_kummer_dimension"] for row in local_reports]
    expected_generic = [row["generic_image_rank"] for row in strict["local"]]
    assert [row["place"] for row in local_reports] == [row["place"] for row in strict["local"]]
    assert local_dimensions == expected_generic
    assert sum(local_dimensions) == 22

    generic_local = generic.transpose().change_ring(GF(2)) * local
    recovered_local = recovered_words.transpose().change_ring(GF(2)) * local
    generic_local_rank = int(generic_local.rank())
    recovered_local_rank = int(recovered_local.rank())
    public_local_rank = int(local.rank())
    assert generic_local_rank == 17
    assert recovered_local_rank == public_local_rank == 21
    for index, report in enumerate(local_reports):
        width_before = sum(row["ambient_signature_width"] for row in local_reports[:index])
        width_after = width_before + report["ambient_signature_width"]
        report["generic_image_dimension"] = int(generic_local[:, width_before:width_after].rank())
        report["recovered_image_dimension"] = int(recovered_local[:, width_before:width_after].rank())
        assert report["generic_image_dimension"] == report["point_kummer_dimension"]
        assert report["recovered_image_dimension"] == report["point_kummer_dimension"]

    ambient = VectorSpace(GF(2), 31)
    generic_space = ambient.subspace(generic.change_ring(GF(2)).columns())
    recovered_space = ambient.subspace(recovered_words.change_ring(GF(2)).columns())
    strict_kernel = ambient.subspace(local.transpose().right_kernel().basis())
    recovered_strict_kernel = recovered_space.intersection(strict_kernel)
    assert (generic_space.intersection(strict_kernel).dimension(), strict_kernel.dimension()) == (0, 10)
    assert recovered_strict_kernel.dimension() == 3
    assert recovered_space + strict_kernel == ambient

    witness = explicit_strict_witness(beta_polynomials, f)
    witness_vector = vector(GF(2), witness["public_word_mod_2"])
    assert witness_vector in strict_kernel
    assert witness_vector not in recovered_space

    character_blocks = complete_split_character_rank(beta_polynomials, f)
    assert character_blocks[-1]["rank_after_block"] == 31

    return {
        "schema": "elliptic-curves.curve302-recovered-quotient-local-filtration.v1",
        "status": "PASS",
        "bindings": {
            relative(path): digest(path)
            for path in (PARENT, SPAN / "input.json", SPAN / "result.json", STRICT, STRICT_INPUT, LOCAL_KUMMER, PUBLIC_SOURCE, Path(__file__))
        },
        "software": {"sage": SAGE_VERSION, "pari": str(pari.version())},
        "displayed_lattices": {
            "D_rank": 31,
            "M17_rank": 17,
            "M24_rank": 24,
            "M24_over_M17_rank": 7,
            "D_over_M17_rank": 14,
            "D_over_M24_rank": 7,
            "M17_smith_in_D": generic_smith,
            "M24_smith_in_D": recovered_smith,
            "M24_integral_public_words": True,
        },
        "kummer_rank_certificate": {
            "complete_split_character_blocks": character_blocks,
            "public_Kummer_dimension": 31,
            "meaning": "Finite quadratic characters at complete-splitting good primes certify that the displayed public group remains 31-dimensional modulo 2.",
        },
        "local_places": local_reports,
        "local_filtration_mod_2": {
            "local_product_dimension": 22,
            "joint_local_dimensions": {"M17": generic_local_rank, "M24": recovered_local_rank, "D": public_local_rank},
            "strict_kernel_dimensions": {"M17": 0, "M24": int(recovered_strict_kernel.dimension()), "D": int(strict_kernel.dimension())},
            "quotient_local_dimension_D_over_M17": public_local_rank - generic_local_rank,
            "recovered_quotient_local_dimension_M24_over_M17": recovered_local_rank - generic_local_rank,
            "recovered_quotient_strict_dimension_M24_over_M17": 7 - (recovered_local_rank - generic_local_rank),
            "strict_residual_dimension_D_over_M24": int(strict_kernel.dimension() - recovered_strict_kernel.dimension()),
            "M24_covers_all_joint_local_patterns_of_D": recovered_local_rank == public_local_rank,
            "every_D_over_M24_mod_2_class_has_a_strict_representative": recovered_space + strict_kernel == ambient,
            "strict_kernel_public_words": binary_words(local.transpose().right_kernel()),
            "M24_strict_kernel_public_words": binary_words(recovered_local.transpose().right_kernel(), recovered_words),
        },
        "explicit_strict_witness": witness,
        "interpretation": [
            "The recovered seven-dimensional quotient maps surjectively onto all four new joint S-local patterns and has a three-dimensional strict kernel.",
            "The remaining seven-dimensional quotient D/M24 is not assigned a noncanonical integral complement.  Modulo two, every one of its classes has a representative in the seven-dimensional quotient of the ten-dimensional strict kernel by its three-dimensional intersection with M24.",
            "This is a local/Kummer filtration of the displayed group only.  It is not a full-Selmer computation, a class-group bound, an Artin/half-ideal certificate, a canonical 7+7 lattice splitting, or a rank upper bound.",
        ],
        "reproducing_command": "sage -python elliptic-curves/cas/verify_curve302_recovered_quotient_local_filtration.sage --check",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", action="store_true", help="write the immutable result once")
    parser.add_argument("--check", action="store_true", help="recompute and compare the stored result")
    args = parser.parse_args()
    assert args.build != args.check
    result = build()
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.build:
        with OUTPUT.open("x") as handle:
            handle.write(rendered)
    else:
        assert OUTPUT.read_text() == rendered
    print(
        "PASS curve302 local quotient filtration: "
        "M17=17, M24=21, D=21; recovered=4 local + 3 strict; residual strict=7",
        flush=True,
    )


if __name__ == "__main__":
    main()
