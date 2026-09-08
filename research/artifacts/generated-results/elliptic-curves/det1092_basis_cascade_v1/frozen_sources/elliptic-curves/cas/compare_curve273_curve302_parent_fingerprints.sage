#!/usr/bin/env sage-python
"""Compare the public 273 and 302 endpoints against the det-1092 parent.

This is a deliberately small, reproducible paired audit:

* exact rational-parameter recognition by ``j_1092(t)=j_273``;
* the same frozen numerical-height/integral-shell recognizer that supplied
  the 302 determinant-1092 lead, now applied to the existing primitive
  rank-17 candidate subspace of the 30 public points of 273; and
* an independent comparison of the *known-public* strict local Kummer
  kernels in the two cubic 2-division fields.

The integral-shell selection is numerical only in its frozen choice of rays.
Everything after that selection (point integrality, quadratic recovery,
root enumeration, discriminant form, and local squareclasses) is exact.
Neither a finite public-point kernel nor a candidate height form is called a
complete strict Selmer group, a generic height lattice, or a parent.
"""

from __future__ import annotations

import argparse
import gzip
import json
import runpy
import signal
from hashlib import sha256
from pathlib import Path

from sage.all import (
    AA,
    GF,
    QQ,
    RealField,
    ZZ,
    EllipticCurve,
    Genus,
    PolynomialRing,
    matrix,
    vector,
)
from sage.quadratic_forms.genera.genus import genera


ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
PARENT_LOADER = ROOT / "elliptic-curves/cas/load_curve302_recovered_parent.sage"
CURVE273 = ROOT / "elliptic-curves/cas/icarm_curve273.py"
CURVE302 = ROOT / "elliptic-curves/cas/icarm_curve302.py"
HEIGHTS = ART / "record_height_lattices_28_29_273_302_v1.json"
CORES = ART / "record_rank17_core_candidates_v1.json"
CERT273 = ART / "icarm_curve273_rank30_v1.json"
CERT302 = ART / "icarm_curve302_rank31_v1.json.gz"
OUT = ART / "curve273_curve302_parent_fingerprints_v1.json"

# These are exactly the rational bad-prime supports already pinned for the
# two public minimal models.  They define the strict-local condition used
# here, rather than a field-class-group computation.
BAD273 = (
    2, 3, 5, 7, 13, 31, 41, 47, 53, 67, 379, 4349, 25721454817,
    97018222656318846556561979214040553412450110580812087282349817173780902099339117104673990259247421230916714670243202937,
)
BAD302 = (
    2, 3, 5, 7, 11, 13, 19, 23, 29, 37, 41, 73, 131, 167, 7547,
    632881, 966509, 18145679437533309132469,
    767028866604834801397681553,
    30580600452196904409276223329355584892025407195996968868775951126238056443210297,
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def rows_as_strings(value):
    return [list(map(str, row)) for row in value.rows()]


def canonical_words(words):
    return json.dumps([list(map(int, word)) for word in words], separators=(",", ":"))


def pairs(n):
    return [(i, j) for i in range(n) for j in range(i, n)]


def gram_from_coordinates(coordinates, n):
    result = matrix(QQ, n)
    for coefficient, (i, j) in zip(coordinates, pairs(n)):
        result[i, j] = result[j, i] = coefficient
    return result


def select_integral_words(curve, basis, numerical_gram):
    """Replay the frozen 302 shell selector, with no search enlargement."""
    from sage.all import pari

    minimum = RealField(192)(pari(numerical_gram).qfminim(None, 1, 2)[1])
    found = pari(numerical_gram).qfminim(QQ(7) / 5 * minimum, 50000, 2)
    vectors = matrix(ZZ, found[2]).transpose()
    assert int(found[0]) == 2 * vectors.nrows()
    chosen = sorted(vectors.rows(), key=lambda word: (word * numerical_gram * word, tuple(word)))[:4000]
    integral = []
    for word in chosen:
        point = sum((coefficient * section for coefficient, section in zip(word, basis) if coefficient), curve(0))
        if point[0].denominator() == 1:
            integral.append(list(map(int, word)))
    return len(chosen), integral


def quadratic_kernel(words, n):
    system = matrix(
        QQ,
        [
            [word[i] * word[j] * (1 if i == j else 2) for i, j in pairs(n)] + [-1]
            for word in words
        ],
    )
    return system.right_kernel().basis_matrix()


def finite_kummer_rank(path: Path, point_count: int):
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            certificate = json.load(handle)
    else:
        certificate = json.loads(path.read_text())
    blocks = certificate["independence_certificate"].get("rows", certificate["independence_certificate"].get("signatures"))
    source_rows = [row for block in blocks for row in block["matrix_rows"]]
    result = matrix(GF(2), source_rows).rank()
    assert result == point_count
    return int(result), len(source_rows)


def two_division_cubic(coefficients):
    """Integral cubic in z=4x for a general integral Weierstrass model."""
    a1, a2, a3, a4, a6 = map(QQ, coefficients)
    assert all(entry in ZZ for entry in (a1, a2, a3, a4, a6))
    b2 = a1 * a1 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3 * a3 + 4 * a6
    ring = PolynomialRing(QQ, "z")
    return ring([16 * b6, 8 * b4, b2, 1])


def strict_local_kernel(label, source, bad_primes, certificate):
    """Compute G_public intersect ker(localization at S and infinity)."""
    from sage.all import pari

    import sys

    sys.path.insert(0, str(ROOT / "elliptic-curves/cas"))
    from research_runtime.local_kummer import LocalSquareclasses

    data = runpy.run_path(str(source))
    coefficients = tuple(map(QQ, data["GENERAL_WEIERSTRASS_COEFFICIENTS"]))
    points = [(QQ(x), QQ(y)) for x, y in data["POINTS"]]
    curve = EllipticCurve(QQ, coefficients)
    assert all(curve.is_on_curve(*point) for point in points)
    assert tuple(sorted(set(map(int, bad_primes)))) == tuple(bad_primes)
    assert all(ZZ(prime).is_prime(proof=True) for prime in bad_primes)
    global_rank, finite_row_count = finite_kummer_rank(certificate, len(points))

    cubic = two_division_cubic(coefficients)
    support = list(map(int, bad_primes))
    nf = pari.nfinit([pari(cubic), support])
    theta = pari.Mod(pari(cubic.parent().gen()), pari(cubic))
    kummer = [pari(4 * x) - theta for x, _y in points]
    beta_polynomials = [cubic.parent()([4 * x, -1]) for x, _y in points]
    for beta in kummer:
        norm = QQ(pari.nfeltnorm(nf, beta))
        assert norm.numerator().is_square() and norm.denominator().is_square()

    local_rows = [[] for _ in points]
    local_records = []
    for prime in support:
        local = LocalSquareclasses(nf, prime)
        signatures = [list(map(int, local.signature(beta))) for beta in kummer]
        assert len({len(signature) for signature in signatures}) == 1
        for row, signature in zip(local_rows, signatures):
            row.extend(signature)
        local_records.append(
            {
                "prime": prime,
                "prime_ideal_count": len(local.primes),
                "local_squareclass_dimension": len(signatures[0]),
                "cumulative_image_rank": int(matrix(GF(2), local_rows).rank()),
            }
        )

    real_roots = cubic.roots(AA, multiplicities=False)
    real_signatures = [
        [int(beta(root) < 0) for root in real_roots]
        for beta in beta_polynomials
    ]
    for row, signature in zip(local_rows, real_signatures):
        row.extend(signature)
    localization = matrix(GF(2), local_rows)
    image_rank = int(localization.rank())
    strict_basis = localization.left_kernel().basis_matrix()
    assert strict_basis.nrows() == len(points) - image_rank
    assert strict_basis * localization == 0
    return {
        "label": label,
        "public_point_count": len(points),
        "public_mod_2_rank": global_rank,
        "finite_certificate_row_count": finite_row_count,
        "two_division_cubic_ascending": list(map(str, cubic.list())),
        "two_division_cubic_irreducible": bool(cubic.is_irreducible()),
        "defining_order_index": str(ZZ(nf[3])),
        "field_discriminant": str(ZZ(nf.disc())),
        "field_discriminant_bit_length": abs(ZZ(nf.disc())).nbits(),
        "strict_bad_rational_primes": support,
        "finite_local_records": local_records,
        "real_embedding_count": len(real_roots),
        "localization_width": localization.ncols(),
        "known_public_local_image_rank": image_rank,
        "known_public_strict_local_kernel_dimension": len(points) - image_rank,
        "known_public_strict_local_kernel_basis": rows(strict_basis),
    }


def det1092_j_test():
    loader = runpy.run_path(str(PARENT_LOADER))["load_curve302_recovered_parent"]
    parent, _basis, specialization = loader()
    assert specialization == 0
    target = runpy.run_path(str(CURVE273))
    curve273 = EllipticCurve(QQ, target["GENERAL_WEIERSTRASS_COEFFICIENTS"])
    numerator = parent.base_ring().ring()((parent.j_invariant() - curve273.j_invariant()).numerator())
    assert numerator.degree() == 24
    roots = numerator.roots(QQ)
    assert not roots and numerator.is_irreducible()
    encoded = "\n".join(map(str, numerator.list())) + "\n"
    return {
        "status": "NO_RATIONAL_PARAMETER",
        "equation": "j_1092(t) = j_273",
        "numerator_degree": numerator.degree(),
        "numerator_irreducible_over_Q": True,
        "rational_roots": [],
        "numerator_coefficient_sha256": sha256(encoded.encode()).hexdigest(),
        "boundary": "A Q-isomorphic fibre, including any quadratic twist, would give a rational root. This excludes only rational parameters on the explicit recovered det-1092 parent.",
    }


def curve273_core():
    core = next(entry for entry in json.loads(CORES.read_text())["curves"] if entry["label"] == "curve273")
    heights = next(entry for entry in json.loads(HEIGHTS.read_text())["curves"] if entry["label"] == "curve273")
    source = runpy.run_path(str(CURVE273))
    curve = EllipticCurve(QQ, source["GENERAL_WEIERSTRASS_COEFFICIENTS"])
    public = [curve(QQ(x), QQ(y)) for x, y in source["POINTS"]]
    embedding = matrix(ZZ, core["saturated_basis_columns_in_public_point_coordinates"]).transpose()
    assert embedding.nrows() == 17 and embedding.ncols() == 30
    assert list(embedding.smith_form()[0].diagonal()) == [1] * 17
    ambient = matrix(RealField(192), heights["height_gram"])
    numerical_gram = embedding * ambient * embedding.transpose()
    basis = [
        sum((coefficient * point for coefficient, point in zip(word, public) if coefficient), curve(0))
        for word in embedding.rows()
    ]
    selected_count, words = select_integral_words(curve, basis, numerical_gram)
    kernel = quadratic_kernel(words, 17)
    assert selected_count == 1036 and len(words) == 450
    assert matrix(ZZ, words).rank() == 17
    assert kernel.nrows() == 1 and kernel[0, -1]
    recovered = gram_from_coordinates(4 * kernel[0] / kernel[0, -1], 17)
    assert all(entry in ZZ for entry in recovered.list())
    recovered = matrix(ZZ, recovered)
    assert recovered.is_positive_definite() and all(entry == 4 for entry in recovered.diagonal())
    assert recovered.det() == 1020
    assert list(recovered.smith_form()[0].diagonal()) == [1] * 16 + [1020]
    assert all(vector(ZZ, word) * recovered * vector(ZZ, word) == 4 for word in words)

    exact_shell = runpy.run_path(
        str(ROOT / "elkies-k3/scripts/certify_curve302_anchor6_triangle_gate.sage")
    )["exact_shell"]
    norm2, norm2_nodes = exact_shell(recovered, 2)
    norm4, norm4_nodes = exact_shell(recovered, 4)
    assert not norm2 and len(norm4) == 2526

    # A primitive embedding U+(-G) in the K3 lattice would require an even
    # ternary complement of signature (2,1), determinant -1020, and the same
    # discriminant form as G.  Sage's finite genus enumeration finds none.
    target_form = Genus(recovered).discriminant_form().normal_form()
    ternary_genera = genera((2, 1), 1020, even=True)
    compatible = [
        genus for genus in ternary_genera
        if genus.discriminant_form().normal_form() == target_form
    ]
    assert len(ternary_genera) == 8 and not compatible
    normal = target_form
    return {
        "status": "EXACT_FROZEN_SHELL_RECOVERY_NOT_A_K3_LATTICE_CANDIDATE",
        "public_embedding_is_primitive": True,
        "selected_ray_count": selected_count,
        "integral_ray_count": len(words),
        "integral_words_sha256": sha256(canonical_words(words).encode()).hexdigest(),
        "quadratic_constraint_rank": 153,
        "quadratic_kernel_dimension": kernel.nrows(),
        "recovered_gram": rows(recovered),
        "determinant": int(recovered.det()),
        "determinant_factorization": "2^2 * 3 * 5 * 17",
        "smith_factors": list(map(int, recovered.smith_form()[0].diagonal())),
        "minimum": 4,
        "signed_norm4_count": len(norm4),
        "exact_shell_nodes": {"norm2": norm2_nodes, "norm4": norm4_nodes},
        "discriminant_form_normal_key": {
            "invariants": list(map(int, normal.invariants())),
            "quadratic_gram": rows_as_strings(normal.gram_matrix_quadratic()),
            "value_module": str(normal.value_module_qf()),
        },
        "even_ternary_signature_2_1_genera_at_determinant": len(ternary_genera),
        "compatible_ternary_genera": 0,
        "primitive_K3_embedding": False,
        "K3_embedding_reason": "No even signature-(2,1) ternary genus has the required discriminant form, so U plus the negative recovered form cannot primitively embed in the K3 lattice.",
        "boundary": "The candidate 17-space and its ray selection are inherited from a bounded numerical-height search. The exact result concerns this frozen selection only; it does not exclude another 17-space, a different parent, or a non-K3 construction of curve 273.",
    }


def build():
    result = {
        "schema": "elliptic-curves.curve273-curve302-parent-fingerprints.v1",
        "status": "DET1092_PARENT_EXCLUDES_273_AND_PAIR_HAS_NO_MATCHING_TESTED_K3_OR_STRICT_FINGERPRINT",
        "input_sha256": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (
                Path(__file__), PARENT_LOADER, CURVE273, CURVE302, HEIGHTS,
                CORES, CERT273, CERT302,
            )
        },
        "det1092_j_recognition": det1092_j_test(),
        "curve273_rank17_core": curve273_core(),
        "strict_local_kummer_comparison": {
            "definition": "For the certified public mod-2 point span, kernel of localization in the product of K_v^*/K_v^{*2} over every listed bad rational prime and every real embedding. This is an exact strict-local kernel, not a complete strict Selmer group or a class-group computation.",
            "curve273": strict_local_kernel("curve273", CURVE273, BAD273, CERT273),
            "curve302": strict_local_kernel("curve302", CURVE302, BAD302, CERT302),
        },
        "boundary": "The direct parent test is exact and negative. The 273 lattice result is exact after a frozen numerical selection but does not identify a generic height pairing. The strict-local comparison is an exact calculation on the supplied public point spans; it neither computes the full S-class/strict Selmer groups nor proves rational solubility of any cover. No claim about the original discoverers' method is made.",
    }
    strict = result["strict_local_kummer_comparison"]
    assert strict["curve273"]["known_public_strict_local_kernel_dimension"] == 14
    assert strict["curve302"]["known_public_strict_local_kernel_dimension"] == 10
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", action="store_true")
    args = parser.parse_args()
    signal.alarm(120)
    result = build()
    rendered = json.dumps(result, indent=2, sort_keys=True, default=lambda value: int(value)) + "\n"
    if args.build:
        assert not OUT.exists()
        OUT.write_text(rendered)
    assert json.loads(rendered) == json.loads(OUT.read_text())
    print(result["status"], flush=True)
