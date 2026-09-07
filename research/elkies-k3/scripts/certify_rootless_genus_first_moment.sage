#!/usr/bin/env sage-python
"""Certify the first root moment of three positive rank-17 genera.

For an even lattice with quadratic form Q(x)=(x,x)/2, the coefficient of
q in its theta series is the signed root count.  Siegel's weighted genus
theta series therefore makes this coefficient a product of local densities.
This checker evaluates that product exactly for the determinant 78, 948, and
950 control genera.  It also recomputes the determinant-78 weighted average
from the complete 1,549-class census and checks equality with the local
formula.

The test ``average signed root count < 2`` is sufficient for a rootless
class, but is not necessary.  The checker deliberately records an
inconclusive result for all three controls; the explicit determinant-948 and
950 representatives still prove rootless existence in those two genera.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import (
    QQ,
    ZZ,
    QuadraticBernoulliNumber,
    QuadraticForm,
    bernoulli,
    fundamental_discriminant,
    kronecker_symbol,
    matrix,
    pari,
    prime_divisors,
)
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = Path(__file__).resolve()
DET78_GRAM = ROOT / "elkies-k3/data/lattice/e6_rank4_det78_frame.txt"
DET948_GRAM = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
DET950_SOURCE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-ns0024-new-rootless-source-route-v1.json"
)
DET78_CENSUS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-e6-rank4-det78-niemeier-frames-v1.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-rootless-genus-first-moment-v1.json"
)


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rational(value):
    value = QQ(value)
    return {
        "denominator": int(value.denominator()),
        "numerator": int(value.numerator()),
        "text": str(value),
    }


def load_matrix(path):
    rows = [
        [ZZ(entry) for entry in line.split()]
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    return matrix(ZZ, rows)


def quadratic_form(gram):
    form = QuadraticForm(ZZ, gram)
    assert form.Hessian_matrix() == gram
    return form


def signed_root_count(gram):
    data = pari(gram).qfminim(2)
    return int(data[0])


def siegel_average_root_count(gram):
    """Return the exact q^1 coefficient of the weighted genus theta series.

    For positive odd rank n=2k+1 with k even and u=1, write

        D = fund_disc((-1)^k 2 det(L)).

    The generic Euler product is expressed with B_{k,chi_D}; at each prime
    dividing 2 det(L), its generic factor is replaced by the exact local
    density alpha_p(1).  The three controls have square rational real-place
    factors.  We keep all factors explicit in the output.
    """

    rank = gram.nrows()
    determinant = ZZ(gram.det())
    assert rank == 17
    assert determinant > 0
    assert gram == gram.transpose()
    assert all(gram[index, index] % 2 == 0 for index in range(rank))

    k = (rank - 1) // 2
    assert k % 2 == 0
    discriminant = fundamental_discriminant(((-1) ** k) * 2 * determinant)
    conductor = abs(discriminant)
    real_square = QQ(2 * conductor) / determinant
    assert real_square.is_square()
    real_factor = real_square.sqrt()
    generalized_bernoulli = QuadraticBernoulliNumber(k, discriminant)
    ordinary_bernoulli = bernoulli(2 * k)

    generic_product = (
        QQ(2 ** (2 * k) * abs(generalized_bernoulli))
        / (conductor**k * abs(ordinary_bernoulli))
        * real_factor
    )

    form = quadratic_form(gram)
    local_corrections = []
    answer = generic_product
    for prime in prime_divisors(2 * determinant):
        prime = ZZ(prime)
        character = kronecker_symbol(discriminant, prime)
        generic_euler_factor = (1 - prime ** (-2 * k)) / (
            1 - character * prime ** (-k)
        )
        local_density = form.local_density(prime, 1)
        correction = local_density / generic_euler_factor
        answer *= correction
        local_corrections.append(
            {
                "character_value": int(character),
                "correction": rational(correction),
                "generic_euler_factor": rational(generic_euler_factor),
                "local_density_at_one": rational(local_density),
                "prime": int(prime),
            }
        )

    return answer, {
        "bad_prime_corrections": local_corrections,
        "fundamental_discriminant": int(discriminant),
        "generalized_bernoulli": rational(generalized_bernoulli),
        "generic_product_before_bad_prime_corrections": rational(generic_product),
        "half_rank_parameter": int(k),
        "ordinary_bernoulli": rational(ordinary_bernoulli),
        "real_factor": rational(real_factor),
        "real_factor_square": rational(real_square),
    }


def det78_enumeration_crosscheck(expected_average):
    source = json.loads(DET78_CENSUS.read_text())
    frames = source["frames"]
    mass = sum(
        (QQ(1) / ZZ(frame["automorphism_group_order"]) for frame in frames),
        QQ(0),
    )
    first_root_moment = sum(
        (
            QQ(frame["signed_root_count"])
            / ZZ(frame["automorphism_group_order"])
            for frame in frames
        ),
        QQ(0),
    )
    average = first_root_moment / mass
    asserted_mass = QQ(source["accounting"]["target_genus_mass"])
    assert len(frames) == 1549
    assert source["accounting"]["genus_mass_closed"] is True
    assert source["accounting"]["rootless_frame_classes"] == 0
    assert all(int(frame["signed_root_count"]) >= 2 for frame in frames)
    assert mass == asserted_mass
    assert average == expected_average
    return {
        "all_classes_rootful": True,
        "class_count": len(frames),
        "genus_mass": rational(mass),
        "local_formula_equals_enumerated_average": True,
        "weighted_first_root_moment": rational(first_root_moment),
        "weighted_mean_signed_root_count": rational(average),
    }


def case_record(name, gram, gram_source, known_rootless, evidence):
    average, formula = siegel_average_root_count(gram)
    roots = signed_root_count(gram)
    assert bool(roots == 0) == known_rootless
    return {
        "control": name,
        "determinant": int(gram.det()),
        "first_moment_criterion": {
            "conclusion": (
                "rootless_class_exists"
                if average < 2
                else "inconclusive_average_not_below_two"
            ),
            "mean_is_strictly_below_two": bool(average < 2),
            "threshold": 2,
        },
        "gram_source": relative(gram_source),
        "known_rootless_representative": known_rootless,
        "known_status_evidence": evidence,
        "rank": gram.nrows(),
        "representative_signed_root_count": roots,
        "siegel_formula": formula,
        "weighted_mean_signed_root_count": rational(average),
        "weighted_mean_signed_root_count_decimal": f"{float(average):.15f}",
    }, average


def build():
    det78 = load_matrix(DET78_GRAM)
    det948 = load_matrix(DET948_GRAM)
    det950_source = json.loads(DET950_SOURCE.read_text())
    det950 = matrix(ZZ, det950_source["new_rootless_frame"]["gram"])

    case78, average78 = case_record(
        "det78_complete_rootful_genus",
        det78,
        DET78_GRAM,
        False,
        "complete 1,549-class mass-closed census",
    )
    case948, _ = case_record(
        "det948_published_R17_rootless_genus",
        det948,
        DET948_GRAM,
        True,
        "explicit rootless Gram matrix",
    )
    case950, _ = case_record(
        "det950_NS0024_rootless_genus",
        det950,
        DET950_SOURCE,
        True,
        "explicit rootless Gram matrix in the NS0024 route artifact",
    )
    assert [case["determinant"] for case in (case78, case948, case950)] == [
        78,
        948,
        950,
    ]
    assert not any(
        case["first_moment_criterion"]["mean_is_strictly_below_two"]
        for case in (case78, case948, case950)
    )

    return {
        "cases": [case78, case948, case950],
        "det78_independent_enumeration_crosscheck": det78_enumeration_crosscheck(
            average78
        ),
        "formula": {
            "convention": "Q(x)=(x,x)/2 and theta_L(q)=sum_x q^Q(x)",
            "criterion": "[q]Theta_G<2 implies positive rootless mass",
            "formula_rank_2k_plus_1_k_even": (
                "2^(2k)*abs(B_(k,chi_D))/(abs(D)^k*abs(B_(2k)))"
                "*sqrt(2*abs(D)/det(L))*product_(p|2det(L))"
                "(alpha_p(1)/((1-p^(-2k))/(1-chi_D(p)*p^(-k))))"
            ),
            "meaning": (
                "The value is the mass-normalized weighted average of the signed "
                "root count over the positive even genus."
            ),
        },
        "inputs": {
            relative(path): "sha256:" + digest(path)
            for path in (DET78_GRAM, DET948_GRAM, DET950_SOURCE, DET78_CENSUS)
        },
        "proof_boundary": (
            "This artifact certifies only the degree-one Siegel root moment and "
            "the <2 sufficient criterion. It does not perform the higher-degree "
            "ADE representation-mass inversion and therefore does not decide "
            "rootlessness when the displayed mean is at least two."
        ),
        "reproduce": {
            "check": (
                "sage -python elkies-k3/scripts/"
                "certify_rootless_genus_first_moment.sage --check"
            ),
            "generate": (
                "sage -python elkies-k3/scripts/"
                "certify_rootless_genus_first_moment.sage"
            ),
        },
        "sage_version": SAGE_VERSION,
        "schema": "elkies-k3-rootless-genus-first-moment-v1",
        "script": relative(SCRIPT),
        "status": "PASS_EXACT_LOCAL_FIRST_MOMENTS_AND_DET78_CENSUS_CROSSCHECK",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    rendered = json.dumps(build(), indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != rendered:
            raise SystemExit(f"stale artifact: run {relative(SCRIPT)}")
        print(f"OK: {relative(output)}")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(rendered)
    print(f"wrote {relative(output)}")


if __name__ == "__main__":
    main()
