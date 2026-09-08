#!/usr/bin/env sage -python
"""Enumerate the canonical degree-two NS quotient of the curve-302 parent.

The recovered parent has geometric Neron--Severi lattice

    NS(X) = U + (-M),   M = MW(X),   rank(M)=17, det(M)=1092,

and all of these classes are rational.  For a class

    C = d O + b F + phi(w)

section translation by ``x in M`` sends ``w`` to ``w+d*x``.  Thus for a
fixed degree the lattice problem is finite only *after* quotienting by
``M/dM``.  At degree two this script visits all 2^17 quotient classes.

For a coset ``c`` let mu(c) be its minimum M-norm.  The degree-two classes
which are nonnegative on every section are exactly:

* rational candidates: mu(c) == 2 (mod 4), mu(c) >= 10;
* genus-one candidates: mu(c) == 0 (mod 4), mu(c) >= 8, c != 0.

Every rational survivor is an actual geometrically irreducible smooth
rational bisection: the parent has only 24 I1 fibres, so the standard
rootless-fibration argument rules out both vertical and two-section
decompositions.  The genus-one rows are only lattice/divisor candidates;
nefness and irreducibility have not been asserted for them.

This is deliberately the lattice gate, not an equation constructor.  A
later equation stage must construct a divisor in a listed orbit, factor its
fibre at t=0, and map every rational factor point to the displayed
``D/sp(M) = Z^14`` quotient.  It must not call a nonsplit lattice orbit a
miss before that construction exists.

Two independent floating GSO precisions drive complete CVP enumerations;
every returned representative and norm is then checked with the integral
Gram matrix.  The cross-precision equality is a numerical-audit guard, not
a replacement for the integral checks.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from sage.all import ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
PARENT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/curve302_recovered_mw17_parent_v1.json"
)
PARENT_PROOF = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/curve302_recovered_mw17_parent_proof_v1.json"
)
GEOMETRIC_PROOF = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/curve302_parent_geometric_picard19_v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "curve302_parent_degree2_multisection_lattice_v1.json"
)
DEFAULT_ORBITS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "curve302_parent_degree2_multisection_orbits_v1.tsv"
)
DIMENSION = 17
DEGREE = 2


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def norm(row, gram):
    return int(row * gram * row)


def mask(row) -> int:
    return sum((int(value) & 1) << index for index, value in enumerate(row))


def vector_text(row) -> str:
    return " ".join(str(int(value)) for value in row)


class CosetOracle:
    """Complete CVP search with exact integral output checks.

    ``Enumeration.enumerate`` is called with no pruning strategy and a ball
    known to contain the digit representative.  Its candidate coordinates,
    residue and norm are independently checked over the integers.
    """

    def __init__(self, gram, float_type: str, precision: int) -> None:
        from fpylll import FPLLL, GSO, IntegerMatrix

        if float_type == "mpfr":
            FPLLL.set_precision(precision)
        self.gram = gram
        self.float_type = float_type
        self.precision = precision
        self.gso = GSO.Mat(
            IntegerMatrix.from_matrix([[int(value) for value in row] for row in gram.rows()]),
            gram=True,
            float_type=float_type,
            update=True,
        )
        self.mu = [
            [self.gso.get_mu(i, j) if i > j else 0.0 for j in range(DIMENSION)]
            for i in range(DIMENSION)
        ]
        # The digit vector is an available representative in every coset.
        # This deliberately generous ball prevents a false empty CVP search.
        self.distance_bound = (
            sum(abs(int(value)) for row in gram.rows() for value in row) / 4.0
            + 1.0
        )

    def solve(self, residue: tuple[int, ...]):
        from fpylll import Enumeration

        target = [
            -(
                residue[i]
                + sum(residue[j] * self.mu[j][i] for j in range(i + 1, DIMENSION))
            )
            / DEGREE
            for i in range(DIMENSION)
        ]
        solutions = Enumeration(self.gso).enumerate(
            0, DIMENSION, self.distance_bound, 0, target=target
        )
        if not solutions:
            raise ArithmeticError("CVP enumeration returned no solution")
        reported_distance, coordinates = solutions[0]
        closest = tuple(int(round(entry)) for entry in coordinates)
        if any(abs(entry - integer) > 1.0e-7 for entry, integer in zip(coordinates, closest)):
            raise ArithmeticError("CVP returned nonintegral coordinates")
        representative = vector(
            ZZ, [residue[index] + DEGREE * closest[index] for index in range(DIMENSION)]
        )
        if any((representative[index] - residue[index]) % DEGREE for index in range(DIMENSION)):
            raise ArithmeticError("CVP returned the wrong residue class")
        exact = norm(representative, self.gram)
        error = abs(DEGREE * DEGREE * reported_distance - exact)
        if error > 1.0e-7 or exact < 0 or exact % 2:
            raise ArithmeticError(f"invalid CVP output: norm={exact}, error={error}")
        return exact, representative, error


def load_inputs():
    parent = json.loads(PARENT.read_text())
    parent_proof = json.loads(PARENT_PROOF.read_text())
    geometric_proof = json.loads(GEOMETRIC_PROOF.read_text())
    if parent["status"] != "EXPLICIT_MODEL_AND_BASIS_INPUTS":
        raise ArithmeticError("unexpected parent input status")
    if parent_proof["status"] != "PASS_FULL_ARITHMETIC_MW17_PARENT":
        raise ArithmeticError("the arithmetic MW17 proof is not available")
    if geometric_proof["status"] != "PASS":
        raise ArithmeticError("the geometric Picard-19 proof is not available")
    if (
        parent_proof["generic_arithmetic_MW_rank"] != 17
        or parent_proof["height_determinant"] != 1092
        or geometric_proof["geometric_Picard_rank"] != 19
        or geometric_proof["geometric_generic_MW_rank"] != 17
        or not geometric_proof["full_geometric_basis_is_displayed_rational_basis"]
    ):
        raise ArithmeticError("the required full geometric MW17 hypotheses changed")
    gram = matrix(ZZ, parent["generic_height_gram"])
    if (
        gram.nrows() != DIMENSION
        or gram.ncols() != DIMENSION
        or not gram.is_positive_definite()
        or gram.det() != 1092
        or any(gram[index, index] % 2 for index in range(DIMENSION))
    ):
        raise ArithmeticError("unexpected even positive-definite MW Gram matrix")
    return parent, gram


def enumerate_degree_two(gram, *, audit_precision: int):
    # A row-coordinate LLL change makes the CVP enumeration inexpensive.  All
    # exported vectors are transported back to the certified parent basis.
    reduced_to_parent = gram.LLL_gram().transpose()
    reduced = reduced_to_parent * gram * reduced_to_parent.transpose()
    if abs(reduced_to_parent.det()) != 1 or reduced.det() != gram.det():
        raise ArithmeticError("invalid LLL coordinate change")
    primary = CosetOracle(reduced, "dd", 160)
    audit = CosetOracle(reduced, "mpfr", audit_precision)

    minima = Counter()
    representative_by_mask = {}
    multiplicity_by_mask = Counter()
    maximum_primary_error = 0.0
    maximum_audit_error = 0.0
    for identifier in range(1 << DIMENSION):
        residue = tuple((identifier >> index) & 1 for index in range(DIMENSION))
        primary_norm, primary_vector, primary_error = primary.solve(residue)
        audit_norm, audit_vector, audit_error = audit.solve(residue)
        if primary_norm != audit_norm:
            raise ArithmeticError(
                f"cross-precision CVP mismatch at mask {identifier}: "
                f"{primary_norm} != {audit_norm}"
            )
        if norm(audit_vector, reduced) != audit_norm:
            raise ArithmeticError("MPFR representative has an incorrect exact norm")
        if norm(primary_vector, reduced) != primary_norm:
            raise ArithmeticError("DD representative has an incorrect exact norm")
        if mask(primary_vector) != identifier or mask(audit_vector) != identifier:
            raise ArithmeticError("representative parity does not equal its mask")
        minima[primary_norm] += 1
        maximum_primary_error = max(maximum_primary_error, primary_error)
        maximum_audit_error = max(maximum_audit_error, audit_error)
        # The two minimizers can differ.  Retain a deterministic exact
        # representative in the reduced basis, then map it back only once.
        choices = sorted((tuple(map(int, primary_vector)), tuple(map(int, audit_vector))))
        representative_by_mask[identifier] = vector(ZZ, choices[0])
        multiplicity_by_mask[identifier] = int(primary_norm)
        if identifier and identifier % 16384 == 0:
            print(
                "CURVE302D2LATTICE|covered={}|total={}|status=RUNNING".format(
                    identifier, 1 << DIMENSION
                ),
                flush=True,
            )

    if sum(minima.values()) != 1 << DIMENSION:
        raise ArithmeticError("incomplete degree-two coset coverage")
    return {
        "reduced_gram": reduced,
        "reduced_to_parent": reduced_to_parent,
        "minima": minima,
        "representative_by_mask": representative_by_mask,
        "minimum_by_mask": multiplicity_by_mask,
        "maximum_primary_error": maximum_primary_error,
        "maximum_audit_error": maximum_audit_error,
    }


def classify(records):
    minimum_by_mask = records["minimum_by_mask"]
    rational = [
        item
        for item, minimum in minimum_by_mask.items()
        if minimum % 4 == 2 and minimum >= 10
    ]
    genus_one = [
        item
        for item, minimum in minimum_by_mask.items()
        if item != 0 and minimum % 4 == 0 and minimum >= 8
    ]
    if any(minimum_by_mask[item] < 10 for item in rational):
        raise ArithmeticError("rational section-nonnegativity filter failed")
    if any(minimum_by_mask[item] < 8 for item in genus_one):
        raise ArithmeticError("genus-one section-nonnegativity filter failed")
    return rational, genus_one


def write_orbits(path, records, rational, genus_one):
    reduced_to_parent = records["reduced_to_parent"]
    parent_gram = (
        reduced_to_parent.inverse()
        * records["reduced_gram"]
        * reduced_to_parent.inverse().transpose()
    )
    representatives = records["representative_by_mask"]
    minima = records["minimum_by_mask"]
    category = {item: "rational" for item in rational}
    category.update({item: "genus_one" for item in genus_one})
    with path.open("w") as stream:
        stream.write(
            "orbit_mask\tcategory\tminimum_norm\treduced_basis_w\t"
            "parent_MW17_w\tlattice_class_at_minimum\n"
        )
        for identifier in sorted(category):
            reduced_vector = representatives[identifier]
            parent_vector = reduced_vector * reduced_to_parent
            if norm(parent_vector, parent_gram) != minima[identifier]:
                raise ArithmeticError("transported representative norm changed")
            genus = 0 if category[identifier] == "rational" else 1
            # C = d O + b F + phi(w), with C.O=0 for the chosen minimum
            # representative only at the threshold shell.  Deep classes use
            # the exact formula shown, rather than being silently dropped.
            b_numerator = (
                2 * DEGREE * DEGREE + minima[identifier] + 2 * genus - 2
            )
            if b_numerator % (2 * DEGREE):
                raise ArithmeticError("lattice class has nonintegral F coefficient")
            b = b_numerator // (2 * DEGREE)
            stream.write(
                "{}\t{}\t{}\t{}\t{}\t{}O+{}F+phi(w)\n".format(
                    identifier,
                    category[identifier],
                    minima[identifier],
                    vector_text(reduced_vector),
                    vector_text(parent_vector),
                    DEGREE,
                    b,
                )
            )


def build(arguments):
    parent, gram = load_inputs()
    records = enumerate_degree_two(gram, audit_precision=arguments.audit_precision)
    rational, genus_one = classify(records)
    arguments.orbits_output.parent.mkdir(parents=True, exist_ok=True)
    write_orbits(arguments.orbits_output, records, rational, genus_one)
    minima = records["minima"]
    if min(minima) != 0 or minima[0] != 1:
        raise ArithmeticError("zero coset was not unique")
    payload = {
        "schema": "elliptic-curves.curve302-parent-degree2-multisection-lattice.v1",
        "status": "PASS_COMPLETE_DEGREE2_TRANSLATION_QUOTIENT",
        "scope": (
            "All M/2M translation orbits for the full geometric MW17 lattice of the "
            "explicit curve-302 parent. Rational survivors are geometrically irreducible "
            "smooth bisections. Genus-one rows are lattice/divisor candidates only. No "
            "bisection equation, t=0 splitting statement, quotient direction, or rank "
            "gain is inferred until a separate equation-and-specialization certificate."
        ),
        "inputs": {
            relative(PARENT): digest(PARENT),
            relative(PARENT_PROOF): digest(PARENT_PROOF),
            relative(GEOMETRIC_PROOF): digest(GEOMETRIC_PROOF),
        },
        "parent": {
            "geometric_NS": "U + (-MW17)",
            "MW_rank": 17,
            "MW_determinant": 1092,
            "geometric_Picard_rank": 19,
            "all_geometric_MW_classes_rational": True,
            "geometric_fibres": "24 I1",
            "specialization_parameter": parent["specialization_parameter"],
            "specialization_model": "literal curve302",
        },
        "normalization": {
            "class": "C=dO+bF+phi(w)",
            "intersection_form": "C^2=-2d^2+2db-<w,w>",
            "section_translation": "w -> w+d*x",
            "degree": DEGREE,
            "finite_quotient": "M/2M",
            "total_translation_orbits": 1 << DIMENSION,
            "minimum_norm": "mu_2(c)=min{<w,w>:w in c}",
            "all_section_intersection_formula": "min_x C.S_x=mu_2(c)/4-5/2+g/2",
        },
        "enumeration": {
            "method": (
                "complete CVP over all 2^17 binary residue vectors in an integral LLL basis; "
                "each result is repeated with double-double and MPFR GSO arithmetic and "
                "its representative/norm are checked integrally"
            ),
            "row_LLL_change_to_certified_parent_basis": [
                list(map(int, row)) for row in records["reduced_to_parent"].rows()
            ],
            "reduced_gram": [list(map(int, row)) for row in records["reduced_gram"].rows()],
            "minimum_norm_histogram": {
                str(value): int(count) for value, count in sorted(minima.items())
            },
            "maximum_dd_distance_rounding_error": records["maximum_primary_error"],
            "maximum_mpfr_distance_rounding_error": records["maximum_audit_error"],
            "mpfr_precision_bits": arguments.audit_precision,
            "complete": True,
        },
        "rational_bisections": {
            "arithmetic_genus": 0,
            "norm_condition": "mu_2(c) == 2 (mod 4)",
            "section_nonnegative_threshold": 10,
            "translation_orbits": len(rational),
            "geometric_conclusion": (
                "Every listed class is effective over QQ and is a geometrically irreducible "
                "smooth rational bisection. This uses the full rational NS lattice, 24 I1 "
                "fibres, and nonnegative intersection with every section; it does not give "
                "an explicit equation in the original Weierstrass chart."
            ),
        },
        "genus_one_bisection_candidates": {
            "arithmetic_genus": 1,
            "norm_condition": "mu_2(c) == 0 (mod 4), c != 0",
            "section_nonnegative_threshold": 8,
            "translation_orbits": len(genus_one),
            "geometric_boundary": (
                "Riemann--Roch gives an effective rational divisor class, but the present "
                "lattice filter does not certify nefness or irreducibility of every genus-one row."
            ),
        },
        "next_equation_gate": {
            "required_per_orbit": [
                "construct a divisor in the exported NS orbit over QQ",
                "verify its generic irreducibility and its degree/genus",
                "factor its fibre at t=0 exactly",
                "transport each rational fibre point into D/sp(MW17)=Z^14",
                "certify independence from the specialized generic image before any rank claim",
            ],
            "prior_bespoke_control": (
                "The 178 conic/three-line twisted-cubic bisections were all nonsplit at t=0; "
                "they are not the complete M/2M quotient."
            ),
        },
        "orbits_tsv": relative(arguments.orbits_output),
        "orbits_tsv_sha256": digest(arguments.orbits_output),
        "reproducing_command": (
            "sage -python elliptic-curves/cas/"
            "enumerate_curve302_parent_degree2_multisections.py --check"
        ),
    }
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--orbits-output", type=Path, default=DEFAULT_ORBITS)
    parser.add_argument("--audit-precision", type=int, default=160)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    if arguments.audit_precision < 80:
        parser.error("--audit-precision must be at least 80")
    if not arguments.check and (
        arguments.output.exists() or arguments.orbits_output.exists()
    ):
        raise FileExistsError("refusing to overwrite an existing generated certificate")
    payload = build(arguments)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    if arguments.check:
        if not arguments.output.exists():
            raise FileNotFoundError(arguments.output)
        stored = json.loads(arguments.output.read_text())
        if stored != payload:
            raise ArithmeticError("the stored degree-two certificate does not replay")
    else:
        arguments.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        "CURVE302D2LATTICE|rational_orbits={}|genus1_candidates={}|"
        "cosets={}|status=PASS_COMPLETE_DEGREE2_TRANSLATION_QUOTIENT".format(
            payload["rational_bisections"]["translation_orbits"],
            payload["genus_one_bisection_candidates"]["translation_orbits"],
            1 << DIMENSION,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
