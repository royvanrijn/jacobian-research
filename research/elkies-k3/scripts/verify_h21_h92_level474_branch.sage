#!/usr/bin/env sage -python
"""Identify the level-474 component of H_21 pulled back to the H_92 chart.

This is an exact characteristic-11 computation.  It consumes the pinned
Gruenewald H_21 Satake polynomial and the Elkies--Kumar H_92 model, expands
the pullback by a sparse Horner scheme, factors it, labels the six components
through the CM-24 point by their tangent directions, and compares the target
component's normalized point counts with the published level-474 genus-two
model.

The computation identifies the correct component modulo 11.  It does not by
itself construct a characteristic-zero birational map.
"""

from sage.all import GF, Infinity, Integer, PolynomialRing, QQ, ZZ, factor, gcd, sage_eval

import argparse
import hashlib
import json
from pathlib import Path
import re
import time


PRIME = 11
H21_SHA256 = "927a5e5c6dffa0f3fa5c386c40fdd6a1895cd3f62b8181a75d8049f1695f6d34"
H92_SHA256 = "427559ecd4c2c19d0a4ed7df1019c8a351ed34f454691e9ef1080a8834e74ea1"
H21_URL = "https://www.maths.usyd.edu.au/u/davidg/ThesisData/SatakeHumbert/level1_21.txt"
EK_SOURCE_URL = "https://export.arxiv.org/e-print/1209.3527"

BRANCH_SLOPES_QQ = (
    "(3*w+5)/7",
    "(-21*w+5)/19",
    "w-1",
    "(11*w-35)/31",
    "(19*w-99)/91",
    "(-7*w-33)/37",
)
EXPECTED_BRANCH_DEGREES = {1: 25, 2: 30, 3: 2, 4: 8, 5: 20, 6: 21}
EXPECTED_FACTOR_DEGREES = (
    2, 8, 10, 11, 13, 14, 18, 20, 21, 23, 24, 25, 28, 30, 31, 33,
    34, 39, 44, 44, 55, 57,
)


def stage(name, **values):
    payload = "|".join(f"{key}={value}" for key, value in values.items())
    print(f"H21H92|stage={name}" + (f"|{payload}" if payload else ""), flush=True)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_input(path, expected, label):
    actual = digest(path)
    if actual != expected:
        raise ValueError(
            f"{label} SHA-256 mismatch: expected {expected}, got {actual} ({path})"
        )
    return actual


def extract_h21(path):
    ring = PolynomialRing(ZZ, names=("s2", "s3", "s5", "s6"))
    text = path.read_text()
    if "symm21:=" not in text:
        raise ValueError("H21 input has no symm21 assignment")
    expression = " ".join(text.split("symm21:=", 1)[1].rsplit(";", 1)[0].split())
    polynomial = ring(expression)
    content = polynomial.content()
    if not content:
        raise ValueError("H21 polynomial is zero")
    primitive = polynomial // content
    return primitive, content


def extract_h92(path):
    text = path.read_text()
    ring = PolynomialRing(QQ, names=("r", "s"))
    environment = {str(generator): generator for generator in ring.gens()}
    values = []
    for name in ("A1", "A", "B1", "B", "B2"):
        match = re.search(rf"\b{name}\s*=\s*(.*?);", text, flags=re.S)
        if match is None:
            raise ValueError(f"H92 input has no {name} assignment")
        expression = re.sub(r"\s+", " ", match.group(1).replace("^", "**")).strip()
        values.append(ring(sage_eval(expression, locals=environment)))
    return ring, tuple(values)


def sparse_horner(terms, power_tables, target_ring, order, level=0):
    if level == len(order):
        return target_ring(sum(terms.values()))
    variable = order[level]
    groups = {}
    for exponent, coefficient in terms.items():
        exponent = tuple(exponent)
        groups.setdefault(exponent[variable], {})[exponent] = coefficient
    exponents = sorted(groups, reverse=True)
    previous = exponents[0]
    result = sparse_horner(groups[previous], power_tables, target_ring, order, level + 1)
    for exponent in exponents[1:]:
        result = (
            result * power_tables[variable][previous - exponent]
            + sparse_horner(groups[exponent], power_tables, target_ring, order, level + 1)
        )
        previous = exponent
    if previous:
        result *= power_tables[variable][previous]
    return result


def build_pullback(h21_primitive, h92_coefficients):
    field = GF(PRIME)
    satake_ring = PolynomialRing(field, names=("s2", "s3", "s5", "s6"))
    h21 = satake_ring(h21_primitive)
    target_ring = PolynomialRing(field, names=("r", "s"), order="lex")
    A1, A, B1, B, B2 = (target_ring(value) for value in h92_coefficients)
    values = (
        -36 * A,
        162 * B,
        -2430 * (A * B + 2 * A1 * B2),
        -2916 * A**3 + 4374 * B**2 + 17496 * B1 * B2,
    )
    order = (0, 2, 1, 3)
    power_tables = []
    for variable, value in enumerate(values):
        maximum = max(tuple(exponent)[variable] for exponent in h21.dict())
        power_tables.append(tuple(value**exponent for exponent in range(maximum + 1)))
    start = time.monotonic()
    pullback = sparse_horner(h21.dict(), power_tables, target_ring, order)
    stage(
        "pullback",
        seconds=f"{time.monotonic() - start:.2f}",
        degree=pullback.total_degree(),
        terms=len(pullback.dict()),
    )
    return pullback


def monomial_valuations(polynomial):
    exponents = tuple(tuple(exponent) for exponent in polynomial.dict())
    return tuple(min(exponent[index] for exponent in exponents) for index in range(2))


def label_cm24_components(factorization):
    base = GF(PRIME)
    univariate = PolynomialRing(base, "Z")
    Z = univariate.gen()
    extension = GF(PRIME**2, "w", modulus=Z**2 - Z + 1)
    w = extension.gen()
    local_ring = PolynomialRing(extension, names=("x", "y"))
    x, y = local_ring.gens()
    slope_ring = PolynomialRing(extension, "z")
    z = slope_ring.gen()
    roots = (
        (3 * w + 5) / 7,
        (-21 * w + 5) / 19,
        w - 1,
        (11 * w - 35) / 31,
        (19 * w - 99) / 91,
        (-7 * w - 33) / 37,
    )
    labels = {}
    for component, component_multiplicity in factorization:
        over_extension = component.change_ring(extension)
        if over_extension(w, -w):
            continue
        translated = local_ring(over_extension(w + x, -w + y))
        order = min(sum(exponent) for exponent in translated.dict())
        initial = local_ring(
            {
                tuple(exponent): coefficient
                for exponent, coefficient in translated.dict().items()
                if sum(exponent) == order
            }
        )
        tangent = slope_ring(initial(z, 1))
        matches = tuple(index for index, root in enumerate(roots, 1) if tangent(root) == 0)
        if len(matches) != 1:
            raise AssertionError(
                f"CM24 component of degree {component.total_degree()} has matches {matches}"
            )
        branch = matches[0]
        if branch in labels:
            raise AssertionError(f"branch {branch} occurs more than once")
        labels[branch] = {
            "polynomial": component,
            "degree": int(component.total_degree()),
            "terms": len(component.dict()),
            "factor_multiplicity": int(component_multiplicity),
            "intersection_multiplicity": int(order),
            "tangent_polynomial_mod_11": str(factor(tangent)),
            "slope_characteristic_zero": BRANCH_SLOPES_QQ[branch - 1],
            "slope_mod_11": str(roots[branch - 1]),
        }
    actual_degrees = {branch: data["degree"] for branch, data in labels.items()}
    if actual_degrees != EXPECTED_BRANCH_DEGREES:
        raise AssertionError(
            f"unexpected branch degrees: {actual_degrees}"
        )
    return labels


def map_polynomial(polynomial, target_field, embedding=None):
    source_ring = polynomial.parent()
    target_ring = PolynomialRing(target_field, names=source_ring.variable_names())
    convert = embedding if embedding is not None else target_field
    return target_ring(
        {
            tuple(exponent): convert(coefficient)
            for exponent, coefficient in polynomial.dict().items()
        }
    )


def local_multiplicity(polynomial):
    return min(sum(exponent) for exponent in polynomial.dict()) if polynomial else Infinity


def initial_form(polynomial, order):
    ring = polynomial.parent()
    return ring(
        {
            tuple(exponent): coefficient
            for exponent, coefficient in polynomial.dict().items()
            if sum(exponent) == order
        }
    )


def split_roots(polynomial, name):
    """Yield (field, embedding_from_base, geometric_root)."""
    base = polynomial.base_ring()
    for irreducible, unused_multiplicity in polynomial.factor():
        if irreducible.degree() == 1:
            for root, unused_root_multiplicity in irreducible.roots():
                yield base, None, root
        else:
            split_field, embedding = irreducible.splitting_field(name, map=True)
            split_ring = PolynomialRing(split_field, polynomial.variable_name())
            split_polynomial = split_ring(
                {
                    (exponent if isinstance(exponent, (int, Integer)) else exponent[0]):
                    embedding(coefficient)
                    for exponent, coefficient in irreducible.dict().items()
                }
            )
            for root, unused_root_multiplicity in split_polynomial.roots():
                yield split_field, embedding, root


def resolve_local(polynomial, depth=0):
    ring = polynomial.parent()
    field = ring.base_ring()
    x, y = ring.gens()
    order = local_multiplicity(polynomial)
    if order < 2:
        return ZZ(0)
    contribution = order * (order - 1) // 2
    tangent = initial_form(polynomial, order)
    univariate = PolynomialRing(field, "u")
    u = univariate.gen()
    finite = univariate(tangent(1, u))

    for irreducible, tangent_multiplicity in finite.factor():
        if tangent_multiplicity < 2:
            continue
        if irreducible.degree() == 1:
            roots = irreducible.roots()
            polynomial_over_split = polynomial
        else:
            split_field, embedding = irreducible.splitting_field(f"a{depth}", map=True)
            split_ring = PolynomialRing(split_field, names=ring.variable_names())
            polynomial_over_split = split_ring(
                {
                    tuple(exponent): embedding(coefficient)
                    for exponent, coefficient in polynomial.dict().items()
                }
            )
            split_univariate = PolynomialRing(split_field, "u")
            irreducible_split = split_univariate(
                {
                    (exponent if isinstance(exponent, (int, Integer)) else exponent[0]):
                    embedding(coefficient)
                    for exponent, coefficient in irreducible.dict().items()
                }
            )
            roots = irreducible_split.roots()
        xx, yy = polynomial_over_split.parent().gens()
        for root, unused_root_multiplicity in roots:
            transformed = polynomial_over_split(xx, xx * (root + yy)) // xx**order
            if local_multiplicity(transformed) >= 2:
                contribution += resolve_local(transformed, depth + 1)

    vertical = univariate(tangent(u, 1))
    if vertical.valuation() >= 2:
        transformed = polynomial(x * y, y) // y**order
        if local_multiplicity(transformed) >= 2:
            contribution += resolve_local(transformed, depth + 1)
    return ZZ(contribution)


def normalization_genus(polynomial):
    field = polynomial.base_ring()
    ring = polynomial.parent()
    r, s = ring.gens()
    degree = polynomial.total_degree()
    homogeneous_ring = PolynomialRing(field, names=("R", "S", "Z"))
    R, S, Z = homogeneous_ring.gens()
    homogeneous = homogeneous_ring(polynomial).homogenize(Z)

    singular_ideal = ring.ideal(
        (polynomial, polynomial.derivative(r), polynomial.derivative(s))
    )
    groebner = singular_ideal.groebner_basis()
    eliminants = [candidate for candidate in groebner if candidate.degree(r) == 0]
    if not eliminants:
        raise AssertionError("affine singular ideal has no s-eliminant")
    univariate_s = PolynomialRing(field, "s")
    s_eliminant = univariate_s(max(eliminants, key=lambda candidate: candidate.degree(s)))

    delta = ZZ(0)
    singularities = []
    for s_field, s_embedding, s_value in split_roots(s_eliminant, "bs"):
        extension_ring = PolynomialRing(s_field, names=("r", "s"))
        rr, ss = extension_ring.gens()
        f_extension = map_polynomial(polynomial, s_field, s_embedding)
        univariate_r = PolynomialRing(s_field, "r")
        rv = univariate_r.gen()
        common = gcd(
            gcd(
                univariate_r(f_extension(rv, s_value)),
                univariate_r(f_extension.derivative(rr)(rv, s_value)),
            ),
            univariate_r(f_extension.derivative(ss)(rv, s_value)),
        )
        for r_field, r_embedding, r_value in split_roots(common, "br"):
            if r_embedding is None:
                final_field = s_field
                final_r = r_value
                final_s = s_value
                final_polynomial = f_extension
            else:
                final_field = r_field
                final_r = r_value
                final_s = r_embedding(s_value)
                final_polynomial = map_polynomial(f_extension, r_field, r_embedding)
            local_ring = PolynomialRing(final_field, names=("x", "y"))
            x, y = local_ring.gens()
            local = local_ring(final_polynomial(final_r + x, final_s + y))
            point_delta = resolve_local(local)
            delta += point_delta
            singularities.append(
                {
                    "chart": "affine",
                    "residue_degree": int(final_field.degree()),
                    "multiplicity": int(local_multiplicity(local)),
                    "delta": int(point_delta),
                }
            )

    infinity_univariate = PolynomialRing(field, "a")
    a = infinity_univariate.gen()
    infinity_polynomials = [
        infinity_univariate(candidate(a, 1, 0))
        for candidate in (
            homogeneous,
            homogeneous.derivative(R),
            homogeneous.derivative(S),
            homogeneous.derivative(Z),
        )
    ]
    infinity_common = infinity_polynomials[0]
    for candidate in infinity_polynomials[1:]:
        infinity_common = gcd(infinity_common, candidate)
    for point_field, embedding, a_value in split_roots(infinity_common, "bi"):
        homogeneous_extension = map_polynomial(homogeneous, point_field, embedding)
        local_ring = PolynomialRing(point_field, names=("x", "y"))
        x, y = local_ring.gens()
        local = local_ring(homogeneous_extension(a_value + x, 1, y))
        point_delta = resolve_local(local)
        delta += point_delta
        singularities.append(
            {
                "chart": "infinity_S_1",
                "residue_degree": int(point_field.degree()),
                "multiplicity": int(local_multiplicity(local)),
                "delta": int(point_delta),
            }
        )

    derivatives = tuple(
        homogeneous.derivative(variable) for variable in (R, S, Z)
    )
    if homogeneous(1, 0, 0) == 0 and all(
        derivative(1, 0, 0) == 0 for derivative in derivatives
    ):
        local_ring = PolynomialRing(field, names=("x", "y"))
        x, y = local_ring.gens()
        local = local_ring(homogeneous(1, x, y))
        point_delta = resolve_local(local)
        delta += point_delta
        singularities.append(
            {
                "chart": "infinity_R_1",
                "residue_degree": 1,
                "multiplicity": int(local_multiplicity(local)),
                "delta": int(point_delta),
            }
        )

    arithmetic_genus = (degree - 1) * (degree - 2) // 2
    genus = arithmetic_genus - delta
    return {
        "degree": int(degree),
        "arithmetic_genus": int(arithmetic_genus),
        "delta": int(delta),
        "geometric_genus": int(genus),
        "geometric_singularities": len(singularities),
        "singularities": singularities,
    }


def rational_branch_count(polynomial):
    ring = polynomial.parent()
    field = ring.base_ring()
    x, y = ring.gens()
    order = local_multiplicity(polynomial)
    if order < 2:
        return 1
    tangent = initial_form(polynomial, order)
    univariate = PolynomialRing(field, "u")
    u = univariate.gen()
    count = 0
    finite = univariate(tangent(1, u))
    for irreducible, unused_tangent_multiplicity in finite.factor():
        if irreducible.degree() != 1:
            continue
        root = irreducible.roots()[0][0]
        transformed = polynomial(x, x * (root + y)) // x**order
        count += rational_branch_count(transformed)
    vertical = univariate(tangent(u, 1))
    if vertical.valuation() > 0:
        transformed = polynomial(x * y, y) // y**order
        count += rational_branch_count(transformed)
    return count


def normalized_point_count(polynomial, extension_degree):
    field = GF(PRIME) if extension_degree == 1 else GF(PRIME**extension_degree, "c")
    elements = tuple(field)
    homogeneous_ring = PolynomialRing(field, names=("R", "S", "Z"))
    R, S, Z = homogeneous_ring.gens()
    homogeneous = homogeneous_ring(polynomial.change_ring(field)).homogenize(Z)
    derivatives = tuple(homogeneous.derivative(variable) for variable in (R, S, Z))
    points = [(first, second, field.one()) for first in elements for second in elements]
    points += [(first, field.one(), field.zero()) for first in elements]
    points += [(field.one(), field.zero(), field.zero())]
    curve_points = [point for point in points if homogeneous(*point) == 0]
    singular = [
        point
        for point in curve_points
        if all(derivative(*point) == 0 for derivative in derivatives)
    ]
    branch_total = 0
    branch_counts = []
    for point in singular:
        pivot = next(index for index, value in enumerate(point) if value)
        free = [index for index in range(3) if index != pivot]
        local_ring = PolynomialRing(field, names=("x", "y"))
        x, y = local_ring.gens()
        values = [None, None, None]
        values[pivot] = field.one()
        values[free[0]] = point[free[0]] / point[pivot] + x
        values[free[1]] = point[free[1]] / point[pivot] + y
        local = local_ring(homogeneous(*values))
        branches = rational_branch_count(local)
        branch_total += branches
        branch_counts.append(branches)
    normalized = len(curve_points) - len(singular) + branch_total
    return {
        "extension_degree": extension_degree,
        "plane_points": len(curve_points),
        "rational_singular_points": len(singular),
        "rational_branches_over_singularities": branch_total,
        "branch_counts": branch_counts,
        "normalized_points": normalized,
    }


def hyperelliptic_point_count(coefficients, extension_degree):
    field = GF(PRIME) if extension_degree == 1 else GF(PRIME**extension_degree, "d")
    ring = PolynomialRing(field, "x")
    polynomial = ring(coefficients)
    count = 0
    for value in field:
        rhs = polynomial(value)
        count += 1 if rhs == 0 else (2 if rhs.is_square() else 0)
    count += 2 if polynomial.leading_coefficient().is_square() else 0
    return count


def frobenius_data(counts):
    q = PRIME
    N1, N2 = counts
    a1 = q + 1 - N1
    a2 = (N2 - q**2 - 1 + a1**2) // 2
    return {
        "a1": int(a1),
        "a2": int(a2),
        "L_polynomial_coefficients_ascending": [1, int(-a1), int(a2), int(-q * a1), q**2],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h21", required=True, type=Path, help="Gruenewald level1_21.txt")
    parser.add_argument("--h92", required=True, type=Path, help="Elkies--Kumar igusa92.txt")
    parser.add_argument("--output", type=Path, help="optional generated JSON result")
    arguments = parser.parse_args()

    stage("inputs")
    hashes = {
        "h21": verify_input(arguments.h21, H21_SHA256, "H21"),
        "h92": verify_input(arguments.h92, H92_SHA256, "H92"),
    }
    h21, h21_content = extract_h21(arguments.h21)
    unused_h92_ring, h92 = extract_h92(arguments.h92)
    stage("parsed", h21_terms=len(h21.dict()), h21_content=h21_content)

    pullback = build_pullback(h21, h92)
    valuations = monomial_valuations(pullback)
    if valuations != (88, 28):
        raise AssertionError(f"unexpected base monomial valuations {valuations}")
    r, s = pullback.parent().gens()
    residual = pullback // (r**valuations[0] * s**valuations[1])
    stage("base_factor_removed", r=valuations[0], s=valuations[1], degree=residual.total_degree())

    factor_start = time.monotonic()
    factorization = residual.factor()
    factor_degrees = tuple(sorted(int(component.total_degree()) for component, unused in factorization))
    stage(
        "factored",
        seconds=f"{time.monotonic() - factor_start:.2f}",
        factors=len(factor_degrees),
        degrees=",".join(map(str, factor_degrees)),
    )
    if factor_degrees != EXPECTED_FACTOR_DEGREES:
        raise AssertionError(f"unexpected factor degrees {factor_degrees}")
    if any(multiplicity != 1 for component, multiplicity in factorization):
        raise AssertionError("residual factorization is not squarefree")

    labels = label_cm24_components(factorization)
    stage(
        "cm24_labeled",
        branches=",".join(f"{branch}:{labels[branch]['degree']}" for branch in sorted(labels)),
    )
    target = labels[6]["polynomial"]

    genus_start = time.monotonic()
    genus = normalization_genus(target)
    stage(
        "target_genus",
        seconds=f"{time.monotonic() - genus_start:.2f}",
        arithmetic=genus["arithmetic_genus"],
        delta=genus["delta"],
        genus=genus["geometric_genus"],
    )
    if genus["geometric_genus"] != 2:
        raise AssertionError(f"target normalization has genus {genus['geometric_genus']}, not 2")

    target_count_records = [normalized_point_count(target, degree) for degree in (1, 2)]
    target_counts = tuple(record["normalized_points"] for record in target_count_records)
    published_models = {
        "level_402": (36, 0, 369, 0, 1170, 0, 729),
        "level_474": (576, 0, -171, 0, 198, 0, -27),
    }
    published_counts = {
        label: tuple(
            hyperelliptic_point_count(coefficients, degree) for degree in (1, 2)
        )
        for label, coefficients in published_models.items()
    }
    stage(
        "point_counts",
        target=target_counts,
        level402=published_counts["level_402"],
        level474=published_counts["level_474"],
    )
    if target_counts != (22, 116):
        raise AssertionError(f"unexpected target normalized counts {target_counts}")
    if published_counts["level_474"] != target_counts:
        raise AssertionError("target component does not match level-474 point counts")
    if published_counts["level_402"] == target_counts:
        raise AssertionError("level-402 control unexpectedly has the target point counts")

    frobenius = frobenius_data(target_counts)
    if frobenius["L_polynomial_coefficients_ascending"] != [1, 10, 47, 110, 121]:
        raise AssertionError(f"unexpected Frobenius polynomial {frobenius}")

    branch_payload = {
        str(branch): {key: value for key, value in data.items() if key != "polynomial"}
        for branch, data in sorted(labels.items())
    }
    output = {
        "schema": "elkies-k3.h21-h92-level474-branch-mod11.v1",
        "status": "PASS_COMPUTATIONAL_BRANCH_IDENTIFICATION",
        "proof_boundary": (
            "Exact modulo-11 component identification and Frobenius match; "
            "no characteristic-zero birational map is asserted."
        ),
        "prime": PRIME,
        "inputs": {
            "h21": {"path": str(arguments.h21), "sha256": hashes["h21"], "url": H21_URL},
            "h92": {"path": str(arguments.h92), "sha256": hashes["h92"], "source_url": EK_SOURCE_URL},
        },
        "h21": {"primitive_terms": len(h21.dict()), "content": str(h21_content)},
        "pullback": {
            "degree": int(pullback.total_degree()),
            "terms": len(pullback.dict()),
            "base_monomial": {"r": valuations[0], "s": valuations[1]},
            "residual_degree": int(residual.total_degree()),
            "residual_terms": len(residual.dict()),
            "factor_degrees": list(factor_degrees),
        },
        "cm24": {
            "minimal_polynomial": "w^2-w+1",
            "point": ["w", "-w"],
            "branches": branch_payload,
            "target_branch": 6,
        },
        "target_component": {
            "degree": labels[6]["degree"],
            "slope": BRANCH_SLOPES_QQ[5],
            "normalization": genus,
            "point_counts": target_count_records,
            "counts_F11_F121": list(target_counts),
            "frobenius": frobenius,
        },
        "published_controls": {
            "level_402": {
                "equation": "y^2=729*x^6+1170*x^4+369*x^2+36",
                "counts_F11_F121": list(published_counts["level_402"]),
            },
            "level_474": {
                "equation": "y^2=-27*x^6+198*x^4-171*x^2+576",
                "counts_F11_F121": list(published_counts["level_474"]),
            },
        },
        "conclusion": {
            "component": "branch_6",
            "height_form_twice": [[21, 6], [6, 92]],
            "determinant": 1896,
            "level": 474,
            "identification": "level_474_candidate_pinned_modulo_11",
        },
    }

    if arguments.output is not None:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        encoded = json.dumps(output, indent=2, sort_keys=True, default=int) + "\n"
        arguments.output.write_text(encoded)
        stage("artifact", path=arguments.output, sha256=hashlib.sha256(encoded.encode()).hexdigest())
    stage("complete", status=output["status"])


if __name__ == "__main__":
    main()
