#!/usr/bin/env python3
"""Exact characteristic-zero cyclic split of the cubic rank-two fiber.

At the first integral exact-rank-two point, this checker decomposes the
18-dimensional relative logarithmic critical algebra into its
14-dimensional interior algebra and two length-two endpoint algebras.  It
then eliminates the Laurent period function P=Q/u^3 from each summand.

The eliminant degrees equal the summand lengths and are pairwise coprime.
Consequently P is a primitive element on every summand; in particular,
1,P,...,P^13 is a cyclic basis of the interior algebra.  This is an exact
characteristic-zero seed for a discrete Gauss--Manin/Picard--Fuchs
calculation, not the connection or telescoping certificate itself.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
from math import factorial, gcd
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from verify_two_pair_sic_bidegree33_rank_two_holonomic_probe import (  # noqa: E402
    POINTS,
    matrix_product,
    pure_moment,
)
from verify_two_pair_sic_bidegree33_rank_two_relative_jacobian import (  # noqa: E402
    q_expression,
    run_singular,
)


OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_interior_cyclic_split.json"
)
EXPECTED_RELATIVE_BASIS = (
    "1",
    "u",
    "u2",
    "u3",
    "u4",
    "t",
    "ut",
    "u2t",
    "u3t",
    "t2",
    "ut2",
    "u2t2",
    "t3",
    "ut3",
    "u2t3",
    "t4",
    "ut4",
    "t5",
)
EXPECTED_INTERIOR_BASIS = (
    "1",
    "u",
    "u2",
    "u3",
    "t",
    "ut",
    "u2t",
    "u3t",
    "t2",
    "ut2",
    "u2t2",
    "t3",
    "ut3",
    "t4",
)
EXPECTED_REMOVED_MONOMIALS = ("u4", "u2t3", "ut4", "t5")


def run_singular_long(code: str, timeout: int = 180) -> str:
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required")
    completed = subprocess.run(
        [singular, "-q"],
        input=code,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=True,
    )
    if "?" in completed.stdout or completed.stderr:
        raise RuntimeError(completed.stdout + completed.stderr)
    return completed.stdout


def singular_calculation() -> str:
    q = q_expression(0)
    return run_singular(
        f"""
LIB "elim.lib";
ring base=0,(u,t),dp;
proc contained(ideal source, ideal target)
{{
  ideal remainder=reduce(source,target);
  int index;
  for(index=1;index<=size(remainder);index++)
  {{
    if(remainder[index]!=0){{return(0);}}
  }}
  return(1);
}}
poly Q={q};
poly A=u*diff(Q,u)-3*Q;
poly C=t*(1-t)*diff(Q,t);

ideal relative_base=std(sat(ideal(A,C),ideal(u)));
ideal interior_base=std(
  sat(
    sat(ideal(A,diff(Q,t)),ideal(u)),
    ideal(t*(1-t))
  )
);
ideal endpoint0_base=std(sat(ideal(A,t),ideal(u)));
ideal endpoint1_base=std(sat(ideal(A,t-1),ideal(u)));
ideal toric_base=std(
  sat(ideal(A,t*diff(Q,t)),ideal(u*t))
);
if(vdim(relative_base)!=18){{ERROR("relative length");}}
if(vdim(interior_base)!=14){{ERROR("interior length");}}
if(vdim(endpoint0_base)!=2){{ERROR("endpoint zero length");}}
if(vdim(endpoint1_base)!=2){{ERROR("endpoint one length");}}
if(vdim(toric_base)!=14){{ERROR("toric logarithmic length");}}
if(!contained(toric_base,interior_base))
  {{ERROR("toric subset interior");}}
if(!contained(interior_base,toric_base))
  {{ERROR("interior subset toric");}}

ideal localized=std(sat(relative_base,ideal(t*(1-t))));
if(!contained(localized,interior_base)){{ERROR("localized subset");}}
if(!contained(interior_base,localized)){{ERROR("interior subset");}}

ideal decomposition=std(
  intersect(intersect(interior_base,endpoint0_base),endpoint1_base)
);
if(!contained(relative_base,decomposition)){{ERROR("relative decomposition");}}
if(!contained(decomposition,relative_base)){{ERROR("decomposition relative");}}
if(vdim(std(interior_base+endpoint0_base))!=0)
  {{ERROR("interior endpoint zero overlap");}}
if(vdim(std(interior_base+endpoint1_base))!=0)
  {{ERROR("interior endpoint one overlap");}}
if(vdim(std(endpoint0_base+endpoint1_base))!=0)
  {{ERROR("endpoint overlap");}}

print("PASS exact relative/interior/endpoint ideal decomposition");
print("PASS pairwise comaximal summands of lengths 14+2+2");
print("PASS toric logarithmic Jacobian equals the interior algebra");
print("RELATIVE_BASIS");
print(kbase(relative_base));
print("INTERIOR_BASIS");
print(kbase(interior_base));

ring r=0,(u,t,z),dp;
ideal relative=imap(base,relative_base);
ideal interior=imap(base,interior_base);
ideal endpoint0=imap(base,endpoint0_base);
ideal endpoint1=imap(base,endpoint1_base);
poly Q=imap(base,Q);
ideal relative_elimination=std(
  eliminate(relative+ideal(z*u3-Q),u*t)
);
ideal interior_elimination=std(
  eliminate(interior+ideal(z*u3-Q),u*t)
);
ideal endpoint0_elimination=std(
  eliminate(endpoint0+ideal(z*u3-Q),u*t)
);
ideal endpoint1_elimination=std(
  eliminate(endpoint1+ideal(z*u3-Q),u*t)
);
if(size(relative_elimination)!=1){{ERROR("relative eliminant count");}}
if(size(interior_elimination)!=1){{ERROR("interior eliminant count");}}
if(size(endpoint0_elimination)!=1){{ERROR("endpoint zero eliminant count");}}
if(size(endpoint1_elimination)!=1){{ERROR("endpoint one eliminant count");}}
poly prel=relative_elimination[1];
poly pint=interior_elimination[1];
poly p0=endpoint0_elimination[1];
poly p1=endpoint1_elimination[1];
if(deg(prel)!=18){{ERROR("relative eliminant degree");}}
if(deg(pint)!=14){{ERROR("interior eliminant degree");}}
if(deg(p0)!=2){{ERROR("endpoint zero eliminant degree");}}
if(deg(p1)!=2){{ERROR("endpoint one eliminant degree");}}
ideal product_eliminant=std(ideal(pint*p0*p1));
if(reduce(prel,product_eliminant)!=0){{ERROR("relative product divisibility");}}
if(reduce(pint*p0*p1,relative_elimination)!=0)
  {{ERROR("product relative divisibility");}}
if(gcd(pint,p0)!=1){{ERROR("interior endpoint zero gcd");}}
if(gcd(pint,p1)!=1){{ERROR("interior endpoint one gcd");}}
if(gcd(p0,p1)!=1){{ERROR("endpoint gcd");}}
if(gcd(pint,diff(pint,z))!=1){{ERROR("interior squarefree");}}
if(gcd(p0,diff(p0,z))!=1){{ERROR("endpoint zero squarefree");}}
if(gcd(p1,diff(p1,z))!=1){{ERROR("endpoint one squarefree");}}

print("PASS period function is primitive on every summand");
print("ENDPOINT0_POLYNOMIAL");
print(p0);
print("ENDPOINT1_POLYNOMIAL");
print(p1);
print("INTERIOR_POLYNOMIAL");
print(pint);
print("RELATIVE_POLYNOMIAL");
print(prel);
print("END");
"""
    )


def marked_block(output: str, marker: str, next_marker: str) -> list[str]:
    lines = output.splitlines()
    start = lines.index(marker) + 1
    end = lines.index(next_marker, start)
    return [line.strip() for line in lines[start:end] if line.strip()]


def parse_basis(lines: list[str]) -> tuple[str, ...]:
    return tuple(
        line.rstrip(",").replace("^", "").replace("*", "").replace(" ", "")
        for line in lines
    )


def parse_univariate_singular_polynomial(raw: str) -> list[int]:
    """Parse Singular's compact integer univariate output, low to high."""

    coefficients: dict[int, int] = {}
    for token in re.findall(r"[+-]?[^+-]+", raw.strip()):
        if not token:
            continue
        if "z" in token:
            scalar, exponent = token.split("z", maxsplit=1)
            if scalar in ("", "+"):
                coefficient = 1
            elif scalar == "-":
                coefficient = -1
            else:
                coefficient = int(scalar)
            power = int(exponent) if exponent else 1
        else:
            coefficient = int(token)
            power = 0
        coefficients[power] = coefficient
    degree = max(coefficients)
    result = [coefficients.get(power, 0) for power in range(degree + 1)]
    content = 0
    for coefficient in result:
        content = gcd(content, abs(coefficient))
    assert content
    result = [coefficient // content for coefficient in result]
    if result[-1] < 0:
        result = [-coefficient for coefficient in result]
    return result


def scalar_multiple(
    left: list[int],
    right: list[int],
) -> Fraction | None:
    if len(left) != len(right):
        return None
    ratio = None
    for left_value, right_value in zip(left, right, strict=True):
        if right_value == 0:
            if left_value != 0:
                return None
            continue
        candidate = Fraction(left_value, right_value)
        if ratio is None:
            ratio = candidate
        elif ratio != candidate:
            return None
    return ratio


def polynomial_multiply(left: list[int], right: list[int]) -> list[int]:
    answer = [0] * (len(left) + len(right) - 1)
    for left_degree, left_value in enumerate(left):
        for right_degree, right_value in enumerate(right):
            answer[left_degree + right_degree] += left_value * right_value
    return answer


def sympy_polynomial(coefficients: list[int], variable: sp.Symbol) -> sp.Poly:
    return sp.Poly(
        sum(
            sp.Integer(coefficient) * variable**degree
            for degree, coefficient in enumerate(coefficients)
        ),
        variable,
        domain=sp.QQ,
    ).monic()


def exact_normalized_moments(maximum: int) -> list[Fraction]:
    raw_u, raw_w = POINTS[0]
    u = [[Fraction(value) for value in row] for row in raw_u]
    w = [[Fraction(value) for value in row] for row in raw_w]
    matrix = matrix_product(u, w)
    moments = [Fraction(1)]
    for order in range(1, maximum + 1):
        moments.append(
            pure_moment(matrix, order)
            / Fraction(factorial(3 * order + 1))
        )
    return moments


def rational_fingerprint(value: Fraction) -> dict[str, object]:
    encoded = f"{value.numerator}/{value.denominator}".encode()
    return {
        "nonzero": value != 0,
        "numerator_bits": abs(value.numerator).bit_length(),
        "denominator_bits": value.denominator.bit_length(),
        "sha256": hashlib.sha256(encoded).hexdigest(),
    }


def idempotent_period_audit(
    interior: list[int],
    endpoint0: list[int],
    endpoint1: list[int],
) -> dict[str, object]:
    """Check whether the raw period functional kills Jacobian summands."""

    z = sp.symbols("z")
    factors = {
        "interior": sympy_polynomial(interior, z),
        "endpoint_t0": sympy_polynomial(endpoint0, z),
        "endpoint_t1": sympy_polynomial(endpoint1, z),
    }
    relative = factors["interior"] * factors["endpoint_t0"] * factors[
        "endpoint_t1"
    ]
    idempotents: dict[str, sp.Poly] = {}
    for name, factor in factors.items():
        complementary = sp.quo(relative, factor)
        inverse = sp.invert(complementary, factor)
        idempotent = sp.rem(complementary * inverse, relative)
        assert sp.rem(idempotent**2 - idempotent, relative).is_zero
        idempotents[name] = idempotent
    assert sp.rem(
        sum(idempotents.values(), sp.Poly(0, z)) - sp.Poly(1, z),
        relative,
    ).is_zero

    moments = exact_normalized_moments(18)

    def period_value(polynomial: sp.Poly, shift: int) -> Fraction:
        return sum(
            Fraction(polynomial.nth(degree))
            * moments[degree + shift]
            for degree in range(polynomial.degree() + 1)
        )

    values = {
        name: [
            period_value(idempotent, 0),
            period_value(idempotent, 1),
        ]
        for name, idempotent in idempotents.items()
    }
    assert all(
        value != 0
        for name in ("endpoint_t0", "endpoint_t1")
        for value in values[name]
    )
    assert sum(value[0] for value in values.values()) == moments[0]
    assert sum(value[1] for value in values.values()) == moments[1]
    return {
        "moments_computed": 19,
        "tests_per_summand": ["period(e)", "period(P*e)"],
        "values": {
            name: [rational_fingerprint(value) for value in pair]
            for name, pair in values.items()
        },
        "endpoint_values_are_nonzero": True,
        "conclusion": (
            "the raw period functional does not descend through the "
            "ordinary Jacobian quotient; an m-dependent divergence "
            "reduction is required before the cyclic split can yield a "
            "recurrence"
        ),
    }


def gradient_lift_seed(interior: list[int]) -> dict[str, object]:
    """Lift the interior eliminant through the raw gradient generators."""

    q = q_expression(0)
    terms = "+".join(
        f"({coefficient})*Q^{power}*u^{42 - 3 * power}"
        for power, coefficient in enumerate(interior)
        if coefficient
    )
    output = run_singular_long(
        f"""
LIB "elim.lib";
ring r=0,(u,t),dp;
poly Q={q};
poly A=u*diff(Q,u)-3*Q;
poly B=diff(Q,t);
list saturation=sat_with_exp(ideal(A,B),ideal(u));
if(saturation[2]!=5){{ERROR("interior u-saturation exponent");}}
ideal interior=std(saturation[1]);
poly F={terms};
poly target=u5*F;
matrix lift_matrix=lift(ideal(A,B),ideal(target));
poly X=lift_matrix[1,1];
poly Y=lift_matrix[2,1];
if(target-X*A-Y*B!=0){{ERROR("gradient lift identity");}}
poly divergence=u*diff(X,u)+diff(Y,t);
poly boundary0=subst(Y,t,0);
poly boundary1=subst(Y,t,1);
if(boundary0==0){{ERROR("unexpected zero t=0 boundary");}}
if(boundary1==0){{ERROR("unexpected zero t=1 boundary");}}
ideal endpoint0=std(sat(ideal(A,t),ideal(u)));
ideal endpoint1=std(sat(ideal(A,t-1),ideal(u)));
poly reduced_divergence=reduce(divergence,interior);
poly reduced_boundary0=reduce(boundary0,endpoint0);
poly reduced_boundary1=reduce(boundary1,endpoint1);
if(size(reduced_divergence)!=14)
  {{ERROR("interior divergence coordinate count");}}
if(size(reduced_boundary0)!=2)
  {{ERROR("endpoint zero coordinate count");}}
if(size(reduced_boundary1)!=2)
  {{ERROR("endpoint one coordinate count");}}
print("PASS exact interior gradient lift");
print("PROFILE");
print(deg(F));
print(size(F));
print(deg(X));
print(size(X));
print(deg(Y));
print(size(Y));
print(deg(divergence));
print(size(divergence));
print(size(boundary0));
print(size(boundary1));
print(deg(reduced_divergence));
print(size(reduced_divergence));
print(deg(reduced_boundary0));
print(size(reduced_boundary0));
print(deg(reduced_boundary1));
print(size(reduced_boundary1));
print("END_PROFILE");
"""
    )
    assert "PASS exact interior gradient lift" in output
    profile = [
        int(value)
        for value in marked_block(output, "PROFILE", "END_PROFILE")
    ]
    assert len(profile) == 16
    (
        eliminant_degree,
        eliminant_terms,
        x_degree,
        x_terms,
        y_degree,
        y_terms,
        divergence_degree,
        divergence_terms,
        boundary0_terms,
        boundary1_terms,
        reduced_divergence_degree,
        reduced_divergence_terms,
        reduced_boundary0_degree,
        reduced_boundary0_terms,
        reduced_boundary1_degree,
        reduced_boundary1_terms,
    ) = profile
    return {
        "saturation_exponent": 5,
        "identity": (
            "u^47*p_int(P)=X*(u*Q_u-3Q)+Y*Q_t, "
            "equivalently u^44*p_int(P)=X*D_u(P)+Y*partial_t(P)"
        ),
        "cleared_eliminant": {
            "total_degree": eliminant_degree,
            "terms": eliminant_terms,
        },
        "X": {"total_degree": x_degree, "terms": x_terms},
        "Y": {"total_degree": y_degree, "terms": y_terms},
        "divergence_DuX_plus_dY": {
            "total_degree": divergence_degree,
            "terms": divergence_terms,
        },
        "endpoint_boundary_terms": {
            "t0": boundary0_terms,
            "t1": boundary1_terms,
        },
        "ordinary_jacobian_normal_forms": {
            "interior_divergence": {
                "total_degree": reduced_divergence_degree,
                "coordinates": reduced_divergence_terms,
            },
            "endpoint_t0": {
                "total_degree": reduced_boundary0_degree,
                "coordinates": reduced_boundary0_terms,
            },
            "endpoint_t1": {
                "total_degree": reduced_boundary1_degree,
                "coordinates": reduced_boundary1_terms,
            },
        },
        "integration_by_parts_identity": (
            "(m+1)*Integral[u^44*p_int(P)*P^m] = "
            "-Integral[(D_u X+partial_t Y)*P^(m+1)] "
            "+Boundary[Y*P^(m+1)]_(t=0)^(t=1)"
        ),
        "next_reduction": (
            "recursively lift the gradient parts discarded by the "
            "14+2+2 ordinary normal forms to construct the full "
            "m-dependent connection"
        ),
    }


def main() -> None:
    output = singular_calculation()
    for marker in (
        "PASS exact relative/interior/endpoint ideal decomposition",
        "PASS pairwise comaximal summands of lengths 14+2+2",
        "PASS toric logarithmic Jacobian equals the interior algebra",
        "PASS period function is primitive on every summand",
    ):
        assert marker in output

    relative_basis = parse_basis(
        marked_block(output, "RELATIVE_BASIS", "INTERIOR_BASIS")
    )
    interior_basis = parse_basis(
        marked_block(
            output,
            "INTERIOR_BASIS",
            "PASS period function is primitive on every summand",
        )
    )
    assert set(relative_basis) == set(EXPECTED_RELATIVE_BASIS)
    assert set(interior_basis) == set(EXPECTED_INTERIOR_BASIS)
    assert (
        set(relative_basis) - set(interior_basis)
        == set(EXPECTED_REMOVED_MONOMIALS)
    )

    endpoint0 = parse_univariate_singular_polynomial(
        marked_block(
            output,
            "ENDPOINT0_POLYNOMIAL",
            "ENDPOINT1_POLYNOMIAL",
        )[0]
    )
    endpoint1 = parse_univariate_singular_polynomial(
        marked_block(
            output,
            "ENDPOINT1_POLYNOMIAL",
            "INTERIOR_POLYNOMIAL",
        )[0]
    )
    interior = parse_univariate_singular_polynomial(
        marked_block(
            output,
            "INTERIOR_POLYNOMIAL",
            "RELATIVE_POLYNOMIAL",
        )[0]
    )
    relative = parse_univariate_singular_polynomial(
        marked_block(output, "RELATIVE_POLYNOMIAL", "END")[0]
    )
    assert len(endpoint0) == 3
    assert len(endpoint1) == 3
    assert len(interior) == 15
    assert len(relative) == 19
    product = polynomial_multiply(
        polynomial_multiply(interior, endpoint0),
        endpoint1,
    )
    assert scalar_multiple(relative, product) is not None
    period_audit = idempotent_period_audit(
        interior,
        endpoint0,
        endpoint1,
    )
    lift_seed = gradient_lift_seed(interior)

    artifact = {
        "format": (
            "two-pair-sic-bidegree33-rank-two-"
            "interior-cyclic-split-v1"
        ),
        "status": (
            "exact characteristic-zero critical-algebra calculation at "
            "one integral rank-two point; not a universal connection"
        ),
        "point": {
            "U": POINTS[0][0],
            "W": POINTS[0][1],
            "coefficient_rank": 2,
        },
        "relative_critical_algebra": {
            "length": 18,
            "standard_monomial_basis": list(EXPECTED_RELATIVE_BASIS),
            "chinese_remainder_lengths": {
                "interior": 14,
                "endpoint_t0": 2,
                "endpoint_t1": 2,
            },
            "interior_standard_monomial_basis": list(
                EXPECTED_INTERIOR_BASIS
            ),
            "relative_basis_monomials_removed_by_interior_localization": (
                list(EXPECTED_REMOVED_MONOMIALS)
            ),
        },
        "toric_logarithmic_critical_algebra": {
            "ideal": (
                "sat((u*Q_u-3Q, t*Q_t), u*t)"
            ),
            "length": 14,
            "equals_interior_algebra": True,
            "standard_monomial_basis": list(EXPECTED_INTERIOR_BASIS),
            "picard_fuchs_target": (
                "the exact toric logarithmic critical rank equals the "
                "sampled minimal Ore-gcd order 14"
            ),
        },
        "period_function_eliminants": {
            "variable": "z=P=Q/u^3",
            "coefficient_order": "low_to_high",
            "endpoint_t0": endpoint0,
            "endpoint_t1": endpoint1,
            "interior": interior,
            "relative": relative,
            "relative_is_scalar_multiple_of_product": True,
            "pairwise_coprime": True,
            "all_three_factors_squarefree": True,
        },
        "cyclic_consequence": {
            "relative_basis": "1,P,...,P^17",
            "interior_basis": "1,P,...,P^13",
            "endpoint_bases": ["1,P", "1,P"],
            "reason": (
                "each eliminant degree equals the corresponding algebra "
                "length, so Q[P] has full dimension on every summand"
            ),
        },
        "jacobian_idempotent_period_audit": period_audit,
        "gradient_lift_seed": lift_seed,
        "next_gate": (
            "construct the m-dependent reduction/connection in the "
            "interior cyclic basis and recover the sampled order-14 "
            "shift operator with divergence certificates"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print("PASS exact 18=14+2+2 Chinese-remainder decomposition")
    print("PASS pairwise-coprime degree 14,2,2 period eliminants")
    print("PASS exact interior cyclic basis 1,P,...,P^13")
    print("PASS exact toric logarithmic critical rank is 14")
    print("PASS endpoint idempotent periods are nonzero")
    print("PASS m-dependent divergence reduction remains necessary")
    print("PASS exact gradient lift of the interior eliminant")
    print("PASS nonzero endpoint boundary terms retained")
    print("PASS characteristic-zero connection seed, not full certificate")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
