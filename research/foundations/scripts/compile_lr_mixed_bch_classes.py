#!/usr/bin/env python3
"""Compile the balanced linear-in-X BCH classes for the F_2 LR problem.

Let D_B and D_C be the lifts of the constant target fields partial_B and
partial_C.  They commute as Lie fields.  For the normalized weight-zero
deformation direction X=N(x,0,-3z), N=v^6*S^4, define

    W_k = (ad(D_B) ad(D_C))^k X
        = {{...{{X,D_B},D_C}...,D_B},D_C}.

The balanced weight-zero part of the term linear in X in
BCH_map(-X,-(D_B+D_C)) is a nonzero Bernoulli scalar times W_k.

On the gamma-boundary face the mixed operator is an exact two-component
polynomial transfer.  Its third-summand constant term is triangular:

    c_(k+1) = -73440 (k+3) (2k+7) c_k,  k>=1.

This proves nonvanishing in every odd BCH order 2k+1.  The class survives the
saturated linear target quotient.  It is not universal over the lower-jet
scheme by itself: with target amplitudes s*partial_B+t*partial_C it is
multiplied by s^k*t^k and vanishes on s*t=0.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp

from compile_lr_rooted_tree_classes import (
    SemiInvariantField,
    constant_source_lift,
    gamma,
    map_bracket,
    u,
    v,
    S,
)


w = sp.Symbol("w")
d_symbol = sp.Symbol("d", integer=True, nonnegative=True)
SEPARATOR = sp.Rational(1, 6)


def normalized_deformation_direction() -> SemiInvariantField:
    """The leading field N(x,0,-3z), normalized by removing 435/7."""
    N = sp.expand(v**6 * S**4)
    return SemiInvariantField(
        0,
        (
            N,
            sp.expand(v * N),
            sp.expand((-sp.Rational(8, 7) * v - S) * N),
        ),
    )


def mixed_pair_operator(field: SemiInvariantField) -> SemiInvariantField:
    """Apply {{field,D_B},D_C}; this equals ad(D_B)ad(D_C)(field)."""
    D_B = constant_source_lift("B")
    D_C = constant_source_lift("C")
    return map_bracket(map_bracket(field, D_B), D_C)


def boundary_face_polynomial(expression: sp.Expr, base: int) -> sp.Expr:
    """Extract sum c_(a,b) w^b over monomials u^a gamma^b with a-b=base."""
    terms = sp.Poly(sp.expand(expression), u, gamma).terms()
    assert max(u_power - gamma_power for (u_power, gamma_power), _ in terms) == base
    return sp.expand(
        sum(
            coefficient * w**gamma_power
            for (u_power, gamma_power), coefficient in terms
            if u_power - gamma_power == base
        )
    )


def boundary_profile(
    field: SemiInvariantField, d_value: int
) -> tuple[sp.Expr, sp.Expr]:
    """Return F,H from the face (u^(d+1)F(w),u^(d+2)F(w),u^d H(w))."""
    first = boundary_face_polynomial(field.logarithmic[0], d_value + 1)
    second = boundary_face_polynomial(field.logarithmic[1], d_value + 2)
    third = boundary_face_polynomial(field.logarithmic[2], d_value)
    assert sp.expand(first - second) == 0
    return sp.factor(first), sp.factor(third)


def symbolic_boundary_operator() -> tuple[sp.Expr, sp.Expr]:
    """Derive the exact leading boundary-face operator for symbolic d,F,H."""
    F_function = sp.Function("F")
    H_function = sp.Function("H")
    F_value = F_function(w)
    H_value = H_function(w)
    generic = SemiInvariantField(
        0,
        (
            u ** (d_symbol + 1) * F_function(u * gamma),
            u ** (d_symbol + 2) * F_function(u * gamma),
            u**d_symbol * H_function(u * gamma),
        ),
    )
    image = mixed_pair_operator(generic)
    offsets = (5, 6, 4)
    faces = []
    for component, offset in zip(image.logarithmic, offsets, strict=True):
        scaled = component.subs(gamma, w / u) / u ** (d_symbol + offset)
        face = sp.simplify(sp.limit(scaled, u, sp.oo))
        assert not face.has(sp.Derivative)
        faces.append(face)
    assert sp.expand(faces[0] - faces[1]) == 0

    F_plain, H_plain = sp.symbols("F H")
    replacements = {F_value: F_plain, H_value: H_plain}
    F_new = sp.collect(sp.expand(faces[0].xreplace(replacements)), (F_plain, H_plain))
    H_new = sp.collect(sp.expand(faces[2].xreplace(replacements)), (F_plain, H_plain))
    assert not F_new.has(F_value, H_value)
    assert not H_new.has(F_value, H_value)

    F_at_zero = sp.factor(F_new.subs(w, 0))
    H_at_zero = sp.factor(H_new.subs(w, 0))
    assert sp.expand(
        F_at_zero
        + 540
        * (d_symbol + 1)
        * (17 * (d_symbol - 1) * F_plain + 317 * H_plain)
    ) == 0
    assert sp.expand(
        H_at_zero + 9180 * (d_symbol + 1) * (d_symbol + 3) * H_plain
    ) == 0
    return F_new, H_new


def bch_balanced_scalar(k: int) -> sp.Expr:
    """Coefficient of W_k in the balanced linear-in-X BCH component."""
    return sp.factor(
        -sp.bernoulli(2 * k)
        * sp.binomial(2 * k, k)
        / sp.factorial(2 * k)
    )


def field_is_zero(field: SemiInvariantField) -> bool:
    return all(entry == 0 for entry in field.logarithmic)


def compile_certificate(max_k: int) -> dict[str, object]:
    if max_k < 2:
        raise ValueError("max_k must be at least two")

    D_B = constant_source_lift("B")
    D_C = constant_source_lift("C")
    assert field_is_zero(map_bracket(D_B, D_C))

    F_new, H_new = symbolic_boundary_operator()
    F_plain, H_plain = sp.symbols("F H")

    field = normalized_deformation_direction()
    rows = []
    previous_leading = None
    previous_F = None
    previous_H = None
    for k in range(1, max_k + 1):
        field = mixed_pair_operator(field)
        d_value = 4 * k + 11
        F_face, H_face = boundary_profile(field, d_value)
        residue = sp.expand(field.logarithmic[2].subs(gamma, 0))
        residue_polynomial = sp.Poly(residue, u)
        assert residue_polynomial.degree() == d_value
        leading = sp.factor(residue_polynomial.LC())
        assert leading == sp.factor(H_face.subs(w, 0))

        if k == 1:
            assert leading == sp.Rational(14438891520, 2401)
        else:
            predicted_F = sp.factor(
                F_new.subs(
                    {
                        d_symbol: d_value - 4,
                        F_plain: previous_F,
                        H_plain: previous_H,
                    }
                )
            )
            predicted_H = sp.factor(
                H_new.subs(
                    {
                        d_symbol: d_value - 4,
                        F_plain: previous_F,
                        H_plain: previous_H,
                    }
                )
            )
            assert sp.expand(F_face - predicted_F) == 0
            assert sp.expand(H_face - predicted_H) == 0
            assert leading == -73440 * (k + 2) * (2 * k + 5) * previous_leading

        bch_scalar = bch_balanced_scalar(k)
        assert bch_scalar != 0
        rows.append(
            {
                "k": k,
                "bch_order": 2 * k + 1,
                "word": "(ad(D_B)*ad(D_C))^k X",
                "residue_degree_u": residue_polynomial.degree(),
                "residue_terms": len(residue_polynomial.terms()),
                "residue_leading_coefficient": sp.sstr(leading),
                "residue_at_u_1_over_6": sp.sstr(
                    sp.factor(residue.subs(u, SEPARATOR))
                ),
                "balanced_bch_scalar": sp.sstr(bch_scalar),
                "balanced_bch_leading_coefficient": sp.sstr(
                    sp.factor(bch_scalar * leading)
                ),
            }
        )
        previous_leading = leading
        previous_F = F_face
        previous_H = H_face

    H_constant = sp.factor(H_new.subs({w: 0, F_plain: 0, H_plain: 1}))
    return {
        "description": "Balanced linear-in-X BCH normal classes for the F_2 target lift",
        "conventions": {
            "map_prelie": "P ▷ Q = (DP)Q",
            "map_bracket": "{P,Q} = P▷Q-Q▷P",
            "bch": "BCH_map(-X,-D)",
            "D": "D_B+D_C",
            "target_directions": {
                "D_B": "ell_F(partial_B), weight 1",
                "D_C": "ell_F(partial_C), weight -1",
            },
            "X": "N(x,0,-3z), N=v^6*S^4, weight 0",
        },
        "commutation": {
            "D_B_D_C_bracket_zero": True,
            "balanced_component": "binomial(2k,k)*(ad(D_B)*ad(D_C))^k X",
            "reason": "[ad(D_B),ad(D_C)]=ad([D_B,D_C])=0",
        },
        "boundary_face": {
            "coordinate": "w=u*gamma",
            "input": "(u^(d+1)F(w),u^(d+2)F(w),u^d H(w))",
            "output": "(u^(d+5)F_new(w),u^(d+6)F_new(w),u^(d+4)H_new(w))",
            "F_new": sp.sstr(F_new),
            "H_new": sp.sstr(H_new),
            "constant_transfer": {
                "F_new_at_w_0": "-540*(d+1)*(17*(d-1)*F(0)+317*H(0))",
                "H_new_at_w_0": sp.sstr(H_constant) + "*H(0)",
            },
        },
        "all_order_theorem": {
            "seed": "c_1=14438891520/2401",
            "degree": "deg_u rho(W_k)=4*k+11",
            "recurrence": "c_(k+1)=-73440*(k+3)*(2*k+7)*c_k for k>=1",
            "bch_scalar": "-B_(2k)*binomial(2k,k)/(2k)!",
            "conclusion": (
                "The balanced coefficient linear in X and of bidegree (k,k) "
                "in (D_B,D_C) has nonzero third saturated associated-graded "
                "normal residue in every odd BCH order 2k+1."
            ),
        },
        "computed_regression": rows,
        "descent_audit": {
            "linear_target_module": (
                "PASS: the residue is taken in the saturated third normal "
                "summand R/(gamma), so every linear lifted-target correction "
                "maps to zero while the displayed residue remains nonzero."
            ),
            "lower_target_amplitudes": (
                "FAILS UNIVERSALITY: replacing D_B,D_C by s*D_B,t*D_C "
                "multiplies W_k by s^k*t^k. The sector vanishes on s*t=0, "
                "including the valid lower-jet choice D_1=0."
            ),
            "status": (
                "This proves actual BCH noncancellation in one multihomogeneous "
                "mixed sector and gives an all-order penalty when both "
                "opposite weights are active. It is not a universal obstruction "
                "over the full lower-jet solution scheme."
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-k", type=int, default=3)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("artifacts/generated-results/lr_mixed_bch_classes.json"),
    )
    parser.add_argument("--no-write", action="store_true")
    arguments = parser.parse_args()

    certificate = compile_certificate(arguments.max_k)
    payload = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()
    if not arguments.no_write:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(payload, encoding="utf-8")

    print("PASS: D_B and D_C commute and the balanced BCH sum collapses to W_k")
    print("PASS: exact boundary-face operator and all-order leading recurrence")
    print("PASS: nonzero mixed BCH normal class in every odd order 2*k+1")
    print("PASS: class survives the saturated linear target quotient")
    print("LIMIT: sector vanishes on the lower-jet locus s*t=0")
    print(f"certificate_sha256={digest}")


if __name__ == "__main__":
    main()
