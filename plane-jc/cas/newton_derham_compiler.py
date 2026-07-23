#!/usr/bin/env python3
"""A small compiler IR for Newton-chain weighted-Wronskian blocks.

The front end deliberately starts *after* Laurent polygons and coefficient
bands have been derived.  It records that missing provenance instead of
inventing equations from corner data.  A complete block is compiled to:

    superelliptic curve -> character -> compact de Rham data
        -> supported primitive -> tail obstruction certificate.

The normalized (72,108) first block is the golden fixture.  The repeated-tail
(96,144) row is represented as source-excluded after its complete-chain
audit, while the (75,125) F2 row remains incomplete because the published
table supplies corners, not Laurent coefficient bands.
"""

from dataclasses import dataclass
from functools import reduce
from operator import mul
from typing import Iterable

import sympy as sp

from superelliptic_derham import (
    SuperellipticDeRham,
    weighted_wronskian_compatibility,
)


Corner = tuple[sp.Rational, sp.Rational]


@dataclass(frozen=True)
class NewtonChainIR:
    """Source-level chain metadata, whether or not bands are yet available."""

    name: str
    corners: tuple[Corner, ...]
    multiplicities: tuple[int, int]
    enumeration_source: str
    status: str
    missing_frontend_data: tuple[str, ...] = ()
    source_reconciliation: str = ""

    @property
    def frontend_complete(self) -> bool:
        return not self.missing_frontend_data


@dataclass(frozen=True)
class WeightedWronskianIR:
    """A fully derived leading block ready for geometric compilation."""

    chain: NewtonChainIR
    t: sp.Symbol
    A: sp.Expr
    R: sp.Expr
    covering_exponent: int
    primitive_weight: int
    primitive_exponents: tuple[int, ...]
    full_primitive_bounds: tuple[int, int]
    normalization: tuple[str, ...] = ()
    scaling_weights: tuple[tuple[sp.Symbol, int], ...] = ()


@dataclass(frozen=True)
class TailObstructionCertificate:
    """Exact representative-level bridge from support equations to H^1."""

    primitive: sp.Expr
    solved_coefficients: tuple[tuple[sp.Symbol, sp.Expr], ...]
    compatibility: tuple[sp.Expr, ...]
    residual: sp.Expr
    residual_degrees: tuple[int, ...]
    tail_differentials: tuple[sp.Expr, ...]
    tail_infinity_orders: tuple[int, ...]
    tail_is_second_kind: bool
    de_rham_coordinates: tuple[sp.Expr, ...]
    low_operator_determinant: sp.Expr
    identity_holds: bool
    tail_basis_certified: bool


@dataclass(frozen=True)
class CompiledDeRhamBlock:
    """Canonical compiler output for one weighted-Wronskian block."""

    ir: WeightedWronskianIR
    genus: int
    character: int
    affine_dimension: int
    compact_dimension: int
    residue_functional: tuple[sp.Expr, ...]
    compact_basis: tuple[sp.Expr, ...]
    certificate: TailObstructionCertificate

    @property
    def local_system_fingerprint(self) -> tuple[object, ...]:
        polynomial = sp.Poly(self.ir.A, self.ir.t)
        support = tuple(i for i in range(polynomial.degree() + 1) if polynomial.nth(i))
        return (
            self.ir.covering_exponent,
            polynomial.degree(),
            self.character,
            sp.gcd(self.ir.covering_exponent, polynomial.degree()),
            self.compact_dimension,
            support,
        )

    @property
    def supported_problem_fingerprint(self) -> tuple[object, ...]:
        return (
            self.local_system_fingerprint,
            tuple(self.ir.primitive_exponents),
            tuple(sp.Poly(self.ir.R, self.ir.t).monoms()),
            tuple(weight for _, weight in self.ir.scaling_weights),
        )


@dataclass(frozen=True)
class LocalSystemReuseCertificate:
    """Exact cyclic-curve isomorphism underlying Gauss--Manin reuse."""

    reference_chain: str
    candidate_chain: str
    base_substitution: tuple[tuple[sp.Symbol, sp.Expr], ...]
    t_scale: sp.Expr
    y_scale: sp.Expr
    covering_exponent: int
    character: int
    curve_identity: sp.Expr


def certify_local_system_reuse(
    reference: WeightedWronskianIR,
    candidate: WeightedWronskianIR,
    *,
    base_substitution: tuple[tuple[sp.Symbol, sp.Expr], ...] = (),
    t_scale: sp.Expr = sp.S.One,
    y_scale: sp.Expr = sp.S.One,
) -> LocalSystemReuseCertificate:
    """Certify ``t_candidate=t_scale*t`` and ``y_candidate=y_scale*y``.

    The base substitution is applied to the candidate coefficients first.
    A successful result proves

        A_candidate(t_scale*t) = y_scale**a * A_reference(t),

    together with equality of the cyclic exponent and character.  It says
    nothing about the forcing term or supported primitive, which belong to
    the stronger plane-JC obstruction problem rather than its local system.
    """

    a = int(reference.covering_exponent)
    if int(candidate.covering_exponent) != a:
        raise ValueError("covering exponents differ")
    character = int(reference.primitive_weight) % a
    if int(candidate.primitive_weight) % a != character:
        raise ValueError("character eigenspaces differ")
    t_scale = sp.sympify(t_scale)
    y_scale = sp.sympify(y_scale)
    if t_scale.is_zero is True or y_scale.is_zero is True:
        raise ValueError("curve-coordinate scales must be nonzero")
    substitution = dict(base_substitution)
    if len(substitution) != len(base_substitution):
        raise ValueError("base substitution contains a duplicate source parameter")
    if candidate.t in substitution:
        raise ValueError("put the candidate t-coordinate in t_scale, not the base map")
    if any(
        value.has(reference.t, candidate.t)
        for value in map(sp.sympify, substitution.values())
    ):
        raise ValueError("a base substitution may not depend on the fiber coordinate")
    pulled_back = candidate.A.subs(substitution, simultaneous=True)
    pulled_back = pulled_back.subs(candidate.t, t_scale * reference.t)
    identity = sp.factor(
        sp.expand(pulled_back - y_scale**a * reference.A)
    )
    if identity != 0:
        raise ValueError(
            "the proposed base/coordinate change does not identify the curve families"
        )
    return LocalSystemReuseCertificate(
        reference_chain=reference.chain.name,
        candidate_chain=candidate.chain.name,
        base_substitution=base_substitution,
        t_scale=t_scale,
        y_scale=y_scale,
        covering_exponent=a,
        character=character,
        curve_identity=identity,
    )


def _substitute_in_order(expression: sp.Expr, substitutions: dict[sp.Symbol, sp.Expr]) -> sp.Expr:
    result = expression
    for variable, value in substitutions.items():
        result = result.subs(variable, value)
    return sp.expand(result)


def _low_operator_determinant(ir: WeightedWronskianIR) -> tuple[sp.Expr, bool]:
    """Determinant of the low-coefficient map on all possible primitives."""

    lower, upper = ir.full_primitive_bounds
    if lower != 0:
        raise ValueError("the current tail certificate expects a zero lower bound")
    t = ir.t
    a = ir.covering_exponent
    b = ir.primitive_weight
    A = ir.A
    exponents = tuple(range(upper + 1))
    columns = []
    for exponent in exponents:
        monomial = t**exponent
        image = sp.Poly(
            sp.expand(a * A * sp.diff(monomial, t) - b * sp.diff(A, t) * monomial),
            t,
        )
        columns.append([image.nth(degree) for degree in exponents])
    matrix = sp.Matrix(len(exponents), len(exponents), lambda i, j: columns[j][i])
    triangular = all(matrix[i, j] == 0 for i in range(len(exponents)) for j in range(len(exponents)) if i < j)
    diagonal_product = reduce(mul, (matrix[i, i] for i in range(len(exponents))), sp.S.One)
    if triangular:
        return sp.factor(diagonal_product), True
    return sp.factor(matrix.det()), False


def compile_weighted_wronskian(ir: WeightedWronskianIR) -> CompiledDeRhamBlock:
    """Compile a derived band equation to deterministic obstruction data."""

    if not ir.chain.frontend_complete:
        raise ValueError(
            f"chain {ir.chain.name} lacks front-end data: "
            + ", ".join(ir.chain.missing_frontend_data)
        )
    t = ir.t
    a = int(ir.covering_exponent)
    b = int(ir.primitive_weight)
    curve = SuperellipticDeRham(t, ir.A, a, check_squarefree=False)
    variables = sp.symbols(f"d0:{len(ir.primitive_exponents)}")
    solved, compatibility = weighted_wronskian_compatibility(
        t,
        ir.A,
        ir.R,
        a,
        b,
        variables,
        ir.primitive_exponents,
    )
    formal_primitive = sum(
        variable * t**exponent
        for variable, exponent in zip(variables, ir.primitive_exponents)
    )
    primitive = _substitute_in_order(formal_primitive, solved)
    residual = sp.Poly(
        sp.expand(
            a * ir.A * sp.diff(primitive, t)
            - b * sp.diff(ir.A, t) * primitive
            - ir.R
        ),
        t,
    )
    residual_degrees = tuple(
        degree for degree in range(max(0, residual.degree()) + 1) if residual.nth(degree) != 0
    )
    residual_coordinates = tuple(sp.factor(residual.nth(degree)) for degree in residual_degrees)
    if residual_coordinates != compatibility:
        raise ArithmeticError("supported solver residual and compatibility vector disagree")

    low_determinant, triangular = _low_operator_determinant(ir)
    lower, upper = ir.full_primitive_bounds
    expected_tail = tuple(range(upper + 1, upper + 1 + curve.character_dimension(b)))
    delta = curve.delta
    infinity_orders = tuple(
        (curve.n * (a + b) - a * (degree + 1)) // delta - 1
        for degree in residual_degrees
    )
    # At a finite branch point the local exponents differ by a.  A dy/y term
    # can occur only in the invariant character b == 0 mod a.  The displayed
    # infinity orders then decide whether the tail representatives are
    # regular there.
    tail_is_second_kind = b % a != 0 and all(order >= 0 for order in infinity_orders)
    tail_basis_certified = (
        triangular
        and low_determinant != 0
        and tail_is_second_kind
        and residual_degrees == expected_tail
        and len(residual_degrees) == curve.character_dimension(b)
    )
    # Since L_A(D)/(a*y^(a+b)) is exact, the class of R/(a*y^(a+b))
    # has coordinates -residual/a in the certified tail basis.
    de_rham_coordinates = tuple(sp.cancel(-value / a) for value in residual_coordinates)
    identity = sp.expand(
        a * ir.A * sp.diff(primitive, t)
        - b * sp.diff(ir.A, t) * primitive
        - ir.R
        - residual.as_expr()
    )
    certificate = TailObstructionCertificate(
        primitive=primitive,
        solved_coefficients=tuple(solved.items()),
        compatibility=compatibility,
        residual=residual.as_expr(),
        residual_degrees=residual_degrees,
        tail_differentials=tuple(t**degree / sp.Symbol("y") ** (a + b) for degree in residual_degrees),
        tail_infinity_orders=infinity_orders,
        tail_is_second_kind=tail_is_second_kind,
        de_rham_coordinates=de_rham_coordinates,
        low_operator_determinant=low_determinant,
        identity_holds=identity == 0,
        tail_basis_certified=tail_basis_certified,
    )
    character = b % a
    return CompiledDeRhamBlock(
        ir=ir,
        genus=curve.genus,
        character=character,
        affine_dimension=curve.character_dimension(character, compact=False),
        compact_dimension=curve.character_dimension(character, compact=True),
        residue_functional=curve.residue_coefficients(character),
        compact_basis=curve.compact_basis(character),
        certificate=certificate,
    )


def normalized_72_108_block() -> WeightedWronskianIR:
    """Golden compiler fixture from the audited first Laurent block."""

    t = sp.symbols("t")
    a2, a3, a4, a5, a6, a7 = sp.symbols("a2:8")
    A = t + a2 * t**2 + a3 * t**3 + a4 * t**4 + a5 * t**5 + a6 * t**6 + a7 * t**7 + t**8
    chain = NewtonChainIR(
        name="72_108_old_tail",
        corners=((sp.Rational(8), sp.Rational(28)), (sp.Rational(11, 4), sp.Rational(7))),
        multiplicities=(3, 2),
        enumeration_source="GGHV 2022 Proposition 4.3 and audited Laurent transcription",
        status="fully derived and locally replayed",
    )
    return WeightedWronskianIR(
        chain=chain,
        t=t,
        A=A,
        R=t**2,
        covering_exponent=2,
        primitive_weight=3,
        primitive_exponents=tuple(range(2, 13)),
        full_primitive_bounds=(0, 12),
        normalization=("A(0)=0", "A'(0)=1", "A monic of degree 8"),
        scaling_weights=tuple(
            (symbol, index) for index, symbol in enumerate((a2, a3, a4, a5, a6, a7), start=1)
        ),
    )


def repeated_tail_96_144_record() -> NewtonChainIR:
    """The raw table row after its source-level complete-chain exclusion."""

    return NewtonChainIR(
        name="96_144_repeated_tail",
        corners=(
            (sp.Rational(8), sp.Rational(40)),
            (sp.Rational(8), sp.Rational(28)),
            (sp.Rational(11, 4), sp.Rational(7)),
        ),
        multiplicities=(3, 2),
        enumeration_source="GGHV 2017 Section 7 length-two complete-chain table",
        status=(
            "source-excluded before Laurent compilation: simple-root "
            "partitions yield (8,4), and the triple-root translated edge "
            "(8,40)->(8,12) has no complete chain"
        ),
        missing_frontend_data=(
            "no Laurent IR: this repeated-tail row is source-excluded",
        ),
        source_reconciliation=(
            "The 2017 table is a necessary over-approximation whose source "
            "comments out the Proposition-3.29 lower-corner filter. The 2016 "
            "remark claims the whole row leads to the impossible last lower "
            "corner (8,4), while the 2022 proof establishes that transition "
            "for the degree-one companion row (8,32). Reusing its common-tail "
            "calculation gives q1=d0=4 here. The residual vertical factor has "
            "degree three; partitions (2,1) and (1,1,1) contain a simple root "
            "and yield (8,4). The remaining triple-root partition translates "
            "to edge (8,40)->(8,12); the permissive complete-chain enumeration "
            "has open counts 1,6,3,0 and no final corner within maximum length "
            "three. See cas/frontier_96_144_source_audit.py and "
            "cas/complete_chain_no_escape.py"
        ),
    )


def frontier_75_125_record() -> NewtonChainIR:
    """The first frontier row, with forced edges but no exhaustive polygons."""

    return NewtonChainIR(
        name="75_125_F2",
        corners=(
            (sp.Rational(5), sp.Rational(20)),
            (sp.Rational(7, 5), sp.Rational(2)),
        ),
        multiplicities=(3, 5),
        enumeration_source="GGHV 2017 family F2 with j=1",
        status=(
            "forced chain, Puiseux chart, bracket monomial, and terminal edge "
            "are derived; lower Laurent boundary remains unclassified"
        ),
        missing_frontend_data=(
            "proof of support control after y -> y+lambda*x^(-1/5)",
            "exhaustive gamma branches and their complete Laurent polygons",
            "all bands from the common-power layer to bracket layer 4",
        ),
        source_reconciliation=(
            "GGV 2014 Section 5 treats the F2 j=0 degree-(50,75) member "
            "with ratio 2:3; its modified coefficient systems do not supply "
            "the missing 3:5 bands for the j=1 degree-(75,125) member. "
            "The exact facts that do follow from GGHV 2017 are encoded in "
            "cas/f2_75_125_frontend.py"
        ),
    )


if __name__ == "__main__":
    compiled = compile_weighted_wronskian(normalized_72_108_block())
    repeated = repeated_tail_96_144_record()
    frontier = frontier_75_125_record()
    print("genus/character dimension", compiled.genus, compiled.compact_dimension)
    print("tail degrees", compiled.certificate.residual_degrees)
    print("tail basis certified", compiled.certificate.tail_basis_certified)
    print("(96,144) front end complete", repeated.frontend_complete)
    print("(75,125) front end complete", frontier.frontend_complete)
