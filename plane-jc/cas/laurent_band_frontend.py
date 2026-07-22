#!/usr/bin/env python3
"""Exact Laurent-polygon front end for Newton/de Rham blocks.

This module compiles *published Laurent polygons* through a monomial chart to
band supports and bracket layers.  It does not attempt to derive those
polygons from an admissible corner chain; that remains the theorem-heavy
front-end boundary.
"""

from dataclasses import dataclass
from typing import Mapping

import sympy as sp

from newton_derham_compiler import NewtonChainIR, WeightedWronskianIR


LatticePoint = tuple[int, int]


@dataclass(frozen=True)
class ConvexLatticePolygon:
    """A boundary-inclusive convex lattice polygon."""

    vertices: tuple[LatticePoint, ...]

    def contains(self, point: LatticePoint) -> bool:
        signs: list[int] = []
        for first, second in zip(self.vertices, self.vertices[1:] + self.vertices[:1]):
            cross = (second[0] - first[0]) * (point[1] - first[1]) - (
                second[1] - first[1]
            ) * (point[0] - first[0])
            if cross:
                signs.append(1 if cross > 0 else -1)
        return not signs or all(sign == signs[0] for sign in signs)

    def lattice_points(self) -> tuple[LatticePoint, ...]:
        min_x = min(point[0] for point in self.vertices)
        max_x = max(point[0] for point in self.vertices)
        min_y = min(point[1] for point in self.vertices)
        max_y = max(point[1] for point in self.vertices)
        return tuple(
            (x, y)
            for x in range(min_x, max_x + 1)
            for y in range(min_y, max_y + 1)
            if self.contains((x, y))
        )


@dataclass(frozen=True)
class MonomialChart:
    """Images of x and y as monomials in t,z, plus [t,z]_(x,y)."""

    x_image: LatticePoint
    y_image: LatticePoint
    coordinate_jacobian: sp.Expr

    @property
    def exponent_determinant(self) -> int:
        return (
            self.x_image[0] * self.y_image[1]
            - self.x_image[1] * self.y_image[0]
        )

    def transform_exponent(self, point: LatticePoint) -> LatticePoint:
        x_power, y_power = point
        return (
            x_power * self.x_image[0] + y_power * self.y_image[0],
            x_power * self.x_image[1] + y_power * self.y_image[1],
        )

    def bands(self, polygon: ConvexLatticePolygon) -> dict[int, tuple[int, ...]]:
        result: dict[int, list[int]] = {}
        for point in polygon.lattice_points():
            t_power, z_power = self.transform_exponent(point)
            result.setdefault(z_power, []).append(t_power)
        return {
            z_power: tuple(sorted(set(t_powers)))
            for z_power, t_powers in sorted(result.items(), reverse=True)
        }


@dataclass(frozen=True)
class LaurentPairIR:
    """A published Laurent pair with enough data to compile its top block."""

    chain: NewtonChainIR
    P_polygon: ConvexLatticePolygon
    Q_polygon: ConvexLatticePolygon
    chart: MonomialChart
    bracket_rhs_xy_exponent: LatticePoint
    top_P_normalization: tuple[tuple[int, sp.Expr], ...]
    full_primitive_bounds: tuple[int, int]
    scaling_weights: tuple[tuple[sp.Symbol, int], ...] = ()


@dataclass(frozen=True)
class LaurentBandCompilation:
    """Band supports and the compiled weighted-Wronskian leading block."""

    source: LaurentPairIR
    P_bands: tuple[tuple[int, tuple[int, ...]], ...]
    Q_bands: tuple[tuple[int, tuple[int, ...]], ...]
    rhs_tz_exponent: LatticePoint
    weighted_wronskian: WeightedWronskianIR


def bracket_layers(
    t: sp.Symbol,
    P_bands: Mapping[int, sp.Expr],
    Q_bands: Mapping[int, sp.Expr],
    coordinate_jacobian: sp.Expr,
) -> dict[int, sp.Expr]:
    """Compile all z-layers of [P,Q]_(x,y) from band polynomials."""

    layers: dict[int, sp.Expr] = {}
    for i, P_i in P_bands.items():
        for j, Q_j in Q_bands.items():
            layer = i + j - 1
            contribution = coordinate_jacobian * (
                j * sp.diff(P_i, t) * Q_j - i * P_i * sp.diff(Q_j, t)
            )
            layers[layer] = sp.expand(layers.get(layer, 0) + contribution)
    return {
        layer: expression
        for layer, expression in sorted(layers.items(), reverse=True)
        if expression != 0
    }


def _band_polynomial(
    t: sp.Symbol,
    prefix: str,
    exponents: tuple[int, ...],
    fixed: Mapping[int, sp.Expr],
) -> sp.Expr:
    return sp.expand(
        sum(
            fixed.get(exponent, sp.Symbol(f"{prefix}{exponent}")) * t**exponent
            for exponent in exponents
        )
    )


def compile_laurent_leading_block(ir: LaurentPairIR) -> LaurentBandCompilation:
    """Compile polygon supports and the highest bracket layer to de Rham IR."""

    P_bands = ir.chart.bands(ir.P_polygon)
    Q_bands = ir.chart.bands(ir.Q_polygon)
    top_P = max(P_bands)
    top_Q = max(Q_bands)
    top_layer = top_P + top_Q - 1
    rhs_tz = ir.chart.transform_exponent(ir.bracket_rhs_xy_exponent)
    if rhs_tz[1] != top_layer:
        raise ValueError(
            f"bracket RHS lies in z-layer {rhs_tz[1]}, expected top layer {top_layer}"
        )
    if ir.chart.coordinate_jacobian == 0:
        raise ValueError("the monomial chart has zero coordinate Jacobian")
    if ir.chart.coordinate_jacobian.free_symbols:
        raise ValueError("the current band compiler requires a constant coordinate Jacobian")
    if abs(ir.chart.exponent_determinant) != 1:
        raise ValueError("the current band compiler requires a unimodular Laurent chart")

    t = sp.symbols("t")
    fixed = dict(ir.top_P_normalization)
    A = _band_polynomial(t, "a", P_bands[top_P], fixed)
    D = _band_polynomial(t, "d", Q_bands[top_Q], {})
    # [A z^p,D z^q]_(x,y)=(-J)(p A D'-q A'D)z^(p+q-1).
    R = sp.cancel(t**rhs_tz[0] / (-ir.chart.coordinate_jacobian))
    d_symbols = tuple(
        sp.Symbol(f"d{exponent}") for exponent in Q_bands[top_Q]
    )
    compiled = WeightedWronskianIR(
        chain=ir.chain,
        t=t,
        A=A,
        R=R,
        covering_exponent=top_P,
        primitive_weight=top_Q,
        primitive_exponents=Q_bands[top_Q],
        full_primitive_bounds=ir.full_primitive_bounds,
        normalization=tuple(
            f"coefficient of t^{exponent} in top P band is {value}"
            for exponent, value in ir.top_P_normalization
        ),
        scaling_weights=ir.scaling_weights,
    )
    if set(D.free_symbols) - {t} != set(d_symbols):
        raise ArithmeticError("top Q-band symbols do not match its lattice support")
    return LaurentBandCompilation(
        source=ir,
        P_bands=tuple(P_bands.items()),
        Q_bands=tuple(Q_bands.items()),
        rhs_tz_exponent=rhs_tz,
        weighted_wronskian=compiled,
    )


def audited_72_108_laurent_case(case: int) -> LaurentPairIR:
    """Return either published Proposition 4.3 polygon pair."""

    if case not in (1, 2):
        raise ValueError("the audited Laurent case must be 1 or 2")
    P_vertices = [(0, 0), (1, 0), (8, 14), (8, 16)]
    Q_vertices = [(0, 0), (2, 1), (12, 21), (12, 24)]
    if case == 1:
        P_vertices.append((0, 8))
        Q_vertices.append((0, 12))
    chain = NewtonChainIR(
        name=f"72_108_laurent_case_{case}",
        corners=((sp.Rational(8), sp.Rational(28)), (sp.Rational(11, 4), sp.Rational(7))),
        multiplicities=(3, 2),
        enumeration_source="GGHV 2022 Proposition 4.3",
        status="published Laurent polygons and locally audited transcription",
    )
    a2, a3, a4, a5, a6, a7 = sp.symbols("a2:8")
    return LaurentPairIR(
        chain=chain,
        P_polygon=ConvexLatticePolygon(tuple(P_vertices)),
        Q_polygon=ConvexLatticePolygon(tuple(Q_vertices)),
        chart=MonomialChart((1, 2), (0, -1), sp.S.NegativeOne),
        bracket_rhs_xy_exponent=(2, 0),
        top_P_normalization=((1, sp.S.One), (8, sp.S.One)),
        full_primitive_bounds=(0, 12),
        scaling_weights=tuple(
            (symbol, index)
            for index, symbol in enumerate((a2, a3, a4, a5, a6, a7), start=1)
        ),
    )


if __name__ == "__main__":
    for case in (1, 2):
        result = compile_laurent_leading_block(audited_72_108_laurent_case(case))
        print(
            "case",
            case,
            "lattice points",
            len(result.source.P_polygon.lattice_points()),
            len(result.source.Q_polygon.lattice_points()),
            "top bands",
            result.P_bands[0],
            result.Q_bands[0],
        )
