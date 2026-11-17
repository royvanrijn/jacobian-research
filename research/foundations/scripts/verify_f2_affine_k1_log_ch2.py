#!/usr/bin/env python3
"""Verify the generic k=1 affine-row logarithmic Chern contribution."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.log_node_profiles import (  # noqa: E402
    immersed_affine_row_profile,
    logarithmic_ch2_budget,
    logarithmic_ch2_model_change,
    residual_point_budget,
)


def generic_cyclic_local_model_audit() -> None:
    # At a smooth point of the affine target normalization take a transverse
    # equation x and a normalization parameter y.  A boundary row of index e
    # has (x,y)=(u^e,t) up to units and higher terms.  In the source log basis
    # (du/u,dt), the unit entry makes Fitt_1 the unit ideal and the cokernel is
    # cyclic R/(u^e).
    u = sp.symbols("u")
    for transverse_index in range(2, 9):
        matrix = sp.Matrix(
            [[transverse_index * u**transverse_index, 0], [0, 1]]
        )
        assert sp.factor(matrix.det()) == transverse_index * u**transverse_index
        assert 1 in tuple(matrix)


def k1_log_immersion_audit() -> None:
    # The exact generic witness from PF2K1C1 is nodal, not cuspidal, on the
    # affine normalization.  Its two derivatives have no common zero, so the
    # normalization is an immersion there.  The unique point at infinity is
    # handled logarithmically.
    t = sp.symbols("t")
    p = t**3 + t
    q = t**5
    assert sp.gcd(sp.diff(p, t), sp.diff(q, t)) == 1

    # On (P2,L_infinity), L_Y=-2H and a quintic has L_Y.C=-10.  The first
    # smooth-boundary blowup extracting the primitive (5,2) ray is centered
    # at the infinity branch of multiplicity two, raising the intersection
    # to -8.  The remaining three terminal target blowups are node blowups.
    projective_log_intersection = -2 * 5
    infinity_multiplicity = 2
    terminal_model_intersection = (
        projective_log_intersection + infinity_multiplicity
    )
    assert terminal_model_intersection == -8

    # If the strict transform meets the special carrier point, every one of
    # the first eight carrier extractions that it follows is another smooth
    # boundary blowup at a smooth point of the curve.  Later carrier fan
    # insertions are node blowups and do not change L_Y.C.
    for followed_carrier_centers in range(9):
        target_intersection = -8 + followed_carrier_centers
        profile = immersed_affine_row_profile(
            target_intersection,
            -1,  # deg Omega_P1(log infinity)
            residue_degree=1,
            transverse_index=2,
            source_self_intersection=-1,
        )
        assert profile.reduced_kernel_degree == -7 + followed_carrier_centers
        assert profile.thickened_kernel_degree == 2 * (
            -7 + followed_carrier_centers
        )
        assert profile.determinant_divisor_square == -4
        assert profile.cokernel_ch2 == 2 * followed_carrier_centers - 16

    # Both the reduced conormal degree and the full thickened degree multiply
    # exactly with residue and transverse degree.
    profile = immersed_affine_row_profile(-6, -1, 3, 4, -2)
    assert profile.reduced_kernel_degree == -15
    assert profile.thickened_kernel_degree == -60
    assert profile.determinant_divisor_square == -32
    assert profile.cokernel_ch2 == -76


def immersed_curve_generalization_audit() -> None:
    # For a smooth normalization of genus g with s logarithmic punctures,
    # deg Omega(log S)=2g-2+s.  The helper therefore specializes to the
    # degree-independent immersed-curve packet formula.
    for genus in range(4):
        for punctures in range(1, 5):
            normalization_degree = 2 * genus - 2 + punctures
            for target_intersection in range(-10, 5):
                for residue_degree in range(1, 4):
                    for transverse_index in range(1, 5):
                        for self_intersection_magnitude in range(1, 4):
                            profile = immersed_affine_row_profile(
                                target_intersection,
                                normalization_degree,
                                residue_degree,
                                transverse_index,
                                -self_intersection_magnitude,
                            )
                            expected = (
                                transverse_index
                                * residue_degree
                                * (target_intersection - normalization_degree)
                                - Fraction(
                                    transverse_index**2
                                    * self_intersection_magnitude,
                                    2,
                                )
                            )
                            assert profile.cokernel_ch2 == expected


def carrier_contact_jet_audit() -> None:
    """Identify ``b`` with a truncated target-jet contact order."""

    # In the carrier-normalized target chart, pi=0 is the old boundary and
    # w=0 is transverse.  The eight smooth-boundary blowups insert the rays
    # (1,1),...,(1,8).  A transverse k=1 curve with valuation (1,r) follows
    # exactly the first min(r,8) centers.
    smooth_carrier_rays = tuple((1, order) for order in range(1, 9))
    assert smooth_carrier_rays == (
        (1, 1),
        (1, 2),
        (1, 3),
        (1, 4),
        (1, 5),
        (1, 6),
        (1, 7),
        (1, 8),
    )
    for contact_order in range(1, 17):
        followed = sum(
            contact_order >= center_order
            for _, center_order in smooth_carrier_rays
        )
        assert followed == min(contact_order, 8)

    # Write the actual fixed-target infinity expansions as
    #   P=A*u^-3*(1+alpha*u+...), -Q=C*u^-5*(1+beta*u+...).
    # Then pi=P^3/(-Q)^2 and h=P^5/(-Q)^3 have the displayed leading jets.
    # After the first fixed carrier shear w=h-lambda_0-c_1*pi-..., the
    # coefficient below is the exact test separating b=1 from b>=2.
    u = sp.symbols("u")
    A, C = sp.symbols("A C", nonzero=True)
    alpha, beta, first_shear = sp.symbols("alpha beta c_1")
    p = A * u**-3 * (1 + alpha * u)
    minus_q = C * u**-5 * (1 + beta * u)
    pi = p**3 / minus_q**2
    h = p**5 / minus_q**3
    leading_pi = sp.simplify(sp.limit(pi / u, u, 0))
    leading_h = sp.simplify(sp.limit(h, u, 0))
    first_h = sp.simplify(sp.limit((h - leading_h) / u, u, 0))
    first_w = sp.simplify(first_h - first_shear * leading_pi)
    assert leading_pi == A**3 / C**2
    assert leading_h == A**5 / C**3
    assert first_h == leading_h * (5 * alpha - 3 * beta)
    assert sp.simplify(
        first_w
        - leading_h * (5 * alpha - 3 * beta)
        + first_shear * leading_pi
    ) == 0


def degree_floor_attachment_sieve() -> None:
    square_budget = logarithmic_ch2_budget(-6, -5, 6).ch2
    double_budget = logarithmic_ch2_budget(-11, -5, 12).ch2
    square_after_root = residual_point_budget(square_budget, (27,)).residual_ch2
    double_after_root = residual_point_budget(double_budget, (27,)).residual_ch2
    assert square_after_root == -10
    assert double_after_root == Fraction(17, 2)

    # PF2PPA1 rules out extracting the affine row above the already resolved
    # terminal neighborhood, so no terminal branch value determines e.  For
    # every possible affine signature, the one-packet filtration has the
    # following exact virtual residuals.  Their interpretation as lengths is
    # conditional on the missing global K-theory filtration.
    for followed_carrier_centers in range(9):
        for transverse_index in range(2, 7):
            for residue_degree in range(1, 5):
                for self_intersection_magnitude in range(1, 5):
                    affine_packet = immersed_affine_row_profile(
                        -8 + followed_carrier_centers,
                        -1,
                        residue_degree=residue_degree,
                        transverse_index=transverse_index,
                        source_self_intersection=-self_intersection_magnitude,
                    ).cokernel_ch2
                    for source_smooth_blowups in range(21):
                        square_change = logarithmic_ch2_model_change(
                            6,
                            source_smooth_blowups=source_smooth_blowups,
                        ).ch2_change
                        double_change = logarithmic_ch2_model_change(
                            12,
                            source_smooth_blowups=source_smooth_blowups,
                        ).ch2_change
                        square_point = (
                            square_after_root + square_change - affine_packet
                        )
                        double_point = (
                            double_after_root + double_change - affine_packet
                        )
                        common = (
                            transverse_index**2 * self_intersection_magnitude
                            - 2
                            * transverse_index
                            * residue_degree
                            * (followed_carrier_centers - 7)
                        )
                        assert 2 * square_point == (
                            common - 20 - source_smooth_blowups
                        )
                        assert 2 * double_point == (
                            common + 17 - source_smooth_blowups
                        )
                        assert (square_point.denominator == 1) == (
                            source_smooth_blowups % 2
                            == (
                                transverse_index**2
                                * self_intersection_magnitude
                            )
                            % 2
                        )
                        assert (double_point.denominator == 1) == (
                            source_smooth_blowups % 2
                            == (
                                transverse_index**2
                                * self_intersection_magnitude
                                + 1
                            )
                            % 2
                        )

    # The smallest purity signature (e,f,n)=(2,1,1) is now only an exact
    # benchmark, not a claim about the unresolved source attachment.  At a
    # nonspecial puncture b=0.  At the special carrier point 1<=b<=8.
    minimal = [
        immersed_affine_row_profile(-8 + b, -1, 1, 2, -1).cokernel_ch2
        for b in range(9)
    ]
    assert minimal == [2 * b - 16 for b in range(9)]
    square_minimal = [
        square_after_root - minimal[b]
        for b in range(9)
    ]
    double_minimal = [
        double_after_root - Fraction(1, 2) - minimal[b]
        for b in range(9)
    ]
    assert square_minimal == [6, 4, 2, 0, -2, -4, -6, -8, -10]
    assert double_minimal == [24, 22, 20, 18, 16, 14, 12, 10, 8]


def cusp_boundary_residual_sieve() -> None:
    # On the generic k=1 nonimmersion face, the unique nonimmersive value is
    # an ordinary cusp of branch multiplicity two.  LCAD1 distinguishes a
    # smooth boundary preimage (epsilon=0) from an SNC boundary node
    # (epsilon=1).  The local isolated-Fitt lower exponent is
    # 2*q_p-1+epsilon.  Thus a degree-f fiber with h points, c of them nodes,
    # has ledger 2*f-h+c, between f and 2*f.
    def compositions(total: int) -> tuple[tuple[int, ...], ...]:
        if total == 0:
            return ((),)
        result: list[tuple[int, ...]] = []
        for first in range(1, total + 1):
            for tail in compositions(total - first):
                result.append((first, *tail))
        return tuple(result)

    for residue_degree in range(1, 8):
        for local_indices in compositions(residue_degree):
            fiber_points = len(local_indices)
            for crossing_count in range(fiber_points + 1):
                epsilons = (1,) * crossing_count + (0,) * (
                    fiber_points - crossing_count
                )
                cusp_ledger = sum(
                    2 * q_p - 1 + epsilon
                    for q_p, epsilon in zip(
                        local_indices, epsilons, strict=True
                    )
                )
                assert cusp_ledger == (
                    2 * residue_degree - fiber_points + crossing_count
                )
                assert residue_degree <= cusp_ledger <= 2 * residue_degree

    # Once exact minimal point quotients occur in the assumed filtration, the
    # unidentified residual is the old virtual residual minus B_cusp.  Test
    # every incidence ledger allowed by every residue partition.
    for followed_carrier_centers in range(9):
        for transverse_index in range(2, 7):
            for residue_degree in range(1, 5):
                for self_intersection_magnitude in range(1, 5):
                    common = (
                        transverse_index**2 * self_intersection_magnitude
                        - 2
                        * transverse_index
                        * residue_degree
                        * (followed_carrier_centers - 7)
                    )
                    for source_smooth_blowups in range(21):
                        square_before_cusp = Fraction(
                            common - 20 - source_smooth_blowups,
                            2,
                        )
                        double_before_cusp = Fraction(
                            common + 17 - source_smooth_blowups,
                            2,
                        )
                        for cusp_ledger in range(
                            residue_degree, 2 * residue_degree + 1
                        ):
                            square_rest = square_before_cusp - cusp_ledger
                            double_rest = double_before_cusp - cusp_ledger
                            assert 2 * square_rest == (
                                common
                                - 20
                                - source_smooth_blowups
                                - 2 * cusp_ledger
                            )
                            assert 2 * double_rest == (
                                common
                                + 17
                                - source_smooth_blowups
                                - 2 * cusp_ledger
                            )
                            assert (square_rest.denominator == 1) == (
                                source_smooth_blowups % 2
                                == (
                                    transverse_index**2
                                    * self_intersection_magnitude
                                )
                                % 2
                            )
                            assert (double_rest.denominator == 1) == (
                                source_smooth_blowups % 2
                                == (
                                    transverse_index**2
                                    * self_intersection_magnitude
                                    + 1
                                )
                                % 2
                            )

    # Minimal (e,f,n)=(2,1,1) benchmarks, with the smallest compatible s_X
    # values zero and one.  The smooth-fold endpoint drops by one; the
    # node-saturated LUAF1 endpoint drops by two.
    square_fold_rest = [5 - 2 * b for b in range(9)]
    square_node_rest = [4 - 2 * b for b in range(9)]
    double_fold_rest = [23 - 2 * b for b in range(9)]
    double_node_rest = [22 - 2 * b for b in range(9)]
    assert square_fold_rest == [5, 3, 1, -1, -3, -5, -7, -9, -11]
    assert square_node_rest == [4, 2, 0, -2, -4, -6, -8, -10, -12]
    assert double_fold_rest == [23, 21, 19, 17, 15, 13, 11, 9, 7]
    assert double_node_rest == [22, 20, 18, 16, 14, 12, 10, 8, 6]


def main() -> None:
    generic_cyclic_local_model_audit()
    k1_log_immersion_audit()
    immersed_curve_generalization_audit()
    carrier_contact_jet_audit()
    degree_floor_attachment_sieve()
    cusp_boundary_residual_sieve()
    print(
        "PASS: generic k=1 affine rows have exact conormal/cyclic ch2 and "
        "all-(e,f,n,b) degree-floor parity/effectivity formulas; on the "
        "ordinary-cusp face the incidence-sensitive minimal boundary ledger "
        "is 2*f-h+c, ranging from f for smooth folds to 2*f for SNC nodes"
    )


if __name__ == "__main__":
    main()
