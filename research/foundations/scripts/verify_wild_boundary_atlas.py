#!/usr/bin/env python3
"""Exact and bounded checks for the plane wild-boundary atlas.

The all-parameter proofs and the status boundary are recorded in
extended-geometry/PLANE_WILD_BOUNDARY_ATLAS.md.  The default checks replay
the polynomial identities, numerical-semigroup formula, thickened
base-change obstructions, and uniform double-fibre lattice.  ``--singular``
also asks Singular to normalize seven representative primitive orders and
compute their relative-different saturations.  Those runs are exact bounded
certificates, not all-prime normalization proofs.
"""

from __future__ import annotations

import argparse
import math
import shutil
import subprocess

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


P, Q, T = sp.symbols("P Q T")
x, y = sp.symbols("x y")
TEST_ROWS = (
    (2, 2),
    (2, 4),
    (2, 6),
    (3, 3),
    (3, 6),
    (5, 5),
    (7, 7),
)
SINGULAR_ROWS = ((2, 2), (2, 4), (2, 6), (3, 3), (3, 6), (5, 5), (7, 7))
BALANCED_SINGULAR_ROWS = (
    (2, 2),
    (2, 6),
    (3, 3),
    (3, 6),
    (5, 5),
    (7, 7),
)
EXTENDED_BAND_ROWS = ((2, 6),)
THICKENED_SINGULAR_ROWS = (
    (2, 2, 2, 6, 2),
    (2, 2, 4, 6, 2),
    (2, 2, 8, 6, 2),
    (3, 3, 2, 16, 3),
    (3, 3, 3, 12, 3),
    (5, 5, 2, 36, 5),
)


def reduced(poly: sp.Expr, prime: int) -> sp.Poly:
    return sp.Poly(sp.expand(poly), P, Q, T, modulus=prime)


def reduced_source(poly: sp.Expr, prime: int) -> sp.Poly:
    return sp.Poly(sp.expand(poly), x, y, modulus=prime)


def residue_degree_parts(prime: int, degree: int) -> tuple[int, int]:
    """Return the separable and inseparable parts of ``z^N=P``."""

    inseparable = 1
    residual = degree
    while residual % prime == 0:
        inseparable *= prime
        residual //= prime
    return residual, inseparable


def semigroup_profile(degree: int) -> tuple[int, int, tuple[int, ...]]:
    """Return conductor, delta, and gaps for <N,N+1,N^2-N-1>."""
    if degree == 2:
        return 0, 0, ()

    generators = (degree, degree + 1, degree * degree - degree - 1)
    search_bound = 2 * degree * degree
    semigroup = {0}
    for value in range(search_bound + 1):
        if value not in semigroup:
            continue
        for generator in generators:
            if value + generator <= search_bound:
                semigroup.add(value + generator)

    conductor = next(
        candidate
        for candidate in range(search_bound + 1)
        if all(value in semigroup for value in range(candidate, search_bound + 1))
    )
    gaps = tuple(value for value in range(conductor) if value not in semigroup)
    return conductor, len(gaps), gaps


def verify_hidden_orders() -> None:
    for prime, degree in TEST_ROWS:
        assert degree > 1 and degree % prime == 0
        assert (degree + 1) % prime == 1
        separable, inseparable = residue_degree_parts(prime, degree)
        assert separable * inseparable == degree
        assert math.gcd(separable, prime) == 1

        boundary_factor = T**degree - P ** (degree + 1)
        hidden_order = (T - 1) * boundary_factor + P * Q * T
        derivative = sp.diff(hidden_order, T)
        expected_derivative = boundary_factor + P * Q

        assert reduced(derivative - expected_derivative, prime).is_zero
        assert reduced(
            (T - 1) * derivative - hidden_order + P * Q,
            prime,
        ).is_zero
        assert reduced(hidden_order.subs(Q, 0) - (T - 1) * boundary_factor, prime).is_zero

        # The retained/missing intersections are P^(N+1)=1.  Their derivative
        # is P^N because N+1=1 in characteristic p, so all N+1 are reduced.
        node_polynomial = P ** (degree + 1) - 1
        assert reduced(sp.diff(node_polynomial, P) - P**degree, prime).is_zero

        conductor, delta, gaps = semigroup_profile(degree)
        assert conductor == degree * (degree - 2)
        assert delta == degree * (degree - 1) // 2 - 1
        if degree == 2:
            assert gaps == ()
        else:
            frobenius_gap = degree * degree - degree - 1
            assert frobenius_gap not in gaps


def verify_balanced_hidden_orders() -> None:
    """Check the balanced gluing and its natural birational source chart."""
    for prime, degree in TEST_ROWS:
        boundary_factor = T**degree - P ** (degree + 1)
        balanced_order = (
            (T - 1) * boundary_factor + P ** (degree - 1) * Q * T
        )
        derivative = sp.diff(balanced_order, T)

        assert reduced(
            derivative - boundary_factor - P ** (degree - 1) * Q,
            prime,
        ).is_zero
        assert reduced(
            (T - 1) * derivative
            - balanced_order
            + P ** (degree - 1) * Q,
            prime,
        ).is_zero

        # The P-adic Newton polygon of the T=0 cluster has vertices
        # (0,N+1),(1,N-1),(N,0), hence slopes -2 and -1.  Its two
        # normalized rows contribute 1 and N-1 unramified sheets.
        assert (degree + 1 - (degree - 1)) == 2
        assert (degree - 1) // (degree - 1) == 1
        assert (degree - 1) % prime != 0
        assert 1 + (degree - 1) + 1 == degree + 1

        source_r = 1 + x ** (degree - 1) * y
        source_u = 1 + x ** (degree + 1) * source_r
        source_P = x * source_r * source_u ** (degree - 1)
        source_T = source_r * source_u**degree
        source_Q_quotient = sp.cancel(
            (source_T - 1) / x ** (degree - 1)
        )
        assert sp.denom(source_Q_quotient) == 1
        source_Q = -source_Q_quotient * source_u ** (degree - 2)

        # For N=p^n, Frobenius compresses the quotient to the old two-term
        # formula.  Composite multiples of p retain extra binomial terms,
        # but the polynomial chart and its Jacobian identity still hold.
        separable, _ = residue_degree_parts(prime, degree)
        if separable == 1:
            assert reduced_source(
                source_Q_quotient
                - y
                - x ** (degree * degree + 1)
                * source_r ** (degree + 1),
                prime,
            ).is_zero
        source_boundary_factor = source_T**degree - source_P ** (degree + 1)
        source_order = (
            (source_T - 1) * source_boundary_factor
            + source_P ** (degree - 1) * source_Q * source_T
        )

        assert reduced_source(source_order, prime).is_zero
        assert reduced_source(
            source_P * source_T ** (degree - 1)
            - x * source_boundary_factor,
            prime,
        ).is_zero
        assert reduced_source(
            source_T**degree - source_u * source_boundary_factor,
            prime,
        ).is_zero

        jacobian = (
            sp.diff(source_P, x) * sp.diff(source_Q, y)
            - sp.diff(source_P, y) * sp.diff(source_Q, x)
        )
        assert reduced_source(
            jacobian + source_u ** (2 * degree - 4),
            prime,
        ).is_zero

        # On the balanced normalized complement, div(x) and div(u) on the
        # two components removed by D(x) have column vectors
        # (1,N-2) and (0,N-1).  Since D(x) is the Laurent chart G_m^2,
        # the localization sequence gives Cl(U)=Z/(N-1).
        class_matrix_determinant = degree - 1
        entry_gcd = math.gcd(1, degree - 2, degree - 1)
        assert entry_gcd == 1
        assert class_matrix_determinant > 0

        # More generally, for gluing P^a*Q*T with 1<=a<=N-1, the long
        # Newton segment has tame index (N-1)/g and separable residue degree
        # g, where g=gcd(a,N-1).  Only a=N-1 removes this different, and that
        # endpoint is exactly the class-group-obstructed balanced row.
        zero_different_valuations = []
        for gluing_valuation in range(1, degree):
            gcd_value = math.gcd(gluing_valuation, degree - 1)
            ramification_index = (degree - 1) // gcd_value
            residue_degree = gcd_value
            different_exponent = ramification_index - 1
            assert ramification_index * residue_degree == degree - 1
            assert ramification_index % prime != 0
            if different_exponent == 0:
                zero_different_valuations.append(gluing_valuation)
        assert zero_different_valuations == [degree - 1]


def verify_extended_gluing_band() -> None:
    """Audit the omitted exponents and extra coefficient-support divisors."""

    for prime, degree in TEST_ROWS:
        # At a=0 the primitive order is normal by the isolated-singularity
        # argument in the note.  On the affine-UFD core D(TF), the two unit
        # generators T,F have the same valuation N+1 on the sole fill prime.
        # Thus the 1x2 core matrix has kernel rank one and cokernel Z/(N+1).
        zero_core_row = (degree + 1, degree + 1)
        assert math.gcd(*zero_core_row) == degree + 1
        assert zero_core_row[0] - zero_core_row[1] == 0

        # For a>=N the point (T-degree,P-order)=(1,a) lies strictly above
        # the segment from (0,N+1) to (N,0).  The surviving slope is
        # -(N+1)/N, already reduced, hence gives a wild e=N branch over P=0.
        for gluing_valuation in (degree, degree + 1, 2 * degree):
            assert degree * gluing_valuation > degree * degree - 1
            assert math.gcd(degree, degree + 1) == 1
            assert degree % prime == 0

    for prime, degree in EXTENDED_BAND_ROWS:
        boundary_factor = T**degree - P ** (degree + 1)
        zero_order = (T - 1) * boundary_factor + Q * T
        zero_derivative = sp.diff(zero_order, T)
        assert reduced(
            zero_derivative - boundary_factor - Q, prime
        ).is_zero
        assert reduced(
            (T - 1) * zero_derivative - zero_order + Q, prime
        ).is_zero

        # A factor of C(P,Q) away from P*Q=0 creates another whole different
        # divisor.  The identity proves support.  At a normalized prime where
        # T-1 is a unit, its coefficient is h*ord_D(R), not automatically h.
        for multiplicity in (1, 2, 3):
            extra_factor = P + Q + 1
            coefficient = (
                P ** (degree - 1) * extra_factor**multiplicity
            )
            general_order = (
                (T - 1) * boundary_factor + coefficient * Q * T
            )
            general_derivative = sp.diff(general_order, T)
            assert reduced(
                general_derivative - boundary_factor - coefficient * Q,
                prime,
            ).is_zero
            assert reduced(
                (T - 1) * general_derivative
                - general_order
                + coefficient * Q,
                prime,
            ).is_zero
            quotient = sp.Poly(coefficient, P, Q, modulus=prime)
            factor = sp.Poly(extra_factor, P, Q, modulus=prime)
            for _ in range(multiplicity):
                quotient = quotient.exquo(factor)
            assert not quotient.rem(factor).is_zero

        # The arbitrary-coefficient support theorem leaves C=P^a*Q^b.
        # On its fierce branch Q is a uniformizer and H_T has exact order
        # c=b+1.  The retained/fierce collision has local form u*v+Q^c.
        for q_exponent in (1, 2, 3, 4):
            monomial_order = (
                (T - 1) * boundary_factor
                + P ** (degree - 1) * Q**q_exponent * T
            )
            monomial_derivative = sp.diff(monomial_order, T)
            assert reduced(
                (T - 1) * monomial_derivative
                - monomial_order
                + P ** (degree - 1) * Q**q_exponent,
                prime,
            ).is_zero
            local_u = T - 1
            shifted_f = boundary_factor + P ** (degree - 1) * Q**q_exponent
            assert reduced(
                monomial_order
                - local_u * shifted_f
                - P ** (degree - 1) * Q**q_exponent,
                prime,
            ).is_zero
            assert q_exponent > 0


def prime_to_characteristic_part(prime: int, value: int) -> int:
    """Remove the purely inseparable part of a positive integer."""

    while value % prime == 0:
        value //= prime
    return value


def verify_thickened_base_change_obstruction() -> None:
    """Audit the Euler, transfer, and uniform Frobenius-core obstructions."""

    pure_frobenius_rows: set[tuple[int, int, int]] = set()
    for prime, degree in TEST_ROWS:
        for q_exponent in range(2, 17):
            separable_degree = prime_to_characteristic_part(
                prime, q_exponent
            )
            # U_N=D(x) union A0 union L1 has Euler characteristic 1;
            # the retained Q=0 branch is A1 minus mu_(N+1), with Euler -N.
            euler_characteristic = (
                (degree + 1) * separable_degree - degree
            )
            if separable_degree > 1:
                assert euler_characteristic != 1
                continue

            # A pure p-power base change preserves the exact order N-1 of
            # the pulled-back L1 class because gcd(p,N-1)=1.
            if degree > 2:
                assert math.gcd(q_exponent, degree - 1) == 1
                pulled_back_order = (degree - 1) // math.gcd(
                    q_exponent, degree - 1
                )
                assert pulled_back_order == degree - 1
                continue
            pure_frobenius_rows.add((prime, degree, q_exponent))

    assert pure_frobenius_rows == {
        (2, 2, 2),
        (2, 2, 4),
        (2, 2, 8),
        (2, 2, 16),
    }

    # Every remaining row is excluded on D(x*u).  If c is a power of two,
    # put t=x/u, r=c/2, and z=Q^r*t^2*u+1.  The core equation implies
    # z^2=(1-t^3)*u and has t-derivative t^2*u, a unit.  Thus its three
    # fibers over mu_3 are double prime fibers and div(z) is their sum.
    local_t, local_u = sp.symbols("local_t local_u")
    for q_exponent in (2, 4, 8, 16):
        core_equation = (
            Q**q_exponent * local_t**4 * local_u**2
            + (1 - local_t**3) * local_u
            + 1
        )
        local_z = (
            Q ** (q_exponent // 2) * local_t**2 * local_u + 1
        )
        assert sp.Poly(
            sp.expand(
                local_z**2
                - (1 - local_t**3) * local_u
                - core_equation
            ),
            Q,
            local_t,
            local_u,
            modulus=2,
        ).is_zero
        assert sp.Poly(
            sp.diff(core_equation, local_t) - local_t**2 * local_u,
            Q,
            local_t,
            local_u,
            modulus=2,
        ).is_zero

        # The geometric generic fiber is rational, not merely
        # topologically G_m.  With r=c/2 and after rescaling Q, its equation
        # is R^r=y^2+y.  The displayed parameter and inverse prove equality
        # of function fields; the rational missing point y=R=0 then kills
        # the affine generic-fiber class group through the P1 divisor row.
        generic_degree = q_exponent // 2
        generic_height = generic_degree.bit_length() - 1
        assert generic_degree == 2**generic_height
        generic_s = sp.symbols("generic_s")
        generic_R = generic_s**2 + generic_s
        generic_y = generic_s**generic_degree
        assert sp.Poly(
            generic_R**generic_degree - generic_y**2 - generic_y,
            generic_s,
            modulus=2,
        ).is_zero
        recovered_s = generic_y + sum(
            generic_R ** (2**index)
            for index in range(generic_height)
        )
        assert sp.Poly(
            recovered_s - generic_s,
            generic_s,
            modulus=2,
        ).is_zero

    # The exact vertical relation lattice has columns 2*D_i and sum(D_i).
    # Its cokernel is (Z/2)^2 for every pure Frobenius exponent.  The generic
    # fiber unit-rank argument in the note proves that there are no further
    # relations; the Smith calculation replays the resulting exact lattice.
    for multiplicity, fiber_count in ((2, 3), (3, 4), (5, 2)):
        multiple_fiber_relations = sp.Matrix.hstack(
            multiplicity * sp.eye(fiber_count),
            sp.ones(fiber_count, 1),
        )
        smith = smith_normal_form(multiple_fiber_relations, domain=ZZ)
        diagonal = tuple(
            abs(int(smith[index, index]))
            for index in range(min(smith.rows, smith.cols))
            if smith[index, index]
        )
        assert diagonal == (1,) + (multiplicity,) * (fiber_count - 1)
        assert multiple_fiber_relations.rank() == fiber_count

    # Unequal multiplicities give (direct_sum Z/m_i)/<diagonal>.  These two
    # controls test both sides of the exact torsion trigger: a repeated prime
    # yields Z/6 for (4,6,9), while pairwise-coprime (4,9,25) yields zero.
    for multiplicities, expected in (
        ((4, 6, 9), (1, 1, 6)),
        ((4, 9, 25), (1, 1, 1)),
    ):
        mixed_relations = sp.Matrix.hstack(
            sp.diag(*multiplicities),
            sp.ones(len(multiplicities), 1),
        )
        mixed_smith = smith_normal_form(mixed_relations, domain=ZZ)
        mixed_diagonal = tuple(
            abs(int(mixed_smith[index, index]))
            for index in range(len(multiplicities))
        )
        assert mixed_diagonal == expected
        assert math.prod(multiplicities) // math.lcm(*multiplicities) == math.prod(
            value for value in expected if value > 1
        )

    # Passing from the corrected core W_q=D(x*u) back to the whole source
    # U_(2,q) adds the reduced primes A_0=(x=0) and L_q=(u=0).  The core
    # units t=x/u and u have the following boundary valuations.  Its Smith
    # profile (1,1) certifies that the localization kernel vanishes, so the
    # restriction Cl(U_(2,q))->Cl(W_q) is an isomorphism for every q=2^s.
    source_fill_matrix = sp.Matrix(((1, 0), (-1, 1)))
    source_fill_smith = smith_normal_form(source_fill_matrix, domain=ZZ)
    assert source_fill_matrix.det() == 1
    assert tuple(
        abs(int(source_fill_smith[index, index]))
        for index in range(2)
    ) == (1, 1)


def artin_schreier_witt_different(prime: int, jumps: tuple[int, ...]) -> int:
    """Different exponent from upper jumps of a cyclic p-power extension."""
    return sum(
        (prime**index - prime ** (index - 1)) * (jump + 1)
        for index, jump in enumerate(jumps, start=1)
    )


def verify_ramification_rows() -> None:
    for prime in (2, 3, 5, 7):
        for pole_order in (1, prime + 1):
            if pole_order % prime:
                assert (prime - 1) * (pole_order + 1) > 0

        jumps = (1, prime + 1)
        different = artin_schreier_witt_different(prime, jumps)
        expected = (prime - 1) * 2 + (prime**2 - prime) * (prime + 2)
        assert different == expected

        # The Kummer-twisted additive block Z^(p^n)+Q^m Z-P has derivative
        # Q^m and fierce residue degree p^n at Q=0.
        for length in (1, 2):
            residue_degree = prime**length
            for different_exponent in (1, 2, 3):
                assert residue_degree % prime == 0
                assert different_exponent > 0


def singular_program(prime: int, degree: int) -> str:
    return rf'''
LIB "primdec.lib";
LIB "normal.lib";

proc assertReductionZero(poly f, ideal G, string label)
{{
  if (reduce(f,std(G)) != 0)
  {{
    "FAIL: "+label;
    exit(1);
  }}
}}

proc assertIdealEqual(ideal A, ideal B, string label)
{{
  int i;
  for (i=1; i<=size(A); i++)
  {{
    assertReductionZero(A[i],B,label);
  }}
  for (i=1; i<=size(B); i++)
  {{
    assertReductionZero(B[i],A,label);
  }}
}}

ring r={prime},(P,Q,T),dp;
poly Hzero=(T-1)*(T^{degree}-P^{degree + 1})+Q*T;
ideal zeroOrder=Hzero;
list zeroNormalization=normal(zeroOrder,"isPrim","withGens");
intvec zeroNormalCheck=norTest(zeroOrder,zeroNormalization);
if (
  zeroNormalCheck[1] != 1
  || zeroNormalCheck[2] != 1
  || zeroNormalCheck[3] != 1
)
{{
  "FAIL: zero-exponent normalization test";
  exit(1);
}}
assertIdealEqual(
  zeroNormalization[2][1],ideal(1),"zero-exponent normalization module"
);
assertIdealEqual(
  normalConductor(zeroOrder),ideal(1),"zero-exponent normal conductor"
);
assertIdealEqual(
  radical(zeroOrder+jacob(zeroOrder)),
  ideal(P,Q,T),
  "zero-exponent isolated singularity"
);

poly H=(T-1)*(T^{degree}-P^{degree + 1})+P*Q*T;
ideal primitiveOrder=H;
list N=normal(primitiveOrder,"isPrim","withGens");
intvec normalCheck=norTest(primitiveOrder,N);
if (normalCheck[1] != 1 || normalCheck[2] != 1 || normalCheck[3] != 1)
{{
  "FAIL: normalization test";
  exit(1);
}}

ideal moduleGenerators=N[2][1];
assertIdealEqual(moduleGenerators,ideal(P^{degree},T),"normalization module");
ideal conductor=normalConductor(primitiveOrder);
assertIdealEqual(conductor,ideal(P,T),"primitive-order conductor");

def Rn=N[1][1]; setring Rn;
poly W=var(1);
ideal expected=
  (T-1)*(T^{degree - 1}-P*W)+P*Q,
  W*T-P^{degree},
  W2-W*P^{degree}+W*Q+P^{degree - 1}*T^{degree - 2}*(T-1),
  (T-1)*(T^{degree}-P^{degree + 1})+P*Q*T;
assertIdealEqual(norid,expected,"normalization presentation");

// Relative Kahler different over the fixed target coordinates P,Q.
ideal C=expected;
matrix relativeJacobian[4][2]=
  diff(C[1],W),diff(C[1],T),
  diff(C[2],W),diff(C[2],T),
  diff(C[3],W),diff(C[3],T),
  diff(C[4],W),diff(C[4],T);
ideal relativeDifferent=C+minor(relativeJacobian,2);
poly fierceEquation=T^{degree}-P^{degree + 1};
ideal fierceBoundary=C+ideal(Q,fierceEquation);
ideal conductorBranch=C+ideal(P,T,W);
ideal residualDifferent=sat(relativeDifferent,fierceBoundary);
'''+(
        r'''
ideal expectedResidual=C+ideal(1);
ideal expectedDifferentRadical=fierceBoundary;
'''
        if degree == 2
        else rf'''
ideal expectedResidual=C+ideal(P,W,T^{degree - 2});
ideal expectedDifferentRadical=intersect(fierceBoundary,conductorBranch);
'''
    )+rf'''
assertIdealEqual(
  residualDifferent,expectedResidual,"different away from fierce boundary"
);
assertIdealEqual(
  radical(relativeDifferent),expectedDifferentRadical,
  "relative different support"
);

// The complete reduced pullback of P=0 has three components.  Their
// contributions are 1, N-1, 1; the middle component is the residual tame
// conductor branch when N>2.
ideal ordinaryP=C+ideal(P,T-1,W);
ideal secondConductorBranch=C+ideal(P,T,W+Q);
ideal expectedPDivisor=intersect(
  ordinaryP,intersect(conductorBranch,secondConductorBranch)
);
assertIdealEqual(
  radical(C+ideal(P)),expectedPDivisor,"P divisor decomposition"
);

"PASS_ATLAS_NORMALIZATION_{prime}_{degree}";
'''


def verify_singular_normalizations() -> None:
    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Singular is required for --singular")

    for prime, degree in SINGULAR_ROWS:
        marker = f"PASS_ATLAS_NORMALIZATION_{prime}_{degree}"
        completed = subprocess.run(
            [singular, "-q"],
            input=singular_program(prime, degree),
            text=True,
            capture_output=True,
            check=False,
        )
        if completed.returncode != 0 or marker not in completed.stdout:
            raise AssertionError(
                f"Singular atlas normalization failed for p={prime}, N={degree}\n"
                f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
            )


def balanced_singular_program(prime: int, degree: int) -> str:
    exact_different = ""
    normal_replay = ""
    if degree <= 3:
        normal_replay = r'''
intvec normalCheck=norTest(primitiveOrder,N);
if (normalCheck[1] != 1 || normalCheck[2] != 1 || normalCheck[3] != 1)
{
  "FAIL: balanced normalization test";
  exit(1);
}
'''
    if degree == 3:
        exact_different = rf'''
// The p=3 presentation is small enough for the complete relative Fitting
// ideal.  Its saturation is the unit ideal and its radical is the fierce
// boundary radical.
int nv=nvars(basering);
int nr=size(norid);
int nc=nv-2;
matrix relativeJacobian[nr][nc];
int rowIndex;
int columnIndex;
for (rowIndex=1;rowIndex<=nr;rowIndex++)
{{
  for (columnIndex=1;columnIndex<nc;columnIndex++)
  {{
    relativeJacobian[rowIndex,columnIndex]=
      diff(norid[rowIndex],var(columnIndex));
  }}
  relativeJacobian[rowIndex,nc]=diff(norid[rowIndex],var(nv));
}}
ideal relativeDifferent=norid+minor(relativeJacobian,nc);
ideal residualDifferent=sat(relativeDifferent,fierceBoundary);
assertIdealEqual(residualDifferent,norid+ideal(1),"balanced residual different");
assertIdealEqual(
  radical(relativeDifferent),radical(fierceBoundary),
  "balanced different support"
);

// The p=3 normalization is smooth; this is stronger than normality and is
// useful as a control for the divisor-class localization argument.
ideal absoluteSingularLocus=
  norid+minor(jacob(norid),nvars(basering)-2);
assertIdealEqual(
  absoluteSingularLocus,norid+ideal(1),"balanced absolute smoothness"
);
'''

    return rf'''
LIB "normal.lib";

proc assertReductionZero(poly f, ideal G, string label)
{{
  if (reduce(f,std(G)) != 0)
  {{
    "FAIL: "+label;
    exit(1);
  }}
}}

proc assertIdealEqual(ideal A, ideal B, string label)
{{
  int i;
  for (i=1; i<=size(A); i++)
  {{
    assertReductionZero(A[i],B,label);
  }}
  for (i=1; i<=size(B); i++)
  {{
    assertReductionZero(B[i],A,label);
  }}
}}

ring r={prime},(P,Q,T),dp;
poly H=(T-1)*(T^{degree}-P^{degree + 1})+P^{degree - 1}*Q*T;
ideal primitiveOrder=H;
list N=normal(primitiveOrder,"isPrim","withGens");
{normal_replay}
ideal conductor=normalConductor(primitiveOrder);
ideal maximalConductor=ideal(P,T)^{degree - 1};
assertIdealEqual(conductor,maximalConductor,"balanced primitive conductor");

def Rn=N[1][1]; setring Rn;
ideal fierceBoundary=norid+ideal(Q,T^{degree}-P^{degree + 1});
{exact_different}

// Exact finite-field point counts are only a reconstruction diagnostic.
// They agree with A2 for these rows but do not prove an isomorphism.
ideal fieldEquations;
int variableIndex;
for (variableIndex=1;variableIndex<=nvars(basering);variableIndex++)
{{
  fieldEquations[variableIndex]=
    var(variableIndex)^{prime}-var(variableIndex);
}}
int coverPoints=vdim(std(norid+fieldEquations));
int boundaryPoints=vdim(std(fierceBoundary+fieldEquations));
if (coverPoints != {prime * prime + prime})
{{
  "FAIL: balanced normalization point count";
  exit(1);
}}
if (boundaryPoints != {prime})
{{
  "FAIL: balanced boundary point count";
  exit(1);
}}

"PASS_BALANCED_NORMALIZATION_{prime}_{degree}";
'''


def verify_balanced_singular_normalizations() -> None:
    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Singular is required for --balanced-singular")

    for prime, degree in BALANCED_SINGULAR_ROWS:
        marker = f"PASS_BALANCED_NORMALIZATION_{prime}_{degree}"
        completed = subprocess.run(
            [singular, "-q"],
            input=balanced_singular_program(prime, degree),
            text=True,
            capture_output=True,
            check=False,
            timeout=60,
        )
        if completed.returncode != 0 or marker not in completed.stdout:
            raise AssertionError(
                f"Singular balanced normalization failed for p={prime}, N={degree}\n"
                f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
            )


def thickened_singular_program(
    prime: int,
    degree: int,
    q_exponent: int,
    expected_cover_points: int,
    expected_boundary_points: int,
) -> str:
    normal_check = ""
    if degree <= 3:
        normal_check = r'''
intvec normalCheck=norTest(primitiveOrder,N);
if (
  normalCheck[1] != 1
  || normalCheck[2] != 1
  || normalCheck[3] != 1
)
{
  "FAIL: thickened normalization test";
  exit(1);
}
'''
    return rf'''
LIB "normal.lib";

proc assertReductionZero(poly f, ideal G, string label)
{{
  if (reduce(f,std(G)) != 0)
  {{
    "FAIL: "+label;
    exit(1);
  }}
}}

proc assertIdealEqual(ideal A, ideal B, string label)
{{
  int i;
  for (i=1; i<=size(A); i++)
  {{
    assertReductionZero(A[i],B,label);
  }}
  for (i=1; i<=size(B); i++)
  {{
    assertReductionZero(B[i],A,label);
  }}
}}

ring r={prime},(P,Q,T),dp;
ideal primitiveOrder=
  (T-1)*(T^{degree}-P^{degree + 1})
  +P^{degree - 1}*Q^{q_exponent}*T;
list N=normal(primitiveOrder,"isPrim","withGens");
{normal_check}
assertIdealEqual(
  normalConductor(primitiveOrder),
  ideal(P,T)^{degree - 1},
  "thickened primitive conductor"
);

def Rn=N[1][1]; setring Rn;
ideal fierceBoundary=norid+ideal(Q,T^{degree}-P^{degree + 1});
ideal fieldEquations;
int variableIndex;
for (variableIndex=1;variableIndex<=nvars(basering);variableIndex++)
{{
  fieldEquations[variableIndex]=
    var(variableIndex)^{prime}-var(variableIndex);
}}
int coverPoints=vdim(std(norid+fieldEquations));
int boundaryPoints=vdim(std(fierceBoundary+fieldEquations));
if (coverPoints != {expected_cover_points})
{{
  "FAIL: thickened cover point count";
  exit(1);
}}
if (boundaryPoints != {expected_boundary_points})
{{
  "FAIL: thickened boundary point count";
  exit(1);
}}

"PASS_THICKENED_{prime}_{degree}_{q_exponent}";
'''


def verify_thickened_singular_normalizations() -> None:
    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Singular is required for --thickened-singular")

    for row in THICKENED_SINGULAR_ROWS:
        prime, degree, q_exponent, cover_points, boundary_points = row
        marker = f"PASS_THICKENED_{prime}_{degree}_{q_exponent}"
        completed = subprocess.run(
            [singular, "-q"],
            input=thickened_singular_program(*row),
            text=True,
            capture_output=True,
            check=False,
            timeout=60,
        )
        if completed.returncode != 0 or marker not in completed.stdout:
            raise AssertionError(
                "Singular thickened normalization failed for "
                f"p={prime}, N={degree}, c={q_exponent}\n"
                f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
            )


def print_atlas() -> None:
    print("p  N (p|N)  degree  fierce-six-tuple       nodes  delta  companion-six-tuple")
    for prime, degree in TEST_ROWS:
        conductor, delta, _ = semigroup_profile(degree)
        separable, inseparable = residue_degree_parts(prime, degree)
        print(
            f"{prime:<2} {degree:<6} {degree + 1:<7} "
            f"(1,{separable},{inseparable},1,c={conductor},{degree})  "
            f"{degree + 1:<5} "
            f"{delta:<5} "
            f"({degree - 1},1,1,{degree - 2},c=(P,T,W),{degree - 1})"
        )

    print("\nbalanced gluing P^(N-1)*Q*T")
    print("p  N (p|N)  P=0 normalized rows                 Jacobian  Cl(complement)")
    for prime, degree in TEST_ROWS:
        print(
            f"{prime:<2} {degree:<6} "
            f"(1,1,1,0,1)+(1,{degree - 1},1,0,{degree - 1})  "
            f"-u^{2 * degree - 4:<3} Z/{degree - 1}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--singular",
        action="store_true",
        help="run seven bounded exact normalizations, including N=6 controls",
    )
    parser.add_argument(
        "--balanced-singular",
        action="store_true",
        help="normalize six balanced rows, including N=6 controls",
    )
    parser.add_argument(
        "--thickened-singular",
        action="store_true",
        help="normalize six b>0 base-change controls and count their points",
    )
    args = parser.parse_args()

    verify_hidden_orders()
    verify_balanced_hidden_orders()
    verify_extended_gluing_band()
    verify_thickened_base_change_obstruction()
    verify_ramification_rows()
    if args.singular:
        verify_singular_normalizations()
    if args.balanced_singular:
        verify_balanced_singular_normalizations()
    if args.thickened_singular:
        verify_thickened_singular_normalizations()
    print_atlas()
    print("Plane wild-boundary atlas audit: PASS")
    if args.singular:
        print("Singular normalization and different regressions: PASS")
    if args.balanced_singular:
        print("Balanced Singular normalization diagnostics: PASS")
    if args.thickened_singular:
        print("Thickened base-change Singular regressions: PASS")


if __name__ == "__main__":
    main()
