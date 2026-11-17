# Tame original-period rigidity and the pure boundary defect

16 November 2026. Internal written deduction from precisely identified
classical theorems. Independent replay, formal verification and external
review are not claimed; review remains deferred. Both COMPLETE exceptional
Mathieu--Zhao cases remain **OPEN**.

The result changes the remaining global route: **a nonconstant eventually
pure phase with a nonzero normalized compact volume period cannot be
cohomologically tame**. More generally, every original marked period of
a tame phase is either zero or has a nonrational ordinary moment resolvent.
Thus tameness cannot provide the missing hypothesis for applying the earlier
semisimplicity criterion to a nonzero exceptional pure phase.

This is a theorem about every phase and every original form under the stated
tameness condition, in every degree. It does not assert that all wild phases
are safe, nor that every wild phase has the required unit defect. The
[intrinsic Haar-source theorem](INTRINSIC_HAAR_STATES_AND_NONSPLIT_PERIOD_EXTENSIONS.md)
and the [separated-core criterion](SEPARATED_HAAR_CORES_AND_GLOBAL_KUMMER_EXTENSIONS.md)
retain the full all-phase/all-original-marker finish line. See the
[exact controls](verify_tame_original_period.py),
[receipt](tame_original_period_verification.json) and
[source audit](tame_original_period_source_audit.json).

## 1. The original-form lattice and its precise classical input

Let X be a smooth complex affine variety of dimension n>=2 and P:X->A1
nonconstant. Write s for the target coordinate, u=partial_s, and use the
FULL ordinary top Gauss--Manin module

    G_P = Omega^n(X)[u]/(d-u*dP wedge)Omega^(n-1)(X)[u],
    s[omega]=[P omega] for u-independent omega,
    s=P-partial_u on general representatives.              (1)

Let M_P be the image of the ORIGINAL top forms Omega^n(X) in G_P.
It is a C[s]-submodule and generates G_P over D_s. Neither arbitrary
u-dependent representatives nor a localized Brieskorn lattice are substituted
for M_P.

Use the decreasing Kashiwara--Malgrange filtration at s=c, normalized so
that (s-c)partial_s acts with eigenvalue beta on gr_V^beta. In this
convention a point module delta_c has basis partial_s^j delta_c and

    (s-c)partial_s(partial_s^j delta_c)
       =-(j+1)partial_s^j delta_c,
    V_c^(>-1)(delta_c)=0.                                (2)

Two standard facts are the only filtration inputs:

1. D-module morphisms are strict for this filtration, including inclusions
   of submodules. Thus if L is point-supported at c in a holonomic B,
   then L intersect V_c^(>-1)B=0.
2. If P is cohomologically tame, the original-form module satisfies

       (M_P)_c subset V_c^(>-1)(G_P)_c                  (3)

   at every finite c. M_P is also finite free over C[s].

For (3), use [Sabbah, Hypergeometric periods for a tame polynomial](https://perso.pages.math.cnrs.fr/users/claude.sabbah/articles/sabbah_bbases.pdf),
Section 8 for the hypothesis, Theorem 10.1 and the full freeness argument
with Lemma 10.7 on printed pages 29--30. His increasing V_(<1) is (3)
in the decreasing convention. The source treats arbitrary smooth affine X,
not only affine space. The 2008
[Nice notes](https://perso.pages.math.cnrs.fr/users/claude.sabbah/livres/sabbah_nice08.pdf),
Proposition 3.2.3, printed page 15, state the same lattice containment and
also give it for M-tame functions. We use that M-tame extension only with
the hypothesis as defined there. No assertion that the two tame notions
are equivalent is needed.

For strictness and (2), see
[Sabbah's Hodge-module chapter 6](https://perso.pages.math.cnrs.fr/users/claude.sabbah/MHMProjectV1/mhm_chap6.pdf),
Proposition 6.2.7, Exercises 6.2.8, 6.2.11 and 6.2.12, printed pages
101--102. Strictness follows by inducing the coherent filtration on an
image and applying uniqueness; (2) follows directly from the displayed
point-module basis. These are classical inputs, not machine-verified
claims or a use of localized semisimplicity.

Here cohomological tameness means that a proper partial compactification
has vanishing cycles supported at finitely many points, all inside X,
as in Sabbah Section 8. We do not extend this import to arbitrary
nonisolated singularities or replace it by weak tameness.

## 2. Original-period rigidity, with all markers

Let H_P be the exponential top source, with differential d+t*dP wedge.
Fourier identification is

    u -> -t,       s -> partial_t,       G_P -> H_P.       (4)

Take ANY D_t-linear formal state lambda:H_P->C[[t]]. In particular, a
compact real n-cycle gamma and an original top form omega give the state
integrating exp(tP)omega. No condition dim H^n(X)=1 is needed here.
Put V=image(lambda). It is a holonomic quotient of H_P; its inverse
Fourier transform B is a regular holonomic quotient of G_P. For every
original omega, let F_omega=lambda([omega]).

**Theorem.** If P is cohomologically tame (or M-tame with (3)), then
F_omega is either zero or is NOT a finite exponential-polynomial sum.
Equivalently, the original moment series

    C_omega(s)=sum_(m>=0) a_m s^(-m-1),
    F_omega(t)=sum_(m>=0) a_m t^m/m!                    (5)

is either zero or nonrational. In particular a nonzero original period
cannot be a polynomial, even if some other elements of V are polynomials.

Proof. Suppose F_omega=sum_c exp(ct)p_c(t) is nonzero, with finitely
many polynomial p_c. Polynomial Chinese-remainder projectors in partial_t
isolate any nonzero summand F_c. They act in (1) as polynomials in s,
so F_c is still the image of an ORIGINAL top form q_c(P)omega in M_P.

The D_t-submodule generated by exp(ct)p_c(t), p_c!=0, is
C[t]exp(ct). Under (4) it is exactly a copy of delta_c in B. By (3)
and functoriality, the image of q_c(P)omega belongs to V_c^(>-1)B.
By strictness of the inclusion delta_c->B and (2), the intersection is
zero. This contradicts F_c!=0. No splitting of B or semisimplicity is
assumed. In particular the argument retains the parameter-zero lattice.

Partial fractions prove the equivalence in (5): a rational ordinary
resolvent gives a sequence satisfying a constant-coefficient recurrence
after finitely many indices, and its exponential generating function is
a finite sum exp(ct)p_c(t), including the finite initial sequence as the
c=0 polynomial term. Conversely these functions give rational ordinary
resolvents. This is a formal identity and requires no convergence.

As a source-level consequence, Hom_Ds(G_P,delta_c)=0 for every c:
a morphism kills M_P by (2)--(3), and M_P generates G_P. This does NOT
say that G_P or an actual quotient B has no point-supported submodule.

## 3. The compulsory origin defect of EVERY eventually pure phase

Return to a normalized compact volume state with lambda([Omega])(0)=1,
including both exceptional quotients. If P is eventually pure, then

    p(t)=lambda([Omega]) in C[t], p(0)=1.                (6)

Therefore P cannot be cohomologically tame or M-tame in the preceding
sense. This applies to pure SAFE phases and pure FAILURES alike. It is
not an additional condition distinguishing the two.

There is an unconditional source obstruction, without a tame hypothesis.
For EVERY eventually pure phase, the ORIGINAL volume class satisfies

    [Omega]_0 notin V_0^(>-1)(G_P)_0.                    (7)

More precisely, if deg p=N, then

    [Omega]_0 notin V_0^(>-N-1)(G_P)_0.                  (8)

Indeed its nonzero image under G_P->B lies in delta_0. In the basis
(2), its highest nonzero derivative order is N, whose filtration index
is -N-1. Strictness of the inclusion and functoriality of the quotient
prove (7)--(8). Thus the original unit has a nonzero class in the exact
local defect

    (G_P)_0 / V_0^(>-1)(G_P)_0.                          (9)

The tame theorem puts every original form on the other side of this
boundary. A general wild phase may have (9) nonzero without being pure;
its existence does not prove purity or failure. No general equality of
(9) with boundary vanishing cycles is asserted here. In particular we
do not silently infer a localized-in-value tameness theorem, or infer
cohomological tameness from Brieskorn finiteness without its hypotheses.

For either exceptional quotient, retain K_P from the separated-core
theorem. The FULL remaining implication is still

    [P Omega] in K_P  => H_P=C[t][Omega]+K_P             (10)

for EVERY actual phase and ALL original fixed markers. It must now be
proved on the non-tame locus with the compulsory origin defect (7).
The earlier unlocalized-semisimplicity sufficient theorem remains valid;
obtaining its hypothesis merely from tame-function theorems cannot reach
a nonconstant pure phase. Neither (7) nor wildness closes (10).

## 4. A genuine tame geometric control: preserve the ORIGINAL form

The distinction between an original form and a general element of the
Fourier module is essential even in an actual volume-period example.
Let

    X={uv+z^2=1},       P=z-1,
    Omega=du wedge dz/u                               (11)

with the residue extension of Omega across u=0. This is a smooth affine
quadric surface, and Omega is nowhere zero. The compact real sphere has
u=x+iy, v=x-iy, x^2+y^2+z^2=1. Its normalized volume makes z uniform on
[-1,1], up to the irrelevant common complex normalization of Omega.

The phase is cohomologically tame. An explicit proper partial
compactification is

    U_1 V_1=(1-z^2) U_0 V_0
             in P1_U x P1_V x A1_z.                    (12)

The total surface is smooth. The boundary consists of the two disjoint
sections (u=infinity,v=0) and (u=0,v=infinity). Near the first section
the equation is V_1=(1-z^2)U_0, so the map and its boundary are a product
in coordinates (U_0,z); likewise at the other section. Hence Rj_*C has
no boundary vanishing cycles. The only affine critical points have
u=v=0,z=+-1 and are isolated Morse points. This verifies the precise
cohomological-tameness hypothesis, not only a generic smoothness test.

Its ORIGINAL unit period is

    F(t)=h(exp(tP))=(1-exp(-2t))/(2t).                   (13)

It is entire and nonpolynomial, as the theorem requires. Nevertheless
the actual period module contains 1, because

    (2t+1)F+tF'=1.                                    (14)

In (14), F' is the original marked period of P, but the coefficients
involve t. Thus 1 is represented by the t-dependent amplitude
(2t+1)+tP, not by the original volume Omega.

For completeness, the full actual marked module here is D_t F.
Angular averaging on the compact sphere turns any original polynomial
Q into a polynomial in z: in the coordinate ring expand Q into powers
of u or v with coefficients in z, using uv=1-z^2, and only equal
angular weights survive. Its period is therefore a polynomial in
partial_t applied to F. Conversely each derivative is an original
marker P^j. In particular V/tV=C by the intrinsic source theorem.

The submodule C[t] supplied by (14) does NOT split off. Its quotient
is nonzero because F is not polynomial. Since V consists of entire
functions, V/C[t] is torsion-free over C[t]. It is t-surjective by
V/tV=C and 1 in V; it is therefore t-invertible. A section would embed
a nonzero t-surjective submodule into C[[t]], which is impossible.
This is a genuine tame geometric nonsplit polynomial-period extension,
but its ORIGINAL unit is not pure. It is not either exceptional quotient,
and it does not realize the earlier nonintegral Kummer family there.

## 5. Evidence and remaining task

The checker verifies the point-module Euler indices, a nontrivial
Chinese-remainder projector, the quadric boundary charts, the compact
height moment formula and identity (14). The universal V-filtration and
tame-lattice arguments remain written proofs with pinned primary inputs.
No degree, support, prefix or singularity enumeration is performed.

The next task is the full actual-source implication (10) on the non-tame
exceptional geometry, with the original unit defect retained. Both full
MZ cases remain open. External review is deferred, not a next step.
