# Cohomological pullback and a uniform original marker

16 November 2026. Internal written proofs with named primary imports.
Independent replay, formal verification and external review are not claimed.
Both COMPLETE exceptional Mathieu--Zhao cases remain **OPEN**; review deferred.

Nonzero intrinsic component classes survive EVERY regular pullback with
nonzero action on ambient H3. In particular, actual invariant phases with
nonzero Neron component at0 exist on BOTH exceptional quotients, for every
parameter in the retained regular Bezout family. Exceptional invariance,
including perfectness of2I, therefore does NOT alone force component torsion.
The eventual-purity hypothesis in the whole-quotient target is essential.

A stronger constructive statement holds for these pullbacks of Long's
phase: ONE fixed regular original marker detects EVERY nonzero polynomial
in that phase, in EVERY degree. Its marked functional is exactly uniform
integration over a real interval. This does not produce an unweighted
pure phase; it makes the missing step precise and removes any marker or
polynomial-recomposition search from that sufficient construction route.
The route is not a replacement for the full exceptional pure-locus target.

See the [checker](verify_haar_pullback.py),
[receipt](haar_pullback_verification.json) and
[source audit](haar_pullback_source_audit.json).

## 1. The reduced intrinsic source is functorial for nonzero H3 maps

Let X=SL2/Gamma and Y=SL2/Delta for ANY finite groups. Let
Phi:X -> Y be a regular algebraic map inducing a NONZERO scalar d on
the primitive rational ambient classes H3=Q(-2). Choose orientations
by the compact cores. Then d is an integer in primitive integral
generators. It is not defined as the generic field degree or a motivic
degree label. Let f be a nonconstant regular phase on Y and g=f o Phi.

Pullback of the complement cohomology, on a common full regular base,
is a morphism of geometric variations

    B_f -> B_g,  u_f |-> d u_g.                         (1)

No properness or etaleness of Phi is needed for ordinary cohomological
pullback. It is functorial for MHS and Gauss--Manin transport. The
canonical lowest-weight reduction is

    W4 B_f=L_f direct_sum C_f, E_f=B_f/C_f, W4 E_f=L_f,

where C_f is pure semisimple with NO globally trivial constituent.
The composite C_f -> B_g -> E_g lands in L_g by weight strictness.
Its underlying flat map must vanish, since Hom(C_f,L_g)=0. Therefore
(1) descends canonically to

    alpha:E_f -> E_g, alpha(u_f)=d u_g.                 (2)

Deleting extra regular values changes neither monodromy image nor this
argument. This is why the canonical compact-middle reduction must be
used before transporting the local component test.

At ANY P1 puncture a, if p_g is a local invariant retraction onto the
source unit, then d^-1*p_g*alpha is one onto the target unit. The same
argument works globally. Consequently, BEFORE any purity assumption,

    xi_boundary,f,a !=0 => xi_boundary,g,a !=0,
    Bad_boundary(f) subset Bad_boundary(g),
    E_f globally unit-nonsplit => E_g globally
       unit-nonsplit.                                  (3)

Infinity is included in the local implication. Additional bad values
can occur on X; no reverse implication or equality of bad sets is claimed.
The hypothesis d!=0 is essential.

The [Neron component criterion](LOCAL_TATE_NORMAL_FORM_AND_NERON_COMPONENTS.md)
translates (3) into

    nu_a(sigma_f)!=0 => nu_a(sigma_g)!=0.               (4)

There is also a direct rational 1-motive description. Dualizing (2)
gives a morphism M_g -> M_f with lattice map d. The retained geometric
1-motive comparison, point full faithfulness and Hom rigidity realize
it algebraically after shrinking and up to isogeny. Its semiabelian
map h satisfies h(sigma_g)=[d]sigma_f in that rational category.
After clearing denominators there is an honest homomorphism and a
nonzero integer multiple on the lattice. Neron functoriality then
gives the same relation on rational components. Isogenies give
isomorphisms on rational component groups because a dual isogeny
composes to integer multiplication. No arbitrary relative essential-
surjectivity or Neron-model producer is required.

If E_f is globally unit-nonsplit, the exact all-degree source reflection
under polynomial recomposition makes E_H(g) globally unit-nonsplit for
EVERY nonconstant H. Thus EVERY eventually pure such H(g), IF one exists,
is a full actual failure by the intrinsic safety theorem. This uses
purity only for the safety conclusion, not for source transport.

## 2. Top-period pullback supplies an exact weighted functional

Normalize nowhere-vanishing regular volume forms Omega_X,Omega_Y to
have integral1 on their compact fundamental cycles, so their periods
are the normalized Haar functionals h_X,h_Y. Write

    Phi^*Omega_Y=J_Phi Omega_X, a=J_Phi/d.              (5)

J_Phi is regular. Pullback of any closed top form R Omega_Y and (1)
give, for EVERY regular R on Y,

    h_X(a*(R o Phi))=h_Y(R), h_X(a)=1.                 (6)

This is a topological period identity, not change of variables for a
finite covering. Phi may be ramified and nonproper. It is dominant:
otherwise its differential has rank at most2 everywhere, so J_Phi=0,
contradicting its nonzero H3 action. In particular g is nonconstant.

The weight a in (6) is an actual fixed original regular function. It
cannot be dropped from the UNWEIGHTED pure moments h_X(H(g)^m).
The unweighted moment functional is still the one in the MZ question.

## 3. ONE original marker gives the interval functional for ALL polynomials

Now take Y=SL2 with xw-yz=1 and Long's pinned phase and marker

    F=(1+y)(1+z+yz), Q=-y,
    h_Y(F^n)=0,
    h_Y(Q F^n)=(-1)^(n-1)/(n+1), n>=1,
    h_Y(Q)=0.                                         (7)

These are retained actual all-order identities, not abstract formal
moments. For EVERY polynomial R in C[s], linearity gives

    h_Y(R(F))=R(0),
    h_Y(Q R(F))=R(0)-integral_0^1 R(-t)dt,
    h_Y((1-Q)R(F))=integral_0^1 R(-t)dt.                (8)

For ANY regular Phi:X -> SL2 with d!=0, put

    g=F o Phi,
    C_Phi=a*((1-Q) o Phi)=a*(1+y o Phi).                (9)

This is ONE fixed original regular marker, independent of R,H and m,
with h_X(C_Phi)=1. Equations(6),(8) give the exact identity

    h_X(C_Phi R(g))=integral_0^1 R(-t)dt                (10)

on the ENTIRE polynomial algebra C[g]. In particular for EVERY nonzero
polynomial H of ANY degree,

    h_X(C_Phi H(g)^m)=integral_0^1 H(-t)^m dt.          (11)

The published theorem of Muger--Tuset, Indag.Math.31(2020),147--151,
Theorem2.1, asserts that the right side of (11) has

    limsup_(m->infinity) |integral_0^1 H(-t)^m dt|^(1/m)>0
       for EVERY nonzero complex polynomial H.         (12)

The primary preprint is [arXiv:1910.12148v1](https://arxiv.org/pdf/1910.12148v1).
The complete proof on PDF2--4 was read, including repeated roots and
endpoint zeros. This is not a positivity argument for complex H or a
finite-moment inference. Nonzero constants are allowed and elementary.

For clarity, the part of that proof needed to rule out eventual zero
has a short explicit form. For a nonconstant polynomial p, form

    A(z)=integral_0^1 dt/(1-z p(t)).                    (13)

Away from the reciprocal critical and endpoint values, partial fractions
express A as -1/z times a sum of logarithms divided by p'(r_i), where
p(r_i)=1/z. Continue along a path to z=infinity. Near a zero of p of
multiplicity e, the corresponding roots separate on the scale z^-1/e,
so 1/p'(r_i)=O(z^((e-1)/e)). The logarithm is bounded, except for an
endpoint zero where it is O(log z). Every term in (13) therefore tends
to0 as O(z^-1/e log z). If all sufficiently high moments vanished,
A would be a polynomial with A(0)=1; its continuation could not tend
to0. This proves non-eventual vanishing. Muger--Tuset additionally use
the same estimate and Liouville to obtain (12).

Therefore EVERY nonzero H(g) fails all-original-marker safety, with the
SAME fixed marker C_Phi. No eventual purity is assumed for this statement.
If ANY nonconstant H(g) is eventually pure for the UNWEIGHTED original
Haar functional, it is an actual MZ counterexample, with (9) already
an original regular invariant bad marker. The primitive g need not be
pure. There is no marker search, no dependence of C_Phi on deg H and
no assumption that only a zero-constant H is considered.

Equation(12) gives infinitely many nonzero marked orders, not necessarily
EVERY order. For example H(s)=s(s+1)(s+1/2) makes H(-t) antisymmetric
under t -> 1-t, so all odd moments vanish while all even moments are
strictly positive. This exact control prevents an all-orders overclaim.

## 4. Actual exceptional components and the integral H3 degree

For BOTH Gamma=2O,2I use EVERY retained
[regular Bezout map](REGULAR_BEZOUT_MAPS_AND_WEIGHTED_HAAR_TAILS.md)

    Phi_b=[[U+F_0 b,V+G_0 b],[F_0,G_0]],
    b in A_Gamma=C[SL2]^Gamma.                         (14)

Here F_0=f(z,w),G_0=g(z,w) are the coprime Klein forms of degrees
(n,m)=(8,12) or(12,20); these f,g are binary forms, not the phase F.
The determinant identity holds globally for EVERY regular b.

The bottom-row torsor maps give a commutative diagram over

    k:(C2 minus0)/Gamma -> C2 minus0, [z,w] |->(f,g).

Each additive Ga-torsor projection is a homotopy equivalence. The finite
map (f,g):C2 -> C2 has degree nm, and its quotient map has degree |Gamma|.
Thus k has degree nm/|Gamma|=2. This is also its degree on oriented
link H3, by the local complex topological degree. The diagram proves

    d_(Phi_b)=2 for EVERY b.                           (15)

This computes the actual primitive integral H3 action; it does NOT
identify the generic field degree of Phi_b, whose fibre polynomial
can depend on b. In the previous normalized weighted-transfer notation,

    W_b=c(z,w)*(H_0+D b), D=z partial_x+w partial_y,
    J_(Phi_b)=2 W_b, a=W_b.                            (16)

The factor2 follows from (15) and normalization of both compact periods.
The earlier raw algebraic volume formula had the factor nm and a source
core volume1/|Gamma|; it is consistent with (16), not an alternative
degree convention. W_b always vanishes on the branch divisor c(z,w)=0
and h_Gamma(W_b)=1, so it cannot be replaced by the constant1.

Write Y_b=V+G_0 b. Equations(9)--(11) specialize to the explicit
original invariant functions

    g_b=(1+Y_b)*(1+F_0+Y_b F_0),
    C_b=W_b*(1+Y_b),
    h_Gamma(C_b R(g_b))=integral_0^1 R(-t)dt.           (17)

This holds for EVERY b and EVERY polynomial R. In particular the SAME
C_b detects every nonzero H(g_b) in every degree by (12).

The primitive target F in (7) has canonical reduced rank2 source with
Kummer representative s/(s+1). Its reduced intrinsic classes at0 and-1
are nonzero. The pullback theorem gives, for EVERY b on BOTH quotients,

    nu_0(sigma_g_b)!=0, nu_(-1)(sigma_g_b)!=0,
    rank N_E(g_b) - rank N_P(g_b)=1 at BOTH values.     (18)

These are ACTUAL invariant algebraic phases and their genuine geometric
sources, not arbitrarily chosen local representations or section models.
Their primitive sources and those of EVERY H(g_b) are globally unit-
nonsplit. Consequently exceptional invariance, finite fundamental group,
or perfectness of2I cannot by themselves force the local component zero.
The family keeps a nonzero local obstruction even at the minimum
possible rank excess1; requiring bounded graph rank would not fix it.

Nothing here proves that any g_b or H(g_b) is unweighted eventually pure.
The earlier Hopf pullback P_b is a DIFFERENT target phase from g_b; its
canonical first-moment exclusions are not asserted for (17). No b or H
search is proposed or performed. The fixed uniform marker closes the
marker/recomposition side, not that missing purity condition.

## 5. The whole-quotient target still needs exceptional purity

The FULL remaining implication is

    f in A_Gamma, h_Gamma(f^m)=0 eventually
       => nu_0(sigma_f)=0, Gamma=2O or2I.              (19)

The new actual controls show why the premise in (19) must enter a proof.
Purely geometric exceptional invariance has genuine nonzero components.
For the sufficient family (17), any eventual-pure nonzero polynomial
would refute (19) and full MZ using ONE marker C_b. Conversely, if the
full MZ statement holds, the unweighted eventual-pure radical of C[g_b]
must be zero for EVERY b. Neither existence nor that entire-family
exclusion is proved here. Proving exclusion only for this family would
still not prove (19) for the whole quotient.

The next whole-quotient milestone remains a proof using BOTH exceptional
invariance and the actual unweighted purity identities to force torsion
of the Picard/NS/gluing section component, or an actual invariant pure
phase violating it. No hidden actual-kernel/genus-semiplicity branch,
scalar norm, low-degree campaign, coefficient-ring promotion or review
gate is introduced. Both COMPLETE exceptional MZ cases remain OPEN.
