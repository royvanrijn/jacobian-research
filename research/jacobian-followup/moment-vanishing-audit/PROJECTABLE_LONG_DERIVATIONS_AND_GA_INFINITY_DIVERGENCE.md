# All projectable Long purity derivations have non-volume-preserving exceptional lifts

16 November 2026. Internal written proof; independent replay, formal
verification and external review are not claimed. External review remains
deferred. Both FULL exceptional Mathieu--Zhao cases remain **OPEN**.
The completed family-transfer theorem and all earlier goals/results remain
in force.

This strengthens the [previous no-lift theorem](LONG_DERIVATION_LIFT_AND_CUSP_BRANCH_OBSTRUCTION.md).
For EVERY original Bezout parameter on BOTH exceptional quotients, EVERY
regular divergence-free target derivation with delta(g)=g that projects to
the (F,Y)-plane has a rational function-field lift with NONZERO divergence.
Regularity of the lift is not needed for this obstruction. This includes
arbitrary polynomial vertical corrections of every degree, not just Long's
original Hamiltonian field. The positive invariant-degree gap rules out
the only possible leading cancellation.

For Long's original field, more is true: the lift's divergence is not even
a rational function of the phase over the bottom invariant field. Thus no
nonzero rational rescaling depending only on the phase produces a scalar
first-order Stokes identity. This closes those transport mechanisms. It
neither excludes nonprojectable target fields or independent original
certificates, nor proves nonpurity of all Bezout phases, full pure-to-safe
inclusion, or an actual pure counterexample.

## 1. Setup and the original Ga coordinate

Let A=C[SL2]^Gamma, Gamma=2O or2I. Use the retained full
[Klein/Ga-torsor model](QUADRATIC_TRACE_MODELS_AND_HAAR_RESIDUES.md), with
bottom invariants F,G,C, respective binary degrees

    (n,m,r)=(8,12,18) or (12,20,30), r=n+m-2,
    B=ker(D:A->A)=C[F,G,C]/(C^2-K(F,G)),
    D=z partial_x+w partial_y,
    K_O=G(F^3-108G^2)/36,
    K_I=12F^5+G^3/144.

Put K0=Frac(B). The rational Ga-torsor is trivial: there is a rational
coordinate v with Dv=1 and Frac(A)=K0(v). Each original regular function
is a polynomial in v over K0, since D is locally nilpotent. If P has
positive v-degree d and leading coefficient ell, then

    ell=D^d(P)/d! belongs to B.                         (1)

In particular ell is REGULAR in the entire positively graded binary
invariant ring, not an arbitrary rational function of F,G,C. This is
the key restriction at infinity in the Ga coordinate.

The original volume has the exact shape, up to a nonzero constant,

    Omega_Gamma=dF wedge dG wedge dv / C.              (2)

For example on an upstairs chart with v0=x/z the SL2 form is a constant
multiple of dv0 wedge dz wedge dw. Any rational slice v differs from
v0 by a bottom function, which does not change this three-form.
The binary identity dF wedge dG=nm C dz wedge dw then gives (2).

For EVERY regular original b, use the retained
[regular Bezout map](REGULAR_BEZOUT_MAPS_AND_WEIGHTED_HAAR_TAILS.md)

    Phi_b=[[X_b,Y_b],[F,G]],
    X_b=U+Fb, Y_b=V+Gb, UG-VF=1.

Write P=Y_b as a polynomial in v. Its degree d is positive. Otherwise
D(Y_b)=0 and the determinant also gives D(X_b)=0, so X_b,Y_b belong to
B, contradicting 1=X_bG-Y_bF in the proper positive-degree ideal (F,G)
of B. We do not need the stronger retained depth bound19/31.

Let ell!=0 be the leading coefficient of P; (1) applies. In the target
coordinates [[X,Y],[F,G]], normalized volume is a constant multiple of

    Omega_*=dX/X wedge dY wedge dF
           =-dF wedge dG wedge dY/G.

Therefore the nonzero Jacobian is

    Phi_b^*Omega_*=J Omega_Gamma,
    J=c C P_v/G, c in C*, P_v=partial_v P.            (3)

This is the earlier J=c' C(H+Db). It is nonzero for every b. All
function-field extensions below are finite separable, so a target
derivation has a UNIQUE rational lift. No regularity of that lift is
silently assumed.

## 2. Exactly which target fields are projectable

Put t=1+Y and

    g=t+F t^2,
    delta X=0, delta Y=-Yt,
    delta F=1+F+2YF, delta G=YG,
    T=X partial_X-G partial_G.

Both delta and T are regular and divergence-free. They satisfy

    delta(g)=g, T(g)=0.

A target derivation is PROJECTABLE here precisely when its values on F
and Y belong to C[F,Y]. This definition does not restrict their degrees.

**Classification.** Every globally regular divergence-free target field
E_* with E_*(g)=g that is projectable has the form

    E_*=delta+S(F,Y)T, S in C[F,Y].                    (4)

Here is a direct proof. Write its induced plane field as V, and set
N=F V(Y)+Y V(F). Differentiating the determinant at X=G=0 proves N
vanishes on 1+YF=0, so N=(1+YF)N0 for a polynomial N0. The plane field
extends regularly with X-value0 and G-value GN0. Its divergence on the
X-chart is div_plane(V). The difference from E_* is S0 T for a regular
S0: modulo the prime ideal (X), G E_*(X)=0 forces E_*(X) divisible by X,
and the determinant then determines the opposite G-component.

Hence 0=div(E_*)=div_plane(V)+T(S0). The first term has T-weight0,
whereas T(S0) has zero weight0 part in the integer weight decomposition
of the coordinate ring. Both vanish. Thus V is a polynomial Hamiltonian
plane field, and S0 is T-invariant, hence belongs to C[F,Y].

The plane projection of delta has Hamiltonian

    H0=Y(1+F+FY), V0(Y)=-partial_F H0,
    V0(F)=partial_Y H0.

The difference of any two Hamiltonians giving V(g)=g commutes with g.
In C(g)(t), the Hamiltonian field of g acts on t as -t^2 partial_t;
its constants are exactly C(g). A rational function of g regular in
C[F,Y] is a polynomial in g: any finite pole has a nonempty fibre and
cannot disappear. The possible plane fields are consequently

    V=V0+R'(g)Z_g, Z_g(Y)=-partial_F g,
    Z_g(F)=partial_Y g, R in C[g].

But on 1+YF=0,

    g=-1-1/Y, F Z_g(Y)+Y Z_g(F)=1/Y.

Tangency of V to this hyperbola forces R'(-1-1/Y)=0 identically, so
R'=0. This proves (4). It is a classification of the specified target
fields, not a classification of all pure phases or all original fields.

## 3. Two leading terms that control every polynomial correction

Let L and M be the rational lifts of delta and T to K0(v). On K0 define
derivations A0,A1,B1 by

    A0(F)=1+F, A0(G)=0,
    A1(F)=2F, A1(G)=G,
    B1(F)=0, B1(G)=G.

Each extends uniquely to C using C^2=K. These derivations can have poles
at C=0; they are not asserted regular on the Klein surface. Since Y=P,

    L|K0=A0+P A1,
    L(v)=[-P(1+P)-A0(P)-P A1(P)]/P_v,
    M|K0=-B1, M(v)=B1(P)/P_v.                         (5)

Logarithmic derivative notation below means E(log a)=E(a)/a in the
function field, with no analytic logarithm chosen. At v=infinity,
P=ell v^d+O(v^(d-1)). Equations (3),(5) give

    L(log J)=u_d v^d+O(v^(d-1)),
    u_d=(ell/d)[A1(log(ell C^d))-(2d-1)],             (6)

and

    M(log J)=w0+O(v^-1),
    w0=-(1/d) B1(log(ell C^d/G^d)).                  (7)

For clarity, the coefficient in (6) comes from the three terms

    ell A1(log(C/G)), A1(ell),
    -(d-1)(ell+A1(ell))/d.

Their sum is exactly (6), because A1(G)/G=1. Formula (7) follows in the
same way from M(v)'s leading linear term. These are rational-function
expansions, not moment asymptotics or a finite degree calculation.

## 4. The branch valuation makes the leading coefficients nonzero

Use the SAME genuine cusp component as in the previous proof:

    G^q=lambda F^p, (p,q)=(3,2) on2O or(5,3) on2I,
    C=0, F,G units at its generic point.              (8)

The base is smooth there, with local coordinates F,C: K_G!=0.
If a base derivation A has A(F) regular and A(C) with a genuine simple
pole, then for any function h with positive valuation e on this divisor,

    A(log h) has a genuine pole of order2.            (9)

Indeed write h=C^e u, u a local unit. The term e A(C)/C has order-2;
A(u)/u has order at worst-1 in these coordinates, so cannot cancel it.
This applies because ell is REGULAR in B and ell C^d/G^d has valuation
at least d>0. A rational leading coefficient with a pole would not give
the same conclusion; (1) is essential.

Both B1 and A1 fail to be tangent to the cusp: their weights on F,G
are (0,1) and (2,1), whereas cusp tangency requires weights proportional
to (q,p). Thus their C-components have genuine simple poles. Equations
(6)-(9) prove

    w0!=0 and u_d!=0.                                (10)

In particular L(log J) has EXACT pole order d at v=infinity.

## 5. All polynomial corrections, including the exceptional slope

The rational lift of (4) is L+S(F,P)M. Since the target field has zero
divergence, Lie differentiation of (3) gives

    div_Omega_Gamma(L+S(F,P)M)
       =-(L+S(F,P)M)(log J).                          (11)

If S has Y-degree k>=2 with nonzero leading coefficient s_k(F), (6),(7)
give a term s_k(F)ell^k w0 v^(kd). It cannot cancel the pole of order d
from L(log J). Thus (11) is NONZERO.

It remains to allow S=s1(F)Y+s0(F), including constant and zero S. Put

    A=A1-s1(F)B1,
    A(F)=2F, A(G)=(1-s1(F))G.

The coefficient of v^d in (L+S(F,P)M)(log J) is

    (ell/d)[A(log(ell C^d))
                -d(1-s1(F))-(d-1)].                 (12)

If A is not tangent to the cusp, (9) makes this nonzero. Tangency can
hold identically only when

    1-s1(F)=2p/q=2m/n.                               (13)

At this exceptional slope, A=(2/n)E_bin is the POSITIVE binary-degree
Euler derivation on B. In particular A(C)/C=2r/n. Vanishing of (12)
would now require

    A(ell)=-[d(1-4/n)+1]ell.                         (14)

For n=8 and n=12 this eigenvalue is strictly NEGATIVE. Every nonzero
homogeneous component of the regular binary invariant ell has
NONNEGATIVE Euler eigenvalue, including a possible constant component.
Equation (14) is impossible. No degree bound on ell or d was used.

This proves (11) nonzero for EVERY polynomial S, EVERY regular b and
EVERY degree, on BOTH quotients. Even the unique RATIONAL lift fails
volume preservation. A fortiori it cannot be a regular divergence-free
original certificate. The degree gap is used in (14), not inferred
from any collection of tested degrees.

## 6. Scalar rescaling of Long's original field cannot close Stokes

Retain S=0. The phase g_b=g(Phi_b) has v-degree2d with leading term
F ell^2 v^(2d). Every nonzero element of K0(g_b) therefore has valuation
at v=infinity divisible by2d. But L(log J) has pole order EXACTLY d,
by (6),(10). Consequently

    div L=-L(log J) does NOT belong to K0(g_b).        (15)

For every nonzero rational function R in C(s), the rational lift of
R(g)delta is R(g_b)L, and its divergence is

    R(g_b)div L+g_b R'(g_b).                         (16)

This still does not belong to K0(g_b). Thus multiplying the target
field by a phase-dependent rational scalar cannot produce an original
first-order Stokes operator with both its phase derivative and divergence
functions only of g_b. It does not repair the missing unweighted purity
identity. No claim like (15) is made for arbitrary polynomial corrections
S of Section5, where the leading order can be a different multiple of d.

As a boundary control, on the trivial quotient with the identity map and
F=z,G=w,C=1, one may take d=1 and ell=G. The leading coefficient in (6)
is zero and J is constant. There is no cusp divisor (8). The argument
therefore does not exclude Long's actual SL2 counterexample.

## 7. Scope relative to the full obligation

All these conclusions survive precomposition by ANY algebraic source
automorphism, since it scales the original volume by a nonzero constant
and conjugates rational derivations without changing zero divergence.
The earlier [intrinsic automorphism theorem](INTRINSIC_REYNOLDS_AND_FINITE_QUOTIENT_DESCENT.md)
is retained; no completeness or exhaustion assertion is added.

The fixed original marker

    C_b=W_b(1+Y_b),
    h(C_b R(g_b))=integral_0^1 R(-u)du

from the retained [uniform-marker theorem](COHOMOLOGICAL_PULLBACK_AND_A_UNIFORM_ORIGINAL_MARKER.md)
still detects every nonzero polynomial in g_b infinitely often. It is
unchanged. A different, independent original derivation E with
E(g_b)=g_b could still prove purity and yield that actual counterexample.
It need not descend through Phi_b. A nonprojectable target field is also
outside this theorem. Neither existence nor exclusion of those fields is
claimed, and eigenfunction certificates are not necessary for purity:
the retained primitive PURE SAFE F8+G12 and F12+G20 controls have nonzero
critical values and cannot satisfy E(f)=f for any regular derivation.

The FULL target remains every pure original phase safe in ALL degrees,
equivalently intrinsic boundary extension splitting, or ONE actual pure
failure with a fixed original bad marker. Both full cases remain OPEN.
Retain all capacity and automorphism constraints, including failure
support at least6 on2O/8 on2I after every allowed automorphism/torus,
and all earlier results/goals. No bounded-degree campaign, additional
decidability result, new arbitrary-ring promotion or review gate is used.
