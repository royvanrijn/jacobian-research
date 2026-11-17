# Full mixed cubic exclusion by the rescaled period at infinity

Current quantitative refinement, 15 November 2026: the [joint-degree actual-seed proof](JOINT_DEGREE_GLOBAL_353.md) gives global pure353, odd352, full W-linear336, rank-one140, mixed-parent255 and cusp255/tangent284. Current least-prefix intervals are global[20,353], odd[20,352] and all-homogeneous[14,336]. Earlier590/888 bounds and rank-one425/mixed-parent320 estimates below are historical valid coarser evidence. Leastness,20/14 sufficiency, expanded I353 cofactors and independent review remain open.

14 November 2026. Internal written proof, with exact algebraic checks and
explicit imported period-comparison theorem. No external review, formalization,
new priority claim, or optimality claim for the finite gate is asserted.

## Result and scope

**Theorem.** For every complex a and every G(Z,T) of total degree at most
three, the polynomial

    P=W*(a+T+Z^2)+G(Z,T),       E(ZW)=E(T^2)=1

has a nonzero positive pure Gaussian moment. Consequently the entire
fourteen-coefficient cubic architecture

    W*(a+bZ+dZ^2+cT)+G(Z,T),    deg G<=3,

satisfies GMC: its c*d!=0 interior has no pure-null members, and its c=0
and d=0 faces are already closed by the preceding cubic-height and
linear-W/T theorems. This does not normalize arbitrary ternary cubics:
WZT, WT^2 and higher W powers remain outside the architecture.

This closes the full family left open in the [previous mixed-chart note](MIXED_WT_WZ2_CHART.md).
The earlier sixty-second standard-basis run remains a terminal **UNKNOWN
computation**; it is not promoted to a certificate. The new proof uses a
different argument. The optimal eight-moment spine theorem and its
seven-moment witness remain useful finite-prefix results. A new sufficient
320-moment full-family gate is proved below; no twelve-moment gate or
optimality is inferred.

## Source reconciliation

A newly inspected source, [Fixed-Seed Cyclicity Loci for Polynomial
Exponential Periods, candidate 0.2.1](https://github.com/ipitchford/cyclicity-loci-exponential-periods/blob/517b980c5e2944b5aa2a46f7325e65c199b9a985/MANUSCRIPT.md),
organizes polynomial twisted reduction and warns that a seed need not
be cyclic and that a constrained family requires its actual tangent
connection. Those are useful methodological inputs, not a theorem
excluding our Gaussian chart. Its quartic/quintic tables are not imported
or claimed independently reproduced. The inspected commit is
517b980c5e2944b5aa2a46f7325e65c199b9a985, dated 8 August 2026; the manuscript
is dated 8 August and its author field is Anonymous. Its SHA256 is
b161d4bff3ebd2d29c275f57b08a4352de11dacb7441a99a8993bb0265b661ee.

The analytic input actually used is Hien--Roucairol,
[Integral representations for solutions of exponential Gauss--Manin
systems, arXiv:0704.1739v1](https://arxiv.org/html/0704.1739v1),
Theorem 1.1 (=3.7), together with its perfect period pairing (Theorem 2.4).
It represents local horizontal functionals by flat rapid-decay cycles.
Apply it to the smooth affine surface C* x C, projection to t and the
regular function Phi below. The phase has degree at least four in z over
C(t), so projection and phase are algebraically independent. Remove
finitely many base values at which a chosen rational frame has poles or
its degree changes. Multiplication by 1/t is an invertible rational gauge
there. No assertion that the original divergent formal Gaussian integral
is an analytic integral is needed.

## Exact phase and the missing prefactor

The preceding formal Gaussian/Lagrange substitution, followed by
z=t*(a+Y), gives the polynomial-in-z rational phase

    F=z/t-a-z^2,
    Phi(t,z)=t*G(z,F)-F^2/2.                        (1)

The corresponding analytic period has normalization

    I_Gamma(t)=(1/t)*integral_Gamma exp(Phi(t,z)) dz. (2)

A fixed nonzero normalization constant is irrelevant. The factor 1/t is
essential: it comes from dY=dz/t. It must not be discarded when forming
the tangent connection.

Put

    V(z)=G(z,-a-z^2),
    L(z)=G_T(z,-a-z^2),
    M(z)=G_TT(z,-a-z^2)/2,
    s=coefficient of T^3 in G.

Taylor expansion, exact because deg_T G<=3, yields

    Phi=t*V+B0+B1/t+B2/t^2,
    B0=z*L-(a+z^2)^2/2,
    B1=z^2*M+z*(a+z^2),
    B2=s*z^3-z^2/2.                               (3)

There are only three degree regimes. This classification includes every
a and every lower-degree coefficient of G.

* If n=deg V is 4,5 or 6, the degree of Phi over C(t) is n. Every Bj
  has degree at most n. After division by t, Phi tends coefficientwise
  to V, with the same nonzero leading coefficient. All critical points
  remain bounded as t tends to infinity.
* If 1<=n<=3, comparison of the degree-six, five and four coefficients
  of V gives s=v1=0 and u2=v0 in the notation of the previous note.
  Thus B0=-z^4/2 plus a polynomial of degree at most three, while
  deg B1<=3 and B2=-z^2/2. The phase has degree four. The same statements
  hold when V is constant.
* If V=k is constant (including V=0), the phase is tk plus a quartic
  polynomial whose coefficients are holomorphic in 1/t at infinity
  and whose limiting quartic coefficient is -1/2.

For example the successive high coefficients of V are

    [z^6]V=-s,    [z^5]V=v1,
    [z^4]V=v0-u2-3a*s.

The checker verifies the entire identity (3), these implications and all
degree bounds universally, rather than sampling a.

## Why formal pure vanishing would give a constant analytic period

Let K=C(t), d=deg_z Phi and

    H=K[z]/D_Phi(K[z]),    D_Phi(q)=q_z+Phi_z*q.

This is a K-vector-space quotient, not a ring quotient. Leading-term
reduction gives the basis 1,z,...,z^(d-2); in particular its rank is d-1
and [1] is nonzero. The normalized tangent connection is

    nabla[q]=[q_t+(Phi_t-1/t)*q].                   (4)

It is well defined because the untwisted derivative commutes with
D_Phi modulo a D_Phi term; subtracting the scalar 1/t preserves its image.

Define a formal functional on K[z] by

    J_t(q)=E_Y[q(t*(a+Y))*exp(Psi(t,Y))],

where Psi is the previous polynomial formal phase. The formal identity
there gives J_t(1)=E(exp(tP)). Gaussian integration by parts proves
J_t(D_Phi q)=0. Differentiating at fixed Y, with z=t*(a+Y), gives

    (J_t(q))'=J_t(q_t+Phi_t*q-q/t),

since the omitted expression is J_t(D_Phi(z*q/t))=0. All expressions
are well-defined in C((t)): each formal t-coefficient is a polynomial
in Y and rational t-coefficients have finite pole order. Thus J_t is a
formal horizontal functional for (4).

Suppose every positive pure moment is zero, so J_t(1)=1. Let U be the
K-submodule generated by [1],nabla[1],nabla^2[1],... . Finite rank makes
U connection-stable. On its generators J_t has values 1,0,0,... . Hence
its restriction takes values in K, is a well-defined rational horizontal
functional, and takes [1] to 1. Any relation among the generators has
zero coefficient of [1] after applying J_t; this justifies well-definedness
without assuming the seed cyclic in H.

At an ordinary complex base point, extend this functional arbitrarily
from the fibre of U to the fibre of H and solve the full dual analytic
connection. Horizontal uniqueness makes its restriction agree with the
functional on U. The Hien--Roucairol period comparison therefore produces
a local flat rapid-decay cycle (a complex linear combination of cycles)
whose normalized period (2) is identically 1. Continue it to a sufficiently
small generic sector at t=infinity, avoiding the finite singular set.
It remains 1. We now show no such period exists.

## All normalized periods at infinity lack a constant branch

The following argument is ordinary one-dimensional steepest descent,
including degenerate bounded critical points. Distinct exponential factors
are separated on a generic sector; coincident factors are grouped before
comparing expansions. A zero scalar period causes no difficulty.

### Degree n>=4

Choose a sector avoiding the finitely many Stokes directions of V.
Rapid-decay paths for V also give a basis of rapid-decay paths for Phi
for sufficiently large t: the degree is constant and its leading
coefficient is t*V_n+O(1). Work with fixed paths of V and expand the
bounded-coefficient factor exp(B0+B1/t+B2/t^2).

Near a critical point z0 of V, write

    V(z)=v+u^k,           k>=2,

in a local analytic coordinate u. The k inverse branches and their
adjacent differences give k-1 local paths, including multiple critical
points. Expanding the analytic amplitude gives the elementary Watson
integrals of u^j*exp(t*u^k). Their powers are t^(-(j+1)/k), j>=0.
The additional coefficients from B1/t,B2/t^2 contribute only further
nonpositive powers. Therefore every local scalar period has expansion

    exp(t*v) * t^(-1-1/k) * series in t^(-1/k),     (5)

or vanishes to higher order; the outer -1 is the prefactor in (2).
The sum of the local dimensions k-1 is n-1, the entire period rank.
The standard local paths can be joined to decaying sectors; if several
critical values coincide, keep their local paths and group that common
exponential. A zero-exponential block (v=0) strictly decays. Every other
block has nonzero exponential factor. Neither supplies a nonzero constant.

This is also a direct treatment of the degenerate-critical-point case:
it does not assume a nonzero Hessian and does not infer an all-coefficient
claim from a generic Morse specialization.

### Degree 1<=n<=3

The n-1 bounded critical points, counted with multiplicity, have precisely
the same local expansions (5). There are 4-n further critical points.
Write V_n for the leading coefficient and set

    z=t^(1/(4-n))*zeta.

After scaling by t^(4/(4-n)), the leading escape phase is

    H(zeta)=V_n*zeta^n-zeta^4/2.

Its nonzero critical points satisfy

    zeta^(4-n)=n*V_n/2.

There are 4-n of them, all simple. At each,

    H(zeta)=(1-n/4)*V_n*zeta^n !=0,
    H''(zeta)=2*(n-4)*zeta^2 !=0.                 (6)

Terms omitted from the leading escape phase have smaller powers of t;
they cannot cancel its nonzero leading critical action. Holomorphic Morse
reduction and a Gaussian expansion give an exponential with nonzero
leading action of degree 4/(4-n) in t, times a prefactor starting at

    t^(-1-1/(4-n)).                               (7)

Coincident leading actions may require grouping or smaller action terms;
each escaping block still has a nonzero highest action. Its projection
can vanish but cannot become a nonzero constant. The n-1 bounded and
4-n escaping dimensions total three, the full quartic period rank.
Thus again no constant branch remains.

### V constant

In this case choose fixed quartic rapid-decay contours of the limiting
phase. The integral without its factor 1/t is holomorphic in epsilon=1/t
near zero, by quartic decay uniform in epsilon. Every normalized period
is therefore

    exp(t*k)*t^(-1)*f(1/t),    f holomorphic at zero. (8)

For k=0 it tends to zero; for k!=0 it has a nonzero exponential factor.
Neither is a nonzero constant. Persistent multiple finite critical points
of the limiting quartic do not affect the contour domination argument.

These regimes exhaust all G and a, contradicting the constant analytic
period forced by pure vanishing. This proves the theorem.

## A sufficient 320-moment gate for the full architecture

**Corollary.** The first 320 zero positive pure moments suffice for
all-order pure vanishing in the fourteen-coefficient architecture. In its
c*d!=0 interior at least one of those moments is nonzero. This is a
sufficient bound, not asserted optimal.

Here is an explicit rational-connection estimate, including all rank-drop
strata. Write d=deg_z Phi, so d is 4,5 or 6 and r=d-1. From (3), every
coefficient of t^2*Phi is polynomial in t of degree at most three. Write

    Phi_z=t^(-2)*sum_i A_i(t)*z^i,
    A=A_(d-1) !=0,       deg_t A_i<=3.

Leading-term reduction uses

    z^(k+d-1) = -(t^2*k/A)*z^(k-1)
                -sum_(i<d-1) (A_i/A)*z^(k+i)

in the vector-space quotient. To reduce any monomial appearing in
z^j*(Phi_t-1/t), 0<=j<=d-2, at most d leading cancellations are needed:
its largest possible z-degree is 2d-2 and the final degree is at most d-2.
Each cancellation introduces at most one A denominator and numerator
t-degree at most three (the derivative term contributes t^2).
Also Phi_t-1/t has denominator t^3 and numerator t-degree at most three.
Consequently the connection matrix is N/H with

    H=t^3*A^d,    deg H<=B=3d+3,    deg N_ij<=B.    (9)

This is allowed to be an overestimate. It remains valid at every fixed
coefficient point after choosing its actual d; poles of H at t=0 cause
no problem when polynomial denominators are cleared.

The derivative columns v_j=nabla^j[1] have the form

    v_j=N_j/H^j,    deg N_j<=j*B,
    N_(j+1)=H*N_j' - j*H'*N_j + N*N_j.            (10)

Let q be the least index at which v_0,...,v_q are dependent, so q<=r.
The first q columns are independent even if the seed is not cyclic in H.
Choose q rows on which they are independent. Cofactors of the resulting
q by (q+1) polynomial-column matrix N_0,...,N_q give a nonzero scalar
relation

    sum_(j=0)^q a_j(t)*v_j=0,
    deg a_j <= D=B*q*(q+1)/2,       a_q!=0.        (11)

Indeed the cofactor omitting column j has degree at most
B*(q*(q+1)/2-j); multiplication by H^j converts N_j to v_j and restores
the displayed uniform bound. The relation extends from the chosen rows
to all rows because the entire q+1-column matrix has rank q.

Its coefficient a_0 is nonzero. Otherwise its minimal-order scalar
equation would admit the constant solution 1; evaluating v_0 at 1 and
all later v_j at 0 would define a rational horizontal functional on the
seed submodule. The same analytic extension and period comparison used
above would produce the impossible constant normalized period.

Apply (11) to the always-defined formal I=E(exp(tP)). If the first
D+q moments vanished, I-1 would have order at least D+q+1. Hence

    a_0=-sum_j a_j*(I-1)^(j)

would have order at least D+1. This contradicts a_0!=0 and deg a_0<=D.
Taking q<=d-1 gives the uniform table

| Actual phase degree d | B | Maximal D | Sufficient nonzero-moment prefix D+q |
|---|---:|---:|---:|
| 4 | 15 | 90 | 93 |
| 5 | 18 | 180 | 184 |
| 6 | 21 | 315 | 320 |

Normalization preserves each pure-moment zero condition, so 320 excludes
the original mixed interior as well. On c=0 the previous 35-moment gate
implies all-order pure vanishing, and on d=0 the previous 92-moment gate
does so. Both are at most 320. Therefore 320 zero moments in the entire
fourteen-coefficient family force the candidate onto a completed face
and imply all-order vanishing and GMC safety there.

This is a degree proof, not a 320-step numerical experiment or an ideal
certificate. It closes the previously missing explicit full-parent gate
without rerunning the timed-out standard-basis calculation.

## Verification and remaining obligations

The [checker](verify_mixed_period_infinity.py) verifies the universal phase
expansion, degree stratification, formal connection bridge, exact escape
actions and Hessians, dimension counts, negative powers and the finite-gate
degree table. Exact sample connection matrices verify the denominator
conventions of (9); the general degree proof is written above. It compares
the formal change with direct Wick moments on each degree regime. These
checks support the explicit algebra, not a computational proof of analytic
period comparison or all steepest-descent assertions. Those are supplied
by the written argument and the stated imported theorem. No producer code
from the newly inspected preprint is executed.

The finite parent equations and their timeout are retained as evidence
for the earlier algorithm, not an outstanding route to full-parent
existence. The explicit sufficient 320-moment gate above is not asserted optimal.
The global ternary cubic problem, omitted WZT/WT^2/higher-W couplings,
external specialist review and formal verification remain open.

The subsequent [logarithmic conic argument](QUADRATIC_W_LOG_PERIODS.md)
closes WZT and all nondegenerate quadratic W coefficients. The subsequent [rank-one classification](RANK_ONE_WT2_CLASSIFICATION.md)
completes all W-linear cubics, with sufficient gate888 and sharp3*deg-Q
mixed cutoff. Higher W powers and the global frame criterion remain open.
