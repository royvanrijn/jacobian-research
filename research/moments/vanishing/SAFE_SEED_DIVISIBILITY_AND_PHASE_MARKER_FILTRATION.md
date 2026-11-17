# Safe Gaussian seeds, original phase markers and the primitive-vacuum problem

17 November 2026. The constant-period-marker deduction in Section 1 was
supplied by the user. This note records its consequences and an elementary
reversible Gaussian construction. These are written proofs and reductions;
general field safety implies purity remains OPEN. Independent replay,
formal verification and external review remain false; priority unassessed.

## 1. Safety is already a rank-one problem

Let K be a characteristic-zero field, C a fixed nonsingular rational
symmetric Gaussian covariance, P in K[X], and use the common normalization

    Phi[A](t)=sum_(m>=0) E_C[A P^m] t^m/m!, f=Phi[1].

An ORIGINAL marker means a polynomial A in K[X], independent of t. Set
O=K[t], V=Phi[K[X]], M=O V and

    F_d=O-span{Phi[A]:deg_X A<=d}.

Thus F_d is the original O-module marker filtration; V need not be an
O-module. Safety means EVERY original Phi[A] is a polynomial. Purity
means f=1. Time differentiation is the actual original-marker identity

    partial_t^j Phi[A]=Phi[P^j A].                       (1)

Suppose first only that f is polynomial, of degree r. Put
mu_r=E_C[P^r]=f^(r), a nonzero scalar. Then the FIXED ORIGINAL marker

    U=P^r/mu_r has Phi[U]=1.                            (2)

For r=0 this says U=1. For r>0 it has spatial degree r deg P.
Equation (2) needs polynomiality of f, not safety. If safety holds as
well, M is contained in O; (2) gives the reverse inclusion. Therefore

    safe => M=O, rank_(K(t)) M=1.                       (3)

This explicitly connects the first safety/purity gap to the existing
[rank-one classification](STEIN_ACTION_AND_RANK_ONE_LATTICES.md). The
rank-at-least-two regular-realization problem is a separate obligation.

## 2. The proposed divisibility condition is exactly the remaining gap

For a SAFE phase the following conditions are equivalent:

* f=1.
* f divides Phi[A] in O for every original marker A.
* f divides Phi[P]=f' in O.
* M=fO.
* F_0=fO is stable under partial_t.
* The intrinsic original marked Stein action is locally nilpotent.
* The finite invertible chosen-seed Stein summand is zero.

The first implication uses safety. The second implies the third by
taking A=P. A nonconstant polynomial in characteristic zero cannot
divide its lower-degree derivative, so the third forces f=1, since
f(0)=1. Alternatively the user's marker (2) gives f divides 1. The
fourth follows from M=O; the fifth is the same derivative test. The
[intrinsic Stein theorem](STEIN_ACTION_AND_RANK_ONE_LATTICES.md) proves
the last two equivalences: on M its action is

    B_f(h)=(h-h(0)f)/t,
    O=fO direct_sum O_(degree<r),
    charpoly(B_f on O_(degree<r))=xi^r f(1/xi).          (4)

For r=0 the second summand is zero. For r>0 it is invertible.
Local FINITENESS already follows from rank one; it cannot replace local
NILPOTENCE here. Nor does regular rank-one arithmetic close the gap:
a nonconstant polynomial seed is already algebraic and globally regular.

The divisibility proposal is a rigidity TARGET, not a proved property
of every safe Gaussian phase. Even its single instance A=P is equivalent
to purity under safety. The original degree-zero piece F_0 and the unit
inside M must therefore be retained. The abstract lattice O alone does
not distinguish seeds 1 and 1+t.

## 3. What phase powers do, and do not, force

Continue to assume only that f is polynomial of degree r. For a formal
phase variable S, the map

    K[S] -> K[t], H(S) |-> H(partial_t) f

has kernel (S^(r+1)) and image K[t]_(degree<=r). Indeed
f,f',...,f^(r) have distinct degrees r,r-1,...,0. If P is nonconstant,
K[S] identifies with the original subalgebra K[P]. Within this subalgebra,
a nonzero constant period first occurs at phase-polynomial degree r.
This is not a lower bound for all spatial markers outside K[P].

The O-span of the first j+1 phase-power periods is the ideal

    J_j=(f,f',...,f^(j))O=(q_j)O.                       (5)

Normalize q_j(0)=1. Over a splitting field, let a run over the roots
of f, with multiplicities m_a; all a are nonzero. Then

    q_j=product_a (1-t/a)^max(m_a-j,0).                 (6)

In particular J_j=O once j reaches m_max=max_a m_a, with m_max=0
when f=1. For a squarefree seed J_1=O by Bezout. This latter recovery
allows TIME-DEPENDENT O-coefficients. It must not be confused with
the fixed original marker (2), whose phase-polynomial degree is r.

If P is safe, every F_d is an ideal (g_d) of O, normalized g_d(0)=1,
and

    g_0=f, g_(d+1) divides g_d, g_d divides f.

Write D=deg P when P is nonconstant. Equation (1), including the
Leibniz rule on O-coefficients, gives

    partial_t F_d subset F_(d+D),
    g_(d+D) divides gcd(g_d,g_d'),
    F_(m_max D)=O=M.                                   (7)

The last assertion follows either by iteration or by (5)-(6), because
P^j has spatial degree jD. A constant phase with polynomial f is P=0;
then F_0=O already. These statements preserve the original filtration,
but finite-degree access to 1 does not turn the original unit f into 1.

In the faithful rank-one classification, safety is exactly d=1,
p=f, where

    M=(f/p)K[t,d^-1].                                  (8)

At a root a of f of multiplicity m_a, all normalized phase-power
periods f^(j)/f have pole order at most m_a. The marker (2) gives 1/f,
which attains the ALLOWED finite bound m_a. It does not exceed that
bound or restart an infinite pole ladder. Thus the constant-period
marker is consistent with the hypothetical nonpure safe lattice; it
does not by itself refute it. A marker above the bound would instead
force d nonconstant and UNSAFETY.

The retained actual unweighted polynomial-seed example with seed 1-t
has Phi[-P]=1 but full image K[t,(1-t)^-1]. It is rank one and UNSAFE:
another original marker exceeds the finite pole bound. See the
[polynomial-seed proof](POLYNOMIAL_SEEDS_AND_PURE_MARKER_DELAY.md).
The positive polynomially WEIGHTED example with seed 1+t/2 has full
image O and is safe. Its original unweighted phase has seed 1; it is
not a counterexample to unweighted purity. See the
[amplitude-content boundary](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md).

## 4. An exact auxiliary Gaussian exponential cancellation

Let psi in tK[t]. Introduce an independent Witt pair Z,W with
E[Z^a W^b]=delta_(a,b) a!, corresponding to rational nonsingular
covariance [[0,1],[1,0]]. Put

    H(Z)=-psi(Z)/Z, Q_psi(Z,W)=W+H(Z).                  (9)

This is an ORIGINAL polynomial phase. The Wick contraction identity
E[F(Z)exp(vW)]=F(v) gives, with an auxiliary source u,

    Phi_(Q_psi)[Z^a W^b]
      =partial_u^b[(t+u)^a exp(t H(t+u))] at u=0
      =exp(-psi(t)) R_ab(t), R_ab in K[t], R_00=1.    (10)

The identity holds coefficientwise formally; it makes no convergence
assumption. Derivatives are finite, so the polynomial assertion holds
for EVERY a,b. In particular its WHOLE actual image is exp(-psi)O,
not merely a supplied scalar solution of that form. This is a direct
special case of the retained linear-momentum Gaussian contraction.

For ANY original phase P in independent variables, form

    P_tilde(X,Z,W)=P(X)+Q_psi(Z,W).

Every product marker factors. Equation (10) proves both the full-image
identity and the exact O-module ORIGINAL marker-filtration identity

    M_tilde=exp(-psi)M,
    F_tilde_d=exp(-psi)F_d for every d>=0,
    f_tilde=exp(-psi)f.                                (11)

For the forward filtered inclusion expand a marker of total degree
at most d into A(X)Z^a W^b with deg A+a+b<=d. Its period belongs
to exp(-psi)F_d by (10). The reverse inclusion uses the very same
old markers A(X) with auxiliary marker 1. This equality concerns
O-module filtrations; it does not identify the spatial vector spaces
V or discard the original time-independent marker maps in (10).

This construction uses two extra variables and preserves a fixed
nonsingular rational covariance by orthogonal direct sum. Over a
coefficient field already split for the standard Gaussian it is also
a linear change of ordinary Gaussian coordinates. No fixed-dimension
or fixed-degree equivalence is asserted. Equations (9)-(11) themselves
work over every Q-algebra; the classification in Section 5 is field-only.

## 5. A reversible primitive-vacuum reformulation

Assume an ACTUAL rank-one Gaussian image M is finite as an O-module.
By (8), this is exactly d=1, and M=gO with g=f/p and g(0)=1. Since
M is stable under time differentiation,

    g'/g in O, g=exp(psi), psi in tK[t],
    f=p exp(psi), p in O, p(0)=1.                     (12)

Alternatively M is a finite torsion-free rank-one module over the PID
O, so has a generator g. Because every element is formal at zero and
M contains f with f(0)=1, that generator is a formal unit. D-stability
then proves (12). The normalized generator is unique. The polynomial
psi here is the GENERATOR EXPONENT; no identification with a chosen
critical-action construction is needed for this argument.

Apply (9) with this psi. Then (11) gives the actual identities

    M_tilde=O, f_tilde=p, F_tilde_d=exp(-psi)F_d.        (13)

Thus P_tilde is SAFE with ALL original markers controlled. If p is
nonconstant, it is a genuine unweighted safe nonpure phase, conditional
on the existence of the original nonprimitive finite image.

Across all dimensions and degrees, the following two open assertions
are therefore equivalent:

* Every safe Gaussian phase over characteristic-zero fields is pure.
* In every actual rank-one Gaussian image finite over K[t], the
  distinguished original unit generates that K[t]-lattice.

The second assertion means p=1 in (12), or M=fO. If the first holds,
(13) forces p=1. Conversely a safe phase already has M=O by (3),
so the second makes f an O-unit, hence f=1. The equivalence is in
the nonsingular rational-covariance Gaussian programme; it does not
silently restrict auxiliary variables to a positive real covariance.
Scalar extension to a splitting field converts the construction to
the ordinary algebraically closed Gaussian setting if desired.

This is a faithful Gaussian realization of exponential cancellation,
not just a formal gauge of an abstract connection. It strengthens the
rank-one reformulation while retaining F_0=pO after cancellation.
It does not establish either open assertion. The finite nonprimitive
vacuum p remains precisely the unresolved obstruction.

## 6. The all-ring boundary of the fixed-marker normalization

Let R be any Q-algebra and suppose f=Phi[1] is polynomial of degree r.
Taylor's identity at zero gives the literal actual-image identity

    1=sum_(j=0)^r (-t)^j f^(j)(t)/j!
     =Phi[sum_(j=0)^r (-tP)^j/j!].                    (14)

Thus safety implies M=R[t] over EVERY R as well. This uses a polynomial
TIME-DEPENDENT amplitude, not division by a possibly nonunit mu_r.
There is no generic field-rank assertion over an arbitrary ring.

For example take R=Q[e]/(e^2), standard covariance, P=eX^2. All
markers are safe since P^2=0, but f=1+et and mu_1=e is not a unit.
No element a+bP of R[P] has period 1: its period is
a+be+ae t; the constant equation a+be=1 implies ae=e, contradicting
the vanishing t coefficient. The time-dependent amplitude 1-tP in
(14) does have period 1. An original marker outside R[P], namely
(3-X^2)/2, also has period 1, since its expectation is 1 and its
product with X^2 has expectation 0. This is an all-ring normalization
boundary and a safe nilpotent nonpure example, not a field counterexample.

## 7. Proof boundaries and the open assertion

The derivative ideals, PID argument and Wick contraction above are written
algebraic derivations. The imported source results are the
[Stein/full-lattice theorem](STEIN_ACTION_AND_RANK_ONE_LATTICES.md) and
[polynomial-seed constructions](POLYNOMIAL_SEEDS_AND_PURE_MARKER_DELAY.md).

The primitive original-vacuum assertion remains open: exclude a nonconstant
f in a genuinely safe unweighted image, or construct an all-marker-safe
counterexample. The original marker filtration and its finite allowed pole
bounds must be preserved. Higher-rank regular transcendental realization
is a separate question.

Classical rational-slope inputs and the current polynomial implementation
are described in the [credit note](BERNSTEIN_CLASSICAL_CREDIT_AND_FORMAL_EXTENSIONS.md)
and [graph-certificate proof](OAKU_GRAPH_CERTIFICATES_AND_ACTUAL_POLYNOMIAL_DECISION.md).
No independent mathematical review or formal verification is claimed.
