# F7 triple-root incidence: a primitive section and rational-point density

14 November 2026. Registered claim **GAOF7MW1**.

The elliptic completion in [GAOF7CUSP1](../gao-f7-three-point/PROOF.md) has
trivial geometric Mordell-Weil torsion and an explicit primitive section
P defined over Q(s), with Shioda self-pairing 58/3. Consequently its rank
is at least one, and the rational points of the auxiliary triple-root
incidence surface are Zariski dense despite its nonrationality.

This is a theorem about the incidence surface. It does not construct a
three-point F7 fibre, exclude any remaining three-point multiplicity type,
compute the full Mordell-Weil rank or settle JC2. No independent end-to-end
review or novelty priority is claimed.

## Previously established geometric input

GAOF7CUSP1 defines the irreducible affine surface

    f(s,a,b)=R''(s;X(s,a,b))/2=0,

where derivatives hold the target fixed before substitution, and proves
that it is birational over Q to the genus-one fibration

    Z^2=H(s,a)=A(s)a^4+B(s)a^3+C(s)a^2+D(s)a+E(s).

The leading coefficient is A=L^2, with

    L=4s(3s^3-2)(3s^3-1).

Thus the generic fibre has rational infinity points. Its elliptic curve
is isomorphic over Q(s) to

    E_s: y^2=x^3-27I(s)x-27J(s),
    I=12AE-3BD+C^2,
    J=72ACE+9BCD-27AD^2-27B^2E-2C^3.

The exact model and its full geometric fibre configuration are retained
in the earlier [model input](../gao-f7-three-point/elliptic-cusp-input.json).
The relatively minimal surface has four I4 fibres at s=0 and the roots
of 3s^3-2, 64 I1 fibres, and one IV fibre at infinity. Its chi(O) is seven.
These statements, including minimality in every chart, are inputs here.

## The explicit section

For any quartic with leading coefficient A=L^2, the formulas

    x_P=3(3B^2-8AC)/(4A),
    y_P=27(B^3-4ABC+8A^2D)/(8L^3)

satisfy y_P^2=x_P^3-27Ix_P-27J identically. This identity follows simply
by substituting the displayed I and J; the checker verifies it with
independent symbolic coefficients before specializing. For our quartic,

    x_P=-3(6237s^12-6156s^9+180s^6+288s^3-40)/(3s^3-1)^2.

The full y_P, whose numerator has degree 29 and denominator (3s^3-1)^3,
is saved in [elliptic-section-input.json](elliptic-section-input.json).
Both SymPy and independent FLINT cross multiplication verify membership
in the exact model. The point is finite in the generic Weierstrass chart,
so it is distinct from the zero section O.

The Jacobian model and invariants are the standard quartic construction;
see Section 2 of [Fisher, On binary quartics and the Cassels-Tate pairing](https://www.dpmms.cam.ac.uk/~taf1000/papers/bq-ctp.pdf).
The new point assertion here uses the explicit identity, rather than
assuming a particular covariant map or its degree.

## Torsion is impossible and all nonzero heights are bounded below

Use the Shioda self-pairing convention h(Q)=<Q,Q>, not half of it. The
height formula is

    h(Q)=2chi(O)+2(Q.O)-sum_v contr_v(Q).

Its positivity, torsion kernel and correction table are in Theorem 11.5,
Section 11.8 and Table 4 of
[Schuett and Shioda, Elliptic surfaces](https://arxiv.org/pdf/0907.0298).
For an I4 fibre the possible self-corrections are 0,3/4,1,3/4. For a IV
fibre they are 0,2/3,2/3. The I1 fibres contribute zero. These numbers
also follow from the diagonal entries of the inverse A3 and A2 Cartan
matrices, which the checker reconstructs directly.

Every section Q distinct from O is a distinct irreducible curve, so
Q.O is nonnegative. Therefore every such section over C(s) satisfies

    h(Q) >= 14-4*1-2/3 = 28/3.                     (1)

A nonzero torsion section would have height zero, contradicting (1).
Thus E(C(s)) has trivial torsion; in particular so does E(Q(s)).

## Exact intersection data for P

Write q=3s^3-1. The numerators of x_P and y_P are coprime to q. At each
of its three distinct roots their orders are -2 and -3. These are good
fibres: gcd(q,disc(E_s))=1. The local parameter -x/y at O has order one,
so each contributes exactly one to P.O.

There are no other finite poles. At the four I4 places the denominator
is a unit and the numerator of y_P is nonzero. On the short nodal
Weierstrass cubic a singular point has y=0, so P specializes to a smooth
point there. It meets the identity component of each I4 fibre, giving
zero correction at all four.

At infinity use T=1/s and the minimal coordinates of GAOF7CUSP1,

    x_infinity=T^14*x_P(1/T),
    y_infinity=T^21*y_P(1/T).

Their orders are 8 and 1. The integral Weierstrass coefficients have
orders 5 and 2. Thus, writing a4=T^5*u(T) and a6=T^2*v(T) with u,v units,
the model is

    y_infinity^2=x_infinity^3+T^5*u*x_infinity+T^2*v.

Blow up its singular point in the T chart, setting x_infinity=T*x1 and
y_infinity=T*y1. The proper transform is

    y1^2=T*x1^3+T^4*u*x1+v.

The lifted section has x1(0)=0 and y1(0)=-839808, with
v(0)=(-839808)^2. It is a smooth point of one of the two exceptional
components, since the derivative with respect to y1 is nonzero. These
are the nonidentity components of the type IV fibre; the zero section
is away from this blowup. The section therefore has correction 2/3 at
infinity. It does not meet O there. Further resolution does not alter
this already smooth neighbourhood of the section.

Combining all places gives

    P.O=3,
    sum_v contr_v(P)=2/3,
    h(P)=14+6-2/3=58/3.                            (2)

## Primitivity and density

Suppose P=nQ for a section Q over C(s) and an integer |n|>=2. Bilinearity
of the height and (1) imply

    h(P)=n^2*h(Q) >= 4*(28/3)=112/3,

contradicting (2). Thus P is primitive even in the geometric group
E(C(s)), not just in E(Q(s)). Its infinite order proves rank at least
one; no larger lower bound or exact rank is asserted.

The multiples nP give infinitely many distinct sections over Q. Each
section is a rational curve over Q with a Zariski-dense set of Q-rational
points, obtained by choosing rational base parameters. A proper closed
subset of an irreducible surface contains only finitely many distinct
irreducible curves. The union of these sections is therefore Zariski
dense in the elliptic surface, and so is its set of Q-rational points.

The isomorphism of the generic genus-one curve with its Jacobian, after
choosing a rational infinity point, and the square-completion map to
f=0 are defined over Q(s). They induce a birational map of surfaces over
Q. Restrict to nonempty open sets on which it is an isomorphism. Removing
the complements can discard only finitely many whole sections. Rational
points on the original incidence surface are consequently Zariski dense
as well. This existence argument does not require a claimed formula
for the inverse quartic-to-Weierstrass map.

In particular, nonrationality does not justify an attempt to exclude
rational triple-root incidences. Their abundance says nothing by itself
about the additional constraints required for exactly three simple roots.

## Verification and the retained three-point attempts

```sh
timeout 60s .venv/bin/python research/jacobian/gao-f7-three-point-structure/verify.py
```

The checker performs 33 exact checks: the generic covariant identity,
model specialization and independent FLINT point identity; pole and
singular-place gcds; local infinity orders and the smooth exceptional
intersection; Cartan corrections and rational height inequalities; and
source hashes. The written section-extension, height and density arguments
are sourced geometric implications, not a formally verified proof.

The same directory retains unfinished work on the actual three-point gate:

- Two distinct Singular transformation-matrix algorithms, liftstd and its
  slimgb variant, each timed out at 30 seconds for the raw cubic-cube
  system modulo 32003. Both had 1536MB memory and 5MB output caps. No lift
  was obtained and no characteristic-zero conclusion follows.
- A balanced coefficient split Q=C-M/4, S=C+3M/4 gives the exact identity
  Q^3*S=C^4-3C^2M^2/8+CM^3/8-3M^4/256. The square-root expansion again
  exposes the earlier Delta denominator. Its necessary four-variable
  equations are exploratory inputs, without a unit certificate.
- The six coefficient equations in degrees 7,3,5,4,1,2 already yield a
  modular unit ideal. This suggests the degree-zero equation may be
  unnecessary, but proves no characteristic-zero exclusion. The five-
  equation basis output hit its 1MB cap, so no complete basis is claimed.

All inputs, partial outputs and receipts are retained; default verification
reruns none of these searches. A targeted literature refresh on 14 November
found no newer F7 result adopted here; this is not an exhaustive priority
search. These attempts themselves give no coefficient exclusion. The later
[GAOF7CUBE1 proof](../gao-f7-three-point-certificate/PROOF.md) excludes the
full cubic-cube family through exact local identities and a full-support
integer determinant. The later
[GAOF7MULT1 proof](../gao-f7-linear-cube/HIGH_ROOT_PROOF.md) proves the
sharp maximum-five bound and excludes(7,2). The subsequent
[GAOF7FIVE1 proof](../gao-f7-quintuple/PROOF.md) classifies all118 fivefold
targets and leaves only(4,3,2),(3,2,2,2) and global three-point existence open. MATH_STATUS.json is the mathematical-status
authority.
