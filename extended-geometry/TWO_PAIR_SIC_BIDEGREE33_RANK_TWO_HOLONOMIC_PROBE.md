# Holonomic probe on the rank-two bidegree-\((3,3)\) stratum

## 1. Outcome and status

There is now a concrete recurrence target for the first rank-two
bidegree-\((3,3)\) experiment.  Put

\[
 C=UW,\qquad
 U\in\operatorname {Mat}_{4\times2},\quad
 W\in\operatorname {Mat}_{2\times4},
 \tag{1.1}
\]

and normalize the moments by

\[
 \nu_m(C)=\frac{\mu_m(C)}{(3m+1)!}.                       \tag{1.2}
\]

At two unrelated integral exact-rank-two points and over each of the three
primes \(1000003,1000033,1000037\), 501 exact values of \(\nu_m\) admit a
scalar recurrence of order \(27\) whose coefficients have degree at most
11 in \(m\).  Each fit uses 335 equations and passes the remaining 139
equations.  At order 27, a degree-10 ansatz fails at both points modulo
\(1000003\).

A second calculation identifies the natural relative-cohomology rank.
The logarithmic Jacobian quotient has exact length \(18\), split as two
length-two endpoint pieces and a length-fourteen interior piece.  At the
same six point/prime combinations, an order-18 recurrence of
\(m\)-degree 18 passes 83 unused equations; degree 17 fails at both
points modulo \(1000003\).  Thus order 27 is a low-\(m\)-degree,
desingularized tradeoff rather than the minimal cohomological order.

After making the coefficient of \(m^{11}\nu_{m+27}\) monic, all six fits
have the same forward coefficient:

\[
\boxed{
 L(m)=
 \frac{
 \displaystyle
 \prod_{k\in\{71,73,74,76,77,79,80,82\}}(3m+k)
 \left(9m^3-1836m^2+210730m-17920380\right)}
 {9\cdot3^8}.}
 \tag{1.3}
\]

This is strong modular evidence for a compact universal recurrence shape.
It is **not** a creative-telescoping certificate.  The other recurrence
coefficients have not been reconstructed in the rank-two parameter ring,
and their universal common parameter denominator has not been found.  A
single generic factor pencil now has an exact modular border-basis
denominator calculation, described in Section 4.4, but its exceptional
fibers have not yet been classified.  In particular, this note does not
prove that bidegree \((3,3)\) is safe and does not produce an all-order
counterexample.

## 2. Exact relative-period realization

Write \(U=(u_{iq})\), \(W=(w_{qj})\), with \(0\leq i,j\leq3\) and
\(q=1,2\).  Define

\[
\begin{aligned}
 A_q(u)&=\sum_{i=0}^3u_{iq}u^{3-i},\\
 B_q(u,t)&=\sum_{j=0}^3
 w_{qj}u^{j-3}t^j(1-t)^{3-j},\\
 P_{U,W}(u,t)&=A_1(u)B_1(u,t)+A_2(u)B_2(u,t).
\end{aligned}                                             \tag{2.1}
\]

If \(\Phi_C(x,y)=\sum c_{ij}x^iy^j\), direct multiplication gives

\[
 P_{U,W}(u,t)=
 \Phi_C\left(1,u,t,\frac{1-t}{u}\right).                  \tag{2.2}
\]

The beta identity applied to the diagonal contraction formula gives

\[
\boxed{
 \nu_m(C)
 =\operatorname {CT}_u\int_0^1P_{U,W}(u,t)^m\,dt.}
 \tag{2.3}
\]

Thus

\[
 \sum_{m\geq0}\nu_m(C)s^m
 =\operatorname {CT}_u\int_0^1
 \frac{dt}{1-sP_{U,W}(u,t)}.                              \tag{2.4}
\]

The checker verifies (2.2) exactly and verifies (2.3) through order four
at both integral factor points.  The factorial normalization in (1.2) is
essential: (2.4) is not the ordinary generating function of the raw
\(\mu_m\).

## 3. The toric number and the forced endpoint discriminant

The generic exponent polygon of \(P_{U,W}\) is

\[
 \operatorname {conv}\{(-3,0),(0,0),(3,3),(-3,3)\}.
 \tag{3.1}
\]

Its Euclidean area is \(27/2\), hence its normalized volume is \(27\).
The low-\(m\)-degree recurrence order agrees with this number, while the
relative logarithmic quotient and the lower-order recurrence both have
rank 18.  Neither number follows automatically from the other.

Indeed, the beta substitution forces the \(u=-3\) face to be

\[
 c_{30}u^{-3}(1-t)^3.                                    \tag{3.2}
\]

For \(c_{30}\ne0\), this face has a triple torus root at \(t=1\).
Consequently the constrained beta family lies on the face discriminant of
the unconstrained Laurent-polynomial family.  A generic stable-GKZ rank
statement therefore cannot be applied to (2.4) without a relative
endpoint analysis.  The number 27 does not by itself prove the recurrence
or bound initial data.

## 4. The modular recurrence probe

For either factor point, the adjacent C++ helper computes the moments
without expanding a sixteen-variable polynomial.  It uses

\[
 F^m=\sum_{k=0}^m\binom mk
 A_1^kA_2^{m-k}P_1^kP_2^{m-k}                            \tag{4.1}
\]

and updates all binary-cubic products by multiplication by one cubic.
This gives the exact sequence modulo the chosen prime.

For order \(R=27\) and coefficient degree \(D=11\), the ansatz

\[
 \sum_{j=0}^{27}
 \left(\sum_{e=0}^{11}a_{j,e}m^e\right)\nu_{m+j}=0
 \tag{4.2}
\]

has \(28\cdot12=336\) unknowns.  The probe uses 335 equations to select a
null vector and then tests all remaining available equations through
\(m+27=500\).  The 139 holdout equations vanish in all six
point/prime combinations.

The monic \(j=27\) coefficient reconstructed from the three primes is
(1.3).  Its eight linear factors are positive for every \(m\geq0\).
The cubic has no integer root because it has no root modulo \(29\).
Therefore

\[
 L(m)\ne0\qquad(m\in\mathbb Z_{\geq0}).                  \tag{4.3}
\]

This eliminates integer singular steps from the common \(m\)-part of the
observed forward coefficient.  It does not eliminate a parameter factor
introduced when a universal recurrence is cleared of denominators.

### 4.1 The exact relative logarithmic quotient

Put \(Q=u^3P\).  The two logarithmic critical polynomials are

\[
 A=uQ_u-3Q,\qquad C=t(1-t)Q_t.                            \tag{4.4}
\]

At the first integral exact-rank-two point, exact characteristic-zero
calculation gives

\[
\begin{aligned}
 \dim_\mathbb Q\mathbb Q[u,t]/(A,C)&=30,\\
 (A,C):u^\infty&=(A,C):u^6,\\
 \dim_\mathbb Q\mathbb Q[u,t]/((A,C):u^\infty)&=18.
\end{aligned}                                             \tag{4.5}
\]

The saturated quotient has standard monomial basis

\[
\begin{gathered}
1,u,u^2,u^3,u^4,\quad
t,ut,u^2t,u^3t,\quad
t^2,ut^2,u^2t^2,\\
t^3,ut^3,u^2t^3,\quad t^4,ut^4,\quad t^5.
\end{gathered}                                             \tag{4.6}
\]

Saturating further by the two endpoints decomposes the length as

\[
 18=2_{t=0}+2_{t=1}+14_{\mathrm{interior}}.               \tag{4.7}
\]

This is the first exact finite relative-cohomology model for the proposed
recurrence.

### 4.2 The natural order-18 recurrence shape

At both integral rank-two points and all three primes, the ansatz

\[
 \sum_{j=0}^{18}
 \left(\sum_{e=0}^{18}b_{j,e}m^e\right)\nu_{m+j}=0
 \tag{4.8}
\]

uses 360 fitting equations and passes 83 unused equations through order
460.  An \(m\)-degree-17 ansatz fails at both points modulo \(1000003\).
The monic forward coefficient factors as

\[
 \prod_{k\in\{44,46,47,49,50,52,53,55\}}(3m+k)
 \cdot H_{U,W}(m),                                       \tag{4.9}
\]

where \(H_{U,W}\) is a point-dependent monic decic.  The degree-eight
factor is common to all six probes.  The decic records
parameter-dependent apparent singularities of the lower-order scalar
operator.  Passing to order 27 removes them from the observed forward
\(m\)-factor, which is why the higher-order relation is preferable for
uniform forward propagation.

A naive certificate polynomial of \(m\)-degree 17 cannot prove (4.8).
After shifting the recurrence to the exponent used in integration by
parts, its leading \(m^{18}\) coefficient has nonzero remainder in every
one of the eighteen basis coordinates (4.6), even after saturation by
\(u\).  The certificate must therefore use a rational-in-\(m\) relative
connection, a higher-degree syzygy cancellation, or an equivalent
desingularized operator.

### 4.3 Why expanded universal interpolation is not the next step

On the scaling family \(C(s)=(1+s)C_0\), the modular rational
reconstructor recovers the predicted coefficient degrees \(27-j\) for
the order-27 recurrence, with denominator one and independent holdouts.
This validates both the recurrence normalization and the interpolation
engine.

On the generic quadratic factor pencil

\[
 (U(s),W(s))=(U_0+sU_1,W_0+sW_1),                        \tag{4.10}
\]

256 samples with twelve holdouts find no rational interpolant for eight
representative recurrence coefficients within the resulting degree
window.  In particular, no candidate appears with combined
numerator/denominator degree at most 243 on this one pencil.  Therefore a
fully expanded twelve-parameter recurrence is the wrong immediate
representation.  The relative connection on the eighteen-element basis
(4.6), followed by determinant and desingularization operations, is the
compact target.

### 4.4 Exact border-basis denominator on the generic pencil

There is nevertheless an exact compact calculation available on the
pencil (4.10).  Over each of
\(\mathbb F_{1000003}(s)\), \(\mathbb F_{1000033}(s)\), and
\(\mathbb F_{1000037}(s)\), saturating the logarithmic Jacobian ideal

\[
 (uQ_u-3Q,\ t(1-t)Q_t):u^\infty
\]

has exponent six, quotient length eighteen, and the same six leading
border monomials

\[
 u^5,\quad u^4t,\quad u^3t^2,\quad u^2t^4,\quad ut^5,\quad t^6.
 \tag{4.11}
\]

After reduction and monic normalization, every border relation has
nineteen terms.  Their denominator degrees, in Singular basis order, are

\[
 74,\ 74,\ 88,\ 94,\ 74,\ 74.                           \tag{4.12}
\]

At all three primes the degree-\(74\) denominator is the common gcd of
the three distinct denominators.  Dividing it from the degree-\(88\) and
degree-\(94\) denominators leaves coprime factors of degrees \(14\) and
\(20\).  Hence this monomial chart is regular away from one squarefree
degree-\(108\) polynomial in \(s\).

This is the first explicit exceptional polynomial in the programme, but
only on one modular pencil and only for the logarithmic border basis.  Its
roots can include monomial-chart artifacts, and it is not the universal
telescoping denominator \(D(U,W)\).  The stable profile at three primes is
evidence for a characteristic-zero pencil determinant of the same
degree, not a reconstruction of that determinant.

## 5. The actual finite-determination statement to prove

Let

\[
 R_2=\mathbb Q[u_{iq},w_{qj}]
\]

be the factor parameter ring.  The modular data suggest a recurrence over
\(\operatorname {Frac}(R_2)\) of the form

\[
 \sum_{j=0}^{27}p_j(m;U,W)\nu_{m+j}=0,\qquad
 \deg_m p_j\leq11,                                      \tag{5.1}
\]

whose monic forward coefficient is \(L(m)\).  After clearing a primitive
common denominator \(D(U,W)\), the forward coefficient should have the
form

\[
 D(U,W)L(m),                                             \tag{5.2}
\]

possibly after removing a common parameter content.

An exact telescoping identity proving (5.1), together with (5.2), would
give the desired first open stratum:

> On \(D\ne0\), 27 consecutive zero moments beginning at any positive
> order propagate to all later orders.

The same statement transfers to the raw \(\mu_m\) after multiplying (5.1)
by \((3(m+27)+1)!\).  This only introduces explicit nonzero factorial
ratios; it does not create characteristic-zero singular steps.

The exceptional analysis is then finite and algebraic only after
\(D(U,W)\) is known.  On each component of \(D=0\), one must specialize
the telescoping matrix, remove its parameter content, and recompute its
generic rank and forward coefficient.  Merely requiring the leading
coefficient to be nonzero as a polynomial in \(m\) is insufficient:
integer-step zeros must also be audited.  Formula (4.3) performs that
audit only for the common open-stratum factor observed here.

## 6. Interaction with the corrected moment system

The corrected ambient system is

\[
 \mu_1,\ldots,\mu_{12},\mu_{14}=0.                       \tag{6.1}
\]

The order-27 recurrence shape does not make (6.1) an all-order condition
by itself.  Even on \(D\ne0\), forward propagation currently asks for
\(\mu_1,\ldots,\mu_{27}\).  The corrected equations omit \(\mu_{13}\) and
do not include orders \(15,\ldots,27\).

The natural order-18 relation narrows the missing bridge set to

\[
 \mu_{13},\mu_{15},\mu_{16},\mu_{17},\mu_{18}.            \tag{6.2}
\]

It does not yet remove those five values, and its parameter-dependent
decic requires an integer-step audit or desingularization.

There are three exact ways the recurrence could still settle (6.1):

1. reduce \(\mu_{13},\mu_{15},\ldots,\mu_{27}\) to zero in the radical of
   the corrected rank-two ideal on \(D\ne0\);
2. derive a lower-order recurrence after quotienting by that ideal; or
3. show directly that (6.1) has only nullcone points on \(D\ne0\), then
   recurse on \(D=0\).

Until one of these calculations is certified, (6.1) remains undecided on
the rank-two stratum.

## 7. Next certificate-producing computation

The immediate target is no longer an unspecified Picard--Fuchs search.
It is the following fixed-size reconstruction problem.

1. Construct the relative connection on the basis (4.6) over the generic
   rank-two factor field, retaining the six border-pivot determinants
   instead of expanding them.  The degree-\(74+14+20\) pencil
   factorization in Section 4.4 is the specialization target for this
   calculation.
2. Take a cyclic determinant to recover the order-18 operator and the
   parameter-dependent decic \(H_{U,W}(m)\).
3. Desingularize the apparent decic steps to the observed order-27,
   \(m\)-degree-11 operator and normalize its forward coefficient to
   (1.3).
4. Store the relative divergence certificates, including the
   \(t=0,1\) endpoint terms, and verify their polynomial identities
   independently.
5. Reduce the five bridge moments in (6.2) and the border-pivot
   determinant against the corrected rank-two moment ideal, then repeat
   on every rank-drop component.

The modular probes have now separated the natural relative rank from its
uniform desingularized form and exposed one exact modular pencil
denominator.  Connection construction, certificate reconstruction, and
recursive analysis of the degree-\(108\) pencil fibers remain the
proof-producing gates.

## 8. Reproduction

Run

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_holonomic_probe.py
```

The command compiles the adjacent C++ helper in a temporary directory,
checks the exact period identities, runs the six modular recurrence
probes, checks the two degree-10 failures, and writes
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_holonomic_probe.json`.
The artifact labels the recurrence as modular evidence rather than a
certificate.

The relative-Jacobian and natural-order calculation is replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_relative_jacobian.py
```

It writes
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_relative_jacobian.json`.
The exact modular pencil border basis is replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_pencil_border.py
```

It writes
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_pencil_border.json`.
The three finite rational-function-field calculations certify the stable
denominator profile \(74,\ 74+14,\ 74+20\) and the squarefree
degree-\(108\) chart polynomial.  They do not reconstruct a universal
parameter determinant or classify its exceptional fibers.

The generic and scaling pencil experiments are implemented in
`scripts/explore_two_pair_sic_bidegree33_rank_two_recurrence_line.py`;
their two artifacts are
`two_pair_sic_bidegree33_rank_two_recurrence_line.json` and
`two_pair_sic_bidegree33_rank_two_recurrence_scaling_line.json` in
`artifacts/generated-results/`.  They remain interpolation experiments
rather than recurrence certificates.
