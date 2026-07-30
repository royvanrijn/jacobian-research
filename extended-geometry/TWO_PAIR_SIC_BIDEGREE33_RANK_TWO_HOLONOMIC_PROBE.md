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
their common parameter denominator has not been found, and no exceptional
locus has yet been classified.  In particular, this note does not prove
that bidegree \((3,3)\) is safe and does not produce an all-order
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
The observed compact recurrence order agrees with this number.
That agreement is evidence, not a rank theorem in this family.

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

1. Build the relative creative-telescoping linear system for (2.4) with
   recurrence order 27 and \(m\)-degree 11.
2. Normalize its forward coefficient to (1.3) over several generic
   modular parameter points.
3. Reconstruct the common parameter denominator \(D(U,W)\) and the
   remaining \(p_j\), preferably in gauge-invariant determinantal
   coordinates.
4. Store the divergence certificates, including the \(t=0,1\) endpoint
   terms, and verify the polynomial identities independently.
5. Reduce \(D\) and the missing bridge moments against the corrected
   rank-two moment ideal, then repeat on every component of \(D=0\).

The modular probe has fixed a plausible compact order and exposed a
forward coefficient with no integer singular steps.  Universal
certificate reconstruction and exceptional-locus stratification remain
the proof-producing gates.

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
