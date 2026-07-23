# The two-real-variable Gaussian frontier

## 1. Result

Let \(X,Y\) be independent standard real Gaussians and put

\[
 Z=\frac{X+iY}{\sqrt2},\qquad
 W=\frac{X-iY}{\sqrt2},\qquad U=ZW.
\]

Then

\[
 \mathbb E(Z^aW^b)=\delta_{a,b}a!.                 \tag{1.1}
\]

[Long's explicit counterexample](https://arxiv.org/abs/2607.18186) proves
that GMC fails in every real Gaussian dimension at least three.  GMC in
dimension one is known, so dimension two is the unique unresolved dimension.
This note does not settle it.  It gives
six positive theorems that force a possible witness beyond the natural
low-complexity collapses:

1. GMC holds for every polynomial of total degree at most two, in every
   Gaussian dimension.
2. In two real variables, a pure-moment-zero polynomial supported on at most
   two rotational weights reduces to one-sided weight support and satisfies
   GMC.
3. In two real variables, the complete ansatz
   \[
     P=W h(Z)+v(Z)
   \]
   satisfies GMC; if all its pure moments vanish, it is either a scalar
   multiple of \(W\) or has strictly positive weight.
4. No total-degree-three polynomial with exactly three nonzero rotational
   weights of both signs can have all pure moments zero.  This is an exact
   finite-chart elimination, not a coefficient search.
5. Of the 33 mixed-sign four-weight supports available to a cubic, 29 are
   impossible.  Only four explicit supports, comprising 24 coefficient
   charts before the final chart-level test, survive the support census.
6. Every mixed-sign four-weight cubic is impossible.  Three exact
   good-prime quotient certificates cover the four charts of the symmetric
   support, and seven rational unit-ideal certificates cover the other
   three supports and 20 charts.  Thus a cubic counterexample needs at least
   five rotational weights.

Consequently a counterexample to \(\operatorname{GMC}(2)\), if one exists,
must have all of the following:

\[
\boxed{
\deg P\geq3,\quad
|\operatorname{wt}(P)|\geq3,\quad
\deg_ZP\geq2,\quad
\deg_WP\geq2.}                                      \tag{1.2}
\]

The last two inequalities are in circular coordinates.  They explain
precisely why collapsing Long's third real variable into the existing
complex pair fails: the resulting polynomial becomes affine in one Wick
source or reduces to two rotational weights.

There is one further degree-sensitive constraint:

\[
 \deg P=3\quad\Longrightarrow\quad
 |\operatorname{wt}(P)|\geq4                       \tag{1.3}
\]

for any counterexample.

If \(\deg P=3\), then either the support has at least five weights or its
four-weight support is one of

\[
\boxed{
\{-3,-1,1,3\},\quad
\{-2,-1,0,1\},\quad
\{-2,-1,1,2\},\quad
\{-1,0,1,2\}.}                                    \tag{1.4}
\]

Theorem 7.2 eliminates the first set in (1.4), leaving the intermediate
list

\[
\boxed{
\{-2,-1,0,1\},\quad
\{-2,-1,1,2\},\quad
\{-1,0,1,2\}.}                                    \tag{1.5}
\]

Theorem 7.3 eliminates all three sets in (1.5).  Consequently

\[
\boxed{\deg P=3\quad\Longrightarrow\quad
|\operatorname{wt}(P)|\geq5}                     \tag{1.6}
\]

for every possible GMC(2) counterexample.

The exact regression and bounded support search are in
[`verify_two_real_gmc_frontier.py`](../scripts/verify_two_real_gmc_frontier.py).
Its machine-readable record is
[`two_real_gmc_frontier.json`](../artifacts/generated-results/two_real_gmc_frontier.json).
The heavier exact quotient-algebra certificate for item 6 is
[`verify_two_real_gmc_symmetric_chart.py`](../scripts/verify_two_real_gmc_symmetric_chart.py),
with record
[`two_real_gmc_symmetric_chart.json`](../artifacts/generated-results/two_real_gmc_symmetric_chart.json).
The seven remaining exact chart certificates are in
[`verify_two_real_gmc_remaining_four_weight.py`](../scripts/verify_two_real_gmc_remaining_four_weight.py),
with record
[`two_real_gmc_remaining_four_weight.json`](../artifacts/generated-results/two_real_gmc_remaining_four_weight.json).

## 2. A rational-exponential lemma

We use one elementary lemma twice.

### Lemma 2.1

Let \(r(t)\in\mathbb C(t)\) be regular at zero with \(r(0)=0\).  If the
formal germ \(\exp(r(t))\) is a rational function, then \(r=0\).

### Proof

Write \(\exp(r)=S\in\mathbb C(t)^\times\) and differentiate:

\[
 r'=\frac{S'}S.                                      \tag{2.1}
\]

At a finite pole of \(r\) of order \(e\), the left side has a pole of order
\(e+1\), whereas a rational logarithmic derivative has only simple poles.
Thus \(r\) has no finite poles.  If the resulting polynomial \(r\) is
nonconstant, then \(r'\) is polynomial and nonzero, while \(S'/S=O(t^{-1})\)
at infinity.  Hence \(r\) is constant, and \(r(0)=0\) gives \(r=0\).
\(\square\)

## 3. GMC for all quadratic Gaussian polynomials

This statement is dimension independent.

### Theorem 3.1

Let \(G\) be a standard real Gaussian vector in \(\mathbb R^n\), and let
\(P\in\mathbb C[x_1,\ldots,x_n]\) have total degree at most two.  If

\[
 \mathbb E(P(G)^m)=0\qquad(m\geq1),                 \tag{3.1}
\]

then, for every polynomial \(Q\),

\[
 \mathbb E(Q(G)P(G)^m)=0
\]

for all sufficiently large \(m\).

### Proof

Write

\[
 P(x)=x^TAx+b^Tx+c
\]

with \(A=A^T\).  Put \(R(t)=(I-2tA)^{-1}\).  Formal Gaussian completion of
the square gives

\[
 M(t):=\mathbb E(e^{tP(G)})
 =\frac{\exp\!\left(tc+\frac{t^2}{2}b^TR(t)b\right)}
        {\sqrt{\det(I-2tA)}}.                       \tag{3.2}
\]

Hypothesis (3.1) is exactly \(M(t)=1\).  Squaring (3.2) gives

\[
 \det(I-2tA)
 =\exp\!\left(2tc+t^2b^TR(t)b\right).              \tag{3.3}
\]

The exponent is rational and vanishes at zero.  Lemma 2.1 makes it zero and
therefore

\[
 \det(I-2tA)=1.                                    \tag{3.4}
\]

Thus \(A\) is nilpotent, so \(R(t)\) is a polynomial matrix in \(t\).

Introduce a source \(s\).  Dividing the completed-square identity with source
by (3.2) yields

\[
 \mathbb E\!\left(e^{s^TG+tP(G)}\right)
 =\exp\!\left(\frac12s^TR(t)s+t\,s^TR(t)b\right).  \tag{3.5}
\]

Apply the finite-order differential operator \(Q(\partial_s)\) and set
\(s=0\).  Because \(R(t)\) is polynomial, the result is a polynomial in
\(t\).  But its coefficient of \(t^m/m!\) is
\(\mathbb E(Q(G)P(G)^m)\).  All sufficiently high coefficients vanish.
\(\square\)

In particular the least possible total degree of any Gaussian counterexample
is three.  Long's four-real-variable cubic attains that bound, so the global
minimum counterexample degree is exactly three even though the minimum
dimension is still either two or three.

## 4. Rotational-weight obstruction

Give \(Z\) weight \(1\) and \(W\) weight \(-1\).  Every polynomial has a
unique finite decomposition

\[
 P=\sum_{k\in\mathbb Z}P_k,\qquad
 P_k=
 \begin{cases}
 Z^kA_k(U),&k\geq0,\\
 W^{-k}A_k(U),&k<0.
 \end{cases}                                       \tag{4.1}
\]

Expectation kills every nonzero weight and restricts on
\(\mathbb C[U]\) to the factorial functional

\[
 \mathcal L(U^j)=j!.                               \tag{4.2}
\]

### Theorem 4.1

If \(P\in\mathbb C[Z,W]\) has at most two rotational weights and all its pure
moments vanish, then \(P\) has one-sided nonzero weight support.  In
particular \(P\) satisfies GMC.

### Proof

Same-sign nonzero weights already give the conclusion by weight: multiplying
by a fixed \(Q\) cannot restore weight zero for arbitrarily large powers.
If one weight is zero and the other is nonzero, the zero-weight part of
\(P^m\) is \(A_0(U)^m\).  The one-variable Factorial Theorem forces
\(A_0=0\).

It remains to take weights \(a>0\) and \(-b<0\).  Write

\[
 a=ga_0,\qquad b=gb_0,\qquad \gcd(a_0,b_0)=1
\]

and

\[
 P=Z^aA(U)+W^bB(U).
\]

A product of \(m\) factors has weight zero only when

\[
 m=(a_0+b_0)r
\]

and exactly \(b_0r\) factors come from the first summand.  Hence

\[
\mathbb E(P^{(a_0+b_0)r})
=\binom{(a_0+b_0)r}{b_0r}
 \mathcal L\!\left[
 \left(U^{ga_0b_0}A(U)^{b_0}B(U)^{a_0}\right)^r
 \right].                                          \tag{4.3}
\]

All other pure moments vanish by weight.  The Factorial Theorem applied to
the bracketed polynomial gives \(A=0\) or \(B=0\).  The support is therefore
one-sided, and fixed multipliers can cancel its growing weight only finitely
often. \(\square\)

This extends the weights \(1,-1\) calculation in Long's dimension-two
remark to arbitrary opposite rotational weights and arbitrary radial
polynomial coefficients.

## 5. Complete obstruction for one affine Wick source

### Theorem 5.1

Let

\[
 P=W h(Z)+v(Z),\qquad h,v\in\mathbb C[Z].           \tag{5.1}
\]

If every pure moment of \(P\) vanishes, then either

1. \(h\) is constant and \(v=0\), so \(P\) is a scalar multiple of \(W\); or
2. \(h\) is divisible by \(Z^2\) and \(v\) by \(Z\), so every monomial of
   \(P\) has strictly positive weight.

In both cases \(P\) satisfies GMC.

### Proof

Let \(g=t h(g)\).  Circular contraction followed by univariate Lagrange
inversion gives, coefficientwise in \(t\),

\[
 \mathbb E\!\left(A(Z)e^{tP}\right)
 =\frac{A(g)\exp(t v(g))}{1-t h'(g)}.              \tag{5.2}
\]

First suppose \(h(0)\ne0\).  Then \(g\in t\mathbb C[[t]]\) has inverse
\(t=z/h(z)\).  Taking \(A=1\) in (5.2), pure-moment vanishing is equivalent
to

\[
 \exp\!\left(\frac{zv(z)}{h(z)}\right)
 =\frac{h(z)-zh'(z)}{h(z)}.                        \tag{5.3}
\]

Lemma 2.1 forces both rational functions in the exponent to be zero and the
right side to be one.  Thus \(v=0\) and \(h'=0\).

If \(h(0)=0\), the unique solution of \(g=t h(g)\) is \(g=0\).  Equation
(5.2) becomes

\[
 \mathbb E(e^{tP})
 =\frac{\exp(t v(0))}{1-t h'(0)}.                 \tag{5.4}
\]

Lemma 2.1, or direct coefficient comparison, gives
\(v(0)=h'(0)=0\).  Therefore \(Z\mid v\) and \(Z^2\mid h\); every term in
(5.1) has positive weight. \(\square\)

Swapping \(Z\) and \(W\) gives the same theorem for polynomials affine in
\(Z\).  A two-real-variable witness must consequently have degree at least
two in each circular coordinate.

## 6. The cubic three-weight stratum is empty

### Theorem 6.1

Let \(P\in\mathbb C[Z,W]\) have total degree at most three and exactly three
nonzero rotational-weight components.  If its support contains both a
positive and a negative weight, then

\[
 \mathbb E(P^m)\ne0
\]

for some \(m\leq8\).

Consequently every cubic polynomial with at most three rotational weights
and vanishing pure moments has one-sided nonzero weight support and satisfies
GMC.

### Proof

The possible weights are

\[
 -3,-2,-1,0,1,2,3.
\]

For weight \(k\), the complete degree-three component is

\[
 P_k=
 \begin{cases}
 Z^k\displaystyle\sum_{j=0}^{\lfloor(3-k)/2\rfloor}a_{k,j}U^j,
     & k\geq0,\\[6pt]
 W^{-k}\displaystyle\sum_{j=0}^{\lfloor(3+k)/2\rfloor}a_{k,j}U^j,
     & k<0.
 \end{cases}                                       \tag{6.1}
\]

There are 27 three-element supports containing both signs.  Cover the locus
where all three components are nonzero by choosing one nonzero coefficient
in each component.  This gives 72 affine charts.

On a chart with support \(k_1<k_2<k_3\), use

\[
 P(Z,W)\longmapsto \lambda P(rZ,r^{-1}W).          \tag{6.2}
\]

Because \(k_1\ne k_2\), suitable
\(\lambda,r\in\mathbb C^\times\) set the selected coefficients in the first
two components to one.  Transformation (6.2) preserves pure-moment
vanishing: scalar multiplication is harmless, and every surviving
contraction has total rotational weight zero.  If \(q\) is the selected
coefficient in the third component, adjoin \(\rho q-1\) to keep that
component nonzero.

For each chart form the exact rational ideal

\[
 I_M=\left\langle
 \rho q-1,\mathbb E(P),\ldots,\mathbb E(P^M)
 \right\rangle .                                  \tag{6.3}
\]

Using (1.1), reduced grevlex Groebner calculation gives \(I_M=(1)\) on every
chart for some \(M\leq8\).  Thus the mixed-sign locus is empty.  If the
support does not contain both signs, its zero-weight component—if
present—vanishes by the one-variable Factorial Theorem, and the remaining
support is one-sided. \(\square\)

The verifier constructs all 27 supports from (6.1), checks all 72
nonvanishing charts, and records the moment-order histogram for each support.
The worst charts are exactly certified at order eight.

## 7. Almost every cubic four-weight support is empty

### Theorem 7.1

Let \(P\in\mathbb C[Z,W]\) have total degree at most three and exactly four
nonzero rotational-weight components containing both signs.  If all pure
moments of \(P\) vanish, then its support is one of the four sets in (1.4).

### Proof

There are 33 mixed-sign four-element subsets of
\(\{-3,-2,-1,0,1,2,3\}\).  Formula (6.1) gives 121 nonvanishing coefficient
charts.  Swapping \(Z\) and \(W\) negates the support and reduces the 33
supports to 18 sign-symmetry orbits.

On every chart use the same two-parameter normalization (6.2) on selected
coefficients in the first two components.  If \(q_3,q_4\) are selected
nonzero coefficients in the remaining components, localize with
\(\rho q_3q_4-1\).  The exact ideals

\[
 \left\langle
 \rho q_3q_4-1,\mathbb E(P),\ldots,\mathbb E(P^M)
 \right\rangle                                    \tag{7.1}
\]

have reduced grevlex basis \([1]\), for some \(M\leq6\), on all charts of 15
of the 18 sign orbits.  Reflecting those certificates excludes 29 supports
and 97 charts.  The three unexcluded sign orbits give exactly the four
supports in (1.4), with 24 charts. \(\square\)

This theorem does not infer nonexistence from a timeout.  The verifier
constructs all 33 supports, explicitly skips the three named exceptional
orbits, and requires a unit basis on every chart of every other orbit.

### Theorem 7.2

There is no pure-moment-zero polynomial with rotational-weight support
\(\{-3,-1,1,3\}\).

### Proof

Normalize the chart as

\[
 P=W^3+W+aZW^2+cZ+dZ^2W+eZ^3,\qquad ce\ne0.       \tag{7.2}
\]

All odd moments vanish by parity.  The second moment is

\[
 2(2ac+6ad+c+2d+6e),
\]

so

\[
 e=-\frac{ac}{3}-ad-\frac c6-\frac d3.            \tag{7.3}
\]

Let \(\widetilde M_j\in\mathbb Z[a,c,d]\) be the primitive numerator of
the \(j\)-th moment after (7.3), and put

\[
 q=c(2ac+6ad+c+2d)=-6ce.
\]

The localized order-eight ideal

\[
 I_8=\langle\widetilde M_4,\widetilde M_6,\widetilde M_8,
                  \rho q-1\rangle
       \subset\mathbb Q[a,c,d,\rho]                \tag{7.4}
\]

has an exact reduced grevlex basis with 53 elements.  Its quotient is
zero-dimensional of vector-space dimension 84.

It remains to decide whether \(\widetilde M_{10}\) vanishes anywhere on
this finite scheme.  Reduce the rational basis modulo
\(p=1\,000\,003\), and independently compute the reduced basis of (7.4)
over \(\mathbb F_p\).  The two bases have the same 53 leading monomials,
and every reduced rational generator has zero remainder modulo the
independent modular basis.  Hence this is a good reduction of (7.4), with
the same 84 standard monomials in characteristic zero and characteristic
\(p\).

In that standard-monomial basis, the exact \(84\times84\) matrix of
multiplication by \(\widetilde M_{10}\) has rank 84 over \(\mathbb F_p\).
Its determinant is therefore nonzero modulo a good prime and consequently
nonzero over \(\mathbb Q\).  Multiplication by
\(\widetilde M_{10}\) is an automorphism of the rational quotient, so
\(\widetilde M_{10}\) is a unit there.  Thus (7.2) cannot satisfy the even
moments through order ten.

For the full support, write

\[
 P=W^3+bW+aZW^2+cZ+dZ^2W+eZ^3.                   \tag{7.5}
\]

The four nonvanishing charts select one of \(b,a\) and one of \(c,d\).
Besides (7.2), it is enough to check the constant/radial and radial/radial
representatives: swapping \(Z,W\) exchanges the constant/radial and
radial/constant charts.  On each representative, eliminate \(e\) with the
second moment and localize at the selected positive-weight coefficient
times \(e\).  The exact results are identical:

| representative | charts covered | rational basis | quotient length | \(M_{10}\) rank mod \(p\) |
|---|---:|---:|---:|---:|
| constant/constant | 1 | 53 | 84 | 84 |
| constant/radial | 2 | 53 | 84 | 84 |
| radial/radial | 1 | 53 | 84 | 84 |

For every row, the independently computed modular basis has the same
leading monomials as the reduction of the rational basis and zero comparison
remainders.  Hence the determinant-lifting argument above applies to all
four charts. \(\square\)

The good-prime comparison is part of the certificate; a modular unit-ideal
calculation by itself would not justify the characteristic-zero conclusion.

### Theorem 7.3

There is no pure-moment-zero cubic with exactly four nonzero rotational
weights containing both signs.

### Proof

Theorem 7.1 reduces the problem to (1.4), and Theorem 7.2 removes its
symmetric support.  It remains to exclude \(E_2,-E_2,E_3\) from (1.5).

For

\[
E_2=\{-2,-1,0,1\},
\]

the first moment centers the zero-weight component:

\[
g+hU\longmapsto h(U-1).
\]

Thus the two coefficient choices in weight zero define the same
nonvanishing locus.  Normalize

\[
 P=W^2+bW+aZW^2+h(U-1)+cZ+dZ^2W                 \tag{7.6}
\]

by setting either \(b=1\) or \(a=1\), and select either \(c\) or \(d\) in
the positive component.  On the resulting four representatives, localize
at \(hc\) or \(hd\).  In every case the exact ideal

\[
\langle \mathbb E(P^2),\ldots,\mathbb E(P^6),
        \rho h(c\text{ or }d)-1\rangle            \tag{7.7}
\]

has reduced standard basis \([1]\) over \(\mathbb Q\).  Swapping \(Z,W\)
also excludes \(-E_2\).  These four calculations cover the 16 original
charts on the reflected pair.

For the balanced support

\[
E_3=\{-2,-1,1,2\},
\]

write

\[
 P=W^2+bW+aZW^2+cZ+dZ^2W+eZ^2.                  \tag{7.8}
\]

The second moment is linear in \(e\), with nonzero constant coefficient, so
eliminate \(e\).  Constant/constant, constant/radial, and radial/radial
representatives cover all four charts up to \(Z,W\) reflection.  After
localizing at the selected positive coefficient times \(e\), each exact
ideal generated by moments through order six again has reduced rational
basis \([1]\).

Hence the last three supports and 20 charts are empty.  Together with
Theorem 7.2, all 33 mixed-sign four-weight supports and all 121
nonvanishing charts are excluded. \(\square\)

## 8. Direct Long-style collapse is already inconsistent

The most literal attempt to absorb the third real Gaussian in Long's
three-variable witness is

\[
 P=W(1+Z)+W^2k(Z).                                 \tag{8.1}
\]

For

\[
 k(Z)=\sum_{j=0}^d c_jZ^j,
\]

the exact rational moment ideals have unit Groebner basis in the following
bounded cases:

| \(d\) | coefficients | moments imposed | reduced basis |
|---:|---:|---:|---|
| 2 | 3 | \(m\leq2\) | \(1\) |
| 3 | 4 | \(m\leq5\) | \(1\) |
| 4 | 5 | \(m\leq6\) | \(1\) |

These are finite-support nonexistence results: an all-order witness in one of
the displayed supports would solve each finite system, so none exists.
They do not exclude unbounded \(Z\)-degree or a more general polynomial
quadratic in \(W\).

## 9. Updated frontier

The two-real problem is no longer an undifferentiated coefficient search.
Any useful next ansatz must begin beyond all three proven walls:

1. at least cubic total degree;
2. at least three rotational weights containing both signs;
3. nonlinear dependence on both \(Z\) and \(W\).

If the witness is cubic, item 2 strengthens to at least five weights.  The
entire four-weight stratum is now empty: Theorem 7.2 excludes the symmetric
support through moment ten, and Theorem 7.3 excludes the other three
supports through moment six.

The five-weight census now excludes 18 of the 21 supports and 92 of the 102
ordinary charts exactly.  All 35 centered/reflected chart presentations
are modular unit ideals through moment eight; only three chart orbits
(ten ordinary charts) still lack characteristic-zero promotion.  They are
the constant/constant chart on \(\{-2,-1,0,1,2\}\) and two charts on
\(\{-3,-1,0,1,2\}\), together with the latter support's reflection.
Thus the smallest unexplored structural classes are these three cubic
five-weight chart orbits and a three-weight quartic.  Their moment
generating functions are not one-coordinate Lagrange determinants.
Progress there requires either finite quotient-algebra elimination or a
new identity for shared Wick self-contractions, exactly the collapse not
supplied by the four-real fixed-point architecture.

The executable next attacks, including the 154 five-to-seven-weight cubic
charts and the cubic invariant-null-cone formulation, are in
[`TWO_REAL_GMC_CLOSURE_PROGRAM.md`](TWO_REAL_GMC_CLOSURE_PROGRAM.md).

## Reproduction

```bash
.venv/bin/python scripts/verify_two_real_gmc_frontier.py
.venv/bin/python scripts/verify_two_real_gmc_symmetric_chart.py
.venv/bin/python scripts/verify_two_real_gmc_remaining_four_weight.py
.venv/bin/python scripts/verify_two_real_gmc_five_weight.py --discovery-only
.venv/bin/python scripts/verify_two_real_gmc_five_weight_frontier.py
```

The script checks (4.3), checks (5.2) through order seven for a nonlinear
test, checks the quadratic tilted-Gaussian identity, proves the cubic
three-weight theorem on all 72 nonvanishing charts, excludes 29 of the 33
mixed-sign four-weight supports on 97 charts, and recomputes the three
Long-style unit Groebner bases over \(\mathbb Q\).  The second command
computes three 53-element rational order-eight bases covering all four
charts of the symmetric support, proves good reduction at
\(p=1\,000\,003\), and verifies full rank of tenth-moment multiplication on
each 84-dimensional modular quotient.  The third command verifies seven
exact rational unit ideals covering the remaining three supports and 20
charts through moment order six.
