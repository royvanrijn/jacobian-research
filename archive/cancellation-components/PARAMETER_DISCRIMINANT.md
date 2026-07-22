# Closed discriminant of the cancellation parameter polynomial

Put `n=mr`, `L=n+r+1`, and

\[
 \mathcal M_{m,r}(q)=
 \sum_{j=0}^{n}(-1)^j\binom Ljq^{n-j}.
\]

The discriminant of this polynomial has a uniform closed form.

## Theorem

For every `m,r>=1`,

\[
 \operatorname{disc}(\mathcal M_{m,r})
 =(-1)^{n(n-1)/2}(r+1)L^{n-1}
   \binom{n+r}{r}^{n-2}.                              \tag{1}
\]

Consequently it is nonzero, so every parameter polynomial is separable over
`Q`, independently of whether its irreducibility is known.

## Proof

Set

\[
 B_{n,r}(x)=\sum_{j=0}^n\binom{j+r}{r}x^j,
 \qquad a=\binom{n+r}{r}.
\]

The binomial transform gives

\[
 x^n\mathcal M_{m,r}(1-x^{-1})=(-1)^nB_{n,r}(x).      \tag{2}
\]

The fractional-linear root change in (2) preserves the discriminant: its
root-difference denominators multiply to `M_(m,r)(1)^(2n-2)=a^(2n-2)`,
exactly cancelling the leading-coefficient factor of `B`.  The sign in (2)
also has no effect because a scalar multiplies a degree-`n` discriminant to
the power `2n-2`.  Hence

\[
 \operatorname{disc}(\mathcal M_{m,r})
 =\operatorname{disc}(B_{n,r}).                       \tag{3}
\]

Consecutive binomial coefficients give the differential identity

\[
 (1-x)B'_{n,r}(x)-(r+1)B_{n,r}(x)=-La x^n.            \tag{4}
\]

If `alpha_1,...,alpha_n` are the roots of `B`, then

\[
 B'(\alpha_i)={-La\alpha_i^n\over1-\alpha_i}.
\]

Moreover,

\[
 \prod_i\alpha_i={(-1)^n\over a},\qquad
 \prod_i(1-\alpha_i)={B(1)\over a}
 ={\binom{n+r+1}{r+1}\over\binom{n+r}{r}}
 ={L\over r+1}.
\]

Multiplying (4) over the roots therefore yields

\[
 \prod_iB'(\alpha_i)=(r+1)L^{n-1}.
\]

Thus

\[
 \operatorname{Res}(B,B')=(r+1)L^{n-1}a^{n-1}.
\]

Using
`disc(B)=(-1)^(n(n-1)/2) Res(B,B')/a` proves (1).  QED

## Rational-square criterion

Formula (1) immediately gives an exact test:

- if `n` is congruent to `2` or `3` modulo `4`, the discriminant is negative
  and is not a square in `Q`;
- if `n` is congruent to `0` modulo `4`, it is a square exactly when
  `(r+1)(n+r+1)` is an integer square;
- if `n` is congruent to `1` modulo `4`, it is a square exactly when
  `(r+1)binom(n+r,r)` is an integer square.

This criterion explains all alternating groups in the exact
[degree-thirty Galois table](PARAMETER_GALOIS_GROUPS.md).  Through `mr<=30`,
the nontrivial square-discriminant pairs are

\[
 (4,3),(2,8),(16,1),(17,1),(1,24),(12,2).
\]

## Complete even-degree square locus

The even-degree part of the square criterion has a simple complete
parametrization.  Write

\[
 r+1=da^2
\]

with `d` squarefree.  If `n=mr` is divisible by four, then (1) is a rational
square if and only if

\[
 n+r+1=db^2                                      \tag{5}
\]

for an integer `b>a`.  Equivalently,

\[
 m={d(b^2-a^2)\over r},                           \tag{6}
\]

where the right side is integral and `d(b^2-a^2)` is divisible by four.
Indeed, two positive integers have square product exactly when they have the
same squarefree part.  Applying this to `(r+1)(n+r+1)` proves both directions.

In particular, square discriminants occur infinitely often for every fixed
`r`.  For any `k>=1`, take

\[
 b=a+2rk,\qquad m=4dk(a+rk).                       \tag{7}
\]

Then `mr=d(b^2-a^2)` is divisible by four and (5) holds.  Thus the parameter
Galois group is contained in `A_(mr)` whenever the corresponding polynomial
is irreducible.  This shows that the alternating cases in the finite table
belong to an infinite arithmetic phenomenon, not a sporadic low-degree list.

For fixed `r`, (6) can also be viewed as a finite union of quadratic
progressions.  Since `gcd(d,r)=1`, integrality is equivalent to

\[
 b^2\equiv a^2\pmod r,
\]

and the condition `4|d(b^2-a^2)` is a further congruence modulo four.  Thus
the admissible `b` occupy finitely many residue classes modulo `4r`, and (6)
maps each class to a quadratic sequence of square-discriminant parameters.

## An odd-degree square family

The odd part of the square locus is governed by the binomial coefficient in
the criterion and need not admit the same two-square parametrization.  One
infinite family is nevertheless immediate.  Set `r=1`.  For `n=m=1 mod 4`,
the criterion becomes

\[
 2\binom{m+1}{1}=2(m+1)\text{ is a square}.
\]

It follows that

\[
 m=2a^2-1,\qquad r=1,\qquad a\text{ odd}               \tag{8}
\]

always gives a rational-square discriminant.  Conversely, every square-
discriminant pair with `r=1` and odd degree has this form: positivity forces
`m=1 mod 4`, and writing `2(m+1)=c^2` forces `c=2a` and hence (8), with `a`
odd.  The first values are

\[
 m=1,17,49,97,161,241,337,449,\ldots.
\]

Thus both parities contain infinite alternating-group candidates.
The exact Frobenius--Jordan regression proves that the next two cases after
`m=17` are not merely candidates:

\[
 \operatorname{Gal}(\mathcal M_{49,1})=A_{49},\qquad
 \operatorname{Gal}(\mathcal M_{97,1})=A_{97}.          \tag{9}
\]

For degree 49 an isolated 13-cycle rules out the only possible nontrivial
block size 7; for prime degree 97, irreducibility already implies
primitivity.  Jordan's theorem and the square discriminants then give (9).

## Fixed-row density of the square locus

The two parity analyses give a useful quantitative simplification.  For
fixed `r`, the even-degree solutions are a finite union of the quadratic
sequences (6), so there are `O_r(sqrt(X))` of them with `m<=X`.

For the odd-degree branch, `r` and `m` must both be odd and the criterion is

\[
 Y^2=(r+1)\binom{r(m+1)}r
 =\frac{r+1}{r!}\prod_{i=0}^{r-1}(r(m+1)-i).           \tag{10}
\]

As a polynomial in `m`, the right side has degree `r` and `r` distinct
roots.  For each fixed odd `r>=3`, (10) is therefore a squarefree
hyperelliptic curve of genus `(r-1)/2>=1`.  Siegel's theorem makes its
integral points finite.  This is the same integral-point mechanism used in
[Filaseta--Moy](https://doi.org/10.4064/cm7474-3-2018), Lemma 4, for their
fixed-degree square-discriminant analysis.  When `r=1`, the curve has genus
zero and its integral square values are exactly the family (8).

Consequently, for every fixed `r`,

\[
 \#\{m\le X:\operatorname{disc}(\mathcal M_{m,r})
             \text{ is a square}\}=O_r(\sqrt X).       \tag{11}
\]

In particular, alternating containment has density zero on each fixed-`r`
row even though it occurs infinitely often there.

The exact regression verifies the geometric-derivative identity, (2), (4),
the resultant formula, and the square criterion against direct discriminant
calculations.  It also checks the converse parametrization (5)--(6), the
infinite even family (7), and the odd family (8), including its converse for
`r=1`, as well as the degree and squarefreeness of the curve (10), in bounded
ranges:

```bash
.venv/bin/python scripts/verify_parameter_discriminant.py
```
