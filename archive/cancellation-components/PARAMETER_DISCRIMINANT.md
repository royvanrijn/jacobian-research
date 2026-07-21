# Closed discriminant of the C24 parameter polynomial

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

The exact regression verifies (2), (4), the resultant formula, and the square
criterion against direct discriminant calculations.  It also checks the
converse parametrization (5)--(6) and the infinite family (7) in bounded
ranges:

```bash
.venv/bin/python scripts/verify_parameter_discriminant.py
```
