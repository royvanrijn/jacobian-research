# Explicit corner cutoffs from Gaussian exclusion of lifts

The [four-column lifting theorem](SFC_HOMOGENEOUS_CENTERED_HIGHER_LIFT.md#6-a-four-column-family-with-nonliftable-gaussian-zeros)
has an explicit numerical cutoff. At scale one a sufficient value is
\[
 H(1)=2^{373029362779174993920000137}.                 \tag{1}
\]
At least one factorial moment of orders \(1,2,3,5\) is nonzero for
the translated offsets \((0,9),(3,6),(4,5),(9,0)\) whenever both
common integer heights exceed (1). This also holds for real gamma
heights. The bound is deliberately enormous. It does not enumerate
exceptional heights or settle the unrestricted odd-start Gaussian problem.

## 1. Quantitative transfer and its hypothesis

Fix \(s\ge2\) distinct offsets \((\lambda_j,D-\lambda_j)\), with
\(0\le\lambda_j\le D\), and \(s\) positive orders of maximum \(M\).
Use \(\sigma^2=L_x^{-1}+L_y^{-1}\) and \(r=L_x^{-1}/\sigma^2\)
in the exact Gaussian normalization. Endpoint values \(r=0,1\)
mean polynomial continuation, not finite positive heights.

Assume **compactified corner exclusion**: the exact polynomial
system has no sequence of nonzero common coefficient directions with
\(\sigma\downarrow0\) and \(r\in[0,1]\). This is a qualitative proof
input; the bound below does not establish it. Gaussian zero-freeness
is sufficient. The four-column family satisfies it through its
lifting obstructions despite its two Gaussian zeros.

For positive integers write \(\operatorname{bit}(u)=\lfloor\log_2u\rfloor+1\).
Define
\[
 B=(s-1)!D^{(s-1)^2},\quad
 L=(2s^2B)^M\bigl(M(1+2D)\bigr)^{MD},\quad
 \tau=1+\operatorname{bit}(L),
\]
\[
 n=2s,\quad b=3s+2,\quad d=M(3D+1),\quad
 N=(2d+1)(2d)^{n-1},\quad \Delta=n(2d-1)+2,
\]
\[
 A=2\tau+\operatorname{bit}(N)+(n+1)\operatorname{bit}(d+1)
       +\operatorname{bit}(b)+3,
\]
\[
 E=4N\Delta A+2\operatorname{bit}(N+1)+\operatorname{bit}(2n).
                                                               \tag{2}
\]

**Explicit cutoff theorem.** Under compactified corner exclusion,
\(H=2^E\) is a sufficient common-height cutoff. At integer heights
the conclusion holds over every characteristic-zero field. Formula
(2) has no unevaluated constants or quantifier-elimination search.

## 2. Integer normalization and size bounds

Set \(V_{ij}=\binom{\lambda_j}{i}\), \(0\le i<s\), and
\(W=\operatorname{adj}(V)\). Distinct exponents give \(\det V\ne0\).
Each cofactor is a sum of \((s-1)!\) products of \(s-1\) entries,
each at most \(D^{s-1}\). Hence \(\lvert W_{ji}\rvert\le B\).
Put \(\ell=s-1\), \(C_j=\sum_{i=0}^{\ell}W_{ji}a_i\sigma^{\ell-i}\).
An integer multiple of the exact moment normalization is
\[
 Q_m=\sigma^{-\ell m}\sum_{|\alpha|=m}\binom m\alpha C^\alpha
 \prod_{\nu=1}^{U}(m+\nu\sigma^2r)
 \prod_{\nu=1}^{mD-U}(m+\nu\sigma^2(1-r)),\qquad
 U=\sum_j\lambda_j\alpha_j.                                  \tag{3}
\]
All numerator powers below \(\sigma^{\ell m}\) cancel. To justify
division at every support and order, use the confluent binomial
basis in the [exact graded normalization](SFC_BIVARIATE_TRANSLATION_GAUSSIAN_JETS.md#1-exact-graded-normalization).
After dividing its \(i\)-th member by \(\sigma^i\), the phase is
polynomial in \(\sigma\) and centered gamma coordinates. Their
moments are polynomial in \(\sigma,r\) by the finite cumulant
formula. Multiplication by \(m^{mD}(\det V)^m\) gives (3).

The separate degrees in \(\sigma,r,a\) are at most \(2mD,mD,m\).
Substitute \(\sigma=1/w\) and multiply by \(w^{2mD}\). The total
degree is at most \(m(3D+1)\le d\).

The sum of the coefficient one-norms of the \(C_j\)'s is at most
\(s^2B\), and each rising-factorial factor has one-norm at most
\(m(1+2D)\). Thus (3) has one-norm at most
\((s^2B)^m(m(1+2D))^{mD}\). Cancellation and division by a monomial
do not increase it. Splitting coefficients into real and imaginary
parts costs at most \(2^m\). The bound \(L\) in (2) therefore
bounds every integer coefficient in absolute value.

## 3. The bounded set and the imported radius theorem

Cover projective coefficient space by \(s\) charts: choose a coefficient
of maximum modulus and divide by it. It becomes \(1\); the other
coefficients satisfy \(\lvert a_i\rvert\le1\). Each chart has
\(2(s-1)\) real coefficient variables, plus \(w,r\), giving \(n=2s\).

Impose the \(2s\) real and imaginary moment equations, the \(s-1\)
closed unit-disk inequalities, and \(w-1\ge0,r\ge0,1-r\ge0\).
This is one weak sign condition involving \(b=3s+2\) polynomials
of degree at most \(d\) and coefficient bit size at most \(\tau\).
It defines a closed set. Coefficients and \(r\) are bounded, and
compactified exclusion rules out \(w\to\infty\). The whole set
is therefore bounded, including its boundary at \(r=0,1\).

[Basu–Roy, Theorem 3, author version of 5 June 2010, p. 5](https://www.math.purdue.edu/~sbasu/jsc_final-06-05-10.pdf),
retained as [the exact primary PDF](sources/basu-roy-explicit-radius-2010.pdf),
now supplies the containing radius
\[
 R=\sqrt n\,(N+1)2^{2N\Delta A}.                            \tag{4}
\]
Our \(A\) has one extra \(\operatorname{bit}(d+1)\), safely
enlarging the printed theorem statement. Section 5 of the source,
pp. 13–14, was also inspected: extrema are reduced to sums of
squares of active equations and the algebraic radius bound.
That theorem and its real algebraic geometry inputs remain imported.

Every common zero with \(w\ge1\) has \(w\le R\). Since
\(w^2=L_xL_y/(L_x+L_y)\ge\min(L_x,L_y)/2\), a failure must satisfy
\(\min(L_x,L_y)\le2R^2\le2^E\). This proves the theorem. At integer
heights the equations are rational, so projective Nullstellensatz
certificates transfer nonvanishing to all characteristic-zero fields.

## 4. Four columns, odd starts and replay

For \((0,9k),(3k,6k),(4k,5k),(9k,0)\), use \(s=4,D=9k,M=5\).
The fifth-jet theorem supplies compactified exclusion: every limiting
root is one of its two simple cubic zeros; the first projected
equation forces \(r\to1/2\), excluding boundary-only sequences.
Both real heights then diverge, and the nonreal next height correction
excludes the sequence. At \(k=1\), formula (2) gives
\[
 (n,b,d,\tau)=(8,14,140,478),\quad
 N=37915129118720000000,\quad\Delta=2234,\quad A=1101,
\]
and (1). Formula (2) gives explicit \(H(k)\) at every positive scale.

For a homogeneous bivariate trinomial of exponent span \(D\), first
absorb its common monomial into the two heights. A verified odd
Gaussian start \(q\) then gives (2) with \(s=3,M=q+2\).
Gaussian zero-freeness makes compactified exclusion immediate. This
applies to every odd start through \(12001\) and the infinite families
in the [odd-start note](SFC_ODD_GAUSSIAN_LOCAL_WINDOW_TESTS.md).
In particular, it applies at every \(q=(2k+1)p^e\), with
\(k\ge0,e\ge1\) and odd prime \(p>4k+2\), without a
character or resultant restriction. Thus every \(q=9\cdot89^e\)
inherits the same explicit formula, with its own \(M=q+2\).
Thus explicit bounds are now available throughout that proved domain.
The bound still depends on the support span and the chosen start;
no cutoff uniform over starts or supports is asserted.
It remains unproved that the joint prime-divisor tests succeed at
every odd start, or that every remaining Gaussian direction cannot lift.
The finite prefix is not an all-order theorem.

The [checker](verify_sfc_explicit_gaussian_cutoff.py) evaluates (2):

```sh
.venv/bin/python research/factorial/image-rigidity/verify_sfc_explicit_gaussian_cutoff.py --audit
```

The [receipt](sfc_explicit_gaussian_cutoff_verification.json) records
an internal pass. At scale one, it expands all 213, 1854, 7945 and
59374 retained terms at orders \(1,2,3,5\), checks exact monomial
divisibility and size bounds, and recovers the Gaussian constant
terms by a separate count-vector calculation. This checks arithmetic
and normalization, not the imported radius theorem or the complete
limiting argument. There is no independent replay, formal verification
or external review. A useful small cutoff and the actual exceptional
height set remain uncomputed.
