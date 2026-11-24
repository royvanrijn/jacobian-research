# Overlapping swap pairs and the limit of reciprocal stability under row merging

Internal research note, 21 November 2026. This supplies an all-order
overlapping two-monomial family and a restricted class of integer row
mergers. It does not prove either unrestricted preservation obligation in
[UNRESTRICTED_SFC_MECHANISMS.tex](UNRESTRICTED_SFC_MECHANISMS.tex), or unrestricted SFC. Earlier theorem sources and their pins are unchanged.

## 1. An overlapping family with a private-degree threshold

Work over any characteristic-zero field `K`. Let `b_i>=0` and `d_i>=1`
be integers, for `1<=i<=r`, and put

\[
 U=\prod_{i=1}^r X_i^{b_i+d_i}Y_i^{b_i},\qquad
 V=\prod_{i=1}^r X_i^{b_i}Y_i^{b_i+d_i},\qquad
 D_0=\sum_{i:b_i=0}d_i.
 \tag{1}
\]

All displayed variables are independent. Additional rows with equal
exponents in the two columns are allowed: for arbitrary integers
`h_l>=0`, multiply both `U,V` by `C(Z)=prod_l Z_l^(h_l)` on additional
variables. Those rows contribute the nonzero scalar
`prod_l(h_l n)!` to the order-`n` moment and impose no degree restriction.
Equivalently, arbitrary equal-exponent rows can be removed before (1).

**Theorem 1.** If `D_0>=2`, every nonzero

\[
                         F=C(Z)(cU+dV),\qquad c,d\in K,
\]

satisfies SFC(2) at every positive starting order. More strongly, the
ordinary moment polynomials for `U+tV` have simple unit-circle roots,
and consecutive orders strictly interlace. Exponents, overlaps, and
starting orders are unrestricted within the stated inequality.

The condition measures the total degree of the part exclusive to each
monomial. It permits arbitrarily large overlaps. For example,

\[
        cX^2 A^2B+dY^2 AB^2
\]

is included. Its shared factor `AB` is not independent of the remaining
terms, so the earlier rule for disjoint monomial multipliers does not
already prove this example.

### Proof

Remove `C`, whose moment factor is nonzero. The exact moment polynomial is

\[
 P_n(t)=\sum_{k=0}^n p_{n,k}t^k,\qquad
 p_{n,k}=\binom nk\prod_i
            (b_i n+d_i k)!\,(b_i n+d_i(n-k))!.
 \tag{2}
\]

Its coefficients are strictly positive rational numbers and palindromic.
Let `Q_n=sum_k t^k/p_(n,k)` be its coefficientwise reciprocal. Define

\[
 \delta_{n,k}=\frac{(k!)^2}{\prod_i(b_i n+d_i k)!}.
\]

Then the required binomial normalization is exactly

\[
 \frac{[t^k]Q_n}{\binom nk}
     =\frac{\delta_{n,k}\delta_{n,n-k}}{(n!)^2}.
 \tag{3}
\]

For each fixed `n`, the consecutive ratio of `delta_(n,k)` is the
restriction to integer `x=k` of

\[
 R_n(x)=\frac{(x+1)^2}
        {\prod_i\prod_{j=1}^{d_i}(b_i n+d_i x+j)},\qquad x\ge0.
 \tag{4}
\]

Every denominator factor is positive, and

\[
 \frac{R_n'(x)}{R_n(x)}
  =\frac2{x+1}
      -\sum_i\sum_{j=1}^{d_i}\frac{d_i}{b_i n+d_i x+j}
  \le\frac{2-D_0}{x+1}\le0.
 \tag{5}
\]

For `b_i=0`, each of its `d_i` summands is at least `1/(x+1)`.
The other summands are positive and can only strengthen the inequality.
Thus `delta_(n,k)` is log-concave for each `n`. Its reversal and their
pointwise product are log-concave, proving the hypothesis for (3).

Apply Mao--Wang--Zheng, *Interlacing on the Unit Circle via
Coefficientwise Reciprocals*, [arXiv:2609.34551v1](https://arxiv.org/abs/2609.34551v1),
Remark 3.4. It permits independently chosen positive palindromic polynomials
of degrees `n,n+1`, each with log-concave binomial-normalized coefficients.
Therefore the dependence of `delta_(n,k)` on `n` is harmless: no fixed
coefficient sequence across the two orders is a hypothesis. The theorem
and its Sections 2--3 proof were inspected in the canonical
[SFC_IDENTICAL_MONOMIAL_PAIR.md](SFC_IDENTICAL_MONOMIAL_PAIR.md); this note uses exactly that input.

Consequently `P_n,P_(n+1)` have no common complex root. Their rational
polynomial Bezout identity extends to every characteristic-zero field.
For `c!=0`, the actual moment is `c^n P_n(d/c)`, multiplied by the
nonzero factor from `C`. This proves the assertion. If one coefficient
is zero, all moments of the remaining monomial are nonzero. QED.

This is a new application in this investigation of a credited preprint,
not a new proof of its interlacing theorem or a literature-priority claim.

## 2. A finite-order refinement and an exact limitation

For every `x>=0`,

\[
 b_i n+d_i x+j\le (b_i n+d_i)(x+1),
\]

so the stronger sufficient test at order `n` is

\[
             W_n:=\sum_i\frac{d_i^2}{b_i n+d_i}\ge2.
 \tag{6}
\]

It implies `R_n'/R_n<=0` by the same calculation. Since `W_n` decreases
with `n`, the single condition `W_(n+1)>=2` proves consecutive
coprimality at orders `n,n+1`. This can apply when `D_0<2`; for example,
`r=1,b_1=1,d_1=4` satisfies it for all `n+1<=4`. This is a proved
parameter inequality, not a finite search assertion.

An exact, potentially sharper finite-order test is also available.
Write `eta_(n,k)=delta_(n,k)delta_(n,n-k)`. For `1<=k<n`,

\[
\begin{aligned}
 E_{n,k}:=\frac{\eta_{n,k}^2}{\eta_{n,k-1}\eta_{n,k+1}}
  &=\left(\frac{k(n-k)}{(k+1)(n-k+1)}\right)^2\\
  &\quad\cdot\prod_i\prod_{j=1}^{d_i}
     \frac{b_i n+d_i k+j}{b_i n+d_i(k-1)+j}
     \frac{b_i n+d_i(n-k)+j}{b_i n+d_i(n-k-1)+j}.
\end{aligned}
 \tag{7}
\]

The exact criterion is `E_(n,k)>=1` at every internal index, and
likewise at order `n+1`. Formula (7) is an explicit test of the imported
hypothesis, not a proof that it always holds.

The threshold in Theorem 1 is **necessary and sufficient for this
normalized-log-concavity condition to hold at every order**. Sufficiency
is (5). For necessity, use the first internal index in (7). Its limit is

\[
             \lim_{n\to\infty}E_{n,1}
                =\frac14\prod_{i:b_i=0}\binom{2d_i}{d_i}.
 \tag{8}
\]

Indeed the factors with `k=1` converge to one for `b_i>0`, while those
with `b_i=0` multiply to `(2d_i)!/(d_i!)^2`. Every factor from the
opposite endpoint tends to one, and the leading factor tends to `1/4`.
If `D_0=0`, the limit is `1/4`; if `D_0=1`, it is `1/2`. Both are
strictly less than one, so the coefficient test fails for all sufficiently
large orders. Thus (6) cannot yield a genuinely larger all-order
parameter class through this same criterion.

This is **not an SFC counterexample**. For instance `XY(cX+dY)` has
`D_0=0`, belongs to the already proved common-power monomial--linear
family, and satisfies SFC. At order two its moment coefficients are
`(48,72,48)`, and its normalized reciprocal coefficients are
`(1/48,1/144,1/48)`, giving `E_(2,1)=1/9<1`. Likewise the elementary
disjoint linear pair has `D_0=1` and is already known to satisfy SFC.

## 3. Actual integer row operations covered by the theorem

Rows in (1) occur as a mirrored pair

\[
                   (b+d,b),\qquad (b,b+d).
\]

The following operations stay within a proved class whenever the endpoint
retains `D_0>=2`:

1. Merge these two rows into the equal row `(2b+d,2b+d)`. This is one
   genuine integer row merger; its resulting common-monomial moment
   factor is nonzero and can be stripped.
2. Merge two first rows together and the corresponding second rows
   together. The resulting mirrored pair has parameters
   `(b+b',d+d')`. This is a matched pair of mergers. The theorem is
   asserted for the completed matched operation, not for an asymmetric
   intermediate matrix.
3. Adjoin `e_1,e_2` as a matched pair of unit rows. This raises `D_0`
   by one. It does not prove preservation for adjoining only `e_1`.

There is also a direct construction from disjoint identical monomials.
The four single-column rows
`(b+d,0),(b,0),(0,b+d),(0,b)` can be merged in crossed pairs to obtain
the two mirrored rows above. Reserve private degree at least two and
perform this construction for any number of desired overlaps.
All intermediate columns remain distinct because the reserved private
variables distinguish them. Only the symmetric endpoint is certified
by this argument.

To obtain arbitrary support sizes, append any number `l` of independent
linear terms to the proved two-term core. A window of `l+2` consecutive
moments is then nonzero. This is the existing linear-padding theorem of
Rocks, author preprint (2015), Theorem 3.12 and Corollary 3.13, inspected
at printed pages 12--13 and already reproduced in
[POWER_DIAGONAL_LINEAR_SFC_TRANSFER.md](POWER_DIAGONAL_LINEAR_SFC_TRANSFER.md), Section 4. Its short proof is:
if `f_m=Lambda(F^m)/m!` and `g_m` are the padded values, then

\[
       \sum g_m z^m=\frac{\sum f_m z^m}
                            {\prod_{j=1}^l(1-\lambda_jz)}.
\]

A block `g_n,...,g_(n+l+1)=0` forces `f_(n+l)=f_(n+l+1)=0` after
multiplication by the degree-`l` denominator. Boundary coefficients only
reduce the active term count. The row operations above may be confined
to the two core columns, leaving these linear rows untouched. This gives
uniform support-size and starting-order preservation within this class,
using the credited padding rule rather than a new addition theorem.

## 4. Reciprocal stability fails under one legal row merger

For an integer support `A`, let `P_(A,m)(c)` be its actual homogeneous
order-`m` moment form, and let `Q_(A,m)` denote its coefficientwise
reciprocal. All its coefficients are positive rationals.

For every `s>=2`, start with the distinct nonconstant features

\[
                        X_1Y_1,\ldots,X_sY_s.
\]

At multiindex `k` of total degree `m`, their moment coefficient is
`m! prod_j k_j!`. Hence

\[
                  Q_{A,m}(c)=\frac{(c_1+\cdots+c_s)^m}{(m!)^2},
 \tag{9}
\]

which is real stable for every `m`: a sum of upper-half-plane numbers
has positive imaginary part.

Merge just the rows for `X_1,X_2`. The new features are
`TY_1,TY_2,X_3Y_3,...,X_sY_s`; the surviving `Y_j` rows make all
columns distinct. At degree two the reciprocal becomes

\[
                  Q_{A',2}(c)
                     =\frac{(c_1+\cdots+c_s)^2-c_1c_2}{4}.
 \tag{10}
\]

It is not real stable. An explicit zero with every variable in the
open upper half plane is

\[
 c_1=\frac{\sqrt{s^2-1}+i}{s},\qquad
 c_2=\frac{-\sqrt{s^2-1}+i}{s},\qquad
 c_j=\frac{i}{s}\ (j\ge3).
 \tag{11}
\]

Here `sum c_j=i` and `c_1c_2=-1`, so (10) vanishes.
Restricting to the first two coefficient variables gives
`(c_1^2+c_1c_2+c_2^2)/4`; its binomial-normalized coefficient list
`(1/4,1/8,1/4)` is not log-concave. Its Hessian is positive definite
on that two-dimensional coordinate subspace, so it also fails the
quadratic Lorentzian signature condition.

This is an actual integer-toric factorial frame, rather than a
general positive-measure or fractional-multiplicity example. It rules
out universal row-merger preservation of reciprocal real stability,
normalized log-concavity, or the Lorentzian condition. It does not
refute regularity preservation or SFC: when `s=2`, the source is a
proved identical pair and the target is `T(c_1Y_1+c_2Y_2)`, which is
also proved. The auxiliary condition fails between two valid endpoints.

## 5. Verification and remaining obligations

`verify_sfc_overlapping_swap_pairs.py` independently constructs actual factorial
coefficients, reconstructs their normalization by rising-product ratios,
checks exact rational Bezout identities for consecutive orders, tests
equal-row factors and linear padding, and verifies (7), (10), and (11).
It retains finite examples satisfying (6) outside `D_0>=2`, known SFC
examples where the criterion fails, and malformed/coefficient-corruption
controls. The finite cases check the algebra; Theorem 1 is the uniform
proof above with the credited interlacing input.

Run from the repository root:

```text
python research/factorial/image-rigidity/verify_sfc_overlapping_swap_pairs.py
```

The unrestricted discrete obligations remain: adjoining one unit row
to an arbitrary nonconstant column, and merging arbitrary integer rows
while keeping distinct columns. Neither follows from this restricted
symmetry class, and neither is claimed. The threshold failure (8)
also shows exactly why the present coefficient method cannot settle
all swap-invariant overlapping pairs.
