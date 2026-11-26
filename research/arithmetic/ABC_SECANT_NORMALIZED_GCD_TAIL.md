# A normalized-scale tail for common secant gcds

For a near-full balanced abc family with small complete
multiplicity, the gcd of two coordinate differences is
typically smaller than any prescribed positive power of
\(n=X^L\). The exceptional fraction is power-small at each
fixed budget. Its exponent can vanish as the budget tends
to zero.

The proof uses shallow residue occupancy to remove the
overlap of each secant gcd with the union of the previous
ones. The remaining lcm divisibility forces a polynomial
through the genuine complete values.
[Shallow residue-list counts (SMI1–SMI2)](ABC_SHALLOW_MOMENT_INTERPOLATION.md#a-residue-list-bound-from-interpolated-moments)
· [Random interpolation (JSI1)](ABC_SHALLOW_JOINT_INTERPOLATION.md#interpolation-against-proper-congruence-hits)

## A radical-scale tail at each fixed budget

Put \(d=1/3\), and use the absolute constants
\(\sigma_*,L_*>0\) from the complete-letter curve bound
(ACG), with \(0<\sigma_*<d\).
For each fixed \(0<L<L_*\), it applies to any finite
degree cap fixed before \(X\to\infty\), with a count
exponent \(d-\sigma_*\) in base \(n=X^L\) up to the
inner allowance \(X^\varepsilon\). It is uniform in
the polynomial at that fixed degree cap. The degree can
depend on \(L\).

Take a family \(S\) of \(N\) distinct positive primitive
triples \((A,B,C)\), \(A+B=C\le X\), with
\[
 \operatorname{rad}(ABC)<C^L,\qquad
 \operatorname{rad}(T)\le n^{d+\kappa}\quad(T=A,B,C),
 \qquad N\ge n^{d-\kappa}.
\]
Assume each complete value on each role occurs in at
most \(n^u\) rows. Identify a row \(x\) with its pair
\((A,B)\). For distinct rows set
\[
           g(x,y)=\gcd(|A(x)-A(y)|,\ |B(x)-B(y)|).
\]
This gcd is positive, with the convention
\(\gcd(0,a)=a\) for \(a>0\).
The same \(g\) is obtained from either other pair of
complete roles.

**NGT1.** Fix \(0<\zeta\le d\), and suppose
\[
   u\le\zeta/32,\qquad
   0<\kappa\le\min(\zeta/128,\sigma_*/8).
\]
Set
\[
 r=\left\lceil\frac{128}{\zeta L}\right\rceil+4,\qquad
 m=\left\lfloor\frac14\binom{r+2}{2}\right\rfloor,\qquad
 \nu=\min(\zeta/32,\sigma_*/2).
\]
At sufficiently large heights at each fixed \(L\),
\[
  \#\{(x,y)\in S^2:x\ne y,\ g(x,y)\ge n^\zeta\}
                    \ll_{L,\zeta} N^2 n^{-\nu/m}.
 \tag{NGT1}
\]
Any original eligibility restriction may be retained in
choosing \(S\). The exponent is positive at the fixed
budget and is of order \(L^2\) as \(L\to0\) at fixed
\(\zeta\). A threshold \(\zeta>d\) can use the result for \(d\).
[Absolute curve saving at fixed-budget degrees (ACG)](ABC_WIDE_COFACTOR_CONGRUENCES.md#an-absolute-curve-saving-with-arbitrary-fixed-budget-degrees)

### Occupancy at the scale of one overlap

A class for the pair \((A,B)\) modulo an integer
\(D\ge n^{\zeta/2}\) contains \(O_{L,\zeta}(n^{d-\zeta/16})\)
rows of \(S\) if it contains a primitive representative
\(y\) of polynomial height. In the use below \(D\le X\).
To see this, split \(D\) into its full prime powers
according to a choice of coordinate of \(y\) that is
a unit at the corresponding prime. Such a choice exists
at every prime by primitivity. One of the resulting
two factors \(D'\) is at least \(n^{\zeta/4}\).
The \(A\)- or \(B\)-coordinate assigned to that factor,
denoted \(W\), is in a specified unit class modulo \(D'\).

For a uniform one-letter bound, split off the part of
\(W\) on valuation orders above an integer cutoff
chosen at the fixed budget so that its radical is at
most \(n^\kappa\). Enumerating this factor costs
\(n^\kappa X^{o_X(1)}\), by radical enumeration and
uniform supported-part counting. It is a unit and can
enter the known coefficient. The remaining exact-order
root boxes have total normalized cap
\(\rho\le d+\kappa+o_X(1)\), with height-subpower many
dyadic choices. If \(\rho\ge d/2\), the half-depth
part of (SMI2), and monotonicity of its gain in the
modulus depth, give gain at least \(\zeta/4\) at \(D'\).
The list count with one class then has exponent at
most \(\rho-\zeta/8\). If \(\rho<d/2\), the raw
dictionary bound is already smaller than
\(n^{d-\zeta/8}\). Restore the factor record and the
complete degree \(n^u\). The resulting exponent
\(d+2\kappa+u-\zeta/8\) is below \(d-\zeta/16\)
with a fixed margin for the counting allowance.
This proves the occupancy estimate uniformly in
the original class and the modulus.
[Monotonicity and half-depth saturation (SMI1–SMI2)](ABC_SHALLOW_MOMENT_INTERPOLATION.md#a-residue-list-bound-from-interpolated-moments)
· [Uniform supported-factor counting](ABC_EFFECTIVE_RADICAL_BUDGET.md#counting-by-a-sublogarithmic-prime-partition-sum)

### Ordered anchors and fresh numerical depth

Let \(\rho_0 N^2\) be the number of ordered pairs
in (NGT1). Jensen's inequality gives at least
\(\rho_0^m N^{m+1}\) ordered stars with centre \(x\)
and \(m\) neighbors \(y_i\), allowing repeated anchors.
The number with repetitions is \(O_m(N^m)\).
Fix any distinct ordered anchors and put
\[
 g_i=g(x,y_i),\qquad
 P_i=\operatorname{lcm}_{j<i} g(y_i,y_j)\quad(i\ge2).
\]
Then
\[
    \gcd\!\left(g_i,\operatorname{lcm}_{j<i}g_j\right)
        =\operatorname{lcm}_{j<i}\gcd(g_i,g_j)
        \ \mid\ P_i .
\]
The equality follows prime by prime by distribution of
minimum over maximum. For divisibility, simultaneous
congruence of \(x\) to \(y_i,y_j\) at a prime-power depth
implies congruence of those two anchors at that depth.

Call a star exceptional for overlap when some
\(\gcd(g_i,\operatorname{lcm}_{j<i}g_j)\ge n^{\zeta/2}\).
The factor causing this exception is a divisor of the
fixed integer \(P_i\le X^{i-1}\), hence has
\(X^{o_X(1)}\) possible values at fixed \(m\). For every
such candidate the centre \(x\) lies in the pair
residue class of \(y_i\). The preceding occupancy
bound applies. Summation over \(i\) and the anchors
gives
\[
    O_{L,\zeta}\bigl(N^{m+1}n^{-\zeta/32}\bigr)
\]
overlap-exceptional stars, using \(N\ge n^{d-\kappa}\)
and taking the divisor allowance smaller than the
available gap.

On any remaining star the first \(g_i\) is at least
\(n^\zeta\), and each subsequent \(g_i\) increases
the lcm by at least \(n^{\zeta/2}\). Thus
\[
                 \operatorname{lcm}(g_1,\ldots,g_m)
                        \ge n^{m\zeta/2}\ge X^{4r}.
\]
Here \(m\ge r^2/16\) and \(r\ge128/(\zeta L)\).
The ordered overlap test loses at most \(\zeta/2\)
of normalized depth per new anchor, even if many
earlier anchors share different parts of its gcd.

For the fixed anchors, form a nonzero integer
polynomial \(F(A,B)\) of degree at most \(r\)
vanishing on all \(y_i\) and having
\(|F(A,B)|\le X^{3r}\) on the height box. To construct
it, map coefficient vectors in
\(\{0,\ldots,\lfloor X^r\rfloor\}^{V_r}\), where
\(V_r=\binom{r+2}{2}\), to their \(m\) anchor
evaluations. There are at most
\((1+V_r X^{2r})^m<X^{rV_r}\) outputs for large
\(X\), and more inputs. A difference of distinct
equal-evaluation vectors is the required polynomial.
For a centre \(x\) on the star, each \(g_i\mid F(x)\).
The lcm exceeds the bound on \(|F(x)|\), so \(F(x)=0\).
The curve bound with inner allowance
\(\varepsilon=L\sigma_*/4\) leaves at most
\(O_{L,r}(n^{d-3\sigma_*/4})\) centres for each
tuple of anchors. Summing again costs at most
\(O_{L,\zeta}(N^{m+1}n^{-\sigma_*/2})\).
Combine these bounds with the repeated-anchor term
and take \(m\)-th roots to prove (NGT1).

## One-coordinate differences with bounded outside radical

For a row pair with distinct values in a role \(T\), define
\[
      v_T=\frac{T(y)-T(x)}{g(x,y)},\qquad
      a_T=\gcd(T(x),T(y)),\qquad
      Z_T=\max(T(x),T(y))/a_T .
\]
Use the hypotheses of (NGT1), and fix \(B_0,K,h_0>0\).
For each first row fix a pool integer \(Q_x\le X^K\).
Only count pairs with \(Z_T\ge X^{h_0}\) and with
\(\operatorname{rad}_{Q_x^c}(|v_T|)\le n^{B_0}\), where
the outside radical omits primes dividing \(Q_x\).

**NGT2.** For any prescribed fixed allowance \(\eta>0\)
and sufficiently small fixed budgets, this count is at most
\[
 O_{B_0,K,h_0,L,\zeta,\eta,\varepsilon}\!\left(
   N^2 n^{-\nu/m}+
   N n^{\,u+\min(d+\kappa,B_0+\zeta)/2+\eta}
       X^\varepsilon\right).
 \tag{NGT2}
\]
The raw radical of the pool is unrestricted. The output
height condition can be supplied on the usual retained
near-balanced secant sector by (UC3).

Indeed discard the exceptional pairs in (NGT1).
Fix \(x\) and \(a_T\mid T(x)\) at divisor cost.
The three positive terms
\[
     T(x)/a_T,\qquad T(y)/a_T,\qquad
     |T(y)-T(x)|/a_T
\]
give a primitive additive equation after reordering.
The first is fixed and the second has radical at most
\(n^{d+\kappa}\). The outside radical of the third
is at most \(n^{B_0+\zeta}\), because that term
divides \(g(x,y)|v_T|\). Apply (FL1) with this third
pool on each shell of the derived height \(Y\ge X^{h_0}\).
The normalized caps there are multiplied by
\(\log X/\log Y\), a bounded factor, and
\(Q_x\le Y^{K/h_0}\). Converting the resulting
exponent to base \(n\) gives the displayed minimum.
The recovered complete value \(T(y)\) specifies
at most \(n^u\) rows. Sum over the shells, the
divisor choices and the first rows, with the
inner allowances chosen inside \(\eta,\varepsilon\).
[Fixed-letter counts with prescribed pools (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)
· [Derived scalar height in secants (UC3)](ABC_UNRECORDED_COFACTOR_SECANT_TRANSFER.md#primitive-secant-coordinates-with-an-arbitrary-scaling-gcd)

## Counting complete primitive directions below a fixed radical cap

Assume also the original height shell \(X/2<C\le X\),
and impose the small cap tolerance and budget threshold of
the rank-one feature bound (CF1).
For distinct rows write
\[
         \boldsymbol v=\frac{(A(y)-A(x),B(y)-B(x))}{g(x,y)},
         \qquad v_C=v_A+v_B .
\]
For nonzero \(v_Av_Bv_C\), the three absolute values,
ordered with the largest last, form a primitive
additive triple. Let \(H_{\boldsymbol v}\) be their maximum.

**NGT3.** Fix \(B_1,\eta>0\). The number of pairs with
nonzero \(v_Av_Bv_C\) and
\(\operatorname{rad}(|v_Av_Bv_C|)\le n^{B_1}\) is at most
\[
 O_{B_1,L,\zeta,\eta,\varepsilon}\left(
       N^2 n^{-\delta_L}
       + n^{\,B_1/3+\zeta+u+(d+\kappa)/2+\eta}
                          X^\varepsilon\right)
 \tag{NGT3}
\]
for a \(\delta_L>0\) at each sufficiently small fixed
budget. In particular, for every fixed \(B_1<3/2\),
the bound is power-smaller than the original pair count
at sufficiently small \(\zeta,u,\kappa,\eta\),
subject to \(u\le\zeta/32\).

Discard the exception from (NGT1). Use rank one in
(CF1) with the complete target \(C(y)\) as the affine
observable, fixing the first row \(x\).
If \(H_{\boldsymbol v}\le X^\theta\), then
\(|C(y)-C(x)|\le n^\zeta X^\theta\).
Use the feature \(C(y)-C(x)\), with \(\theta\) and
then \(L\) small enough for the feature differences
to fit the (CF1) threshold. Likewise, if
\(\max(C(x),C(y))/\gcd(C(x),C(y))\le X^\theta\),
fix the gcd \(a\mid C(x)\) at divisor cost and use
the feature \(C(y)/a\le X^\theta\).
Each step removes only a power-small fraction of
the original pairs, taking \(\kappa\) and the inner
allowances smaller than the (CF1) gap.
The bound concerns the original primitive triples
in both applications of (CF1).
[Small affine features on one complete role (CF1)](ABC_COMPLETE_LETTER_SMALL_FEATURES.md#a-strict-saving-for-small-affine-minors)

On the rest, \(H_{\boldsymbol v}\ge X^\theta\).
The limiting \(1/3\) theorem bounds the number of
oriented primitive directions of radical at most \(n^{B_1}\)
by \(n^{B_1/3+\eta_1}X^\varepsilon\), for arbitrarily small
fixed \(\eta_1>0\) and sufficiently small fixed budgets.
For the height conversion, cover \([\theta,1]\) by a
finite fine mesh in \(\log_X H_{\boldsymbol v}\).
For a bin with endpoints \(h_-,h_+\), apply the theorem
at height \(X^{h_+}\) with the fixed budget
\((B_1+\eta_0)L/h_-\), where \(\eta_0>0\) is small.
The mesh and the small-budget counting allowance
give the claimed exponent after choosing them in terms
of \(B_1,\theta,\eta_1\). There are boundedly many
signs and orientations.
[Limiting small-budget count](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)

Now fix a direction and a positive scale \(g<n^\zeta\).
The complete target difference \(g v_C\) is fixed
and nonzero. Enumerate \(a=\gcd(C(x),C(y))\) among
divisors of its absolute value. The primitive shifted
target quotients have individual radical caps at most
\(n^{d+\kappa}\), their shift is fixed, and their
maximum is at least \(X^\theta\) on the retained
class. Apply (MCB2) with the shift in the fixed
coefficient and use exact-order singleton layers
on the two moving terms. This counts at most
\(n^{(d+\kappa)/2+\eta_1}X^\varepsilon\)
ordered complete target pairs. For each such target
pair there are at most \(n^u\) choices for the first
row. The fixed displacement \(g\boldsymbol v\)
then gives the second row. Sum over \(g\), the directions
and the divisor choices, with the internal allowances
smaller than \(\eta\). This proves (NGT3).
[Fixed-coefficient half bound (MCB2)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#bounded-group-dictionaries-with-arbitrary-fixed-coefficients)

This is a direct count of all three primitive direction
coordinates with an explicit total radical threshold. The
arbitrary fixed single-coordinate thresholds available from
(UC2) are stronger qualitative exclusions on these
small-multiplicity charts. The explicit direction count
does not use that unrecorded-cofactor estimate.

## Relation to numerical peeling and the endpoint

The fixed-numerical-threshold consequence of (NGT1) also
follows from the earlier peeling pools and unit-class
occupancy on intermediate multiplicity charts. The
normalized result charges the full gcd by the arbitrarily
small test \(n^\zeta\) for correspondingly small complete
degree. The previous normalized peeling bound (BC5)
instead gives \(g/\gcd(g,P_x)\le n^{C_*}\) for an
absolute \(C_*\), allowing an anchor-dependent pool.
One could also use the absolute fixed-budget degree
bound (ACG) in that greedy peeling argument with
degree depending on \(\zeta,L\), and then use the
shallow occupancy estimate to control the pooled part.
The ordered-anchor proof avoids an assigned pool for each
row and makes the dependence of the exceptional exponent
explicit.
[Numerical and normalized peeling (BC3–BC5)](ABC_BOUNDED_DEGREE_CURVES_AND_SECANT_GCD.md#peeling-down-to-the-normalized-radical-scale)
· [Numerically large gcds and minimum dictionaries](ABC_UNRECORDED_COFACTOR_SECANT_TRANSFER.md#numerically-large-scaling-gcds-without-a-direction-cap)

When \(B_0+\zeta<d+\kappa\), (NGT2) prices the
short-radical coordinate rather than the other complete
value. The earlier (UC2) already excludes arbitrary
fixed directional caps on full pair-scale intermediate
multiplicity families, with an absolute small-budget
gap for a fixed sector. The bound here is useful when
one needs the explicit smaller conditional price or a
full gcd bounded at the normalized scale.
The factor \(1/m\) in (NGT1) records the loss from
building a star out of a sparse pair population.
Since \(m\) grows as \(L^{-2}\) at fixed \(\zeta\),
this estimate still lacks an absolute exponent gap
as the budget tends to zero. It also needs a positive
small-multiplicity margin measured against \(\zeta\).
