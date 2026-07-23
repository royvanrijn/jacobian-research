# Higher-dimensional GMC on restricted supports

## 1. Scope and conclusion

The global Gaussian Moments Conjecture is false in every real dimension at
least three.  Nevertheless, the lower-face/prime-isolation proof in two real
variables has a genuine multiradial extension on a large, explicitly
checkable class of supports.

The correct condition is slightly stronger and cleaner than uniqueness on
one exposed face:

> the balanced radial polytope must have a **least vector** in the
> componentwise order.

Equivalently, its coordinatewise infimum must be attained.  Under this
condition, pure Gaussian moment-nullness forces instability for the angular
torus, and torus instability gives the eventual mixed-moment conclusion.

For a fixed polynomial, only its full balanced radial polytope is needed.
Requiring the same property on every balanced sub-support is a convenient
uniform hypothesis for a class of allowed supports.  Requiring it on every
balanced face is sufficient but stronger than the pointwise theorem.

The proof below is for \(2r\) real Gaussian variables, arranged as \(r\)
circular complex pairs.  A real spectator can also be handled after adding
parity data and replacing one factorial by a double factorial, but that
variant is not needed for the main support theorem and is recorded only as a
remark.

## 2. Circular support data

Let

\[
 Z_i=\frac{X_i+iY_i}{\sqrt2},\qquad
 W_i=\frac{X_i-iY_i}{\sqrt2}
 \quad(1\leq i\leq r).
\]

Then

\[
 \mathbb E(Z^\alpha W^\beta)
 =\mathbf 1_{\alpha=\beta}\,\alpha!,
 \qquad
 \alpha!=\prod_{i=1}^r\alpha_i!.
\]

Make the injective substitution

\[
 Z_i\longmapsto T_i,\qquad
 W_i\longmapsto U_iT_i^{-1}.
\]

After combining equal monomials, write

\[
 P=\sum_{a\in A}c_aT^{w_a}U^{n_a},
 \qquad
 c_a\in\mathbb C^\times,\quad
 w_a\in\mathbb Z^r,\quad n_a\in\mathbb N^r.
 \tag{2.1}
\]

The angular weight is \(w_a\), and the radial degree is \(n_a\).  If

\[
 F_m(U)=\operatorname{CT}_{T}(P^m),
 \]

then

\[
 \mathbb E(P^m)=\mathcal L(F_m),\qquad
 \mathcal L(U^n)=n!:=\prod_i n_i!.
 \tag{2.2}
\]

Assume \(0\in\operatorname{conv}\{w_a:a\in A\}\).  Define the balanced
radial polytope

\[
 \mathcal B(A)=
 \left\{
 \sum_{a\in A}\lambda_an_a:
 \lambda_a\geq0,\quad
 \sum_a\lambda_a=1,\quad
 \sum_a\lambda_aw_a=0
 \right\}.
 \tag{2.3}
\]

It is the set of average radial degrees of fractional weight-zero products.

We say that \(A\) has the **least-balanced-degree property** if there is
\(\rho\in\mathcal B(A)\) such that

\[
 \rho\leq x\quad\text{componentwise for every }x\in\mathcal B(A).
 \tag{2.4}
\]

Thus \(\rho\) is a least element, not merely a minimal element selected by
one scalar objective.  A practical equivalent test is

\[
 \rho_i=\min_{x\in\mathcal B(A)}x_i\quad(1\leq i\leq r),
 \qquad
 (\rho_1,\ldots,\rho_r)\in\mathcal B(A).
 \tag{2.5}
\]

All these tests are rational linear programs.

## 3. Least-balanced-degree theorem

> **Theorem 3.1 (multiradial prime isolation).**  
> Let \(P\) be as in (2.1), and suppose that its support \(A\) has the
> least-balanced-degree property whenever \(0\) lies in its angular weight
> convex hull.  If
> \[
> \mathbb E(P^m)=0\qquad(m\geq1),
> \tag{3.1}
> \]
> then
> \[
> 0\notin\operatorname{conv}\{w_a:a\in A\}.
> \tag{3.2}
> \]
> Consequently \(P\) is unstable for the angular torus
> \((\mathbb C^\times)^r\), and for every polynomial \(Q\),
> \[
> \mathbb E(QP^m)=0
> \qquad\text{for all sufficiently large }m.
> \tag{3.3}
> \]

### Proof

Suppose instead that \(0\) lies in the angular weight convex hull, and let
\(\rho\) be the least point of \(\mathcal B(A)\).  Minimize the strictly
positive scalar functional

\[
 |x|_1=x_1+\cdots+x_r
 \]

on \(\mathcal B(A)\).  The least-point condition makes \(\rho\) its unique
minimizer.  Linear-programming duality supplies
\(\theta\in\mathbb Q^r\) and \(\mu=|\rho|_1\) such that

\[
 |n_a|_1\geq \mu+\theta\mathbin{\cdot}w_a
 \qquad(a\in A).
 \tag{3.4}
\]

Let \(A_0\) be the equality set and

\[
 P_0=\sum_{a\in A_0}c_aT^{w_a}U^{n_a}.
 \tag{3.5}
\]

An optimal balanced distribution is supported on \(A_0\), so zero lies in
the convex hull of the angular weights occurring in \(P_0\).  Apply the
Duistermaat--van der Kallen constant-term theorem to \(P_0\), viewed as a
Laurent polynomial in \(T\) with coefficients in \(\mathbb C(U)\).  One
may equivalently specialize \(U\) generically and apply the theorem over
\(\mathbb C\).  For some \(s\geq1\),

\[
 \operatorname{CT}_{T}(P_0^s)\neq0.
 \tag{3.6}
\]

Every balanced product from \(A_0\) has average radial degree \(x\) in
\(\mathcal B(A)\) with \(|x|_1=|\rho|_1\).  Since \(x\geq\rho\)
componentwise, necessarily \(x=\rho\).  Hence (3.6) is one monomial:

\[
 \operatorname{CT}_{T}(P_0^s)=cU^n,
 \qquad n=s\rho\in\mathbb N^r,\quad c\neq0.
 \tag{3.7}
\]

It is the coefficient at the unique lowest balanced radial vector in
\(F_s=\operatorname{CT}_{T}(P^s)\).

More importantly, every monomial \(U^d\) in \(F_{sp}\) comes from a
balanced product of length \(sp\).  Its average degree belongs to
\(\mathcal B(A)\), so (2.4) gives

\[
 d\geq sp\rho=pn
 \quad\text{componentwise}.
 \tag{3.8}
\]

Adjoin \(c^{-1}\) to the finitely generated coefficient domain and reduce
at a sufficiently good rational prime \(p\).  The assumed moment identity
at exponent \(sp\) can be divided integrally by

\[
 (pn)!:=\prod_i(pn_i)!.
 \tag{3.9}
\]

In characteristic \(p\), Frobenius and constant-term extraction give

\[
 F_{sp}
 =\operatorname{CT}_{T}\bigl((P^s)^p\bigr)
 =F_s^p.
 \tag{3.10}
\]

Thus every surviving exponent is \(pj\), and the coefficient at \(pn\) is
\(c^p\).  If \(j\neq n\), (3.8) gives \(j\geq n\), so some coordinate
satisfies \(j_i>n_i\).  The normalized factorial multiplier contains

\[
 \frac{(pj_i)!}{(pn_i)!},
 \]

which is divisible by \(p\).  All non-\(p\)-dilated coefficients vanish by
Frobenius.  Reducing the normalized moment identity modulo \(p\) therefore
leaves \(c^p=0\), a contradiction.

This proves (3.2).  By strict separation, some integral one-parameter
subgroup \(\eta\in\mathbb Z^r\) has
\(\eta\cdot w_a>0\) for every \(a\).  The \(\eta\)-weights of \(P^m\)
then escape linearly to \(+\infty\), while those of a fixed \(Q\) remain
bounded.  Hence \(QP^m\) has no angular weight zero for large \(m\), which
proves (3.3). \(\square\)

## 4. What “unique componentwise-minimal” must mean

For a finite partially ordered set, a unique minimal element is
automatically least.  The balanced object (2.3), however, is a convex
polytope with infinitely many points.  The safest formulation is the
attained-coordinatewise-infimum condition (2.5).

A unique minimizer of one strictly positive functional is not enough.
For example,

\[
 \mathcal B=\operatorname{conv}\{(0,2),(1,0)\}
\]

has a unique minimizer for many positive linear functionals, but neither
endpoint is below the other.  An off-face term may decrease one coordinate,
so its factorial is not divisible by the factorial attached to the exposed
term.  The normalized reduction modulo \(p\) then has no reason to be
integral.

This is the exact point at which scalar lower-face exposure in two real
variables differs from multiradial exposure.

## 5. Support classes

### 5.1 Affine radial degree as a function of angular weight

Suppose there are \(b\in\mathbb Q^r\) and a rational linear map \(L\) such
that

\[
 n_a=b+Lw_a\qquad(a\in A).
 \tag{5.1}
\]

Then every balanced distribution has average radial degree \(b\), so

\[
 \mathcal B(A)=\{b\}.
\]

The theorem applies.  This includes polynomials that are separately
homogeneous in each real Gaussian pair: if the degree in
\((X_i,Y_i)\) is \(d_i\), then in circular coordinates

\[
 n_i=\frac{d_i-w_i}{2}.
\]

This recovers, and conceptually extends, the multihomogeneous special case
proved by Derksen--van den Essen--Zhao.

### 5.2 Monotone affine-line and rank-one radial supports

Suppose that on the balanced fiber

\[
 n=b+qv
\]

with \(v\in\mathbb R_{\geq0}^r\) (or
\(v\in\mathbb R_{\leq0}^r\)).  The endpoint with smallest \(q\) (or largest
\(q\)) is a least vector.  Thus affine-line and rank-one radial supports
work when their radial direction is contained in one closed orthant.

The sign condition is essential.  Long's four-real-variable counterexample

\[
 P_4=(1+Z_2)\bigl(W_1(1-Z_1)+W_2\bigr)
 \tag{5.2}
\]

has two circular radial vectors

\[
 (1,0)\quad\text{and}\quad(0,1).
\]

Both occur on angular weight-zero monomials, so its balanced radial
polytope contains the two incomparable endpoints and lies on the affine
line \(x+y=1\).  Thus “radial degrees lie on one affine line” and
“rank-one radial support” are both false without monotonicity; (5.2) is
already a counterexample in the literal multifactorial setting.

### 5.3 A distinguished minimal balanced circuit

Let

\[
 \Delta_A=\left\{
 \lambda\geq0:\sum_a\lambda_a=1,\ 
 \sum_a\lambda_aw_a=0
 \right\}.
\]

Its vertices are the normalized positive circuits, including zero-weight
atoms.  The polytope \(\mathcal B(A)\) is their image under
\(\lambda\mapsto\sum_a\lambda_an_a\).

Consequently it is enough to find one circuit vector \(\rho\) such that

\[
 \rho\leq
 \frac{\sum_a\gamma_an_a}{\sum_a\gamma_a}
 \quad\text{componentwise}
 \tag{5.3}
\]

for every positive circuit \(\gamma\).  Then every point of
\(\mathcal B(A)\), being a convex combination of circuit images, lies above
\(\rho\).  This validates the “one distinguished minimal balanced circuit”
idea with two qualifications:

1. circuit radial degrees must be divided by circuit length;
2. domination must be checked against every vertex of the balanced
   relation polytope, not only against an arbitrarily chosen Hilbert-basis
   list.

Several circuits may have the same least radial vector.  Uniqueness of the
circuit is unnecessary; the Duistermaat--van der Kallen step prevents
coefficient cancellation from persisting through every power.

### 5.4 Strictly positive separating functionals

A positive functional is useful for exposing the initial form, but one
functional alone does not control coordinate tradeoffs.  The sufficient
condition is the cone containment

\[
 \mathcal B(A)\subseteq \rho+\mathbb R_{\geq0}^r,
 \tag{5.4}
\]

with \(\rho\in\mathcal B(A)\).  Equivalently, all coordinate minima are
simultaneously attained.  Thus a “strictly positive separating functional”
subclass is valid only when separation is strengthened to (5.4), or when
additional support geometry makes scalar separation imply (5.4).

### 5.5 One complex Gaussian linear form and spectators

If \(P=F(L)\) depends on only one complex Gaussian linear form \(L\), while
the spectator variables occur only in \(Q\), GMC reduces to the one-form
case and holds.  If \(\mathbb E(L^2)\neq0\), a complex orthogonal change of
Gaussian coordinates makes \(L\) a nonzero multiple of one real Gaussian
coordinate, and the univariate theorem applies.  If
\(\mathbb E(L^2)=0\), then \(L\) can be completed to a circular pair
\((Z,W)\); \(F(Z)\) has only nonnegative angular weights, and pure
moment-nullness forces its weight-zero coefficient to vanish.  The
remaining support is strictly positive.  Spectators in \(Q\) do not alter
the escape of angular weight.

More generally, spectator dependence in \(P\) is safe when it preserves
the least-balanced-degree property, for example when every term has the
same spectator monomial factor or when the spectator radial degrees form a
monotone chain.  Arbitrary spectator dependence is not safe: Long's
three-real-variable polynomial

\[
 W+WZ-T^2-\frac32ZT^2-\frac12Z^2T^2
\]

uses one circular pair and one real spectator precisely to create the
incomparable tradeoff \((1,0)\leftrightarrow(0,1)\).

## 6. A finite support-class criterion

Let \(S\) be an allowed monomial support, and permit arbitrary polynomials
whose actual nonzero support is a subset \(A\subseteq S\).  Then GMC holds
for every such polynomial if:

\[
 \text{for every }A\subseteq S\text{ with }
 0\in\operatorname{conv}\{w_a:a\in A\},
 \quad\mathcal B(A)\text{ has a least vector}.
 \tag{6.1}
\]

This is a finite combinatorial criterion.  For any fixed \(A\), it can be
checked by:

1. solving the \(r\) coordinate linear programs in (2.5);
2. testing whether the vector of coordinate minima belongs to
   \(\mathcal B(A)\).

For structured families, it is usually better to prove (5.1), monotone
rank one, or circuit domination than to enumerate all sub-supports.

## 7. Relation to the proposed balanced-face theorem

The proposed statement

> if every balanced face has a unique componentwise-minimal radial degree,
> moment-nullness implies torus instability

is correct after making “unique componentwise-minimal” mean “the
coordinatewise infimum is attained,” and it is stronger than necessary.
The proof only needs that property for the full balanced radial polytope of
the actual support.  A hereditary theorem for all polynomials supported in
a fixed allowed set should quantify over all balanced sub-supports as in
(6.1), rather than only over geometrically exposed faces.

The surviving piece of the two-dimensional proof is therefore exact:

\[
\boxed{
\begin{array}{c}
\text{balanced radial fiber has a least vector}\\
\Downarrow\\
\text{one exposed initial form contributes }cU^n\\
\Downarrow\ \text{prime dilation}\\
\text{every other surviving vector raises a coordinate}\\
\Downarrow\\
\text{a factorial ratio is divisible by }p.
\end{array}}
\]

The higher-dimensional obstruction is equally exact: incomparable radial
vectors permit factorial tradeoffs, and Long's counterexample realizes the
smallest such tradeoff.

## 8. Odd real spectators

For an unpaired real Gaussian \(T\),

\[
 \mathbb E(T^{2k})=(2k-1)!!,\qquad
 \mathbb E(T^{2k+1})=0.
\]

After parity is incorporated into the balancing data, the arithmetic
isolation still works at odd primes: if \(j>n\), then

\[
 \frac{(2pj-1)!!}{(2pn-1)!!}
\]

contains the odd multiple \(p(2n+1)\).  What must be added is a
\(\mathbb Z/2\)-character bookkeeping step ensuring that the exposed
balanced term has even spectator degree.  The three-variable counterexample
shows that this bookkeeping does not restore the theorem when the
circular and spectator radial degrees are incomparable.

## 9. Literature position

The ingredients are classical or already present in GMC:

- Duistermaat and van der Kallen prove the Laurent-polynomial
  constant-term theorem used to produce a nonzero power of the exposed
  initial form.
- Derksen, van den Essen, and Zhao formulate GMC, prove the univariate and
  multihomogeneous special cases, and use essentially the same prime
  divisibility of Gaussian factorials in their Propositions 4.6--4.7.
- Long gives explicit failures in every dimension at least three; his
  three-variable witness is also the sharp warning against nonmonotone
  rank-one radial support.

The least-balanced-degree formulation packages these ingredients into a
support criterion.  No claim of priority is made here; a literature search
did not locate this exact polyhedral statement as a named theorem.

Primary references:

- H. Derksen, A. van den Essen, and W. Zhao,
  [*The Gaussian Moments Conjecture and the Jacobian
  Conjecture*](https://arxiv.org/abs/1506.05192), Israel J. Math. 219
  (2017), 917--928.
- J. J. Duistermaat and W. van der Kallen,
  [*Constant terms in powers of a Laurent
  polynomial*](https://doi.org/10.1016/S0019-3577(98)80020-7),
  Indag. Math. 9 (1998), 221--231.
- C. D. Long,
  [*Small Counterexamples to the Gaussian Moments
  Conjecture*](https://arxiv.org/abs/2607.18186), 2026.
