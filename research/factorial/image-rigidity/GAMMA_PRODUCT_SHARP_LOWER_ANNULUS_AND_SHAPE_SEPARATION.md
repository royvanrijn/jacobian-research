# A classical Gamma-factor circle bound and sharper product separation

For \(\gamma>0\) and integers \(n,d\geq1\), write
\[
 G_{n,\gamma}(z)=\sum_{k=0}^n\binom nk(\gamma)_kz^{n-k},\qquad
 T_{n,\gamma,d}(t)=\sum_{k=0}^n\binom nk(\gamma)_{dk}t^k.
\]
The circle bound for the individual Gamma factor below is classical. Its
combination with circular composition gives an all-shape product annulus
and sharper sufficient conditions for separation of actual moment orders.
No unrestricted all-order Strong Factorial assertion follows for a fixed
lower exponent.

## 1. The classical factor bound, with an elementary proof

**Lemma GammaFactorInnerCircle.** Every zero of \(G_{n,\gamma}\) is
simple and has modulus at least \(\gamma\). The modulus inequality is
strict if \(n\geq2\). For \(n=1\) the sole zero is \(-\gamma\).

Coefficient comparison gives the exact differential equation
\[
 zG''+(1-\gamma-n-z)G'+nG=0.                         \tag{A1}
\]
For example, at the coefficient of \(z^{n-k}\), the cancellation is
\(k\binom nk(\gamma)_k=(n-k+1)(\gamma+k-1)
\binom n{k-1}(\gamma)_{k-1}\).
The constant coefficient is \((\gamma)_n>0\), so every zero is nonzero.
If a zero \(z_0\ne0\) had multiplicity \(l\geq2\), the term
\(zG''\) in(A1) would have a nonzero coefficient in degree \(l-2\)
of its expansion at \(z_0\), while the other terms have higher orders.
This proves simplicity without invoking analytic ODE uniqueness.

Choose a zero \(z_i\) of minimum modulus, and label the other zeros
\(z_j\). Evaluating(A1) at \(z_i\) gives
\[
 2\sum_{j\ne i}\frac{z_i}{z_i-z_j}=z_i+\gamma+n-1.    \tag{A2}
\]
For \(u=z_j/z_i\), one has \(|u|\geq1\), \(u\ne1\), and
\[
 \Re\frac1{1-u}
  =\frac12+\frac{1-|u|^2}{2|1-u|^2}\leq\frac12.       \tag{A3}
\]
Thus \(\Re z_i\leq-\gamma\), whence \(|z_i|\geq\gamma\).
Minimum modulus proves the lower bound for every zero. Notice that this
argument bounds the real part of the chosen minimum-modulus zero only;
it does not assert that all zeros have real part at most \(-\gamma\).

If \(n\geq2\) and \(|z_i|=\gamma\), then \(z_i=-\gamma\),
and equality holds in each use of(A3). Consequently every zero has
modulus \(\gamma\). Vieta gives \(\sum_jz_j=-n\gamma\).
Since each \(\Re z_j\geq-\gamma\), every zero must be \(-\gamma\),
contrary to simplicity. This proves strictness.

The upper factor bound used below is the retained
**Gamma matrix disk lemma**, `lem:gamma-disk` in
`research/factorial/image-rigidity/GAMMA_PRODUCT_ROOT_MODULUS.tex`:
\[
 |z|\leq\gamma+n-1,                                  \tag{A4}
\]
strict for \(n\geq2\). Its proof uses the exact Gamma tridiagonal
matrix and a positive comparison vector for every Hermitian part;
it holds for every real \(\gamma>0\), without \(\gamma\geq n\).
This is an explicit imported proof input, not a consequence of(A1)--(A3).

**Classical comparison.** In the inspected primary scan of M. G. de Bruin,
E. B. Saff and R. S. Varga, *On the zeros of generalized Bessel
polynomials. I*, Proceedings A84(1), 20 March1981, printed2, definition(2.1)
and assumption(2.3), the convention is
\[
 y_n(v;b)=\sum_{k=0}^n\binom nk(n+b-1)_k(v/2)^k,
 \quad b\in\mathbb R,\quad n+b-1>0.
\]
Taking \(b=\gamma-n+1\) gives
\(G_{n,\gamma}(z)=z^ny_n(2/z;b)\). Printed3, Theorem3.1, gives the
Saff--Varga cardioid
\[
 \{re^{i\theta}:0<r<(1-\cos\theta)/(n+b-1)\}
       \cup\{-2/(n+b-1)\}.
\]
The following paragraph credits Dočev's earlier containing circle
\(|v|\leq2/(n+b-1)\). Under inversion this circle is exactly the
lower bound \(|z|\geq\gamma\). Thus(A1)--(A3) provide a self-contained
proof of a classical bound, not a new priority claim. The stronger
cardioid theorem is context only and is not used as a proof premise.
The [cached scan](sources/varga-generalized-bessel-1981-I-II.pdf)
contains papers I and II together; the locators above refer to paper I.

## 2. Every positive shape: a sharper product annulus

**Theorem GammaProductSharpInnerAnnulus.** For every root \(t\) of
\(T_{n,\gamma,d}\),
\[
 (\gamma)_d\ \leq\ |1/t|\ \leq\ (\gamma+d(n-1))_d.     \tag{A5}
\]
Both inequalities are strict if \(n\geq2\). If \(n=1\), both are
equalities. In particular the lower inequality requires only
\(\gamma>0\), not \(\gamma>d(n-1)\).

Here are the composition conventions and complete transfer argument.
For degree-\(n\) polynomials in the binomial basis put
\[
 p(z)=\sum_{k=0}^n\binom nk a_kz^k,\quad
 q(z)=\sum_{k=0}^n\binom nk b_kz^k,\quad
 p\diamond q=\sum_{k=0}^n\binom nk a_kb_kz^k.
\]
Use the retained self-contained **ClosedDiscGraceComposition** lemma,
Section1 of `GAP_TWO_GRACE_GLOBAL_LOCALIZATION.md` in the same repository
directory: if all roots of \(p\) lie in \(D(c,R)\) and \(q(0)\ne0\),
every root of \(p\diamond q\) is \(-\alpha\beta\), for some
\(\alpha\in D(c,R)\) and some root \(\beta\) of \(q\).
For discs centered at zero this multiplies upper root-modulus bounds.
The minus sign is essential but does not affect these modulus bounds.

For lower bounds use reversal
\(\mathcal R_np(z)=z^np(1/z)\). Direct coefficient comparison gives
\[
 \mathcal R_n(p\diamond q)
    =(\mathcal R_np)\diamond(\mathcal R_nq).           \tag{A6}
\]
When both the leading and constant coefficients are nonzero, reversal
preserves the degree and replaces roots by their reciprocals. Thus if
the roots of the two factors have moduli at least \(L_1,L_2>0\),
apply the same closed-disc lemma to their reversals and reverse back:
the composed polynomial has lower bound \(L_1L_2\). Its leading and
constant coefficients remain nonzero. This avoids any exterior-domain
extension of Grace's theorem.

Set \(\alpha_j=(\gamma+j)/d>0\), \(0\leq j<d\), and let
\(C=G_{n,\alpha_0}\diamond\cdots\diamond G_{n,\alpha_{d-1}}\).
The rising-factorial identity gives, coefficient by coefficient,
\[
 (\gamma)_{dk}=d^{dk}\prod_{j=0}^{d-1}(\alpha_j)_k,
 \quad
 A_{n,\gamma,d}(x):=x^nT_{n,\gamma,d}(1/x)
             =d^{dn}C(x/d^d).                         \tag{A7}
\]
The factor bounds(A3)--(A4), closed-disc composition, and reversal give
\[
 \prod_j\alpha_j\leq |z|\leq\prod_j(\alpha_j+n-1)
          \quad(C(z)=0).
\]
Every factor has finitely many roots. For \(n\geq2\), choose upper
radii strictly smaller than \(\alpha_j+n-1\), and reversed upper
radii strictly smaller than \(1/\alpha_j\), before composing. Both
product bounds are consequently strict. Multiplying moduli by \(d^d\)
in(A7) proves(A5). The degree-one case is \(T=1+(\gamma)_dt\).

**Corollary OrdinaryNormalizedGlobalAnnulus.** For every real \(a>0\),
integers \(m,d\geq1\), and every root \(x\) of the actual ordinary
reversed polynomial
\[
 A_m(x)=\sum_{k=0}^m\binom mk(am+1)_{dk}x^{m-k},
 \qquad
 a^d<|x|/m^d\leq(a+d)^d.                             \tag{A8}
\]
Indeed the lower endpoint in(A5), divided by \(m^d\), is
\(\prod_{j=0}^{d-1}[a+(1+j)/m]>a^d\). The upper endpoint is
\(\prod_{j=0}^{d-1}[a+d+(1-d+j)/m]\leq(a+d)^d\).
These arguments also cover \(m=1\). Ordinary polynomial phases use
integer \(a\); (A8) as a polynomial statement allows real \(a>0\).

## 3. Shape separation and a sharper soft complement

**Corollary ImprovedShapeShiftSeparation.** Let \(d\geq1\),
\(1\leq r<s\) be integers, \(\gamma>0\), and
\(\Delta\geq d(r-1)\). Then \(T_{r,\gamma,d}\) and
\(T_{s,\gamma+\Delta,d}\) are coprime over \(\mathbb C\).

Every reciprocal root of the first has modulus at most
\(U=(\gamma+d(r-1))_d\). Every reciprocal root of the second has
modulus strictly greater than \(L=(\gamma+\Delta)_d\), because
\(s\geq2\). The rising factorial is strictly increasing on positive
arguments, so \(L\geq U\). The intervals are disjoint, including at
equality in the shift hypothesis. No condition on \(r+s\) is required.

**Corollary ImprovedLinearLowerExponentCriterion.** Let \(a\geq0\)
and \(\beta>0\) be real, let \(d\geq1\) be an integer, and put
\(P_j=T_{j,aj+\beta,d}\). Then
\[
 a(s-r)\geq d(r-1)\quad\Longrightarrow\quad
                    \gcd(P_r,P_s)=1.                 \tag{A9}
\]
In particular, for every integer \(N\geq2\),
\[
 a\geq d(N-2)\quad\Longrightarrow\quad
             P_1,\ldots,P_N\text{ pairwise coprime}.  \tag{A10}
\]
Use the preceding corollary with \(\gamma=ar+\beta\) and
\(\Delta=a(s-r)\). For(A10),
\((r-1)/(s-r)\leq N-2\), with the maximum attained at
\((r,s)=(N-1,N)\). Equality at the lower-exponent threshold is allowed.
These are sufficient criteria, not necessary criteria.

For integer \(a\geq0\), the Gamma-\(\beta\) moment functional
\(\mathcal L_\beta(X^k)=(\beta)_k\) gives the exact identity
\[
 \mathcal L_\beta\bigl([X^a(A+BX^d)]^j\bigr)
       =A^j(\beta)_{aj}P_j(B/A).                     \tag{A11}
\]
Thus no nonzero complex phase satisfying(A10) has two zero moments
among orders1 throughN. If one coefficient vanishes, the phase is a
monomial and all its positive moments are nonzero. For rational positive
\(\beta\) and integer \(a,d\), the coprime polynomials have rational
coefficients, so their Bezout identity persists over every
characteristic-zero field. More generally, for fixed real parameters
the complex coprimality identity descends to the real coefficient field
and persists under extensions of that field. No common-root converse
over a specified non-algebraically-closed field is asserted.

For the ordinary critical family \(a=(d-1)^2/4\), the strict complement
of the adjacent sufficient condition is
\[
 a<d(N-2)\quad\Longrightarrow\quad d<4N-6,
 \qquad d/(N-1)<4.                                   \tag{A12}
\]
Indeed \(d^2-(4N-6)d+1<0\), so \(d<4N-6\). This is a necessary
condition for a possible failure at terminal orderN, not an existence
claim or an all-order cusp theorem. The sharper criterion changes
none of the currently retained all-order or activation questions.

## Exact reproduction and limits

Run the standard-library checker from this package directory:

    python3 verify_gamma_product_sharp_lower_annulus.py \
      --repo-root /path/to/jacobian-research --output review.json

The mandatory ledger pins this proof, the checker, exact fixture input,
the inspected1981 primary scan, and both precise retained proof imports.
The checker verifies full formal-parameter ODE identities, the cleared
complex inequality, exact binomial composition/reversal and dilation,
actual moment normalization, rational coprimality controls at the new
threshold including equality, and guards against several false
extensions. These finite controls support the written universal proof;
they do not replace it. No assertion of simplicity for the product
polynomial, an optimal separation threshold, a fixed-a all-order result,
external review, or formal verification is made.
