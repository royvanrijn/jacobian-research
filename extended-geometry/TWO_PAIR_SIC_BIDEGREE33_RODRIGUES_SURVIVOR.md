# The bidegree-\((3,3)\) Rodrigues survivor

## 1. Statement and status

Use the contraction pairs \((W,Z),(V,Y)\), and put

\[
\begin{aligned}
F={}&W^3Z^2Y+W^3Y^3-W^2VZ^3-2WV^2Z^2Y\\
   &\qquad-WV^2Y^3-V^3ZY^2.                         \tag{1.1}
\end{aligned}
\]

Its coefficient matrix in
\[
M_{ij}=W^{3-i}V^iZ^{3-j}Y^j,\qquad 0\leq i,j\leq3,
\]
is
\[
C=
\begin{pmatrix}
0&1&0&1\\
-1&0&0&0\\
0&-2&0&-1\\
0&0&-1&0
\end{pmatrix},
\qquad \det C=1.                                      \tag{1.2}
\]

> **Theorem 1.1 (full-rank moment--nullcone survivor).** Over every
> characteristic-zero field,
> \[
> \mathcal E_2(F^m)=0\qquad(m\geq1).                  \tag{1.3}
> \]
> Nevertheless \(F\) is not in the one-sided nullcone, since its
> coefficient matrix is invertible.  Consequently the
> moment--nullcone assertion \(\mathrm{MN}_3\) is false.

This is not an SIC counterexample.

> **Theorem 1.2 (explicit SIC safety).** Let \(Q\) be any polynomial and
> let \(e\) be the largest coordinate degree in \(Z,Y\) of one of its
> bihomogeneous components.  Then
> \[
> \boxed{\mathcal E_2(QF^m)=0\quad\text{for every }m>3e.} \tag{1.4}
> \]

Thus the first semistable all-order pure-moment survivor in balanced
degree three is a genuine false positive for the stronger
moment--nullcone strategy: it satisfies the SIC conclusion with an
explicit linear bound.

The dehomogenized coefficient polynomial factors as
\[
\begin{aligned}
F_C(x,y)
 &=y+y^3-x-2x^2y-x^2y^3-x^3y^2\\
 &=-(xy+y^2+1)(x^2y+x-y).                              \tag{1.5}
\end{aligned}
\]
Equivalently,
\[
F=-(WZ^2+VZY+WY^2)(V^2Y+WVZ-W^2Y).                    \tag{1.6}
\]
In the Clebsch--Gordan coordinates of the
[bidegree-\((3,3)\) frontier](TWO_PAIR_SIC_BIDEGREE33_FRONTIER.md),
the only nonzero coordinates are
\[
s_6=-1,\qquad t_3=\frac12,\qquad r_0=-1.               \tag{1.7}
\]

## 2. Beta and constant-term realization

For a balanced form \(H\) of bidegree \((D,D)\), dehomogenize with
\[
x=\frac{ut}{1-t},\qquad y=u^{-1}.
\]
The beta identity gives
\[
\frac{\mathcal E_2(H)}{(D+1)!}
=\int_0^1\operatorname{CT}_u
 \left((1-t)^D H_C\left(\frac{ut}{1-t},u^{-1}\right)\right)\,dt.
\tag{2.1}
\]

For (1.1), write \(a=1-t\).  Direct substitution and (1.5) give the
two-factor Laurent polynomial
\[
\boxed{
p(t,u)=a^3F_C\left(\frac{ut}{a},u^{-1}\right)
=u^{-3}(u^2+a)(a^2-tu^2).}                             \tag{2.2}
\]
Therefore
\[
\frac{\mathcal E_2(F^m)}{(3m+1)!}
=\int_0^1\operatorname{CT}_u p(t,u)^m\,dt.             \tag{2.3}
\]

If \(m\) is odd, every exponent of \(u\) in \(p^m\) is odd, so the
constant term is zero.  Let \(m=2n\).  Setting \(z=u^2\), the required
constant term is the coefficient of \(z^{3n}\) in
\[
(z+a)^{2n}(a^2-tz)^{2n}.
\]
The beta integral then gives
\[
\begin{aligned}
\int_0^1\operatorname{CT}_u p^{2n}\,dt
&=\frac{(-1)^n(2n)!^2}{(3n+1)!}
  \sum_{r=0}^{n}\frac{(-1)^r}{r!(n-r)!}\\
&=\frac{(-1)^n(2n)!^2}{(3n+1)!\,n!}(1-1)^n=0.
\end{aligned}                                          \tag{2.4}
\]
This proves Theorem 1.1 without a finite-moment cutoff.

## 3. The Rodrigues identity

It suffices first to consider a balanced monomial multiplier
\[
Q_{eij}=W^{e-i}V^iZ^{e-j}Y^j,\qquad 0\leq i,j\leq e.
\tag{3.1}
\]
Its angular restriction is
\[
t^i(1-t)^{e-i}u^{i-j}.                                 \tag{3.2}
\]
Put
\[
d=i-j,\qquad r=\frac{m-d}{2},\qquad s=\frac{m+d}{2}.
\tag{3.3}
\]
If \(r,s\) are not nonnegative integers, phase parity already makes the
mixed contraction zero.  Otherwise the relevant phase coefficient in
(2.2) is
\[
\begin{aligned}
\operatorname{CT}_u(u^dp^m)
&=[z^{2r+s}](z+a)^m(a^2-tz)^m\\
&=\boxed{
\frac{(-1)^r a^d}{s!}
\frac{d^s}{dt^s}\left(t^m a^m\right).}                 \tag{3.4}
\end{aligned}
\]
The second equality follows directly from the Leibniz rule.  Multiplying
by (3.2) and using \(d=i-j\) gives
\[
\frac{\mathcal E_2(Q_{eij}F^m)}{(3m+e+1)!}
=\frac{(-1)^r}{s!}
\int_0^1t^i(1-t)^{e-j}
\frac{d^s}{dt^s}\left(t^m(1-t)^m\right)\,dt.            \tag{3.5}
\]

Every boundary term in \(s\) integrations by parts vanishes: before the
last integration the differentiated factor \(t^m(1-t)^m\) still vanishes
at both endpoints, since \(s\leq m\).  Hence (3.5) is, up to a nonzero
sign,
\[
\frac1{s!}\int_0^1
\frac{d^s}{dt^s}\left(t^i(1-t)^{e-j}\right)
t^m(1-t)^m\,dt.                                        \tag{3.6}
\]
The first factor has degree
\[
i+e-j=e+d.
\]
It is therefore zero when \(s>e+d\), equivalently when
\[
m>2e+d.                                                 \tag{3.7}
\]
Since \(d\leq e\), the uniform bound \(m>3e\) works for every balanced
monomial of bidegree \((e,e)\), and hence for every balanced multiplier
of that bidegree.

For completeness, let a bihomogeneous multiplier have dual degree
\(\alpha\) and coordinate degree \(\beta\).  If \(\alpha>\beta\), total
degree makes its contraction with \(F^m\) zero.  If
\(\alpha\leq\beta\), every coefficient of the output polynomial is
detected by multiplying once more by a dual monomial of degree
\(\beta-\alpha\).  The resulting scalar contraction has a balanced
multiplier of bidegree \((\beta,\beta)\), so it vanishes for
\(m>3\beta\).  Bihomogeneous decomposition proves (1.4).

## 4. Exact sparse discovery and sharpness

The form was extracted by a complete coefficient-torus census, not by a
floating-point solve.

* Supports of size at most five were already excluded exactly outside the
  one-sided cases.
* Among all \(7{,}588\) mixed six-entry supports, moments through order
  twelve exclude \(7{,}586\).  The remaining two normalized systems each
  have one rational point.  One is (1.1); the other is its Weyl/torus
  transform
  \[
  x+y+2xy^2+x^2y^3-x^3-x^3y^2
  =-(x^2-xy-1)(xy^2+x+y).                              \tag{4.1}
  \]
* Every one of the \(11{,}200\) mixed seven-entry coefficient tori is
  excluded over \(\mathbb Q\) by moments through order twelve.

Thus any degree-three SIC counterexample has coefficient support at least
eight.  The statement is about support in the displayed monomial basis;
it is not a coordinate-invariant rank lower bound.

## 5. Exact local geometry on the null-quadratic chart

Normalize the null quadratic and the opposite sextic coefficient by
\[
r_0=s_6=-1,
\]
and use the remaining unipotent stabilizer to set \(s_5=0\).  The chart
variables are
\[
s_0,s_1,s_2,s_3,s_4,t_0,t_1,t_2,t_3,t_4.              \tag{5.1}
\]
At the point (1.7), the \(9\times9\) Jacobian of
\(\mu_2,\ldots,\mu_{10}\), using every column except \(t_4\), has
determinant
\[
3445505947738252325099075904000\ne0.                   \tag{5.2}
\]
The formal implicit-function theorem therefore gives a unique formal
curve through the point with parameter \(\varepsilon=t_4\).  Exact
coefficient lifting shows that \(\mu_{11}\) restricts to
\[
\mu_{11}
=\frac{558209902860000}{44871740771}\varepsilon^5
+O(\varepsilon^6).                                     \tag{5.3}
\]
Hence the Rodrigues point is an isolated length-five local component of
\((\mu_2,\ldots,\mu_{11})\) on this normalized chart.

On the five-variable tangent slice retaining only
\((s_0,s_3,t_0,t_3,t_4)\), the reduced Gröbner basis of moments through
order ten is
\[
\begin{gathered}
17t_0-70t_4,\ s_0,\ t_4^2,\ 7t_3t_4+17s_3,\ s_3t_4,\\
4t_3^2-1,\ 68s_3t_3+7t_4,\ s_3^2.                     \tag{5.4}
\end{gathered}
\]
The quotient has length four and reduced support \(t_3=\pm\frac12\),
with every other retained coordinate zero.  These are the two normalized
copies of the same Weyl/torus orbit.  Equations (5.2)--(5.4) prove local
isolation; they do not classify the whole ten-variable chart.

## 6. Relation to the current literature

[Long's Laurent counterexample](https://arxiv.org/abs/2607.19012) has the
same interval--torus flavor, but its detecting multiplier is a bare
negative phase.  In the beta model (3.2), a regular balanced polynomial
multiplier of negative phase necessarily carries an endpoint factor.
Thus Long's Laurent multiplier does not directly descend to a regular
two-pair SIC multiplier in bidegree three.  The
[Müger--Tuset reduction](https://arxiv.org/abs/2410.11622) explains the
broader interval--torus setting, while the earlier
[Dings--Koelink results](https://arxiv.org/abs/1404.4215) treat restricted
families of \(SU(2)\) matrix elements rather than this full coefficient
space.

[Wilson's face-isolation theorem](https://arxiv.org/abs/2607.23887)
proves the real two-variable Gaussian Moments Conjecture by isolating a
Newton face with \(p\)-adic factorial valuations.  Its moment weight has
one radial factorial.  Formula (2.4) for the present balanced problem has
the complementary pair \(I!(3m-I)!\); moving off a face raises one factor
while lowering the other, so the valuation separation used there does not
transfer directly.  None of these results classifies the remaining
bidegree-\((3,3)\) two-pair strata.

## 7. Reproduction

Run

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rodrigues_survivor.py

.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_null_quadratic_s6.py \
  --orders 2,3,4,5,6,7,8,9,10,11 --skip-solver \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_null_quadratic_s6_local.json
```

The first checker verifies the determinant, factorization, direct pure
moments through order thirty, the beta formula against direct mixed
contraction, the Rodrigues identity on a finite symbolic audit range, and
the uniform mixed cutoff.  The second checker constructs the exact chart
moments, verifies (5.2)--(5.3), and asks Singular to reproduce (5.4).

The complete size-six and size-seven sparse censuses are in
[`two_pair_sic_bidegree33_sparse_six_support_screen.json`](../artifacts/generated-results/two_pair_sic_bidegree33_sparse_six_support_screen.json)
and
[`two_pair_sic_bidegree33_sparse_support7_screen.json`](../artifacts/generated-results/two_pair_sic_bidegree33_sparse_support7_screen.json).
Their sharded reproduction commands are recorded in
[`REPRODUCE.md`](../REPRODUCE.md).
