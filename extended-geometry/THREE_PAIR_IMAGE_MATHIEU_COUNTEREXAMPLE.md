# A four-term counterexample in three contraction pairs and four real Gaussians

## 1. Explicit statement

Work over a characteristic-zero field and use the three contraction pairs

\[
 (\tau,t),\qquad (w,z),\qquad (v,y).
\]

Put

\[
 \mathcal M_3
 =(\partial_t-\tau)\mathbb C[\tau,w,v,t,z,y]
  +(\partial_z-w)\mathbb C[\tau,w,v,t,z,y]
  +(\partial_y-v)\mathbb C[\tau,w,v,t,z,y].
                                                               \tag{1.1}
\]

Define

\[
 \boxed{
 \begin{aligned}
 f&=\tau(t-y)(wz+vt),\\
 g&=y.
 \end{aligned}}                                               \tag{1.2}
\]

> **Theorem 1.1.** For every \(m\geq1\),
> \[
>  f^m\in\mathcal M_3,\qquad
>  gf^m\notin\mathcal M_3.
> \]
> More precisely, for the contraction map \(\mathcal E_3\),
> \[
>  \boxed{
>  \mathcal E_3(f^m)=0,\qquad
>  [t]\mathcal E_3(gf^m)=(-1)^{m-1}(m+1)!\,m!.}              \tag{1.3}
> \]

Thus \(\mathcal M_3\) is not a Mathieu--Zhao subspace.  The polynomial \(f\)
has four expanded terms, coefficients in \(\{1,-1\}\), zeta-degree two,
\(z\)-degree two, and ordinary total degree four.  The multiplier \(g\) has
bidegree \((0,1)\).

Since the one-pair Special Image Conjecture is known, the minimum pair
dimension satisfies

\[
 \boxed{2\leq r_{\rm SIC}\leq3.}                              \tag{1.4}
\]

No assertion about \(\operatorname{SIC}(2)\) or term-count minimality is made.

## 2. Image membership as contraction

For multi-indices \(\alpha,\beta\in\mathbb N^3\), set

\[
 \mathcal E_3(\zeta^\alpha z^\beta)
 =\partial_z^\alpha z^\beta.
                                                               \tag{2.1}
\]

Zhao's image-kernel identity gives

\[
 \mathcal M_3=\ker\mathcal E_3.                               \tag{2.2}
\]

Every monomial of \(f\) has zeta-degree two and \(z\)-degree two.
Consequently every monomial of \(\mathcal E_3(f^m)\) has total residual
degree zero; the contraction is a scalar.  Since \(g\) has bidegree
\((0,1)\), the mixed contraction has total residual degree one.  It is enough
to compute its coefficient of \(t\).

## 3. Long's two-pair circular seed

First omit the pair \((\tau,t)\) and consider

\[
 P=(1-Y)(WZ+V),\qquad Q=Y.                                   \tag{3.1}
\]

Let the centered circular contraction functional be

\[
 \mathcal F(W^aZ^bV^cY^d)
 =\delta_{a,b}\delta_{c,d}\,a!\,c!.                           \tag{3.2}
\]

Expanding the second factor gives

\[
 P^m=\sum_{k=0}^m\binom mk
 W^{m-k}Z^{m-k}V^k(1-Y)^m.                                  \tag{3.3}
\]

For a nonzero contraction, the \(V,Y\) pair selects
\([Y^k](1-Y)^m=(-1)^k\binom mk\).  Hence

\[
\begin{aligned}
 \mathcal F(P^m)
 &=\sum_{k=0}^m(-1)^k\binom mk^2(m-k)!k!\\
 &=m!\sum_{k=0}^m(-1)^k\binom mk
 =0.                                                         \tag{3.4}
\end{aligned}
\]

For \(QP^m=YP^m\), the same pair instead selects
\([Y^k]Y(1-Y)^m=(-1)^{k-1}\binom m{k-1}\).  Therefore

\[
\begin{aligned}
 \mathcal F(YP^m)
 &=m!\sum_{k=1}^m(-1)^{k-1}\binom m{k-1}\\
 &=(-1)^{m-1}m!.                                             \tag{3.5}
\end{aligned}
\]

These are all-order binomial identities.

### 3.1 Four-real Gaussian consequence

Let \(X_1,Y_1,X_2,Y_2\) be independent standard real Gaussian variables and
put

\[
 Z_j=\frac{X_j+iY_j}{\sqrt2},\qquad
 W_j=\frac{X_j-iY_j}{\sqrt2}=\overline{Z_j}
 \qquad(j=1,2).
                                                               \tag{3.6}
\]

Define

\[
 \boxed{
 \begin{aligned}
 P_4^{\rm small}
   &=(1-Z_2)(W_1Z_1+W_2)\\
   &=W_1Z_1+W_2-W_1Z_1Z_2-W_2Z_2,\\
 Q_4^{\rm small}&=Z_2.
 \end{aligned}}                                               \tag{3.7}
\]

> **Corollary 3.1.** For every \(m\geq1\),
> \[
>  \boxed{
>  \mathbb E\!\left((P_4^{\rm small})^m\right)=0,\qquad
>  \mathbb E\!\left(Q_4^{\rm small}(P_4^{\rm small})^m\right)
>    =(-1)^{m-1}m!\ne0.}                                     \tag{3.8}
> \]
> Hence \(P_4^{\rm small},Q_4^{\rm small}\) form a four-term cubic
> counterexample to \(\operatorname{GMC}(4)\).

Indeed, independence and circular symmetry give

\[
 \mathbb E(W_1^aZ_1^bW_2^cZ_2^d)
 =\delta_{a,b}\delta_{c,d}\,a!\,c!,
                                                               \tag{3.9}
\]

so Gaussian expectation is exactly the functional \(\mathcal F\) in
(3.2).  Equations (3.4) and (3.5) therefore prove (3.8).  More generally,
for \(A\in\mathbb C[z]\),

\[
 \mathbb E\!\left(A(Z_2)(P_4^{\rm small})^m\right)
 =m!\,[z^m]A(z)(1-z)^{m-1}.                                 \tag{3.10}
\]

The term count in (3.7) is in the same natural circular coordinates used in
Christopher D. Long's Gaussian-moment paper.  His displayed four-real
cubic

\[
 (1+Z_2)\bigl(W_1(1-Z_1)+W_2\bigr)
\]

has six terms and mixed moments \(m!\).  Thus (3.7) reduces the displayed
support from six terms to four at the same real dimension and total degree;
its alternating mixed moments have generating function \(t/(1+t)\) instead
of \(t/(1-t)\).  Long's five-term quartic in three real variables remains
lower-dimensional.  No global term-minimality or literature-priority claim
is made here.

Long's Gaussian paper is
[*Small Counterexamples to the Gaussian Moments
Conjecture*](https://arxiv.org/abs/2607.18186), arXiv:2607.18186v1.

## 4. One-pair bihomogenization

The polynomial \(P\) has \(W,V\)-degree one and \(Z,Y\)-degree at most two.
Its total-degree-two homogenization in a new variable \(t\) is

\[
 \widetilde P
 =t^2P\left(w,\frac zt,v,\frac yt\right)
 =(t-y)(wz+vt).                                               \tag{4.1}
\]

Multiplication by \(\tau\) makes both sides have degree two:

\[
 f=\tau\widetilde P.
\]

In a term of \(P^m\) that survives circular contraction, the total
\((Z,Y)\)-degree is \(m\), matching the total \((W,V)\)-degree \(m\).
After homogenization, its \(t\)-degree is \(2m-m=m\), exactly matching the
exponent of \(\tau\).  The new contraction pair contributes \(m!\).
Equation (3.4) gives

\[
 \mathcal E_3(f^m)=m!\mathcal F(P^m)=0.                       \tag{4.2}
\]

For \(YP^m\), a surviving term of \(P^m\) has \((Z,Y)\)-degree \(m-1\).
The homogenized term has \(t\)-degree \(2m-(m-1)=m+1\).  Contracting the
\(\tau^m\) factor leaves one power of \(t\).  Thus (3.5) gives

\[
 [t]\mathcal E_3(gf^m)
 =(m+1)!\mathcal F(YP^m)
 =(-1)^{m-1}(m+1)!\,m!\ne0.                                 \tag{4.3}
\]

Combining (2.2), (4.2), and (4.3) proves Theorem 1.1.

## 5. Relation with the \(xz\)- and \(SU(2)\) counterexamples

Write the coordinate functions on \(SU(2)\) as
\(\left(\begin{smallmatrix}a&c\\b&d\end{smallmatrix}\right)\).  Under

\[
 W=a,\qquad Z=d,\qquad V=b,\qquad Y=-c,
\]

the seed and multiplier become

\[
 P\longmapsto(1+c)(ad+b),\qquad Q\longmapsto-c.
\]

This is exactly Christopher D. Long's \(SU(2)\) counterexample.  His
Müger--Tuset substitution further maps it to

\[
 (1-z^{-1})\bigl((1-x)+xz\bigr),
\]

the three-term counterexample to \(xz(1,1)\).  Thus the seed and its
\(xz\)- and \(SU(2)\)-interpretations are due to Long.  The degree-two
one-pair bihomogenization (4.1)--(4.3), its three-pair SIC consequence, and
the compact formula (1.2) are repository-derived.

Long's paper is
[*Counterexamples to the \(xz\)-Conjecture and the Mathieu Conjecture for
\(SU(2)\)*](https://arxiv.org/abs/2607.19012), arXiv:2607.19012v1.

## Reproduction

Run

```bash
python3 scripts/verify_three_pair_image_mathieu_counterexample.py
```

The dependency-free checker expands \(f\), performs exact sparse
contractions through \(m=10\), verifies (1.3) on that prefix, and separately
checks the two all-order binomial identities through \(m=99\).  It also
performs an independent four-real Gaussian Wick replay of (3.8) through
\(m=10\).  The proofs above, rather than either finite cutoff, establish the
identities for every \(m\).
