# An eight-term counterexample in three contraction pairs

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
 f&=\tau^3(t+z)\bigl(wt^3-vy(t+y)^2\bigr),\\
 g&=z.
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
>  [t]\mathcal E_3(gf^m)=(3m+1)!\,m!.}                        \tag{1.3}
> \]

Thus \(\mathcal M_3\) is not a Mathieu--Zhao subspace.  The polynomial \(f\)
has eight expanded terms, coefficients in
\(\{1,-1,-2\}\), zeta-degree four, \(z\)-degree four, and ordinary total
degree eight.  The multiplier \(g\) is one \(z\)-variable of bidegree
\((0,1)\).

Since the one-pair Special Image Conjecture is known, the minimum pair
dimension now satisfies

\[
 \boxed{2\leq r_{\rm SIC}\leq3.}                              \tag{1.4}
\]

No assertion about \(\operatorname{SIC}(2)\) is made.

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

Every monomial of \(f\) has zeta-degree four and \(z\)-degree four.
Consequently every monomial of \(\mathcal E_3(f^m)\) has total residual
degree zero; the contraction is a scalar.  Since \(g\) has bidegree
\((0,1)\), the mixed contraction has total residual degree one.  We only
need its coefficient of \(t\).

## 3. The two-pair circular identity

First omit the pair \((\tau,t)\) and consider

\[
 P=(1+Z)\bigl(W-VY(1+Y)^2\bigr),\qquad Q=Z.                   \tag{3.1}
\]

Let the centered circular contraction functional be

\[
 \mathcal F(W^aZ^bV^cY^d)
 =\delta_{a,b}\delta_{c,d}\,a!\,c!.                           \tag{3.2}
\]

Expanding according to the number \(k\) of
\(-VY(1+Y)^2\) factors gives

\[
 P^m=\sum_{k=0}^m(-1)^k\binom mk
 W^{m-k}V^k(1+Z)^mY^k(1+Y)^{2k}.                             \tag{3.3}
\]

For a nonzero contraction, the first circular pair selects the coefficient
of \(Z^{m-k}\), namely \(\binom mk\).  The second pair selects the constant
term of \((1+Y)^{2k}\).  Hence

\[
\begin{aligned}
 \mathcal F(P^m)
 &=\sum_{k=0}^m(-1)^k\binom mk^2(m-k)!k!\\
 &=m!\sum_{k=0}^m(-1)^k\binom mk
 =0.                                                         \tag{3.4}
\end{aligned}
\]

For \(QP^m=ZP^m\), the first pair instead selects
\([Z^{m-k-1}](1+Z)^m=\binom m{k+1}\).  Therefore

\[
\begin{aligned}
 \mathcal F(ZP^m)
 &=\sum_{k=0}^{m-1}(-1)^k
   \binom mk\binom m{k+1}(m-k)!k!\\
 &=m!\sum_{k=0}^{m-1}(-1)^k\binom m{k+1}
 =m!.                                                        \tag{3.5}
\end{aligned}
\]

These are all-order binomial identities.

## 4. One-pair bihomogenization

The polynomial \(P\) has \(W,V\)-degree one and \(Z,Y\)-degree at most four.
Its total-degree-four homogenization in a new variable \(t\) is

\[
 \widetilde P
 =(t+z)\bigl(wt^3-vy(t+y)^2\bigr).                            \tag{4.1}
\]

Multiplication by \(\tau^3\) makes both sides have degree four:

\[
 f=\tau^3\widetilde P.
\]

In a term of \(P^m\) that survives circular contraction, the total
\((Z,Y)\)-degree is \(m\), matching the total \((W,V)\)-degree \(m\).
After homogenization, its \(t\)-degree is therefore \(4m-m=3m\), exactly
matching the exponent of \(\tau\).  The new contraction pair contributes
\((3m)!\).  Equation (3.4) gives

\[
 \mathcal E_3(f^m)=(3m)!\mathcal F(P^m)=0.                   \tag{4.2}
\]

For \(ZP^m\), a surviving term of \(P^m\) has \((Z,Y)\)-degree \(m-1\).
The homogenized term has \(t\)-degree \(4m-(m-1)=3m+1\).  Contracting the
\(\tau^{3m}\) factor leaves one power of \(t\).  Thus (3.5) gives

\[
 [t]\mathcal E_3(gf^m)
 =(3m+1)!\mathcal F(ZP^m)
 =(3m+1)!\,m!\ne0.                                          \tag{4.3}
\]

Combining (2.2), (4.2), and (4.3) proves Theorem 1.1.

## 5. Provenance and novelty status

The centered two-pair polynomial (3.1) is the \(H(z)=z,\lambda=1\)
specialization of the repository's
[weighted-seed Gaussian bridge](WEIGHTED_GAUSSIAN_BRIDGE.md).  The
one-pair bihomogenization (4.1)--(4.3), its three-pair SIC consequence, and
the compact formula (1.2) are repository-derived.

No matching three-pair formula was found in the targeted public searches
performed when this note was written.  This is not a priority claim:
external review and a dedicated literature search remain necessary.

## Reproduction

Run

```bash
python3 scripts/verify_three_pair_image_mathieu_counterexample.py
```

The dependency-free checker expands \(f\), performs exact sparse
contractions through \(m=10\), verifies (1.3) on that prefix, and separately
checks the two all-order binomial identities through \(m=99\).  The proof
above, rather than either finite cutoff, establishes the theorem for every
\(m\).
