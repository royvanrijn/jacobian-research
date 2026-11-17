# Split and connected full fibers in one Keller map

## 1. The fixed-map theorem

For every integer \(N\ge3\), there is one explicit polynomial Keller map

\[
F_N:\mathbb A^3_{\mathbb Q}\longrightarrow\mathbb A^3_{\mathbb Q},
\qquad \det DF_N=1,
\]

with the following two properties.

1. One rational target has full fiber
   \[
   F_N^{-1}(y_{\mathrm{split}})
   \simeq\operatorname{Spec}\mathbb Q^N.
   \]
2. Infinitely many rational targets have connected full fibers of degree
   \(N\):
   \[
   F_N^{-1}(y_v)
   \simeq\operatorname{Spec}L_v,
   \qquad [L_v:\mathbb Q]=N.
   \]

Thus one fixed Keller cover in every geometric degree at least three carries
both a completely split fiber and infinitely many field fibers.  The
connected targets may all be taken on one rational affine line through the
split target.

The theorem concerns the simultaneous occurrence of sharply different
fiber types.  It does **not** assert that two arbitrary rank-\(N\) finite
étale algebras can be prescribed in one map, nor that the fields \(L_v\) are
pairwise nonisomorphic.

## 2. The explicit all-degree map

Put

\[
G_N(S)=\prod_{j=0}^{N-1}(S-j)
      =g_1S+g_2S^2+\cdots+g_NS^N.
\]

Its linear coefficient is

\[
g_1=(-1)^{N-1}(N-1)!\ne0.
\]

The cubic coefficient is, up to sign, the positive elementary symmetric
sum \(e_{N-3}(1,\ldots,N-1)\), with \(e_0=1\) when \(N=3\).  Hence
\(g_3\ne0\), and \(G_N\) is an admissible root-engineered quadratic-gauge
seed.

For completeness, the map is given directly as follows.  Set

\[
t=1+xy,\qquad
q=t^2z+\frac{g_1}{g_3}y^2(1+3t).
\]

Define \(F_N=(\Pi,U,C)\) by

\[
\begin{aligned}
\Pi={}&tq,\\
U={}&-\frac12\left(
y+3\frac{g_3}{g_1}xq
+2\frac{g_2}{g_1}tq
+\sum_{k=4}^Nk\frac{g_k}{g_1}t^2x^{k-2}q^k
\right),\\
C={}&x(5-3t)-\frac{g_3}{g_1}x^3z
-\sum_{k=4}^N(k-2)\frac{g_k}{g_1}(xq)^k.
\end{aligned}
\]

The root-engineered quadratic-gauge identity gives determinant \(-2\)
before the output replacement \(U=-B/2\).  Therefore

\[
\boxed{\det DF_N=1.}
\]

In these displayed target coordinates, the inverse polynomial is

\[
\boxed{
E_{\pi,u,c}(S)
=g_1S+\pi(g_2S^2+g_3S^3)
+\sum_{k=4}^Ng_k\pi^kS^k
+g_1uS^2-\frac{g_1}{2}c.
}
\]

## 3. The split fiber

At

\[
y_{\mathrm{split}}=(1,0,0)
\]

the inverse identity is

\[
E_{1,0,0}(S)=G_N(S)=\prod_{j=0}^{N-1}(S-j).
\]

All roots are simple and rational.  Since \(\pi=1\), the
scheme-theoretic reconstruction theorem applies and identifies the entire
fiber:

\[
\boxed{
F_N^{-1}(1,0,0)
\simeq\operatorname{Spec}\mathbb Q[S]/(G_N)
\simeq\operatorname{Spec}\mathbb Q^N.
}
\]

## 4. Infinitely many field fibers on one line

Restrict the target to

\[
\ell=\{(1,0,v):v\in\mathbb A^1\}.
\]

The inverse polynomial becomes

\[
\boxed{
P_{N,v}(S)=G_N(S)-\frac{g_1}{2}v.
}
\]

This polynomial is irreducible over \(\mathbb Q(v)\).  Indeed, work first
in \(K[S,v]\), where \(K=\mathbb Q\).  If

\[
G_N(S)-\frac{g_1}{2}v=A(S,v)B(S,v),
\]

then its degree one in \(v\) forces one factor, say \(A\), to have
\(v\)-degree zero.  Comparing the coefficient of \(v\) shows that \(A\)
divides the nonzero constant \(-g_1/2\) in \(K[S]\), so \(A\) is a unit.
Gauss's lemma gives irreducibility in \(\mathbb Q(v)[S]\).  The same proof
after extending the constant field shows geometric irreducibility.

Hilbert irreducibility now supplies an infinite Hilbert subset
\(\mathcal H_N\subset\mathbb Q\) such that

\[
P_{N,v}(S)\ \text{is irreducible over }\mathbb Q
\qquad(v\in\mathcal H_N).
\]

Every such polynomial is separable because \(\mathbb Q\) has characteristic
zero.  Again \(\pi=1\), so scheme-theoretic reconstruction gives

\[
\boxed{
F_N^{-1}(1,0,v)
\simeq\operatorname{Spec}\mathbb Q[S]/(P_{N,v}),
\qquad v\in\mathcal H_N.
}
\]

The quotient is a degree-\(N\) field.  These are full fibers because the
inverse polynomial has degree \(N=\operatorname{gdeg}(F_N)\), and the
reconstruction theorem accounts for the complete fiber scheme.

## 5. A tiny explicit quintic certificate

For \(N=5\),

\[
G_5(S)=S(S-1)(S-2)(S-3)(S-4)
=S^5-10S^4+35S^3-50S^2+24S.
\]

The map specializes to

\[
t=1+xy,\qquad
q=t^2z+\frac{24}{35}y^2(1+3t),
\]

\[
\begin{aligned}
\Pi={}&tq,\\
U={}&-\frac12y-\frac{35}{16}xq+\frac{25}{12}tq
       +\frac56t^2x^2q^4-\frac5{48}t^2x^3q^5,\\
C={}&x(5-3t)-\frac{35}{24}x^3z
       +\frac56(xq)^4-\frac18(xq)^5.
\end{aligned}
\]

The two targets can be taken on the theorem's line:

\[
y_{\mathrm{split}}=(1,0,0),\qquad
y_{\mathrm{field}}=(1,0,1).
\]

Their inverse polynomials are

\[
\boxed{
E_{1,0,0}(S)=S(S-1)(S-2)(S-3)(S-4)
}
\]

and

\[
\boxed{
E_{1,0,1}(S)
=S^5-10S^4+35S^3-50S^2+24S-12.
}
\]

Modulo \(5\), the second polynomial is

\[
S^5-S-2.
\]

This Artin--Schreier polynomial is irreducible over \(\mathbb F_5\).
If \(\alpha^5-\alpha=2\), then Frobenius sends
\(\alpha\) to \(\alpha+2\), so its Frobenius orbit has exactly five elements.
Thus its minimal polynomial has degree five.  Equivalently, the exact checker
verifies the degree-five Rabin criterion.  Gauss's lemma makes the displayed
quintic irreducible over \(\mathbb Q\).

Consequently

\[
F_5^{-1}(1,0,0)\simeq\operatorname{Spec}\mathbb Q^5,
\]

whereas

\[
F_5^{-1}(1,0,1)
\simeq
\operatorname{Spec}
\mathbb Q[S]/(S^5-10S^4+35S^3-50S^2+24S-12)
\]

is the spectrum of a degree-five field.

## 6. Constraint on arbitrary simultaneous prescription

The same inverse formula explains why this theorem does not yet solve the
arbitrary pair problem.  Suppose two fixed displayed polynomials are

\[
P_i(S)=\sum_{k=0}^Na_{i,k}S^k
\]

and both are required to equal inverse polynomials for one fixed seed.
For nonzero \(\pi_1,\pi_2\), put \(\rho=\pi_2/\pi_1\).  Exact coefficient
comparison forces

\[
a_{1,1}=a_{2,1}=g_1,
\]

\[
a_{2,3}=\rho a_{1,3},\qquad
a_{2,k}=\rho^k a_{1,k}\quad(4\le k\le N).
\]

Only the constant and quadratic coefficients move freely through \(c\) and
\(u\).  Increasing the seed degree therefore does not make the two-identity
system underdetermined inside this gauge: every additional seed coefficient
brings two equations but only one new coefficient.

Arbitrary simultaneous prescription would require extra target-controlled
coefficient directions, a systematic use of nonlinear primitive generators
for the two étale algebras, or a different Keller construction.

## 7. Reproduction

Run

```bash
.venv/bin/python scripts/verify_split_and_connected_full_fibers.py
```

The checker uses exact rational arithmetic and verifies:

- the displayed quintic map has Jacobian one;
- its inverse polynomial vanishes identically on the source chart;
- the split and field inverse identities;
- squarefreeness and irreducibility modulo \(5\);
- the all-degree seed, target-line, and generic-irreducibility formulas
  through degree twelve as a bounded regression.

The general scheme-level reconstruction is independently replayed by

```bash
.venv/bin/python scripts/verify_finite_etale_keller_fibers.py
```

The bounded all-degree loop is a regression, not the proof of the universal
statement.  The proof is the coefficient argument and Hilbert
irreducibility in Section 4.
