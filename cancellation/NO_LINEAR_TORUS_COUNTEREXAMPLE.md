# A Keller counterexample with no linear torus symmetry

## 1. Statement

Over any characteristic-zero field, put

\[
 t=1+xy,\qquad
 q=t^2z-y^2(1+3t)
\]

and define \(F:\mathbb A^3\to\mathbb A^3\) by

\[
\boxed{
\begin{aligned}
F_1={}&-\frac12tq,\\
F_2={}&y-3xq-tq+2t^2x^2q^4,\\
F_3={}&x(5-3t)+x^3z-(xq)^4.
\end{aligned}}
\tag{1}
\]

Then:

\[
\boxed{\det JF=1,}
\tag{2}
\]

\(F\) is not injective, and the only matrices \(A,B\in M_3(k)\)
satisfying

\[
\boxed{B\,F(x)=JF(x)\,A x}
\tag{3}
\]

are \(A=B=0\).  In particular, \(F\) admits no nontrivial linear
\(\mathbb G_m\)-equivariance.  The same remains true after arbitrary
independent linear changes of source and target.

For the explicit map, a stronger affine-linear calculation also holds.  If
\(a,b\in k^3\) and

\[
BF(x)+b=JF(x)(Ax+a),
\tag{3b}
\]

then \(A=B=0\) and \(a=b=0\).  Consequently the example has no nontrivial
affine-linear \(\mathbb G_m\)-equivariance, even after independent affine
changes of source and target.

More generally, let

\[
G(S)=g_1S+g_2S^2+g_3S^3+g_4S^4,\qquad g_1g_3g_4\ne0,
\]

and normalize the associated quadratic-gauge map to determinant one by
multiplying its first target coordinate by \(-1/2\).  The same conclusion
holds for every such quartic seed.  For the 18 rational coefficient rows
selected below, before rowwise denominator clearing, the universal
determinant is

\[
\boxed{\frac{10935}{4}\frac{g_4^6}{g_1^6}\ne0.}
\tag{3a}
\]

Thus the loss of linear torus symmetry is uniform on the admissible quartic
family, not merely generic.

This is deliberately a linear statement.  It does not claim that \(F\) is
not polynomially left--right equivalent to some other torus-equivariant
map.

## 2. Construction and collision

Formula (1) is the normalized root-engineered quadratic-gauge map associated
with

\[
G(S)=S(S-1)(S+1)(S-2)
    =2S-S^2-2S^3+S^4.
\tag{4}
\]

The unscaled quadratic-gauge formula has determinant \(-2\).  Multiplying
its first target coordinate by \(-1/2\) gives (2).

The complete fiber over \((-1/2,0,0)\) contains the four distinct rational
points

\[
\begin{aligned}
p_0&=(0,1,5),\\
p_1&=(-1,2,-9),\\
p_{-1}&=(1/3,-4,-27),\\
p_2&=(2/3,-1,45).
\end{aligned}
\tag{5}
\]

Direct substitution gives

\[
F(p_0)=F(p_1)=F(p_{-1})=F(p_2)=(-1/2,0,0).
\tag{6}
\]

Thus (1) is a determinant-one Keller counterexample.

## 3. The finite linear-symmetry certificate

Write the unknown matrices in row-major order and use the column order

\[
(a_{11},a_{12},a_{13},a_{21},a_{22},a_{23},a_{31},a_{32},a_{33},
 b_{11},b_{12},b_{13},b_{21},b_{22},b_{23},b_{31},b_{32},b_{33}).
\tag{7}
\]

Expand the three components of

\[
BF-JF\,A(x,y,z)^T
\tag{8}
\]

and equate every monomial coefficient to zero.  The complete coefficient
matrix has size \(734\times18\).  The following 18 coefficient rows already
have full rank.  A row label \((i;\alpha,\beta,\gamma)\) means the
coefficient of \(x^\alpha y^\beta z^\gamma\) in component \(i\) of (8).

```text
(1;12,10,4)  0  0  0   0  0  0   0  0  0   0   1    0    0    0    0   0   0    0
(1;12, 8,4)  0  0  0   0  0  0   0  0  0   0   0    1    0    0    0   0   0    0
(1; 4, 3,0)  0  0  0   0  0  0   1  0  0   0   0    0    0    0    0   0   0    0
(1; 4, 2,1)  0  0  0   1  0  0   0  0  0   0   0    0    0    0    0   0   0    0
(1; 3, 4,0)  0  0  0   0  0  0   0  1  0   0   0    0    0    0    0   0   0    0
(1; 3, 3,1)  3  0  0   0  3  0   0  0  1  -1   -2   0    0    0    0   0   0    0
(1; 3, 2,2)  0  0  0   0  0  1   0  0  0   0   0    0    0    0    0   0   0    0
(1; 2, 4,1)  0  1  0   0  0  0   0  0  0   0   0    0    0    0    0   0   0    0
(1; 2, 4,0)  2  0  0   0  4  0   0  0  0  -1   -2   0    0    0    0   0   0    0
(1; 2, 3,2)  0  0  1   0  0  0   0  0  0   0   0    0    0    0    0   0   0    0
(1; 2, 2,1)  2  0  0   0  2  0   0  0  1  -1   -2   0    0    0    0   0   0    0
(2;12,10,4) 12  0  0   0 10  0   0  0  4   0   0    0    0   -1    0   0   0    0
(2;12, 8,4)  0  0  0 180  0  0   0  0  0   0   0    0    0    0    1   0   0    0
(2; 3, 3,1)  6  0  0   0  6  0   0  0  2   0   0    0   -1   -2    0   0   0    0
(2; 3, 2,1)  3  0  0   0  2  0   0  0  1   0   0    0    0   -1    0   0   0    0
(3;12,10,4)  0  0  0   0  0  0   0  0  0   0   0    0    0    0    0   0   1    0
(3;12, 8,4) 12  0  0   0  8  0   0  0  4   0   0    0    0    0    0   0   0   -1
(3; 3, 3,1)  0  0  0   0  0  0   0  0  0   0   0    0    0    0    0   1   2    0
```

Call this integer matrix \(M\).  Its exact determinant is

\[
\boxed{\det M=-5\ne0.}
\tag{9}
\]

Consequently (3) has only the zero solution.  This is stronger than
excluding semisimple solutions: there are no nonzero nilpotent or mixed
linear infinitesimal symmetries either.

There is also a short elimination proof hidden in the matrix.  Nine unit
rows first give

\[
a_{12}=a_{13}=a_{21}=a_{23}=a_{31}=a_{32}
=b_{12}=b_{13}=b_{32}=0.
\]

The rows labelled \((2;12,8,4)\) and \((3;3,3,1)\) then give
\(b_{23}=b_{31}=0\).  Put
\(\alpha=a_{11},\beta=a_{22},\gamma=a_{33}\).  Three rows give

\[
\begin{aligned}
b_{11}&=3\alpha+3\beta+\gamma,\\
b_{11}&=2\alpha+4\beta,\\
b_{11}&=2\alpha+2\beta+\gamma.
\end{aligned}
\]

Hence \(\beta=-\alpha\) and \(\gamma=-2\alpha\).  Two further rows give
\(b_{22}=-\alpha\) and \(b_{22}=-6\alpha\), so \(5\alpha=0\).  In
characteristic zero this kills \(\alpha,\beta,\gamma,b_{11},b_{22}\), and
the last two remaining rows kill \(b_{21},b_{33}\).  Thus every one of the
18 unknowns vanishes without performing a large determinant calculation.

For a general quartic seed, extracting these same labelled rows without
primitive integer rescaling gives (3a).  This proves the uniform
quartic-family statement in Section 1.  Notice that \(g_2\) and \(g_3\)
cancel from the determinant; \(g_3\ne0\) is still required for the displayed
polynomial quadratic-gauge formula.

## 4. Why independent linear coordinates do not help

Let \(S,T\in\operatorname{GL}_3(k)\) and

\[
\widetilde F(x)=T F(Sx).
\]

If

\[
\widetilde B\,\widetilde F(x)
 =J\widetilde F(x)\,\widetilde A x,
\]

then, after multiplying by \(T^{-1}\) and setting \(u=Sx\),

\[
(T^{-1}\widetilde B T)F(u)
 =JF(u)(S\widetilde A S^{-1})u.
\tag{10}
\]

Thus linear-symmetry pairs for \(F\) and \(\widetilde F\) correspond by
conjugation.  Since (9) kills every pair for \(F\), it kills every pair
after independent linear source and target changes.

Finally, differentiating any linear \(\mathbb G_m\)-equivariance gives
(3), with semisimple integral-weight generators \(A,B\).  Equation (9)
therefore excludes every nontrivial linear torus action.

Allowing constant terms gives the 24-variable equation (3b).  Its complete
coefficient matrix has size \(785\times24\) and full column rank.  The
clean-room replay verifies this rank independently over \(\mathbb Q\).
Affine coordinate changes conjugate affine-linear vector fields to
affine-linear vector fields, so the vanishing persists after independent
affine source and target changes.

## 5. Exact reproduction

Run

```bash
make verify-linear-torus-free
```

This runs two independent exact checks.

1. The SymPy checker expands (1) over \(\mathbb Q\), verifies (2) and (6),
   extracts only the displayed 18 rows, clears them to primitive integer
   rows, and verifies (9).  It avoids solving the full linear system because
   the fixed minor is already a complete certificate.  The same extraction
   over \(\mathbb Q(g_1,g_2,g_3,g_4)\) verifies the universal determinant
   (3a), while the affine specialization verifies the \(785\times24\) rank.
2. A clean-room replay uses only the Python standard library.  It implements
   sparse polynomial arithmetic over `Fraction`, reconstructs the map,
   Jacobian, collision, and all 734 coefficient rows independently, compares
   the displayed integer matrix entry by entry, and computes its determinant
   by fraction-free Bareiss elimination.  It also reconstructs the
   \(785\times24\) affine-linear system and proves its exact rank is 24.
