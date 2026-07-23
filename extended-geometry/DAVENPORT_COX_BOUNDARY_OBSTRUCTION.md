# Boundary lattice and affine-space obstruction for the Davenport Sunada pair

The global Davenport--Sunada covers admit determinant-one Cox suspensions
over one common target, but those morphisms live on affine boundary
complements.  This note computes the complete height-one boundary pullback
and proves that neither affine stabilization nor a coordinate-preserving
polynomial suspension can turn that construction into a polynomial Keller
map of affine space.

Work over

\[
 K=\mathbb Q(a),\qquad a^2+a+2=0.
\]

Use the notation \(g_T(Y),h_T(Z),\Delta(T,u)\) from the
[global Sunada theorem](GLOBAL_SUNADA_KELLER_COVERS.md).

## 1. Exact boundary factorization

On the point cover \(u=g_T(Y)\), exact factorization in \(K(T)[Y]\) gives

\[
 \boxed{\Delta(T,g_T(Y))=E_{3,g}(T,Y)E_{6,g}(T,Y)J_g(T,Y)^2,} \tag{1.1}
\]

where

\[
 J_g(T,Y)=g_T'(Y),                                    \tag{1.2}
\]

\[
 E_{3,g}(T,Y)
 =Y^3+(1+3a)TY+(7+5a)T,                              \tag{1.3}
\]

and

\[
\begin{aligned}
E_{6,g}(T,Y)={}&Y^6+(10+8a)TY^4+(6+8a)TY^3\\
&+(-19+28a)T^2Y^2+(-18+36a)T^2Y\\
&-27T^2+(-88-36a)T^3.
\end{aligned}                                         \tag{1.4}
\]

All three displayed factors are monic and irreducible in \(K(T)[Y]\).
Their coefficients are primitive over \(K[T]\), so Gauss's lemma makes them
irreducible in \(K[T,Y]\).  Conjugating \(a\) and replacing \(Y\) by \(Z\)
gives

\[
 \Delta(T,h_T(Z))=E_{3,h}E_{6,h}J_h^2,\qquad J_h=h_T'. \tag{1.5}
\]

The multiplicity vector

\[
 \boxed{(1,1,2)}                                      \tag{1.6}
\]

has the expected branch-cycle meaning.  Over each of the three finite
branch values, two critical points have ramification index two.  The
derivative therefore supplies the doubled height-one component; the other
degree-three and degree-six components consist of unramified sheets lying
over the branch divisor.

## 2. Units and Picard groups

The common target surface and the point-cover source are

\[
\begin{aligned}
B&=\operatorname {Spec}K[T,u,\Delta^{-1}],\\
X_g&=\operatorname {Spec}
K[T,Y,(E_{3,g}E_{6,g}J_g)^{-1}].
\end{aligned}                                         \tag{2.1}
\]

Both coordinate rings are localizations of polynomial UFDs.  Hence

\[
\operatorname {Pic}(B)=\operatorname {Pic}(X_g)=0,     \tag{2.2}
\]

but their unit lattices are different:

\[
\mathcal O(B)^*/K^*=\mathbb Z[\Delta],                \tag{2.3}
\]

\[
\mathcal O(X_g)^*/K^*
=\mathbb Z[E_{3,g}]\oplus
 \mathbb Z[E_{6,g}]\oplus
 \mathbb Z[J_g].                                      \tag{2.4}
\]

The pullback of the target generator is

\[
\pi_g^*[\Delta]=(1,1,2),                              \tag{2.5}
\]

while the plane-cover Jacobian has unit vector

\[
[J_g]=(0,0,1).                                        \tag{2.6}
\]

The primitive Cox coordinate \(Z=z/J_g\) cancels (2.6), exactly as in the
global Sunada construction.  It does not remove the boundary unit lattice.

There is, however, no finite Smith obstruction.  The target-pullback and
Jacobian rows form the matrix

\[
\begin{pmatrix}
1&1&2\\
0&0&1
\end{pmatrix}.                                        \tag{2.7}
\]

Its maximal minors have gcd one, so its row lattice is saturated.  Appending
one primitive boundary character gives, for example,

\[
\boxed{
\begin{pmatrix}
1&1&2\\
0&0&1\\
1&0&0
\end{pmatrix}
\in\operatorname {GL}_3(\mathbb Z).
}                                                     \tag{2.8}
\]

Thus one further primitive character beyond the derivative-unit suspension
is sufficient at the integral-lattice level.  Unlike a general
three-factor degree ledger with nontrivial Smith index, this construction
does not require a finite root cover or stack merely to complete its
exponent matrix.  What remains is to realize the last row of (2.8) by a
polynomial affine modification.

The most obvious realization does not work.  The primitive row in (2.8)
selects \(E_{3,g}\), and (1.3) is linear in \(T\):

\[
E_{3,g}=Y^3+T L(Y),\qquad
L(Y)=(1+3a)Y+(7+5a).                                  \tag{2.9}
\]

Hence the tentative coordinate change

\[
(T,Y)\longmapsto(E_{3,g},Y)
\]

has Jacobian \(L(Y)\).  Exact gcd calculations give

\[
\gcd(L,E_{3,g})=\gcd(L,E_{6,g})=\gcd(L,J_g)=1.         \tag{2.10}
\]

Thus \(L=0\) is a genuinely new boundary divisor, not a combination of the
three ledger columns.  The linear-in-\(T\) chart realizes the desired
primitive character only by recreating the extra-divisor obstruction.  A
successful affine modification must absorb \(L\) simultaneously or use a
nonlinear chart.

## 3. Stable no-straightening theorem

For every reduced domain \(R\),

\[
R[t_1,\ldots,t_m]^*=R^*.                              \tag{3.1}
\]

Consequently affine stabilization preserves the unit ranks in (2.3)--(2.4).
In particular, for every \(m,n\ge0\),

\[
B\times\mathbb A^m\not\simeq\mathbb A^{m+2},\qquad
X_g\times\mathbb A^n\not\simeq\mathbb A^{n+2}.         \tag{3.2}
\]

The same holds for \(X_h\).  More strongly, the source and target
complements cannot become isomorphic after stabilization because their unit
ranks are three and one.

Thus the direct Cox--Sunada morphisms cannot be converted into affine-space
Keller maps by polynomial coordinate changes, even after adding arbitrarily
many identity variables.  Any absolute construction must change the
varieties by filling or modifying boundary divisors; stable straightening is
not an available shortcut.

## 4. Coordinate-preserving polynomial suspension no-go

There is a second obstruction which does not use units.  Suppose a
polynomial map of affine spaces retains the Davenport core as its first two
target coordinates:

\[
(T,Y,\mathbf z)\longmapsto
\bigl(T,g_T(Y),R_1(T,Y,\mathbf z),\ldots,
R_m(T,Y,\mathbf z)\bigr).                             \tag{4.1}
\]

Its Jacobian matrix is block triangular, so

\[
\det DF
=g_T'(Y)\det\left(\frac{\partial(R_1,\ldots,R_m)}
                        {\partial(z_1,\ldots,z_m)}\right). \tag{4.2}
\]

The right side cannot be a nonzero constant because \(g_T'\) is a
nonconstant nonunit in the polynomial ring.  Therefore:

> No polynomial suspension which preserves the Davenport cover coordinates
> can be a Keller map of affine space, in any number of added variables.

The rational Cox coordinate \(z/J_g\) is not an avoidable presentation
artifact; division by the derivative is forced in this entire
coordinate-preserving class.

## 5. The remaining construction problem

The global group theory, common cover, zeta equality, and Jacobian-unit
cancellation are complete.  Equations (1.1) and (2.5) now isolate the
absolute problem to a concrete three-column boundary ledger:

\[
\begin{array}{c|ccc}
 &E_3&E_6&J\\ \hline
\pi^*\Delta&1&1&2\\
\operatorname {Jac}(\pi)&0&0&1.
\end{array}                                           \tag{5.1}
\]

An affine-space construction must use a nontrivial source and target
modification which:

1. fills the target divisor \(\Delta=0\);
2. controls all three source divisors in (1.1);
3. cancels the Jacobian row without inverting a polynomial coordinate; and
4. preserves the point/line \(GL_3(\mathbb F_2)\) Gassmann closure on the
   surviving common open.

The unimodular completion (2.8) proves that this is a geometric realization
problem, not an integral-lattice obstruction.  At most one additional
primitive boundary character is missing after the direct derivative-unit
suspension, but (2.9)--(2.10) show that its simplest coordinate realization
adds a fourth divisor.

The existing weighted lift achieves an analogous cancellation only after
specializing to a tangent pencil and adjoining tangent marks.  The next
absolute search should therefore be organized around affine modifications
of the explicit ledger (5.1), not around further stable coordinate changes
or prime-interval arguments for the cancellation parameter polynomial.

Adjoining square roots on both sides does solve the polynomial determinant
ledger without a reciprocal coordinate: the
[oriented Davenport Cox maps](ORIENTED_DAVENPORT_COX_MAPS.md) use
`W^2=E_3E_6`, `D^2=Delta`, and `D=JW` to obtain residue Jacobian one and the
reduced boundary vector `(1,1,1)`.  Their source and target remain
non-affine-space Cox charts.

## 6. Verification

Run

```bash
.venv/bin/python scripts/verify_davenport_cox_boundary.py
```

The checker proves the factor pattern \(3\cdot6\cdot6^2\), identifies the
doubled factor with the derivative, records the unit ranks and unimodular
completion (2.8), certifies the new coprime divisor (2.10), and verifies the
block-determinant identity (4.2) with three arbitrary suspension variables.
