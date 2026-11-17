# Degree-four exclusion at the foundational mixed weights

## 1. Result

Let \(k\) be a characteristic-zero field.  Give \(\mathbb A^3\) the split
\(\mathbb G_m\)-action with weights

\[
\operatorname{wt}(x,y,z)=(1,-1,-2).
\]

The following closes the full coordinate-degree-four locus for this action.

> **Theorem.**  Every \(\mathbb G_m\)-equivariant Keller map
> \(F:\mathbb A^3_k\to\mathbb A^3_k\) of coordinate degree at most four,
> whose target weights are a permutation of \((1,-1,-2)\), is a polynomial
> automorphism.  After permuting target coordinates and equivariantly
> normalizing the linear part, it belongs to one of three explicit tame
> families.

This is a restricted exclusion, not a proof that every degree-four Keller
map in dimension three is invertible.  It eliminates the smallest complete
degree-four stratum carrying the same mixed-sign torus action as the
foundational degree-seven counterexample.

## 2. Complete normalized support

The three weights are distinct.  A Keller map has an invertible linear part,
so after a target permutation and diagonal rescaling its linear part is the
identity.  Equivariance forbids constant terms.  Enumerating all monomials
of total degree at most four with the required weights gives

\[
\begin{aligned}
F_1={}&x+a x^2y+b x^3z,\\
F_2={}&y+c xz+dxy^2+e x^2yz,\\
F_3={}&z+fy^2+gxyz+hxy^3+i x^2z^2.                 \tag{2.1}
\end{aligned}
\]

Thus the complete normalized problem has nine coefficients.  Expanding
\(\det JF-1\) gives eleven source monomials.  Their coefficients generate an
ideal

\[
I\subset k[a,b,c,d,e,f,g,h,i].
\]

The exact radical calculation is

\[
\sqrt I=P_1\cap P_2\cap P_3,                         \tag{2.2}
\]

where

\[
\begin{aligned}
P_1&=(i,h,g,e,b,a,cf-d),\\
P_2&=(i,g,e,d,c,b,a),\\
P_3&=(h,e,d,b,a,g^2-4fi,cg-2i,2cf-g).
\end{aligned}
\]

Equivalently, a compact Gröbner basis for the radical is

\[
\begin{gathered}
a,\ b,\ e,\ hi,\ di,\ gh,\ dh,\ ch,\ g^2-4fi,\\
dg,\ cg-2i,\ 2cf-2d-g.                              \tag{2.3}
\end{gathered}
\]

Equation (2.2), applied after extending scalars to an algebraic closure, is
a set-theoretic classification of every characteristic-zero Keller point in
the normalized coefficient space; embedded or nilpotent scheme structure
cannot create another map.

## 3. The three components are tame

On \(V(P_1)\), formula (2.1) becomes

\[
F=(x,\ y+cx(z+fy^2),\ z+fy^2).
\]

Writing the target as \((X,Y,Z)\), its inverse is

\[
x=X,\qquad y=Y-cXZ,\qquad z=Z-f(Y-cXZ)^2.            \tag{3.1}
\]

On \(V(P_2)\),

\[
F=(x,\ y,\ z+fy^2+hxy^3),
\]

with the evident triangular inverse

\[
(X,Y,Z)\longmapsto(X,Y,Z-fY^2-hXY^3).               \tag{3.2}
\]

On \(V(P_3)\), the relations give \(g=2cf\) and \(i=c^2f\), so

\[
F=(x,\ y+cxz,\ z+f(y+cxz)^2).
\]

Its inverse is

\[
z=Z-fY^2,\qquad y=Y-cX(Z-fY^2),\qquad x=X.           \tag{3.3}
\]

Each family is a composition of elementary shears.  This proves the theorem.

## 4. Meaning for the minimum-degree frontier

Vistoli's theorem gives the unrestricted lower bound

\[
\delta_3\ge4,
\]

while the foundational map gives \(\delta_3\le7\).  The present theorem does
not raise the lower bound.  It proves that a degree-four counterexample
cannot preserve the foundational \((1,-1,-2)\) grading.  Such a
counterexample must therefore do at least one of the following:

1. have no nontrivial linear torus symmetry;
2. use a different mixed weight signature;
3. acquire equivariance only after nonlinear rather than affine changes;
4. lie in a degree-four leading-form stratum not visible to a torus action.

The next finite extension is to enumerate primitive mixed weight signatures
whose degree-four support contains a coupled nonlinear dependency cycle.
Signatures supporting only acyclic shears are automatically tame and should
be removed combinatorially before coefficient elimination.

## 5. Reproduction

Run

```bash
.venv/bin/python scripts/verify_degree_four_foundational_weight.py
```

The checker independently:

1. enumerates the complete equivariant monomial support;
2. derives all eleven determinant equations;
3. asks Singular to compute \(\sqrt I\) and compares it with
   \(P_1\cap P_2\cap P_3\) in both directions;
4. substitutes all three components into the determinant;
5. composes each displayed inverse on both sides.

## Literature boundary

- Angelo Vistoli, *The Jacobian conjecture in dimension 3 and degree 3*,
  J. Pure Appl. Algebra **142** (1999), 79--89,
  <https://doi.org/10.1016/S0022-4049(98)00040-1>.
- T. Shaska, *Graded Keller maps and the Jacobian Conjecture*,
  arXiv:2607.20210 (2026), <https://arxiv.org/abs/2607.20210>.
  That paper proves the all-one-sign theorem and the dimension-two theorem
  for every signature, and derives the quotient-Jacobian identity for the
  foundational mixed weights.  The degree-four component calculation above
  is the bounded three-dimensional continuation.
