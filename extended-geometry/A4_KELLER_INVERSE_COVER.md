# The oriented quartic \(A_4\) Keller inverse cover

## 1. First proper-monodromy checkpoint

The smallest useful absolute pilot for the
[Keller inverse-Galois program](../KELLER_INVERSE_GALOIS_PROGRAM.md)
is not a regular cyclic cover. If the degree-\(N\) function-field extension
of a characteristic-zero polynomial Keller map is normal, the
Campbell--Razar--Wright Galois-case theorem makes the map invertible. Thus a
nontrivial regular action cannot be the generic inverse action of an
absolute polynomial Keller map.

The natural four-point action

\[
 A_4\leq S_4
\]

avoids this obstruction: it is transitive, but a point stabilizer has order
three and is not normal. This note constructs the exact \(A_4\) inverse cover
as a determinant-one morphism of smooth affine boundary complements. It
therefore reaches the **chart output** of the program, with a controlled full
fiber and all three unramified \(A_4\) cycle types. It does not yet give a
polynomial self-map of affine space.

## 2. The universal oriented quartic

Work over \(\mathbb Q\). Put

\[
 P(T)=T^4+pT^2+qT+r
\]

and

\[
\begin{aligned}
\Delta(p,q,r)={}&256r^3-128p^2r^2+144pq^2r-27q^4\\
                &\quad+16p^4r-4p^3q^2.
\end{aligned}                                       \tag{2.1}
\]

Direct calculation gives

\[
 \operatorname{Disc}_T(P)=\Delta.
\]

Define the oriented squarefree base

\[
 B^\circ=
 \operatorname{Spec}
 \frac{\mathbb Q[p,q,r,D,D^{-1}]}
      {(D^2-\Delta(p,q,r))}.                         \tag{2.2}
\]

It is a smooth affine threefold because
\(\partial(D^2-\Delta)/\partial D=2D\) is a unit. The root incidence

\[
 X^\circ=
 \operatorname{Spec}_{B^\circ}
 \mathcal O_{B^\circ}[T]/(P(T))                     \tag{2.3}
\]

is finite étale of degree four. Indeed \(D\) is a unit, so the discriminant
is a unit and \(P'(T)\) is a unit in the root algebra.

## 3. Exact generic monodromy

Over

\[
 K_0=\mathbb Q(p,q,r),
\]

the generic depressed quartic has geometric group \(S_4\). Indeed its
\(p=0\) slice is the universal inverse pencil

\[
 T^4+qT+r,
\]

whose geometric monodromy is \(S_4\) by the
[universal symmetric-monodromy
theorem](../verified/UNIVERSAL_SYMMETRIC_MONODROMY.md). The ambient generic
group contains the slice group and is therefore also \(S_4\).

As an independent arithmetic certificate, one exact specialization is

\[
 T^4-T-1.
\]

It is irreducible modulo \(2\), giving a 4-cycle, and has factorization type
\((3,1)\) modulo \(7\), giving a 3-cycle. A transitive subgroup of \(S_4\)
containing both cycles is \(S_4\).

Inside the generic splitting field, the unique quadratic sign subfield is

\[
 K_0(\sqrt{\Delta}).
\]

Base change to the oriented equation \(D^2=\Delta\) therefore replaces
\(S_4\) by its index-two subgroup \(A_4\). Since \(A_4\) remains transitive in
its natural four-point action, \(P\) stays irreducible over the oriented
function field. Consequently

\[
 \boxed{
 G_{\mathrm{geom}}=G_{\mathrm{arith}}=A_4
 }
                                                            \tag{3.1}
\]

for \(X^\circ\to B^\circ\). The permutation representation is the natural
degree-four action, not the regular degree-twelve action.

## 4. Derivative-unit Kellerization

The derivative

\[
 j=P'(T)=4T^3+2pT+q
\]

is a unit on \(X^\circ\). Its exact inverse modulo \(P\) has denominator
\(\Delta=D^2\):

\[
 \frac1{P'(T)}
 =
 -\frac{N(T,p,q,r)}{\Delta},                         \tag{4.1}
\]

where

\[
\begin{aligned}
N={}&8T^3p^3-32T^3pr+36T^3q^2-4T^2p^2q-48T^2qr\\
   &+8Tp^4-48Tp^2r+42Tpq^2+64Tr^2\\
   &+4p^3q-48pqr+27q^3.
\end{aligned}
\]

Thus (4.1) is a regular function on the affine chart (2.3).

Choose the hypersurface residue form \(\omega_B\) on \(B^\circ\) and the
root-incidence residue form \(\omega_X\) on \(X^\circ\). With compatible
orientations,

\[
 f^*\omega_B=P'(T)\omega_X.                          \tag{4.2}
\]

Introduce one source coordinate \(z\) and one target coordinate \(Z\). The
one-unit suspension

\[
\begin{aligned}
\widehat f:
 X^\circ\times\mathbb A^1_z
 &\longrightarrow B^\circ\times\mathbb A^1_Z,\\
(p,q,r,D,T,z)
&\longmapsto
\left(p,q,r,D,\frac{z}{P'(T)}\right)
\end{aligned}                                       \tag{4.3}
\]

is a finite étale degree-four morphism and

\[
 \widehat f^*(\omega_B\wedge dZ)
 =\omega_X\wedge dz.                                \tag{4.4}
\]

Hence its residue Jacobian is one.

There is also a direct local-coordinate check. On the root incidence solve

\[
 r=-T^4-pT^2-qT.
\]

In coordinates \((p,q,T,z)\), the nontrivial part of (4.3) is

\[
 (r,Z)=
 \left(-T^4-pT^2-qT,\frac z{4T^3+2pT+q}\right).
\]

Therefore

\[
 \det
 \frac{\partial(p,q,r,Z)}
      {\partial(p,q,T,z)}
 =-1.                                                \tag{4.5}
\]

The sign depends only on the chosen ordering of the residue forms.

Adjoining \(z\) does not change the inverse cover: for a fixed target value
\(Z=Z_0\), every root \(T\) has the unique lift

\[
 z=P'(T)Z_0.
\]

Thus the generic monodromy of (4.3) is still exactly \(A_4\), rather than an
ambient \(S_4\) which merely contains special \(A_4\)-fibers.

## 5. One controlled full \(A_4\) fiber

Take

\[
 P_0(T)=T^4-7T^2-3T+1.                              \tag{5.1}
\]

Its discriminant is

\[
 \operatorname{Disc}(P_0)=33489=183^2.              \tag{5.2}
\]

The quartic is irreducible, and its cubic resolvent

\[
 R_0(Y)=Y^3+7Y^2-4Y-37
\]

is irreducible. A transitive subgroup of \(A_4\) is \(V_4\) or \(A_4\);
irreducibility of the resolvent excludes \(V_4\). Hence

\[
 \operatorname{Gal}(P_0/\mathbb Q)=A_4.             \tag{5.3}
\]

It has four real roots, so this particular connected fiber is totally real.

At the rational oriented target

\[
 (p,q,r,D)=(-7,-3,1,183)
\]

and any \(Z_0\in\mathbb Q\), the complete fiber of (4.3) is

\[
 \operatorname{Spec}\mathbb Q[T]/(P_0),
\qquad z=P_0'(T)Z_0.                                \tag{5.4}
\]

It has length four, equal to the generic degree, so no sheet is lost at the
boundary.

The good-prime factorization types

\[
\begin{array}{c|c|c}
\text{prime}&\text{factor degrees of }P_0&\text{\(A_4\) class}\\
\hline
2&(3,1)&\text{3-cycle}\\
11&(2,2)&\text{double transposition}\\
233&(1,1,1,1)&\text{identity}
\end{array}                                         \tag{5.5}
\]

realize every cycle partition occurring in the natural \(A_4\)-action.

## 6. What remains for an absolute map

Construction (4.3) is already a prescribed-monodromy Keller inverse cover,
but on smooth affine boundary complements:

- the target remembers the orientation \(D\) and inverts it;
- the source is the root-incidence cover;
- the coordinate \(z/P'(T)\) uses the unit \(D^{-2}\).

The direct one-coordinate pole clearing does not solve this. Put
\(z=D^2w\). Equation (4.1) makes the new output

\[
 Z=-wN(T,p,q,r),
\]

which is polynomial, but the local Jacobian becomes

\[
 \det
 \frac{\partial(p,q,r,Z)}
      {\partial(p,q,T,w)}
 =-D^2.                                               \tag{6.1}
\]

Thus the denominator disappears only by reintroducing the discriminant as a
critical factor. Polynomial automorphisms or identity stabilization cannot
turn this factor into a constant; a successful completion has to couple the
orientation modification to at least one additional boundary coordinate.

Forgetting \(D\) takes the quotient by the sign involution and changes the
generic group from \(A_4\) back to \(S_4\). This is exactly the
section-transfer trap in geometric form: the simplest affine quotient loses
the datum that enforces proper monodromy.

The next absolute problem is therefore precise:

> Complete the oriented target and root-incidence source to affine spaces,
> absorbing the \(D^{-2}\) derivative unit into polynomial coordinates,
> without quotienting the orientation or changing the four-sheet function
> field.

The [oriented cubic Cox
charts](ORIENTED_CUBIC_COX_CHART.md) identify the likely source obstruction.
After one root is selected, the remaining three roots carry a cyclic
orientation; the natural oriented cubic chart is not affine space, and its
known affine quotient forgets the orientation. More strongly, the
[linear-hyperplane
classification](LINEAR_HYPERPLANE_COX_CLASSIFICATION.md) proves that every
linear coefficient slice of this normalized three-factor source is
non-affine even after arbitrary polynomial stabilization. An \(A_4\)
absolute construction must therefore use a genuinely
orientation-preserving nonlinear affine modification or change the boundary
model itself.

The follow-up
[affine \(A_4\) frontier](A4_AFFINE_KELLER_FRONTIER.md) removes the
orientation denominator completely. It constructs a polynomial
constant-residue-Jacobian \(A_4\) map on the oriented cubic quotient and a
literal polynomial map \(\mathbb A^3\to\mathbb A^3\) with generic \(A_4\)
monodromy and determinant \(4W^2K^3L\). Thus prescribed monodromy and
affine-space algebraization are now simultaneous; cancelling this explicit
divisor without changing the cover is the remaining Keller gate.

## 7. Reproduction

Run

```bash
.venv/bin/python scripts/verify_a4_keller_inverse_cover.py
```

The checker verifies the quartic discriminant, the exact quotient-ring
inverse of \(P'\), the constant suspension Jacobian, the \(D^2\) obstruction
to direct pole clearing, an \(S_4\) certificate before orientation, the
explicit totally real \(A_4\) fiber, and its three good-prime cycle types.
