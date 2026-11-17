# Davenport normalized-boundary involution

The Davenport derivative boundary has normalized affine curve
\(\mathbb P^1\) minus three points, while the split determinant
modification is centered on \(\mathbb G_m\).  There is a canonical
boundary-level quotient which merges exactly two punctures and produces
\(\mathbb G_m\).

However, this involution does not preserve the two branches above the node.
It therefore fails to descend from the normalization to the singular
derivative curve and cannot yet be lifted to the ambient Davenport cover.

Work over

\[
K=\mathbb Q(a),\qquad a^2+a+2=0.
\]

## 1. Explicit rational normalization

Let \(J(T,Y)=g_T'(Y)\).  On the torus put

\[
z=\frac{T}{Y^2}.
\]

Then \(J/Y^4=0\) is quadratic in \(Y\).  Its discriminant is

\[
-4z^2(z+1)^2\bigl((4a+16)z+2-3a\bigr).
\]

Introduce \(v\) by

\[
v^2=-\bigl((4a+16)z+2-3a\bigr),
\qquad
z=\frac{3a-2-v^2}{4a+16}.                            \tag{1}
\]

The quadratic formula expresses \(Y\), and hence \(T=zY^2\), rationally in
\(v\).  After cancellation over \(K\), the pole divisor is

\[
\boxed{
(v+a)
\left(v^2+(-2+3a)v-8-2a\right).
}                                                     \tag{2}
\]

The linear factor gives one \(K\)-rational puncture.  The quadratic factor
is irreducible over \(K\) and gives one conjugate pair.

The ordinary node \((T,Y)=(-1/4,-1/2)\) has \(z=-1\).  Its two points in
the normalization are

\[
v_+=a+4,\qquad v_-=-a-4.                             \tag{3}
\]

## 2. The unique puncture-swapping involution

Move the rational puncture \(v=-a\) to infinity:

\[
\xi=\frac1{v+a}.
\]

In the \(\xi\)-coordinate, the quadratic puncture pair has trace
\(-1/2\).  Therefore

\[
\iota_\xi(\xi)=-\frac12-\xi                          \tag{4}
\]

is the unique nontrivial \(K\)-automorphism preserving the three-puncture
set.  It fixes the rational puncture and swaps the conjugate pair.

The invariant

\[
w=\left(\xi+\frac14\right)^2                         \tag{5}
\]

sends the quadratic pair to one finite \(K\)-point and the rational
puncture to infinity.  Consequently

\[
\left(
\mathbb P^1\setminus\{P_0,P_1,P_2\}
\right)/\langle\iota\rangle
\simeq\mathbb P^1\setminus\{\infty,w_0\}
\simeq\mathbb G_m.                                  \tag{6}
\]

Thus the unit-rank mismatch from the
[derivative-center audit](DAVENPORT_DERIVATIVE_CENTER_MISMATCH.md) can be
repaired on the normalized boundary itself.

## 3. Conductor obstruction

For an automorphism of the normalization to descend to the nodal curve, it
must preserve the unordered pair \(\{v_+,v_-\}\) above the node.
Exact substitution in (4) gives

\[
\iota(v_-) = v_-,
\qquad
\iota(v_+)=-\frac{5a+6}{4}.
\]

The second value is neither \(v_+\) nor \(v_-\).  Hence

\[
\boxed{
\iota\{v_+,v_-\}\ne\{v_+,v_-\}.
}                                                     \tag{7}
\]

The involution does not preserve the conductor and does not descend to
\(J=0\).

Because the rational puncture is the unique degree-one orbit in (2), every
\(K\)-automorphism preserving the puncture set fixes it.  The only
possibilities are the identity and the transposition (4).  Therefore the
singular affine derivative curve has no nontrivial \(K\)-automorphism
induced by a puncture-preserving automorphism of its normalization.

## 4. Minimal equivariant conductor closure

The conductor cannot be repaired by adding only one new identification.
The original node relation is

\[
v_+\sim v_-.
\]

Applying \(\iota\) gives

\[
\iota(v_+)\sim\iota(v_-)=v_-.
\]

The three points

\[
v_+,\qquad v_-,\qquad \iota(v_+)
\]

are distinct.  Therefore the smallest \(\iota\)-invariant equivalence
relation containing the original conductor glues all three points.  Its
pushout has a three-branch seminormal point rather than a node.

In the quotient normalization, \(v_+\) and \(\iota(v_+)\) form one orbit,
while \(v_-\) is fixed.  The descended conductor identifies these two
quotient points, so the quotient curve is still nodal.  It is not the
smooth \(\mathbb G_m\) obtained from the fully normalized boundary.

Thus an equivariant finite re-gluing does not solve the center problem: the
minimal repair introduces a third branch, echoing the extra-boundary
phenomenon in controlled suspensions.

## 5. Consequence

The desired \(\mathbb G_m\) center is now available as a canonical quotient
of the normalized boundary, but not as a quotient of the actual nodal
boundary.  A surviving descent must first modify the conductor—for example
by fully separating the node before quotienting without a finite
re-gluing—or realize a different involution in a larger oriented
neighborhood.

This is narrower than the original affine-space problem: the puncture
geometry is solved, and the obstruction is the explicit two-point
conductor (3).

## 6. Verification

Run

```bash
.venv/bin/python scripts/verify_davenport_boundary_involution.py
```

The checker constructs the rational normalization, factors its puncture
divisor, verifies the involution and its \(\mathbb G_m\) quotient, locates
the two node branches, proves the conductor failure (7), and verifies that
the minimal invariant conductor closure has three branches and nodal
quotient.
