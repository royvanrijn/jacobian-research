# The \(A_4\) determinant-ledger reduction

## 1. Outcome

For the polynomial \(A_4\) cone

\[
 \Phi(U,V,W)=(WN_1,WN_2,WH)
\]

from the [affine \(A_4\) frontier](A4_AFFINE_KELLER_FRONTIER.md), the
Jacobian divisor

\[
 \det D\Phi=4W^2K^3L
\]

admits an exact target-ledger reduction.  If

\[
\begin{aligned}
\mathcal B(P,Q,R)
={}&P^3-3PQ^2+2Q^3-9PQR+9Q^2R\\
   &-27PR^2+27QR^2+27R^3,
\end{aligned}
\]

then

\[
 \boxed{\mathcal B(\Phi)=W^3K^3L^2}                 \tag{1.1}
\]

and consequently

\[
 \boxed{\frac{\mathcal B(\Phi)}{\det D\Phi}
       =\frac{WL}{4}.}                              \tag{1.2}
\]

Thus the target divisor \(\mathcal B=0\) accounts for the whole \(K^3\)
column and one copy of \(L\).  The residual source ledger is \(WL\), not
the original four-factor expression.

Two natural completion routes can now be closed exactly.

1. Every correction of the oriented-surface map by a multiple of its
   defining defect still has zero Jacobian at one fixed rational point.
   This excludes corrections of every polynomial degree, not merely
   affine-linear corrections.
2. Polynomial reparametrization inside the present common-denominator
   target model cannot remove the cone factor \(W^2\).  A rational target
   modification or a genuinely coupled stabilization is necessary.

This does not construct an absolute Keller map.  It reduces the surviving
problem to a coupled affine modification for the two-factor residual
boundary \(WL\).

## 2. The target \(B\)-divisor

Retain the polynomials

\[
\begin{aligned}
H={}&8U^3-6UV^2-18UV-54U
      -2V^3-9V^2-27V-27,\\
K={}&4U^2+4UV+6U+V^2+3V+9,\\
L={}&U^3-3UV^2-9UV-27U
      +2V^3+9V^2+27V+27,\\
M={}&U^2+2V^2+6V+18,
\end{aligned}
\]

with \(N_1=MK\) and

\[
\begin{aligned}
N_2={}&8U^3V+12U^2V^2+36U^2V+108U^2\\
     &+6UV^3+36UV^2+108UV+162U\\
     &+V^4+9V^3+27V^2+54V.
\end{aligned}
\]

For

\[
 \alpha=\frac{N_1}{H},\qquad \beta=\frac{N_2}{H},
\]

the polynomial called \(B\) in the Jensen--Ledet--Yui generic \(A_4\)
polynomial satisfies the stronger pullback identity

\[
\boxed{
B(\alpha,\beta)
=\alpha^3-3\alpha\beta^2+2\beta^3-9\alpha\beta
9\beta^2-27\alpha+27\beta+27
=\frac{K^3L^2}{H^3}.
}                                                     \tag{2.1}
\]

Equation (1.1) is its homogeneous cone form.  On the rational plane cover,

\[
 \det\frac{\partial(\alpha,\beta)}{\partial(U,V)}
 =\frac{4K^3L}{H^3},
\]

so already there

\[
 \frac{B(\alpha,\beta)}
 {\det\partial(\alpha,\beta)/\partial(U,V)}
 =\frac L4.                                          \tag{2.2}
\]

The extra \(W\) in (1.2) is exactly the common-denominator cone boundary.

## 3. Geometry of the residual divisor

The two plane factors have elementary normalizations.

First,

\[
 K=(2U+V)^2+3(2U+V)+9.                               \tag{3.1}
\]

Thus \(K=0\) is irreducible over \(\mathbb Q\) and becomes two parallel
affine lines over \(\mathbb Q(\sqrt{-3})\).  Identity (1.2) shows that this
entire geometrically reducible multiplicity-three column belongs to the
target ledger.

Second, put \(t=U-V\).  The normalization of \(L=0\) is

\[
\begin{aligned}
V&=\frac{-t^3+27t-27}{3t(t-3)},\\
U&=\frac{(2t-3)(t^2-3t+9)}{3t(t-3)}.
\end{aligned}                                       \tag{3.2}
\]

The inverse is \(t=U-V\).  The values \(t=0,3,\infty\) are precisely the
three points above the projective boundary, hence

\[
 \widetilde{V(L)}\simeq
 \mathbb P^1\setminus\{0,3,\infty\}.                 \tag{3.3}
\]

The residual factor \(L\) therefore has three places at infinity.  It is
not one of the \(\mathbb A^1\)- or \(\mathbb G_m\)-normalized one-boundary
cores handled by the established weighted, cancellation, or quadratic-gauge
suspensions.

## 4. All-degree obstruction to the defect lift

Let

\[
\begin{aligned}
\delta(x,y)&=x^2y^2-4x^3-4y^3+18xy-27,\\
R&=d^2-\delta(x,y),\\
F_0(x,y,d)&=(x^2-2y,\ y^2-2x,\ d(xy-1)).
\end{aligned}
\]

Consider arbitrary polynomials \(f,g,h\) and the most general correction
which preserves the oriented-surface map modulo its defining equation:

\[
 \widetilde F=F_0+R(f,g,h).                           \tag{4.1}
\]

At the rational point

\[
 p=(-1,-1,0)
\]

one has \(R(p)=0\) and

\[
 DF_0(p)=
\begin{pmatrix}
-2&-2&0\\
-2&-2&0\\
0&0&0
\end{pmatrix},
\qquad
dR(p)=(32,32,0).                                    \tag{4.2}
\]

The derivative of the correction at \(p\) is

\[
 D(R(f,g,h))(p)=
 \begin{pmatrix}f(p)\\g(p)\\h(p)\end{pmatrix}dR(p),
                                                               \tag{4.3}
\]

which has rank at most one.  Its row space is contained in the same
one-dimensional row direction as \(DF_0(p)\).  Hence

\[
 \operatorname{rank}D\widetilde F(p)\le1
\quad\text{and}\quad
 \boxed{\det D\widetilde F(p)=0}.                    \tag{4.4}
\]

Therefore no correction (4.1), of any degree, is Keller.  Nonlinear
searches inside this defect-multiple ansatz cannot succeed.

## 5. Rigidity of polynomial cone reparametrization

There is also an intrinsic obstruction inside the present affine target
ring.  The polynomial \(H\) is irreducible over \(\mathbb Q\).  Along its
generic point, both \(\alpha\) and \(\beta\) have a simple pole and the
ratio

\[
 [N_1:N_2]\big|_{H=0}
\]

is nonconstant.  One exact certificate is that the tangent derivative of
\(N_1/N_2\) along \(H=0\) is nonzero modulo \(H\).

It follows that

\[
 q(\alpha,\beta)\in\mathbb Q[U,V],
 \quad q\in\mathbb Q[\alpha,\beta]
 \quad\Longrightarrow\quad q\in\mathbb Q.            \tag{5.1}
\]

Indeed, if \(q\) has positive degree \(m\), its top homogeneous part gives
a pole of order \(m\) along \(H=0\).  Cancellation would force a nonzero
homogeneous binary form to vanish on the nonconstant ratio
\([N_1:N_2]\), which is impossible.

Put \(T=WH\).  Let

\[
 q(\alpha,\beta,T)\in\mathbb Q[\alpha,\beta,T]
\]

have polynomial pullback to \(\mathbb Q[U,V,W]\).  Its restriction to
\(W=0\) is constant by (5.1): every positive power of \(T\) vanishes there,
and the \(T^0\)-coefficient can be polynomial only when it is constant.
Consequently, any three polynomial coordinates obtained from this affine
target ring and pulled back to the cone are constant on \(W=0\).

After subtracting those constants, all three coordinates are divisible by
\(W\).  Their \(U\)- and \(V\)-derivative columns are therefore divisible
by \(W\), proving

\[
 \boxed{W^2\mid\det D(q_1,q_2,q_3)}.                 \tag{5.2}
\]

This is a no-go theorem for polynomial target reparametrizations of the
current common-denominator model.  It does not exclude rational target
affine modifications: those are exactly where new boundary valuations can
enter.

## 6. Stabilization must feed back into the cover coordinates

There is a second elementary restriction.  If a polynomial stabilization
retains \(\Phi\) as its first three outputs,

\[
 G(U,V,W,z)=\bigl(\Phi(U,V,W),G_{\mathrm{aux}}(U,V,W,z)\bigr),
\]

then its derivative is block triangular and

\[
 \det DG=(\det D\Phi)
          \det\frac{\partial G_{\mathrm{aux}}}{\partial z}.     \tag{6.1}
\]

Thus \(4W^2K^3L\) still divides \(\det DG\).  Adding primitive variables
only below the existing cover block cannot cancel the divisor.  Any
successful stabilization must let the new variables enter at least one of
the three cover outputs, while still providing a target-field certificate
that the generic extension remains the natural \(A_4\)-extension.

## 7. Surviving construction problem

The [pure-target ledger lift](A4_PURE_TARGET_LEDGER_LIFT.md) completes the
first step below polynomially: after adjoining \(z\), the fourth output
\((WL/4)z\) makes the full determinant exactly
\(\mathcal B(\Phi)=W^3K^3L^2\).  Thus the remaining defect is a pure target
pullback and the generic \(A_4\) extension is unchanged.

The remaining attack is now specific:

1. use the target divisor \(\mathcal B=0\) to absorb the \(K^3L\) portion
   certified by (1.2);
2. construct a coupled source/target affine modification for \(WL\);
3. allow the modification variables to feed back into the cover outputs,
   as required by (6.1);
4. prove that eliminating the auxiliary variables recovers the original
   degree-four \(A_4\) function-field extension;
5. pass a complete polar-ledger test across the three punctures of the
   normalization (3.3).

An absolute solution of these five steps would be a degree-four
counterexample to the Jacobian conjecture.  The present reduction should
therefore be read as a sharp construction frontier, not as evidence that a
routine suspension remains.

## 8. Reproduction

Run

```bash
.venv/bin/python scripts/verify_a4_ledger_reduction.py
```

The checker verifies (1.1), (2.1), the residual ratio (1.2), the
normalization of \(L\), the nonconstant \(H\)-boundary ratio, and the
all-degree rank obstruction at \((-1,-1,0)\).
