# Fermat-symmetry normal forms for the smooth-quartic line fibers

## Status

This note proves `HC4NHM23`.  It continues `HC4NHM17` and `HC4NHM20`
invariantly: the automorphism group of the fixed binary Fermat cubic reduces
the fifteen degenerate polar fibers to three orbit types.  The rational
fiber `tau=-1` and the two roots of `tau^2-tau+1` form one orbit, so the
exact generic-component and first-secondary-stratum exclusions already
proved at `tau=-1` transport to all three fibers.

Replay the covariance and orbit identities with

```bash
.venv/bin/python scripts/verify_hc4_smooth_quartic_fermat_symmetry_orbits.py
```

This is a transport theorem, not a new sampled calculation.  The standard
bases being transported are the exact characteristic-zero certificates of
`HC4NHM17`.

## 1. The normalized reciprocal packet is covariant

Retain

\[
s_3=\frac{x^3+y^3}{3}+z^2(ux+vy)+wz^3,
\qquad \ell=y+\tau x+\sigma z,
\tag{1.1}
\]

and the boundary matrix

\[
A_0=\begin{pmatrix}
0&0&-y^2\\
0&0&x^2\\
-y^2&x^2&px^2+qxy+ry^2
\end{pmatrix}.
\tag{1.2}
\]

For an invertible coordinate matrix (T), put

\[
s_3'(X)=s_3(TX),\qquad
A'(X)=\det(T)T^{-1}A(TX)T^{-\mathsf T}.
\tag{1.3}
\]

If (d=\nabla s_3), then

\[
d'=T^{\mathsf T}d(TX),\qquad
\operatorname{adj}(A')+d'd'^{\mathsf T}
=T^{\mathsf T}\bigl(\operatorname{adj}(A)+dd^{\mathsf T}\bigr)(TX)T.
\tag{1.4}
\]

Consequently Hessian curl of
(C=(\operatorname{adj}(A)+dd^{\mathsf T})/z) is preserved.  If the old
(z)-coordinate is (c z), set

\[
\ell'(X)=\ell(TX)/c.
\tag{1.5}
\]

Then divisibility of (R=\det(A)/z) by (ell), and the scalar equation

\[
\ell(R+e^{\mathsf T}d)=zR,\qquad e=Ad/z,
\tag{1.6}
\]

are also preserved.  Thus (1.3) transports the complete 81-equation packet,
not merely its visible pivot.

Two transformations preserve the normal form (1.1)--(1.2).

1. For (lambda^2+lambda+1=0), take
   (T=\operatorname{diag}(\lambda,1,1)).  Then

   \[
   \tau'=\lambda\tau,quad
   (p',q',r')=(p,\lambda^2q,\lambda r),quad
   (u',v',w')=(\lambda u,v,w),quad \sigma'=\sigma.
   \tag{1.7}
   \]

2. To normalize the reflection, take

   \[
   T=\begin{pmatrix}0&1&0\\1&0&0\\0&0&\tau\end{pmatrix}.
   \tag{1.8}
   \]

   It sends

   \[
   \tau'=\tau^{-1},\quad
   (p',q',r')=-\tau^{-1}(r,q,p),quad
   (u',v',w')=(\tau^2v,\tau^2u,\tau^3w),quad \sigma'=\sigma.
   \tag{1.9}
   \]

The general (zB) term remains a general (zB') term under both invertible
linear transformations.  Hence no chart coefficient is lost.

## 2. Equivariance of the polar divisor

For the pivot (Delta(\tau;p,q,r)) of `HC4NHM16`, direct substitution gives

\[
\Delta(\lambda\tau;p,\lambda^2q,\lambda r)
=\lambda^2\Delta(\tau;p,q,r),
\tag{2.1}
\]

and

\[
\Delta\left(\tau^{-1};-\frac r\tau,-\frac q\tau,-\frac p\tau\right)
=\tau^{-7}\Delta(\tau;p,q,r).
\tag{2.2}
\]

Thus the two generators act on the degenerate fibers exactly as the
dihedral presentation of the (S_3)-automorphism group of
(s^3+t^3).

## 3. Three, not fifteen, slope normal forms

The degree-fifteen degeneration resultant from `HC4NHM20` reorganizes as

\[
\begin{aligned}
D(\tau)
={}&(\tau^3+1)
\cdot(\tau^{12}+44\tau^9+586\tau^6+44\tau^3+1).
\end{aligned}
\tag{3.1}
\]

It is invariant under (	au\mapsto\lambda\tau), and reflection gives
(D(\tau^{-1})=\tau^{-15}D(\tau)).  Put (s=\tau^3).  Away from the first
orbit (s=-1), the reciprocal quartic is

\[
s^4+44s^3+586s^2+44s+1
=s^2\bigl(j^2+44j+584\bigr),
\qquad j=s+s^{-1}.
\tag{3.2}
\]

The two remaining quotient values are

\[
j=-22+10i,qquad j=-22-10i.
\tag{3.3}
\]

Therefore the fifteen slopes form exactly three (S_3)-orbits:

\[
\boxed{3+6+6}.
\tag{3.4}
\]

The size-three orbit is (	au^3=-1).  The other twelve slopes require only
two algebraic normal-form calculations, one for each value in (3.3).  The
degree-four and degree-eight rational factors in `HC4NHM20` are arithmetic
factorizations; they are not the geometric orbit decomposition.

## 4. Transport of the first orbit

At (	au=-1), `HC4NHM17` treats the resultant line

\[
3p+q+3r=0
\tag{4.1}
\]

and the residual polar line

\[
p-r=0.
\tag{4.2}
\]

Apply (1.7).  The two images have
(	au^2-	au+1=0), and their equations can be written invariantly over
(mathbf Q(\tau)) as

\[
-3\tau p+(\tau-1)q+3r=0
\tag{4.3}
\]

and

\[
\tau p+r=0,
\tag{4.4}
\]

respectively.  These are precisely the resultant and residual-polar factors
of the split conic.

The first secondary strata from `HC4NHM17` transport as well.  On (4.4),

\[
\tau^2q^2-3\tau pq+8p^2=0,
\tag{4.5}
\]

is the image of (q^2+3pq+8p^2=0).  On (4.3),

\[
7p-33\tau^2r=0
\tag{4.6}
\]

is the image of (7p-33r=0).  Since (1.3)--(1.6) preserve the full
reciprocal-Hessian ideal and determinant-zero support, the exact bases at
(	au=-1) prove the corresponding statements at both quadratic roots.

## 5. Result and remaining line fibers

> **Theorem `HC4NHM23` -- Fermat-symmetry orbit reduction.**  The fifteen
> two-line fibers of the smooth-quartic polar divisor form exactly three
> automorphism orbits of sizes (3,6,6).  The complete normalized
> reciprocal-Hessian packet is covariant under the orbit generators.  Hence
> the generic loci of both line components, and their first registered
> secondary strata, are determinant-zero at every slope in
> (	au^3=-1): the rational slope treated by `HC4NHM17` and both roots of
> (	au^2-	au+1).

This replaces four additional exceptional-line specializations by transport
of one exact certificate.  Among non-generic polar fibers, only the two
six-point normal-form types (3.3) remain genuinely new.  Further lower
secondary strata inside the first orbit, and the parameterization-denominator
strata of the generic polar conic, remain open and are not claimed here.
