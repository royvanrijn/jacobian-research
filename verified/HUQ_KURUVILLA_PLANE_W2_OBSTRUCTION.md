# The plane characteristic-two map has no Keller lift modulo four

This is a repository extension of Mondello's external plane theorem,
[*A Dimension-Two Counterexample to the Separable Jacobian Conjecture in
Characteristic Two*, arXiv:2608.02634v1](https://arxiv.org/abs/2608.02634).
The lifting obstruction, stabilization, and Witt-tower statements below are
not claims from Mondello's paper.

## 1. Statement

Let

\[
\begin{aligned}
 P&=x+x^2y+x^4+x^6y^2,\\
 Q&=y+x^5+x^6y+x^7y^2+x^8y^3
\end{aligned}
\]

over \(\mathbb F_2\).  This is the determinant-one, noninjective plane map
from the
[characteristic-two audit](HUQ_KURUVILLA_CHARACTERISTIC_TWO_AUDIT.md).

**Theorem (HKM2W1).**

1. For every determinant-one plane map over \(\mathbb F_2\), the cokernel of
   its first Jacobian-variation operator is canonically the top algebraic de
   Rham cohomology of the source:
   \[
   \operatorname{coker}\mathcal D_F
   \simeq H^2_{\mathrm{dR}}(\mathbb F_2[x,y])
   \simeq xy\,\mathbb F_2[x^2,y^2].
   \]
2. For the displayed map, the integral Jacobian error has class
   \[
   [K]=xy\bigl(1+x^6+x^8y^2\bigr)=xyu^2\ne0,
   \qquad u=1+x^3(1+xy).
   \]
   Hence no polynomials
   \(\widetilde P,\widetilde Q\in(\mathbb Z/4)[x,y]\) reducing to
   \(P,Q\) have constant Jacobian determinant.
3. Nonvanishing of this class is invariant under arbitrary polynomial
   left--right equivalence of plane maps over \(\mathbb F_2\).  Thus no
   polynomially left--right equivalent representative has a Keller lift
   through \(W_2(\mathbb F_2)\).
4. Every first-order Jacobian obstruction becomes exact after adjoining one
   identity coordinate.  For this map the exact stabilization extends to an
   explicit compatible Keller lift through every finite Witt ring
   \(W_n(\mathbb F_2)=\mathbb Z/2^n\).

Consequently the plane map has no compatible lift through higher Witt
vectors, through \(\mathbb Z_2\), or to characteristic zero.  These are
unrestricted polynomial statements: no correction-degree bound is imposed.

## 2. First-order Jacobian correction

Use the displayed zero-one formulas as integral representatives \(P_0,Q_0\).
Every lift modulo four is uniquely of the form

\[
 \widetilde P=P_0+2A,\qquad \widetilde Q=Q_0+2B
\]

for polynomials \(A,B\) modulo two.  Direct differentiation gives

\[
 \det D(P_0,Q_0)=1+2K,
\]

where, modulo two,

\[
\begin{aligned}
K={}&xy+x^7y+x^{10}y+x^5y^2+x^{11}y^2\\
   &+x^9y^3+x^{12}y^3+x^{13}y^4.                 \tag{2.1}
\end{aligned}
\]

In particular

\[
 [xy]K=1.                                         \tag{2.2}
\]

Modulo four, the determinant of an arbitrary lift is

\[
 \det D(\widetilde P,\widetilde Q)
 =1+2\bigl(K+\mathcal D_F(A,B)\bigr),            \tag{2.3}
\]

with the first Jacobian-variation operator

\[
 \boxed{\mathcal D_F(A,B)
 =A_xQ_y+P_xB_y+A_yQ_x+P_yB_x}.                   \tag{2.4}
\]

over \(\mathbb F_2\).  The signs become plus signs in characteristic two.

## 3. The full cokernel is top de Rham cohomology

The following argument applies to every determinant-one polynomial plane map
\(G=(G_1,G_2)\) over \(\mathbb F_2\).  Put
\(R=\mathbb F_2[x,y]\) and \(\omega=dx\wedge dy\).  The first variation
satisfies

\[
 \mathcal D_G(A,B)\,\omega
 =dA\wedge dG_2+dB\wedge dG_1
 =d\bigl(A\,dG_2+B\,dG_1\bigr).                \tag{3.1}
\]

Because \(dG_1\wedge dG_2=\omega\), the forms \(dG_1,dG_2\) are an
\(R\)-basis of \(\Omega_R^1\).  Hence every one-form is uniquely of the
form \(A\,dG_2+B\,dG_1\), and (3.1) identifies the image of
\(\mathcal D_G\) with the coefficients of all exact two-forms.  Therefore

\[
 \boxed{\operatorname{coker}\mathcal D_G
 \simeq H^2_{\mathrm{dR}}(R/\mathbb F_2).}       \tag{3.2}
\]

This quotient is elementary.  A monomial \(x^iy^j\) is an \(x\)-derivative
when \(i\) is even and a \(y\)-derivative when \(j\) is even.  Conversely,
a nonzero \(x\)-derivative has even \(x\)-exponent, and a nonzero
\(y\)-derivative has even \(y\)-exponent.  Thus every class has a unique
representative supported on monomials odd in both variables:

\[
 \boxed{H^2_{\mathrm{dR}}(R/\mathbb F_2)
 =xy\,\mathbb F_2[x^2,y^2]\,\omega.}           \tag{3.3}
\]

In particular the cokernel is rank one over the Frobenius subring
\(\mathbb F_2[x^2,y^2]\), but infinite-dimensional over \(\mathbb F_2\).
The old \(xy\)-coefficient test detects only the constant coefficient of
this Frobenius-linear obstruction.

## 4. The exact Cartier class of the Jacobian error

Projecting (2.1) to the odd--odd monomials from (3.3) gives

\[
\begin{aligned}
 [K]
 &=xy+x^7y+x^9y^3\\
 &=xy\bigl(1+x^6+x^8y^2\bigr).                  \tag{4.1}
\end{aligned}
\]

With \(r=1+xy\) and \(u=1+x^3r\),

\[
 u^2=1+x^6r^2=1+x^6+x^8y^2,
\]

so

\[
 \boxed{[K]=xyu^2\ne0.}                          \tag{4.2}
\]

Equivalently, the Cartier isomorphism sends this class to

\[
 C\bigl([K\omega]\bigr)=u\,dx\wedge dy.        \tag{4.3}
\]

This computes the entire obstruction, not only one coefficient.  In
particular \([xy]K=1\) is the lowest-coordinate witness for the nonzero
Cartier class.

For completeness, choose any integral representatives of a determinant-one
map over \(\mathbb F_2\) and write their Jacobian as \(1+2K_G\) modulo
four.  Changing the representatives by \(2(A,B)\) changes \(K_G\) by
\(\mathcal D_G(A,B)\).  Hence

\[
 \mathfrak o(G):=[K_G]\in H^2_{\mathrm{dR}}(R/\mathbb F_2) \tag{4.4}
\]

is independent of every choice.  A constant polynomial represents zero in
top de Rham cohomology, so \(G\) has a constant-Jacobian lift through
\(\mathbb Z/4\) if and only if \(\mathfrak o(G)=0\).  Equation (4.2)
therefore proves the plane lifting obstruction.

## 5. Polynomial left--right invariance

Let \(\sigma,\tau\) be polynomial automorphisms of the source and target
plane over \(\mathbb F_2\), and set

\[
 G=\tau\circ F\circ\sigma.
\]

By the Jung--van der Kulk theorem, every plane polynomial automorphism is a
composition of affine and triangular automorphisms.  Lifting the coefficients
of those generators gives polynomial automorphisms over \(\mathbb Z/4\)
with constant unit Jacobian.  We may therefore compute the obstruction of
\(G\) using compositional lifts.  The chain rule gives, up to a constant
unit whose first-order error is de Rham exact,

\[
 \mathfrak o(G)=\sigma^*\mathfrak o(F).           \tag{5.1}
\]

Pullback by \(\sigma\) is an automorphism of algebraic de Rham cohomology.
Consequently

\[
 \boxed{\mathfrak o(G)=0\iff\mathfrak o(F)=0.}   \tag{5.2}
\]

Since \(\mathfrak o(F)=xyu^2\ne0\), no polynomially left--right equivalent
plane representative over \(\mathbb F_2\) has a constant-Jacobian lift
through \(\mathbb Z/4\).  The same proof works after a perfect extension of
\(\mathbb F_2\), using the corresponding Witt-vector lifts of affine and
triangular generators.

## 6. General stabilization theorem

There is an all-dimensional form of (3.1).  Let
\(R_n=k[x_1,\ldots,x_n]\), let \(G=(G_1,\ldots,G_n)\) have determinant
one, and let \(\omega_n=dx_1\wedge\cdots\wedge dx_n\).  For a correction
\(A=(A_1,\ldots,A_n)\), variation of
\(dG_1\wedge\cdots\wedge dG_n\) gives

\[
 \mathcal D_G(A)\omega_n
 =d\left(
 \sum_{i=1}^n(-1)^{i-1}A_i\,
 dG_1\wedge\cdots\widehat{dG_i}\cdots\wedge dG_n
 \right).                                         \tag{6.1}
\]

Since the \(dG_i\) form a basis, this identifies

\[
 \operatorname{coker}\mathcal D_G
 \simeq H^n_{\mathrm{dR}}(R_n/k).                \tag{6.2}
\]

Now adjoin an identity coordinate \(z\).  If \(K\omega_n\) represents any
first obstruction for \(G\), then

\[
 K\omega_n\wedge dz
 =d\bigl((-1)^n zK\omega_n\bigr).               \tag{6.3}
\]

Thus every first Jacobian obstruction becomes exact after one identity
stabilization.  Explicitly, correcting the new coordinate by \(zK\) kills
the first-order error.  This is a general cohomological explanation of the
instability seen here; it does not by itself construct compatible lifts at
all higher Witt levels.

## 7. One stabilization gives a full finite-Witt tower

The same calculation gives a sharp contrast.  Since

\[
 \det D(P_0,Q_0)=1+2K,
\]

put \(h=2K\) and, for every \(n\ge2\), define

\[
 S_n=\sum_{j=0}^{n-1}(-h)^j,
 \qquad
 \widetilde F_n(x,y,z)=\bigl(P_0,Q_0,zS_n\bigr) \pmod {2^n}.  \tag{7.1}
\]

This reduces to \((P,Q,z)\).  Its Jacobian matrix is block lower triangular.
The finite geometric-series identity gives

\[
 \det D\widetilde F_n
  =(1+h)S_n=1-(-h)^n=1\pmod {2^n}.               \tag{7.2}
\]

Moreover \(S_{n+1}\equiv S_n\pmod {2^n}\), so these polynomial maps form a
compatible tower.  Thus the obstruction class is nonzero in the plane
correction complex but becomes exact after one identity stabilization.  This
is not merely a failure of the particular \([xy]\) detector: (7.1) cancels
the entire Jacobian error at every finite Witt level.

The degrees of \(S_n\) grow with \(n\).  Their inverse limit is the restricted
two-adic power series

\[
 (1+2K)^{-1}=\sum_{j\ge0}(-2K)^j,                              \tag{7.3}
\]

not a polynomial in \(\mathbb Z_2[x,y]\).  Accordingly, (7.1) is a compatible
formal/Witt lift; it is not a finite-degree characteristic-zero polynomial
Keller map.

## 8. Consequences and boundary of the result

The obstruction closes the direct mixed-characteristic route for this exact
plane map.  It is stronger than the earlier observation that the displayed
integer formulas themselves are not Keller: allowing arbitrary higher-degree
corrections divisible by two does not help.

It now also closes the proposed escape through polynomial plane
left--right equivalence: the obstruction is the functorial de Rham class
\(xyu^2\), not a coordinate-specific coefficient accident.  It does not
obstruct unrelated positive-characteristic Keller maps.

Stabilization at all finite Witt levels is decided by (7.1), but polynomial
algebraization over \(\mathbb Z_2\) is not.  The remaining quantitative
questions are the exact minimum degree at Witt level \(n\), lower bounds on
degree growth, possible corrections with smaller support, and impossibility
of a uniformly bounded-degree or algebraic inverse-limit family.

## 9. Exact reproduction

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_plane_w2_obstruction.py
```

The checker computes the integral Jacobian error, its full odd--odd de Rham
representative \(xyu^2\), and the old \(xy\)-coefficient witness.  It checks
the first-variation/exact-form identity, the monomial description of the full
cokernel, representative-independence on symbolic corrections, and the
explicit stabilized primitive \(zK\).  It also verifies the universal
geometric-series induction step and replays the first finite levels of the
stable Witt tower.  The all-level claim is the exact identity (7.2), not an
inference from that regression.  The checker performs no bounded correction
or left--right search.
