# A bidegree-\((4,4)\) counterexample in two contraction pairs

## 1. Statement

Work over a characteristic-zero field and use the contraction pairs

\[
 (\xi _1,z_1),\qquad(\xi _2,z_2).
\]

Put

\[
 \mathcal M_2
 =(\partial_{z_1}-\xi _1)k[\xi _1,\xi _2,z_1,z_2]
  +(\partial_{z_2}-\xi _2)k[\xi _1,\xi _2,z_1,z_2],
 \tag{1.1}
\]

and define the four bilinears

\[
\begin{aligned}
 R&=\xi _1z_1+\xi _2z_2,&
 Z&=\xi _1z_2,\\
 W&=2\xi _2z_1,&
 T&=\xi _1z_1-\xi _2z_2.
\end{aligned}
\tag{1.2}
\]

They satisfy the rank-one-quadric identity

\[
 T^2=R^2-2ZW.
\tag{1.3}
\]

Set

\[
 \boxed{
 F=(R+Z)\left(R^2W-\frac12(2R+Z)T^2\right),
 \qquad Q=Z.}
\tag{1.4}
\]

> **Theorem 1.1.** For every \(m\geq1\),
> \[
>  \boxed{
>  \mathcal E_2(F^m)=0,\qquad
>  \mathcal E_2(QF^m)
>   =\frac{(4m+2)!\,m!}{(2m+1)!!}\ne0.}
> \tag{1.5}
> \]
> Consequently \(F^m\in\mathcal M_2\) but
> \(QF^m\notin\mathcal M_2\) for every \(m\geq1\).
> Thus \(\operatorname{SIC}(2)\) is false.

The one-pair Special Image Conjecture is known, so the minimum failing pair
dimension is exactly two.  The witness \(F\) is bihomogeneous of bidegree
\((4,4)\), has ordinary total degree eight, and has sixteen expanded terms.
The multiplier \(Q\) has bidegree \((1,1)\).

In the monomial bases
\(\xi _1^i\xi _2^{4-i}\) and \(z_1^jz_2^{4-j}\), the coefficient matrix of
\(F\), with \(i,j=0,\ldots,4\), is
\[
 \begin{pmatrix}
 -1&2&0&0&0\\
 -3/2&2&6&0&0\\
 -1/2&3/2&6&6&0\\
 0&1&3/2&2&2\\
 0&0&-1/2&-3/2&-1
 \end{pmatrix},
\qquad \det=48.
\tag{1.6}
\]
In particular \(F\) has tensor rank five, so this is a genuinely
nonseparable two-pair witness rather than a constant-coefficient GVC point
on the rank-one Segre cone.

## 2. Contraction as circular Gaussian expectation

For a monomial, let

\[
 \mathcal E_2(\xi^\alpha z^\beta)
 =\partial_z^\alpha z^\beta.
\tag{2.1}
\]

Zhao's image-kernel identity gives

\[
 \mathcal M_2=\ker\mathcal E_2.
\tag{2.2}
\]

Both \(F^m\) and \(QF^m\) have equal dual and coordinate degrees.  Their
contractions are therefore scalars.  If \(G_1,G_2\) are independent
standard circular complex Gaussians, the Wick rule gives, for every
balanced polynomial \(H\),

\[
 \mathcal E_2(H)
 =\mathbb E\,
 H(\overline G_1,\overline G_2,G_1,G_2).
\tag{2.3}
\]

All coefficients in (1.4) are rational, and both sides of (1.5) are
rational scalars obtained by formal differentiation.  It therefore
suffices to prove the identities over \(\mathbb C\) using (2.3); the
resulting rational identities then extend to every characteristic-zero
field.

Write

\[
 S=|G_1|^2+|G_2|^2,\qquad
 (G_1,G_2)=\sqrt S\,(U_1,U_2),
\tag{2.4}
\]

where \(U\) is uniform on the unit sphere in \(\mathbb C^2\) and is
independent of \(S\).  The variable \(S\) has the gamma distribution of
shape two, hence

\[
 \mathbb E(S^n)=(n+1)!.
\tag{2.5}
\]

On the sphere put

\[
 x=\overline U_1U_2,\qquad
 y=2\overline U_2U_1,\qquad
 t=|U_1|^2-|U_2|^2.
\tag{2.6}
\]

Then \(t^2+2xy=1\), and radial homogeneity gives

\[
 F=S^4p,\qquad Q=Sx,
\quad
 p=(1+x)\left(y-\frac12(2+x)t^2\right).
\tag{2.7}
\]

It remains to calculate two angular moments.

## 3. The all-order angular identity

The Hopf coordinates may be chosen so that \(t\) is uniform on
\([-1,1]\), the phase of \(x\) is uniform, and

\[
 xy=\frac{1-t^2}{2}.
\tag{3.1}
\]

Phase averaging is therefore constant-term extraction after substituting
\(y=(1-t^2)/(2x)\).  Equation (2.7) becomes

\[
 p=\frac{1+x}{2x}
 \left(1-t^2(1+x)^2\right).
\tag{3.2}
\]

Since the integrand is even in \(t\), define

\[
 H_m(X)
 =X^m\int_0^1(1-s^2X^2)^m\,ds.
\tag{3.3}
\]

The phase coefficient and the \(t\)-integral give

\[
\begin{aligned}
 \mathbb E_U(p^m)
 &=2^{-m}[u^m]H_m(1+u),\\
 \mathbb E_U(xp^m)
 &=2^{-m}[u^{m-1}]H_m(1+u).
\end{aligned}
\tag{3.4}
\]

After the change of variables \(v=sX\),

\[
 H_m(X)
 =X^{m-1}J_m(X),
\qquad
 J_m(X)=\int_0^X(1-v^2)^m\,dv.
\tag{3.5}
\]

Now \(J_m'(X)=(1-X^2)^m\) has a zero of order \(m\) at \(X=1\).
Consequently the Taylor coefficients through order \(m\) of
\(H_m(X)\) at \(X=1\) are the same as those of
\(J_m(1)X^{m-1}\).  Equation (3.4) immediately yields

\[
 \mathbb E_U(p^m)=0,\qquad
 \mathbb E_U(xp^m)=2^{-m}J_m(1).
\tag{3.6}
\]

The remaining beta integral is

\[
\begin{aligned}
 J_m(1)
 &=\int_0^1(1-v^2)^m\,dv\\
 &=\frac12B\left(\frac12,m+1\right)
 =\frac{2^m m!}{(2m+1)!!}.
\end{aligned}
\tag{3.7}
\]

Thus

\[
 \mathbb E_U(p^m)=0,\qquad
 \mathbb E_U(xp^m)=\frac{m!}{(2m+1)!!}.
\tag{3.8}
\]

Combining (2.3), (2.5), (2.7), and (3.8) gives

\[
\begin{aligned}
 \mathcal E_2(F^m)
 &=(4m+1)!\,\mathbb E_U(p^m)=0,\\
 \mathcal E_2(QF^m)
 &=(4m+2)!\,\mathbb E_U(xp^m)
 =\frac{(4m+2)!\,m!}{(2m+1)!!},
\end{aligned}
\]

which proves Theorem 1.1.

## 4. Frontier consequences

The counterexample lies in balanced bidegree \((4,4)\).  It is therefore
compatible with the complete positive theorem in bidegree \((2,2)\) and
does not decide the still-open full bidegree-\((3,3)\) stratum.  It
falsifies the proposed all-degree two-pair moment--nullcone equality at
\(d=4\): the pure moments of \(F\) all vanish, while the fixed multiplier
\(Q\) proves that \(F\) cannot be one-sided.

Because the coefficient matrix has full rank, the witness does not
contradict the split-symbol GVC theorem.  That theorem concerns rank-one
forms \(A(\xi)P(z)\), whereas (1.4) is nonseparable.  The two-variable GVC
frontier and the ordinary-Laplacian polarization problem therefore remain
open in exactly their previously stated nonhomogeneous or degree-raising
forms.

## Reproduction

Run

```bash
python3 scripts/verify_two_pair_image_mathieu_counterexample.py
```

The dependency-free checker builds \(F\) from (1.2), verifies (1.3), checks
the sixteen-term expansion and coefficient-matrix determinant, and performs
exact sparse contraction through \(m=8\).  It also replays through \(m=99\)
the two finite sums obtained by expanding (3.3).  The argument in Sections
2--3, rather than either finite cutoff, proves (1.5) for every \(m\).

No term-count minimality or literature-priority claim is made here.

## Source for the Image framework

A. van den Essen, D. Wright, and W. Zhao,
[*On the Image Conjecture*](https://arxiv.org/abs/1008.3962),
J. Algebra 340 (2011), 211--224, supplies the contraction
image-kernel identity and the one-pair positive theorem used for the sharp
dimension conclusion.
