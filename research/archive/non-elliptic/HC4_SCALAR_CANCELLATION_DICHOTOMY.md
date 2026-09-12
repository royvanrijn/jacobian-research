# Scalar cancellation dichotomy for reverse HC4 descent

## Status

This note separates scalar exact cancellation into its genuine
quadratic-pivot strata and closes every zero-corner case. In the nonzero
corner it closes ranks three and four and reduces ranks two and one to HC2
or the exact JC2 cotangent packet.

> **Theorem HC4RSD11 (nonzero-corner pencil equivalence).** Let \(K\) have
> characteristic zero and
>
> \[
> \Phi(t,x)=\frac{\lambda}{2}t^2+tA(x)+B(x),
> \qquad \lambda\in K^\times,
> \quad x\in K^4.
> \tag{0.1}
> \]
>
> Put
>
> \[
> \psi(x)=B(x)-\frac{A(x)^2}{2\lambda}.
> \tag{0.2}
> \]
>
> Then
>
> \[
> \det\operatorname{Hess}_{t,x}\Phi=c
> \quad\Longleftrightarrow\quad
> \det\operatorname{Hess}_x(\psi+sA)=c/\lambda
> \quad\text{for every }s.
> \tag{0.3}
> \]
>
> Equal-gradient pairs of \(\Phi\) are in bijection with equal-gradient
> pairs of the appropriate pencil member \(\psi+sA\). Conversely, every
> constant-Hessian pencil \(\psi+sA\) has the suspension
>
> \[
> \Phi=\frac{\lambda}{2}
> \left(t+\frac{A}{\lambda}\right)^2+\psi.
> \tag{0.4}
> \]
>
> Thus nonzero-corner determinant-term cancellation is exactly the problem
> of classifying four-variable constant-Hessian pencils. Every direct HC4
> candidate embeds by taking \(A\) affine, while nonlinear \(A\) defines a
> smaller structured pencil locus.

> **Theorem HC4RSD12 (zero-corner graph-coordinate obstruction).** Let
>
> \[
> \Phi(t,x)=tA(x)+B(x)
> \tag{0.5}
> \]
>
> have nonzero constant five-variable Hessian determinant. Suppose there is
> a constant direction \(v\in K^4\) with
>
> \[
> D_vA\in K^\times.
> \tag{0.6}
> \]
>
> Then \(\nabla\Phi\) is a polynomial automorphism. In particular \(\Phi\)
> has no gradient collision, independently of whether
> \(\operatorname{Hess}(B+sA)\) is singular or survives by an exact
> specialized determinant remainder.

> Every zero-corner scalar parent with \(\deg A\le2\) satisfies (0.6).
> Consequently no five-variable Hessian collision can pass through a
> quadratic zero-corner scalar pivot.

> **Theorem HC4RSD13 (rank-three nonlinear pencil classification).** Let
> \(\psi,A\in K[x_1,\ldots,x_4]\), with \(\deg A\le2\),
> \(\operatorname{rank}\operatorname{Hess}A=3\), and
>
> \[
> \det\operatorname{Hess}(\psi+sA)=\delta\in K^\times
> \quad\text{for every }s.
> \tag{0.7}
> \]
>
> Then every member of the pencil has a polynomially invertible gradient.
> More precisely, after scalar extension and constant affine changes,
>
> \[
> A=xy+\frac12z^2,
> \qquad
> \psi=wx+y\bigl(\alpha z+\beta(x)\bigr)+G(x,z),
> \tag{0.8}
> \]
>
> where \(\alpha\in K^\times\). Its Hessian determinant is
> \(\alpha^2\), independently of \(s,\beta,G\).

> **Theorem HC4RSD14 (rank-two pencil reduction to JC2).** Under (0.7),
> suppose \(\operatorname{rank}\operatorname{Hess}A=2\). After scalar
> extension and constant affine changes, put \(A=xy\) and let
> \(E=\operatorname{Hess}_{z,w}\psi\). Then \(\det E=0\).
>
> If \(\operatorname{rank}E=1\), every member of the pencil has injective
> gradient. After normalization it has the form
>
> \[
> \psi=xz+C(x,y,w),
> \qquad
> \det\operatorname{Hess}_{y,w}C\in K^\times.
> \tag{0.9}
> \]
>
> If \(E=0\), the complete form is
>
> \[
> \psi=z\,b(x,y)+w\,c(x,y)+D(x,y),
> \qquad
> \det J(b,c)\in K^\times.
> \tag{0.10}
> \]
>
> The gradient of \(\psi+sA\) is a polynomial automorphism if and only if
> the plane Keller map \((b,c)\) is a polynomial automorphism. Thus the sole
> rank-two survivor is exactly the cotangent lift of JC2.

> **Theorem HC4RSD15 (rank-one constant-chart reduction).** Under (0.7),
> suppose \(\operatorname{rank}\operatorname{Hess}A=1\). Normalize
> \(A=x^2/2\) and put \(E=\operatorname{Hess}_{y,z,w}\psi\). Then \(E\)
> has generic rank two. If its ternary singular-Hessian normal form can be
> chosen in a passive chart independent of \(x\), exactly two packets occur:
>
> 1. the constant-kernel packet
>    \[
>    \psi=xw+C(x,y,z),
>    \qquad
>    \det\operatorname{Hess}_{y,z}C\in K^\times,
>    \tag{0.11}
>    \]
>    whose pencil gradients are injective by HC2;
> 2. the exceptional packet
>    \[
>    \psi=zP(x,y)+wQ(x,y)+R(x,y),
>    \qquad
>    \det J(P,Q)\in K^\times,
>    \tag{0.12}
>    \]
>    which is exactly the cotangent lift of a plane Keller map.
>
> Therefore the only rank-one packet not reduced to HC2 or JC2 is the
> globalization locus where the ternary de Bondt--van den Essen normal-form
> chart moves rationally with the active variable \(x\).

> **Theorem HC4RSD16 (unit-transverse globalization obstruction).** In the
> setting of HC4RSD15, the rationally \(x\)-moving chart is empty. In the
> constant-kernel type, the primitive passive kernel \(v(x)\) satisfies
>
> \[
> \nabla_u(v^{\mathsf T}d)=E_xv=-Ev',
> \tag{0.13}
> \]
>
> so the unit identity (7.4) forces \(v'\in\ker E\) and freezes the
> projective kernel. In the exceptional type, write invariantly over \(K(x)\)
>
> \[
> \psi=u^{\mathsf T}g(x,\ell)+R(x,\ell),
> \qquad \ell=p(x)^{\mathsf T}u.
> \tag{0.14}
> \]
>
> With \(q=g_\ell\) and \(v=p\mathbin{\times}q\), one has
>
> \[
> \operatorname{adj}E=-vv^{\mathsf T},
> \qquad
> v^{\mathsf T}d=v^{\mathsf T}g_x+
> (u^{\mathsf T}q+R_\ell)v^{\mathsf T}p'.
> \tag{0.15}
> \]
>
> Constancy along every affine level of \(\ell\) gives
> \(v^{\mathsf T}p'=0\). Either \(q\bmod p\) has fixed projective direction,
> which is already the constant-kernel type, or \(p'\parallel p\). Hence the
> distinguished exceptional covector is projectively constant and (0.12)
> applies in constant passive coordinates. Thus every rank-one quadratic
> pencil reduces to HC2 or exactly the JC2 cotangent packet.

Together with HC4RSD8--HC4RSD10, the zero-corner conclusion is complete for
quadratic scalar pivots. The singular-pencil descendants are triangular and
every zero-corner quadratic parent is injective. The remaining scalar
quadratic problem is the nonzero-corner constant-Hessian pencil identified
by HC4RSD11. HC4RSD13 closes its pivot-Hessian ranks three and four;
HC4RSD14 reduces rank two exactly to JC2, HC4RSD15 identifies the two
rank-one fiberwise charts, and HC4RSD16 proves that unit transversality
freezes both. Thus every nonlinear quadratic scalar pencil reduces to HC2
or exactly the JC2 cotangent endpoint.

The universal determinant identities are replayed by
[scripts/verify_hc4_scalar_cancellation_dichotomy.py](scripts/verify_hc4_scalar_cancellation_dichotomy.py),
which writes
[artifacts/generated-results/hc4_scalar_cancellation_dichotomy.json](artifacts/generated-results/hc4_scalar_cancellation_dichotomy.json).
That artifact records the `HC4RSD11--16` stage.  Its higher-degree pencil
`open_frontier` is historical and is superseded inside the auxiliary
relative-nilpotent branch by `HC4MR1`.

## 1. The nonzero corner is a completed-square gauge

Put

\[
g=\nabla A,
\qquad
M(t)=\operatorname{Hess}(B+tA).
\]

The parent Hessian is

\[
\begin{pmatrix}
\lambda&g^{\mathsf T}\\
g&M(t)
\end{pmatrix}.
\tag{1.1}
\]

Fix the pivot-gradient output

\[
y=\partial_t\Phi=\lambda t+A.
\]

Set

\[
s=\frac{y}{\lambda}=t+\frac{A}{\lambda}.
\]

Then direct differentiation of (0.2) gives

\[
\nabla(\psi+sA)=\nabla B+t\nabla A,
\qquad
\operatorname{Hess}(\psi+sA)=M(t)-\lambda^{-1}gg^{\mathsf T}.
\tag{1.2}
\]

The exact block identity

\[
\det
\begin{pmatrix}
\lambda&g^{\mathsf T}\\
g&M
\end{pmatrix}
=\lambda\det(M-\lambda^{-1}gg^{\mathsf T})
\tag{1.3}
\]

proves (0.3), because \((t,x)\mapsto(s,x)\) is a polynomial automorphism.
It also proves collision equivalence: equality of parent gradients first
fixes \(y\), and the other four parent-gradient coordinates are exactly
\(\nabla(\psi+sA)\). A collision of that pencil member lifts by taking

\[
t_i=\frac{y-A(x_i)}{\lambda}.
\tag{1.4}
\]

Conversely, if the whole pencil has constant Hessian determinant, (0.4)
gives

\[
\partial_t\Phi=\lambda t+A,
\qquad
\nabla_x\Phi=\nabla\psi+
\left(t+\frac{A}{\lambda}\right)\nabla A.
\]

A collision of \(\nabla(\psi+s_0A)\) lifts on the level
\(t=s_0-A/\lambda\). Taking \(A\) affine makes the Hessian pencil constant
for every \(\psi\), so this branch contains direct HC4. For nonlinear \(A\),
the all-\(s\) determinant identity is additional structure and remains a
classification target.

## 2. Exact graph-coordinate factorization

Normalize (0.6) by constant affine coordinates and scaling:

\[
x=(u,w),
\qquad
u=(u_1,u_2,u_3),
\qquad
A=w+q(u).
\tag{2.1}
\]

Use the polynomial coordinate

\[
r=A=w+q(u)
\]

and write

\[
C(u,r)=B(u,r-q(u)).
\tag{2.2}
\]

Finally set

\[
\tau=t+C_r(u,r).
\tag{2.3}
\]

In the original variable order \((t,u,w)\), the five parent-gradient
coordinates become

\[
\boxed{
\partial_t\Phi=r,
\qquad
\nabla_u^{\rm old}\Phi
=\nabla_u C(u,r)+\tau\nabla q(u),
\qquad
\partial_w\Phi=\tau.
}
\tag{2.4}
\]

Both changes

\[
(t,u,w)\longleftrightarrow(r,u,\tau)
\]

are triangular polynomial changes with constant Jacobian up to sign. Thus
the Jacobian determinant of (2.4), equivalently the parent Hessian
determinant, is

\[
\boxed{
\det\operatorname{Hess}\Phi
=-\det\operatorname{Hess}_u(C(u,r)+\tau q(u)).
}
\tag{2.5}
\]

This can also be checked directly from the universal Hessian block. If

\[
p=\nabla q,
\quad Q=\operatorname{Hess}q,
\quad H=\operatorname{Hess}_uC,
\quad b=\nabla_uC_r,
\quad d=C_{rr},
\]

then the parent Hessian is

\[
\begin{pmatrix}
0&p^{\mathsf T}&1\\
p&H+bp^{\mathsf T}+pb^{\mathsf T}+dpp^{\mathsf T}+\tau Q&b+dp\\
1&(b+dp)^{\mathsf T}&d
\end{pmatrix},
\tag{2.6}
\]

whose determinant is exactly \(-\det(H+\tau Q)\).

## 3. HC3 makes every graph-coordinate parent injective

If the parent determinant is \(c\in K^\times\), (2.5) says that for every
fixed \(r,\tau\), the ternary potential

\[
C_{r,\tau}(u)=C(u,r)+\tau q(u)
\tag{3.1}
\]

has Hessian determinant \(-c\). By HC3 its gradient is injective over the
algebraic closure.

Suppose two parent gradients agree. Equation (2.4) first gives equal \(r\)
and equal \(\tau\). The middle three equations are then equal gradients of
the same \(C_{r,\tau}\), so HC3 gives equal \(u\). Finally

\[
w=r-q(u),
\qquad
t=\tau-C_r(u,r)
\tag{3.2}
\]

recover the remaining variables. Thus the full five-variable gradient is
injective. Ax--Grothendieck makes it a polynomial automorphism, and faithful
flatness descends the conclusion from the algebraic closure to \(K\).

## 4. Why every quadratic zero-corner pivot is covered

Let

\[
A=\frac12x^{\mathsf T}Qx+a^{\mathsf T}x+a_0,
\qquad Q=Q^{\mathsf T}.
\tag{4.1}
\]

For the zero-corner parent, with \(M(t)=\operatorname{Hess}(B+tA)\),

\[
c=-\nabla A^{\mathsf T}\operatorname{adj}(M(t))\nabla A.
\tag{4.2}
\]

The right side writes the nonzero constant \(c\) as a polynomial
combination of the entries of \(\nabla A=Qx+a\). Hence these affine entries
generate the unit ideal, so \(Qx+a=0\) has no solution over the algebraic
closure. Therefore \(Q\) is singular and \(a\notin\operatorname{im}Q\).
Since \(Q\) is symmetric, there is \(v\in\ker Q\) with

\[
D_vA=v^{\mathsf T}a\ne0.
\tag{4.3}
\]

This is precisely (0.6). No singularity, corank, or exact-remainder
assumption on the reduced pencil was used.

## 5. Rank-three nonlinear constant-Hessian pencils

Write \(Q=\operatorname{Hess}A\). Rank four is impossible in (0.7), since
the coefficient of \(s^4\) is \(\det Q\). Suppose \(\operatorname{rank}Q=3\).
After scalar extension and a constant congruence, split the variables as
\((u,w)\), with \(u\in K^3\), and write

\[
\operatorname{Hess}\psi=
\begin{pmatrix}
K&d\\
d^{\mathsf T}&e
\end{pmatrix},
\qquad
\operatorname{Hess}A=
\begin{pmatrix}
Q&0\\
0&0
\end{pmatrix},
\quad \det Q\ne0.
\tag{5.1}
\]

The coefficient of \(s^3\) in (0.7) is \((\det Q)e\), so \(e=0\).
Integration gives

\[
\psi=w\,b(u)+C(u),
\qquad
d=\nabla b,
\qquad
K=\operatorname{Hess}C+w\operatorname{Hess}b.
\tag{5.2}
\]

Now

\[
\det\operatorname{Hess}(\psi+sA)
=-d^{\mathsf T}\operatorname{adj}(K+sQ)d.
\tag{5.3}
\]

The coefficient of \(s^2\) gives the eikonal equation

\[
d^{\mathsf T}\operatorname{adj}(Q)d=0.
\tag{5.4}
\]

The constant term in (5.3) is the nonzero unit \(\delta\), so the three
components of \(d\) generate the unit ideal.

We use the following elementary unimodular eikonal lemma.

> If a polynomial \(b\) in three variables has unimodular gradient and
> \(\nabla b\) is null for a nondegenerate constant quadratic form, then
> \(b\) is affine.

To prove it, use null coordinates \((\xi,\eta,z)\), so the equation is

\[
b_\xi b_\eta+b_z^2=0.
\tag{5.5}
\]

If none of the derivatives vanishes, primitive factorization in the UFD
\(K[\xi,\eta,z]\) gives

\[
b_\xi=a^2,\qquad b_\eta=-h^2,\qquad b_z=ah,
\qquad \gcd(a,h)=1.
\tag{5.6}
\]

The three mixed-partial equations imply successively

\[
a_\xi=0,\qquad h_\eta=0,\qquad
h_\xi=2a_z,\qquad a_\eta=-2h_z,\qquad
(ah)_z=0.
\tag{5.7}
\]

Coprimality in the last equation gives \(a_z=h_z=0\), hence also
\(h_\xi=a_\eta=0\). Thus \(a,h\) are constants. The cases with a zero
derivative are even more immediate because unimodularity makes the
remaining derivative a unit. This proves the lemma and makes \(d\) a
nonzero constant isotropic vector.

The orthogonal group of \(Q\) is transitive on its nonzero null vectors.
Normalize

\[
A=xy+\frac12z^2,
\qquad
b=x.
\tag{5.8}
\]

Writing \(\psi=wx+C(x,y,z)\), equation (5.3) becomes

\[
\det\operatorname{Hess}(\psi+sA)
=C_{yz}^2-C_{yy}C_{zz}-sC_{yy}.
\tag{5.9}
\]

Therefore \(C_{yy}=0\) and \(C_{yz}^2=\delta\). Exact integration gives

\[
C=y\bigl(\alpha z+\beta(x)\bigr)+G(x,z),
\qquad \alpha^2=\delta,
\tag{5.10}
\]

which is (0.8). For \(F=\nabla(\psi+sA)\), recover

\[
\begin{aligned}
x&=F_w,\\
z&=\alpha^{-1}\bigl(F_y-\beta(x)-sx\bigr),\\
y&=\alpha^{-1}\bigl(F_z-G_z(x,z)-sz\bigr),\\
w&=F_x-y\beta'(x)-G_x(x,z)-sy.
\end{aligned}
\tag{5.11}
\]

Thus every pencil member is a polynomial automorphism.

## 6. Rank-two pencils and the plane-Jacobian endpoint

Normalize \(A=xy\), with active variables \((x,y)\) and passive variables
\((z,w)\). Write

\[
\operatorname{Hess}(\psi+sA)
=
\begin{pmatrix}
K+sQ&D\\
D^{\mathsf T}&E
\end{pmatrix},
\qquad
Q=\begin{pmatrix}0&1\\1&0\end{pmatrix}.
\tag{6.1}
\]

The coefficient of \(s^2\) is \(-\det E\), so \(\det E=0\).

### 6.1 Passive rank zero

If \(E=0\), exact integration gives

\[
\psi=z\,b(x,y)+w\,c(x,y)+D_0(x,y).
\tag{6.2}
\]

With \(J=J(b,c)\), direct block expansion gives

\[
\det\operatorname{Hess}(\psi+sA)=(\det J)^2.
\tag{6.3}
\]

Thus \(G=(b,c)\) is an arbitrary plane Keller map. The last two gradient
coordinates of \(\psi+sA\) are exactly \(G(x,y)\). If \(G\) is invertible,
they recover \((x,y)\), and the first two coordinates then recover
\((z,w)\) through the invertible matrix \(J^{\mathsf T}\).

Conversely, if \(G(p)=G(q)\) with \(p\ne q\), choose the two fiber vectors
\((z,w)\) over \(p\) and \(q\) so that the first two gradient coordinates
have one common prescribed value. This is possible because both Jacobian
matrices are invertible. The cotangent lift then has a gradient collision.
Therefore

\[
\nabla(\psi+sA)\text{ is invertible}
\quad\Longleftrightarrow\quad
G\text{ is invertible}.
\tag{6.4}
\]

This is an exact equivalence with JC2, not an exclusion.

### 6.2 Passive rank one

Apply the rank-one polynomial-Hessian normal form over \(K(x,y)\), followed
by primitive polynomial factorization:

\[
E=\rho\,\ell\ell^{\mathsf T},
\qquad
\ell=(p(x,y),r(x,y)).
\tag{6.5}
\]

Write \(d_x,d_y\) for the two rows of the active/passive block and put

\[
\Delta_i=\det(d_i,\ell).
\]

The coefficient of \(s\) in the pencil determinant is

\[
[s]\det\operatorname{Hess}(\psi+sA)
=2\rho\,\Delta_x\Delta_y.
\tag{6.6}
\]

Choose the branch \(\Delta_y=0\). If the passive gradient is written as

\[
\nabla_{z,w}\psi=\varphi\,\ell+a(x,y),
\tag{6.7}
\]

then \(\varphi\) genuinely depends on the passive variables and

\[
\Delta_i
=\varphi\det(\partial_i\ell,\ell)
+ \det(\partial_i a,\ell).
\tag{6.8}
\]

Equation \(\Delta_y=0\) first fixes the projective direction of \(\ell\)
along \(y\). Since \(d_y=h\ell\) for a polynomial \(h\), the constant pencil
term factors exactly as

\[
\delta
=\bigl(h^2-\rho\psi_{yy}\bigr)\Delta_x^2.
\tag{6.9}
\]

Both factors are polynomials, so \(\Delta_x\) is a unit. Equation (6.8)
then fixes the projective direction of \(\ell\) along \(x\) as well, and
the remaining common factor is a unit. Hence \(\ell\) is constant.

Normalize \(\ell=(0,1)\). Equations (6.8)--(6.9) now give

\[
\psi=z\,b(x)+C(x,y,w),
\qquad b'(x)\in K^\times.
\]

After rescaling and removing affine terms,

\[
\psi=xz+C(x,y,w),
\tag{6.10}
\]

and the pencil identity becomes

\[
\det\operatorname{Hess}(\psi+sxy)
=C_{yw}^2-C_{yy}C_{ww}
=\delta.
\tag{6.11}
\]

For two equal gradients, \(F_z=x\) first fixes \(x\). The pair
\((F_y,F_w)\), after subtracting the known constant \(sx\), is the gradient
of the binary potential \(C(x,\cdot,\cdot)\), whose Hessian determinant is
\(-\delta\). HC2 fixes \((y,w)\), and \(F_x=z+C_x+sy\) finally fixes \(z\).
Thus the passive-rank-one branch is injective.

## 7. Rank-one pencils and ternary singular-Hessian charts

Normalize \(A=x^2/2\) and write

\[
\operatorname{Hess}\psi=
\begin{pmatrix}
k&d^{\mathsf T}\\
d&E
\end{pmatrix},
\qquad E=\operatorname{Hess}_{y,z,w}\psi.
\tag{7.1}
\]

The pencil determinant is

\[
\det\operatorname{Hess}(\psi+sA)
=\det\operatorname{Hess}\psi+s\det E.
\tag{7.2}
\]

Hence \(\det E=0\). Since the full determinant is nonzero, \(E\) cannot
have rank at most one: one border raises rank by at most two. Thus
\(\operatorname{rank}E=2\). Over the fraction field, write

\[
\operatorname{adj}(E)=\rho\,vv^{\mathsf T}.
\tag{7.3}
\]

The constant determinant is

\[
\boxed{\delta=-\rho\,(v^{\mathsf T}d)^2.}
\tag{7.4}
\]

After primitive polynomial factorization, both \(\rho\) and
\(v^{\mathsf T}d\) are polynomial units. This is the exact globalization
gate on the ternary singular-Hessian classification.

The characteristic-zero de Bondt--van den Essen theorem, in the form
recorded in
[de Bondt's three-variable Hessian paper](https://arxiv.org/abs/1203.6605),
applied over \(K(x)\), has two normal-form types.

### 7.1 Constant kernel chart

If the kernel direction can be chosen constant in the passive variables
and in \(x\), normalize it to \(\partial_w\). Integration gives

\[
\psi=w\,b(x)+C(x,y,z).
\tag{7.5}
\]

Equation (7.4) becomes

\[
\delta=-b'(x)^2\det\operatorname{Hess}_{y,z}C.
\tag{7.6}
\]

Both factors are units. Normalize \(b=x\). The last gradient coordinate
recovers \(x\); HC2 applied to \(C(x,\cdot,\cdot)\) recovers \(y,z\); and
the first coordinate recovers \(w\). This proves (0.11).

### 7.2 Constant exceptional chart

The exceptional ternary zero-Hessian form is linear in two passive
variables after choosing one passive linear form. If that chart is
independent of \(x\), its complete four-variable form is

\[
\psi=zP(x,y)+wQ(x,y)+R(x,y).
\tag{7.7}
\]

Direct expansion gives

\[
\det\operatorname{Hess}\psi=\det J(P,Q)^2.
\tag{7.8}
\]

This is the same cotangent-lift equivalence as (6.2)--(6.4). It is
invertible exactly when the plane Keller map \((P,Q)\) is invertible.

Over \(K(x)\), every rank-two ternary singular Hessian has one of these two
forms. What remains is not a third fiberwise normal form: it is the
globalization problem in which the required passive linear chart has
rational \(x\)-dependence. Equation (7.4), Hessian integrability, and the
Piola identity are the exact constraints on that moving chart.

### 7.3 The unit gate freezes the rational chart

We now solve that globalization problem. First suppose that the passive
kernel is independent of \(u=(y,z,w)\), but is represented over \(K(x)\) by
a vector \(v(x)\). Clear denominators and divide by the gcd of its entries.
Because \(K[x]\) is a PID, the resulting primitive vector is unimodular and

\[
\operatorname{adj}E=\rho vv^{\mathsf T}
\tag{7.9}
\]

with \(\rho\in K[x,u]\). Equation (7.4) is a factorization of a unit in the
UFD \(K[x,u]\). Hence

\[
\rho\in K^\times,
\qquad
v^{\mathsf T}d\in K^\times.
\tag{7.10}
\]

Since \(v\) is independent of \(u\), differentiating the second identity in
the passive variables and differentiating \(Ev=0\) in \(x\) give

\[
0=\nabla_u(v^{\mathsf T}d)=E_xv=-Ev'.
\tag{7.11}
\]

Thus \(v'\) lies in the one-dimensional kernel of \(E\), so the projective
class of \(v\) is constant. Primitivity removes the remaining scalar factor.
This reduces the apparently moving constant-kernel packet to Section 7.1.

For the exceptional type, the ternary classification over \(K(x)\) can be
written without choosing its two transverse coordinates. There are
\(p(x)\in K(x)^3\), \(g(x,\ell)\in K(x)[\ell]^3\), and
\(R(x,\ell)\in K(x)[\ell]\) such that

\[
\psi=u^{\mathsf T}g(x,\ell)+R(x,\ell),
\qquad \ell=p(x)^{\mathsf T}u.
\tag{7.12}
\]

A component of \(g\) parallel to \(p\) may be absorbed in \(R\). Put
\(q=g_\ell\). Direct differentiation gives

\[
E=pq^{\mathsf T}+qp^{\mathsf T}
 +(u^{\mathsf T}g_{\ell\ell}+R_{\ell\ell})pp^{\mathsf T}.
\tag{7.13}
\]

The \(3\)-by-\(3\) adjugate and the bordered factor are therefore

\[
\operatorname{adj}E
=-(p\mathbin{\times}q)(p\mathbin{\times}q)^{\mathsf T},
\tag{7.14}
\]

\[
(p\mathbin{\times}q)^{\mathsf T}d
=(p\mathbin{\times}q)^{\mathsf T}g_x
 +(u^{\mathsf T}q+R_\ell)
   (p\mathbin{\times}q)^{\mathsf T}p'.
\tag{7.15}
\]

Here \(g_x\) means the partial derivative with \(\ell\) held fixed. The
square of (7.15) is the nonzero constant four-variable Hessian determinant.
Vary \(u\) by \(r\in\ker p\), keeping \(\ell\) fixed. The change in (7.15)
is exactly

\[
(r^{\mathsf T}q)(p\mathbin{\times}q)^{\mathsf T}p'.
\tag{7.16}
\]

Rank two supplies an \(r\in\ker p\) with \(r^{\mathsf T}q\ne0\). Hence
\((p\mathbin{\times}q)^{\mathsf T}p'=0\), or

\[
p'\in\operatorname{span}_{K(x,\ell)}\{p,q\}.
\tag{7.17}
\]

If \(p'\) is not parallel to \(p\), (7.17) forces every \(q(x,\ell)\) into
the plane spanned by \(p,p'\). The passive kernel
\(p\mathbin{\times}q\) then has direction independent of \(u\), so this is
the constant-kernel type already closed above. On the genuinely exceptional
type, therefore, \(p'\parallel p\); its projective class is constant. Absorb
the rational scalar into \(\ell\) and take \(p\) constant. Equation (7.12)
is then affine in two fixed transverse passive coordinates. Since the
original \(\psi\) is polynomial, coefficient extraction gives polynomial
\(P,Q,R\), and the fixed form (7.7) follows.

This proves HC4RSD16. There is no moving rank-one packet: rank one is HC2
or the exact JC2 cotangent lift.

## 8. Revised scalar frontier

The zero-corner quadratic scalar-pivot mechanism is exhausted:

- singular zero-corner descendants are explicitly triangular by
  HC4RSD8--HC4RSD10;
- every zero-corner quadratic parent is collision-free by HC4RSD12, even
  when its reduced pencil is nonsingular;
- every nonzero-corner cancellation is a four-variable constant-Hessian
  pencil by HC4RSD11. Its affine directions are exactly direct HC4;
  rank-four quadratic directions are impossible and rank three is
  triangular by HC4RSD13. In rank two, HC4RSD14 closes passive rank one and
  identifies passive rank zero exactly with the cotangent lift of JC2.
  In rank one, HC4RSD15 reduces both constant ternary normal-form charts to
  HC2 or the same JC2 cotangent packet, and HC4RSD16 proves that the bordered
  unit freezes every rationally \(x\)-moving chart.

The live reverse-descent mechanisms are therefore nonlinear
constant-Hessian pencils with pivot direction of degree at least three in
the nonzero-corner scalar branch, moving matrix-pivot planes, and genuinely
mixed/coisotropic canonical transformations. For quadratic scalar pivots,
the only unresolved endpoint is exactly JC2, not an HC4-specific component.
Direct degree-five HC4 classification remains the parallel route.

## 9. Reproduction

Run:

~~~bash
.venv/bin/python scripts/verify_hc4_scalar_cancellation_dichotomy.py
# committed `HC4RSD11--16` stage artifact only, without symbolic replay:
.venv/bin/python scripts/verify_hc4_scalar_cancellation_dichotomy.py --audit-existing-only
~~~

The maintenance-only mode verifies the committed artifact and explicitly
reports that its higher-degree pencil handoff is stage-local; it neither
recomputes nor rewrites the identities.

The checker verifies the universal \(5\)-by-\(5\) block determinant, the
graph-coordinate Hessian factorization, exact gradient coordinates, a
degree-unbounded nonlinear calibration, and the nonzero-corner suspension
identity. It also verifies the rank-three determinant faces and triangular
normal form, the rank-two channel factorization, the safe passive-rank-one
normal form, the rank-one constant-chart factors, the exact cotangent-lift
determinant, the exceptional moving-chart adjugate and border factor, and a
moving-kernel calibration exposing the forbidden nonconstant gate.
