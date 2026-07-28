# Rank-three sextic reduction with simultaneous cubic and quartic layers

## Status

This note treats

\[
 \psi=q_2+h_3+h_4+h_6
\]

in four Hessian variables when
\(\operatorname{Hess}(h_6)\) has generic rank three.  The quadratic form
\(q_2\) is nondegenerate and the other forms are arbitrary homogeneous
forms of the indicated degrees.  There is no support restriction.

> **Theorem `HC4T31`.**  If
> \(\operatorname{rank}\operatorname{Hess}(h_6)=3\) generically and
> \(\det\operatorname{Hess}(\psi)\) is a nonzero constant, then
> \(\nabla\psi\) cannot identify a nonzero antipodal pair.

This is the first support-free result with cubic, quartic, and sextic
layers simultaneously present.  It does not treat sextic Hessian rank at
most two or a non-coordinate coisotropic embedding.

We work after scalar extension to an algebraic closure.  This preserves
both the collision and the constant-Hessian identity.

## 1. The two top determinant layers

Put

\[
 H_0=\operatorname{Hess}(q_2),\quad
 A=\operatorname{Hess}(h_3),\quad
 B=\operatorname{Hess}(h_4),\quad
 C=\operatorname{Hess}(h_6).
\]

Spatial scaling gives

\[
 \det(H_0+\lambda A+\lambda^2B+\lambda^4C)=\det H_0. \tag{1.1}
\]

Gordan--Noether gives a constant kernel direction \(t\) of the rank-three
matrix \(C\).  In coordinates with
\(C=0\oplus C_3\), the coefficients of degrees fourteen and thirteen in
(1.1) are

\[
 \det(C_3)B_{tt},\qquad \det(C_3)A_{tt}.
\]

Consequently

\[
 D_th_6=0,\qquad D_t^2h_4=0,\qquad D_t^2h_3=0.       \tag{1.2}
\]

After subtracting the common gradient value of the antipodal pair, write

\[
 \psi=\frac{\kappa}{2}t^2+t\,s(u)+\phi(u),\qquad
 \deg s\le3.                                         \tag{1.3}
\]

If \(\kappa\ne0\), polynomial Schur descent eliminates \(t\) and carries
the two distinct critical points to a three-variable constant-Hessian
potential.  This contradicts \(\mathrm{HC}_3\).

## 2. The cubic bordered lemma

It remains to take \(\kappa=0\).  The coefficient of \(t^2\) in the
bordered Hessian determinant is

\[
 -J(s),\qquad
 J(s)=(\nabla s)^{\mathsf T}
       \operatorname{adj}(\operatorname{Hess}s)\nabla s. \tag{2.1}
\]

Write

\[
 s=c+\ell+a_2+a_3
\]

with homogeneous pieces.  We need the following small classification.

> **Cubic bordered lemma.**  If \(s\) is a polynomial of degree at most
> three in three variables and \(J(s)=0\), then \(s\) is independent of
> some nonzero constant direction.

Here is an exact proof.  The degree-six part of (2.1) is

\[
 \frac32a_3\det\operatorname{Hess}(a_3).
\]

Thus the ternary cubic \(a_3\) has singular Hessian.  Gordan--Noether
makes it binary.  Over the algebraic closure, its nonzero binary orbit is
represented by exactly one of

\[
 xy(x+y),\qquad x^2y,\qquad x^3.                    \tag{2.2}
\]

Write the remaining variable as \(m\), and write the quadratic and linear
coefficients involving it as

\[
 q_{xm}xm+q_{ym}ym+\ell_m m.
\]

For either rank-two orbit in (2.2), the complete coefficient ideal of
\(J(s)\) is

\[
 (q_{xm},q_{ym},\ell_m)^2.                           \tag{2.3}
\]

Hence \(s\) is independent of \(m\).

For the rank-one orbit \(x^3\), put

\[
 a=q_{xy},\quad b=q_{xm},\quad c=q_{yy},\quad
 d=q_{ym},\quad p=\ell_y,\quad q=\ell_m.
\]

The radical of the coefficient ideal is the intersection of

\[
 P_1=(d,c,bp-aq),\qquad P_2=(q,d,b).                 \tag{2.4}
\]

On \(P_2\), the polynomial is independent of \(m\).  On \(P_1\), its
\((y,m)\)-dependence is

\[
 (p,q)\binom ym+x(a,b)\binom ym,
\]

and \(bp=aq\) makes the two coefficient rows proportional.  A nonzero
constant direction annihilates both.  The checker proves (2.4) without
trusting a primary decomposition: it verifies both ideal inclusions and
the explicit radical powers \(2,3,2,3\).

Finally, if \(a_3=0\), the leading part of (2.1) gives
\(\det\operatorname{Hess}(a_2)=0\).  In rank two, the constant part
\(\ell^{\mathsf T}\operatorname{adj}(\operatorname{Hess}a_2)\ell=0\)
makes the one-dimensional quadratic kernel orthogonal to \(\ell\).
In ranks one and zero, intersect the quadratic kernel with
\(\ker\ell\).  This proves the lemma in every case.

## 3. The second descent

Choose the direction supplied by the cubic bordered lemma as \(m\), so

\[
 \psi=t\,s(x,y)+\phi(x,y,m).
\]

The coefficient of \(t\) factors as

\[
 -\phi_{mm}
 (\nabla s)^{\mathsf T}
 \operatorname{adj}(\operatorname{Hess}_{x,y}s)\nabla s. \tag{3.1}
\]

If \(\phi_{mm}=0\), then

\[
 \psi=t\,s(x,y)+m\,g(x,y)+h(x,y)
\]

and

\[
 \det\operatorname{Hess}(\psi)=\operatorname{Jac}(s,g)^2. \tag{3.2}
\]

Here \(\deg s\le3\) and \(\deg g\le5\).  Moh's plane degree bound makes
\((s,g)\) a polynomial automorphism, and the cotangent lift (3.2) is also
an automorphism.  It has no collision.

Otherwise the binary bordered invariant in (3.1) vanishes.  The
two-variable singular-polynomial-Hessian classification applied to
\(t\,s(x,y)\) says that \(s\) is a polynomial in one linear form.  If its
degree were at least two, its gradient would vanish somewhere over the
algebraic closure, making the first bordered Hessian row zero.  Therefore
\(s\) is linear.  Normalize it to \(s=c+x\).

Both antipodal critical points have \(s=0\), hence \(c=x=0\).  Expansion
along the \(t,x\) block gives

\[
 \det\operatorname{Hess}(\psi)
 =-\det\operatorname{Hess}_{y,m}\phi.
\]

At \(x=0\), \(\mathrm{HC}_2\) makes the remaining gradient injective.
Antipodality forces \(y=m=0\), and the last critical equation fixes \(t\)
uniquely.  The two points coincide, which is the final contradiction.

## Reproduction

Run:

```bash
.venv/bin/python \
  scripts/verify_hc4_meng_triple_rank_three_reduction.py
```

The checker verifies the top determinant layers, the bordered identity,
the two rank-two ideals (2.3), both radical inclusions and explicit powers
in (2.4), and the quadratic-leading boundary.  It also replays the shared
cotangent-lift and terminal \(\mathrm{HC}_2\) identities.

The external structural inputs are the
[Gordan--Noether classification](https://arxiv.org/abs/1501.05168),
the de Bondt--van den Essen
[singular-Hessian classification](https://www.sciencedirect.com/science/article/pii/S0021869304004867),
the Hessian conjecture in dimensions at most three, and
[Moh's plane degree bound](https://www.math.purdue.edu/~ttm/jacobian.pdf).
