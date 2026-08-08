# Polynomial-graph closure for the final rank-three `[4]` HC4 stratum

## Status

This note continues `HC4RSD71`.  It treats the case where the smooth generic
part of the gradient-image hypersurface is an affine polynomial graph

\[
Y:\quad y_4=\mathcal H(y_1,y_2,y_3),\qquad \mathcal H\in K[p_1,p_2,p_3].
\]

> **Theorem HC4RSD72 — polynomial graph closure.**  In the final rank-three
> `[4]` relative-nilpotent HC4 packet, if after a constant affine target change
> the gradient-image hypersurface is a polynomial graph, then the top kernel
> has a constant linear relation.  Hence the packet reduces to the already
> closed linearly-dependent/fixed-kernel branch.
>
> In particular, a genuinely unresolved HC4 packet cannot have a polynomial
> graph as its gradient image in any affine target direction.

The proof is degree-free and does not assume that the quotient map to the
gradient image is birational.

## 1. Slice supplied by a polynomial graph

Write

\[
F=\nabla A=(p_1,p_2,p_3,\mathcal H(p)).
\]

For

\[
g(y)=y_4-\mathcal H(y_1,y_2,y_3)
\]

the associated kernel field is

\[
k=\nabla g(F)=(-\nabla\mathcal H(p),1).
\]

Because \(Tk=0\), where \(T=\operatorname{Hess}A\),

\[
D:=k\cdot\nabla,\qquad Dk=0,
\]

and \(Dx_4=1\).  Hence the final orbit coordinate is an explicit polynomial
slice.  Put

\[
u_i=x_i-x_4k_i(x),\qquad i=1,2,3.
\]

Then \(Du_i=0\), and invariance of \(k\) along its own flow gives the exact
polynomial inverse

\[
x_i=u_i+x_4k_i(u,0),\qquad x_4=t.
\]

On the slice \(t=0\), set

\[
A_0(u)=A(u,0),\qquad p=\nabla A_0(u),\qquad
K=\operatorname{Hess}A_0.
\]

The horizontal orbit velocity is

\[
\kappa(u)=-\nabla\mathcal H(p(u)).
\]

Since \((u,t)\mapsto(u+t\kappa(u),t)\) is a polynomial automorphism, each
horizontal map \(u\mapsto u+t\kappa(u)\) is a polynomial automorphism and

\[
\det(I+tJ\kappa)=1.
\]

Therefore \(J\kappa\) is nilpotent.  As

\[
J\kappa=-\mathcal H''(p)K,
\]

on a local Legendre chart with \(L''=K^{-1}\) we obtain

\[
\boxed{\det(L''-z\mathcal H'')=\det L''}\tag{1.1}
\]

for every \(z\).

## 2. Ternary singular-Hessian normal form

Assume the projective Gauss rank is two, so \(\operatorname{rank}\mathcal H''=2\)
generically.  The associated kernel of \(\mathcal H''\) is a three-variable
quasi-translation.  By the dimension-three quasi-translation classification,
after a constant linear change of the \(p\)-coordinates its projective kernel
has a constant linear invariant.  Integrating the Hessian-kernel equations
leaves exactly two cases.

1. The kernel line is constant.  Then \(\mathcal H\), up to an affine term,
   depends on two constant linear forms, so \((-\nabla\mathcal H,1)\) has a
   constant linear relation.  This is already closed by the previous
   linearly-dependent-kernel theorem.
2. The only genuinely moving case is
   \[
   \boxed{\mathcal H(t,r,s)=a(t)r+b(t)s+c(t).}\tag{2.1}
   \]

We therefore assume (2.1) and derive a contradiction.

On the open set \(a'(t)\ne0\), put

\[
q(t)=\frac{b'(t)}{a'(t)},\qquad h=r+q(t)s.
\]

The kernel direction is proportional to

\[
d=q(t)\partial_r-\partial_s.
\]

The coefficient of \(z^2\) in (1.1) is the rank-one adjugate obstruction and
gives

\[
d^{\mathsf T}L''d=0.
\]

Hence the local Legendre potential is affine along the \(d\)-lines:

\[
L=sP(t,h)+Q(t,h).\tag{2.2}
\]

## 3. The remaining determinant coefficient kills the mixed tail

Substitute (2.1)--(2.2) into the coefficient of \(z\) in (1.1), using
\(b'=a'q\).  Exact differentiation gives

\[
[z]\det(L''-z\mathcal H'')
=P_h\Bigl(
 (a''h+c'')P_h
 -3s a'q'P_h
 -2a'q'Q_h
 -2a'P_t
\Bigr).\tag{3.1}
\]

If the kernel genuinely moves, \(q'\ne0\), and comparison of the coefficient
of \(s\) in (3.1) forces

\[
\boxed{P_h=0.}\tag{3.2}
\]

Thus

\[
L=sP(t)+Q(t,h).\tag{3.3}
\]

The transverse Hessian determinant now factorizes exactly as

\[
\boxed{
\det L''
=-\bigl(P'(t)+q'(t)Q_h\bigr)^2Q_{hh}.
}\tag{3.4}
\]

## 4. Polynomiality closes the moving scroll

The slice coordinates are \(u=\nabla L\).  From (3.3),

\[
u_2=Q_h,\qquad
u_3=P(t)+q(t)u_2.\tag{4.1}
\]

But \(t=p_1(u)=\partial_{u_1}A_0(u)\) is a polynomial in the slice coordinates.
Differentiating (4.1) with respect to \(u_1\), and using the nonvanishing
factor in (3.4), gives

\[
\partial_{u_1}t=0.
\]

Hence

\[
t=T(u_2,u_3)\in K[u_2,u_3].\tag{4.2}
\]

No birationality assumption is needed.  Let \(h_0\in K[u_2,u_3]\) be a
closed/generative polynomial for \(T\):

\[
T=R(h_0),\qquad R\in K[z],
\]

and \(K(h_0)\) is the relative algebraic closure of \(K(T)\) in
\(K(u_2,u_3)\).  Such a generative polynomial exists and is unique up to an
affine change.

The element

\[
P(T)=u_3-q(T)u_2
\]

belongs to that relative algebraic closure, hence is a rational function of
\(h_0\).  Write

\[
P(T)=\widetilde P(h_0),\qquad q(T)=\widetilde q(h_0).
\]

Then

\[
u_3=\widetilde P(h_0)+\widetilde q(h_0)u_2,
\]

so

\[
K(u_2,u_3)=K(u_2,h_0).
\]

Consequently the polynomial \(h_0(u_2,u_3)\) has degree one in \(u_3\):

\[
h_0=\alpha(u_2)u_3+\beta(u_2).\tag{4.3}
\]

Inverting (4.3),

\[
\frac{h_0-\beta(u_2)}{\alpha(u_2)}
=\widetilde P(h_0)+\widetilde q(h_0)u_2.\tag{4.4}
\]

The left side is linear in \(h_0\).  Therefore both
\(\widetilde P\) and \(\widetilde q\) are affine in \(h_0\).  Comparing the
coefficient of \(h_0\) gives

\[
\frac1{\alpha(u_2)}=c+d u_2.
\]

Since \(\alpha\in K[u_2]\), both factors in
\(\alpha(u_2)(c+d u_2)=1\) are units.  Hence \(d=0\), so
\(\widetilde q\) is constant.  Therefore \(q\) itself is constant, contrary
to \(q'\ne0\).

The moving scroll is impossible.  Thus the ternary kernel line is constant,
and the four-component source kernel has a constant linear relation.  The
previous HC4 closure applies.

## 5. Consequence

A final linearly-independent rank-three `[4]` HC4 obstruction cannot have a
polynomial graph gradient image in any affine target direction.  Any remaining
candidate must use genuinely non-graph developable geometry; equivalently,
normalizing one component of the projective normal requires a nonconstant
rational denominator.

This is substantially narrower than the smooth graph-chart reduction of
`HC4RSD71`: the polynomial graph case is now closed globally, including finite
quotient degree.

## 6. External inputs

The only external structural inputs are:

- M. de Bondt, *Quasi-translations and singular Hessians*, for the
  three-dimensional quasi-translation classification;
- I. Arzhantsev and A. Petravchuk, *Closed and Irreducible Polynomials in
  Several Variables*, for existence and uniqueness (up to affine change) of a
  generative/closed polynomial.

The determinant identities (3.1) and (3.4) are checked by the companion
symbolic verifier.
