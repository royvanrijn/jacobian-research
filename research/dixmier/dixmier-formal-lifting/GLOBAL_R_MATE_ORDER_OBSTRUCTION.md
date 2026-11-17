# A global first-order R-mate obstruction

Let K have characteristic zero, L=K((h)), and use the rank-two Weyl
algebra with [P,X]=[Z,Q]=h. Fix R=X(2-3XQ). For the globally admitted
mixed image S_lambda in the
[global admission note](GLOBAL_ADMISSION_AND_QUARTIC_CRITICAL_FACE.md),
there is no globally admitted D of momentum differential order at
most one satisfying

\[
                         [D,R]=h,\qquad [D,S_\lambda]=0.  \tag{1}
\]

This remains true after the actual global first-order perturbation
S_{lambda,nu}=S_lambda+nu*h*X^2*C described below. The statement
allows arbitrary coefficients in L, including poles in h. It concerns
first-order D in the original momentum variables P,Z. It makes no
nonexistence claim about D of higher momentum order.

## The global normal form forces polynomial r-dependence

[QEP1](EXACT_POSITION_CENTRALIZER_AND_COMPLETION.md) gives

\[
 \operatorname{Cent}(R)=L[X,Q][B;h\Delta],\qquad
 B=3X^2P+2(1-3XQ)Z,
\]

and the global frame

\[
 D_0=(1+3XQ)P/2-3Q^2Z,\qquad [D_0,R]=h.
\]

Thus any global first-order solution of [D,R]=h has exactly the form
D=D_0+f_1(X,Q)B+f_0(X,Q), with f_i in L[X,Q]. The frame change
preserves momentum order, since its coefficient matrix has unit
determinant.

On localizing only at X, QEP1's exact finite chart is

\[
 A_2(L)[X^{-1}]
 =L[r,v,v^{-1}][a;h\delta][d;h\partial_r],\qquad
 r=R,\quad v=X^{-1},\quad \delta=v^{-1}\partial_v.
                                                                  \tag{2}
\]

Here [d,r]=h, [a,v]=h/v and all other chart-generator brackets vanish.
Both directions of this first-order frame change have coefficients
Laurent in v and polynomial in r; no r-denominator is introduced.
Consequently a global first-order D satisfying [D,r]=h must become

\[
                   D=d+f(r,v)a+q(r,v),\qquad
                          f,q\in L[r,v,v^{-1}].          \tag{3}
\]

Its coefficient of d is exactly one by [D,r]=h. Polynomial dependence
on r in (3), rather than an arbitrary rational coefficient, is the
essential global admission condition used below.

## A divisor obstruction in the highest momentum coefficient

Let any finite S commuting with r have positive a-order m and highest
coefficient p_m(r,v). The coefficient of a^m in [D,S]/h, for (3), is

\[
            \partial_r p_m+f\delta(p_m)-m p_m\delta(f).   \tag{4}
\]

Multiplication by q contributes only at orders below m. Suppose p_m
has positive r-divisor order N: p_m=r^N u with r not dividing u.
The first term in (4) has exact r-divisor order N-1, whereas the two
vertical terms are divisible by r^N. Their sum cannot vanish in
L[r,v,v^{-1}]. Thus no global first-order D can satisfy (1).
This general local-divisor observation needs no mate for S and no
bound on S's lower coefficients.

For S_lambda, direct expansion of its exact centered square gives

\[
                     p_4=-3\lambda^2h^2r^5v^{-6}.        \tag{5}
\]

Therefore N=5 and the obstruction applies. Equivalently, over the
larger rational coefficient field, (4) forces

\[
 \frac5r-\frac{6f}{v^2}-\frac4v\partial_v f=0,
               \qquad f=\frac{5v^2}{14r}.               \tag{6}
\]

The homogeneous branch is c(r)*v^{-3/2}, which is not rational in v
unless c=0. The unique rational solution displayed in (6) is not
polynomial in r and so cannot come from a global D. It is only a
necessary leading coefficient over that larger field; no localized
first-order solution D is asserted.

## A different admitted perturbation keeps this R-mate obstruction

Use the exact global operators

\[
 \begin{split}
 A_0&=Q+3XQ^2-XB/3=a-v,\qquad C=RA_0-2/3,\\
 W&=X^3C^2,\\
 S_b&=-3RA_0^2+(18XQ-8)A_0+9Q-3hRX,\\
 S_\lambda&=S_b-3\lambda h(RE+CW+WC)-3R\lambda^2h^2W^2,\\
 E&=2X^2C^2+2RhX^4C-R^2h^2X^6.
 \end{split}
\]

Every factor is globally admitted and commutes with R. For h-independent
lambda!=0 and nu!=0 define

\[
                    S_{\lambda,\nu}=S_\lambda+\nu hX^2C.
                                                                  \tag{7}
\]

The addition is global, finite, first order in momentum, and vanishes
classically. Hence it preserves Long's second classical image and
leaves the generic order-four coefficient (5) unchanged. In particular
(7) has no global first-order R-mate by (4).

For the critical face, write z=a-2/(3r), so C=r(z-v). Under
h=t^5, v=t^2V and z=tA, the extra term divided by t^2 becomes
nu*r*V^{-2}(A-tV). The critical reduction is therefore

\[
 P_{\beta,\sigma}
       =g\bigl(A+\beta V^{-3}A^2\bigr)^2+\mu V
                                             +\sigma V^{-2}A,
 \quad [A,V]=\eta/V,                                    \tag{8}
\]

with g=-3r, beta=lambda*r^2, mu=2, sigma=nu*r and eta=1.
The following checks describe this face without constructing a mate.

For weights wt(V)=1 and wt(A)=3, the added first-order term has weight
one. The top symbol remains g*(xi+beta*V^{-3}xi^2)^2, with a nonzero
cubic term of weight six. Thus it remains outside QOI1's strict
quartic cone. For the canonical weight wt(V)=1, wt(A)=-2, the new
term has weight minus four and does not change the bottom weight
minus fourteen. With H=V^2A and C_2(H)=H(H-2eta), the bottom
polynomial remains g*beta^2*C_2(H-7eta)*C_2(H). Its translation-divisor
character is (1+t^2)(1+t^7), zero at a primitive fourteenth root.

The extra term also breaks the q-fixing twisted adjoint used for
the unperturbed mixed square, where q=A+beta*V^{-3}A^2.
For D=eta*V^{-1}partial_V, its unique coefficient-fixing adjoint
q^dagger=q has

\[
 A^\dagger=-A+u,\qquad
 u=\frac{D(\beta V^{-3})-1}{\beta V^{-3}}
                         =-3\eta V^{-2}-\beta^{-1}V^3.
\]

The adjoint of sigma*V^{-2}A has A-coefficient -sigma*V^{-2}, rather
than +sigma*V^{-2}. Thus (8) is not fixed by this adjoint when sigma!=0.
This statement only identifies the failure of that particular
selfadjoint-square argument; it asserts neither a mate nor the absence
of every possible obstruction.

## Evidence and limitations

[The compact exact controls](check_global_r_mate_order.py) check the
highest-coefficient equation, its unique rational solution, the global
perturbation's admission, the unchanged quartic coefficient, the
critical reduction, both weights and bottom divisor character, and the
failure of the q-fixing adjoint. The all-polynomial-r exclusion follows
from the written divisor argument, not a coefficient-degree search.
No higher-order R-mate, finite S-mate, bounded Long lift, formal
verification, independent whole-proof replay or external review is
claimed.
