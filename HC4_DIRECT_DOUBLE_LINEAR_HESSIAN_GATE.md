# Repeated-linear Hessian-factor gates for direct `HC4`

## Status

This note continues `HC4-DIR1--2` in the direct homogeneous filtration of an
arbitrary four-variable constant-Hessian polynomial.  It treats the first
non-squarefree Hessian-discriminant strata, closes exact linear multiplicities
three and four, and closes both boundary ranks of exact linear multiplicity
five.  At exact multiplicity six it closes the complete generic-corank-one
boundary and reduces the lower-rank boundary to three explicit packets.

The identifiers `HC4-DIR1--27` have not yet been promoted to
`MATH_STATUS.json`.  The proof and exact identity checker below form a proof
package awaiting registry audit; they do not change authoritative repository
status.

Let `K` be a characteristic-zero field and let

\[
\Psi\in K[x_1,x_2,x_3,w],\qquad
\det\operatorname{Hess}\Psi=\delta\in K^\times
\]

have degree `D>2`.  Suppose its top homogeneous form has generic Hessian rank
three.  After a constant linear change, write

\[
\Psi_D=f(x_1,x_2,x_3),\qquad
A_0=\operatorname{Hess}f,
\qquad \Delta=\det A_0\ne0.
\]

> **Theorem HC4-DIR3 — double-linear Hessian-factor gate.**  Suppose
> 
> \[
> \Delta=\ell^2R,
> \tag{0.1}
> \]
> 
> where `ell` is linear, `R` is squarefree, and `gcd(ell,R)=1`.  If `D>=5`,
> no such constant-Hessian polynomial `Psi` exists.
>
> For `D=4`, the only possible first transverse motion forces the split top
> form
> 
> \[
> f=c\ell^4+h(y,z),
> \]
> 
> after constant coordinates.  This lies in the already closed
> degree-at-most-four packet.

> **Corollary HC4-DIR3a — triple-linear gate.**  If instead
>
> \[
> \Delta=\ell^3R,
> \qquad R\text{ squarefree},
> \qquad \gcd(\ell,R)=1,
> \]
>
> then no such `Psi` exists for `D>=6`.  For `D=5`, the only possible first
> transverse motion forces
>
> \[
> f=c\ell^5+h(y,z).
> \]
>
> Thus the triple-linear degree-five boundary is an additive split top, not a
> genuinely ternary top cone.

> **Theorem HC4-DIR3b — exact triple-linear closure.**  The degree-five split
> packet in `HC4-DIR3a` is empty.  Consequently exact linear multiplicity
> three admits no four-variable constant-Hessian completion for `D>=5`.

> **Theorem HC4-DIR4 — lower-rank quadruple-linear gate.**  Suppose
>
> \[
> \Delta=\ell^4R,
> \qquad R\text{ squarefree},
> \qquad \gcd(\ell,R)=1,
> \]
>
> and suppose that `Hess(f)` has generic rank at most one on `ell=0`.  Then no
> constant-Hessian completion exists for `D>=5`.
>
> Thus an exact quadruple linear component can survive only on the generic
> rank-two boundary, where its first motion has order one and satisfies the
> linear-field system (5.8) below.

> **Corollary HC4-DIR4a — order-two quadruple gate.**  Under the same exact
> quadruple-factor hypothesis, a first transverse motion of order two is
> impossible for every `D>=5`, regardless of the generic boundary rank.
> The sole normal-form candidate is a degree-six split top, and the quadratic
> suspension descent in Section 5.3 contradicts constant Hessian determinant.

> **Theorem HC4-DIR5 — rank-two order-one collapse.**  In the remaining
> generic rank-two, order-one system, the linear-field matrix has rank one.
> Moreover
>
> \[
> D=6,
> \qquad L=\ell B,
> \qquad B(\ell)\ne0,
> \qquad f=C\ell^6+h_6(y,z).
> \]
>
> Thus the entire exact quadruple-linear frontier is reduced to one
> degree-six additive-top packet with order-one transverse motion.

> **Theorem HC4-DIR6 — exact quadruple-linear closure.**  The degree-six
> packet in `HC4-DIR5` is empty.  Consequently, if
>
> \[
> \det\operatorname{Hess}f=\ell^4R,
> \qquad R\text{ squarefree},
> \qquad \gcd(\ell,R)=1,
> \]
>
> then no four-variable constant-Hessian completion exists for `D>=5`.

> **Theorem HC4-DIR7 — corank-one quintuple-linear gate.**  Suppose
>
> \[
> \det\operatorname{Hess}f=\ell^5R,
> \qquad R\text{ squarefree},
> \qquad \gcd(\ell,R)=1,
> \]
>
> and `Hess(f) mod ell` has generic rank two.  Then no four-variable
> constant-Hessian completion exists for `D>=5`.

> **Theorem HC4-DIR8 — lower-rank quintuple-linear gate.**  Under the exact
> quintuple-factor hypothesis of `HC4-DIR7`, suppose instead that
> `Hess(f) mod ell` has generic rank at most one.  Then no four-variable
> constant-Hessian completion exists for `D>=5`.
>
> Consequently exact linear multiplicity five is empty in every boundary
> rank.

> **Theorem HC4-DIR9 — generic sextuple high-order gate.**  Suppose
>
> \[
> \det\operatorname{Hess}f=\ell^6R,
> \qquad R\text{ squarefree},
> \qquad \gcd(\ell,R)=1,
> \]
>
> and `Hess(f) mod ell` has generic rank two.  Then the first transverse
> order satisfies `j<=3`, and the cases `j=3` and `j=2` are impossible.
> Consequently every surviving generic-corank-one sextuple packet has `j=1`
> and the rank-at-most-two quadratic-vector system (5.64a)--(5.64c) below.

> **Theorem HC4-DIR10 — lower-rank sextuple reduction.**  Under the exact
> sextuple-factor hypothesis, generic boundary rank zero is impossible and
> rank one forces `j<=2`.  At `j=2`, the only top packet is the degree-five
> tangent form (5.67).  At `j=1`, the top lies in one of the two explicit
> order-six boundary jets (5.68)--(5.69).

> **Theorem HC4-DIR11 — sextuple quadratic-Jacobian rank-zero gate.**  In the
> surviving generic-corank-one system of `HC4-DIR9`, put
> `N=Jac(Q)`.  The case `rank(N mod ell)=0` is impossible.  Hence its exact
> remaining boundary ranks are one and two.

> **Theorem HC4-DIR12 — sextuple quadratic-Jacobian rank-one reduction.**
> In the rank-one case left by `HC4-DIR11`, the axial and normal image
> orientations are impossible.  After coordinates preserving `ell=x`, every
> survivor has tangent image and can be written
>
> \[
> Q=q\,\partial_z+x^2v,
> \qquad q_0=q|_{x=0}\ne0,
> \tag{0.2}
> \]
>
> with `v` constant.  If
>
> \[
> t=\operatorname{ord}_x(\partial_zf),
> \]
>
> then `3<=t<=6`.  Thus boundary rank one is reduced to four finite tangent
> packets.  For `t=3,4,5`, one may further normalize `q_0=y^2`; for `t=6`,
> `q_0` is an arbitrary nonzero binary quadratic, with `m` odd whenever
> `q_0` is not a square.  The precise leading equations are
> (5.64k)--(5.64m) below.

> **Theorem HC4-DIR13 — rank-one functional-composition split.**  Each of
> the four tangent packets in `HC4-DIR12` has exactly two types.  In the
> composite type, `q`, `F`, and `D_u f` are binary in `x,l` for a linear
> `l`.  For `t<6`, `D_u(l)=0` and this becomes
> `f=zG(x,y)+h(x,y)`; at `t=6` there is one additional transverse composite
> orientation.  In the primitive type,
>
> \[
> F=\sum_k c_k x^{m-2k}q^k,
> \tag{0.3}
> \]
>
> and its largest index is `k=(m+3-t)/2`.  Thus the primitive packets have
> parity `(t,m)=(3,even),(4,odd),(5,even),(6,odd)`.  For `t<6`, their
> quadratic is uniquely normalized to `q=y^2+xz`.

> **Theorem HC4-DIR14 — invariant composite sextuple gate.**  In the
> composite orientation `D_u(l)=0` of `HC4-DIR13`, exact sextuple order
> eliminates `t=3,5,6`.  The sole survivor has
>
> \[
> t=4,
> \qquad
> f=zG(x,y)+h(x,y),
> \qquad
> \operatorname{ord}_xG=4,
> \qquad
> h|_{x=0}=A y^{m+2},\ A\ne0.
> \tag{0.4}
> \]
>
> Its exact order-six coefficient is given in (5.64v).  Thus the composite
> frontier is one invariant `t=4` pencil plus the transverse `t=6`
> orientation.

> **Theorem HC4-DIR15 — invariant composite closure.**  The `t=4` pencil in
> `HC4-DIR14` is impossible.  Consequently the entire invariant composite
> orientation `D_u(l)=0` is empty.  The only composite rank-one survivor is
> the transverse orientation `D_u(l)!=0` at `t=6`.

> **Theorem HC4-DIR16 — transverse composite closure.**  The transverse
> composite `t=6` orientation is impossible.  Hence every composite case in
> `HC4-DIR13` is empty, and the rank-one quadratic-Jacobian frontier consists
> only of the four primitive conic rows (5.64r).

> **Theorem HC4-DIR17 — sextuple quadratic-Jacobian rank-one closure.**  The
> four primitive conic rows are impossible.  Consequently
>
> \[
> \operatorname{rank}(\operatorname{Jac}(Q)\bmod x)=1
> \]
>
> never occurs in the surviving generic-corank-one sextuple system.  Together
> with `HC4-DIR11`, this leaves boundary Jacobian rank exactly two.

> **Theorem HC4-DIR18 — rank-two boundary root-count gate.**  In that sole
> rank-two system, either
>
> \[
> f|_{x=0}=0,
> \qquad
> Q_x\in(x^2),
> \tag{0.5}
> \]
>
> or the binary form `f(0,y,z)` has at most three distinct projective roots.
> After algebraic closure the nonzero boundary is therefore one of
>
> \[
> y^D,
> \qquad
> y^a z^{D-a},
> \qquad
> y^a z^b(y-z)^c,
> \quad a+b+c=D,
> \tag{0.6}
> \]
>
> with positive displayed exponents.  This root-count reduction is valid for
> any homogeneous order-one system whose boundary matrix is linear of rank
> two; it is not specific to sextuple multiplicity.

> **Theorem HC4-DIR19 — rank-two one-root jet gate.**  In the one-root
> profile of `HC4-DIR18`, boundary-preserving coordinates force
>
> \[
> \begin{aligned}
> f&=y^{m+2}+C x^2y^m+O(x^3),\\
> Q_x&=\lambda xy+\alpha x^2,\\
> Q_y&=-\frac{2C\lambda}{m+2}x^2,\\
> Q_z&=q_2(y,z)+x r_1(y,z)+\gamma x^2,
> \end{aligned}
> \tag{0.7}
> \]
>
> with `C*lambda!=0`.  Thus its arbitrary degree-`m` second boundary jet
> collapses to the single monomial `C*y^m`.

> **Theorem HC4-DIR20 — rank-two one-root closure.**  The normal form (0.7)
> is incompatible with the `dx` component of the order-one matrix equation.
> Hence the one-root rank-two profile is empty.

> **Theorem HC4-DIR21 — remaining root-profile first-jet gate.**  For the
> two-root boundary `f_0=y^a z^b`, a boundary-preserving shear reduces
> `f_1=[x]f` to exactly one of
>
> \[
> 0,
> \qquad
> \kappa y^{a-2}z^{b+1}\ (a\ge2),
> \qquad
> \kappa y^{a+1}z^{b-2}\ (b\ge2).
> \tag{0.8}
> \]
>
> For the three-root boundary `f_0=y^a z^b(y-z)^c`, it reduces to
>
> \[
> f_1=\kappa y^a z^b(y-z)^{c-1}.
> \tag{0.9}
> \]
>
> Thus no arbitrary degree-`D-1` first boundary jet remains in either root
> profile.

> **Theorem HC4-DIR22 — remaining root-profile Schur gate.**  In the
> two-root zero-first-jet packet and in every three-root packet, one has
>
> \[
> f=f_0+x^8h_{D-8}+O(x^9),
> \qquad D\ge8.
> \tag{0.10}
> \]
>
> In particular the three-root parameter `kappa` in (0.9) vanishes.  The two
> exceptional outer jets in (0.8) can survive only with `a>=4` or `b>=4`,
> respectively, and their second jets are the explicit monomials
> (5.64au)--(5.64av) below.

> **Theorem HC4-DIR23 — normal rank-two closure.**  The normal packet
> (0.5) is empty.  More precisely, its scalar and matrix equations force,
> after boundary-preserving coordinates,
>
> \[
> f=xCy^{m+1}+x^2\frac{\alpha C}{\rho}y^{m-1}z+O(x^3),
> \qquad
> Q^x=\alpha x^2,
> \qquad
> Q_{\rm tan}\bmod x=\rho y^2\partial_z,
> \tag{0.11}
> \]
>
> with `alpha*C*rho!=0`.  The boundary value of `Q(F) in (x^3)` then
> kills the remaining diagonal tangent weight, while the next tangent
> matrix coefficient is a nonzero immutable multiple of `m+3`.  Hence no
> normal rank-two packet survives.

> **Theorem HC4-DIR24 — delayed-jet rank-collapse gate.**  In any rank-two
> quadratic-field packet with invertible binary boundary Hessian, the
> simultaneous conditions
>
> \[
> f_1=f_2=0
> \tag{0.12}
> \]
>
> force `rank(Jac(Q) mod x)<=1`.  Indeed the first tangent field `l`
> satisfies both `l(f_0)=0` and
> `Jac(l)^T*grad(f_0)=0`, and differentiation gives
> `Hess(f_0)l=0`.  Thus every order-eight two- and three-root packet in
> (0.10) is empty.  This gate depends only on the two missing jets and the
> invertible boundary Hessian, not on root multiplicities or exact sextuple
> order.

> **Theorem HC4-DIR25 — outer-jet closure.**  Both exceptional two-root
> outer jets in (0.8) are empty.  For the first orientation, polynomiality
> of the boundary Hessian kernel fixes
>
> \[
> Q_0=\left(y^2,-\frac{2\kappa}{a}yz,
>                 \frac{\kappa}{b}z^2\right)
> \tag{0.13}
> \]
>
> after rescaling.  The coefficient of `x` in `Q(f)=x^3F` fixes the four
> possible entries of `Q_1`; the two tangent components of the matrix
> equation at the same order then contain the unavoidable coefficient
> `-(2a+4b)*kappa^2`.  The mirrored orientation is symmetric.  Consequently
> the complete generic-corank-one exact-sextuple system is empty.

> **Theorem HC4-DIR26 — lower-rank sextuple synchronization.**  The
> `x^2g` order-one packet in (5.68) is empty.  The pure-power order-one
> packet (5.69) either has zero exact-sextuple coefficient or is the unique
> resonance `m=3`, where
>
> \[
> f=y^5+cx^3y^2+x^4(uy+vz)+dx^5,
> \qquad
> L=x\partial_z,
> \qquad cv\ne0,
> \tag{0.14}
> \]
>
> and
>
> \[
> \det\operatorname{Hess}f
> =-32v^2x^6(cx^3+10y^3).
> \tag{0.15}
> \]
>
> In particular this is already a subfamily of the degree-five tangent top
> (5.67).  Thus all lower-rank exact-sextuple survivors have the single top
> geometry `f=C*z*x^4+h_5(x,y)`, with either the order-two constant direction
> `partial_z` or the resonant order-one direction `x*partial_z`.

> **Theorem HC4-DIR27 — lower-rank order-two scalar-pivot closure.**  No
> order-two completion of the common top family (5.67) can be an HC4
> counterexample.  Every such completion is affine in the top-kernel
> variable and has quadratic pivot
>
> \[
> P=2Cx^2+\ell(x,y,z).
> \tag{0.16}
> \]
>
> If the tangent part of `ell` vanishes, its bordered Hessian determinant is
> a nonconstant square times a polynomial and cannot be a unit.  Otherwise
> `P` has a constant unit direction, and the registered quadratic
> zero-corner scalar-parent theorem `HC4RSD12` reduces every collision fiber
> to `HC3`.  Thus the sole exact-sextuple lower-rank counterexample frontier
> is the order-one resonance (0.14).

Consequently a genuinely ternary rank-three top cone of degree at least five
can survive `HC4-DIR2--27` only if its Hessian determinant has a repeated
nonlinear factor, a linear factor of multiplicity at least six, or at least
two distinct repeated factors.  Exact linear multiplicities three, four, and
five are completely closed by `HC4-DIR3b--8`; at exact multiplicity six the
generic rank-one quadratic-Jacobian branch is empty.  The only generic
sextuple order-one branch has quadratic-Jacobian boundary rank two and is
reduced by `HC4-DIR18` to the normal packet (0.5) and the three root profiles
(0.6); `HC4-DIR20` closes the one-root profile, `HC4-DIR23` closes the normal
packet, `HC4-DIR24` closes every delayed two-/three-root row, and `HC4-DIR25`
closes the two outer jets.  Thus the complete generic-corank-one exact
sextuple stratum is empty.  `HC4-DIR26` further synchronizes the three
lower-rank Hessian handoffs of `HC4-DIR10` to the single degree-five tangent
top family (5.67), with order two or one along its same tangent direction;
`HC4-DIR27` closes the complete order-two counterexample branch.  Only the
degree-five order-one resonance (0.14) remains at exact linear multiplicity
six.

## 1. The first moving off-diagonal block

Put

\[
m=D-2
\]

and use the scaled Hessian polynomial

\[
M(t)=H_D+tH_{D-1}+\cdots+t^mH_2,
\qquad
\det M(t)=\delta t^{4m}.
\tag{1.1}
\]

Split off the constant top-kernel direction:

\[
M(t)=
\begin{pmatrix}
A(t)&b(t)\\
b(t)^{\mathsf T}&c(t)
\end{pmatrix},
\qquad A(0)=A_0,\qquad b(0)=c(0)=0.
\tag{1.2}
\]

If `b(t)=0`, then

\[
\det M(t)=\det A(t)c(t)
\]

has `t`-adic order at most `m`, contradicting (1.1).  Hence let `j>=1` be
minimal with

\[
b(t)=t^jb_j+O(t^{j+1}),\qquad b_j\ne0.
\tag{1.3}
\]

Every coefficient of `c(t)` below order `2j` vanishes.  At order `2j`, the
block determinant gives

\[
c_{2j}\Delta=b_j^{\mathsf T}\operatorname{adj}(A_0)b_j,
\tag{1.4}
\]

where `c_(2j)=0` if `2j>m`.  In particular, the right side is divisible by
`Delta`.

## 2. The radical divisibility gate

Let `pi` be a simple irreducible factor of `R`.  Modulo `pi`, the matrix
`A0` has generic rank two and

\[
\operatorname{adj}(A_0)=\rho vv^{\mathsf T}.
\]

Equation (1.4) gives `(v^T b_j)^2=0`, hence `pi` divides every component of
`adj(A0)b_j`.

The same conclusion holds for `ell`.  If `A0 mod ell` has rank two, reduce
(1.4) modulo `ell` as above.  If it has rank at most one, every two-by-two
minor, and therefore every entry of `adj(A0)`, is already divisible by `ell`.
Since the factors are coprime,

\[
\ell R\mid\operatorname{adj}(A_0)b_j.
\tag{2.1}
\]

All entries of `A0` have degree `m`, while `b_j` has degree `m-j`.  Therefore

\[
\deg\bigl(\operatorname{adj}(A_0)b_j\bigr)=3m-j,
\qquad
\deg(\ell R)=3m-1.
\tag{2.2}
\]

If `j>1`, (2.1)--(2.2) force `adj(A0)b_j=0`; multiplying by `A0` then gives
`Delta*b_j=0`, a contradiction.  Thus

\[
\boxed{j=1.}
\tag{2.3}
\]

Equality of the two degrees in (2.2) now gives a nonzero constant vector `B`
such that

\[
\operatorname{adj}(A_0)b_1=\ell R B.
\tag{2.4}
\]

Multiplying by `A0` and using (0.1) yields

\[
\boxed{A_0B=\ell b_1.}
\tag{2.5}
\]

## 3. Curl-freeness forces a pure power

The first lower homogeneous layer is affine in `w`, so

\[
\Psi_{D-1}=w a(x_1,x_2,x_3)+h(x_1,x_2,x_3),
\qquad b_1=\nabla a,
\]

with `deg a=m`.  Because `B` is constant,

\[
A_0B=\nabla(D_Bf).
\]

Equation (2.5) becomes

\[
\nabla(D_Bf)=\ell\nabla a.
\tag{3.1}
\]

Taking exterior derivative gives

\[
d\ell\wedge da=0.
\]

Choose `x=ell`.  Then `a_y=a_z=0`, and homogeneity gives

\[
a=cx^m,\qquad c\ne0.
\tag{3.2}
\]

Integrating (3.1),

\[
D_Bf=c' x^{m+1},\qquad c'\ne0.
\tag{3.3}
\]

## 4. The two normal forms contradict multiplicity two

There are two constant-linear cases.

If `B(x)!=0`, choose two linear coordinates `y,z` killed by `B` and rescale
so `B=partial_x`.  Equation (3.3) integrates to

\[
f=Cx^{m+2}+h(y,z).
\tag{4.1}
\]

Its Hessian determinant is

\[
\det\operatorname{Hess}f
=C(m+2)(m+1)x^m\det\operatorname{Hess}_{y,z}h.
\tag{4.2}
\]

Since the binary factor is independent of `x`, the exact `x`-multiplicity is
`m`.  Assumption (0.1) therefore forces `m=2`, equivalently `D=4`.  This is
the displayed split top form in the theorem.

If `B(x)=0`, choose constant coordinates with `B=partial_z`.  Then

\[
f=Czx^{m+1}+h(x,y),
\tag{4.3}
\]

and direct expansion gives

\[
\det\operatorname{Hess}f
=-C^2(m+1)^2x^{2m}h_{yy}.
\tag{4.4}
\]

Thus the `x`-multiplicity is at least `2m`, which is larger than two for
`D>=4`.

For `D>=5`, both cases contradict (0.1).  This proves `HC4-DIR3`.

For the triple-linear corollary, the same argument improves the divisibility
by one power:

\[
\ell^2R\mid\operatorname{adj}(A_0)b_j.
\tag{4.5}
\]

To see this, localize at `ell` and diagonalize the symmetric matrix over the
DVR.  The valuation partitions of a determinant of order three are

\[
(0,0,3),\qquad(0,1,2),\qquad(1,1,1).
\]

In the first case, integrality of `b_j^T A0^{-1}b_j` forces valuation at least
two in the unique singular component.  In the second, the order-two component
first has valuation at least one and then the order-one component does too.
In the third, every adjugate entry already has valuation at least two.  Thus
(4.5) holds in all cases.

Now `deg(ell^2*R)=3m-1`, so again `j=1` and the quotient is a constant vector.
Because `Delta/(ell^2*R)=ell`, equations (2.5)--(4.4) are unchanged.  The
transverse normal form has exact `ell`-multiplicity `m`, forcing `m=3` and
`D=5`; the tangent normal form has multiplicity at least `2m` and is
impossible for `D>=4`.  This proves `HC4-DIR3a`.

### 4.1 The terminal degree-five triple packet

It remains to exclude the split packet left by `HC4-DIR3a`.  Write its
forced layers as

\[
\begin{aligned}
h_5&=Cx^5+H_5(y,z),\\
h_4&=\alpha wx^3+r_4(x,y,z),\\
h_3&=Bw^2x+w r_2(x,y,z)+\phi_3(x,y,z),
\end{aligned}
\tag{4.6}
\]

where `C*alpha!=0`.  The scalar Schur coefficient gives

\[
2Bx=\frac{(3\alpha x^2)^2}{20Cx^3},
\qquad
\boxed{B=\frac{9\alpha^2}{40C}}.
\tag{4.7}
\]

Let

\[
\mathcal H(\lambda)
=H_2+\lambda H_3+\lambda^2H_4+\lambda^3H_5,
\qquad
R=\det\operatorname{Hess}_{y,z}(H_5)\ne0.
\]

The two-top-passive part of its determinant is controlled by

\[
\Phi
=\lambda^3Cx^5+\lambda^2\alpha wx^3+\lambda Bw^2x,
\]

whose binary Hessian determinant is

\[
(40BC-9\alpha^2)\lambda^4x^4
-4B^2\lambda^2w^2.
\tag{4.8}
\]

The first coefficient is exactly (4.7).  For the second, project the full
four-by-four determinant to passive degree six and `x`-degree zero.  The two
passive rows and columns must use `Hess_(y,z)(H_5)`, leaving binary weight
two.  Its partitions are `1+1` and `2+0`.  The latter cannot retain two
powers of `w` without also retaining a power of `x`.  In (4.6), the arbitrary
linear-`w` term has zero `ww` derivative and contributes no `w` to its `xw`
derivative.  Thus the complete `3+3+1+1` channel is the immutable term

\[
-4B^2\lambda^8w^2R
=-\frac{81\alpha^4}{400C^2}\lambda^8w^2R\ne0.
\tag{4.9}
\]

This contradicts the constant determinant identity, proves `HC4-DIR3b`, and
closes exact triple-linear multiplicity.

## 5. The lower-rank quadruple component

Assume now

\[
\Delta=x^4R,
\qquad R\text{ squarefree},
\qquad x\nmid R,
\tag{5.1}
\]

and that `A0 mod x` has rank at most one over `K(y/z)`.  The condition
`Delta!=0` means that the generic rank is exactly one.

### 5.1 The first motion has order at most two

At `x`, diagonalize `A0` by congruence over the DVR.  Its invariant-factor
valuations add to four.  Generic boundary rank at most one leaves

\[
(0,1,3),\qquad (0,2,2),\qquad (1,1,2).
\tag{5.2}
\]

The last partition would make every entry of `A0` divisible by `x`; Hessian
integrability then makes `f` divisible by `x^3`, forcing
`ord_x(det Hess(f))>=7`, contrary to (5.1).  We retain it in the valuation
estimate because it only strengthens the conclusion.

For `(0,2,2)` and `(1,1,2)`, every entry of the adjugate is divisible by
`x^2`.  For `(0,1,3)`, polynomiality of

\[
b_j^{\mathsf T}A_0^{-1}b_j
\]

forces the component of `b_j` in the order-three direction to be divisible
by `x`; hence `x^2` again divides every component of
`adj(A0)b_j`.  At each simple factor of `R`, the radical argument of Section
2 applies.  Consequently

\[
x^2R\mid\operatorname{adj}(A_0)b_j.
\tag{5.3}
\]

Because `deg(x^2R)=3m-2` and
`deg(adj(A0)b_j)=3m-j`, equation (5.3) gives

\[
j\le2.
\tag{5.4}
\]

If `j=2`, the quotient is a constant vector `B` and

\[
A_0B=x^2\nabla a.
\tag{5.5}
\]

Curl-freeness gives `a=c*x^(m-1)` and
`D_B f=c'*x^(m+1)`.  The two normal forms of Section 4 then have
`x`-multiplicity `m` and at least `2m`.  Exact multiplicity four leaves only
`m=4` in the transverse form

\[
f=Cx^6+h(y,z).
\tag{5.6}
\]

But `Delta!=0` makes the binary Hessian of `h` generically nonsingular, so
`Hess(f) mod x` has rank two.  This contradicts the lower-rank hypothesis.
Thus `j=2` is impossible here.

### 5.2 The two rank-one boundary jets

It remains to exclude `j=1`.  The quotient in (5.3) is now a linear vector

\[
L=M(x,y,z)^{\mathsf T},
\]

and multiplication by `A0` gives

\[
A_0L=x^2\nabla a.
\tag{5.7}
\]

The scalar order-two determinant coefficient is also essential.  It gives

\[
\boxed{
\operatorname{Hess}(f)L=x^2\nabla a,
\qquad
L(a)=x^2c,
}
\tag{5.8}
\]

where `a` is homogeneous of degree `m` and `c` of degree `m-2`.  These are
identities, not a bounded ansatz.

Write `D=m+2` and expand `f` in powers of `x`.  Since the binary Hessian of
`f|_(x=0)` has rank at most one, the zero-Hessian lemma for binary forms and
the remaining two-by-two minors give, after a constant change preserving
`x`, exactly one of

\[
f=x^2g_m(x,y,z),
\tag{5.9}
\]

or

\[
f=y^{m+2}+x^3g_{m-1}(x,y,z).
\tag{5.10}
\]

Here one may first extend `K` to its algebraic closure; nonexistence after
that scalar extension implies nonexistence over `K`.

Indeed, if the boundary value is nonzero, it is a pure power.  Normalize it
to `y^(m+2)`.  Rank one then forces the next two `x`-coefficients to be the
first two Taylor coefficients of `(y+alpha*x)^(m+2)`; shifting `y` removes
them and yields (5.10).  If the boundary value is zero, the arrow-shaped
boundary Hessian forces the coefficient of `x` to vanish, giving (5.9).

For (5.9), put `g_0=g|_(x=0)`.  Exact order four says that the bordered
Hessian

\[
E(g_0)=
\det\begin{pmatrix}
2g_0&2(g_0)_y&2(g_0)_z\\
2(g_0)_y&(g_0)_{yy}&(g_0)_{yz}\\
2(g_0)_z&(g_0)_{yz}&(g_0)_{zz}
\end{pmatrix}
\tag{5.11}
\]

is nonzero.  Contracting the first equation in (5.8) with the Euler vector
gives

\[
a=\frac{m+1}{m}x^{-2}L(f).
\tag{5.12}
\]

Divisibility by `x^2` first forces `L(x)=lambda*x`.  Let `N` be the vector
field induced by `L` on `x=0`.  The first equation of (5.8), reduced at the
boundary, gives

\[
N(g_0)=-\lambda g_0.
\tag{5.13}
\]

The scalar equation `x^2|L(a)` then gives `lambda=0`.  Write
`L(g)=xq`.  The tangential boundary part of (5.8) is

\[
N^{\mathsf T}\nabla g_0=0.
\tag{5.14}
\]

If `N!=0`, a nonzero constant direction annihilates `g_0`, so `g_0` is a
power of one binary linear form; then (5.11) vanishes.  Hence `N=0` and
`L=xB` for a constant tangent vector `B`.  Equation (5.8) now forces

\[
D_Bg=cx^{m-1}.
\]

After taking `B=partial_z`,

\[
f=cx^{m+1}z+x^2h(x,y),
\qquad
\det\operatorname{Hess}f
=-c^2(m+1)^2x^{2m+2}h_{yy}.
\tag{5.15}
\]

This has `x`-order at least `2m+2>4`.

For (5.10), exact order four gives

\[
(g_0)_z\ne0,
\qquad
6g_0(g_0)_{zz}-9(g_0)_z^2\ne0.
\tag{5.16}
\]

The order-zero and order-one parts of (5.8) force

\[
L(x)=\mu x,
\qquad L(y)=0,
\qquad L(z)=px+ry+sz.
\tag{5.17}
\]

Put `h=3*mu*g+L(g)`.  The boundary part of (5.8) is

\[
3m\mu g_0=(2m-1)h_0,
\qquad
\nabla h_0=-m(r,s)(g_0)_z,
\qquad
(\mu+(ry+sz)\partial_z)h_0=0.
\tag{5.18}
\]

If `mu!=0`, the first two identities make `g_0` a power of a linear form
`l` with `l_z!=0`.  They give

\[
s=-\frac{3\mu}{2m-1},
\]

while the last identity gives

\[
s=-\frac{\mu}{m-1}.
\]

These are compatible only for `m=2`, whereas `D>=5` means `m>=3`.
Therefore `mu=0`, and (5.18) gives `r=s=0`.  Now `L=px*partial_z`.
Comparing the leading `dx` coefficient in (5.8) gives

\[
3=2\frac{m+1}{m},
\]

again forcing `m=2`, contrary to `m>=3`.  This excludes (5.10), proves
`HC4-DIR4`, and completes the lower-rank quadruple component.

### 5.3 The order-two split descends and contradicts itself

When `A0 mod x` has generic rank two, the invariant factors are `(0,0,4)`.
The same estimate (5.3)--(5.4) applies, but the order-two transverse form
(5.6) now has the correct boundary rank, so the boundary-rank argument alone
cannot discard it.  The full polynomial does.

Write its degree-six top and degree-three `w`-coefficient as

\[
P_6=Cx^6+h_6(y,z),
\qquad A_3=\alpha x^3,
\qquad C\alpha\ne0.
\tag{5.19}
\]

All coefficients of the lower-right Hessian block before order four vanish.
Since `D=6`, the order-four coefficient comes from the quadratic layer and is
a constant `q`.  Thus the entire polynomial is, up to irrelevant affine
terms,

\[
\Psi=P(x,y,z)+wA(x,y,z)+\frac q2w^2.
\tag{5.20}
\]

The order-four determinant coefficient is

\[
\det\operatorname{Hess}_{y,z}(h_6)
\bigl(30Cq-9\alpha^2\bigr)x^4.
\]

It must vanish, so

\[
q=\frac{3\alpha^2}{10C}\ne0.
\tag{5.21}
\]

Take the Schur complement of the constant entry `q` and put

\[
Q=P-\frac{A^2}{2q},
\qquad s=w+\frac Aq.
\]

The constant-Hessian identity becomes

\[
\det\operatorname{Hess}\Psi
=q\det\bigl(\operatorname{Hess}Q+s\operatorname{Hess}A\bigr)
=\delta.
\tag{5.22}
\]

Because `(x,y,z,w) -> (x,y,z,s)` is a polynomial automorphism, (5.22) is a
polynomial identity in `s`.  Setting `s=0` makes
`det Hess(Q)=delta/q` constant.  But (5.19)--(5.21) give

\[
Q_6
=\left(C-\frac{\alpha^2}{2q}\right)x^6+h_6(y,z)
=-\frac23Cx^6+h_6(y,z),
\]

whose ternary Hessian determinant is nonzero.  It is the top homogeneous
part of `det Hess(Q)`, a contradiction.  This proves `HC4-DIR4a`.

Hence the exact unresolved quadruple-linear input initially becomes

\[
\boxed{
j=1,\qquad
\operatorname{Hess}(f)L=x^2\nabla a,\qquad
L(a)=x^2c,\qquad
\operatorname{rank}(\operatorname{Hess}(f)\bmod x)=2.
}
\tag{5.23}
\]

This is the next direct `HC4` target.

### 5.4 The order-one field has rank one

We now solve the linear-field part of (5.23).  Put

\[
F=\frac{m}{m+1}a=x^{-2}L(f),
\qquad L=M(x,y,z)^{\mathsf T}.
\]

The Hessian equation, its scalar coefficient, and the identity
`grad(L(f))=Hess(f)L+M^T*grad(f)` are equivalent to

\[
L(f)=x^2F,
\qquad
M^{\mathsf T}\nabla f
=2xF\,dx-\frac1m x^2\nabla F,
\qquad
L(F)\in(x^2).
\tag{5.24}
\]

Set `x=0`.  The middle identity gives

\[
M^{\mathsf T}\nabla f|_{x=0}=0.
\tag{5.25}
\]

If `M` were invertible, then both the boundary value of `f` and its first
`x`-coefficient would vanish.  The boundary Hessian would have rank at most
one, contrary to (5.23).  Hence `rank(M)<=2`.

Suppose `rank(M)=2`, and let `k` span `ker(M^T)`.  Equation (5.25) writes the
boundary gradient as a polynomial multiple of the fixed vector `k`.  There
are two orientations relative to the boundary plane.

If the tangent part of `k` is nonzero, a constant change preserving `x`,
followed by a tangent shift by a multiple of `x`, gives

\[
k=\partial_y,
\qquad
f=y^{m+2}+x^2g_m,
\qquad
L(x)=\mu x,
\qquad L(y)=0,
\qquad L(z)=px+ry+sz.
\tag{5.26}
\]

Boundary rank two forces `g_0=g|_(x=0)` to be nonzero, and rank two of `M`
forces `mu!=0`.  If `N=(ry+sz)*partial_z`, the first two identities in (5.24)
give

\[
N(g_0)=-\mu g_0,
\qquad F_0=\mu g_0.
\]

But the boundary value of the scalar identity is then

\[
N(F_0)=-\mu^2g_0\ne0,
\]

contradicting `L(F) in (x^2)`.

If the tangent part of `k` is zero, normalize `k=partial_x`.  Then

\[
f=xp_{m+1}(y,z)+x^2g_m(x,y,z),
\qquad L(x)=0.
\tag{5.27}
\]

Let `N` be the tangent boundary field induced by `L`.  The boundary Hessian
equation and the order-one part of (5.24) give

\[
N(p)=0,
\qquad N^{\mathsf T}\nabla p=0.
\tag{5.28}
\]

The field `N` is nonzero, since otherwise `M` would have rank at most one.
A nonzero constant column of `N` therefore annihilates `p`, so after tangent
coordinates `p=y^(m+1)`.  Equations (5.27)--(5.28) put the field in the form

\[
L=bx\partial_y+(dx+ry+sz)\partial_z,
\qquad b\ne0.
\tag{5.29}
\]

The first component of (5.24) gives

\[
F_0=\frac b2p_y.
\]

Its `y`- and `z`-components at the next order give

\[
rg_{0,z}=-\frac1m(F_0)_y\ne0,
\qquad
sg_{0,z}=-\frac1m(F_0)_z=0.
\]

Thus `s=0`.  The last component of the full middle identity in (5.24) now
forces `F_z=0`.  Consequently

\[
L(F)=bxF_y.
\]

Divisibility by `x^2` would make `(F_0)_y=0`, contradicting
`F_0=(b/2)p_y`.  This excludes `rank(M)=2`.

It remains that `rank(M)=1`.  Write

\[
L=\lambda u,
\]

where `u` is a nonzero constant vector and `lambda` is linear.  The middle
identity in (5.24), together with
`lambda*D_u(f)=x^2F`, gives in the fraction field

\[
\frac{dF}{F}=2m\frac{dx}{x}-m\frac{d\lambda}{\lambda}.
\tag{5.30}
\]

Therefore

\[
F=C\frac{x^{2m}}{\lambda^m}.
\tag{5.31}
\]

Since `F` is a nonzero polynomial, `lambda` is proportional to `x`.  After
rescaling,

\[
L=xu,
\qquad F=Cx^m,
\qquad D_u(f)=C'x^{m+1}.
\tag{5.32}
\]

If `u(x)=0`, the tangent normal form of Section 4 gives `x`-multiplicity at
least `2m>4`.  If `u(x)!=0`, take `u=partial_x` while keeping two invariant
coordinates.  Then

\[
f=C''x^{m+2}+h(y,z),
\]

whose exact `x`-multiplicity is `m`.  Hence `m=4`, or `D=6`.  This proves
`HC4-DIR5` and sharpens (5.23) to the sole packet

\[
\boxed{
D=6,
\qquad
f=Cx^6+h_6(y,z),
\qquad
L=x\partial_x
}
\tag{5.33}
\]

up to constant coordinates and rescaling.

### 5.5 The terminal degree-six face

It remains to exclude (5.33).  Write the homogeneous layers of the full
degree-six potential, using `w` for the top-kernel variable, as

\[
\begin{aligned}
h_6&=Cx^6+H_6(y,z),\\
h_5&=\alpha wx^4+r_5(x,y,z),\\
h_4&=Bw^2x^2+w r_3(x,y,z)+\phi_4(x,y,z),\\
h_3&=\gamma w^3+w^2L_1(x,y,z)+w g_2(x,y,z)+\phi_3(x,y,z).
\end{aligned}
\tag{5.34}
\]

Here `C*alpha!=0`.  The first order-one Schur coefficient is exactly the
scalar equation already used in (5.8).  Since

\[
(h_6)_{xx}=30Cx^4,
\qquad
\partial_x\partial_w h_5=4\alpha x^3,
\qquad
\partial_w^2h_4=2Bx^2,
\]

it gives

\[
2Bx^2=\frac{(4\alpha x^3)^2}{30Cx^4},
\qquad
\boxed{B=\frac{4\alpha^2}{15C}}.
\tag{5.35}
\]

Now form the homogeneous Hessian pencil

\[
\mathcal H(\lambda)
=H_2+\lambda H_3+\lambda^2H_4+\lambda^3H_5+\lambda^4H_6.
\tag{5.36}
\]

Put

\[
R(y,z)=\det\operatorname{Hess}_{y,z}H_6\ne0.
\]

In the determinant of (5.36), consider first the
`lambda^13*w*x^4*(y,z)-degree-8` channel.  The only nonzero determinant
permutations in this channel use both top passive Hessian entries and the
`(x,w)` entries displayed in (5.34).  The alternative layer partition
`4+3+3+3` cannot contribute: retaining `w` from `w*x^4` uses the `xx`
entry, while its `w`-row entry also uses the already occupied `x` direction.
The channel is therefore

\[
\bigl(180C\gamma-8\alpha B\bigr)wx^4R.
\tag{5.37}
\]

It must vanish, so

\[
\boxed{
\gamma=\frac{2\alpha B}{45C}
=\frac{8\alpha^3}{675C^2}.
}
\tag{5.38}
\]

Next take the `lambda^11*w^3*(y,z)-degree-8` channel.  The possible layer
partitions are

\[
4+4+2+1,
\quad 4+3+3+1,
\quad 4+3+2+2,
\quad 3+3+3+2.
\]

Only the first can retain three powers of `w` while supplying the `w` row
and column: it uses the two top passive Hessian entries, the `xx` derivative
of `B*w^2*x^2`, and the `ww` derivative of `gamma*w^3`.  Every arbitrary
form in (5.34) has too little `w`-degree in the required rows.  Hence this
entire channel is the immutable coefficient

\[
12B\gamma w^3R
=\frac{128}{3375}\frac{\alpha^5}{C^3}w^3R.
\tag{5.39}
\]

It is nonzero in characteristic zero because `C*alpha*R!=0`.  This
contradicts the constant determinant identity, proves `HC4-DIR6`, and closes
the exact quadruple-linear Hessian-factor stratum.

### 5.6 The generic quintuple-linear component

Assume now

\[
\Delta=x^5R,
\qquad R\text{ squarefree},
\qquad x\nmid R,
\qquad \operatorname{rank}(A_0\bmod x)=2.
\tag{5.40}
\]

The DVR invariant factors at `x` are `(0,0,5)`.  Polynomiality of
`b_j^T*A_0^(-1)*b_j` forces order at least three in its singular component.
Together with the simple factors of `R`, this gives

\[
x^3R\mid\operatorname{adj}(A_0)b_j.
\tag{5.41}
\]

Since `deg(x^3R)=3m-2`, degree comparison gives `j<=2`.

If `j=2`, the quotient is constant and

\[
A_0B=x^2\nabla a.
\]

The normal forms of Section 4 force `m=5`, `D=7`, and

\[
f=Cx^7+H_7(y,z).
\tag{5.42}
\]

Write the first moving layer as `h_5=alpha*w*x^4+...`.  The scalar Schur
coefficient puts `B_2*w^2*x` in `h_3`, with

\[
B_2=\frac{4\alpha^2}{21C}.
\]

The two-top-passive-Hessian channel is governed by

\[
\Phi_2
=\lambda^5Cx^7+\lambda^3\alpha wx^4+\lambda B_2w^2x.
\]

Its binary Hessian determinant is

\[
(84B_2C-16\alpha^2)\lambda^6x^6
+8B_2\alpha\lambda^4wx^3
-4B_2^2\lambda^2w^2.
\tag{5.43}
\]

The first coefficient vanishes by the displayed value of `B_2`, while the
next is

\[
\frac{32\alpha^3}{21C}\lambda^4wx^3.
\]

To see that this coefficient is immutable in the full four-by-four
determinant, take its passive-degree-ten component.  Both passive rows and
columns must use the two weight-five entries of
`Hess_(y,z)(H_7)`; replacing either by a lower layer loses passive degree and
cannot retain the prescribed `w*x^3` monomial.  The remaining `(x,w)` minor
has weight four.  Its possible partitions are `3+1` and `2+2`.  The first is
the displayed `h_5`--`h_3` cross term.  At weight two there is no quadratic
`w` term, because the first off-diagonal order is two and the first scalar
Schur term occurs at order four; hence `2+2` has `w`-degree zero.  The
coefficient is therefore nonzero in the full determinant, and `j=2` is
impossible.

If `j=1`, multiplication by `A_0` gives exactly the order-one system (5.23).
The rank analysis of Section 5.4 applies unchanged and again forces (5.42),
now with `L=x*partial_x`.  Write the forced weighted layers as

\[
\Phi_1
=\lambda^5Cx^7
+\lambda^4\alpha wx^5
+\lambda^3B_1w^2x^3
+\lambda^2\gamma w^3x
+\lambda\eta w^3.
\tag{5.44}
\]

The successive two-passive-Hessian channels give

\[
B_1=\frac{25\alpha^2}{84C},
\qquad
\gamma=\frac{5B_1\alpha}{63C}
=\frac{125\alpha^3}{5292C^2},
\qquad
\eta=0.
\tag{5.45}
\]

The intervening coefficient `6*(-4*B_1^2+15*alpha*gamma)` vanishes
identically after (5.45).  The terminal channel is

\[
-9\gamma^2w^4R,
\tag{5.46}
\]

which is nonzero.  Here too the passive-degree-ten projection forces both
weight-five passive Hessian entries.  The remaining binary weight is four.
For the `w^4` monomial, the only surviving partition is `2+2`, using twice
the `xw` derivative of `gamma*w^3*x`; it contributes `-9*gamma^2*w^4`.
The alternatives `3+1` and `4+0` have `w`-degree at most three after the two
required derivatives.  Any determinant permutation using a passive row or
column outside `Hess_(y,z)(H_7)` has passive degree below ten in this
`x`-degree-zero channel.  Thus arbitrary remainders cannot enter (5.46).
This excludes `j=1` and proves `HC4-DIR7`.

### 5.7 The lower-rank quintuple component

It remains to remove the complementary boundary rank in exact multiplicity
five.  If `rank(A_0 mod x)=0`, Hessian integrability makes `f` divisible by
`x^3`, and then `ord_x(det Hess(f))>=7`.  Thus exact order five has boundary
rank one.  Its DVR invariant factors are

\[
(0,1,4)\qquad\text{or}\qquad(0,2,3).
\tag{5.47}
\]

In the first case polynomiality of `b_j^T*A_0^(-1)*b_j` puts two powers of
`x` in the order-four component of `b_j`; in the second it puts one and two
powers in the order-two and order-three components.  Inspection of the
adjugate valuations therefore gives, in both cases,

\[
x^3R\mid\operatorname{adj}(A_0)b_j.
\tag{5.48}
\]

Again `j<=2`.  At `j=2`, the constant quotient satisfies
`A_0B=x^2*grad(a)`.  The Section 4 normal forms have `x`-multiplicity `m` and
at least `2m`.  Exact multiplicity five leaves only `m=5` in the transverse
form `f=C*x^7+H_7(y,z)`, whose boundary Hessian has rank two.  This
contradicts the present hypothesis.

For `j=1`, multiplication by `A_0` gives the same system (5.8).  The rank-one
boundary classification (5.9)--(5.10) applies, but exact order five asks for
the coefficient immediately after the one used in Section 5.2.

First let

\[
f=x^2g_m,
\qquad
g=g_0+xg_1+O(x^2).
\tag{5.49}
\]

The order-four coefficient is `E(g_0)` from (5.11).  Euler's identities for
the binary degree-`m` form `g_0` give

\[
E(g_0)=-\frac{2(m+1)}{m-1}
g_0\det\operatorname{Hess}_{y,z}(g_0).
\tag{5.50}
\]

If `g_0=0`, the determinant starts in order at least seven.  Otherwise the
vanishing required before exact order five makes `g_0` a pure power; normalize
`g_0=y^m`.  The next determinant coefficient is

\[
[x^5]\det\operatorname{Hess}f
=-2m(m+1)y^{2m-2}(g_1)_{zz}.
\tag{5.51}
\]

Hence exact order five requires `(g_1)_{zz}!=0`.

The boundary equations used in (5.12)--(5.14) first force `L(x)=0`.  With
`N` the boundary field and `g_0=y^m`, they put

\[
L=bx\partial_y+(dx+ry+sz)\partial_z.
\tag{5.52}
\]

Writing `F=x^(-2)L(f)=L(g)`, the `x`-component of (5.24) in order two and
the tangent components in order three give

\[
F_1=\frac{m^2b}{2m-1}y^{m-1},
\qquad
r(g_1)_z=-\frac{m(m-1)b}{2m-1}y^{m-2},
\qquad
s(g_1)_z=0.
\tag{5.53}
\]

If `b!=0`, these identities make `(g_1)_z` independent of `z`, contrary to
(5.51).  If `b=0`, exactness gives `r=s=0`; then
`L=dx*partial_z` with `d!=0`, and the leading `dx` coefficient of (5.24) is

\[
1=2-\frac2m.
\]

It forces `m=2`, outside `D>=5`.  Thus (5.49) is impossible.

Now let

\[
f=y^{m+2}+x^3g_{m-1},
\qquad
g=g_0+xg_1+O(x^2).
\tag{5.54}
\]

The order-four coefficient is the expression in (5.16).  Over `K(y)`, the
identity `2*g_0*(g_0)_{zz}-3*(g_0)_z^2=0` has no nonconstant polynomial
solution in `z`: at a root of multiplicity `r>0`, its lowest coefficient is
`-r*(r+2)`.  Hence its vanishing forces `(g_0)_z=0`.  Exact order five then
requires `g_0!=0`, so homogeneity gives `g_0=c*y^(m-1)`, and

\[
[x^5]\det\operatorname{Hess}f
=6(m+2)(m+1)y^m g_0(g_1)_{zz}.
\tag{5.55}
\]

In particular `(g_1)_{zz}!=0`.  Equations (5.18) now give `h_0=mu=0` and

\[
L=(px+ry+sz)\partial_z.
\]

The next `dx` coefficient of (5.24) first forces
`(ry+sz)*(g_1)_z=0`.  Exactness makes `(g_1)_z` nonzero, so `r=s=0`.
For `L=px*partial_z`, the following coefficient is

\[
1=2-\frac3m,
\]

and hence `m=3`.  But `g_1` then has degree `m-2=1`, so
`(g_1)_{zz}=0`, the final contradiction.  This proves `HC4-DIR8` and closes
the exact quintuple-linear stratum.

### 5.8 The generic sextuple component above first order

Assume

\[
\Delta=x^6R,
\qquad R\text{ squarefree},
\qquad x\nmid R,
\qquad \operatorname{rank}(A_0\bmod x)=2.
\tag{5.56}
\]

The invariant factors are `(0,0,6)`.  Polynomiality of
`b_j^T*A_0^(-1)*b_j` puts three powers of `x` in its singular component, and
the simple factors of `R` contribute their radical.  Therefore

\[
x^3R\mid\operatorname{adj}(A_0)b_j.
\tag{5.57}
\]

The quotient has degree `3-j`, so `j<=3`.

At `j=3`, the quotient is a constant vector `B` and
`A_0B=x^3*grad(a)`.  Curl-freeness makes `a` a power of `x`, and the Section
4 normal forms, exact multiplicity six, and boundary rank two leave only

\[
D=8,
\qquad
f=Cx^8+H_8(y,z).
\tag{5.58}
\]

Write the forced moving terms as `h_5=alpha*w*x^4+...` and
`h_2=(q/2)*w^2+...`.  The scalar Schur coefficient fixes

\[
q=\frac{2\alpha^2}{7C}.
\]

The maximal passive sector is the binary Hessian of

\[
\lambda^6Cx^8+\lambda^3\alpha wx^4+\frac q2w^2,
\]

namely

\[
(56Cq-16\alpha^2)\lambda^6x^6
+12\alpha q\lambda^3wx^2.
\tag{5.59}
\]

After the first coefficient vanishes, the second is
`(24/7)*(alpha^3/C)*lambda^3*w*x^2`.  In the full determinant it is multiplied
by `lambda^12*R`.  Passive degree twelve forces both top passive Hessian
entries.  The remaining binary weight three has partitions `3+0` and `2+1`;
the latter has no `ww` entry because the first scalar Schur term occurs only
at order six.  Thus the terminal channel is immutable and nonzero, excluding
`j=3`.

At `j=2`, the quotient is a linear field `L=M(x,y,z)^T`.  Put

\[
F=\frac{m-1}{m+1}a=x^{-3}L(f).
\]

The Hessian equation, the gradient identity, and the scalar Schur coefficient
give

\[
\boxed{
L(f)=x^3F,
\qquad
M^T\nabla f=3x^2F\,dx-\frac{2}{m-1}x^3\nabla F,
\qquad
L(F)\in(x^3).
}
\tag{5.60}
\]

The boundary value of the middle identity again gives
`M^T*grad(f)|_(x=0)=0`.  If `M` has rank three, the boundary Hessian has rank
at most one.  Suppose it has rank two.  If `ker(M^T)` has a tangent
component, use the normal form `f=y^(m+2)+x^2g`; the order-one `dx` component
of (5.60) forces the coefficient of `x*partial_x` in `L` to vanish, leaving
`rank(M)<=1`.  If the kernel is normal to `x=0`, write
`f=x*p_(m+1)(y,z)+x^2g`.  The order-one part of (5.60) forces both the
boundary field and the constant tangent part of `L` to annihilate `p` and
its gradient.  Taking `p=y^(m+1)` leaves only one nonzero row in `M`, again a
contradiction.  Hence `rank(M)=1`.

Write `L=lambda*u`.  Equation (5.60) gives

\[
\frac{dF}{F}
=\frac{3(m-1)}2\frac{dx}{x}
-\frac{m-1}2\frac{d\lambda}{\lambda}.
\tag{5.61}
\]

Polynomiality of `F` forces `lambda` proportional to `x`; otherwise the
residue along `lambda=0` is negative.  Thus `F=C*x^(m-1)` and
`D_u(f)=C'*x^(m+1)`.  Exact multiplicity six and boundary rank two again
force

\[
D=8,
\qquad
f=Cx^8+H_8(y,z),
\qquad
L=x\partial_x.
\tag{5.62}
\]

Now the forced weighted sector is

\[
\lambda^6Cx^8
+\lambda^4\alpha wx^5
+\lambda^2Bw^2x^2,
\qquad
B=\frac{25\alpha^2}{112C}.
\]

Its binary Hessian determinant is

\[
(112BC-25\alpha^2)\lambda^8x^8
-12B^2\lambda^4w^2x^2.
\tag{5.63}
\]

The terminal multiplier is
`-(1875/3136)*(alpha^4/C^2)`.  Passive degree twelve again forces the two top
passive entries.  The alternative binary partition `3+1` cannot retain two
powers of `w`: the order-five scalar coefficient is independent of `w`, so
the weight-one layer has no cubic `w` term.  Hence this terminal channel is
also immutable and `j=2` is impossible.

The sole generic-corank-one sextuple input is therefore `j=1`.  The quotient
in (5.57) is a homogeneous quadratic vector field `Q`, and multiplication by
`A_0` together with the scalar coefficient gives the exact remaining system

\[
\boxed{
A_0Q=x^3\nabla a,
\qquad
x^3\mid Q(a),
\qquad
\deg Q=2,
\qquad
\deg a=m.
}
\tag{5.64a}
\]

This has an order-one matrix form.  Put

\[
F=\frac{m}{m+1}a=x^{-3}Q(f),
\qquad
N=\operatorname{Jac}(Q).
\]

Since `Q` is homogeneous quadratic, `N` is linear and `N(x,y,z)^T=2Q`.
The gradient identity and the scalar coefficient turn (5.64a) into

\[
\boxed{
Q(f)=x^3F,
\qquad
N^T\nabla f=3x^2F\,dx-\frac1m x^3\nabla F,
\qquad
Q(F)\in(x^3).
}
\tag{5.64b}
\]

On `x=0`, the middle identity gives `N^T*grad(f)=0`.  If `N mod x` were
invertible over `K(y/z)`, then the boundary value and first `x`-coefficient
of `f` would vanish, making the boundary Hessian rank at most one.  Hence

\[
\operatorname{rank}(N\bmod x)\le2.
\tag{5.64c}
\]

Thus the surviving generic sextuple packet is precisely a rank-at-most-two
order-one linear-matrix boundary system attached to a homogeneous quadratic
field.

The rank-zero subcase closes immediately.  If `N mod x=0`, then every entry
of the linear matrix `N` is divisible by `x`.  Jacobian integrability of `N`
forces

\[
Q=x^2u
\]

for a constant vector `u`.  Equation (5.64a) reduces to
`A_0u=x*grad(a)`.  The constant-vector normal forms, exact multiplicity six,
and boundary rank two give

\[
D=8,
\qquad
f=Cx^8+H_8(y,z),
\qquad
u=\partial_x.
\tag{5.64d}
\]

Write the forced maximal-passive ladder as

\[
\Phi_0
=\lambda^6Cx^8
+\lambda^5\alpha wx^6
+\lambda^4Bw^2x^4
+\lambda^3\gamma w^3x^2
+\lambda^2\delta w^4.
\]

Its binary Hessian determinant has successive coefficients

\[
\begin{aligned}
&112BC-36\alpha^2,\\
&-36B\alpha+336C\gamma,\\
&-40B^2+672C\delta+108\alpha\gamma,\\
&-20B\gamma+360\alpha\delta,\\
&144B\delta-24\gamma^2,\\
&24\delta\gamma,
\end{aligned}
\tag{5.64e}
\]

in the monomials `lambda^10*x^10`, `lambda^9*w*x^8`, through
`lambda^5*w^5`.  The first three determine

\[
B=\frac{9\alpha^2}{28C},
\qquad
\gamma=\frac{27\alpha^3}{784C^2},
\qquad
\delta=\frac{27\alpha^4}{43904C^3}.
\]

The next two coefficients then vanish identically, while the terminal one is

\[
24\delta\gamma
=\frac{2187\alpha^7}{4302592C^5}\ne0.
\tag{5.64f}
\]

In the full determinant this is multiplied by `lambda^12*R`.  Passive degree
twelve and `x`-degree zero force the two top passive Hessian entries.  The
remaining binary weight five has partitions `3+2`, `4+1`, and `5+0`; only
`3+2` can retain five powers of `w` while supplying the `w` row and column.
Thus arbitrary remainders cannot enter (5.64f).  This proves `HC4-DIR11` and
sharpens (5.64c) to boundary rank one or two.

We next classify boundary rank one.  Write `Q_0=Q|_(x=0)`.  If `Q_0=0`,
then `Q=xL` for a homogeneous linear field `L`.  Substitution in (5.64b),
with `M=Jac(L)`, gives

\[
L(f)=x^2F,
\qquad
M^T\nabla f=2xF\,dx-\frac1m x^2\nabla F,
\qquad
L(F)\in(x^2).
\tag{5.64g}
\]

This is exactly the linear-field system used in `HC4-DIR5`.  Its rank
collapse and rank-one integration force `L=x*partial_x` and a split top.
Exact sextuple multiplicity then gives `m=6` and
`f=C*x^8+H_8(y,z)`.  But `Q=xL=x^2*partial_x` has
`rank(Jac(Q) mod x)=0`, already excluded by `HC4-DIR11`.  Thus the axial
rank-one case is empty.

Suppose now that `Q_0` is nonzero.  Choose a nonzero binary-quadratic
component `q_0` of `Q_0`.  The rank-one minors give
`d(Q_0)_i wedge dq_0=0`.  Euler's identity then shows that every component
of `Q_0` is a constant multiple of `q_0`.  The remaining minors say the
coefficient vector of the terms linear in `x` lies in the same constant
line.  Absorbing their common linear form into a quadratic lift `q` gives

\[
Q=u q+x^2v,
\qquad
N\bmod x=u(\nabla q|_{x=0})^T,
\tag{5.64h}
\]

where `u,v` are constant, `q` is quadratic, and
`q_0=q|_(x=0)` is nonzero.  Put

\[
U=D_u(f),\qquad V=D_v(f).
\]

The first two equations in (5.64b) become

\[
qU+x^2V=x^3F,
\qquad
U\,dq+2xV\,dx
=3x^2F\,dx-\frac1m x^3dF.
\tag{5.64i}
\]

There are two image orientations.  If `u(x)` is nonzero, a change preserving
`x=0` makes `u=partial_x`.  The boundary value of (5.64i) gives
`partial_x(f)|_(x=0)=0`.  If `x^s f_s` is the first positive `x`-term of
`f`, the tangent components of the second identity force `s>=4`.  The first
identity then forces `D_v(f|_(x=0))=0` and

\[
F=x^{s-4}(s q_0f_s)+O(x^{s-3}).
\]

The left side of the `dx` component of the second identity has order at
least `s-1`, whereas its right side has the nonzero term

\[
\left(3-\frac{s-4}{m}\right)sq_0f_s\,x^{s-2}dx.
\]

Here `s<=m+2`, so the displayed scalar cannot vanish.  This is a
contradiction.  If there is no positive `x`-term at all, (5.64i) instead
gives `F=0`, hence `a=0` and then `Q=0`, also a contradiction.  The normal
image orientation is therefore empty.

It remains that `u` is tangent.  Normalize `u=partial_z` and retain
`q_0!=0`.  Set

\[
t=\operatorname{ord}_x(\partial_zf).
\tag{5.64j}
\]

The tangent components of (5.64i) give `t>=3`.  Moreover `t` is finite:
otherwise `f` is independent of `z` and `det Hess(f)=0`.  Since the boundary
Hessian has rank two and its `z` row and column vanish, its `(x,y)` block is
invertible over `K(y)`.  The Schur complement has order at least
`min(t,2t-2)=t`.  Exact order six therefore gives `t<=6`.

Write `U_t`, `V_(t-2)`, and `F_(t-3)` for the first nonzero coefficients of
`partial_zf`, `D_vf`, and `F`.  Degree gives `t<=m+1`, and comparison of the
first nonzero terms in (5.64i) gives

\[
\begin{aligned}
q_0U_t+V_{t-2}&=F_{t-3},\\
2V_{t-2}&=\left(3-\frac{t-3}{m}\right)F_{t-3},\\
U_t\,dq_0&=-\frac1m dF_{t-3}.
\end{aligned}
\tag{5.64k}
\]

Consequently, with `e=(m+3-t)/2`,

\[
\boxed{
d\log F_{t-3}=e\,d\log q_0,
\qquad
U_t=-\frac em\frac{F_{t-3}}{q_0},
\qquad
V_{t-2}=\frac12\left(3-\frac{t-3}{m}\right)F_{t-3}.
}
\tag{5.64l}
\]

For `t<=5`, the leading term of `Q(F) in (x^3)` is
`x^(t-3) q_0*partial_z(F_(t-3))`.  Hence `F_(t-3)` is independent of `z`,
and (5.64l) forces `q_0` to be independent of `z`.  After scale,
`q_0=y^2`.  The four packets are therefore

\[
\begin{array}{c|c|c|c|c}
t & \text{range} & q_0 & F_{t-3} & U_t\\ \hline
3 & m\ge3 & y^2 & c y^m & c' y^{m-2}\\
4 & m\ge3 & y^2 & c y^{m-1} & c' y^{m-3}\\
5 & m\ge4 & y^2 & c y^{m-2} & c' y^{m-4}\\
6 & m\ge5 & \text{arbitrary} & c q_0^{(m-3)/2}
  & c' q_0^{(m-5)/2}
\end{array}
\tag{5.64m}
\]

The last row is interpreted divisor-wise through (5.64l).  If `q_0` is a
square, the displayed half-powers are powers of its linear square root and
every `m>=5` is allowed.  If `q_0` is not a square, polynomiality requires
`m` odd.  This proves `HC4-DIR12`.  It is an exact reduction, not yet a
nonexistence result for the four rows of (5.64m).

There is a global refinement of these leading equations.  Eliminate `V`
from (5.64i) and put

\[
S=\frac q{x^2},
\qquad
T=\frac F{x^m}.
\]

The resulting one-form identity is

\[
U\,dS=-\frac{x^{m+1}}m\,dT,
\qquad
dS\wedge dT=0.
\tag{5.64n}
\]

On the affine chart `x=1`, write `S=hat(q)(Y,Z)` and
`T=hat(F)(Y,Z)`.  The elementary quadratic centralizer lemma says that
`d hat(q) wedge d hat(F)=0` has exactly two forms.  If `hat(q)` is a
composite quadratic, an affine change makes `hat(q)` univariate in a linear
coordinate `Y`, and then `hat(F)` is also in `K[Y]`.  Otherwise `hat(q)` is
primitive and `hat(F) in K[hat(q)]`.  For completeness, after algebraic
closure the primitive quadratic normal forms are a nondegenerate quadratic
or `Y^2+Z`; the kernels of their Hamiltonian derivations are respectively
`K[hat(q)]`.  This also descends to `K`.

In the composite case, homogenization and `q_0!=0` give a homogeneous linear
form `l`, with nonzero boundary value, such that

\[
q,F,U\in K[x,l].
\tag{5.64o}
\]

In particular `q_0` is a square.  If `D_u(l)=0`, choose `l=y` and the
complementary tangent coordinate along `u` to obtain

\[
f=zG(x,y)+h(x,y),
\qquad G=U,
\qquad \operatorname{ord}_xG=t.
\tag{5.64o'}
\]

This orientation is automatic for `t<6`, because (5.64m) says
`D_u(q_0)=0`.  At `t=6` there is one other composite orientation:
`D_u(l)!=0`.  Taking a tangent coordinate `r` invariant under `u` then gives
`f=P(x,l)+h(x,r)`, with `D_u(P)=U`.  These are the composite pencils that
remain in the four rows of (5.64m).

In the primitive case the centralizer lemma and homogeneity give

\[
F=\sum_{0\le k\le\lfloor m/2\rfloor}c_kx^{m-2k}q^k.
\tag{5.64p}
\]

Equation (5.64n) shows that a term of index `k` contributes order
`m+3-2k` to `U`.  Therefore the largest nonzero index is exactly

\[
k_{\max}=\frac{m+3-t}{2}.
\tag{5.64q}
\]

This must be an integer.  The primitive parity rows are consequently

\[
\begin{array}{c|c|c}
t& m& k_{\max}\\ \hline
3&\text{even},\ m\ge4&m/2\\
4&\text{odd},\ m\ge3&(m-1)/2\\
5&\text{even},\ m\ge4&(m-2)/2\\
6&\text{odd},\ m\ge5&(m-3)/2.
\end{array}
\tag{5.64r}
\]

For `t<6`, (5.64m) already gives `q_0=y^2`.  Write
`q=y^2+x(ay+bz)+cx^2`.  Primitivity is precisely `b!=0`; shifts and scales
preserving `x=0` and `partial_z` then give

\[
q=y^2+xz.
\tag{5.64s}
\]

This proves `HC4-DIR13`.

The invariant composite pencil admits a further exact gate.  For
`f=zG(x,y)+h(x,y)`, direct expansion gives

\[
\begin{aligned}
\Delta={}&-G_x^2(h_{yy}+zG_{yy})
+2G_xG_y(h_{xy}+zG_{xy})\\
&-G_y^2(h_{xx}+zG_{xx}),\\
[z]\Delta={}&-\frac{m+1}{m}G\det\operatorname{Hess}_{x,y}(G).
\end{aligned}
\tag{5.64t}
\]

If `ord_x(G)=t`, the `h` part has order at least `2t-2`; the `z` part has
order at least `3t-2` (and possibly higher when the leading binary monomial
is a pure `x` power).  Exact order six therefore excludes `t=5,6`.

For `t=3`, write the first binary coefficients as

\[
\begin{aligned}
G&=c_0x^3y^{m-2}+c_1x^4y^{m-3}+O(x^5),\\
h&=a_0y^{m+2}+a_1xy^{m+1}+a_2x^2y^m+O(x^3),
\end{aligned}
\qquad c_0\ne0.
\]

The first two determinant coefficients are

\[
\begin{aligned}
[x^4]\Delta
&=-9a_0c_0^2(m+1)(m+2)y^{3m-4},\\
[x^5]\Delta\big|_{a_0=0}
&=-3a_1c_0^2(m+1)(m+4)y^{3m-5}.
\end{aligned}
\tag{5.64u}
\]

Exact order six first forces `a_0=0`.  Boundary Hessian rank two then forces
`a_1!=0`, so the second line of (5.64u) is nonzero in characteristic zero.
Thus `t=3` is impossible.

At `t=4`, write
`G=c_0x^4y^(m-3)+O(x^5)` and retain the notation
`h|_(x=0)=a_0y^(m+2)`.  Now

\[
[x^6]\Delta
=-16a_0c_0^2(m+1)(m+2)y^{3m-6}.
\tag{5.64v}
\]

Exact order six is therefore equivalent at the leading boundary to
`a_0!=0`.  This proves `HC4-DIR14`.

It remains to close that `t=4` pencil.  Return to

\[
Q=q\partial_z+x^2(v_x\partial_x+v_y\partial_y+v_z\partial_z),
\qquad
f=zG+h.
\]

The coefficient of `z` in `Q(f)=x^3F` is

\[
x^2(v_xG_x+v_yG_y),
\]

so

\[
v_xG_x+v_yG_y=0.
\tag{5.64w}
\]

If `v_x=v_y=0`, absorb `x^2v_z` into `q`; then `Q=q*partial_z` and
`F=x^(-3)qG`.  Substitution in the middle equation of (5.64b) integrates to

\[
Gq^{m+1}=C x^{3(m+1)}.
\tag{5.64x}
\]

This is impossible because `q_0!=0`.

Hence `(v_x,v_y)` is nonzero.  Equation (5.64w) makes the homogeneous binary
form `G` a power of a linear form.  Since `ord_x(G)=4`, necessarily

\[
m=3,
\qquad
G=Cx^4,
\qquad
v_x=0,
\qquad
v_y\ne0.
\]

After again absorbing `v_z`, the coefficient of `x^2` in `Q(f)=x^3F` is
`v_y*h_y|_(x=0)`.  But (5.64v) requires
`h|_(x=0)=a_0y^5` with `a_0!=0`, so this coefficient is nonzero.  This final
contradiction proves `HC4-DIR15`.

Finally consider the transverse composite orientation.  Take `l=y`, a
tangent coordinate `z` invariant under `u`, and normalize `u=partial_y`.
Then

\[
q,F,U\in K[x,y],
\qquad
f=P(x,y)+h(x,z),
\qquad
U=P_y,
\qquad
\operatorname{ord}_xU=6.
\tag{5.64y}
\]

Write `v=v_x*partial_x+v_y*partial_y+v_z*partial_z`.  All terms of
`Q(f)=x^3F` except

\[
x^2(v_xh_x+v_zh_z)
\]

belong to `K[x,y]`.  Hence

\[
D_{(v_x,v_z)}h\in K[x].
\tag{5.64z}
\]

By homogeneity it equals `c*x^(m+1)`.  Differentiating (5.64z) gives

\[
\operatorname{Hess}_{x,z}(h)
\binom{v_x}{v_z}
=c(m+1)x^m\binom10.
\tag{5.64aa}
\]

On `x=0`, the `(x,z)` Hessian is invertible: the `y` row and column of
`Hess(f)` vanish there, while its boundary rank is two.  Equation (5.64aa)
therefore forces `v_x=v_z=0`.  Absorb the remaining `x^2v_y` into `q`.
Then `Q=q*partial_y`, and exactly the calculation (5.64x), with `G=U`, gives

\[
Uq^{m+1}=C x^{3(m+1)},
\]

again contradicting `q|_(x=0)!=0`.  This proves `HC4-DIR16`.

There is a uniform closure for the primitive rows with `t>=4`.  Equation
(5.64k) and the nonzero scalar in its second line give

\[
\operatorname{ord}_x(V)=t-2\ge2.
\]

Therefore

\[
\operatorname{Hess}(f)v=\nabla V=0\pmod{x}.
\tag{5.64ab}
\]

Likewise `U=D_u(f) in (x^t)` gives `Hess(f)u=0 mod x`.  The boundary Hessian
has rank two, so its kernel is the line `K*u`; hence `v=c*u`.  Absorbing
`c*x^2` into `q` reduces to `Q=q*u`, and (5.64x) again gives

\[
Uq^{m+1}=C x^{3(m+1)},
\]

contrary to `q_0!=0`.  Thus the primitive rows `t=4,5,6` are empty.

The argument just used is degree-free and will be useful beyond sextuple
multiplicity:

> **Kernel-locking lemma.**  Let a polynomial Hessian in `n` variables have
> boundary rank `n-1` on `x=0`.  If constant vectors `u,v` satisfy
> `D_u(f),D_v(f) in (x^2)`, then `u` and `v` are proportional.  Indeed their
> gradients put both vectors in the one-dimensional kernel of
> `Hess(f) mod x`.

No homogeneity or quadratic-field hypothesis is needed for this lemma.

At `t=3`, the composite row was already eliminated by (5.64u).  In the
primitive row, (5.64r) gives even `m`, `e=m/2`, and

\[
q=y^2+xz,
\qquad
F=cq^e+O(x^2),
\qquad c\ne0.
\]

The coefficient of `x` in `Q(F)` is then

\[
[x]Q(F)=ce\,y^{2e}=\frac{cm}{2}y^m\ne0.
\tag{5.64ac}
\]

The `x^2v` part and all lower powers of `q` start at order at least two, so
they cannot alter (5.64ac).  This contradicts `Q(F) in (x^3)` and proves
`HC4-DIR17`.

We now reduce that rank-two target.  Put

\[
N_0=N\bmod x,
\qquad
g_0=\nabla f\bmod x
=\bigl(f_1,(f_0)_y,(f_0)_z\bigr),
\tag{5.64ad}
\]

where `f_0=f|_(x=0)` and `f_1=[x]f`.  The boundary value of (5.64b) is

\[
N_0^Tg_0=0.
\tag{5.64ae}
\]

The matrix `N_0` is linear and has rank two.  Its two-by-two minors give a
homogeneous polynomial generator `chi` of `ker(N_0^T)`.  After removing their
common factor, `chi` is primitive and

\[
k=\deg(\chi)\le2.
\]

Since boundary Hessian rank two prevents `g_0=0`, unique factorization in
`K[y,z]` gives

\[
g_0=H\chi,
\qquad
\deg H=m+1-k.
\tag{5.64af}
\]

Write `chi=(chi_x,chi_y,chi_z)`.  If `chi_y=chi_z=0`, primitivity makes
`chi` a constant normal vector.  Then `(f_0)_y=(f_0)_z=0`, hence `f_0=0`.
Moreover `chi^T N_0=0` says the normal row of `N_0` vanishes.  Because
`N=Jac(Q)`, this is exactly

\[
Q_x\in(x^2).
\tag{5.64ag}
\]

This is the normal packet (0.5).

Otherwise put

\[
s=\deg\gcd(\chi_y,\chi_z)\le k.
\]

Equation (5.64af) gives

\[
\deg\gcd((f_0)_y,(f_0)_z)=m+1-k+s.
\tag{5.64ah}
\]

If a nonzero binary form of degree `D=m+2` has `r` distinct projective roots,
characteristic zero gives

\[
\deg\gcd((f_0)_y,(f_0)_z)=D-r.
\]

Comparison with (5.64ah) yields the exact count

\[
\boxed{r=1+k-s\le3.}
\tag{5.64ai}
\]

After algebraic closure these are precisely the three profiles (0.6).  This
proves `HC4-DIR18`.  Notice that the proof used only a linear rank-two
boundary matrix and a polynomial gradient in its left kernel; the root bound
is independent of sextuple multiplicity.

We can sharpen the one-root profile.  Normalize `f_0=y^(m+2)`.  Before any
further change, the boundary Hessian determinant is

\[
-(m+2)(m+1)y^m((f_1)_z)^2.
\]

Boundary rank two therefore gives `(f_1)_z=0`.  Since `f_1` is homogeneous,
it is a multiple of `y^(m+1)`, and a shear of `y` by `x` removes it.  We may
write

\[
f=y^{m+2}+x^2g_m(y,z)+O(x^3),
\qquad g_m\ne0.
\tag{5.64aj}
\]

The boundary Hessian kernel is now `partial_z`.  Since `A_0Q_0=0`, and since
the second row of `N_0` annihilates the boundary gradient, the quadratic
field has the form

\[
\begin{aligned}
Q_x&=x\ell(y,z)+\alpha x^2,\\
Q_y&=\beta x^2,\\
Q_z&=q_2(y,z)+xr_1(y,z)+\gamma x^2,
\end{aligned}
\qquad \ell\ne0.
\tag{5.64ak}
\]

The coefficient of `x^2` in `Q(f)=x^3F` is

\[
q_2(g_m)_z+2\ell g_m+\beta(m+2)y^{m+1}=0.
\tag{5.64al}
\]

Independently, the first possible Hessian-determinant coefficient is

\[
[x^2]\Delta
=2(m+2)(m+1)y^m
\left(g_m(g_m)_{zz}-2(g_m)_z^2\right).
\tag{5.64am}
\]

Exact order six forces the parenthesis to vanish.  Over `K(y)`,

\[
\partial_z\left(\frac{(g_m)_z}{g_m^2}\right)
=\frac{g_m(g_m)_{zz}-2(g_m)_z^2}{g_m^3}=0.
\]

If `(g_m)_z` were nonzero, `1/g_m` would be a nonconstant affine polynomial
in `z`, impossible for the reciprocal of a polynomial.  Hence

\[
g_m=C y^m,
\qquad C\ne0.
\]

Equation (5.64al) then gives

\[
\ell=\lambda y,
\qquad
\beta=-\frac{2C\lambda}{m+2},
\qquad
\lambda\ne0.
\tag{5.64an}
\]

This proves `HC4-DIR19`.

But (0.7) is already inconsistent with the middle equation in (5.64b).  Its
`dx` component at order `x` is

\[
[x]\bigl(N^T\nabla f\bigr)_x
=2C\lambda y^{m+1}+2\beta(m+2)y^{m+1}.
\tag{5.64ao}
\]

The right side of (5.64b) starts at order `x^2`, so (5.64ao) must vanish.
Substitution of (5.64an) leaves

\[
-2C\lambda y^{m+1}\ne0.
\]

This contradiction proves `HC4-DIR20` and closes the one-root profile.

It remains to normalize the first jets of the two other root profiles.  For
`f_0=y^a z^b`, put

\[
G=y^{a-1}z^{b-1}=gcd((f_0)_y,(f_0)_z).
\]

The root-count identity (5.64ai) reads `k-s=1`.  Since `k<=2`, there are two
possibilities.  If `(k,s)=(1,0)`, then `H=G` and the tangent components of
the primitive kernel generator are

\[
\chi_y=az,
\qquad
\chi_z=by.
\tag{5.64ap}
\]

They span all linear forms, so a shear
`y -> y+p*x`, `z -> z+q*x` removes `chi_x` and hence `f_1`.

If `(k,s)=(2,1)`, the common linear factor of `chi_y,chi_z` must divide
`G`.  For the factor `y` (which requires `a>=2`), one has

\[
H=y^{a-2}z^{b-1},
\qquad
(\chi_y,\chi_z)=(ayz,by^2).
\tag{5.64aq}
\]

The same shears remove the `yz` and `y^2` components of `chi_x`; the sole
complement is `z^2`.  Thus

\[
f_1=\kappa y^{a-2}z^{b+1}.
\]

The common factor `z` gives symmetrically
`f_1=\kappa y^{a+1}z^{b-2}`.  This proves (0.8).

For `f_0=y^a z^b(y-z)^c`, put

\[
G=y^{a-1}z^{b-1}(y-z)^{c-1}.
\]

Now (5.64ai) forces `(k,s)=(2,0)` and `H=G`.  Direct differentiation gives

\[
\begin{aligned}
\chi_y&=z\bigl((a+c)y-az\bigr),\\
\chi_z&=y\bigl(by-(b+c)z\bigr).
\end{aligned}
\tag{5.64ar}
\]

These coprime quadratics span a two-plane in `K[y,z]_2`; together with `yz`
they form a basis.  The two shear parameters remove the first two basis
directions from `chi_x`, leaving

\[
f_1=\kappa Gyz
=\kappa y^a z^b(y-z)^{c-1}.
\tag{5.64as}
\]

This proves `HC4-DIR21`.

The boundary Hessian now fixes the second jet.  Both the two- and three-root
binary Hessians

\[
B_0=\operatorname{Hess}_{y,z}(f_0)
\]

are invertible over `K(y/z)`.  Boundary rank two and the block Schur
complement therefore give

\[
2f_2=(\nabla f_1)^T B_0^{-1}\nabla f_1.
\tag{5.64at}
\]

For the first exceptional two-root jet in (0.8), direct substitution gives

\[
f_2
=\frac{\kappa^2(ab-a-4b)}{2ab}
y^{a-4}z^{b+2}.
\tag{5.64au}
\]

The numerator coefficient is nonzero for `a=2,3`, so polynomiality requires
`a>=4`.  The symmetric exceptional jet gives

\[
f_2
=\frac{\kappa^2(ab-4a-b)}{2ab}
y^{a+2}z^{b-4},
\qquad b\ge4.
\tag{5.64av}
\]

If the two-root first jet is zero, (5.64at) gives `f_2=0`.

For the three-root profile, set

\[
\begin{aligned}
R_2&=ab(y-z)^2+ac z^2+bc y^2,\\
P_2&=ab(y-z)^2+a(c-1)z^2+b(c-1)y^2.
\end{aligned}
\tag{5.64aw}
\]

The same Schur calculation gives

\[
f_2
=\frac{\kappa^2}{2}
y^a z^b(y-z)^{c-2}\frac{P_2}{R_2}.
\tag{5.64ax}
\]

The quadratic `R_2` is coprime to `y*z*(y-z)`, since its restrictions to the
three root lines are respectively nonzero multiples of `z^2`, `y^2`, and
`y^2`.  If `R_2` divided `P_2`, equal `yz` coefficients would force
`P_2=R_2`, but
`R_2-P_2=a z^2+b y^2` is nonzero.  Thus (5.64ax) is polynomial only when
`kappa=0`.  Equation (5.64at) then also gives `f_2=0`.

Finally suppose `f_1=f_2=0`, in either the zero-first-jet two-root packet or
the three-root packet.  If `x^s h_{D-s}` is the first positive `x`-term,
then `s>=3` and the Schur complement of `B_0` begins with

\[
s(s-1)x^{s-2}h_{D-s}.
\]

Hence

\[
\operatorname{ord}_x\Delta=s-2.
\]

Exact order six forces `s=8`, proving

\[
f=f_0+x^8h_{D-8}+O(x^9),
\qquad D\ge8.
\tag{5.64ay}
\]

This proves `HC4-DIR22`.

We next close the normal packet.  To avoid confusing a vector component with
a derivative, write `Q^x,Q^y,Q^z` for the components of `Q`.  Since `f_0=0`,
put

\[
f=xp_n(y,z)+x^2g_{n-1}(y,z)+x^3h_{n-2}(y,z)+O(x^4),
\qquad n=m+1.
\tag{5.64az}
\]

Boundary Hessian rank two makes `p_n` nonconstant.  Equation (5.64ag) and
homogeneity give

\[
Q^x=\alpha x^2,
\qquad
(Q^y,Q^z)=q(y,z)+x l(y,z)+x^2c,
\tag{5.64ba}
\]

where `q` is quadratic, `l` is linear, and `c` is constant.  The coefficients
of order `x` and `x^2` in the first equation of (5.64b), and the order-`x`
normal and tangent components of its middle equation, give

\[
q(p)=0,
\qquad
\alpha p+q(g)+l(p)=0,
\qquad
2\alpha p+l(p)=0,
\qquad
\operatorname{Jac}(q)^{\mathsf T}\nabla p=0.
\tag{5.64bb}
\]

The last identity and `rank(N mod x)=2` force `Jac(q)` to have rank one.
Indeed rank two would kill `grad(p)`, while rank zero would leave only the
single `x`-column of `N mod x`.  A homogeneous binary quadratic map of
Jacobian rank one has constant image line: Euler contraction of
`dq_1 wedge dq_2=0` makes its two components proportional.  Thus
`q=q_2u` for a constant tangent vector `u`.  Now `q(p)=0` says `D_u p=0`, so

\[
p=CM^n
\tag{5.64bc}
\]

for a linear form `M` with `D_uM=0`.  The third identity in (5.64bb) gives
`l(M)=-(2*alpha/n)M`.  If `alpha=0`, both `q` and `l` take values in the line
`Ku`, again making `rank(N mod x)<=1`.  Hence `alpha!=0`.  The other two
identities in (5.64bb) reduce to

\[
q_2D_ug=\alpha CM^n.
\tag{5.64bd}
\]

Unique factorization forces `q_2=rho*M^2`, with `rho!=0`.  Normalize
`M=y`, `u=partial_z`.  Then, for constants `A,r,s`,

\[
\begin{aligned}
p&=Cy^n,\\
q&=\rho y^2\partial_z,\\
g&=B y^{n-2}z+A y^{n-1},
\qquad B=\frac{\alpha C}{\rho},\\
l&=-\frac{2\alpha}{n}y\partial_y+(ry+sz)\partial_z.
\end{aligned}
\tag{5.64be}
\]

Let `F_0=F mod x`.  The coefficient of `x^3` in `Q(f)=x^3F` is

\[
F_0=\rho y^2h_z+U y^{n-1}
       +B\left(s+\frac{4\alpha}{n}\right)y^{n-2}z
\tag{5.64bf}
\]

for a scalar `U` whose value is irrelevant.  The `z`-component at order
`x^3` in the middle equation of (5.64b) gives

\[
\rho y^2h_{zz}
=-B\left(ns+\frac{4\alpha}{n}\right)y^{n-2}.
\tag{5.64bg}
\]

Substitution into (5.64bf) leaves

\[
F_0=V y^{n-1}-B(n-1)s y^{n-2}z.
\tag{5.64bh}
\]

The boundary value of the scalar condition `Q(F) in (x^3)` is therefore

\[
q(F_0)=-\rho B(n-1)s y^n,
\tag{5.64bi}
\]

so `s=0`.  Return to the `y`-component of the same order-`x^3` middle
equation.  Its coefficient of `y^(n-3)z` is now

\[
-\frac{2\alpha B(n+2)}{n},
\tag{5.64bj}
\]

whereas the right side `-(1/(n-1))*partial_y(F_0)` has no such monomial.
This is nonzero because `alpha*B!=0`.  The contradiction proves
`HC4-DIR23`.

There is also a short closure that applies simultaneously to all the delayed
root packets.  Suppose the binary boundary Hessian
`B_0=Hess_(y,z)(f_0)` is invertible and `f_1=f_2=0`.  The boundary value of
`Hess(f)Q=x^3*grad(a)` first forces

\[
Q_0=(q_2(y,z),0,0).
\tag{5.64bk}
\]

Write the first tangent coefficient of `Q` as `x*l(y,z)`.  The coefficient
of `x` in `Q(f)=x^3F` and the tangent part at order `x` of the middle equation
in (5.64b) say respectively

\[
l(f_0)=0,
\qquad
\operatorname{Jac}(l)^{\mathsf T}\nabla f_0=0.
\tag{5.64bl}
\]

Differentiating the first identity and using the second gives

\[
B_0l=0.
\tag{5.64bm}
\]

Thus `l=0` over `K(y/z)`.  But then `N mod x` has only its normal row and has
rank at most one, contradicting the rank-two packet.  This proves
`HC4-DIR24`.  Notice that neither the value eight in (5.64ay) nor the number
of boundary roots entered the argument.  The reusable statement is that two
missing normal jets plus an invertible tangent Hessian collapse a quadratic
boundary Jacobian from rank two to rank at most one.

It remains only to close the two outer jets.  Treat

\[
f_0=y^az^b,
\qquad
f_1=\kappa y^{a-2}z^{b+1},
\qquad
f_2=\frac{\kappa^2(ab-a-4b)}{2ab}y^{a-4}z^{b+2},
\tag{5.64bn}
\]

where `kappa!=0`; the zero value is already covered by `HC4-DIR24`.
The kernel of the boundary Hessian is generated over `K(y/z)` by

\[
\left(1,-\frac{2\kappa z}{ay},
          \frac{\kappa z^2}{by^2}\right).
\tag{5.64bo}
\]

Since `Q_0` is a nonzero homogeneous quadratic vector, polynomiality fixes,
after rescaling `Q` and `F`,

\[
Q_0=\left(y^2,-\frac{2\kappa}{a}yz,
                 \frac{\kappa}{b}z^2\right).
\tag{5.64bp}
\]

Write

\[
Q_1=(Ay+Bz,\ Cy+Dz,\ Ey+Gz).
\tag{5.64bq}
\]

The coefficient of `x` in `Q(f)=x^3F` is zero.  Its four monomial
coefficients give

\[
B=E=0,
\qquad
\kappa A+aD=0,
\qquad
aC+bG=0.
\tag{5.64br}
\]

After these substitutions, the two tangent components at order `x` of the
middle equation in (5.64b) are a monomial multiple of the same polynomial

\[
a^2bC y^2-a^2bD yz-(2a+4b)\kappa^2z^2.
\tag{5.64bs}
\]

The right side has no tangent term before order `x^3`, so this polynomial
must vanish.  Its first two coefficients give `C=D=0`, while its last is
nonzero for positive `a,b` and `kappa!=0`.  This contradiction closes the
first outer jet.  Interchanging `y,a` with `z,b` closes the second and proves
`HC4-DIR25`.

Thus every generic-corank-one sextuple packet is empty.  The three lower-rank
Hessian packets in `HC4-DIR10` remain separate.  This strengthens the
reduction asserted in `HC4-DIR9` to a complete generic-boundary closure.

### 5.9 The lower-rank sextuple reduction

If the boundary Hessian has rank zero, integrability gives `x^3|f` and hence
`ord_x(Delta)>=7`.  Thus exact order six has boundary rank one.  Its invariant
factors are

\[
(0,1,5),\qquad(0,2,4),\qquad(0,3,3).
\]

The same DVR calculation now gives the stronger divisibility

\[
x^4R\mid\operatorname{adj}(A_0)b_j,
\qquad
j\le2.
\tag{5.65}
\]

At `j=2`, the quotient is constant and

\[
A_0B=x^2\nabla a.
\tag{5.66}
\]

The transverse normal form has exact multiplicity `m` but boundary rank two,
while the tangent form has multiplicity at least `2m`.  Exact multiplicity six
therefore leaves only `m=3`, `D=5`, and

\[
f=Czx^4+h_5(x,y),
\qquad
(h_5)_{yy}|_{x=0}\ne0,
\qquad
B=\partial_z.
\tag{5.67}
\]

This packet is genuinely isotropic: `B(a)=0`, so its first scalar Schur norm
vanishes.  It is not excluded by the high-order channels above.

At `j=1`, the quotient is linear and gives the order-one system (5.24).
The two rank-one boundary forms are still (5.9)--(5.10).  Exact order six
sharpens their first unused jets as follows.  In the first form, after
normalizing the boundary pure power,

\[
\begin{aligned}
f&=x^2\bigl(y^m+xg_1+x^2g_2+O(x^3)\bigr),\\
g_1&=u y^{m-1}+v y^{m-2}z,\\
[x^6]\det\operatorname{Hess}f
&=(m+1)\left((m-8)v^2y^{3m-6}
-2m y^{2m-2}(g_2)_{zz}\right)\ne0.
\end{aligned}
\tag{5.68}
\]

In the pure-power-boundary form,

\[
\begin{aligned}
f&=y^{m+2}+x^3\bigl(c y^{m-1}+xg_1+x^2g_2+O(x^3)\bigr),\\
g_1&=u y^{m-2}+v y^{m-3}z,\\
[x^6]\det\operatorname{Hess}f
&=2(m+1)(m+2)\left(3c y^{2m-1}(g_2)_{zz}
-8v^2y^{3m-6}\right)\ne0.
\end{aligned}
\tag{5.69}
\]

Equations (5.67)--(5.69) are the complete lower-rank sextuple handoff.  This
proves `HC4-DIR10`; no nonexistence claim is made there for the residual
packets.  We now synchronize them.

First take (5.68) and use the order-one system (5.24).  The same boundary
argument as in (5.12)--(5.14) first gives

\[
L=\alpha x\partial_y+(\beta x+ry+sz)\partial_z.
\tag{5.70}
\]

Indeed divisibility makes `L(x)=lambda*x`; the normal component of (5.24)
and `L(F) in (x^2)` then force `lambda=0`, while
`L mod x` must annihilate `g_0=y^m`.  Write

\[
F=xF_1+x^2F_2+O(x^3),
\qquad
P=m\alpha+rv.
\]

The normal coefficient at order `x^2` and the two tangent coefficients at
order `x^3` give

\[
F_1=P y^{m-1}+svy^{m-2}z,
\qquad
m\alpha=\left(2-\frac1m\right)P,
\qquad
sv=0,
\qquad
rv=-\frac{m-1}{m}P.
\tag{5.71}
\]

Put `H=alpha*(g_1)_y+beta*(g_1)_z`.  At the next order, the normal equation,
the definition of `F`, and the `z`-tangent equation are

\[
F_2=\frac{m}{2(m-1)}H,
\qquad
F_2=H+(ry+sz)(g_2)_z,
\qquad
s(g_2)_z=-\frac1m(F_2)_z.
\tag{5.72}
\]

If `v!=0`, (5.71) gives `s=0`; the last equation in (5.72) then gives
`alpha=0`, after which (5.71) gives `r=0`.  The two formulas for `F_2` force
`beta=0` because `m>=3`, contradicting `L!=0`.  Hence `v=0`.  Equations
(5.71) then give `alpha=0`, and (5.72) gives

\[
(ry+sz)(g_2)_z=0.
\tag{5.73}
\]

If the boundary field is nonzero, `(g_2)_z=0`.  Otherwise
`L=beta*x*partial_z`; comparison of the next normal coefficient gives either
`(g_2)_z=0` or the sole weight resonance `m=3`, where `g_2` is linear and
still `(g_2)_{zz}=0`.  Thus `v=(g_2)_{zz}=0` in (5.68), making its displayed
exact-sextuple coefficient zero.  The first order-one packet is empty.

Now take (5.69).  The boundary value of the middle equation in (5.24) kills
the complete `y`-component of `L`.  Writing its other components generally,
the normal coefficients of orders `x` and `x^2` successively kill the
boundary part and then the `x`-coefficient of `L(x)`.  Hence

\[
L=(\delta x+ey+sz)\partial_z.
\tag{5.74}
\]

If `v!=0`, then `F` starts at order two with
`v(ey+sz)y^(m-3)`.  The next normal equation kills `e=s=0`.  Thus
`L=delta*x*partial_z`, and the following normal coefficient reads

\[
1=2-\frac3m.
\tag{5.75}
\]

It forces `m=3`.  If `v=0`, the identical comparison one order later makes
`(g_2)_{zz}=0`; its only apparent resonance is `m=4`, where `g_2` is linear
and the same conclusion holds.  Thus the exact coefficient in (5.69) is
nonzero only at `m=3`, with `v!=0`.  After rescaling `L`, that packet is

\[
f=y^5+cx^3y^2+x^4(uy+vz)+dx^5,
\qquad
L=x\partial_z,
\qquad
F=vx^3.
\tag{5.76}
\]

Direct differentiation gives

\[
\det\operatorname{Hess}f
=-32v^2x^6(cx^3+10y^3).
\tag{5.77}
\]

Squarefreeness of the cofactor requires `c!=0`.  Formula (5.76) already has
the form (5.67), with its `z*x^4` coefficient equal to `v`.  Therefore all
three lower-rank handoffs have one top geometry,

\[
\boxed{f=Czx^4+h_5(x,y),}
\tag{5.78}
\]

carrying either the order-two direction `partial_z` or the resonant order-one
direction `x*partial_z`.  This proves `HC4-DIR26`.  It is a synchronization
theorem, not yet a nonexistence theorem for the common family (5.78).

The order-two side of (5.78) nevertheless closes at the level relevant to
`HC4`.  Here `D=5`, so `m=3`, and the first motion occurs at `j=2`.  Since
every bottom-right coefficient below order `2j=4` vanishes in (1.4), while
the scaled Hessian polynomial stops at order three, its bottom-right block is
identically zero.  The complete four-variable potential is therefore

\[
\Psi=wP(x,y,z)+H(x,y,z).
\tag{5.79}
\]

Moreover (5.66), (5.67), and `B=partial_z` give

\[
\nabla f_z=4Cx^3dx=x^2\nabla a,
\qquad
a=2Cx^2.
\tag{5.80}
\]

The only later `w`-linear layer has degree one, hence

\[
P=2Cx^2+\ell(x,y,z)
\tag{5.81}
\]

for an affine-linear `ell`.  Put `g=grad(P)`.  If the `(y,z)` part of
`ell` is zero, then `g=(4Cx+ell_x,0,0)`.  The bordered determinant of
`Hess(Psi)` is exactly

\[
\det\operatorname{Hess}\Psi
=-(4Cx+\ell_x)^2
  \det\operatorname{Hess}_{y,z}H.
\tag{5.82}
\]

This cannot be a nonzero constant because its first factor is nonconstant.
If the tangent part of `ell` is nonzero, a constant tangent vector `u`
satisfies `D_uP in K^times`.  Thus (5.79) is precisely a quadratic
zero-corner scalar parent with a unit pivot direction.  The registered
theorem `HC4RSD12` puts `P` in a graph coordinate, factors the bordered unit
into the Hessian determinant on each ternary fiber, and uses `HC3` to exclude
a gradient collision on that fiber.  Therefore no order-two member of
(5.78) is an HC4 counterexample.  This proves `HC4-DIR27` and leaves only the
order-one resonance (5.76) on the exact-sextuple lower-rank frontier.

That last resonance already identifies the broader next target.  Its first
`w`-linear coefficient is fixed by

\[
\operatorname{Hess}(f)(x\partial_z)
=4vx^4dx=x^2\nabla a_3,
\qquad
a_3=\frac{4v}{3}x^3.
\tag{5.83}
\]

The scalar norm is zero, so the complete potential has the scalar-parent
form

\[
\Psi=H(x,y,z)+wP(x,y,z)+\frac{\eta}{2}w^2,
\qquad
P=\frac{4v}{3}x^3+P_{\le2}.
\tag{5.84}
\]

Thus the next problem is no longer a list of sextuple boundary jets: it is
the classification of scalar parents whose pivot has pure-cube leading form.
The zero-corner case `eta=0` and the nonzero-corner reverse-Schur pencil
`eta!=0` are the two natural branches.  A uniform closure of (5.84) would
also provide the right template for higher repeated-factor packets whose
extremal motion is forced to a one-polynomial composite by (6.4)--(6.5).

## 6. Verification and next boundary

The same argument gives a useful budget on every repeated-factor stratum.
Factor

\[
\Delta=\prod_i\pi_i^{e_i}
\]

and assume `A0` has generic corank one along each `pi_i`.  Put

\[
H=\prod_i\pi_i^{\lceil e_i/2\rceil},\qquad
G=\prod_i\pi_i^{\lfloor e_i/2\rfloor},\qquad
\kappa=\deg G.
\tag{6.1}
\]

Over the DVR at `pi_i`, a symmetric Schur complement reduces `A0` to a
two-by-two unit block plus one entry of valuation `e_i`.  Polynomiality of
`b_j^T A0^{-1}b_j` forces the singular component of `b_j` to have valuation
at least `ceil(e_i/2)`.  Therefore

\[
H\mid\operatorname{adj}(A_0)b_j.
\tag{6.2}
\]

Since `deg H=3m-kappa`, degree comparison gives the **multiplicity budget**

\[
\boxed{j\le\kappa.}
\tag{6.3}
\]

The squarefree theorem is the case `kappa=0`.  The double-linear theorem has
`kappa=1`, so `j=1`.

At the extremal order `j=kappa`, the quotient in (6.2) is a nonzero constant
vector `B`, and

\[
A_0B=G\nabla a,
\qquad
dG\wedge da=0,
\qquad
\deg a=m+1-\kappa.
\tag{6.4}
\]

If `G` is a closed/generative polynomial, then `a` belongs to `K[G]`.
Homogeneity consequently forces

\[
\kappa\mid m+1,
\qquad
a=cG^{(m+1-\kappa)/\kappa},
\qquad
D_Bf=c'G^{(m+1)/\kappa}.
\tag{6.5}
\]

Thus every repeated Hessian factor has a quantitative cost: the first flag
rotation must occur no later than the half-multiplicity degree `kappa`, and an
extremal rotation is already a one-polynomial composite.

### Reproduction

Run

```bash
.venv/bin/python scripts/verify_hc4_direct_double_linear_hessian_gate.py
```

The checker verifies the all-degree determinant identities (4.2)--(4.4), the
terminal triple packet (4.7)--(4.9), and (5.15), the boundary coefficients
(5.11) and (5.16), the degree comparisons,
the incompatible weights in (5.18), and the rank-one integration
(5.30)--(5.32).  It also verifies the four weighted-binary channels ending
in (5.39), both quintuple-linear channel ladders (5.43)--(5.46), and the two
order-five boundary coefficients (5.51) and (5.55).  It further verifies the
two sextuple high-order channel pairs (5.59) and (5.63), the logarithmic
residues (5.61), the lower-rank order-six coefficients (5.68)--(5.69), the
normal and outer rank-two contradictions (5.64az)--(5.64bs), the delayed-jet
product-rule gate, the lower-rank recurrence synchronization (5.70)--(5.77),
and the bordered scalar-pivot determinant (5.82).
The UFD
and DVR divisibility arguments, the two rank-two boundary-gradient
exclusions, and the determinant-channel uniqueness are written proof steps
and are not replaced by a bounded search.

The next repeated-factor strata are now explicit:

1. one nonlinear repeated Hessian factor;
2. the sole exact-sextuple order-one resonance (5.76), or a linear factor of
   multiplicity at least seven;
3. two or more distinct repeated factors.

These are the remaining rank-three top-cone inputs to the ternary Schur
dichotomy.  Rank-at-most-two top Hessians remain a separate synchronization
problem.
