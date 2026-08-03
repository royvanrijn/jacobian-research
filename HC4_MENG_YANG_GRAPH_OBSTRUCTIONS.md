# Affine and polynomial-graph obstructions for the Meng--Yang `HC(5)` family

## Status

No `HC(4)` counterexample is constructed here.  This note proves the
following exact results for attempts to retain the rational collision of the
five-variable Meng--Yang potential:

> **Theorem `HC4MYA1` -- affine-hyperplane obstruction.**  Every affine
> four-plane restriction of the scaled Meng--Yang family is either
> Hessian-degenerate or has nonconstant Hessian determinant.

> **Theorem `HC4MYG3` -- cubic-graph obstruction.**  No graph
> \(y_3=R(x_1,x_2,y_1,y_2)\) with \(\deg R\le3\) gives a four-variable
> potential with constant Hessian determinant.

> **Theorem `HC4MYG4` -- quartic-graph obstruction.**  No graph
> \(y_3=R(x_1,x_2,y_1,y_2)\) with \(\deg R\le4\) gives a four-variable
> potential with constant Hessian determinant over a characteristic-zero
> field.

> **Theorem `HC4MYGJ1` -- plane normal-jet theorem.**  On the plane
> \(x_1=0\), the graph determinant is affine-linear with unit coefficient in
> the normal jet \(\partial_{x_1}R\).  In particular, every correction
> vanishing to second order on that plane is invisible there, and no graph
> \(R=R_{\le4}+x_1^2S\), with \(S\) of arbitrary degree, can have constant
> Hessian determinant.

> **Theorem `HC4MYG5J` -- quintic leading-jet rigidity.**  If a graph of
> degree at most five has constant Hessian determinant, its quintic trace on
> \(x_1=0\) is \(\kappa x_2^5\), its quartic trace is independent of
> \(y_2\), and the next determinant faces satisfy the explicit conic system
> in Section 5.  Thus only one of the 21 \(x_1\)-free quintic coefficients
> survives the leading faces.

> **Theorem `HC4MYG5S` -- sparse quintic-trace obstruction.**  For the v2
> potential, no degree-at-most-five graph with plane trace
> \(\kappa x_2^5+d x_2^3y_1+\rho x_2^2y_2\) has constant Hessian
> determinant.  This includes arbitrary \(x_1\)-linear normal jet and every
> remaining \(x_1^2\)-divisible term allowed by the degree bound.

None of the exclusions uses a collision equation.  For the v2 member, a graph
that contains the two marked points transfers their gradient collision
automatically.  Thus the determinant obstructions are strictly stronger than
their collision-preserving versions.

The full degree-five graph problem remains open.  The new normal-jet identity
explains both why degree five is the first unclassified single-graph degree
and why further restrictions to the same plane cannot close it by themselves.

## 1. The scaled Meng--Yang family

Work over a characteristic-zero field \(K\).  Put

\[
 x=x_1,\qquad y=x_2,\qquad p=y_1,\qquad q=y_2,
 \qquad r=y_3,\qquad u=1+xy,
\]

and define

\[
\begin{aligned}
 A={}&u^3p+3xu^2q-x^3r,\\
 B={}&y^2u(4+3xy)p
      +\bigl(y+3xy^2(4+3xy)\bigr)q
      +(2x-3x^2y)r.                                  \tag{1.1}
\end{aligned}
\]

The full scaled family is

\[
 \Psi_{L,M,N}=LA^2+MA+NB,
 \qquad LN\ne0.                                      \tag{1.2}
\]

The [Meng--Yang Schur-descent bridge](MENG_YANG_SCHUR_DESCENT_BRIDGE.md)
gives

\[
 \det\operatorname{Hess}\Psi_{L,M,N}=8LN^4.          \tag{1.3}
\]

Indeed, factor (1.2) as
\(N(B+(L/N)A^2+(M/N)A)\), use
\(\lambda=2L/N\) and bordered determinant \(c=-4\) in the bridge, and then
scale the five-variable Hessian by \(N\).  The parameter \(M\) does not
enter the determinant.

The v2 representative is

\[
 (L,M,N)=(1,13,2),
 \qquad \Psi=A^2+13A+2B,                              \tag{1.4}
\]

and its two marked points and common gradient are

\[
\begin{aligned}
 P_+&=(1,-3/2,0,0,0),&
 P_-&=(-1,3/2,0,0,0),\\
 \nabla\Psi(P_+)&=(0,0,-1/2,0,0)
                 =\nabla\Psi(P_-).                  \tag{1.5}
\end{aligned}
\]

## 2. Affine four-planes

Let \(H\subset\mathbb A^5\) be an affine hyperplane with normal

\[
 n=(\alpha,\beta,\gamma,\delta,\epsilon),             \tag{2.1}
\]

where \((\gamma,\delta,\epsilon)\) is its dual part.  If
\(i(w)=z_0+Tw\) parameterizes \(H\), then

\[
 \operatorname{Hess}(\Psi\circ i)
 =T^{\mathsf T}(\operatorname{Hess}\Psi\circ i)T.     \tag{2.2}
\]

This is equivalently the bordered-Hessian calculation.  Choose a vector
\(v\) with \(n^{\mathsf T}v=1\) and put \(Q=(T\ v)\).  For any symmetric
matrix \(G\),

\[
 \det(T^{\mathsf T}GT)
 =-\det(Q)^2
   \det\begin{pmatrix}G&n\\n^{\mathsf T}&0\end{pmatrix}. \tag{2.3}
\]

Thus changing affine coordinates on the same four-plane only composes the
restricted determinant with an affine automorphism and multiplies it by a
nonzero square.

### 2.1 Nonzero dual normal

Split the projective dual normal into three disjoint charts.  Normalize the
first nonzero dual component and solve the hyperplane equation for its
variable.  Direct coefficient extraction from (2.2) gives:

| dual-normal chart | eliminated coordinate | retained coordinates | unavoidable coefficient |
|---|---|---|---|
| \(\gamma\ne0\) | \(p=h-ax-by-cq-dr\) | \((x,y,q,r)\) | \([x^{12}q^4]=793152L^4\) |
| \(\gamma=0,\ \delta\ne0\) | \(q=h-ax-by-dr\) | \((x,y,p,r)\) | \([x^6p^4]=2160L^4\) |
| \(\gamma=\delta=0,\ \epsilon\ne0\) | \(r=h-ax-by\) | \((x,y,p,q)\) | \([x^2p^4]=2160L^4\) |

Every source-normal, remaining dual-normal, and translation parameter is
retained in these extractions.  The displayed coefficients are independent
of all of them, as well as of \(M,N\).  Since \(L\ne0\), the restricted
Hessian determinant is nonconstant on every chart.

### 2.2 Zero dual normal

If \((\gamma,\delta,\epsilon)=0\), the hyperplane cuts the two source
variables to one affine parameter \(s\) and retains all three dual
variables \(Y=(p,q,r)\).  The restriction has the form

\[
 f(s,Y)=L\bigl(a(s)^{\mathsf T}Y\bigr)^2+c(s)^{\mathsf T}Y. \tag{2.4}
\]

Its dual--dual Hessian block is

\[
 \operatorname{Hess}_{Y}f=2L\,a(s)a(s)^{\mathsf T}.   \tag{2.5}
\]

Consequently the three dual rows of the full Hessian lie in the span of the
single source column and the row \((0,a(s)^{\mathsf T})\).  Those three rows
have rank at most two; after adjoining the source row, the full rank is at
most three.  Hence

\[
 \det\operatorname{Hess}f\equiv0.                    \tag{2.6}
\]

Sections 2.1 and 2.2 exhaust every affine four-plane and prove `HC4MYA1`.

## 3. Polynomial graphs through degree three

Let

\[
 i_R(x,y,p,q)=(x,y,p,q,R(x,y,p,q)),
 \qquad \deg R\le3,                                  \tag{3.1}
\]

and put

\[
 D_R=\det\operatorname{Hess}(\Psi_{L,M,N}\circ i_R). \tag{3.2}
\]

Only three cubic coefficients are needed:

\[
 \rho=[y^2q]R,\qquad
 \sigma=[ypq]R,\qquad
 \tau=[p^2q]R.                                       \tag{3.3}
\]

Do not expand the generic determinant.  Restrict its Hessian entries to the
two lines

\[
 \ell_0(t)=(0,t,0,0),\qquad
 \ell_1(t)=(0,t,t,0),                                \tag{3.4}
\]

and extract only the indicated univariate coefficients.  Exact determinant
expansion gives

\[
\begin{aligned}
 [t^6]D_R(\ell_0(t))
   &=N^4(16\rho+89)^2,                                \tag{3.5}\\
 [t^6]D_R(\ell_1(t))
   &=N^4\bigl(16(\rho+\sigma+\tau)+89\bigr)^2.        \tag{3.6}
\end{aligned}
\]

If \(D_R\) is constant, (3.5) first forces

\[
 \rho=-\frac{89}{16}.                                \tag{3.7}
\]

Equation (3.6) then forces

\[
 \sigma=-\tau.                                       \tag{3.8}
\]

After only these two substitutions, the next line coefficient loses every
remaining coefficient of \(R\), as well as \(M\):

\[
 [t^5]D_R(\ell_1(t))
 =\frac{197}{4}LN^3.                                 \tag{3.9}
\]

This is nonzero because \(LN\ne0\), contradicting constancy.  This proves
`HC4MYG3`.  Notice that the proof excludes even a constant zero determinant;
the collision conditions never enter.

## 4. The quartic two-slope jet

Write \(R=R_0+\cdots+R_4\), with \(R_j\) homogeneous of degree \(j\), and
retain every coefficient of every piece.  Use the two-slope pencil

\[
 \ell_{c,d}(t)=(0,t,ct,dt).                         \tag{4.1}
\]

The highest relevant line coefficient is the square

\[
 [t^8]D_R(\ell_{c,d}(t))
 =256N^4\bigl(\partial_qR_4(0,1,c,d)\bigr)^2.       \tag{4.2}
\]

Thus constancy forces

\[
 \partial_qR_4(0,y,p,q)=0.                          \tag{4.3}
\]

This kills all ten \(x\)-free quartic coefficients containing \(q\),
including the proposed escape

\[
 [y^3q]R=[x_2^3y_2]R=0.                             \tag{4.4}
\]

For the next square, put

\[
\begin{gathered}
 \rho=[y^2q]R_3,\quad \sigma=[ypq]R_3,\quad
 \tau=[p^2q]R_3,\\
 a=[p^4]R_4,\quad b=[yp^3]R_4,\quad
 e=[y^2p^2]R_4,\quad f=[y^3p]R_4,\quad g=[y^4]R_4,\\
 \eta=[q^3]R_3,\quad \theta=[pq^2]R_3,\quad
 \iota=[yq^2]R_3.
\end{gathered}                                      \tag{4.5}
\]

After (4.3), exact truncated determinant expansion gives

\[
\begin{aligned}
 [t^6]D_R(\ell_{c,d}(t))=N^4\bigl(&-8c^3a
 +(16\tau-6b)c^2+32cd\theta +(16\sigma-4e)c\\
 &+48d^2\eta+32d\iota+16\rho-2f+89\bigr)^2.
                                                               \tag{4.6}
\end{aligned}
\]

It follows triangularly that

\[
 a=\eta=\theta=\iota=0,\qquad
 b=\frac83\tau,\qquad e=4\sigma,\qquad
 f=8\rho+\frac{89}{2}.                              \tag{4.7}
\]

The \(t^7\) coefficient vanishes after (4.3).  Substitute (4.7) into the
\(t^5\) coefficient.  Four coefficients of its polynomial in \(c\) are
enough:

\[
\begin{array}{c|c}
\text{coefficient, after earlier rows vanish}&\text{value}\\ \hline
[c^5]&-\dfrac{64}{3}LN^3\tau^2\\[2mm]
[c^3]\ \text{after }\tau=0&64LN^3\sigma^2\\[1mm]
[c]\ \text{after }\tau=\sigma=0
  &2LN^3(160\rho^2+1968\rho+6021)\\[1mm]
[1]\ \text{after }\tau=\sigma=0
  &4LN^3(16\rho+99)g.
\end{array}                                         \tag{4.8}
\]

Because \(LN\ne0\) and the field has characteristic zero, constancy forces

\[
 \tau=\sigma=g=0,\qquad
 Q(\rho):=160\rho^2+1968\rho+6021=0.                \tag{4.9}
\]

For the last implication, \(Q(-99/16)=-243/8\), so the last factor in
(4.8) cannot vanish through \(16\rho+99\).  Equations (4.3), (4.7), and
(4.9) leave the exact slice

\[
 R_4(0,y,p,q)=\left(8\rho+\frac{89}{2}\right)y^3p. \tag{4.10}
\]

The discriminant and roots of \(Q\) are

\[
 \operatorname{disc}(Q)=19584=576\cdot34,qquad
 \rho=\frac{3(-82\pm\sqrt{34})}{40}.                \tag{4.11}
\]

At this stage (4.11) leaves two apparent branches over
\(K(\sqrt{34})\).  The next coefficient eliminates both.

### 4.1 The terminal \(t^4\) coefficient

Set

\[
 \mu=[p^3]R_3,qquad \nu=[q^2]R_2.                  \tag{4.12}
\]

After only (4.3), (4.7), and (4.9), three coefficients of
\([t^4]D_R(\ell_{c,d}(t))\) are

\[
\begin{aligned}
 [c^4t^4]D_R(\ell_{c,d}(t))&=36N^4\mu^2,\\
 [d^2t^4]D_R(\ell_{c,d}(t))&=1024N^4\nu^2,\\
 [dt^4]D_R(\ell_{c,d}(t))\big|_{\mu=\nu=0}
   &=4LN^3P(\rho),                                  \tag{4.13}
\end{aligned}
\]

where

\[
 P(\rho)=8\rho^2+99\rho+279.                       \tag{4.14}
\]

Constancy first forces \(\mu=\nu=0\), and then \(P(\rho)=0\).  But

\[
 Q(\rho)-20P(\rho)=-3(4\rho-147),\qquad
 P(147/4)=\frac{58887}{4}\ne0.                     \tag{4.15}
\]

Thus \(P\) and \(Q\) have no common root in any characteristic-zero field
(their resultant is \(16959456\)).  This contradiction proves `HC4MYG4`.
It uses neither the collision equations nor an external classification
theorem.

### 4.2 Collision transfer and the remaining frontier

Suppose now that (1.4) is used and the graph contains the two marked points:

\[
 R(1,-3/2,0,0)=R(-1,3/2,0,0)=0.                     \tag{4.16}
\]

The graph chain rule is

\[
 \nabla(\Psi\circ i_R)
 =D i_R^{\mathsf T}\nabla\Psi.                       \tag{4.17}
\]

The two ambient gradients in (1.5) agree, and their omitted \(r\)-component
is zero.  Therefore the terms involving \(\nabla R\) vanish at both points,
and (4.17) transfers the collision without any equality between the two graph
tangent spaces.

Degree four is now closed even before (4.16) is imposed.  For single-graph
searches, the next line-jet problem is therefore \(\deg R=5\).  This does
not close nonlinear generating families that are not globally one graph,
nonconstant mixed source--dual Schur pivots, or higher-degree nonlinear
symplectic transformations.

This inherited-collision graph problem is separate from the strongest
parallel coordinate theorem route.  The complete coordinate
\(q_2+h_3+h_4+h_6\) chart is already closed by
[`HC4TC1`](HC4_MENG_TRIPLE_RANK_ONE.md), while the quintic layer
\(q_2+h_3+h_4+h_5+h_6\) reaches the unresolved Hessian-metric divisibility

\[
 \det(\bar C)\mid
 (\nabla s_4)^{\mathsf T}\operatorname{adj}(\bar C)\nabla s_4             \tag{4.18}
\]

in [`HC4CD5`](HC4_QUINTIC_COMMON_DIRECTION.md).  Repeating cubic--sextic
coordinate enumeration or the already excluded short Hamiltonian words does
not address either frontier.

## 5. The plane normal-jet equation

Write

\[
 T(y,p,q)=R(0,y,p,q),\qquad
 S(y,p,q)=\partial_xR(0,y,p,q).                     \tag{5.1}
\]

At \(x=0\), direct differentiation of (1.2) gives

\[
 \Psi_r=\Psi_{rr}=0,
 \qquad
 (\Psi_{xr},\Psi_{yr},\Psi_{pr},\Psi_{qr})
 =(2N,0,0,0).                                      \tag{5.2}
\]

Let \(w=(x,y,p,q)\), \(e_x=(1,0,0,0)^{\mathsf T}\), and

\[
 g=(S,T_y,T_p,T_q)^{\mathsf T}.
\]

The graph chain rule therefore reduces exactly to

\[
 \operatorname{Hess}(\Psi\circ i_R)|_{x=0}
 =\Psi_{ww}|_{x=0,r=T}
  +2N(e_xg^{\mathsf T}+ge_x^{\mathsf T}).           \tag{5.3}
\]

The lower-right \((y,p,q)\)-minor of (5.3) is independent of \(R\) and has
determinant

\[
 -2LN^2.                                            \tag{5.4}
\]

Only the \((x,x)\)-entry contains \(S\), with coefficient \(4N\).  Hence

\[
 \boxed{
 D_R(0,y,p,q)=\mathcal F(T,T_y,T_p,T_q)-8LN^3S,
 }                                                   \tag{5.5}
\]

where \(D_R=\det\operatorname{Hess}(\Psi\circ i_R)\) and \(\mathcal F\)
is the expression obtained by setting \(S=0\).  Since \(LN\ne0\), for every
trace \(T\) and every target constant \(C\), there is a unique polynomial
normal jet

\[
 S=\frac{\mathcal F(T,T_y,T_p,T_q)-C}{8LN^3}         \tag{5.6}
\]

that makes \(D_R=C\) on the plane, provided the right side satisfies the
desired degree bound.

This has two opposite consequences.

First, if a correction lies in \((x^2)\), then its value and first jet vanish
on \(x=0\).  Equation (5.3) is unchanged.  Therefore

\[
 D_{R_{\le4}+x^2U}|_{x=0}=D_{R_{\le4}}|_{x=0}.       \tag{5.7}
\]

The quartic contradiction in Section 4 occurs entirely on this plane, so
(5.7) excludes \(R_{\le4}+x^2U\) for arbitrary polynomial \(U\).  This proves
`HC4MYGJ1`.

Second, an \(x\)-linear quintic term is exactly \(xS_4(y,p,q)\), with
\(S_4\) quartic.  Once all determinant faces above degree four vanish,
(5.6) can cancel the complete degree-four residual.  This is the mechanism
that was unavailable through degree four.

### 5.1 The degree-five leading faces

Write the plane trace as

\[
 T=T_0+T_1+\cdots+T_5,
\]

with \(T_j\) homogeneous of degree \(j\) in \((y,p,q)\).  On the pencil
\((y,p,q)=(t,ct,dt)\), the top coefficient is

\[
 [t^{10}]D_R(0,t,ct,dt)
 =256N^4\bigl(\partial_qT_5(1,c,d)\bigr)^2.          \tag{5.8}
\]

Thus \(T_{5,q}=0\).  The \(t^9\) coefficient then vanishes, and the next
square is

\[
 [t^8]D_R(0,t,ct,dt)
 =4N^4\bigl(T_{5,p}-8yT_{4,q}\bigr)^2(1,c,d).        \tag{5.9}
\]

Put \(Q=T_{4,q}\).  Equation (5.9) makes \(Q\) independent of \(q\) and
gives \(T_{5,p}=8yQ\).  The intervening odd face is

\[
 [t^7]D_R(0,t,ct,dt)
 =16LN^3Q\bigl(T_{5,y}-4pQ\bigr)(1,c,d).             \tag{5.10}
\]

The polynomial ring is a domain.  If the first factor in (5.10) vanishes,
then \(Q=0\).  If the second factor vanishes, differentiate it in \(p\), use
\(T_{5,p}=8yQ\), and use that \(Q\) is a homogeneous cubic.  This gives

\[
 7Q-3pQ_p=0.                                        \tag{5.11}
\]

A monomial \(y^{3-k}p^k\) in \(Q\) would require \(7=3k\), which has no
integer solution.  Hence this branch also has \(Q=0\), and consequently

\[
 \boxed{T_5=\kappa y^5,\qquad T_{4,q}=0.}            \tag{5.12}
\]

More generally the same compatibility calculation has a resonance only
when the graph degree is divisible by three.  Degree five is nonresonant.

For the next faces set

\[
\begin{gathered}
 \rho=[y^2q]T_3,\qquad \sigma=[ypq]T_3,
 \qquad \tau=[p^2q]T_3,\\
 C=[y^2p^2]T_4,\qquad D=[y^3p]T_4,\\
 U=16\sigma-4C,\qquad
 V=16\rho-2D+89,\qquad
 \lambda=\frac{4L\kappa}{N}.
\end{gathered}                                      \tag{5.13}
\]

The highest slope terms of the \(t^6\) face first force

\[
 [q^3]T_3=[pq^2]T_3=[yq^2]T_3=[p^4]T_4=0,
 \qquad [yp^3]T_4=\frac83\tau.                     \tag{5.14}
\]

After (5.14), the complete residual face is

\[
 \frac{[t^6]D_R(0,t,ct,dt)}{N^4}
 =(Uc+V)^2
  +\lambda(20\tau c^2+20\sigma c+20\rho+123).      \tag{5.15}
\]

Two coefficients of the next face are triangular:

\[
 [c^5t^5]D_R=-\frac{64}{3}LN^3\tau^2,
 \qquad
 [c^3t^5]D_R\big|_{\tau=U=0}=64LN^3\sigma^2.       \tag{5.16}
\]

It follows that

\[
 \tau=U=\sigma=0,qquad
 \boxed{V^2+\lambda(20\rho+123)=0.}                \tag{5.17}
\]

Equations (5.12)--(5.17) prove `HC4MYG5J`.  They use no collision condition.
Of the 56 quintic coefficients of \(R_5\), the plane trace initially contains
21; only \(\kappa=[y^5]R_5\) survives.  The 15 coefficients in
\(x\,K[y,p,q]_4\) form the first normal repair, and the remaining 20
coefficients are divisible by \(x^2\).  By (5.7), a degree-five candidate
must have either \(\kappa\ne0\) or a nonzero \(x\)-linear quartic normal jet.

### 5.2 The first transverse variation

There is another useful unit coefficient.  Replace a plane-flat graph by

\[
 R\longmapsto R+x^2U.
\]

The plane Hessian is unchanged.  At first order in \(x\), only its
\((x,x)\)-entry changes, by \(12NU(0,y,p,q)\).  Multiplying by the cofactor
(5.4) gives

\[
 [x]D_{R+x^2U}-[x]D_R
 =-24LN^3U(0,y,p,q).                                \tag{5.18}
\]

For \(\deg R\le5\), this remaining normal jet has degree at most three.
Therefore every homogeneous component of \([x]D_R\) above degree three is
an invariant obstruction: the \(x^2U\) freedom cannot cancel it.

## 6. A sparse degree-five trace is impossible

Now specialize to the v2 parameters \((L,M,N)=(1,13,2)\), and assume the
complete plane trace is

\[
 T=\kappa y^5+d y^3p+\rho y^2q.                    \tag{6.1}
\]

All other coefficients of the graph remain free subject to total degree at
most five.  Put

\[
 V=16\rho-2d+89,\qquad \lambda=2\kappa,
 \qquad \mathcal Q(\rho)=160\rho^2+1968\rho+6021.  \tag{6.2}
\]

Equations (5.15) and (5.16) reduce to

\[
\begin{aligned}
 E_6&=V^2+\lambda(20\rho+123)=0,\\
 E_5&=\mathcal Q(\rho)-3V(4\rho+13)+30\lambda=0.
\end{aligned}                                      \tag{6.3}
\]

Eliminate \(\lambda\).  The residual plane curve is

\[
\begin{aligned}
 G(\rho,V)={}&30V^2+240V\rho^2+2256V\rho+4797V\\
 &-3200\rho^3-59040\rho^2-362484\rho-740583=0.     \tag{6.4}
\end{aligned}
\]

Modulo (6.4), \(\mathcal F(T)\) has degree four.  Therefore (5.6) gives the
unique normal jet

\[
 S=\frac{\mathcal F(T)-C}{64}                       \tag{6.5}
\]

for which the determinant equals the arbitrary target constant \(C\) on
\(x=0\).

Compute the first transverse determinant on the source axis \(p=q=0\),
reduce its coefficients modulo \(G\), and retain the two terms immune to
(5.18):

\[
 [x]D_R(x,y,0,0)\equiv
 \frac85H_7(\rho,V)y^7+
 \frac{104}{5}H_5(\rho,V)y^5+O(y^3),               \tag{6.6}
\]

where

\[
\begin{aligned}
H_7={}&1920V\rho^3+36240V\rho^2+223344V\rho+448443V\\
     &-25600\rho^4-650880\rho^3-6194880\rho^2
       -26158356\rho-41346207,\\
H_5={}&1212V\rho-7716V+51040\rho^2+627792\rho+1920699.
\end{aligned}                                      \tag{6.7}
\]

The exact eliminations are

\[
\begin{aligned}
 \operatorname{Res}_V(G,H_7)
   &=-120\,\mathcal Q(\rho)^2 A_4(\rho),\\
 \operatorname{Res}_V(H_7,H_5)
   &=15\,\mathcal Q(\rho) A_3(\rho),               \tag{6.8}
\end{aligned}
\]

with

\[
\begin{aligned}
A_4={}&3200\rho^4+72240\rho^3+579240\rho^2
       +1899693\rho+2000700,\\
A_3={}&53760\rho^3+858080\rho^2+4224396\rho+6004503.
\end{aligned}
\]

The monic gcd of the two resultants is exactly \(\mathcal Q\).  On that
remaining locus the exact Groebner basis is

\[
 (G,H_7,H_5,\mathcal Q)=(V,\mathcal Q).             \tag{6.9}
\]

Thus simultaneous transverse cancellation first forces

\[
 V=0,\qquad \mathcal Q(\rho)=0,
 \qquad \lambda=\kappa=0.                          \tag{6.10}
\]

This is precisely the old apparent quartic branch, now repaired on the plane
by the \(xS_4\) term.  It still fails transversely.  Reducing the full first
transverse polynomial modulo \(\mathcal Q\) gives

\[
 [xy^4q]D_R=216(5\rho+12).                          \tag{6.11}
\]

The correction (5.18) has degree at most three and cannot alter (6.11).
But

\[
 \mathcal Q(-12/5)=\frac{11097}{5}\ne0.            \tag{6.12}
\]

Equations (6.10)--(6.12) are a contradiction and prove `HC4MYG5S` over
every characteristic-zero field.

### 6.1 A collision-containing plane-flat near miss

The plane equations themselves genuinely have rational degree-five
solutions.  One particularly small point is

\[
 \rho=-\frac{123}{20},\qquad V=0,\qquad
 \lambda=\frac{51}{50},\qquad
 \kappa=\frac{51}{100},\qquad d=-\frac{47}{10}.     \tag{6.13}
\]

Thus

\[
 T_*=\frac{51}{100}y^5-\frac{47}{10}y^3p
     -\frac{123}{20}y^2q.                           \tag{6.14}
\]

Let \(S_*=(\mathcal F(T_*)-C)/64\).  At the two marked points, the graph
values are

\[
 R(P_+)=-\frac{25C-17165601}{1600},\qquad
 R(P_-)= \frac{25C-17165601}{1600}.                 \tag{6.15}
\]

Choosing

\[
 C=\frac{17165601}{25}                              \tag{6.16}
\]

makes \(R=T_*+xS_*\) contain both collision points and gives
\(D_R=C\) identically on \(x=0\).  The collision therefore transfers.  This
is not an `HC(4)` example: its first transverse coefficient already has

\[
 [xy^7]D_R=\frac{22032}{125}\ne0.                   \tag{6.17}
\]

This exact near miss is useful evidence that the plane gate is now solved in
the positive direction and that the real obstruction has moved to transverse
compatibility.

## 7. Consequences for the degree-five search

The full degree-five graph problem is not closed.  The surviving generic
shape is

\[
 R_5=\kappa y^5+xS_4(y,p,q)+x^2U_3(x,y,p,q),        \tag{7.1}
\]

with at least one of \(\kappa,S_4\) nonzero.  Lower trace terms must satisfy
(5.14)--(5.17), and the plane-flat equation fixes \(S_4\) rather than leaving
it as a search parameter.

Over \(\mathbb Q\), the complementary top-cone theorem HC4MYG5N in
[the quintic Schur-frontier note](HC4_MENG_YANG_QUINTIC_SCHUR_FRONTIER.md)
also forces \(\kappa\ne0\) and a constant kernel direction in the
\((p,q)\)-plane.  Its two charts are

\[
 R_5=\kappa y^5+xT_4(x,y,p-aq),
 \qquad
 R_5=\kappa y^5+xU_4(x,y,q).                        \tag{7.2}
\]

In particular, the quartic part of the uniquely forced normal jet from
(5.6) must satisfy respectively

\[
 (a\partial_p+\partial_q)S_4=0,
 \qquad\text{or}\qquad
 \partial_pS_4=0.                                  \tag{7.3}
\]

Thus the external cone theorem and the exact plane identity meet in a small
differential ideal before any full determinant is formed.

The next efficient calculation is therefore not another line on \(x=0\).
It is the joint ideal formed by:

1. the high plane faces (5.12)--(5.17);
2. the unique normal solution (5.6);
3. all homogeneous components of \([x]D_R\) above degree three, which are
   immune to \(x^2U_3\) by (5.18);
4. only after that elimination, the collision-containment equations and the
   remaining degree-at-most-three transverse repair.

The sparse theorem shows that this transverse ideal can close a family even
when the plane determinant and the inherited collision have both been solved
exactly.  The highest-value next calculation is to impose (7.3) on the
general trace satisfying (5.14)--(5.17), eliminate its quartic normal jet,
and only then enlarge the sparse trace by the surviving monomials.  This
intersects the plane, top-cone, and first-transverse gates rather than
repeating any one of them.

The first such intersection is `HC4MYG5K` in the quintic Schur-frontier
note.  On the \(\partial_q\)-kernel chart it retains the complete
degree-at-most-two lower trace.  Three immutable transverse equations generate
the generic unit ideal over \(\mathbb Q\), a quartic-normal square kills the
\(q^2\) coefficient, and the exceptional nonzero \(yq\) branch ends in two
coprime transverse polynomials.  Broader cubic and quartic traces on this
chart and the other kernel charts remain open.

## 8. Reproduction

Run

```bash
.venv/bin/python scripts/verify_hc4_meng_yang_graph_obstructions.py
```

The checker constructs all affine-normal parameters symbolically, extracts
the three chart coefficients without forming the full determinant, builds a
generic degree-at-most-three graph, verifies (3.5)--(3.9), and checks the v2
collision gradient.  It then retains every quartic coefficient, computes the
two-slope determinant only through degree eight, and verifies
(4.2), (4.6), (4.8), (4.11), and the terminal contradiction
(4.13)--(4.15).  Finally it derives the exact plane graph Hessian, checks the
unit normal coefficient (5.5), verifies the degree-five leading faces
(5.8)--(5.17), constructs the sparse plane-flat quotient, reproduces the
resultants (6.8), and checks the transverse contradiction (6.11)--(6.12) and
the collision-containing near miss (6.13)--(6.17).
