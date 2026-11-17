# Meng-Yang parent and normal-jet identities

These formulas define the fixed parent and the all-normal triangular
recurrence used by the current quintic graph proof. The complete bounded
graph exclusion is in the
[residual algebra proof](../jacobian/affine-complements/FULL_COMPLEX_QUINTIC_GRAPH_EXCLUSION_AND_CUBIC_NILPOTENCE.md).
Formal local solvability does not imply polynomial termination or full HC4.

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

### 5.3 All normal orders

The same coefficient is unit-triangular at every normal order.  Let
\(U\in K[y,p,q]\) and put \(\widetilde R=R+x^kU\), with \(k\ge1\).  Along
\(x=0\),

\[
 \Psi_r=0,\qquad \frac{d}{dx}\Psi_r(x,y,p,q,R)=\Psi_{xr}+\Psi_{rr}R_x=2N.
                                                               \tag{5.19}
\]

There is no hidden analytic expansion here.  From (1.1)--(1.2), exactly

\[
 \Psi_r=-2Lx^3A-Mx^3+N(2x-3x^2y),
 \qquad \Psi_{rr}=2Lx^6,
\]

and hence, with \(\delta=x^kU\),

\[
 \Psi(r+\delta)-\Psi(r)=\Psi_r(r)\delta+Lx^6\delta^2.
\]

Thus in the \(x\)-adic polynomial ring

\[
 \Psi(x,y,p,q,\widetilde R)-\Psi(x,y,p,q,R)
 =2Nx^{k+1}U+O(x^{k+2}).                            \tag{5.20}
\]

Its Hessian difference first appears in the \((x,x)\)-entry, where it is

\[
 2Nk(k+1)x^{k-1}U+O(x^k).                           \tag{5.21}
\]

Every other entry has \(x\)-order at least \(k\).  Multiplication by the
universal complementary cofactor (5.4) now gives

\[
 \boxed{
 D_{R+x^kU}-D_R
 \equiv-4LN^3k(k+1)x^{k-1}U\pmod{x^k}.
 }                                                   \tag{5.22}
\]

This proves `HC4MYGJ2`.  In characteristic zero with \(LN\ne0\), (5.22)
recursively solves the complete normal jet.  Equation (5.6) is the \(k=1\)
step; after coefficients through \(x^{k-2}\) have been killed, there is a
unique \(U\) that kills the coefficient of \(x^{k-1}\).  Therefore every
trace \(T\in K[y,p,q]\) and target constant \(C\) determine a unique

\[
 R=T+xS+x^2U_2+x^3U_3+\cdots\in K[y,p,q][[x]]       \tag{5.23}
\]

with \(D_R=C\).  The theorem supplies formal local existence and uniqueness;
it does not say that (5.23) terminates, converges away from the plane, or
contains the two marked points.  The all-degree graph problem is therefore a
polynomial-termination and global-collision problem, not a failure of normal
jet solvability.
