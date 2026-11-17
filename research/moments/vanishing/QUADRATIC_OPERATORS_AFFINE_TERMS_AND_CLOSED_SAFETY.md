# Quadratic operators, affine terms and closed safety

15 November 2026. This extends the [normalized Laplacian result](LAPLACIAN_INVERSE_CERTIFICATES_AND_CLOSED_SAFE_LOCI.md).
It is an internal written proof with exact small controls, not an independently
reviewed theorem or a claim of novelty. The [registry](../../MATH_STATUS.json)
is authoritative. Previous source-bound proof notes are retained unchanged.

## Scope

Let \(K\) be a characteristic-zero field, \(A\) any constant symmetric
\(n\times n\) matrix, and
\[
 \Lambda_A=\sum_{i,j}A_{ij}\partial_i\partial_j.
\]
The matrix may be singular or zero, and may vary as part of the coefficient
family. Let \(P\in K[x]\) be arbitrary, including constant and linear terms.
Pure identities mean \(\Lambda_A^mP^m=0\) **as polynomials**, not merely at
the origin. Safety means eventual vanishing for every polynomial multiplier.
No first-order or constant term is added to the differential operator.

Classical inputs are the inverse-degree bound of Gabber and Wang's quadratic
map theorem, retained with their explicit import boundary in the preceding
note. The commutator lemma is [de Bondt1206.2836v2, Theorem1](https://arxiv.org/html/1206.2836v2).
The new formal argument below removes the normalization and invertibility
restrictions instead of silently removing hypotheses from Zhao's theorems.
For comparison, [Zhao0704.1689v2, Proposition1.2 and Theorem1.3](https://arxiv.org/html/0704.1689v2)
explicitly retain order at least two. Its arXiv version is 25 February 2008;
the HTML body date of August 2026 is not evidence of a new version.

## 1. A Laurent heat logarithm without normalization

For every polynomial \(P\), the t-adic map
\[
 F_t=x-tA\nabla P
\]
has a unique inverse \(G_t=x+O(t)\) with coefficients in \(K[x]\). Define
\[
 U_t=P(G_t)-\frac t2\nabla P(G_t)^T A\nabla P(G_t).                 \tag{1}
\]
Differentiating the inverse identity gives, without inverting \(A\),
\[
 \nabla U_t=\nabla P(G_t),\quad
 G_t=x+tA\nabla U_t,\quad
 \partial_tU_t=\tfrac12\nabla U_t^TA\nabla U_t,\quad U_0=P.       \tag{2}
\]
The first equality follows from
\(DG_t=(I-tA\operatorname{Hess}P(G_t))^{-1}\) and the cancellation of
\(I-t\operatorname{Hess}P(G_t)A\) in the derivative of (1).

Introduce the heat series
\[
 H(t,s)=\exp\!\left(\frac{t\Lambda_A}{2s}\right)e^{sP},\qquad
 L(t,s)=sP+\log(e^{-sP}H(t,s)).                                  \tag{3}
\]
At each fixed t order the second factor in the logarithm is a finite Laurent
polynomial in \(s\). Indeed conjugating \(\Lambda_A\) by \(e^{sP}\)
is a differential operator polynomial in \(s\), of degree at most two.
Thus (3) is a well-defined t-adic logarithm with
\(L\in K[x][s,s^{-1}][[t]]\), initially \(L(0,s)=sP\). Its equation is
\[
 \partial_tL=\frac{\Lambda_A L+\nabla L^TA\nabla L}{2s}.           \tag{4}
\]
**Every coefficient of L has s degree at most one.** Induct on t order:
\(\Lambda_A L/s\) has s degree at most zero, and
\(\nabla L^TA\nabla L/s\) has degree at most one. The coefficient of
\(s^1\) consequently solves the Hamilton--Jacobi equation in (2), and is
exactly \(U_t\). No connected-graph assertion or contour comparison is needed.

If the first \(n\) pure identities vanish, the coefficient
\(\Lambda_A^jP^r\) is zero whenever \(1\le r\le j\le n\): apply further
powers to \(\Lambda_A^rP^r=0\). Hence, modulo \(t^{n+1}\),
\(H\) has constant s coefficient one and no negative s powers. Taking its
logarithm in s agrees with (3), and (4)'s degree bound gives
\[
 L=sU_t\pmod{t^{n+1}},\qquad \Lambda_AU_t=0\pmod{t^n}.           \tag{5}
\]
Writing \(J=A\operatorname{Hess}P\), (2) gives the exact trace identity
\[
 (\Lambda_AU_t)(F_t(x))
   =\operatorname{Tr}\big(J(I-tJ)^{-1}\big)
   =\sum_{j\ge0}t^j\operatorname{Tr}(J^{j+1}).                   \tag{6}
\]
Thus (5) forces the first n traces to vanish. Newton identities and
Cayley--Hamilton imply \(J^n=0\) over \(K[x]\).

Conversely, if \(J\) is nilpotent, (6) makes \(U_t\) harmonic for
\(\Lambda_A\). Then \(e^{sU_t}\) solves the heat equation of (3), and
t-adic uniqueness gives \(H=e^{sU_t}\). Its constant s coefficient is one,
proving every pure identity. We have proved, with all affine terms retained,
\[
 \boxed{\ \Lambda_A^mP^m=0\ (1\le m\le n)
 \iff (A\operatorname{Hess}P)^n=0
 \iff\Lambda_A^mP^m=0\ (m\ge1).\ }                             \tag{7}
\]
These are polynomial identities. They supply no degree-free Gaussian origin
moment test. They are field statements, not a full ordinary coefficient-ideal
identification over nonreduced algebras.

## 2. Full inverse criterion, including singular A

For a pure pair, comparing positive powers of s in \(H=e^{sU_t}\) gives
\[
 U_t^k=k!\sum_{j\ge0}
 \frac{t^j\Lambda_A^jP^{j+k}}{2^j j!(j+k)!},\qquad k\ge1.        \tag{8}
\]
Therefore safety is equivalent to the single multiplier \(Q=P\) having an
eventually zero tail, and to \(U_t\) being polynomial in \(x,t\). It is also
equivalent to \(F_t\) having a polynomial inverse over \(K[t]\).

One direction uses \(G_t=x+tA\nabla U_t\). The reverse direction uses the
**full scalar formula (1)**, which is polynomial whenever \(G_t\) is.
For singular A, recovering only \(A\nabla U_t\) does not determine \(U_t\).
For example, with \(A=\operatorname{diag}(1,0)\) and \(P=xy\),
\[
 G_t=(x+ty,y),\quad U_t=xy+ty^2/2,
\]
while \(A\nabla U_t=(y,0)\) misses the time-dependent kernel potential.

As in the preceding proof, if \(N=\deg_tU_t\) (take N=0 when U=0), equation
(8) and the commutator lemma imply the finite-certificate bound
\[
 \Lambda_A^m(QP^m)=0\quad\text{whenever }m>(N+1)\deg Q.           \tag{9}
\]
For constant Q, purity handles all positive m. All these statements work for
zero A as well. No inverse of A appears in the certificate or its verification.

## 3. Adding an affine potential preserves purity and safety exactly

Write \(P=P_*+c^Tx+c_0\), where \(P_*(0)=0,\nabla P_*(0)=0\).
The Hessian is unchanged. If \(G_*,U_*\) are the inverse and potential for
\(P_*\), direct substitution gives
\[
 G_t(x)=G_*(x+tAc),\qquad
 U_t(x)=U_*(x+tAc)+c^Tx+c_0+\tfrac t2c^TAc.                     \tag{10}
\]
The potential equality follows from (2), its initial value and uniqueness.
Equations (7), (8), (10) prove invariance of both purity and safety in both
directions. Formula (10), including the last scalar term, is essential: a
gradient-only argument could incorrectly erase it. Safety need not have the
same numerical cutoff before and after the affine addition.

## 4. Closed safety while the quadratic operator varies in rank

Fix \(n\ge1,D\ge2\), let every symmetric A vary, and let every polynomial
P of degree at most D vary. Put
\[
 B=(D-1)^{n-1},\qquad C(n,D)=2nB-n+1.                            \tag{11}
\]
On the pure locus, normalize P to \(P_*\). Set \(V=A\nabla P_*\),
\(M_0=JV(0)\), and
\[
 E_t=(I-tM_0)^{-1}=\sum_{i=0}^{n-1}t^iM_0^i.
\]
Compute the formal inverse of \(x-tV\), by homogeneous x degree, up to B:
\[
 T_1=E_tx,\quad
 T_k=tE_t\left[\sum_{j\ge2}V_j(T_1+\cdots+T_{k-1})\right]_k,
 \quad T=\sum_{k=1}^BT_k.                                       \tag{12}
\]
Every coefficient is polynomial in A and P's coefficients; (7) guarantees
\(M_0^n=0\) at every pure field point, even on rank-drop strata.
Impose all coefficients of
\[
 (x-tV)(T)=x,\qquad T(x-tV)=x.                                  \tag{13}
\]
These finite equations define exactly the safe field pairs (A,P). If safe,
Gabber bounds the inverse x degree by B over K(t), so its formal expansion is
T. Conversely (13) supplies an inverse and Section 2 proves safety. Thus the
safe locus is Zariski closed **across all matrix ranks**, not merely separately
on each rank stratum. Every polynomial family pulls it back to a closed set.
This is a field/radical conclusion, not an assertion that these equations give
the full ordinary ideal of every cutoff condition.

For the explicit common cutoff, fix a field point of matrix rank q. A linear
change over the algebraic closure puts the operator into a q-variable
Laplacian with the remaining variables passive. Treat passive variables as
parameters and use their fraction field. For q>=2, the preceding normalized
proof gives, with \(B_q=(D-1)^{q-1}\) and \(L_q=(2q-1)B_q-q\),
\[
 \deg_tU_*\le L_q-1,\qquad \deg_xU_*\le B_q+1.
\]
Here normalization is in the active variables, so its affine coefficients may
belong to the passive fraction field. Formula (10) bounds the full potential's
t degree by \(L_q+B_q=2qB_q-q\). Hence its cutoff is at most
\(2qB_q-q+1\). The q=1 case has a phase affine in its active variable and
obeys m>deg Q directly; q=0 is zero. The expression in (11) increases with q,
so all safe pairs obey the single uniform bound
\[
 \boxed{\ \Lambda_A^m(QP^m)=0\quad(m>C(n,D)\deg Q).\ }           \tag{14}
\]
Clearing the passive fraction field and undoing the linear change yields the
original polynomial identities. No rank-dependent change is used to define
the closed equations (13). Consequently the increasing closed field cutoff
loci stabilize by C(n,D), even when A varies. For D<=1 elementary total degree
gives the common bound m>deg Q. This is not a bound for a general operator of
arbitrary order, or for Gaussian origin moments.

## 5. Cubics and an obstruction to quadraticization

If deg P<=3 and the pair is pure, \(x-tA\nabla P_*\) is a Keller map of
degree at most two over K(t). Wang's theorem makes it invertible there.
Gabber bounds its inverse degree by B. Recursion (12) makes every one of the
finitely many inverse coefficients polynomial in t. Therefore **every pure
cubic pair is safe for every homogeneous quadratic operator, in every
dimension, including affine phases and singular operators**.

In particular, no valid eventual-vanishing transfer can send the Dvorsky
order-three cubic obstruction to a polynomial phase of degree at most three
for any homogeneous quadratic operator, regardless of how many auxiliary
variables are introduced or which polynomial multiplier is used. A successful
quadraticization must increase phase degree or leave these hypotheses.
This complements DVG4: its canonical one-pair obstruction holds in every phase
degree, whereas the present obstruction allows all auxiliary dimensions but
only phase degree<=3. Neither result excludes a higher-degree, different
quadratic completion. Protected DVG1/DVG4 and minimum-dimension records are
unchanged. The reported quartic endpoints retain their unreplayed-input status.

## Verification and remaining work

The [checker](verify_quadratic_operator_safety.py) checks affine translations,
full scalar potentials, both inverse compositions, exact pure/shifted heat
coefficients, a degenerate kernel-potential mutation, a polynomial rank-drop
family and non-HN harmonic controls. The [receipt](quadratic_operator_safety_verification.json)
and [source audit](quadratic_operator_safety_source_audit.json) retain exact
versions and the finite-check boundary. The general equations (12)--(13) are a
written algorithm, not an executed universal producer. Independent review of
the Laurent-logarithm truncation and rank-uniform proof is outstanding.
General Gaussian closedness/uniformity, operators with lower-order terms or
order>=3, ordinary nilpotents, and smaller genuine quartic witnesses remain
separate obligations. The ternary Delta-power frontier is preserved.
