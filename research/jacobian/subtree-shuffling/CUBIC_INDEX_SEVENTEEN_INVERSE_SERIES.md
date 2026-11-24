# A closed inverse slice and exact degrees at cubic Engel index seventeen

Let \(K\) be a field of characteristic zero. The fixed homogeneous cubic
map \(H:K^{18}\to K^{18}\) in the
[literal index-seventeen input](cubic_index_seventeen_shuffle_failure_inputs.json)
has \((JH)^{17}=0\) and \((JH)^{16}\ne0\). Write
\[
 F_t=I+tH,\qquad G_t=F_t^{-1}=\sum_{k\ge0}t^kG_k
\]
for its formal inverse. Each nonzero \(G_k\) is homogeneous of ordinary
degree \(2k+1\).

The previous [jet construction](CUBIC_INDEX_SEVENTEEN_SHUFFLE_FAILURE.md)
used nontermination to obtain a nonzero coefficient in every interval
\([k,3k-2]\). Here a single explicit inverse slice proves \(G_k\ne0\)
for **every** \(k\). It also gives a smaller jet family with exact inverse
degree, rather than only a lower bound.

**Theorem ClosedAlgebraicInverseSlice.** In the coordinate order
\((x_1,\ldots,x_{12},w_1,\ldots,w_5,\tau)\), put
\[
 v=(1,1,0,\ldots,0,1),\qquad
 A(t)=\frac{(1-4t^2)^{-1/2}-1}{2t^2}
     =\sum_{j\ge0}c_jt^{2j},\qquad
 c_j=\frac12\binom{2j+2}{j+1}.
 \tag{1}
\]
The square root is the formal branch with constant term one, and the
quotient defining \(A\) is a power series with constant term one. Then
\[
 (G_t(v))_{x_1}=A(t),\qquad
 (G_t(v))_{x_6}=-tA(t)^2.
 \tag{2}
\]
Consequently, for every \(j\ge0\),
\[
 (G_{2j}(v))_{x_1}=c_j\ne0,\qquad
 (G_{2j+1}(v))_{x_6}
   =-\sum_{i=0}^{j}c_ic_{j-i}\ne0.
 \tag{3}
\]
These are nonzero rational numbers, so the nonvanishing persists over
every characteristic-zero extension.

**Corollary EveryIndexSeventeenWeightedCoefficientSurvives.** In the free
one-generator symmetric ternary Engel algebra \(A_{3,p}\), the weighted
inverse class \(U_k\) is nonzero for every \(k\ge0\) and every \(p\ge17\).
Thus constant-vector membership in the corresponding length-\(p\)
shuffle span fails at each individual grade, not only arbitrarily far
along the sequence.

This corollary concerns weighted inverse sums. It is stronger than the
existence of a nonzero individual tree, and uses a different argument
from the classical index-four one-generator obstruction. It does not
decide the weighted inverse question for \(4\le p\le16\).

## 1. Derivation of the inverse slice

The mixed cubic seed is \(K(X)=X+Q(X)+C(X)\), with \(C=\iota C_*\)
supported in its first five outputs. Its homogeneous companion is
\[
 H(X,W,\tau)=(\tau Q(X)-\tau^2\iota W,\ C_*(X),\ 0).
 \tag{4}
\]
The companion and its exact index are defined by the retained literal
input; no choice of inverse branch is needed to specify \(H\).

For target \((Y,0,1)\), the equation \(F_t(X,W,1)=(Y,0,1)\) gives
\[
 W=-tC_*(X),\qquad X+tQ(X)+t^2C(X)=Y.
 \tag{5}
\]
By homogeneity of \(Q,C\), the latter equation is \(K(tX)=tY\).
Thus \(X=t^{-1}K^{-1}(tY)\) as a formal series. All components of
\(K^{-1}(tY)\) have zero constant term, so this division introduces
no negative powers.

We use the three-variable base map from de Bondt,
arXiv:2610.24889v1, equations (6)--(8) and Section 4.6:
\[
 \begin{aligned}
 g_1&=a-3a^2b-a^3c,\\
 g_2&=3V^2ac+(9V^2-6V-2)b,\\
 g_3&=V^3c+(6V^2+2V)b^2,\qquad V=1+2ab.
 \end{aligned}
 \tag{6}
\]
The nine helper outputs in that cubic reduction are
\[
 \begin{aligned}
 f_4&=x_4+a^2c+3ab,& f_5&=x_5+b^3,& f_6&=x_6+a^2,\\
 f_7&=x_7+ac,& f_8&=x_8+b^2,& f_9&=x_9+ax_5,\\
 f_{10}&=x_{10}+ax_8,&
 f_{11}&=x_{11}+2ax_9,&
 f_{12}&=x_{12}+ab+ax_{10}.
 \end{aligned}
 \tag{7}
\]
The coefficient two in \(f_{11}\) is the retained rationally conjugate
scale. Setting these nine outputs to zero leaves the first three outputs
equal to (6): all correction products contain a vanishing helper output.

On the formal branch of \(g_3=0\) at the origin, \(V\) is a unit and
\[
 c=-\frac{(6V+2)b^2}{V^2}.
 \tag{8}
\]
Substitution into the other two coordinates yields, exactly,
\[
 g_2=b,\qquad g_1=\frac{a(1+ab)}{(1+2ab)^2}.
 \tag{9}
\]
For target \((t,t,0)\), set \(b=t\) and \(R=1+2at\). Equation (9) gives
\[
 t^2=\frac{R^2-1}{4R^2},\qquad
 R=(1-4t^2)^{-1/2},\qquad a=\frac{R-1}{2t}.
 \tag{10}
\]
The formal branch is uniquely fixed by \(R(0)=1\). Equations (7) give
\(x_6=-a^2\). Dividing these two coordinates by \(t\) in (5) proves (2).
The binomial expansion proves (1); every \(c_j\) is a positive rational
number. Hence the sums in (3) are nonzero rational numbers. This proves
the theorem without an asymptotic or finite-prefix argument.

For the corollary, polarize \(-H\). Its multiplication operator is
\(-JH(a)/3\), so the resulting algebra satisfies every Engel identity
of length \(p\ge17\). Evaluation of the universal inverse coefficient
at \(v\) gives \(G_k(v)\), with the same fixed sign convention.
Equation (3) therefore witnesses \(U_k\ne0\) at every grade. The
[shuffle-orbit dictionary](SUBTREE_SHUFFLE_AND_INVERSE_DEGREE.md)
then gives the stated failure of membership.

## 2. A compressed triangular jet family

Keep four independent variables \(u_1,u_2,u_3,\tau\). For each
\(1\le j\le M\), introduce a vector \(z_j\in K^{17}\). Define
\[
 Y(t)=
 (u_1,u_2,u_3,0,\ldots,0,\tau)
 +\sum_{j=1}^{M}(z_j,0)t^j.
 \tag{11}
\]
There are \(17M+4\) coordinates. Define \(\mathcal F_M\) by leaving
\(u_1,u_2,u_3,\tau\) fixed and setting
\[
 (\mathcal F_M)_j
 =z_j+[t^{j-1}]H_{\mathrm{first}\ 17}(Y(t)),
 \qquad 1\le j\le M.
 \tag{12}
\]
In particular \(\tau\) remains a genuine independent variable. Only its
positive-order jet coordinates are omitted.

**Theorem SharpCompressedIndexSeventeenJets.** For every \(M\ge0\),
\(\mathcal F_M\) is a triangular polynomial automorphism in \(17M+4\)
variables and
\[
 \deg\mathcal F_M^{-1}=2M+1.
 \tag{13}
\]
For \(M\ge1\), its nonlinear part is homogeneous cubic. For \(M\ge16\),
its nonlinear Jacobian has exact nilpotency index seventeen.

### Triangularity and the upper degree bound

The coefficient in (12) uses only \(z_1,\ldots,z_{j-1}\) and the four
initial variables. Thus recursive substitution supplies a polynomial
inverse. Give the initial variables weight one and each coordinate of
\(z_j\) weight \(2j+1\). A term contributing to level \(j\) has three
input levels summing to \(j-1\), hence the same weight \(2j+1\).
Every inverse coordinate retains its weight, and all input weights are
positive. Its ordinary degree is therefore at most \(2M+1\).

### The exact lower bound

Set the initial target to \(u_1=u_2=\tau=1,\ u_3=0\), and all positive
target levels to zero. Its inverse series is \(G_t(v)\) modulo
\(t^{M+1}\). Equation (3) gives a nonzero coefficient in level \(M\):
use \(x_1\) when \(M\) is even and \(x_6\) when \(M\) is odd.
Scaling the three nonzero initial target variables together by \(s\)
multiplies this coefficient by \(s^{2M+1}\), since \(G_M\) has that
ordinary degree. The inverse therefore contains a term of degree
\(2M+1\). This proves (13), including \(M=0\).

### The exact Jacobian index

Regard a tangent vector as a truncated series of the same form as (11).
The derivative of the nonlinear part of (12) is the operator
\[
 \delta Y(t)\longmapsto t\,JH(Y(t))\,\delta Y(t)
                  \pmod {t^{M+1}}.
 \tag{14}
\]
Its output has zero constant term and zero \(\tau\)-coordinate, so it
lies in the permitted tangent subspace. Repeated application is ordinary
matrix multiplication in \(K[\mathbf u,\mathbf z,t]/(t^{M+1})\);
the coefficient matrix is not differentiated again. The retained
identity \((JH)^{17}=0\) therefore proves seventeenth-power vanishing.

For \(M\ge16\), take input tangent \(d\tau\) in level zero and inspect
output coordinate \(x_4\) in level sixteen. Its coefficient in the
sixteenth power of (14) is
\[
 ((JH(Y(0)))^{16})_{x_4,\tau}
   =288u_1^{10}u_2^4u_3^3\tau^{15}\ne0.
 \tag{15}
\]
This is exactly the retained seed witness after the allowed constant
coordinates are set to zero. Its inputs all remain independent in
(11). Thus the index is seventeen. Freezing \(\tau\) itself would
remove this witness; that different restriction is not used.

For every fixed \(p>17\), direct sum with the cubic chain on
\(q_0,\ldots,q_{p-1}\),
\[
 T_0=0,\quad T_1=q_0^3,\quad
 T_j=q_0^2q_{j-1}\quad(2\le j<p),
\]
gives exact index \(p\) for the nonlinear part of \(I+T\).
The inverse of \(I+T\) has degree \(2p-1\). Consequently,
for every \(M\ge p-1\), the direct sum has dimension \(17M+4+p\) and
inverse degree exactly \(2M+1\).

No minimum dimension or least failing Engel length is asserted.
The construction strengthens the explicit upper part of the boundary;
the unresolved interval \(4,\ldots,16\) remains unchanged.

## 3. Exact verification and primary context

The source map is Alpöge's, in de Bondt's normalization and cubic
stabilization. The source statements inspected here are
[arXiv:2610.24889v1](https://arxiv.org/html/2610.24889v1),
equations (6)--(8) and Sections 4.6--4.7. The present claims are the
explicit slice, all-grade weighted nonvanishing, and the exact degrees
of the compressed jet family. No historical-priority claim is made.

The standard-library checker uses the pinned literal 45-term cubic map.
It verifies all eighteen inverse equations as polynomial identities
modulo
\[
 (1-4t^2)A(1+t^2A)-1=0,
 \tag{16}
\]
not by sampling finitely many values of \(t\). It also replays the full
seed Jacobian powers, verifies the restricted witness (15), and checks
the compressed jet construction against its coefficientwise derivative
at several small orders. A separate direct recurrence supplies finite
controls for the displayed binomial and convolution coefficients.
Incorrect branch relations, omitted helpers and perturbed cubic data
are rejected.

Run from the repository root:

    python3 research/jacobian/subtree-shuffling/verify_cubic_index_seventeen_inverse_series.py \
      --output /tmp/cubic_index_seventeen_inverse_series_verification.json

The accompanying source ledger pins the proof, checker, literal seed
and reused exact arithmetic. The receipt records the actual replay.
Universal conclusions follow from the written identities, induction and
operator argument above. Internal checks and review do not assert
formal verification or external review.
