# Laplacian inverse certificates and closed safe loci

15 November 2026. This is an all-degree synthesis and a finite-certificate
deduction, not a novelty claim. The new low-dimensional quartic inputs at the
end are **source reports, not independently replayed witnesses**. Mathematical
state and programme archiving remain separate in the [registry](../../MATH_STATUS.json).

## Scope and sources

Work over a characteristic-zero field, extending scalars to its algebraic
closure when necessary. Put \(\Delta=\sum_{i=1}^n\partial_i^2\), and assume
\(P(0)=0,\nabla P(0)=0\). These normalization hypotheses are retained throughout.
“Pure” means the **polynomial identities** \(\Delta^mP^m=0\), not just their
values at the origin. “Safe” means that for every polynomial \(Q\),
\(\Delta^m(QP^m)=0\) eventually. The operator here is \(\Delta\) itself.

Retained primary reading:

* [Zhao, math/0409534v2](https://arxiv.org/html/math/0409534v2), Theorems
  3.4, 4.3, 6.2 and Section 7: inversion/heat identities, the Hessian criterion,
  and the inverse-degree bound used in Proposition 7.4.
* [van den Essen--Zhao, 0704.1690v1](https://arxiv.org/html/0704.1690v1),
  Theorem 1.5 and its complete proof: the fixed-potential arbitrary-multiplier
  equivalence. [de Bondt, 1206.2836v2](https://arxiv.org/html/1206.2836v2),
  Theorem 1 and Corollary 2, give the operator commutator formulation.
* [de Bondt, 2610.24889v1](https://arxiv.org/html/2610.24889v1), submitted
  17 October 2026: definitions, Sections 3.1--3.3 and 4.7, with selected base
  and cubic-reduction formulas. Its ancillary Maxima reductions were not run.

Gabber's inverse-degree theorem, published by Bass--Connell--Wright, is an
imported classical input: an automorphism of degree at most \(r\) in dimension
\(n\) has inverse degree at most \(r^{n-1}\). Its explicit use and statement in
Zhao's Proposition 7.4 were read; the original AMS PDF returned HTTP 403.
No claim of reading that PDF is made. Wang's degree-two theorem is likewise
imported through Zhao's explicit degree-three consequence.

## 1. A single inverse controls every multiplier

By Zhao's criterion, purity is equivalent to the Hessian \(A=\operatorname{Hess}P\)
being nilpotent, and is already equivalent to the first \(n\) pure identities.
This finite test is independent of \(\deg P\) in this **polynomial-identity,
ordinary-Laplacian** setting. It gives no degree-free Gaussian moment test.

Write
\[
 F_t=x-t\nabla P,\qquad G_t=F_t^{-1}=x+t\nabla U_t,
 \qquad U_t(0)=0.
\]
The inverse initially belongs to \(K[x][[t]]\). The normalized potential is
uniquely determined; for example
\[
 U_t=P(G_t)-\frac t2\langle\nabla P(G_t),\nabla P(G_t)\rangle.
\]
Zhao's identity is
\[
 U_t^k=k!\sum_{j\ge0}
 \frac{t^j\Delta^jP^{j+k}}{2^j j!(j+k)!}\quad(k\ge1).                 \tag{1}
\]
It follows that the following are equivalent for this fixed pure \(P\):

1. Every polynomial multiplier is eventually killed.
2. The single multiplier \(Q=P\) is eventually killed.
3. \(U_t\in K[x,t]\).
4. \(F_t\) is a polynomial automorphism over \(K[t]\).

Here (1) of the list implies (2); equation (1) with \(k=1\) equates (2) and
(3); the inverse formula and normalization equate (3) and (4). If
\(N=\deg_tU_t\), equation (1) gives
\[
 \Delta^jP^{j+k}=0\quad\text{for }j>kN.                             \tag{2}
\]
For a multiplier of degree \(d\ge1\), take \(j=m-d,k=d\). The right side of
(2) applies whenever \(m>(N+1)d\). The constant-coefficient commutator lemma
then gives
\[
 \boxed{\ \Delta^m(QP^m)=0\quad\text{for }m>(N+1)\deg Q.\ }        \tag{3}
\]
Indeed \(\Lambda^{m-d}f=0,\deg Q\le d\) implies
\(\Lambda^m(Qf)=0\): commuting through a polynomial lowers its degree on
each commutator, and after at most \(d\) commutators a remaining power of
\(\Lambda\) annihilates \(f\). For constant \(Q\), purity handles every
\(m\ge1\). Set \(N=0\) when \(U=0\). Bound (3) need not be sharp.

Thus a supplied polynomial inverse, checked by both composition identities,
is a finite certificate for **all** multipliers and all orders. No general
holonomic integration, period-kernel solver, or infinite marked replay is
required to verify such a certificate. Inverse production is a separate task.

## 2. A uniform bound, including nonhomogeneous normalized potentials

Fix \(n\ge2,D\ge2\), and allow all normalized pure potentials of degree at
most \(D\). Set
\[
 r=D-1,\quad B=r^{n-1},\quad L=(2n-1)B-n.                           \tag{4}
\]
Every **safe field point** has
\[
 \deg_xG_t\le B,\quad \deg_tG_t\le L,\quad\deg_tU_t\le L-1,
 \qquad \Delta^m(QP^m)=0\quad(m>L\deg Q).                          \tag{5}
\]
The first inequality is Gabber's theorem over \(K(t)\). The time bound does
not follow from it alone; the following recursion supplies the missing step.

Let \(H=\nabla P=A_0x+\sum_{j=2}^rH_j(x)\), with \(H_j\) homogeneous and
\(A_0^n=0\). Let
\[
 C_t=(I-tA_0)^{-1}=\sum_{i=0}^{n-1}t^iA_0^i.
\]
Writing the formal inverse by its homogeneous **x** degree gives
\[
 G_1=C_tx,\qquad
 G_k=tC_t\left[\sum_{j=2}^rH_j(G_1+\cdots+G_{k-1})\right]_k
 \quad(k\ge2).                                                    \tag{6}
\]
Each coefficient is a polynomial in \(t\) and the coefficients of \(P\).
Induction proves
\[
 \deg_tG_k\le(2n-1)k-n.                                           \tag{7}
\]
For \(k=1\) this is \(n-1\). In a term with \(j\ge2\) children of degrees
\(k_i\), the right side is bounded by
\(n+\sum_i((2n-1)k_i-n)\le(2n-1)k-n\).
Apply (7) for \(k\le B\). Since \(G-x=t\nabla U\), differentiation and
the normalization \(U(0)=0\) yield the last time-degree bound in (5).
Equation (3) completes the proof. In dimension one the pure normalized
potential is zero, so handle that case separately.

## 3. Closed safety and a finite family decision

The normalized pure coefficient locus is algebraic: impose the coefficients
of \(\Delta^mP^m\) for \(1\le m\le n\), or the characteristic coefficients
of \(\operatorname{Hess}P\). On that locus use (6) to compute
\[
 T=G_1+\cdots+G_B.
\]
Then impose every coefficient of the two finite polynomial identities
\[
 F_t(T)=x,\qquad T(F_t)=x.                                        \tag{8}
\]
These equations define **exactly the safe field points**. A safe point has
inverse degree at most \(B\), so its formal inverse equals \(T\). Conversely
(8) supplies a polynomial inverse and Section 1 supplies safety.

Consequently the safe locus in this whole normalized fixed-\(n,D\)
standard-Laplacian family is Zariski closed and has the common cutoff (5).
The pullback to any polynomial coefficient family is also closed. All safe
field points lie in the closed cutoff locus for the explicit integer \(L\),
and every point in that cutoff locus is safe, so the increasing closed cutoff
chain stabilizes there by \(L\). This is a finite polynomial decision and a
finite stratification certificate; its potentially enormous generic producer
has not been implemented here. It is a field/radical assertion. It neither
computes full ordinary moment ideals nor extends the inverse-degree theorem
to arbitrary nonreduced coefficient algebras.

For a homogeneous safe potential of degree \(d\ge3\), homogeneity improves
the bound to
\[
 N\le\frac{(d-1)^{n-1}-(d-1)}{d-2},\qquad
 m>\frac{(d-1)^{n-1}-1}{d-2}\deg Q.                              \tag{9}
\]
The \(t^j\) term of \(\nabla U\) has x degree
\(d-1+j(d-2)\); Gabber bounds that degree. This is the conditional-on-safety
bound, not an assertion that every homogeneous HN potential is safe.

## 4. The degree boundary and the new report

Wang's theorem makes every normalized pure potential of degree at most three
safe in every dimension: \(F_t\) has degree at most two over \(K(t)\), so its
inverse has x degree at most \(2^{n-1}\); recursion (6) makes its finitely many
x coefficients polynomial in t. For homogeneous cubics, (9) gives the common
cutoff \(m>(2^{n-1}-1)\deg Q\). Homogeneous quadratics follow from the finite
nilpotent geometric inverse, and degrees zero and one are elementary by total
degree. Therefore a homogeneous ordinary-Laplacian failure has degree at least
four, uniformly in dimension.

The October 2026 de Bondt paper reports nonhomogeneous degree-four failures
in dimension 34 and homogeneous degree-four failures in dimension 36; it
also reports degree-eight failures in dimensions 28 and 30 respectively.
Its input map dimensions for the quartic path are 17 and 18.
The polynomial configurations for this path are supplied through ancillary
Maxima files. These files have **not** been replayed here. Hence degree four
is the reported sharp boundary, with that external assurance boundary retained.

The universal transfer itself is short. Assume \(H(0)=0\). In independent split
coordinates \((u,v)\), let \(h=v\cdot H(u)\) and
\(\Delta_B=2\sum\partial_{u_i}\partial_{v_i}\). For the split metric
\(B=\left(\begin{smallmatrix}0&I\\I&0\end{smallmatrix}\right)\),
\[
 B\operatorname{Hess}h=
 \begin{pmatrix}JH&0\\\sum_i v_i\operatorname{Hess}H_i&JH^T\end{pmatrix}.
\]
If \(JH\) is nilpotent, this block matrix is nilpotent. If \(u-tH(u)\)
has a collision, the gradient map for \(h\) inherits it on \(v=0\), and
cannot have a polynomial inverse. Section 1 then forces VC failure with
the single multiplier \(h\). A complex linear isometry converts
\(\Delta_B\) to the standard Laplacian. Homogeneous \(H\) of degree \(d\)
gives homogeneous \(h\) of degree \(d+1\) in twice as many variables.
This argument validates the transfer's scope, but not the unexecuted input
producer. No reduced-ring specialization theorem from Section 3.2 is needed;
in particular its displayed A versus A-inverse normalization typo is not used.

## 5. Retained boundaries and verification

This closes field uniformity and closedness for the normalized
ordinary-Laplacian HN family in every fixed dimension and degree. General
Gaussian safe-locus uniformity remains open. The Dvorsky five-variable
operator is still order three. The retained ternary \(\Delta^6\), degree-12
witness and the search for lower powers are unchanged. Doubling to 34 or 36
variables does not compress a witness into three variables. Ordinary-Laplacian
small-dimension endpoints and minimum quartic witness dimension remain open.

The [checker](verify_laplacian_inverse_certificates.py) checks exact small
polynomial inverse/Hamilton--Jacobi identities, both compositions, pure and
shifted contractions, multiplier tails, and a truncated false inverse rejection.
It also verifies the displayed three-variable base-map determinant and collision
from de Bondt's Section 2.4.1; that is not the 17/18-variable reduction replay.
See the [receipt](laplacian_inverse_verification.json) and
[source audit](laplacian_inverse_source_audit.json).
These finite controls do not prove the universal recursion, Gabber or Wang,
or certify any reported high-dimensional counterexample. The written synthesis
and new cutoff/closed-locus deduction still require independent review.
