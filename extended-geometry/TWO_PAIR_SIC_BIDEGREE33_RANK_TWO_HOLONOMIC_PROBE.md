# Holonomic probe on the rank-two bidegree-\((3,3)\) stratum

## 1. Outcome and status

There is now a concrete recurrence target for the first rank-two
bidegree-\((3,3)\) experiment.  Put

\[
 C=UW,\qquad
 U\in\operatorname {Mat}_{4\times2},\quad
 W\in\operatorname {Mat}_{2\times4},
 \tag{1.1}
\]

and normalize the moments by

\[
 \nu_m(C)=\frac{\mu_m(C)}{(3m+1)!}.                       \tag{1.2}
\]

At two unrelated integral exact-rank-two points and over each of the three
primes \(1000003,1000033,1000037\), 501 exact values of \(\nu_m\) admit a
scalar recurrence of order \(27\) whose coefficients have degree at most
11 in \(m\).  Each fit uses 335 equations and passes the remaining 139
equations.  At order 27, a degree-10 ansatz fails at both points modulo
\(1000003\).

A second calculation identifies the natural relative-cohomology rank.
The logarithmic Jacobian quotient has exact length \(18\), split as two
length-two endpoint pieces and a length-fourteen interior piece.  At the
same six point/prime combinations, an order-18 recurrence of
\(m\)-degree 18 passes 83 unused equations; degree 17 fails at both
points modulo \(1000003\).

Exact shift-Ore comparison sharpens this picture.  The order-27 operator
is not a left multiple of the order-18 operator.  Instead, at all six
samples,

\[
 R_{18}=Q_4G_{14},\qquad R_{27}=Q_{13}G_{14},             \tag{1.3}
\]

where \(G_{14}\) is their greatest common right divisor and the subscripts
denote shift order.  After clearing rational content, every coefficient of
\(G_{14}\) has \(m\)-degree 58, and the resulting operator annihilates all
487 available moment rows.  Its order 14 matches the interior length in
the exact \(2+2+14\) logarithmic-Jacobian decomposition.  At the first
point modulo \(1000003\), Section 4.4.2 upgrades these 487 rows to an
exact all-order divergence certificate with a separately verified zero
endpoint trace.  The other five samples remain bounded modular evidence.
At the first characteristic-zero point, multiplication by \(P\) is now
proved cyclic on the degree-14 interior critical algebra.  The sampled
order-14 operator also has a stable exact characteristic-zero
reconstruction and its leading descent identity closes over \(\mathbb Q\);
see Section 4.4.3.  However, the raw period is nonzero on both endpoint
idempotent pairs, so it does not descend through that ordinary Jacobian
algebra.  The \(m\)-dependent
divergence connection is essential and is now explicit at one modular
fiber; no characteristic-zero or universal factor is identified.

After making the coefficient of \(m^{11}\nu_{m+27}\) monic, all six fits
have the same forward coefficient:

\[
\boxed{
 L(m)=
 \frac{
 \displaystyle
 \prod_{k\in\{71,73,74,76,77,79,80,82\}}(3m+k)
 \left(9m^3-1836m^2+210730m-17920380\right)}
 {9\cdot3^8}.}
 \tag{1.4}
\]

This is strong modular evidence for a compact universal recurrence shape.
It is **not** a creative-telescoping certificate.  The other recurrence
coefficients have not been reconstructed in the rank-two parameter ring,
and their universal common parameter denominator has not been found.  A
single generic factor pencil now has an exact modular border-basis
denominator calculation, described in Section 4.6.  All fifteen roots of
that denominator lying in the three tested prime fields have now been
specialized exactly: eleven lower the relative length from \(18\) to \(17\),
while four preserve the complete \(2+2+14\) profile even though the chosen
border chart fails.  The non-linear exceptional closed points and the
characteristic-zero determinant remain unclassified.  In particular, this
note does not prove that bidegree \((3,3)\) is safe and does not produce an
all-order counterexample.

## 2. Exact relative-period realization

Write \(U=(u_{iq})\), \(W=(w_{qj})\), with \(0\leq i,j\leq3\) and
\(q=1,2\).  Define

\[
\begin{aligned}
 A_q(u)&=\sum_{i=0}^3u_{iq}u^{3-i},\\
 B_q(u,t)&=\sum_{j=0}^3
 w_{qj}u^{j-3}t^j(1-t)^{3-j},\\
 P_{U,W}(u,t)&=A_1(u)B_1(u,t)+A_2(u)B_2(u,t).
\end{aligned}                                             \tag{2.1}
\]

If \(\Phi_C(x,y)=\sum c_{ij}x^iy^j\), direct multiplication gives

\[
 P_{U,W}(u,t)=
 \Phi_C\left(1,u,t,\frac{1-t}{u}\right).                  \tag{2.2}
\]

The beta identity applied to the diagonal contraction formula gives

\[
\boxed{
 \nu_m(C)
 =\operatorname {CT}_u\int_0^1P_{U,W}(u,t)^m\,dt.}
 \tag{2.3}
\]

Thus

\[
 \sum_{m\geq0}\nu_m(C)s^m
 =\operatorname {CT}_u\int_0^1
 \frac{dt}{1-sP_{U,W}(u,t)}.                              \tag{2.4}
\]

The checker verifies (2.2) exactly and verifies (2.3) through order four
at both integral factor points.  The factorial normalization in (1.2) is
essential: (2.4) is not the ordinary generating function of the raw
\(\mu_m\).

## 3. The toric number and the forced endpoint discriminant

The generic exponent polygon of \(P_{U,W}\) is

\[
 \operatorname {conv}\{(-3,0),(0,0),(3,3),(-3,3)\}.
 \tag{3.1}
\]

Its Euclidean area is \(27/2\), hence its normalized volume is \(27\).
The low-\(m\)-degree recurrence order agrees with this number, while the
relative logarithmic quotient and the lower-order recurrence both have
rank 18.  Neither number follows automatically from the other.

Indeed, the beta substitution forces the \(u=-3\) face to be

\[
 c_{30}u^{-3}(1-t)^3.                                    \tag{3.2}
\]

For \(c_{30}\ne0\), this face has a triple torus root at \(t=1\).
Consequently the constrained beta family lies on the face discriminant of
the unconstrained Laurent-polynomial family.  A generic stable-GKZ rank
statement therefore cannot be applied to (2.4) without a relative
endpoint analysis.  The number 27 does not by itself prove the recurrence
or bound initial data.

## 4. The modular recurrence probe

For either factor point, the adjacent C++ helper computes the moments
without expanding a sixteen-variable polynomial.  It uses

\[
 F^m=\sum_{k=0}^m\binom mk
 A_1^kA_2^{m-k}P_1^kP_2^{m-k}                            \tag{4.1}
\]

and updates all binary-cubic products by multiplication by one cubic.
This gives the exact sequence modulo the chosen prime.

For order \(R=27\) and coefficient degree \(D=11\), the ansatz

\[
 \sum_{j=0}^{27}
 \left(\sum_{e=0}^{11}a_{j,e}m^e\right)\nu_{m+j}=0
 \tag{4.2}
\]

has \(28\cdot12=336\) unknowns.  The probe uses 335 equations to select a
null vector and then tests all remaining available equations through
\(m+27=500\).  The 139 holdout equations vanish in all six
point/prime combinations.

The monic \(j=27\) coefficient reconstructed from the three primes is
(1.4).  Its eight linear factors are positive for every \(m\geq0\).
The cubic has no integer root because it has no root modulo \(29\).
Therefore

\[
 L(m)\ne0\qquad(m\in\mathbb Z_{\geq0}).                  \tag{4.3}
\]

This eliminates integer singular steps from the common \(m\)-part of the
observed forward coefficient.  It does not eliminate a parameter factor
introduced when a universal recurrence is cleared of denominators.

### 4.1 The exact relative logarithmic quotient

Put \(Q=u^3P\).  The two logarithmic critical polynomials are

\[
 A=uQ_u-3Q,\qquad C=t(1-t)Q_t.                            \tag{4.4}
\]

At the first integral exact-rank-two point, exact characteristic-zero
calculation gives

\[
\begin{aligned}
 \dim_\mathbb Q\mathbb Q[u,t]/(A,C)&=30,\\
 (A,C):u^\infty&=(A,C):u^6,\\
 \dim_\mathbb Q\mathbb Q[u,t]/((A,C):u^\infty)&=18.
\end{aligned}                                             \tag{4.5}
\]

The saturated quotient has standard monomial basis

\[
\begin{gathered}
1,u,u^2,u^3,u^4,\quad
t,ut,u^2t,u^3t,\quad
t^2,ut^2,u^2t^2,\\
t^3,ut^3,u^2t^3,\quad t^4,ut^4,\quad t^5.
\end{gathered}                                             \tag{4.6}
\]

Saturating further by the two endpoints decomposes the length as

\[
 18=2_{t=0}+2_{t=1}+14_{\mathrm{interior}}.               \tag{4.7}
\]

This is the first exact finite relative-cohomology model for the proposed
recurrence.

### 4.2 Exact characteristic-zero cyclic splitting

At the first integral rank-two point, put

\[
\begin{aligned}
 I_{\rm int}&=((A,Q_t):u^\infty):(t(1-t))^\infty,\\
 I_0&=(A,t):u^\infty,\qquad
 I_1=(A,t-1):u^\infty.
\end{aligned}
\]

Exact characteristic-zero ideal arithmetic gives the pairwise-comaximal
decomposition

\[
 (A,C):u^\infty=I_{\rm int}\cap I_0\cap I_1,\qquad
 18=14+2+2.                                             \tag{4.7a}
\]

The interior standard monomials are

\[
 1,u,u^2,u^3,\ t,ut,u^2t,u^3t,\ t^2,ut^2,u^2t^2,\
 t^3,ut^3,t^4.                                         \tag{4.7b}
\]

Thus interior localization removes precisely
\(u^4,u^2t^3,ut^4,t^5\) from (4.6).

This interior algebra is also the exact logarithmic critical algebra of
the Laurent polynomial \(P=Q/u^3\).  Indeed, the same calculation proves

\[
 \bigl((uQ_u-3Q,tQ_t):(ut)^\infty\bigr)=I_{\rm int}.
\]

In particular, the logarithmic critical rank at this rank-two fiber is
exactly \(14\), with the standard monomials (4.7b).  This identifies the
correct fixed-fiber Picard--Fuchs target rank and explains why the
sampled common Ore factor below has order \(14\).  It does not by itself
prove the scalar recurrence: a filtered twisted-de-Rham reduction must
still control classes at infinity and the interval endpoints.

There is a more useful cyclic description.  Eliminate \(u,t\) after
adjoining \(z u^3-Q\), so that \(z=P=Q/u^3\).  The eliminant degrees on
the interior and two endpoint summands are \(14,2,2\), respectively.  The
endpoint polynomials are

\[
\begin{aligned}
 p_0(z)&=1259712z^2-22590364z+102200491,\\
 p_1(z)&=9747z^2-5218042z+846742771.
\end{aligned}                                           \tag{4.7c}
\]

The degree-14 interior polynomial is stored exactly in the generated
artifact.  All three factors are squarefree and pairwise coprime, and the
degree-18 relative eliminant is their product up to a rational scalar.
Since every eliminant degree equals the length of its algebra,

\[
 A_{\rm int}\simeq\mathbb Q[z]/(p_{\rm int}(z)),
 \qquad 1,P,\ldots,P^{13}\ \text{is a basis}.            \tag{4.7d}
\]

This supplies an exact characteristic-zero Krylov basis for the interior
critical algebra.  It does **not** yet supply the period connection.  The
Chinese-remainder idempotents \(e_{\rm int},e_0,e_1\) were evaluated
against the exact moments through \(\nu_{18}\).  Both
\((\nu(e_0),\nu(Pe_0))\) and
\((\nu(e_1),\nu(Pe_1))\) are nonzero.  Hence the raw period functional
does not descend through the ordinary Jacobian quotient: gradient terms
vanish only after the \(m\)-dependent divergence correction.  Any
derivation of the order-14 recurrence must construct that correction
explicitly.

The first divergence certificate is now explicit.  The interior
saturation exponent with respect to \(u\) is five.  If
\(p_{\rm int}\) is the degree-14 eliminant, exact lifting through
\((A,Q_t)\) gives polynomials \(X,Y\) satisfying

\[
 \boxed{
 u^{47}p_{\rm int}(P)=X(uQ_u-3Q)+YQ_t,
 }
 \qquad\text{or}\qquad
 u^{44}p_{\rm int}(P)=X D_uP+Y\partial_tP.              \tag{4.7e}
\]

The cleared eliminant has 2752 terms; \(X\) and \(Y\) have respectively
6750 and 6791 terms.  Their divergence

\[
 D_uX+\partial_tY
\]

has 6749 terms.  Integration by parts therefore gives the exact seed

\[
\begin{aligned}
 (m+1)\int u^{44}p_{\rm int}(P)P^m
 &=-\int(D_uX+\partial_tY)P^{m+1}\\
 &\quad+[YP^{m+1}]_{t=0}^{t=1}.                         \tag{4.7f}
\end{aligned}
\]

The two restrictions of \(Y\) are nonzero, with 45 terms at \(t=0\) and
85 terms at \(t=1\).  Ordinary Jacobian reduction puts the divergence in
all fourteen interior coordinates and the two restrictions in both
coordinates of their endpoint algebras.  This is not yet the connection:
the gradient parts discarded by these normal forms must themselves be
lifted recursively with their \(m\)-dependent integration-by-parts
weights.  Thus the endpoint extension is not bookkeeping; it appears in
the first exact connection identity and survives its first reduction.

#### 4.2.1 Exact rational \(D\)-module seed

The ordinary generating function of the normalized moments has the
rational-period presentation

\[
 \sum_{m\geq 0}\nu_m z^m
 =\operatorname{CT}_u\int_0^1\frac{1}{1-zQ/u^3}\,dt
 =\oint\int_0^1\frac{u^2}{u^3-zQ(u,t)}\,dt\,du.
\]

Put \(H=u^3-zQ\).  Macaulay2's exact first-order annihilator calculation
for \(H^{-1}\) returns 34 generators.  The resulting left ideal is
holonomic of rank one.  Taking the kernel of multiplication by \(u^2\)
gives 76 exact differential operators annihilating the specific
integrand \(u^2/H\); this ideal is again holonomic of rank one.  Thus the
input to a \(D\)-module pushforward is now an exact all-order object, not
a fitted recurrence.

This does not yet prove an operator for the interval period.
The sequential pushforward in \(t\) and \(u\) remains to be completed,
and any resulting ambient operator must retain or discharge the
\(t=0,1\) certificate boundaries.  The exact annihilator counts,
holonomicity tests, ranks, and software versions are stored in
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_dmodule_picard_fuchs_research.json`.

### 4.3 The natural order-18 recurrence shape

At both integral rank-two points and all three primes, the ansatz

\[
 \sum_{j=0}^{18}
 \left(\sum_{e=0}^{18}b_{j,e}m^e\right)\nu_{m+j}=0
 \tag{4.8}
\]

uses 360 fitting equations and passes 83 unused equations through order
460.  An \(m\)-degree-17 ansatz fails at both points modulo \(1000003\).
The monic forward coefficient factors as

\[
 \prod_{k\in\{44,46,47,49,50,52,53,55\}}(3m+k)
 \cdot H_{U,W}(m),                                       \tag{4.9}
\]

where \(H_{U,W}\) is a point-dependent monic decic.  The degree-eight
factor is common to all six probes.  The decic is absent from the
order-27 forward coefficient, but the Ore calculation below shows that
this is not explained by a direct left multiplication
\(R_{27}=Q R_{18}\).

A naive certificate polynomial of \(m\)-degree 17 cannot prove (4.8).
After shifting the recurrence to the exponent used in integration by
parts, its leading \(m^{18}\) coefficient has nonzero remainder in every
one of the eighteen basis coordinates (4.6), even after saturation by
\(u\).  The certificate must therefore use a rational-in-\(m\) relative
connection, a higher-degree syzygy cancellation, or an equivalent
desingularized operator.

### 4.4 The sampled order-14 interior factor

Work in the shift Ore ring

\[
 \mathbb F_p(m)[S;\sigma],\qquad
 \sigma(f(m))=f(m+1),\qquad Sf(m)=f(m+1)S.               \tag{4.10}
\]

Exact left Euclidean division of \(R_{27}\) by \(R_{18}\) first leaves an
order-17 remainder.  At every one of the six point/prime samples, the
successive remainder orders are

\[
 17,\ 16,\ 15,\ 14,\ -1.                                \tag{4.11}
\]

Consequently their monic greatest common right divisor is \(G_{14}\), and
exact division gives

\[
 R_{18}=Q_4G_{14},\qquad R_{27}=Q_{13}G_{14}.             \tag{4.12}
\]

After taking the least common denominator of the rational coefficients
of \(G_{14}\) and removing their polynomial gcd, all fifteen primitive
coefficients have degree 58.  The forward coefficient contains

\[
 \prod_{k\in\{32,34,35,37,38,40,41,43\}}(3m+k)           \tag{4.13}
\]

and a residual factor of degree 50.  Direct substitution in the computed
sequence verifies the primitive order-14 relation for every
\(0\leq m\leq486\), giving 487 exact modular identities per sample.

The equality

\[
 \operatorname {ord}(G_{14})=14
 =\dim(\text{interior logarithmic quotient})             \tag{4.14}
\]

is the structural source of the Picard--Fuchs description.  It
agrees with the exact characteristic-zero cyclic algebra (4.7d), while
the two endpoint pieces account for the four extra dimensions in the
relative critical model.  The nonzero idempotent-period audit above rules
out the shortcut of simply projecting the raw period to that algebra.
Section 4.4.2 constructs and certifies the \(m\)-dependent divergence
at the first point modulo \(1000003\).  Section 4.4.3 gives a stable
exact characteristic-zero operator lift and verifies its leading
divergence level, but not the remaining descent or endpoint identities;
the generic-parameter operator also remains unproved.  The exact Krylov
construction implemented in
[`SUPERELLIPTIC_DERHAM_ENGINE.md`](../plane-jc/SUPERELLIPTIC_DERHAM_ENGINE.md)
remains the model for lifting this fixed-fiber calculation.

#### 4.4.1 Two fixed-fiber scalar shortcuts fail

Two bounded calculations at the first point modulo \(1000003\) sharpen
what the connection must retain.  First, the leading \(m^{58}\)
coefficient of the sampled \(G_{14}\), interpreted as a Laurent
integrand, has nonzero class in the length-eighteen relative critical
algebra.  Hence a direct polynomial-in-\(m\), zero-boundary divergence
certificate whose certificate degree is at most 57 cannot start.

Allowing a degree-58 leading Koszul syzygy does not repair this
obstruction.  Put
\[
 A=uQ_u-3Q,\qquad C=t(1-t)Q_t,\qquad T=t(1-t).
\]
On the relative critical quotient, the linear map
\[
 R\longmapsto
 \left[
 Q\bigl(D_u(CR)+\partial_t(-TAR)\bigr)
 \right]
\]
has rank zero, while the leading obstruction is nonzero.  Thus the
obvious leading correction \((X_{58},Y_{58})=R(C,-A)\) also cannot
produce the sampled scalar relation.

These are exact modular no-go calculations for the displayed ansätze,
not a disproof of the sampled recurrence.  They confirm that the
\(14+2+2\) endpoint-extended connection must be constructed before
scalarization; prescribing \(G_{14}\) inside a zero-boundary scalar
ansatz discards essential endpoint data.

#### 4.4.2 An all-order modular fixed-fiber certificate

The first point modulo \(p=1000003\) now has an exact certificate for
the sampled \(G_{14}\).  Write

\[
 G_{14}=\sum_{j=0}^{14}g_j(m)S^j,\qquad
 {\cal S}(m,P)=\sum_{j=0}^{14}g_j(m)P^j.
\]

In the Laurent presentation
\(\mathbb F_p[u,U,t]/(uU-1)\), descending the 58 powers of \(m\)
produces

\[
 X(m)=\sum_{r=0}^{57}X_r m^r,\qquad
 Y(m)=\sum_{r=0}^{57}Y_r m^r
\]

with

\[
 Q{\cal S}(m,P)=
 Q\bigl(D_uX+\partial_tY\bigr)+m(AX+Q_tY).              \tag{4.14a}
\]

Canonical reduction modulo \(uU-1\) is essential here: it turns the
previous multi-hour expanded lift into three restartable chunks.  The
uncorrected descent reaches \(m^0\) with a nonzero 307276-term residual
\(T\).  This is not a failed recurrence.  The final Koszul freedom

\[
 (X_0,Y_0)\longmapsto(X_0+Q_tR,\;Y_0-AR)
\]

changes the divergence by

\[
 H(R)=Q_t(D_u+3)R-A\partial_tR.
\]

A descending \(t\)-block solve gives a 298606-term Laurent polynomial
\(R\) satisfying

\[
 \boxed{T=QH(R)}.                                      \tag{4.14b}
\]

An independent Singular replay verifies all 58 coefficient identities,
the three restart residuals, and (4.14b).  After applying this correction,
(4.14a) is equivalent to the exact telescoping identity

\[
 {\cal S}(m,P)P^m
 =D_u\!\left(X(m)P^m\right)
  +\partial_t\!\left(Y(m)P^m\right).                   \tag{4.14c}
\]

The endpoints are not discarded.  At this fiber

\[
\begin{aligned}
 P(u,0)&=11+23u^{-1}+91u^{-2}+216u^{-3},\\
 P(u,1)&=354+149u+37u^2+19u^3.
\end{aligned}
\]

Consequently each endpoint constant term is an
exponential-polynomial, respectively \(11^mE_0(m)\) and
\(354^mE_1(m)\).  Exact extraction from the corrected \(Y\) gives

\[
 E_0=E_1=0\quad\hbox{in }\mathbb F_{1000003}[m].        \tag{4.14d}
\]

Taking \(\operatorname{CT}_u\int_0^1\) in (4.14c) therefore proves

\[
 \boxed{
   \sum_{j=0}^{14}g_j(m)\nu_{m+j}=0
   \quad\text{for every }m\geq0
   \text{ over }\mathbb F_{1000003}.
 }                                                       \tag{4.14e}
\]

This replaces the 487-row observation by an all-order modular theorem
at one fixed fiber.  By itself it does **not** reconstruct the operator
or certificates over characteristic zero, prove a generic rank-two
parameter identity, cover the other five point/prime samples, or
classify the exceptional locus.  The finite characteristic-zero lift
obtained subsequently is recorded next.

#### 4.4.3 Characteristic-zero operator lift and leading descent

The fixed-fiber operator itself can now be lifted exactly.  We computed
its normalized image at 205 consecutive primes above \(10^6\).  A
simultaneous projective reconstruction from the first 200 images uses one
common denominator rather than 885 independent rational reconstructions.
LLL first succeeds in dimension 24.  After primitive integer
normalization the largest coefficient has 2397 bits; all 885 coefficients
then agree at each of five fresh holdout primes.

This is not merely a modular consistency check.  Exact rational replay
at the original integral point gives all 27 available recurrence
identities through moment \(40\).  The forward coefficient factors
exactly as

\[
 \prod_{k\in\{32,34,35,37,38,40,41,43\}}(3m+k)
 \cdot h_{50}(m),                                      \tag{4.14f}
\]

with \(h_{50}\in\mathbb Z[m]\) of degree 50.  Thus the stored primitive
integer operator is a stable exact characteristic-zero lift of the
sampled \(G_{14}\), within this finite reconstruction and replay scope.
These checks alone still do not prove that it annihilates every moment.

The first coefficient of the characteristic-zero telescoping descent is
also exact.  Over
\(\mathbb Q[u,U,t]/(uU-1)\), the coefficient of \(m^{58}\) admits the
required lift with a 5722-term \(X_{57}\) and a 5769-term \(Y_{57}\).
An independent ambient Singular calculation replays that identity and
matches the 6354-term restart residual at \(m^{57}\).  Hence the top
descent equation is closed over \(\mathbb Q\); this is stronger than
reduction modulo many primes, but it leaves levels \(m^{57},\ldots,m^1\),
the terminal \(m^0\) syzygy, and both endpoint identities open.

A direct attempt to expand five rational descent levels exceeded 1800
seconds and reached 4.3 GB without completing.  This is a performance
failure, not a mathematical obstruction.  It identifies the remaining
problem as representation growth in the certificates.

The appropriate next representation is reduction-based creative
telescoping.  Bostan--Lairez--Salvy's
[Griffiths--Dwork algorithm](https://arxiv.org/abs/1301.4313) is
specifically designed to derive a Picard--Fuchs operator without
expanding certificates.  Lairez's
[extended reduction](https://arxiv.org/abs/1404.5069) covers singular
hypersurfaces, and the newer
[mixed \(D\)-module reduction](https://arxiv.org/abs/2504.12724) has an
experimental
[Julia implementation](https://github.com/HBrochet/MultivariateCreativeTelescoping.jl).
The generic projective interface of that implementation did not finish
our homogenized seed within 900 seconds.  In this problem it also loses
the already proved toric \(14+2+2\) decomposition, so it is not the final
engine.

The proof-producing route is therefore a specialized relative reduction
map: reduce the successive \(z\)-derivatives of

\[
 \frac{u^2}{u^3-zQ(u,t)}
\]

directly in the known 14-dimensional interior quotient while carrying
the two two-dimensional endpoint trace blocks.  An exact dependence in
this \(14+2+2\) state space supplies the all-order operator without
materializing every \(X_r,Y_r\); the retained endpoint blocks supply the
boundary audit.  Only after this compact reduction is certified should
one reconstruct a full expanded certificate, if one is still desired.

The combined checker
`scripts/verify_two_pair_sic_bidegree33_rank_two_characteristic_zero_lift.py`
replays the 205 prime images, the 27 exact rational moment identities,
and the independent \(m^{58}\) descent.  Its status is deliberately
finite: no characteristic-zero all-order theorem is claimed here.

#### 4.4.4 Compact relative Picard--Fuchs bridge

The beta chart admits a useful birational compression.  Put

\[
 x=u,\qquad y=\frac{ut}{1-t},\qquad t=\frac{y}{x+y}.
\]

An exact symbolic calculation gives

\[
 P(u,t)=\frac{\Phi(x,y)}{(x+y)^3},\qquad
 \frac{du}{u}\,dt=\frac{dx\,dy}{(x+y)^2},              \tag{4.14g}
\]

where

\[
\begin{aligned}
\Phi={}&19x^3y^3+17x^3y^2+13x^3y+11x^3\\
 &+37x^2y^3+31x^2y^2+29x^2y+23x^2\\
 &+149xy^3+127xy^2+113xy+91x\\
 &+354y^3+302y^2+268y+216.
\end{aligned}
\]

Thus the generating form becomes

\[
 \frac{(x+y)\,dx\,dy}{(x+y)^3-z\Phi(x,y)},             \tag{4.14h}
\]

with a sixteen-term denominator of total degree six.  The projective
Griffiths--Dwork calculation of Brochet--Chyzak--Lairez now finishes: it
returns a closed-cycle operator of differential order eight, with coefficient
degrees \(82,81,\ldots,74\).  That operator does not annihilate the interval
period, because its discarded exact forms retain boundary traces.

The interval calculation has a smaller inhomogeneous relation

\[
 L_8F=B_{55},\qquad
 L_8=\sum_{k=0}^8p_k(z)\partial_z^k,                   \tag{4.14i}
\]

where every \(p_k\) has degree at most 72 and \(B_{55}\) has degree exactly
55.  At \(p=1000003\), the normalized relation is the unique vector in the
tested order-eight, degree-72 box; the adjacent boxes \((8,71)\), \((7,72)\),
and \((7,80)\) have zero nullity.  Its support satisfies
\(z^k\mid p_k\) and \(\deg p_k=k+64\), except that \(p_3\) starts at \(z^4\).
The coefficient conversion \(n=m+64\) gives

\[
 R_{64,8}=\sum_{j=0}^{64}r_j(m)S^j,\qquad
 z^e\partial_z^k\longmapsto
 S^{64-e+k}(m+64-e+k)_{\underline{k}}.                 \tag{4.14j}
\]

Every \(r_j\) has degree eight.  A complete independent modular replay
verifies all eight descending divergence identities, removes the 138651-term
terminal residual with a 132615-term Koszul correction, and proves that both
corrected endpoint exponential-polynomials vanish.  Therefore

\[
 R_{64,8}\nu=0\quad\text{for every }m\geq0
 \quad\text{over }\mathbb F_{1000003}.                 \tag{4.14k}
\]

This is a second all-order modular certificate, now with only eight descent
levels.  Its purpose is the characteristic-zero lift, not a stronger modular
statement than (4.14e).

The compact relation itself has a stable exact lift.  Simultaneous projective
reconstruction from 95 prime images, followed by five fresh holdouts,
recovers all 657 coefficients; the maximum primitive coefficient height is
1527 bits.  Exact rational moments confirm that the residual in (4.14i) is
supported precisely in degrees zero through 55 on all 93 available rows.
Converting the lifted \(L_8\) gives a primitive characteristic-zero
\(R_{64,8}\).  Exact rational shift-Ore division, taking about eight minutes
and 136 MB on the reference machine, proves

\[
 \boxed{R_{64,8}=Q_{50}G_{14,58}}                      \tag{4.14l}
\]

with zero remainder.  The first and last coefficients of \(Q_{50}\) have
degree pair \((0,50)\); the other 49 have degree pair \((50,100)\).  All 51
coefficients of the forward denominator are positive, hence that coefficient
never vanishes at an integer \(m\geq0\).  A separate rank-two formula computes
exact moments through the required range in under a second and proves

\[
 (G_{14,58}\nu)_0=\cdots=(G_{14,58}\nu)_{49}=0.        \tag{4.14m}
\]

Consequently an all-order characteristic-zero certificate for \(R_{64,8}\)
would force \(G_{14,58}\nu=0\) for every \(m\geq0\): the quotient recurrence
and its complete initial data now have no escape.

The remaining gate is specifically the characteristic-zero relative
certificate for \(R\).  Its first expanded original-coordinate descent level
reached the 900-second guard at 8.6 GB; the same modular level takes 47 seconds
below 1 GB, and all eight modular levels take under four minutes.  Ordinary
top-pole Griffiths reduction is insufficient: the \(H^{-9}\) numerator has an
18-term remainder in the Laurent Jacobian quotient.  The next proof-producing
engine is therefore extended relative reduction in the compact chart, or
termwise reconstruction of the eight-level modular certificate.  Operator
reconstruction, exact factorization, forward nonvanishing, and the 50 initial
values are no longer open.

A second-prime scout at \(p=1000033\) makes the reconstruction alternative
concrete.  All eight \(Y_r\) supports, five of the eight \(X_r\) supports, and
the 132615-term terminal \(R\)-support are identical at the two primes.  Each
of the other three \(X_r\) supports differs by exactly one monomial, consistent
with an isolated coefficient vanishing at one prime.  Thus union-support CRT
is the next bounded experiment.  This two-prime support agreement is not
itself a characteristic-zero lift.

The modular certificate is replayed by
`scripts/verify_two_pair_sic_bidegree33_rank_two_compact_relative_modular_all_order.py`.
The characteristic-zero reconstruction and exact factorization are replayed
by `scripts/verify_two_pair_sic_bidegree33_rank_two_compact_relative_pf_lift.py`.
The latter deliberately remains a finite bridge, not an all-order
characteristic-zero theorem.

### 4.5 Why expanded universal interpolation is not the next step

On the scaling family \(C(s)=(1+s)C_0\), the modular rational
reconstructor recovers the predicted coefficient degrees \(27-j\) for
the order-27 recurrence, with denominator one and independent holdouts.
This validates both the recurrence normalization and the interpolation
engine.

On the generic quadratic factor pencil

\[
 (U(s),W(s))=(U_0+sU_1,W_0+sW_1),                        \tag{4.15}
\]

256 samples with twelve holdouts find no rational interpolant for eight
representative recurrence coefficients within the resulting degree
window.  In particular, no candidate appears with combined
numerator/denominator degree at most 243 on this one pencil.  Therefore a
fully expanded twelve-parameter recurrence is the wrong immediate
representation.  The interior connection suggested by (4.14), with the
endpoint extension retained for certificates, is the compact target.

### 4.6 Exact border-basis denominator on the generic pencil

There is nevertheless an exact compact calculation available on the
pencil (4.15).  Over each of
\(\mathbb F_{1000003}(s)\), \(\mathbb F_{1000033}(s)\), and
\(\mathbb F_{1000037}(s)\), saturating the logarithmic Jacobian ideal

\[
 (uQ_u-3Q,\ t(1-t)Q_t):u^\infty
\]

has exponent six, quotient length eighteen, and the same six leading
border monomials

\[
 u^5,\quad u^4t,\quad u^3t^2,\quad u^2t^4,\quad ut^5,\quad t^6.
 \tag{4.16}
\]

After reduction and monic normalization, every border relation has
nineteen terms.  Their denominator degrees, in Singular basis order, are

\[
 74,\ 74,\ 88,\ 94,\ 74,\ 74.                           \tag{4.17}
\]

At all three primes the degree-\(74\) denominator is the common gcd of
the three distinct denominators.  Dividing it from the degree-\(88\) and
degree-\(94\) denominators leaves coprime factors of degrees \(14\) and
\(20\).  Hence this monomial chart is regular away from one squarefree
degree-\(108\) polynomial in \(s\).

This is the first explicit exceptional polynomial in the programme, but
only on one modular pencil and only for the logarithmic border basis.  Its
roots can include monomial-chart artifacts, and it is not the universal
telescoping denominator \(D(U,W)\).  The stable profile at three primes is
evidence for a characteristic-zero pencil determinant of the same
degree, not a reconstruction of that determinant.

Every linear factor of the degree-\(108\) polynomial over the three base
fields can nevertheless be specialized exactly.  There are respectively
\(4,6,5\) such roots.  All fifteen specializations keep the coefficient
matrix at exact rank two and retain saturation exponent six.  Their
relative-length profiles are

\[
\begin{array}{c|c|c}
\text{number}&\text{profile}&\text{denominator component}\\ \hline
4&2+2+14&74\\
4&1+2+14&14\\
4&2+1+14&14\\
2&2+2+13&14\\
1&2+2+13&20.
\end{array}                                             \tag{4.18}
\]

Thus all ten accessible roots of the extra degree-\(14\) factor and the one
accessible root of the extra degree-\(20\) factor are genuine
relative-length changes on these reductions.  The four accessible roots of
the common degree-\(74\) factor preserve the full endpoint/interior length
profile and are border-chart failures at this level.  This is an exact
classification over the stated finite fields, not evidence that every
geometric root of the three factors has the same type in characteristic
zero.

## 5. The actual finite-determination statement to prove

Let

\[
 R_2=\mathbb Q[u_{iq},w_{qj}]
\]

be the factor parameter ring.  The modular data suggest a recurrence over
\(\operatorname {Frac}(R_2)\) of the form

\[
 G_{14}\nu=
 \sum_{j=0}^{14}g_j(m;U,W)\nu_{m+j}=0,\qquad
 \deg_m g_j\leq58,                                      \tag{5.1}
\]

whose sampled forward coefficient contains (4.13) and a
parameter-dependent degree-50 residual factor.  The order-27,
\(m\)-degree-11 operator with forward coefficient (1.4) is a sampled left
multiple of (5.1), not a left multiple of the order-18 operator.  After
clearing a primitive common parameter denominator \(D_{14}(U,W)\), the
forward coefficient of (5.1) should have the form

\[
 D_{14}(U,W)
 \prod_{k\in\{32,34,35,37,38,40,41,43\}}(3m+k)
 H_{50}(m;U,W),                                         \tag{5.2}
\]

possibly after removing a common parameter content.

An exact telescoping identity proving (5.1), together with (5.2), would
give the sharper desired first open stratum after auditing the integer
steps of \(H_{50}\):

> On \(D_{14}\ne0\), 14 consecutive zero moments beginning at any positive
> order propagate to all later orders.

The same statement transfers to the raw \(\mu_m\) after multiplying (5.1)
by \((3(m+14)+1)!\).  This only introduces explicit nonzero factorial
ratios; it does not create characteristic-zero singular steps.  It does
change the integral recurrence lattice after reduction modulo \(p\).
Consequently the normalized-period operator and the raw-moment operator
must be stored separately in any arithmetic continuation.

The exceptional analysis is then finite and algebraic only after
\(D_{14}(U,W)\) is known.  On each component of \(D_{14}=0\), one must
specialize the telescoping matrix, remove its parameter content, and
recompute its generic rank and forward coefficient.  Merely requiring the
leading coefficient to be nonzero as a polynomial in \(m\) is
insufficient: integer-step zeros must also be audited.  Formula (4.3)
performs that audit only for the order-27 common factor; no
characteristic-zero audit of \(H_{50}\) is available.

## 6. Interaction with the corrected moment system

The corrected ambient system is

\[
 \mu_1,\ldots,\mu_{12},\mu_{14}=0.                       \tag{6.1}
\]

The order-27 recurrence shape does not make (6.1) an all-order condition
by itself.  Even on its open parameter locus, forward propagation asks
for \(\mu_1,\ldots,\mu_{27}\).  The corrected equations omit
\(\mu_{13}\) and do not include orders \(15,\ldots,27\).

The sampled order-14 factor narrows the candidate missing bridge set to

\[
 \boxed{\mu_{13}}.                                       \tag{6.2}
\]

Indeed, (6.1) together with \(\mu_{13}=0\) supplies the consecutive block
\(\mu_1,\ldots,\mu_{14}\).  This reduction is conditional on a universal
certificate for (5.1), a parameter-locus analysis, and the integer-step
audit of its degree-50 forward residual.

There are three exact ways the recurrence could still settle (6.1):

1. reduce \(\mu_{13}\) to zero in the radical of the corrected rank-two
   ideal on \(D_{14}\ne0\);
2. derive a lower-order recurrence after quotienting by that ideal; or
3. show directly that (6.1) has only nullcone points on
   \(D_{14}\ne0\), then recurse on \(D_{14}=0\).

Until one of these calculations is certified, (6.1) remains undecided on
the rank-two stratum.

## 7. Next certificate-producing computation

The immediate target is no longer an unspecified Picard--Fuchs search.
It is the following fixed-size reconstruction problem.

1. Starting from the exact one-fiber cyclic basis (4.7d), construct the
   \(m\)-dependent divergence reduction at that characteristic-zero
   point.  Retain the two endpoint extensions, because the idempotent
   period audit proves that ordinary Jacobian projection is insufficient.
2. Apply the exact Krylov calculation to this discrete connection and
   recover the order-14, \(m\)-degree-58 target (5.1), together with its
   polynomial divergence certificates.
3. Lift the construction over the generic rank-two factor field while
   retaining the six border-pivot determinants.  The
   degree-\(74+14+20\) pencil factorization in Section 4.6 is the
   specialization target.
4. Recover the order-4 and order-13 left quotients in (4.12), thereby
   explaining both displayed low-degree recurrences and normalizing the
   order-27 forward coefficient to (1.4).
5. Store the relative divergence certificates, including the
   \(t=0,1\) endpoint terms, and verify their polynomial identities
   independently.
6. Reduce the single bridge moment \(\mu_{13}\) and the border-pivot
   determinant against the corrected rank-two moment ideal.  On the
   exceptional locus, separate chart changes from the length-\(17\)
   relative modules and recompute the scalar operator on each component.
7. After the characteristic-zero operator is certified, primitively
   normalize both its period and raw-moment forms over the parameter ring.
   At selected good primes, compute the companion \(p\)-curvature together
   with the uncancelled step-matrix Smith/singularity ledger.  The latter,
   not the generic curvature characteristic polynomial, is the input for
   any prime-power moment phase claim.

The modular probes have now separated the full relative rank 18 from a
sampled scalar interior factor of order 14 and from its order-18 and
order-27 left multiples.  They have also exposed one exact modular pencil
denominator.  The exact one-fiber cyclic algebra is now constructed, and
it explicitly shows why the divergence connection cannot be skipped.
Connection construction, certificate reconstruction, and
recursive analysis of the non-linear degree-\(108\) pencil fibers remain
the proof-producing gates.  The arithmetic postprocessing in step 7 is
specified now, but it does not precede those gates; see
[`TWO_PAIR_SIC_FROBENIUS_CURVATURE_BRIDGE.md`](TWO_PAIR_SIC_FROBENIUS_CURVATURE_BRIDGE.md)
for the order-one calibration.

## 8. Reproduction

Run

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_holonomic_probe.py
```

The command compiles the adjacent C++ helper in a temporary directory,
checks the exact period identities, runs the six modular recurrence
probes, checks the two degree-10 failures, and writes
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_holonomic_probe.json`.
The artifact labels the recurrence as modular evidence rather than a
certificate.

The relative-Jacobian and natural-order calculation is replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_relative_jacobian.py
```

It writes
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_relative_jacobian.json`.
The exact characteristic-zero interior split is replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_interior_cyclic_split.py
```

At the first integral rank-two point it certifies the pairwise-comaximal
critical-algebra decomposition \(18=14+2+2\), eliminates \(P=Q/u^3\) on
all three summands, and proves that
\(1,P,\ldots,P^{13}\) is a basis of the interior algebra.  It also
proves that the saturated toric logarithmic ideal
\(((uQ_u-3Q,tQ_t):(ut)^\infty)\) equals the interior ideal, giving the
exact logarithmic Picard--Fuchs target rank \(14\).  It then
computes the first nineteen exact moments and verifies that both endpoint
idempotent pairs have nonzero period values.  Thus the ordinary Jacobian
split is a connection seed, not a substitute for the \(m\)-dependent
divergence calculation.  The artifact is
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_interior_cyclic_split.json`.

The exact rational \(D\)-module seed is replayed by

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_dmodule_picard_fuchs.py \
  --annihilator-only
```

This requires Macaulay2, its `BernsteinSato` package, and Normaliz.  It
certifies the 34-generator rank-one holonomic annihilator of \(H^{-1}\)
and the 76-generator rank-one holonomic annihilator of \(u^2/H\), then
writes
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_dmodule_picard_fuchs_research.json`.
Omitting `--annihilator-only` requests the sequential \(t,u\)
pushforward.  That long computation is not part of the retained
certificate until it finishes and its endpoint terms are audited.

The exact Ore-factor comparison is replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_ore_gcd.py
```

It recomputes both recurrence operators and the 501 moments at all six
samples, performs exact left Euclidean division over
\(\mathbb F_p(m)[S;\sigma]\), and writes
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_ore_gcd.json`.
It verifies the order-14 common right factor directly on 487 moment rows
per sample; it does not certify a universal operator.

The exact modular pencil border basis is replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_pencil_border.py
```

It writes
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_pencil_border.json`.
The three finite rational-function-field calculations certify the stable
denominator profile \(74,\ 74+14,\ 74+20\) and the squarefree
degree-\(108\) chart polynomial.  They also classify all fifteen
base-field roots by coefficient rank, saturation exponent, total relative
length, and endpoint/interior profile.  They do not reconstruct a universal
parameter determinant or classify the non-linear exceptional closed points.

The generic and scaling pencil experiments are implemented in
`scripts/explore_two_pair_sic_bidegree33_rank_two_recurrence_line.py`;
their two artifacts are
`two_pair_sic_bidegree33_rank_two_recurrence_line.json` and
`two_pair_sic_bidegree33_rank_two_recurrence_scaling_line.json` in
`artifacts/generated-results/`.  They remain interpolation experiments
rather than recurrence certificates.
