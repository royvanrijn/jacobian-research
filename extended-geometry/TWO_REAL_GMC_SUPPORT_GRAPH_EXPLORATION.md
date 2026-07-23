# Support graphs for the two-real Gaussian moment problem

## 1. Updated starting point

For

\[
Z=(X+iY)/\sqrt2,\qquad W=(X-iY)/\sqrt2,\qquad U=ZW,
\]

one has

\[
\mathbb E(Z^aW^b)=\delta_{a,b}a!.
\]

Give \(Z,W\) rotational weights \(1,-1\).  The repository now proves more
than the preliminary support census: every polynomial of total degree at
most three satisfies \(\operatorname{GMC}(2)\).  Thus a possible
counterexample has degree at least four.

The proposed support-classification direction is still natural, but
“weight-support graph” needs a fixed definition.  Two plausible definitions
give different smallest cyclic supports.

## 2. Canonical circuit-incidence graph

Embed

\[
\mathbb C[Z,W]\hookrightarrow\mathbb C[U][T,T^{-1}],
\qquad Z\mapsto T,\quad W\mapsto UT^{-1}.
\]

Write

\[
P=\sum_{k\in S}T^kB_k(U).
\]

The expectation is

\[
\mathbb E(P^m)=
\mathcal L\!\left(\operatorname{CT}_T P^m\right),
\qquad \mathcal L(U^j)=j!.                         \tag{2.1}
\]

A canonical combinatorial object is the incidence graph between:

- the nonzero weights \(k\in S\setminus\{0\}\); and
- the primitive support-minimal nonnegative relations
  \(\alpha\) satisfying \(\sum_k\alpha_k k=0\).

For a rank-one weight set, every such circuit has one positive and one
negative weight.  If \(p\) positive and \(q\) negative levels occur, the
incidence graph is the subdivision of \(K_{p,q}\).  Consequently it is a
forest exactly when

\[
\min(p,q)\leq1.                                    \tag{2.2}
\]

The zero weight contributes a singleton invariant and does not create a
cycle.  Under this definition, the support \(\{-1,0,1\}\) is a forest.  The
first cycle is \(K_{2,2}\), so it requires at least two positive and two
negative nonzero weight levels.

This gives a precise restricted conjecture:

> **Circuit-forest conjecture.** If the nonzero rotational support of \(P\)
> has at most one positive level or at most one negative level, then \(P\)
> satisfies \(\operatorname{GMC}(2)\).

The known two-weight theorem is only the \(p=q=1\) base case.  The conjecture
also contains genuinely new stars \(K_{1,q}\) and \(K_{p,1}\).

## 3. Expanded coefficient graph

One may instead expand every \(B_k(U)\) into radial monomials and retain
parallel invariant products as separate edges.  This can create cycles even
for weights \(\{-1,0,1\}\).  But contracting the \(U=ZW\) direction merges
those parallel edges and returns the circuit-incidence graph above.

Therefore the two claims

1. “contract invariant \(ZW\)-coefficients,” and
2. “the first cyclic case is \(\{-1,0,1\}\),”

cannot both hold without extra graph structure.  Any theorem or search
should specify whether radial coefficient polynomials are vertices, edges,
labels, or contracted data.

## 4. The exact three-level functional equation

The smallest interaction between two different invariant polynomials is
nevertheless the three-level ansatz

\[
P=W A(U)+C(U)+ZB(U).
\]

Put

\[
D(U)=U A(U)B(U).
\]

Constant-term extraction gives

\[
\begin{aligned}
\mathbb E(e^{tP})
&=\mathcal L\!\left[
 e^{tC(U)}
 \sum_{r\geq0}\frac{t^{2r}D(U)^r}{(r!)^2}
\right]\\
&=\mathcal L\!\left[
 e^{tC(U)}I_0\!\left(2t\sqrt{D(U)}\right)
\right].                                           \tag{4.1}
\end{aligned}
\]

Thus all pure moments vanish exactly when (4.1) equals \(1\).  Equivalently,
for every \(m\geq1\),

\[
\sum_{r=0}^{\lfloor m/2\rfloor}
\frac{m!}{(r!)^2(m-2r)!}
\mathcal L\!\left(C^{m-2r}D^r\right)=0.            \tag{4.2}
\]

This is the desired finite-alphabet functional equation: it depends only on
the two invariant polynomials \(C,D\).  The prime-endpoint theorem below
closes both axes and rules out every coupled cancellation between them.

The ordinary moment series is even simpler:

\[
\boxed{
\sum_{m\geq0}\mathbb E(P^m)t^m
=
\mathcal L\!\left(
\frac1{\sqrt{(1-tC)^2-4t^2D}}
\right).}                                         \tag{4.3}
\]

Indeed, before applying \(\mathcal L\), the right side is the constant term
of

\[
\frac1{1-t(TB+C+T^{-1}UA)}.
\]

This algebraic resolvent is preferable to the Bessel series for
classification.  It also shows that the three-level problem depends only on
\(C\) and \(D\), not on the individual factorization \(D=UAB\).

### 4.1 All-degree Bessel--factorial rigidity

The functional equation admits a direct prime-endpoint proof.  No
differential equation, degree bound, algebraic-coefficient specialization,
or external one-variable theorem is needed.

> **Bessel--factorial rigidity theorem.** Let \(C,D\in\mathbb C[U]\), and put
> \[
> H_m(C,D)=
> \sum_{r=0}^{\lfloor m/2\rfloor}
> \frac{m!}{(r!)^2(m-2r)!}C^{m-2r}D^r.
> \]
> If
> \[
> \mathcal L(H_m(C,D))=0\qquad(m\geq1),
> \]
> then \(C=D=0\).

The essential congruence is the following reusable statement:

> **Prime-endpoint lemma.** For every odd prime \(p\),
> \[
> H_p(C,D)\equiv C^p,\qquad
> H_{2p}(C,D)\equiv C^{2p}+2D^p\pmod p.           \tag{4.5}
> \]

Put the coefficients of \(C,D\), together with the inverses of their lowest
nonzero coefficients, in their literal finite-type subring
\(R\subset\mathbb C\).  For every sufficiently large prime \(p\),
\(\operatorname{Spec}R\) has a characteristic-\(p\) fiber, and the selected
coefficients remain units after reduction.  If
\(e=\operatorname{ord}_U D\geq2\operatorname{ord}_U C=2a\), cancel
\((ap)!\) from the \(p\)-th factorial moment and reduce (4.5); only the
lowest coefficient of \(C^p\) survives.  If \(e<2a\), cancel \((ep)!\) from
the \(2p\)-th moment; only twice the lowest coefficient of \(D^p\) survives.
Both are units, a contradiction.  The same two arguments directly cover
the axes.

The standalone note
[Prime-endpoint rigidity for the Bessel--factorial transform](PRIME_ENDPOINT_BESSEL_FACTORIAL_RIGIDITY.md)
gives the finite-type reduction lemma, all coefficient valuations, the
factorial-normalization audit, and the complete proof.

For the original three-level polynomial, \(D=UAB\).  Consequently:

> **All-degree three-level support theorem.** If
> \[
> P=WA(U)+C(U)+ZB(U)
> \]
> has \(\mathbb E(P^m)=0\) for every \(m\geq1\), then
> \[
> C=0,\qquad A=0\ \text{or}\ B=0.
> \]
> Thus its rotational support is one-sided, and \(P\) satisfies
> \(\operatorname{GMC}(2)\).

Two independent arithmetic regressions check the endpoint congruences and
both order cases:
[`verify_two_real_gmc_three_level_rigidity.py`](../scripts/verify_two_real_gmc_three_level_rigidity.py)
uses SymPy, while
[`audit_prime_endpoint_rigidity_independent.py`](../scripts/audit_prime_endpoint_rigidity_independent.py)
uses a separate pure-Python polynomial implementation.

The mixed GMC moments have an equally useful marked form.  For
\(Q=T^\ell R(U)\), only indices with the appropriate weight imbalance
survive in

\[
\mathcal L\!\left(
\operatorname{CT}_T
T^\ell R(U)
\exp(t(TB+C+T^{-1}UA))
\right).                                           \tag{4.4}
\]

Once a solution of (4.1) is found, (4.4) is a direct finite search for a
multiplier witnessing failure.

## 5. Exact low-degree theorem

In total degree at most \(n\), put

\[
r=\lfloor n/2\rfloor,\qquad s=\lfloor(n-1)/2\rfloor.
\]

After the first moment centers \(C\), the invariant coordinates are

\[
\begin{aligned}
C&=\sum_{j=1}^{r}c_j(U^j-j!),\\
D&=\sum_{j=1}^{2s+1}d_jU^j.                       \tag{5.1}
\end{aligned}
\]

Conversely, every nonzero \(D\) in (5.1) factors over \(\mathbb C\) as
\(UAB\) with \(\deg A,\deg B\leq s\): split the at most \(2s\) roots of
\(D/U\) into two groups.  Thus (5.1) exactly parameterizes this support
after quotienting the relative \(A/B\) scaling.

Coefficient charts cover \(C\ne0,D\ne0\).  Weighted scalar multiplication
normalizes one selected \(c_j\), a selected \(d_j\) is localized, and the
second moment eliminates the highest \(D\)-coefficient.  Exact
characteristic-zero calculations give:

| degree | invariant charts | moment cutoff | result |
|---:|---:|---:|---|
| 4 | 6 | 6 | reduced basis \([1]\) on every chart |
| 5 | 10 | 8 | reduced basis \([1]\) on every chart |
| 6 | 15 | 9 | reduced basis \([1]\) on every chart |

Therefore:

> **Bounded certificate regression.** In total degrees four, five, and six,
> the first \(6,8,9\) moments, respectively, already exclude every chart
> with all three weight components \(-1,0,1\) nonzero.

The exact replay is
[`verify_two_real_gmc_three_weight_low_degree.py`](../scripts/verify_two_real_gmc_three_weight_low_degree.py).

These calculations are now finite-cutoff certificates and independent
regressions for the all-degree theorem.  They remain useful because the
prime-endpoint proof uses arbitrarily large prime-indexed moments and does
not give the displayed small cutoffs.

## 6. Unit-star rigidity and the forest boundary

The prime method proves the smallest genuine star, and every star whose
unique weight on one side is \(1\).  For
\[
P=TB_1+C+T^{-1}UA_1+T^{-2}U^2A_2,
\]
the primitive invariants are
\[
X=UA_1B_1,\qquad Y=U^2A_2B_1^2.
\]
Their first prime-indexed endpoints are
\[
\begin{aligned}
H_p&\equiv C^p,\\
H_{2p}&\equiv C^{2p}+2X^p,\\
H_{3p}&\equiv C^{3p}+6C^pX^p+3Y^p\pmod p.
\end{aligned}
\]
Comparing
\[
\operatorname{ord}_U(C),\qquad
\frac{\operatorname{ord}_U(X)}2,\qquad
\frac{\operatorname{ord}_U(Y)}3
\]
isolates one invariant; ties belong to the shorter circuit.

For general support \(\{0,1,-d_1,\ldots,-d_q\}\), the relation monoid is
free on \(X_d=U^dA_dB^d\), of circuit length \(\ell_d=d+1\).  Choose a
nonzero invariant of minimum
\(\operatorname{ord}_U(X_d)/\ell_d\), breaking ties by minimum length.
The \(\ell_dp\)-th moment isolates it after prime reduction and the
factorial functional.  Hence every unit star and its reflection has
one-sided support and satisfies GMC(2).

The literal graph-only leaf-removal statement needs additional semigroup
data.  For weights \(\{5,-2,-3\}\), the pairwise circuit graph is a star,
but the mixed relation \(5-2-3=0\) produces a length-three invariant
invisible in that graph.

There is, however, a stronger way around this obstruction.  For arbitrary
\[
P=\sum_kT^kB_k(U),
\]
plot the points
\[
(k,\operatorname{ord}_U B_k).
\]
The lower convex envelope over weight zero has a supporting line
\[
\operatorname{ord}_U B_k\geq\rho+\theta k.
\]
Every zero-weight monomial of length \(m\) has radial order at least
\(\rho m\).  The contact-face Laurent polynomial has weight zero in its
Newton interval, so the Duistermaat--van der Kallen theorem supplies a
power \(r\) with nonzero constant term \(c\).  At moment \(rp\), Frobenius
and the factorial functional isolate \(c^p\), a contradiction.

This lower-face argument proves that every pure-moment-zero polynomial has
strictly one-sided rotational support.  Thus it proves all of GMC(2), not
only the circuit-forest conjecture.

The proof and precise boundary are in
[`UNIT_STAR_GAUSSIAN_RIGIDITY.md`](UNIT_STAR_GAUSSIAN_RIGIDITY.md); the
general theorem is
[`TWO_REAL_GMC_LOWER_FACE_THEOREM.md`](TWO_REAL_GMC_LOWER_FACE_THEOREM.md).

## 7. Literature connection

The support program sits at the intersection of three existing theorems,
but is not currently implied by any of them.

1. [Derksen--van den Essen--Zhao](https://arxiv.org/abs/1506.05192)
   introduced GMC, proved \(\operatorname{GMC}(1)\), and proved
   \(\operatorname{GMC}(2)\) for homogeneous polynomials.  Their homogeneous
   proof separates radius from angle and applies the torus Mathieu theorem.
   Equation (4.1) identifies why the inhomogeneous case is different: the
   factorial functional couples radial degrees after angular
   constant-term extraction.
2. [Duistermaat--van der Kallen](https://doi.org/10.1016/S0019-3577(98)80020-7)
   classified Laurent polynomials whose positive powers have zero constant
   term; geometrically, vanishing forces the Newton support away from the
   origin.  This is the direct ancestor of the weight-null-cone and
   support-graph viewpoint.  It cannot be applied after (2.1), because
   \(\mathcal L\) can cancel different radial coefficients even when the
   constant term itself is nonzero.
3. [van den Essen--Wright--Zhao](https://arxiv.org/abs/1008.3962)
   introduced the Factorial Conjecture and proved its one-variable case
   (Theorem 4.9).  That result closes the axes \(C=0\) and \(D=0\) in
   (4.1); the prime-endpoint argument in Section 4.1 closes their
   Bessel-weighted coupling.
4. [Long](https://arxiv.org/abs/2607.18186) gives explicit failures in
   every dimension at least three and leaves dimension two as the only
   unresolved dimension.  His examples show that Lagrange/determinant
   cancellation can defeat pure moments, but their auxiliary Gaussian
   source has no analogue in (4.1).  The Bessel-factorial equation is the
   exact two-real-variable replacement problem.

The theorem in Section 4.1 was the first rigidity result interpolating
between the Duistermaat--van der Kallen constant-term theorem and the
factorial functional.  Section 6 upgrades the same prime mechanism to the
complete two-real-variable theorem.

## 8. Alternative structure behind the theorem

The prime proof settles the three-level family, so the following
reformulations are no longer required for rigidity.  They remain useful as
possible templates for genuine stars and circuits, and as explanations of
the analytic structure hidden by the elementary valuation argument.

### 8.1 Laplace--Bessel saddles and the formal-series caveat

On polynomials,

\[
 \mathcal L(f)=\int_0^\infty e^{-u}f(u)\,du.
\]

If

\[
 c=\deg C,\qquad d=\deg D,\qquad q=\max(c,d/2),
\]

the two large-argument Bessel phases are formally

\[
 \Phi_\pm(u,t)=-u+tC(u)\pm2t\sqrt{D(u)}.           \tag{8.1}
\]

For \(q>1\), a generic saddle has

\[
 u\asymp t^{-1/(q-1)},\qquad
 \Phi_\pm(u,t)\asymp t^{-1/(q-1)}.                \tag{8.2}
\]

Thus the leading saddle slopes are separated unless \(d=2c\).  On that
balanced line the first collision occurs precisely when one of the leading
coefficients

\[
 c_c\pm2\sqrt{d_{2c}}
\]

vanishes, equivalently \(c_c^2=4d_{2c}\).  Subsequent collisions are read
from the Newton polygon of

\[
 C(U)^2-4D(U).                                    \tag{8.3}
\]

There is an important analytic qualification.  Equations (4.1)--(4.3) are
identities of formal moment series; in radial degree at least two those
series generally have zero radius of convergence.  All moments vanishing
therefore says that a sectorial Laplace integral has asymptotic expansion
\(1\), not automatically that the analytic integral is identically \(1\).
Nontrivial saddle terms may be exponentially flat at \(t=0\).  A complete
asymptotic proof must consequently establish that every nonzero pair
\((C,D)\) has a nonzero Stokes contribution on some admissible \(t\)-ray.
The Newton classification above gives the finite list of leading collisions
that such a proof must resolve.

### 8.2 Bessel ODE and creative telescoping

Put

\[
 F(t,U)=e^{tC(U)}I_0(2t\sqrt{D(U)}).
\]

The Bessel equation gives the exact differential identity

\[
 \left[
 t^2\partial_t^2+(t-2t^2C)\partial_t
 +t^2(C^2-4D)-tC
 \right]F=0.                                      \tag{8.4}
\]

The factorial functional has the adjoint relation

\[
 \mathcal L(Uf)=\mathcal L((U\partial_U+1)f).      \tag{8.5}
\]

More generally, multiplication by \(U^k\) can be transferred one factor at
a time, producing the rising Euler operator

\[
 (U\partial_U+1)(U\partial_U+2)\cdots
 (U\partial_U+k).                                 \tag{8.6}
\]

For symbolic elimination it is preferable to use the algebraic resolvent

\[
 G(t,U)=Q(t,U)^{-1/2},\qquad
 Q=(1-tC)^2-4t^2D,
\]

which satisfies

\[
 2Q\,\partial_tG+Q_tG=0,\qquad
 2Q\,\partial_UG+Q_UG=0.                          \tag{8.7}
\]

The concrete creative-telescoping target is an operator
\(T(t,\partial_t)\) and an algebraic certificate \(R(t,U)G\) such that

\[
 T(G)=(\partial_U-1)(R\,G)+B(t,U).                \tag{8.8}
\]

After integration against \(e^{-U}dU\), the first term is a boundary term.
Substituting the proposed solution \(\mathcal L(G)=1\) into the resulting
ODE gives algebraic conditions on the coefficients of \(C,D\).  This route
also exposes the irregular singularities where exponentially flat solutions
can occur.  Formal substitution alone does not rule those solutions out;
boundary conditions and Stokes data must still be tracked.  A useful next
computation is to derive the minimal telescoper for generic symbolic degrees
\((c,d)=(2,3)\), then \((3,5)\), and identify the coefficient pattern before
attempting arbitrary degree.

There is already a finite first-order closure before scalar telescoping.
Write

\[
 Q(t,U)=\sum_{j=0}^s q_j(t)U^j,\qquad
 M_k(t)=\mathcal L(U^kG(t,U)).
\]

Integration by parts in the \(U\)-equation of (8.7) gives, for every
\(k\geq0\),

\[
\boxed{
 2\sum_{j=0}^s q_jM_{k+j}
 -\sum_{j=0}^s(2k+j)q_jM_{k+j-1}
 =2\delta_{k,0}Q(t,0)G(t,0),
}                                                  \tag{8.9}
\]

where terms with negative subscripts are omitted.  The coefficient of
\(M_{k+s}\) is \(2q_s\).  Hence, after localizing at the leading coefficient,
every radial moment reduces recursively to

\[
 M_0,\ldots,M_{s-1}.                              \tag{8.10}
\]

Multiplying the \(t\)-equation of (8.7) by \(U^k\) and applying
\(\mathcal L\) gives

\[
 2\sum_{j=0}^s q_jM'_{k+j}
 +\sum_{j=0}^s q'_jM_{k+j}=0.                    \tag{8.11}
\]

Differentiate the reductions (8.9), substitute them into (8.11) for
\(0\leq k<s\), and solve the resulting linear system.  This produces an
explicit meromorphic Pfaffian system

\[
 \partial_t
 \begin{pmatrix}M_0\\ \vdots\\ M_{s-1}\end{pmatrix}
 =
 A(t)
 \begin{pmatrix}M_0\\ \vdots\\ M_{s-1}\end{pmatrix}
 +b(t).                                           \tag{8.12}
\]

The initial data are \(M_k(0)=k!\).  Equivalently, the prime theorem says
that this distinguished formal horizontal section cannot remain in the
affine hyperplane \(M_0=1\) unless \(C=D=0\).
Eliminating the other \(s-1\) coordinates yields the scalar telescoper in
(8.8).  The poles introduced by \(q_s(t)\), which normally vanishes at
\(t=0\), make the expected irregular singularity explicit rather than
hiding it in moment growth.

The exact regression
[`verify_two_real_gmc_resolvent_system.py`](../scripts/verify_two_real_gmc_resolvent_system.py)
constructs the four-dimensional system for a centered degree-\((2,3)\)
pair and checks it against the factorial moment series.

### 8.3 Laguerre operator form

Let \(L_n(U)\) denote the ordinary Laguerre polynomials.  They satisfy

\[
 \int_0^\infty e^{-U}L_n(U)L_m(U)\,dU=\delta_{n,m},
\]

so \(\mathcal L\) selects the \(L_0\)-coefficient.  Multiplication by \(U\)
is the tridiagonal Jacobi operator

\[
 U L_n=(2n+1)L_n-(n+1)L_{n+1}-nL_{n-1}.           \tag{8.13}
\]

If \(J\) denotes this operator and

\[
 H_m(C,D)=
 \sum_{r=0}^{\lfloor m/2\rfloor}
 \frac{m!}{(r!)^2(m-2r)!}C^{m-2r}D^r,
\]

then the moment equations become

\[
 \langle e_0,H_m(C(J),D(J))e_0\rangle=0
 \quad(m\geq1).                                  \tag{8.14}
\]

Here \(C(J)\) and \(D(J)\) are commuting finite-band operators.  This turns
all-degree rigidity into a boundary spectral-moment problem for one fixed
Jacobi matrix.  It also supplies a sparse exact computational basis: the
outermost Laguerre bands and their unique extremal paths can be studied
before any coefficient-chart localization.

For related factorial transforms in higher radial dimension,
(8.7)--(8.8) remains a concrete theorem engine.
The saddle analysis identifies its irregular singularities and exceptional
Newton polygons; the Laguerre form may expose a triangular or
positivity-free band argument.

## 9. The first genuine circuit cycle

Before the general lower-face theorem, the first test beyond the
three-level forest was the \(K_{2,2}\) support

\[
 P=T^2B_2+TB_1+T^{-1}UA_1+T^{-2}U^2A_2.
\]

Its four primitive circuit invariants are

\[
\begin{aligned}
 X_{11}&=UA_1B_1,&
 X_{12}&=U^2A_2B_1^2,\\
 X_{21}&=U^2B_2A_1^2,&
 X_{22}&=U^2A_2B_2.
\end{aligned}                                    \tag{9.1}
\]

They obey the toric cycle relation

\[
 \boxed{X_{12}X_{21}=X_{11}^2X_{22}.}             \tag{9.2}
\]

The first constant terms are already sparse:

\[
\begin{aligned}
 \operatorname{CT}_T(P^2)&=2(X_{11}+X_{22}),\\
 \operatorname{CT}_T(P^3)&=3(X_{12}+X_{21}).
\end{aligned}                                    \tag{9.3}
\]

The prime-valuation method in fact survives this first cycle.

### 9.1 Frobenius before coefficient bounds

Put
\[
 K_m=\operatorname{CT}_T(P^m).
\]
For every prime \(p\) and fixed \(k\), Frobenius in
\(\mathbb F_p[U][T,T^{-1}]\) gives the universal identity
\[
 \boxed{K_{kp}\equiv K_k^p\pmod p.}               \tag{9.4}
\]
This is the correct reduction to make before introducing any radial
coefficient bounds.

Give the four invariants their circuit lengths
\[
 d_{11}=d_{22}=2,\qquad d_{12}=d_{21}=3.
\]
For a nonzero invariant put
\[
 \alpha_{ij}=\operatorname{ord}_U X_{ij},\qquad
 \rho_{ij}=\frac{\alpha_{ij}}{d_{ij}}.             \tag{9.5}
\]
Every monomial of \(K_m\) has weighted degree \(m\), so its radial order is
at least \(m\rho\), where \(\rho=\min\rho_{ij}\).

The same finite-type coefficient reduction and prime-valuation argument as
in Section 4.1 turns (9.4) into the following initial-form rule.  If
\(k\rho\) is an integer, then
\[
 \mathcal L(K_{kp})=0\quad\text{for all sufficiently large primes }p
 \quad\Longrightarrow\quad
 [U^{k\rho}]K_k=0.                                \tag{9.6}
\]
Indeed, every term of \(K_{kp}\) has order at least \(kp\rho\).  After
division by \((kp\rho)!\), the \(p\)-divisible error in (9.4) disappears
modulo \(p\), while every higher term of \(K_k^p\) crosses a new multiple
of \(p\).  Only the \(p\)-th power of the coefficient in (9.6) remains.

### 9.2 Classification of valuation ties

When all four invariants are nonzero, (9.2) gives
\[
 \alpha_{12}+\alpha_{21}=2\alpha_{11}+\alpha_{22},
\]
or, in normalized slopes,
\[
 3\rho_{12}+3\rho_{21}
 =4\rho_{11}+2\rho_{22}.                          \tag{9.7}
\]
Equation (9.7) leaves only three types of minimum:

1. one circuit has a unique minimum;
2. exactly one length-two circuit and one adjacent length-three circuit
   tie for the minimum;
3. all four circuits tie.

Two length-two minima force both length-three slopes to be minimal, and
two length-three minima force both length-two slopes to be minimal.  Any
three-way tie similarly forces the fourth.  Thus the list is exhaustive.

A unique minimum is impossible by (9.6), taking \(k=2\) or \(3\): its
lowest coefficient survives alone.

For an adjacent two-way tie, write \(x\) for the lowest coefficient of the
length-two invariant and \(y\) for that of the length-three invariant.
The initial equations at \(k=6\) and \(k=12\) are
\[
 20x^3+15y^2=0,                                   \tag{9.8}
\]
and one of
\[
\begin{aligned}
 924x^6+7920x^3y^2+495y^4&=0,\\
 924x^6+27720x^3y^2+495y^4&=0.                   \tag{9.9}
\end{aligned}
\]
The first line occurs for \(X_{11}\), the second for \(X_{22}\).
Substituting \(y^2=-4x^3/3\) into (9.9) gives respectively
\[
 -8756x^6=0,\qquad -35156x^6=0.                  \tag{9.10}
\]
Hence \(x=y=0\), contradicting the definition of the tied face.

For the four-way tie, let \(a,b,c,d\) be the lowest coefficients of
\(X_{11},X_{12},X_{21},X_{22}\).  Equations (9.3), the fourth moment, and
the initial form of (9.2) give
\[
\begin{aligned}
 a+d&=0,& b+c&=0,\\
 a^2+4ad+d^2&=0,& bc&=a^2d.                      \tag{9.11}
\end{aligned}
\]
The first and third equations give \(a=d=0\); the remaining two give
\(b=c=0\), again impossible.

Boundary faces of the toric hypersurface cause no exception.  Their
nonzero supports consist of single circuits or pairs.  Adjacent pairs are
covered by (9.8)--(9.10).  A tied pair \(X_{11},X_{22}\) is killed by
\(K_2,K_4\), and a tied pair \(X_{12},X_{21}\) by \(K_3,K_6\).

Consequently:

> **First-cycle rigidity theorem.** Let
> \[
> P=T^2B_2(U)+TB_1(U)+T^{-1}UA_1(U)+T^{-2}U^2A_2(U).
> \]
> If \(\mathcal L(\operatorname{CT}_T P^m)=0\) for every \(m\geq1\), then
> \[
> X_{11}=X_{12}=X_{21}=X_{22}=0.
> \]
> Hence all nonzero rotational weights of \(P\) have the same sign, and
> \(P\) satisfies \(\operatorname{GMC}(2)\).

Thus the first toric cycle produces no surviving valuation cone.  The
prime method extends beyond forests; the next obstruction must involve a
larger circuit semigroup whose minimal face is not separated by finitely
many low invariant moments.

The exact regression
[`verify_two_real_gmc_first_cycle_rigidity.py`](../scripts/verify_two_real_gmc_first_cycle_rigidity.py)
enumerates the toric moments, verifies (9.4), and checks every tied-face
elimination above.
