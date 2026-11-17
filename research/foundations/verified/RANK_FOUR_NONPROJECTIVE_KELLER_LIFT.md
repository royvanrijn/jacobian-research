# Rank-four nonprojective Keller descent: Kummer class, formal lift, and target-symmetry obstructions

## 0. Result and scope

The rank-four projective-descent theorem leaves the first genuinely
nonprojective primitive change open.  This note applies that change to the
quartic quadratic-gauge Keller family and proves six exact statements.

1. Over an arbitrary characteristic-zero field, the rank-four
   quadratic-gauge orbit has a precise Kummer descent class in
   \(K^\times/K^{\times5}\).  This separates arithmetic failure of a
   ground-field equivalence from geometric failure of projective root
   transport.
2. There is a rational primitive nonprojective witness for which the Kummer
   class vanishes.  After one explicit rational left--right scaling, the
   problem becomes a comparison of two regular fibers of the **same fixed
   map** \(F_{-124416}\).
3. The desired marked-root motion has an exact divergence-free polynomial
   first-order lift.  The repository's formal-orbit theorem then supplies a
   unique lift to every finite order, and every finite jet has a polynomial
   `SAut(A^3)` representative.
4. The straight target translation cannot have a polynomial source lift:
   its \(-4\) iterate takes a four-point fiber to an exact two-point fiber.
5. More generally, every target component of degree at most twelve is one
   of the residual \(\mu_5\) symmetries; none carries the two endpoint
   targets.
6. The necessary logarithmic-boundary and constant-Jacobian equations have
   no endpoint solution through target degree eighteen.

This does **not** construct a polynomial automorphism at the endpoint.  The
remaining question is whether a different, nontranslation target symmetry
of degree at least nineteen can realize the prescribed permutation of the
four collision-frame sheets.
Formal deformation theory, the finite-etale conormal module, and the coarse
quartic stable invariant cannot decide that question.

## 1. Credit and provenance

The collision algebra, diagonal/off-diagonal terminology, and its
finite-etale base-change interface used here are credited to Chloe van der
Vlugt's *Collision Ideals and Off-Diagonal Sheets*, with public Lean source
published by the GitHub account `what-social-construct`.  The exact
manuscript/account distinction, pinned commit, licenses, disclosed Codex
assistance, and literature credits are recorded in the
[external audit](COLLISION_IDEALS_EXTERNAL_AUDIT.md).

The all-rank frame criterion, the Kummer calculation below, the neutral
witness, and the explicit Keller lifts are repository derivations.  They are
not attributed to the external paper.  In particular, that paper supplies
the collision interface but does not claim the rank-four cross-ratio theorem
or the nonprojective lift developed here.

The mathematical inputs internal to this repository are:

- [rank-four collision framing and cross ratio](RANK_FOUR_COLLISION_CROSS_RATIO.md);
- [all-rank projective Tschirnhaus descent](ALL_RANK_COLLISION_PROJECTIVE_DESCENT.md);
- [quadratic-gauge stable moduli](QUADRATIC_GAUGE_STABLE_MODULI.md);
- [formal polynomial-orbit triviality](../extended-geometry/FORMAL_ORBIT_TRIVIALITY.md);
- [intrinsic quartic torus exclusion](../cancellation/NO_ALGEBRAIC_TORUS_EQUIVARIANCE.md).

## 2. The fixed quartic Keller family

Put

\[
 d=1+xy,\qquad
 Q=d^2z+y^2(1+3d).
\]

For \(U\ne0\), define

\[
\begin{aligned}
 P_U&=dQ,\\
 B_U&=-\frac12\left(y+3xQ+4U d^2x^2Q^4\right),\\
 C_U&=x(5-3d)-x^3z-2U(xQ)^4,
\end{aligned}
\]

and write

\[
 F_U=(P_U,B_U,C_U):\mathbb A^3\longrightarrow\mathbb A^3.
 \tag{2.1}
\]

The exact Jacobian identity is

\[
 \boxed{\det DF_U=1.}                                  \tag{2.2}
\]

At a target \(y=(\pi,b,c)\), its inverse-root equation is

\[
 E_{U,y}(S)
 =U\pi^4S^4+\pi S^3+bS^2+S-\frac c2.                 \tag{2.3}
\]

If (2.3) is separable, the fiber is scheme-theoretically

\[
 F_U^{-1}(y)\simeq
 \operatorname{Spec}K[S]/(E_{U,y}).
 \tag{2.4}
\]

Thus a normalized quartic relation

\[
 E(S)=a_0+S+a_2S^2+a_3S^3+a_4S^4,\qquad a_3a_4\ne0,  \tag{2.5}
\]

is compiled by

\[
 \boxed{
 U=\frac{a_4}{a_3^4},\qquad
 y=(a_3,a_2,-2a_0).
 }                                                       \tag{2.6}
\]

## 3. General rank-four ground-field reduction

Let \(K\) be a characteristic-zero field.  For coefficient pairs
\((a_3,a_4),(b_3,b_4)\in(K^\times)^2\), the exact two-torus action from the
stable-moduli theorem is

\[
 b_3=\alpha^{-2}\beta^{-1}a_3,\qquad
 b_4=\alpha^{-3}\beta^{-4}a_4.                        \tag{3.1}
\]

Eliminating \(\beta\) gives

\[
 \alpha^5
 =
 \frac{b_4/b_3^4}{a_4/a_3^4}.                         \tag{3.2}
\]

Consequently two quartic quadratic-gauge maps are stably polynomially
left--right equivalent **over \(K\)** exactly when

\[
 \boxed{
 \kappa(a,b):=
 \frac{b_4/b_3^4}{a_4/a_3^4}
 \in K^{\times5}.
 }                                                       \tag{3.3}
\]

The sufficiency is the displayed polynomial scaling.  Necessity follows
from the intrinsic normalization/Fitting reconstruction in the
stable-moduli theorem: a \(K\)-defined equivalence induces the ordered torus
coordinate change over \(K\), so its \(\alpha,\beta\) lie in \(K^\times\).
Equivalently, (3.3) is the \(K^\times/K^{\times5}\) form of the residual
\(\mu_5\) stack stabilizer.  Over an algebraically closed field it vanishes,
which explains why the coarse geometric quartic quotient is one point.

For the compressed maps (2.1), \(a_3=b_3=1\).  If

\[
 U'=\gamma^5U,
\]

then

\[
 \sigma_\gamma(x,y,z)
 =(\gamma x,\gamma^{-1}y,\gamma^{-2}z),               \tag{3.4}
\]

\[
 \tau_\gamma(P,B,C)
 =(\gamma^{-2}P,\gamma^{-1}B,\gamma C)                \tag{3.5}
\]

satisfy the polynomial identity

\[
 \boxed{
 F_{U'}\circ\sigma_\gamma
 =\tau_\gamma\circ F_U.
 }                                                       \tag{3.6}
\]

This gives a general decision procedure for two normalized quartic
presentations:

1. compute \(U=a_4/a_3^4\) and \(U'=a'_4/(a'_3)^4\);
2. test the Kummer class \(U'/U\);
3. if it is a fifth power, use (3.4)--(3.6) and reduce the comparison to two
   marked fibers of one fixed map;
4. only then test the genuinely geometric, framed fiber transport.

The Kummer class is an arithmetic obstruction to a ground-field
quadratic-gauge equivalence.  It is not a new geometric obstruction to the
Jacobian problem.

## 4. Why the first quadratic witness mixes two effects

Take the ordered roots

\[
 r=(1,2,3,4),\qquad u=r+r^2=(2,6,12,20).
\]

The two normalized relations are

\[
 E_r(S)
 =-\frac{S^4}{50}+\frac{S^3}{5}
  -\frac{7S^2}{10}+S-\frac{12}{25},
\]

\[
 E_u(S)
 =-\frac{S^4}{2304}+\frac{5S^3}{288}
  -\frac{127S^2}{576}+S-\frac54.
\]

Their compressed data are

\[
\begin{array}{c|c|c}
 &U&(\pi,b,c)\\ \hline
 r&-25/2&(1/5,-7/10,24/25)\\
 u&-2985984/625&(5/288,-127/576,5/2).
\end{array}
\]

The four columns \(1,r,u,ru\) have rank four, so the change is genuinely
nonprojective.  But

\[
 \frac{U_u}{U_r}
 =\frac{5971968}{15625}
 =2^{13}3^6 5^{-6}                                   \tag{4.1}
\]

is not a rational fifth power.  This witness therefore combines:

- the geometric rank-four projective defect; and
- a separate rational \(\mu_5\)-twist between the two ambient maps.

It is unsuitable as the cleanest test of nonprojective Keller descent over
\(\mathbb Q\), because a negative result could be caused by (4.1) alone.

## 5. Arithmetic-neutral nonprojective witness

Instead take

\[
 r=(-6,-3,0,3),\qquad
 u=r+r^2-18=(12,-12,-18,-6).                          \tag{5.1}
\]

The values of \(u\) are distinct, so it is a primitive coordinate on the
split rank-four algebra.  The matrix with columns \(1,r,u,ru\) again has
rank four.  Hence (5.1) is primitive and nonprojective.

The normalized relations are

\[
 E_0(S)
 =-\frac{S^4}{54}-\frac{S^3}{9}+\frac{S^2}{6}+S,
 \tag{5.2}
\]

\[
 E_1(S)
 =-\frac{S^4}{3456}-\frac{S^3}{144}
  +\frac{S^2}{96}+S+\frac92.                          \tag{5.3}
\]

Their compressed data are

\[
\begin{array}{c|c|c}
 &U&(\pi,b,c)\\ \hline
 E_0&-243/2&(-1/9,1/6,0)\\
 E_1&-124416&(-1/144,1/96,-9).
\end{array}                                           \tag{5.4}
\]

Now

\[
 \frac{-124416}{-243/2}=1024=4^5.                    \tag{5.5}
\]

With

\[
 \sigma(x,y,z)=(4x,y/4,z/16),\qquad
 \tau(P,B,C)=(P/16,B/4,4C),                           \tag{5.6}
\]

equation (3.6) becomes

\[
 \boxed{
 F_{-124416}\circ\sigma
 =\tau\circ F_{-243/2}.
 }                                                       \tag{5.7}
\]

The old target moves to

\[
 y_{\rm lin}
 =\tau(-1/9,1/6,0)
 =(-1/144,1/24,0),                                    \tag{5.8}
\]

while the desired target is

\[
 y_1=(-1/144,1/96,-9).                                \tag{5.9}
\]

Thus the residual is

\[
 \boxed{\delta=y_1-y_{\rm lin}=(0,-1/32,-9).}         \tag{5.10}
\]

The original roots scale to

\[
 s=4r=(-24,-12,0,12),
\]

and the requested collision-frame coordinate is

\[
 \boxed{
 u=\frac{s^2}{16}+\frac s4-18.
 }                                                       \tag{5.11}
\]

Equations (5.7)--(5.11) are the promised fixed-map reduction.  The
arithmetic ambiguity and the unmarked ambient-map moduli have been removed:
the remaining problem lies entirely inside \(F_{-124416}\).

## 6. A finite-etale target line, and why it is not yet the framed lift

Join the two fixed-map targets by

\[
 y_\lambda
 =y_{\rm lin}+\lambda\delta
 =\left(
 -\frac1{144},
 \frac1{24}-\frac{\lambda}{32},
 -9\lambda
 \right).
 \tag{6.1}
\]

Substitution into (2.3) gives the exact factorization

\[
 \boxed{
 E_\lambda(S)
 =-\frac{(S-12)(S+12)(S^2+24S+108\lambda)}{3456}.
 }                                                       \tag{6.2}
\]

Its discriminant is

\[
 \boxed{
 \operatorname{disc}_S(E_\lambda)
 =-\frac{(\lambda+4)^2(3\lambda-4)^3}{1358954496}.
 }                                                       \tag{6.3}
\]

Therefore (6.2) is a rank-four finite-etale algebra over

\[
 \mathbb Q\left[
 \lambda,\frac1{(\lambda+4)(3\lambda-4)}
 \right],
\]

and both \(\lambda=0\) and \(\lambda=1\) are regular.  This proves that the
two unmarked endpoint fibers occur in one explicit finite-etale family of
the fixed Keller map.

For the constant target velocity \(\delta\), put

\[
 W=(DF_{-124416})^{-1}\delta.
 \tag{6.4}
\]

Because the Jacobian determinant is one, \(W\) is polynomial.  Exact
expansion gives

\[
 \operatorname{div}W=0,\qquad
 \deg W=(31,29,31),\qquad
 \#\operatorname{supp}W=(118,115,149).                \tag{6.5}
\]

The formal target-lift theorem gives a unique formal source automorphism
\(\widehat A_\lambda\), based at the identity, such that

\[
 F_{-124416}\circ\widehat A_\lambda
 =
 T_{\lambda\delta}\circ F_{-124416}.                  \tag{6.6}
\]

However, (6.2) has two constant sections \(S=\pm12\) and one quadratic
sheet.  The desired map (5.11) sends the constant-start set

\[
 \{-12,12\}\longmapsto\{-12,-6\},
\]

so it crosses the constant/quadratic decomposition.  An algebraization of
the based straight-line lift (6.6) would preserve that relative
decomposition and would not, by itself, realize the required collision
frame.  This is the exact distinction between unmarked fiber motion and
nonprojective framed descent.

### 6.1 A global fiber-orbit obstruction

There is a stronger obstruction to algebraizing the straight translation.
For any morphism \(F:X\to Y\), if automorphisms \(A\) and \(B\) satisfy

\[
 F\circ A=B\circ F,                                   \tag{6.7}
\]

then for every \(m\in\mathbb Z\), \(A^m\) restricts to a
scheme-theoretic isomorphism

\[
 F^{-1}(y)\simeq F^{-1}(B^m y).                       \tag{6.8}
\]

Thus geometric point count, fiber length when finite, and every
scheme-theoretic fiber invariant are constant along the full
\(\mathbb Z\)-orbit of \(B\).  This elementary observation is a useful
general obstruction to lifting a target automorphism through a Keller map.

Apply it to \(B=T_\delta\), \(y=y_{\rm lin}\), and \(m=-4\).  Formula (6.2)
becomes

\[
 E_{-4}(S)
 =-\frac{(S-12)^2(S+12)(S+36)}{3456}.                \tag{6.9}
\]

For \(\pi\ne0\), every affine source point has

\[
 E'(S)=\frac1{1+xy}\ne0,                              \tag{6.10}
\]

and every simple inverse root reconstructs one source point.  Hence the
repeated root \(S=12\) in (6.9) is missing from the affine fiber, while the
simple roots \(-36,-12\) reconstruct exactly

\[
 \left(-\frac94,\frac5{12},-\frac{731}{9}\right),
 \qquad
 \left(3,-\frac5{12},-\frac14\right).                 \tag{6.11}
\]

Therefore

\[
 \#F_{-124416}^{-1}(y_{\rm lin})=4,\qquad
 \#F_{-124416}^{-1}(y_{-4})=2.                        \tag{6.12}
\]

Equations (6.8) and (6.12) are incompatible.  Consequently

\[
 \boxed{
 \nexists A\in\operatorname{Aut}(\mathbb A^3):
 F_{-124416}\circ A=T_\delta\circ F_{-124416}.
 }                                                       \tag{6.13}
\]

This proves that the formal target-translation lift (6.6) does not
algebraize, even without imposing the desired frame permutation.  It does
not exclude a different polynomial target automorphism taking
\(y_{\rm lin}\) to \(y_1\).

### 6.2 Exact low-degree target self-equivalence group

The intrinsic decorated-normalization theorem gives a second general
reduction.  For every compressed quartic \(F_U\) over an algebraically
closed characteristic-zero field, the residual symmetries

\[
\begin{aligned}
 \sigma_\zeta(x,y,z)
   &=(\zeta x,\zeta^{-1}y,\zeta^{-2}z),\\
 \tau_\zeta(P,B,C)
   &=(\zeta^{-2}P,\zeta^{-1}B,\zeta C),
\end{aligned}
\qquad \zeta^5=1,                                    \tag{6.14}
\]

satisfy

\[
 F_U\circ\sigma_\zeta=\tau_\zeta\circ F_U.            \tag{6.15}
\]

They are not merely examples among affine-target symmetries.

**Affine-target theorem.**  If polynomial automorphisms \(A,T\) satisfy

\[
 F_U\circ A=T\circ F_U
\]

and \(T\) is affine, then

\[
 (A,T)=(\sigma_\zeta,\tau_\zeta)
\quad\text{for a unique }\zeta\in\mu_5.               \tag{6.16}
\]

Indeed, functoriality sends every left--right self-equivalence to an
automorphism of the intrinsic decorated ramified-normalization stratum.
That automorphism group is exactly \(\mu_5\).  Compose with the inverse of
the corresponding pair in (6.14).  The resulting affine target
automorphism restricts to the identity on the normalization and hence fixes
the dense prime discriminant hypersurface pointwise.  Each coordinate
difference from the identity is affine-linear and vanishes on that
nonlinear hypersurface, so every difference is zero.  The remaining source
map is target-fixed.  Generic \(S_4\) monodromy makes its deck group

\[
 N_{S_4}(S_3)/S_3=1,
\]

so the source map is also the identity.  This proves (6.16).  Over a
nonclosed field \(K\), the \(K\)-defined affine-target group is
\(\mu_5(K)\).

The same argument gives much more than the affine case.  Write the prime
ramified discriminant component as

\[
 \operatorname{disc}_S(E_{U,(P,B,C)})
 =-\frac{P^2}{4}H_U(P,B,C).                           \tag{6.17}
\]

Its defining polynomial has ordinary total degree thirteen, with unique
top-degree term

\[
 \operatorname{in}_{13}(H_U)
 =128U^3P^{10}C^3.                                   \tag{6.18}
\]

For an arbitrary polynomial target component \(T\), compose off its
decorated \(\tau_\zeta\) action as above.  The resulting map fixes the prime
hypersurface \(H_U=0\) pointwise.  Hence every coordinate satisfies

\[
 T_i-(\tau_\zeta)_i\in(H_U).                          \tag{6.19}
\]

If one difference is nonzero, its degree is at least
\(\deg H_U=13\).  If all are zero, deck rigidity recovers
\((A,T)=(\sigma_\zeta,\tau_\zeta)\).  Therefore

\[
 \boxed{
 T=\tau_\zeta\text{ for some }\zeta^5=1
 \quad\text{or}\quad
 \deg T\ge13.
 }                                                       \tag{6.20}
\]

For the present map,

\[
 \frac{-124416}{1/2}=(-12)^5,
\]

so it is rationally in the same quartic class as the small reference map in
the decorated-normalization theorem.  If
\(\tau_\zeta(y_{\rm lin})=y_1\), equality of their common nonzero first
coordinate forces \(\zeta^2=1\).  Together with \(\zeta^5=1\), this gives
\(\zeta=1\), but then the second and third coordinates do not move.
Therefore:

\[
 \boxed{\text{any target symmetry carrying }y_{\rm lin}\text{ to }y_1
 \text{ has polynomial degree at least thirteen.}}     \tag{6.21}
\]

This does not classify target automorphisms of degree at least thirteen in
the kernel of the decorated-boundary action.

### 6.3 Exact exclusion through target degree eighteen

The first degrees allowed by (6.21) reduce to a finite logarithmic system.
After composing off the residual \(\mu_5\) action, write

\[
 T=\operatorname{id}+H_UV.                            \tag{6.22}
\]

If \(\deg T\le18\), then \(\deg V\le5\).  Since an automorphism preserving
the prime divisor has \(T^*H_U=\rho H_U\), Taylor expansion modulo \(H_U^2\)
gives the necessary identity

\[
 V(H_U)-QH_U-\kappa=0,\qquad
 \deg Q\le\deg V-1.                                   \tag{6.23}
\]

Exact rational coefficient matrices for (6.23) have nullities

\[
\begin{array}{c|rrrrrr}
\deg V\le m&0&1&2&3&4&5\\ \hline
\dim\mathcal L_m&0&0&0&0&1&7.
\end{array}                                           \tag{6.24}
\]

Every solution has \(\kappa=0\).  The endpoint equation is

\[
 H_U(y_{\rm lin})V(y_{\rm lin})=\delta,\qquad
 H_U(y_{\rm lin})=-\frac1{16}.                        \tag{6.25}
\]

The evaluation map from the seven-dimensional space in (6.24) has rank
three, so (6.25) leaves exactly four affine parameters.

For such a candidate,

\[
 DT
 =I+V(\nabla H_U)^T+H_UDV.                           \tag{6.26}
\]

The equality \(\kappa=0\) makes the normal action along \(H_U=0\) the
identity.  Thus a polynomial automorphism would have
\(\det DT=1\).  Evaluate \(\det DT-1\) at the ten target points

\[
\begin{gathered}
(1,1,0),(1,0,1),(1,1,1),(-1,1,0),(2,0,0),\\
(1,-1,2),(2,1,-1),(-1,2,1),(2,-1,1),(-2,1,2).
\end{gathered}                                        \tag{6.27}
\]

After primitive denominator clearing, these are ten cubic polynomials in
the four remaining parameters.  Their exact characteristic-zero Gröbner
basis in Singular is

\[
 \boxed{[1].}                                         \tag{6.28}
\]

Every genuine constant-Jacobian target automorphism would vanish on all ten
evaluations, so (6.28) is an exact contradiction.  Combining this with the
\(\mu_5\) endpoint failure proves

\[
 \boxed{\text{any target symmetry realizing the endpoint frame has }
 \deg T\ge19.}                                        \tag{6.29}
\]

This is an exact finite computer-algebra proof over \(\mathbb Q\), not a
random or bounded-coefficient search.  It uses only necessary equations, so
their inconsistency is sufficient.  No claim is made in target degree
nineteen or above.

## 7. The framed path and its polynomial first-order lift

To retain the desired labels, use

\[
 u_\theta(r)=r+\theta(r^2-18),\qquad
 r\in\{-6,-3,0,3\}.                                   \tag{7.1}
\]

The four values are

\[
 (18\theta-6,-9\theta-3,-18\theta,3-9\theta).
\]

Their Vandermonde is

\[
 8748(3\theta-1)^2(3\theta+1)(6\theta-1)(9\theta-1).
 \tag{7.2}
\]

After normalizing the linear coefficient to one, the relation coefficients
are

\[
\begin{aligned}
 a_0(\theta)
 &=\frac{18\theta(3\theta-1)(3\theta+1)}
         {36\theta^2-3\theta-1},\\
 a_2(\theta)
 &=\frac{27\theta^2-24\theta+1}
 {6(3\theta-1)(36\theta^2-3\theta-1)},\\
 a_3(\theta)
 &=-\frac{3\theta+1}
 {9(3\theta-1)(36\theta^2-3\theta-1)},\\
 a_4(\theta)
 &=-\frac1
 {54(3\theta-1)(36\theta^2-3\theta-1)}.
\end{aligned}                                         \tag{7.3}
\]

The compressed seed is

\[
 U(\theta)
 =-\frac{
 243(3\theta-1)^3(36\theta^2-3\theta-1)^3
 }{2(3\theta+1)^4}.                                   \tag{7.4}
\]

The open obtained by removing the factors in (7.2), the linear
normalization factor \(36\theta^2-3\theta-1\), and \(3\theta+1\) contains
both endpoints.  At those endpoints,

\[
\begin{aligned}
 U(0)&=-243/2,&
 y(0)&=(-1/9,1/6,0),\\
 U(1)&=-124416,&
 y(1)&=(-1/144,1/96,-9).
\end{aligned}                                         \tag{7.5}
\]

The first derivatives at zero are

\[
 \dot U(0)=1458,\qquad
 \dot y(0)=(-1/3,-4,-36).                             \tag{7.6}
\]

Put \(F_0=F_{-243/2}\), and define

\[
 X=(DF_0)^{-1}
 \left(\dot y(0)-1458\,\partial_UF_U|_{U=-243/2}\right).
 \tag{7.7}
\]

The determinant-one identity makes \(X\) polynomial.  Exact expansion gives

\[
 \boxed{
 \operatorname{div}X=0,\quad
 \deg X=(55,53,55),\quad
 \#\operatorname{supp}X=(512,510,612).
 }                                                       \tag{7.8}
\]

The four source points reconstructed from (5.2) are

\[
\begin{array}{c|c|c}
 r&p_r&X(p_r)\\ \hline
 -6&(-2,1/3,-5)&(-30,-8,303)\\
 -3&(3,-2/3,1)&(45,6,-25)\\
 0&(0,-1/3,-5/9)&(-18,2,31/3)\\
 3&(-1,4/3,3)&(3,4,9).
\end{array}                                           \tag{7.9}
\]

For the intrinsic root coordinate \(S=x/(1+xy)\),

\[
 dS_{p_r}(X)=r^2-18.                                  \tag{7.10}
\]

Thus \(X\) induces exactly the derivative of the nonprojective labeled
motion (7.1), not merely an unmarked motion of the four-point set.

For the translated family

\[
 \mathcal G_\theta
 =F_{U(\theta)}-y(\theta)+y(0),
\]

the formal-orbit theorem gives a unique formal determinant-one source
automorphism \(\widehat\alpha_\theta\) satisfying

\[
 \mathcal G_\theta=F_0\circ\widehat\alpha_\theta.
 \tag{7.11}
\]

The inverse formal automorphism transports the marked fiber, and its first
coefficient is \(X\).  In dimension three, every finite jet of this formal
transport is represented by a polynomial curve in `SAut(A^3)`.  This proves
all finite orders separately.  It does not prove convergence, termination,
rationality in \(\theta\), or a polynomial automorphism at \(\theta=1\).

## 8. What this buys us

The former rank-four question mixed four issues.  They are now separated:

\[
\begin{array}{c|c}
\text{issue}&\text{status}\\ \hline
S_4\text{ labeling}&\text{resolved by the collision frame}\\
\text{projective root transport}&\text{classified by rank}(1,r,u,ru)\le3\\
\mathbb Q\text{-descent of the ambient quartic map}
  &\text{classified by }U'/U\in\mathbb Q^{\times5}\\
\text{finite-order source transport}
  &\text{exists uniquely and is polynomial at every jet}\\
\text{straight target-translation lift}
  &\text{globally impossible by the }-4\text{ fiber drop}\\
\text{target degree at most }12
  &\mu_5\text{ exactly; none carries the endpoints}\\
\text{endpoint target degrees }13\text{ through }18
  &\text{excluded by logarithmic/Jacobian equations}\\
\text{target degree at least }19\text{ with the }q\text{-frame}
  &\text{open}.
\end{array}
\]

In particular:

- the quartic coarse stable quotient cannot obstruct the neutral witness,
  because its Kummer class is trivial;
- finite-etale local and conormal data cannot obstruct it;
- unrestricted infinitesimal or Artin deformation theory cannot obstruct
  it, by formal orbit triviality;
- the straight target line connects the unmarked fibers but has the wrong
  relative sheet decomposition;
- its target translation has no polynomial lift at all, because an integer
  iterate changes the fiber cardinality from four to two;
- every target self-equivalence through degree twelve is in \(\mu_5\), and
  its endpoint orbit test fails;
- the exact logarithmic-boundary frontier and ten-point Jacobian unit ideal
  exclude endpoint target degrees thirteen through eighteen;
- the obstruction, if one exists, must be global: boundary algebraization,
  filtered polynomial complexity, or a discrete decorated-normalization
  constraint.

The intrinsic torus theorem says that the canonical boundary decoration of
the small quartic map has automorphism group scheme \(\mu_5\) and no
infinitesimal automorphisms.  Section 6.2 upgrades this to a complete
classification through target degree twelve.  It still does not classify
degree-at-least-thirteen discrete or unipotent polynomial
self-equivalences in general, but Section 6.3 excludes their endpoint
problem through degree eighteen.  The degree-at-least-nineteen
possibilities are exactly what remains relevant here.

The earlier
[root-changing suspension audit](../cancellation/ESCAPE_THREE_SUSPENSION_FAMILIES.md)
eliminates root-only Möbius recharts by an unavoidable order-two boundary
pole.  The transformation (5.11) is generically quadratic, not birational on
the root line, so any global lift must use the finite collision frame, mix
the root with \(P\), or introduce an additional primitive reconstruction
coordinate.  Enlarging the root-only Möbius ansatz cannot reach this case.

## 9. Exact next problem

The straight translation is now excluded by (6.13).  The minimal remaining
endpoint problem is:

> Determine whether there are polynomial automorphisms
> \(A,T\in\operatorname{Aut}_{\mathbb Q}(\mathbb A^3)\) such that
> \[
> F_{-124416}\circ A=T\circ F_{-124416},\qquad
> T(y_{\rm lin})=y_1,
> \]
> and the induced fiber isomorphism realizes (5.11), or prove that every
> such degree-at-least-nineteen target symmetry violates an intrinsic boundary
> divisor.

Any candidate \(T\) must preserve every fiber invariant along its complete
integer orbit by (6.8).  It must also induce an allowed automorphism of the
decorated normalization, while (6.29) forces \(\deg T\ge19\).  A useful
bounded search should impose those orbit and boundary constraints, together
with the endpoint permutation, before expanding a source automorphism
ansatz.

## 10. Exact regression

Run

```bash
.venv/bin/python scripts/verify_rank_four_nonprojective_keller_lift.py
```

The checker verifies:

- the fifth-power orbit elimination;
- both explicit quartic presentations and their projective defects;
- the nontrivial and neutral Kummer ratios;
- the exact source/target scaling (5.7);
- both endpoint fibers and the residual target translation;
- the factorization and discriminant of the fixed-map target line;
- the two-point fiber at \(\lambda=-4\) and the resulting global
  translation-lift obstruction;
- the degree-thirteen prime discriminant and residual \(\mu_5\) orbit test
  used by the low-degree target theorem;
- the framed rational path and its endpoint data;
- \(\det DF_U=1\);
- exact reconstruction of all four source points;
- both polynomial first-order lifts, their divergences, degrees, term
  counts, and induced root velocities.

This first checker deliberately ends with a scope line stating that no
degree-at-least-thirteen endpoint self-equivalence is claimed.

The exact degree-eighteen continuation requires Singular:

```bash
.venv/bin/python scripts/verify_rank_four_degree_eighteen_target_obstruction.py
```

It reconstructs the logarithmic systems through multiplier degree five,
checks the nullities `(0,0,0,0,1,7)`, imposes the endpoint condition, and
generates ten exact determinant evaluations in the four surviving
parameters.  Singular returns the reduced basis `[1]` over `QQ`.  Its scope
line leaves target degree nineteen and above open.
