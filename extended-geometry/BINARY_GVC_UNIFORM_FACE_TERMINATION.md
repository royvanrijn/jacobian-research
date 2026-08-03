# Uniform Hall localization and terminal-face rigidity for binary GVC

## 1. Scope

Let \(k\) be a characteristic-zero field,
\[
 \Lambda=\lambda(\partial_x,\partial_y),\qquad P,Q\in k[x,y],
\]
and suppose
\[
 \Lambda^m(P^m)=0\qquad(m\geq1).
 \tag{1.1}
\]
Write \(r\) for the lowest positive order of \(\Lambda\) and \(d=\deg P\).
The results below are uniform in \(r\) and \(d\).

They prove three uniform pieces of the proposed binary termination
argument:

1. every leading pure-zero component is localized by one root
   multiplicity, with no enumeration of global binary root partitions; and
2. every genuine unequal-weight common-threshold face is automatically
   terminal; and
3. finite-character digit separation splits every genuinely
   scale-compatible torsion--torus trace.

Section 7 first bounds the circuit support, closes every terminal packet,
and records the small-prime obstructions to simpler exposure arguments.
The final proved step uses repeated equal base-\(p\) digits at completely
split primes.  Those digit rows recover every power sum of the finitely
many character-component moments, so Newton identities split a
scale-compatible trace componentwise.

The later consecutive-residue incidence theorem further proves that, once
such a fixed marked packet exists, every nonfree factorization ambiguity
reduces to already-safe beta circuits in every radial span.  Thus there is no
remaining post-promotion semigroup-classification hypothesis.

This does **not** yet prove unrestricted \(\operatorname{GVC}(2)\).
The Hall--jet valuation produces affine, prime-dependent carry shells
before it produces a trace family.  What remains is a promotion theorem
showing that every proportional-depth shell either becomes a fixed
scale-compatible family, is Hall-terminal, or loses support.  Section 7
gives an exact weighted-trace classification and a small factorial
counterexample which show why that promotion is a genuine extra step.

## 2. Uniform Hall localization

Let \(\lambda_r\) and \(P_d\) be the indicated homogeneous pieces.  After
scalar extension, write
\[
 \lambda_r=\prod_{i=1}^r D_{v_i},\qquad
 P_d=\prod_{j=1}^d L_j.
 \tag{2.1}
\]
The translated split-symbol construction and the
Duistermaat--van der Kallen theorem say that the \(r\) derivative copies
cannot be matched to \(r\) distinct polynomial factors on which they act
nontrivially.

> **Theorem 2.1 (uniform binary Hall localization).**  Assume \(d>r\).
> There is a root direction of \(\lambda_r\), of multiplicity
> \(1\leq\mu\leq r\), such that after a linear change of coordinates
> \[
>  M_\mu=X^\mu Y^{r-\mu},\qquad
>  P_d=y^{\,d-\mu+1}C_{\mu-1}(x,y).
>  \tag{2.2}
> \]

Indeed, a Hall-deficient subset containing two nonparallel derivative
directions sees every one of the \(d\) polynomial factors: a nonzero binary
linear form cannot annihilate two independent vectors.  Hence a deficient
subset consists of the \(\mu\) copies of one direction.  If \(c\) polynomial
factors annihilate it, Hall failure is
\[
 d-c<\mu,\qquad\text{equivalently}\qquad c\geq d-\mu+1,
 \tag{2.3}
\]
which is exactly (2.2).  Conversely, (2.2) gives the displayed Hall
deficiency.

Local division by \(M_\mu\) gives a second uniform restriction.  A
normalized operator monomial of excess order \(h\),
\[
 X^aY^{r+h-a},
 \tag{2.4}
\]
is supported in one of the two wings
\[
 a<\mu\qquad\text{or}\qquad a>\mu+h.
 \tag{2.5}
\]
The removed differential unit is invertible and locally finite on
polynomials.  Formula (2.5), rather than the global root partition of
\(\lambda_r\), is the finite support input for every local Newton face.

## 3. Unequal-weight face rigidity

Fix distinct positive integers \(u,v\), and put
\[
 \operatorname{wt}(x)=\operatorname{wt}(X)=u,\qquad
 \operatorname{wt}(y)=\operatorname{wt}(Y)=v.
 \tag{3.1}
\]

> **Theorem 3.1 (unequal-weight terminal-face theorem).**  Let
> \(A(X,Y)\) and \(B(x,y)\) be nonzero weighted-homogeneous polynomials of
> the same positive weight \(W\).  If
> \[
>  A(\partial)^m(B^m)=0\qquad(m\geq1),
>  \tag{3.2}
> \]
> then their Newton segments are disjoint.
>
> Consequently, after possibly exchanging \(x\) and \(y\), there is an
> integer \(\epsilon>0\) such that every monomial of \(A\) has
> \(x\)-exponent at least \(\epsilon\) larger than every monomial of \(B\).
> For every fixed \(Q\),
> \[
>  A(\partial)^m(QB^m)=0
>  \qquad
>  \left(m>\frac{\deg_x Q}{\epsilon}\right).
>  \tag{3.3}
> \]

Because \(A\) and \(B\) have the same positive weight, the expression in
(3.2) has weight zero and is a scalar.  Their supports lie on the same
affine lattice line
\[
 ua+vb=W.
 \tag{3.4}
\]
Since \(u\ne v\), ordinary degree \(a+b\) is strictly monotone along this
line.

Assume that the Newton segments meet.  Let \(\rho\) be the unique point of
their intersection having least ordinary degree.  Clear its lattice
denominator.  Applying the one-variable constant-term theorem to the two
independently shifted Laurent polynomials gives an integer \(k>0\) for
which
\[
 c=[X^{k\rho}]A^k\ne0,\qquad
 e=[x^{k\rho}]B^k\ne0.
 \tag{3.5}
\]
The use of two independent Laurent variables makes the constant term of
the product equal to the product of the two desired coefficients, so the
same \(k\) works on both sides.

Put \(A_0=A^k\), \(B_0=B^k\), \(\alpha=k\rho\), and
\(s=\alpha_x+\alpha_y\).  Specialize the finitely generated coefficient
domain to a number field while retaining \(ce\ne0\), and choose a
sufficiently large unramified rational prime \(p\) at which \(c\) and \(e\)
are units.  Expand
\[
 A_0(\partial)^p(B_0^p)
 =\sum_\gamma
   \gamma!\,[X^\gamma]A_0^p\,[x^\gamma]B_0^p.
 \tag{3.6}
\]
Every contributing \(\gamma\) lies in
\[
 p\operatorname{Newt}(A_0)\cap p\operatorname{Newt}(B_0),
 \tag{3.7}
\]
so \(\gamma_x+\gamma_y\geq ps\).

Choose \(p\) larger than every coordinate occurring in either Newton
segment.  If \(\gamma\) is not divisible by \(p\) coordinatewise,
Frobenius makes each of the two coefficients in (3.6) divisible by \(p\),
while
\[
 v_p(\gamma!)=
 \left\lfloor\frac{\gamma_x}{p}\right\rfloor+
 \left\lfloor\frac{\gamma_y}{p}\right\rfloor
 \geq s-1.
 \tag{3.8}
\]
Such a summand therefore has valuation at least \(s+1\).

If \(\gamma=p\beta\), its factorial valuation is
\(\beta_x+\beta_y\).  The unique point of the intersection with ordinary
degree \(s\) is \(\alpha\).  At every other \(\beta\), either the factorial
valuation is at least \(s+1\), or one of the two Frobenius coefficients is
divisible by \(p\).  The term \(\gamma=p\alpha\), on the other hand, has
valuation exactly \(s\) and unit coefficient congruent to \(c^pe^p\).
It is the unique term of minimum valuation.  Hence (3.6) is nonzero,
contradicting (3.2) at \(m=kp\).

The Newton segments are therefore disjoint.  Along (3.4), disjoint
segments are strictly ordered in the \(x\)-coordinate and oppositely
ordered in the \(y\)-coordinate.  One coordinate supplies the linear
derivative deficit in (3.3).

When \(u=v\), the face is ordinarily homogeneous instead.  The
split-symbol theorem already proves GVC for the homogeneous binary
operator \(A(\partial)\); disjoint Newton segments need not follow and are
not claimed.

The same prime argument applies on a moving output ray.

> **Theorem 3.2 (shifted-ray endpoint theorem).**  Keep \(u\ne v\), let
> \(\delta=(\delta_x,\delta_y)\in\mathbb N^2\), and suppose \(A\) has
> weight \(W\), while \(B\) has weight
> \(W+\operatorname{wt}(\delta)\).  If
> \[
>  [x^{m\delta_x}y^{m\delta_y}]
>  A(\partial)^m(B^m)=0\qquad(m\geq1),
>  \tag{3.9}
> \]
> then
> \[
>  \operatorname{Newt}(A)\cap
>  \bigl(\operatorname{Newt}(B)-\delta\bigr)=\varnothing.
>  \tag{3.10}
> \]

If the segments in (3.10) met, choose their unique
least-ordinary-degree endpoint \(\rho\), and expose nonzero coefficients at
\(k\rho\) and \(k(\rho+\delta)\) in one common power as above.  At the
prime dilation, choose \(p\) larger than every coordinate occurring after
the carrier shift.  The factorial in (3.6) is replaced by
\[
 \prod_{i=x,y}
 \frac{(\gamma_i+pk\delta_i)!}{(pk\delta_i)!}.
 \tag{3.11}
\]
For \(p\)-divisible \(\gamma=p\beta\), its \(p\)-adic valuation is
\(\beta_x+\beta_y\); for a non-\(p\)-divisible \(\gamma\), its valuation is
at least \(s-1\), and the two non-Frobenius coefficients again supply two
additional factors of \(p\).  Thus the endpoint is still the unique term
of minimum valuation.

After clearing the denominator of a rational ray, Theorem 3.2 applies to
every fixed rational normalized output direction.  Therefore an
asymptotic defect path cannot be supported on one unequal-weight endpoint
fiber.  Any surviving convolution must repeatedly have an
ordinary-homogeneous lower endpoint or mix several fibers before a single
endpoint is exposed.

## 4. Common-threshold consequence

> **Corollary 4.1 (automatic terminal closure).**  Suppose a positive
> unequal weight and a threshold \(W\) satisfy
> \[
>  \operatorname{wt}(\text{every monomial of }\lambda)\geq W,\qquad
>  \operatorname{wt}(\text{every monomial of }P)\leq W.
>  \tag{4.1}
> \]
> Then (1.1) implies the GVC conclusion for every fixed \(Q\).

Let \(A\) and \(B\) be the equality faces in (4.1).  The highest-weight
part of (1.1) is \(A(\partial)^m(B^m)=0\), so Theorem 3.1 gives a linear
coordinate deficit on the equality face.  Every strict operator or
polynomial selection consumes a positive integral amount of weight
defect.  A nonzero term of \(\Lambda^m(QP^m)\) can contain only
\(O(\operatorname{wt}Q)\) strict selections.  They cannot repair the
linear equality-face deficit, and eventual vanishing follows.

Thus no exact moment radical is needed after an unequal common threshold
has been found.

There is also a density version that treats the ordinarily homogeneous
alternative.

> **Corollary 4.2 (single-face asymptotic termination).**  Fix one
> exposed face pair.  Suppose every selection pattern contributing to a
> proposed mixed coefficient uses only \(o(m)\) factors off that pair.
> If the pair is unequal-weight, Theorems 3.1--3.2 kill the coefficient.
> If the operator face is ordinarily homogeneous, the split-symbol
> separator kills it.  Hence a nonzero asymptotic defect must use at least
> two distinct faces with positive limiting selection density.

For an unequal-weight face, the coordinate gap is \(cm\) for some
\(c>0\).  For an ordinary homogeneous face, translated complete
polarization supplies an integral Laurent weight with the same linear
gap.  The finitely many original operator and polynomial monomials give a
uniform bound \(C\) on the amount by which one off-face selection can
change either gap.  Thus \(o(m)\) off-face selections and the fixed
multiplier contribute only \(o(m)+O_Q(1)\), which cannot repair the linear
deficit.

There is a global Newton criterion for the remaining multifactorial
convolution.  Define
\[
 \mathcal N_\Lambda=\operatorname{Newt}(\lambda),\qquad
 \mathcal D_P=
 \operatorname{conv}\{\beta\in\mathbb N^2:\partial^\beta P\ne0\}.
 \tag{4.2}
\]
The second polygon is the convex hull of the downward monomial support of
\(P\).

> **Theorem 4.3 (least-intersection termination criterion).**  If
> \(\mathcal N_\Lambda\cap\mathcal D_P\) is empty or has a componentwise
> least point, then (1.1) implies the GVC conclusion.  Consequently, a
> binary counterexample requires this intersection to have at least two
> incomparable Pareto-minimal points.

Introduce a generic translation \(z=(z_1,z_2)\), angular Laurent
variables \(T=(T_1,T_2)\), and radial variables \(U=(U_1,U_2)\).  Put
\[
 F_z(T,U)=
 \lambda(U_1T_1^{-1},U_2T_2^{-1})
 P(z_1+T_1,z_2+T_2).
 \tag{4.3}
\]
For the multifactorial constant-term functional
\[
 \Gamma(U^\alpha T^\nu)=
 \begin{cases}
  \alpha!,&\nu=0,\\
  0,&\nu\ne0,
 \end{cases}
 \tag{4.4}
\]
one has the exact identity
\[
 \Gamma(F_z^m)=\Lambda^m(P^m)(z).
 \tag{4.5}
\]

The support of \(F_z\) is Cartesian over the coefficient field \(k(z)\).
A term indexed by an operator exponent \(\alpha\) and a translated
polynomial exponent \(\beta\) has angular weight \(\beta-\alpha\) and
radial degree \(\alpha\).  Hence its balanced radial polytope is exactly
\[
 \mathcal B(F_z)=\mathcal N_\Lambda\cap\mathcal D_P.
 \tag{4.6}
\]
Indeed, an angularly balanced distribution gives equal average operator
and polynomial exponents; conversely, independent convex
representations of any point in the intersection give such a
distribution.

If the intersection is nonempty and has a least point, the multiradial
prime-isolation argument applies to (4.5): pure vanishing forces zero
outside the convex hull of the angular weights of \(F_z\).  If the
intersection is empty, this exclusion is immediate.  Strict angular
separation then has a linear gap on \(F_z^m\), while the translated support
of a fixed \(Q\) is bounded.  Thus
\(\Gamma(Q(z+T)F_z^m)=0\) eventually, which is precisely the mixed GVC
identity at generic \(z\), and therefore identically.

## 5. The ordinary-homogeneous Pareto tie

The unequal-weight argument leaves one geometrically distinguished
possibility.  If an incomparable Pareto edge is not seen by unequal
positive weights, then it is parallel to
\[
 (1,-1);
\]
all of its radial exponent vectors have the same ordinary degree.

For a single multiplicative radial channel this case is rigid in every
degree.

> **Theorem 5.1 (homogeneous binary factorial rigidity).**  Let
> \(C\in\mathbb C[U_1,U_2]\) be homogeneous.  If
> \[
>  \mathcal L(C^m)=0\qquad(m\geq1),\qquad
>  \mathcal L(U_1^aU_2^b)=a!b!,
>  \tag{5.1}
> \]
> then \(C=0\).

This is the homogeneous binary case of the Factorial Conjecture proved by
Liu and Sun.  There is also a short interval-moment formulation.  If
\(e=\deg C\), then the beta integral gives
\[
 \mathcal L(C^m)
  =(em+1)!\int_0^1 C(t,1-t)^m\,dt.             \tag{5.2}
\]
The one-variable polynomial-moment theorem applied to
\(C(t,1-t)\) makes (5.1) imply \(C(t,1-t)=0\), and homogeneity then gives
\(C=0\).

Consequently an ordinary-homogeneous Pareto contribution which really has
the form \(C^m\) cannot support a defect: its radial support disappears.
Together with Theorem 3.1, this closes every one-channel Pareto edge,
including the slope \(-1\) edge not covered by unequal weights.

The qualification “one-channel” is essential.  Let \(S\) mark a Hall or
jet selection and put
\[
 G(S,U)=S\,U_1+S^{-1}U_2.
 \tag{5.3}
\]
Then
\[
 \operatorname{CT}_S G=0,\qquad
 \operatorname{CT}_S G^2=2U_1U_2.
 \tag{5.4}
\]
Thus
\[
 \operatorname{CT}_S(G^m)
 \ne\bigl(\operatorname{CT}_S G\bigr)^m
 \tag{5.5}
\]
even in the smallest two-channel example.  A moving rational coefficient
can be made multiplicative before extraction,
\[
 [S^{qm}]H(S)^m=\operatorname{CT}_S(S^{-q}H(S))^m,
 \tag{5.6}
\]
after clearing denominators and restricting to a subsequence.  But (5.6)
produces a factorially weighted **beta--torus** moment, not the factorial
moment of one fixed binary polynomial.  When the radial degree is \(e\),
the exact identity is
\[
 \mathcal L_U\operatorname{CT}_S G(S,U)^m
 =(em+1)!\int_0^1
   \operatorname{CT}_S G(S,t,1-t)^m\,dt.
 \tag{5.7}
\]
Neither the interval-moment theorem nor the
Duistermaat--van der Kallen theorem permits the two operations on the
right of (5.7) to be interchanged.

There is a second compatibility failure.  A toric blow-up can make two
incomparable lattice vectors comparable by the Euclidean algorithm, but
it does not preserve \(\mathcal L\).  Under
\[
 (a,b)\longmapsto(a+b,b)
 \tag{5.8}
\]
the factorial weight is multiplied by
\[
 \frac{(a+b)!b!}{a!b!}=\frac{(a+b)!}{a!},
 \tag{5.9}
\]
which depends on the exponent vector.  Hence ordinary resolution of the
Pareto edge is not a valid Hall reduction of the contraction unless the
resulting binomial factors are also controlled.

At this stage of the argument, the exact target can be stated without
filtration language:

> **Beta--torus Hall problem.**  For the rank-one Cartesian polynomials
> \(G\) produced by a binary operator symbol and a translated polynomial,
> prove that all moments in (5.7) can vanish only if the angular support is
> strictly separated or the active Hall/jet support drops.

Theorem 5.1 proves this when there is no genuine \(S\)-convolution.
Prime-endpoint Bessel rigidity proves the adjacent three-level special
case.  The general positive-density multi-level case is not a consequence
of either theorem separately.  Theorem 7.4 quaterdecies supplies the
finite-trace separation after a scale-compatible trace has been exposed;
the promotion from the prime-dependent Hall shell to that trace remains
open.

The conclusion cannot be replaced by the assertion that the beta--torus
kernel itself is Mathieu--Zhao.  Long's three-term circuit is the sharp
warning.  The following calculation also shows why it does not obstruct
the Hall conclusion.

> **Theorem 5.2 (minimal Bernstein circuit is Hall-terminal).**  Put
> \[
>  \Phi(H)=\int_0^1\operatorname{CT}_Z H(U,Z)\,dU
> \]
> and
> \[
>  F=(c+dZ^{-1})\bigl(a(1-U)+bUZ\bigr).
>  \tag{5.10}
> \]
> Then, for every \(m\geq1\),
> \[
>  \Phi(F^m)=\frac{(ac+bd)^m}{m+1}.
>  \tag{5.11}
> \]
> Hence all pure moments vanish exactly when \(ac+bd=0\).  In that case
> the associated binary differential pair is a linear Hall annihilator
> and becomes one-sided after a linear coordinate change.

Indeed, choose \(k\) copies of \(dZ^{-1}\) and \(k\) copies of \(bUZ\).
Angular balance forces the two selection counts to agree.  The beta
identity
\[
 \int_0^1 U^k(1-U)^{m-k}\,dU
 =\frac1{(m+1)\binom mk}
 \tag{5.12}
\]
cancels one of the two binomial coefficients, leaving
\[
 \frac1{m+1}
 \sum_{k=0}^m\binom mk(ac)^{m-k}(bd)^k,
\]
which is (5.11).

To recover the differential pair, take
\[
 A(X,Y)=bX+aY,\qquad P(x,y)=dx+cy.
 \tag{5.13}
\]
With \(Z=T_2/T_1\) and \(U_1=U,\ U_2=1-U\),
\[
 A(U_1/T_1,U_2/T_2)P(T_1,T_2)=F.
 \tag{5.14}
\]
The equation \(ac+bd=0\) is precisely \(A(\partial)P=0\).
Thus suitable linear coordinates give \(A(\partial)=\partial_u\) and
\(P=v\), after which every polynomial multiplier has an immediate
derivative-count cutoff.

For \(a=b=c=1,d=-1\), (5.10) is Long's polynomial
\[
 (1-Z^{-1})\bigl((1-U)+UZ\bigr).
 \tag{5.15}
\]
It has zero pure moments but
\[
 \Phi(Z^{-1}F^m)=\frac{(-1)^{m-1}}{m+1}.
 \tag{5.16}
\]
This disproves the unrestricted beta--torus Mathieu assertion.  It is not
a GVC counterexample: \(Z^{-1}=T_1/T_2\) is not the translate of a
polynomial multiplier, and (5.13) is already Hall-terminal.  Therefore a
valid induction must retain the polynomial-multiplier cone and allow a
split-symbol/Hall separator as a terminal outcome; radial support loss in
the original torus coordinates is too strong.

Theorem 5.2 closes the smallest two-by-two beta circuit exactly.  Any
unresolved Hall convolution must contain a genuine higher-jet selection
level after all homogeneous split faces and linear annihilator circuits
have been removed.  This is the first case in which the Bernstein
finite-difference cancellation is not merely a disguised Hall
annihilator.

The first such circuit can also be closed without a degree bound.  It is
the affine-rank-two parallelogram obtained by compressing the balance
relations of the Dvorsky--Long circuit to two exponent coordinates.

> **Theorem 5.3 (primitive cusp parallelogram is impossible).**  Let
> \(r,s\geq1\), and let \(u,v\in\mathbb Z^2\) be linearly independent
> with
> \[
>  u+v=(-r,s).
> \tag{5.17}
> \]
> Assume that the four exponents
> \[
>  (r,0),\quad (r,0)+u,\quad (r,0)+v,\quad (0,s)
> \tag{5.18}
> \]
> are nonnegative.  Give all four corresponding polynomial monomials
> nonzero coefficients, and give both monomials of
> \[
>  A(X,Y)=aX^r+bY^s,\qquad ab\ne0,
> \tag{5.19}
> \]
> nonzero coefficients.  If \(r\ne s\), then the three scalar equations
> \[
>  [A(\partial)^m(P^m)]_{(0,0)}=0,\qquad m=1,2,3,
> \tag{5.20}
> \]
> are inconsistent.  If \(r=s\), the operator is homogeneous and the
> split-symbol theorem is terminal.  Thus no nondegenerate primitive
> cusp parallelogram supports the remaining beta--torus Hall
> convolution.

To prove the theorem, assume \(r\ne s\) and use the affine lattice
coordinates \(U,V\) along
\(u,v\).  Scaling the two lattice characters and the overall polynomial
coefficient puts the four polynomial coefficients in the form
\[
 p(U,V)=1+U+V+tUV,\qquad t\ne0.
\tag{5.21}
\]
Absorb all remaining coefficient ratios into \(q\), the coefficient ratio
of \(Y^s\) to \(X^r\).  In the \(m\)-th power, choosing \(k\) copies of
\(Y^s\) forces the polynomial selection to have affine exponent
\((k,k)\).  Hence
\[
 M_m=\sum_{k=0}^m
  \binom mk q^k(r(m-k))!(sk)!
  [U^kV^k]p(U,V)^m.                            \tag{5.22}
\]
Put
\[
 C_n=\binom{2n}{n},\qquad
 T_n=\frac{(3n)!}{(n!)^3}.
\tag{5.23}
\]
The first two equations in (5.20) give, successively,
\[
 q=-\frac{r!}{s!\,t},\qquad
 t=\frac4{C_r+C_s-4}.                           \tag{5.24}
\]
The denominator in (5.24) is positive.  The diagonal coefficients needed
at the third moment are
\[
 [U^kV^k]p^3=
 1,\ 3(t+2),\ 3t(t+2),\ t^3
 \quad(k=0,1,2,3).
\tag{5.25}
\]
After (5.24), division by \((r!)^3\) reduces the third moment to
\[
 E_{r,s}
 =T_r-T_s+\frac92(C_s-C_r)(C_r+C_s-2).          \tag{5.26}
\]

It remains only to check that \(E_{r,s}\ne0\) for \(r\ne s\).  By
antisymmetry it is enough to take \(r<s\).  Directly,
\[
 E_{1,2}=24,\qquad E_{1,3}=-54,\qquad E_{2,3}=-78.
\tag{5.27}
\]
For \(s\geq4\), set \(R_n=T_n/C_n^2\).  The ratios
\[
 \frac{T_n}{T_{n-1}}
 =\frac{3(3n-1)(3n-2)}{n^2}\geq6
\tag{5.28}
\]
and
\[
 \frac{R_{n+1}}{R_n}
 =\frac{3(3n+1)(3n+2)}{4(2n+1)^2}>1             \tag{5.29}
\]
give
\[
 T_s-T_{s-1}\geq\frac56T_s
 \geq\frac56R_4C_s^2
 =\frac{165}{28}C_s^2
 >\frac92C_s^2.                                  \tag{5.30}
\]
Since
\[
 (C_s-C_r)(C_r+C_s-2)<C_s^2,
\tag{5.31}
\]
equations (5.30)--(5.31) imply \(E_{r,s}<0\).  This proves the
claim.

Theorem 5.3 is the required nonzero Bernstein pivot for the saturated
primitive parallelogram.  It is stronger than a bounded degree-seven
regression: it excludes the entire two-endpoint/four-channel cusp family
in all orders and all degrees using only the first three scalar moments.
It does **not** by itself prove that every sparse multi-level Hall
relation admits such a parallelogram inside its original support.
Introducing a missing lattice point during toric saturation changes a
zero support coefficient into a new variable, so quadratic generation of
the saturated toric ideal is not yet a support-decreasing induction.

The same calculation closes every four-channel configuration, including
the sparse dilations whose saturation is missing the parallelogram point.

> **Corollary 5.4 (all primitive four-channel circuits).**  Keep the
> two primitive operator endpoints
> \[
>  A=(r,0),\qquad B=(0,s),
> \]
> and suppose the active polynomial support consists of the four
> distinct points \(A,B,C,D\), all with nonzero coefficients.  Then this
> support is terminal.

At moment two an additional polynomial pair can affect the scalar
contraction only when its exponent sum lies in
\[
 \{2A,A+B,2B\}.                                    \tag{5.32}
\]
Apart from the three endpoint pairs, such an equality is one of three
types: a collinear midpoint, a centered adjacent three-level circuit, or a
noncollinear parallelogram.  The first is a one-channel face, the second is
closed by prime-endpoint Bessel rigidity, and the third is Theorem 5.3.
In the centered case the unused fourth point is strictly one-sided in the
normal angular weight; unless a second return occurs, it therefore has a
separator and cannot participate with positive density.
If two distinct additional quadratic returns occur, the unique affine
dependences force all four points to be collinear, returning to the first
case.  If no additional return occurs, the off-endpoint channels do not
enter the first two scalar moments.  The first moment fixes the product of
the two endpoint coefficient ratios, and the second moment, after this
substitution, is
\[
 (2r)!-4(r!)^2+\frac{(r!)^2}{(s!)^2}(2s)!
 =(r!)^2(C_r+C_s-4).                               \tag{5.34}
\]
This vanishes only when \(r=s=1\).  That exceptional operator is
homogeneous and split-symbol terminal.  Thus every nonhomogeneous
primitive four-channel circuit either reduces to a previously closed
quadratic return or dies before its higher toric relation can enter.

For example, in affine support coordinates the sparse configuration
\[
 \{(0,0),(1,0),(0,1),(h,h)\},\qquad h\geq2,
\tag{5.35}
\]
has no supported parallelogram and first relation
\[
 h(1,0)+h(0,1)=(h,h)+(2h-1)(0,0),
\tag{5.36}
\]
but (5.34), rather than quadratic generation after saturation, already
closes it.

Finite truncations do not inherit this four-channel conclusion.  The
smallest exact warning already occurs for five polynomial channels.

> **Example 5.5 (finite-prefix minor inheritance fails).**  Put
> \[
>  A(X,Y)=X^2-\frac13Y^3
> \]
> and
> \[
>  P(x,y)=x^2+y^3+c+dx+exy^3.
> \tag{5.37}
> \]
> If \(M_m=[A(\partial)^m(P^m)]_{(0,0)}\), then
> \[
> \begin{aligned}
>  M_1&=0,\\
>  M_2&=-8(2de-11),\\
>  M_3&=96(15ce^2+21de-122),\\
>  M_4&=-5760(200ce^2-12d^2e^2+188de-967).
> \end{aligned}
> \tag{5.38}
> \]
> Thus
> \[
>  (c,d,e)=\left(\frac{13}{30},\frac{11}{2},1\right)
> \tag{5.39}
> \]
> gives
> \[
>  M_1=M_2=M_3=0,\qquad M_4=1\,205\,760.
> \tag{5.40}
> \]
> Nevertheless no four-channel restriction of (5.37) has its first
> three moments zero.  Omitting either endpoint makes \(M_1=\pm2\);
> omitting \(x\) or \(xy^3\) makes \(M_2=88\); and omitting the constant
> channel makes \(M_3=-624\).

This is not a pure-moment-zero pair and hence is neither a counterexample
to the beta--torus Hall problem nor a GVC counterexample.  It does rule
out any proof which tries to infer minor inheritance from a fixed initial
block of moments: the first genuine five-channel cancellation can hide
every four-channel obstruction through moment three, and moment four is
a higher-minor Bernstein pivot.

The mixed fourth-pivot support visible in Example 5.5 extends to an
all-degree family.

> **Theorem 5.6 (unit-line half-bridge pivot).**  Let \(n\geq2\) be even,
> and put
> \[
>  A(X,Y)=X-\frac1{n!}Y^n,\qquad
>  P(x,y)=x+y^n+c+d y^{n/2}+e x y^{n/2},
> \tag{5.41}
> \]
> with \(cde\ne0\).  Then the four scalar equations
> \[
>  [A(\partial)^m(P^m)]_{(0,0)}=0,\qquad 1\leq m\leq4,
> \tag{5.42}
> \]
> are inconsistent.

Set
\[
 C=C_n=\binom{2n}{n},\qquad
 T=T_n=\frac{(3n)!}{(n!)^3},\qquad
 Q=Q_n=\frac{(4n)!}{(n!)^4},
\]
and use the two coefficient invariants
\[
 u=de,\qquad v=ce^2.
\]
The first moment vanishes by the normalization in (5.41), while direct
selection gives
\[
\begin{aligned}
 M_2={}&C-2-4u,\\
 M_3={}&-6\left(
  3v-3(C-2)u+\frac{T-9C+12}{6}\right),\\
 M_4={}&144\left(
  (C-2)v+\frac C2u^2-\frac{T-6C+6}{3}u
  +\frac{Q-16T+72C-72}{144}\right).
\end{aligned}
\tag{5.43}
\]
Thus \(M_2=M_3=0\) determines \(u,v\), and substitution into the last
line gives
\[
 M_4=\frac12H_n,
\qquad
 H_n=
 2Q-40CT+48T+81C^3-180C^2+132C-48.
\tag{5.44}
\]
Now \(H_2=-480\).  For \(n\geq4\), put
\[
 R_n=\frac{Q_n}{C_nT_n}
 =\frac{\binom{4n}{n}}{\binom{2n}{n}}.
\]
Here \(R_4=26\), and
\[
 \frac{R_{n+1}}{R_n}
 =\frac{4(4n+1)(4n+3)}{3(3n+1)(3n+2)}>1.
\tag{5.45}
\]
Consequently \(2Q-40CT>0\).  The remaining terms are positive because
\[
 81C^3-180C^2+132C-48
 =3(3C-4)(9C^2-8C+4)>0.
\]
Hence \(H_n>0\) for \(n\geq4\), proving the theorem.

An exact finite-moment search complements the example.  For each endpoint
pair below it enumerates zero through three additional support points in
\[
 [0,4r]\mathbin{\times}[0,4s]
\]
and saturates each successive moment ideal by the product of all active
channel coefficients.

| \((r,s)\) | supports | first unit at \(M_2\) | at \(M_3\) | at \(M_4\) |
|---:|---:|---:|---:|---:|
| \((1,2)\) | \(13\,288\) | \(13\,082\) | \(173\) | \(33\) |
| \((1,3)\) | \(41\,728\) | \(41\,299\) | \(394\) | \(35\) |
| \((1,4)\) | \(95\,368\) | \(94\,636\) | \(679\) | \(53\) |
| \((2,3)\) | \(253\,576\) | \(252\,442\) | \(1\,074\) | \(60\) |

There is no torus survivor.  A point outside the displayed box cannot
occur in a balanced return through order four, so the calculation covers
arbitrary nonnegative polynomial support with at most five channels for
these four endpoint pairs: some \(M_m\), \(m\leq4\), is nonzero.  Only 181
of the \(403\,960\) support ideals require the genuine fourth-moment
pivot.  This is an exact rational finite-moment computation, not a proof
for arbitrary endpoint orders or for three operator endpoints.  The
remaining proof target is correspondingly sharper: classify those
fourth-pivot return graphs uniformly in \(r,s\), or find the first endpoint
pair at which one survives.

For that purpose, record a balanced selection by
\[
 (k,n_A,n_B,n_C,n_D,n_E),
\]
where \(k\) is the number of \(Y^s\)-operator selections and the remaining
entries are the polynomial-channel multiplicities.  Canonicalizing the
complete order-two-through-four row sets under endpoint exchange and the
six permutations of \(C,D,E\) collapses the 181 fourth-pivot supports to
exactly **14 return-matrix types**.  All 14 already occur for \((1,2)\);
the other three cusps introduce no new type.  Their realization counts
are
\[
 32,23,17,17,17,14,14,12,8,8,6,6,4,3.
\]
Thus the observed fourth-pivot problem is a finite template problem, not
403,960 unrelated coefficient eliminations.  This classification remains
a computation at this point; Theorem 5.9 below explains why the same
14-type list is exhaustive in the reduced early-entry regime.

Late channel entry cannot repair an earlier core obstruction.  The
following elementary observation turns that statement into a finite
reduction when combined with the prime-endpoint witnesses.

> **Lemma 5.7 (late-entry transfer).**  Let \(A\) be a
> constant-coefficient operator symbol, let \(P_0\) be a polynomial, and
> let \(\gamma\in\mathbb N^2\).  Define
> \[
> \tau_A(P_0;\gamma)=
> \min\left\{m\geq1:
> \begin{array}{l}
> \text{there are }1\leq j\leq m,\quad
> \beta\in\operatorname{supp}(P_0^{m-j}),\\
> \alpha\in\operatorname{supp}(A^m)
> \text{ with }\beta+j\gamma=\alpha
> \end{array}\right\},
> \tag{5.46}
> \]
> with minimum \(\infty\) when the set is empty.  Then, for every scalar
> \(c\) and every \(m<\tau_A(P_0;\gamma)\),
> \[
> [A(\partial)^m(P_0+c z^\gamma)^m]_{z=0}
> =
> [A(\partial)^m(P_0^m)]_{z=0}.
> \tag{5.47}
> \]
> In particular, a nonzero core moment at an order below
> \(\tau_A(P_0;\gamma)\) excludes the enlarged support.

Indeed, every term in the difference of the two sides uses \(j\geq1\)
copies of \(z^\gamma\).  It contributes at the origin only when its total
polynomial exponent equals an exponent of \(A^m\), which is exactly the
condition excluded in (5.46).  The same statement holds for several new
channels when \(\tau\) is defined as the first return using any of them.

> **Corollary 5.8 (prime witnesses bound entry order).**  Suppose a fixed
> core has the following property after finite-type specialization: there
> is \(p_0\) such that, for every prime \(p\geq p_0\), at least one of its
> moments of order \(p\) and \(2p\) is nonzero.  Then a pure-zero
> extension by a new channel must satisfy
> \[
>  \tau_A(P_0;\gamma)\leq4(p_0+1).
> \tag{5.48}
> \]

Otherwise choose, by Bertrand's postulate, a prime
\[
 p\in\left(\left\lfloor\frac{\tau}{4}\right\rfloor,
                 2\left\lfloor\frac{\tau}{4}\right\rfloor\right).
\]
Then \(p\geq p_0\) and \(2p<\tau\), so Lemma 5.7 transfers both candidate
prime moments from the core to the extension, a contradiction.  Thus the
one-channel and three-level prime-endpoint cases reduce every sufficiently
late fifth channel to a finite list of early-entry return matrices.  The
remaining uniform issue is to control \(p_0\) on each core type, not to
search unbounded filtration depth.

> **Theorem 5.9 (fourteen reduced early-entry return types).**  Keep two
> operator endpoints
> \[
>  A=(r,0),\qquad B=(0,s),
> \]
> and three distinct off-endpoint polynomial channels \(C,D,E\).  Call a
> balanced return primitive if its selection row is not a sum of
> lower-order balanced rows.  Assume:
>
> 1. there is a primitive non-endpoint return at order two;
> 2. every one of \(C,D,E\) occurs in a primitive return by order three;
> 3. single-extra-channel returns and coincident channels have already
>    been removed by the one-channel, centered three-level, and
>    support-loss reductions.
>
> Then, up to endpoint exchange and permutation of \(C,D,E\), the
> primitive return data has one of exactly 14 forms: two
> double-quadratic types and twelve quadratic--cubic types.

Write
\[
 T_m(k)=(r(m-k),sk).
\]
A reduced primitive quadratic return must use two distinct extra
channels, so normalize it to
\[
 D+E=T_2(q).
\tag{5.51}
\]
If \(C\) also enters at order two, its return can be written
\[
 C+E=T_2(q').
\tag{5.52}
\]
The two target levels are distinct, since otherwise \(C=D\).  The three
two-element subsets of \(\{0,1,2\}\) are reduced by endpoint exchange to
\[
 \{q,q'\}=\{0,1\}\quad\text{or}\quad\{0,2\},
\tag{5.53}
\]
giving the two double-quadratic types.

Otherwise \(C\) first enters through a primitive cubic return.  That
return contains no copy of \(A\) or \(B\), since subtracting such a copy
would leave a lower-order balanced return.  It cannot contain one copy
each of \(C,D,E\): subtracting (5.51) would make \(C\) equal to \(A\) or
\(B\).  The single-channel case \(3C=T_3(\ell)\) was removed by
assumption.  Hence, after exchanging \(D,E\), the cubic is
\[
 aC+(3-a)E=T_3(\ell),
\qquad a\in\{1,2\}.
\tag{5.54}
\]
Endpoint exchange puts \(q\in\{0,1\}\).  For \(q=0\), all four levels
\(\ell=0,1,2,3\) remain, producing \(2\cdot4=8\) types.  For \(q=1\),
endpoint exchange fixes the quadratic level and sends
\(\ell\mapsto3-\ell\), leaving \(\ell=0,1\), hence \(2\cdot2=4\) types.
Together with (5.53), this gives \(2+8+4=14\).  The exact census realizes
every type already at \((r,s)=(1,2)\).

> **Lemma 5.10 (no new primitive return at order four).**  Under the
> hypotheses of Theorem 5.9, every balanced order-four row is a sum of
> the primitive order-two and order-three rows and endpoint rows.  Thus
> the special primitive order-four branch is empty.

Project orthogonally to the endpoint line and write the nonzero transverse
coordinate of \(D\) as \(t\).  In a quadratic--cubic type, (5.51) and
(5.54) give
\[
 \bar D=t,\qquad \bar E=-t,\qquad
 \bar C=\frac{3-a}{a}t.
\tag{5.56}
\]
If \(a=1\), a transverse-zero selection with multiplicities
\((n_C,n_D,n_E)\) satisfies
\[
 n_E=n_D+2n_C,\qquad
 n_C+n_D+n_E=3n_C+2n_D\leq4.
\]
Its nonzero possibilities are exactly the quadratic row, twice that row,
and the cubic row.  If \(a=2\), the equation is
\[
 n_C=2(n_E-n_D),\qquad
 n_C+n_D+n_E=3n_E-n_D\leq4,
\]
and gives the same list.  Adding endpoint copies supplies all remaining
order-four rows.  In a double-quadratic type,
\(\bar C=\bar D=t,\ \bar E=-t\), so
\[
 n_E=n_C+n_D,\qquad
 n_C+n_D+n_E=2(n_C+n_D)\leq4.
\]
The solutions are the two quadratic rows and their pairwise sums.  Finally,
\(t=0\) is the ordinary-homogeneous case already closed by Theorem 5.1.

The coefficient elimination for these types can be written once and for
all.  Define, with out-of-range subscripts interpreted as zero,
\[
 W_{m,k}=
 \binom mk\left(-\frac{r!}{s!}\right)^k
 (r(m-k))!(sk)!,
 \qquad
 E_m=\sum_{k=0}^m\binom mk W_{m,k},
\tag{5.55}
\]
and put
\[
 S_q=W_{3,q}+W_{3,q+1},\qquad
 A_q=W_{4,q}+2W_{4,q+1}+W_{4,q+2}.
\tag{5.57}
\]

> **Corollary 5.11 (eight factorial obstructions).**  The twelve
> quadratic--cubic types reduce in coefficient-paired
> copies to the six scalar obstructions
> \[
>  \mathcal H_{q,\ell},
>  \qquad
>  (q,\ell)\in
>  \{(0,0),(0,1),(0,2),(0,3),(1,0),(1,1)\},
> \tag{5.58}
> \]
> and the two double-quadratic types reduce to
> \(\mathcal D_{0,1}\) and \(\mathcal D_{0,2}\).

For a quadratic--cubic type, set
\[
 u=-\frac{E_2}{2W_{2,q}},\qquad
 v=-\frac{E_3+6S_qu}{3W_{3,\ell}}.
\tag{5.59}
\]
The fourth moment after the first two equations is exactly
\[
 \mathcal H_{q,\ell}
 =
 E_4+12A_qu+6W_{4,2q}u^2
 +12(W_{4,\ell}+W_{4,\ell+1})v.
\tag{5.60}
\]
The two choices \(a=1,2\) in (5.54) give the same expression because the
relevant multinomial coefficients are both \(3\) at order three and
\(12\) at order four.

For a double-quadratic type \(\{q,q'\}\), let \(u,v\) be the solution of
\[
 \begin{pmatrix}
  2W_{2,q}&2W_{2,q'}\\
  6S_q&6S_{q'}
 \end{pmatrix}
 \binom uv
 =
 -\binom{E_2}{E_3}.
\tag{5.61}
\]
When the determinant is nonzero, the remaining fourth moment is
\[
 \begin{aligned}
 \mathcal D_{q,q'}={}&E_4
 +12A_qu+6W_{4,2q}u^2
 +12A_{q'}v+6W_{4,2q'}v^2\\
 &+12W_{4,q+q'}uv.
 \end{aligned}
\tag{5.62}
\]
Thus the reduced early-entry problem has become eight explicit
factorial-ratio inequalities, plus the degenerate determinant branches
of (5.61).  Exact evaluation finds no zero of any of the eight
obstructions for \(1\leq r,s\leq100,\ r\ne s\); this is a regression, not
an all-\((r,s)\) proof for the final wedge not closed below.

Several parts of that arithmetic residue admit uniform certificates.  Put
\[
 C_n=\binom{2n}{n},\qquad
 L_n=\binom{3n}{n},\qquad
 M_n=\binom{4n}{n}.
\tag{5.63}
\]

> **Lemma 5.12 (determinants and the corner obstruction).**  For
> \(r\ne s\), both determinants in (5.61) are nonzero.  Moreover
> \[
>  \mathcal H_{0,0}(r,s)>0.
> \tag{5.64}
> \]

After division by the positive common factorial powers, the determinants
are
\[
\begin{aligned}
 \Delta_{0,1}
 &=12C_r(2L_r-3C_r+3C_s-6)>0,\\
 \Delta_{0,2}
 &=-12C_rC_s(L_r+L_s-6)<0.
\end{aligned}
\tag{5.65}
\]
Indeed \(L_n\geq3\), with equality only at \(n=1\), and
\(2L_n\geq3C_n\), again with equality only at \(n=1\).  These facts follow
directly by taking successive ratios.

For (5.64), also note
\[
 C_n\geq2,\qquad
 L_n\geq\frac32C_n,\qquad
 M_n\geq\frac43L_n.
\tag{5.66}
\]
The sequences \(C_n\), \(L_n/C_n\), and \(M_n/L_n\) are increasing.
Substitute the exact
successor ratios
\[
\frac{C_{n+1}}{C_n}=\frac{2(2n+1)}{n+1},\quad
\frac{T_{n+1}}{T_n}
=\frac{(3n+1)(3n+2)(3n+3)}{(n+1)^3},\quad
\frac{Q_{n+1}}{Q_n}
=\frac{\prod_{j=1}^4(4n+j)}{(n+1)^4}
\tag{5.67}
\]
in
\(\mathcal H_{0,0}(r,s+1)-\mathcal H_{0,0}(r,s)\).  On writing
\[
\begin{array}{lll}
C_r=2+x_r,&L_r=\frac32C_r+y_r,&M_r=\frac43L_r+z_r,\\
C_s=2+x_s,&L_s=\frac32C_s+y_s,&M_s=\frac43L_s+z_s,
\end{array}
\]
and \(s=1+n_0\), its numerator expands into 932 monomials with
nonnegative rational coefficients and positive constant coefficient.
Thus \(\mathcal H_{0,0}\) is strictly increasing in \(s\).  At \(s=1\),
the stronger \(r\geq2\) bounds
\[
 C_r\geq6,\qquad L_r\geq\frac52C_r,\qquad
 M_r\geq\frac{28}{15}L_r
\]
give a 19-term coefficient-positive base expansion.  Finally
\(\mathcal H_{0,0}(1,1)=0\), so the first positive difference handles
\(r=1,s\geq2\).  The exact expansions are checked by
`scripts/research_binary_gvc_eight_obstructions.py --prove-h00`.

> **Lemma 5.13 (the two negative corner obstructions).**  For every
> \(r\ne s\),
> \[
>  \mathcal H_{0,3}(r,s)<0,\qquad
>  \mathcal D_{0,2}(r,s)<0.
> \tag{5.68}
> \]

Write \(T_n=C_nL_n\) and \(Q_n=C_nL_nM_n\), as above.  In addition to
(5.66), the exact one-step ratios satisfy
\[
\begin{aligned}
 \frac{C_{n+1}}{C_n}-3
 &=\frac{n-1}{n+1},\\
 \frac{L_{n+1}/L_n}{C_{n+1}/C_n}-\frac53
 &=\frac{(n-1)(n+2)}{12(2n+1)^2},\\
 \frac{M_{n+1}/M_n}{L_{n+1}/L_n}-\frac75
 &=\frac{(n-1)(n+2)(17n^2+17n+6)}
 {45(3n+1)^2(3n+2)^2}.
\end{aligned}
\tag{5.69}
\]
Consequently, if \(s>r\), there are nonnegative real numbers
\(x,y,z,a,b,d\) such that
\[
\begin{array}{lll}
 C_r=2+x,&L_r=\frac32C_r+y,&M_r=\frac43L_r+z,\\
 C_s=C_r(3+a),&
 L_s=L_r\left(\frac53(3+a)+b\right),&
 M_s=M_r\left(\frac75\left(\frac53(3+a)+b\right)+d\right).
\end{array}
\tag{5.70}
\]
Indeed, each ratio across the gap \(s-r\) is a product of the one-step
ratios in (5.69).  The reverse ordering has the same parametrization with
\(r,s\) exchanged.

After the substitutions \(T_i=C_iL_i,\ Q_i=C_iL_iM_i\), the positive
denominators of \(-\mathcal H_{0,3}\) and
\(-\mathcal D_{0,2}\) are respectively
\[
 2C_r,\qquad 6C_rC_s(L_r+L_s-6)^2.
\tag{5.71}
\]
The numerator of \(-\mathcal H_{0,3}\), expanded in the six slacks of
(5.70), has 266 terms when \(s>r\) and 361 terms when \(r>s\).
The corresponding numerator of \(-\mathcal D_{0,2}\) has 2236 terms in
either orientation.  Every coefficient in all four expansions is a
nonnegative rational number and every constant coefficient is positive.
This proves (5.68); the \(\mathcal D_{0,2}\) expression is also invariant
under endpoint exchange.  The exact expansions and the ratio identities
are checked by
`scripts/research_binary_gvc_eight_obstructions.py
--prove-negative-corners`.

> **Lemma 5.14 (ordered-tail closure).**  For every \(r\ne s\),
> \[
>  \mathcal H_{0,1}(r,s)\ne0,\qquad
>  \mathcal H_{0,2}(r,s)\ne0,\qquad
>  \mathcal H_{1,0}(r,s)\ne0,\qquad
>  \mathcal H_{1,1}(r,s)\ne0,\qquad
>  \mathcal D_{0,1}(r,s)\ne0.
> \tag{5.72}
> \]

For the proof, put
\[
\begin{aligned}
 c_n&=\frac{C_{n+1}}{C_n},&
 \rho_n&=\frac{L_{n+1}/L_n}{C_{n+1}/C_n},&
 \sigma_n&=\frac{M_{n+1}/M_n}{L_{n+1}/L_n}.
\end{aligned}
\]
The three sequences are strictly increasing, since
\[
\begin{aligned}
 c_{n+1}-c_n
 &=\frac2{(n+1)(n+2)},\\
 \rho_{n+1}-\rho_n
 &=\frac{3(n+1)}{2(2n+1)^2(2n+3)^2},\\
 \sigma_{n+1}-\sigma_n
 &=\frac{64(n+1)(18n^4+72n^3+103n^2+62n+15)}
 {9(3n+1)^2(3n+2)^2(3n+4)^2(3n+5)^2}.
\end{aligned}
\tag{5.73}
\]
Together with the monotonicity in (5.66), this gives a sharper cone for
any \(u>v\geq N\):
\[
\begin{array}{lll}
 C_v=C_N+x,&
 L_v=\dfrac{L_N}{C_N}C_v+y,&
 M_v=\dfrac{M_N}{L_N}L_v+z,\\
 C_u=C_v(c_N+a),&
 L_u=L_v(c_N+a)(\rho_N+b),&
 M_u=M_v(c_N+a)(\rho_N+b)(\sigma_N+d),
\end{array}
\tag{5.74}
\]
with all six slacks nonnegative.  Products across a gap of more than one
step only increase the three ratio lower bounds.

After substituting (5.74), coefficient-positive numerator expansions
give the following ordered tails:
\[
\begin{array}{c|c|c|c}
\text{expression}&\text{region}&\text{sign}&\text{terms}\\ \hline
\mathcal H_{0,1}&s>r\ge2&+&456\\
\mathcal H_{0,1}&r>s\ge4&-&570\\
\mathcal H_{0,2}&s>r\ge6&+&456\\
\mathcal H_{0,2}&r>s&-&361\\
\mathcal H_{1,1}&s>r\ge2&+&340\\
\mathcal H_{1,1}&r>s\ge6&+&340\\
\mathcal H_{1,0}&s>r\ge2&+&340\\
\mathcal D_{0,1}&s>r\ge2&+&1500.
\end{array}
\]
On a fixed endpoint ray, the first line of (5.74) alone gives
coefficient-positive expansions of 14 or 19 terms for the remaining
quadratic--cubic strips.  The corresponding \(\mathcal D_{0,1}\) rays
have 24 or 43 terms.  These fixed rays, followed by exact evaluation of
their finite complements, prove the first, second, and fourth claims in
(5.72), and reduce the other two to
\[
 \mathcal H_{1,0}(r,s),\quad \mathcal D_{0,1}(r,s),
 \qquad r>s\ge4.
\tag{5.75}
\]

The broad linear cone misses one genuine factorial correlation on
(5.75).  Define \(K_n=L_nM_n/C_n^3\).  Its exact successor ratio is
\[
 \frac{K_{n+1}}{K_n}
 =\frac{(n+1)(4n+1)(4n+3)}{2(2n+1)^3}
 =1+\frac{8n^2+7n+1}{2(2n+1)^3}>1.
\tag{5.76}
\]
Hence \(K_s\geq K_4\) on the wedge.  Replace the lower-endpoint \(M\)
slack in (5.74) by the multiplicative slack
\[
 L_sM_s=K_4C_s^3+w,\qquad w\geq0,
\tag{5.77}
\]
and retain the three ordered gap-ratio slacks.  The numerator of
\(-\mathcal H_{1,0}\) then has 408 terms and that of
\(-\mathcal D_{0,1}\) has 1692 terms.  Every coefficient is nonnegative
and each constant coefficient is positive.  Their cleared denominators
are positive on the same cone, proving the remaining two claims in
(5.72).

The ratio identities, every coefficient and denominator sign, the
positive constant coefficients, and the finite complements are checked by
`scripts/research_binary_gvc_eight_obstructions.py --prove-three-more`.

Primary references for Theorem 5.1 and its moment input are:

- D. Liu and X. Sun,
  [*The Factorial Conjecture and Images of Locally Nilpotent
  Derivations*](https://doi.org/10.1017/S0004972719000546),
  Bull. Aust. Math. Soc. 101 (2020), 71--79;
- J. P. Françoise, F. Pakovich, Y. Yomdin, and W. Zhao,
  *Moment vanishing problem and positivity: some examples*,
  Bull. Sci. Math. 135 (2011), 10--32.
- C. D. Long,
  [*Counterexamples to the \(xz\)-Conjecture and the Mathieu Conjecture
  for \(SU(2)\)*](https://arxiv.org/abs/2607.19012).

## 6. The pre-separation termination frontier

The general fixed-depth theorem kills every bounded filtration defect.
Theorem 2.1 and (2.5) reduce every leading component to one root
multiplicity and two normalized jet wings.  Corollary 4.1 kills every
defect path that reaches one fixed common threshold.

Before applying finite-trace separation, the only surviving architecture
is:

1. the defect depth is proportional to \(m\);
2. no single positive weight is a common threshold for all terms used with
   positive limiting frequency;
3. no fixed rational output ray has one unequal-weight exposed endpoint;
4. \(\mathcal N_\Lambda\cap\mathcal D_P\) has no componentwise least
   point; and
5. coefficient extraction at that moving depth uses at least two distinct
   faces with positive limiting selection density, attached to
   incomparable Pareto endpoints.

The last point is essential.  At depth \(\lfloor\rho m\rfloor\), the
relevant coefficient is a convolution coefficient of a power.  It need
not be the power of one fixed rank-one face, so Theorem 3.1 cannot simply
be applied to it.  This is the same multiplicativity loss that appears in
attempts to turn a finite-channel SIC tensor into a GVC pair.

A complete binary termination theorem must solve the **restricted**
beta--torus Hall problem above: the positive-density coupled convolution
between incomparable Pareto endpoints must either acquire a least radial
vector after a factorial-compatible Hall/jet reduction, have smaller
Hall/jet support after cancellation, or acquire a split-symbol separator
whose gap dominates every polynomial multiplier.  The last alternative
is required by Theorem 5.2 and cannot be replaced by unrestricted
beta--torus Mathieu--Zhao behavior.

The degree-at-most-six calculations verify this three-way alternative in
their bounded support ranges.  Theorems 5.1--5.6 remove the one-channel,
minimal Bernstein, primitive four-channel, and unit-line half-bridge
branches.  Theorems 5.9--5.14 now classify the reduced two-endpoint
early-entry regime uniformly and prove all eight resulting factorial
obstructions nonzero.  Thus no fourth-pivot arithmetic inequality remains
in that regime.

Equations (5.4), (5.9), and (5.16) prevent an automatic
arbitrary-support straightening.  A complete proof must use genuinely
all-moment information to expose a proper primitive-return packet; simple
deletion of channels is invalid.  Example 5.5 shows that no fixed
three-moment truncation can do this: the first three moments of the whole
five-channel convolution can vanish while every four-channel minor
already fails.  Section 7 gives the exact return-cone formulation,
proves that its extreme circuits have uniformly bounded support, and
isolates the all-moment exposure statement.  The repeated-digit argument
in Theorem 7.4 quaterdecies proves the component-separation part
simultaneously for two and for at least three active operator endpoints.
It does not by itself prove that the Hall--jet initial shell is
scale-compatible.

## 7. The factorial circuit-exposure theorem

The last descent can be stated directly on the two independent selection
multisets, without introducing artificial pairings between operator and
polynomial channels.  Write
\[
 \lambda=\sum_{i\in I}a_iX^{\alpha_i},\qquad
 P=\sum_{j\in J}b_jx^{\beta_j},
\]
where \(\alpha_i,\beta_j\in\mathbb N^2\).  Define the balanced return
semigroup
\[
 \mathsf R=
 \left\{(p,q)\in\mathbb N^I\times\mathbb N^J:
 |p|=|q|,
 \sum_i p_i\alpha_i=\sum_jq_j\beta_j
 \right\}.
\tag{7.1}
\]
Its order and radial vector are
\[
 d(p,q)=|p|=|q|,\qquad
 \rho(p,q)=\sum_i p_i\alpha_i.
\]

> **Proposition 7.1 (exact return-semigroup expansion).**  The scalar
> pure moment is
> \[
> \left[\Lambda^m(P^m)\right]_{x=0}
> =
> \sum_{\substack{(p,q)\in\mathsf R\\d(p,q)=m}}
> \binom mp\binom mq\,a^pb^q\,\rho(p,q)!.
> \tag{7.2}
> \]
> A channel absent from every nonzero element of \(\mathsf R\) never
> contributes to any scalar moment and may be deleted.

Formula (7.2) is the two multinomial expansions followed by equality of
the total operator and polynomial exponents.  The deletion statement is
therefore genuine support loss, not the invalid operation of discarding a
channel which participates in later returns.

Let \(\mathsf K\) be the real cone obtained from (7.1).  Its defining
matrix has the signed columns
\[
 (1,\alpha_i)\quad(i\in I),\qquad
 (-1,-\beta_j)\quad(j\in J)
\tag{7.3}
\]
in a vector space of dimension three.

> **Proposition 7.2 (uniform circuit bound).**  Every extreme ray of
> \(\mathsf K\) has support on at most four variables among
> \(p_i,q_j\).  Equivalently, every extremal balanced return is a
> positive circuit involving at most four operator and polynomial
> channels in total.

An extreme nonnegative kernel ray has inclusion-minimal support.  Its
supported columns in (7.3) form a minimally dependent set in
\(\mathbb R^3\), and such a circuit has at most four elements.  This
bound is independent of all operator orders and polynomial degrees.

Normalize a nonzero return by \(d(p,q)=1\).  The radial image of the
resulting compact return polytope is exactly
\(\mathcal N_\Lambda\cap\mathcal D_P\), as in (4.6).  Every vertex of
that image lifts to a vertex of the return polytope.  Hence:

> **Corollary 7.3 (two-circuit Pareto skeleton).**  If the balanced
> radial polytope has no componentwise least point, then it has two
> incomparable Pareto-minimal vertices lifted by two circuit returns,
> each supported on at most four channel variables.
>
> If only two operator channels are active, the union of the two circuits
> uses at most six polynomial channels.  Equality can occur only when
> each circuit uses exactly one operator channel and three polynomial
> channels; the two polynomial triples are disjoint and are centered at
> the two different operator endpoints.

Indeed, each circuit contains at least one operator and one polynomial
variable, so it uses at most three polynomial variables.  If two circuits
use six distinct polynomial variables, both attain that maximum and
therefore each has only one operator variable.  Incomparability forces
those operator variables to be the two different endpoints.  Thus the
only packet beyond the closed at-most-five-channel regime is an opposite
three-by-three packet, not an arbitrary large support.

Prime dilation gives one further exact restriction.  For a primitive
integer circuit \(h=(p,q)\), put \(d=d(h)\) and \(\rho=\rho(h)\).  Its
pure \(N\)-fold repeat inside the moment of order \(Nd\) is
\[
 \mathcal T_h(N)=
 \binom{Nd}{Np}\binom{Nd}{Nq}
 a^{Np}b^{Nq}(N\rho)!.
\tag{7.4}
\]

> **Lemma 7.4 (prime profile of a repeated circuit).**  After
> specialization to a number field, let \(\ell\) be a sufficiently large
> unramified prime at which the active coefficients and the two base
> multinomial coefficients are units.  Then
> \[
>  v_\ell\!\left(\mathcal T_h(\ell)\right)
>  =\rho_x+\rho_y.
> \tag{7.5}
> \]

Lucas congruences make the two dilated multinomial coefficients congruent
to their undilated values modulo \(\ell\).  The coefficient powers are
units, while Legendre's formula gives
\(v_\ell((\ell\rho_i)!)=\rho_i\) when \(\ell>\rho_i\).
Thus unequal ordinary radial degrees are visible at a single good-prime
dilation.  The genuinely coupled remainder consists of equal-profile
circuits and mixed sums of them, exactly the ordinary-homogeneous
beta--torus tie.

There is a stronger probe inside the equal-total layer which is useful,
but not by itself an exposure theorem.

> **Lemma 7.4 bis (radial digit spectrum).**  Fix
> \(a+b=c+d=S\).  If
> \[
>  v_p\!\left((p^ka)!(p^kb)!\right)
>  =
>  v_p\!\left((p^kc)!(p^kd)!\right)
> \]
> for every prime \(p\) and one (equivalently every) \(k\geq0\), then
> \(\{a,b\}=\{c,d\}\).

Indeed Legendre's digit-sum formula gives
\[
 v_p((p^ka)!(p^kb)!)
 =\frac{p^kS-s_p(a)-s_p(b)}{p-1}.
\]
Equality for every prime is equivalent to
\(v_p(a!b!)=v_p(c!d!)\) for every prime, hence to
\(a!b!=c!d!\).  On \(0\leq a<S/2\), the successive ratio
\[
 \frac{(a+1)!(S-a-1)!}{a!(S-a)!}
 =\frac{a+1}{S-a}
\]
is strictly less than one, so the factorial product determines the
unordered pair.  Thus higher prime powers separate every radial profile
except coordinate reversal.  This still does not isolate a partial
convolution: at a small distinguishing prime, the circuit multinomials
and coefficient valuations can contribute at the same order.  The
remaining arithmetic issue is therefore packet exposure in the presence
of those extra valuations, not failure of the radial factorials to
distinguish profiles.

The same-prime calculation also sees the nonhomogeneous operator jet
exactly.  Fix an integer \(s\) no larger than the ordinary degree of any
active channel, and put
\[
 h_i=|\alpha_i|-s,\qquad k_j=|\beta_j|-s.
\]
For a balanced return, the common jet excess is
\[
 J(u,v)=\sum_i u_ih_i=\sum_jv_jk_j.
\]

> **Lemma 7.4 ter (binary jet--carry score).**  Let \(D\) be the largest
> ordinary degree of an active channel.  Fix a base order \(d\), and let
> \(\ell>dD\) be a prime at which the active coefficients are units.
> For a balanced return \((u,v)\) of order \(\ell d\), write
> \[
>  e_u=\frac1\ell\sum_i(u_i\bmod\ell),\qquad
>  e_v=\frac1\ell\sum_j(v_j\bmod\ell),
> \]
> let \(t_x,t_y\) be the least residues of \(\rho_x,\rho_y\), and define
> \[
>  c_\rho=\frac{t_x+t_y-(J\bmod\ell)}{\ell}\in\{0,1\}.
> \]
> Then
> \[
> v_\ell\!\left(
>  \binom{\ell d}{u}\binom{\ell d}{v}\rho(u,v)!
> \right)
> =sd+\left\lfloor\frac J\ell\right\rfloor
>   +e_u+e_v-c_\rho.                     \tag{7.5a}
> \]
> Moreover \(\kappa(u,v)=e_u+e_v-c_\rho\geq0\), with equality if and
> only if every entry of both \(u\) and \(v\) is divisible by \(\ell\).
> Thus a Frobenius dilation of a base return of jet excess \(j\) has
> score \(sd+j\), while every non-Frobenius return has a strict carry
> penalty \(\kappa\geq1\).

Because \(\ell>d\), every selection entry is less than \(\ell^2\).
Legendre--Kummer therefore gives
\[
 v_\ell\binom{\ell d}{u}=e_u,\qquad
 v_\ell\binom{\ell d}{v}=e_v.
\]
The radial coordinates are less than \(\ell^2\), and their sum is
\(s\ell d+J\).  Hence
\[
 v_\ell(\rho_x!\rho_y!)
 =sd+\left\lfloor\frac J\ell\right\rfloor-c_\rho.
\]
This proves (7.5a).  If \(c_\rho=0\), equality forces
\(e_u=e_v=0\).  If \(c_\rho=1\), equality would force exactly one of
\(e_u,e_v\) to vanish.  That selection vector would be divisible by
\(\ell\), making its radial vector—and hence the common radial
vector—divisible by \(\ell\), a contradiction.

The corresponding unit is also explicit.  Write
\[
 u_i=\ell\bar u_i+u_i^{(0)},\qquad 0\leq u_i^{(0)}<\ell,
\]
and similarly for \(v\).  Removing the multiples of \(\ell\) block by
block and using Wilson's theorem gives
\[
 \ell^{-e_u}\binom{\ell d}{u}
 \equiv
 (-1)^{e_u}
 \frac{d!}{\prod_i\bar u_i!\prod_i u_i^{(0)}!}
 \pmod\ell,                                      \tag{7.5b}
\]
with the analogous formula for \(v\).  If
\(\rho_i=\ell H_i+t_i\), then
\[
 \ell^{-(sd+\lfloor J/\ell\rfloor-c_\rho)}\rho_x!\rho_y!
 \equiv
 (-1)^{H_x+H_y}H_x!H_y!t_x!t_y!
 \pmod\ell.                                      \tag{7.5c}
\]

The zero-low-digit restriction is unnecessary for the valuation and unit
calculation.  Its affine extension is the useful input for the next
promotion test.

> **Lemma 7.4 ter bis (two-digit affine jet--carry transform).**  Keep the
> hypotheses of Lemma 7.4 ter, strengthen the good-prime bound to
> \(\ell>D(d+1)\), let \(0\leq r<\ell\), and consider order
> \(\ell d+r\).  Write
> \[
>  u_i=\ell\bar u_i+u_i^{(0)},\qquad
>  \sum_i u_i^{(0)}=r+\ell e_u,
> \]
> and similarly for \(v\).  If
> \(\rho_i=\ell H_i+t_i\), put
> \[
>  c_{\rho,r}=
>  \frac{t_x+t_y-((sr+J)\bmod\ell)}{\ell}\in\{0,1\}.
> \]
> Then
> \[
> \begin{aligned}
>  v_\ell\binom{\ell d+r}{u}&=e_u,\\
>  v_\ell\!\left(\binom{\ell d+r}{u}
>                  \binom{\ell d+r}{v}\rho!\right)
>  &=sd+\left\lfloor\frac{sr+J}{\ell}\right\rfloor
>    +e_u+e_v-c_{\rho,r}.                       \tag{7.5c'}
> \end{aligned}
> \]
> Moreover
> \[
>  \ell^{-e_u}\binom{\ell d+r}{u}
>  \equiv
>  (-1)^{e_u}
>  \frac{d!r!}
>       {\prod_i\bar u_i!\prod_i u_i^{(0)}!}
>  \pmod\ell.                                   \tag{7.5c''}
> \]

The digit sum of the numerator is \(d+r\), while the sum of the digit
sums in the denominator is

\[
 \sum_i\bar u_i+\sum_i u_i^{(0)}
 =(d-e_u)+(r+\ell e_u)
 =d+r+(\ell-1)e_u.
\]

Legendre--Kummer proves the first line.  Since

\[
 H_x+H_y
 =sd+\left\lfloor\frac{sr+J}{\ell}\right\rfloor-c_{\rho,r},
\]

the radial factorial gives the second.  Removing multiples of \(\ell\)
from numerator and denominator and applying Wilson proves (7.5c'').  Unlike
the \(r=0\) penalty in Lemma 7.4 ter, the relative affine score
\(e_u+e_v-c_{\rho,r}\) can equal \(-1\): two no-multinomial-carry low
digits may themselves cross one radial digit boundary.  This disappears
when \(\ell\) is chosen larger than the maximum radial size of the bounded
low digit.  The formula, rather than nonnegativity at one prescribed prime,
is the invariant needed for joint low-digit tomography.

> **Corollary 7.4 ter ter (one-high-digit quotients are free).**  At orders
> \(\ell+r\), every selection side in the two-digit transform has high-digit
> count either zero or one.  Hence its high-digit quotient is empty or is one
> marked channel; it cannot contain two inequivalent semigroup
> decompositions.

Indeed \(\sum_i\bar u_i=1-e_u\), and nonnegativity forces
\(e_u\in\{0,1\}\), with the same argument on the other side.  Thus the
nonfree quotient semigroup which obstructs circuit peeling at orders
\(\ell d\) with \(d>1\) disappears completely in the one-high-digit probes.
The remaining compatibility problem is narrower: low-digit bridges carrying
different singleton marks can still cancel in the same translated Hall
shell.  Proving that the translation-Hasse rows triangularize those
singleton marks, or producing a nonflat kernel, is the next exposure lemma.

The later
[affine singleton-localization theorem](BINARY_GVC_AFFINE_SINGLETON_LOCALIZATION.md)
removes the state-combinatorial part of this ambiguity.  In the homogenized
radial/source configuration, every state sharing a singleton fibre differs
from it by the fixed integer kernel and conformally decomposes into its fixed
Graver basis.  Exposed singleton columns admit only finitely many
prime-independent signed low corrections, with explicit carry and unit
formulas; every unbounded branch repeats one fixed primitive move with
positive density.  Its radial-carry Hasse formula further compresses every
long binary carry interval with residue (R) to the single row
(-F^{(R+1)}(-1)).  This does not by itself triangularize the **linear**
specialized packet sum.  The later translation-tangent theorem shows that a
coefficient-blind module-inheritance implication is false; the valid missing
statement must use the complete Cartesian derivative tower to expose a flat
twist or a conformal curvature block with common high quotient.

For a one-carry multinomial, \(d!=d(d-1)!\); hence (7.5b)--(7.5c)
factor the normalized weight into a high-digit quotient weight, a
low-digit bridge weight, and an explicit scalar.  The high-digit
quotient need not itself be balanced: its balance defect is exactly
the radial vector carried by the low-digit atom.

There is a canonical algebraic package for these terms.  For a symbol
\(A=\sum_i a_i z^{\gamma_i}\) over an integral model, define
\[
 \Phi_\ell(A)=\sum_i a_i^\ell z^{\ell\gamma_i},\qquad
 G_\ell(A)=\frac{A^\ell-\Phi_\ell(A)}{\ell}.       \tag{7.5d}
\]
The quotient is integral, and its reduction modulo \(\ell\) is exactly
the generating polynomial of the one-carry selections, with the unit
weights in (7.5b).  The binomial expansion gives the exact first-ghost
identity
\[
 A^{\ell d}\equiv
 \Phi_\ell(A)^d+
 \ell d\,\Phi_\ell(A)^{d-1}G_\ell(A)
 \pmod{\ell^2}.                                   \tag{7.5e}
\]
Thus the one-sided score correction is a contraction with one ghost.
When both selections carry once and \(c_\rho=1\), the product of the
operator and polynomial ghosts loses one power of \(\ell\) through the
radial factorial and enters the same score shell.  Its residual radial
kernel is explicitly \(t!(\ell-t)!\equiv(-1)^t t\pmod\ell\).

The ghost is an additive \(p\)-derivation cocycle:
\[
 G_\ell(A+B)=G_\ell(A)+G_\ell(B)
 +\frac{(A+B)^\ell-A^\ell-B^\ell}{\ell}.          \tag{7.5f}
\]
At this stage the shell problem is finite and algebraic:
triangularize the one-ghost terms and the bilinear ghost pairing under
successive support splits, including the next-order correction of the
pure Frobenius terms.  A triangular identity with nonzero diagonal would
give factorial-compatible support descent.  Theorem 7.4 quaterdecies
supersedes the subsequent componentwise triangularization once (SC) is
known, but not this shell-promotion step.

The two primitive diagonal blocks can be written in closed form.

> **Proposition 7.4 quater (centered and beta ghost diagonals).**
> Let \(\ell>3\) be prime.
>
> 1. For a centered triple \(A(z)=a+bz+cz^2\), the central Cartier
>    coefficient of \(G_\ell(A)\) is
>    \[
>     ac\,b^{\ell-2}H_\ell(ac/b^2),                \tag{7.5g}
>    \]
>    where
>    \[
>     H_\ell(X)=
>     \sum_{t=1}^{(\ell-1)/2}
>     \frac{(\ell-1)!}{t!^2(\ell-2t)!}X^{t-1}
>     \in\mathbb F_\ell[X].
>    \tag{7.5h}
>    \]
>    The terminal geometric relation is a universal factor:
>    \(H_\ell(1)=0\).
> 2. For a two-by-two beta atom, put \(X\) for the ratio of the
>    level-one operator--polynomial coefficient product to the
>    level-zero product.  Its two-ghost radial diagonal is
>    \[
>     K_\ell(X)=
>     \sum_{n=1}^{\ell-1}\frac{(-1)^n}{n}X^n.
>    \tag{7.5i}
>    \]
>    It has the universal support and Hall factors \(X(X+1)\).

For (7.5g), choosing \(t\) copies of each endpoint and
\(\ell-2t\) copies of the middle level gives (7.5h) by (7.5b).
At \(X=1\), its sum is
\[
 \frac{[z^\ell](1+z+z^2)^\ell-1}{\ell}.
\]
The numerator is divisible by \(\ell^2\).  Indeed
\[
 (1+z+z^2)^\ell=(1-z^3)^\ell(1-z)^{-\ell}.
\]
The term with no \(z^3\)-selection contributes
\(\binom{2\ell-1}{\ell-1}\equiv1\pmod{\ell^2}\);
every term with a positive number of \(z^3\)-selections contains one
factor \(\ell\) from \(\binom{\ell}{j}\) and another from the remaining
binomial coefficient by Lucas.  This proves \(H_\ell(1)=0\).

For (7.5i),
\(\ell^{-1}\binom{\ell}{n}\equiv(-1)^{n-1}/n\), while the binary radial
unit is \(n!(\ell-n)!\equiv(-1)^n n\).  Combining the two ghost
multinomials and the radial unit gives the displayed coefficient.
The factor \(X\) is support loss, and
\[
 K_\ell(-1)=\sum_{n=1}^{\ell-1}\frac1n=0
\]
by pairing \(n\) with \(\ell-n\).  This is exactly the linear Hall
annihilator from Theorem 5.2.

The universal factors do not exhaust a block at one fixed prime.  For
example,
\[
\begin{aligned}
 H_{11}(X)&=-(X-4)(X-1)(X^2+3),\\
 K_7(X)&=-X(X-2)^2(X+1)(X+3)^2
\end{aligned}
\qquad\text{over the indicated prime fields}.      \tag{7.5j}
\]
Thus the strongest fixed-prime triangularity claim is false.  Even
cross-prime avoidance fails for algebraic beta ratios.  Before reduction,
the beta diagonal is
\[
 D_\ell(X)=\frac{(1+X)^\ell-1-X^\ell}{\ell}.
\tag{7.5k}
\]
For every prime \(\ell\geq5\),
\[
 X(X+1)(X^2+X+1)\mid D_\ell(X).                  \tag{7.5l}
\]
The first two factors are immediate.  If \(X^2+X+1=0\), then
\(1+X=-X^2\); separating \(\ell\equiv1,2\pmod3\) proves (7.5l).
Conversely the gcd of \(D_5\) and \(D_7\) is exactly the polynomial in
(7.5l), so these are all beta roots persistent for every
\(\ell\geq5\).

The primitive cube-root branch is not Hall-terminal, but it is not a
pure-moment survivor: the original beta atom already has first moment
proportional to \(1+X\), which is nonzero there.  Hence a successful
triangular argument must use a **multi-row atom block** containing the
ordinary atom pivot together with its ghost diagonal.  The centered
block analogously requires its Bessel endpoint row when the Cartier
diagonal vanishes.  The exact factorizations, the persistent beta gcd,
and a bounded rational cross-prime search are replayed by
`scripts/research_binary_gvc_ghost_shell.py`.

> **Corollary 7.4 quinquies (augmented ghost atom blocks are terminal).**
> Adjoin to each primitive ghost diagonal its ordinary atom endpoint
> rows.  The beta block then has only the Hall zero \(X=-1\), while the
> centered block forces support loss.

For the beta atom, Theorem 5.2 supplies the ordinary row \(1+X\).
Consequently every prime-dependent root and the persistent primitive
cube-root branch are removed; the only common zero is \(X=-1\), already
the linear Hall annihilator.

For the centered atom, use the Bessel endpoint variables
\[
 U=C^\ell,\qquad V=D^\ell.
\]
The prime-endpoint rows are
\[
 U,\qquad U^2+2V.
\]
Their Jacobian with respect to \(U,V\), restricted after the first row,
is triangular with diagonal \(1,2\).  Since \(\ell\) is odd, the block
forces \(U=V=0\), hence \(C=D=0\), which is support loss.  This closes the
finite atom arithmetic completely.  The compatibility problem left here
is to put these rows over the same high-digit quotient.  Finite-trace
separation avoids the later componentwise triangularization after (SC),
but does not put the rows over that quotient.

Formula (7.5a) is the exact migration filtration.  Its level \(L\)
consists of the returns satisfying
\[
 \left\lfloor J/\ell\right\rfloor+\kappa=L.
\]
It therefore mixes Frobenius dilations of base returns of jet excess
\(L\) with carry bridges having smaller high-digit jet excess.

On an ordinary-homogeneous radial layer, choose \(s=r\); then \(J=0\).
The score-\(rd\) initial sum is therefore exactly the Frobenius image,
up to the common sign \((-1)^{rd}\), of the base moment of order \(d\).
Under the pure premise it vanishes again and supplies no proper packet.
The first potentially new information is the score-\(rd+1\) shell.
Equation (7.5a) gives only two possibilities:

1. \(c_\rho=0\) and \(\{e_u,e_v\}=\{0,1\}\);
2. \(c_\rho=1\) and \(e_u=e_v=1\).

These are respectively the one-sided and two-sided **one-carry
bridges**.  Neither automatically loses support: the smallest
full-support examples are the two-level selections
\((1,1)\leftrightarrow(1,1)\) at \((r,d,\ell)=(1,1,2)\), and
\((3)\leftrightarrow(1,1,1)\) on levels
\((1)\leftrightarrow(0,1,2)\) at \((2,1,3)\).  They are the already
closed beta and centered-three-level circuits, but larger full-support
bridges exist.  Thus the next valid target is a primitive decomposition
of the one-carry shell, not a claim that every first correction already
loses support.  The exact bounded regression is
`scripts/research_binary_gvc_frobenius_carry.py`.

> **Corollary 7.4 sexies (primitive one-carry atoms).**  On the
> homogeneous line of levels \(0,\ldots,r<\ell\), every
> inclusion-minimal one-sided bridge uses at most three low-digit
> channels.  Every inclusion-minimal two-sided bridge uses exactly two
> low-digit channels on each side.

For three distinct levels \(a,b,c\), the kernel vector
\[
 (b-c,\ c-a,\ a-b)\pmod\ell
\]
has three nonzero entries.  Taking their representatives in
\(\{1,\ldots,\ell-1\}\), their sum is either \(\ell\) or \(2\ell\);
negating the vector interchanges those two sums.  Thus every triple
supports a positive one-sided bridge of total \(\ell\), and any
one-sided support of size at least four contains a proper bridge.

For a pair \(a\ne b\), the selections
\[
 n\text{ copies of }a,\qquad \ell-n\text{ copies of }b
\]
realize the nonzero residue \(n(a-b)\), and hence every nonzero residue
as \(n=1,\ldots,\ell-1\).  Consequently any two operator levels and any
two polynomial levels support a common nonzero radial residue.  A
two-sided bridge needs at least two channels on each side, so its
minimal support is exactly \(2+2\).

The primitive low-digit atoms are therefore precisely the support sizes
of the centered three-level and beta/parallelogram circuits already
closed above.  The step not proved by this atom analysis alone is that the
normalized
score-\(rd+1\) sum is compatible with peeling such an atom from its
high-digit quotient.  Theorem 7.4 quaterdecies closes the resulting
finite trace after scale-compatible promotion; proving that promotion
remains Open Problem 7.8.

The quotient itself has a second, independent toric obstruction.

> **Proposition 7.4 septies (first projected-scroll quotient
> obstruction).**  Circuit-only peeling of the high-digit quotient is
> false after inactive levels have been projected away.  The first
> obstruction, up to reversing the level order, is
> \[
>  R_3B_1B_2=R_0B_3^2.                           \tag{7.5m}
> \]
> It is a primitive color-homogeneous partition identity of degree
> three and support five.  On
> \[
>  \operatorname{supp}R=\{0,3\},\qquad
>  \operatorname{supp}B=\{1,2,3\},
> \]
> the two sides of (7.5m) are the only points with one \(R\)-selection,
> two \(B\)-selections, and total level six.  Hence no circuit of support
> at most four joins them.  After restoring level \(R_2\), the
> obstruction decomposes nonconformally as
> \[
>  R_3B_1B_2\longrightarrow R_2B_2^2
>  \longrightarrow R_0B_3^2.                    \tag{7.5n}
> \]

For the projected support, the assertion follows by direct enumeration.
If the \(R\)-selection is \(R_3\), the two \(B\)-levels must sum to
three, giving \(B_1B_2\).  If it is \(R_0\), they must sum to six,
giving \(B_3^2\).  Any proper color-homogeneous subidentity would have
to retain the unique \(R\)-selection on each side.  No one-element
\(B\)-subselection repairs the resulting level gap three, so (7.5m) is
primitive.  Its difference has five nonzero coordinates, whereas a
circuit of the rank-three colored configuration has support at most
four.  The two arrows in (7.5n) preserve both color counts and total
level and have support four.

This does **not** reopen the isolated low-digit blocks of Corollary 7.4
quinquies.  It shows instead that tensoring those blocks with a quotient
and then reducing only by quotient circuits is not a valid global
argument on a projected support.  The replacement remains finite for
every fixed symbol support: the quotient configuration is a coordinate
projection of a rational normal scroll, hence has a finite reduced
Gröbner basis.  Petrović's scroll theorem bounds the degree of every
Graver element of the complete scroll by its degree and bounds every
reduced-Gröbner-basis element after coordinate projection by the degree
of the projected variety.  Thus the refinement is a finite list of
**projected quotient Gröbner blocks**, not another unbounded
filtration-depth search.

There is also an exact threshold for when primitive quotient moves cease
to be visible in some term order.  Bogart--Hemmecke--Petrović prove that
the universal Gröbner basis equals the Graver basis for a rational
normal scroll \(S(a_1,\ldots,a_c)\) unless the scroll dominates
\[
 S(6),\qquad S(5,4),\qquad\text{or}\qquad S(4,3,2).
\]
For the present two-color quotient only the first two alternatives
occur.  Equality is inherited by coordinate subconfigurations.  Hence
every primitive projected quotient move whose completed two-color
scroll has spans below those thresholds is Gröbner-visible; the first
new non-Gröbner primitive architectures require one span at least six,
or two spans at least five and four.  This matches the point at which
the degree-five/six finite closures stop being automatic, and gives two
specific completed-scroll architectures for the next attack.

The exact first obstruction and its two completed paths are replayed by
`scripts/research_binary_gvc_quotient_graver.py`.

> **Corollary 7.4 octies (the first projected quotient block is
> terminal).**  In the reduced fiber of Proposition 7.4 septies,
> (7.5m) either reduces to the already closed circuit blocks or is
> isolated by the radial digit spectrum.

If \(R_2\) is active, use (7.5n).  If \(B_0\) is active, use the second
support-four path
\[
 R_3B_1B_2\longrightarrow R_3B_0B_3
 \longrightarrow R_0B_3^2.                     \tag{7.5o}
\]
Suppose both are absent.  Reversing the levels \(i\mapsto3-i\) sends
(7.5m) to
\[
 R_0B_1B_2=R_3B_0^2,
\]
which is unavailable because \(B_0\) is absent.  Lemma 7.4 bis says
that coordinate reversal is the only other radial profile with the
same prime digit spectrum.  Hence, once the reduced fiber contains no
additional state of the same profile, good-prime dilation exposes the
two-state block itself.  Its ordinary row and ghost row form the
terminal augmented beta block of Corollary 7.4 quinquies.  Thus the
first non-circuit quotient identity does not survive.  The case not
covered by this block alone starts with higher projected Gröbner blocks
for which several same-profile states and their reversals coexist; the
finite-trace theorem handles them simultaneously after scale-compatible
promotion.

There is a useful distinction between **exposing** such a packet and
proving it terminal after exposure.

> **Proposition 7.4 nonies (every exposed scroll packet is
> terminal).**  Fix quotient colors \(1,\ldots,c\), positive color
> counts \(d_1,\ldots,d_c\), and one radial level \(w\).  If the whole
> color-homogeneous quotient fiber with these data is
> factorial-compatibly exposed, then it either has a nonzero ordinary
> endpoint row or its surviving levels have a one-dimensional
> split-symbol separator.  In particular, the dominance-minimal
> non-universal-Gröbner scroll packets \(S(6)\) and \(S(5,4)\) are
> terminal once their whole radial profile is exposed.

For color \(j\), put
\[
 C_j(z)=\sum_i c_{j,i}z^i .
\]
After removing the radial factorial common to the fiber, its order-\(N\)
ordinary endpoint row is
\[
 [z^{Nw}]\prod_{j=1}^c C_j(z)^{Nd_j}
 =
 \operatorname {CT}_z
 \left(
   z^{-w}\prod_{j=1}^c C_j(z)^{d_j}
 \right)^N .                                      \tag{7.5p}
\]
This is exactly the product of the colorwise multinomial expansions;
no pairing of operator and polynomial selections has been introduced.
If (7.5p) vanished for every \(N>0\), the
Duistermaat--van der Kallen constant-term theorem would put \(0\)
outside the Newton interval of
\[
 z^{-w}\prod_j C_j(z)^{d_j}.
\]
In one Laurent variable all its surviving exponents then have one
strict sign.  That is precisely a split-symbol separator; cancellation
of an entire \(C_j\) is support loss.  Otherwise an ordinary endpoint
row is nonzero and the augmented ghost block is terminal by
Corollary 7.4 quinquies.

For completeness, the two minimal exceptional fibers in the
Bogart--Hemmecke--Petrović classification can be seen explicitly.  The
\(S(6)\) witness is the one-color primitive identity
\[
 0+2+6=1+3+4 .
\]
It already lies in the ordinary-homogeneous one-color sector.  On the
support of the genuinely bichromatic \(S(5,4)\) witness,
\[
 R_0R_4B_0B_4=R_1R_5B_1^2,                       \tag{7.5q}
\]
the fiber with color counts \((2,2)\) and total level \(8\) has exactly
four states:
\[
 R_4^2B_0^2,\quad
 R_1R_5B_1^2,\quad
 R_0R_4B_0B_4,\quad
 R_0^2B_4^2 .
\]
The middle positive state is the midpoint, in exponent space, of the
first and last states, while \(R_1R_5B_1^2\) is isolated from the other
three by support-at-most-four circuits.  Thus (7.5q) is the first
genuinely bichromatic failure of circuit or universal-Gröbner
exposure, but Proposition 7.4 nonies proves that its *whole* four-state
profile is terminal.  The repeated-ray factorial partitions
\[
 (2,2),\qquad(2,1,1),\qquad(1,1,1,1)
\]
have respective Stirling bases \(16,4,1\), an independent asymptotic
check that no cancellation can remain confined to distinct pure rays.
The exact fibers, convex certificates, circuit components, and
factorial signatures are replayed by
`scripts/research_binary_gvc_quotient_graver.py`.

Proposition 7.4 nonies means that one does **not** need to prove
terminality separately for every higher projected Graver block.  The
issue left at that point was solely inheritance: the Hall/jet filtration
must
expose one complete color-count/radial-profile packet, rather than an
arbitrary circuit or a proper selection of its states.  Several
profiles can still occur in the same leading quotient convolution, so
the proposition does not by itself prove Conditional Theorem 7.6.
Finite-trace digit separation handles the profiles after (SC), but does
not establish (SC).

The first two prime-power corrections of a complete profile are not an
additional obstruction.

> **Proposition 7.4 decies (exposed packets have an integral Witt
> recursion).**  Keep the colors and counts of Proposition 7.4 nonies.
> Suppose color \(j\) has completed span \(s_j\), and put
> \[
>  T=\sum_jd_js_j,\qquad q=T-w,\qquad
>  F(z)=z^{-w}\prod_jC_j(z)^{d_j}.
> \]
> Let \(p>\max(w,q,d_1,\ldots,d_c)\) be an odd good prime and
> \(S_k=1+p+\cdots+p^{k-1}\).  Write \(R\) for the integral
> coefficient ring of the \(C_j\), and let
> \(\phi:R\longrightarrow R\) fix the integers and send every
> coefficient variable to its \(p\)-th power.  Then
> \[
>  \mathcal G_k=
>  (-1)^{TS_k}
>  \frac{(p^kw)!(p^kq)!}{p^{TS_k}}\,
>  \operatorname {CT}F^{p^k}                    \tag{7.5r}
> \]
> is integral and satisfies
> \[
>  \mathcal G_k\equiv\phi(\mathcal G_{k-1})\pmod {p^k}
>  \qquad(k\ge1).                                \tag{7.5s}
> \]
> Consequently \((\mathcal G_k)_{k\ge0}\) is the ghost sequence of a
> unique \(p\)-typical Witt vector over the \(p\)-torsion-free integral
> coefficient ring.  In particular
> \[
> \begin{aligned}
>  X_1&=\frac{\mathcal G_1-\mathcal G_0^p}{p},\\
>  X_2&=\frac{\mathcal G_2-\mathcal G_0^{p^2}-pX_1^p}{p^2}
> \end{aligned}                                  \tag{7.5t}
> \]
> are integral.  Thus no new denominator or residual correction first
> appears at the \(p^2\)-ghost layer.

There are two independent Gauss congruences behind (7.5s).  First,
\[
 a_N=\operatorname {CT}F^N
\]
is an integral Gauss sequence relative to \(\phi\).  One proof groups
the closed words in the Laurent monomials of \(F\) into cyclic orbits:
primitive necklaces give the integral Witt coordinates, and repetition
gives
\[
 a_{p^k}\equiv\phi(a_{p^{k-1}})\pmod {p^k}.
\tag{7.5u}
\]
This is the elementary closed-word form of the constant-term
congruences studied by Mellit--Vlasenko.

Second, for \(0\leq a<p\), set
\[
 U_k(a)=(-1)^{aS_k}\frac{(p^ka)!}{p^{aS_k}}.
\]
After removing the multiples of \(p\),
\[
 \frac{U_k(a)}{U_{k-1}(a)}
 =
 (-1)^a
 \prod_{\substack{1\leq n\leq ap^k\\p\nmid n}}n .
\]
Modulo \(p^k\), each of the \(a\) complete blocks contains every unit
class modulo \(p^k\); their product is \(-1\) for odd \(p\).
Therefore
\[
 U_k(a)\equiv U_{k-1}(a)\pmod {p^k}.              \tag{7.5v}
\]
The factorial factor in (7.5r) is \(U_k(w)U_k(q)\), and \(\phi\)
fixes this integer.  Multiplying (7.5u) and (7.5v) proves (7.5s).
The Dwork recursive criterion for \(p\)-typical ghost components over
a \(p\)-torsion-free ring with Frobenius lift then proves (7.5t) and
all higher coordinates.  After an integer specialization, \(\phi\)
becomes the identity, which is the simpler consecutive congruence
used in the regression below.

The exact calculation
`scripts/research_binary_gvc_witt_rees.py` checks (7.5t) at \(p=13\)
for two unit coefficient specializations of each of the support-five,
\(S(6)\), \(S(5,4)\), and first larger reversal-symmetric packets.  In
all eight cases the first residual has valuation exactly one and both
\(\mathcal G_2-\mathcal G_1\) and the second Witt residual have
valuation exactly two.

Proposition 7.4 decies also identifies the precise limitation of the
Witt shortcut.  If several profiles \(\pi\) coexist, each gives a Witt
vector \(X_\pi\), but
\[
 \sum_\pi \operatorname {ghost}(X_\pi)=0
 \quad\Longrightarrow\quad
 \sum_\pi X_\pi=0
\]
does **not** imply \(X_\pi=0\) separately.  The ghost functor is
injective, not a source of profile idempotents.  Before repeated-digit
separation, a Witt--Rees closure would have needed one of two stronger
statements:

1. the Hall/jet Rees filtration has a profile-separating initial
   idempotent, so its least Witt coordinate is one complete profile; or
2. cancellation between two least profile Witt vectors forces a
   split-symbol separator or support loss.

This identified a strictness/separation problem, not a missing higher
prime-power correction.  Repeated equal digits provide the required
separation without constructing a Witt idempotent.

There is also no need to separate color counts once an entire oriented
radial vector has been exposed.

> **Proposition 7.4 undecies (complete equal-radial unions are
> terminal).**  Fix \(d\geq1\) and
> \(\rho\in\mathbb N^2\).  At scale \(N\), sum every state in (7.2)
> with
> \[
>  d(p,q)=Nd,\qquad \rho(p,q)=N\rho,
> \]
> without fixing its operator or polynomial color counts, and remove
> the common radial factorial \((N\rho)!\).  The resulting complete
> union is
> \[
> \begin{aligned}
>  E_N(d,\rho)
>  &=[X^{N\rho}]\lambda(X)^{Nd}\,
>    [Y^{N\rho}]P(Y)^{Nd}\\
>  &=\operatorname {CT}_{X,Y}
>    \left(
>      X^{-\rho}Y^{-\rho}\lambda(X)^dP(Y)^d
>    \right)^N .                                  \tag{7.5w}
> \end{aligned}
> \]
> If the radial fiber is nonempty on the active support, the sequence
> \(E_N(d,\rho)\) cannot vanish for every \(N\geq1\).  Consequently a
> complete oriented radial-vector union is terminal after exposure,
> even when its achievable color-count set is not a face or has
> persistent lattice holes.

The first equality in (7.5w) simply performs the operator and
polynomial multinomial sums independently.  The second converts the two
moving coefficient extractions into one fixed Laurent polynomial.  If
the radial fiber is nonempty, then
\[
 \rho\in d\operatorname {Newt}(\lambda)
 \quad\hbox{and}\quad
 \rho\in d\operatorname {Newt}(P).
\tag{7.5x}
\]
Newton polytopes are additive under products, so the origin belongs to
the Newton polytope of the Laurent polynomial in (7.5w).  The
Duistermaat--van der Kallen theorem therefore forbids all of its
positive-power constant terms from vanishing.  If (7.5x) fails after
specialization, the corresponding radial return has disappeared from
the actual support, which is precisely support loss or a separating
face.

Face saturation of the projected color-count set would be too strong.
The smallest warning already has
\[
 C(z)=a(1+z^2)+bz.
\tag{7.5y}
\]
At order \(2N\) and level \(2N\), the achievable number of selections
from \(a(1+z^2)\) is
\[
 0,2,\ldots,2N,
\]
so every odd count is missing at every scale.  Nevertheless the whole
union is \([z^{2N}]C(z)^{2N}\).  Its first two rows are
\[
 b^2+2a^2,\qquad b^4+12a^2b^2+6a^4,
\]
and reduction of the second modulo the first gives \(-14a^4\); hence
even this first nonsaturated family is terminal in two rows.

Proposition 7.4 undecies improves the inheritance target: one need not
produce an idempotent for every color-count profile.  It is enough to
expose all states at one **oriented** radial vector.  What it does not
yet handle is a proper subset selected by multinomial carries, or a
leading sum in which \(\rho\) and its coordinate reversal coexist with
the same factorial weight.  Those are the only places where the
complete-union recombination (7.5w) is unavailable.

The exact selection sums, the persistent parity hole, and the two-row
terminal calculation are replayed by
`scripts/research_binary_gvc_equal_radial_union.py`.

The first coordinate-reversed width also closes exactly.

> **Proposition 7.4 duodecies (width-two reversal packets are
> terminal).**  Let
> \[
>  G(z)=\sum_{k=-2}^{2}c_kz^k
> \]
> and suppose both \(1\) and \(-1\) belong to
> \(\operatorname {Newt}(G)\).  Put
> \[
>  S_N=[z^N]G(z)^N+[z^{-N}]G(z)^N.                \tag{7.5z}
> \]
> Then \(S_1,\ldots,S_8\) cannot all vanish.  More precisely, after
> fixing the least and greatest active exponent, the four saturated
> endpoint charts
> \[
> (-1,1),\quad(-1,2),\quad(-2,1),\quad(-2,2)
> \]
> close after respectively \(2,4,4,8\) rows.

For the first chart the argument is visible without elimination:
\[
 S_1=c_{-1}+c_1,\qquad
 S_2=c_{-1}^2+c_1^2,
\]
so \(S_1=S_2=0\) forces both endpoints to vanish.  On each of the other
three charts, adjoin an inverse \(t\) to the endpoint product and compute
the exact rational Gröbner basis of
\[
 (S_1,\ldots,S_h,\ 1-tc_{\min}c_{\max}).
\tag{7.5aa}
\]
The unit first occurs at \(h=4,4,8\), respectively.  Thus every common
zero loses an endpoint and no longer contains both target slopes in its
Newton interval.

The calculation
`scripts/research_binary_gvc_reversal_union.py --prove-width-two`
replays the four saturated ideals.  It also exhausts the projective
five-coefficient space modulo \(5,7,11\); every point dies by row eight.
The modular census is only a regression, while the rational saturated
ideals prove the displayed width-two statement.

Both formerly unresolved mechanisms have the same algebraic shape.  A
coordinate reversal is a two-component trace, while a carry-selected
proper subset is a finite-character projection of a complete radial
union.  They can be encoded in a group algebra
\[
 \mathbb C[\mathbb Z^r\times F]
\tag{7.5ab}
\]
for a finite abelian group \(F\), followed by the coefficient of the
identity in both the free and torsion factors.  Duistermaat--van der
Kallen handles \(F=0\), and the characteristic-zero finite-group
constant-term kernel is a Mathieu subspace, but neither result alone
proves their tensor-product composite.  The next uniform target is
therefore a **torsion--torus trace lemma for the rank-one Cartesian
elements arising from (7.2)**.  Its \(F=C_2\), width-two case is
Proposition 7.4 duodecies.  An unrestricted statement for arbitrary
group-algebra elements is not asserted.

The two separate inputs are J. J. Duistermaat and W. van der Kallen,
[*Constant terms in powers of a Laurent
polynomial*](https://webspace.science.uu.nl/~kalle101/powers.pdf),
Indag. Math. 9 (1998), 221--231, and W. Zhao and R. Willems,
[*Analogue of the Duistermaat--van der Kallen theorem for group
algebras*](https://arxiv.org/abs/1009.5794),
Cent. Eur. J. Math. 10 (2012), 974--986.

The torsion--torus problem has an exact determinantal reduction.  Let
\[
 \mathcal A_F=\mathbb C[\mathbb Z^r\times F]
\]
for a finite abelian group \(F\), let
\(\operatorname {Reg}_F\) be the regular representation in the \(F\)
factor, and write \(\operatorname {CT}_{\mathbb Z^r\times F}\) for the
coefficient of the identity in both factors.

> **Proposition 7.4 terdecies (regular trace and logarithmic determinant
> reduction).**  For every \(u\in\mathcal A_F\) and \(N\geq1\),
> \[
>  \operatorname {CT}_{\mathbb Z^r\times F}(u^N)
>  =
>  \frac1{|F|}
>  \operatorname {CT}_{\mathbb Z^r}
>  \operatorname {Tr}\bigl(\operatorname {Reg}_F(u)^N\bigr).
> \tag{7.5ac}
> \]
> If
> \[
>  D_u(t)=\det\!\left(I-t\operatorname {Reg}_F(u)\right),
> \]
> then, as a formal power series,
> \[
>  \sum_{N\geq1}
>  \operatorname {CT}_{\mathbb Z^r\times F}(u^N)t^N
>  =
>  -\frac1{|F|}
>  \operatorname {CT}_{\mathbb Z^r}
>  \left(t\frac{\partial}{\partial t}\log D_u(t)\right).
> \tag{7.5ad}
> \]
> Moreover, if \(p\equiv1\pmod{\exp(F)}\), reduction modulo \(p\)
> satisfies
> \[
>  u^p=\phi_p(u),                                  \tag{7.5ae}
> \]
> where \(\phi_p\) multiplies every free exponent by \(p\), raises
> coefficients to their \(p\)-th powers, and fixes the torsion
> coordinate.

Left multiplication by an element of \(\mathbb C[F]\) has trace
\(|F|\) times its identity coefficient, which proves (7.5ac)
coefficientwise in the free lattice.  The standard identity
\[
 -t\partial_t\log\det(I-tA)
 =\sum_{N\geq1}\operatorname {Tr}(A^N)t^N
\]
gives (7.5ad).  Equation (7.5ae) is the freshman's-dream identity in the
commutative group algebra, together with \(f^p=f\) for
\(f\in F\) when \(p\equiv1\pmod{\exp(F)}\).

For the Hall/jet application, let \(e_\chi\) be the character idempotents
of \(\mathbb C[F]\).  The structured element has the form
\[
 u=
 \sum_{\chi\in\widehat F}
 e_\chi\,
 X^{-\rho_\chi}Y^{-\rho_\chi}
 \lambda(X)^dP(Y)^d.                              \tag{7.5af}
\]
Its identity coefficient is exactly the average of the complete radial
rows indexed by the finite carry or reversal characters.  Thus (7.5ad)
packages the previously inseparable profile sum into one scalar
log-determinant.

The trace can in fact be split without a new determinantal theorem.

> **Theorem 7.4 quaterdecies (finite-trace digit separation).**  Let
> \(K\) be a characteristic-zero field and let
> \(f_1,\ldots,f_q\) be Laurent polynomials in any finite number of free
> variables.  If
> \[
>  \sum_{i=1}^q\operatorname {CT}(f_i^N)=0
>  \qquad(N\geq1),                                 \tag{7.5ag}
> \]
> then
> \[
>  \operatorname {CT}(f_i^N)=0
>  \qquad(1\leq i\leq q,\ N\geq1).                 \tag{7.5ah}
> \]

Put \(a_i(n)=\operatorname {CT}(f_i^n)\).  It is enough to work over a
number field.  Indeed all coefficients and the equations (7.5ag) lie in
a finitely generated \(\mathbb Q\)-algebra.  If some \(a_i(m)\) were
nonzero, localize at it and take a closed point.  Zariski's lemma gives a
number-field specialization which preserves both (7.5ag) and that
nonzero value.  Clear the finitely many coefficient denominators.

Fix \(m\geq1\) and \(1\leq k\leq q\).  Choose an arbitrarily large
rational prime \(p\) which splits completely in the coefficient number
field and avoids the cleared denominators.  Such primes exist by
Chebotarev.  Require also that \(p\) exceed \(m\) times the absolute value
of every Laurent exponent coordinate occurring in the \(f_i\).  Set
\[
 N_{p,k}=m(1+p+\cdots+p^{k-1}).                   \tag{7.5ai}
\]
In the residue field, the freshman's dream gives
\[
 f_i^{N_{p,k}}
 =
 \prod_{j=0}^{k-1} f_i(z^{p^j})^m.               \tag{7.5aj}
\]
If exponents \(\beta_j\) selected from the \(j\)-th factor contribute
to the constant term, then
\[
 \beta_0+p\beta_1+\cdots+p^{k-1}\beta_{k-1}=0.
\]
Reducing coordinatewise modulo \(p\), the bound on \(p\) forces
\(\beta_0=0\).  Divide by \(p\) and repeat.  Thus every \(\beta_j=0\);
there are no signed base-\(p\) carries between the blocks.  Since the
prime splits completely, Frobenius fixes every residue coefficient, and
therefore
\[
 a_i(N_{p,k})\equiv a_i(m)^k\pmod p.             \tag{7.5ak}
\]
Equation (7.5ag) at \(N_{p,k}\) shows that
\[
 \sum_{i=1}^q a_i(m)^k\equiv0\pmod p.
\]
This holds at infinitely many completely split primes, so the algebraic
integer on the left is zero.  We have obtained the first \(q\) power
sums of \(a_1(m),\ldots,a_q(m)\).  Newton's identities force every
elementary symmetric function to vanish; hence
\[
 \prod_{i=1}^q(T-a_i(m))=T^q.
\]
Thus every \(a_i(m)=0\).  Since \(m\) was arbitrary, (7.5ah) follows,
contradicting the chosen specialization if the original conclusion had
failed.

This theorem is not restricted to Cartesian components.  Applying it
after diagonalizing the finite abelian factor proves the formerly
targeted statement.

> **Corollary 7.4 quindecies (torsion--torus trace separation).**  Let
> \(F\) be a finite abelian group and
> \(u\in K[\mathbb Z^r\times F]\).  If
> \[
>  \operatorname {CT}_{\mathbb Z^r\times F}(u^N)=0
>  \qquad(N\geq1),
> \]
> then every character component \(u_\chi\) satisfies
> \[
>  \operatorname {CT}_{\mathbb Z^r}(u_\chi^N)=0
>  \qquad(N\geq1).
> \]
> Consequently the coefficient-of-the-identity kernel in
> \(K[\mathbb Z^r\times F]\) is a Mathieu--Zhao subspace.

After adjoining the character values, Proposition 7.4 terdecies turns
the premise into the sum in (7.5ag).  Theorem 7.4 quaterdecies gives
componentwise pure vanishing.  Duistermaat--van der Kallen puts the
origin outside every \(\operatorname {Newt}(u_\chi)\).  Strict
separation then makes
\(\operatorname {CT}(v_\chi u_\chi^N)=0\) for every fixed multiplier
\(v_\chi\) and all large \(N\), uniformly over the finite character
set.  Averaging proves the Mathieu--Zhao assertion.

For the structured element (7.5af),
\[
 0\in\operatorname {Newt}(u_\chi)
 \quad\Longleftrightarrow\quad
 \rho_\chi\in d\operatorname {Newt}(\lambda)
 \cap d\operatorname {Newt}(P).                  \tag{7.5al}
\]
Therefore an active origin-containing component cannot participate in
an all-order zero trace.  A coefficient twist which removes a Newton
endpoint is support loss; otherwise strict separation is exactly the
free Laurent separator required by the Hall descent.

For clarity, a **scale-compatible** finite carry projection really is a
power trace.  Label the operator and polynomial channels by a homomorphism
\(\pi:\mathbb Z^{I}\oplus\mathbb Z^{J}\to F\), fix \(c\in F\), and let
\(\lambda_\chi,P_\chi\) denote the polynomials obtained by multiplying
each channel coefficient by the character of its label.  Character
orthogonality gives the exact restricted complete-union identity
\[
\begin{aligned}
 &\sum_{\substack{|p|=|q|=Nd,\ \rho(p,q)=N\rho\\
                   \pi(p,q)=Nc}}
   \binom{Nd}{p}\binom{Nd}{q}a^pb^q\\
 &\qquad =
 \frac1{|F|}\sum_{\chi\in\widehat F}
 \operatorname {CT}_{X,Y}
 \left(
  \chi(-c)X^{-\rho}Y^{-\rho}
  \lambda_\chi(X)^dP_\chi(Y)^d
 \right)^N.                                      \tag{7.5am}
\end{aligned}
\]
Finite unions of carry classes are obtained by adding copies of the
right side; coordinate reversal adds the two reversed radial copies.
Thus every packet already known to satisfy the displayed \(N c\)
congruence is a finite trace of the precise form covered by Theorem 7.4
quaterdecies.  Formula (7.5am) does not turn a fixed affine congruence
\(\pi(p,q)=c\) into a power trace: its Fourier factor
\(\chi(-c)\) then lies outside the \(N\)-th power.

Repeated digits give an exact description of what can cancel in this
affine case as well.

> **Theorem 7.4 sexdecies (weighted finite-trace classification).**
> Let \(w_1,\ldots,w_q\in K\) and let \(f_1,\ldots,f_q\) be Laurent
> polynomials over a characteristic-zero field.  Put
> \(a_i(n)=\operatorname {CT}(f_i^n)\).  If
> \[
>  \sum_{i=1}^q w_i a_i(N)=0\qquad(N\geq1),       \tag{7.5an}
> \]
> then, for every positive tuple \(m_1,\ldots,m_s\),
> \[
>  \sum_{i=1}^q w_i\prod_{j=1}^s a_i(m_j)=0.     \tag{7.5ao}
> \]
> Declare \(i\sim i'\) when
> \(a_i(n)=a_{i'}(n)\) for every \(n\geq1\).  For every equivalence
> class \(C\) whose common moment sequence is not identically zero,
> \[
>  \sum_{i\in C}w_i=0.                            \tag{7.5ap}
> \]
> Conversely, (7.5ap) for every nonzero class implies (7.5an).

For (7.5ao), specialize to a number field as in Theorem 7.4
quaterdecies and choose arbitrarily large completely split good primes
\(p\).  At
\[
 N=m_1+m_2p+\cdots+m_sp^{s-1},
\]
bounded signed digit uniqueness gives
\[
 a_i(N)\equiv\prod_{j=1}^s a_i(m_j)\pmod p.
\]
Equation (7.5an) and infinitely many such primes prove (7.5ao) in
characteristic zero.  There are only finitely many distinct vectors
\((a_i(1),a_i(2),\ldots)\).  Multivariate interpolation on finitely many
coordinates gives, for each nonzero vector, a polynomial with zero
constant term which is one on that vector and zero on all the others.
Expanding it into monomials and applying (7.5ao) proves (7.5ap).  The
converse is immediate.  Taking every \(w_i=1\) recovers Theorem 7.4
quaterdecies because a nonempty class has nonzero total weight in
characteristic zero.

The weighted alternative is real, even when every Newton interval
contains the origin.  Put
\[
 f_1=z+z^{-1},\qquad f_2=z^2+z^{-2}.
\tag{7.5aq}
\]
Dilation \(z\mapsto z^2\) gives
\[
 \operatorname {CT}(f_1^N)-\operatorname {CT}(f_2^N)=0
 \qquad(N\geq1),
\]
but for every odd \(N\),
\[
 \operatorname {CT}\!\left(z^{-1}f_1^N\right)
 -\operatorname {CT}\!\left(z^{-1}f_2^N\right)
 =\binom N{(N-1)/2}\ne0.                         \tag{7.5ar}
\]
Thus the coefficient kernel at a fixed nonidentity torsion class is not
Mathieu--Zhao.  Theorem 7.4 sexdecies says that its only all-order
mechanism is cancellation among character components with the same
complete Laurent period sequence.  Such isoperiodic components need not
be equal; the period-preserving Laurent mutations of M. Akhtar,
T. Coates, S. Galkin, and A. Kasprzyk,
[*Minkowski Polynomials and
Mutations*](https://arxiv.org/abs/1212.1785), give much larger families
in two or more free variables.

There is also a direct factorial warning.  Let
\[
 C(x,y)=y^2+4xy+2x^2,\qquad
 \mathcal L(x^ay^b)=a!b!.
\]
Then \(\mathcal L(C)=10\).  At \(p=11\) and repeated-digit index
\(N=1+11=12\), exact calculation gives
\[
 v_{11}\!\left(\mathcal L(C^{12})\right)=3,\qquad
 \frac{\mathcal L(C^{12})}{11^2}\equiv0\pmod {11},
\tag{7.5as}
\]
whereas naive Laurent digit factorization predicts
\(\mathcal L(C)^2\equiv1\pmod {11}\) after the common valuation two.
Equal-score radial and multinomial carries supply the extra
cancellation.  Hence one cannot apply Theorem 7.4 quaterdecies directly
to a Hall shell before proving that the radial factorial has become
common and that the remaining class is scale-compatible.

The exact regular trace, log determinant, \(C_2,C_3\) Frobenius checks,
and reversal diagonalization are replayed by
`scripts/research_binary_gvc_torsion_torus_trace.py`.  The signed
base-\(p\) digit factorization and the Newton endpoint are replayed by
`scripts/verify_binary_gvc_torsion_torus_digit_separation.py`.  The
weighted theorem, the affine \(C_2\) obstruction, and (7.5as) are
replayed by
`scripts/verify_binary_gvc_weighted_trace_obstruction.py`.

The maximal two-endpoint packet from Corollary 7.3 is nevertheless
terminal once it is exposed at all orders.

> **Theorem 7.5 (opposite three-by-three packet is terminal).**  Keep
> operator endpoints \(A=(r,0)\), \(B=(0,s)\).  Let \(h_A,h_B\) be
> disjoint primitive circuits, where \(h_A\) uses \(d_A\) copies of \(A\)
> and three polynomial channels with positive multiplicities
> \(q^A_1,q^A_2,q^A_3\), and \(h_B\) is defined similarly.  If the return
> semigroup packet is generated by \(h_A,h_B\), its factorial-weighted
> partial moments cannot vanish at every positive order when both
> circuits are active.

Absorb the active coefficients of each circuit into nonzero scalars
\(u,v\).  Since the two polynomial triples are disjoint, the normalized
packet generating series factors:
\[
 \sum_{m\geq0}\frac{\mathcal M_m}{(m!)^2}t^m
 =F_A(t)F_B(t),
\tag{7.6}
\]
where
\[
\begin{aligned}
 F_A(t)&=\sum_{N\geq0}
 \frac{(Nd_Ar)!}
 {(Nd_A)!\prod_{j=1}^3(Nq^A_j)!}\,
 u^Nt^{Nd_A},\\
 F_B(t)&=\sum_{N\geq0}
 \frac{(Nd_Bs)!}
 {(Nd_B)!\prod_{j=1}^3(Nq^B_j)!}\,
 v^Nt^{Nd_B}.
\end{aligned}
\tag{7.7}
\]
Here \(\sum_jq^A_j=d_A\) and \(\sum_jq^B_j=d_B\).
If every positive packet moment vanished, then \(F_AF_B=1\).
When \(d_A\ne d_B\), the coefficient in the smaller positive degree is
already nonzero.  Hence \(d_A=d_B=d\).

Write \(A_N,B_N\) for the positive rational coefficients in (7.7) before
\(u^N,v^N\).  The coefficient of \(t^d\) gives
\(A_1u+B_1v=0\).  (If either \(u\) or \(v\) were zero, this equation
would force both to be zero.)  After this substitution, the coefficient
of \(t^{2d}\) is
\[
 A_1^2u^2(R_A+R_B-1),
\qquad
 R_A=\frac{A_2}{A_1^2}
 =\frac{\binom{2dr}{dr}}
 {\binom{2d}{d}\prod_{j=1}^3\binom{2q^A_j}{q^A_j}},
\tag{7.8}
\]
with the analogous formula for \(R_B\).

If \(r=1\), then
\[
 R_A=\frac1{\prod_j\binom{2q^A_j}{q^A_j}}\leq\frac18.
\]
If \(r\geq2\), then \(R_A>1\).  Indeed the central binomial coefficients
satisfy
\[
 \binom{2(a+b)}{a+b}>
 \binom{2a}{a}\binom{2b}{b}
\qquad(a,b>0)
\tag{7.9}
\]
by Vandermonde's identity, so
\(\binom{2d}{d}>\prod_j\binom{2q^A_j}{q^A_j}\), and therefore
\[
 \binom{2dr}{dr}\geq\binom{4d}{2d}>
 \binom{2d}{d}^{\,2}>
 \binom{2d}{d}\prod_j\binom{2q^A_j}{q^A_j}.
\]
Thus \(R_A+R_B\leq1/4<1\) when \(r=s=1\), while it is \(>1\) if either
endpoint order is at least two.  The coefficient (7.8) never vanishes.
The exact coefficient identity and bounded central-binomial regression
are replayed by
`scripts/research_binary_gvc_eight_obstructions.py
--verify-opposite-packet`.

The finite-trace separation theorem finishes the descent once the
following promotion hypothesis has been established.

> **Conditional Theorem 7.6 (scale-compatible binary circuit
> exposure).**  Suppose all pure moments (7.2) vanish.  Assume at each
> Hall--jet descent stage the following property:
>
> **(SC)** A tied initial shell which does not already lose support or
> have a separator exposes, at every scale \(N\), a nonempty packet with
> one common radial factorial.  After removing that factorial, the
> packet is a fixed finite union of complete radial rows and congruence
> fibers
> \[
>   \pi_\nu(p,q)=Nc_\nu
> \]
> with \(F_\nu,\pi_\nu,c_\nu\) independent of \(N\) and of the good
> prime used to expose the shell.  Pure vanishing makes the normalized
> packet sum zero in characteristic zero for every \(N\).
>
> Then, after deleting channels absent from \(\mathsf R\), iteration
> produces one of:
>
> 1. a componentwise least radial vector;
> 2. a split-symbol separator with a linear multiplier gap;
> 3. a reduced two-operator packet with at most five polynomial channels;
> 4. a strict loss of operator or polynomial support.

Under (SC), formula (7.5am) writes the entire scale family as
\[
 \sum_{\nu=1}^s\operatorname {CT}(f_\nu^N)=0
 \qquad(N\geq1)                                  \tag{7.10}
\]
for finitely many Laurent polynomials \(f_\nu\).  Here the index \(\nu\)
includes the finite carry class, its Fourier character, and the possible
coordinate reversal.  Multinomial weights are already the coefficients
of the powers in (7.5am), while the radial factorial is common to the
tied class and has been divided out.

Theorem 7.4 quaterdecies splits (7.10) componentwise.  Duistermaat--van
der Kallen then gives one of two outcomes for every component: its
origin-containing Newton support has disappeared, which is strict
support loss, or its support has a strict free Laurent separator.  If no
carry twist changes the support, summing the characters back gives the
complete oriented radial union of Proposition 7.4 undecies as the
factorial-compatible pure-zero face.  Coordinate reversal is split in
the same application, without a width bound.  Each iteration removes a
channel, produces a separator, or passes to a proper face of a finite
return cone, so it terminates.

Every terminal outcome in Conditional Theorem 7.6 is already closed:
Theorem 4.3 handles the first, Theorems 3.1--3.2 and the homogeneous
split-symbol theorem handle the second, and Theorems 5.9--5.14 handle the
third; Theorem 7.5 handles the maximal six-channel packet which can occur
immediately before it.  The fourth permits induction.

> **Corollary 7.7 (conditional binary GVC).**  If (SC) holds for every
> Hall-reduced binary return system, then for every characteristic-zero
> field, every
> \(\Lambda\in k[\partial_x,\partial_y]\), and all \(P,Q\in k[x,y]\),
> \[
>  \Lambda^m(P^m)=0\quad(m\geq1)
>  \quad\Longrightarrow\quad
>  \Lambda^m(QP^m)=0\quad(m\gg0).
> \]

The standard translation/polarization reduction converts every
coefficient of the polynomial pure moments into the scalar return
expansion (7.2).  Conditional Theorem 7.6 terminates every such return
expansion under (SC).
The separator outcomes have a linear gap, so the bounded support of the
fixed multiplier \(Q\) cannot repair them for large \(m\).  Taking the
maximum threshold over the finitely many coefficients of \(Q\) proves
the conditional conclusion.

> **Open Problem 7.8 (affine-carry promotion).**  Prove (SC), or produce
> a rank-one Cartesian Hall packet in which an affine carry class
> survives as a nonterminal isoperiodic character cluster.

There is a sharper attack on the isoperiodic alternative which uses the
translation relations discarded by an abstract Laurent trace.  Write
\[
 P(z+Y)=\sum_{\beta\in\mathcal D_P}p_\beta(z)Y^\beta,
 \qquad
 p_\beta(z)=\frac{\partial^\beta P(z)}{\beta!}.
\tag{7.11}
\]
The derivative support \(\mathcal D_P\) is a finite down-set, and
\[
 \partial_{z_i}p_\beta=(\beta_i+1)p_{\beta+e_i}.
\tag{7.12}
\]
For two nonvanishing character twists \(c_\beta,c'_\beta\), put
\(r_\beta=c'_\beta/c_\beta\).  Call the relative twist **flat** when
there are roots of unity \(\eta_1,\eta_2\) such that
\[
 r_{\beta+e_i}=\eta_i r_\beta
\tag{7.13}
\]
on every active edge of \(\mathcal D_P\).

> **Lemma 7.9 (flat translation twists are scale-compatible).**  On
> each connected derivative support, (7.13) implies
> \(r_\beta=C\eta^\beta\).  Consequently, if
> \[
>  H(Y)=\sum_\beta c_\beta p_\beta(z)Y^\beta,\qquad
>  H'(Y)=\sum_\beta c'_\beta p_\beta(z)Y^\beta,
> \]
> then \(H'(Y)=C H(\eta_1Y_1,\eta_2Y_2)\), and
> \[
>  [Y^{N\rho}]H'(Y)^N
>   =(C\eta^\rho)^N[Y^{N\rho}]H(Y)^N.            \tag{7.14}
> \]
> The same conclusion, with \(\rho\) reversed, holds after a coordinate
> reversal.  Thus scalar, torus, and reversal collisions contribute
> only \(N\)-dependent character phases and belong to the
> scale-compatible fibers of (7.5am).

The proof is path propagation on the connected down-set; (7.14) is
coefficient extraction.  This elementary lemma identifies the right
rigidity statement.  A nonflat ratio on a two-dimensional down-set has a
minimal witness on a collinear three-level path or a unit square: the
corresponding multiplicative second difference or square curvature is
not one.  These are exactly the Bessel triple and beta bridge supports
closed in Corollary 7.4 quinquies when isolated.

> **Translation-curvature target.**  Let an isoperiodic character class
> of a rank-one Cartesian Hall packet remain after support loss and
> Laurent separation.  Prove that either its relative twists are
> scalar--torus--reversal equivalent on the active Taylor down-set, or
> some derivative in \(z\), applied to the all-order period identity,
> exposes a nonflat three-level or square curvature block with the same
> high-digit quotient.

This target would close Open Problem 7.8.  In the flat case Lemma 7.9
gives (SC).  In the nonflat case the demanded common-quotient exposure
feeds the already proved augmented atom determinant.  Differentiating
the period identity is the natural mechanism because (7.12) replaces
one Taylor channel by its adjacent channel without changing the
high-digit factor.  What is not yet proved is that the resulting
inserted identities can be triangularized before different insertions
cancel.

The subsequent
[translation-tangent theorem](BINARY_GVC_TRANSLATION_TANGENT_RIGIDITY.md)
proves the primitive one-direction linearized case exactly.  For
\(H(t)=(1+t)^d\) and \(\gcd(d,r)=1\), every polynomial
\(B(t)=\sum_jb_j\binom dj t^j\) satisfying
\[
 [t^{rN}]B(t)H(t)^{N-1}=0\qquad(N\geq1)
\]
has \(b_j=c(j-r)\); this is precisely the tangent to the flat scalar--torus
action.  An iterated cyclotomic-neighbourhood argument then proves that every
\(q^a\)-order pair of Taylor twists with identical complete moving rows is
flat once the underlying prime \(q\) is sufficiently large.  The same note
also proves that no
coefficient-blind Hilbert-module inheritance theorem can replace the
Cartesian argument: two free module translates carrying a common factorial
weight can cancel isoperiodically at every pure order while a fixed multiplier
survives.  Thus the remaining triangularization must use higher
cyclotomic-adic rows for exceptional small primes, mixed-prime torsion, and
the genuinely two-dimensional Taylor down-set, always over one common high
quotient.

An exact bounded search tests the smallest version of this target.  For
one translated monomial and a \(C_2\) twist, set
\[
 h_\epsilon(t)=\sum_{j=0}^d
   \epsilon_j\binom dj t^j,\qquad
 a_{\epsilon,r}(N)=[t^{rN}]h_\epsilon(t)^N .
\tag{7.15}
\]
After quotienting scalar sign, \(t\mapsto-t\), and coefficient reversal,
all sign twists for \(2\leq d\leq12\), all
\(1\leq r<d\), and all rows \(N\leq2d+4\) give no unexplained equal
period pair.  This is 8,188 twists and 81,924 exact moving rows.
The first apparent collisions recover the expected symmetries:
\[
\begin{aligned}
 (1,-1,-1,1)&\longleftrightarrow(1,1,-1,-1)
 &&(t\mapsto-t,\ d=3,r=2),\\
 (1,-1,-1,-1,1,-1,1)&\longleftrightarrow
 (1,-1,1,-1,-1,-1,1)
 &&(\text{reversal},\ d=6,r=3).
\end{aligned}
\]
The command is
`scripts/search_binary_gvc_translation_isoperiodic_twists.py
--max-degree 12 --extra-depth 4`.  It is evidence for the target, not a
proof.  The optional two-direction run
`--max-degree 3 --extra-depth 4 --rectangles` also checks every \(C_2\)
twist on the \((2,2)\) and \((3,2)\) binomial Taylor rectangles: 2,304
twists and 4,352 moving rows through depths 12 and 14 give no collision
outside scalar, two-torus, reversal, and coordinate-exchange symmetry.
This remains bounded evidence.  General Laurent moment rigidity cannot
be assumed: the Laurent
moment problem already has non-compositional phenomena, as shown by
F. Pakovich, C. Pech, and A. Zvonkin,
[*Laurent polynomial moment problem: a case
study*](https://arxiv.org/abs/0910.2691), and period-preserving mutations
give further isoperiodic families.

This is now the exact binary gap.  Lemma 7.4 ter controls one prime and
one order \(\ell d\); its carry rows factor into a low-digit insertion
and a generally unbalanced high-digit quotient.  It does not show that
these rows, as \(N\) varies, have congruence \(\pi(p,q)=Nc\), nor that
the finite group and the packet are independent of the exposing prime.
For a fixed affine residue the character weights remain outside the
power, and Theorem 7.4 sexdecies permits cancellation precisely along
isoperiodic clusters.  Equation (7.5as) independently shows that the
Laurent repeated-digit congruence cannot simply be moved through the
radial factorial.

For context, the still-unproved word is
**factorial-compatible**.  Polyhedral decomposition
alone supplies the circuits in Corollary 7.3, but does not imply that the
corresponding partial sum of (7.2) vanishes.  Example 5.5 rules out
deriving that inheritance from the first three moments.  Lemma 7.4 shows
that good-prime valuation separates unequal radial profiles, confining
the gap to equal-profile circuit mixtures.
For two operator endpoints, Theorem 7.5 closes the maximal opposite
three-by-three packet once it is exposed.  With more operator endpoints,
the same four-variable circuit bound remains valid.  Theorem 7.4
quaterdecies supplies component separation once such a packet has been
promoted to (7.10); it does not supply the promotion itself.

Lemma 7.4 ter rules out one tempting shortcut: the lowest same-prime
initial form is always the Frobenius image of a lower moment, so it
contains no new exposure information.  It replaces that failed route by
the exact jet--carry filtration (7.5a).  At its first new level,
Corollary 7.4 sexies leaves only centered triples and two-by-two bridge
atoms, while (7.5d)--(7.5f) package their full sum as one- and two-ghost
cocycles.  Proposition 7.4 quater shows that their universal diagonal
zeros are exactly the known centered and Hall relations, but also rules
out a fixed-prime converse.  Formula (7.5l) further rules out
single-row cross-prime avoidance.  The immediate finite lemma is
therefore block triangularity after adjoining the ordinary beta moment
row and the centered Bessel endpoint row to the ghost diagonals.  One
must then incorporate the next unit correction of the Frobenius terms
under the same Hall support split.  Corollary 7.4 quinquies proves the
isolated atom blocks terminal.  Their compatibility with the common
high-digit quotient remains the promotion problem.  Proposition
7.4 septies rules out
reducing that compatibility to circuit moves on arbitrary projected
supports: the first missing block is the support-five identity (7.5m).
The projected-scroll Gröbner degree bound proves that only finitely many
quotient blocks occur for a fixed symbol support.  Proposition 7.4
nonies closes every *whole* scroll fiber after exposure, including the
minimal non-universal-Gröbner \(S(6)\) and \(S(5,4)\) fibers.  The
finite-trace theorem bypasses the block-by-block census only after the
shell has become a fixed scale-compatible trace: it then separates the
complete color-count/radial-profile components simultaneously.

The subsequent
[prime-power tomography census](BINARY_GVC_PRIME_POWER_TOMOGRAPHY.md)
computes this projected hierarchy through radial span seven.  At
\(p=5,7,11,13\), \(e\leq3\), its finite adelic scalar collisions are exactly
the all-scale factorial-partition collisions; there are no accidental
finite-window collisions.  Its only fully decorated \(C_2,C_3\)-blind moves
are the span-seven orbits of the exact all-span family
\[
 R_{s+6}B_aB_{a+1}=R_sB_{a+3}B_{a+4}.
\]
Every such move is Graver primitive and has a two-state exact projected fibre.
Its factorial weights, low digits, and Kummer data agree at every scale, and
its \(C_2,C_3\) marked-character traces also agree.  A \(C_4\) character
always separates its primitive endpoints.

The same note now closes the family after every fixed finite-character
promotion, even if the character groups the endpoints.  At scale \(N\), the
complete \(C_2,C_3\)-blind fibre is
\[
 z_t=(N-t,t,N-t,N-t,t,t),\qquad0\leq t\leq N,
\]
and its normalized row is
\[
 \binom{2N}{N}\sum_t\binom Nt^3U^{N-t}V^t.
\]
If a further character has relative order \(h\), its endpoint class at
scales \(h,2h\) gives respectively
\(U^h+V^h\) and
\(U^{2h}+\binom{2h}{h}^3U^hV^h+V^{2h}\).  Their common zero in
characteristic zero is \(U=V=0\), hence support loss.  Thus the six-step
packet is character-separated, terminal, or loses support after fixed
promotion.

The subsequent
[nonfree-factorization theorem](BINARY_GVC_NONFREE_FACTORIZATION_TOMOGRAPHY.md)
removes the remaining abstract mixed-semigroup ambiguity at every fixed
radial span.  If \(q=\lceil s/2\rceil\), the two complete marked histograms
for \(C_q,C_{q+1}\) are injective on levels \(0,\ldots,s\) when \(s\) is
odd.  When \(s\) is even, their kernel is one cycle, and the corresponding
two-colour return is a conformal sum of the already-safe beta swaps
\[
 R_iB_{i+q+1}=R_{i+q+1}B_i,\qquad0\leq i<q.
\]
It follows atomwise that every nonfree factorization collision with equal
marked signatures is safe.  The exact Hilbert census through span six
independently confirms the smaller \(C_2,C_3,C_4\) signatures and exhibits
the first factorial-only square.

What remains is therefore only the Hall-specific problem of
factorial-compatibly inheriting one fixed marked packet from the
prime-dependent affine shell.  A module-only version is false; the corrected
statement must use the complete Cartesian translation tower.  Its primitive
one-direction tangent and every large-underlying-prime, prime-power-character
case are now proved, leaving exceptional small primes, mixed-prime torsion,
and two-dimensional curvature triangularization over a common high quotient.
Neither the censuses nor the
new tangent theorem supplies that final promotion or a GVC counterexample.

The toric input used here is S. Petrović,
[*On the universal Gröbner bases of varieties of minimal
degree*](https://arxiv.org/abs/0711.2714),
Math. Res. Lett. 15 (2008), 1211--1221, together with
T. Bogart, R. Hemmecke, and S. Petrović,
[*Equality of Graver bases and universal Gröbner bases of colored
partition identities*](https://arxiv.org/abs/1004.0840).

A new preprint of M. Wilson,
[*A face-isolation proof of the two-variable Gaussian Moments
Conjecture*](https://arxiv.org/abs/2607.23887), uses the same
good-prime/Frobenius-to-Laurent strategy in one complex weight
coordinate.  It does not prove the binary GVC statement here.  Under
the Gaussian realization of a binary differential pairing there are
two complex weight coordinates and the radial factor is
\(\rho_1!\rho_2!\), not one factorial.  Wilson's own higher-pair
discussion identifies exactly the resulting obstruction:
Duistermaat--van der Kallen only puts \(0\) outside a
two-dimensional Laurent Newton polygon, which is weaker than
one-sidedness, and the weight-zero part returns to the multivariate
Factorial Conjecture.  The paper therefore validates the face-isolation
mechanism and the algebraic-specialization step.  The repeated-digit
trace separation above supplies the additional finite-profile
idempotence which that one-radial argument does not contain.
