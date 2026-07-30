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

They prove two pieces of the proposed binary termination argument:

1. every leading pure-zero component is localized by one root
   multiplicity, with no enumeration of global binary root partitions; and
2. every genuine unequal-weight common-threshold face is automatically
   terminal.

The remaining step is an existence theorem for a common threshold along a
defect path whose depth is proportional to \(m\).  Thus this note does not
prove unrestricted \(\operatorname{GVC}(2)\).

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

The exact remaining target can therefore be stated without filtration
language:

> **Beta--torus Hall problem.**  For the rank-one Cartesian polynomials
> \(G\) produced by a binary operator symbol and a translated polynomial,
> prove that all moments in (5.7) can vanish only if the angular support is
> strictly separated or the active Hall/jet support drops.

Theorem 5.1 proves this when there is no genuine \(S\)-convolution.
Prime-endpoint Bessel rigidity proves the adjacent three-level special
case.  The general positive-density multi-level case is not a consequence
of either theorem.

Primary references for Theorem 5.1 and its moment input are:

- D. Liu and X. Sun,
  [*The Factorial Conjecture and Images of Locally Nilpotent
  Derivations*](https://doi.org/10.1017/S0004972719000546),
  Bull. Aust. Math. Soc. 101 (2020), 71--79;
- J. P. Françoise, F. Pakovich, Y. Yomdin, and W. Zhao,
  *Moment vanishing problem and positivity: some examples*,
  Bull. Sci. Math. 135 (2011), 10--32.

## 6. What remains of termination

The general fixed-depth theorem kills every bounded filtration defect.
Theorem 2.1 and (2.5) reduce every leading component to one root
multiplicity and two normalized jet wings.  Corollary 4.1 kills every
defect path that reaches one fixed common threshold.

The only remaining architecture is therefore:

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

A complete binary termination theorem must solve the beta--torus Hall
problem above: the positive-density coupled convolution between
incomparable Pareto endpoints must either acquire a least radial vector
after a **factorial-compatible** Hall/jet reduction, where the multiradial
prime theorem applies, or have smaller Hall/jet support after
cancellation.  The degree-at-most-six calculations verify exactly this
alternative in their bounded support ranges.  Theorem 5.1 removes every
one-channel ordinary-homogeneous tie in all degrees, but (5.4) and (5.9)
show why it does not supply the missing multi-channel induction.  That
uniform beta--torus statement remains open, so unrestricted
\(\operatorname{GVC}(2)\) is not proved here.
