# Radial cores and one-sided Gaussian phases: finite images and scalar rigidity

Let K have characteristic zero. Use one Witt Gaussian pair Z,W with
E[ZW]=1 and an independent ordinary Gaussian Y with E[Y^2]=1; all other
basic pairings vanish. Put S=ZW and consider

    P=P0+R, P0=g(S)+B(S)Y+aY^2,                     (R1)

where a is constant, g,B are arbitrary polynomials, and every monomial
of R has strictly positive Z-minus-W angular weight. Strictly negative
support is included by reversing the weights. There is no degree bound
on R, including its W or Y degree. The case R=0 is included.

For ORIGINAL fixed spatial markers define

    Phi_P[Q]=E[Q exp(tP)], V_P=Phi_P[K[Z,W,Y]],
    M_P=K[t]V_P.

The following are equivalent:

1. The actual vacuum is p(t)exp(psi(t)), with p in K[t], p(0)=1 and
   psi in tK[t].
2. The full original image M_P is finite over K[t].
3. a=0 and g,B are constant.

In every such case, writing g=g0 and B=b,

    p=1, psi=g0t+b^2t^2/2, M_P=exp(psi)K[t].       (R2)

Thus the original vacuum is primitive. For nonconstant P one also has
literally V_P=M_P. The necessary constant-phase exception is R=0,b=0:
then V_P=K exp(g0t), whereas M_P=exp(g0t)K[t].

Within this WHOLE one-sided family, purity and safety are equivalent.
Both hold exactly when

    a=0, g=0, B=0, so P=R.                         (R3)

The result is scoped to (R1), with constant a and a core linear in Y
apart from its constant quadratic coefficient. It is not a universal
Gaussian safety-to-purity or finite-vacuum-primitivity theorem.

## 1. The actual radial scalar source retains its endpoint

Give Z,W,Y weights 1,-1,0. Since P0 has weight zero and R is strictly
one-sided, every vacuum term containing R has nonzero angular weight.
Consequently the ACTUAL scalar periods agree:

    Phi_P[1]=Phi_P0[1].                             (R4)

This uses original Wick contraction and no localization or alternate
marker map. The stronger equality M_P=M_P0 is proved in the
[angular-core full-image theorem](GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md).
The finite-image argument below also supplies its sufficient all-marker
direction directly.

Work first over C. The original radial moments are E_S[S^m]=m!. Put
d=1-2at. Coefficientwise ordinary-Y completion gives

    Phi_P0[1]=H E_S exp(phi+S),
    H=d^(-1/2), phi=-S+t g(S)+t^2 B(S)^2/(2d).       (R5)

The branch H(0)=1 is selected. No convergence of an integral over the
positive real axis is presumed. Its Gamma-integral notation represents
the factorial moments with endpoint S=0.

Assume statement 1, and normalize

    phi*=phi-psi, H*=d^(-1/2)/p,
    L(q)=H* E_S[q exp(phi*+S)].                     (R6)

Then L(1)=1. For every q in C(t)[S] its value is a Laurent-t series
with finite principal part: phi*+S is t-divisible and every coefficient
uses only finitely many original factorial moments. The exact identity

    E_S[f'-f]=-f(0)                                (R7)

implies L((partial_S+phi*_S)f)=0 for f in S C(t)[S]. Time differentiation
gives

    partial_t L(q)
      =L(partial_t q+(phi*_t+H*'/H*)q).             (R8)

Thus L is an ACTUAL horizontal unit functional on the relative module

    Mrel=C(t)[S]/(partial_S+phi*_S)(S C(t)[S]).      (R9)

The condition f(0)=0 is essential. On all polynomials the endpoint
term in (R7) survives, so the absolute de Rham quotient would be wrong.

## 2. An endpoint-sensitive leading-coefficient bound

We give the relative argument explicitly. Let varphi in C(t)[x] have
degree D>=2 and top coefficient a_D. Let h be a nonzero algebraic scalar
with h'/h rational. Normalize v(r)=1 at r=1/t=0. Then

    D v(h)>v(a_D)                                  (R10)

excludes a rational horizontal quotient of the differential seed span
in the pair (A1,{0}) taking the seed h dx to one. Only convergent
meromorphic algebraic parameter coefficients are used in this argument;
no arbitrary nonconvergent formal extension is asserted.

### The relative rank and its boundary class

Write A=C(t)[x] and Dx=partial_x+varphi_x. The relative complex is
xA --Dx--> A dx, and its quotient is A/Dx(xA). Its connection is
q -> q_t+(varphi_t+h'/h)q. Since Dx(xq) has degree deg q+D, the
relative rank is D, with basis 1,x,...,x^(D-1). The absolute quotient
A/Dx A has rank D-1. Injectivity of Dx on polynomials gives the exact
sequence of differential modules

    0 -> C(t)b -> Mrel -> Mabs -> 0,
    b=[varphi_x],
    nabla b=(varphi_t(t,0)+h'/h)b.                  (R11)

The last identity follows because Dx(varphi_t-varphi_t(t,0)) is a
relative exact term. The boundary class and its parameter action are
therefore present throughout the argument.

At an ordinary parameter, add one path from x=0 to a decay sector at
infinity to a perfect absolute rapid-decay frame. Its value on b is
-h exp(varphi(t,0)), by the fundamental theorem of calculus, while all
closed absolute cycles vanish on b. These paths give a full rank-D
relative perfect pairing and it is horizontal under flat transport.
The absolute pairing is the retained ordinary exponential comparison;
the extra path establishes the relative pairing without discarding
the endpoint or assuming the extension (R11) splits.

### A polar endpoint action

Put a0=varphi(t,0). If a0 has a nonzero polar part, the rank-one boundary
module in (R11) has no zero-exponential part. After finite ramification,
exponential decomposition is exact, so the zero-exponential part of the
relative module identifies with that of the absolute module, with the
same seed image. The retained absolute all-period leading bound gives
seed order at least v(h)-v(a_D)/D>0 there. This argument is applied
after the scalar gauge; a gauge making a0 regular belongs to the next
case instead.

### A regular endpoint action and its additional period

Suppose a0 has no polar part. Write varphi=a0+sum_(j>=1) a_j x^j and set

    rho=max_(a_j!=0)(-v(a_j)/j), delta=r^rho,
    F_r(u)=varphi(t,delta u)-a0.

After finite ramification, every coefficient of F_r is integral and
at least one is a unit. Its limit F(u) is nonconstant and F(0)=0.
Moreover

    v(h)+rho>=v(h)-v(a_D)/D>0.                     (R12)

The endpoint stays u=0 and the volume is h delta du.

Choose a sufficiently narrow parameter sector and a descending ray in
the phase plane avoiding the finitely many limiting directions of the
algebraic critical values. Identically zero critical values only meet
the starting endpoint. A local inverse branch from u=0, allowing finite
ramification in phase, continues along that ray and reaches an actual
decay sector without crossing a nonzero critical value. Branch choices
are finite and continuous on the sector; this is a flat relative path
with nonzero endpoint boundary.

On a fixed u-disc, its phase curve has uniformly bounded algebraic
degree. The retained bounded-degree curve-length estimate gives bounded
arc length; monotonic phase prevents repeated traversal. Each bounded
piece therefore contributes O(h delta), also when critical points
collide with one another or with the endpoint.

Let q=deg F and choose a fixed large radius A0 so the leading monomial
of F dominates outside that radius. Every additional higher coefficient
of F_r has positive order. Choose sigma>0 smaller than each such order
divided by its degree excess over q. Then throughout

    A0<=|u|<=|r|^(-sigma),

one has |F_r(u)| comparable to |u|^q. The bounded-degree length estimate
on dyadic annuli, multiplied by exp(-c|u|^q) on the descending path, is
uniformly summable. Thus this portion is also O(h delta). A branch
cannot escape to a far root without first crossing these annuli, where
its phase already has this large size.

At |u|=|r|^(-sigma), the phase size is at least c|r|^(-sigma q).
Beyond that cutoff, polynomial root bounds and the fixed-degree inverse
phase curves give at most algebraic parameter and polynomial phase
growth. The remaining exponential tail dominates both. The complete
endpoint path consequently satisfies

    Iendpoint=exp(a0) O(h delta).                  (R13)

Finite logarithmic factors preserve the strict positive power in (R12).
Further finite local zooms around endpoint collisions insert positive
displacement factors; the relative endpoint section is allowed because
the volume is regular there.

### The full zero-exponential frame and the imported comparisons

The closed absolute frame plus this endpoint path is the full relative
frame. The absolute zero-exponential seed periods have positive order.
On a generic parameter ray, (R13) rules out any growing nonzero-action
component of the endpoint period. Remaining nonzero-action components
are exponentially small; removing them leaves its zero-exponential
part with the same positive power, allowing finite logarithmic factors.
Coincident zero actions and regular extensions do not remove the
positive powers. If a0 is polar use the previous action-block argument.

Equivalently, include the endpoint section in the full proper marked
model and retain the relative chain cone. Nearby-cycle comparison and
proper direct image commute with this cone, keeping every endpoint,
annulus and gluing map. Finite/finite necks carry a positive power times
a logarithm, and the positive finite-side factor survives descending
integration at finite/pole nodes. Thus no crossover or endpoint term
is omitted from the zero-exponential calculation.

A rational unit quotient on the seed span extends at an ordinary
parameter to a full analytic horizontal functional and is represented
by the perfect relative period pairing. Its seed value would be one,
contradicting the positive-order bound. This proves (R10). The boundary
class b itself can have a constant period; that does not give a unit
period on the seed [1].

The imports are precise: the absolute all-period bound is
[Formal stationary balls](../../../results/gaussian-boundary/polynomial_stationary_balls.tex),
Theorem `thm:stationaryballs`; ordinary seed extension is
[seed-period comparison](../../../results/gaussian-boundary/seed_period_comparison.tex),
Lemma `lem:seedfunctional`; bounded-degree arc length is
[the analytic core](../../../results/gaussian-boundary/analytic_core.tex),
Lemma `lem:curvelength`. The three exponential comparison facts and their
primary-source audit, the complete compact/neck/tail estimates, and the
relative chain-cone formulation are retained in the canonical
[H21 fixed-puncture comparison](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#fixed-puncture-logarithmic-seed-bound),
[all-centres and crossover estimates](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#all-centres-punctures-and-crossover-estimates),
and [relative leading-coefficient proof](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#a-relative-leading-coefficient-bound).
The argument above supplies the endpoint frame; it does not relabel
an absolute theorem as a relative theorem without that construction.
These comparison dependencies retain their existing review limits.

## 3. Higher radial degree contradicts the scalar premise

Suppose D=max(deg g,2deg B)>=2, interpreting the degree of zero as
-infinity. This is the actual spatial degree of phi* in (R6), since

    a_D=t g_D+t^2 C_D/(2d), C_D=[S^D]B(S)^2          (R14)

cannot vanish identically when g_D,C_D are not both zero: multiplication
by d gives t g_D+t^2(C_D/2-2a g_D). The scalar psi changes no spatial
coefficient.

If a=0, v(a_D) is -1 or -2, whereas v(H*)=deg p>=0. If a!=0, then

    a_D=[g_D(r-2a)+C_D/2]/[r(r-2a)].

Its order is -1, or zero if the constant numerator cancels; that case
leaves the nonzero numerator g_D r. Meanwhile v(H*)=1/2+deg p. Thus
(R10) is strict in every higher-degree case, including this leading
cancellation and every polynomial p. The ACTUAL relative functional
(R6)-(R9) supplies precisely the forbidden unit seed quotient. Hence
D>=2 is impossible under statement 1.

## 4. The remaining degree-one phase is primitive

The term -S in phi* makes its actual degree at least one. The remaining
case is exactly B=b constant and g=g0+g1S. Its exact original vacuum is

    f0=exp(A)/(sqrt(d)(1-g1t)),
    A=g0t+b^2t^2/(2d).                              (R15)

Therefore f0=p exp(psi) implies

    exp(2(A-psi))=p^2 d(1-g1t)^2.

The exponent is rational. At a finite pole its derivative has pole
order at least two, whereas the logarithmic derivative of a rational
function has at most simple poles. Thus the exponent has no finite
poles and is polynomial. Comparison at infinity makes it constant;
its initial value is zero. Consequently A=psi and the polynomial
identity p^2 d(1-g1t)^2=1 forces

    p=1, a=g1=0, psi=g0t+b^2t^2/2.                  (R16)

This proves statement 1 implies statement 3. Statement2 implies statement 1
by the [general finite actual-image theorem](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#7-every-finite-actual-image-and-its-original-fixed-marker-image).
All comparison arguments work over C. Embed the finitely generated
field of phase coefficients and hypothetical p,psi coefficients into C
to obtain necessity over every characteristic-zero field. The original
factorial source and the sufficient formulas below are algebraic over K.

## 5. Every original marker in the finite branch and its cutoff

Assume a=0,g=g0,B=b. If R!=0 let gamma>=1 be its least positive angular
weight, after reversing negative support. For R=0 set N_D=0; otherwise
put N_D=floor(D/gamma). For an ORIGINAL marker Q of degree at most D,
angular balance and ordinary Gaussian translation give the exact FINITE
formula

    Phi_P[Q]=exp(g0t+b^2t^2/2)
       *sum_(n=0)^N_D t^n/n!
          E[Q(Z,W,Y+bt) R(Z,W,Y+bt)^n].             (R17)

Every summand is polynomial. The vacuum has only n=0 and is the
primitive exponential in (R2). Therefore the full actual image is
exp(psi)K[t], proving statement 3 implies statement 2. The
[nonconstant-phase original-image theorem](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#7-every-finite-actual-image-and-its-original-fixed-marker-image)
gives literal V_P=M_P whenever P is nonconstant; the constant-phase
exception is stated after (R2).

When b=0 the normalized period in (R17) has time degree at most N_D,
independently of the degree or Y support of R. In particular every safe
member has the original cutoff

    E[Q P^m]=0 for m>floor(deg Q/gamma),             (R18)

with bound zero if R=0. This is an upper bound, not a claimed exact
filtered image for arbitrary R. Purity forces a=0 and constant g,B
by statement 1 with p=1,psi=0, and then (R16) gives g0=b=0. Conversely
(R3) makes every original marked period polynomial by (R17). Hence
purity and safety are equivalent in the entire stated family.

For b!=0 and R!=0 put beta=deg_Y R. The normalized time degree in (R17)
is at most

    D+(1+beta)floor(D/gamma).                       (R19)

More precisely for a monomial Q=Z^iW^jY^k with positive R, terms can
occur only for 0<=n<=floor((j-i)/gamma), and their degree is at most
k+(1+beta)n. If j<i the period is zero. Reverse i,j for negative R.
For R=0 the normalized degree is at most deg_Y Q<=D, and it is zero
when b=0. These bounds retain ORIGINAL marker degrees.

At a safe phase, (R18) also gives moment r-jet budget
floor((1+d/gamma)r) for perturbations of original degree at most d:
the coefficient involving j perturbation factors has marker degree at
most dj and vanishes once m-j>floor(dj/gamma). This is a jet bound at
an already known safe phase. The
[cutoff-to-jet and certified-tail theorem](MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md#1-the-finite-jet-theorem)
then says that an actual coefficient nilpotence certificate of index nu reduces the remaining intervening moments to
floor((1+d/gamma)(nu-1)); no automatic ideal-prefix theorem is asserted.
For the zero phase use gamma=infinity in these formulas.

## 6. Exact safe recompositions

Let S_P={f in K[T]: f(P) is safe}. The
[general safe-recomposition theorem](GAUSSIAN_SAFE_RECOMPOSITION_CLASSIFICATION.md)
and the scalar criterion above give the exact ideal

    S_P=(T-g0) if a=0, g=g0 constant, B=0;
    S_P=(0) otherwise.                             (R20)

Indeed a safe nonconstant recomposition forces M_P finite. Its primitive
exponent in (R2) must then be linear, which is equivalent to b=0. The
unique safe shift is g0. Conversely P-g0=R is safe whenever the displayed
condition holds. For a constant phase, f(P)=0 is exactly f(g0)=0 and gives
the same ideal. This classifies safe polynomial recompositions even when
the recomposed phase has higher W or Y degree than (R1). No common bad
marker for all rejected recompositions is asserted.

## 7. Evidence and scope

The [focused exact checker](verify_gaussian_radial_one_sided_finite_images.py)
and [receipt](gaussian_radial_one_sided_finite_images_verification.json)
retain bounded radial Wick/Y-completion comparisons, relative spatial
and parameter identities under a nontrivial scalar gauge, explicit
nonzero endpoint terms, the cancelling leading-coefficient case, and
original marked one-sided controls with W^3 and Y^4 terms. Those checks
verify their recorded finite algebraic scope. They do not prove the
unrestricted relative comparison or the theorem by moment extrapolation.

The scalar pure-radial conclusion and the corresponding H21 one-sided
purity statement were already retained in the
[all-degree radial proof](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#all-degree-radial-obstruction-and-angular-support)
and [weakly one-sided support corollary](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#weakly-one-sided-angular-support-in-every-degree).
Here the scalar premise is strengthened to an arbitrary polynomial factor
and polynomial exponential gauge, and the conclusion classifies full
actual finite images with primitive vacuum. Formula (R17) and its bounds
keep every original marker. The strict remainder also permits arbitrary
higher W and Y degrees while its radial core keeps a constant quadratic
ordinary fibre. It makes no assertion about arbitrary cores
in K[S,Y] or phases with both positive and negative angular support.
The endpoint is explicitly retained, and the analytic comparison imports
keep their existing primary-source audits and pending independent-review
status. No independent replay, formal verification, external mathematical
review, publication or priority is asserted.
