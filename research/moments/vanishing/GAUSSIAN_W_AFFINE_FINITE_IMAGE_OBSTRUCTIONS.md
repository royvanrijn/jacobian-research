# Finite actual images and primitive vacua in W-affine Gaussian families

Internal written proof, 17 November 2026. This upgrades a unit-period
obstruction to a finite WHOLE actual-image obstruction and gives complete
finite-image criteria in the branches stated below. It uses the original
mixed-source map, not an alternate scalar-seed trace. Internal source
audits checked every scalar-shift face and the global closure. Independent
mathematical review, formal verification and external review are not
claimed; priority is unassessed.

The exact source and complete marked-curve comparison are retained in
the [H5 section of the affine-isotropic note](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#common-roots-at-which-h-vanishes-at-least-as-strongly-as-c),
in particular its [original-source construction](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#the-original-source-and-exact-localization),
[exact leading-face charts](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#exact-scales-and-all-possible-leading-faces)
and [complete component/neck/tail estimates](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#all-centres-punctures-and-crossover-estimates).
The constant mixed-source branch uses the precisely stated
[ordinary leading-coefficient obstruction](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#constant-mixed-source-purity-forces-a-quadratic-phase)
and its already declared imported proof dependencies. No stronger external
assurance is asserted here. The finite-image scalar normalization uses
the [general original-image theorem](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#7-every-finite-actual-image-and-its-original-fixed-marker-image).
The other branch inputs are the multibase note's
[finite-image theorem, Section 5](GAUSSIAN_MULTIBASE_QUADRATIC_FIBRE_PRIMITIVITY.md#5-every-finite-actual-image-has-primitive-vacuum)
and [separate mixed zero-root theorem, Section 9](GAUSSIAN_MULTIBASE_QUADRATIC_FIBRE_PRIMITIVITY.md#9-a-separate-mixed-source-extension-with-zero-root).

For the comparison argument below set K=C, O=K[t], and
M=O Phi_P[K[Z,W,Y]]. All images retain the original unit.

Let

    P=W(h(Z)+c(Z)Y)+aY^2+u(Z)Y+k(Z)

where a in C is CONSTANT and h,c,u,k are complex polynomials, with
c nonconstant and

    (h(0),c(0))!=(0,0).

Assume that at every common root alpha of h,c,

    ord_alpha h>=ord_alpha c.                        (F1)

Then the full actual Gaussian image M is NOT finite over O. There are
no degree or multiplicity bounds. In particular,
no nonconstant squarefree c has finite M for any h,a,u,k: outside the
zero-root case this follows here, while at a common origin its nonzero
c'(0) is rejected by the separate exact zero-root finite-image theorem.

## 1. Exact original source and the finite-image normalization

Set r=1/t,

    s=(rz-h(z))/c(z),
    phi=-s^2/2+t[a s^2+u(z)s+k(z)],
    eta=r dz/c(z).

The original circular Lagrange root is

    zeta=t[h(zeta)+c(zeta)Y].

If c(0)!=0, c(zeta) is an ordinary t-adic unit. If c(0)=0,h(0)!=0,
then zeta=t h(0)+O(t^2), and for m=ord_0 c,

    c(zeta)=c_m h(0)^m t^m(1+O(t)),

whose leading coefficient is nonzero and independent of Y. Thus every
required c-localized marker has a coefficientwise Laurent-t expansion
with polynomial Y coefficients. The retained integration-by-parts
identities make its ORIGINAL actual evaluation a horizontal functional
on the normalized rational-phase system with seed eta. Its seed value
is exactly Phi_P[1]. No rational marker with an undefined original
evaluation has been introduced.

Suppose M were finite. The retained all-dimension finite-image theorem
then gives

    Phi_P[1]=exp(psi(t))p(t),
    psi in tC[t], p in C[t], p(0)=1.

Subtract the spatially constant polynomial psi from phi and multiply
the actual functional by exp(-psi). The new seed value is p. Equivalently
use normalized volume eta/p with phase phi-psi; its seed value is one.
The scalar changes of connection are exact and preserve the original
map up to these explicit scalars. At infinity, dividing volume by p
only INCREASES every normalized volume order by deg p>=0. It remains
to prove positive order for eta at all finite-phase components of the
shifted phase phi-psi.

## 2. Four original parameter levels and one scalar shift

Expand

    phi=r^2 F2+r F1+F0+r^-1 Fm1,
    F2=-z^2/(2c^2),
    F1=zh/c^2+a z^2/c^2,
    F0=-h^2/(2c^2)-2a zh/c^2+zu/c,
    Fm1=a h^2/c^2-uh/c+k.                           (F2)

If deg psi<=1, the shift simply changes the regular polynomial k by
a constant, or does nothing. If Dpsi=deg psi>=2, the scalar shift adds
a leading monomial with parameter weight -Dpsi and spatial exponent
zero. Its lower parameter levels have strictly larger weights. In
every leading face, the minimum is at most -Dpsi<=-2, so any original
group of weight at least -1 is absent.

At a c-root alpha of multiplicity m, put x=(z-alpha)^-1 and v(x)=-N<0.
An exact final disc chart has x=xi+r^rho X, epsilon=rho+N>=0. The seed
order is

    v(eta)=1-(m-1)N+epsilon.                        (F3)

If m=1 or N<1/(m-1), it is positive already. Hence only m>=2 and
N>=1/(m-1) require a leading-face argument. As in the retained exact-
chart proof, a leading Laurent face f(tau), tau=r^N x, must have a
constant nonzero root of multiplicity mu, and

    epsilon>=-L/mu,                                (F4)

where L<0 is its leading parameter weight. The coefficient of X^mu
has exact valuation L+mu epsilon at the final unit Puiseux centre;
fixed derivatives preserve the strictly higher original weights.
The scalar shift contributes no additional derivatives but is included
in the leading Laurent polynomial before counting its root multiplicity.

## 3. Shared H5 roots

A shared root is nonzero under the source hypothesis. Write h=cv
locally, with v regular, and let m=ord_alpha c.

First suppose a=0. The r-level-1 group has pole degree at most m and
is strictly above the r-level-2 group in the bad range, since their
weight difference is at least mN-1>0. The level-zero group has actual
pole degree A<=m. The level-minus-one group k-uv is holomorphic, so
its weight is at least -1.

For deg psi<=1 the possible groups are levels 2,0,-1, with exponents
2m,A<=m,B<=0. A leading root has mu<=2; if its highest level is zero
it has only two monomials and mu=1, while a highest-level-minus-one
face is a single monomial and has no nonzero root. The bounds below
therefore cover this case as well. For deg psi>=2, the level-minus-one
group cannot lead, and the only
possible dangerous groups are

    level 2, exponent 2m;
    level 0, exponent A<=m;
    scalar level -Dpsi, exponent 0.

Thus every nonzero leading root has mu<=2. A face containing level 2
has L=2-2mN. Equations (F3)-(F4) give order at least

    -1+(m+1)N>0 if mu=1,   or N>0 if mu=2.

A face whose highest level is zero has only that group and the scalar
group, hence mu=1. Its A is positive and comparison with level 2 gives
(2m-A)N<=2. The resulting volume lower bound is

    1+(A-m+1)N.

It is at least one if A>=m-1; otherwise the endpoint bound gives
(A+2)/(2m-A)>0. A single scalar monomial has no nonzero root.
This exhausts all shifted faces and proves positivity.

Now suppose a!=0. F1 has exact pole degree 2m, with nonzero leading
coefficient a alpha^2/c_m^2. F2 has the same pole degree but parameter
level two, so is higher. F0 has pole degree at most m and Fm1 is
holomorphic. In the bad range they are strictly above F1:

    L0-L1>=-1+mN>0,
    Lm1-L1>=-2+2mN>0.

For deg psi<=1 the leading face is a single F1 monomial and cannot
be finite. For deg psi>=2 the only possible additional leading term
is the scalar shift, so mu=1. With L=1-2mN, (F3)-(F4) give the
strictly positive order (m+1)N. This covers every shared H5 root for
arbitrary a.

## 4. Unshared roots of arbitrary multiplicity

If h(alpha)!=0 and a=0, F0 has exact pole degree 2m, while F1,F2
have degree at most 2m but higher parameter levels. Fm1 has degree
at most m. In the bad range,

    Lm1-L0>=-1+mN>0.

Hence only the F0 monomial and, if present at the same weight, the
scalar shift can lead. For deg psi<=1 there is no finite bad-range
component; for deg psi>=2 a possible leading root has mu=1. Its
weight L=-2mN yields volume order at least 1+(m+1)N>0.

If a!=0, Fm1 instead has exact pole degree 2m with coefficient
a h(alpha)^2/c_m^2. It is strictly below all other original groups,
whose pole degrees are at most 2m and parameter levels are higher.
Again only a scalar shift can join it. Its possible root has mu=1
and L=-1-2mN, giving order at least 2+(m+1)N>0. The degree-one
scalar shift changes only a regular term and cannot cancel this pole.
These estimates also hold at alpha=0 when that root is unshared.

## 5. Complete comparison and the polynomial-seed contradiction

At ordinary finite centres eta has positive r order. At spatial infinity,
deg c>=1 gives order at least 1+(deg c-1)N>=1. Scalar phase shifts do
not change those volume estimates. Every higher-order finite horizontal
volume pole remains phase-irregular: its leading pole coefficient is
proportional to (r alpha-h(alpha))^2 times -1/2+a/r, nonzero in C(r)
outside the excluded shared origin. The infinity volume is logarithmic
or holomorphic.

The shifted phase remains spatially nonconstant. Its r^2 coefficient
-z^2/(2c^2) is nonconstant unless c is proportional to z. In that
exceptional shape the source hypothesis forces h(0)!=0, and the
r-level-1 coefficient has the nonconstant term h/(kappa^2 z).
A spatially constant shift cannot remove it.

Consequently the retained COMPLETE marked-curve comparison applies
to phi-psi, with all compact components, crossover annuli, descending
tails, regular punctures and gluing maps retained. All zero-exponential
eta seed periods tend to zero; eta/p has still more positive order.
The normalized actual source would provide a unit seed quotient for
eta/p. Equivalently, it would supply the nonzero polynomial p as a
zero-exponential eta period, which cannot tend to zero. This is the
required contradiction. No local stationary period has been substituted
for the original mixed source.

## 6. A companion constant-c primitive-vacuum result

If c=lambda!=0 is constant, the retained ordinary polynomial phase
has normalization H=(lambda t)^-1. Under finite M write its actual
seed exp(psi)p as above, shift phi by psi and normalize volume H/p.
For actual spatial degree D>=3, the top coefficient still has
v_infinity<=1, whereas v_infinity(H/p)=1+deg p. The retained ordinary
leading-coefficient criterion excludes the unit quotient. Thus the
original phase is quadratic, by the same successive coefficient
comparison as the retained constant-c proof: h,u are affine and k
has degree at most two.

Write a quadratic phase as X^T A_Q X+l^T X+k0, with A_Q symmetric and
Gaussian covariance Gamma. Its exact period is

    f=det(I-2t Gamma A_Q)^(-1/2)
       exp(tk0+t^2 l^T(I-2t Gamma A_Q)^(-1)Gamma l/2).

If f=exp(psi)p with polynomial psi,p, the rational-exponential
argument forces the rational exponent to equal psi, and then
det(I-2t Gamma A_Q)*p^2=1. Hence this determinant is one and p=1.
Conversely nilpotent Gamma A_Q makes the precision inverse, mean and covariance polynomial in t,
so the whole actual image is finite and its vacuum is primitive.
This proves finite-image vacuum primitivity for the entire constant
nonzero-c mixed branch, without a preassigned degree bound. More
explicitly, its finite-image criterion is that h,u are affine, k has
degree at most two, and Gamma A_Q is nilpotent.
Safety additionally requires its displayed polynomial exponent to vanish.

## 7. Complete finite-image criterion in the nonconstant-c H5 family

Drop the origin hypothesis, retain nonconstant c and the common-root
condition (F1), and keep a constant. Then

    M is finite over O
      iff a=0 and h,c belong to z^2 C[z].           (F5)

The polynomials u,k are unrestricted. Necessity follows because Sections
1-5 exclude finite M unless h(0)=c(0)=0; the exact mixed zero-root
theorem in the linked multibase note then forces a=h'(0)=c'(0)=0.
Conversely those conditions are sufficient by that zero-root theorem,
without needing any root-multiplicity condition. In every finite case,

    f=exp(k(0)t+u(0)^2t^2/2), M=fO.                 (F6)

Since c is nonconstant, P is nonconstant and the literal original image
V also equals fO. Thus the vacuum is primitive. Safety is equivalent
to u(0)=k(0)=0 in addition to (F5), and then P is pure. In particular
no nonconstant squarefree c has a finite actual image: (F5) would
require a double root at zero.

## 8. Complete finite-image criterion when c=0

There is no root-order condition on h in this separate branch. The
multibase finite-image theorem specialized to one pair and one ordinary
variable gives

    M is finite over O
      iff a=0 and either
        h is a nonzero constant, or h belongs to z^2 C[z].  (F7)

The polynomials u,k again are arbitrary. For clarity, the exact scalar
source root and Jacobians are

    r=t h(r),   J=1-t h'(r),   D=1-2at.

The multibase criterion says precisely that r is polynomial and J=D=1.
Thus a=0. If h(0)!=0, the root r is nonconstant (it begins h(0)t), and
J=1 gives h'(r)=0. Composition with a nonconstant polynomial is
injective, so h'=0: h is constant. If h(0)=0, the unique source root
is r=0, and J=1 is exactly h'(0)=0. These two alternatives give (F7),
and each makes r polynomial and J=D=1, proving sufficiency as well.

For h=h0!=0 the exact primitive unit and image are

    r=h0t,
    f=exp(t k(h0t)+t^2 u(h0t)^2/2), M=fO.           (F8)

For h in z^2 C[z] they instead are

    r=0, f=exp(k(0)t+u(0)^2t^2/2), M=fO.            (F9)

The literal original image V equals M for every nonconstant P; for a
constant P=k0 the exact exception is V=C exp(k0t). These finite-image
criteria and formulas also hold over every characteristic-zero field,
by the purely algebraic multibase theorem.

Safety is exactly vanishing of the exponent in (F8) or (F9). Explicitly,
for h=h0!=0 it is k(z)=-z u(z)^2/(2h0); for h in z^2 C[z] it is
u(0)=k(0)=0. All these safe phases are pure. This recovers the known
shear and positive-weight forms while proving the stronger primitive
finite-image statement.

## 9. Constant quadratic fibre and assurance boundaries

The hypothesis that a is constant is essential in Section 6. For any
polynomial f(Z), the Gaussian-preserving shear from the
[exact automorphism-kernel theorem](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#8-exact-original-kernels-for-gaussian-automorphism-orbits)
sends the safe phase YW to

    (Y-Zf)(W+fY-Zf^2/2)
      =W(Y-Zf)+fY^2-(3/2)Zf^2Y+(1/2)Z^2f^3.

Its original unit is one and every original marker is safe by exact
equivariance with the inverse shear. For nonconstant f this is a
nonquadratic safe phase with constant mixed coefficient c=1 but
variable quadratic fibre A(Z)=f(Z); f=Z gives degree five. Thus the
constant-c quadratic conclusion cannot be extended by silently allowing
variable a.

The unrestricted proofs above are written source, valuation and
finite-image arguments. No finite checker is needed or claimed to
certify their universal range. Their imported analytic comparison and
leading-coefficient obligations remain exactly those linked at the
start. These families have primitive finite vacua and satisfy safety
implies purity; general field safety implies purity remains open.
