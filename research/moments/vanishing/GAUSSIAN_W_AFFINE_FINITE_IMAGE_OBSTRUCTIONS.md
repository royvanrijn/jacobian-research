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

**Complete finite-image conclusion.** For CONSTANT a and arbitrary
NONCONSTANT c, with NO common-root or multiplicity restriction,

    M finite over O iff a=0 and h,c belong to z^2 C[z].

In every finite case f=exp(k(0)t+u(0)^2t^2/2) is primitive and V=M=fO.
Safety additionally requires u(0)=k(0)=0, and every safe member is pure.
Section 14 proves this universal statement using a fixed ORIGINAL marker
that clears the volume poles. The earlier face proofs are retained for
their stronger SCALAR-only conclusions on the stated subclasses. The
scalar criterion of Section 12 is NOT automatically extended to arbitrary
root orders by the new finite-image argument.

First retain the H5 scalar-source obstruction used in Sections 1-5.

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

## 10. Simple and double common roots in the finite-image criterion

Keep a constant and let c be NONCONSTANT. At each nonzero common root
alpha of c,h, assume either

    ord_alpha h>=ord_alpha c,

or

    ord_alpha c=2, ord_alpha h=1.                   (F10)

Unshared roots have no multiplicity restriction, and no root-order
condition at the common origin is needed: that case has its separate
exact zero-root formula. Then the complete criterion is

    M finite over K[t]
      iff a=0 and h,c belong to z^2 K[z].            (F11)

The polynomials u,k are arbitrary. Every finite image has primitive unit
f=exp(k(0)t+u(0)^2t^2/2), and since P is nonconstant its literal original
image is V=M=fK[t]. Safety additionally requires u(0)=k(0)=0 and is then
equivalent to purity. In particular (F11) classifies every nonconstant c
whose roots all have multiplicity at most two. Constant c is explicitly
excluded here and has the separate quadratic branches of Section 6.

Only the double/simple-order alternative in (F10) needs a new local
estimate. Work outside h(0)=c(0)=0, so this shared alpha is nonzero.
Assume finite M, normalize its actual unit exp(psi)p as in Section 1,
shift the phase by psi and use the seed eta/p. Put r=1/t,
x=(z-alpha)^(-1), v(x)=-N<0. In an exact final chart the extra scale
epsilon>=0 gives

    v(eta)=1-N+epsilon.

For N<1 this is positive. At m=2,ell=1, the pole bounds for the four
original groups (F2) are

    a=0:  j2=4, j1=3, j0<=2, jm1<=1;
    a!=0: j2=4, j1=4, j0<=3, jm1<=2.

The fourth-power leading coefficient is nonzero: for a=0 it occurs in
F2, and for a!=0 in F1 with coefficient a alpha^2/c_2^2.

At N=1 the minimum original weight is at most -2 if a=0 and at most
-3 if a!=0. The scalar shift can add only one leading exponent-zero
monomial. The nonzero leading face therefore has at most five distinct
monomials; a nonzero root has finite multiplicity mu<=4. Finite phase
requires such a root and the exact-chart estimate epsilon>=-L/mu>0.
Thus v(eta)=epsilon>0. If the scalar term has strictly smaller weight,
it is the sole leading monomial and admits no finite-phase component.

For N>1 and a=0, level two uniquely dominates the other ORIGINAL groups:

    L1-L2>=N-1>0,
    L0-L2>=2N-2>0,
    Lm1-L2>=3N-3>0,    L2=2-4N.

Only the scalar shift can join it. A tied binomial has a simple nonzero
root, so mu=1 and the seed order is at least

    1-N-L2=3N-1>0.

For N>1 and a!=0, level one instead uniquely dominates:

    L2-L1=1,
    L0-L1>=N-1>0,
    Lm1-L1>=2N-2>0,    L1=1-4N.

Again only the scalar term can join. Its possible root is simple, giving
seed order at least 1-N-L1=3N>0. These alternatives include every final
cancellation refinement; dividing by p only increases the infinity
volume order.

All remaining source, ordinary-centre, infinity, spatial-puncture,
crossover-neck and descending-tail arguments are precisely Sections 1-5.
Thus no finite image exists outside the common-zero-origin case. The
same estimates exclude a pure unit directly by taking psi=0,p=1,
without a finite-image premise. Inside the common-zero-origin case,
Section 9 of the multibase note proves exactly a=h'(0)=c'(0)=0
and the stated primitive exponential image. This proves (F11), including
all original markers in its sufficient direction.

The necessity proof above is complex analytic through the retained
comparison inputs. It descends to every characteristic-zero field K:
under a hypothetical finite image, adjoin the finite phase,p,psi
coefficients to Q and embed that finitely generated field in C. The
entire actual scalar identity specializes, while polynomial degrees,
geometric common-root orders and the nonzero source leading coefficients
are preserved under field extension. The complex contradiction applies.
Sufficiency is the algebraic all-field zero-root formula.

The older sharp PURE-unit classification through common-root
multiplicity eight does not automatically extend to finite images.
The scalar polynomial psi adds a new constant face; already at
(m,ell)=(4,2), leading weight data can permit a quartic degeneracy.
No finite-image example or further multiplicity theorem is claimed by
that observation. The finite-image assertion here has exactly (F10)'s
scope and the complete comparison dependencies of Sections 1-5.
Section 14 resolves the finite actual-image question for all root orders
by changing to a fixed original marker. The quartic warning remains
valid for this unit-volume/scalar-only face method.


## 11. Triple common roots

Keep c NONCONSTANT. The composite criterion extends further: at each
nonzero common root allow either ord h>=ord c, or c-multiplicity at
most THREE. Unshared roots still have arbitrary multiplicity, and the
common origin still uses its separate exact source formula. Then

    M finite over K[t]
      iff a=0 and h,c belong to z^2 K[z],            (F12)

with arbitrary u,k, primitive unit exp(k(0)t+u(0)^2t^2/2), and literal
V=M. Safety additionally requires u(0)=k(0)=0 and is equivalent to purity.
In particular this classifies every nonconstant c with root multiplicities
at most three. The source and field-descent arguments are those of
Section 10; only the new local triple-root estimate follows.

Outside h(0)=c(0)=0, let alpha be a common root with ord c=3 and
ell=ord h in {1,2}. Thus alpha!=0. Write c=C(z-alpha)^3+... and
h=H(z-alpha)^ell+..., with C,H nonzero. After subtracting the scalar
psi as before, put x=(z-alpha)^(-1), v(x)=-N<0, and use the exact final
extra scale epsilon. The normalized volume has order 1-2N+epsilon.
For N<1/2 it is positive already. At N=1/2 the original F2 weight is
-1, so the minimum L is negative. Any nonzero leading face has finite
multiplicity (at most four after the scalar term); finite phase forces
epsilon>=-L/mu>0. Thus the boundary N=1/2 is positive as well.

For ell=2, N>1/2 and a=0, F1 has pole order four and lies strictly
above F2, by L1-L2=-1+2N>0. If u(alpha)!=0, F0 has exact pole order
three and Fm1 order one; Lm1-L0=-1+2N>0. If u(alpha)=0, Fm1 is
holomorphic, hence above the F2 weight 2-6N<-1. Therefore only F2,
F0 and the scalar shift can lead, giving mu<=2. A face containing F2
has L=2-6N and volume bound -1+4N>0 for mu=1, or N>0 for mu=2.
If F0 is highest, its actual pole exponent A<=3 obeys (6-A)N<=2.
Since N>1/2, A=3. Its possible scalar binomial has mu=1 and volume
at least 1+N>0. A single monomial has no nonzero root.

For ell=2 and a!=0, F1 has exact pole order six, F0 order four and
Fm1 order two. Above F1 their weight differences are respectively
1, -1+2N and -2+4N, all positive. Only a scalar can join the leading
monomial; its possible root is simple and volume is at least 4N>0.

For ell=1 and 1/2<N<1, if a=0 then F0 has exact pole order four.
The weight differences above it are 1-N, 2-2N and at least -1+2N,
all positive. Its scalar binomial has mu=1 and volume at least 1+2N.
If a!=0, Fm1 has exact pole order four and uniquely dominates the
other original groups: their differences are 1-N, 2-2N and 3-2N.
Its possible scalar binomial gives volume at least 2+2N>0.

For ell=1 and N>1, a=0 makes F2 unique among originals; the differences
are N-1, 2N-2 and at least 4N-3. A scalar binomial has mu=1 and volume
at least 4N-1>0. If a!=0, F1 is unique, with differences 1, N-1 and
2N-2; its possible scalar binomial gives volume at least 4N>0.

It remains to check ell=1,N=1, where several original levels coincide.
Put tau=rx and B(tau)=tau^2(alpha tau-H). The exact leading faces are

    -B(tau)^2/(2C^2), L=-4, if a=0;
     a B(tau)^2/C^2,  L=-5, if a!=0.               (F13)

All u,k terms, and all other original terms, are strictly higher weight.
The scalar shift can add only a constant at this weight; if it has smaller
weight it is alone and admits no finite-phase component.

Every nonzero root of B^2 plus ANY constant has multiplicity at most two.
For zero constant, the only nonzero B root H/alpha is simple. For nonzero
constant a multiple root requires B'=tau(3alpha tau-2H)=0. Its nonzero
critical point 2H/(3alpha) is simple, with

    B=-4H^3/(27alpha^2)!=0, B''=2H!=0.

Hence mu<=2 in both cases of (F13). The base volume is -1, so epsilon
is at least two when a=0 and at least 5/2 when a!=0. The resulting
orders are at least one and 3/2. This covers the final critical scale,
including every coefficient cancellation and scalar-shift degree.

Thus every finite-phase triple-root component has positive order.
The unchanged complete source, puncture, neck and descending-tail proof
excludes a finite image outside the common-zero-origin case. Taking
psi=0,p=1 also excludes a pure unit directly. The exact zero-root result
then gives (F12) and its sufficient all-original-marker direction.
No finite-image extension of the older multiplicity-eight PURE-unit
result is inferred; the assertion here has exactly the stated composite
H5-or-multiplicity-at-most-three scope.

## 12. A scalar criterion on the covered mixed-source branches

There is a stronger consequence of the proofs. Keep a constant and
restrict either to nonzero constant c, or to the NONCONSTANT c family
of Section 11. For the ACTUAL Gaussian unit f, the following are equivalent:

    f=p exp(psi) for some p in K[t], psi in tK[t];
    the full actual image M is finite over K[t].

Whenever these hold, p=1 and the original vacuum generates M. Consequently
polynomiality of f, purity and all-original-marker safety are equivalent
in these covered branches. The statement does not include c=0.

To check the extra direction, inspect the source arguments. Finiteness
was used only to obtain f=p exp(psi), with p(0)=1. The actual c-localized
functional exists for every phase in the nonzero-root case. After the
phase shift it takes eta/p to one, so its restriction to that seed's
cyclic differential module is the same rational unit quotient used in
Sections 1-5, 10 and 11. All volume and complete-comparison estimates
therefore apply directly to this scalar identity. In the zero-root case
the scalar proof of Section 9 in the multibase note likewise starts from
f=p exp(psi), forces a=h'(0)=c'(0)=0 and p=1, and then controls every
original marker. For nonzero constant c, Section 6 uses the scalar identity
to force a quadratic phase with nilpotent precision, again giving p=1
and all-marker finiteness. The converse is the general finite-image theorem.

These implications also descend to characteristic-zero coefficient fields
by the finite-coefficient embedding argument of Section 10: the assumed
scalar identity alone is preserved, and the sufficient Gaussian formulas
are algebraic. Finally, a polynomial exp(psi) forces psi=0, proving the
stated polynomial-unit, purity and safety equivalence. No claim about
arbitrary Gaussian phases follows from this scalar criterion.


## 13. Common roots strictly below half the c-order

Keep a constant and c NONCONSTANT. At every nonzero common root alpha
of c,h, allow any of the following conditions:

    ord_alpha h>=ord_alpha c;
    ord_alpha c<=3;
    1<=ord_alpha h<(ord_alpha c)/2.                  (F14)

Unshared roots have arbitrary multiplicity, and the common origin still
uses its separate exact source formula. Then the complete criterion is

    M finite over K[t]
      iff a=0 and h,c belong to z^2 K[z],            (F15)

with arbitrary u,k. In every finite case the actual unit is primitive,

    f=exp(k(0)t+u(0)^2t^2/2), V=M=fK[t].

The scalar criterion of Section 12 extends to exactly this composite
class: an exponential-polynomial actual unit exists if and only if M
is finite, and polynomiality of f, purity and safety are equivalent.
The new alternative is STRICT: equality ord h=(ord c)/2 is not covered
by the LOCAL estimate below. The (m,ell)=(4,2) warning in Section 10
remains for this scalar-only method; the independent full-image argument
of Section 14 now covers finite images at that equality and all orders.

Only the third alternative in (F14) needs proof. Outside the common-zero
origin case, fix a shared alpha!=0 and write

    m=ord_alpha c, ell=ord_alpha h, d=m-ell>ell>=1.

Assume an actual scalar identity f=p exp(psi), p(0)=1, and use the
unconditional source functional, shifted phase phi-psi and seed eta/p
from Sections 1-5 and 12. Put r=1/t, x=(z-alpha)^(-1), v(x)=-N<0.
In an exact final chart the extra relative scale epsilon>=0 gives

    v(eta)=1-(m-1)N+epsilon,
    epsilon>=-L/mu,                                 (F16)

where L<0 is the leading phase weight and mu is the multiplicity of
its selected nonzero root in tau=r^N x. For N<1/(m-1) the volume is
positive already. A sole leading monomial has no nonzero root; k
leading monomials give mu<=k-1 by the retained Vandermonde argument.
Division by the polynomial p only increases the infinity order.
Absorb the degree-at-most-one part of psi into k. Any remaining scalar
leader has INTEGER degree D>=2 and contributes just one exponent-zero
monomial of weight -D. The four original groups are precisely (F2).

Suppose first a=0. Their exact pole orders at alpha are

    j2=2m, j1=2m-ell, j0=2d, jm1=B<=d.

Here j0 is exact: its square term has order 2d>m, so zu/c cannot
cancel it. If Fm1 vanishes omit it; otherwise B is its ACTUAL highest
Laurent exponent after absorbing the linear shift, and can be negative.
For N<1/ell, F1 and F2 lie strictly above F0 by 1-ell N and 2-2ell N.
Only F0,Fm1 and the scalar can lead. If F0 belongs to a binomial face,
mu=1, L=-2dN, and (F16) gives

    v(eta)>=1+(m-2ell+1)N>0.

If Fm1 and the scalar tie with F0 higher, then mu=1 and
(2d-B)N<=1. Since m-1-B>=ell-1>=0, the volume obeys

    v(eta)>=2-(m-1-B)N
            >=(3m-4ell-B+1)/(2d-B)>0.               (F17)

The numerator is at least 2m-3ell+1>0. If all three tie, then
N=1/(2d-B) and D=2d/(2d-B)<=2. The integer condition D>=2 forces
D=2, B=d, N=1/d. Now mu<=2 and the volume is at least
2-(m-1)/d=(m-2ell+1)/d>0. These are all possible nontrivial faces.

For a=0 and N>1/ell, F2 uniquely dominates the other original groups.
The differences above it are ell N-1, 2ell N-2, and

    -3+(2m-B)N>=-3+(m+ell)N>m/ell-2>0.

Only the scalar can join. Its possible binomial root is simple, giving
v(eta)>=-1+(m+1)N>0.

At N=1/ell, write c=C(z-alpha)^m+... and h=H(z-alpha)^ell+...,
where C,H are nonzero. Define

    B0(tau)=tau^d(alpha tau^ell-H).

Because d>ell, all u,k terms are strictly higher than the square terms.
For a=0 the exact leading face is -B0^2/(2C^2), of weight L=-2d/ell.
The scalar may add a constant; if it has smaller weight it is the sole
leader and admits no finite-phase component. Every nonzero root of
B0^2 plus any constant has multiplicity at most two. Indeed B0 has
simple nonzero roots, and

    B0'=tau^(d-1)(m alpha tau^ell-dH).

Every nonzero critical point is simple, and there
B0=-(ell/m)H tau^d!=0. Thus no constant shift creates multiplicity
above two. Equation (F16) gives epsilon>=d/ell and

    v(eta)>=1-(m-1)/ell+d/ell=1/ell>0.              (F18)

Finally let a!=0. The exact pole orders are

    j2=2m, j1=2m, j0=2m-ell, jm1=2d.

Their a-dependent leaders cannot cancel against lower-order u,k terms.
For N<1/ell, Fm1 uniquely dominates the other original groups, with
differences 1-ell N, 2-2ell N and 3-2ell N. A possible scalar binomial
has mu=1 and volume at least 2+(m-2ell+1)N>0. For N>1/ell, F1
uniquely dominates, with differences 1, ell N-1 and 2ell N-2; its
possible scalar binomial gives volume at least (m+1)N>0. At resonance
the leading face is a B0^2/C^2 of weight -1-2d/ell. The same shifted-
square argument gives mu<=2 and

    v(eta)>=1/2+1/ell>0.                            (F19)

These alternatives cover all leading faces and every subsequent exact
cancellation refinement. All orders are strict before division by p.
The unchanged complete source, unshared-root, ordinary-centre, infinity,
spatial-puncture, crossover-neck and descending-tail arguments of
Sections 1-5 therefore exclude the assumed scalar identity outside the
common-zero-origin case. The exact zero-root result then gives (F15)
and its all-original-marker sufficient direction. No all-marker
finiteness assumption was used after the initial scalar normalization,
so Section 12's stronger scalar criterion extends as stated.

The same finite-coefficient field embedding as in Section 10 gives the
necessity over every characteristic-zero field K; the source leading
coefficients and geometric root orders survive extension. Sufficiency
is the algebraic zero-root formula. This proves the composite class
(F14), without adding the equality ell=m/2 or any unrestricted
finite-image consequence of the older multiplicity-eight pure theorem.

## 14. A fixed original logarithmic marker removes every root-order restriction

Let c be ANY nonconstant polynomial and keep a constant. Then, without
any common-root, multiplicity or degree bound on h,c,u,k,

    M finite over K[t]
      iff a=0 and h,c belong to z^2 K[z].            (F20)

This holds over every characteristic-zero field K. In every finite case,

    psi=k(0)t+u(0)^2t^2/2,
    f=exp(psi), V=M=exp(psi)K[t].                   (F21)

The phase is nonconstant because c!=0, so V=M has no constant-phase
exception here. The original vacuum is primitive. Safety is equivalent
to (F20) together with u(0)=k(0)=0, and consequently every safe member
is pure. This is a FULL-image theorem; it does not assert that an
exponential-polynomial UNIT alone forces these conditions for arbitrary
root orders.

### The exact logarithmic marker and its nonzero original period

First work over C and outside the common-zero-origin case:

    (h(0),c(0))!=(0,0).

Take rho=gcd(c,c') with monic gcd convention, and put c_red=c/rho.
Thus c_red is squarefree, has exactly the roots of c, and has positive
degree. Equivalently rho=c/c_red is a fixed nonzero ORIGINAL polynomial
marker. The normalized differential becomes

    rho*eta=r dz/c_red,   r=1/t.                    (F22)

The differential omega=dz/c_red is fixed and logarithmic on P1: it has
at most simple poles at the roots of c, and its infinity pole is simple
when deg c_red=1 and absent when deg c_red>=2.

Section 1's original c-localized source construction needs only the
displayed origin condition, not H5. For clarity, its exact pullback is

    zeta=t[h(zeta)+c(zeta)Y],
    J=1-t[h'(zeta)+c'(zeta)Y],
    zeta_Y=t c(zeta)/J,
    eta -> dY/J,
    phi(t,zeta)=-Y^2/2+t[aY^2+u(zeta)Y+k(zeta)].

Hence the actual normalized period of rho*eta is precisely the ORIGINAL
marked period

    Phi_P[rho]=E_Y[rho(zeta)
        exp(t[aY^2+u(zeta)Y+k(zeta)])/J].            (F23)

This period is NONZERO. If c(0)!=0, then rho(0)!=0 and its value at
t=0 is rho(0). If c(0)=0 and h(0)!=0, let m0=ord_0 c>=1 and
s=ord_0 rho=m0-1. Write rho(z)=rho_s z^s+O(z^(s+1)), rho_s!=0.
Since zeta=t h(0)+O(t^2), its leading coefficient is independent of Y,
and (F23) gives the exact leading term

    Phi_P[rho]=rho_s h(0)^s t^s+O(t^(s+1))!=0.      (F24)

Both J^-1 and the fibre exponential have constant coefficient one.
This proof covers s=0 as well. The nonzero period may start in positive
t order; no nonzero initial mass is required.

### Full-image finiteness supplies the needed marked polynomial

Suppose M were finite. The general finite-image theorem gives its
normalized generator exp(psi), psi in tC[t], without assuming primitive
vacuum. Since rho is an ORIGINAL marker, there is a polynomial

    p_rho=exp(-psi)Phi_P[rho] in C[t],   p_rho!=0.    (F25)

This is the decisive use of the WHOLE finite-image hypothesis. A scalar
identity for Phi_P[1] alone does not imply (F25).

Shift the rational phase to phi-psi. Multiplying the original horizontal
functional by exp(-psi) is the corresponding exact scalar connection
gauge. Choose the distinguished seed to be

    H(t)omega,   H(t)=r/p_rho(t).

By (F22)-(F25) its actual seed value is one. Its scalar normalization
has strictly positive infinity valuation

    v_infinity(H)=1+deg p_rho>0.                    (F26)

Even if p_rho vanishes at zero, 1/p_rho is a legitimate element of
C(t), with a Laurent series having finite principal part at zero. The
source identities
allow these Laurent coefficients, and the seed value is still exactly
one. No arbitrary rational spatial marker has been evaluated outside
the original source domain.

### Direct application of the retained complete logarithmic seed lemma

Apply the [fixed-puncture logarithmic seed bound](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#fixed-puncture-logarithmic-seed-bound)
with

    S={roots of c} union {infinity},
    U=P1 minus S, phase phi-psi, fixed form omega=dz/c_red.

Its statement allows EVERY rational phase with fixed spatial pole set
S, including an arbitrary spatially constant polynomial shift psi.
Its proof retains the complete finite-phase components, actual annulus
complexes and gluing, descending tails and regular punctures. Every
zero-exponential omega seed period is O(1+log|t|); multiplication by
H with (F26) makes it tend to zero. The same lemma explicitly excludes
a unit horizontal formal seed functional supplied at another parameter
point, using its ordinary-parameter seed-extension and full period
pairing. Thus the original formal unit just constructed is impossible.

This invocation has exactly that lemma's existing comparison dependencies
and pending independent-review status. No local quartic estimate or
additional stationary-phase assertion is needed. All roots of c are
already fixed punctures, and the logarithmic differential remains
logarithmic after any marked-model refinement. The scalar psi may change
finite balls and their action values but cannot invalidate the lemma.

We have proved NO finite actual image outside h(0)=c(0)=0 for ANY
nonconstant c. More precisely, the same argument shows that the nonzero
original period Phi_P[rho] cannot equal exp(psi)q for ANY polynomial
psi and nonzero polynomial q. This identifies a fixed original marked
obstruction, rather than replacing the Gaussian unit by an untracked
localized period.

### Complete criterion, the (4,2) case and field descent

Inside h(0)=c(0)=0, the [exact mixed zero-root theorem](GAUSSIAN_MULTIBASE_QUADRATIC_FIBRE_PRIMITIVITY.md#9-a-separate-mixed-source-extension-with-zero-root)
proves precisely a=h'(0)=c'(0)=0 for finite M. It supplies every original
marker in the sufficient direction and gives (F21). Combined with the
preceding obstruction this proves (F20), primitive vacuum, and the safety
criterion. For example a nonzero shared root with
ord_alpha c=4, ord_alpha h=2 requires no extra hypothesis: if the origin
has (h(0),c(0))!=(0,0), there is no finite image; at a common origin,
finiteness additionally requires a=h'(0)=c'(0)=0 as in (F20),
independently of all other roots. At the quartically degenerate (4,2)
component N=epsilon=1/2, rho has local order three and raises the
original zero volume order by 3N=3/2, exactly the logarithmic order
1+epsilon.
The zero-volume quartic face discussed in Section 10 is therefore no
obstruction to this FULL-image classification.

For example h=z^2(z-1)^2, c=z^2(z-1)^4 and a=u=k=0 has a nonzero
shared root with orders (4,2), while its common origin satisfies (F20).
It is pure and safe by strict positive angular support. Here
rho=z(z-1)^3 and Phi_P[rho]=0 identically by the same weight argument;
the nonzero marked-period premise of the logarithmic obstruction is
absent, exactly as required.

The necessity descends to every characteristic-zero field. Under finite
M its generator exp(psi) is defined over that field. Adjoin to Q the
finitely many phase, rho, c_red, psi and p_rho coefficients and embed
that finitely generated field in C. The exact original marked identity
(F25) and the nonzero leading term (F24), or rho(0)!=0, are preserved.
The squarefree/gcd identities and the origin condition survive the
embedding, so the complex logarithmic contradiction applies. The
zero-root sufficient direction is already algebraic over every field.

For the broader constant-c branches, Sections 6 and 8 remain the complete
finite criteria. Together with (F20), every finite actual image in this
entire constant-a W-affine family has a primitive vacuum; hence safety
implies purity throughout this architecture. General field safety implies
purity is not proved outside it.

Finally, Section 12's stronger SCALAR-only criterion remains restricted
to its stated classes. For an unclassified scalar exponential-polynomial
unit outside those classes, (F25) is unavailable; the logarithmic marker
would instead demonstrate failure of the common finite lattice. No new
purity classification or scalar-only exclusion is claimed here.

The [focused checker](verify_gaussian_w_affine_logarithmic_marker.py) and
[receipt](gaussian_w_affine_logarithmic_marker_verification.json) replay
finite exact radical/volume identities and original leading-marker
controls, including mixed root orders (4,2) and c(0)=0,h(0)!=0. They do
NOT certify unrestricted multiplicity or the imported comparison theorem.
The written proof above closes the root-order restriction.
