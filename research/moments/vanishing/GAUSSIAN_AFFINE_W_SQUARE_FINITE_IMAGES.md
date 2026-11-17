# Affine W-square coefficient: scalar rigidity and finite actual images

All periods are ORIGINAL unweighted Gaussian Wick periods. The argument
keeps a possible polynomial vacuum factor and the full scalar exponential
gauge. Its exact scalar criterion also classifies the full actual image.

Let lambda!=0, mu,c,a be constants over a characteristic-zero field K, and
let h,u,k be arbitrary one-variable polynomials. Use the original Witt
Gaussian pair Z,W, with E[ZW]=1, and independent ordinary Gaussian Y,
with E[Y^2]=1; all other basic pairings vanish. Write

    Phi_P[Q]=E[Q exp(tP)], V_P=Phi_P[K[Z,W,Y]],
    M_P=K[t]V_P.

Put

    P=(lambda Z+mu)W^2+W(h(Z)+cY)+aY^2+u(Z)Y+k(Z).

Then the following are equivalent:

1. The full original image M_P is finite over K[t].
2. The actual vacuum is p(t)exp(psi(t)), with p polynomial, p(0)=1,
   and psi in tK[t].
3. a=0 and h,u,k are constant.

In every such case p=1, psi=k0 t+u0^2 t^2/2, and
V_P=M_P=exp(psi)K[t]. Safety is exactly u0=k0=0; every safe member is
pure. This includes c=0. The proof uses the retained complete exponential
curve comparison with its existing review limitations.

## 1. Exact curve and actual scalar source

Write ell=lambda(z-alpha), alpha=-mu/lambda, and

    d=1-2at, A=t ell+t^2 c^2/(2d),
    B=t h+t^2 c u/d, Psi0=t k+t^2 u^2/(2d),
    q=c^2+2d ell/t,
    phi=Psi0-(z-B)^2/(4A), eta=dz/(tU), U^2=q.

The following algebraic expression fixes all signs and normalizations:

    phi=N/q,
    N=t k q-d(z-th)^2/(2t^2)+c u(z-th)+t ell u^2.       (1)

In powers of t, the numerator is

    N=-z^2/(2t^2)+(a z^2+z h)/t
       +[2ell k-2a z h-h^2/2+cuz]
       +t[(c^2-4aell)k+a h^2-cuh+ell u^2].             (2)

Use U as the rational coordinate and put

    beta=alpha-tc^2/(2lambda d), s=t/(2lambda d),
    z=beta+sU^2, eta=dU/(lambda d).                    (3)

Consequently phi is an even Laurent polynomial in U, with only the
possible negative power U^-2. Its polynomial-part degree D is even.
If M=deg_z N, then D=2(M-1) whenever M>=1.

There is an actual horizontal scalar functional on this full punctured
U-line, with the exact original vacuum as its eta-value. This assertion
has no lower degree requirement on ell. For c!=0 put

    g=gcd(ell,c^2), e=ell/g, C=c^2/g,
    Rmov=C[z(4ae-C)-2eh]+4e^2cu.

The exact source cases are as follows:

* If ell(0)!=0, use the
  [completed-z original-source construction](GAUSSIAN_UNIT_ELL_FORMAL_SOURCE_FINITE_IMAGE_EXCLUSION.md#1-completed-z-wick-contraction-preserves-the-full-original-source).
* If c!=0 and ell(0)=0, g=gcd(ell,c^2)=1. The moving determinant
  restriction is acyclic unless Rmov=0, by the retained
  [exact source criterion and localization](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#exact-source-criterion-and-localization-of-the-original-seed).
  In that exceptional case e=ell=lambda z and the
  [moving-common-root bridge](GAUSSIAN_MOVING_COMMON_ROOT_SCALAR_EXCLUSION.md#5-algebraic-verification-on-the-residual-rmov0-branch)
  applies with gamma=0, nu=1.
* If c=0, the determinant's sole spatial divisor is z=alpha. Its
  W-source th(alpha)-alpha is nonzero unless alpha=h(0)=0. Clean
  restriction gives the full original-source localization. In the
  remaining simple common origin, the
  [fixed-origin scalar bridge](GAUSSIAN_SIMPLE_ORIGIN_ALTERNATE_SCALAR_EXCLUSION.md)
  gives the actual vacuum and every covariant iterate.

The exact Gaussian reduction, normalized determinant factor and tracked
seed are given in the
[unit-ell elimination](GAUSSIAN_UNIT_ELL_FORMAL_SOURCE_FINITE_IMAGE_EXCLUSION.md#2-gaussian-fibre-elimination-and-the-tracked-cover-seed)
and the retained
[normalized-cover construction](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#normalized-cover-phase-and-horizontal-poles).
These are source assertions from the unit-ell, moving-root,
fixed-simple-origin and clean localization proofs. Their constructions
apply at degree one; only their old volume conclusions had a degree
restriction. No original W-marker identification is needed for the
alternate scalar bridges.

The
[general finite actual-image theorem](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#7-every-finite-actual-image-and-its-original-fixed-marker-image)
gives statement 2 from statement 1. No primitive vacuum is assumed here.
Assume statement 2. Gauge the curve phase by -psi and divide eta by p.
Its genuine horizontal cyclic functional has value one. The comparison
contradiction below needs only this scalar premise. Division by p adds
nonnegative infinity order and does not change spatial poles. It is
therefore enough to prove strict positivity for the unnormalized eta.

## 2. Excluding polynomial-part degree at least four when a!=0

Suppose D=2K>=4, so M=K+1>=3. At t=infinity, beta and s are bounded,
and s is a nonzero unit. The scalar in eta has order one. The highest
U coefficient equals s^M times the highest z coefficient of N.
By (2), that latter nonzero coefficient is a Laurent polynomial in t
with exponents among -1,0,1: the sole t^-2 term has z-degree two and
cannot occur at degree M>=3. Therefore its infinity order vtop<=1.

The U^-2 coefficient is b=N(t,beta). Since beta is bounded and the
coefficients of (2) have infinity order at least -1, v(b)>=-1.

At finite U-centres the volume has strictly positive order. At infinity
write an exact final disc chart

    U=xi(r)+r^(-N0+epsilon)X,
    v(xi)=-N0<0, epsilon>=0, r=1/t.

Its volume order is 1-N0+epsilon. If N0<1 it is already positive.
For N0>=1 the pole term b/U^2 is integral, including after every further
zoom. Meanwhile the top polynomial coefficient has order at most one
and D>=4. The coefficient of X^D in a finite-phase chart therefore
forces

    epsilon>=N0-vtop/D,
    v(eta)>=1-vtop/D>0.                              (4)

The pole term cannot cancel a negative-order top polynomial Taylor
coefficient: its corresponding coefficient is integral. Lower
polynomial degrees have zero D-th derivative. This is the retained
exact-final-chart argument, allowing arbitrary lower coefficients and
any scalar phase shift.

## 3. Excluding polynomial-part degree at least four when a=0

Now beta=alpha-tc^2/(2lambda) is polynomial, s=t/(2lambda), and eta is
constant dU/lambda. Write

    phi=b(t)U^-2+a0(t)+sum_(j=1)^K a_(2j)(t)U^(2j), K>=2.

For EVERY j>=2, Taylor expansion of (1) gives

    a_(2j)=s^(j+1)/(j+1)! * partial_z^(j+1)N(t,beta). (5)

The derivative order is at least three and kills the sole t^-2 term
-z^2/2. Every remaining coefficient of N has t-exponent at least -1.
Since beta is polynomial, (5) is a polynomial divisible by t^j.
If it is nonzero, its pole order B_j at time infinity satisfies

    B_j=deg_t a_(2j)>=j>0.                           (6)

The pole b(t) may grow arbitrarily and is not assigned an artificial
bound. Likewise the scalar phase shift psi is unrestricted.

At a component with finite nonconstant U residue, some high coefficient
is polar by (6). The leading phase cannot be constant as a rational
function of that residue. Thus it cannot be a finite-phase component.
Every component specializing to a finite point, including U=0, has
strictly positive volume order equal to its positive displacement order.

It remains to examine U=infinity. Put v(U)=-N0<0, and let R0>0 be the
pole order of the leading spatial/scalar action after scaling
U=r^-N0 tau. The highest high-degree term ensures

    R0>=B_K+2K N0>4N0.                              (7)

If the leading face is a nonzero constant, no finite-phase chart occurs.
Otherwise at a finite nonzero centre tau0 it is an even Laurent
polynomial with powers among -2,0,2,...,2J. Let mu be its root
multiplicity at tau0. Multiplication by tau^2 makes it a polynomial in
tau^2 of degree at most J+1. Since tau0!=0, this substitution is
unramified, and

    mu<=J+1.                                        (8)

For J>=2, the leading coefficient a_(2J) has action
R0=B_J+2J N0 and B_J>=J. Therefore

    R0/mu>=R0/(J+1)>N0.

If J=1, (7) and mu<=2 imply the same strict inequality. If there is
no positive power, the nonconstant face is b*tau^-2 plus a scalar,
so mu<=1 and (7) again suffices.

In a final Puiseux disc around tau0 with displacement r^epsilon X,
the mu-th Taylor coefficient has exact leading order -R0+mu epsilon.
All lower-action groups have strictly larger order in this coefficient;
the full Puiseux centre does not cancel the nonzero leading mu-th
derivative. Finite phase consequently forces

    epsilon>=R0/mu>N0.

Its volume order is -N0+epsilon>0. This includes faces containing an
arbitrarily large pole coefficient b and an arbitrary scalar psi.

## 4. The complete comparison consequence

In both cases every finite-phase generic component of every final marked
model has strictly positive seed-volume order. The volume is regular at
finite horizontal punctures. Its pole at infinity is a DOUBLE pole,
rather than a logarithmic pole; infinity is phase-irregular for D>=4.
Thus a logarithmic-seed-only comparison is not sufficient.

Here the required far-infinity estimate is explicit. Let H=1/(lambda d)
be the scalar in eta and a_D its top phase coefficient. In the coordinate
v=a_D^(1/D)U, the volume prefactor is H a_D^(-1/D). Its infinity order is

    1-vtop/D>0                     if a!=0,
    B_K/(2K)>=1/2                  if a=0.             (9)

Division by p adds deg p. Descending exponential integration at this
irregular puncture absorbs the fixed double pole. The remaining rational
phase coefficients and algebraic inverse branches have only algebraic
parameter growth; a phase cutoff at a negative power of the parameter
therefore gives exponentially smaller continuation tails.

Use the retained
[complete compact-piece, crossover-annulus and descending-tail estimates](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#all-centres-punctures-and-crossover-estimates),
with the full proper marked model and every gluing map retained. On a
finite/finite node the positive orders on both branches give a positive
power times a logarithmic neck length. At a finite/pole node the positive
finite-side factor survives substitution into the descending tail; fixed
extra pole orders are absorbed by that exponential integration. Regular
finite punctures cause no higher-order amplitude pole. The positive
finite-component estimates in Sections 2-3 and the far-tail bound (9)
therefore cover all pieces, including successive cancellation scales.

The ordinary-parameter perfect pairing and zero-exponential comparison
are the retained inputs recorded in the
[fixed-puncture comparison proof](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#fixed-puncture-logarithmic-seed-bound)
and its
[complete period comparison](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#infinity-and-complete-period-comparison).
The explicit treatment of an irregular infinity pole is also given in
the retained
[relative leading-coefficient proof](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#a-relative-leading-coefficient-bound).
At an ordinary parameter, the normalized scalar functional extends from
its cyclic seed span to the full fibre and parallel-transports to an
analytic curve period. Its constant seed value belongs to the
zero-exponential part, whereas all such seed periods decay by the
complete estimates. This is the contradiction.

Hence the scalar premise forces D<=2, equivalently deg_z N<=2.
These are the same explicit complete-comparison dependencies as in the
retained curve proofs. Bounded calculations do not certify the imported
all-component statements, and no new review assurance is inferred.

## 5. Low degree forces a cubic phase of the retained architecture

Compare (2) coefficientwise under deg_z N<=2. The t^-1 coefficient
az^2+zh forces deg h<=1. The t^0 coefficient then gives

    deg k<=max(1,deg u),

because ell has degree one and c is constant. If deg u=m>=1, the
term ell u^2 in the t^1 coefficient has degree 2m+1, whereas every
other term has degree at most max(m+1,2). This is impossible. Hence

    h=h0+h1z, u=u0, k=k0+k1z.                       (10)

Swap the ORIGINAL Witt pair Z,W. This is now the exact multibase
quadratic-fibre architecture with one base coordinate x. Its exact source
formula follows from the
[self-contained multibase Wick contraction](GAUSSIAN_MULTIBASE_QUADRATIC_FIBRE_PRIMITIVITY.md#1-a-self-contained-multibase-wick-contraction):

    P=Wnew H(x)+aY^2+(cx+u0)Y+G(x),
    H=lambda x^2+h1x+k1, G=mu x^2+h0x+k0.

Let r(t) be the selected formal root r=tH(r), and put

    J=1-t(2lambda r+h1),
    J^2=(1-h1t)^2-4lambda k1t^2,
    R=tG(r)+t^2(cr+u0)^2/(2d).

The exact ORIGINAL vacuum is

    f=exp(R)/(J sqrt(d)).                           (11)

R is algebraic over C(t), not presumed rational. From f=p exp(psi),

    exp(2(R-psi))=p^2 d J^2.                         (12)

An algebraic function with algebraic exponential is constant. For a
self-contained proof, pass to the smooth projective algebraic curve
on which the function and its nonzero algebraic exponential live.
At a pole, the differential of the function has pole order at least
two, while the logarithmic differential of a rational function on
that curve has only simple poles. The identity d(exp S)/(exp S)=dS
rules out all poles of S. A regular function on the projective curve
is constant. This argument applies to (12) and is valid also if the
right side is already rational.

Since R(0)=psi(0)=0, (12) gives R=psi and

    p^2 d [(1-h1t)^2-4lambda k1t^2]=1.

Every factor is polynomial with initial value one. Thus p=1, a=0,
h1=0 and k1=0. The selected root is r=0, so

    psi=k0 t+u0^2 t^2/2.

## 6. Sufficiency and original markers

If a=0 and h,u,k are constant, write

    P=u0Y+k0+Pminus,
    Pminus=(lambda Z+mu)W^2+cWY+h0W.

Every monomial of Pminus has strictly negative Z-minus-W angular
weight. For an original marker Q of degree D, angular balance truncates
the exponential of Pminus after at most D factors. Translate the
ordinary Gaussian Y by tu0. It follows directly that every original
marked period is exp(k0t+u0^2t^2/2) times a polynomial. The unit is
exactly the exponential factor. Thus the vacuum is primitive and M
is finite. The retained
[nonconstant-phase original-image theorem](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#7-every-finite-actual-image-and-its-original-fixed-marker-image)
gives literal V=M. Safety is equivalent to this exponential factor being a
polynomial, hence u0=k0=0; then the phase is strictly negative angular
and pure.

All deductions are algebraic except the explicitly retained complete
curve comparison. Necessity descends to any characteristic-zero field
by embedding the finitely generated field of phase coefficients and
hypothetical p,psi into C. Sufficiency is algebraic over the original
field. No finite-jet inference, external review or formal verification
is asserted.

## 7. Verification and scope

The [focused exact checker](verify_gaussian_affine_w_square_finite_images.py)
and [receipt](gaussian_affine_w_square_finite_images_verification.json)
retain generic Schur/numerator identities, eight exact high-coefficient
controls with c=0 and c!=0, three independently computed original Wick
vacuum comparisons through moment four, and the generic selected-root
discriminant identity. They check bounded algebra and normalizations.
The all-degree coefficient bounds, all-final-chart multiplicity argument,
complete comparison and scalar algebraic-exponential rigidity are the
written proof above.

The result covers every constant c, including zero, with lambda!=0 and
arbitrary-degree h,u,k. It does not classify arbitrary nonlinear ell or
nonconstant c. The alternate scalar bridges retain the actual scalar
seed and all its parameter-covariant derivatives; they do not identify
original W-marked periods with alternate localized values. The direct
angular proof supplies every original marker in the finite branch.
Imported primary-source audits and pending independent-review status
remain unchanged. Independent replay, formal verification, external
mathematical review, publication and priority are not asserted.
