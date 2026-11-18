# Finite actual images for every variable W-square coefficient

Use one Witt Gaussian pair Z,W with E[ZW]=1 and an independent ordinary
Gaussian Y with E[Y^2]=1. All other basic pairings are zero. Let K be
a characteristic-zero field and consider

    P=ell(Z)W^2+W(h(Z)+c(Z)Y)+aY^2+u(Z)Y+k(Z),

where a is constant and ell,h,c,u,k are polynomials with deg ell>=2.
There is NO restriction on c, including c=0, and no root-multiplicity
bound. The coefficients u,k are arbitrary.

For the ORIGINAL fixed-marker periods and their polynomial-time image,

    Phi_P[Q]=E[Q exp(tP)], V_P=Phi_P[K[Z,W,Y]],
    M_P=K[t]V_P,

the complete criterion is

    M_P finite over K[t]
      iff a=0, ell in Z^3K[Z], h,c in Z^2K[Z].      (V1)

Every finite image has primitive vacuum and literal original image

    psi=k(0)t+u(0)^2t^2/2,
    Phi_P[1]=exp(psi), V_P=M_P=exp(psi)K[t].         (V2)

Safety is equivalent to (V1) plus u(0)=k(0)=0. In particular every
safe member of this architecture is pure. Section5 gives the original
mixed cutoff C=1 for safe phases and the normalized degree bound C=2
for every finite image.

This is a WHOLE-actual-image theorem. It does not infer a general
scalar-unit-only or purity-only classification, and makes no claim
for deg ell<2. Independent review and formal verification are not asserted.

For c!=0, (V1)-(V2) are the
[complete mixed-source theorem](GAUSSIAN_QUADRATIC_W_MIXED_FINITE_IMAGES.md),
whose degree condition is automatic here. The body below proves the
remaining c=0 branch. Thus henceforth

    P=ell(Z)W^2+h(Z)W+aY^2+u(Z)Y+k(Z).

For the comparison argument work first over C; Section4 gives the
characteristic-zero descent. The usual polynomial logarithmic marker reduces ell to fixed-cover
multiplicities one or two. Its only lower-degree exception is

    ell=kappa(z-alpha)^m, kappa!=0, m>=3 odd.        (O1)

For this exception a marker one power lower moves the double volume
pole to a finite point where the phase is irregular. Outside the
common-zero-origin case alpha=h(0)=0, this rules out finite M_P.

**Exact comparison inputs.** We use the three underlying comparison
facts listed in the retained
[fixed-puncture discussion](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#fixed-puncture-logarithmic-seed-bound),
together with the complete marked-model and
[finite/pole node estimate](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#a-double-root-the-remaining-quadratic-mixed-coefficient)
that integrates exp(-s) times any fixed algebraic power.
The transition-annulus bound is
[analytic_core.tex](../../../results/gaussian-boundary/analytic_core.tex),
Lemma `lem:curvelength`; the rational-unit restriction and ordinary-
parameter analytic extension are
[seed_period_comparison.tex](../../../results/gaussian-boundary/seed_period_comparison.tex),
Lemma `lem:seedfunctional`.
The double pole below does NOT satisfy the simple-pole hypothesis of
the fixed-logarithmic seed bound itself. We use its underlying complete
comparison facts, with the higher-pole estimate explicitly proved in
Section3. Their existing primary-source audit and pending independent-
review status remain unchanged.

## 1. Original marker, normalized curve and exact phase

Set r=1/t and normalize the fixed cover by

    z=alpha+w^2,
    U=sqrt(2kappa(r-2a)) w^m.

The original seed eta=r dz/U and the fixed ORIGINAL polynomial marker

    sigma(z)=(z-alpha)^((m-3)/2)

give exactly

    sigma eta=H(r) dw/w^2,
    H(r)=2r/sqrt(2kappa(r-2a)),
    nu=v(H)=1/2 if a=0, nu=1 if a!=0.               (O2)

The marked volume is holomorphic at w=infinity and at every finite
w!=0. Its only higher horizontal pole is the DOUBLE pole w=0.

The reduced phase after a scalar polynomial shift psi(t) is

    varphi=-(r(alpha+w^2)-h(alpha+w^2))^2/(4kappa r w^(2m))
          +k(alpha+w^2)/r
          +u(alpha+w^2)^2/[2r(r-2a)]-psi(1/r).       (O3)

At w=0 its actual leading spatial pole coefficient is

    b(r)=-(r alpha-h(alpha))^2/(4kappa r).

It is nonzero in C(r) precisely when (alpha,h(alpha))!=(0,0).
Thus under the stated outside-origin condition the double volume pole
is ALWAYS phase-irregular; subtracting psi cannot turn it into a
regular puncture with an uncontrolled double-pole circle.

The ORIGINAL source survives localization at z=alpha. The restricted
W-source is t h(alpha)-alpha, a nonzero scalar in C(t). The operator
partial_W+t h(alpha)-alpha is bijective on polynomials by its finite
geometric derivative inverse, so the full restricted exponential
de Rham complex is acyclic. Smooth-divisor residue and higher-pole
comparison therefore give a horizontal localization isomorphism.
The remaining determinant factor is a scalar in C(t). Standard Gaussian
fibre elimination transports the ORIGINAL marker sigma to (O2).

This period is nonzero. If alpha!=0 then sigma(0)!=0. If alpha=0 and
h(0)!=0, write s=(m-3)/2. Only h(0)W has negative angular weight;
ell W^2 has positive weight and every other monomial has nonnegative
weight. Consequently

    Phi_P[sigma]=h(0)^s t^s+O(t^(s+1))!=0.           (O4)

For m=3 this simply says the unit starts with1.

Assume finite M. The general actual-image theorem gives exp(psi)K[t],
so p_sigma=exp(-psi)Phi_P[sigma] is a NONZERO polynomial. Dividing
(O2) by p_sigma adds deg p_sigma>=0 to its infinity valuation, even
if p_sigma(0)=0. The shifted actual functional takes that normalized
seed to one and supplies a rational cyclic unit quotient. It suffices
to prove positive volume for (O2), before that additional division.

## 2. Every finite-phase component has positive order

Away from w=0 the fixed form dw/w^2 is holomorphic, including w=infinity.
On every model component centred there, its relative order is nonnegative,
and multiplication by H gives strictly positive order. This covers
ordinary finite components, shrinking ordinary discs and spatial infinity,
irrespective of psi. Only components approaching w=0 need phase information.

Put N=v(w)>0. In an exact final chart

    w=r^N(tau_c+r^epsilon X), epsilon>=0,

the volume order is nu-N+epsilon. For N<1/2 it is already positive.

A leading phase face of weight -R<0 must have a nonzero root tau0.
If that root has multiplicity mu, the exact-chart argument gives

    epsilon>=R/mu.                                  (O5)

Indeed its X^mu coefficient has exact valuation -R+mu epsilon at
the final unit Puiseux centre. The coefficient remains a unit; every
strictly higher original level remains higher after finitely many
tau derivatives. Thus all subsequent cancellations/refinements are covered.

The spatially regular terms in (O3) have action at most2 when a=0,
and at most1 when a!=0. At N=1/2 only a constant from the u term can
have action2; its w-dependent corrections have strictly higher valuation.
The linear part of psi can be absorbed into k, and every remaining
scalar leader contributes just a constant to the selected face.

### Case h(alpha)!=0

For N>=1/2 the unique original leading monomial is

    -h(alpha)^2/(4kappa r w^(2m)),
    R=1+2mN>=m+1>2.

All other square terms are strictly higher, as are the regular u,k terms.
Only a scalar psi can join it. If that scalar is strictly lower it is
alone and has no finite-phase component; if it ties, the face is a
binomial in tau^(-2m) and a constant, whose nonzero roots are simple.
The same is true if no scalar leads: a sole monomial has no nonzero root.

Thus every possible finite component has epsilon>=R and

    nu-N+epsilon>=1/2-N+1+2mN>0.                     (O6)

This includes N=1/2.

### Case h(alpha)=0, hence alpha!=0, and N>1/2

Now the unique original leader is

    -r alpha^2/(4kappa w^(2m)), R=2mN-1>m-1>=2.

Strictness R>2 holds even when m=3. If h has a simple z-root at alpha,
the next square weights differ by2N-1 and4N-2; higher root orders only
increase those differences. Terms from r w^2 are also higher.
Again only scalar psi can join, with simple binomial roots. Hence

    nu-N+epsilon>=1/2-N+2mN-1
                 =-1/2+(2m-1)N>0.                 (O7)

### The boundary h(alpha)=0, alpha!=0, N=1/2

Write h1=h'(alpha). The exact leading rational face of the square term is

    F(tau)=-(alpha-h1 tau^2)^2/(4kappa tau^(2m)),
    R=m-1>=2.                                      (O8)

It is nonconstant because alpha!=0. The scalar psi can add a constant.
When m=3,a=0, the leading u(alpha)^2/(2r^2) term can also add a constant;
every other regular-fibre term is higher. A scalar of larger action
would be the only leader and give no finite-phase component.

Every nonzero root of F plus ANY constant has multiplicity at most two.
Indeed

    F'=(alpha-h1 tau^2)
        *(m alpha-(m-2)h1 tau^2)/(2kappa tau^(2m+1)).

If h1=0 there is no nonzero critical point. Otherwise the first factor's
roots are simple and give double zero-level roots. The second factor's
roots are simple, distinct from the first, and satisfy
alpha-h1 tau^2=-2alpha/(m-2)!=0. Thus F'' is nonzero there.
This proves the stated multiplicity bound at every shifted level.

Consequently epsilon>=(m-1)/2 and

    nu-1/2+epsilon>=(m-1)/2>0.                      (O9)

There is no missing quartic or higher degeneracy at this boundary.

These alternatives exhaust all N and every scalar degree. Components
with nonconstant leading tau residue at N=1/2 would carry the nonzero
negative phase F(tau)+constant and cannot be finite-phase components.
The exact root charts above treat every remaining refinement.

## 3. Complete neck and high-pole tail estimates

Use the retained complete marked model: separate the fixed punctures,
the zero, pole and value-one sections of varphi, resolve the parameter
boundary and retain the full annulus/gluing complex. The finite-phase
component bounds just proved are strict. On finitely many compact pieces
they therefore give a common positive power of r.

At a finite/finite node uv=r (after the usual ramification and semistable
resolution), the marked puncture sections have been separated from the
node. Write the relative form as f(u,v)du/u. Its horizontal poles are
absent there; positive generic order on both adjacent components makes
f a positive monomial in u,v times a bounded unit. The annulus integral
is therefore a positive power of r times O(1+|log r|).

At a finite/pole node, let the finite component have positive order a0.
A local phase coordinate and u=r/v write the form as

    r^a0 v^b dv/v

times a bounded unit, for a fixed integer b. A descending phase coordinate
T proportional to v^(-p), p>=1, changes the tail integral to r^a0 times
an exponentially convergent integral of e^(-T) times a fixed power of T.
Thus the positive finite-component order survives despite a higher
horizontal pole.

For the actual puncture w=0, the deepest tail has
varphi~b(r)w^(-2m). On a descending inverse branch,

    dw/w^2=constant*b(r)^(-1/(2m))
             T^(1/(2m)-1)(1+o(1))dT.

This gives only algebraic parameter growth and a fixed power of T,
absorbed by exponential integration. In fact the scalar prefactor
H b(r)^(-1/(2m)) already has positive valuation:
nu-1/(2m)>0 when h(alpha)=0 and nu+1/(2m)>0 otherwise.
All transition components where the coefficient vanishes or changes
dominance were included in Section2.

Truncate each node tail before leaving its chart, where the phase is
already bounded above by -|r|^(-delta). On the remaining finitely many
transition annuli, inverse-phase paths have bounded algebraic degree.
The retained dyadic curve-length bound gives at most fixed algebraic
parameter growth and logarithmic length; the exponential cutoff defeats
it. Regular puncture circles occur only where the form is holomorphic
or already has the positive component order. The double pole w=0 remains
phase-irregular, so it supplies no omitted regular-circle contribution.

These are the same underlying complete component/annulus/tail estimates
used in the retained H5 and cover arguments, with the actual fixed order2
pole explicitly integrated. The fixed-log seed lemma is NOT invoked
outside its simple-pole hypothesis.

Therefore every zero-exponential normalized sigma-eta period tends to
zero. The rational unit quotient supplied by finite M would have value
one. The retained ordinary-parameter extension and full perfect period
pairing give the contradiction. This proves no finite M in (O1) outside
the common-zero-origin case.

## 4. Consequence for every c=0 coefficient ell of degree at least two

For arbitrary ell with deg ell>=2 set g=ell up to a nonzero scalar and
rho=product_(ell-root alpha)
       (z-alpha)^floor((ord_alpha ell-1)/2).

The reduced ell/rho^2 has root multiplicities1 or2. Its degree is at
least two EXCEPT exactly when ell has a SINGLE root of ODD multiplicity.
Thus for all other ell the weighted-cover VOLUME ARGUMENT applies
with its ordinary marker rho; the c!=0 classification is not being
used as a black box here. The exceptional ell are exactly (O1), now proved
using sigma. The original fixed-clean, unit-ell and common-simple-origin
source constructions are unchanged.

More explicitly, on the fixed normalization v^2=ell/rho^2 the marked
form is a nonzero constant times r/sqrt(r-2a) dz/v. The fixed form dz/v
is regular at simple roots, logarithmic at double roots, and logarithmic
or holomorphic at infinity when deg(ell/rho^2)>=2. Its pullback has
nonnegative relative order on every marked-model component. The scalar
has positive order1/2 or1, giving the complete volume estimate with all
positive-genus cycles retained.

The source case split is exhaustive. If ell(0)!=0, use the
[full formal original source](GAUSSIAN_UNIT_ELL_FORMAL_SOURCE_FINITE_IMAGE_EXCLUSION.md#1-completed-z-wick-contraction-preserves-the-full-original-source);
rho(0)!=0. If ell(0)=0,h(0)!=0, every fixed root has nonzero W-source,
so the same exact fixed-divisor localization used in Section1 retains
Phi_P[rho]; if rho vanishes to order s at zero, its first term is
rho_s h(0)^s t^s!=0 by angular balance. If ell(0)=h(0)=0 and the raw
ell order is one, the
[simple-origin scalar/base-Z bridge](GAUSSIAN_SIMPLE_ORIGIN_ALTERNATE_SCALAR_EXCLUSION.md)
retains Phi_P[rho] by inserting rho(z) into its original residue formula;
here rho(0)!=0. These source constructions require no multiplicity bound.
At a common origin with higher raw ell order, use the angular argument
that follows. No original W-marker localization is asserted in the
simple-origin scalar construction.

The comparison phases in these cases are spatially nonconstant. A clean
fixed root supplies an actual pole of -(z-th)^2/(4t ell), which the
spatially regular ordinary-fibre terms cannot cancel. At a common raw
simple origin, the t^-1 coefficient -z^2/(4ell) has a simple zero at
zero and is not constant. A spatially constant phase can therefore
occur only in the higher common-origin case handled algebraically below.

At a common origin with raw ell order at least two, the weakly one-sided
quadratic angular core forces, under finite M,

    a=0, ell_2=h'(0)=0.

Equivalently ell belongs to Z^3K[Z] and h belongs to Z^2K[Z]. Conversely
these conditions leave the zero-weight part u(0)Y+k(0) and a strictly
positive angular remainder. The same exact finite expansion after the
Gaussian shift Y->Y+t u(0) gives ALL marked periods as exp(psi) times
polynomials, with primitive

    psi=k(0)t+u(0)^2t^2/2.

The original nonconstant-phase image theorem then gives literal V=M.
Therefore the COMPLETE c=0, deg ell>=2 criterion is

    M finite iff a=0, ell in Z^3K[Z], h in Z^2K[Z]. (O10)

Safety additionally requires u(0)=k(0)=0, and every safe member is pure.
Together with the c!=0 theorem, (O10) settles every constant-a phase of
the displayed quadratic-W architecture with deg ell>=2, without a gcd
multiplicity restriction.

Necessity descends to every characteristic-zero field by the finite
phase, marker, psi and p_sigma (or p_rho) coefficient embedding. The
half-power normalization requires at most a finite algebraic extension.
Sufficiency is the algebraic original Gaussian tilt/weight formula.
The imported comparison inputs keep their existing internal assurance
and pending independent-review status; this argument adds no external
review or formal verification claim.


## 5. Primitive finite image and original cutoff bounds

Under (V1) the phase is

    P=u(0)Y+k(0)+Pplus,

where EVERY monomial of Pplus has strictly positive Z-minus-W weight.
For an original marker Q of degree D, angular balance truncates the
positive-part expansion after D factors. Ordinary Gaussian translation
gives the exact ALL-marker formula

    Phi_P[Q]=exp(psi)
       *sum_(n=0)^D t^n/n!
          E[Q(Z,W,Y+t u(0))
            Pplus(Z,W,Y+t u(0))^n].                (V3)

In particular the unit is exactly exp(psi), which generates M_P.
The phase is nonconstant since deg ell>=2, so the
[finite original-image theorem](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#7-every-finite-actual-image-and-its-original-fixed-marker-image)
gives V_P=M_P literally.

The ordinary Y degree of Pplus is at most one. For a marker monomial
Z^iW^jY^b, balance forces n<=j-i<=j. The Gaussian shift contributes
t degree at most b+n; including t^n yields at most b+2n<=2D.
Consequently

    deg_t(exp(-psi)Phi_P[Q])<=2deg Q.               (V4)

If u(0)=0 there is no shift and the bound improves to deg Q.
Safety makes exp(psi) polynomial, forcing u(0)=k(0)=0; conversely
these conditions make every original period polynomial. Then the whole
phase has positive angular support, so

    E[Q P^n]=0 for n>deg Q.

Thus C=1 is a sufficient ORIGINAL mixed cutoff throughout the safe
family. The
[cutoff-to-jet theorem](MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md#1-the-finite-jet-theorem)
gives moment r-jets through (d+1)r for perturbations of original
spatial degree at most d. This is a jet statement at a known safe
phase, not an arbitrary-point finite purity test.

## 6. Verification and assurance

The [focused exact checker](verify_gaussian_variable_w_square_finite_images.py)
and [receipt](gaussian_variable_w_square_finite_images_verification.json)
record bounded original-moment, fixed-cover, leading-face and tail
identity controls. The unrestricted proof is the original-source,
exact-scale, complete-comparison and angular argument above; a finite
moment prefix does not establish it. The higher horizontal volume pole
is treated by its phase-irregular tail, not by silently applying a
logarithmic-only bound. Imported comparison statements retain their
existing internal audit and pending independent-review status. No
external review, formal verification, publication or priority claim
is added.
