# Polynomial Gaussian time changes and exact ordinary cutoff regions

16 November2026. Internal all-degree/all-ring proof, with finite exact
Wick and ordinary-ideal controls. Independent replay, formal verification
and external review remain false; priority unassessed. This construction
applies to EVERY polynomial Gaussian phase, without a weight or purity
premise. It is a reduction/transport theorem, not an exhaustive safe
normal-form classification. See the [checker](verify_gaussian_time_changes.py),
[certificate](gaussian_time_change_certificate.json) and
[source audit](gaussian_time_change_source_audit.json).

## 1. An actual polynomial phase realizes polynomial time substitution

Let B(X) be any polynomial over a commutative Q-algebra R, with the
original nonsingular rational Gaussian covariance. Write

    F_Q(s)=E_X[Q exp(sB)].

Adjoin an independent circular Gaussian pair U,A, with E[U^i A^j]
equal to i! if i=j and0 otherwise. The joint precision has hyperbolic
block UA. In standard real coordinates one may take U=x+iy,
A=(x-iy)/2; adjoining i is faithfully flat when necessary.

For ANY polynomial a(s) in sR[s], put C(s)=a(s)/s and

    B_a(X,U,A)=A+C(U)B(X).                              (1)

The phase and every marker are original finite polynomials. For integers
alpha,beta>=0 one has the exact ORIGINAL-amplitude identity

    Phi_a[U^alpha A^beta Q]
      = [partial_z^beta {z^alpha F_Q(t C(z))}]_(z=t).    (2)

Indeed expand exp(tA) and exp(t C(U)B) coefficientwise and use
E[U^i A^j]=delta_ij i!. Equivalently E[A^beta h(U)exp(tA)]
=h^(beta)(t) for a polynomial h, then extend t-adically. At each time
order only finitely many polynomial terms occur. There is no analytic
interchange or convergence assumption. Chain rule in(2) expresses every
image as a FINITE R[t]-linear combination of F_(B^k Q)(a(t)); all
coefficients are polynomial, including when a has ramification.
Conversely beta=0 gives t^alpha F_Q(a(t)). Hence the FULL actual global
polynomial-Weyl image is exactly

    M_(B_a)=sum_Q R[t] F_Q(a(t)).                       (3)

The sum uses all original X markers. Closure under partial_t follows
from partial_t F_Q(a)=a'(t)F_(BQ)(a). Formula(3) proves equality of
actual images, not merely a map from a larger annihilating source.
No image-kernel computation is required. When a(t)=t, the global image
is unchanged, though the original marker-degree filtration is enlarged.

## 2. General finite polynomial clocks reflect full safety

First work over a characteristic-zero field K, with a nonconstant.
Forward safety follows immediately from(2). For the converse suppose
EVERY marked image of B_a is polynomial. Then F_Q(a(t)) is polynomial
for every original Q. Since t is integral over K[a(t)], every F_Q(s)
is algebraic and integral over K[s]. The same holds for ALL derivatives

    partial_s^j F_Q=F_(B^j Q).                         (4)

There is already an elementary polynomial-degree proof of descent.
Put q_j(t)=F_Q^(j)(a(t)), all polynomial by hypothesis, and d=deg a.
The chain rule gives q_j'=a' q_(j+1). If q_j is nonconstant of degree N,
then q_(j+1) has degree N-d. This cannot continue indefinitely. Thus
some q_j is constant and q_(j+1)=0. Formal substitution by a nonconstant
zero-origin polynomial is injective over a field, so F_Q^(j+1)=0 and
F_Q is polynomial. No regularity theorem is needed for this argument.
The following geometric lemma gives the complementary finite-pole
interpretation and does not assume the existence of a polynomial clock.

The following classical algebraic-function argument now applies.

**Integral derivative lemma.** If f in K[[s]] is algebraic over K(s)
and every derivative partial_s^j f is integral over K[s], then f lies
in K[s]. To prove this, pass to geometric places of K(s,f). At a finite
point b choose a local parameter z with s-b=z^e. Integrality of f gives
a nonnegative Puiseux expansion. If a fractional exponent n/e occurs,
take the least such exponent after the integral powers. More than n/e
derivatives produce a nonzero negative-power term; no distinct exponent
can cancel it. This contradicts integrality of that derivative. Thus
f has only integer powers at EVERY finite place, and K(s,f)/K(s) is
geometrically unramified over the affine line.

For completeness, characteristic-zero Riemann--Hurwitz rules out a
nontrivial connected degree d cover ramified only over infinity:

    2g-2=-2d+sum_(v|infinity)(e_v-1)<=-d-1.

For d>1 this is less than-2. Thus d=1 geometrically, and f is rational
with no finite poles, hence polynomial. Constants descend because f has
coefficients in K (or by the finite-column linear relation argument).
The local derivative and characteristic-zero/tame hypotheses are essential.
The exact classical formula and tame multiplicities are given in
[Stacks,53.12](https://stacks.math.columbia.edu/tag/0C1B).

Applying the lemma separately to every Q proves, WITHOUT PURE PREMISE,

    B_a safe iff B safe.                               (5)

This uses the ENTIRE marked derivative family, not just one algebraic
series. The countercontrol f(s)=sqrt(1+s), a(t)=t^2+2t has f(a)=1+t,
but f'(a)=1/[2(1+t)] is not polynomial. A single selected series can
be polynomialized by a finite polynomial substitution.

Over EVERY Q-algebra assume the leading coefficient of a is a unit, so
every prime-residue clock has the same positive degree. The retained
[full mixed all-ring descent](GENERAL_GAUSSIAN_FAMILY_LOCUS_AND_NORMALIZATION.md)
then extends(5): safety is equivalent to safety over every prime residue
field, for BOTH original polynomial Gaussian phases. No purity,
Noetherian hypothesis or bounded whole nilradical is added.
Applying this equivalence over R/Ann(c) gives the EXACT normalization
identity

    N_(B_a)=N_B,
    N_B={c: c F_Q is polynomial for EVERY Q}.            (6)

If also the lowest nonzero coefficient of a is a unit, substitution is
formally injective over every R. The unit period is F_1(a), so purity is
preserved and reflected. Triangular comparison at degrees ord(a)*n
also shows equality of the FULL ORDINARY pure ideals. For a general
clock no same-prefix or prescribed mixed-cutoff ideal equality is claimed.
The leading-unit condition in the all-ring descent is explicit, not
inferred from a nonzero coefficient in a ring with zero divisors.
It cannot simply be dropped: over R=Q[e]/e^2, take B=X and a(t)=et.
The original real linear Gaussian phase is nonpure and unsafe, whereas
B_a=A+eX is pure and safe. Its unit period is1; every original marked
period is polynomial by(2) and the finite expansion in e. Its clock
vanishes on the reduced ring. This is a genuine Gaussian countercontrol,
separate from the scalar sqrt example.

## 3. Monomial clocks retain exact individual moments

Now a(t)=t^r, with integer r>=1, so B_r=A+U^(r-1)B. Put
M_n(Q)=E_X[QB^n]. For each original monomial marker U^alpha A^beta Q,

    Phi_r[U^alpha A^beta Q]
      =sum_(n>=0) (alpha+(r-1)n)_beta * M_n(Q)/n!
                         *t^(alpha+rn-beta).            (7)

Here (x)_beta is the falling factorial. A term is zero when
alpha+(r-1)n<beta; these are exactly the would-be negative powers
of exp(tA)'s index. Thus(7) is an ordinary formal power series even
when the displayed exponent could otherwise be negative.
Equivalently its ORIGINAL m-th moment is zero unless

    n=(m-alpha+beta)/r is an integer>=0,
    alpha+(r-1)n>=beta;

when these conditions hold it is

    (m!/n!)*(alpha+(r-1)n)_beta *M_n(Q).                (8)

The nonzero scalar is a rational UNIT over every Q-algebra. These are
identities of individual ordinary moments, with all nilpotents retained.
In particular the first N pure moments of B_r generate EXACTLY the
ideal of the first floor(N/r) pure moments of B. No new pure test is
inferred from a finite numerical prefix.

## 4. Exact ordinary mixed ideals and sharp real cutoff regions

For C,E>=0 let I_B(C,E) be the ordinary ideal of all original marked
moments M_m(X^gamma) with m>=1 and m>C|gamma|+E. It is an ideal in the
ACTUAL coefficient ring, before reducing or taking radicals. Then

    I_(B_r)(C,E)=R,                         0<=C<1,
    I_(B_r)(C,E)=I_B(C/r,E/r),              C>=1.       (9)

For C<1 choose the marker U^alpha, with alpha>C alpha+E. Its alpha-th
moment is alpha!, a unit. For C>=1 a lifted forbidden term in(8), with
q=deg Q, satisfies

    rn>Cq+E+(C-1)alpha+(C+1)beta>=Cq+E.

The index n cannot be0. Thus every lifted generator is a unit multiple
of a base forbidden moment. Conversely take alpha=beta=0 and m=rn
to lift EVERY base forbidden moment, proving both ordinary inclusions.
This proves scheme-level equality for all parameter rings, primary
thickenings and embedded components, not just equality of safe fields.

Consequently the entire real affine sufficient-cutoff region is exactly

    Region(B_r)={(C,E): C>=1,(C/r,E/r) in Region(B)}.    (10)

Degree filtrations and nonempty-region descriptions assume R is nonzero;
the zero ring has all moments zero and is trivial. The ordinary ideal
identity(9) itself also holds over the zero ring.

If B has sharp asymptotic affine slope c, its infimum for B_r is
max(1,rc); the same formula holds for the infimum pure-linear coefficient
when the linear region is nonempty. Infima and their attainment are
distinguished by(10). If r*c<1 the boundary C=1 corresponds to a base
coefficient strictly above its infimum; no endpoint is invented. Safety,
purity and normalization are therefore preserved/reflected over EVERY R.

There is an equally exact original-marker filtration formula. Write
h_B(D)=max_(deg Q<=D) deg_t F_Q, allowing infinity for an infinite
tail. Finite-dimensionality makes a finite maximum attained. Then

    h_(B_r)(D)=D+max_(0<=q<=D){r h_B(q)-q}.             (11)

The upper bound follows from(7): momentum degree beta can only lower
the time degree, giving r h_B(q)+alpha-beta. The lower bound uses
U^(D-q)Q, beta=0. Formula(11) is valid over arbitrary coefficient
rings too; nonzero integer factors are units. It proves sharpness using
ORIGINAL markers, without identifying a generic lattice with a filtration.

For B=eta Z^d in an original circular pair Z,W, the nonzero moments
have q=dn and are (dn)! eta^n. If C>=1 put delta=r-Cd. For delta<=0
the ideal in(9) is0; for delta>0 it is EXACTLY

    (eta^(floor(E/delta)+1)).                           (12)

For nilpotence order L and r>d, (11) becomes
D+(r-d)min(floor(D/d),L-1); when r<=d it is D. This retains literal
ordinary powers, not only their common radical. With eta a unit the
sharp coefficient is max(1,r/d). The base family is already known safe;
the new result is the all-phase transport and exact cutoff/ideal law.

### The ordinary cutoff equality extends to EVERY finite polynomial clock

In fact(9)-(10), with r replaced by d=deg a, hold for EVERY polynomial
clock a with zero constant term and UNIT leading coefficient over EVERY
Q-algebra. The monomial formula(8) is the stronger individual-moment
description available only for monomial clocks.

To prove the extension, first suppose B is safe. Each F_Q is a finite
polynomial. The term of index n in(2) has time degree at most
alpha+d*n-beta; it vanishes if alpha+(d-1)n<beta. With beta=0 its
top coefficient is multiplied by the UNIT leading coefficient of a
to the power n. Hence the SAME exact marker-filtration formula is

    h_(B_a)(D)=D+max_(q<=D){d h_B(q)-q}.               (13)

It follows directly that, for C>=1,E>=0, the clock cutoff holds iff
the base cutoff with C/d,E/d holds. For C<1 the marker U^alpha at
time alpha again has moment alpha!, irrespective of B or lower clock
coefficients, so the clock cutoff is impossible over a nonzero ring.

Now let J and I be the ordinary clock and base cutoff ideals in the
given coefficient ring. Modulo I the base bound holds, so B is safe;
(13) proves the clock bound and J subset I. Modulo J the clock bound
holds, so the clock is safe; the already proved ALL-RING descent(5)
makes B safe in this quotient. Formula(13) then gives the base bound,
so I subset J. The leading coefficient of a remains a unit in both
quotients. Thus

    I_(B_a)(C,E)=R,                       0<=C<1,
    I_(B_a)(C,E)=I_B(C/d,E/d),            C>=1,         (14)

as literal ORDINARY ideals, without purity, reduction or Noetherianity.
This quotient argument uses the retained full mixed all-ring descent;
it does not infer ideal equality from a fieldwise radical equality.
It yields the complete real affine/linear cutoff-region transport by
d for arbitrary polynomial clocks. The pure-prefix dilation still uses
ord(a), under its lowest-coefficient-unit hypothesis, rather than deg a.
Global safety/cutoff transport is therefore stronger than a claim about
the generic rank: selected ranks can collapse for nonmonomial clocks.

## 5. General generic transport and faithful monomial pullback

For a general degree-d polynomial clock over a field, let V and V_a
be the original and new ACTUAL generic connections. Equation(3) gives
a surjection from the finite pullback of V to V_a. Formal substitution
is also an injection V --> Res_(K(t)/K(s)) V_a, where s=a(t) and the
derivation on the restriction is (1/a'(t))*partial_t. The second map
is K(s)-linear and differential; it uses only injectivity of substitution
on formal Laurent functions, not K(t)-linear independence. Thus

    rank(V_a)<=rank(V)<=d*rank(V_a).                   (15)

Finite tame pullback and restriction of scalars preserve regularity,
as do subquotients: formally they multiply/divide slopes by the positive
ramification index. These TWO actual-module maps consequently prove

    V regular on the whole P1 iff V_a regular there.   (16)

Locally V is regular at a base point iff V_a is regular at EVERY
preimage. No equality of individual positive slopes is asserted for
arbitrary clocks. Algebraicity of all marked periods is likewise
preserved and reflected, since K(t)/K(a(t)) is algebraic. Therefore
the existence of regular transcendental marked periods is preserved
and reflected by EVERY nonconstant polynomial clock, even when the
selected generic rank drops. This entire transport cannot create a
regular transcendental core from one without it.

Monomial clocks give a stronger statement, with no generic quotient loss.

Over a characteristic-zero field put u=t^r. The field extension
K(t)/K(u) has basis1,t,...,t^(r-1). If a rational relation between the
F_Q(t^r) is expanded in this basis, its r residue classes of Laurent
exponents give r relations between the ORIGINAL F_Q(u). Therefore

    K(t) tensor_(K(u)) M_B --> K((t)),
    b tensor f |--> b f(t^r)

is INJECTIVE. This uses the Laurent-series residue decomposition, not
an unsupported assertion that every selected algebraic pullback is
faithful. Combining this with(3) proves that the full actual generic
connection of B_r is the genuine finite pullback of that of B, with
EXACTLY THE SAME generic rank. In particular this construction cannot
hide an irregular source behind a smaller actual quotient.

Every geometric formal slope at0 or infinity is multiplied by r,
with its multiplicity unchanged. Away0/infinity the map is unramified;
the slopes occur unchanged over each preimage. Thus global regularity
is preserved AND reflected. Over algebraically closed constants the
classical Picard--Vessiot correspondence also preserves the differential
Galois identity component under this finite algebraic base extension:
the new group has finite index in the old one. Neither statement
identifies the full global boundary lattice with an abstract inverse image.
It is (3) that specifies the actual global lattice here.

For arbitrary polynomial a of degree larger than its order at0, the
Laurent residue proof no longer applies. Selected algebraic branches may
collapse after substitution, as the sqrt countercontrol shows. No general
actual-rank preservation or faithfulness is claimed for those clocks.

This boundary also has an ACTUAL Gaussian realization by the retained
[universal algebraic encoder](ALGEBRAIC_FIXED_POINT_LIFTS_AND_UNIVERSAL_GAUSSIAN_IMAGES.md).
Realize sqrt(1+s)-1 with full generic field K(s,sqrt(1+s)), rank2 and
unit period1. Applying the clock a(t)=t^2+2t makes every original field
element rational in t. Formula(3) makes the full new actual rank1, since
the unit remains1. Thus safety reflection does not justify rank preservation
for general clocks. No smaller dimension or newly computed phase is claimed.

## 6. Published knowledge and remaining general work

The integral-derivative descent uses classical characteristic-zero
Puiseux ramification and Riemann--Hurwitz, with the exact tame formula
checked in the [Stacks source](https://stacks.math.columbia.edu/tag/0C1B).
The Wick/source substitution and scalar-extension calculations above are
explicit proofs; no new general classification is imported.

[Hien--Roucairol,BSMF136(2008),505-532](https://www.numdam.org/articles/10.24033/bsmf.2564/)
gives rapid-decay integral representations for exponential Gauss--Manin
systems. Its selected [v1 Section3.1](https://arxiv.org/html/0704.1739v1#S3.SS1)
constructs the twisted direct image and distinguishes the regularity
question. It supplies no theorem identifying every Gaussian selected
image with a regular compact period. No such transfer is used here.
The original regular-transcendental realization/classification question
remains open beyond the filtered class; polynomial clocks and their
monomial subfamily are now understood transports, not a solution of it.

Earlier formal/all-ring,ordinary/numerical,family,origin and filtered
results stay closed. All Oaku operators,torsion witnesses and thresholds
remain sealed;Sabbah/Maisonobe retain the convergent rational-slope credit.
Nonconvergent formal/all-ring priority is unassessed. No capped-job restart,
delegation,outreach,assurance promotion or broad-goal completion.
