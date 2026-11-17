# Complete finite actual images with quadratic W and mixed WY terms

Use one Witt Gaussian pair Z,W with E[ZW]=1 and an independent ordinary
Gaussian Y with E[Y^2]=1. All other basic pairings are zero. Over a
characteristic-zero field K consider

    P=ell(Z)W^2+W(h(Z)+c(Z)Y)+aY^2+u(Z)Y+k(Z),

where a is constant and ell,h,c,u,k are polynomials. Put

    Phi_P[Q]=E[Q exp(tP)], V_P=Phi_P[K[Z,W,Y]],
    M_P=K[t]V_P.

The ORIGINAL fixed spatial markers define V_P; M_P allows polynomial
time coefficients. A phase is safe if all those original marked periods
are polynomial, and pure if Phi_P[1]=1.

**Complete mixed-source theorem.** Assume

    c!=0, E=max(deg ell,2deg c)>=2,                 (Q1)

with deg0=-infinity. There are no root-multiplicity or degree bounds
beyond (Q1). Then

    M_P is finite over K[t]
      iff a=0, ell belongs to Z^3K[Z],
               h,c belong to Z^2K[Z].             (Q2)

The polynomials u,k are arbitrary. In every finite case

    psi=k(0)t+u(0)^2t^2/2,
    Phi_P[1]=exp(psi), V_P=M_P=exp(psi)K[t].         (Q3)

Thus the vacuum is primitive. Safety is exactly (Q2) together with
u(0)=k(0)=0; every safe member is pure. Every NONCONSTANT c satisfies
(Q1). If c is constant nonzero, (Q1) requires deg ell>=2, and (Q2) is
impossible, so no finite image exists. The lower-E branches are outside
this theorem.

This is a FULL-actual-image classification. The new weighted-marker
argument does not assert that an exponential-polynomial scalar unit,
or purity alone, forces (Q2). Independent mathematical review and formal
verification are not claimed.

## 1. The scalar normalization and the retained curve system

The
[general finite actual-image theorem](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#7-every-finite-actual-image-and-its-original-fixed-marker-image)
gives, under finite M_P,

    M_P=exp(psi)K[t], psi in tK[t].                 (Q4)

It does NOT initially give primitive vacuum. In particular EVERY original
polynomial marker rho with Phi_P[rho]!=0 supplies a nonzero polynomial

    p_rho=exp(-psi)Phi_P[rho] in K[t].              (Q5)

Work first over C. If ell=0, (Q1) makes c nonconstant, and the
[complete W-affine theorem, Section14](GAUSSIAN_W_AFFINE_FINITE_IMAGE_OBSTRUCTIONS.md#14-a-fixed-original-logarithmic-marker-removes-every-root-order-restriction)
already proves (Q2)-(Q3) in that branch. Hence suppose ell!=0.
Throughout the complex curve argument write K_t=C(t) for the rational
time field; K in the theorem denotes the constant coefficient field.

Put d=1-2at, r=1/t and

    A=t ell+t^2c^2/(2d), B=th+t^2cu/d,
    Psi0=tk+t^2u^2/(2d), D=exp(Psi0),
    q=c^2+2d ell/t=G+2r ell, G=c^2-4aell.

The normalized curve and its phase and seed are

    U^2=q, U!=0,
    phi=Psi0-(z-B)^2/(4A), eta=r dz/U.              (Q6)

Gaussian W,Y elimination identifies its odd sector with the determinant
half-density of the localized original source. The precision determinant
is -t^2q, so its normalized volume differs from eta only by a fixed
nonzero complex constant. The original-source or scalar bridges below
fix that constant by the actual marked period.

The complete comparison inputs are those retained in the
[varying-quadratic-W proof](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#a-varying-quadratic-w-coefficient-clean-localization),
especially its
[finite-centre estimates](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#positive-volume-order-over-finite-centres)
and
[infinity and complete period comparison](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#infinity-and-complete-period-comparison).
They keep the full marked semistable curve, positive-genus cohomology,
annulus complexes and gluing maps, descending tails, and the ordinary-
parameter perfect period pairing. The proof here supplies the needed
source in every case and removes the old multiplicity restriction on
the volume. Its external comparison dependencies and pending independent-
review status remain unchanged.

## 2. An original polynomial marker reduces every common multiplicity

Take the monic gcd and its squarefree decomposition over K:

    g=gcd(ell,c^2)=product_(m>=1) g_m^m,
    rho=product_(m>=1) g_m^floor((m-1)/2).           (Q7)

The factors g_m are pairwise coprime and squarefree. Thus rho is a
fixed ORIGINAL polynomial marker over K, including when its roots are
not individually K-rational. At a common root of multiplicity m, the
rho exponent is floor((m-1)/2).

One has rho^2|ell, rho^2|c^2 and rho|c. Set

    ell*=ell/rho^2, c*=c/rho, U*=U/rho.
    q*=q/rho^2=c*^2+2(r-2a)ell*.

Then

    rho eta=r dz/U*,
    gcd(ell*,c*^2)=g/rho^2,                         (Q8)

whose root multiplicities are exactly ONE for odd m and TWO for even m.
No common-root support disappears. On the U!=0 open curve rho is
invertible. Conversely every rho root remains a root of g/rho^2, so
rho is invertible on the U*!=0 open curve as well. Hence this is an
isomorphism of the punctured normalized covers, not an extra untracked
localization of the source.

If a c-root has multiplicity n, its rho exponent is at most n-1.
Consequently

    rho divides gcd(c,dc/dz).

If c is nonconstant, c* is still nonconstant. If c is constant nonzero,
rho=1 and deg ell*=deg ell>=2 by (Q1). Therefore in every case

    E*=max(deg ell*,2deg c*)>=2.                    (Q9)

When G*!=0, the retained finite-centre volume inequalities apply to (Q8): their
hypothesis is precisely common multiplicity at most two, and their
finite-centre proof imposes no lower bound on deg ell*. At a finite
centre, if v(z-alpha)=N>0, l=ord ell*, j=ord G*, the minimum of j,l
is at most two. G*-dominance gives
1+(1-j/2)N; if j>2 then l<=2 and N<1/(j-l), giving the positive
endpoint at least1/2. r ell*-dominance gives
1/2+(1-l/2)N>=1/2. At cancellation the leading coefficient j-l
prevents derivative cancellation and

    rho eta=2r dU*/(G*'+2r ell*')

gives strictly more than the same positive endpoint. Ordinary finite
centres have order at least one. These formulas include minimum zero.

At infinity, if deg ell*>=2, use the retained complete infinity
inequalities. If deg ell*<=1, c* is nonconstant and

    deg G*=2deg c*>=2>deg ell*.

The G* term dominates at every infinity scale and gives order
1+(deg c*-1)N>=1. Thus the possible drop in ell-degree during (Q8)
does not create an omitted chart.

If G*=0 identically, ell* is proportional to c*^2 and has degree at
least two. The cover is fixed, with scalar normalization proportional
to sqrt(r). Its differential is logarithmic at every double root and
regular at every simple root, and logarithmic or holomorphic at infinity.
The same argument applies componentwise if the geometric normalization
disconnects.

It follows that rho eta has strictly positive order on EVERY generic
component of EVERY sufficiently refined marked model, independently of
the phase. Its generic horizontal poles are logarithmic or holomorphic.
All these statements survive a scalar polynomial phase shift psi.

## 3. Nonzero actual marker and the weighted finite-image contradiction

We will construct a horizontal functional on the full original curve
whose value on rho eta is precisely Phi_P[rho]. It may be an alternate
scalar/base-Z bridge; no original W-marker identification is required
in that case.

Outside h(0)=c(0)=0, this marked period is NONZERO. If rho(0)!=0 it has
that nonzero initial value. Otherwise write rho=rho_k z^k+..., k>0.
A positive rho exponent means that the origin multiplicity of g is at
least three. Thus ord0 ell>=3 and ord0 c>=2, and the displayed outside-
origin condition forces h(0)!=0.

All monomials of P except h(0)W then have nonnegative Z-minus-W angular
weight. No power smaller than k can balance rho, and in power k the
only balancing term is rho_k Z^k(h(0)W)^k. Therefore exactly

    Phi_P[rho]=rho_k h(0)^k t^k+O(t^(k+1))!=0.       (Q10)

At a common origin with ord0 ell=1, the gcd multiplicity there is at
most one, so rho(0)!=0 again.

Under finite M_P, (Q5) therefore holds wherever these conditions apply.
Gauge to phase phi-psi and use the seed rho eta/p_rho. Its actual
functional value is one. Even if p_rho(0)=0, its inverse is a valid
rational time coefficient with a finite Laurent principal part. It
changes no spatial pole and increases every infinity order by deg p_rho.

Horizontality on the entire curve module implies, on the differential
span of this seed,

    L'(nabla_t^n(rho eta/p_rho))=partial_t^n(1).

Its restriction is thus a genuine rational horizontal quotient onto
K_t, not merely a scalar formal coincidence. Section2 and the retained
complete component, annulus, tail and gluing comparison exclude this
unit quotient. This is the contradiction.

The WHOLE finite-image premise is essential to obtain the polynomial
p_rho. A scalar identity for Phi_P[1] alone supplies no such information
about the marked period. All source cases below are used only to
justify the actual identity on rho eta and its covariant derivatives.

## 4. Full source outside the exceptional moving residue

Suppose first ell(0)!=0. The
[completed-z original-source theorem](GAUSSIAN_UNIT_ELL_FORMAL_SOURCE_FINITE_IMAGE_EXCLUSION.md#1-completed-z-wick-contraction-preserves-the-full-original-source)
extends the ORIGINAL functional after inverting q. Its construction
has no multiplicity premise: q^-1 has Laurent-t coefficients in C[[z]],
and circular contraction uses only the finite z Taylor jet selected by
each finite W degree. All spatial Stein identities and the parameter
connection remain valid. Gaussian fibre reduction then retains every
original polynomial marker, in particular rho. Since g divides ell,
ell(0)!=0 also gives rho(0)!=0. Thus Section3 applies.

Now assume ell(0)=0 and (h(0),c(0))!=(0,0). Define

    e=ell/g, C=c^2/g, f=4ae-C,
    J0=-2e+tf, Rmov=C[zf-2eh]+4e^2cu.              (Q11)

Then gcd(e,C)=1. The determinant is t g J0; its fixed and moving
spatial divisors are disjoint over K_t. At a fixed g-root alpha, J0(alpha)
could vanish identically only if e(alpha)=C(alpha)=0, which is impossible.

Every fixed root is clean. Its restricted full phase has W-source

    b_alpha=t h(alpha)-alpha,

nonzero for alpha!=0. If alpha=0 is a g-root, then c(0)=0 and the
outside-origin condition gives h(0)!=0, so b_alpha is again nonzero.
On that plane the operator partial_W+b_alpha is bijective on polynomials,
with its finite geometric derivative inverse. Its full restricted
exponential de Rham complex is therefore acyclic.

When J0 has positive spatial degree, it is irreducible and separable
over K_t: it is primitive and t-linear with gcd(e,f)=1, and e/f is not
constant. On its finite residual field, complete the nondegenerate
Y-square. Its remaining W-linear source is a nonzero scalar multiple
of S=d(th-z)+t^2cu. Substitution t=2e/f gives

    f^2 S=Rmov.

A nonzero fixed polynomial Rmov cannot vanish on the moving divisor:
a fixed root would force a common zero of e,f. Thus if Rmov!=0 that
source is a unit at every conjugate moving root, and the full moving
restriction is acyclic by the same W-linear contraction.

The retained
[source-localization proof](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#exact-source-criterion-and-localization-of-the-original-seed)
now transports the original source to the curve. Its smooth-divisor
residue sequence and higher-pole comparison give the localization
isomorphism; the canonical map commutes with the parameter connection.
If J0 is spatially constant, no moving divisor or Rmov condition is
present. This source proof is independent of the multiplicity of g:
one localizes at its reduced support.

If Rmov=0, reduction modulo e gives

    Rmov=-z C^2 modulo e,

so coprimality forces e|z. Only e constant or e=kappa z remains.
The latter is Section5. The constant case has a direct two-stage
source proof, with no moving-nullsource assumption.

Indeed for e constant nonzero,

    J0^-1=(-2e)^-1
        *sum_(j>=0)(t(4ae-C)/(2e))^j

has POLYNOMIAL-z t coefficients. First extend the actual original
functional by this series; its spatial and parameter identities remain
coefficientwise valid. Then invert g. At every fixed root, J0 restricts
to a nonzero scalar in K_t, leaving the entire W,Y plane and the same
nonzero source b_alpha. The restriction is acyclic by partial_W+b_alpha,
so fixed-divisor localization is a cohomology isomorphism. The regular
connection pole filtration contracts higher poles by their nonzero
normal residue, and the first residue term is the acyclic restricted
complex. The canonical localization map is horizontal even if that
intermediate filtration is not parameter-stable.

This transports the full original source to the determinant-localized
curve. It proves the needed rho identity without assigning an arbitrary
original value to 1/g. Hence Section3 applies in the constant-e case too.

## 5. The moving-common-root and fixed-simple-origin scalar bridges

For completeness the remaining source constructions retain only the
needed ORIGINAL base-Z marker and all its parameter-covariant iterates.

**Moving residual e=kappa z, Rmov=0.** Coprimality gives C(0)!=0.
The coefficient of z in (Q11) forces

    h(0)=-C(0)/(2kappa)!=0.

The equation 2de+tC=0 has an algebraic germ beta in tC[[t]], with
beta=t h(0)+O(t^2). At this root the exact identity is

    t^2Rmov(beta)=-4d e(beta)^2[beta-B(beta)].

Both e(beta) and g(beta) are nonzero in C((t)), so A(beta)=0 and
beta-B(beta)=0. Also A_z(beta)!=0. If gamma=ord0g, then c^2=gC
makes gamma even and

    A_z(beta)=kappa g_gamma h(0)^gamma
                t^(gamma+1)+O(t^(gamma+2)).         (Q12)

Thus the
[moving-common-root bridge](GAUSSIAN_MOVING_COMMON_ROOT_SCALAR_EXCLUSION.md)
applies for ARBITRARY gamma, not only zero or two. To recall its exact
mechanism, translate z=beta+x and factor

    A=xA1, z-B=xF, H=2A1/F^2, F=1+O(t).

The actual circular Wick-residue identity is translation-invariant
because beta is t-divisible. It gives the root x=H S^2 and

    Phi_P[rho]=d^(-1/2)
        E_S[rho(beta+x)D(beta+x)/(FJ)],
    J=1-H_x S^2.                                   (Q13)

The same formula with rho=1 is the retained actual scalar identity;
linearity of its residue calculation proves (Q13). If A_z(beta) starts
a0t^nu, then x starts 2a0t^nu S^2 and

    U starts 2a0t^(nu-1)S.

Every allowed curve function and every covariant iterate therefore
evaluates in the Laurent-t ring with coefficients in C[S,S^-1]. The fixed Laurent Gaussian prescription
G(S^(-2n))=(-1)^n/(2n-1)!!, with odd moments zero, kills partial_S-S
at every integer exponent. The resulting horizontal curve functional
has value (Q13) on rho eta. Its moving-coordinate chain rule includes
beta'(t). No original W/Gysin-marker identification is asserted.
The construction requires no bound on nu or gamma; Section2 supplies
the NEW weighted volume bound. Thus Section3 applies.

**Common origin with ord0 ell=1.** Write ell=zL, h=zH0,c=zC0,L(0)!=0.
The
[fixed simple-origin construction](GAUSSIAN_SIMPLE_ORIGIN_ALTERNATE_SCALAR_EXCLUSION.md)
is the same residue argument with beta=0:
A=t z Ltilde, z-B=zF and H=2tLtilde/F^2 give x=z=H S^2.
Here U starts 2L(0)S. Every curve function and finite covariant iterate
again lies in the fixed Laurent-S ring, and its prescribed Gaussian
functional has the exact actual base-Z value

    Phi_P[rho]=d^(-1/2)E_S[rho(z)D(z)/(FJ)].

There is no multiplicity premise in this SOURCE construction. Since
rho(0)!=0, the period is nonzero and Section3 supplies the weighted
contradiction.

**Spatial nonconstancy in every comparison case.** If ell is not
proportional to z^2, the t^-1 coefficient -z^2/(4ell) of phi is
nonconstant. If ell=kappa z^2, all comparison cases above lie outside
h(0)=c(0)=0. The coefficient of t^0 is

    h/(2kappa z)+c^2/(8kappa^2 z^2).

If c(0)!=0 it has an unavoidable double pole; if c(0)=0,h(0)!=0 it
has an unavoidable simple pole. Thus phi-psi is spatially nonconstant.
The common-origin order-two case instead uses the angular argument
below, not this comparison assertion.

The preceding cases are exhaustive outside the remaining common origin:
ell(0)!=0 is Section4; otherwise fixed roots are clean outside the
common origin and Rmov is either nonzero, or e is constant or kappa z.
At a common origin the raw ell order is one or at least two.

## 6. The common-origin angular core and all-marker sufficiency

It remains to treat ell(0)=h(0)=c(0)=0 with ord0 ell>=2. Set S=ZW and
split P into its zero angular-weight part and its positive remainder:

    Pzero=ell_2 S^2+h'(0)S+c'(0)SY
              +aY^2+u(0)Y+k(0),
    P=Pzero+Pplus.

Every monomial of Pplus has strictly positive Z-minus-W weight.
The
[scalar quadratic angular-core theorem](GAUSSIAN_MULTIBASE_QUADRATIC_FIBRE_PRIMITIVITY.md#10-a-quadratic-zero-weight-part-with-arbitrary-positive-remainder)
therefore says that finite M_P forces

    ell_2=h'(0)=c'(0)=a=0.

In its scalar ordinary fibre, the null-chain condition c'(0)^2=0
forces c'(0)=0 over a field. Thus precisely (Q2) is necessary. No
root-multiplicity condition is used here.

Conversely (Q2) gives

    P=u(0)Y+k(0)+Pplus

with strictly positive angular Pplus. For an ORIGINAL marker Q of
degree D, angular balance truncates the expansion after at most D
positive factors. Ordinary Gaussian translation gives exactly

    Phi_P[Q]=exp(psi)
       *sum_(n=0)^D t^n/n!
        E[Q(Z,W,Y+t u(0))
          Pplus(Z,W,Y+t u(0))^n].                  (Q14)

The sum is polynomial. Its unit is exactly exp(psi), since every n>=1
term for Q=1 has positive angular weight. Hence M_P=exp(psi)K[t]
with primitive vacuum. The phase is nonconstant because c!=0, so the
retained original-image theorem gives V_P=M_P. This proves (Q2)-(Q3).

In every finite member Pplus has Y degree at most one. For a marker
monomial Z^iW^jY^alpha, angular balance gives n<=j-i and translation
in (Q14) gives time degree at most alpha+2n. Thus the normalized
original period has degree at most 2deg Q; when u(0)=0, the bound is
deg Q because there is no Y translation. These are sufficient bounds
for the whole family, not a general exact filtered-image assertion.

Safety forces exp(psi) to be polynomial, hence u(0)=k(0)=0.
Conversely those conditions make every marked period polynomial by
(Q14). In that case P has strictly positive angular support, and

    E[Q P^n]=0 for n>deg Q.                         (Q15)

Thus every safe member has original mixed cutoff C=1. The
[cutoff-to-jet theorem](MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md#1-the-finite-jet-theorem)
then gives moment r-jets through (d+1)r for perturbations of original
spatial degree at most d. This concerns jets at a known safe phase,
not a universal arbitrary-point finite purity test.

Over any characteristic-zero field, the squarefree factors in (Q7)
and all finite-image data are defined over that field. A hypothetical
finite image gives psi and p_rho there; embed the finitely generated
field of these coefficients and the phase into C. The nonzero marked
leading term, source relations and gcd data survive. This descends
necessity. Formula (Q14) is algebraic and proves sufficiency directly.

## 7. An exact degree-five filtered image and sharp moment-jet budget

Consider the explicit phase

    Pdag=Z+Z^2W+Z^3W^2+Z^2WY
        =Z(1+S+S^2+SY), S=ZW.                      (Q16)

Every monomial has angular weight exactly one. Thus a marker of degree
at most D has polynomial period of degree at most D. For W^n, only
the moment with power n can survive, and

    Phi_Pdag[W^n]=A_n t^n,
    A_n=E[S^n(1+S+S^2+SY)^n]/n!.                   (Q17)

The rational integer A_n is positive: all surviving circular and ordinary
Gaussian contractions have nonnegative coefficients, and the constant
term of the parenthesized polynomial contributes n!/n!=1. In particular
A_n is nonzero over every characteristic-zero field. Consequently

    {Phi_Pdag[Q]:deg Q<=D}=K[t]_(<=D)               (Q18)

for EVERY D. The original marker A_n^(-1)W^n lifts t^n and has the least
possible degree n. Direct initial diagonal values are

    A_0=1, A_1=9, A_2=535, A_3=89593, A_4=29527341.

There is an exact finite formula for EVERY original marker, and hence
for the full original kernel. Write R(S,Y)=1+S+S^2+SY and use the
product functional E_S,Y[S^jY^k]=j! E[Y^k]. For a monomial,

    Phi_Pdag[Z^iW^jY^k]=0                          if j<i,
    Phi_Pdag[Z^iW^jY^k]
      =t^(j-i)/(j-i)! E_S,Y[S^jY^k R^(j-i)]       if j>=i.  (Q19)

For Q=sum q_(i,j,k)Z^iW^jY^k, its period is zero exactly when

    sum_(i,k) q_(i,i+n,k) E_S,Y[S^(i+n)Y^k R^n]=0
                           for every n>=0.         (Q20)

Only finitely many such equations occur for each Q. This is a complete
linear-space kernel criterion in the ORIGINAL coordinates. It is not
an ordinary ideal: B=2W-9WY has zero period, whereas
Phi_Pdag[WB]=-280t^2.

For degree-at-most-d perturbations, the general jet bound (d+1)r is
sharp for the full perturbation space. Indeed take H=vW^d and moment
m=(d+1)r. Its coefficient of v^r is

    binom(m,r) E[W^(dr)Pdag^(dr)]!=0.              (Q21)

For any j<r the corresponding monomial has total angular weight
m-(d+1)j=(d+1)(r-j)>0, so these lower perturbation coefficients vanish.
Hence this moment has perturbation order exactly r, for every r>=1
and d>=0. The cutoff concerns individual moment jets; no least ideal
prefix is asserted. An actual coefficient nilpotence certificate of
index nu leaves only intervening moments through (d+1)(nu-1) to check.
In fact angular balance makes the entire one-parameter moment zero
unless (d+1) divides m; at m=(d+1)r it is exactly the nonzero coefficient
in (Q21) times v^r. This does not make the jet budget a least ideal prefix.

## 8. A no-finite-image corollary including c=0

For EVERY c, including c=0, suppose deg ell>=2 and every root of
gcd(ell,c^2) has multiplicity at most two. Then M_P is not finite.

For c!=0, this follows from (Q2): finiteness would force ord0 ell>=3
and ord0 c>=2, so their gcd would have origin multiplicity at least
three. For c=0, g is proportional to ell and the unweighted volume
already satisfies the retained multiplicity-at-most-two bounds.
If ell(0)!=0 use the full formal source. If ell(0)=0,h(0)!=0 use the
two-stage constant-e proof, now with C=0 and scalar J0. If ell(0)=h(0)=0,
raw order one uses the simple-origin bridge; raw order at least two
uses the angular core, which forces order at least three and contradicts
the hypothesis on g.

This is a NO-FINITE-IMAGE corollary. In particular there are no safe
members, but absence of a finite full image does not by itself imply
impurity or classify the scalar unit.

## 9. Scope and assurance

The proof combines exact original-source extensions, explicitly tracked
base-Z scalar bridges, the polynomial logarithmic marker, retained
complete curve comparison, and the known angular-core finite criterion.
The two-stage constant-e argument is proved here; no separate unproved
moving-source exception is hidden in the case split.

The multiplicity restriction is removed by the ORIGINAL marker rho,
not by extrapolating an older pure-unit multiplicity bound. The
nonzero polynomial p_rho uses the WHOLE actual image. Scalar-unit-only
classification, purity implies safety outside previously proved
subclasses, and universal Gaussian safety implies purity are not
asserted. The present theorem proves safety implies purity within
(Q1), with primitive vacuum for every finite actual image there.

All unrestricted steps are written arguments. Any finite exact checker
has only its recorded bounded algebraic scope. Imported comparison
statements retain their existing primary-source audit and pending
independent-review status; no external review, formal verification,
publication or priority claim is added by this note.

The [focused exact checker](verify_gaussian_quadratic_w_mixed_finite_images.py)
and [receipt](gaussian_quadratic_w_mixed_finite_images_verification.json)
retain explicit high-multiplicity cover reductions, the constant-e
geometric inverse and a zero original escaping-kernel period, nonzero
marked periods with zero initial mass, a weighted moving-source identity,
original finite-branch Wick controls through moment six, and the first
five diagonal entries in (Q17). These are bounded exact controls; the
all-degree filtration, original kernel and sharp jet budget follow from
(Q16)-(Q21).
