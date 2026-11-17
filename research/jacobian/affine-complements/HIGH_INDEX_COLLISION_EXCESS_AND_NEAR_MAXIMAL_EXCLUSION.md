# Image-mass collision, high-index excess and the complete near-maximal exclusion

PF2RULED112, 16 November 2026. Keep the ORIGINAL finite normalization
pi:N->A2 of a nonautomorphic complex plane Keller map of generic field
degree d, its WHOLE marked U=A2, and ALL r omitted primes. N may be
singular. No knot, ruling, Laurent support, residue or finite-image
singularity hypothesis is imposed in the main theorem.

**All-degree image-mass theorem.** For a nonproper image component C,
let w_C=sum_(D maps to C)(e_D*f_D), and let b_C be the number of those
primes with residue degree one. If 2w_C>d, then the normalization of
C is injective and C is homeomorphic to A1. At least one omitted prime
has a DIFFERENT image, each of the b_C residue-one primes has positive
own finite excess, and

    w_C+b_C <= d-2.                                  (0)

If the WHOLE nonproper-value set is irreducible, its TOTAL generic
omitted mass and residue-normalized excess satisfy

    sum_D e_D*f_D <= floor(d/2),
    sum_(D,a) eta_a >= d-1-floor(d/2).                (0a)

This covers arbitrary residue degrees and shared image components.
It does not infer irreducibility of the nonproper set from that of the
finite branch curve.

**All-degree high-index theorem.** If an omitted prime R has transverse
index e with 2e>d, then:

* its residue degree is one and its image normalization is injective;
  the image C_R is homeomorphic to the complex line;
* its own finite excess is POSITIVE: E_R=sum_a(ell_a-e)>=1;
* another omitted prime has a DIFFERENT image component. In particular
  the whole nonproper-value set is reducible and r>=2;
* the exact global budget satisfies e+r<=d-1.

Consequently, in EVERY generic degree d>=5, every omitted index satisfies

    e <= d-3.                                        (1)

This closes the ENTIRE e=d-2 boundary class, in ALL degrees d>=5,
including singular normalizations and arbitrary other omitted primes.
The e=d-1 exclusion was already imported from Orevkov in PF46; it is
not claimed as a new exclusion. The inference (1) broadens the retained
Chau/Orevkov inputs beyond PF32's particular index-two Wright model.
No historical priority or independent proof of those inputs is claimed.

If the WHOLE nonproper-value set is irreducible, the stronger bound is

    e <= floor(d/2) for EVERY omitted prime.          (2)

This includes any sole-omitted-prime source, smooth or singular. It is
the irreducibility of the NONPROPER-VALUE set, not merely of the finite
branch curve, that is required in (2).

**The next boundary strip.** For d>=7, a prime with e=d-3 forces exactly
two omitted primes: R with (e,f)=(d-3,1), and one unramified D with
(e,f)=(1,1). There is exactly one finite excess unit, on R. Its image
C_R has injective normalization and exactly one nonimmersive parameter;
the other image C_D is different, is everywhere immersively normalized,
and identifies at least two distinct normalization parameters. This is
a complete necessary profile, NOT an exclusion or a realization of it.

**Quotient corollary.** Combining (1) with PF106, every finite cyclic
quotient singularity of an actual normalization in degree d>=5 has
order n<=d-4. Thus degree-five actual normalizations are smooth;
degree-six ones can have only A1 quotient singularities. Neither all
degree-five nor all degree-six Keller maps are excluded here.

## 1. Retained global inputs and their exact scopes

The finite-flat boundary ledger of
[PF2OMIT6](../jc2-omitted-divisors/RESIDUE_LEDGER.md) holds for ALL actual
normal N. Each omitted prime has unibranch support and normalization
A1. Its image has normalization A1. For its finite residue polynomial
of degree f, write q(a) for local degree, ell_a for ambient finite fibre
length, and eta_a=ell_a-e*q(a)>=0. Then

    d-1=sum_D e_D*f_D + sum_(D,a) eta_a,
    E_D=e_D*(f_D-1)+sum_a eta_a.                      (3)

Orevkov, *On three-sheeted polynomial mappings of C2*, Lemma5.2,
has GENERAL quantifiers, despite the title: at a finite actual boundary
point with ell_a=e_D, the induced image-normalization germ is a
nonsingular embedding. The source comparison with the original normal
finite Stein factor is the one already audited in PF105/PF2OMIT6.
No assumption that a generic index determines every local length is
made. The lemma is imported, not independently reproved here.

Use TWO separate statements from the PUBLISHED Chau2003 paper,
*Two remarks on non-zero constant Jacobian polynomial maps of C2*:

1. The Remark on printed43/fullPDF5 says the exceptional-value set of
   a nonsingular polynomial map has NO irreducible component isomorphic
   to the affine line. Other components may be present.
2. Corollary3 on printed40/fullPDF2 says the WHOLE exceptional-value
   set of a Keller map cannot be a simply connected curve.

For a Keller map the critical-value set is empty, so this exceptional
set is exactly the nonproper-value set S_F. The second statement is
NOT a prohibition on a contractible component of a larger union.
The first applies only to an algebraically smooth embedded A1
component; it does not itself exclude a cuspidal contractible component.
These distinctions are essential to the two different deductions below.
Both statements and their proof context are explicit imports, already
retained in the repository; all six published PDF pages were freshly
rendered and visually inspected. The withdrawn stronger Jelonek2020
claim is not an input. No new C* action or target-field descent is assumed.

Finally the ORIGINAL finite normalization gives exactly

    S_F = union_D pi(D).                              (4)

The retained boundary theorem gives N minus U as this union of omitted
curves; no isolated missing component is discarded.

To check both directions, a sequence escaping U whose images stay
bounded has a cluster point in the proper finite preimage in N; this
point lies outside U. Conversely a point of an omitted divisor is
approached by U, and such a sequence cannot stay bounded in the actual
source A2, or its limit would also belong to U. Thus every image point
of an omitted divisor is nonproper. No component-open deletion retains
this proof or its original finite-map budget automatically.

## 2. High index gives injectivity, then two independent costs

Fix R with 2e>d. Equation (3) implies e*f<=d-1<2e, hence f=1.
The residue map R^nu=A1->C_R^nu=A1 is therefore a degree-one polynomial
and is an isomorphism. Write gamma:A1->C_R for the finite normalization
parameterization. If two distinct parameters mapped to the same target
point, they would correspond to distinct points of R, by unibranchness.
Each ambient local fibre length is at least e. Their sum would be at
least 2e>d, contradicting finite flat degree d. Therefore gamma is
INJECTIVE, without an immersion or zero-excess hypothesis.

Since gamma is finite it is proper. A proper continuous bijection from
A1(C) onto C_R(C) is a homeomorphism. Thus C_R is topologically a line,
even though it can have a cusp and need not be isomorphic to A1.

Suppose E_R=0. Because f=1, every eta_a=0 and ell_a=e at every point
of R. Orevkov Lemma5.2 makes gamma immersive everywhere. Its already
proved global injectivity now makes C_R a smooth embedded affine line.
But C_R is an irreducible component of S_F in (4), contradicting
Chau's no-A1-component Remark. Consequently

    E_R=sum_a eta_a >= 1.                            (5)

This uses no assertion that the whole branch locus is C_R. Other
ramified primes and all unramified omitted images remain permitted.

Separately, if every omitted prime had image C_R, (4) would give
S_F=C_R, homeomorphic to the simply connected complex line. Chau's
Corollary3 forbids this WHOLE set. Hence there is another omitted
prime with a different image component, and r>=2. This proves (2).
In particular an irreducible nonproper-value set has no prime of
index strictly greater than half the field degree.

Every other prime costs at least one in (3). Add these r-1 costs to
the cost e and the excess unit (5):

    d-1 >= e+(r-1)+1 = e+r,  r>=2.                  (6)

This is the all-degree theorem. If d>=5 and e>=d-2, then 2e>d, so (6)
would force d-1>=e+2>=d, a contradiction. Values e>=d were already
forbidden by (3), and e=d-1 was Orevkov's retained maximal-index gate.
Thus every actual omitted index satisfies (1), with no bound on the
number of other primes or any finite image singularities. Small d=4
is not inserted into the strict-half argument at e=2; the published
Domrina quartic exclusion remains the separate retained theorem.

## 3. The total mass above one image component

Let C be any component of S_F, and sum over ALL omitted primes whose
image is C. The retained boundary theorem makes different omitted
primes pairwise disjoint in the ORIGINAL N; each is unibranch. For any
parameter c of C^nu=A1, the finite polynomial D^nu->C^nu has total
local degree f_D above c. By ell_a>=e_D*q_D(a), the omitted points
over this one normalized target branch consume ambient fibre length
at least w_C. If two distinct normalization parameters represented
the same target point, the corresponding sets of omitted points would
be disjoint and would consume at least 2w_C. Thus 2w_C>d forces the
normalization C^nu->C to be injective. Its finiteness again makes C
homeomorphic to A1. This argument does not require f_D=1.

Chau's Corollary3 now requires a different image component in S_F,
and hence an outside omitted prime costing at least one. Chau's
no-A1-component Remark separately says C cannot be algebraically
smooth. Choose a singular parameter of its injective normalization.
For EVERY residue-one prime above C, its normalized map to C^nu is an
isomorphism, so it has a distinct omitted point above this parameter.
Orevkov's Lemma5.2 prevents ell=e at that point. Since q=1 there,
eta>=1 on each of these b_C primes. Distinct primes have separate
ambient points, so the costs add. The global ledger gives

    d-1 >= w_C+b_C+1,

which is (0). A known SUBLIST of same-image primes can already supply
a necessary version of this filter: a lower bound on w_C exceeding
d/2 proves injectivity, and all its known residue-one primes still
cost a unit each. No unknown prime is deleted from the finite map.

If S_F is irreducible, every omitted prime maps onto this C. Then
2w_C>d is impossible, since it would make the WHOLE S_F simply
connected. Hence w_C<=floor(d/2), and subtracting from (3) proves
(0a). The excess here is sum eta, not sum E_D: the latter also includes
the residue ramification terms e_D*(f_D-1).

An image component carrying an UNRAMIFIED omitted prime also has
w_C<=floor(d/2), even when the whole S_F is reducible. PF2OMIT2 gives
f=1 and zero own excess everywhere on such a prime; the positive unit
forced by 2w_C>d would contradict that retained fact.

For example the same-image numerical list (e,f)=(3,1),(3,1),(2,1)
at d=12 has no individually high index, but w=8>d/2 and b=3.
It fails w+b<=10. The entire such common-image class is excluded by
the written theorem, irrespective of its target equation or knot.
Conversely w=7,b=3,d=12 meets that necessary budget; no realization
is inferred. For smooth actual N, PF105 already imports residue-one
for every omitted prime, so b_C is just their number above C.

## 4. The e=d-3 profile and quotient-order consequences

Let d>=7 and e=d-3. Then 2e>d. The remaining budget after paying R
is two. By Section2 it must pay for both positive excess on R and
another prime with a distinct image. Both costs are positive integers,
so equality forces exactly one other prime D with e_D*f_D=1, exactly
one eta_R=1, and no other eta. Thus e_D=f_D=f_R=1 and

    ell_R(a*)=e+1,   ell_R(a)=e for a!=a*,
    ell_D(b)=1 for every b.                          (7)

Orevkov's lemma makes C_R immersive away from a*. Since gamma_R is
injective and its image cannot be a smooth A1 component, a* is indeed
nonimmersive, and is its only such parameter. The same lemma makes
gamma_D immersive everywhere. The two image curves are different;
if gamma_D were also injective, its proper immersed image would be
an A1 component, again forbidden by Chau. Therefore gamma_D must
identify distinct parameters. No transversality of those collisions,
source smoothness at a*, target equations, or finite-cover realization
of this profile is asserted.

At every finite singular point of an actual normalization, PF106
gives marked cyclic quotient type 1/n(1,-e), with gcd(n,e)=1 and
2<=n<=e-1. Combining with (1) yields n<=d-4. For d=5 no n>=2 remains,
so N is smooth. For d=6 a singular point has n=2. Its marked type is
1/2(1,1), an A1 singularity. The primary quotient and character results
are credited to [PF106](CYCLIC_QUOTIENT_ORDER_AND_MORSE_RESIDUE_BUDGET.md),
not independently reclassified here. Lower-degree inverses are NOT
inferred from smoothness or these necessary boundary restrictions.

## 5. Stronger use of a partial linking computation

For any actual end link containing K_R, let S be its core subgroup in
H_1(W), and let g be the minimum number of generators of H_1(W)/S.
[PF111](WHOLE_PLANE_THOM_LATTICE_AND_COMPLETE_PERIPHERAL_OBSTRUCTION.md)
forces at least g omitted primes outside that link. None is R.
Section2 also forces at least one other omitted prime regardless of
the link choice. Thus r-1>=max(1,g), and (5) strengthens the budget to

    d-1 >= e+max(1,g)+1,
    equivalently d>=e+max(1,g)+2.                    (8)

The positive excess is an independent cost; it cannot be paid by
declaring an additional unramified omitted prime. An integral linking
matrix or zero individual row defects still does not certify that
S=H_1(W); PF111's complete generation hypothesis is retained.
This is an all-degree filter, not a realization criterion.

## 6. Seventeen sheets: application and exact comparison control

The exploratory control from PF111 is now settled by the broader
theorem (1). Work over F16=F2[z]/(z^4+z+1), with primitive ell=z,
t=ell+ell^(-1)+1. The matrices

    A=[[1,t],[0,1]],  B=[[0,1],[1,1]]

act on P1(F16). Their image has order4080, namely PSL_2(F16), and
the trefoil relation A^2=B^3=1 holds. The actual meridian action has
cycle type (15,1,1), and lambda=mu^(-6). Boundary orbit types are
(e,f,a)=(15,1,6),(1,1,0),(1,1,0). The exact seed has

    H_1(W)=Z/51,     H_1(E)=Z^3 plus Z/3,
    F=(6/17)*J_3+diag(0,-6,-6),
    every core order=17, every row denominator=17,
    every individual cleared-row defect=1,
    H_1(W)/<all cores>=Z/3.                          (9)

Thus it is a concrete genuine cover whose individual defects vanish
while the full linking matrix misses a Z/3 quotient. It does not
contradict PF111's conditional completeness theorem.

The PF111 b=2, M=15 sector has gcd(15,2-6)=1. For ANY finite list of
coprime ordinary cables r_j,k_j>1 with

    r_j=1 modulo15, k_j=2 modulo15,
    rho(z_j)=rho(mu_j)^2*rho(lambda_j),               (10)

the full peripheral action is preserved. Writing R=product r_j,

    H_1(W_final)=Z/51,
    F_final=R^2 F+diag(6*(1-R^2)/15,0,0).

The cores generate a subgroup of order17 if 17 does not divide R,
and zero otherwise. The hidden quotient is Z/3 or Z/51. This is the
unbounded written cable theorem from PF111, not inference from the
sample cable checks. The parameters are unbounded: for any r=15h+1>1,
k=r*(2+15ell)+15 is coprime to r and meets (10) for every ell>=0.

If ANY of these were the actual FULL branched end cover of a plane
Keller normalization of degree17, the index15/f1 orbit would be an
omitted prime. But (1) requires e<=14. Thus the ENTIRE stated colored
iterated family is excluded, including all possible extra unramified
omitted images and arbitrary finite image singularities. The empty
cable list is allowed. This is not merely Orevkov's old e=d-1 case:
here e=d-2=15. The application follows from the new all-degree bound;
its knot computations are comparison evidence, not an additional
unproved source hypothesis. The underlying uncolored knots with other
representations, and all degree17 maps, are not excluded.

## 7. Evidence, genuine controls and the remaining goal

[The manifest](sources/high-index-budget/manifest.json) pins the
PUBLISHED Chau2003 PDF (all six complete pages freshly viewed), Orevkov
fullPDF7-10 freshly viewed, and unchanged predecessor contexts. The
Chau and Orevkov theorems are IMPORTED, not independently reproved.
The new written arguments are the arbitrary-index fibre collision
bound, its image-mass form (0)-(0a), two separate global costs, (1)-(2), the complete next-strip
profile, quotient consequences and the strengthened partial-link budget.
Formal verification and external review remain false.

The genuine PF107 finite A1-quotient cover has degree12, index9, a sole
omitted-axis image (t^2,t^3), and wrong open C* times C. It is etale
on that open, but it is NOT a polynomial Keller map on a WHOLE plane.
It violates the r>=2 conclusion and illustrates exactly why finite
covers with an arbitrary etale open cannot replace the source hypothesis.
PF110's all-p affine whole-plane exact-volume sources still supply no
finite compatible Keller target, and are not excluded as surfaces.

The strict inequality 2e>d is essential in Section2. The boundary
case d=6,e=3 is not eliminated by that argument; the separate retained
triple-torus exclusion has its own hypotheses. The e=d-3 profile for
d>=7 remains a necessary open class. Other indices, nonproper image
arrangements and the ACTUAL Laurent finite algebra remain OPEN.
Published Domrina2000 quartic retirement and shared trace priorities
are preserved; the full all-degree research objective stays active.

Run [the checker](verify_high_index_budget.py); its exact controls and
source seals are in [the receipt](high_index_budget_receipt.json).
Finite sample checks establish neither the unbounded index theorem
nor a Keller realization of any numerical profile.
