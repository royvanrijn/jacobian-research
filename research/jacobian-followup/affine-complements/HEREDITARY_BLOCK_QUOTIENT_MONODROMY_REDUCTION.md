# The primitive reduction holds for every block quotient

PF2RULED122, 16 November 2026. Let G act on the d sheets of the
ORIGINAL finite normalization of an ACTUAL nonautomorphic complex
PLANE Keller map. Every nontrivial transitive BLOCK QUOTIENT of this
action has two nonidentity branch meridians with DISJOINT supports.
Every PRIMITIVE block quotient is necessarily natural alternating/
symmetric or binary affine. This constrains ALL hypothetical maps,
including imprimitive monodromy, without recognizing an intermediate
source as A2. It does not prove an inverse theorem or close the full
finite-algebra problem. The full goal remains ACTIVE.

## 1. Exact support budget on a genuine block quotient

Retain the whole marked U=A2 and ALL omitted primes of the original
normal finite pi:N->A2; N may be singular. Let Omega->Gamma be a
G-equivariant surjection onto b>1 blocks of size t, so d=tb, and let
Q be the faithful image on Gamma. This includes t=1. It does NOT
include an arbitrary abstract homomorphic image with an unrelated
permutation action. The finite intermediate normalization M->A2 exists
by the block/field correspondence; no plane, smoothness, or whole-plane
open in M is assumed.

Index all irreducible ORIGINAL full branch images by C_j. Write
sigma_j for a meridian, w_j for its original support, a_j for its
support on Gamma, and delta_j=w_j-t*a_j. A moved block has ALL its
letters moved, hence delta_j>=0. The exact original support ledger is

    tb-1=sum_j w_j+t0+Eta,
    t*(b-sum_j a_j)=1+t0+Eta+sum_j delta_j.           (1)

In particular sum_j a_j<=b-1. Terms with a_j=0 are retained in (1);
the full branch B_Q of M is exactly the union of C_j with a_j>0.
It is nonempty: otherwise purity makes M a connected finite etale
cover of A2, forcing b=1.
The original PF120 cap also gives, for EVERY quotient branch meridian,

    a_j<5b/6.                                       (2)

An involution on blocks need NOT lift to an involution on Omega. We
do not import PF120's stronger involution cap for such a lift.

## 2. The full weighted-branch exclusion descends without a plane source

Suppose B_Q is polynomially rectifiable to a reduced curve homogeneous
for strictly positive weights. At its weighted origin y, radial
contraction identifies the global and small local complements of B_Q.
Thus the LOCAL quotient monodromy at y is transitive on all b blocks,
as in PF113. Extra original branch curves not in B_Q may be present;
loops can be perturbed off them, so the original local monodromy group
K still projects onto this transitive quotient group.

Each original finite local factor O at p over y is a K-orbit of sheets.
Its nonempty block image is K-invariant, hence equals all Gamma. In
particular |O|>=b>1: no factor is unramified or in U. Each p lies on
one ramified omitted curve and its local branch image is ONE analytic
branch, by the disjoint/unibranch original boundary theorem.

Moreover that image belongs to B_Q. If it did not, meridians around
every B_Q branch would act trivially on O, since this factor extends
unramified across those branches. Their conjugates normally generate
the quotient local image; the other original branch meridians project
trivially. The block action induced from O would be trivial, contrary
to its transitivity on b>1 blocks. This uses normal meridian generation
and invariance of O, not a chosen global generating list.

Every irreducible positive-weighted plane component has injective
normalization A1 and passes through y. Explicitly its reduced equation
is a coordinate axis or x^u-c*y^v with coprime u,v, in the chosen
polynomial coordinates. Thus it has just ONE normalization parameter
over y. For each original D_i mapping to such a component, all points
p on D_i over y have normalized residue degrees summing to f_i.
Using the exact normalized local excess eta_p=ell_p-e_i*q_p,

    sum_(p on D_i over y) ell_p=e_i*f_i+sum_(p over y) eta_p.

All points of the complete original fibre are on these curves, and
the curves are disjoint. Finite flatness of the original N gives

    d=sum_(all p over y) ell_p
      <=sum_(ALL D_i) e_i*f_i+sum_(ALL p on D_i) eta_p=d-1,

a contradiction. Other primes were not deleted; they only increase
the upper bound. This works for f_i>1, quotient points and arbitrary
indices. It assumes neither an original singleton fibre nor a smooth
intermediate model.

Therefore EVERY block quotient has no positively weighted FULL branch
union. In particular its FULL branch curve cannot be irreducible and
contractible, by the precise Palka/Lin--Zaidenberg import retained in
PF113. Extra unramified images and original branch images invisible
in Q do not change this conclusion.

## 3. Disjoint branch supports descend to EVERY block quotient

If B_Q had no finite multibranch point, its components would be
disjoint with injective A1 normalizations, hence contractible. One
component contradicts Section2. For several components, the argument
in PF121 Section2 applies to these ORIGINAL image curves: Chau rules
out a smooth A1 component of the original S_F; a singular contractible
component can be rectified to x^u-y^v=0, u,v>=2 coprime. A second
disjoint A1-parametrized component would map nonconstantly into the
positive-genus shifted fibre x^u-y^v=c!=0, impossible by
Riemann--Hurwitz. Only ONE component is rectified.

Take two distinct analytic B_Q branches at such a multibranch point.
PF121 gives original meridians sigma,tau with disjoint nonempty sheet
supports. Their quotient meridians are both nonidentity. If sigma
moves a block, it moves EVERY letter in it, so tau fixes that whole
block pointwise. Thus their quotient supports are also disjoint.
Consequently

    mu(Q)<=floor(b/2).                              (3)

This gives the half-support obstruction to EVERY block quotient,
without transferring the original boundary hypotheses to M.
It immediately excludes quotient degrees2 and3, and every faithful
interior k-subset action by PF121's elementary oriented-star proof.
The PF118 natural nonbinary affine and stated projective actions
are likewise excluded as block quotients by their strict-majority
minimum support. No finite group census is needed.

## 4. Product actions on a primitive quotient are excluded

Apply the published Burness--Guralnick primitive classification to
(3). Retain PF120's exact table filter: apart from alternating/symmetric
subset, binary affine and the three binary classical actions, only
product actions still require consideration. In that product case
Gamma=Delta^s, s>=2, with an almost simple primitive base of size c>=5.

PF120's elementary fixed-point dichotomy says that a nontrivial top
coordinate permutation is either conjugate to a pure coordinate swap
(fixed ratio1/c), or has fixed ratio<=1/8. For c>=6 both violate (2).
The only remaining base size is c=5, the natural A5/S5 base. Here
only a pure coordinate swap is possible as a branch meridian with
nontrivial top; its moved support is 4b/5.

Every nonidentity element of Sym(5) wr Sym(s) in product action moves
at least 2b/5: a nonidentity base coordinate fixes at most3/5, and
a nontrivial top has the stronger dichotomy above. A second quotient
branch component would therefore make sum a_j>=6b/5>b-1, contrary
to (1). Hence this swap would require B_Q irreducible. Its ORIGINAL
branch support is w_j>=t*(4b/5)=4d/5>d/2; PF112 then makes C_j
contractible. Section2 forbids this irreducible FULL B_Q.

Thus all quotient branch inertia has trivial top. Its normal closure
is Q (purity and trivial finite etale covers of the target A2), so
the top group is trivial. The resulting invariant coordinate-fibre
partition contradicts primitivity. ALL primitive product-action block
quotients are excluded, without assuming a lift of a swap is an
involution or applying any plane theorem to M.

## 5. A transvection's support has trivial pointwise stabilizer

For each of the three binary classical actions retained in PF120,
we prove that ANY group element fixing the support of a transvection
pointwise is the identity. The second element need not be a transvection.

For Sp_(2a)(2) on minus quadratic forms, a>=3, let
S_v={q:Arf(q)=1,q(v)=0}, the support of t_v. Differences of pairs
of forms in S_v include B(u,-) for every u in v-perp. Indeed one
can choose a minus form q with q(v)=q(u)=0: for independent orthogonal
v,u extend to two symplectic pairs and put the minus type on the
remaining dimension>=2; u=v is immediate. Then q+B(u,-) is also
minus, and still vanishes on v. An element fixing S_v pointwise
therefore fixes v-perp pointwise. Its symplectic possibilities are
1 and t_v. The latter moves S_v, so only1 remains.

For O_(2a)^+(2) on nonsingular vectors, a>=4, and O_(2a)^-(2) on
nonzero singular vectors, a>=3, the support of t_v is
S_v={x:B(x,v)=1,q(x)=epsilon}, with q(v)=1 and epsilon1 or0.
This set spans V. To check it, fix u in v-perp. The case u=v uses
any x in S_v, since x+v is also in S_v. If u,v are independent,
impose B(x,v)=1 and B(x,u)=q(u). This affine codimension2 subspace
has direction of alternating rank at least 2a-4>=2, so its quadratic
restriction is nonconstant and takes both F2-values. Choose q(x)=epsilon;
then x,x+u belong to S_v. Their differences span v-perp, and S_v
also contains a point outside that hyperplane, so spans V. A linear
element fixing S_v pointwise is therefore the identity. Triality and
the degree27 small identification preserve this permutation property.

Now apply the explicit Guralnick--Magaard ALL-ELEMENT inclusive-half
import recorded in PF121: of the disjoint quotient branch pair, one
moves at most b/2, and in a retained binary classical action it is a
transvection. The other fixes its support pointwise, hence is the
identity by the preceding proof, a contradiction. This excludes all
three classical primitive block quotients, without asserting that
EVERY quotient branch inertia fixes half the blocks.

The GM primary original-paper web text was read; its original PDF
was NOT obtained or visually audited, as PF121 honestly records.
No OCR subset-parameter bound or independent CFSG proof is used.

## 6. Current all-map consequence and precise remaining boundary

EVERY primitive block quotient of ANY actual plane Keller monodromy
is necessarily

    natural A_b/S_b; OR binary affine F2^a:H, b=2^a.

For the affine quotient, mu=b/2 and (1) force exactly ONE quotient
branch image C, whose quotient inertia has half-support and is a
transvection (PF118's elementary affine fixed-hyperplane calculation).
Section2 makes C noncontractible. Its ORIGINAL support w_C satisfies
w_C>=d/2, while w_C>d/2 would make C contractible by PF112. Hence

    w_C=d/2=t*a_C.

There are no extra original moved sheets inside blocks fixed by this
meridian. Its original ramification indices above C are even, but
they need NOT all be2: a single original cycle may wind repeatedly
around a two-block orbit. Other original branch images may project
trivially. The linear part H is irreducible and normally generated by
linear transvections, by projecting affine inertia normal generation.

The original PF119 two-cusp-profile theorem is retained when t=1.
It is NOT transferred to a general intermediate model or claimed
for all original indices above such an affine quotient. Likewise a
degree4 QUOTIENT is not ruled out by published degree4 PLANE inversion:
M has not been recognized as A2. Degree2/3 quotients are excluded by
Section3, and natural A5 is excluded by mu>5/2. Prime quotient degree
therefore forces natural alternating/symmetric, but not inversion.

This is a restriction on G-equivariant block quotients, NOT arbitrary
quotient representations, a theorem about all subgroup composition
factors, or an induction treating intermediate sources as planes.
Natural alternating/symmetric quotients, binary affine quotients and
the actual remaining imprimitive/global finite-algebra compatibility
remain OPEN. Local cusp S3 orbit profile 3^(d/4)1^(d/4) is compatible
with binary affine action; McLaughlin's finer linear classification
and a general cyclic-complement claim were not imported as closures.
All original primes, sealed predecessors, shared analytic/multiplication
priorities, PF117 local trace closure and published quartic retirement
are preserved. No historical priority, formal verification, external
review or full goal completion is claimed.
