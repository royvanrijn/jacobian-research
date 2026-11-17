# Local monodromy capacity and the weighted full-branch exclusion

PF2RULED113, 16 November 2026. Work with a nonautomorphic complex plane
Keller map F, its ORIGINAL finite normalization pi:N->A2 of field degree
d>1, the WHOLE marked U=A2, and ALL r omitted primes. N may be singular.
Write B for the reduced FULL finite-map branch-image union, and S_F for
the whole nonproper-value set. These two curves need not be equal.

**Universal finite-fibre and local-monodromy theorem.** Every ambient
local finite fibre length ell_p satisfies

    ell_p <= d-r.

Consequently every COMPLETE normalization fibre has at least
ceil(d/(d-r))>=2 points. The local monodromy at every finite target
point has a nontrivial invariant partition by these local factors,
with each block of size at most d-r. In particular it is NEVER
transitive on all d sheets. This is about pi:N->A2, not the fibre of
the nonproper source restriction F. It does not assert that F is
surjective or forbid a singleton SOURCE fibre.

**Complete full-branch class exclusion.** The FULL B is never carried
by a polynomial target automorphism to the zero set of a nonconstant
polynomial homogeneous for STRICTLY POSITIVE integer weights. This
includes reducible weighted branch unions, arbitrary finite local
excesses, singular normalizations and every degree. The proof uses
global weighted radial contraction on the TARGET complement; it needs
no lift of a C* action to N or U and no source-action descent.

**Contractible branch corollary.** A nonempty irreducible FULL branch
curve B can NEVER be topologically contractible. The Lin-Zaidenberg/
Abhyankar-Moh-Suzuki classification, in Palka Theorem A, is an explicit
import for this corollary. Additional UNRAMIFIED omitted images outside
B are allowed. Thus this is stronger than applying Chau's obstruction
only when the WHOLE nonproper-value set is contractible.

**All-degree index improvement.** An omitted index e with 2e>d requires
another RAMIFIED prime with a DIFFERENT image, and

    e+r <= d-2.

In particular EVERY omitted index satisfies

    e<=d-4 for d>=7.                                (1)

The ENTIRE e=d-3 class left open in PF112 is now CLOSED for d>=7,
in every normal source and with every other prime retained. Its
unramified-companion necessary profile is not a remaining frontier.
PF112's e<=d-3 for d=5,6 remains valid; this does not prove inverses
in either degree or in the full remaining range.

**Shared-image improvement.** With w_C=sum_(D maps to C)e_D*f_D and
b_C the number of residue-one primes above C, the PF112 high-mass
condition now gives

    2w_C>d  =>  w_C+b_C<=d-3.                        (2)

For an irreducible FULL branch curve, its generic omitted mass satisfies
w_B<=floor(d/2), even when S_F has other unramified-image components.
This is not a bound on the total mass outside B. The PF112 bound for an
irreducible WHOLE S_F remains separate.

## 1. The Euler ledger bounds every local cluster

Use the all-normal prime-indexed boundary ledger of PF2OMIT2:

    d-1=sum_D e_D + sum_D E_D,
    E_D=sum_(p in D)(ell_p-e_D)>=0.                 (3)

For an omitted point p on D, isolate its summand and keep every other
prime and excess term. Then

    d-1 >= ell_p + sum_(D'!=D)e_D' >= ell_p+r-1,

so ell_p<=d-r. On U, pi is the original etale Keller map and ell_p=1.
The same bound holds there because (3) implies r<=d-1. Finite flatness
over the smooth target gives sum_(p in pi^-1(y))ell_p=d. Therefore

    #pi^-1(y) >= ceil(d/(d-r)) >= 2.                 (4)

The weaker ell_p<=d-1 already excludes a complete singleton fibre.
These are elementary consequences of the retained exact ledger, not
claims of historical novelty. The fibre can consist entirely of
several omitted points, so (4) gives no surjectivity theorem for F.
Dropping other primes and retaining the original d would invalidate
the bound used here.

## 2. Local factors give invariant monodromy blocks

Let pi^-1(y)={p_1,...,p_s}. Choose pairwise disjoint neighborhoods V_i
of these points. For a sufficiently small target ball Q about y,
proper finiteness gives pi^-1(Q) contained in their union: otherwise a
sequence of counterexamples over targets tending to y has a limit
over y outside that union, by properness. Shrinking Q if needed gives
the usual finite analytic local factors, of constant ranks ell_p_i.
This is exactly the local-degree conservation used in the boundary
ledger. At a regular basepoint of Q minus B, the factor in V_i contains
ell_p_i sheets. Path lifting cannot move between these disjoint open
and closed factors. Hence they form an invariant partition for

    rho(image(pi_1(Q minus B)->pi_1(A2 minus B))).

There are s>=2 nonempty blocks by (4), each of size at most d-r. Thus
the action is intransitive. No claim that a generic ramification index
equals an exceptional local degree has entered this proof.

More generally, any finite target neighborhood for which the induced
local complement group acts transitively on all d sheets is an
obstruction to the ACTUAL Keller normalization. Surjectivity onto the
global monodromy group is a sufficient way to force that transitivity.
A finite topological or algebraic cover without the WHOLE-plane Keller
hypothesis need not satisfy (3)-(4), and can have transitive local
monodromy. No exclusion based only on finite-order eigenvalues or
abstract local monodromy is being asserted.

The global cover pi^-1(A2 minus B)->A2 minus B is connected: its source
is a nonempty Zariski open subset of the irreducible normal N. Thus
its global sheet action is transitive. Purity gives B as the union of
images of the ramified omitted primes, with no extra isolated branch
component off that union. If B were empty, the original finite cover
would be etale over simply connected A2 and d=1. These standard
finite-cover facts are the retained PF46/PF98 context, not a new
assumption that unramified omitted images are absent.

## 3. Positive weighted radial flow makes the full local group global

Suppose first B={h=0}, with h a nonconstant polynomial and

    h(t^a*x,t^b*y)=t^L*h(x,y),  a,b,L positive integers.

Reducedness can be imposed on h, but only its zero set is used. The
origin lies in B. On its complement define the positive real flow

    T_t(x,y)=(t^a*x,t^b*y),  t>0.

It preserves the complement. For z=(x,y)!=0 the function

    ||T_t(z)||^2=|x|^2*t^(2a)+|y|^2*t^(2b)

is strictly increasing from zero to infinity: its derivative is
2a|x|^2*t^(2a-1)+2b|y|^2*t^(2b-1)>0. Therefore every orbit meets
each Euclidean sphere of positive radius exactly once. The meeting
time t_epsilon(z) is continuous (indeed smooth off the origin).
The homotopy T_(t_epsilon(z)^s)(z), 0<=s<=1, is a strong deformation
retraction of the global complement onto S_epsilon minus B.

Choose 0<epsilon<delta. The SAME homotopy stays in the ball of radius
delta whenever its initial point is there: the norm varies monotonically
between its initial value and epsilon. Thus both the local complement
Q_delta minus B and the global complement retract onto the SAME
S_epsilon minus B, fixing its basepoint. The inclusion induces an
isomorphism on fundamental groups. The local monodromy therefore
equals the globally transitive sheet action, contradicting Section2.

If B is only polynomially equivalent to such a curve, transport this
argument by that global target automorphism. Its small inverse-image
neighborhoods may be used instead of round balls in Section2.

This proves the complete positive-weighted FULL branch exclusion.
It does not follow from one weighted component, from a weighted
analytic germ alone, or from weights containing zero or negative
entries. None supplies the required global contracting flow on the
FULL complement. Individual local branch singularities can remain
weighted homogeneous in an allowed global arrangement.
An uncolored infinity knot alone does not certify global weighted
homogeneity or contractibility of the full branch curve.

## 4. An irreducible contractible full branch curve is impossible

Use the primary author paper Karol Palka,
*A new proof of the theorems of Lin-Zaidenberg and Abhyankar-Moh-Suzuki*,
arXiv1405.5391v2, 9 February2015, Theorem A. Its explicit convention
is that curves are irreducible. It states that an irreducible
topologically contractible complex algebraic plane curve has equation
x^n=y^m in polynomial plane coordinates, for coprime n>m>0. The
singular Lin-Zaidenberg case has m>=2; m=1 is the smooth embedding
case. This deep classification is IMPORTED, not independently reproved.
It was already retained in the HC4 contractible-curve work.

The resulting polynomial is homogeneous with positive weights (m,n).
Section3 excludes it as the FULL B. This proof continues to apply if
S_F has extra unramified omitted images elsewhere: they do not change
the finite-map branch complement. It makes no claim that F itself is
a covering outside B, or that the source open is C*-stable. The cover
to which the topological proof applies is the ORIGINAL finite pi.
Different contractible components cannot be rectified separately and
then declared simultaneously weighted homogeneous; a full union
requires its own global equation in one coordinate system.

## 5. Close the e=d-3 strip and sharpen shared-image budgets

PF112 proves that 2e>d gives f=1, injective image normalization,
contractible image C_R, and own excess E_R>=1. If every ramified prime
had image C_R, this would be the entire FULL branch curve. Section4
forbids it. Hence some RAMIFIED prime has a DIFFERENT image, costing
at least two in the global ledger. Counting all r-1 other primes at
cost one, then paying one extra for this other ramified prime and one
for R's own excess, gives

    d-1 >= e+(r-1)+1+1 = e+r+1,
    e+r <= d-2,   r>=2.                             (5)

Therefore e<=d-4 in the strict-high-index range. In the complementary
range e<=floor(d/2). For d>=7, floor(d/2)<=d-4. This proves (1) for
EVERY prime, not just a special knot or finite cover degree. It closes
the ENTIRE e=d-3 class for all d>=7. The sealed PF112 necessary
profile is now an excluded predecessor profile, not a current open
queue; its proof bytes are preserved.

For the image-mass version use the exact PF112 definitions w_C and
b_C. If 2w_C>d, that theorem makes C contractible and forces an own
excess unit on each residue-one prime above it. Such an image cannot
carry an unramified omitted prime, since PF2OMIT2 gives zero excess
and residue one on every unramified prime. Thus C is a branch
component. If it were the whole B, Section4 would contradict it.
Another RAMIFIED image component costs at least two outside C, so

    d-1 >= w_C+b_C+2,

which is (2). As before, known same-image sublists already give
necessary lower-bound filters without deleting any unknown prime.
If B itself is irreducible, 2w_B>d would make it contractible and
impossible; hence w_B<=floor(d/2). Omitted images outside B are
unramified and are not included in this last mass.

The local-cluster cap also strengthens PF111's partial-link budget:
for a link containing K_R, if g is the minimum generator count of its
unseen end-homology quotient, then r-1>=g. For strict-high index,
(5) implies d>=e+max(1,g)+3. The extra ramified prime and own excess
are separate costs. There is no automatic generation certificate
from zero rational-linking row defects.

## 6. The new remaining profile and quotient corollaries

For d>=9 and e=d-4, strict high index holds. The remaining three
budget units must pay for R's positive excess and another ramified
prime with a distinct image. Equality forces EXACTLY two omitted
primes, BOTH RAMIFIED:

    R: (e,f)=(d-4,1), one eta=1,
    S: (e,f)=(2,1), all eta=0.                       (6)

There are no unramified omitted primes or other excesses. C_R has
injective normalization and exactly one nonimmersive parameter, with
local length e+1=d-3 there and e elsewhere. The other image C_S is
different and immersively normalized by Orevkov5.2. It must identify
distinct parameters, or it would be a smooth A1 component forbidden
by Chau's component Remark. This is a COMPLETE NECESSARY profile,
not a realized finite Keller map or an exclusion of this next class.
Source smoothness at the R parameter is not inferred.

Combining (1) with PF106 gives cyclic quotient order n<=d-5 for d>=7.
In degree7 the only possible nontrivial quotient is A1. In degree8
one has n<=3; marked type1/n(1,-e), gcd(n,e)=1 and e<=4 leave A1
or A2. These are local source singularity restrictions, not an
identification of the end with an ADE link or an inverse theorem.
PF112's degree-five smoothness and degree-six A1 restriction remain.

## 7. Genuine finite controls and verification boundary

For EVERY m>=2, the genuine finite map (z,t)->(z^m,t) has a single
point of local length m over the origin, and transitive local cyclic
monodromy over its weighted branch line. Its restriction away from
z=0 is etale, but that open is C* times C, NOT A2. Its Jacobian on
the whole source is m*z^(m-1), and its putative Keller Euler ledger
would assert m-1>=m. The finite cover itself is not excluded.

The genuine cubic cover (s,t)->(s,t^3+s*t) has smooth source,
weighted cusp branch 4u^3+27v^2=0, generic index2 and a unique fibre
point of length3 at the cusp. Its critical parabola is s=-3t^2;
the etale open is again C* times C. Its Jacobian t-derivative is
3t^2+s, not a nonzero constant. It demonstrates why even a smooth
finite cover of a contractible full branch curve is allowed without
the marked Keller-plane hypothesis.

The retained PF107 degree12/index9 finite A1 quotient has full branch
the cusp (t^2,t^3), unique fibre point over its vertex, and wrong
open C* times C. Its homogeneous invariant fibre algebra has length12.
It is an ACTUAL singular finite cover, not an abstract lattice, and
is not disproved by this theorem. PF110's whole-plane affine exact-
volume sources are likewise retained; no finite compatible target
is supplied or excluded as a surface by these statements alone.

[The checker](verify_local_monodromy_branch.py) seals predecessor
proofs and the exact Palka source version, checks the unbounded
cluster/index/mass arithmetic and positive-weight support identities,
and verifies the genuine finite controls including the singular
length12 fibre. The universal weighted-flow and local-factor proofs
are written above; they do not rest on finite samples or a group census.
See [the receipt](local_monodromy_branch_receipt.json) and
[source manifest](sources/local-monodromy-branch/manifest.json).

The Palka classification and PF112's Chau/Orevkov inputs remain
declared imports. No historical priority, formal verification or
external review is claimed. Branch arrangements with no local group
giving a global transitive action, the new two-ramified-prime profile
and the ACTUAL Laurent finite algebra remain OPEN. Published
Domrina2000 quartic retirement and shared trace priorities persist;
the full research objective remains active.
