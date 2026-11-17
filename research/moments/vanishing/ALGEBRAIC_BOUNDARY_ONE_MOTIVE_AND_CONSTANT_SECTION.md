# Algebraic boundary 1-motive and the constant-section criterion

16 November 2026. Internal written proofs with named deep primary imports.
Independent replay, formal verification and external review are not claimed.
Both COMPLETE exceptional Mathieu--Zhao cases remain **OPEN**. Review deferred.

The remaining intrinsic source condition has an algebraic geometric model
for EVERY phase: a rank-one lattice mapping to a semiabelian family. Unit
splitting is equivalent to this ONE section coming from a constant point
of a constant semiabelian family, after an isogeny. It does not require
the section to be torsion or the mixed Hodge extension to split.

Under EVERY polynomial recomposition the semiabelian family undergoes
Weil restriction and its section retains ALL branches. A scalar norm
is insufficient, even on an actual eventually pure failing phase.

See the [checker](verify_haar_boundary_motive.py),
[receipt](haar_boundary_motive_verification.json) and
[source audit](haar_boundary_motive_source_audit.json).

## 1. The full geometric source has a rank-one 1-motivic lattice

Let f be ANY nonconstant regular phase on X=SL2/Gamma, ANY finite Gamma.
Use the FULL regular algebraic base S and canonical boundary source E
from the [intrinsic safety theorem](INTRINSIC_SOURCE_SPLITTING_AND_FULL_HAAR_SAFETY.md).
Its original types are (2,2),(3,2),(2,3),(3,3), and W4 E=L=Q(-2)
has rank1. Therefore

    H=(E(2))^dual                                      (1)

has types (0,0),(-1,0),(0,-1),(-1,-1), with GrW0 H=Q of rank1.
All the genus and Tate graph remain in W-1 H. If their ranks are 2g
and t respectively, rank H=1+2g+t. No genus or graph vanishing is imposed.

We construct, up to isogeny and after possibly deleting extra REGULAR
parameter values, an ALGEBRAIC 1-motive

    M_f=[Z -> G_f],  1 |-> sigma_f,
    0 -> T_f -> G_f -> A_f -> 0,                       (2)

whose Hodge realization is H. Here A_f is an abelian scheme of dimension
g, T_f is an algebraic torus of rank t with finite character monodromy,
and sigma_f is an algebraic section. The lattice is globally constant.
This is a geometric construction from compactification/Picard data,
not an unrestricted equivalence for every abstract relative VMHS.

Choose a smooth proper relative compactification Y of U_s=X minus F_s,
with SNC boundary D, over a smaller algebraic base. Resolution and
spreading supply this for every f; no tameness is used. Poincare duality
and the localization sequence give

    B^dual(-2)=Hc3(U_s)(1)=H3(Y_s,D_s)(1),
    H2(Y_s) -> H2(D_s) -> Hc3(U_s) -> H3(Y_s).          (3)

The last group is pure of weight3; the entire compact-support group
in (3) has weights0,1,2, because B has original weights4,5,6. Hence the
map to H3(Y_s) is zero and the boundary map onto Hc3 is surjective.
After twist1, (1) is a sub-MHS of the 1-motivic part of this relative
cohomology, with the SAME entire W-1 part and one constant Tate lattice.

The precise deep import is Barbieri-Viale--Rosenschon--Saito (BRS),
Annals158(2003), Theorem0.1, Construction3.1 and comparison5.4--5.6.
It constructs the 1-motivic part of RELATIVE cohomology algebraically,
with its extension class, not just its graded groups. Only the rational
form is needed. No integral strengthening from Theorem0.2 is assumed.

Here are the finite algebraic inputs to that construction. Resolve the
pair by a complex of smooth proper varieties and form its derived Picard
groups P>=i. A point of P>=i is a line bundle on its i-th term together
with a trivialization of its alternating restriction to the next term;
that trivialization must satisfy the cocycle identity on the following
term. Construction3.1 forms the discrete lattice from kernels/cokernels
of Neron--Severi maps, and the semiabelian part from

    G_r=coker(P_(r-2)^0 -> P_>=(r-1)^0),  r=3.         (4)

The boundary morphism from that lattice to G_r is sigma before the
rank-one pullback. BRS3.2--3.3 prove representability and semiabelianity;
BRS5.4 identifies the ACTUAL extension class. Triple-intersection gluing
is essential: Jacobians alone omit its torus part.

Apply this construction to the generic geometric pair and pull back its
lattice along the ambient rank-one subspace GrW0 H. This subspace is
algebraic: its type(0,0) classes lift, after finite cover and integer
rescaling, through the Neron--Severi complex by the pure semisimplicity
and the degree-two Lefschetz(1,1) theorem on the smooth proper terms.
The map of that rank-one lattice into (4) gives (2).

For precision about families: perform the construction over a finitely
generated field of definition and its generic parameter, spread the
FINITE chosen Picard/NS data and comparison maps, and take a finite
extension of the parameter field when needed. After shrinking, its
normalization is a finite etale cover. BRS's natural Picard, cycle and
relative-cohomology maps commute with the geometric comparison and
Gauss--Manin transport. Their comparison on very general fibres gives
the required isomorphism of variations on the open base. This spreads
the CHOSEN ambient lattice, not every fibre's entire Picard group:
Noether--Lefschetz jumps are not declared absent or removed by deleting
only finitely many points. Galois descent up to isogeny is effective:
use Weil restriction and the rational averaging idempotent, whose image
exists in the isogeny category. Thus (2) can be taken on an open subset
of the ORIGINAL parameter line. No general relative essential-
surjectivity theorem or zero-fixed-part assumption is required.

Deleting additional regular values changes no monodromy image. Moreover
unit splitting is preserved AND reflected by ANY finite base cover:
on its normal closure average a fixed ambient lift over the finite deck
group. A rational average still maps to1. Consequently these construction
choices do not weaken the original GLOBAL source-unit test.

This supplies a finite geometric construction from a supplied resolved
pair, not a software implementation of resolution, Picard schemes,
fixed-part computation or the full exceptional pure locus.

## 2. The exact condition is a constant point, not Hodge splitting

For a section sigma of a semiabelian scheme G/S, say CONSTANT TYPE UP TO
ISOGENY if there exist a constant semiabelian variety G0 over C, a point
a in G0(C), an algebraic homomorphism h:G0 times S -> G (on a nonempty
open base), and a nonzero integer m such that

    [m] sigma = h(a).                                  (5)

All statements use rational realizations. This property is invariant
under isogeny, finite base cover and deletion of extra regular values.
The finite-cover assertion follows either from the proof below and
averaging, or by descending the sum of the conjugate constant images.

The general rank-one 1-motive lemma is

    underlying local system of [Z -> G] has a GLOBAL
    invariant lift of its lattice generator
       iff sigma has constant type up to isogeny.       (6)

To prove sufficiency, [Z -> G0] defined by a is a CONSTANT mixed Hodge
structure/local system. The morphism to [Z -> G] with lattice map m and
semiabelian map h sends a fixed lift to one mapping to m. Divide by m.

For necessity take the FULL invariant part F of the Hodge realization
of [Z -> G]. The mixed fixed-part theorem makes F a constant MHS. A
fixed lift means F -> GrW0=Q is onto; strictness gives GrW0 F=Q. F is
itself of 1-motivic type, so the POINT equivalence of Deligne, in the
precise rational form BRS Proposition1.5, gives a constant 1-motive
[Z -> G0]. The inclusion F -> H is a horizontal Hodge morphism.

Its semiabelian part is algebraic, after shrinking S if necessary.
Indeed fibrewise full faithfulness of Proposition1.5 gives the unique
semiabelian homomorphism. A horizontal Hodge map varies holomorphically;
the relative Hom functor of semiabelian schemes is unramified and locally
of finite type, by rigidity of abelian homomorphisms and the discrete
character lattice for tori. Its chosen horizontal branch lies generically
in a finite type, quasi-finite Hom component. On a dense open it is a
finite etale branch; the analytic section selects a degree-one component,
so it is algebraic there. This suffices because monodromy images do not
change on shrinking. Clearing rational denominators makes the lattice
component a nonzero integer m. Commutativity of the 1-motive morphism
then gives exactly (5). This uses full faithfulness on points and Hom
rigidity, not relative essential surjectivity for arbitrary variations.

Duality in (1) identifies the fixed lift in (6) with the intrinsic
source retraction p:E -> Q, p(u)=1. Combining (6) with the intrinsic
safety theorem yields the ALL-PHASE criterion

    under EVENTUAL PURITY,
    ALL original markers are safe iff sigma_f has
    constant type up to isogeny.                       (7)

Constant type is NOT the assertion sigma is torsion or that [Z -> G]
splits as a mixed Hodge structure. The constant example G=Gm, sigma=2
has trivial monodromy and passes (6). Its rational Hodge extension is
nonsplit: log(2) is not a rational multiple of 2*pi*i. Conversely simply
quotienting by a constant semiabelian subgroup is too weak: G=Gm over
C*, sigma(s)=s, has zero quotient by the entire constant group, but its
logarithm has nonzero winding and (6) fails. The POINT in (5) must be
constant; an arbitrary varying section of a constant group is insufficient.

## 3. EVERY recomposition is Weil restriction with a diagonal lattice

Let Hpoly be ANY nonconstant polynomial. On full compatible regular
bases let pi:T -> S' be the finite etale map r=Hpoly(tau). Then, up to
isogeny, the boundary 1-motive for Hpoly(f) is

    [Z -> Res_(T/S') G_f],
    sigma_Hpoly(f) = the ALL-BRANCH tuple of sigma_f.    (8)

Its abelian dimension is deg(Hpoly)*g and torus rank deg(Hpoly)*t.
Its lattice STILL has rank1. To see (8), the established boundary
source relation is

    E_Hpoly(f)=pi_*E_f / (trace-zero ambient L).

Dualizing makes the lattice the DIAGONAL copy of Q in pi_*Q, while
leaving ALL negative-weight branch systems. The Betti realization of
Weil restriction is exactly that induced system. The algebraic 1-motive
in (8) therefore has the same geometric comparison; fibrewise full
faithfulness and Hom rigidity identify them on an open base. The dual
of ambient SUM is diagonal1, with no division by the number of branches.

This is NOT a scalar norm or a product discarding the norm-one torus.
Constant type of the full tuple is equivalent to primitive source-unit
splitting by finite-cover restriction/averaging, before purity. Thus (7)
decides EVERY existing eventual-pure Hpoly(f) through ONE primitive
constant-type condition, without assuming f itself pure or constructing
such Hpoly. The earlier actual-rank exclusion and safe-centre theorem
remain: source-split actual rank>1 permits NO eventual-pure composite;
rank1 has a unique safe centre c and exactly Hpoly(c)=0.

## 4. A genuine pure failure has constant scalar norm

Retain Long's known nonexceptional phase and marker

    F=(1+y)(1+z+yz), Q=-y, xw-yz=1,
    h(F^k)=0, h(Q F^k)=(-1)^(k-1)/(k+1), k>=1.

Its boundary source has rank2, with genus0 and torus rank1. The Kummer
representative of its section, modulo a constant multiplier and rational
isogeny normalization, is g(s)=s/(s+1), since

    C_Q=1/s+log g(s), dlog g=1/[s(s+1)].               (9)

The geometric 1-motive construction above justifies an algebraic section;
the residue divisor fixes this representative modulo constants. No claim
about the constant mixed-Hodge extension value is needed.

Put tau=F+1/2 and Hpoly(tau)=tau^2-1/4. Then the ACTUAL composite is

    P=F(F+1).

It is pure in EVERY positive order by binomial expansion. The SAME fixed
original marker has the exact NONZERO tail

    h(Q P^n)=(-1)^(n-1) integral_0^1 x^n(1-x)^n dx
            =(-1)^(n-1)*(n!)^2/(2n+1)!, n>=1.         (10)

But the two branch representatives in (8) obey

    g(tau)=(tau-1/2)/(tau+1/2),
    g(-tau)=1/g(tau),  Norm(g)=1.                     (11)

An unnormalized constant multiplier only changes this to a CONSTANT
norm. The full tuple lies in the norm-one torus but moves; its section
is not of constant type. The actual Cauchy function retains the difference
of logarithms with the trace weights, not their unweighted sum:

    C_Q^P(r)=1/r+(1/tau)log((tau-1/2)/(tau+1/2)),
    tau^2=r+1/4.                                      (12)

The physical determinations agree near infinity. At r=0 its logarithm
has nonzero coefficient. This is an ACTUAL pure example, not a free
choice of an obstructing relation kernel. It remains on SL2 and is NOT
right2O/2I invariant; it supplies no exceptional counterexample.

## 5. The whole-quotient target remains full exceptional constancy

For EVERY actual exceptional eventually pure phase, prove (5) for its
geometrically constructed section sigma_f, or realize a genuine admissible
nonconstant-type section with an original pure phase and bad marker.
This is the SAME necessary-and-sufficient full target, now algebraic.
It covers all degrees, markers, genus, graph ranks and polynomial branches.
No constancy theorem on the entire exceptional pure locus is proved here.
No general numerical NS/Picard/trace computation is implemented, no
coefficient-ring topological promotion or universal machine proof is claimed.

The recent relative-equivalence preprint arXiv:2603.16545v3 is NOT an
import. Its no-fixed-part-after-every-finite-cover assumption fails for
our allowed nonzero toric part: on a splitting cover, H1 of that torus
is a nonzero invariant subspace of H1 G. Its displayed toric forgetting
map cannot serve as an unrestricted injectivity assertion because of
the constant sigma=2 example. The published BRS geometric construction,
point equivalence and actual mixed fixed part are the precise inputs used.

Both COMPLETE exceptional cases remain OPEN; review remains deferred.
