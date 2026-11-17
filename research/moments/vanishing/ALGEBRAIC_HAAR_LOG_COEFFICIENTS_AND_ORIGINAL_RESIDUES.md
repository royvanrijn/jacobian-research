# Algebraic Haar log coefficients and original residues

16 November 2026. Internal written proofs; independent replay, formal
verification and external review are not claimed. Both COMPLETE exceptional
Mathieu--Zhao cases remain **OPEN**. Review remains deferred.

For a fixed phase, ONE finite algebraic function field contains every
integer-log coefficient of every original marked Cauchy period at any
puncture with semisimple eigenvalue-1 middle-fibre monodromy. The field
does not grow with multiplier degree. Those coefficients are explicit
pairings of the ACTUAL projected Haar cycle with algebraic double residues
at boundary nodes. Finitely many ORIGINAL markers separate every possible
Tate graph class.

The original unit can replace all those tests only under an additional
cyclicity condition. A Wronskian gives that condition exactly. It is not
universal: a genuine SL2 square composite hides an integer logarithm from
the unit while a fixed marker sees it, and the known actual pure safe
exceptional family has an unbounded unit-invisible graph space.

See the [checker](verify_haar_residue_log.py),
[receipt](haar_residue_log_verification.json) and
[source audit](haar_residue_log_source_audit.json).

## 1. One finite field for the entire Tate system

Let X=SL2/Gamma for ANY finite Gamma, let Omega be its algebraic volume,
and put Z=integral_M Omega!=0 on the compact core. Fix a nonconstant
regular f. Let S be the parameter line with the FULL finite bifurcation
set removed, including finite bifurcations caused by infinity. The
nearby fibres F_s={f=s} and complements U_s form algebraic families on S.

By the retained [weight theorem](HODGE_WEIGHTS_OF_HAAR_VARIATION_AND_BOUNDARY_GRAPHS.md),

    T_s=W_(-4)H_2(F_s;Q)
       =(Gr_W^4 H2(F_s;Q))^*                           (1)

is dual to a pure Tate variation of type (2,2). Its monodromy is FINITE
globally on S. Indeed the integral cohomology lattice induces an invariant
torsion-free lattice on this rational weight quotient. Its flat pure
polarization is positive definite, since the only Hodge type is (2,2).
The integral isometries of a positive definite lattice form a finite group:
the image of every fixed basis vector has fixed norm and finitely many
possible integral coordinates. This is also Proposition3.5(1) of
Doran--Harder--Vanhove; no Feynman-integral result is used here.

There is consequently a finite algebraic cover of S trivializing (1).
We may enlarge its function field to one finite extension

    L / C(s)                                          (2)

with a smooth projective compactification Y_L of the GENERIC fibre and
SNC boundary D_L, all components and nodes geometrically labelled.
To see existence directly, compactify and resolve the generic surface
over an algebraic closure of C(s). This finite-type construction and its
finitely many component/node data are defined over a finite extension.
Spread out and shrink S to retain a relative SNC compactification, with
its nodes split after a further finite cover. Enlarging the finite
extension or shrinking this open set does not alter the generic field
statement or a local continuation to an original puncture.

In the labelled model, let E,V be the graph edges and vertices, choose
orientations, and let partial:Q^E->Q^V be its incidence map. Then

    Gr_W^4 H2(F)=ker(partial)(-2),
    T=Q^E/image(partial^t), with Tate twist (2).        (3)

These are flat identifications after (2). The graph may be disconnected
if the generic fibre is. No connected-fibre or tameness assumption is made.

## 2. Residues of EVERY original top form

For Q in C[X], define the regular relative two-form

    eta_Q(s)=Res_(f=s) (Q Omega/(s-f)).                 (4)

For example, if Omega=df wedge beta locally, eta_Q=-Q beta. It is
meromorphic on Y_L with poles confined to D_L, of finite but possibly
HIGH order depending on Q. At an oriented node e with local boundary
coordinates u,v, set

    r_(Q,e)(s)=[u^-1 v^-1] eta_Q(s)/(du wedge dv).     (5)

This lies in L. It is a finite Laurent-coefficient calculation in algebraic
local parameters, not an assumption of simple poles. Repeated coefficients
are obtained from finite algebraic jets after clearing the actual pole
orders. Coordinate changes preserving the labelled boundary branches
preserve the oriented small torus period, so give the same residue.
For the small normal two-torus tau_e,

    integral_(tau_e) eta_Q=(2*pi*i)^2 r_(Q,e).          (6)

The logarithmic residue comparison identifies these tori with the edge
generators of (3). Their relations are the vertex coboundaries. Directly,
the normal circle bundle over a boundary curve with little discs removed
at its nodes is a compact three-chain whose boundary is the signed sum
of the incident node tori. Stokes therefore gives

    r_Q=(r_(Q,e))_e in ker(partial) tensor L.           (7)

Exact meromorphic two-forms have zero double Laurent residues: each term
is a partial derivative and cannot supply a nonzero u^-1 v^-1 coefficient.
Thus (5) depends only on the top-weight cohomology class. Deligne's
logarithmic comparison and the top-weight graph identification show that
it is precisely the projection to (3), with normalization (6). Genus
classes and lower weights are killed; no splitting of their representatives
is needed.

For c in Q^E/image(partial^t), represented by edge coefficients c_e,

    integral_c eta_Q=(2*pi*i)^2 sum_e c_e r_(Q,e).      (8)

The result is independent of the representative by (7). Projected Haar
classes are RATIONAL here; a spectral projector need not preserve the
integral lattice of the unprojected tube sequence.

## 3. The exact integer-log coefficient and the physical cycle

At a finite puncture a, assume middle-fibre monodromy T_F is semisimple
on generalized eigenvalue 1. Continue the ACTUAL compact Haar cycle
from infinity, and write

    (T_U-I)gamma_s=Tube(d_f),
    c_a=d_(f,1)=pi_1(T_F)d_f.                          (9)

The retained weight theorem puts c_a in T_s, fixed by local T_F. Its flat
continuation on the cover (2) is a constant rational graph class. The
generalized eigenvalue-1 part of complement monodromy has square-zero
logarithm: its quotient is semisimple and its kernel is the ambient
constant line. Consequently the INTEGER exponent class of the original
Cauchy period has at most ONE logarithm:

    C_(Q,1)(s)=A_Q(s)+B_(Q,a)(s)log(s-a),              (10)

where A_Q and B_(Q,a) are single-valued meromorphic germs at a.
Other exponent classes, including their possible logarithms, remain
separate. Neither A_Q nor the full C_Q is claimed algebraic.

With positive normal-circle orientation and (4), one loop changes (10)
by 2*pi*i B_(Q,a). The tube residue identity and (8) therefore give

    B_(Q,a)(s)=(2*pi*i)^2/Z * <c_a,r_Q(s)>.           (11)

This is the actual continued coefficient, not an arbitrary solution of
a Picard--Fuchs equation. It belongs to the SAME finite field L for
EVERY original Q and every a satisfying the stated semisimplicity
condition. The constant normalizations are complex constants, so are
in its ground field. A choice of sheet and its actual c_a determine the
germ; local T_F-invariance makes that particular algebraic germ
meromorphic at a even if other sheets of L ramify there.

Equation (11) needs no purity assumption. It gives a uniform algebraic
degree bound [L:C(s)] for all these coefficients for the FIXED phase,
not a numerical bound uniform over phases or matrix degrees. Finite
monodromy of the Tate system does not split the full complement extension;
the logarithm in (10) is exactly the remaining extension data.

## 4. A finite ORIGINAL-marker test, with no degree campaign

Let r=dim_Q T. There exist r ORIGINAL markers Q_1,...,Q_r such that
their residue vectors form an L-basis of ker(partial) tensor L.
Indeed eta_1 trivializes the canonical module of the smooth affine
generic fibre. Every regular relative top form is a finite L-linear
combination of the eta_Q with ORIGINAL Q in C[X]. These top forms
represent all H2_dR of that fibre. Projecting to its top weight shows
that all r_Q span (3), so an r-element basis can be selected.

These markers depend only on f and a chosen generic residue model,
not on the puncture a, the actual continuation, or the degree of a later
marker. In particular, at every puncture as in Section3,

    c_a=0 iff B_(Q_i,a)=0 for all 1<=i<=r.             (12)

If f is EVENTUALLY PURE and a=0, these equivalent conditions are also
equivalent to safety for EVERY original multiplier. Any nonzero c_0
gives an original Q_i with a genuine integer logarithm and an infinite
marked tail. When r=0 the test is empty and safety follows under the
same hypotheses. No original-marker degree bound, source/compactification
producer or implemented all-phase decision algorithm is supplied.

This result applies to every polynomial composite on its full fibre,
including all disconnected primitive branches. Alternatively combine
(11) branchwise with the exact
[primitive transfer theorem](PRIMITIVE_HAAR_VARIATION_AND_ALL_POLYNOMIAL_COMPOSITES.md).
Do not discard root branches or infer purity of the primitive factor.

## 5. Precisely when the original unit suffices

In a constant flat basis of (3), write the ORIGINAL UNIT residue vector
as rho(s)=(rho_1,...,rho_r) in L^r. The constant field of the connected
function field L with derivation partial_s is C. Thus

    K_unit={c in C^r : sum_j c_j rho_j(s)=0 in L}       (13)

is an exact finite-dimensional space of graph directions invisible to
the unit. It is not silently identified with the actual Haar variation.

The Wronskian

    W_unit=det(partial_s^(i-1)rho_j)_(1<=i,j<=r)       (14)

is nonzero iff K_unit=0. Here the flat basis must be CONSTANT; replacing
it by a parameter-dependent algebraic basis requires the corresponding
connection derivatives. For r=0 the determinant is 1.

For completeness, the characteristic-zero differential-field criterion
follows by induction. Divide by a nonzero first function; the Wronskian
reduces to that of the derivatives of the remaining ratios. A zero
Wronskian gives a nontrivial CONSTANT relation among those derivatives
by induction. Integrating the relation adds a constant and gives a
constant relation among the original ratios. The reverse implication
is immediate. The constants of an algebraic function field over C(s)
are C: on its connected smooth curve a meromorphic function with zero
differential is constant.

If f is eventually pure, its unit Cauchy period is a finite Laurent
polynomial, so (11) for Q=1 forces c_0 in K_unit. Therefore, under the
semisimple fibre1 hypothesis,

    W_unit!=0 => c_0=0 => ALL original markers safe.    (15)

A full failure in this case must have W_unit=0 and a NONZERO ACTUAL
rational vector c_0 in K_unit. Nonzero K_unit alone proves no failure.
This is an all-degree sufficient criterion, not a universal theorem of
unit cyclicity on either exceptional algebra.

## 6. Genuine SL2 control: the unit hides a marked logarithm

On SL2 take f=xw=1+yz and Omega=dx/x wedge dy wedge dz on the x-chart.
For f=s outside {0,1}, x and y are nonzero independent coordinates;
F_s=(C*)^2, with w=s/x and z=(s-1)/y. The relative residue in (4) is

    eta_1=-dx/x wedge dy/y.                            (16)

Its standard torus period is -(2*pi*i)^2. In Hopf coordinates on SU2,
x=sqrt(u)e^(i theta), y=sqrt(1-u)e^(i phi), 0<=u<=1,
direct substitution gives Omega=-du wedge dtheta wedge dphi. Hence
Z=(2*pi*i)^2 for the corresponding compact orientation, and Haar makes
u=f uniformly distributed on [0,1]. Thus

    C_1^f(s)=log(s/(s-1)).                             (17)

The primitive Tate rank is 1 and its unit residue is nonzero. Formula
(11) gives projected c_0=-[torus], c_1=+[torus], with the same orientation.

Now put P=f^2 and r=v^2. Its full generic fibre is TWO tori. On its two
branches f=+v,-v the residue vectors for the original markers 1 and f are

    rho_1=(-1/(2v), +1/(2v)),
    rho_f=(-1/2, -1/2).                               (18)

The two rows in (18) separate the full graph space, but the unit
Wronskian is zero. Exact compact periods are

    C_1^P(r)=1/(2v) log((v+1)/(v-1)),
    C_f^P(r)=1/2 log(r/(r-1)).                         (19)

Near r=0, the first expression is a half-power term plus a holomorphic
series and has NO integer-exponent logarithm. The second has integer-log
coefficient +1/2. Fibre monodromy exchanges the two tori and is semisimple;
the projected actual class is (-1/2,-1/2) in this torus basis, exactly as
(11),(18) predict. The rational projection is not an integral class.

This is an ACTUAL three-dimensional Haar control, not an abstract module.
It is NOT pure: h(P^m)=1/(2m+1), and h(f P^m)=1/(2m+2). It therefore
does not contradict MZ or assert an exceptional failure. It does prove
that unit log-freeness alone is insufficient, and that primitive unit
cyclicity is NOT preserved automatically by polynomial recomposition.

## 7. Unit cyclicity is not universal even on the actual pure safe locus

The [known exceptional biform family](UNBOUNDED_SAFE_BOUNDARY_HOMOLOGY_ON_EXCEPTIONAL_QUOTIENTS.md)
has effective left G_m weight 1, is pure and all-marker safe, and is
algebraically trivial over C* by that action. Its Tate rank satisfies
r>=2k/N-1, for all the previously proved admissible k. The action
preserves Omega and multiplies f by s. Consequently, in its induced
constant flat graph frame,

    rho_1(s)=s^-1 rho_1(1).                            (20)

Its coordinate functions span at most ONE dimension over C. Therefore

    dim_C K_unit >= r-1 >= 2k/N-2.                    (21)

This is UNBOUNDED on BOTH actual PURE SAFE exceptional loci, while the
actual c_0 is zero. Thus universal nonzero unit Wronskian on the entire
pure locus is false. This is an obstruction to a proposed shortcut, not
a new one-weight safety-sector closure or a counterexample on a narrower
unresolved candidate locus.

The complete remaining task stays: control the ACTUAL vector in (11),
now through one finite algebraic residue system, and handle every
nonsemisimple fibre1 case. No large degree/support/prefix campaign,
automatic purity reflection, arbitrary-ring topological promotion or
external-review gate is introduced. Both whole exceptional MZ cases
remain OPEN. The full marked resolvents may still be transcendental;
algebraic logarithm coefficients do not meet the retained hypothesis
that the resolvents THEMSELVES are algebraic.
