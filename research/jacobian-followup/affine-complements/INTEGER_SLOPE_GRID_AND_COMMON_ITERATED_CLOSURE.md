# Integer slope grids and common-pencil iterated closure

15 November 2026. PF2RULED93 excludes the ENTIRE common-pencil class
of genuine positive iterated-torus quartic ends covered by PF82's
complete peripheral-state reduction, at every finite length. The
whole smooth source, entire affine-line boundary, actual preferred
peripherals, genuine stages and simple-meridian full-S4 final
monodromy hypotheses of that reduction remain in force. This is not a full quartic or
Jacobian-conjecture theorem. Complementary-pencil tails remain open.

## 1. A degree-free grid of missing interior slopes

Let P be an oriented irreducible manifold with one incompressible
torus boundary and actual primitive basis(m,L). Suppose, INTEGRALLY,

    H1(P)=Z*L plus T, [m]=0.                         (1)

T may have torsion or free summands. For EVERY integer a, put
gamma_a=a*m+L. In the gamma_a filling, the core is represented,
up to sign and multiples of gamma_a, by m, because their boundary
intersection has absolute value1. Its class is thus integrally0.
PF90's null-homologous-core lemma makes gamma_a lie OUTSIDE the
L-space INTERIOR of P. If the filling has positive b1, it is already
not an L-space. Otherwise that lemma uses the retained primary
Floer-simple-core implication and Ni--Wu Proposition5.1: an
integrally null-homologous Floer-simple knot in an L-space has
genus0, contrary to boundary incompressibility of its exterior.
No assertion about all null-homologous knots is made.

PF87 proves(1) for ANY filled connected cyclic n-pattern of winding
coprime to n, with T the actual cyclic branched-companion homology.
Thus the grid assertion is degree-free whenever the additional
irreducibility/incompressibility hypotheses hold. In particular an
actual C2 tail of odd winding w>=3 has P disk(w,w) for even k, or
disk(w,2) for odd k. It satisfies those hypotheses, with T=Z/w or0.

Consequently, if a nonsolid incompressible graph exterior E with
b1(E)=1 has ANY gamma_a outside its L-space interior, gluing ANY
such C2 tail P to E is not an L-space: the matched gamma_a is
missing from both interiors. This uses Rasmussen Corollary6.1,
without assuming either interior is nonempty. The closed union is
irreducible graph. Its next drilled R exterior is incompressible
graph by the actual annulus construction of PF82/92. A positive-b1
obstruction persists by the integral transfer; otherwise the next
normal filling is non-L-space and PF82/92 propagate exclusion
through EVERY further finite tail. PF74 excludes the common source.

## 2. Finite denominator certificates for a graph root

We newly import Rasmussen arXiv:1511.04413v3 Theorem4.6, with the
FULL Section4.2 conventions and Definition3.7. The pinned HTML and
manifest record this deeper graph-interval input; it is not an
independent proof of the published theorem. Integer Seifert
coefficients are permitted and e0=0 is absorbed as in Section4.2.

For an orientable nonsolid graph root with cone coefficients a_i
and daughter intervals [[g_i-,g_i+]], define

    F_-(N)=-[1+sum floor(N*a_i)
                   +sum(ceil(N*g_i+)-1)]/N,
    F_+(N)=-[-1+sum ceil(N*a_i)
                    +sum(floor(N*g_i-)+1)]/N.         (2)

If a finite target slope x were in its L-space INTERIOR, then the
manifold would be Floer simple with nonempty L-space set. Theorem4.6
would give endpoints y_-=max F_-, y_+=min F_+. Therefore witnesses

    F_-(N1)>=x>=F_+(N2)                              (3)

contradict interior membership by Definition3.7. Finite input
coefficients keep these endpoints finite. When y_-=y_+=x, the
interval is the complement of x, so equality also suffices. If the
L-space set is empty or a singleton, its interior was already empty.
We need no separate classification of when that set is nonempty,
no search for all denominators and no claimed Floer computation
from a homology matrix.

## 3. The actual drilled S3 root, with its normal section fixed

PF90 already excludes ALL nonsolid incoming S3 companion histories
through all common tails. The only solid histories start with the
trefoil and have PF91's odd cyclic2 loops. Write their actual cap
as A*m+L, with A=W^2, and let J be the ordinary positive iterated
torus companion. PF91 gives

    c=2g(J)-1=2A-D, D>=1 odd,
    L(E(J))=[c,infinity],
    B_companion=2A+w_last for all-plus nonempty histories,
    B_companion=2A-w_last after a last-minus loop.     (4)

These are PF91's actual recurrences B_new=w*k. There can be no
solid continuation after a minus loop.

At first simple S3 arrival, r=3s odd, k even, gcd(r,k)=1. Put

    J0=A*r-k/2 !=0, d=k-r*c=eta_arrival+1.             (5)

PF92 excludes every common arrival eta<=0 through all later tails.
We can therefore assume eta_arrival>0: it is even, so d>=3 is odd.
There is at least one later tail, since PF91 closes final S3 arrivals.

Drill the ACTUAL R from the main arrival piece before capping the
solid companion. Its base has three boundaries and cone s. The
fibre has degree1 and the normal R filling is a SECTION: the
branched meridian cycle lifts a primitive section on its boundary,
and R is a regular fibre. This follows from PF81's actual cable
sections, not from a homology recognition of R.

Choose r*u-k*v=1, 1<=v<r. Use PF81's same X,Y,Z sections on the
other boundaries, and U the actual normal R meridian, with a
coherent choice of signs. Its presentation is

    X^s=h^(-v), XYZU=h^j,
    n=h^(-v)Y^r, L1=h^u Y^(-k),
    m=h^(-2v)Z^r, L=h^u Z^(-k/2).                    (6)

Filling U restores EXACTLY PF81's retained relative presentation
XYZ=1. In that filled piece h has infinite order (H1 free rank2),
so j=0. This fixes the integer section shift; it is not merely a
choice of an arbitrary Seifert zero framing. The solid cap gives

    Z^J0=h^beta, beta=2A*v-u.                        (7)

The root Q now has annulus base and orders(s,abs(J0)), with order1
entries absorbed. In boundary basis(h,-U), the normal R slope is0,
the cone coefficients and daughter endpoints are

    a1=-3v/r, a2=beta/J0=2v/r-1/(r*J0),
    x=v/r, z=x+1/(r*d), daughter [[z,x]].              (8)

Signed J0 is essential in a2, even when its order is abs(J0).
For the daughter this follows directly from(6): t*n+L1 has slope
(u-t*v)/(k-r*t). The map is increasing with derivative
1/(k-r*t)^2 and sends c to z, infinity to x. It sends the FULL
ordinary positive L-space interval[c,infinity] to [[z,x]]. The
pole inside that interval accounts for the wrap. Since d>=3,
0<x<z<1, including when v=r-1.

When either cone is genuine, Q has negative orbifold Euler
characteristic, and E_R=Q union E(J) is an actual nonsolid,
irreducible, boundary-incompressible graph manifold with b1=1.
The gluing fibre intersection is r*B_companion-k>0 by the genuine
infinity inequality. Thus(2) applies to THIS rooted decomposition.

## 4. Two witnesses close every nonproduct root

There is one daughter, so(2) simplifies to

    F_-(N)=-[floor(N*a1)+floor(N*a2)+ceil(N*x)]/N,
    F_+(N)=-[ceil(N*a1)+ceil(N*a2)+floor(N*z)]/N.      (9)

The following table supplies(3) for the ACTUAL normal slope0:

| Case | N1 | F_-(N1) | N2 | F_+(N2) |
| --- | --- | --- | --- | --- |
| J0>=2 | r | 1/r | r | 0 |
| J0<=-2 | r | 0 | r | -1/r |
| J0=1, r>=9 | r | 1/r | 1 | 0 |
| J0=-1, r>=9 | 1 | 0 | r | -1/r |

Indeed r*a1=-3v, r*a2=2v-1/J0, r*x=v,
r*z=v+1/d with 0<1/d<1. These prove the first two rows,
and the r entries in the last two. If J0=1, k=2Ar-2 and
v=(r+1)/2; if J0=-1, k=2Ar+2 and v=(r-1)/2.
In both unit cases a2=1 and, at r>=9, -2<a1<-1,
0<x<z<1. Thus F_+(1)=0 in the positive unit case and
F_-(1)=0 in the negative unit case, proving the table.

Section2 puts the normal R meridian outside the L-space INTERIOR
for every nonproduct root, at every solid-history depth and all
remaining arrival parameters. PF92's endpoint tail theorem now
excludes EVERY nonempty later simple tail. We have imposed no
affine-source hypothesis on a proper prefix, and no sign condition
on any earlier prefix degree.

## 5. The product root is closed by the integer grid

The only unhandled roots have s=1,abs(J0)=1, namely r=3,
k=6A +/-2. A last-minus history has B=2A-w, w>=3, so
k<3B forbids BOTH unit-cap values. Hence the history is all-plus,
c=2A-1. The value k=6A-2 has eta=0 and was closed through all
tails in PF92. The remaining value is

    r=3, k=6A+2, eta_arrival=4, J0=-1.               (10)

Both cones have order1, so the actual root is a product collar
and E_R is the ordinary E(J). PF81's product map is
m=-n, L=(k/2)n+L1. Capping A*m+L consequently identifies
the normal R slope projectively with

    gamma=(2A+1)n+L1.                               (11)

For an explicit check directly in the DRILLED presentation, s=1
gives X=h^(-v), and J0=-1 gives Z=h^(-beta). Equation(6) then
gives U^(-1)=h^(-v-beta)Y. Meanwhile(2A+1)n+L1 has exactly
that word, since 3(2A+1)-k=1 and
u-(2A+1)v=-v-beta. Thus(11) tracks the actual normal section;
it is not inferred from a homeomorphism type of the filled end.

The ordinary meridian n has geometric intersection1 with gamma.
The actual normal/preferred pair(m_R,L_R) is an integral basis;
therefore, in that basis, n is projectively a*m_R+L_R for SOME
integer a. We need neither its precise value nor an identification
of the preferred framing with a Seifert zero framing.

The slope n lies outside the L-space interior of E(J): it is
infinity in[c,infinity], or directly use PF90's null-core lemma
for its S3 filling. Section1 puts its matching integer-grid slope
outside the interior of EVERY actual C2 filled tail. Corollary6.1
therefore excludes the first tail, and PF82/92 exclude all further
tails. The proof does not transport the genus endpoint through an
extra annulus and does not assume square-invariance of linking.

## 6. Exhaustion and exact scope

PF82's complete ten-state tree/planar reduction has simple full-S4
state3. Every trajectory either starts there or first enters it
from the following list. State numbers refer to its sealed receipt.

| First simple arrival | All-length common exclusion |
| --- | --- |
| Initial state3 | PF88 |
| Incoming states4/5 (S4) | PF86 |
| Incoming states6/9 (cyclic C3/C4) | PF85 |
| Incoming state7 (A4) | PF87 |
| Incoming state2 (S3) | PF90 nonsolid; PF92 eta<=0; Sections3--5 here |

After arrival every passing stage is exactly an odd-winding cyclic2
tail, by the same complete state reduction. Thus ALL genuine
positive iterated-torus quartic COMMON-pencil end matches in that
scope are excluded, for unbounded parameters and arbitrary finite
stage count. The common source obstruction is PF74's no nontrivial
left-orderable image, with its whole smooth affine-source hypotheses.

COMPLEMENTARY S3 survivors remain: arrivalr3, or all-plus r>=9
with exact k=2Ar +/-2 and later tails. PF89's later A4 ordered
terminal pairs(2,4),(2,8),(4,10), including its order/tail restrictions,
remain open. This theorem does not settle those complementary
sources, degenerate/non-genuine stages, higher-degree packets,
all quartic ends, finite Keller maps, or the full research goal.
KIG28 concerns FULL generic A4 monodromy of a Keller self-map;
it does not identify local/peripheral/intermediate A4 states with
that full action and is not used to close this frontier.

The [checker](verify_integer_grid_common.py),
[receipt](integer_grid_common_receipt.json) and
[primary manifest](sources/integer-grid-common/manifest.json)
retain exact symbolic identities, signed root certificates, actual
integer boundary controls and predecessor seals. The written proofs
establish the unbounded assertions; finite controls do not. Deep
Floer/graph/LO/peripheral-rigidity inputs are imported. No formal
verification, independent full reproof, external review, or novelty
claim is made.
