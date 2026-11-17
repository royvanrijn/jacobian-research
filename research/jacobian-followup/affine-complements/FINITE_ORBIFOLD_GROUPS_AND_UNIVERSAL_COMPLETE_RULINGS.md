# Finite orbifold groups and the complete ruling behind noncyclic deck symmetry

PF2RULED134. All actual complex polynomial plane Keller degrees and all
genuine finite intermediates. This uses the classical affine-line ruling
and sphere-orbifold inputs already adopted in PF132 at their GENERAL
finite-fundamental-group scope. No Euler<=2 hypothesis is used to obtain
the group containers. The full Wright/HC4 finite-target, trace-energy and
EMPTY ORIGINAL fibre problems remain OPEN.

The main reduction is global: every genuine intermediate's finite pi1,
and every subgroup of the actual generic deck group, is cyclic, dihedral,
A4, S4 or A5. Noncyclic cases have a universal complete ruling and must
retain a fibre whose component multiplicities, AFTER division by their
gcd, are ALL at least2. This supplies structure for the remaining rigid
branch; it does not exclude that branch or supply a finite Keller target.

## 1. Retained geometric and classical inputs

For an actual original plane source X=A2, field extension L/F, and
genuine intermediate N->M->A2, let V be the ENTIRE etale locus of M.
PF128/129 and KU0/1 give V smooth affine, negative logarithmic Kodaira
dimension, with finite pi1 of order m dividing t=[L:C(M)]. Its target
coframe gives the exact nowhere-zero form dP wedge dQ.

For ANY subgroup G<=Aut_F(L), put T=U_N/G. HC4MYLA0/1 gives the same
surface hypotheses for T and pi1(T)=G; U_N is its simply connected
finite etale cover. T need not be the ENTIRE etale locus of the finite
normalization in L^G. All original omitted primes are retained.

We use the following precisely credited retained inputs:

* Dubouloz--Kishimoto, arXiv1408.1328v1, Section1.2/Theorem2, PDF3:
  negative log Kodaira dimension supplies an A1-fibration over a smooth
  curve. Over C no constant-field extension is needed.
* Miyanishi, *Affine pseudo-planes and affine pseudo-coverings*,
  Oberwolfach Report19/2005, printed1111/PDF39, Lemmas2.3--2.4:
  pi1 is the base orbifold group, using WHOLE-FIBRE gcd multiplicities.
  Finite groups over A1 are cyclic; over P1 the three-multiple-fibre
  signatures are spherical. Lemma2.4 does NOT require irreducible fibres.
  That extra condition begins only in the subsequent Definition2.5.
  The full relevant page was reread and visually checked for this use.
* Thurston, Chapter13, Theorem13.3.6 and adjacent table, printed312/PDF18:
  good/bad orientable sphere orbifolds and their spherical covers with
  exact local orders, as explicitly adopted in PF132.
* HC4MYML0's retained Dubouloz2002.12838v1 local-structure package:
  reduced components are disjoint A1s, and a multiplicity-one selection
  gives a Zariski A1-bundle. Its source/proof-sketch boundary is retained.

These are adopted classical results, not independent classifications or
proofs supplied by the checker. The orbifold meridian/gcd argument is
also written in the credited coupling packet and PF132. Here their
consequences are applied beyond its low-Euler specialization.

## 2. All finite fundamental groups lie in the Möbius list

Let Y be either V or T above, and choose ANY A1-fibration rho:Y->B.
Finite pi1(Y) surjects onto pi1(B), so B=A1 or P1. Let a special fibre
have reduced components C_ij with multiplicities mu_ij, and put

    m_i=gcd_j(mu_ij).

The contractible generic fibre gives the exact orbifold presentation.
Over A1 it is a free product of cyclic groups of orders m_i; finiteness
leaves at most one nontrivial factor. Over P1 it is

    <gamma_1,...,gamma_r | product gamma_i=1, gamma_i^m_i=1>.

Zero or one cone gives the trivial group; two cones give C_gcd(m1,m2),
including bad unequal spindles. In the noncyclic finite case the base
orbifold is good spherical, with exactly one of these signatures:

| Signature | Group | Order |
|---|---|---:|
| (2,2,n), n>=2 | D_(2n), meaning dihedral of order2n | 2n |
| (2,3,3) | A4 | 12 |
| (2,3,4) | S4 | 24 |
| (2,3,5) | A5 | 60 |

These are the orientation-preserving spherical symmetry groups. Thus all
finite groups here embed ABSTRACTLY in PGL2(C), including cyclic groups.
For noncyclic groups Section3 provides the actual faithful action on the
base of a universal ruling. No such base action is asserted from a bad
cyclic orbifold's nonexistent manifold cover.

Consequences in every original degree include:

* An ABELIAN generic deck subgroup is cyclic or C2 x C2. In particular
  C2^q, q>=3, is impossible as a DECK subgroup for ANY monodromy action.
* Every odd-order finite pi1 here, and every odd-order deck subgroup, is
  cyclic. Thus every genuine intermediate with odd relative degree t has
  cyclic pi1, including all intermediates of an odd-degree original map.
* Quaternion groups, noncyclic odd p-groups and groups outside the list
  cannot occur in these roles. These are NOT exclusions of the full
  monodromy, its arbitrary subgroups, or its composition factors.

Combine with PF132/129: a noncyclic pi1(V) forces chi(V)>=3 and a genuine
branch costs support at least2, so b=deg(M->A2)>=5. Since m divides t,
the original degree obeys d=b*t>=5*m. For deck subgroups the same bound
is PF133. The A4, S4 and A5 cases therefore require original degrees at
least60,120 and300 respectively; dihedral order2n requires at least10n.
These are necessary degree bounds, not inversion in the smaller degrees.

## 3. The universal cover has a complete ruling, with exact fibre data

Suppose Gamma=pi1(Y) is NONCYCLIC, and let g=|Gamma|. Every A1-fibration
of Y must have base P1 and a good spherical orbifold. Let

    tau:P1'->P1

be its universal orbifold cover, of degree g and local orders m_i.
Normalize the base change to obtain W->Y. The generic A1 field is
geometrically integral, so W is connected. At a component with
mu_ij=m_i*a_ij, the exact-order base change is divisorially unramified.
Purity over smooth Y makes W->Y finite etale everywhere. It has degree
g, hence is the UNIVERSAL cover of Y. Its ruling W->P1' is complete-type.

There is an exact, generally NONREDUCED fibre formula. Over each of the
g/m_i base points above the i-th special value, the fibre has

    m_i copies of EACH original component,
    with multiplicity a_ij=mu_ij/m_i on each copy.       (1)

Indeed W->Y is etale and its restriction above C_ij=A1 is the trivial
g-sheet cover, so it has g separate A1 lifts in total. The base deck
group permutes the g/m_i base points transitively, giving m_i lifts at
each, and pullback of a local base parameter divides the multiplicity
by m_i. All points are retained; the roots do not force the a_ij to1.

In particular a whole multiple irreducible fibre produces m_i REDUCED
components over each of those base points, while a mixed fibre
m_i*(2C_1+3C_2) produces m_i copies of multiplicities2 and3. The latter
survives this cover. Do not replace (1) by either universal irreducibility
or universal reducedness.

Let s_i be the number of original reduced components. Euler integration
and Riemann--Hurwitz give the useful exact check

    chi(Y)=2+sum_i(s_i-1),
    chi(W)=2+sum_i (g/m_i)*(m_i*s_i-1)=g*chi(Y),
    2=g*(2-sum_(cone i)(1-1/m_i)).                    (2)

The finite set in the middle sum contains all degenerate fibres, with
m_i=1 permitted. The last sum contains only actual cone points. This
keeps reducible ordinary fibres as well as reducible multiple fibres.

For Y=T=U_N/G, uniqueness of its connected degree-g universal cover
identifies W with the ACTUAL U_N. The ruling

    U_N -> P1'

is G-equivariant, with faithful base action the corresponding dihedral
or platonic action. Thus a noncyclic generic deck group forces an
ACTUAL complete-type A1-ruling of the full etale completion, not merely
a conjectural ruling of some substituted birational surface. HC4MYLM1
still makes U_N additively rigid: complete-base rulings are not Ga quotients
over A1. This is additional structure, not a contradiction to rigidity.

## 4. Exact volume forces a missing gcd component in every such ruling

The universal cover W is smooth affine and carries the pulled-back
nowhere-zero exact form. HC4MYML0 therefore requires some fibre of
W->P1' to have NO multiplicity-one component. By (1), the corresponding
fibre of Y has

    mu_ij/m_i >=2 for EVERY component j.               (3)

Equivalently, NO component of that fibre has multiplicity equal to the
whole-fibre gcd. This holds for EVERY chosen ruling of a noncyclic Y,
not just for a conveniently selected ruling. This is the precise
obstruction to applying the tempting gcd-component chart shortcut.

Such a fibre is necessarily reducible: a single multiplicity equals its
own gcd. If chi(Y)=3, the Euler formula in (2) leaves exactly ONE reducible
fibre, with precisely TWO components. It must have the form

    m*(a*C_1+b*C_2),  gcd(a,b)=1, a,b>=2.              (4)

All other fibres are irreducible. The bad fibre can be a cone fibre
(m>1) or an ordinary orbifold value (m=1). We do not assume the three
whole multiple fibres to be irreducible. This sharpens PF133's unresolved
index-five profile to an explicit mixed-fibre escape, uniformly rather
than by enumerating multiplicities. The profile itself remains OPEN.

## 5. Checks and remaining scope

The checker pins the retained primary PDFs and complete relevant source
records, verifies the spherical order/Euler identities, generalized
split-fibre formula and the gcd-component obstruction, and supplies tiny
permutation controls distinguishing dihedral, quaternion and platonic
groups. Written arguments and the explicit classical imports establish
the unbounded conclusions; finite samples do not classify arbitrary groups.
Registered predecessor mathematical checkers are not rerun.

The new containers concern finite surface fundamental groups and generic
DECK groups. Natural A/S and full-AGL MONODROMY families remain open.
The surviving cyclic, dihedral and platonic deck cases are not constructed
or excluded here beyond PF133's index/support bounds. No automatic common
quadratic refinement, affine ruling, source-plane action, polynomial target,
trace integrality or full HC4/JC2 closure is asserted. Exact-volume mixed
fibres remain the geometric escape; the actual finite-target/coframe data
are still essential. Formal verification, independent replay, external
review, priority and full-goal completion remain unclaimed.
