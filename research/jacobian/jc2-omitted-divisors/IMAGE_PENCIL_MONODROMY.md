# Geometric monodromy and every omitted point of the nonlinear image pencil

14 November 2026. PF2OMIT18 applies to the four oriented residual
degree-six gap21 rows of PF2OMIT13. It determines the geometric
generic monodromy over the pencil defined by the known image equation
h=H_C, and proves that h composed with the original Keller map has
geometrically integral generic fibre. It also gives the exact genus,
puncture and sheet-loss ledger for that fibre. Actual B,H equations,
the global cover and the original source identification remain absent.

## 1. The function-field setup and horizontal second image

Let F be the hypothetical Keller map, L=C(U,V) the target field,
K=C(x,y) its degree-six source field, and M/L its Galois closure.
By PF2OMIT8, and its nonimmersed extension in PF2OMIT13,

    Gal(M/L)=S6, K=M^S5.                              (P1)

Put kappa=C(h), h=H_C. PF2OMIT16--17 show that the projective
normalization of a generic h-level has genus12 and one infinity point.
The infinity-root proof works also over the algebraic closure of
C(h), so this level is geometrically integral. Thus L/kappa is regular:
kappa is algebraically closed in L, and L is linearly disjoint from
every algebraic extension of kappa.

The distinct second ramified image B meets C by PF2OMIT17. The function
h restricted to B's A1 normalization is therefore nonconstant: a
constant zero would give B=C, while a nonzero constant would make
B disjoint from C. Write

    b=deg(h restricted to B^nu)>0.                    (P2)

It is the total finite C/B intersection multiplicity. In particular
b>=6; if deg B=5j with j odd then b>=9j, and for j1 the sharper bound
is10. The upper bound b<=30j remains valid. Generic h-levels meet B
in exactly b distinct smooth points, avoiding all special image values.

If the optional unramified omitted image H is distinct from C, the
same argument makes h restricted to H nonconstant. Its degree counts
its generic level intersections as well. If H=C, it lies over h0 and
is absent from a generic nonzero level. This is a vertical omitted
prime, not a retained generic point.

## 2. Inertia fixes the algebraic constants over the pencil

Let a be the algebraic closure of kappa in M, finite over kappa.
It is invariant under S6 because S6 fixes kappa. Hence

    G_geom=Gal(M/L*a) is a normal subgroup of S6.      (P3)

Here it is the geometric generic monodromy group, after extending
kappa to its algebraic closure and choosing a component of the
Galois closure. This description uses the regularity of L/kappa.

At the generic point of B the target divisorial valuation is trivial
on kappa. Indeed h restricted to B is nonconstant, so no nonzero
rational function of h vanishes identically there. Any extension of
this valuation to M is also trivial on the algebraic field a: a
valuation extending the trivial valuation across a finite algebraic
extension is trivial, as follows at once from a minimal polynomial.
Consequently a injects into the residue field at that prime.

An inertia element fixes every residue, so it fixes every element of a.
If it moved c in a, its difference from c would have positive valuation
while being a nonzero element of the trivially valued field a. This
is impossible. Thus the B inertia belongs to G_geom. By normality,
G_geom contains its full S6 normal closure.

In rows A,B,N this inertia is a single transposition. Its normal
closure is S6, so G_geom=S6, L*a=L, and regularity gives a=kappa.

In row C it is a single3-cycle. Its normal closure is A6, so
A6 is contained in G_geom. The sign field from PF2OMIT8 is

    L(sqrt(H_C))=L(sqrt(h)).                          (P4)

The square root is algebraic over kappa and lies in a. Therefore
G_geom is contained in A6, giving equality. Since [L*a:L]=2 and
L/kappa is regular, a=kappa(sqrt(h)). We have proved the exact table:

| Oriented row | Algebraic constants of M over kappa | Geometric generic group | Geometric Galois-closure components |
|---|---|---|---:|
| A,B,N | kappa | S6 | 1 |
| C | kappa(sqrt(h)) | A6 | 2 |

The optional unramified omitted prime contributes no inertia to this
argument. Its effect on the source curve is recorded in Section4.

## 3. The source generic fibre stays integral in every row

The algebraic closure of kappa in K is K intersect a. In rows A,B,N
it is kappa. In row C, S5 contains odd permutations, so it acts
nontrivially on sqrt(h); its fixed part in a is again just kappa.
Thus K/kappa is regular in all four rows. Equivalently the source
polynomial h composed with F is a primitive fibration with geometrically
integral generic fibre. It does not acquire a quadratic source factor
from the sign field.

Over an algebraic closure of kappa, the degree-six curve cover has
stabilizer S5 in S6, or A5 in A6. Both actions are2-transitive, hence
primitive: a block containing two letters must contain all letters
by the transitivity of a point stabilizer on the other five.
Consequently this geometric generic curve cover has no proper
intermediate degree-two/degree-three field either. This is a statement
about the nonlinear pencil's generic curve fields, not just the
already established surface-field primitivity.

The two Galois-closure components in row C must not be mistaken for
two source components. The source generic fibre is integral; its
geometric monodromy is the transitive degree-six action of A6.

## 4. Exact genus, punctures and omission on a generic h-level

Work over an algebraic closure of kappa. Let Y be the smooth affine
target generic level of h, whose compactification has genus12 and one
point at infinity, so chi(Y)=-23. Let Gamma be the smooth projective
normalization of its degree-six pullback. Its connectedness is now
proved in Section3. Write g_h for its genus and n_h for its number of
geometric points over the unique infinity point. Then1<=n_h<=6.

This n_h is specific to the nonlinear h-pencil. It is not the
infinity count of a general target line, nor automatically the number
of original boundary divisors on a compact surface.

Put e=e_S, equal to2 in A,B,N and3 in C. The first ramified prime R
has image C=h0 and does not meet this generic nonzero level. The
second prime contributes b omitted points, each with length and
ramification index e. The isolated cusp in B or N is avoided. Put
u for the number of additional unramified omitted points on the level:

| Row/grouping | u |
|---|---:|
| A, H distinct | deg(h restricted to H^nu), a positive integer>=6 |
| A, H=C | 0 |
| A, H=B | b |
| B,N,C | 0 |

If H=B, its points are distinct from the omitted ramified points,
despite sharing their target values. The full generic affine ledger is

    removed points=b+u,
    lost sheets=e*b+u,
    finite different=(e-1)*b.                        (P5)

Infinity contributes different6-n_h. Riemann--Hurwitz, including the
genus of the target curve, therefore gives

    2g_h-2=6*(2*12-2)+(e-1)*b+(6-n_h),
    2g_h=140+(e-1)*b-n_h.                            (P6)

The actual source generic fibre has total puncture count and Euler
characteristic

    r_h=n_h+b+u,
    chi_h=2-2g_h-r_h=-138-e*b-u
         =6*chi(Y)-(lost sheets).                    (P7)

The normal affine pullback before deleting its omitted points has
Euler characteristic-138-(e-1)*b. Thus identity inertia of the extra
u points does not remove them from the puncture or Euler ledger.

The complete necessary numerical genus possibilities are

    e2: g_h in ceil(b/2)+{67,68,69}, n_h congruent to b modulo2;
    e3: g_h in b+{67,68,69}, n_h even.               (P8)

Every such generic source normalization has genus at least70 in e2
rows and73 in row C, using b>=6. For a quintic second image b>=10,
the bounds strengthen to72 and77. At the generic high-contact quintic
incidence b12 they give73--75 and79--81. These are necessary numerical
options, not constructed covers or assertions that all options occur.
The genus3--5 and5--7 values for a general linear target pencil are
different and remain correct; the nonlinear h-pencil cannot replace it.

For a numerical control, b12,n_h2 gives g_h75 in A and g_h81 in C.
With H=B the A puncture count is26; in C it is14. Both Euler values
are-174. With H=C, A instead has14 punctures and Euler-162.
These controls are not global Keller realizations.

The infinity permutation has sign(-1)^((e-1)*b), since the finite
branch cycles and the product of genus commutators determine its sign.
The [replay](verify_image_pencil.py) checks all partitions of6 against
this parity and the group, while retaining them as numeric possibilities.

## 5. What this supplies and what remains

Any actual joint image/cover proposal must satisfy the geometric group,
constant-field, integral-source-fibre, primitive-curve-cover and full
omission constraints above. The constant quadratic field occurs only
in row C's Galois closure; it never supplies a quadratic source field.
The target h-level itself is known from C, while its degree-six cover,
its actual infinity partition, and the b,u values still require B,H
and the global normalization.

The special level h0 contains C and needs its own fibre analysis. In
particular the whole omitted first prime and any H=C prime lie over
that level. The generic formulas do not account for their central
fibre factors, or for collisions at other exceptional h-values.
No actual global cover or source A2 identification, no unrestricted
singleton-zero properness or whole-r7 conclusion follows here.
