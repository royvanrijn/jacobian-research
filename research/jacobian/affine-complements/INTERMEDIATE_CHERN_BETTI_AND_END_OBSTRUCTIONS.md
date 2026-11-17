# Intermediate Chern classes, exact affine Betti profiles and end exclusions

PF2RULED124, 16 November 2026. Let N be the ORIGINAL finite normalization
of an ACTUAL nonautomorphic complex PLANE Keller map, U=A2 its WHOLE
marked source, and retain ALL omitted primes. For EVERY nontrivial
finite intermediate N->M->A2, the rational canonical Chern class of M
is nonzero. Critical divisor classes are rationally independent in
H2-cohomology. The minimal resolution has nonzero rational Chern class
and its second Betti number is at least the number of critical curves
plus the number of resolution curves.

Consequently EVERY such intermediate, singular or smooth, has no S3 or
negative ADE oriented end. An L(p,1) end forces M smooth with exactly
one critical curve and b2=1. In a binary affine quotient of degree b,
the exact b2 is 3b/4-1 (immersive image) or b/2-1 (one ordinary cusp).
Thus an L(p,1) end is excluded for ALL b>=8 and for the degree4
immersive case. The degree4 one-cusp case and the general finite algebra
remain OPEN. No intermediate is assumed to be a plane.

## 1. Retained inputs and what is added here

HC4MYKN1 already proves NON-TORSION Weil canonical class for EVERY
nontrivial finite intermediate, in arbitrary source dimension. Its
torsion-canonical, factorial and global complete-intersection exclusions
are RETAINED RESULTS, not new claims here.

UADANVOL1 already constructs exact volume on generalized Danielewski
surfaces. HC4MYLA1 already pins the known cyclic pseudo-plane controls
W_k={x^2*y=z^k-1}, their actual open A2 charts and their free cyclic
quotients T_k. They give surjective etale A2->T_k, scalar units, finite
fundamental group, rational acyclicity and exact nonvanishing volume.
These controls are NOT new constructions here. They are not full finite
intermediate normalizations over the required plane target. In particular
the full canonical obstruction does not apply to an OPEN etale quotient
as though it were that entire finite normalization.

The additional work here is a TOPOLOGICAL strengthening for the current
plane intermediates, the exact binary affine Betti calculation, and the
resulting transfer of complete end exclusions. It uses PF123's newly
proved local geometry before transferring any filling argument.

## 2. Independent critical classes and nonzero rational Chern class

Write h:N->M for the finite degree-t map and pi_M:M->A2 for degree b>1.
By PF123 the critical curves E_1,...,E_s in M are disjoint smooth A1s;
singularities are cyclic quotient-axis points. Thus M is locally
Q-factorial and Q-Gorenstein. Choose a common positive m making all
m E_j Cartier. Every E_j has a preimage original omitted prime D_i,
and all its prime preimages are omitted: transverse indices multiply,
so ramification downstairs cannot meet the etale whole source U.

Let rho:Y_N->N be the minimal resolution. PF111 gives

    Pic(Y_N) --c1--> H^2(Y_N;Z)
       = direct_sum Z*[strict D_i] + direct_sum Z*[exceptionals].

The divisor of (h*rho)^*O_M(m E_j) is supported on this boundary, with
positive coefficient m*e(D_i/E_j) at every strict D_i mapping to E_j.
For different j, those D_i-coordinate sets are disjoint and nonempty.
Therefore these pulled-back classes are rationally independent. By
naturality of c1, the original classes c1(O_M(m E_j)) are independent
in H^2(M;Q) as well. No general injectivity of h on all Picard groups
or an intermediate divisor-basis theorem has been assumed.

The pulled-back target volume represents the canonical Weil divisor

    K_M = sum_j (e_j-1) E_j.

Purity and trivial connected finite etale covers of A2 make s>=1.
Every coefficient is positive. Defining c1(K_M) rationally by any
Cartier multiple, the independence proves

    c1(K_M)!=0 in H^2(M;Q),   b2(M)>=s.             (1)

In particular M is not a rational homology plane. The older algebraic
torsion-canonical exclusion is credited above. It also excludes every
normal finite factorial cover of M: if f:X->M is finite surjective and
X factorial, f^*(m K_M) is principal. Pushing the principal divisor
down by the function-field norm makes deg(f)*m K_M principal, contrary
to HC4MYKN1. Planes and finite quotients of factorial affine varieties
are therefore already excluded as full intermediates; an etale dominant
map from A2 is not a finite factorial cover.

## 3. Fundamental group and resolution Chern class

Every loop in N_reg moves off its boundary curves into U=A2 and
contracts there. Filling the isolated normal singular points only
quotients the fundamental group, so pi1(N)=0. Remove from M its
singular points and the branch/exceptional-value set of the finite h.
Over this dense open V, h is a connected finite etale cover of degree t.
Its fundamental subgroup H has index t in pi1(V). The map pi1(V)->pi1(M)
is surjective, and H maps trivially because its map factors through N.
Hence

    pi1(M) is finite,   |pi1(M)| divides t.          (2)

Let rho_M:Y_M->M be the minimal resolution, with a exceptional curves.
Cyclic quotient links are rational homology spheres. Excision for
removing the finite singular set therefore makes

    H^2(M;Q) -> H^2(M_reg;Q)

an isomorphism (the adjacent relative groups are H1 and H2 of those
links, both zero over Q). If c1(K_Y_M) were zero, its restriction to
M_reg would make (1) zero, impossible. Thus

    c1(TY_M)!=0 in H^2(Y_M;Q).                     (3)

The pulled-back critical classes and the a exceptional classes are
independent: restriction to M_reg first kills all critical coefficients
by (1), and negative definiteness of each exceptional chain then kills
the remaining coefficients. Consequently

    b2(Y_M)>=s+a.                                 (4)

Replacing a quotient cone by its simply connected resolution plumbing
kills the same link subgroup in van Kampen. Thus pi1(Y_M)=pi1(M), so
b1(Y_M;Q)=0. None of these assertions identify Y_M with a plane.

## 4. Transfer the complete ADE and circle-bundle end arguments

A large core of Y_M is a minimal relatively compact strictly pseudoconvex
surface: every compact complex curve maps to a point in affine M, hence
is exceptional, and minimal resolution has no holomorphic (-1) curve.
Use EXACTLY the Bogomolov--de Oliveira/ Epstein/Gray imports audited in
PF109. They give a small Stein deformation V of this core preserving
its topology, rational Chern class and boundary contact isotopy class.
Thus V has nonzero rational c1 by (3), b1=0 and b2>=s+a. The resulting
Stein symplectic form is exact and the filling minimal. The original
resolved complex structure is not assumed exact, and the finite map
is not claimed holomorphic after deformation.

For a negative ADE oriented link, retain the precise contact uniqueness
and Ohta--Ono filling imports of PF108--109. They apply to ANY minimal
Stein filling of that canonical contact link; simple connectedness of
our starting core is not an extra hypothesis. The filling has the
negative Cartan topology of rank L, while the contact theta is L-2.
Gompf's formula then reads

    theta=c1(TV)^2+L-2.

It forces c1(TV)^2=0, so negative definiteness forces c1(TV)=0 rationally,
contradicting (3). This excludes EVERY negative ADE end of EVERY
nontrivial finite intermediate, including singular M. Oriented end
identification is an explicit premise; determinant or homology alone
does not supply it.

PF110's OSS planar obstruction likewise transfers: if this contact is
planar, its INTEGRAL boundary c1 must be nonzero. Rational boundary c1
zero is insufficient. The unique positive tight S3 contact is planar
with integral c1=0; hence S3 is also excluded. All negative finite-SU(2)
quotient ends are thus excluded in the stated orientation.

For L(p,1), retain the published PVHM rank bound of PF110. Its p=4
factorization calculation needs b1=0 to turn Euler characteristic into
b2=length-3; (2) supplies exactly this rational vanishing. It does not
require pi1(V)=0. Thus b2(V)<=1 in every case. Equations (3)--(4) force

    s=1, a=0, M smooth, b2(M)=1.                   (5)

This counts CRITICAL curves in M, not all omitted primes of the original
N. No primitive marking of E_1 in Pic(M) has been proved; the stronger
PF110 index/parity conclusions for an original sole-boundary plane
completion are not transferred without that missing marking.

## 5. Exact binary affine intermediate Betti profiles

Let Q be a binary affine block quotient of degree b=4r. PF123 gives
smooth M, r critical curves, one branch C with normalization A1, and
kappa=0 or1 ordinary cusps. Let v count its double normalization values,
including arbitrary tangency orders. Then

    chi(C)=1-v,
    chi(C minus its singular values)=1-2v-kappa.

The finite fibre cardinalities of M->A2 are b off C, b-r at a smooth
C value, and b-2r at a double value OR the ordinary cusp. At a double
value the critical factors exhaust all sheets; at a cusp they have
length3 with r remaining etale points. Euler integration therefore gives

    chi(M)=b*(1-chi(C))+(b-r)*(1-2v-kappa)
                         +(b-2r)*(v+kappa)
          =b-r*(1+kappa)=r*(3-kappa).               (6)

M is a smooth affine Stein surface, so has the homotopy type of a real
two-dimensional CW complex. Equation (2) gives b1=0 rationally. Hence

    b2(M)=r*(3-kappa)-1
      =3b/4-1 (immersive), or b/2-1 (one cusp).      (7)

Finite H1 torsion is allowed; no integral acyclicity is inferred. The
formula is independent of the number and orders of the double contacts.
It uses proper FINITE pi_M, not the nonproper etale source restriction.

For an L(p,1) end, (5) and (7) now force r=1 and kappa=1. Thus ONLY
b=4 with one ordinary cusp survives this end test. ALL b>=8 binary
affine intermediates and ALL immersive binary affine intermediates
have NO L(p,1) end. This is an all-degree end-class exclusion, not
exclusion of every binary affine quotient or of every quotient degree4.

## 6. Remaining task and the HC4 missing-point interface

The full normalization M is different from its open etale locus. The
retained pseudo-plane controls have exact unit volume and surjective
etale plane charts; their finite completions cannot satisfy (1) with
the required target plane. Rebuilding their known formulas supplies no
new closure. A further unmarked exact-volume argument cannot replace
the actual finite map, global canonical class and embedded multiplication.

PF123's double values give EMPTY ORIGINAL source fibres even for
imprimitive monodromy and higher original indices. There is a further
precise pencil restriction for the shared HC4 route. Translate such a
value to0 and take a line through0 transverse to both branch tangents.
On the normalization of its inverse image in N, the radial scaling
coordinate has a zero of order e_i at EACH of the two points belonging
to every visible original D_i. Its complete zero profile is therefore

    (e_1,e_1,e_2,e_2,...),  all e_i EVEN, sum 2e_i=d. (8)

Indeed the transverse local factors are ordinary index-e_i folds, and
PF123 says these points exhaust the complete fibre. A generic pencil
line is transverse to the remaining branch locus and avoids its other
singular points; one may work componentwise if its inverse curve is
disconnected. At the block level, the two
branch transvections have complementary supports, so the radial
meridian is a FIXED-POINT-FREE affine INVOLUTION with profile 2^(b/2).
It is NOT necessarily a translation: complementary affine transvections
can have the same covector and different direction vectors. Nor must
its original lift be an involution; its cycles are exactly (8).

For the basepoint-free pencil [P:Q] after translating the missing value,
(8) says the radial function Q on a generic finite-normalized pencil
curve has an entirely even zero divisor over0. This does NOT make Q
a square in the function field: pole orders at infinity and divisor
classes have not been controlled. The shared polynomial zero-free
Euler lift with divergence2 is retained context, not a new theorem or
an automatic completeness assertion here. No proof excluding this
full-sheet empty-fibre/even-zero profile is claimed. Natural
alternating/symmetric quotients, binary affine
global compatibility, the surviving degree4 one-cusp circle-bundle
intermediate, invisible original branch images and the actual full
finite algebra B=C[x,y] remain OPEN. Sealed predecessors, all original
primes, shared analytic/multiplication priorities and published quartic
PLANE retirement are preserved. No historical priority, new independent
proof of the imported filling theorems, fresh PDF visual audit, formal
verification, external review or full goal completion is claimed.
