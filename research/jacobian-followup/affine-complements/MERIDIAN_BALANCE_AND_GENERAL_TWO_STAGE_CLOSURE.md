# Meridian balance and the general two-stage quartic closure

15 November 2026. PF2RULED81 excludes the ENTIRE general genuine
TWO-STAGE QUARTIC boundary-packet class, allowing arbitrary lower-weight
key terms. All coprime first pairs p,q>=2, final winding r>=2 and
ordinary cabling coefficient k>=2 are covered, at every source order,
pole and characteristic depth under PF46/50's actual packet hypotheses.
No representation or source pencil remains in this target class.

PF80 already closed eta0 and odd rk. This note closes its entire
strictly negative EVEN remainder, k<=r*(pq-p-q)-1, by actual marked
end models with k independent of r*(pq-p-q)+1. The new common-source
argument propagates a missing L-space slope through a two-boundary
piece using an integral meridian balance; its trefoil exception has
two explicit interval bounds. No longer coefficient search is needed.

The source lemmas are degree-free under the explicit common or
TWO-line complementary smooth/exact/whole-chart hypotheses. This is
not a classification of arbitrary higher-degree packets. Longer target
sequences, higher degrees and the full goal remain OPEN; no full
quartic Keller exclusion is asserted.

## 1. A general consequence of a meridian balance

Let P be an orientable irreducible graph manifold with two incompressible
torus boundaries T_m,T_n, and primitive slopes m,n. Let A,C be non-solid
graph manifolds with one incompressible torus boundary, glued to those
boundaries by their ACTUAL maps. Suppose:

* alpha*m+beta*n=0 in H1(P;Z), with nonzero integer coefficients
  alpha,beta and the specified boundary orientations;
* filling P at n gives a non-solid irreducible boundary-incompressible
  graph manifold P(n) with b1=1;
* m on A and n on C lie outside their respective L-space INTERIORS;
* Y=A union P is a non-solid boundary-incompressible graph manifold
  with b1=1, and all the gluings below are between incompressible tori.

After filling n, alpha*[m]=0 in H1(P(n)), so the surviving m is torsion.
It is therefore the rational longitude. Filling it cannot produce an L-space,
since H1 remains infinite. In particular m is outside L^o(P(n)).
Rasmussen Corollary6.1 applied to A and P(n), with the actual matched
slope m, gives

    Y(n)=A union P(n) is NOT an L-space.

Thus n is outside L(Y), and hence outside L^o(Y). Apply the SAME
gluing theorem to Y and C, using the matched slope n:

    A union P union C is NOT an L-space.                 (1)

This is a two-step argument about actual fillings and INTERIORS.
The coefficients need not be units; the application below has alpha=beta=1.
P(n) and Y must be nonsolid; no theorem is applied to a solid torus.
In each use below their irreducibility/boundary-incompressibility and
b1 hypotheses are checked. The closed end is an irreducible orientable
graph manifold. Némethi Theorem1.5 (the published graph-manifold
equivalence) makes its group nontrivial and LEFT-ORDERABLE. It cannot
be a PF74 common-source end group, which has no nontrivial homomorphism
to any left-orderable group. No negative-definite target star or
classification of rational Seifert quotients is required by (1).

## 2. Actual S3/S3 sections at arbitrary k

Use PF80's preferred meridian/longitude and actual two-stage knot.
In this sector write e=2b,o=3t,r=3s, with s,t odd and gcd(2b,3t)=1;
k is even and gcd(r,k)=1. Put A=bt and choose 0<v<r,u with

    r*u-k*v=1, h=M^k*L^r=z^r,
    beta=M^u*L^v, w=z^(-v), mu=beta*w.

beta is an actual outer section (its intersection with h is1), and mu
is the true inner meridian section. Their base relation determines the
geometric cone loop w; it is not an arbitrary relabeling of a cone
generator. In the main three-sheet piece P, after ALL true branch
meridian-cycle fillings, take X=lift(w^3), and Y,Z the fixed- and
two-sheet beta boundary sections. Then

    X^s=h^(-v), XYZ=h^j, h central,
    n=h^(-v)*Y^r, L1=h^u*Y^(-k),
    m=h^(-2v)*Z^r, L=h^u*Z^(-k/2).                     (2)

The two pairs(n,L1),(m,L) are genuine integral boundary bases:
their exponent determinants are kv-ru=-1. These formulas retain k;
no substitution k=r*(eo-e-o)+1 is made.

The lifted meridian disk in the actual companion solid torus has
boundary the fixed and two-sheet meridian lifts, so PF73's integral
surface lemma gives m+n=0. This depends on sheet degrees1 and2,
not on k or the genus of that lifted surface. No pushforward to the
companion EXTERIOR is asserted after a meridian filling.
P is an annulus Seifert piece with cone s. Its H1 is free rank2 and
h is nonzero, since gcd(v,s)=1. In additive notation, (2) and m+n=0
give r*j*h=0, hence j=0. Consequently

    X^s=h^(-v), XYZ=1                                   (3)

is the ACTUAL relative section presentation at every k.
The fixed-sheet cable caps are retained. The other pieces are

    E1: the positive T(2b,3t) exterior;
    E': disk cones(b,2b,t), boundary(m,L).

For (b,t)!=(1,1), both companions are nonsolid with negative orbifold
Euler characteristic. PF76's actual companion analysis, independent
of the cable coefficient, shows their meridians m,n are outside their
L-space INTERIORS. The weighted cubic regular-fibre filling used there
depends only on(b,t). This is an application of that retained endpoint
theorem to unchanged companion pieces, not of its old eta0 conclusion.

If s=1, P is T2 times I. Elimination in (2)-(3) gives the actual gluing

    m=-n, L=(k/2)*n+L1.                                 (4)

It matches those missing-interior meridians, so PF76's two-piece
gluing obstruction excludes every non-trefoil COMMON r=3 case at
arbitrary k.

If s>=3, fill P at n. The filling has distance r from h, so P(n)
is disk Seifert with genuine cones(s,r); it is nonsolid, irreducible
and boundary-incompressible. The balance m+n=0 makes m its rational
longitude. More explicitly PF72's integral two-cone transfer, applied
in the remaining basis(m,L), gives

    H1(P(n))=Z*L plus Z/s, m=0, h=rL.                   (5)

The map used for that transfer is to the actual companion SOLID TORUS,
where L generates homology; it is not a map to its knot exterior.
Gluing E' gives H1((E' union P)(n))=Z/s plus Z/b, a finite group.

The actual Y=E' union P has b1=1 by Mayer--Vietoris: P's boundary
maps rationally inject into its rank2 H1, and E' has rank1. All these
pieces have incompressible torus boundaries. The fibre intersections
(r*eo-k) and half that are nonzero, since k<r*eo. Thus Y is nonsolid
and boundary-incompressible; the closed full end is irreducible graph.
Apply Section1 with A=E', C=E1. This excludes EVERY non-trefoil
S3/S3 common source for ALL s>=3 and arbitrary k.

## 3. Two explicit bounds for the trefoil common exception

Now (e,o)=(2,3), r=3s with s odd. In the negative-even remainder,
2<=k<=r-1, k even. At r=3 only k=2 is possible. The full target
homology is then infinite (PF80's ordinary satellite presentation,
or (2)-(3) with the solid-companion meridian below), contradicting
the source's finite end H1. Thus assume s>=3, r>=9.

E' is ACTUALLY solid with meridian m+L. Capping P by it produces
Y with base disk and genuine cones

    (s,J), J=r-k/2 >= (r+1)/2.

Its remaining boundary is the fixed companion torus(n,L1). Equation(3)
and the actual cap give

    X^s=h^(-v), Z^J=h^(2v-u), XYZ=1.

In Rasmussen's basis(h,-Y), the two (unnormalized) cone coefficients
and the actual n slope are

    a1=-v/s=-3v/r,
    a2=(2v-u)/J=2v/r-1/(rJ),
    x= slope(n)=v/r.                                    (6)

The sign and the section are fixed by (2), not inferred from H1.
The rational longitude is -a1-a2=x+1/(rJ), as an independent check.
The section-slope0 filling is a rational-homology Seifert space with
nonzero Euler number a1+a2, making Y a regular-fibre complement for
the purpose of Rasmussen Proposition3.9. Y itself is nonsolid.
That proposition, with the integer parts kept in a1,a2, gives

    y_minus=max_(l>0) [-1-floor(l*a1)-floor(l*a2)]/l,
    y_plus=min_(l>0) [1-ceil(l*a1)-ceil(l*a2)]/l,
    L(Y)=[[y_minus,y_plus]].                            (7)

At l=r, the two floors are -3v and2v-1. Thus y_minus>=x.
For l0=r-k>0, Bezout gives l0*v=N*r+1 for an integer N>=0.
Since r>3 and J>1, the two ceilings are exactly

    ceil(l0*a1)=-3N,
    ceil(l0*a2)=2N+1,

where the second identity uses 2J-l0=r, so its nonintegral part
is precisely1/J. Therefore

    y_plus <= N/l0 = x-1/(r*l0) < x <= y_minus.           (8)

Definition3.7's circular-interval convention now puts x OUTSIDE the
INTERIOR of L(Y), including the possibility x=y_minus. No search of
denominators or assumption that the n filling itself is not an L-space
has been made. The endpoint distinction is essential.
The ordinary trefoil E1 has n outside its L-space interior by PF76.
Rasmussen Corollary6.1 excludes their actual gluing as an L-space;
it is irreducible graph and has nontrivial LO group. PF74 rules out
its common source. This closes the ENTIRE S3/S3 COMMON sector.

## 4. Cyclic companions at arbitrary k

The two cyclic sectors have r,k even product, h image1, and ordinary
R fibre. PF65's cyclic companion replacement keeps the R mark and
adds the genuine cyclic branched-companion H1 summands. The pattern
is now computed WITHOUT any eta0 relation.

For cyclic3 write r=4t,k=3s. Replacing companions by the unknot gives
the degree4 branched cover of T(r,k), with lifted cone orders
(t,s,3s) and Euler number -4/(rk)=-1/(3st).
For cyclic4 write r=3s,k=4t; its orders are(s,3s,t), with the same
Euler number. The actual base degree is4, generic fibre degree1 and
base genus0. Actual gcd(r,k)=1 gives lcm=3st in either case. PF72's
exact fibre-order formula shows h has order1. The full pattern group
H1 is CYCLIC, of order s in either case, by the cone presentation.
R is that ordinary fibre, so its class is ZERO. Replacement preserves
this zero mark in every companion summand. Hence PF67 closes EVERY
cyclic COMMON source, independently of k.

In cyclic3, if 3|pq the added companion group is (Z/other exponent)^2,
contradicting the source's zero-R cyclic-H1 requirement. If 3 does not
divide pq, the three-sheet companion has generic fibre degree3,
base degree1, and is disk(p,q); the fixed-sheet companion is the same
torus exterior. Both are nonsolid. The drilled main cable has THREE
base boundaries and cone r/4. The R boundary is in the MIDDLE with
two genuine companion JSJ edges, since both fibre intersections have
magnitude r*pq-k>0. PF77's source endpoint has at most one neighbour.
Full peripheral recognition therefore excludes all cyclic3 COMPLEMENTARY
sources as well.

For cyclic4, 4|e is excluded by its added (Z/o)^3 homology and zero R.
If e=2 mod4, the companion has fibre degree2 and cones(e/2,o,o).
If both first exponents are odd, its fibre degree is4 and its disk
cones are(p,q). These actual fibre-orbit descriptions are retained.
For s=r/3>=3 the root annulus cones(s,3s) are nonconsecutive, excluding
the complementary source. At s=1 the root has one cone3 and R-distance1.
The rho=2 source collar has only distances0 or3 and rho=1 has cone2.
Thus rho=4,q_src=3 is the only higher-order root option and the D side
must be a SINGLE disk matching the companion cones.

If e=2 mod4 with e>2 the companion already has three genuine cones,
violating PF79's terminal budget. For e=2 its cones(o,o) force smooth
D with sigma=1,p_D=o, so det Q=13*(2o)-(3o)^2<0 for o>=3.
For both first exponents odd, let P=max(p,q), K=min(p,q), P>=5,K>=3.
The target end H1 is0: the pattern order s is1 and the cyclic4
branched-companion H1 vanishes. A matching D disk has p_D=P,
kappa=K,sigma=P+1-K and n_D=P^2+(1-K)P+K. The full source
determinant would have to satisfy

    1=13n_D-9P^2=4P^2+13(1-K)P+13K,
    (13K-4P-17)*(P-1)=16.                              (9)

Since P is odd and P>=5, P-1 must be4,8 or16. The corresponding
P=5,9,17 give K=41/13,55/13,86/13, none integral. This finite
DIVISOR consequence of an unrestricted identity closes every all-odd
cyclic4 r=3 case, with no target-degree cutoff. All cyclic sectors
are now excluded in BOTH source pencils.

## 5. The other noncyclic sectors, without the eta0 substitution

Use the complete general census in PF80. Its planar even-product
sectors and ACTUAL cone lists are unchanged at arbitrary k; only the
peripheral shears change. Their earlier generator budgets remain
valid because they kill all boundary/fibre classes at any slope.
The following transfer is verified using actual bases, not residue
coincidence:

* S4/S4: PF70's repeated-cone quotient forces e=2. Write o=3t,
  r=3s with s,t odd. The actual basis(H,m=M^4) has
  h=rH+j*m, j=(k-r*eo)/4, an integer coprime to s.
  PF72's integral transfer still gives m=0,h=rH and cable torsionZ/s.
  Thus PF71's full H1=Z/(8t) plus Z/s, order4 R mark and self-linking1/2
  hold at ANY such k. Common-source self-linking is strictly below1/2.
  Complementary sources fail the terminal budget because the companion
  disk(2,2,t,3t) always has at least THREE genuine cones.
* A4/S3: the actual basis(m=M^3,L) gives h=(k/3)m+rL, and cones
  (r/2,r). Coprimality gives the unit coefficient in PF72, so m=0,
  h=rL, H1=Z/(eo/3) plus Z/(r/2), R=0 remain valid. All common
  sources fail. Cyclic H1 forces r=2 mod4. For nonsolid companions,
  the complete PF79 terminal/nonconsecutive-root proof is unchanged:
  it uses only these orders, R-distance1 and the full target H1.
  Its r=2,e=2,o=3t residual source determinants are2t(3t+14) and
  12t(3t+1), whereas the target order is2t, independently of k.
* S4/C2: here r is odd,k even, e=4b,o=3t. The basis(m=M^2,L)
  gives h=(k/2)m+rL, cones(r,r) and the same unit transfer. Full
  H1=Z/t plus Z/r and R=0 close all common sources. Every nonsolid
  companion fails the complementary consecutive-root criterion,
  since both root orders are the equal genuine number r.

The nonsolid companion fibres differ from the cable fibre, with
positive determinant proportional to r*eo-k. They remain genuine JSJ
edges, with negative-Euler pieces. All source terminal, marked-slope,
Klein and Waldhausen hypotheses are exactly those audited in PF79/80.

There are TWO small solid-companion complementary exceptions to check.
The exact whole-Seifert source classification from PF79 is stronger
than merely consecutive cones: its three-cone alternative is

    disk(q_src,q_src+1,2), q_src>=2.                      (10)

For A4 at(e,o)=(2,3), k is an odd multiple of3 and r=2 mod4 with
3 not dividing r. Negative eta gives k<=r-1, so r>=10. The actual
solid meridian is m+2L OR7m+2L, both retained as in PF79. Its cap
orders are J=r-2k/3 OR7r-2k/3, both at least4. The target single
disk has three cones(r/2,r,J), all exceeding2. It cannot match (10)
or the two-cone source alternative.

For S4/C2 at(e,o)=(4,3), negative eta gives k<=5r-1. The actual
meridian is3m+L OR9m+L, giving J=3r-k/2 OR9r-k/2, with J>=2.
The target disk has cones(r,r,J), with r odd>=3. The only repeated
order in (10) occurs at q_src=2 and is the EVEN number2. It cannot
match the target's repeated ODD r. Both cap signs are thus excluded.

## 6. Complementary S3/S3 and final exhaustion

At a non-trefoil first pair both companions are nonsolid; after drilling
R the central base has THREE boundaries and cone s=r/3. Both edge
fibre intersections are positive, r*eo-k and half that. The R boundary
has two canonical JSJ neighbours, contradicting PF77's endpoint.
This works for arbitrary k, including s=1 (a pair-of-pants product,
not a collar).

At the trefoil pair, capping the actual solid companion gives root
annulus cones(s,J), J=r-k/2. For s>=3 and negative eta, J>=
(3s+1)/2, so J-s>1. Source root cones must be consecutive, contradiction.
For s=1 only k=2 remains and the full target H1 is infinite as in
Section3. Thus every S3/S3 complementary source also closes.

PF80's full9216-cell certificate exhausts the general two-stage
peripheral representations:448 graph-cycle failures,128 positive-genus
tree failures,256 odd-product classes already closed, and1920 even
planar tree types in the SIX sectors above. Sections1-6 close ALL
1920 remaining types for every actual integer lift with eta<0.
PF80 closed eta0, PF59 eta>0, and PF80 the odd-product sector.
PF50's actual source pencil dichotomy has been retained throughout.

Therefore the ENTIRE general genuine TWO-STAGE QUARTIC target class
is excluded in BOTH pencils, at all p,q,r,k, all source data, and
with arbitrary lower-weight key terms. Polynomial curves and infinity
representations are not declared nonexistent; the contradiction is
with their FULL marked quartic finite normalization/source packet.

## Evidence, limits and next gate

The [checker](verify_meridian_balance_two_stage.py),
[receipt](meridian_balance_two_stage_receipt.json) and
[manifest](sources/meridian-balance-two-stage/manifest.json) seal the
complete generic census, arbitrary-k section/mark identities, meridian
balance, explicit trefoil denominator bounds, actual generic satellite
H1/R controls and the finite divisor16 certificate. Universal claims
come from the written proofs, not interpolation in those controls.
Complete current Rasmussen PDF8/9/18/41 and Némethi PDF2/3 were read
visually. Interval/gluing/orderability and retained JSJ/recognition
theorems are imported, not independently reproved. The old general
Seifert quotient remains a historical exploratory lead; this proof
does not need to classify its rational cases or guess its section.

The next target gate is THREE OR MORE genuine characteristic stages:
propagate actual missing-interior slopes through longer pieces using
the full lifted-surface balance coefficients, and compare their
boundary-rooted JSJ data to the two-line source. Do not assume a longer
splice is two-stage or that its fibre image remains identity.
Degenerate stages, general higher-degree source packets and the full
quartic/higher-degree goal remain OPEN. No finite Keller map, formal
verification, independent full reproof, external review or novelty/
priority claim is asserted; the full goal remains active.
