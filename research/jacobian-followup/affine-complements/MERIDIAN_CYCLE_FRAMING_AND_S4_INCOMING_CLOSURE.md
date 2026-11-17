# Meridian-cycle framing and the S4 incoming closure

15 November 2026. PF2RULED86 excludes incoming states4/5 of PF82
in BOTH source pencils, at arbitrary companion depth and arbitrary
later tail length, without a prefix eta assumption. These are the
noncyclic S4 companion states with a transitive 4-cycle meridian and
longitude1 or its meridian square.

The new all-degree lemma reads the filling-core self-linking directly
from a lifted preferred-longitude surface. It needs no classification
of the companion and no detailed coefficient calculation. For the
quartic states the value is0 or1/2, violating the strict common-source
fold bound. Complementary sources fail a separate persistent JSJ
branching/terminal-cone obstruction.

Together with PF85, all cyclic incoming states6/9 and S4 incoming
states4/5 are closed. Incoming states2(S3),7(A4), and initial state3
trajectories outside closed intrinsic seed ranges remain open.
No full quartic or higher-degree exclusion is asserted.

## 1. An all-degree longitude-surface framing lemma

Let E be a connected n-sheet UNBRANCHED cover of any oriented knot
exterior. The positive meridian permutation M is an n-cycle.
Its commuting preferred longitude permutation is necessarily M^ell,
0<=ell<n, since the centralizer of an n-cycle in S_n is C_n.
There is ONE boundary torus. Its ACTUAL primitive basis is

    m=lift(M^n), Q=lift(L*M^(-ell)).                    (1)

Lift an oriented Seifert surface of the knot, as in PF73. If
l=n/gcd(n,ell) is the order of the longitude permutation, the surface
has n/l boundary lifts. Each is parallel to

    lift(L^l)=l*Q+(ell*l/n)*m.

Their integral sum proves

    n[Q]+ell[m]=0 in H1(E;Z).                          (2)

No regular-cover, Seifert-base, degree, genus or finite-homology
hypothesis is used in this balance.

Fill m by a solid torus to obtain W. Its core K is represented by Q.
Assume H1(W) is finite. The lifted surface has boundary nQ+ell m;
cap its m-boundaries by meridian disks in the filling solid torus
and join its n Q-boundaries to n parallel core copies. A second
parallel core can be kept inside the filling solid torus, disjoint
from the surface in E. Its intersections with the caps give

    lambda_W(K,K)=+/-ell/n mod Z.                      (3)

The overall sign depends on the boundary convention. Integer framing
changes do not affect(3). This is a chain-level calculation; it is
not inferred from the order of K, and does not assert that Q vanishes.
In particular the two orientation-independent cases are

    ell=0: self-linking0;
    n even, ell=n/2: self-linking1/2.                  (4)

If W has positive b1, any finite-end application must instead retain
that obstruction; no torsion self-linking formula is invoked then.

PF73's integral replacement theorem applies whenever an actual filled
pattern has H1=Z*Q plus T, m=0, with boundary maps containing NO T
coordinate. If its mark is c*Q, attaching it preserves W's linking
summand and gives self-linking +/-c^2*ell/n. Thus c odd preserves
both special obstructions in(4). The cover/framing statement is valid
at ALL finite degrees and for arbitrary companions, including those
whose exteriors are not graph manifolds. The pattern hypotheses are
separate and must be proved in every application.

## 2. Quartic longitude images and actual final cable

Take incoming companion image S4, M a4-cycle. Its longitude image
is1 or M^2, as PF82's complete states4/5 say. There is also a direct
group argument: a Seifert surface has all its loops in the knot
commutator subgroup, so its boundary longitude lies in the SECOND
commutator subgroup. S4''=V4, and V4 intersect <M>={1,M^2}.
Therefore ell=0 or2 independently of a length cutoff.

For a genuine final cable(r,k), use PF80's actual words ru-kv=1,

    F=M^k L^r=z^r, beta=M^u L^v,
    mu=beta*z^(-v), lambda=mu^(-rk)*F.

The final simple periphery requires mu=tau a transposition and
lambda=1. F belongs to <M> and also equals tau^(rk). A4-cycle's
cyclic group contains no single transposition, so

    F=1, rk even, k+ell*r=0 mod4.

If r were even then k would be odd, making k+ell*r odd. Hence
r is odd, k even. From z^r=1 in S4, z is identity or a3-cycle;
identity cannot produce a transposition from beta, a generator of
<M>. Thus z is a3-cycle and

    r=3s, s odd, cable image S4.                       (5)

The actual filled four-sheet cable piece P has base disk with cones
(s,3s), orders1 omitted. ALL branch and fixed-sheet meridian-cycle
fills are retained. The ordinary R is its regular fibre F.
In the basis(1),

    F=A*m+r*Q, A=(k+ell*r)/4, gcd(A,s)=1.

The map P to the actual companion SOLID TORUS sends m to0,Q to
a longitude generator and F to r times it. PF72's integral unit
transfer, including s=1, gives

    H1(P)=Z*Q plus Z/s, m=0, [R]=r*Q.                 (6)

This uses no map to the companion exterior AFTER a meridian filling.
PF73 now gives

    H1(full end)=H1(E(m)) plus Z/s,
    [R]=(r[K],0), lambda(R,R)=+/-r^2*ell/4.            (7)

If H1(E(m)) is infinite the full end cannot be a PF61 source end.
Otherwise(7) is0 for ell0 and1/2 for ell2. Both violate PF71's
strict 0<lambda_source(R,R)<1/2 for any actual COMMON source with a
sole fold and the stated smooth/whole-line/common-pencil hypotheses.
That source theorem has no field-degree, omitted-line-count or
characteristic-depth bound. No graph-manifold hypothesis on E is
needed in this common-source application.

Every later passing state3 cable has odd winding by PF82; its actual
mark transfer multiplies self-linking by an odd square and may add
unmarked torsion. Thus the same obstruction persists through ALL tails.

## 3. Complementary sources: the unbranched companion is retained

For the actual iterated torus-knot companion J, E is its full
connected four-sheet unbranched exterior cover. Since M is transitive
on all sheets, there is only ONE boundary. The final cable has no
other companion boundary to cap: its fixed meridian fills take place
inside P, not inside E. Consequently EVERY lifted companion JSJ torus
stays within the ramification exterior.

Unbranched covers of the genuine knot pieces are irreducible,
boundary-incompressible, and have negative-Euler orientable bases.
Different original fibre lines remain different in a finite-index
torus cover. For adjacent ordinary knot cable levels the determinant
is proportional to r*B-k, with integral B=pq or r_old*k_old.
It cannot vanish by gcd(r,k)=1 and r>=2. Thus these are actual
canonical JSJ pieces, with fibre uniqueness and full peripheral
recognition exactly as in PF77/79. No homology-only recognition,
Klein fibre assumption or geometric deletion of cyclic loops is used.

PF82/84's complete all-length state graph gives a finite dichotomy
for a path reaching state4 or5 before final arrival at state3:

1. Its first S4 image was on the INITIAL torus knot in state4/5.
2. It first became S4 at a cable from state1(C2, two fixed sheets)
   or state6(C3, one fixed sheet), into state4/5.

Afterward it has only cyclic4 self loops until the final state3
arrival. There are no other ways to enter4/5 in the complete table.
Paths already in state3 cannot return to them.

## 4. Both branches of the dichotomy violate the source path

In case1 let a,b be the torus generator permutations. Their common
central fibre is identity since the image is S4. Coprime torus
exponents force their orders coprime. A planar connected base with
one boundary requires c(a)+c(b)=5. The ONLY full-S4 possibility
is a transposition and a3-cycle: the4-cycle/3-cycle possibility has
positive genus and is already excluded by PF69.

Write the actual even exponent e>=2 and odd exponent o=3t.
The terminal cover has disk cone orders

    e/2,e,e,t,3t, with orders1 omitted.

The two e cones and the3t cone are always genuine, so there are
AT LEAST THREE. Each later cyclic annulus has genuine cones and
two incompressible boundaries. The terminal disk remains terminal
under all attachments. PF79's complementary source has exactly TWO
genuine terminal cones opposite R whenever any canonical edge exists.
The final drilled cable is annulus with cone3 or cones(s,3s), so a
genuine edge does exist. Contradiction.

In case2 the cable image at the first S4 level is transitive. Its
boundary toward the next level is ONE torus, because its new meridian
is a4-cycle. Its older companion boundary has THREE torus lifts
for state1 (orbit sizes2,1,1), or TWO for state6 (sizes3,1).
They attach genuine unbranched covers of a NONTRIVIAL companion knot.
All are nonsolid and boundary-incompressible.

Hence the first S4 cable vertex has degree4 or3 after the next
level/final R piece is attached. Its base has at least three boundary
circles and negative Euler characteristic. It is an actual branching
canonical JSJ vertex; no caps at the eventual final arrival touch
its old companion tori. Later cyclic4 self loops and state3 tails
cannot remove this vertex, by the actual nonzero fibre determinants.
But PF79's complementary source exterior is a PATH. Contradiction.

This closes ALL incoming states4/5 in the complementary pencil,
at any companion depth and any later tail length. If an earlier
graph cycle or base genus occurs, PF62/69 had already excluded it.

## 5. Evidence and the remaining frontier

The [checker](verify_meridian_cycle_s4.py),
[receipt](meridian_cycle_s4_receipt.json) and
[manifest](sources/meridian-cycle-s4/manifest.json) retain the integral
surface relation on full unbranched CW matrices, true filled-pattern
and marked-group controls, the exact4/5 history dichotomy, planar
initial torus types and actual branching examples in both longitude
states. The written surface/cap proof computes self-linking; CW group
orders alone do not compute it. PF71's independent known half-linking
control is retained. General statements are not extrapolated from
the finite examples.

The relevant primary literature on rational linking and rational
longitudes was consulted; no additional classification theorem is
imported. PF73's oriented surface and torsion-replacement proofs,
PF72's integral transfer and the complete retained Hatcher/Neumann-
Raymond/relative-Waldhausen audits supply the stated deep inputs.

PF85 closes cyclic incoming6/9; this note closes S4 incoming4/5.
PRIORITY: incoming2(S3) and7(A4), with their ACTUAL fixed-sheet caps,
plus initial full-S4 state3 torus trajectories outside already closed
intrinsic seed ranges. The sharp two-noncyclic-level bound remains.
No full longer-stage/quartic/higher-degree exclusion, arbitrary
higher-degree smoothness, finite Keller map, formal verification,
external review or novelty/priority claim is asserted. Full goal OPEN.

