# Full peripheral linking and an all-length longitude-joining recursion

PF2RULED101, 16 November 2026. The user's knot-theory route extends beyond
the longitude-identity and integrally null-core hypotheses of
[PF92](RECURSIVE_LINKING_AND_ENDPOINT_TAIL_CLOSURE.md). This note proves
the full-peripheral surface balance, rational core-linking reconstruction
with its exact integral torsion defect, and a complete all-length recursion
for a family with genuinely nontrivial longitude monodromy. No Keller map
or affine realization of that family is asserted. The actual finite-algebra
gap in [the finish-line note](LAURENT_FINITE_ALGEBRA_FINISH_LINE.md) remains
open; published Domrina2000 quartic retirement remains in force.

The all-parameter part uses a small **unimodular certificate over Z[t]**,
not a finite parameter sweep. The final recursion applies at arbitrary
finite cover degree subject to the stated peripheral orbit types; the
eight- and twenty-four-sheet trefoil covers are concrete seeds.

## 1. Actual boundary tori, filling normals and the shear balance

Let E be any connected finite cover of an oriented knot exterior in S3.
Use the positively linking meridian mu and Seifert longitude lambda.
Their commuting permutations act jointly on the sheets. On boundary
orbit i, let the stabilizer lattice have normalized basis

    (e_i,0), (a_i,f_i),   e_i,f_i>=1, 0<=a_i<e_i.       (1)

The orbit has size e_i f_i; it need not be one meridian cycle.
Longitude can join f_i such cycles. The actual primitive oriented basis
of this covering torus is

    m_i=lift(mu^e_i), L_i=lift(mu^a_i lambda^f_i).

Fill all m_i to obtain W; orient the cores K_i along L_i. They map to
the base knot with longitudinal degree f_i. These are the actual normal
fillings, not a replacement by lifted lambda curves. The integral torus
lattice calculation is the retained [PF98 helper](peripheral_boundary_lattices.py).

The preimage of a Seifert surface gives the **integral** relation

    sum_i (e_i L_i-a_i m_i)=0 in H1(E;Z).              (2)

Indeed put g_i=gcd(e_i,a_i). On orbit i, lambda has g_i cycles of
length f_i e_i/g_i. Each boundary lift represents

    -(a_i/g_i)m_i+(e_i/g_i)L_i.

Adding those g_i curves proves (2). This includes a_i=0, when g_i=e_i.
It also explains why f_i does not multiply the coefficients in (2).
No null-homology assumption on the filling cores enters this proof.

## 2. Rational matrix, core orders and the integral information it loses

Assume every K_i is rationally null-homologous in W. W itself may have
positive first Betti number. The drilling exact sequence over Q and
intersection with rational bounding surfaces make the m_i independent
and give unique rational coefficients

    L_i=sum_j F_ij m_j in H1(E;Q).                    (3)

For i!=j, F_ij is lk_W(K_i,K_j). Diagonal entries keep the chosen lifted
framing L_i. A rational bounding surface for a multiple of K_i, punctured
at every K_j, proves this identification. Rational linking is symmetric;
thus (2) implies the complete shear balance

    F=F^T,  F*e=a.                                   (4)

The conventions agree with the rational linking definition and chain
method in [Cahn et al., Dihedral linking invariants](https://arxiv.org/html/2112.14790v2).
Their odd-dihedral algorithm supplies the comparison on the overlap;
it is not being applied to arbitrary permutation covers. In particular,
the new eight-sheet cover below is not called dihedral.

There is essential integral information beyond F. Let q_i be the exact
order of [K_i] in H1(W;Z), taking q_i=1 for the zero class, and let D_i
be the least common denominator of row i of F. An integral bounding
surface for q_i K_i proves

    D_i divides q_i,
    q_i L_i=sum_j(q_i F_ij)m_j in H1(E;Z).             (5)

The minimal-denominator expression

    T_i=D_i L_i-sum_j(D_i F_ij)m_j

is a torsion class of **exact** order q_i/D_i in H1(E;Z). Its image in
H1(W) is D_i[K_i], of that order; (5) gives the reverse bound, proving
equality. Also the denominator of F_ij divides gcd(q_i,q_j), by using
bounding multiples of both components. These statements are valid
with other free or torsion homology summands retained.

When q_i=1 for every i, all relations (3) hold integrally and F is
integral, recovering PF92. Without that hypothesis one must not
clear denominators and declare the resulting row zero. Nor does F
alone determine H1(W), the core orders or the complete torsion linking
form of W. The full integer presentation and peripheral rows are
retained to compute arbitrary further fillings exactly.

For one boundary, (4) reads F_11=a/e. A nonzero normalized shear
therefore forbids an integrally null core. For two integrally null
cores, put e_1=g u,e_2=g v with gcd(u,v)=1. Necessarily g divides both
a_i. Write a_i=g b_i. The off-diagonal entry ell satisfies

    v ell=b_1 mod u,  u ell=b_2 mod v.                (6)

These two congruences have a unique residue ell_0 modulo u v, and
every integral solution is

    ell=ell_0+u v z,
    F_11=(a_1-e_2 ell)/e_1,
    F_22=(a_2-e_1 ell)/e_2, z in Z.                   (7)

PF92's divisible linking grid is precisely a=0. This is a necessary
arithmetic filter, not a realization of any matrix by a knot cover.

## 3. Full peripheral computation at arbitrary cable length

The [new helper](full_peripheral_linking.py) lifts presentation paths
as pairs (destination permutation, integral oriented edge flow). The
product rule is

    dest_(AB)(s)=dest_B(dest_A(s)),
    flow_(AB)(s)=flow_A(s)+flow_B(dest_A(s)).           (8)

Inverse paths reverse their signed flows. Binary powers give signed
integer exponents without expanding words. For an initial T(p,q) use
<A,B | A^p=B^q> and its actual meridian/longitude. Each cable (r,k)
adjoins z with z^r=mu^k lambda^r. Choose r u-k v=1; then

    mu_new=mu^u lambda^v z^-v,
    lambda_new=mu_new^(-r k) z^r.                    (9)

These are the root-amalgam and both peripheral words from
[Clay--Watson, Section2](https://arxiv.org/html/1103.2358v1).
Changing (u,v) to (u+k h,v+r h) leaves the displayed meridian unchanged:
mu^k lambda^r=z^r commutes with mu,lambda and z inside the cable space.
Thus the signed Bezout choice used in the certificate agrees with the
positive choice in that source. No surgery or orderability theorem is
imported from the paper.

Every lifted relator bounds a lifted 2-cell. Collapse a maximal tree of
the covering presentation graph, then retain those integral relation
rows. This is a cellular presentation of H1(E), independent of whether
the presentation complex is aspherical. Evaluate the closed words
mu^e_i and mu^a_i lambda^f_i at the correct joint-orbit seed to obtain
the actual boundary rows. The helper validates all group relations,
commutation and connectedness before using them. No longitude-identity
shortcut remains.

For rational linking, quotient the relation space over Q and solve (3).
An integer HNF basis separately gives the exact order of each longitude
modulo the relations and normal fillings. Membership of a rational row
in the span alone is insufficient; (5) and the exact torsion defects are
verified in the integer row lattice. Different maximal trees are exact
changes of generators and give the same marked answers. Adding any
chosen primitive filling slope b_i m_i+c_i L_i to these actual rows
gives its full integral H1, including summands invisible in F.

## 4. A complete nontrivial-longitude cable theorem

Suppose every peripheral orbit of a connected companion cover has

    (e_i,f_i,a_i)=(1,2,0) or (3,2,0).                 (10)

Thus mu^3=1, lambda^2=1, and lambda interchanges pairs of meridian
cycles. Take **any** t>=1 and the genuine ordinary cable

    r=36t+1, k=36t-5.                                (11)

Both parameters exceed one and gcd(r,k)=gcd(r,6)=1. Extend the actual
permutation representation by rho(z)=rho(mu)rho(lambda). All peripheral
permutations commute; since r=k=1 modulo6 the root relation holds.
The new meridian and longitude permutations remain mu and lambda.
Connectedness persists since the old transitive image is retained.

Each covered cable-space component has presentation

    <m,l,z | [m,l]=1, z^r=m^k l^r>,

with cyclic sheet action of size 2e, e=1 or3. In cyclic coordinates the
actions of m,l are shifts (0,1) at e=1, or (2,3) at e=3. The action of z
is their sum. Each component has one old and one new boundary torus.
The [exact Z[t] certificate](longitude_joining_transfer_certificate.py)
proves for **every** integer t>=1

    H1(pattern;Z)=Z*m_new direct-sum Z*L_old,
    m_old=r*m_new, L_new=r*L_old.                     (12)

Here is the certificate's universal proof mechanism. The lifted
commutator and root relators give two matrices over Z[t], with four
marked port rows. Powers are computed using a^(b+6q)=a^b(a^6)^q;
a^6 has identity destination, so its flow multiplies by the integer
polynomial q, even for negative exponents. The fixed tree depends
only on the permutations. A finite list of elementary integral row
and column operations, and cancellation of unit pivots, reduces each
matrix to zero relations on two generators. All operations are
unimodular over Z[t], so specialize to unimodular operations at every
integer t; no exceptional integer parameter is discarded.

For e=1 the final port rows (m_old,L_old,m_new,L_new) are

    (36t+1,0), (5-36t,1), (1,0),
    (-1296t^2+144t+5,36t+1).

For e=3 they are

    (-36t-1,36t+1), (36t-1,2-36t), (-1,1),
    (1296t^2-1,-1296t^2+36t+2).

In both cases m_new,L_old have determinant +/-1, and the two identities
in (12) are polynomial identities. The residual e=1 row
(5-36t,36t+1,0) acquires a unit by column1+=column0 followed by
column0+=(6t-1)column1. The e=3 row (-2,36t+1,2) acquires a unit by
column1+=18t column0. These are the only parameter-dependent unit
creation steps; the full fixed operation lists are retained. This
establishes (12) for unbounded parameters, not from sample evaluation.

Mayer--Vietoris now glues these two-port pieces to the ACTUAL old cover.
Each piece attaches to only one old boundary, so there is no extra
cycle in the gluing graph. Eliminating its L_old generator leaves

    H1(E_new)=(H1(E_old) plus direct-sum_i Z*m_i,new)
                    /<m_i,old-r*m_i,new>.             (13)

After normal filling, H1(W_new) is naturally H1(W_old), and

    [K_i,new]=r*[K_i,old],
    q_i,new=q_i,old/gcd(q_i,old,r) for finite q_i,old.  (14)

An infinite-order core remains of infinite order. Rational nullity is
therefore preserved. From (12)-(13), whenever it holds,

    F_new=r^2 F_old.                                 (15)

This proof keeps arbitrary other homology summands and even nonzero
torsion core classes. It does not identify the full torsion linking
forms of W_old and W_new. Iteration is legitimate because the entire
peripheral action, including longitude joining, is unchanged. Hence
for any finite list t_1,...,t_s>=1, with R=product_j(36t_j+1),

    F_s=R^2 F_0, H1(W_s)=H1(W_0),
    [K_i,s]=R*[K_i,0], q_i,s=q_i,0/gcd(q_i,0,R).       (16)

No upper bound on length or parameters is used. The theorem requires
the actual orbit types (10) and the specified root/coloring choice.
It is not asserted for all cable representations or for every knot
at infinity of a polynomial curve.

## 5. Explicit seeds and essential torsion controls

Use T(2,3), with group generators A^2=B^3. In SL(2,F3), take

    A=((0,2),(1,0)), B=((0,2),(1,1)).

Their actions on the eight nonzero vectors give a connected eight-sheet
cover. Meridian has two 3-cycles and two fixed sheets; longitude -I joins
them in pairs. The actual boundary types are (1,2,0),(3,2,0), not four
separate tori. Normal filling gives H1(W)=0, so all cores are integrally
null. In the index1/index3 ordering,

    F_0=[[-9,3],[3,-1]].                              (17)

Thus the entire family (11) has mutual linking 3 R^2 and the matrix
R^2 F_0, at every finite length. On the regular twenty-four-sheet action
there are four boundary tori of type (3,2,0), normal filling H1=0,
and F_0 has diagonal -3 and all off-diagonals1. Its whole matrix also
scales by R^2. These are colored-cover invariants, not a single value
for the underlying uncolored iterated knot.

For a distinct essential control, map A to (01) and B to (123) on four
sheets. This connected trefoil cover has one boundary (e,f,a)=(4,1,2).
Its exact data are

    H1(E)=Z plus Z/2, H1(W)=Z/8,
    F=[1/2], core order q=4, row denominator D=2.

The rational expression 2L-m is a NONZERO class of order2 in H1(E),
while 4L-2m=0 integrally. This simultaneously shows that F*e=0 is false
with shear, rationally null is weaker than integrally null, and clearing
the least denominator can lose torsion. Even the same scalar F=1/2
occurs on an abstract solid-torus peripheral model with H1(W)=Z/2 and
core order2; the matrix alone is not the full marked invariant.

The checker compares the longitude-identity overlap with the sealed PF92
helper and its unmodified-author-code dihedral comparisons at degrees
3,5,7, comparing only off-diagonal entries. It also independently lifts
short literal words, checks both maximal trees, all integer surface
balances and core orders, and checks the symbolic certificates above.
Large parameters use binary flows; there is no expanded cable diagram
or deeper coefficient campaign.

## 6. Scope for the active Wright problem

The algebraic work this turn did not prove that actual finite Keller
normalization forces a polynomial terminal key. The tame codifferent
can make low trace moments regular; Newton identities force integrality
only after enough actual moments are regular. Affineness and exactness
of the retained Laurent controls do not supply that requirement. No
automatic trace regularity or source-action descent is inferred.

The knot fallback now has a complete all-degree full-peripheral
reconstruction and an all-length family beyond the longitude-one
restriction. Before using it as a Wright exclusion, match the ACTUAL
finite target cover and marked source to its orbit types, fillings,
rational matrix and integral core classes. Neither an abstract group
nor a matching uncolored knot supplies that comparison. The complete
global Laurent-key/finite-algebra implication remains open, and no
new higher-degree Keller exclusion is claimed from these examples.
Quartic realizations remain retired by published HC4MYKR0. Closure
precedes external-review expense; no outreach is requested or performed.

Run `.venv/bin/python research/jacobian-followup/affine-complements/verify_full_peripheral_linking.py`
under the retained 512MiB limit. The [receipt](full_peripheral_linking_receipt.json)
pins the proof, helpers, symbolic certificates, selected primary sources
and predecessor receipts. Universal claims rest on (2)-(16), including
the Z[t] unit certificate. Finite controls verify signs, implementation
and comparison; they do not establish arbitrary realizability, a full
torsion linking pairing, formal verification or external review.
