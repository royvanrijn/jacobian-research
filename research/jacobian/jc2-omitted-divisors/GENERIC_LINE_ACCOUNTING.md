# Generic-line genus, punctures and the completed first-image catalogue

14 November 2026. PF2OMIT14 applies under the four oriented degree-six
gap-21 hypotheses of PF2OMIT13. It adds global degree/genus obligations
for the unknown images and reconciles the complete first-image collision
catalogues. It neither constructs those images nor supplies the source map.

## 1. Every other image has degree a positive multiple of five

Use the normalized target coordinates in which the distinguished image
C has normalization degrees(6,10). Chau's
[Theorem3](https://arxiv.org/html/0711.3894v1), with a
[pinned source](sources/generic_line_manifest.json), says that each
nonproper component of a plane Keller map has a polynomial parametrization
whose coordinate-degree ratio equals the source map's degree ratio.
Such a parametrization factors through the normalization A1 by a
nonconstant polynomial self-map. Composition multiplies both degrees
by the same positive integer. The normalization degrees of C therefore
force the ratio3:5, and every other normalized nonproper image has
degrees(3k,5k), k>=1. Its projective degree is5k by birationality.

Write deg B=5k, and in row A write deg H=5l. If H=C then l=2;
if H=B then l=k. These are image degrees; f_R=f_S=f_D=1 remains a
separate residue-field statement. The identities below hold for all
positive k,l, with additional source-degree restrictions as follows.

Write the source coordinate degrees as(3M,5M), forced by the same
ratio. Choose a generic admissible source coordinate T=x+c*y as in
PF2OMIT3. With beta=y, the equation P(T-c*beta,beta)-U is monic of
beta-degree3M after a constant rescaling. Primitivity of T over the
target field says that Q generates the full degree3M extension over
C(U,T). Its norm polynomial is therefore the primitive inverse R_6,
monic of V-degree3M. In particular the V^(3M) coefficient is a
nonzero constant, so the leading T^6 coefficient has V-degree at most
3M-1.

On the other hand the complete pole formula of PF2OMIT3 gives

    a6=constant*H_C^p_R*H_B^p_S*H_H^(t*p_D),
    deg_V a6=3W, deg a6=5W,
    W=2p_R+k*p_S+t*l*p_D<=M-1.                     (G0)

Here the H polynomials are reduced image equations, their V-degrees
are6,3k,3l, and repeated images have their exponents added. All pole
orders are positive. The coefficient2 in W is C's degree parameter,
not its generic ramification index; S's index does not enter this pole
charge. This is why unramified omission must be retained here too.

For F2, M=25. The full inverse leading coefficient has degree at most120,
and the necessary image-degree options are

| Grouping | Degree-parameter constraint | Maximum deg B |
|---|---|---:|
| A, three distinct images | k+l<=22 | 105 |
| A, H=C | k<=20 | 100 |
| A, H=B | k<=11 | 55 |
| B,N,C | k<=22 | 110 |

The general upper bound110 improves the individual120 bound of PF2ATC1
by charging the already known C. The shared-image improvements come
from charging the extra omitted prime even though no new image support
is added. Actual pole orders can only strengthen(G0). These caps require
the source degrees(75,125); they are not imposed on the general hypothesis.

## 2. A generic target line: three distinct counts

Let W be the finite normalization of target P2 in the source field.
For a sufficiently general target line L, the curve Gamma=W times_P2 L
is smooth, connected and projective, of genus g, and maps to L=P1 with
degree6. Generic Bertini avoids the finite singular locus of W and
gives smoothness. Connectedness follows from the ample pullback of a
line: two disjoint nonzero smooth components would have positive
self-intersections and zero mutual intersection on a resolution,
contradicting the Hodge index theorem. Thus Gamma is irreducible.

Over generic target infinity let the full DVR data be(e_i,f_i), with

    sum e_i*f_i=6, n=sum f_i, 1<=n<=6.               (G1)

The integer n is the number of geometric infinity points of Gamma,
not automatically the number of global infinity primes. Their total
different is6-n. No individual f_i is set to one here without its
separate model-specific input.

Put e=e_S, so e=2 in A,B,N and e=3 in C. Put t=1 in A and0 elsewhere.
Choose L to avoid all finite singular values and pairwise intersections
of distinct image supports and to cut each image transversely. The complete
affine omitted ledger on Gamma is:

| Original omitted prime | Number of points on Gamma | Length per point | Total lost sheets | Different |
|---|---:|---:|---:|---:|
| R over C | 10 | 2 | 20 | 10 |
| S over B | 5k | e | 5ek | 5(e-1)k |
| D over H, row A only | 5l | 1 | 5l | 0 |

R,S,D are counted as primes. If H shares C or B, its points on Gamma
are still distinct from the points on the ramified prime in the
original affine normalization. Their deletion is not charged only
once because the image equations agree. Generic L avoids the isolated
cusps, so the global excess1 in B,N does not add a generic-line branch
point or an extra term to this table.

Riemann--Hurwitz, including all infinity primes, now gives

    2g-2=-12+10+5(e-1)k+(6-n),
    2g=6+5(e-1)k-n.                                 (G2)

The affine normal curve N_L has Euler characteristic

    chi(N_L)=2-2g-n=-4-5(e-1)k.                     (G3)

The actual source curve U_L removes10+5k+5tl additional distinct affine
points. Its total puncture count and Euler characteristic are

    r=n+10+5k+5tl,
    chi(U_L)=2-2g-r=-14-5ek-5tl.                    (G4)

Thus the missing-sheet charge20+5ek+5tl, the affine missing-point
count10+5k+5tl and the different10+5(e-1)k are different quantities.
Their difference identity holds because every affine ramified point is
omitted. A ramification-only ledger would miss5l punctures in row A.
The degrees10,5k,5l here belong to a generic linear target projection.
For a specified coordinate direction, use its actual image-projection
degrees and infinity data; the P-coordinate direction has degrees6,3k,3l.
The inverse monicity budget(G0) uses that coordinate direction separately.

## 3. Genus and degree gates for the unknown second image

Since n has the parity required by (G2), the complete numerical genus
possibilities at a fixed image degree are

    e=2: g in {ceil(5k/2),ceil(5k/2)+1,ceil(5k/2)+2},
    e=3: g in {5k,5k+1,5k+2}.                        (G5)

Consequently every residual row has g>=3, with g>=5 in row C. The
cyclic-index3 row additionally requires g mod5 in{0,1,2}. Conversely,
if g is known from an actual source calculation, it constrains B by

    e=2: deg B is a multiple of5 in[2g-5,2g],
    e=3: deg B is a multiple of5 in[g-2,g].          (G6)

These are necessary numerical options, not existence assertions. In the
F2 degree-capped case, e2 rows have3<=g<=57, while row C has5<=g<=112
and the stated congruence. Row A has the sharper upper bounds55,52,30
in the distinct,H=C,H=B groupings respectively. The generic source Euler characteristics are

    A: -14-10k-5l, B,N: -14-10k, C: -14-15k.        (G7)

For example row A with H=B,k=2,n=2 has(g,r)=(7,32), whereas row C
with k=2,n=2 has(g,r)=(12,22). Both have chi(U_L)=-44 and total lost
sheet charge50. This numerical control is not a pair of global Keller
realizations. It demonstrates why Euler characteristic alone cannot
replace the separate ramification and omission columns.

The sum n=2 used in the excluded entire-single-prime r7 proof cannot
be carried into the new rows from that proof. There it follows from
the single affine branch degree10, Riemann--Hurwitz and the genus floor.
If an independent argument did give n=2 in a residual row, (G2) would
instead force g>=7; in e2 rows k would be even. Thus the old genus-two
adjoint-pencil reduction does not apply to the surviving rows.

The same genus applies to a generic target line through **any fixed
affine target value**, including a singleton zero. In the four rows
N is smooth and d(pi) has rank at least one everywhere: the local
blocks are unramified, fold, ordinary cusp or cyclic3. A generic line
through the fixed value therefore has smooth pullback at each of the
finitely many base points. Bertini handles the complement; infinity
is avoided as above, and ampleness gives connectedness. Its divisor
class is still the pullback of a line, so its smooth genus is the same.
This is the target radial pencil, not the source-coordinate fibration
T=x+c*y, whose finite fibres are A1 by PF2OMIT3. Rational exceptional
curves in the radial principalization do not make its general curve
rational. Equations(G3)-(G4) concern a fully generic line; a line through
a fixed exceptional value can have a different affine puncture count.

## 4. Exact Euler jumps of a target projection

Fix a generic linear target function h, and put H=h after F on the
source. For each omitted prime D, let h_D=h restricted to its A1
normalization, m_D=deg h_D and r_D(a)=#h_D^-1(a). Write
eta_D(q)=ell_D(q)-e_D, so eta is zero except at the one ordinary cusp
in rows B,N, where it is one. Euler integration of the finite original
source fibres over the target line h=a gives the exact identity

    chi(H^-1(a))=6-sum_D sum_(q:h_D(q)=a) ell_D(q),
    delta(a):=chi(H^-1(a))-chi(H^-1(generic a))
       =sum_D e_D*(m_D-r_D(a))
          -sum_D sum_(q:h_D(q)=a) eta_D(q).          (G8)

This formula does not require distinct projection-critical values.
It counts the original omitted points once, including when images
share a component. Source H has no critical points because F is etale;
the changes counted here take place at the omitted boundary.

For a generic direction all projection critical points are simple,
including the compulsory cusp point. Their values are distinct between
distinct image supports. The exceptions to simplicity are finitely many
inflection directions. Coincident tangent critical values correspond
to bitangent or common-tangent lines; characteristic-zero plane duality
makes these finite for nonline distinct irreducible supports. Shared
images retain their forced coincidences. Thus a smooth-image critical
point contributes Euler jump e_D; the ordinary cusp contributes2-1=1.
Let N_j count target projection values with nonzero Euler jump j. Then

| Row/grouping | N1 | N2 | N3 |
|---|---:|---:|---:|
| A, three distinct images | 5l-1 | 5k+8 | 0 |
| A, H=C | 0 | 5k-1 | 9 |
| A, H=B | 0 | 9 | 5k-1 |
| B or N | 1 | 5k+7 | 0 |
| C | 0 | 9 | 5k-1 |

These are nonzero Euler-jump counts, not a claimed classification of
all source fibration topology or a construction of its values. Row A
with H=B and row C have identical spectra: jump3 can be a fold plus
an unramified deletion at separate points, or a cyclic-index3 point.
It must not be used as the divisorial ramification index.

Summing(G8), polynomial Riemann--Hurwitz on each h_D gives

    sum_a delta(a)=sum_D e_D*(m_D-1)-sum_D E_D
                 =1-chi(H^-1(generic a)).           (G9)

The final equality is exactly the all-prime budget5=sum e_D+sum E_D.
The genus statement for generic lines through a fixed point above is
compatible with this jump: at the cusp value, the compact genus stays
the same and one affine puncture disappears, giving the Euler jump1.
At other projection-critical values the compact pullback can be singular,
so that fixed-point generic-direction argument is not a claim that
all parallel-line fibres have the same genus.

## 5. The first-image exceptional catalogue is now complete

Import PF2R7COLL1 and PF2R7KCOLL1/2, retaining the residual cube-root
choices. The two normalized coefficient charts have268 distinct
nonnodal points:34 on L plus237 on K minus3 intersection points.
There are163 coefficient points with triple values,99 with tacnodes
and6 with cusps. These are coefficient points, not affine-moduli classes.
At the individual singular values of the known first image, its
complete omitted factors are:

| Image type | Normalization points | Image-parameter fibre length | R ambient factors | Length remaining out of6 |
|---|---:|---:|---|---:|
| Node, tacnode or A5 | 2 | 2 | 2+2 | 2 |
| Ordinary triple | 3 | 3 | 2+2+2 | 0 |
| Ordinary cusp | 1 | 2 | 3 | 3 |

The existing matching formulas of PF2OMIT11/13 determine all allowed
additional-image completions. Every triple value is avoided by every
other omitted image. In the H=C grouping, all163 triple coefficient
points are excluded; the generic twelve-node opens and99 tacnode
points remain subject to global matching. The nonimmersed first-image
case belongs to row N and is treated separately.

Across only the268 nonnodal coefficient points, the singularity
**incidence** totals are2508 nodes,165 ordinary triples,99 tacnodes,
6 cusps and3 A5 values. The triple incidence count is165 rather than163
because L,d=0 has three triple values. Their delta check is

    2508+3*165+2*99+6+3*3=12*268=3216.

These totals are not counts of distinct points in one target plane.
They reconcile all exceptional catalogues and keep tangency conductor,
normalization multiplicity and ambient missing length separate.

The [replay](verify_generic_line.py) checks the source pin, all numeric
DVR profiles of total rank6, branch-cycle parity, genus options, Euler
and puncture conservation, shared-image controls and catalogue sums.
Universal Bertini/Hodge, Riemann--Hurwitz and the named degree-ratio
and collision imports remain written arguments. Actual B,H equations,
infinity data and source realization remain absent; no whole r7,
singleton-zero properness or JC2 theorem follows.
