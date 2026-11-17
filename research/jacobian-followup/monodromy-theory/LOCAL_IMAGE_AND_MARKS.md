# Local monodromy controls the actual affine image

Divisorial inertia tests one group element at a time. At a branch
intersection or singularity, an entire local subgroup must fix the
same sheet if the point is in the image of a Keller map.
For the two exceptional PSL2(11) actions this gives different answers:
one class of S3 subgroups fixes two sheets in one action and none in
the other, although their elementwise permutation characters agree.

These are necessary image and arithmetic constraints. No exceptional
degree-eleven realization or absolute exclusion is supplied.

## 1. Pure boundary and local fixed sheets

Work over algebraically closed characteristic-zero k. For a Keller
map F:X=A^d -> Y=A^d let L/K be its generic extension, M/K its
Galois closure, G its group and Omega=G/H its sheets.
Let nu:Z -> Y be the finite normalization in L.
Zariski main identifies X with an open subset of normal affine Z.

For a geometric target point y, put R=O_(Y,y)^sh and K_y=Frac(R).
A component of M tensor_K K_y has Galois group J_y, defined up
to conjugacy in G. This is the **strict local monodromy group**.
At a geometric generic point of a divisor it is inertia.
It is not the codimension-one decomposition group before closing
the residue field.

**KMFI94.** The whole complement Z-X is a union of prime divisors;
there are no additional omitted components of codimension>=2.
Let R_nu be the reduced support of all divisorial ramification in Z.
There is a union T of omitted unramified prime divisors such that
\[
                   X=Z\setminus(R_\nu\cup T).               \tag{1}
\]
For every geometric y there is an injection
\[
                   F^{-1}(y)\hookrightarrow\Omega^{J_y},
        \qquad \#F^{-1}(y)\le|\Omega^{J_y}|.                 \tag{2}
\]
If T is empty, equality holds for every y.
In particular J_y with no common fixed sheet forces y outside F(X).
Conversely a fixed sheet alone does not guarantee a source point.

Here are the proofs, including the boundary qualifications.
Remove all divisorial components of Z-X and call the remaining
open U. Then U-X has codimension>=2. Normality gives
O(U)=O(X)=k[x_1,...,x_d], by height-one extension.
The corresponding map U -> X restricts to the identity on X.
Composing with X -> U gives the identity on the dense open X,
hence on all U by separatedness. Thus U=X.
This is also the complement-of-affine-open result in
[Stacks, Tag0BCQ](https://stacks.math.columbia.edu/tag/0BCQ).

The nonetale locus of nu is exactly R_nu:
[Stacks, Tag0BMB](https://stacks.math.columbia.edu/tag/0BMB)
applies at any point not lying on a ramified height-one
generalization, since Z is normal and Y regular.
Every ramified divisor is omitted from the etale source.
The pure-boundary assertion therefore gives (1).

Over strictly henselian R, each source point determines a unique
etale section. The henselian decomposition of a quasi-finite algebra
separates its finite part from the part missing the closed fiber
([Stacks, Lemma10.153.3](https://stacks.math.columbia.edu/tag/04GE)).
Different source points give different generic sections.
Meanwhile L tensor_K K_y has factors indexed by the J_y-orbits
on Omega. Its degree-one factors are precisely the fixed sheets.
Their normal integral models are copies of R; they are exactly the
etale normalization points over y. This proves (2).
If T is empty, (1) is the full etale locus, so all those points
are present and equality follows.

KMFI87 gives another useful form of the equality hypothesis:
if rank Cl(Z) equals the number of ramified normalization primes,
there are no additional omitted unramified prime divisors.
Then the exact fiber count is the local fixed-sheet mark.

The image is open and its complement has codimension>=2, by
etaleness and the divisorial surjectivity in KMFI1.
Omission of a higher-codimension target stratum is consequently
allowed. This does not conflict with the purity of the **source
normalization boundary** proved above.

## 2. Complete local subgroup marks for exceptional PSL2(11)

Fix the two nonisomorphic eleven-sheet actions Omega_a,Omega_b.
The suffix a on A5 means the stabilizer class for Omega_a.
The suffix a on S3 means the class that fixes two Omega_a sheets.
D10 and D12 denote groups of orders10 and12.

**KMFI95.** The complete subgroup mark table is:

| Subgroup class J | Number in G | Fixed sheets in Omega_a | Fixed sheets in Omega_b |
|---|---:|---:|---:|
| 1 | 1 | 11 | 11 |
| C2 | 55 | 3 | 3 |
| C3 | 55 | 2 | 2 |
| V4 | 55 | 1 | 1 |
| C5 | 66 | 1 | 1 |
| S3a | 55 | 2 | 0 |
| S3b | 55 | 0 | 2 |
| C6 | 55 | 0 | 0 |
| D10 | 66 | 1 | 1 |
| C11 | 12 | 0 | 0 |
| A4 | 55 | 1 | 1 |
| D12 | 55 | 0 | 0 |
| F55=C11 semidirect C5 | 12 | 0 | 0 |
| A5a | 11 | 1 | 0 |
| A5b | 11 | 0 | 1 |
| G | 1 | 0 | 0 |

There are620 subgroups, in16 conjugacy classes.
For Omega_a a subgroup has no common fixed sheet if and only if
it contains a conjugate of C6, C11 or S3b.
For Omega_b replace S3b by S3a.
These are the complete inclusion-minimal obstructions:55+12+55=122
subgroups in each action. They are local subgroup obstructions,
not exclusions of those subgroups as abstract global monodromy groups.
The table does not assert that every listed subgroup is realizable
as a strict local group: its divisorial inertia must also satisfy
KMFI1 and generate the local group.

The S3 orbit profiles make the difference explicit:
\[
 \begin{array}{c|cc}
 &\Omega_a&\Omega_b\\ \hline
 S3a&(1,1,3,6)&(2,3,3,3)\\
 S3b&(2,3,3,3)&(1,1,3,6).
 \end{array}                                                \tag{3}
\]
Every element of either S3 subgroup has a fixed sheet in both
actions. A common fixed sheet is a stronger requirement.
Averaging the elementwise fixed-point counts gives the number
of orbits, four in both rows, not the number of common fixed
sheets. Equal permutation characters do not determine these marks.

The finite proof is exhaustive, not an appeal to a subgroup database.
Use all660 canonical SL2(F11)/{+I,-I} matrices and their multiplication
table. Element orders are1,2,3,5,6,11.
For subgroups containing neither order6 nor11, an odd-order
subgroup has order1,3,5 or15. The order15 case is cyclic by Sylow
counting and impossible in G.
An even-order subgroup contains an involution; all55 involutions
are conjugate. Starting with one C2, enumerate every one-element
extension, discarding an extension as soon as it contains an
order6 or11 element. The property is inherited by subgroups, so
this pruning cannot lose an eligible overgroup.
There are10244 tested extensions and exactly25 retained overgroups:
one C2, three V4, six S3, six D10, three A4 and six A5.

For subgroups containing an element of order6, the complete
overgroups of one C6 are C6,D12,G (1302 tested extensions).
For those containing order11, the complete overgroups of one C11
are C11,F55,G (1254 tested extensions).
Conjugating these finite lists and adding C3,C5,1 exhausts all
possibilities and yields the620 groups above.
The checker evaluates both fixed-sheet marks for every subgroup
and tests inclusion-minimality over the full list.

## 3. Normal-crossing and cusp strata

**KMFI96.** At a point where the **whole** branch divisor of the
normal closure is simple normal crossing, the strict local group
is abelian and generated by its branch-meridian inertias.
For exceptional11 its possibilities and fiber bounds are
\[
 \begin{array}{c|rrrrrr}
 J_y&1&C2&V4&C3&C5&C6\\ \hline
 \#F^{-1}(y)\text{ at most}&11&3&1&2&1&0 .
 \end{array}                                                \tag{4}
\]
In particular a normal-crossing stratum carrying both inertia2
and inertia3 is omitted from the affine image.
Inertia5 cannot meet inertia2 or3 at such a stratum at all:
commuting elements would give order10 or15, absent from G.
This concerns the full branch support; two components crossing
transversely in the presence of other singular branch components
do not satisfy this hypothesis.

To prove local abelianity, in the strictly henselian regular local
ring write the normal-crossing components as parameters t_i.
Adjoin common multiples of their tame inertia roots t_i^(1/N).
The resulting local ring is regular and strictly henselian.
The DVR Abhyankar lemma kills all height-one ramification;
purity makes the normalized cover etale and hence split.
Thus J_y is a quotient of the abelian Kummer group.
This is the iterated parameter version of
[Stacks, Tag0EYG](https://stacks.math.columbia.edu/tag/0EYG).
Inertia11 is excluded already by KMFI1. The full abelian subgroup
census gives exactly the groups in (4).

At a cusp with local subgroup S3a, Omega_b forces omission while
Omega_a allows at most two points; for S3b the roles reverse.
This statement assumes that the cusp's actual local group is S3;
a cusp shape alone does not specify its finite monodromy quotient.
The local decomposition (3) accounts for all normalization branches.

### An actual Keller map omits a smooth complete-intersection curve

The known cubic S3 map KMFI11,88 is a control against turning
forced omission into an absolute action exclusion.
Write u=1+xy, g=1-(3/2)xy+x^2z and
\[
 A=(u+u^2-2u^3g)/x^2,\quad
 B=(1+2u-3u^2g)/x,\quad C=xg .
\]
These are polynomials with whole Jacobian1.
Its finite normalization has the original A3 as its full etale
locus (KMFI88). The omitted target is **exactly**
\[
             \Sigma:\quad 3BC-1=0,\qquad27AC^2-1=0.          \tag{5}
\]
This smooth curve is isomorphic to G_m, with parameter C:
\(B=1/(3C), A=1/(27C^2)\).

For C!=0 the finite cover is represented by the projective roots
of \(P(R)=CR^3-R^2+BR-A\). It has a simple root unless all three
coincide. The triple-root equations are exactly (5).
At C=0 the root at infinity is always simple.
This proves the exact image assertion using the full etale source.
Translating R by1/(3C) gives a depressed cubic with coefficients
\[
  p=(3BC-1)/(3C^2),\qquad
  q=(3(3BC-1)-(27AC^2-1))/(27C^3).
\]
These are transverse local coordinates at Sigma.
The ordered-root cover has local group S3, acting transitively
on the three roots; it has no common fixed sheet.
Thus KMFI94 detects this actual omitted stratum.

There is also a global polynomial check. Put
P0=3BC-1, Q0=27AC^2-1 after substituting the source map,
S=3ug-1 and L=uS-1. Direct identities give
\[
 P0=3g-S^2,\qquad Q0=S^3+3P0(S+1),\qquad1-uP0=SL,
\]
and hence
\[
 1=P0[3u-3u^2P0+u^3P0^2-3(S+1)L^3]+Q0 L^3.              \tag{6}
\]
The two pulled-back equations are unimodular despite the target
curve being nonempty, smooth and a complete intersection.
A rule forbidding such omissions would contradict this actual map.

## 4. The exceptional action has no larger arithmetic normalizer

**KMFI97.** For either exceptional action,
\[
                         N_{S_{11}}(G)=G.                 \tag{7}
\]
Consequently any characteristic-zero Keller map with this geometric
generic action has the same arithmetic group PSL2(11), and its
generic normal closure has no additional algebraic constants.

This also has a complete small-group proof. Fix generators a,b
of orders2,3 with ab of order11. There are55 choices for an
involution image and24 choices for the order3 image with product
order11. Every automorphism is determined by these images, so
|Aut(G)|<=1320.
The660 inner automorphisms are distinct since the center is trivial.
Conjugation of matrices by diag(2,1), of nonsquare determinant in
F11, is an automorphism theta interchanging the two A5 classes.
It is not inner and supplies another660 automorphisms.
Thus Aut(G) has exactly1320 elements; all outer automorphisms
interchange those classes.

A permutation normalizing an eleven-sheet action must preserve
the conjugacy class of its point stabilizer. It therefore induces
an inner automorphism of G. The permutation centralizer is trivial:
the point stabilizer fixes only its one distinguished point, so
a commuting permutation fixes that point and then every point by
transitivity. This proves (7), without enumerating S11.

Over an arbitrary characteristic-zero coefficient field k,
the geometric group is normal in the arithmetic permutation group.
The latter lies in (7), hence equals G.
Its quotient by the geometric group is the Galois group of the
normal-closure constant extension, so that extension is trivial.
This is an arithmetic necessity for a hypothetical realization,
not a construction.

## Replay and boundaries

Run

    .venv/bin/python research/jacobian-followup/monodromy-theory/verify_local_image.py

[Helper](local_image.py) · [model](local_image_model.json) ·
[receipt](local_image_receipt.json) ·
[pinned primary inputs](sources/local-image-manifest.json).

The checker exhausts every subgroup and both sets of marks,
all minimal bad subgroups, abelian strata and the automorphism
image bound. It replays the cubic full Jacobian, triple-root
locus and polynomial omission certificate.
The all-dimensional boundary and strict-local arguments have
written proofs using the pinned primary inputs.
Formal verification, independent whole-proof replay and external
review are not claimed. Exceptional11 occurrence and the complete
faithful transitive Keller action spectrum remain open.
