# Global collision chords and compulsory defect drops in every pencil

Date: 2026-11-16. Claims HC4MYLB0 and HC4MYLB1.

FULL THEOREM TARGET: prove HC4/JC2 or construct an actual counterexample.
CURRENT LAST UNSOLVED IMPLICATION: the ACTUAL polynomial plane source forces
Tr(Lambda_L*Lambda_L) subset B, equivalently regular inverse Gram entries;
the alternative full marked-energy integrality implication is also OPEN.
This note does not prove either implication. External review/outreach
remains DEFERRED.

Unlike the [LA deck-group restriction](GENERIC_DECK_GROUP_AND_UNRAMIFIED_BOUNDARY_COST.md),
this argument needs NO nontrivial generic deck group. It applies to EVERY
hypothetical noninvertible complex plane Keller map in every field degree.
The [KW source-line conductor theorem](SOURCE_LINE_CONDUCTOR_AND_PRIMITIVE_ORDER_COLLISIONS.md)
already forces positive affine normalization defect on every line. We put
all finite line collisions into one smooth affine surface, retaining their
SCHEME multiplicities, and obtain a stronger global constraint.

## 1. Every line is its image's finite immersive normalization

Let K=(P,Q):X=A^2->Y=A^2 be complex polynomial Keller and noninvertible.
For ANY affine source line ell choose a parameter t. At least one of
p=P|ell,q=Q|ell is nonconstant, so t is integral over C[p,q]. The finite
normal image has normalization A^1: as in the retained KW/PF53 argument,
the polynomial parametrization has rational function field and only one
place at infinity. Factoring through that normalization gives a polynomial
t->r(t). Since the invertible Keller differential makes (p',q') nowhere
simultaneously zero, r' cannot vanish. Thus r is linear and

    A_ell=C[p,q] subset R_ell=C[t]

is a finite BIRATIONAL normalization, immersive at every t. This repeats
the retained argument for arbitrary line orientations; it does not assume
one fixed monic presentation works in all orientations.

Define the total AFFINE normalization defect

    delta(ell)=dim_C(R_ell/A_ell).                          (1)

Every individual image branch is smooth. Singularities can be nodes,+tangencies, or multiple-branch points; no generic-node assumption has
been made. The existing Gwozdziewicz one-line injectivity theorem, imported
in PF52, gives delta(ell)>=1 for EVERY ell under noninvertibility.

## 2. LB0: the smooth affine collision surface and its chord map

Form W=X times_Y X. It is affine, smooth, and pure dimension2, since
the first projection is a base change of the etale map K. Its diagonal
Delta is OPEN because K is unramified and CLOSED because K is separated.
Thus the OFF-DIAGONAL C=W minus Delta is also affine, smooth and pure
dimension2. Every connected component is an irreducible smooth surface.
It is nonempty under noninvertibility. No finite-cover assumption on
K or on either projection C->X is made.

Swapping the two distinct points is a FREE involution. Its affine
quotient D=C/S2 is smooth affine, and C->D is finite etale of degree2.
These are the same finite free-action quotient foundations retained in
KMFI47; no assertion that D is an affine plane is involved.

Let L_line=P^2_dual minus {[0:0:1]}, the variety of ALL affine source
lines. To a pair p=(x1,y1),q=(x2,y2) attach its unique chord,

    [a:b:c]=[y1-y2:x2-x1:x1*y2-x2*y1].                    (2)

The first two entries are not both zero. They therefore define a
morphism C->L_line, invariant under the swap, and hence a morphism

    chi:D->L_line.                                       (3)

For ell=[a:b:c], its ordered scheme fibre is exactly

    (ell times_(K(ell)) ell) minus its diagonal.           (4)

This is an equality of SCHEMES, not only of closed point sets. For
example, near a vertical chord x=c0 the difference y1-y2 is invertible.
The chord coordinates are m=(x1-x2)/(y1-y2), b=x1-m*y1.
The ideal (m,b-c0) is exactly (x1-c0,x2-c0). The other chart is analogous.
Thus fixing the chord is the same as imposing both line equations on C.

The finite normalization in Section1 has only finitely many nontrivial
identifications, so (4) is a finite zero-dimensional scheme. The map chi
is therefore quasi-finite and separated. It is SURJECTIVE, because EVERY
source line has a collision by the existing one-line theorem.

The EXACT fibre length is

    length_C chi^-1(ell)=delta(ell).                       (5)

To verify the multiplicities, complete at an image singular point with
normalization preimages t1,...,tr. The branches are smooth. For i!=j the
ordered local fibre ring is

    C[[s,t]]/(p_i(s)-p_j(t),q_i(s)-q_j(t)),

whose length is the plane-branch intersection multiplicity I_ij. The
diagonal contributes none, since the normalization map is unramified.
For a reduced plane curve with smooth branches its delta invariant is
sum_(i<j) I_ij. One can derive this directly from
0->O_(A union B)->O_A direct_sum O_B->O_(A intersect B)->0 and induction;
intersection orders add when another smooth branch is compared with a
product of branch equations. Thus the ordered fibre has length 2*delta.
The finite etale swap quotient halves that length, proving (5).
In particular a TACNODE contributes length2 to (5), not1, while an
ordinary triple point contributes3 unordered pairs. Finite singularities
have not been confused with omitted points or points at infinity.

Every component of D dominates L_line: quasi-finiteness and dimension2
force a dense image in that irreducible surface. Miracle flatness (the
ALREADY IMPORTED Stacks tag00R4, used in KU1) makes chi FLAT, because
both surfaces are smooth and its fibres are zero-dimensional. It need
NOT be etale. Write N for its TOTAL generic degree, summing all components.
On a dense open it is finite etale, so N is the generic value of (1).

There is also a SECOND, different cover carried by the ACTUAL collision
data. The common value K(p)=K(q) defines e:D->Y, which is etale by descent
from C->Y through the finite etale swap quotient. Its total generic degree
is binomial(d,2), where d is the Keller field degree. Thus D has the exact
coordinate coframe e^*dP,e^*dQ and exact nonvanishing volume. The two cover
degrees N and binomial(d,2) must NOT be equated. If generic monodromy is
transitive on unordered pairs, D is irreducible, by the generic-cover orbit
decomposition. This includes the natural symmetric/alternating packets.

Finally, chi^*O_(P2_dual)(2) is TRIVIAL on D. Indeed (2) is a nowhere-zero
trivializing vector for the pulled-back tautological line on C, and the
swap negates that vector. Its square is invariant and descends to a
nowhere-zero section. Only order DIVIDING2 is asserted, not exact order2.
For the pencil through a=(a_x,a_y), the regular swap-odd function on C

    f_a=(x1-a_x)*(y2-a_y)-(x2-a_x)*(y1-a_y)

has zero divisor the pulled-back pencil. Therefore its invariant square
descends to g_a in O(D) with EXACT divisor

    div(g_a)=2*D_a.                                      (5a)

This includes all Cartier multiplicities and supplies an additional
actual-collision constraint beyond an abstract affine covering of L_line.
It is not an assertion that D itself is A^2 or that its Picard group is
torsion-free.

Zariski Main, already used in KMFI94 and the normalization programme,
embeds D openly in the full normalization Z_ch of L_line in its component
function fields. This normalization is finite over L_line. Its normal
surfaces are Cohen-Macaulay, so the same miracle-flatness argument makes
Z_ch finite LOCALLY FREE of total rank N. We have, at EVERY line,

    1<=delta(ell)<=N,
    N-delta(ell)=sum of lengths of the omitted local
                 Artinian factors in (Z_ch)_ell.          (6)

The exact factors are those of the FULL normalization fibre; no primitive
order or ordinary branch-count substitute is used. The source D is normal,
which is essential to this normalization factorization.

The omitted set Z_ch minus D is purely divisorial. To see why an isolated
omission is impossible, take a normal neighbourhood of such a point
avoiding all divisorial omissions. All coordinate functions of the AFFINE
variety D extend across the puncture by normality. They define a map of
the neighbourhood into D whose composition with D->Z_ch agrees with the
identity on a dense open and hence everywhere. The omitted point would
belong to D, a contradiction. This is the same affine-open/Hartogs mechanism
retained in KMFI94, now on the collision normalization; no affine-plane
source hypothesis is transferred to D.

Consequently the set S_ch of lines with delta(ell)<N is the closed
nonproperness locus of chi and, when nonempty, is a union of curves.
Neither chi nor its restriction to a pencil is assumed proper.

## 3. LB1: every complete source pencil has a strict defect drop

For a finite source point a, ALL affine lines through a form a closed
PENCIL P_a=P^1 inside L_line: its dual linear equation does not contain
the removed point for the line at infinity. Set D_a=chi^-1(P_a), with
its WHOLE scheme structure. It is a closed affine curve, possibly
nonreduced. Since chi is flat, D_a->P_a is flat, quasi-finite and
surjective. Equivalently it is a Cartier divisor in the smooth surface
D, Cohen-Macaulay and with every component dominating P_a. Its fibres
have lengths delta(ell) by (5).

The following elementary consequence of Zariski Main is useful here.
If an AFFINE scheme T of pure dimension1 is flat, quasi-finite and
surjective over P^1, with generic degree m, then

    every fibre has length <=m,
    some fibre has length STRICTLY LESS than m.            (7)

Here are the details, including nonreduced T. Factor T as an open of a
finite scheme Z->P^1 by Zariski Main. Replace Z by the scheme-theoretic
closure of T, then divide its finite O_(P1)-algebra by its base torsion.
This does not change T, which is flat over P^1. The resulting finite
algebra is torsion-free over the smooth curve P^1, hence locally free of
rank m. Thus T is an open of a finite FLAT degree-m cover Z'. At a fibre,
its Artinian local factors are a subset of those of Z', proving the first
inequality. T cannot equal Z': Z' is projective of positive dimension,
whereas T is affine. Some point is omitted, and its positive local fibre
length gives the strict inequality. No normalization of T is taken:
normalizing a nonreduced or singular pencil could lose its required lengths.

For an exact primary reference to the ALREADY USED Zariski Main
foundation, Stacks Lemma37.43.3, tag05K0, states this factorization for a
quasi-finite separated morphism to a quasi-compact quasi-separated base.
The actual statement and proof were read and cached. P^1 and L_line meet
these hypotheses. The theorem and miracle flatness are retained proved
inputs, not newly proved here or inferred from symbolic examples.

Apply (7) to D_a. If m_a is the generic affine line-image defect in the
pencil through a, then

    m_a>=2,
    1<=delta(ell)<=m_a for EVERY ell through a,
    some ell through a has delta(ell)<m_a.                 (8)

The first inequality follows from the last one and delta>=1. These
statements hold for EVERY finite source point a, without a bound on
field degree, line-image degree, genus at infinity, singularity type,
collision multiplicity, or number of components of D. In particular

    a ZARISKI-GENERAL source line has delta>=2.             (9)

Indeed a general complete pencil meets the dense open where delta=N,
so its generic m_a=N, and (8) gives N>=2. Alternatively, N=1 would make
D a single birational component and Zariski Main would identify chi with
an open immersion into the normal L_line; its surjectivity would identify
affine D with NONAFFINE L_line, which contains complete P^1s.

The loss set S_ch in Section2 is therefore NONEMPTY and meets EVERY
complete pencil. Some lines may have defect1; it is the GENERAL member
of EACH complete pencil that has defect at least2. A fixed PARALLEL
family is an A^1, not a complete P^1 in L_line, and is not subject to
the strict-drop argument without a separate infinity analysis.

For a general admissible linear source ruling of KW0, its conductor
polynomial h has general fibre degree 2*delta>=4. KW0's every-special-line
positive-even bound remains2 at minimum; no degree>=4 claim is made for
every line or every fixed exceptional direction. This strengthens the
universal generic finite-collision requirement, not the Keller field-degree
frontier and not the current Laurent-key finite-target compatibility gap.

## 4. Exact hostile controls

The following are PLANE CURVE controls, not whole Keller maps.

The nodal parametrization p=t^2-1,q=t*(t^2-1) has a single collision
between t=1,-1. Its divided-difference equations reduce to s=-t and
t^2-1=0: ordered length2, unordered length1, delta1.

The tacnodal parametrization p=t^2-1,q=t*(t^2-1)^2 has the SAME two
distinct points over its singular image. Its divided-difference equations
reduce to s=-t and (t^2-1)^2=0: ordered length4, unordered length2,
delta2. Its implicit curve is q^2=p^4*(p+1). The full normalization is
immersive, including at those points. A locally etale TWO-DIMENSIONAL
extension along the entire source line y=0 is

    K0(x,y)=(x^2-1-y,
             x*(x^2-1)^2+y*(3*x-(5/2)*x^3)).

Its Jacobian is 1+y*(3-(15/2)*x^2), equal to1 on the whole line, but
NOT a global nonzero constant. The two image tangent vectors at x=1,-1
are parallel. This shows why local invertibility at the collision points
does not justify declaring the chord map etale or counting distinct
pairs instead of their scheme lengths.

The immersive triple-point parametrization p=t^3-t,q=t*(t^3-t) has
three preimages -1,0,1 over its singular image, with distinct tangent
directions. It has six ordered pairs, three unordered pairs, and delta3.
Thus a higher-multiplicity image point must not be mistaken for one pair.

For the necessity of scheme lengths and flatness in (7), a double cover
of curves can merge two distinct points into one retained point while its
fibre LENGTH stays2. This is ramification, not an omitted sheet. Conversely
Gm->A1,t->t^m is quasi-finite flat but not surjective at zero; it cannot
replace the actual surjective pencil map. The affine target A1 itself
admits finite flat covers of constant positive degree, so it cannot be
substituted for the COMPLETE pencil P1 in the nonproperness argument.

The lower bound2 is SHARP for the abstract covering argument alone.
The disjoint union of the two standard affine charts of P^1 maps
surjectively and etale to P^1, has generic degree2, and drops to length1
at each of the two missing opposite-chart points. At surface level, the
disjoint union of the direction charts {a!=0} and {b!=0} of L_line consists
of two affine planes and gives the same generic degree2 behaviour in
EVERY complete pencil. These are not collision surfaces of Keller maps:
they supply neither the two source-point maps nor their common K-value.
Thus the affine-cover argument alone cannot improve2 to an arbitrarily
large bound or a contradiction without additional ACTUAL collision data.

## 5. Scope of the advance

LB0/1 applies even to trivial generic deck group and therefore to the
surviving symmetric/alternating monodromy packets. It closes the candidate
generic-delta1 collision geometry in EVERY field degree and proves a
compulsory loss of finite collision multiplicity somewhere in EVERY complete
source pencil. It does not make chi etale or finite, forbid defect>=2,
bound N or m_a above, or prove that one line is injective. It does not
prove that all these collision losses occur over finite target values;
the common K-value may escape to infinity too. No direct contradiction
with the entire nonproperness curve of K is inferred from such a loss.

The last universal trace-dual regularity/marked-energy implication remains
OPEN. The checker verifies exact chord-ideal identities, the curve controls,
and provenance/scope; it does NOT certify the global algebraic-geometric
proof or HC4. Independent replay, formal verification and external review
remain false.
