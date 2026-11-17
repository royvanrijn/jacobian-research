# Cubic Gaussian frames and a terminating search for the least pure prefix

Current quantitative refinement, 15 November 2026: the [exact Euler representative proof](EULER_REPRESENTATIVES_GLOBAL_297.md) gives global pure297, odd296, conic/full-W-linear258 and all-homogeneous258. Ordinary-node complements are297/271/244, with homogeneous nonisotropic-node101; cusp255/tangent284 are unchanged. Current intervals are global[20,297], odd[20,296] and all-homogeneous[14,258]. Earlier353/336 bounds below are historical valid coarser evidence. Leastness,20/14 sufficiency, expanded I297 cofactors and independent review remain open.

Current quantitative refinement, 15 November 2026: the [joint-degree actual-seed proof](JOINT_DEGREE_GLOBAL_353.md) gives global pure353, odd352, full W-linear336, rank-one140, mixed-parent255 and cusp255/tangent284. Current least-prefix intervals are global[20,353], odd[20,352] and all-homogeneous[14,336]. Earlier590/888 bounds and rank-one425/mixed-parent320 estimates below are historical valid coarser evidence. Leastness,20/14 sufficiency, expanded I353 cofactors and independent review remain open.

Current quantitative refinement, 15 November 2026: the [balanced logarithmic seed proof](BALANCED_LOGARITHMIC_GLOBAL_590.md) gives global pure590, all-homogeneous336, conic336 and full W-linear425. Current intervals are global/odd[20,590] and all-homogeneous[14,336]; smooth[20,94] and homogeneous-smooth[14,56] are unchanged. Earlier888/1855 bounds and frame-plus888 statements below retain their historical scopes. Leastness,20/14 sufficiency, expanded I590 cofactors and independent review remain open.

Current global refinement,15 November2026: the [rescaled singular seed proof](RESCALED_SINGULAR_GLOBAL_888.md)
now proves888 sufficient as a global pure-only prefix, without an extra
frame condition. Current intervals are global[20,888], odd[20,888] and
all-homogeneous[14,888]. Earlier1855 bounds and reservations about888
alone below are historical. The new proof remains internal and
independently unreviewed; leastness and expanded global cofactors are open.

Current lower-bound refinement,15 November2026: the [odd-invariant proof](ODD_CUBIC_CLASSIFICATION_AND_PREFIX_LOWER_BOUND.md)
raises the least global pure-prefix interval to[20,1855] and proves
homogeneous lower14, by algebraic existence without displayed coefficients.
Earlier lower10 statements remain valid explicit-witness evidence.

Current quantitative status,15 November2026: the [global pure-prefix proof](GLOBAL_TERNARY_CUBIC_PURE_PREFIX.md)
now gives a sufficient1855 pure-only gate and least-prefix interval[10,1855].
Earlier absence-of-global-upper or unquantified-node statements below are
historical. Their earlier class/certificate scopes remain unchanged; the
new proof remains internal and independently unreviewed.

15 November 2026. Internal written algebraic deductions. The computability
result depends on the retained internal cubic theorem; it does not upgrade
that theorem's analytic assurance. No independent replay, external review
or formal verification. No full twenty-coefficient prefix search is run here.

## Results and exact boundary

**G3FRAMEDET20261115.** The existence of a Gaussian W-linear frame for
any ternary polynomial of total degree<=3 is equivalent to rank<8 of an
explicit8-by16 matrix with entries linear in its quadratic and cubic
coefficients. Its12870 maximal minors give equations of degree8 over Q
in circular Gaussian coordinates. This is a frame test, not a pure-null test.

**G3PURECOMPUTE20261115.** There is a terminating exact algorithm for
the least integer K>=0 such that the first K zero pure Gaussian moments
force all positive pure moments to vanish for every ternary polynomial
of total degree<=3. It also computes the least threshold on any specified
constructible coefficient family defined over the algebraic numbers.
Every insufficient K has an algebraic-coefficient witness in that family.
The algorithm returns checkable radical-membership certificates for
sufficiency and algebraic points for insufficiency.

The numerical global K is **not computed**. Its retained lower bound is10.
This numerical question is registered separately as **OP-GMC3-PURE-PREFIX**.
The subsequent [weighted binary class gates](SQUAREFREE_SINGULAR_FINITE_PREFIXES.md)
quantify concurrent/cuspidal/tangent leading curves by104/720/935; the
ordinary-node leading types still need quantitative bounds.
Neither888 alone nor10 is asserted sufficient globally. The computability
deduction uses the [full internal cubic theorem](TERNARY_CUBIC_GMC_CLOSURE.md),
its frame-plus888 certificate, and standard exact polynomial algorithms.
It is separate from ideal-generation cutoffs, nilradical indices, primary
decompositions, global term optimality and lower-power GVC.

## 1. A rational parametrization of the null directions

Use circular coordinates (Z,W,T) with E(ZW)=E(T^2)=1 and all other
quadratic contractions zero. The metric on directions is

    B(v,v)=2*vZ*vW+vT^2.

Its projective null conic is parametrized, without irrational constants, by

    v(s,t)=(s^2,-t^2/2,s*t),  [s:t] in P^1.       (1)

This is the usual Veronese conic followed by an invertible linear change.
It has no zero vector for [s:t] in P^1 and covers the entire null conic,
including Z and W directions at the two endpoints.

Write P=F3+F2+F1+F0 in homogeneous pieces. Form four binary quartics

    f0=(partial_Z F3)(v), f1=(partial_W F3)(v),
    f2=(partial_T F3)(v), f3=F2(v).              (2)

Zero quartics are allowed. By the retained exact frame criterion, a frame
exists if and only if these four forms have a common projective root.
Ordinary partial derivatives suffice: changing the metric raises the
gradient by an invertible matrix and does not change its zero set.
Equivalently D_v^2P=0, since

    D_v^2P(x)=2*sum_i (partial_i F3)(v)*x_i+2*F2(v).

## 2. The finite matrix and its complete common-root proof

Let S_d be binary homogeneous forms of degree d. Define

    M: S_3^4 -> S_7, (h0,...,h3) |-> sum_i hi*fi. (3)

Choose descending monomial bases. The columns are fi*s^(3-j)*t^j for
i=0,...,3 and j=0,...,3, so M is8-by16. Its entries are rational
linear forms in the16 coefficients of F3,F2. Lower coefficients do
not enter the frame test.

If the fi share a projective root, evaluation there is a nonzero linear
functional on S_7 annihilating the image of M. Thus rank M<8.

Conversely, suppose they have no common root. Choose a nonzero quartic
a in their span. It has only finitely many projective roots. At each root
some fi is nonzero; a generic linear combination b avoids all these
finitely many proper hyperplanes in the coefficient space. Thus a,b are
coprime degree4 forms in the same span. The map

    S_3^2 -> S_7, (h,k) |-> h*a+k*b

is injective: a divides k*b, hence divides k, which has degree3<4, so
k=0 and then h=0. Both sides have dimension8, so it is surjective.
Its image is contained in the image of M. Hence rank M=8.

This proves the equivalence in every specialization, including zero
quartics and repeated roots. No generic-rank assumption is made. A common
root at infinity is retained because the forms and matrix are homogeneous.
There are binomial(16,8)=12870 maximal minors, each homogeneous of degree8
in the coefficients (some can vanish identically). Clearing rational
constants changes neither the frame locus nor the rank test.

For an individual algebraic input, a homogeneous gcd of the four quartics
is usually cheaper than enumerating minors. Its positive degree is an
exact existence certificate over C. It does not require that a frame
direction lie in the input's own coefficient field.

## 3. A finite target ideal in the full twenty-coefficient space

Let R=Q[c_alpha: |alpha|<=3], and let P=sum c_alpha*Z^alphaZ*
W^alphaW*T^alphaT. There are20 coefficients. Wick contraction gives
M_j=E(P^j) in R effectively, with degree j; the circular normalization
has rational (indeed integral) contractions.

Let J be the ideal generated by the maximal minors of (3), and put

    T=J+(M1,...,M888),    I_K=(M1,...,M_K),
    I_infinity=(M_j: j>=1).                     (4)

The internal cubic theorem says exactly V(T)=V(I_infinity) over C.
Consequently the [Hilbert Nullstellensatz](https://stacks.math.columbia.edu/tag/00FS)
gives

    radical(T)=radical(I_infinity).             (5)

This equality is about reduced zero sets. It asserts neither
T=I_infinity nor that the moment ideal is generated by moment888.

All tests below can be performed over Q: extending a field is faithfully
flat, so an element's power belongs to an extended ideal precisely when
it belongs to the original ideal. For algebraic coefficient families,
replace Q by an explicitly presented number field containing the data.

## 4. A search that halts, and then finds the least threshold

List the finitely many generators g of T in (4). For K>=888, compute
I_K and test g in radical(I_K) for every g, using

    g in radical(I_K) iff 1 in (I_K,1-u*g) in R[u]. (6)

Equation (6) follows directly by localizing (R/I_K) at g: that localization
is the zero ring exactly when some power of g is zero. A Groebner basis
decides unit-ideal membership and can retain its polynomial cofactors.

Increase K until all tests succeed. This process halts. By (5), each of
the finitely many g has g^e in I_infinity for some finite e. A polynomial
ideal membership is a finite sum, so it uses only finitely many moments.
Taking the largest of these finitely many indices gives a finite K at
which every test succeeds. This is a termination proof, not a proposed
test based on a long run of unchanged ideals.

If N is the first sufficient value found, test the same finite target
at K in[0,N], by a monotone binary search. At every K, the inclusion
I_K subset radical(T) holds by the all-order classification. Therefore
all generator tests in(6) succeed exactly when V(I_K)=V(T). The first
successful K is the least global pure threshold. K=0 means no moment
hypothesis is necessary; it is included for constrained families.

For insufficiency, some g is not in radical(I_K). The proper ideal
(I_K,1-u*g) has a point over the algebraic closure of Q by the same
Nullstellensatz. Its coefficients supply a polynomial with its first K
moments zero but g nonzero, hence not pure-null. Exact polynomial-system
algorithms can produce a closed point and present its coordinates by
minimal polynomials and isolating data. The source theorem implies such
a point has some later nonzero moment; enumerating its exact moments
therefore eventually produces an explicit failed order as well.

These are all finite, effective operations. Their cost can be enormous.
No full moment expansion through888, global minor ideal, or global
twenty-variable radical test has been executed in this follow-up.

## 5. Constructible coefficient constraints

A basic constructible family is specified by polynomial equations H=0
and nonvanishing conditions h1*...*hr!=0 over a number field. Adjoin y
and impose y*(h1*...*hr)-1=0. Work in this finite-type coefficient
algebra, add H and this inverse equation to every ideal in(6), and use
the same fixed target generators. Equality (5) and its finite-support
argument survive this restriction. Thus the same search terminates.

A finite union of basic constructible sets is handled stratum by stratum;
the least uniform threshold is the maximum of their least thresholds.
Empty strata have threshold0. Insufficient thresholds again have
algebraic points, now satisfying every specified equation and open condition.
This supplies computable exact finite gates for algebraically specified
cubic supports and subfamilies, without claiming their numerical values.

The distinction is the same one needed in finite-threshold discussions
such as [Liu's v1 two-variable record](https://zenodo.org/records/21674142):
computable does not mean that a sharp numerical bound has been evaluated.
No two-variable theorem is used to establish the ternary result here;
the new application rests on the retained ternary classification.
Only that v1 record's abstract is used for this comparison. Its latest-version
link resolves to the differently titled [quadratic v3 record](https://zenodo.org/records/21688251);
the older abstract is not presented as the current revision's contents.

## 6. A completed small family and the plateau trap

On P=a*(Z^2*T+W^3), the retained sharp cusp formula gives M1..M9=0
identically and M10=1508855040000*a^10. Its only pure-null point is a=0.
The frame matrix at a=1 has rank8, so its minors on this line generate
(a^8), up to nonzero rational multiples. The target radical is(a).

Thus the radical tests fail through K=9 and succeed at K=10. They
recover the sharp threshold10 on this particular line. They also show
why checking I_K=I_(K+1) is not a stopping criterion: I1=...=I9=(0),
yet I10 is nonzero. Further, a^8 is in radical((a^10)) but not in
(a^10). A set-theoretic finite gate cannot establish an ideal-generation
cutoff or scheme thickness.

The [checker](verify_cubic_frame_prefix.py) and
[receipt](cubic_frame_prefix_verification.json) verify the null conic,
directional derivative identity, matrix dimensions and exact ranks/gcds
in stated specializations, rational maximal-minor certificates for the
no-frame examples, and the completed one-parameter radical/plateau tests.
The universal common-root proof and global termination are written
arguments, not extrapolations from these regressions. The
[source audit](cubic_frame_prefix_source_audit.json) retains the dependency
and source boundaries.
