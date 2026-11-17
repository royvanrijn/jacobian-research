# Canonical etale completion: boundary topology and free deck action

For every actual complex polynomial plane Keller map, the canonical entire
etale completion has a pure boundary whose distinct components are disjoint
and whose normalizations are affine lines. Each component is homeomorphic
to the complex affine line. Consequently the usual Euler and free-deck
divisibility formulas require no smooth-boundary or common-pencil input.
The argument does **not** prove that these boundary curves are smooth, that
their normalization maps are algebraic isomorphisms, or that a common source
coordinate extends across them.

## 1. Setup, scope and direct conclusions

Let
```
K=(P,Q): X=Spec R=A2_C -> Y=Spec B=A2_C,
R=C[x,y], B=C[P,Q], F=C(P,Q), L=C(x,y),
c=Jac(P,Q) in C*.
```
Let S be the full normalization of B in L, Z=Spec S, and let U be the
ENTIRE etale locus of Z->Y. Every omitted prime is retained in Z; U is
not obtained by selecting a convenient subcollection. Write
```
A=B[Lambda] subset R,
Lambda={q in L : Tr_(L/F)(qS) subset B},
D=(U minus X)_red,   b0=# irreducible components of D.
```
The canonical trace-dual reconstruction KU0 and KU1, Sections1 and3 of
[the retained proof](CANONICAL_TRACE_DUAL_ETALE_RECONSTRUCTION.md), identify
U=Spec A and make U->Y etale. Those two sections, rather than KU2 or its
differential IC input KQ0, are the algebraic construction used here.

The conclusions proved below are:

* U is a smooth connected affine surface, and the given X is an open subset.
* D is empty or is a pure reduced curve. Its irreducible components are
  pairwise disjoint; each normalization is A1 and is a homeomorphism.
* Each added unramified curve is its entire closure in the full affine Z.
  Every component open M_i=X union D_i is a smooth affine open of Z.
* pi1_top(U)=1, H^2(U,Z)=Z^b0 and chi(U)=1+b0.
* A*=C*, Pic(U)=Cl(U)=Z^b0, with the actual boundary divisors as free basis.
* Omega=c^-1 dP wedge dQ is a regular exact nowhere-zero two-form on U,
  restricting to dx wedge dy on X. The pulled-back target coframe is free.
* Every generic deck subgroup G subset Aut_F(L) acts regularly and freely
  on U. Its order divides both [L:F] and 1+b0.

In particular, an actual nontrivial generic deck group requires b0>0.
This is not an inverse theorem: ramified omitted primes can remain in Z
outside U, and trivial deck group does not give generic degree one.

The boundary and topology assertions in Sections2--4 hold more generally
for **any** smooth complex affine surface containing an open A2. Their
proof does not use an etale target map, a canonical volume form, a source
pencil, whole boundary charts, or smoothness of the boundary.

## 2. The actual source is open, and the boundary is pure

The source map K is etale, hence quasi-finite. Its algebraic normalization
factor X->Z is an open immersion by the affine version of Zariski's Main
Theorem, [Stacks03GT](https://stacks.math.columbia.edu/tag/03GT).
Here normalization in R agrees with the full normalization in L:
an element of L integral over B is integral over R and hence belongs to R,
since R is normal. Every source point belongs to the etale locus of Z->Y
on this open, so X subset U. KU0/1 makes U smooth affine and connected.

We give a purity argument that needs only normal affineness of U and
affineness of the given open X. If p were an isolated irreducible component
of D, then p lies on no positive-dimensional boundary component. Because U
is affine, choose f in A with f(p)=1 and f=0 on all the other boundary
components. Then W=D_U(f) is a normal affine surface and
```
W intersect X = W minus {p} = D_X(f|_X)
```
is affine. Normal Hartogs extension across the codimension-two point p
identifies O(W) with O(W minus {p}); equivalently both rings are the
intersection of the same height-one local rings in L. An open immersion
between these two affine schemes with the identical map of global rings
would be an isomorphism, contradicting the presence of p. Thus no isolated
boundary component exists. Since D is proper closed in an integral
surface, it is empty or a pure curve.

Purity alone does not control the omitted curves at points outside U.
The normalization argument in Section4 will in fact show that the full
closure in Z of every added unramified prime already lies in U. No such
claim is made for the ramified omitted primes, which stay outside U.

## 3. Simple connectivity and the vanishing used for the boundary

A loop in U based in X can be moved, with its basepoint fixed, off D.
Indeed choose a finite smooth stratification of the algebraic curve D;
its strata have real dimension at most2 inside the real four-manifold U.
General position for a one-dimensional loop avoids each stratum, since
1+2<4. Perturbation relative to the basepoint is possible because that
point has a neighbourhood disjoint from D. Consequently
```
pi1(X) -> pi1(U)
```
is surjective. X(C)=C^2 is simply connected, so pi1(U)=1. This argument
does not suppose that D is nonsingular.

The other topological input is the affine Morse theorem:
a smooth complex affine surface has the homotopy type of a CW complex
of real dimension at most2. Thus
```
H^j(U,Z)=0 for j>2.                                      (1)
```
The primary statement used is Milnor, *Morse Theory*, Theorem7.2,
printed p39, with proof on pp39--41; the exact cached source and adopted
foundations are recorded in Section8. In brief, a generic squared-distance
function on a closed complex embedding is a proper Morse exhaustion;
the complex Hessian pairs its relevant eigenvalues with opposite signs,
bounding every index by the complex dimension. Morse theory then attaches
only cells of those indices. We adopt this classical Morse theorem,
rather than claiming to reprove its existence or handle-attachment inputs.

All groups below have integral coefficients. Use Borel--Moore homology
H_i^BM, with its closed/open exact sequence and Poincare duality on
oriented manifolds:
```
H_i^BM(U) = H^(4-i)(U).                                  (2)
```
These standard topological foundations apply to the triangulable
algebraic sets in question. Since X=C^2 has BM homology Z in degree4
and zero in every other degree, the closed/open sequence for D subset U
gives isomorphisms
```
H_i^BM(D) -> H_i^BM(U)   for i=1,2.                      (3)
```
In particular (1)--(3) imply
```
H_1^BM(D)=H^3(U)=0.                                      (4)
```
The indices in (2) are ordinary cohomology against BM homology; no
compact-support cohomology substitution is made.

## 4. Curve normalization forces disjoint topological affine lines

Let nu:Dtilde->D be the finite normalization and let T be a finite set
containing all singular points of D. Set Ttilde=nu^-1(T). The complement
Dtilde minus Ttilde maps isomorphically to D minus T.

Each connected component of Dtilde is a smooth affine curve, because
normalization is finite. It is noncompact: a complete positive-dimensional
curve cannot be a closed subvariety of an affine variety. Thus
H_0^BM(Dtilde)=0; likewise H_0^BM(D)=0. The closed/open sequences for
these two finite subsets have the form
```
0 -> H_1^BM(Dtilde) -> H_1^BM(D minus T) -> Z^Ttilde -> 0,
0 -> H_1^BM(D)      -> H_1^BM(D minus T) -> Z^T      -> 0.
```
The right vertical map is the summation map over each fibre nu^-1(p).
Its kernel is the direct sum of groups Z^(r_p-1), where
r_p=#nu^-1(p), the number of local analytic branches. Comparing the
sequences gives the precise normalization sequence
```
0 -> H_1^BM(Dtilde) -> H_1^BM(D)
  -> direct-sum_(p in T) Z^(r_p-1) -> 0.                 (5)
```
This explains the finite-normalization step directly; no equality between
a residue-field degree and a branch count is assumed.

By (4), both nonnegative contributions in (5) vanish. Hence r_p=1 at
every point: D is locally unibranch. Distinct irreducible components cannot
meet, since an intersection would have at least two normalization
preimages. Now let C be one normalization component, let g be the genus
of its smooth projective completion, and let s>=1 be its number of deleted
points. For this noncompact oriented real surface, duality gives
```
rank H_1^BM(C)=rank H^1(C)=2g+s-1.
```
It vanishes by (5), so g=0 and s=1. Over C this is exactly C=A1.

The finite normalization map of each original component is now bijective,
since all its fibres have one point. A finite morphism of complex curves
is proper analytically, so this continuous bijection is a homeomorphism.
Therefore
```
D = disjoint-union D_i,   D_i topologically C,
normalization(D_i)=A1.                                  (6)
```
No rationality theorem for a resolution or for target nonproperness curves
was used.

There is also complete coverage of each added unramified prime in the
**original full affine normalization Z**. Let barD_i be its closure in Z.
It lies in Z minus X, because that complement is closed and contains its
generic point; hence D_i=barD_i intersect U. Normalize barD_i to the smooth
affine curve C_i. Normalization commutes with open restriction, so the
inverse image of D_i is its normalization A1. A smooth affine curve
containing an open A1 must equal that A1: its smooth projective completion
is P1, and its nonempty set of deleted points is contained in the single
point deleted for the given A1. Therefore C_i has no point above
barD_i minus D_i. Surjectivity of normalization proves
```
barD_i=D_i subset U.                                     (6a)
```
Thus an unramified added prime cannot have exceptional nonetale points
of its closure missing from U. Its possible intrinsic cusp singularities
remain a separate issue.

For each i, the actual component open
```
M_i=X union D_i = U minus union_(j!=i) D_j               (6b)
```
is smooth and affine and is an open subset of the full Z containing the
WHOLE D_i. To see affineness, each D_j is an effective Cartier divisor
on smooth U. Their finite sum is again an effective Cartier divisor.
Its complement immersion into U is affine: on any affine open cover
trivializing its invertible ideal the preimage is a principal affine
open. Affineness of a morphism is local on the target, and U is affine.
For b0=1 this simply says M_i=U. These component opens therefore need no
whole-root-chart theorem for their affineness or their finite-point
coverage. Their boundary still has only the algebraic normalization A1
proved in(6), rather than a proved smooth A1 structure.

This argument allows an algebraic cusp. For example the normalization
t |-> (t^2,t^3) of the affine cusp is a finite bijective homeomorphism from
A1, but its derivative vanishes at0 and the cusp is not a smooth curve.
This is a countercontrol to the implication
"normalization A1 and topological A1 imply algebraic A1"; it is **not**
claimed to be a boundary example for the full Keller or affine-pair setup.

Equation(6) gives H_2^BM(D)=Z^b0. Equations(2)--(3) yield
H^2(U,Z)=Z^b0. Euler additivity gives
```
chi_c(U)=chi_c(X)+chi_c(D)=1+b0.
```
A smooth complex surface has chi_c(U)=chi(U) by oriented real
four-dimensional duality and finite Betti numbers. Thus
```
chi(U)=1+b0.                                             (7)
```
Equivalently, the vanishing above degree2 together with pi1(U)=1 and
H^2(U,Z)=Z^b0 gives the same ordinary Euler count.

## 5. Units, divisor classes and the canonical exact coframe

Restriction embeds A* into R*=C*, so A*=C*. A divisor on smooth U
restricts to a principal divisor on the factorial source X. Subtracting
that rational principal divisor leaves a divisor supported on D.
Conversely, if sum n_i D_i is principal, its defining rational function
restricts to a unit of R and hence is a nonzero constant. All n_i vanish.
Therefore
```
Pic(U)=Cl(U)=direct-sum_i Z[D_i].                         (8)
```
This proof does not use the smoothness or topology of D_i.

Because U->Y is etale, dP,dQ are a global basis of Omega^1_U and
```
Omega=c^-1 dP wedge dQ=d(c^-1 P dQ)                      (9)
```
is regular, nowhere zero and exact. It restricts to dx wedge dy on X,
so K_U is trivial. The target derivations also extend uniquely to
commuting regular vector fields on U. On the source they are
```
D_P=c^-1(Q_y partial_x-Q_x partial_y),
D_Q=c^-1(-P_y partial_x+P_x partial_y).
```
Their regularity is not a claim of local nilpotence or complete flow.
In particular it does not extend a chosen source translation partial_x
across the boundary without further input.

For completeness, each actual source affine line is closed in U.
Its restriction algebra A' subset C[t] contains a nonconstant polynomial
a(t), since the map to affine U is not constant. C[t] is finite over
C[a(t)] and therefore finite over the intermediate A'. Thus the line
maps finitely onto its closed image in U. Its map to U is already a
locally closed immersion, so it is a closed immersion. The same argument
applies to any polynomial A1-parametrized source curve as a closed image;
it does not make arbitrary affine source curves proper in U.

## 6. Generic deck transformations act regularly and freely

Let sigma in Aut_F(L). Integrality over B is intrinsic, so sigma(S)=S.
It induces an automorphism of the full Z over Y, and preserves the entire
etale locus U. Thus the rational generic deck action is regular on U,
without any prior assertion about boundary charts or its action on X.

For an automorphism sigma of U over Y, its fixed-point subscheme is the
inverse image of the diagonal under
```
(id,sigma):U -> U times_Y U.
```
The diagonal is open because U->Y is etale, and closed because it is
separated. Hence the fixed locus is open and closed in connected U.
For a nonidentity generic automorphism it cannot be all of U, and is
therefore empty. This proves freeness at all points, including every
boundary point, rather than only generically. The form Omega in (9)
is deck invariant because P,Q are fixed.

Every subgroup G is finite, since L/F is finite separable. Its free action
has an affine smooth quotient U/G and gives a finite etale covering
U->U/G of degree |G|. Equivalently on complex points it is a free finite
covering of finite-type spaces. Euler multiplicativity and (7) imply
```
|G| divides 1+b0.                                        (10)
```
Artin's fixed-field theorem also gives
[L:F]=|G|[L^G:F], so
```
|G| divides gcd([L:F],1+b0).                              (11)
```
No assertion that G preserves the original open X or permutes the original
D_i is needed or made.

## 7. What remains of the boundary input

The complete canonical construction and Sections2--6 supply smooth
affineness of U, simple connectivity, a pure disjoint topological-A1
boundary, whole closure coverage of every added unramified prime, smooth
affine component opens X union D_i, the exact Euler count, the boundary
divisor basis, the regular free generic deck action and its Euler
divisibility.

They do not establish any of the following:

* smoothness of the reduced curves D_i, or algebraic isomorphisms D_i=A1;
* whole-component polynomial root charts, including all finite points;
* a common **source coordinate** y with nonnegative value at every
  omitted U-divisor, or extension of its source translation;
* triviality of the deck group, trace-dual multiplication regularity,
  or inversion of the actual Keller map.

The first two are stronger algebraic geometry than(6), even after the
complete finite-point coverage in(6a)--(6b). A common-coordinate or Hamiltonian extension argument must
prove its own hypotheses; the exact coframe alone does not supply them.
Thus this note narrows the LA0 input and independently validates the
Euler/free-deck consequences. It does not adopt the unproduced full LA0,
PF46, PF50, a common-pencil assertion or the differential input KQ0.

## 8. Source audit and verification scope

The essential retained construction is KU0/1, specifically Sections1 and3
of CANONICAL_TRACE_DUAL_ETALE_RECONSTRUCTION.md. The algebraic Zariski
main statement [Stacks03GT](https://stacks.math.columbia.edu/tag/03GT)
was inspected with its displayed proof; the affine, finite-type,
quasi-finite source hypotheses above match it.

The primary affine Morse source is J. Milnor, *Morse Theory*, Annals of
Mathematics Studies51, 1963; cached fifth printing1973 of the1969 corrected text:
[sources/milnor-morse-theory-1969.pdf](sources/milnor-morse-theory-1969.pdf),
downloaded from the
[IMPA-hosted primary book](https://luis.impa.br/aulas/topdif/Milnor_MorseTheory.pdf).
Theorem7.2 and its proof were read on printed pp39--41, PDF zero-index
pages47--49. The theorem concerns closed complex analytic manifold
embeddings, so applies to a closed algebraic embedding of the smooth
affine U. The proof explicitly uses the earlier generic squared-distance
and Morse results, especially Lemma6.9 and the handle-attachment
theory. These are adopted classical foundations, not independently
reproved here. The companion source-audit JSON pins the exact PDF and
records these locators.

Normal Hartogs extension, algebraic curve normalization, general position,
Borel--Moore closed/open exactness, proper pushforward and manifold
duality, and finite covering Euler multiplicativity are standard
foundational inputs. Their actual applications, especially the two
sequences leading to(5), are written above. No new surface-classification
theorem is imported.

The exact checker tests the closed/open degree ledger, the normalization
branch-summation maps, genus/puncture formulas, Picard-basis relations,
the cusp sharpness control, maximal affine-A1 curve coverage, Cartier
component-open affine charts, deck Euler divisibility and source/coframe
identities. Its finite controls do not prove the unrestricted statements;
the written argument does. Source/checker/input hashes are recorded in
the receipt. No external review, formal verification or independent
replay is claimed by creating that receipt.
