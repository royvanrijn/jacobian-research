# Symmetry groups of polynomial Keller maps: affine inertia and realization

14 November 2026. This is a synthesis with an additional affine obstruction.
It does not claim a complete classification. The requested classification
remains OP-KMON in [MATH_STATUS.json](../../MATH_STATUS.json). Statements
labelled KMFI1--3 are proved below; the finite replay checks consequences,
not the universal geometric proof. No novelty-priority claim is made.

## 1. The object being classified

Work first over an algebraically closed field k of characteristic zero.
A Keller map is a polynomial morphism F:X=A^m -> Y=A^m with
det DF in k*. It is dominant, generically finite and etale. Write

\[
 K=k(Y),\quad L=k(X),\quad M=\text{Galois closure of }L/K,\qquad
 G=\operatorname{Gal}(M/K),\quad H=\operatorname{Gal}(M/L).
\]

The degree is n=[L:K]=[G:H]. Analytic continuation over C, or the geometric
etale fundamental group in general, acts on the cosets G/H. This is the
monodromy action. It is faithful because M is the Galois closure.

The action must be recorded together with H. The deck group of the generic
cover is N_G(H)/H, equivalently the centralizer of G in S_n. It need not
equal G. Polynomial source symmetries covering target symmetries form yet
another object. These three meanings of "symmetry group" must be kept separate.

For a nonclosed ground field, geometric and arithmetic monodromy are
different groups, G_geom normal in G_arith <= S_n. Their quotient describes
the algebraic constant field of the Galois closure. The F11 example below
makes this distinction necessary, even with an explicit rational formula.

The absolute spectrum consists of actions realized by polynomial Keller
self-maps of affine space, allowing dimension to vary. The chart spectrum
allows smooth affine source and target varieties and a chosen volume form.
A field-versal spectrum additionally asks to realize all torsors of a
specified action. These are different problems.

## 2. Three distinct realization steps

A finite etale algebra in a special fiber does not prescribe the generic
monodromy group. A special target slice can have a smaller group than the
ambient cover. The universal pencil H(W)-sW+t has group S_n for every
characteristic-zero degree-n polynomial H; choosing a Chebyshev or a
monomial seed alone does not preserve its dihedral or cyclic slice group.
See the retained [universal-pencil proof](../../archive/non-elliptic/verified/UNIVERSAL_SYMMETRIC_MONODROMY.md).

There is, however, a general cover-preserving chart construction. If
f:V -> U is a finite etale map of smooth affine curves with coordinates x,t
on open subsets of the projective line, write t=f(x). On a chart where x
is a coordinate and f' is a unit,

\[
 (x,z)\longmapsto (f(x),z/f'(x))
 \tag{2.1}
\]

has determinant one and preserves the function-field extension after adding
one independent transcendental variable. A monic root presentation
P(t,x)=0 gives the analogous derivative-unit construction.

Such charts can have many nonconstant units. Since A^m has only constant
units, adjoining affine coordinates does not make these charts affine
spaces. Changing the boundary through affine modifications is a separate
problem. A candidate absolute realization must exhibit polynomial
coordinates everywhere, a nonzero constant full Jacobian, an exact generic
field identification, and actual affine source and target spaces.

## 3. Every affine divisor is met without ramification

**Lemma.** Let F:A^m -> A^m be Keller over k. For every prime divisor
D=V(q) in the target there is a prime divisor E in the source dominating D,
and its ramification index in L/K is one.

Indeed q(F) is nonconstant, since F is dominant. A nonconstant polynomial
over an algebraically closed field is not a unit and has a prime factor.
Every component E of V(q(F)) has dimension m-1. The restriction of a
quasi-finite morphism is quasi-finite, so its image has dimension m-1;
its image closure is D. At the generic points the local rings are DVRs.
The local map is etale, hence has ramification index one and separable
residue extension. This works for m=1 as well.

In fact q(F) is squarefree: every component just considered has index one.
Distinct target prime divisors have disjoint sets of source prime-divisor
factors, since such a factor cannot dominate two different divisors.
Consequently pullback on the free divisor groups is a primitive injection:
each target prime maps to a nonempty sum of distinct source primes, and
the supports of these sums are disjoint. Choosing one source prime over
each target prime gives a splitting. The cokernel is free.

One elementary consequence, useful for quadratic resolvents, is

\[
 b(F)\in k(X)^{*d}\quad\Longrightarrow\quad
 b\in k(Y)^{*d}\qquad(d\ge2).
 \tag{3.1}
\]

Write the divisor of b in target prime divisors and compare valuations on
one source prime over each. Every exponent must be divisible by d; the
remaining constant is a d-th power in algebraically closed k. This
argument needs neither a prescribed primitive polynomial nor a Jacobian
ansatz. It cannot be applied to a selected target slice instead of the
whole affine target.

## 4. The fixed-sheet inertia theorem

**KMFI1.** For every prime divisor D of the affine target, a geometric
inertia group I_D of M/K fixes at least one sheet of G/H. The normal
subgroup generated by all these inertia groups is G. Consequently

\[
 \boxed{\langle gHg^{-1}:g\in G\rangle=G.}
 \tag{4.1}
\]

Here "geometric inertia" means inertia after strict henselization at D;
it is cyclic in residue characteristic zero. The group is determined up
to conjugacy. The theorem concerns divisors of A^m, not arbitrary
divisors at infinity on a compactification.

To prove the first assertion, use E from Section 3. In the standard
inertia action on embeddings of L into M, the orbits have sizes equal to
ramification indices of geometric branches. Index one gives a fixed
sheet. Equivalently, for a suitable conjugate of H one has I_D <= gHg^-1.
An unramified branch may have a nontrivial residue-field extension before
strict henselization; this does not remove its inertia fixed points.

For the second assertion, let N be the normal subgroup generated by all
the I_D, and normalize Y in M^N. This normalization is finite, since Y
is excellent and M^N/K is finite. All codimension-one inertia is killed,
so the finite morphism is unramified in codimension one. Its source is
normal and its target regular. Purity of the branch locus makes it
finite etale everywhere. Affine space over an algebraically closed
characteristic-zero field has no nontrivial connected finite etale cover.
Thus M^N=K and N=G. Since every I_D lies in a conjugate of H, (4.1) follows.

The purity input is [Stacks, Lemma 58.21.4](https://stacks.math.columbia.edu/tag/0BMB).
For simple connectedness, descend the finite etale cover to a finitely
generated characteristic-zero field, embed it in C and use comparison
with the contractible analytic space C^m. Equivalence between algebraic
finite etale covers and finite analytic covers is the Riemann-existence
comparison of SGA 1, Expose XII, Theorem 5.1 and Corollary 5.2.
Isomorphism and geometric connectedness descend/extend with algebraically
closed constants. These are the only global inputs in this proof.

An equivalent formulation of (4.1) is that there is **no nontrivial
intermediate Galois field** K < E <= L with E/K Galois. Such a field
corresponds to a normal subgroup N of G containing H. Its normalization
over Y would be the finite etale cover just excluded. No polynomial
factorization through E needs to be assumed.

This obstruction also applies to an arithmetic Galois intermediate
field over a characteristic-zero k: finite etale covers of A^m_k are
constant covers after geometric triviality, whereas k is relatively
algebraically closed in k(x_1,...,x_m). We use the geometric formulation
for the action classifications below.

The assumptions matter. Kummer covers of G_m can have regular monodromy,
because an affine divisor has been deleted and the target is not simply
connected. In characteristic p, x -> x^p-x is finite etale on A^1 with
regular C_p monodromy. These are controls against dropping the target
or characteristic hypotheses.

## 5. Absolute exclusions beyond regular actions

**KMFI2.** The following are necessary consequences of KMFI1.

1. No nontrivial regular action occurs: H=1 cannot normally generate G.
   This recovers the classical Galois-case obstruction.
2. No nontrivial nilpotent group occurs in any faithful transitive action.
   For a finite nilpotent group every maximal subgroup is normal. A proper
   H lies in such a subgroup, so its normal closure is proper.
3. No proper normal block stabilizer is allowed. Under the
   block/intermediate-field correspondence, a normal block stabilizer
   would give a nontrivial Galois intermediate field. This does not
   exclude all imprimitive actions or algebraize arbitrary intermediate
   fields into polynomial factors.
4. H must surject onto every quotient in which its normal closure would
   otherwise be proper; in particular H -> G_ab is surjective. This last
   test alone is weaker than (4.1), for example for perfect groups.

For the natural degree-n action of
D_n=<r,s | r^n=s^2=1, srs=r^-1>, of order 2n, the point stabilizer is <s>.
Conjugating s by r gives r^2s, so

\[
 \langle\!\langle s\rangle\!\rangle
   =\langle s,r^2\rangle
   =\begin{cases}D_n& n\text{ odd},\\
                 \text{an index-two subgroup}&n\text{ even}.
     \end{cases}
 \tag{5.1}
\]

Thus **every even-degree natural dihedral action is impossible absolutely,
in every affine dimension**, irrespective of masks, blowdowns or
stabilization. For even n both reflection classes fail. In fact these and
the regular action exhaust faithful transitive dihedral actions: the
rotation part of any stabilizer is a normal subgroup contained in its
core, so a core-free stabilizer is trivial or a reflection. Therefore
D_n for even n has no absolute faithful transitive realization at all.

There is an independent elementary check of both exclusions. Every
quadratic subextension has the form K(sqrt(b)); equation (3.1) forbids
such a subextension inside L. More generally, over algebraically closed
k every cyclic prime-degree extension is Kummer, so (3.1) forbids it.
A finite nilpotent group has a normal maximal subgroup of prime index
containing any given proper stabilizer. The even dihedral action has
the explicit quadratic quotient (5.1). Thus these exclusions do not
depend on the global purity proof of KMFI1. On the even Dickson cover
itself, the identity P_(2m)(a,u)+2u^m=P_m(a,u)^2 displays the quadratic
resolvent directly.

Odd natural dihedral actions pass (4.1); this is not a construction.
D_3=S_3 is already realized by the retained symmetric family. D_5 and the
larger odd dihedral actions still require absolute realizations or new
obstructions. The Dickson derivative-unit charts remain valid in both
parities. Their even-degree absolute search is superseded by (5.1).

The same computation gives a useful affine-group test. For an abelian
group V with a faithful complement A acting on it, in the transitive
action of V semidirect A on V,

\[
 \langle\!\langle A\rangle\!\rangle=[V,A]\rtimes A.
 \tag{5.2}
\]

Indeed a normal subgroup containing A contains every commutator [v,a],
and quotienting by [V,A] makes A normal. Hence the necessary condition is
V/[V,A]=0. It excludes generalized dihedral actions when V/2V is nonzero.
For F_20=C_5 semidirect C_4 this quotient vanishes; its natural action
survives this test. These universal group statements have group-theoretic
proofs; the bounded finite replay illustrates and checks them.

## 6. The PSL2(11) local sieve

**KMFI3.** For the natural degree-twelve action of PSL2(F_11), the complete
element-order and cycle dictionary is

| order | number of elements | cycles | allowed nontrivial affine-divisor inertia? |
|---:|---:|---|---|
| 1 | 1 | 1^12 | unramified |
| 2 | 55 | 2^6 | no |
| 3 | 110 | 3^4 | no |
| 5 | 264 | 5^2 1^2 | yes |
| 6 | 110 | 6^2 | no |
| 11 | 120 | 11 1 | yes |

The dictionary follows either by the split/nonsplit/unipotent
classification or by enumerating SL2(F_11)/{+/-I}; the replay does both
the projective action and its fixed-point test. In general the natural
PSL2(q) action permits only split semisimple and unipotent nontrivial
divisor inertia. Nonsplit semisimple inertia is a derangement.

The modular (2,3,11) cover on X_0(11) has two derangement branch values.
A Mobius rechoice of the target point at infinity can move only one of
them out of the affine line. Thus that product presentation cannot
satisfy the affine-divisor test. This is a statement about that cover and
its affine branch placement, not an exclusion of the natural action.

Suppose now that a degree-twelve cover of P^1 with this natural group
has genus zero and every *finite* branch inertia has a fixed sheet.
Finite inertia can have only orders 5 and 11, of indices 8 and 10.
The index at infinity belongs to {0,6,8,10}. Riemann--Hurwitz gives

\[
 8a+10b+\operatorname{ind}(\sigma_\infty)=22.
 \tag{6.1}
\]

For nonnegative integers a,b the unique solution is a=2, b=0,
ind(sigma_infty)=6. Thus the unique possible genus-zero passport with
this branch-placement condition is

\[
 5^2 1^2\mid5^2 1^2\mid2^6.
 \tag{6.2}
\]

This is a necessary passport classification under the genus-zero
hypothesis, not a uniqueness theorem for dessins or a genus-zero theorem
for arbitrary inverse curves of Keller maps.

The passport really occurs. In PSL2(F_11), matrices modulo scalar sign

\[
 A=\begin{pmatrix}0&1\\10&3\end{pmatrix},\quad
 B=\begin{pmatrix}0&1\\10&8\end{pmatrix},\quad
 C=AB=\begin{pmatrix}1&3\\3&10\end{pmatrix}
 \tag{6.3}
\]

have orders 5,5,2, satisfy ABC=1 projectively and generate all 660
elements. Their natural cycle indices are 8,8,6, giving genus zero.
The checker certifies generation rather than inferring it from a
passport. Riemann existence gives a rational cover over C.

A known explicit formula is van Hoeij--Vidunas's F11. The
[full F11 audit](F11.md) records its exact coefficient field, the target
reparametrization, rational derivative, unit lattice and an arithmetic
discriminant test. The affine follow-up proves that this rational curve
input cannot be completed absolutely by independent stabilization and
birational recharts.

## 7. A global common-index obstruction for curve-cover fields

The [affine-completion proof](F11_AFFINE_COMPLETION.md) adds KMFI5--8.
A multi-pole rational function cannot become polynomial after composition
with a nonconstant rational function on affine space. F11 nevertheless
has an explicit polynomial cone model of degree twelve, and a polynomial
threefold model whose Jacobian is the pullback of U(U-V)V. Those models
are not Keller.

The stronger field obstruction is independent of the chosen affine
target boundary. For a rational curve cover of degree n, let e be the
gcd of the ramification indices in any one fiber. If its field pair,
after adjoining the same independent variables on both sides, admits
any absolute Keller model in birational source and target coordinates,
then

\[
 \gcd(e,n/e)=1 .
\]

The proof compares a forced d-th power, d=n/e, with valuations at the
original curve fiber, even when those valuations have moved to the
boundary of the new affine target. Local fixed-sheet inertia alone
does not detect this condition.

For natural degree-twelve PSL2(11), every genus-zero curve cover has
an involution branch with fiber 2^6. Its e=2 and n/e=6 violate the
condition. Positive-genus curve sources cannot become rational after
independent stabilization. Hence no curve cover in this natural
action can yield an absolute Keller map through independent affine
stabilization and birational recharts. F11 and an explicit rational
target replacement that passes every affine-divisor inertia test are
both excluded by this global obstruction.

This does not exclude the action for a genuinely higher-dimensional
field pair that is not induced from a curve in this way. The two
exceptional degree-eleven actions survive this common-index test.

## 8. Global denominator descent

[KMFI9--11](CURVE_DESCENT_LCM.md) strengthen the curve-generated test.
For every branch fiber let e_b be the gcd of its complete list of
ramification indices. A curve field pair that admits an absolute Keller
model after independent stabilization and arbitrary birational affine
recharts must satisfy

\[
 \operatorname{lcm}_b e_b=n.
\]

The proof uses a denominator class modulo n-th powers and degree zero
of its divisor on the original target curve, rather than only one local
valuation. The degree-eight PSL2(7) control has every fiber gcd equal
to one: it passes the earlier local and single-fiber tests but fails
this global condition.

For an action define C(G,Omega) as the lcm, over all its permutations,
of the gcd of their cycle lengths. Curve-generated Keller realization
requires C=n. Natural PSL2(q) for every odd q has C=(q+1)/2, so all
such curve inputs are excluded. The exceptional degree-seven and
degree-eleven actions pass this group test. In prime degree, the
global condition forces an actual totally ramified fiber; for
degree-eleven PSL2(11), the index budget then forces the retained
(2,3,11) polynomial passport.

The filter has an actual positive control: the cubic weighted Keller
map is explicitly identified with the curve tau=-(v+1)/v^3 after
adjoining two independent parameters. Both field reconstructions and
the global constant Jacobian are verified. Thus curve-generated
Keller maps do exist. These theorems do not assert that every
higher-dimensional Keller field pair descends to a curve.

## 9. Subgroup reduction and the rationality of the target

[KMFI17--18](SUBGROUP_REDUCTION.md) test the tempting step of taking
a PSL2(11) subgroup of a realized symmetric or alternating group.
Galois-closure base change really reduces the group when the subgroup
acts transitively, but its target field must still be rational to
support an affine-space map.

For the known independently stabilized degree-eleven A11 and S11
inputs, the PSL2(11) quotient targets have genera 5678, 7331 and
9343. All three therefore remain nonrational after any independent
stabilization; no birational affine rechart can complete those
specific field pairs to an affine-space map. The smaller group
has been realized as a field action, with no affine-space target.

For the two A12 inputs, the analogous generic target-pencil genera
are 527826 and 524161. Those genera alone do not decide rationality:
rational surfaces can have positive-genus generic fibers.
[KMFI19](SURFACE_FORMS.md) supplies the missing global obstruction.
Explicit PSL2(11)-invariant two-forms extend across every finite
branch collision and the weighted projective boundary, and exact
nonzero certificates show they do not vanish. They descend to
smooth projective target models. Both targets are therefore
not stably rational, so all five audited subgroup reductions
are excluded even after independent stabilization and arbitrary
birational coordinate changes. Independent Klein Shabat completion
and the exceptional degree-eleven PSL2(11) occurrence question remains open;
the natural degree-twelve action is excluded absolutely by KMFI67 below.

[KMFI20--22](ORBIT_FORMS.md) give a general representation criterion
for the two A12 normal closures. If r_m is the subgroup's number of
orbits on m-subsets of inverse labels, the target geometric-genus
bounds are r_6-r_1 for the cone and r_6-r_2 for the shear.
The complete group argument then excludes every proper transitive
subgroup reduction of either family, including M11 and M12.
These are stable-rationality obstructions for the specified fields.
The same form modules obstruct source rationality of many ordered
and unordered root-configuration covers, even while their generic
group stays A12. Thus changing the point stabilizer inside an
actual normal closure does not automatically give another
affine-space realization of its new permutation action.

## 10. Positive constructions and the remaining classification

The retained weighted examples give absolute natural S_n actions for
all n>=3; primitive monodromy makes such maps atomic.
[KMFI14--16](ALL_ALTERNATING.md) prove actual affine-three-space
realization of every natural A_n action, n>=5. Higher-power radicands
supply every odd degree, and the level-two cyclic cones of
Kistner--Shaska supply every even degree at least six. Global
polynomiality, constant Jacobian, degree and geometric group all
have written proofs and exact controls. The odd examples have
explicit independently stabilized curve presentations. In contrast,
every A_(2^a) example, a>=3, is proved to lie outside that class,
even after further independent stabilization and birational changes.
Thus the curve-generated spectrum is strictly smaller than the full
affine-space spectrum, and C=n cannot be imposed on arbitrary
higher-dimensional field pairs. The
[earlier alternating audit](ALTERNATING_REALIZATIONS.md) retains
the attributed quadratic beta family and the alternating-shear A12 map.
Products of maps
in independent variable blocks give product actions: the covers and
their Galois closures are linearly disjoint over independent rational
target coordinates. Identity stabilization preserves the action.
Composition gives a transitive subgroup of an appropriate wreath
product and an intermediate-field block system. Equality with the full
wreath product requires a separate disjointness proof.

For reference, the [atomicity proof](../../archive/non-elliptic/verified/PRIMITIVE_MONODROMY_ATOMICITY.md)
uses the block/intermediate-field correspondence and the fact that
polynomial factors of a Keller map are themselves Keller.
The [Ritt product theorem](../../archive/non-elliptic/verified/KELLER_RITT_PRODUCT_THEOREM.md)
provides a more precise statement for the retained independent-product
family. No such factorization theorem follows for a general action
without affine-space models of the intermediate fields. KMFI30 now gives
an explicit infinite family where the converse fails: the maps are atomic
and their monodromy has a unique proper block system.

| action | absolute conclusion | remaining issue |
|---|---|---|
| nontrivial regular | excluded | none for this action |
| any nontrivial nilpotent G | excluded | none for these actions |
| natural D_n, n even | excluded | chart constructions remain valid |
| natural S_n, n>=3 | realized in retained constructions | classification does not follow |
| natural A_n, n>=5 | realized on A^3 by KMFI14--15 | other transitive actions and the full classification remain separate |
| signed 2^(n-1):A_n on 2n letters, even n>=6 | realized on A^4 by KMFI29; atomic but imprimitive by KMFI30 | other signed and imprimitive actions remain separate |
| natural A_(2^a), a>=3 | actual Keller examples provably outside independent curve stabilization by KMFI16 | C=n is a curve-class condition |
| natural A_3 | excluded as nontrivial nilpotent | none for this action |
| natural D_n, n odd>=5 | survives these tests | absolute map or further obstruction |
| every faithful transitive A_4 action | excluded absolutely by KMFI66 | none for this group |
| natural F_20 | survives the present action tests | actual map or further obstruction |
| exceptional degree-11 PSL2(11) actions | curve-induced compatible inputs must have the polynomial (2,3,11) passport | absolute completion remains open |
| natural PSL2(q), q=3 mod4 | excluded absolutely by KMFI67, including degree12 at q11 | none for these actions |\n| natural odd-q PSL2(q), q=1 mod4 | independently curve-induced inputs excluded by KMFI9--10 | other-field occurrence remains open |

For the exceptional degree-eleven actions an order-six element has
cycles 6,3,2 and no fixed sheet, and an order-eleven element is an
eleven-cycle. Both are forbidden at affine divisors. The Shabat
(3,2,11) placement keeps the eleven-cycle at infinity and therefore
passes this local test. The existing Gassmann, conductor and boundary
calculations remain relevant; see the
[retained benchmark](../../archive/non-elliptic/extended-geometry/PSL2_11_KELLER_ACTION_SPECTRUM.md).

The [fiber-incidence theorem](FIBER_INCIDENCE.md), KMFI23, gives a further
completion constraint for independently stabilized polynomial covers:
the graph joining source-pencil fibers to target prime components
must be connected, with ramification index times source multiplicity
equal to target multiplicity on every edge. In the exceptional
degree-eleven Shabat case, two distinct unramified tail fibers must
be reducible and nonreduced, one mixing multiplicities divisible
and not divisible by three, and one mixing even and odd multiplicities.
The two finite branch fibers need at least fourteen distinct source
prime components in total. An actual S4 map realizes the minimal
four-edge incidence tree for its double branch, giving a positive
control. These conditions exclude simplified completion pencils,
not every possible affine completion.

The [normalized sweep analysis](SWEEP_MONODROMY.md), KMFI24--25,
audits another construction route. Every normalized ruled-surface sweep
of nontrivial degree N has full S_N monodromy, by a birational
coefficient-curve, double-transitivity and simple-tangency proof.
This includes all free ruled members of Gao's four-dimensional
construction. In the different degree-ten branch, Gao F5 is an
actual rational-coefficient affine four-space Keller map with geometric A10,
arithmetic S10 and constant field Q(sqrt(6)).
The geometric sign follows from the complete finite sweep and its
square Jacobian, while an exact prime-seven Frobenius proves the
large arithmetic group. These do not produce a PSL2(11) action.

The [nonlinear sweep frontier](SWEEP_FRONTIER.md), KMFI26--28,
extends this construction analysis. Every normalized four-dimensional
curve-type sweep has full symmetric monodromy, even with arbitrary
nonlinear transverse potential. An explicit nonlinear affine
four-space S4 Keller map supplies a control. Finite normalized sweeps
have uniform-cycle inertia (k+1)^delta with the fixed-sheet condition;
the natural PSL2(11) action first permits such a Gao sweep in
dimension six. Finite tautological sweeps with a polynomial multiplier
depending on one parameter have branch-surface degree at most two:
the degree-one case is alternating and the degree-two case imprimitive.
Thus neither specified PSL2(11) action comes from a primitive sweep
in this whole multiplier class. The formerly marked degree-twelve
imprimitive candidate is now completed by KMFI29 below.

The [signed-alternating family](SIGNED_ALTERNATING.md), KMFI29--31,
adds actual affine-four-space realizations of 2^(n-1):A_n on 2n signed
letters for every even n>=6, with full polynomial stages, constant
Jacobians and exact geometric groups. At n=6 the group has order
11520 in degree twelve; over Q(sqrt(-15)) its arithmetic group is
C2 wr S6 and constants are Q(sqrt(-15))(i,sqrt(5)).
The uniform arithmetic description is given in the family proof.
Every map in this family is absolutely and stably atomic despite its imprimitive monodromy:
its unique proper block field would require a degree-two Keller factor.
Its generic deck involution does not extend to the ordinary affine source.
Their cycle-gcd invariant is n for n a power of two and 2n otherwise.
Thus the power-two subfamily also lies outside every independent
curve stabilization. These are realizations, not just admissible actions.

The [asymptotic-direction theorem](SWEEP_ASYMPTOTICS.md), KMFI32--34,
now excludes the entire finite tautological four-dimensional route to
the exceptional degree-eleven PSL2(11) action, without a multiplier
assumption: the branch-surface degree is at most two, whereas that
action requires three. The proof uses the surface second fundamental form.
Every polynomial multiplier is independent of the second parameter;
a rational multiplier in k(w)[z] has at most one simple pole in its
nonconstant z coefficient. The complete homogeneous rational class
reduces to independently stabilized polynomial covers with group A_N
for odd N and S_N for even N. Its centered presentation cannot supply
the marked origin required by its positive-weight affine stage.
These are construction-class obstructions, not absolute exclusions.

The [rational-sweep classification](RATIONAL_SWEEPS.md), KMFI35--37,
settles all z-dependent rational multipliers in k(w)[z] whose
normalization factor lambda is independent of z. After single-pole
normalization, lambda=w^r and rho=(kappa*z+R(w))/w with R Laurent.
Affine R reduces to the homogeneous class; every non-affine R has
full natural alternating group A_(2(r+j)), where j is its highest
non-affine exponent. The proof checks the missing normalization
branches rather than assuming the sweep finite. It also determines
the exact arithmetic sign class and constant field.
For lambda=w^2, R=alpha*w^D+w^2, every D>=3 has a complete
ordinary-affine-four-space lift with constant Jacobian and degree2D+4.
These are additional actual alternating maps; their actions were
already known to occur by KMFI15. Positive z degree in lambda
and other direction fields still need separate analysis.

The [tautological boundary theorem](TAUTOLOGICAL_BOUNDARY.md),
KMFI38--40, removes finiteness from the local present-cycle gate.
Rank at most one in a tautological sweep forces a linear source
parameter and full symmetric monodromy. At rank two, gamma=0
supplies a three-cycle in the full inertia, even if other sheets
are missing. Fixed-sheet inertia therefore excludes the natural
degree-twelve PSL2(11) action from every tautological
four-dimensional Gao lift, regardless of rank or multiplier chart.
For general divisor-preserving normalized sweeps, the natural
PSL2(11) k=4 or10 dimension gates persist without finiteness.
Every genuine w-pole in a k(w)[z] multiplier collapses a source line.
When that multiplier depends on z, the full Gauss-cover argument
proves X birational onto its branch surface for all lambda.
An exceptional degree-eleven proposal in this rational chart
must therefore exhibit exactly two further ramified index-three
normalization branches above the same branch divisor, outside S.
KMFI42 below now rules out that coincidence in this chart.

The [complete z-dependent rational chart](RATIONAL_CHART_CLASSIFICATION.md),
KMFI41--43, now settles that entire remaining multiplier class.
The triangular wave recursion is necessary and sufficient, with exact
pole bound r+(s+2)m>=2 for a negative least Laurent exponent.
After removing affine R, its degree is
N=(s+2)(r+(s+1)j), always composite.
A constant-leading resultant proves finiteness on Y2!=0, so every
omitted affine branch lies over the separate plane Y2=0.
Its generic inertia is e^g times fixed sheets,
g=gcd(s+1,r-1),e=(s+1)/g. The branch surface has one isolated
three-cycle and no missing ramification above it.
Normal generation and the cyclic complement of the boundary plane
force the full natural group: S_N exactly when r is even and s is odd,
and A_N otherwise. Affine R retains the homogeneous KMFI34 groups.
Neither PSL2(11) field action occurs anywhere in this z-dependent
rational chart, even after independent stabilization or arbitrary
birational recharts of the same field pair.
KMFI44 below closes the polynomial-multiplier frontier; rational
charts outside k(w)[z], other directions and higher architectures
remain separate.

The [polynomial multiplier theorem](POLYNOMIAL_MULTIPLIER_CLASSIFICATION.md),
KMFI44--46, proves that every nonconstant polynomial multiplier h(w)
gives a finite sweep: its weighted heat-polynomial leading model is finite.
The exact degree and full group follow. The quadratic even exception
has geometric SL2(F3) on eight nonzero vectors, arithmetic GL2(F3)
over Q and constants Q(sqrt(-3)); KMFI60 below now completes this
finite sweep field on ordinary affine4 over Q(sqrt(-3)). Higher even exceptions have the signed
V_n:A_n groups already realized in KMFI29.
Together with KMFI43 and the constant-multiplier rank-one case,
the entire k(w)[z] chart now has exactly natural S_n, natural A_n
(n>=3,n!=4), the binary degree-eight action and the even signed family.
Both PSL2(11) field actions and natural A4 are excluded throughout
this whole chart. The formerly proposed nonfinite polynomial h(w)
route is empty under the stated hypotheses.

The [regular-deck rigidity theorem](GENERAL_DECK_RIGIDITY.md), KMFI47,
is independent of the sweep architecture: every polynomial endomorphism
sigma of a characteristic-zero Keller source with F sigma=F is identity.
Its finite-order generic action would be free by etaleness, but a free
finite quotient would force |H| to divide chi_c(A^m)=1.
This does not force N_G(H)/H to be trivial: every nonidentity generic
deck map must instead have a pole somewhere on the affine source.
Simultaneous source/target symmetries are a separate object.
Neither theorem supplies an absolute PSL2(11) exclusion.

The [projective pole-direction theorem](PROJECTIVE_POLE_RIGIDITY.md),
KMFI48--50, constrains rational multipliers outside that chart.
Write the normalization as A*v*v_perp with primitive polynomial
pair v=(b,c). Every homogeneous solution has primitive pair degree
at most one: a partial-fraction fixed-direction identity and
quadratic residue traces force this bound for every scalar degree.
It is linearly equivalent to constant direction or the known
monomial eigen-direction form.
Applied to the highest tensor, this gives an all-scalar obstruction:
the primitive leading ratio has degree at most one; a constant top
direction has a pure transverse common factor, while a linear one
has only its eigenline factors. A passing leading test is not an
existence certificate.
Every affine primitive pair reduces to the completed k(w)[z] chart
and has natural alternating/symmetric group.
For a primitive quadratic pair with nonzero constant scalar,
the complete coefficient reduction leaves only A5,A6 orA8.
The arbitrary-scalar quadratic case is completed by KMFI54--55 below;
higher-degree pairs with these special leading terms remain unresolved.

The [quadratic scalar obstruction](QUADRATIC_SCALAR_OBSTRUCTION.md),
KMFI51--53, supplies a universal constant-leading-direction obstruction
and exact weight0..20 Bezout certificates.

The [complete quadratic reduction](QUADRATIC_CLASSIFICATION.md),
KMFI54--55, now removes every scalar-degree bound.
The bottom weighted minor satisfies a tridiagonal continuant recurrence.
Every edge factor is strictly positive on the real line, and every
leading principal coefficient is positive by a normalized LDL bound
q=i/(4(i+1))<1/4. An elementary induction therefore gives simple real
roots and strict interlacing in every size. The top minor has no real
root, so the two maximal minors are coprime at every weight.
No scalar can repair a curved diagonal denominator.
Together with KMFI51 this reduces every primitive quadratic wave pair
to the completed k(w)[z] chart, for arbitrary polynomial scalar.
Its geometric field action is natural alternating, natural symmetric
or binary SL2(F3) on8; neither PSL2(11) action nor natural A4 occurs.
The complete unbounded action families and exact degree formulas are
given in the new note. The only prime degree in this quadratic class
is5, with group A5. At scalar degree one the exact list remains
A8,A10,SL2(F3) on8.
Higher primitive-pair degrees, other direction fields and architectures
remain open. These field results supply no new binary Keller realization.

The [constant-scalar cubic classification](CUBIC_WAVE_CLASSIFICATION.md),
KMFI56--58, is the first completed class beyond that affine chart.
The primitive pair b=w^2-4z/3,c=w^3-wz+D has a curved denominator
that no affine conjugacy can put into k(w)[z]. Its sweep is finite
of degree21. Its degree-three Gauss field and the branch-degree
bound force a birational ramification surface, so full affine
inertia is one3-cycle and the group is A21.
For D=1 an exact sextic marked point, the opposite scaling
w=gamma^2u,z=gamma v, and one triangular Taylor correction give
an actual polynomial Keller map on ordinary affine four-space.
Its full Jacobian is the nonzero marked gamma0.
This is an affine completion of the curved input, for an action
already in the known spectrum; arbitrary birational inequivalence
to all earlier fields is not asserted.
The complete coefficient exhaustion with constant scalar leaves
exactly A7,A8,A10,A12,A21 for primitive pair degree3.
Both PSL2(11) actions are excluded in this whole class.
Nonconstant cubic scalars and primitive-pair degree>=4, other
directions and architectures remain open.

The [branch-centered completion criterion](BRANCH_CENTERED_COMPLETION.md),
KMFI59--61, now gives an exact affine-completion theorem for every
rank-two normalized polynomial tautological sweep, including
nonfinite ones. For its branch image d(Y)=0 the polynomial pullback
is gamma^3*h. A nonconstant h supplies an actual unramified source
sheet over a good branch value. Affine centering there makes the
weighted orders of X at least1,2,3, permitting the universal scaling
w=gamma u,z=gamma^2v and one quadratic Taylor correction.
If h is constant, gamma supplies a nontrivial cyclic cubic
intermediate field, violating normal generation by the sheet
stabilizer. Consequently normal generation, nonconstant tail and
ordinary-affine4 completion of the pair after one independent
variable are equivalent in this construction class.
For a finite input this is N>3delta, so every degree>6 input qualifies.
This does not assert that every allowed action has such a presentation.

The binary degree8 field now has an actual ordinary-affine4 witness.
Its branch fiber factors as (q-1)^3(q+3); the regular tail point is
defined over Q(sqrt(-3)). Explicit polynomial coordinates have
degrees5,19,22,33 and constant full Jacobian4sqrt(-3)/3.
The entire expanded-coordinate Jacobian is verified exactly,
independently of the determinant-chain proof. The geometric and
arithmetic groups over Q(sqrt(-3)) are SL2(F3) on8, with no extra
normal-closure constants. The map is absolutely and stably atomic
despite its unique pair block and generic deck C2; the deck
involution has no regular source extension.
Its C=4<8 excludes every independent curve-induced Keller route.
The four-sheet A4 block does not give a polynomial Keller factor.
All field-action types in the completed k(w)[z] chart now have
actual action witnesses except the excluded S2 and A3 cases.
KMFI66--67 below now exclude A4 and natural degree12 PSL2(11)
absolutely. Exceptional degree11 PSL2(11) and the full spectrum remain open.

The [homogeneous signed family and forbidden block fields](HOMOGENEOUS_BLOCK_OBSTRUCTION.md),
KMFI62--64, sharpen the intermediate-field distinction.
For every l>=2 the homogeneous h=w^(2l-2),lambda=w sweep has an
explicit ordinary-affine4 completion from a root of a degree4l-6
marked polynomial. Its group is binary8 for l2 and even signed
V_(2l):A_(2l) on4l letters otherwise. The maps are stably atomic.
Their unique pair-block fields have a two-sided reconstruction as
tau=F_l(u)^l/u with independent Y1,e,C. The branch passports
l^2,3*1^(2l-3),(2l-1)*1 give gcd-lcm l<2l.
Thus these exact intermediate fields have no Keller presentation,
even after common independent stabilization and arbitrary recharts.
Both sides remain rational and their stabilizers normally generate
A_(2l): normal generation and field rationality alone are therefore
insufficient outside KMFI59's actual normalized polynomial class.
For l>=3 the same quotient actions A_(2l) occur through other
Keller fields, so this is a field obstruction, not an action exclusion.

In the l6 block A12 normal closure, natural PSL2(11) on12 and
exceptional PSL2(11) on11+1 give exact fixed target curve genera
74038 and74226. The independent spectator reconstruction makes
these intrinsic non-stably-rational targets. The transitive twelve-
sheet reduction and the exceptional eleven-sheet component are
both excluded, without an absolute PSL2(11) conclusion.

The [central prime-cover obstruction](CENTRAL_LIFT_OBSTRUCTION.md),
KMFI65--68, now gives a new absolute action condition:
if a prime l does not divide the sheet stabilizer order, then
H^2(G,C_l)=0. A hypothetical nonsplit central l-cover has an
unramified Brauer obstruction, since every affine inertia has
order prime to l. Brauer purity and Br(A^d)[l]=0 make it lift.
Independent Kummer twists at the finitely many affine prime
divisors remove every new central inertia. The preimage of H
then gives a nontrivial cyclic l-cover of the affine source
unramified at every source prime, contradicting purity and
simple connectedness. The proof controls all affine divisors
and imposes no unproved condition at infinity.

The nonsplit SL2(F_q) double cover and odd natural stabilizer
exclude natural PSL2(q) for every odd prime power q=3 mod4.
Thus the natural degree12 PSL2(11) action is excluded in every
affine dimension and every construction, not just in a chart.
The q3 case excludes natural A4. Its other faithful transitive
actions already fail KMFI1, so A4 is excluded as an abstract
geometric monodromy group altogether.
The actual binary8 group SL2(F3) is preserved: its mod2 H^2
vanishes in the full normalized cochain calculation.
For exceptional degree11 PSL2(11), the stabilizer has order60
and its lifted order120 subgroup is perfect; the coprime
argument does not decide that case.

Combining the A4 exclusion with the nilpotent obstruction and
the actual S3,S4 maps gives the complete action spectrum through
degree4: identity, none, S3, S4 respectively. In degrees3,4
arithmetic monodromy must equal that full symmetric geometric
group, with no additional normal-closure constants.

The [marked-root wreath family](MARKED_WREATH_ACTIONS.md), KMFI69--72,
adds actual Q-defined maps on ordinary affine n-space, for every n>=3,
with degree n(n-2) and determinant(n-2)/2. Harish's construction and
Fonollosa's geometric group theorem are reproduced with direct global
polynomial checks. Put m=n-2. The group is C_m wr S_n for odd m,
and the subgroup sum(colors)=parity(root permutation) for even m.
It has order m^n*n!/gcd(2,m).

Its exact normal-closure constants over Q are Q(zeta_m). Arithmetic
monodromy includes the cyclotomic unit action on colors; when n=2
modulo4 the parity also changes by(u-1)/2. Thus root-fixed scalar
conjugation alone can give the wrong embedded group.
All proper blocks lie inside a root block and have precisely the
divisor sizes d|m. Every corresponding upper extension is cyclic,
so none supplies a nontrivial polynomial Keller inner factor.
These maps are absolutely and stably atomic for every m.

Their geometric generic deck groups are C_m, whereas the generic
deck group over Q has order gcd(m,2), and every polynomial fixed-target
source deck endomorphism is the identity. Independent products give
actual maps with every finite abelian geometric generic deck group.
These are deck realizations; nontrivial abelian monodromy is still
excluded by KMFI2. This family passes C=N via an actual full N-cycle,
so that necessary curve test does not classify its exact fields.

The [polynomial coefficient graphs](COEFFICIENT_GRAPHS.md), KMFI73--76,
extend the marked-root route beyond constant coefficient slices.
For every n>=3 and arbitrary polynomials lambda_j=h_j(b), the full
inverse-image graph has an explicit ordinary-affine3 source chart.
Exact first-jet coordinate exchanges fix x and have constant nonzero
slopes; there is no degree bound on the graph polynomials.
The resulting entire Jacobian is(-1)^(n-2)*(n-2)/2.

More generally, a graph h_j(q,r,b) has full S_n root group if every
term q^a*b^c surviving r=0 satisfies c+j-1>a*(n-1).
Scaling at q=infinity gives exactly Z^n+Z-B and a simple collision
divisor; Kummer independence gives the full marked G_n and exact
cyclotomic arithmetic constants. Every q-independent polynomial
graph passes, so these entire field routes cannot produce PSL roots.
The actual b-graph maps are absolutely stably atomic, with complete
divisor block lattices and generic cyclic deck group C_(n-2).
Every finite cyclic generic deck group thus has an explicit affine3
atomic witness; products give finite abelian deck groups.

The strict inequality cannot be weakened. The exact graph
h2=-q*b^(n-2)-r*b^(n-3), all other h_j=0, has rational root T=b.
Its marked field is cyclic of degree n-2 over the target, despite
both fields being rational. Its actual smooth affine source component
has nonconstant unit p=1+xy, with inverse1-by, so is not even stably
ordinary affine space. A proper root-group graph must violate the
strict weight gate, but that alone supplies no affine realization.
Exceptional degree11 PSL2(11) and the full spectrum remain undecided.

The [linear coefficient pencils and source-unit obstructions](LINEAR_PENCIL_OBSTRUCTION.md),
KMFI77--80, close additional graph routes outside the strict weight gate.
After removing the full common factor of H,J,L in P=qH+rJ+L,
the basepoint-free plane map has degree e and image degree d, with
N=e*d,d>=2. Sectional monodromy is S_d in characteristic zero.
For e=1 the moving group is S_N; otherwise it is imprimitive.
Thus every proper primitive moving group, including exceptional PSL2(11),
is excluded from this entire affine-linear q,r graph class.
Every irreducible prime-degree member has symmetric root group.

Without an offset, L=T-b, and the full root classification in every
degree is S_n or a fixed b plus S_(n-1), according as H(b),J(b)
do not both vanish or both vanish. The transitive derivative-marked
fields have full G_n and the previously proved cyclotomic constants.

A polynomial pencil with a positive-degree common factor has a
rational marked moving field but an explicit nonconstant cofactor unit
on its actual source component. More generally every actual component
selected by a nonlinear monic factor has a nonconstant unit p=1+xy.
This survives polynomial unimodular root frames: their denominator
linear form becomes the nonconstant unit. Therefore no nonlinear
ordinary-affine Keller component comes from a constant nonzero leading
coefficient base change. Monic curve-factor embeddings can have
precisely computed relative C_m wr G_f groups and rational fields
while their actual sources remain non-stably-affine.
These are construction/source exclusions, not absolute action exclusions.
Genuinely nonlinear q,r dependence, nonmonic static factors or a
different affine architecture remain possible next gates for exceptional11.

The [exceptional PSL2(11) residue analysis](PSL11_RESIDUE_OBSTRUCTION.md),
KMFI81--84, adds an absolute constraint beyond constructor exclusions.
Over algebraically closed constants, for every Keller map of degree n
and prime l not dividing n, a nonzero
target l-Brauer class must have some nonzero residue on the actual
affine source. Restriction is injective by corestriction, and source
Brauer purity forbids all residues disappearing at infinity.

For both exceptional11 actions the complete tame table has affine
inertia2,3,5 only. Decomposition centralizes inertia. Their respective
centralizers have orders12,6,5, giving seven rows. The unramified
source residue degrees are(1,1,1),(1,2),(3),(3) for inertia2 with
D=C2,V4,C6,D12; (1,1),(2) for inertia3 with D=C3,C6; and(1) for
inertia5. Nonempty subsets of these orbits may be retained.

The SL2(11) central-cover Brauer residue is the sign/discriminant
character of all three involution-fixed tails, including omitted ones.
It is nonzero exactly for V4 and D12 decomposition, and survives
on a retained degree1 V4 orbit or the degree3 D12 orbit.
Consequently any noncyclic involution decomposition must be accompanied
by some such actual odd-degree source prime. A nonempty family of V4
divisors retaining only quadratic orbits, with no D12 divisor, is
globally impossible despite being locally allowed at each divisor.

If all involution decompositions are cyclic, the Brauer class vanishes
and a proper SL2(11) field lift exists. Its lifted A5 stabilizer is
perfect, so the earlier cyclic etale-source contradiction is unavailable.
The all-cyclic genus0(2,3,11) passport still passes the new test.
Exceptional11 occurrence and the full action spectrum remain open.

The [complete PSL2(11) curve-action reduction](PSL11_CURVE_ACTIONS.md),
KMFI85--86, now covers every faithful transitive action in the
independently stabilized curve field class. A general prime-subgroup
bound controls the p-part of every cycle-gcd. Since a Sylow2 of
PSL2(11) is V4 and all involutions are conjugate, all even degrees fail.
Exhausting all overgroups of one V4 leaves degrees165,55 and11.
The V4 and A4 actions have no genus0 index budget. The D12 degree55
action has only necessary branch-order budgets(3,3,3) and(2,2,2,3);
both have branch-gcd lcm1 instead of55. Thus only the two exceptional
degree11 actions survive, necessarily with polynomial passport(2,3,11).
The exhaustion retains all660 matrices and4448 subgroup extensions.
This closes alternative curve-action routes, without turning the
surviving passport into an affine realization or excluding arbitrary
higher-dimensional field pairs.

The [finite-normalization boundary theorem](BOUNDARY_CLASS_OBSTRUCTION.md),
KMFI87--89, adds an absolute condition on actual affine completions.
The omitted prime divisors form a free basis for Cl(Z); the
canonical class is their nonzero vector of ramification indices minus1.
Hence a nontrivial Keller normalization cannot have torsion canonical
class, cannot be a normal global affine complete intersection, and
cannot be a normal globally monogenic finite model. These fixed-target
obstructions survive independent stabilization. A raw nonnormal
hypersurface is not excluded, nor is smoothness of its normalization.

The retained cubic map illustrates the distinction exactly: its monic
hypersurface is nonnormal, while its integral closure is a smooth
rank-one determinantal threefold with two affine3 charts, Cl=Z and
primitive canonical class. The unique index2 ramification boundary
is an affine plane whose complement is the original ordinary affine3.

For either exceptional11 action, the complete local prime census
forces normalization class rank at least2. Its canonical vector has
exact divisibility1,2 or4 according as inertia2 occurs, only inertia3
and5 occur with some3, or only inertia5 occurs. For the independent
(2,3,11) curve input, decomposition equals each nontrivial prime
inertia, so rank Cl>=4a+3b>=3. These requirements do not supply a
normalization, an affine source chart or a sufficient action criterion.

The [norm-target reduction](NORM_TARGET_REDUCTION.md), KMFI90--91,
forces every dominant degree11 marked descent of the rational m12
cofactor family to contain the same norm double cover. The
[completed finite-fiber proof](NORM_TARGET_FINITE_FIBERS.md), KMFI92--93,
now proves its exact geometric genus3. A finite ordered-root cover
turns the form into a polynomial times a complete-intersection residue.
Both finite critical fibers, zero-root contacts, two crepant blowups
and all projective coefficient charts are resolved. At infinity exactly
the numerators q^3,q^2*r,q*r^2 survive.

Thus every dominant rational-target base change selecting an
eleven-sheet marked PSL2(11) component in this family is impossible,
even after birational coordinate changes and independent stabilization.
This is stronger than the original component-unit obstruction.
Non-dominant specializations, forgetting the marking and unrelated
field architectures remain outside the theorem. The absolute
exceptional action and complete spectrum are still open.

The [strict-local image theory](LOCAL_IMAGE_AND_MARKS.md), KMFI94--97,
now tests an entire local subgroup at a target point. Every actual
source point supplies a common fixed sheet. The normalization boundary
is pure divisorial, and absence of omitted unramified primes makes
the fixed-sheet mark equal the whole fiber cardinality.

A complete620-subgroup census separates the two exceptional11 actions:
their S3 classes have marks(2,0) and(0,2), despite the same elementwise
characters and four orbits in both restrictions. Minimal bad subgroups
are C6,C11 and the opposite S3 class. Normal-crossing inertia2 and3
forces an omitted stratum, while inertia5 cannot meet2 or3 with the
whole branch support normal crossing. Omission itself is allowed:
the actual cubic map misses exactly the smooth complete-intersection
G_m curve3BC=27AC^2=1. Its local group is the natural transitive S3.

The exceptional action is self-normalizing in S11. Its arithmetic
group must therefore equal geometric PSL2(11), with no additional
normal-closure constants. None of these necessary local or arithmetic
conditions supplies an exceptional affine realization or an absolute
action obstruction.

The [complete A_k branch-stratum test](ARTIN_BRANCH_STRATA.md), KMFI98--99,
now imposes the actual two-meridian braid relation. All184041
ordered allowed-meridian pairs form295 conjugacy orbits; their
relations have periods dividing660, exhausting every A_k singularity.
At an ordinary cusp the noncyclic images are S3 for inertia2,
A4 for inertia3 and either A5 class for inertia5. Both A5 and S3
marks distinguish the exceptional actions.

These restrictions have sharp limits. A tacnode A3 already has
full PSL2(11) local representations, from explicit order5 matrices
whose product has order2. Without inertia5 the first full image is
A9. The resulting target points must be omitted, which is compatible
with the existing affine-image constraints. Local quotient existence
is not an actual Keller realization, and no global action exclusion
has been proved.

The [complete sweep-boundary theorem](SWEEP_BOUNDARY_COMPLETENESS.md),
KMFI100--102, now closes the entire normalized three-variable
tautological polynomial field route to exceptional11, including
nonfinite sweeps of every scalar and primitive-pair degree.
A ramification surface covered by actual horizontal lines forces
a polynomial one-coordinate multiplier and hence a finite sweep.
The projective line-family compactification shows that any omitted
normalization branch over that surface would force just such a ruling.
Thus the sweep is finite near its generic ramification value and
full inertia is exactly one or two three-cycles, with all other
sheets fixed. Exceptional11 needs three three-cycles and is impossible.

The natural projective PSL2(q) sweep-field spectrum is exactlyq2,4.
Odd dihedral natural actions of degree>=5 and F20 also fail the
complete inertia test. These are field/construction exclusions,
not absolute Keller nonoccurrence. Generic properness over the
ramification surface does not mean every sweep is globally finite:
the exact nonfinite control has an escaping horizontal boundary
over a separate plane. Future exceptional11 work must leave the
entire normalized rank2 tautological architecture; increasing a
primitive-pair degree or adding a nonconstant scalar cannot evade it.

The [all-rank completion theorem](HIGHER_RANK_COMPLETION.md),
KMFI103--105, supplies the missing higher-dimensional affine step
for full-rank normalized tautological fields. A finite formal
curve jet and right unipotent Gaussian elimination give an actual
triangular polynomial source automorphism with polynomial inverse.
Larger weights then give global polynomial quotients and the whole
constant Jacobian in every rank. Normal generation of the sheet
stabilizer, an actual unramified branch tail, a regular second
preimage and completion after one common variable are equivalent.

Within this full-rank architecture, exceptional11 now reduces to
rank4 alone. A proved input would automatically give an A6 Keller
map. Tame centralizerC5 makes its ramification surface birational,
forces one additional index5 normalization prime outside the sweep,
and makes the unramified tail irreducible and birational.
Thus it must be nonfinite and its branch hypersurface must have
horizontal line rulings. The rank2 rigidity does not automatically
extend to that four-dimensional hypersurface. No such field input
has been supplied. A fully expanded A6 control verifies the new
construction with exact S3 action and whole Jacobian-3.
Rank-deficient higher directions and unrelated architectures remain
outside the full-rank reduction.

The [exceptional-divisor calculation](COMPLETION_EXCEPTIONAL_DIVISOR.md),
KMFI106--108, now determines actual fibers of the higher-weight
completion on its exceptional hyperplane C0. One source component
maps isomorphically onto all C0, and a disjoint second component
maps isomorphically onto Z0!=0 there. Fibers have exactly2 points
offH and1 onH. Local Kummer base change gives full decomposition
and inertia generated by sigma^(r+1), with trivial residue fields.
The exact class-rank increment is cycles(sigma^(r+1))-1.

For the remaining PSL11 rank4 case, that power is the identity.
All11 normalization primes over C0 split, but only2 are retained.
Nine unramified divisors are necessarily omitted, and the completed
normalization has class rankCl(Z_S)+10>=11. This is a property of
the specified representative, not the original sweep normalization
or arbitrary affine models. A fully specified higher-weight binary
control has ClZ7 and canonical vector(2,0,0,0,0,0,0).
No PSL11 field input or absolute action obstruction is obtained.

The [integrated-power family](INTEGRATED_SWEEP_PENCILS.md),
KMFI109--111, supplies full-rank normalized sweeps in every rank.
Their degree is rd+1. Every sweep with unit-triangular first
r coordinates reduces to exactly that integrated field; the proof
places no degree bound on the remaining polynomial coordinates.

In rank4 their geometric groups are C5 for linear h and A_(4d+1)
for every degree d>=2. The fourth derivative of the critical
value with respect to h(t) recovers 24t, giving generic branch
separation. Arithmetic enlarges to S_(4d+1) exactly when 4d+1
is not a square in the coefficient field, with constants
k(sqrt(4d+1)). Thus the whole triangular rank4 class is finite
and cannot supply the required PSL11 input.

The quadratic h=t^2 member supplies a fully specified ordinary-A6
Keller map over a quartic field, with full-rank input, exact degree9,
geometric/arithmetic A9 and geometric normalization ClZ8. Every
finite rectifier coefficient is retained; global polynomiality
and the whole constant Jacobian follow from KMFI103--104.
The final A6 determinant is not claimed to have been expanded.
A9 was already realized: this is a higher-rank construction control,
not a new abstract action. General nontriangular rank4 inputs
and the absolute exceptional11 action remain open.

The [all-rank integrated classification](ALL_RANK_INTEGRATED_ACTIONS.md),
KMFI112--114, now closes this entire triangular architecture.
The critical-value derivative d^r b/da^r=r!t gives generic
separation, and the decomposition proof supplies primitivity
even for composite cycle length. Jones's composite-cycle
Corollary1.3, with its classification input, then gives symmetric
groups for odd rank and alternating groups for even rank.
Linear h gives the absolutely forbidden regular cyclic case.
Even-rank arithmetic constants are
k(sqrt((-1)^(rd/2)*(rd+1))); odd ranks add no constants.

The residual cofactor at a good critical value gives an exact
coefficient extension of degree at most r(d-1) for an actual
A^(r+2) completion. For r>=2 its geometric normalization has
ClZ^(rd) and canonical vector(r,0,...,0) in the specified
higher-weight representative. This is not a class-group claim
over the nonclosed coefficient field.

Every integrated degree11 field, in every rank, has S11,A11
orC11 action, even after common independent stabilization and
birational recharts. Polynomial equivalences preserving the
normalized sweep-line foliation are affine on directions.
Nilpotency does not justify a Nijenhuis normal form: a valid
full-rank rank4 example has torsion component4t-2.
General nontriangular inputs are not reduced to this class;
no exceptional11 map or absolute obstruction follows.

The [differential reconstruction theorem](SWEEP_DIFFERENTIAL_RECONSTRUCTION.md),
KMFI115--117, now applies to every normalized polynomial input,
including nontriangular ones. The first component is reconstructed
by homogeneous Euler inversion from the remaining vector g.
After solving trace and linear integrability, rank4 has exactly
the three polynomial nilpotency equations trP^2,trP^3,trP^4.
A cyclic coframe df,dfP,...,dfP^(r-1) is invertible precisely
on the full-rank locus.

At a full-rank point the derivative of P along its kernel has
the nonzero eigenvalue df(u). This excludes every affineP and
every matrix image confined to a fixed linear nilpotent space;
a constant simultaneous triangularization is impossible.
All quadraticX inputs have rankDX<=r-1, and an explicit family
attains that bound with a regular nilpotent block.
The Jordan frame at a single branch point remains valid.
Exact Piola identities test all powers ofP; the first-order
matrix compatibility needs a separate closedness equation in
rank2 but already forces it in ranks>=3.

These restrictions cover the entire remaining full-rank rank4
normalized architecture. A candidate needs degreeP>=2 and a
g component of degree>=3, as well as the existing nonfinite
boundary and exact PSL11 group requirements. The coefficient
formulation is not a solution of its nonlinear trace equations
or a sufficient abstract-action criterion. Neither the absolute
exceptional action nor the full spectrum is decided.

The [homogeneous axis theorem](HOMOGENEOUS_AXIS_OBSTRUCTION.md),
KMFI118--120, treats every higher-rank matrix field homogeneous
of one positive ordinary degree n. The axis fiber is reduced
and etale, with one distinguished point and free n-orbits.
Every additional axis point lies on the rank-drop locus ofX.
A conditional Euler--Jacobi argument proves that the whole
homogeneous class is globally nonfinite.

For exceptional11, torus residue and all A5 normal cyclic
quotients leave four positive-axis rows: C2/1,V4/C2,C3/1,C5/1
at matrix degrees2,2,3,5. The V4/C2 case omits four geometric
normalization points of total local generic rank8.
Only degree5 with all11 actual points can be finite etale near
the axis. This is a conditional statement: the one-point axis
case survives in every homogeneous degree, and nonhomogeneous
inputs are not covered.

A rational all-rank full-rank norm construction has just that
one-point axis fiber but a regular cyclic field of degree r+1.
Its exact binomial equation and axis ideal certificate show
why neither full rank nor a one-point special fiber supplies
a Keller realization or determines generic degree1.
No exceptional11 field or absolute action obstruction follows.

The [factorized homogeneous theorem](FACTORIZED_HOMOGENEOUS_FIELDS.md),
KMFI121--123, now closes an entire architecture in the one-point
axis case. If H=g-fp=F(p)Mp with M constant, full rank forces M
to have distinct eigenvalues and F to be a positive monomial in
its eigen-coordinates. The exact field reduces to a polynomial
of degree1+degF; all degree11 actions in this class are S11,A11
or forbidden regularC11. Even-rank arithmetic constants are
k(sqrt(-11)). A remaining homogeneous PSL11 candidate must
have primitive projective coordinate degree>=2 after removing
the coordinate gcd.

The explicit exponent(7,1,1,1) member gives a nonfinite full-rank
degree11 A11 input and an actual ordinary-A6 Keller completion.
Every finite coefficient is computed over an irreducible degree24
field; the whole Jacobian is beta-1 by the completion identity.
Its arithmetic group is S11 with exact quadratic constants.
The original normalization has geometric ClZ12, and the completed
representative has ClZ22 and canonical vector(6,6,6,4,0^18).
A retained residue-degree7 prime is not seven separate primes.
No expanded final A6 determinant or class-group equality over
the nonclosed coefficient field is claimed.

This establishes that a one-point axis fiber can coexist with
actual noncyclic A11 monodromy. It does not realize a new abstract
group or settle PSL11. Unfactored homogeneous and nonhomogeneous
inputs, absolute exceptional11 and the complete spectrum remain open.

The [primitive multiplier theorem](PRIMITIVE_MULTIPLIER_CRITERION.md),
KMFI124--126, removes an unbounded scalar search for each fixed
homogeneous primitive direction. For every full-rank normalized
branch, the Piola identities first make qvec-q0*p generically
cyclic for every nonzero constant target vector. The branch Gauss
image spans the entire dual space and has no fixed cylinder
direction; this does not assert a dominant Gauss map.

Write homogeneous H=FV with primitive V. Full rank forces
detDV=delta*h0^r and the Krylov matrix generated by V under DV
to be invertible. For h=lambda*h0 the characteristic-feedback
identity uniquely determines eta=dlogF. Matching its whole
rational form to nonnegative integer logarithmic residues on
the Krylov determinant is necessary and sufficient for a polynomial
F, unique up to a constant. A separate cyclic coframe tests full
rank. Univariate consistency and residue polynomials retain all
possible scalar degrees and isolated lambda values.

The entire rank4 primitive quadratic binary-shear direction
V=(u^2,a*u*v,b*u*w,c*u*t+Q(v,w)), Q nonzero, fails full rank for
every homogeneous scalar multiplier, including all resonant and
repeated parameters. Exact adjugate identities and incompatible
Darboux poles prove the all-parameter obstruction. A complex
resonant example admits the honest polynomial multiplier F=vw
and rank3, so polynomiality alone cannot replace the coframe.
Other primitive quadratic directions, primitive degree>=3 and
nonhomogeneous inputs remain open. The fixed-direction scalar
criterion does not classify generic groups or settle the absolute
exceptional11 action.

The [one-curved-coordinate theorem](ONE_CURVED_DIRECTION.md), KMFI127,
now closes a larger architecture in every rank and homogeneous degree.
Allow arbitrary projective polynomial Q in the last coordinate of
V=(u^d,a_i*u^(d-1)*x_i,c*u^(d-1)*t+u^d*Q(x/u)).
The arrow resolvent and Krylov determinant have identities independent
of d after dehomogenizing, valid also at resonant eigenvalues.
Polynomial logarithmic residues eliminate every mixed curvature term;
the Darboux pole coefficient (k-1)*(a_i-1)/(c-a_i) eliminates every
remaining univariate degree k>=2 term.

Thus full rank forces Q affine and V=u^(d-1)*M*p.
Every survivor has primitive degree1, with exactly the old positive
monomial multipliers and their constant spectral condition.
All degree11 fields in this class are S11,A11 or regularC11.
Exact full-rank affine survivors at displayed degrees2,3,5,10 and
a nonlinear rank3 control show why neither high displayed degree
nor polynomial multiplier existence can replace the primitive and
full-rank tests. Multiple curved outputs, general nontriangular
primitive directions and nonhomogeneous inputs remain open;
no absolute exceptional11 action obstruction is claimed.

The [coupled triangular projective theorem](TRIANGULAR_PROJECTIVE_RIGIDITY.md),
KMFI128--129, extends this closure to arbitrary numbers of curved coordinates.
For a rational projective field v, put B=Dv and W=v-x.
If W is cyclic, a rational eigen-coordinate dq*B=c*dq,
W(q)=(c-1)q must be affine: its Hessian is self-adjoint under B
and kills W, hence kills every Krylov vector.

For homogeneous directions V0=u^d and
Vi=a_i*u^(d-1)*x_i+u^d*Qi(x1/u,...,x(i-1)/u),
the projective derivative is triangular with constant diagonal.
Full-rank normalization forces W cyclic. A descending adjugate
calculation gives a polynomial multiplier graph divisor at the
last curved coordinate. Its residue equations, with an explicit
quotient-cyclicity proof on that divisor, identify its gradient
with the left eigenrow. The lemma makes that coordinate affine.
Descending removes all curvature, through every repeated or resonant
parameter case. Thus every originalQi is affine and every survivor
has primitive degree1, with exactly the KMFI121 multipliers.
No degree11 member gives exceptionalPSL.

Exact affine coupled full-rank and nonlinear coupled rank3 controls
preserve the hypotheses. This theorem covers every rank, direction
degree and homogeneous scalar degree within that architecture.
Nonconstant projective diagonal, arbitrary nontriangular homogeneous
directions and nonhomogeneous fields remain open. It is a construction
classification, not an absolute exceptional-action obstruction.

The [rank-one projective curvature theorem](RANK_ONE_PROJECTIVE_CURVATURE.md),
KMFI130--132, now closes genuinely nontriangular quadratic directions.
First, a homogeneous quadratic presentation withV0=u^2 and triangular
dehomogenization automatically has constant diagonal under full-rank
normalization: its Jacobian-power factor must be a power ofu.
Thus the preceding triangular theorem applies without a diagonal assumption.

For arbitrary projective v=Mx+b+cQ with quadraticQ, a nonzero DcQ
requires an impossible negative polynomial cofactor forFbar.
The proof includes isotropic Q(c)=0 cases withDcQ nonzero.
Univariate longitudinal curvature is likewise excluded in every degree.
When the projective quadratic part has one-dimensional image in rank4,
the constant Jacobian gate putsM^-1*c in the radical ofQ.
All radical cases either have a triangular invariant flag or reduce
affinely to the companion form v=(y,z-sigma*x+tau*y+x^2,delta*x).
Its unique feedback form has the universal polynomial coefficient
[x^2*y*z^4]D^2*(partial_x eta_y-partial_y eta_x)=8.
The denominator gate has a parameter-independent coefficient-2,
so repeated parameters cannot hide a closed multiplier form.

Hence every rank4 homogeneous quadratic direction in this class,
with matchedV0=u^2, is excluded to full rank in every scalar degree.
A normalized nonlinear rank3 radical example remains valid, so the
whole theorem does not exclude rank-deficient normalization.
Quadratic nonlinear images of dimension>=2, directions lacking that
matched pure-square coordinate, higher primitive degrees and
nonhomogeneous fields remain open. No absolute exceptional11
obstruction or new abstract action has been obtained.

The [two-image projective reduction](TWO_IMAGE_PROJECTIVE_REDUCTION.md),
KMFI133--135, now makes the remaining matched pure-square quadratic rank4
direction class exhaustive. Normalize the quadratic part by the invertible
linear matrixM. De Bondt's source-attributed quadratic nilpotent-Jacobian
triangularization and explicit constant basis changes leave exactly
H_A=(0,x^2,y^2) orH_B=(0,x^2,xy).
The full map remainsv=M*(x+H)+b with arbitrary invertibleM andb.
TriangularizingH does not triangularizev.

A coupled companion family meets both nonlinear shapes and is closed
in every parameter and homogeneous scalar degree: its polynomial
xy-curl numerator has coefficients-8*beta^2 and32*gamma^4.
They force both couplings zero, and a negative cofactor rejects that
remaining branch. Singular Krylov parameters use the separate full-rank
rejection gate, with an exact retained checkpoint.

A broader theorem closes every completely split homogeneous scalar
for polynomial automorphic projectivev, in every rank and degree.
Each linear scalar factor must be an affine eigen-coordinate ofv;
otherwise the trace has a pole. Automorphy makes its image irreducible.
Full rank requires those constant covectors span all projective directions,
forcingv affine and primitive degree1. Both nonlinear normal forms have
explicit polynomial inverses, so a remaining primitive quadratic candidate
needs a nonlinear irreducible scalar factor, of degree at least2.
This leaves generalM,b and such nonlinear factors open; generic groups,
actual exceptional11 fields and absolute action obstructions are not supplied.

The [low-degree multiplier divisor theorem](LOW_DEGREE_MULTIPLIER_DIVISORS.md),
KMFI136--138, now rules out the smallest remaining scalar case.
For polynomial automorphic projectivev, the first two feedback residues
at an irreducible quadratic divisor make its finite Taylor expansion
q(v) divisible byq. Automorphy gives q(v)=constant*q, in every multiplier
degree and without needing full rank for this invariance lemma.

In the remaining two-image rank4 quadratic class, all three possible
ranks of the leading quadratic form then contradict full-rank cyclicity.
A nondegenerate ternary form cannot vanish on the dense quadratic image
plane; ranks1 and2 force a fixed affine direction. Thus no quadratic
factor may occur anywhere insideF. A single non-u linear factor likewise
forces a hyperplane flag in the arbitrary matrixM. The triangular branch
is already closed; the other two branches have impossible negative
constant-derivative or Euler eigenvalues. This handles factors mixed
with nonlinear higher-degree divisors, not only completely splitF.

The precise remaining scalar form isF=constant*u^e*product Di^bi,
with e>=0 and every non-u irreducibleDi of degree>=3, at least oneDi.
The factoru is not excluded: it disappears on the affine projective chart.
An exact normalized nonlinear two-image F=x control has rank3 and whole
Jacobiangamma4, showing why full rank remains indispensable.
This necessary-factor reduction is retained as a component result.
KMFI141 below closes its matched quadratic higher-factor frontier.
Other architectures and the absolute exceptional11 action remain open.

The [complete matched pure-square quadratic theorem](QUADRATIC_PURE_SQUARE_RIGIDITY.md),
KMFI139--141, now closes all higher-degree scalar branches left above.
Every multiplier divisor, in any degree, is preserved by the straight-line
homotopy x+t*(v-x). The proof uses the finite Krylov residue sequence and
a formal inverse only to establish a polynomial ideal identity; no global
polynomial inverse of that homotopy is assumed.

For an ambient projective polynomial automorphism of three variables,
the selected ruling-line family on a nonplane divisor is a curve and
generically supplies one selected line through a point.
The automorphism therefore fixes each selected affine line and acts
affinely on it. If the leading base locus lies in finitely many
hyperplanes, no scalar divisor can have degree exceeding degreev:
otherwise a fixed affine covector would kill the full Krylov span.
Both remaining quadratic normal forms have the required hyperplane
base locus. Their factor degrees are therefore at most2, already
excluded above; pureu powers likewise force an affine projective map.

Thus every full-rank homogeneous quadratic direction with matchedV0=u^2
hasv affine andV=u*M*p. Its scalar and exact generic field reduce
to KMFI121; every degree11 action here is S11,A11 or regularC11.
This closes all arbitraryM,b and all cubic-or-higher factors in that
matched class. An affine full-rank rank4 norm control and a nonlinear
normalized rank3 control preserve both sides of the statement.
An unrelated square output does not establish the matched hypothesis.
Other quadratic directions, higher primitive degree, nonhomogeneous
fields and the absolute exceptional11 action remain open.

The [unmatched quadratic chart theorem](UNMATCHED_QUADRATIC_CHARTS.md),
KMFI142--145, now reduces the whole determinant-compatible quadratic
direction class. A midpoint calculation with both signs proves projective
injectivity wherever the critical linear coordinateu is nonzero.
The map is etale there and hence an open immersion; the affine-space
unit and divisor arguments make its image the complement of one target
hyperplanealpha, withalpha(V)=u^2.
The two hyperplanes need not agree. The radial degree2 ofV itself
is not the generic degree or monodromy of the Keller sweep.

In rank4 the projective chart automorphism and the retained nilpotent
quadratic classification give exactly four primitive forms V=A*K:
chained squares, a mixed chain, a single square shear and a product shear.
Every invertible GL4 output matrixA is retained.
The matched subcase is closed, while the independent target matrix
cannot be removed by a target-only change in a normalized sweep.

The full-source divisor homotopy requires no pure-power component:
D(p+tV) is divisible by every scalar factorD.
The chart isomorphism then gives D(V)=constant*u^degD*D off the source
boundary; an unmatched target-boundary factoralpha is impossible.
A projective ruling-line and base-point argument bounds each non-u
factor degree by3, independently of its multiplicity and totaldegF.
Larger factors would force a matched square and the closed KMFI141 case.
Linear, quadratic and cubic factors and powersofu remain in the unmatched
frontier; the old matched lower bound is not transferred to it.
The arbitrary matrix/factor solutions, their actual generic groups,
ordinary affine witnesses and the full action spectrum are still open.

The [pure critical-boundary multiplier theorem](PURE_BOUNDARY_MULTIPLIER.md),
KMFI146--148, closes all pure powers F=c*u^s for primitive quadratic rank4
full-rank normalized inputs, without an exponent or output-matrix bound.
A positive boundary factor forces V_u divisible byu. After conjugation
by diag(u,1,1,1), the nilpotent feedback specializes to a weighted matrix.
If source and target boundaries differ, its high quadratic image lies
in a plane. A binary-Jacobian and trace/determinant argument excludes
nilpotency at every nonnegative weight unless the offending linear form
vanishes, which forces the matched chart already closed byKMFI141.

Mixed factors have an additional rank-one term involving dlogG0.
The exact computation retains it, so this theorem does not close mixed
multipliers. Every remaining primitive quadratic candidate needs at
least one non-u factor of degree1,2 or3; powersofu may accompany it.
Negative weights-1,-2 are genuine matrix exceptions but not polynomial
scalar exponents. Actual normalized primitive quadratic rank3 controls
exist at every positive exponent with regular nilpotent feedback and
whole Jacobiangamma4, retaining the full-rank hypothesis.
Exceptional11 and the complete action spectrum remain open.

The [mixed critical-boundary theorem](MIXED_BOUNDARY_MULTIPLIERS.md),
KMFI149--150, now closes all accompanying u powers as well.
The boundary matrix keeps its entire logarithmic rank-one correction.
After factoring out the exact image-plane coordinate from the scalar,
the tangential binary determinant has coefficient(e+degj+2)/2.
Nilpotency makes the binary quadratics proportional; both branches
then force a negative polynomial logarithmic residue.
Thus the boundary linear form must annihilate the image plane.

Single-shear forms have boundary image dimension1, so the form vanishes
and the charts match. The two chain forms instead acquire a2x2 output
flag. Their entirez-infinity matrix retains the logarithmic derivative
of the highest-z coefficient. The square chain contradicts a polynomial
coefficient, while the mixed chain requires residue-(e+1).
Hence normalization alone forces matched charts wheneveru dividesF.
Primitive full rank is already excluded there byKMFI141.

Every remaining primitive quadratic rank4 full-rank scalar is therefore
nonconstant and coprime tou. Its irreducible factor degrees are1,2 or3,
with the unmatched targetalpha forbidden; total degree and multiplicities
remain unbounded. General nonboundary factor arrangements, actual generic
groups, higher directions, nonhomogeneous fields and the full action
spectrum are still open. No new abstract occurrence group is claimed.

The [complete chained-square theorem](CHAINED_SQUARE_QUADRATIC_RIGIDITY.md),
KMFI151--152, closes the entire arbitrary-matrix family
V=A*(u^2,ux,uy+x^2,uz+y^2), in every homogeneous scalar degree.
A unique projective base point first excludes pointwise-fixed divisor
surfaces: a curve section avoiding the base point would have incompatible
pullback degrees under an identity represented by quadratic sections.
Each nonlinear divisor then has a selected ruling by degree1 preserved
lines, all passing through that point, so its equation is independent ofz.
Weighted preservation separately puts every linear factor inspan(u,x).

The fullz-trace of the feedback reducesF to a binary homogeneous scalar,
hence a product of linear factors. A critical-hyperplane invariant branch
forces a matched square. The other branch has exactly one linear factor,
made into x by a simultaneous constant conjugacy preserving this family.
The first trace forcesb0; the x-square coefficient of the whole second
trace is2*(s+2)*c*d, nonzero for all remaining parameters.
Normalized rank3 chained-square controls still have wholeJacgamma4,
so full rank remains essential.

The remaining primitive quadratic rank4 directions are the mixed chainK_B,
square shearK_C and product shearK_D with arbitrary invertible outputA.
Their scalars must be nonconstant and coprime tou, with irreducible factor
degrees1..3 but unbounded multiplicities. Their actual fields/groups,
higher primitive degrees, nonhomogeneous inputs, exceptional11 and the
full faithful-transitive spectrum are not decided.

The [complete single-square theorem](SINGLE_SQUARE_QUADRATIC_RIGIDITY.md),
KMFI153--154, closes V=A*(u^2,ux,uy,uz+x^2) for arbitrary invertibleA
and every homogeneous multiplier degree. A preliminary polynomial lemma
excludes nonlinear fixed hypersurfaces forV=uMp+a*ell^2 whenMb=a,
u(b)0 andu(a) nonzero, in every dimension>=3.
Elementary exterior coefficient identities force a fixed quadric
to haveV=L*p+D*c; transversality makesM rank1, impossible.

A scalar divisor's selected ruling line has baseq onB=u=x0 and containsMq.
Varyingq gives either a plane or a skew-directrix quadric.
The latter has an explicit limiting line inu0 sent outsideu0, contradicting
the polynomial line-preservation identity. A fixedq is an invariant cone
vertex; quotienting by itsM-eigenline gives a forbidden nonlinear fixed
curve in dimension3. Thus every scalar factor is linear.
Weighted preservation makes them annihilate botha andb.
The derivative of the whole trace equation alongb is thenu(a)*F0,
the final contradiction. A normalized rank3 control over the exact field
beta^2-beta+1=0 has whole sweep Jacobiangamma4 and retains full-rank necessity.

The remaining primitive quadratic rank4 families areK_B mixed chain and
K_D product shear, with arbitrary output matrices. Their nonconstant
multipliers are coprime tou and have irreducible factor degrees1..3,
with total degree and multiplicities unbounded. Their actual generic
actions, higher directions, nonhomogeneous inputs and the absolute
exceptional11 action remain unresolved.

The [complete primitive quadratic theorem](QUADRATIC_PRIMITIVE_RIGIDITY.md),
KMFI155--157, now closes every primitive homogeneous quadratic rank4
direction for arbitrary matrices and every homogeneous scalar degree.
A general polynomial-density lemma reads the highest translation
coefficient: ifpartial_bV=u*a,u(b)0,u(a) nonzero andR is independent
ofb, then div(FV)=RF would require a polynomialu-Euler eigenvalue-1.
It forcesA03 zero in every canonical family without assumingF_z0.

A second lemma scales normal coordinates along a linear stratum where
V,h vanish. If its normal derivativeB is invertible, the nilpotent
feedback would have a normal limiting determinant(m+1)*detB, nonzero
for every polynomial vanishing orderm>=0. The upper-right tangential
block is retained; both lower blocks vanish.
For the remaining mixed chain this forcesA13 zero. Its highest total
(y,z)-degree trace then has another impossible negative Euler eigenvalue.
The single shears reduce directly to the closed matched square case.

Thus all four primitive quadratic rank4 families are excluded.
Nonprimitive quadratic survivors are exactlyV=u*M*p, withG=uF among
the KMFI121 scalar-times-linear data. Their exact generic field is its
retained one-variable polynomial with common spectators; no new
all-degree group classification is inferred. Degree11 at rank4 isA11.
A homogeneous rank4 exceptional11 input must now have primitive
direction degree>=3. Normalized lower-rank quadratic inputs and
nonprimitive full-rank inputs remain valid.
Higher primitive directions, genuinely nonhomogeneous fields, actual
exceptional11 occurrence or an absolute action obstruction and the full
faithful-transitive spectrum remain open.

The [nilpotent normal-stratum theorem](NILPOTENT_NORMAL_STRATA.md),
KMFI158--161, strengthens the all-degree test beyond invertible normals.
For every polynomialJ, B+(Bq)*dlogJ is nilpotent exactly whenB is
nilpotent and alldJ*B^j*q vanish. A nonzero invertible Fitting component
would force a polynomial Euler eigenvalue-1. For a regular nilpotent
block, a homogeneousJ is a pure power of its left-kernel coordinate.
This concerns the initial normal form, without global divisibility ofF.

The full zero-stratum limit therefore forces singular normal derivatives
to be nilpotent too. Every nonzero zeroa of a homogeneous primitive
direction satisfies nilpotentDV(a), by applying the test to its radial
line. A second extension permits any polynomial coefficientP in
partial_bV=P*a whenpartial_bP=partial_bR0 andpartial_aP nonzero.

For the triangular cubic directionA*(u^3,u^2x,u^2y+ux^2,u^2z+xy^2),
these tests forceA02=A03=A13 zero. Both relevant normal determinants
are identically zero, so their nonzero eigenvalues require the stronger
criterion. An explicit swapped-output member hasdetDV=-3u^8 andunit
Krylov determinant1, but base-point derivative diag(1,0,0,0).
It is excluded in every scalar degree. Normalized primitive cubic
rank3 controls still have whole sweep Jacobiangamma4.

The remaining matrices even in that cubic family are not classified.
All quadratic shapes of the fourth-power baseh0 for arbitrary cubicV,
higher primitive degrees, genuinely nonhomogeneous fields and the
absolute exceptional11 action remain open. No occurrence group or
complete action criterion is supplied by these necessary tests.

The [complete triangular cubic-chain theorem](CUBIC_TRIANGULAR_CHAIN_RIGIDITY.md),
KMFI162--163, now closes every arbitrary output matrix and scalar degree
in the retained familyA*(u^3,u^2x,u^2y+ux^2,u^2z+xy^2).
A general last-coordinate test applies toV=(v(q),c(q)z+R(q)).
Cyclicity makes the eigenvalue adjugate componentW=Kz+L nonzero even
at repeated roots. Full rank rejectsc=h, polynomiality rejectsK0,
and degree in z gives the exact positive-integer identitynu*K=-(c-h)^r.
The single-root expression overk(q) does not remove denominator checks.

In the cubic family, a highest-y trace first removesA23.
The whole adjugate coefficient then forcesA12 zero: its other candidate
makesK identically zero, a retained rejection branch.
Highest weighty1,z2 followed by highestz in the trace forcesA01 zero.
The matrix is now lower triangular, and its matched pure-cube
constant-diagonal projective chart has nonzero quadratic curvature,
contradictingKMFI129. Normalized cubic rank3 controls remain valid.

Only this full cubic architecture has been closed. Other cubic directions,
including otherh0=u^2 cases and all other quadratic shapes ofh0,
higher primitive degrees, genuinely nonhomogeneous inputs, actual
exceptional11 or an absolute obstruction and the full action spectrum
remain unresolved. No new occurrence group is claimed.

The [rational linear-fiber theorem](RATIONAL_LINEAR_FIBER_RIGIDITY.md),
KMFI164--165, now closes a family with an irreducible nonlinear critical
factor. A rational homogeneous degree1 directionW=(Cq,mu*z+T(q)) with
polynomial normalized full-rank input has cyclicDW,W and cyclic headq.
Its graph shiftS is fixed by the last-coordinate adjugate.
The entire graph residue sequence forcesdT=dS*(muI-C).
Symmetry and Euler giveH_S*C=C^t*H_S andH_S*q0; head cyclicity makes
the Hessian vanish. ThusT is linear, without diagonalizingC or applying
a nonlinear source change to the normalized input.

ForV=QMp+(Mb)R, where irreducible nonlinearQ has exactly one constant
translation directionb andR isb-independent ofdegreedegQ+1,
the Jacobian is(degQ+1)*detM*Q^r. Translation forcesMb in that kernel,
so division byQ gives the rational fiber case with polynomial scalarFQ.
It forcesQ dividingR; the remaining scalar-times-linear input has the
forbidden nonlinear scalar factorQ. AllM,R and scalar degrees in this
architecture are therefore excluded to full rank.

A concrete primitive rank4 cubic hasJacobian150*(x^2+y^2+w^2)^4,
nonzero exact Krylov determinant andsquare-zero derivatives at every
projective basepoint. It still fails the full graph-gradient condition.
The earlier local gates are necessary, not sufficient.
Other directions with rank3 critical quadrics, nondegenerate rank4
quadrics, higher degrees, nonhomogeneous inputs and the absolute
exceptional11 action remain open. No occurrence group is claimed.

A complete theory must decide every faithful transitive action, or give a
necessary-and-sufficient criterion with a proof of sufficiency. KMFI1--165
give necessary conditions, exclusions and verified occurrence families,
but do not prove a sufficient criterion for every action. Alternating
simple realizations are established. Natural degree12 PSL2(11) is
now excluded absolutely. The remaining PSL2(11) gate is an actual
exceptional degree11 map or a new absolute obstruction for that action. Repeating
completion of the natural F11 curve field after independent stabilization
is excluded. An allowed permutation action or a determinant-one chart
still does not supply the required absolute realization.
