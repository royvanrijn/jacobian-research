# A Keller inverse-Galois program

## 1. The stretch beyond fiber transfer

The
[finite-étale fiber theorem](verified/FINITE_ETALE_KELLER_FIBERS.md)
starts with a squarefree degree-\(N\) polynomial \(P\), transfers

\[
 \operatorname{Spec} K[T]/(P)
\]

into a full fiber of an explicit polynomial Keller map, and preserves the
splitting field and every local property of that one fiber. This settles
occurrence of an already-existing finite étale algebra. It does not prescribe
the monodromy of the generic inverse cover surrounding the fiber. In the
weighted linear-pencil realization, that generic monodromy is always
\(S_N\).

The next program is therefore:

> **Keller inverse-Galois problem.** Let \(K\) be a number field and let
> \(G\leq S_N\) be a prescribed transitive permutation group. Construct an
> explicit Keller family whose generic degree-\(N\) inverse cover has
> monodromy \(G\), while retaining full fibers and allowing compatible real
> and finite-place splitting conditions to be imposed on rational fibers.

The permutation representation is part of the input. Abstractly isomorphic
groups with inequivalent transitive degree-\(N\) actions define different
problems.

The
[Beckmann--Black specialization refinement](KELLER_BECKMANN_BLACK_SPECIALIZATION.md)
makes the uniform arithmetic quantifiers precise.  The literal regular-action
version is impossible for a nontrivial absolute Keller map by the classical
Galois-case theorem.  For a core-free point stabilizer \(H<G\), however, a
generic polynomial for the action on \(G/H\) has a universal
derivative-unit suspension: one determinant-one morphism of smooth affine
charts preserves the generic \(G\)-cover and realizes every point field
\(L^H\) as a complete fiber.  This settles the chart version for
\(A_4,D_5,F_{20},A_5\); absolute affine-space completion remains open in
all four cases.

This is not a claim that the general inverse Galois problem has been solved.
The first target is a **transfer theorem for an already-constructed explicit
regular \(G\)-cover**. If such a cover over \(K\) is not known, Kellerization
does not manufacture it. Conversely, an unconditional construction over
\(\mathbb Q\) for every finite group, together with specialization retaining
the full arithmetic group, would have consequences for the classical inverse
Galois problem and cannot be treated as a merely geometric corollary.

## 2. What “generic inverse monodromy \(G\)” should mean

Let \(B\) be an integral \(K\)-variety and consider a relative polynomial map

\[
 \mathcal F:\mathbb A^m_B\longrightarrow\mathbb A^m_B
\]

whose relative Jacobian determinant is a unit in \(K\). On a nonempty target
open \(U\), suppose that \(\mathcal F^{-1}(U)\to U\) is finite étale of degree
\(N\). Its geometric generic fiber gives a transitive subgroup

\[
 G_{\mathrm{geom}}\leq S_N,
\]

and its arithmetic function-field extension gives
\(G_{\mathrm{arith}}\leq S_N\), with \(G_{\mathrm{geom}}\triangleleft
G_{\mathrm{arith}}\).

The clean target is

\[
 G_{\mathrm{geom}}=G_{\mathrm{arith}}=G. \tag{2.1}
\]

When equality cannot be obtained over \(K\), the geometric and arithmetic
groups must be stated separately. A construction with geometric group \(G\)
and a larger arithmetic normalizer is not a realization of (2.1).

If \(B\) is an affine space and the coefficients of \(\mathcal F\) are
polynomial in its parameters, adjoining the parameters as unchanged
coordinates,

\[
 (x,b)\longmapsto (\mathcal F_b(x),b),
\]

turns the relative construction into one absolute polynomial Keller map.
When \(B\) is a punctured curve, a boundary complement, a stack, or a
nonrational Hurwitz space, this promotion is a genuine affine-completion
problem rather than a formal step.

## 3. The section-transfer trap

It is essential to distinguish the following statements.

1. A \(G\)-extension occurs as one complete fiber of a Keller map.
2. A positive-dimensional family of \(G\)-extensions occurs along a section
   of the target.
3. The generic inverse cover of the Keller map itself has monodromy \(G\).

The finite-étale theorem proves the first statement for every already-given
finite étale algebra of rank at least three. Uniform coefficient formulas may
sometimes prove the second. Neither implies the third. Adding independent
target directions can enlarge the generic group, and the two-parameter
inverse pencil

\[
 H(W)-sW+t
\]

has generic group \(S_N\) for every degree-\(N\) polynomial \(H\).
Consequently a proper subgroup \(G<S_N\) cannot be obtained merely by choosing
a special seed or by imposing arithmetic conditions on isolated targets of
that pencil.

The new geometric problem is **cover-preserving Kellerization**: the
constant-Jacobian completion must retain the function-field extension of the
input \(G\)-cover on a dense target open, not only recover selected vertical
fibers.

## 4. A cover-first transfer theorem to seek

Start with explicit data

\[
 \pi:X\longrightarrow B
\]

such that:

- \(B\) is a geometrically integral \(K\)-variety;
- \(\pi\) is a generically finite degree-\(N\) cover with prescribed
  geometric and arithmetic monodromy \(G\leq S_N\);
- after removing an explicit discriminant divisor \(D\), the restriction
  \(X^\circ\to B^\circ=B\setminus D\) is finite étale;
- a primitive element gives a monogenic presentation
  \[
    X^\circ=\operatorname{Spec}_{B^\circ}
    \mathcal O_{B^\circ}[T]/(P_G(T;b));
  \]
- the branch, pole, and reconstruction divisors are explicitly controlled.

The desired transfer theorem would construct a constant-Jacobian morphism
between smooth affine spaces, or first between controlled affine charts,
whose inverse cover over a dense open is isomorphic to \(X^\circ\to
B^\circ\). It should preserve:

1. the degree-\(N\) permutation representation and its Galois closure;
2. the full scheme fiber
   \(\operatorname{Spec}K[T]/(P_G(T;b))\) at every certified regular target;
3. specified marked sheets or subgroup quotients;
4. the local factorization data at all selected good primes;
5. explicit bad divisors and a verifiable constant Jacobian.

There are three increasingly strong outputs:

- **relative output:** a Keller morphism over the original base or a finite
  marked cover of it;
- **chart output:** a determinant-one morphism of smooth affine boundary
  complements with the exact \(G\)-cover;
- **absolute output:** a polynomial Keller self-map of affine space whose
  generic inverse cover has group \(G\).

The
[global Davenport construction](extended-geometry/GLOBAL_SUNADA_KELLER_COVERS.md)
already reaches the first two kinds of output for the point and line actions of
\(\operatorname{GL}_3(\mathbb F_2)\). Its affine-space descent obstruction
shows why the third output is substantially stronger.

The
[\(\operatorname{PSL}_2(11)\) action-spectrum benchmark](extended-geometry/PSL2_11_KELLER_ACTION_SPECTRUM.md)
adds a second global Gassmann chart.  Two nonconjugate exceptional \(A_5\)
stabilizers give conjugate degree-eleven Shabat polynomials and determinant-one
common-target Keller charts, while the natural degree-twelve action of the
same group and the same \((2,3,11)\) closure is the explicit elliptic modular
curve \(X_0(11)\).  The degree-five/six direct correspondences normalize to
genus-one \(A_4\)- and genus-two \(D_{10}\)-quotients, with reduced
five/eight-node conductors.  The genus-one normalization is the
conductor-\(121\) curve \(v^2+uv=u^3+u^2-2u-7\), has \(j=-121\), and is
not isogenous over \(\mathbb Q(\sqrt{-11})\) to \(X_0(11)\).  The genus-two
normalization is bielliptic: its Jacobian is \((2,2)\)-isogenous to the
product of the conductor-\(11\) class containing \(X_0(11)\) and the CM
conductor-\(121\) class of discriminant \(-11\).  The three positive-genus
boundary-unit lattices are exact: ranks \(3,14,17\), with an index-six
saturation on \(X_0(11)\) and quotient
\(\mathbb Z^3\oplus\mathbb Z/5\) in genus two.  Along both correspondence
projections, the two rank-six unit images meet in the two common base units;
their rank-ten sum is primitive and has free cokernel of rank four/seven in
the genus-one/genus-two unit lattice.  Exact normalization-module
interpolation now constructs all four/seven residual masks.  On the
genus-two quotient, the difference of the two infinity coefficients
descends through both projection images and has values
\((-2,0,0,0,0,1,1)\), proving that three asymmetric pole profiles are
intrinsic.  Thus even the geometry between two Gassmann actions carries new
obstructions.  The action and its quotient geometry, not the abstract group,
are the correct unit of classification.

For the natural prime actions the adjacent genus comparison is now uniform:
the standard \((2,3,p)\) modular triple has
\[
 g_{\mathrm{nat}}(p)
 ={p-6-3\left(\frac{-1}{p}\right)
        -4\left(\frac{-3}{p}\right)\over12}.
\]
Hence the natural genera at \(p=7,11,13\) are \(0,1,0\).  The obstruction at
\(11\) is a split/nonsplit arithmetic event, not a monotone group-size
phenomenon; the next comparison for \(7\) and \(13\) must use boundary and
affine-completion data.

<!-- status-consumer: KAS1 f40f8588d37ade00 -->
<!-- status-consumer: KAS2 f56459cc921661ea -->
<!-- status-consumer: KAS3 b95d888270f98c59 -->
<!-- status-consumer: KAS4 45a513f714702919 -->
<!-- status-consumer: KAS5 2baa200b6712564f -->

## 5. The proposed construction pipeline

### Step A: choose an explicit \(G\)-cover

Use a regular realization, a generic polynomial, a Hurwitz family, a
resolvent construction, a quotient of a generically free linear action, or a
known exceptional cover. Record both geometric and arithmetic monodromy in
the chosen degree-\(N\) action.

### Step B: add a marked primitive sheet

Pass to a marked-root incidence when necessary. The mark supplies a primitive
element, makes reconstruction explicit, and rigidifies deck transformations.
The finite marked cover of coefficient or Hurwitz space should be tracked
rather than silently absorbed into the ground field.

### Step C: solve the Jacobian ledger

Compare the source ramification, target discriminant, primitive derivative,
and boundary-unit divisors. Add the minimum number of Cox, Kummer, or affine
modification coordinates needed to make the determinant divisor principal
with constant residue. This is the point at which cover preservation must be
checked on function fields.

### Step D: control fullness

Exhibit an inverse-polynomial presentation on a reconstruction open and prove
scheme-theoretically that every squarefree specialized root remains in that
open. The fiber length must equal the generic degree; no inverse sheet may
escape through the completion boundary.

### Step E: impose arithmetic conditions

For a finite set of good primes, choose conjugacy classes of \(G\) in its
degree-\(N\) action and realize the corresponding factorization partitions.
Combine the resulting residue conditions with an allowed real chamber and a
Hilbert subset retaining arithmetic monodromy \(G\). Constructive weak
approximation then produces a rational target when the chosen parameter
chart has the required local points.

Ramified local extensions, specified inertia groups, and decomposition groups
belong to a later Grunwald layer. They are not consequences of the current
unramified factorization-type theorem.

## 6. How the existing machinery enters

- [**Finite-étale fiber transfer**](verified/FINITE_ETALE_KELLER_FIBERS.md)
  supplies exact scheme reconstruction, fullness criteria, scalar-extension
  compatibility, and an effective constant-Jacobian model for individual
  algebras.
- [**Marked-root geometry**](MARKED_ROOT_KELLER_MAPS.md) supplies primitive
  sheets, normalized incidence covers, deck rigidity, and explicit control
  of what a finite marking base change does to the cover.
- [**Adelic engineering**](verified/ADELIC_FIBER_ENGINEERING.md) supplies real
  chambers, finite-prime factorization conditions, constructive CRT lifts,
  and height bounds once a suitable rational parameter chart is available.
- [**Sunada-Gassmann
  machinery**](extended-geometry/GLOBAL_SUNADA_KELLER_COVERS.md) supplies
  nonisomorphic subgroup covers with a common Galois closure and identical
  Frobenius statistics. It tests whether Kellerization preserves more than
  the abstract Galois group.
- [**Cox ledgers and affine
  modifications**](extended-geometry/DAVENPORT_COX_BOUNDARY_OBSTRUCTION.md)
  address the main geometric obstruction: converting a determinant-one
  morphism on a boundary complement into a polynomial self-map of affine
  space without changing the inverse cover.

These ingredients divide the program cleanly into arithmetic input,
cover-preserving Jacobian completion, and affine-space descent.

## 7. Programme 3 queue and current checkpoints

The examples now motivate an action-level problem rather than only a
construction queue:

> **Keller monodromy action-spectrum problem.** Classify faithful transitive
> permutation actions which occur as geometric generic monodromy of
> noninvertible polynomial Keller self-maps of affine space.  Determine the
> obstruction to promotion from determinant-one smooth-affine charts to
> absolute affine-space maps.

The first filtration records regularity and the block interval, Hurwitz
inertia and quotient genus, stable birational type, the boundary-unit
pullback and derivative class, arithmetic normalizers, and essential
dimension only at the versal layer.  These fields are developed and tested in
the [action-spectrum note](extended-geometry/PSL2_11_KELLER_ACTION_SPECTRUM.md).

The systematic queue is:

| stage | groups or actions | first model to audit | current repository state |
|---:|---|---|---|
| 1 | cyclic and dihedral | Kummer and Dickson/Chebyshev incidence | cyclic absolute no-go; uniform dihedral ledger; canonical one-normal completion excluded in every degree |
| 2 | \(A_4,S_4,A_5\) | oriented discriminants and low-dimensional generic polynomials | \(A_4\) has the advanced frontier below; \(S_4\) has a six-edge collision core and a decomposable absolute group-only checkpoint, but its ordinary collision-frame algebraization is open; \(A_5\) absolute descent is queued |
| 3 | Frobenius groups | semidirect-product generic polynomials | \(F_{20}\) has a universal smooth-affine chart, a complete finite exceptional-color audit, and a rational connected \(q\)-conductor cover; the complete packet closure compresses to \(3D_d+D_q+D_r=\operatorname{div}(P_X)\), while the natural incidence algebra and its finite normalization fail the value-one gate; thirteen controlled-transform chart types supply primitive exceptional variables for every positive exceptional color and verify literal local \((3,1,1)\) cancellation; seven bivariate corner charts close the exceptional cusp and \(q\)-\(r\) intersections, and weighted Taylor--Cox charts now attach cusp \(E_4\) to strict \(r\), triple \(E_2\) to strict \(d\), and the \(q\)-\(r\) \(A\)-packet to strict \(q\); the \(q\)-node incidence is proved conductor-degenerate, the apparent positive triple \(E_1\)-\(E_2\) edge is deleted by root-center separation, and a compact residue product passes on the punctured conductor; the remaining strict-\(r\)/conductor attachments, global Cox ring, affine completion, and absolute descent remain open |
| 4 | \(\operatorname{PSL}_2(\mathbb F_q)\) | natural projective-line and low-index subgroup torsors | the natural prime-triangle genus is explicit and gives \(0,1,0\) at \(q=7,11,13\); at \(q=11\), the degree-12 quotient is \(X_0(11)\), the two exceptional degree-11 \(A_5\)-coset actions form a Gassmann pair with Shabat and determinant-one chart outputs, their degree-5/6 correspondences normalize to genus 1/2, all residual masks are explicit in the normalization algebras, and a descended infinity imbalance obstructs uniform poles on three genus-two classes; boundary audits at \(7,13\), ambient polynomial assembly, and absolute descent are open |
| 5 | affine linear groups | natural affine actions and additive-polynomial invariants | \(F_{20}\) is the first overlap with stage 3; its \((4,1)\), unramified-crossing, and \((2,2,1)\) boundary colors are compiled, while higher-rank audits are queued |
| 6 | selected Gassmann pairs | common-closure point/subspace actions | the degree-seven \(\operatorname{GL}_3(\mathbb F_2)\) pair has relative and Cox outputs; absolute descent is open |

Every row is to be reported with the common eight-field card: invariant
ring, discriminant/orientation cover, boundary unit lattice, class group,
derivative units, essential dimension, affine-modification candidates, and
the exact geometric/arithmetic monodromy distinction.

The current concrete checkpoints are:

1. [**Cyclic and dihedral
   actions.**](extended-geometry/ABSOLUTE_INVERSE_GALOIS_CYCLIC_DIHEDRAL_AUDIT.md)
   Every faithful transitive cyclic action is regular, so the Galois-case
   theorem excludes nontrivial cyclic absolute inverse monodromy altogether.
   The natural degree-\(n\) dihedral action is nonregular.  Its Dickson
   incidence gives polynomial invariant rings, an exact discriminant and
   boundary-unit ledger, and a derivative-unit determinant-one chart.  The
   odd orientation cover has an \(A_{n-1}\) class-group obstruction, while
   direct one-coordinate polynomial pole clearing is excluded in every
   degree.  The
   [all-degree affine-completion audit](extended-geometry/DIHEDRAL_ALL_DEGREE_AFFINE_COMPLETION_OBSTRUCTIONS.md)
   further excludes the canonical two-mask construction with any single
   normal coordinate, including every nonlinear resonant degree.  The
   remaining canonical search requires two coupled normal directions;
   alternative determinant-\(\Delta_n\) blowdowns remain separate.
2. [**\(A_N\) in its natural
   action.**](extended-geometry/A4_KELLER_INVERSE_COVER.md) The
   square-discriminant condition is conceptually minimal, and the oriented
   discriminant/Cox machinery is already designed for it. In degree four the
   exact \(A_4\) inverse cover now has a determinant-one derivative-unit
   suspension on smooth affine boundary complements, a complete \(A_4\)
   fiber, and all three unramified cycle types. The
   [affine follow-up](extended-geometry/A4_AFFINE_KELLER_FRONTIER.md)
   further gives a polynomial \(\mathbb A^3\)-map with generic \(A_4\)
   monodromy and determinant \(4W^2K^3L\). The
   [exact ledger reduction](extended-geometry/A4_LEDGER_REDUCTION_AND_RIGIDITY.md)
   absorbs \(K^3L\) into a target divisor and leaves the residual source
   boundary \(WL\). It closes all defect-multiple ambient corrections and
   forces any surviving stabilization to couple new variables back into the
   cover coordinates. The
   [pure-target lift](extended-geometry/A4_PURE_TARGET_LEDGER_LIFT.md)
   then adjoins one coordinate and makes the full determinant exactly the
   pullback of that target divisor, without changing the \(A_4\) extension.
   Factoring this log-Keller map into an ordinary Keller map requires at
   least two source-dependent masks entering the cone outputs. The
   [two-mask frontier](extended-geometry/A4_TWO_MASK_FACTORIZATION_FRONTIER.md)
   supplies the puncture-adapted birational target blowdown with determinant
   \(\mathcal B\), excludes every rechart preserving the mask zero section,
   and reduces the next construction to an incidence embedding. The
   [normalized-boundary assembly audit](extended-geometry/A4_NORMALIZED_BOUNDARY_ASSEMBLY_AUDIT.md)
   gives the smooth normalization an explicit determinant-one ambient
   completion, but proves that the resulting inverse blowdown is
   nonpolynomial with nonconstant rational Jacobian. More generally, every
   automorphic incidence embedding is excluded. Homogeneous radial scaling
   is an explicit nonautomorphic solution of the corrected log-Jacobian
   equation, but an all-degree order argument excludes the entire
   `(SP,SQ,SR,A,C)` family for arbitrary polynomial masks. Any surviving
   incidence must therefore couple the three cubic base coordinates
   nonradially.
3. [**Decomposable and collision-frame proper
   groups.**](extended-geometry/S4_COLLISION_FRAME_KELLER_FRONTIER.md) The
   composition \(F\circ F\) of the foundational degree-three map is an
   absolute degree-nine Keller map with determinant four.  Its three blocks
   of three sheets give monodromy inside \(S_3\wr S_3<S_9\), with nonabelian
   quotient \(S_3\).  This is a literal but decomposable group-only
   checkpoint.  Independently, factoring a depressed quartic into two
   quadratics gives the exact \(S_4<S_6\) edge cover, together with its
   discriminant, primitive-conductor, two-normal, rational cotangent, and
   relative logarithmic ledgers.  Complementary edges form a three-block
   system, so no group-theoretic atomicity is claimed.  A one-normal
   zero-section-preserving completion is impossible; an ordinary polynomial
   Keller algebraization with a proved affine-space source remains open.
4. **\(\operatorname{GL}_3(\mathbb F_2)\) in degree seven.** This is the
   existing nonsymmetric benchmark. The relative and Cox realizations are
   known; absolute affine-space descent remains open.
5. [**\(\operatorname{PSL}_2(11)\) in degrees eleven and
   twelve.**](extended-geometry/PSL2_11_KELLER_ACTION_SPECTRUM.md) The two
   nonconjugate exceptional \(A_5\) stabilizers are Gassmann equivalent.  Their
   conjugate genus-zero Shabat polynomials over
   \(\mathbb Q(\sqrt{-11})\) have square discriminant, a corrected exact
   branch factorization, and derivative-unit common-target Keller charts.
   The source/target unit ranks are six and two, the derivative is nonzero
   modulo the target pullback lattice, and stable straightening is excluded.
   The direct \(C_5,C_6\) correspondences normalize to the degree-\(55\)
   \(A_4\) quotient of genus one and degree-\(66\) \(D_{10}\) quotient of
   genus two.  Their affine models have five/eight ordinary nodes and reduced
   conductor, with exact boundary pullback through both projections.  The
   first normalization is the conductor-\(121\) curve with \(j=-121\); exact
   traces above (23) prove it is not isogenous to the natural
   degree-twelve \(X_0(11)\) quotient, whose Weierstrass \(j\)-map is also
   explicit.  The second is bielliptic, with Jacobian \((2,2)\)-isogenous to
   the product of that conductor-\(11\) isogeny type and the CM
   conductor-\(121\) type.  Their positive-genus boundary-unit lattices have
   exact ranks \(3,14,17\); the natural quotient has an index-six support
   saturation, and the genus-two quotient has Smith cokernel
   \(\mathbb Z^3\oplus\mathbb Z/5\).  Both correspondence pullback images
   have rank ten and primitive free cokernels of ranks four and seven.  The
   two pulled-back derivative rows pass the weak integral ledger, while two
   additional mask characters cannot complete either full lattice.
   Projection exchange fixes the four residual \(C_5\) classes and splits
   the \(C_6\) quotient as five fixed classes plus one exchanged pair.
   Effective divisor bases realize every class; the \(C_5\) simple-pole
   lattice has index two, so one double pole is unavoidable.  The first
   one-class-per-monomial supports occur in normal degrees two and three.
   Normalization-module interpolation constructs exact formulas for all
   masks and proves an intrinsic infinity imbalance on \(C_6\).  Turning
   those rational normalization functions into an ambient polynomial output
   block satisfying the Jacobian and inverse-adjugate equations remains
   open; none has an absolute realization.
6. **Gassmann pairs inside a fixed \(G\).** Once one \(G\)-closure has been
   Kellerized, nonconjugate almost-conjugate subgroups should produce
   arithmetically indistinguishable inverse covers. This is the strongest
   test that the construction retains the original cover rather than only
   its generic degree.

## 8. Milestones and failure criteria

### Milestone I: one proper subgroup of \(S_N\)

Produce an absolute polynomial Keller map with certified generic inverse
monodromy \(G<S_N\), a complete connected rational fiber, and at least two
prescribed unramified Frobenius classes.

The group-only clause is now met by the decomposable map \(F\circ F\), but
the milestone as stated is not: this checkpoint does not certify the required
connected rational fiber and Frobenius data, and it is not an atomic family.
The collision-frame route remains the preferred structural target.

### Milestone II: a functorial relative theorem

Give checkable hypotheses on a monogenic \(G\)-cover and its divisor ledger
which imply a relative determinant-one Keller realization preserving the
Galois closure.

### Milestone III: local-global specialization

Prove a Hilbert-with-local-conditions statement inside the admissible Keller
parameter open, retaining the full arithmetic group \(G\), real signature,
and finitely many good-prime splitting conditions.

### Milestone IV: Sunada compatibility

For a Gassmann triple \((G,H_1,H_2)\), construct two nonisomorphic Keller
inverse covers with common \(G\)-closure and identical good-fiber zeta
functions; then determine when both descend to absolute affine-space maps.

The chart/relative clause is now met in two structurally different cases:
Fano point/line stabilizers in \(\operatorname{GL}_3(\mathbb F_2)\), and the
two exceptional \(A_5\) classes in \(\operatorname{PSL}_2(11)\).  Absolute
descent remains open in both.

A proposed construction fails the program if it only places the desired
extension on a proper target section while the ambient generic group is
\(S_N\), if it loses inverse sheets at the completion boundary, or if a
marking base change silently enlarges the arithmetic group.

## 9. Immediate research tasks

1. Write the root-engineered quadratic-gauge construction over a coefficient
   base and compute the generic function-field group before and after each
   added target coordinate.
2. Formulate a cover-preservation lemma for controlled-boundary suspensions:
   the reconstructed root algebra should be the original finite cover after
   restriction to the common open.
3. Search for a nonradial, nonautomorphic incidence map
   \(\alpha:\mathbb A^5\to\mathbb A^5\) satisfying
   \[
   \mathcal B(\operatorname{pr}_{1,2,3}\alpha)
   =u\,\mathcal B(P,Q,R)\det D\alpha,
   \]
   both adjugate divisibilities for the inverse two-mask blowdown, and
   reconstruction of the original \(A_4\) target field. The normalized
   automorphic assembly and every radial base `(SP,SQ,SR)` with arbitrary
   polynomial masks are excluded. Regular cyclic cubics and \(V_4\) quartics
   remain absolute no-go cases because their degree extensions are Galois.
4. Run fixed-support mixed-incidence searches for the quartic collision
   frame in its conductor-adapted parametrization
   \(m=-2a\tau+\sigma\), \(p=-a^2-2\tau^2+n\).  Feed two masks into at least
   one coefficient output, impose the modified adjugate divisibilities and a
   constant full Jacobian, then certify the degree-six \(S_4\) field by
   elimination.  Every surviving source must pass units, class group/UFD,
   ML, Derksen, Hodge--Deligne, and topological gates before an affine-space
   recognition proof is attempted.
5. Combine the finite-field witness search with conjugacy classes of a
   prescribed subgroup \(G\), rather than partitions of \(N\) sampled from
   \(S_N\).
6. Recast the Davenport point/line construction as the first complete
   relative instance of the cover-first theorem and isolate exactly which
   hypothesis fails at affine-space descent.
7. Separate three arithmetic outputs in every computation: connectedness of
   the degree-\(N\) fiber, its splitting-field group, and the geometric
   monodromy of the surrounding inverse cover.
8. The positive-genus pullback/cokernel and support pre-sieves are complete.
   The two
   projection images have rank ten and primitive cokernels
   \(\mathbb Z^4,\mathbb Z^7\); the order-five class on the genus-two curve
   is a prior divisor-class gate, whereas the order six on \(X_0(11)\) is an
   internal saturation defect among evident units.  A two-derivative ledger
   is feasible.  Effective residual bases exist with exchange signatures
   four fixed and five fixed plus one pair; their first separate-monomial
   supports have normal degrees two and three.  The next calculation should
   construct explicit Riemann--Roch representatives and solve the full
   Jacobian and inverse-adjugate system on those finite supports.
9. Audit the natural \(\operatorname{PSL}_2(q)\) actions first at
   \(q=7,11,13\).  Keep the rigid quotient genus, the field of the Nielsen
   component, and the affine-completion ledger separate from essential
   dimension; the latter constrains a versal Keller family, not one rigid
   cover.

The program succeeds when explicit Galois covers can be transported into
constant-Jacobian geometry without forgetting their generic monodromy, while
the existing full-fiber and adelic controls remain available.
