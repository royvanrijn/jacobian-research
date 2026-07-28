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

## 7. First pilot groups

The most informative order is not simply increasing group size.

1. **Nonregular cyclic and dihedral actions.** Kummer and Chebyshev-type
   presentations give explicit covers, but naturally introduce tori and
   quotient boundaries. Regular actions are excluded for absolute polynomial
   Keller maps by the Galois-case invertibility theorem.
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
   automorphic incidence embedding is excluded.
3. **Small imprimitive groups.** The actions of \(V_4\), \(D_4\), and wreath
   products test whether block systems survive marked-root completion rather
   than expanding to \(S_N\).
4. **\(\operatorname{GL}_3(\mathbb F_2)\) in degree seven.** This is the
   existing nonsymmetric benchmark. The relative and Cox realizations are
   known; absolute affine-space descent remains open.
5. **Gassmann pairs inside a fixed \(G\).** Once one \(G\)-closure has been
   Kellerized, nonconjugate almost-conjugate subgroups should produce
   arithmetically indistinguishable inverse covers. This is the strongest
   test that the construction retains the original cover rather than only
   its generic degree.

## 8. Milestones and failure criteria

### Milestone I: one proper subgroup of \(S_N\)

Produce an absolute polynomial Keller map with certified generic inverse
monodromy \(G<S_N\), a complete connected rational fiber, and at least two
prescribed unramified Frobenius classes.

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
3. Search for a nonautomorphic incidence map \(\alpha:\mathbb A^5\to
   \mathbb A^5\) satisfying
   \[
   \mathcal B(\operatorname{pr}_{1,2,3}\alpha)
   =u\,\mathcal B(P,Q,R)\det D\alpha,
   \]
   both adjugate divisibilities for the inverse two-mask blowdown, and
   reconstruction of the original \(A_4\) target field. The normalized
   automorphic assembly is excluded. Regular cyclic cubics and \(V_4\)
   quartics remain absolute no-go cases because their degree extensions are
   Galois.
4. Combine the finite-field witness search with conjugacy classes of a
   prescribed subgroup \(G\), rather than partitions of \(N\) sampled from
   \(S_N\).
5. Recast the Davenport point/line construction as the first complete
   relative instance of the cover-first theorem and isolate exactly which
   hypothesis fails at affine-space descent.
6. Separate three arithmetic outputs in every computation: connectedness of
   the degree-\(N\) fiber, its splitting-field group, and the geometric
   monodromy of the surrounding inverse cover.

The program succeeds when explicit Galois covers can be transported into
constant-Jacobian geometry without forgetting their generic monodromy, while
the existing full-fiber and adelic controls remain available.
