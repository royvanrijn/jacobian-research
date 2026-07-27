# Local Q2 certificates for arithmetic Keller fibers

## Status

Research direction motivated by David Roe and David Turturean's 2026 explicit
presentation of the absolute Galois group of `Q_2` by four generators and two
relations, with two independent Lean 4 formalizations.

This does **not** alter the polynomial Keller constructions or their Jacobian
proofs.  It supplies a new arithmetic description and possible certificate
layer for the finite etale fiber algebras already realized by those maps.

Primary source:

- David Roe and David Turturean, *A Presentation of the Absolute Galois Group
  of Q_2*, project page: <https://roed314.github.io/gq2/>.

## 1. Exact bridge

For a field `K`, finite etale `K`-algebras are contravariantly equivalent to
finite continuous `G_K`-sets, where

\[
 G_K=\operatorname{Gal}(K^{\mathrm{sep}}/K).
\]

If a Keller map over `Q_2` has a complete fiber with coordinate algebra `A`,
then its geometric points form the finite continuous set

\[
 X_A=\operatorname{Hom}_{\mathbb Q_2\text{-alg}}
      (A,\overline{\mathbb Q}_2)
\]

with its natural `G_{Q_2}` action.  Conversely this action determines `A` up
to `Q_2`-algebra isomorphism.

The Roe--Turturean presentation therefore turns the local arithmetic type of a
Keller fiber into finite action data for four named generators satisfying two
relations.  At the level of a degree-`N` fiber, this can be represented by four
permutations in `S_N`, together with checks of the two relators and any marking
data required by the presentation.

The central point is:

> The finite-etale realization theorem constructs the fiber; the explicit
> presentation of `G_{Q_2}` can provide a compact certificate of its complete
> local Galois structure.

## 2. What this can be used for

### 2.1 Local fiber fingerprints

For a rational fiber algebra `A/Q`, base change to `Q_2` and record the
resulting finite `G_{Q_2}`-set.  The orbit decomposition gives the factor
and residue-degree decomposition over `Q_2`; stabilizers retain the full local
extension type rather than only the degrees.

For two global fibers, unequal `Q_2` action certificates immediately prove
that the global `Q`-algebras are nonisomorphic.  Equal local certificates do
not by themselves prove global isomorphism.

### 2.2 Small independently checkable certificates

A candidate certificate format is:

1. a squarefree polynomial `P in Q[T]` defining the fiber;
2. a certified factorization/decomposition of `P` over `Q_2`;
3. four permutations describing the action of the presented generators on
   the geometric roots;
4. exact verification of the two group relations;
5. a comparison certificate identifying this action with the etale algebra
   `Q_2[T]/(P)`.

Items 3--4 are finite combinatorics.  Item 5 is the substantive local-field
bridge and should not be silently inferred from matching orbit sizes.

### 2.3 Search and engineering

The presentation suggests reversing the usual workflow:

1. enumerate small transitive permutation actions satisfying the two
   relations and the presentation's marking constraints;
2. interpret them as candidate finite extensions, or products of extensions,
   of `Q_2`;
3. construct a squarefree polynomial with that local algebra;
4. feed the polynomial into the universal arithmetic-fiber Keller
   construction.

This would produce Keller fibers selected by prescribed local Galois action,
not merely by degree or factorization pattern.

For global `Q`-fibers one can ask for simultaneous local specifications at a
finite set of primes and then use approximation/Hilbert irreducibility to seek
a global polynomial meeting them.  The `Q_2` presentation makes the dyadic
local condition unusually explicit.

### 2.4 Formalization reuse

The two Lean formalizations may provide reusable definitions and lemmas for
profinite presentations, finite continuous actions, and marked generators.
The immediate formal target should remain modest: connect a finite action of
the presented group to a finite etale algebra only after identifying exactly
which local-field and Galois-category infrastructure is available.

This is a downstream formalization project.  It is not required for the
current finite-etale Keller theorem, whose algebraic realization statement is
field-uniform.

## 3. What this does not currently solve

- It does not simplify the Keller map formulas.
- It does not strengthen the constant-Jacobian calculation.
- It does not address injectivity, properness, or the plane Jacobian problem.
- It does not classify global finite etale `Q`-algebras.
- A tuple of permutations satisfying the relators is not automatically a
  certificate for a particular polynomial without the comparison step.
- The explicit presentation describes `G_{Q_2}`; it does not imply that every
  abstract finite permutation representation chosen without continuity or
  marking checks is arithmetically realizable.

## 4. Open problems

### Q2-1. Certificate schema

Extract the exact four generators, two relators, topology, and marking
conventions from Roe--Turturean and define a versioned finite certificate
schema for actions on at most `N` points.

**Deliverable:** a verifier that checks the finite permutation relations and
reports orbit/stabilizer invariants.

### Q2-2. Polynomial-to-action compiler

Given a squarefree `P in Q[T]`, compute and certify the finite
`G_{Q_2}`-set attached to `Q_2[T]/(P)` in the Roe--Turturean generators.

This is substantially stronger than factoring `P` over `Q_2`; it requires
recovering the action of the named generators or a rigorously equivalent
marked description.

### Q2-3. Action-to-polynomial realization

Given a certified finite action of the presented `G_{Q_2}`, construct an
explicit squarefree polynomial `P in Q_2[T]` whose etale algebra has that
action.  Seek integral or rational approximations suitable for input to the
existing Keller-fiber formulas.

### Q2-4. Prescribed local Keller fibers

Prove an explicit theorem of the following shape:

> Every finite continuous `G_{Q_2}`-set, supplied in presentation-certificate
> form, occurs as the geometric fiber action of an explicit Keller map over
> `Q_2`.

Abstractly this follows from the finite-etale realization theorem plus the
Galois-category equivalence.  The open content is an effective compiler with
checkable output.

### Q2-5. Global fibers with dyadic specification

Construct infinite families of `Q`-defined Keller maps/fibers whose base
change to `Q_2` has a prescribed certified action while retaining chosen
behavior at other primes and over `R`.

A first test should use the existing quartic and quintic arithmetic-zoo
polynomials, compute their exact dyadic action certificates, and determine
whether the local data separates examples not already separated by elementary
factorization and discriminant invariants.

### Q2-6. Lean bridge

Audit the two Roe--Turturean formalizations and determine the smallest import
surface needed to state:

\[
 \text{finite action certificate}
 \Longrightarrow
 \text{finite continuous }G_{\mathbb Q_2}\text{-set}
 \Longleftrightarrow
 \text{finite etale }\mathbb Q_2\text{-algebra}.
\]

Do not begin by formalizing all local Galois theory.  First identify which
part can be isolated as a finite-action interface and which theorems must
remain trusted external inputs.

## 5. Priority assessment

This branch is valuable for arithmetic packaging, explicit examples, and
credibility through independently checkable certificates.  It is orthogonal
to the current highest-priority obstruction problems in
`OPEN_PROBLEMS_FOR_MAP_EXTENSIONS.md`.

Recommended order:

1. Q2-1: extract and implement the finite relation verifier;
2. compute coarse `Q_2` decompositions for the existing arithmetic zoo;
3. assess whether named-generator actions can be extracted effectively;
4. only then attempt the Lean Galois-category bridge.

The likely near-term publishable object is a **local arithmetic certificate
attached to an explicit Keller fiber**, not a new Jacobian counterexample or a
simplification of the map construction.
