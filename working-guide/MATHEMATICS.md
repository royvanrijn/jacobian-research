# Mathematical validity

A result is ready to become a paper claim when its statement, proof and dependencies agree, and a reader can check the argument. Fluent prose, an agent's confidence and a passing script do not establish this agreement. Use this procedure for each new result and whenever a proof or hypothesis changes.

## Fix the statement before evaluating the proof

Write the exact claim with its objects, domains and quantifiers. Record the base field or coefficient ring, characteristic, dimension, regularity, degree or support bounds, and any exceptional cases. State whether parameters are fixed or vary. Distinguish existence from an effective construction, generic validity from validity everywhere, and a necessary condition from a sufficient one.

Check the order of quantifiers. A bound that depends on an individual input does not give one bound for an entire family. If a statement says that for every object there exists a cutoff, identify which data the cutoff may depend on. Prove any stronger uniform assertion separately.

Use a small example to confirm that each definition has the intended meaning. Check empty and zero cases, constants, small dimensions and degrees, singular or reducible objects, repeated factors, and changes of characteristic when relevant. State conventions that affect the claim rather than leaving them to a checker implementation.

Classify the statement before drafting its title or abstract:

| Kind of statement | What must accompany it |
|---|---|
| Proved result | A complete argument, including every indispensable imported result. |
| Conditional result | Explicit assumptions and a proof of the implication under those assumptions. |
| Computational result | Exact search domain, algorithm, coverage argument and reproducible evidence. |
| Conjecture or heuristic | A clear label, motivation and limits of the available evidence. |
| Counterexample | An explicit admissible object and verification of every hypothesis and the failed conclusion. |

A result can fit more than one row. For example, an exact computational theorem may depend on a stated finite reduction. Keep every condition in the theorem and in summaries of it.

## Audit the proof as a dependency graph

For each principal claim, keep a compact internal record of the exact statement, proof location, dependencies, evidence type, reviewer and unresolved obligations. Reuse the repository's current result records instead of creating a parallel status system. This record belongs in supporting documentation; the paper itself must present the mathematical argument.

The executable graph is part of [MATH_STATUS.json](../research/MATH_STATUS.json).
Use the [dependency guide](../research/DEPENDENCY_GUIDE.md) to record exact nodes,
proof-backed implications and source pins. Its [generated diagrams](../research/DEPENDENCIES.md)
and [scope audit](../research/SCOPES.md) show current routes and gaps. A graph pass
checks the declared structure and evidence freshness; it does not validate the
mathematics or establish complete dependency coverage.

Start at the main theorem and work backwards. For every lemma or external theorem used, check its precise hypotheses, conclusion and application. Cite the exact theorem, proposition or page and version when numbering can change. A matching title, abstract or search snippet cannot establish that a source proves the needed fact. If the source is unavailable, the dependency remains unverified until the proof is supplied or checked.

Look for circularity: a lemma must not use the main theorem, an equivalent unresolved conjecture or a later result that depends on it. Separate a new reduction from a solution of the reduced problem. An equivalence theorem does not resolve either side without an independent argument for one side.

Expand every step that carries mathematical weight. For an asserted extension, descent, limit, specialization or gluing step, identify the theorem that permits it and verify its assumptions. Check that denominators are units on the claimed domain, coordinate changes are admissible, constructions preserve the required category, and transformations preserve both directions of any claimed equivalence. For infinite sums or limiting arguments, justify convergence or formal validity and any interchange of operations.

A proof gap is an open obligation. Repair it, weaken the statement to what is proved, state the result conditionally, or withdraw the claim. Do not conceal the gap with phrases such as "standard," "clearly" or "by routine arguments." The same rule applies when an agent generated the missing step confidently.

## Challenge scope changes explicitly

Common invalid transitions include finite to unrestricted, local to global, generic to every point, formal to polynomial, numerical to exact, and equality after a quotient or localization to equality in the original object. At each transition, write the additional argument required.

**Finite testing.** For an integer M >= 0, the polynomial

`p_M(t) = product_{j=0}^M (t - j)`

vanishes at every integer from 0 through M and is nonzero over the rationals. Checking a longer initial range cannot prove a general identity without a degree bound or another completeness argument. Finite testing can prove an identity of known degree at most d from d + 1 distinct roots over a field; the bound and distinctness are part of the proof.

**Formal versus polynomial inversion.** Over a field of characteristic zero, `f(x) = x + x^2` has a compositional inverse in the formal power series ring because its linear coefficient is a unit. It has no polynomial compositional inverse: if a nonconstant polynomial g satisfied `f(g(x)) = x`, then `deg(f composed with g) = 2 deg(g)`, which cannot be 1. A formal inverse construction therefore needs an additional termination argument before it yields a polynomial inverse.

**Localization.** In `k[x, x^{-1}]`, multiplication by x is invertible. In `k[x]`, it is not surjective because 1 is not divisible by x. A construction after inverting a factor need not descend to the original ring. Identify the exact descent condition and prove it.

These examples diagnose transitions; they do not replace checking the particular proof. Seek examples that violate the claimed conclusion while satisfying all its hypotheses. A failed attempt at a counterexample is evidence about that attempt, not a proof of the theorem.

### Check simultaneous elimination on the same fibre

When several residuals are imposed after a finite elimination, verify that
the certificate forces them to vanish at the same geometric point. Separate
norms or pairwise resultants can vanish using different witnesses. For a
finite free algebra \(A/R\), the coefficients of the single norm
\(N_{A/R}(\sum_j T^j\rho_j)\) test for one common geometric zero of all the
residuals. Each local factor contributes its length to the norm. The
[conductor-18 calculation](../research/foundations/extended-geometry/CM_NEWTON_SEVENTEEN_ROW_LOCUS.md#a-single-norm-for-the-common-root-condition)
gives the proof and a rank-two example. A complete common-root certificate can
also use ideal membership or subresultants. Keep each leading-coefficient
chart, denominator and excluded slice explicit, and retain nilpotents when
making a statement about scheme length.

## Match computational evidence to the claim

Read the checker, its inputs and its output semantics before relying on a success message. Determine whether it verifies a certificate, independently derives an identity, searches a bounded space, samples examples, or merely checks file consistency.

Exact arithmetic can certify an exact equality for the represented inputs. Symbolic algebra needs the correct coefficient domain and assumptions; simplification under unstated nonzero or branch assumptions can change the problem. Modular checks need a reconstruction or lifting argument to establish a characteristic-zero claim. Floating point residuals and plots support exploration; a proof using numerical computation needs rigorous error bounds, such as certified interval enclosures, and an argument connecting those bounds to the theorem.

A computational proof must account for the entire stated domain. Explain the finite reduction, enumeration completeness, symmetry reductions, exceptional cases and why the algorithm terminates. Testing representatives requires a proof that every admissible case has a representative. Searching all cases under a resource bound proves only the recorded bounded assertion unless a theorem makes that bound sufficient.

For each indispensable computation:

1. Record the mathematical input, coefficient domain, algorithm, pinned dependencies, command and expected invariant. Retain the actual input or certificate, rather than only a checksum or success flag.
2. Check at least one small example by hand or an independently implemented method. Use known positive cases and deliberately invalid inputs or corrupted certificates to check that relevant failures are rejected.
3. Check that every invariant needed for the conclusion is verified. A checker that only recomputes a displayed formula cannot detect an omitted hypothesis or an incomplete case list.
4. Rerun after relevant source or input changes and record the actual outcome. A timeout, crash, skipped case or incomplete run is not a pass. Distinguish an old receipt from a run against current bytes.

Use a simpler certificate verifier when possible. It should check the claimed mathematical properties without trusting the generator's conclusion. A second agent reviewing the same derivation can help find errors, but shared models, code and assumptions can produce correlated mistakes.

Describe verification precisely. An author's rerun is an internal replay. A separately performed replay checks reproducibility. An independent mathematical review examines whether the computation and proof establish the stated theorem. Record who performed each check and what was examined; none of these labels implies the others.

## When to extend a search

Before raising a degree cap, expansion order, precision, or search range,
identify the unresolved proof step and what the next computation could establish.
Choose a run that serves at least one of these purposes.

1. **Close a delimited obligation.** Use a proved reduction, exhaustive cutoff,
   or specific certificate requirement to identify what remains to be checked.
2. **Find a mechanism that extends.** Use selected examples to isolate a
   recurrence, invariant, compatibility obstruction, or uniform estimate.
   Identify where the same argument would have to work beyond the tested range.
3. **Test a proposed argument.** State what a pass or a failure would establish
   for a particular hypothesis or proof step. Choose the smallest case that can
   decide it and retain the exact witness or certificate.

After a batch, record the verified range and any change to the proof gap. If
another run would only enlarge the range while leaving the same gap open,
revisit the argument before extending the search.

### Applications in this repository

- **Weyl lift degree and order.** The canonically graded scalar Long scheme at
  spatial cap \(27\) fails through \(h^7\). At a larger cap, identify which
  new coefficient directions could remove that residual obstruction. Use a
  higher jet to test an obstruction uniform in the cap or a recursion with one
  spatial bound valid at every order. [Cap-27 obstruction](../research/dixmier/dixmier-formal-lifting/GRADED_LONG_CAP27_FULL_OBSTRUCTION.md)
- **Factorial subsum precision.** For \(f(1,5)=\infty\), the missing step is
  to handle valuation jumps of unbounded size. Use another precision check to
  test a proposed growing lookback rule or a carry invariant that reaches the
  target residue \(-1\). Select a jump where that rule could fail.
  [Factorial subsums and valuation jumps](../research/factorial/ERDOS_404_FACTORIAL_SUMS.md)
- **Moment prefixes.** For ternary Gaussian cubics, sixteen moments suffice
  locally along the exceptional orbit. To decide the global sixteen-moment
  question, identify residual components away from that orbit. On each
  component, seek a purity certificate or a point satisfying the first sixteen
  equations with a certified nonzero later moment. Compute later moments to
  decide a specified component. [Local rigidity and global scope](../research/moments/vanishing/GAUSSIAN_EXCEPTIONAL_SIXTEEN_MOMENTS.md)
- **Finite arithmetic coverage.** The thirds-balanced Eisenstein problem has
  a proved reduction to all-inert targets in a fixed finite interval. An exact
  check covering that residual set would settle the stated exception list.
  [Exact reduction](../research/arithmetic/EISENSTEIN_SIFTED_ANGULAR_DISTRIBUTION.md#exact-balanced-exceptions-prime-families-and-reduced-targets)

## Directional-conductor computations

Prioritize structural results in the plane Keller conductor work. Reserve
searches at higher conductors for computations with one of the purposes below.

1. **Extend an infinite family.** The [three-to-two Laurent exclusion](../research/hessian/hc4-triangular-graphs/CONDUCTOR_FOURTEEN_CUBIC_ORDER_ALTERNATE_EMBEDDINGS.md#the-general-three-to-two-laurent-obstruction)
   covers the sequences \((3m,2m,r)\) with \(2\le r\le m\) and
   \(\gcd(m,r)=1\). For \(r>m\), identify which additional Laurent terms
   can enter the polynomial Jacobian identity and what constraints they
   satisfy. This is the obstruction to extending the existing argument.
2. **Extract a uniform lemma from the case proofs.** Compare the
   [quartic–septic parameter argument](../research/hessian/hc4-triangular-graphs/CONDUCTOR_EIGHTEEN_QUARTIC_SEPTIC_EXCLUSION.md)
   and the [singular tangent determinant argument](../research/hessian/hc4-triangular-graphs/HYPERELLIPTIC_CONDUCTORS_AND_SINGULAR_TANGENT_DETERMINANT.md).
   Identify hypotheses on the parameter invariants, finite poles and
   singular-fibre lifts that make the mechanism work independently of
   conductor degree. Include the behavior at finite degree drops in the
   proposed statement.
3. **Use individual cases to test the lemma.** Before calculating a new
   conductor case, state which proposed hypothesis or proof step it tests.
   Retain the case if it exposes a failure of the generalization or supplies
   a specific ingredient needed in its proof.

## State the exact scope of formal verification

Pin the prover, library versions and build command. Build the target theorem and inspect its transitive assumptions. In Lean, inspect `#print axioms` for the relevant declaration and dependencies, check for `sorry` and `sorryAx`, and identify custom axioms. Distinguish ordinary foundational axioms from unproved mathematical inputs; report each material input with its exact statement.

Check the bridge between the formal theorem and the paper: the definitions, domains, quantifiers, hypotheses and conclusion must match. Formal verification of an algebraic identity, a finite certificate or a conditional theorem does not by itself verify an entire paper. Explain any interpretation or reduction that remains outside the prover. A successful build proves the checked statement under its assumptions; it does not establish that the assumptions model the intended claim.

Do not edit an assumption into a definition or introduce an axiom merely to make a build pass while describing the result as proved. If an essential mathematical input has no checked proof or source, state the mathematical result conditionally in the theorem, abstract and verification statement. If the paper proves or correctly sources that input but the formal development assumes it, the paper's theorem can be unconditional while the formalization remains conditional on the imported input. Describe those two scopes separately.

## Review and correct without inflating assurance

Ask a reviewer to reconstruct the main argument from its definitions and cited results, identify its weakest step, and challenge the transitions above. Review the claim separately from the prose: a readable proof can still be false, and a correct computation can support an incorrectly generalized conclusion. Record substantive objections and their resolutions in the existing support record.

Before calling a result ready for review, verify that every principal claim has a complete proof or an explicit conditional status; all computational claims stay within checked coverage; all imported results have been checked; and formalization claims identify their actual assumptions and scope. External review or submission requirements remain separate from these internal checks.

When an error changes a statement, trace every dependent result. Update the proof, abstract, introduction, conclusions, current result catalogue, manuscript metadata and verification description as applicable. Rerun affected checks and mark dependent claims unresolved until repaired. Preserve exact deposited editions; put corrected manuscripts in the established revision area and use the publication's correction process when needed. Never silently replace a deposited file or turn an internal correction into a claim of external review.
