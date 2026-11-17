# Ten term lower bound for normalized degree seven Keller collisions

Over any characteristic-zero field, let F=(F1,F2,F3) satisfy

\[
F(0)=F(e_1)=0,\qquad JF(0)=I,\qquad \det JF=1.
\]

If every displayed coordinate has ordinary total degree at most seven,
then the nonlinear part contains **at least ten monomial occurrences**,
counted with coordinate multiplicity.

This strengthens the [eight-term exclusion](../sparse-eight/PROOF.md).
The proof combines necessary coordinate conditions with an exhaustive
support search. It imports the established plane degree-100 theorem through
the [constant-pairing reduction](../sparse-eight/COORDINATE_GATES.md).
It does not enumerate every raw determinant-balanced nine-term support.

## Necessary row conditions

Every completed support of a collision in degree at most ten must meet
these conditions:

1. Each output row has a nonlinear monomial containing its own variable.
   Otherwise Fi=xi+H(other variables). Straightening this elementary source
   coordinate produces plane fibres of degree at most D^2<=100. Moh's
   plane theorem makes every fibre injective, excluding the collision.
2. Each row has a nonlinear monomial involving another source variable.
   Otherwise Fi is either xi alone, already excluded in (1), or a nonlinear
   univariate polynomial. Its nonconstant derivative would divide det JF.
3. If every nonlinear term in row i contains xi and some has xi-exponent
   exactly one, the final row must acquire a term independent of xi.
   Without it, Fi is divisible by xi and its normal derivative on xi=0
   is nonconstant. Restricting the determinant to that hyperplane contradicts
   the polynomial unit identity.
4. A transverse row with exactly one pure-x term must acquire another
   pure-x term, since F2(e1)=F3(e1)=0 and all support coefficients are nonzero.

These conditions concern completed exact supports. Adding a missing monomial
can repair a partial support, so none is discarded merely because it
currently violates a row condition.

The only external input is the plane theorem used in (1).
[Moh's primary author account](https://www.math.purdue.edu/~ttm/jacobian)
states injectivity when both plane component degrees are at most 100.
The reduction preserves distinctness of any colliding points, applies over
every characteristic-zero coefficient field, and needs no global inverse
for a family of fibres.

## Strongly connected dependence

Draw an edge i to j, i!=j, when Fi contains a nonlinear monomial involving
xj. The dependence graph of a hypothetical collision must be strongly
connected.

If it were not, order its strongly connected components topologically.
The Jacobian becomes block triangular under simultaneous row and column
permutation. Its determinant is the product of the diagonal-block
determinants, so every factor must be a nonzero polynomial unit.
With three vertices, a nontrivial component partition has a singleton
component. Its diagonal block is the scalar derivative dFi/dxi.

Condition (1) makes that derivative nonconstant in characteristic zero:
distinct own-variable monomials have distinct nonzero derivatives. This
contradicts its being a unit. Therefore all three vertices form one
strongly connected component.

Equivalently every nonempty proper vertex subset has an outgoing edge.
If a partial support has none, any final support containing it must acquire
an atom producing such an edge. The independent implementation instead
tests reachability from each vertex and repairs a proper reachability set.

## Complete gated closure

The collision-axis argument in the
[earlier census](../sparse-seven/PROOF.md) requires a pure-x term in F1 and,
after the simultaneous transverse swap, at least two in F2. The 90 seeds
are the smallest such exponents (p,q,r), with

\[
2\le p\le7,\qquad2\le q<r\le7.
\]

The universe is all 348 coordinate-labelled nonlinear atoms of degree two
through seven. No mixed monomial, weighted form or coefficient box is
assumed.

For each partial support S, the [production compiler](gated_search.py):

- closes nonconstant singleton Jacobian determinant buckets by adding every
  possible missing contributing triple;
- repairs a single transverse axis term by every possible second axis term;
- at a determinant-balanced support, repairs one failed row or dependence
  condition by **every** atom that could repair it;
- retains a support only after all necessary gates pass, and, if below the
  size bound, permits arbitrary additional atoms.

Completeness follows by induction inside a hypothetical final support T.
Every chosen singleton repair, axis repair, row repair or missing-edge repair
has a choice entirely contained in T. At an already valid proper subset,
arbitrary atom addition includes an atom of T minus S. The size bound makes
the process finite. The seed and smallest-exponent filters preserve at least
one orientation of every candidate.

An additional budget bound is safe: requirements confined to different
output rows need different new atoms. Within a transverse row, an axis-count-one
repair is disjoint from a missing-own-variable repair. The code counts both.
No graph requirement is added to this bound, avoiding double counting an
atom that could repair several edges.

All 90 shards terminate with empty frontiers after a total of **6,040,767**
partial-support visits and retain **zero supports** through size nine.

## Independent traversal and positive control

The [independent standard-library traversal](verify_coverage.py) imports no
production search code. It reconstructs the monomial universe and determinant
contributions using a separate Leibniz expansion, starts all seeds in one
global traversal without smallest-exponent shard filters, and uses:

- direct row predicates rather than production row masks;
- graph reachability rather than production cut masks;
- per-row intersections of repair-atom sets rather than the production
  lower-bound case split;
- a different tie and traversal order.

It finishes after **7,228,226** visits with an empty frontier and no retained
support. Exact candidate-set equality with the union of all production
shards is required. Interrupted passes retain their actual DFS frontier;
only terminal execution produces the passing comparison receipt.

The [controls](verify_controls.py) independently normalize the actual
foundational counterexample to the same collision frame. Its 49-term support
passes every new gate. A complete production branch chain stays inside that
support from its seed to the full counterexample. This guards against an
implementation that spuriously rejects every support.

A separate degree-three, support-at-most-eight regression compares the gated
search with all 173 supports from the original unseeded determinant closure.
The filtered candidate sets agree. Actual checkpoint resumption from a
one-node incomplete pass is tested in both degree-three shards.

There is no residual coefficient ideal in the nine-term proof: zero cases
reach rational algebra. The manifest's empty algebra list is explicit;
its vacuous all-unit flag is not the reason for the exclusion.

## Reproduction and boundaries

A scoped integrity audit, with optional fresh positive-control replay:

~~~sh
.venv/bin/python research/jacobian-followup/sparse-nine/verify.py --replay-controls
~~~

To perform a fresh independent enumeration, choose an unused checkpoint path:

~~~sh
.venv/bin/python research/jacobian-followup/sparse-nine/verify_coverage.py \
  --seconds 240 --checkpoint tmp/nine-independent-replay.pkl.gz
~~~

If it is incomplete, repeat exactly that command to resume its saved
frontier. A completed checkpoint can be compared again, but that comparison
does not constitute a fresh enumeration. The original independent execution
was resumed twice and took about 555 seconds, excluding compilation and
checkpoint loading.

To regenerate production evidence in a separate directory:

~~~sh
.venv/bin/python research/jacobian-followup/sparse-nine/run_all.py \
  --output tmp/nine-production-replay
~~~

Do not overwrite the predecessor sources or reinterpret old raw-support
counts as counts for this gated search. The proof is a finite combinatorial
exclusion using KPAIR1 and Moh. It is not proof-assistant formalization,
external review or an unrestricted support theorem.

The [49-term normalized example](NORMALIZED_WITNESS.md) gives an upper bound
in this same raw-degree-at-most-seven frame. The later
[48-term transverse completion](../collision-frame-search/PROOF.md) improves
the current normalized interval to [10,48]. Its invariant degree profile is (7,6,4), so the example
does not attain a profile strictly below that flag. The cardinality-unbounded
lower-profile census remains open.

