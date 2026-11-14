# Nine term lower bound for normalized degree seven Keller collisions

Over any characteristic-zero field, suppose F=(F1,F2,F3) satisfies

\[
F(0)=F(e_1)=0,\qquad JF(0)=I,\qquad \det JF=1,
\]

and every displayed coordinate has ordinary total degree at most seven.
Then the nonlinear part of F contains **at least nine monomial occurrences**,
counted with coordinate multiplicity.

This closes eight-term attainment in the same frame as the
[seven-term predecessor](../sparse-seven/PROOF.md). The conclusion is about
the normalized frame; affine changes can change support cardinality.
No whole dense profile above maximum degree two is excluded.

## Complete finite cover

The [predecessor proof](../sparse-seven/PROOF.md) gives the degree-independent
axis necessity: F1 has a pure-x term and, after the simultaneous transverse
swap, F2 has at least two. Let p be the smallest pure-x exponent of F1 and
q<r the two smallest in F2. The 90 triples

\[
2\le p\le7,\qquad2\le q<r\le7
\]

cover every candidate. All 348 nonlinear coordinate-labelled atoms of
degrees two through seven remain available, subject to the specified
smallest-axis-exponent rules. Mixed nonlinear terms, additional axis terms
and dependence on all three variables remain unrestricted.

The checkpointed [production search](../sparse-seven/search.py), run with
support bound eight, uses the same exhaustive determinant-singleton repair,
transverse collision repair and balanced-superset branching proved complete
in the predecessor. That induction does not depend on a bound of seven.
Every shard is terminal with an empty frontier.

The complete counts are:

| Nonlinear occurrences | Oriented supports | Residual swap orbits |
|---:|---:|---:|
| at most five | 0 | 0 |
| six | 450 | 450 |
| seven | 1,876 | 1,876 |
| eight | 20,857 | 17,482 |
| total | 23,183 | 19,808 |

Here oriented means that F2 carries at least two axis terms. A support with
axis terms in both transverse rows can contribute both members of a swap
orbit. Therefore oriented counts must not simply be doubled to obtain
labelled counts.

The [independent standard-library enumerator](verify_coverage.py) uses its own
monomial universe and Leibniz determinant expansion. It starts all seeds in
one global traversal without the production shard filters. It visits
1,947,529 supports and reproduces the exact union of all production supports,
including all 20,857 oriented size-eight supports. Set equality is checked
before quotienting, not inferred from the aggregate counts.

## Rational coefficient exclusion

For each residual support, the exact ideal contains every nonconstant
Jacobian determinant coefficient, all collision equations, and

\[
\rho\prod_{a\in S}c_a-1.
\]

The last equation enforces every declared coefficient to be nonzero. The
constant determinant coefficient is already one from JF(0)=I.

Every one of the **19,808 ideals** has rational Groebner basis {1} in both
Singular and SymPy. This includes all 17,482 new size-eight orbits and all
2,326 smaller regression orbits. Hence no coefficient assignment exists over
any characteristic-zero extension field. All unretained supports were
excluded by the necessary axis or singleton conditions.

The [verifier](verify.py) checks exact shard and orbit coverage and ordered
support IDs in both solver ledgers. It independently builds the derivative
entries and expands the six signed determinant permutations, comparing the
entire coefficient ideals, including collision and saturation equations.
The independent input compiler is retained unchanged in the predecessor's
[verifier](../sparse-seven/verify.py). Both rational backends were rerun after
this input comparison. The [receipt](validation.json) pins all relevant
sources and results.

This proof uses no plane degree theorem. It is an exhaustive exact finite
exclusion, without a formal proof-assistant certificate or external review.

## Structural explanation

The same finite support cover also has a short structural classification:

| Size-eight orbit class, in the listed priority order | Orbits |
|---|---:|
| A coordinate is nonlinear and univariate in its own source variable | 17,435 |
| A coordinate is its own source variable plus a polynomial in the others | 46 |
| The remaining hyperplane-normal-factor obstruction | 1 |

A nonlinear univariate coordinate would contribute its nonconstant derivative
as a factor of det JF, which cannot equal one.

The elementary-shear coordinate is straightened by a triangular source
automorphism. Each resulting plane fibre has degree at most 49 and is
injective by Moh's degree-100 plane theorem. This alternative structural
route uses that external import, unlike the rational coefficient proof above.

The sole remaining support is

\[
\begin{aligned}
F_1&=x+a z^2+b xz+c x^3,\\
F_2&=y+dxy,\\
F_3&=z+e z^2+f xz+g x^2+h x^3,
\end{aligned}
\]

with all eight coefficients nonzero. On y=0 the second Jacobian row is
(0,1+dx,0), so the restricted determinant has the nonconstant factor 1+dx.
It therefore cannot be one. The [classification replay](verify_structure.py)
retains all 48 non-univariate supports, including the one seven-term
regression, and verifies the category membership.

The general [coordinate gates](COORDINATE_GATES.md) explain these obstructions
and give a stronger affine-invariant test for future candidates.

## Reproduction and remaining work

From the repository root:

~~~sh
.venv/bin/python research/jacobian-followup/sparse-eight/verify.py \
  --replay-coverage --replay-algebra
.venv/bin/python research/jacobian-followup/sparse-eight/verify_coordinate_gates.py
~~~

The independent coverage replay takes roughly 90 seconds in the recorded
environment; rational input compilation and replay take roughly 70 seconds.
Both are actual replays, not hash-only audits. Checkpoint existence is never
treated as proof of terminal search status.

Regenerate production evidence in a separate directory with:

~~~sh
.venv/bin/python research/jacobian-followup/sparse-seven/run_all.py \
  --bound 8 --output tmp/sparse-eight-replay --independent
~~~

The schema name in production output is inherited from the seven-term
compiler; the explicit degree=7 and bound=8 parameters control this run.
The previous source and certificates are unchanged.

The first unexcluded normalized support size is nine. No nine-term
counterexample or attained absolute support minimum is claimed.

