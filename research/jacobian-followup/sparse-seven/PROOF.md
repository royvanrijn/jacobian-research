# Eight term lower bound for normalized Keller collisions

Let k be a characteristic-zero field and let
F=(F1,F2,F3) in k[x,y,z]^3 satisfy

\[
F(0)=F(1,0,0)=0,\qquad JF(0)=I,\qquad \det JF=1.
\]

Assume every displayed coordinate has ordinary total degree at most seven.
Then the nonlinear part of F contains at least eight monomial occurrences,
counted separately in each coordinate.

This extends the [six-term exclusion](../../archive/non-elliptic/extended-geometry/GLOBAL_LOW_DEGREE_SUPPORT_CENSUS.md).
It imposes no weighted, equivariant, transverse-linear, coefficient-box or
invariant-degree-profile restriction. The support bound refers to this
normalized frame; affine normalization can change support cardinality.
All 74 invariant profiles below (7,6,4) are covered by their raw degree bound,
but no entire higher-degree dense profile is excluded by this theorem.

## Collision axis reduction

Write h(t)=F(t,0,0) and v=h'. Since JF(t,0,0) has determinant one,
its first column v is unimodular over k[t].

At least one nonlinear pure-x term occurs in F1, because h1(1)=0 while
h1'(0)=1. Neither F2 nor F3 can contain exactly one pure-x term:
its coefficient would have to vanish at the collision point.
If both have no pure-x terms, then v=(h1',0,0). Unimodularity would
force h1' to be a unit. Its value at zero gives h1'=1 and hence h1=t,
contradicting h1(1)=0.

Thus, after simultaneously exchanging y with z and F2 with F3 if needed,
F2 contains at least two pure-x terms. There are at least three pure-axis
occurrences in total. This uses only the elementary part of the
[axis theorem](../../archive/non-elliptic/extended-geometry/COLLISION_AXIS_UNIMODULAR_FRONTEND.md);
no resultant criterion or choice of an SL3 completion is required.

Let p be the smallest pure-x exponent in F1, and q<r the two smallest in
F2. Every oriented candidate belongs to exactly one of

\[
(p,q,r),\qquad 2\le p\le7,\quad 2\le q<r\le7.
\]

There are 6 times 15 = 90 shards. Each starts with x^p in F1 and x^q,x^r
in F2. Further pure-axis terms in any row remain allowed, subject only to
these minimal-exponent definitions. Every mixed monomial of degrees two
through seven is allowed in every row. In particular, no assumption about
linear dependence on y or z enters the search.

## Finite determinant closure

There are 116 nonlinear monomials in three variables of degree at most
seven, and 348 coordinate-labelled support atoms. The fixed linear terms
are x,y,z. For exponents alpha,beta,gamma selected in the three rows,
their contribution to det JF is

\[
\det(\alpha,\beta,\gamma)\,
c_{1,\alpha}c_{2,\beta}c_{3,\gamma}\,
x^{\alpha_1+\beta_1+\gamma_1-1}
y^{\alpha_2+\beta_2+\gamma_2-1}
z^{\alpha_3+\beta_3+\gamma_3-1},
\]

with coefficient one for a fixed linear term. Nonzero contributions with
the same source exponent form a determinant bucket. All nonconstant
buckets must vanish. A bucket with exactly one nonzero contribution cannot
vanish on an exact coefficient support.

For each partial support S containing its axis seed:

1. If its size exceeds seven, discard it. If a transverse row has exactly
   one pure-x term, at least one more term in that row will be needed.
   Discard S if this lower bound already exceeds seven.
2. If S has a singleton determinant bucket, branch over **every** universe
   triple contributing nontrivially to that bucket, adding its missing
   support atoms whenever the new size is at most seven and the shard's
   exponent rules allow them. A bucket with no possible repair closes
   the branch. Choosing the bucket with fewest repairs changes only order.
3. Otherwise, repair a transverse row with a single pure-x term by branching
   over every allowed second pure-x term in that row.
4. If no repair is needed, retain S as balanced. If fewer than seven terms
   occur, branch over every possible additional atom. This last step retains
   balanced supersets rather than stopping at a balanced proper subset.

Completeness follows by induction inside any hypothetical final support T.
Its seed is included. Every singleton in a partial S contained in T must
have a second contributing triple in T, so at least one repair remains
inside T. The collision repair likewise has a choice inside T.
At a balanced proper subset, the arbitrary-addition step includes an atom
of T minus S. The finite size bound ensures eventual arrival at T.
Memoizing a visited support cannot change its descendants.

The [90 completed shard ledgers](results/manifest.json) contain:

| Nonlinear occurrences | Oriented balanced supports | Residual swap orbits |
|---:|---:|---:|
| at most five | 0 | 0 |
| six | 450 | 450 |
| seven | 1,876 | 1,876 |

Every balanced support has axis terms in exactly one transverse row, so its
swap supplies the other orientation. No support is fixed by that swap.
The raw labelled counts including both orientations are 900 and 3,752.

The standard-library-only [independent enumerator](verify_coverage.py)
reconstructs all monomials and determinant contributions using a separate
Leibniz expansion. It starts all 90 seeds in one global traversal, without
the production minimal-exponent filters. Its 266,295 visited supports yield
exactly the same 2,326 oriented balanced supports, by set equality, not just
matching totals. A separate [small-degree regression](regression.json)
compares the new shards against the original unseeded enumerator and tests
resumption from an explicitly incomplete checkpoint.

## Exact coefficient exclusion

For every residual support S, assign a separate variable c_a to each atom.
The ideal consists of:

- every nonconstant coefficient of det JF minus one;
- 1 plus the sum of the pure-x coefficients in F1;
- the sum of the pure-x coefficients in each other row when present;
- the exact-support equation rho times the product of all c_a minus one.

The constant determinant coefficient is already one because JF(0)=I.
The rho equation removes every coefficient-zero degeneration and ensures
that a declared support really is exact. Lower supports are independently
covered by the search.

For **all 2,326 ideals**, the reduced rational Groebner basis is {1} in both
Singular and SymPy. The seven-term contribution is 1,876; the remaining 450
are the six-term regression. Thus no coefficient assignment over any
characteristic-zero extension field can realize a collision on any retained
support. Singleton buckets and the axis proof close every other support of
size at most seven.

The [verifier](verify.py) independently compiles the polynomial maps'
derivative entries, expands the six signed Jacobian permutations, and
compares the resulting complete coefficient ideals with those supplied to
the solvers. It checks the ordered support IDs in both rational result
ledgers and can rerun both backends. These computations establish the stated
finite exclusion; they are not a formal proof-assistant certificate or an
external review.

## Reproduction and remaining boundary

From the repository root, using the versions in [validation](validation.json):

~~~sh
.venv/bin/python research/jacobian-followup/sparse-seven/verify.py --replay-algebra
.venv/bin/python research/jacobian-followup/sparse-seven/verify_regression.py
~~~

The verification includes a fresh independent enumeration of the whole cover
and both rational algebra backends. To regenerate the production shards in
a separate directory:

~~~sh
.venv/bin/python research/jacobian-followup/sparse-seven/run_all.py \
  --output tmp/sparse-seven-replay --independent
~~~

An interrupted shard is labelled incomplete and retains its DFS frontier.
The solver refuses to certify an unfinished shard. Source changes require a
new output directory instead of silently resuming old checkpoints.

This predecessor closes seven-term attainment. The later
[eight-term follow-up](../sparse-eight/PROOF.md) also excludes eight terms
and is followed by the [nine-term exclusion](../sparse-nine/PROOF.md),
which raises the current normalized lower bound to ten. The six- and
seven-term certificates above remain preserved; no attained global
support minimum follows from them.

