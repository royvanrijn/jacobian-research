# Exact cubic target-completion obstruction for all linear \(K_{12}\) graphs

## Status

This note proves an **exact bounded theorem over \(\mathbf Q\)** for the
twelve-variable degree-three Keller collision \(K\).  None of the nine full
normalized linear graph-coordinate families with source pivot

\[
 z_4,z_5,\ldots,z_{12}
\]

admits a one-stage target completion of degree at most three that restores
source degree at most three.  This strengthens the
[all-nine bilinear obstruction](K12_CUBIC_GRAPH_BILINEAR_OBSTRUCTION.md) by
including every cubic target monomial.  It rules out an eleven-variable
descendant in this precise architecture, not every possible compression of
\(K\).

## 1. Schur-reduced discovery

For a selected retained output \(L_i(a)\), let \(M^{(3)}_{j,i}(a)\) contain
the source-degree-above-three coefficient vectors of all nonconstant target
monomials of degrees one, two, and three in the other ten raw retained
outputs.  The formal basis has

\[
 10+55+220=285
\]

members.  A completion requires

\[
 b_{j,i}(a)\in\operatorname{im}M^{(3)}_{j,i}(a),              \tag{1}
\]

where \(b_{j,i}(a)\) is the high-degree part of \(L_i(a)\).

The finite-field scout first constructs the certified linear/bilinear
column space from the preceding theorem and then adds the 220 cubic columns
in its quotient.  It tests support at most one with nonzero values
\(\{-2,-1,1,2\}\), plus 250 deterministic random points per parameter
count, over both \(\mathbf F_{101}\) and \(\mathbf F_{103}\).  Across 2,902
evaluations there is no survivor.  More strongly, every tested point has
the same selected-component signature: the bilinear rank is 56 for \(z_4\)
and 55 for the other four constant-rank families, all 220 cubic columns are
independent in the quotient, and adjoining the defect raises rank by one.

That census is only a discovery experiment.  The rest of this note is an
exact reconstruction over the complete rational parameter spaces.

## 2. Six constant-minor families

Sparse elimination at a rational parameter point selects all nonzero
completion columns and one augmented row.  Reconstructing those row sets
over the full parameter rings gives the following nonzero constants:

| pivot | matrix | nonzero columns by target degree | column determinant | augmented/column |
|---:|---:|---:|---:|---:|
| \(z_4\) | \(12214\times276\) | \(1,55,220\) | \(-2^{32}3^{74}23\) | \(2/3\) |
| \(z_8\) | \(108338\times277\) | \(2,55,220\) | \(2^{23}3^{144}\) | \(1\) |
| \(z_9\) | \(5432\times275\) | \(0,55,220\) | \(2^8 3^{35}7\) | \(-1\) |
| \(z_{10}\) | \(4059\times275\) | \(0,55,220\) | \(2^{23}3^{85}23\) | \(7\) |
| \(z_{11}\) | \(6842\times275\) | \(0,55,220\) | \(-2^{24}3^{84}7\,23\) | \(-1/2\) |
| \(z_{12}\) | \(6433\times275\) | \(0,55,220\) | \(-2^{39}3^{84}7\,23\) | \(-1/2\) |

Thus each column matrix has full rank everywhere, while its augmentation by
\(b_{j,i}(a)\) has strictly larger rank everywhere.  Notice that \(z_8\),
which required a stratified argument at target degree two, acquires a single
constant full-column obstruction minor after the cubic columns are added.

## 3. The \(z_5,z_6,z_7\) determinant covers

The remaining three families retain the same rank-drop strata seen in the
bilinear theorem.  Their exact generic column determinants have the
following parameter factors; every corresponding augmented determinant is
the listed nonzero constant multiple.

| pivot | generic determinant factors | augmented/column | closed stratum |
|---:|---|---:|---|
| \(z_5\) | \(\tau_0^{43}\), \(\tau_1^{42}\) | \(9\) | \(\tau_0=\tau_1=0\) |
| \(z_6\) | \(\tau_0^{51}\), \(\tau_1^{50}\), \(\tau_2^{45}(3\tau_2-7\tau_5)(3\tau_2+7\tau_5)^3(6\tau_2+7\tau_5)\), \(\tau_2^2\tau_5^{56}\) | \(-1\) | \(\tau_0=\tau_1=\tau_2=0\) |
| \(z_7\) | \(\tau_0^7\), \(\tau_1^7\), \(\tau_2^7\) | \(7\) | \(\tau_0=\tau_1=\tau_2=0\) |

For \(z_6\), if \(\tau_2\ne0\) and any of the three linear factors in the
third determinant vanishes, then \(\tau_5\ne0\), so the fourth determinant
is nonzero.  The displayed opens therefore cover the complements of the
closed strata.

On those strata the same exact column relations as at target degree two
remain valid:

\[
\begin{aligned}
 z_5:&\quad c_{46}=\tfrac13c_2+2\tau_2c_{47},\\
 z_6:&\quad c_{45}=\tfrac13c_1,\\
 z_7:&\quad c_0=\tau_3c_{44}-2\tau_4c_{45}.
\end{aligned}
\]

After deleting the dependent column, the \(z_5\) closed stratum is covered
by minors with factors \((7\tau_4+3)^8\) and \(\tau_4^{40}\), both with
augmented ratio 9.  The \(z_6\) and \(z_7\) closed strata have constant
column determinants \(-2^9 3^{59}7\) and \(-2^{22}3^{88}\), with augmented
ratios \(-1\) and 7.  Hence (1) fails on every closed stratum as well.

## 4. The exact bounded theorem

Combining all certificates gives:

> **All-nine cubic-completion obstruction.** None of the nine normalized
> linear graph-coordinate families with source pivot \(z_4,\ldots,z_{12}\)
> admits a one-stage target completion of degree at most three that lowers
> every retained output to source degree at most three.

The remaining routes begin with target degree four on the full cubic graph
families, nonlinear graph coordinates, or ordered target stages.  This
theorem does not establish a general dimension-eleven lower bound.

## 5. Reproduction

Run the finite-field discovery pass with:

```bash
.venv/bin/python scripts/search_k12_schur_cubic_completions.py \
  --support-max 1 --values=-2,-1,1,2 --random-samples 250 \
  --random-seed 20260804 --primes 101,103 --keep-closest 5 \
  --output artifacts/generated-results/k12_schur_cubic_completion_modular_search.json
```

Run the exact reconstruction with:

```bash
make verify-k12-graph-cubic-completion-obstruction
```

The generated records are
[`k12_schur_cubic_completion_modular_search.json`](../artifacts/generated-results/k12_schur_cubic_completion_modular_search.json)
and
[`k12_graph_cubic_completion_obstruction.json`](../artifacts/generated-results/k12_graph_cubic_completion_obstruction.json).
