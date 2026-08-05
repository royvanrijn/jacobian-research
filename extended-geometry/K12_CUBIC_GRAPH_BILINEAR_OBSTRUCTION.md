# Exact bilinear obstruction for all linear \(K_{12}\) graph families

## Status

This note proves an **exact bounded theorem over \(\mathbf Q\)** for the
twelve-variable degree-three Keller collision \(K\).  It upgrades the
finite-field scout in the
[tensor optimization frontier](K12_TENSOR_OPTIMIZATION_FRONTIER.md) to a
characteristic-zero certificate on every normalized linear graph-coordinate
family.

None of the nine families with source pivot

\[
 z_4,z_5,\ldots,z_{12}
\]

admits a one-stage target completion of degree at most two that restores
source degree at most three.  Every parameter value is covered.  This rules
out an eleven-variable descendant in this precise linear-graph/bilinear-
completion architecture; it does not prove that dimension twelve is
minimal.

## 1. The complete graph-completion problem

Let

\[
 g_a(y)=\sum_{i=1}^{12}a_i y_i,
 \qquad a_j=1,
\]

be one of the nine linear target coordinates classified in the
[coordinate-pair frontier](K12_TO_K11_COORDINATE_PAIR_FRONTIER.md).  Its
pullback \(g_a(K)\) is affine in the source pivot \(z_j\).  Restriction to
the common collision level \(g_a(K)=g_a(K(p))\) therefore gives a polynomial
graph

\[
 z_j=R_a(z_1,\ldots,\widehat z_j,\ldots,z_{12}).
\]

After substitution, choose one retained output \(L_i(a)\) of source degree
greater than three.  Let \(M_{j,i}(a)\) be the coefficient matrix of all
source monomials of degree greater than three in every nonconstant target
monomial of degree at most two in the other ten retained raw outputs.  Its
formal basis consists of ten linear and fifty-five quadratic monomials.  Let
\(b_{j,i}(a)\) be the corresponding high-degree coefficient vector of
\(L_i(a)\).  A linear/bilinear target completion exists only if

\[
 b_{j,i}(a)\in\operatorname{im}M_{j,i}(a).                    \tag{1}
\]

All matrices and minors below are reconstructed exactly over their full
rational parameter rings.

## 2. Five constant-minor families

For five pivots, rows selected at the literal parameter point yield a
nonzero constant full-column minor.  One additional augmented row gives a
second nonzero constant minor:

| pivot | parameters | selected output | matrix | column determinant | augmented/column |
|---:|---:|---:|---:|---:|---:|
| \(z_4\) | 9 | 3 | \(999\times56\) | \(3^{23}/2\) | \(-3\) |
| \(z_9\) | 10 | 2 | \(396\times55\) | \(-3^7/2\) | \(-9\) |
| \(z_{10}\) | 10 | 3 | \(323\times55\) | \(3^{25}/2\) | \(-1\) |
| \(z_{11}\) | 10 | 1 | \(537\times55\) | \(-3^{24}\) | \(-1/2\) |
| \(z_{12}\) | 10 | 1 | \(514\times55\) | \(-16\,3^{25}\) | \(-1/2\) |

Thus \(M_{j,i}(a)\) has full column rank everywhere, and adjoining
\(b_{j,i}(a)\) raises its rank everywhere.  Condition (1) fails on each
complete parameter space.  For \(z_4\), obstructing output 3 alone already
prevents simultaneous repair of bad outputs 3 and 10.

## 3. Four stratified families

The other four matrices lose rank on coordinate strata, so they require a
determinant-open cover followed by an exact closed-stratum calculation.

| pivot | matrix | generic determinant opens | closed stratum |
|---:|---:|---|---|
| \(z_5\) | \(2487\times58\) | \(\tau_0^{23}\), \(\tau_1^{23}\) | \(\tau_0=\tau_1=0\) |
| \(z_6\) | \(2951\times57\) | \(\tau_0^{17}\), \(\tau_1^{17}\), \(\tau_2^{10}(3\tau_2+7\tau_5)^3\), \(\tau_2\tau_5^{13}\) | \(\tau_0=\tau_1=\tau_2=0\) |
| \(z_7\) | \(880\times56\) | \(\tau_0^5\), \(\tau_1^5\), \(\tau_2^5\) | \(\tau_0=\tau_1=\tau_2=0\) |
| \(z_8\) | \(6711\times57\) | \(\tau_0^{21}\), \(\tau_1^{16}\), \(\tau_2^{20}\) | \(\tau_0=\tau_1=\tau_2=0\) |

Every listed column determinant has a nonzero constant factor, and the
corresponding augmented determinant is a nonzero constant multiple of it.
The \(z_6\) pair involving \(3\tau_2+7\tau_5\) and \(\tau_5\) covers the
possible cancellation within the \(\tau_2\ne0\) open.  Hence the displayed
opens cover the complement of the stated closed strata.

On the \(z_5\) closed stratum, one completion column satisfies

\[
 c_{46}=\tfrac13c_2+2\tau_2c_{47}.
\]

After removing it, two column/augmented minor pairs have common determinant
factors \((7\tau_4+3)^6\) and \(\tau_4^{14}\), respectively, and augmented-
to-column ratio \(-9\).  These two opens cover the remaining parameter
line, so (1) still fails everywhere.

On the \(z_6\) closed stratum,

\[
 c_{45}=\tfrac13c_1.
\]

Deleting that dependent column leaves constant column and augmented minors
\(-3^{18}/8\) and \(3^{20}/8\), with ratio \(-9\).

On the \(z_7\) closed stratum,

\[
 c_0=\tau_3c_{44}-2\tau_4c_{45}.
\]

Deleting it leaves equal nonzero constant column and augmented minors
\(-3^{25}\).

Finally, the \(z_8\) closed stratum is exactly the complete five-parameter
quadratic \(z_8\) graph family already excluded by
[BCR5](K12_Z8_CUBIC_COMPLETION.md), in fact through target degree three.

## 4. Conclusion and boundary

Combining the constant and stratified certificates proves:

> **All-nine linear-graph bilinear obstruction.** None of the nine full
> normalized linear graph-coordinate families with source pivot
> \(z_4,\ldots,z_{12}\) admits a one-stage target completion of degree at
> most two that lowers every retained output to source degree at most three.

The linear graph-coordinate/bilinear-completion route from \(K_{12}\) to
dimension eleven is therefore closed.  The subsequent
[cubic target-completion theorem](K12_CUBIC_TARGET_COMPLETION_OBSTRUCTION.md)
also closes target degree three for all nine families.  Still outside the
combined results are target degree at least four on the full cubic graph
families, nonlinear target coordinates, and ordered target stages.  Those
gaps are genuine: this is not a general dimension-eleven lower bound.

## 5. Reproduction

Run

```bash
make verify-k12-cubic-graph-bilinear-obstruction
```

The checker reconstructs \(K\), all nine full graph families, every linear
and quadratic completion column, every open-cover determinant, each closed-
stratum column relation, and the final augmented-rank obstructions.  It
writes
[`k12_cubic_graph_bilinear_obstruction.json`](../artifacts/generated-results/k12_cubic_graph_bilinear_obstruction.json),
including the selected source-monomial rows needed to replay every minor.
