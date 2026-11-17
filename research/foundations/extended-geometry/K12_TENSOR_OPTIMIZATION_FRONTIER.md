# Tensor and cubic-graph optimization frontier for (K_{12})

## Status

This note records two different kinds of progress toward compressing the
[twelve-variable degree-three Keller collision](../verified/TWELVE_VARIABLE_DEGREE_THREE_KELLER_COUNTEREXAMPLE.md):

1. an **exact rational tensor-module obstruction** for the displayed
   (K_{12}) tensor and its cubic-homogeneous parent (G_{19});
2. a **bounded finite-field discovery experiment** in the genuinely cubic
   graph-coordinate families left open by the earlier exact completion
   theorems.

The first result is a theorem about these displayed tensors.  The second is
not a proof over (mathbf Q), and neither result proves that dimension twelve
is minimal.

## 1. The two tensor flattenings

Let

\[
 T\in W\otimes \operatorname{Sym}^d(V^*)
\]

encode a homogeneous polynomial map (H:V\to W).  Two elementary
flattenings directly measure the constant row and column modules of (JH):

\[
 \lambda_T:W^*\longrightarrow \operatorname{Sym}^d(V^*),
 \qquad u\longmapsto u\mathbin{\cdot}H,
\]

and

\[
 \delta_T:V\longrightarrow W\otimes\operatorname{Sym}^{d-1}(V^*),
 \qquad v\longmapsto D_vH.
\]

Thus

\[
 \ker\lambda_T=\{u:u^TJH=0\},\qquad
 \ker\delta_T=\{v:JH\,v=0\}.
\]

The second equality means that a nonzero constant source direction can be
removed by a linear tensor quotient only if (delta_T) has a kernel.  If

\[
 H=\sum_{s=1}^R u_s\ell_s^d,
\]

then both flattening ranks are at most (R).  Their maximum is therefore a
rigorous lower bound for a partially symmetric pure-power decomposition.

## 2. Exact ranks for (K_{12}) and (G_{19})

Write the nonlinear part of (K_{12}) as (Q+C), with (Q) quadratic and
(C) cubic.  Exact coefficient row reduction over (mathbf Q) gives:

| tensor | output rank (operatorname{rank}\lambda_T) | input rank (operatorname{rank}\delta_T) | common left kernel | common right kernel | coefficient-matrix span | pure-power lower bound |
|---|---:|---:|---:|---:|---:|---:|
| (Q) | 11 | 12 | 1 | 0 | 12 | 12 |
| (C) | 6 | 12 | 6 | 0 | 19 | 12 |
| (H_{19}) | 18 | 19 | 1 | 0 | 38 | 19 |

For (Q), the left kernel is the sixth output direction because (Q_6=0).
For (C), it is exactly the last six output directions, recovering the
known cubic-output rank six.  For (H_{19}), the sole left annihilator is
the nineteenth, fixed (	au) output.

The important new obstruction is on the input side:

> **Constant-module obstruction.** The cubic tensor (C) of (K_{12})
> uses all twelve source directions, and the cubic tensor (H_{19}) uses all
> nineteen source directions.  Equivalently, the coefficient matrices of
> their Jacobians have zero common right kernel.  No linear change of source
> coordinates can expose a constant direction on which either tensor is
> independent.

Consequently, quotienting a common constant row/column module cannot compress
the displayed (K_{12}) cubic tensor, and it cannot reduce (G_{19}) below
nineteen variables.  Also, a decomposition of the displayed tensors as
(sum u_s\ell_s^3) needs at least twelve and nineteen summands,
respectively.  These statements do not obstruct nonlinear graph coordinates,
nonconstant modules, exact Schur elimination, or a different tensor in the
same Keller-collision problem.

The checker also replays both rational collisions and the exact companion
scaling identity

\[
 E_\tau(z)=\tau^{-1}K_{12}(\tau z)
\]

used in the determinant certificate for (G_{19}).  The direct determinant
of (K_{12}) remains in the canonical BCR2 verifier rather than being
duplicated here.

## 3. Genuinely cubic graph-coordinate scout

The exact
[linear graph-coordinate classification](K12_TO_K11_COORDINATE_PAIR_FRONTIER.md)
has nine normalized families with source pivots (z_4,\ldots,z_{12}).
Subsequent exact theorems close the six subfamilies whose graph correction is
quadratic through target degree two, and much further in several cases.  The
larger families with cubic graph correction were still unscreened even for
bilinear target completion.

The new sparse scout evaluates integer parameter points with:

- support at most two and nonzero values in
  ({-2,-1,1,2});
- 250 additional deterministic random points per parameter count;
- independent reduction modulo (101) and (103);
- all ten linear and all fifty-five quadratic monomials in the other raw
  retained outputs as the completion basis for each bad component.

The same signed integer points are used at both primes.  The census is:

| pivot | parameters | points per prime | genuinely cubic graphs per prime | bilinear survivors mod 101 | bilinear survivors mod 103 |
|---:|---:|---:|---:|---:|---:|
| 4 | 9 | 863 | 863 | 0 | 0 |
| 5 | 7 | 615 | 615 | 0 | 0 |
| 6 | 8 | 731 | 731 | 0 | 0 |
| 7 | 9 | 862 | 749 | 0 | 0 |
| 8 | 8 | 729 | 545 | 0 | 0 |
| 9 | 10 | 1,011 | 830 | 0 | 0 |
| 10 | 10 | 1,011 | 830 | 0 | 0 |
| 11 | 10 | 1,011 | 830 | 0 | 0 |
| 12 | 10 | 1,011 | 830 | 0 | 0 |

There are 15,688 modular point evaluations in total, including 13,646
evaluations of genuinely cubic graphs.  No point makes every high-degree
retained component bilinearly completable at either prime.  Hence there is no
rational reconstruction candidate in this bounded support window.

The subsequent
[exact bilinear obstruction](K12_CUBIC_GRAPH_BILINEAR_OBSTRUCTION.md)
does more than reconstruct isolated candidates: constant minors and
rank-stratified determinant covers now exclude the complete parameter spaces
for all nine pivots (z_4,\ldots,z_{12}) over (mathbf Q).  The modular census
remains the broader experimental record that selected the support patterns
used in the exact calculation.

This negative census is a search result, not an exact exclusion.  In
particular, it does not cover arbitrary rational parameter values, target
degree at least three on the cubic graph families, nonlinear target
generators, or ordered completion stages.  The subsequent
[cubic target-completion theorem](K12_CUBIC_TARGET_COMPLETION_OBSTRUCTION.md)
now supplies the exact degree-three exclusion.  The value of this earlier
census is to remove the
sparsest support patterns from the discovery queue and to provide a reusable
two-prime implementation for broader or structured sampling.

## 4. Next exact targets

The calculations redirect the compression search away from constant tensor
quotients and one-stage target completion through degree three on linear
graph coordinates.  The most economical remaining branches are:

1. raise the target completion degree to four on the full cubic graph
   families, using the cubic obstruction minors as fixed Schur pivots;
2. allow an ordered completion whose first stage changes the
   generators available to the second;
3. search nonlinear graph coordinates, nonconstant common modules, or a
   jointly affine Schur block rather than a constant right kernel;
4. start with a different low-rank collision tensor instead of deforming the
   full-input-rank (K_{12}) cubic tensor.

## 5. Reproduction

Run the exact tensor audit with:

```bash
.venv/bin/python scripts/verify_k12_tensor_module_frontier.py
```

Run the pinned modular scout with:

```bash
.venv/bin/python scripts/search_k12_cubic_graph_bilinear_completions.py \
  --support-max 2 --values=-2,-1,1,2 --random-samples 250 \
  --random-seed 20260804 --primes 101,103 --keep-closest 5 \
  --output artifacts/generated-results/k12_cubic_graph_bilinear_modular_search.json
```

The generated records are
[`k12_tensor_module_frontier.json`](../artifacts/generated-results/k12_tensor_module_frontier.json)
and
[`k12_cubic_graph_bilinear_modular_search.json`](../artifacts/generated-results/k12_cubic_graph_bilinear_modular_search.json).
