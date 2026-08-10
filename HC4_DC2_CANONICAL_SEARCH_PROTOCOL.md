# Shared canonical-transformation search for `HC4` and `DC2`

## Status and purpose

This is a search protocol, not a theorem about `HC4` or `DC2`. It makes the
canonical-word layer common to both programmes while keeping their admission
gates logically separate.

The exact fixed-dimensional Hessian result is `HC4MR2`:

\[
JC2\quad\Longleftrightarrow\quad PHC4,
\]

where `PHC4` is `HC4` restricted to potentials lying on a nontrivial
constant-Hessian pencil. See the
[relative-nilpotent master reduction](HC4_RELATIVE_NILPOTENT_MASTER_REDUCTION.md#exact-restricted-equivalence-with-jc2).
The missing unrestricted step is pencil recognition, not another reduction
inside the pencil class.

At the software level,
[`scripts/canonical_transform_search.py`](scripts/canonical_transform_search.py)
now supplies the convention-explicit Poisson matrix, shared mixed-line
Hamiltonian alphabet, Hamiltonian vector fields, exact invariant shears, word
composition, pullback, and symplectic verification used by both the HC4
mixed-pivot search and the DC2 Hessian-symbol search. The two consumers retain
different mathematical gates.

## 1. What is genuinely shared

For Darboux variables `z=(q,p)` and a constant Poisson matrix `Pi`, a declared
canonical-word letter is retained only when its time-one map is exactly

\[
T_H(z)=z+\Pi\nabla H(z).
\tag{1.1}
\]

The engine verifies that every component of the velocity is invariant under
the Hamiltonian flow; (1.1) is therefore a finite exact flow, not a truncated
Lie series. It records the sign convention for `Pi`, composes words by
simultaneous substitution, and checks

\[
DT\,\Pi\,(DT)^{\mathsf T}=\Pi,
\qquad \det DT=1.
\tag{1.2}
\]

The DC2 Cayley generator also uses `Pi*grad(H)` directly after independently
checking its Jacobian integrability; it does not label every such map a
time-one shear. These identities are common infrastructure. They do not imply
either of the programme-specific statements below.

| gate | `HC4` consumer | `DC2` consumer |
|---|---|---|
| ambient object | a potential pulled back by a canonical chart | a four-dimensional classical symplectic map or symbol |
| indispensable exact test | constant nonzero Hessian determinant after rechart/reduction | all six Poisson relations and determinant one |
| geometric status | collision retained or explicitly transferred | classical nonautomorphism/noninjectivity proved |
| hard rejection | parent or descended Hessian determinant nonconstant | polynomial inverse or known automorphism orbit |
| terminal work | recognize a relative pencil and invoke `HC4MR2` | construct one Weyl endomorphism, then prove it non-surjective |

A nonlinear symplectic change of the independent variables does not generally
preserve a constant Hessian determinant. Conversely, a constant-Hessian
gradient need not be symplectic for the fixed Darboux form. The common engine
therefore shares transformations and orbit data, not the two admission
theorems.

## 2. The automorphism veto for `DC2`

Every Hamiltonian shear in the current canonical alphabet has the polynomial
inverse `T_{-H}`. Every finite word in these letters is consequently a
polynomial symplectic automorphism. In the declared elementary alphabet its
factorwise Moyal/Weyl lift is also an exact Weyl automorphism.

Therefore:

> **Automorphism veto.** A canonical word by itself cannot be a `DC2`
> counterexample. It is an orbit action, polarization, or calibration
> control. Only a four-dimensional symplectic etale map whose polynomial
> inverse is not already known may enter the `DC2` candidate queue.

This is the lesson of the two-pencil false positives in the
[Hessian-symbol optimization](extended-geometry/DC2_HESSIAN_SYMBOL_OPTIMIZATION.md):
a native-support obstruction can disappear after one correction monomial is
added, while factorwise composition already supplies the exact Weyl
automorphism.

## 3. Orbit-level treatment of the moving-symbol problem

The unit of search is a **canonical orbit record**, not one principal symbol.
For every seed `F` and every canonical word `T`, record:

1. exact formulas for `T`, `T^{-1}`, and the transformed seed;
2. the classical automorphism/nonautomorphism status of the seed;
3. the collision or fiber data transported by the word;
4. the factorwise Weyl lift of every declared word letter;
5. native and word-saturated correction supports;
6. the first obstruction that survives every known factorization and
   polarization in the orbit.

Killing one symbol removes only that representative. Killing an entire
declared orbit requires either a coordinate-invariant obstruction or an exact
proof that every permitted word and support saturation fails. Conversely, a
known factorwise Weyl lift overrides a restricted-support rank jump for that
representative.

## 4. Shared queue

The next searches should use one generated word ledger with two consumers.

### `HC4` lane

1. Start from a collision-normalized Hessian or cotangent seed.
2. Apply the common mixed-line and coisotropic word alphabet.
3. Audit the parent Hessian determinant after every nonlinear rechart.
4. Retain only charts that transfer a marked collision.
5. Test affine, oblique, and coisotropic Schur blocks.
6. On any descended four-variable constant-Hessian potential, run the four
   projective constant-null-covector charts of `HC4MR3` first.
7. Only if the rank-one scheme is empty, solve higher-rank pencil-admission
   equations. `HC4MR2` makes any positive admission a reduction to `JC2`
   rather than a new HC4 frontier.

### `DC2` lane

1. Accept only four-dimensional exact symplectic etale seeds.
2. Apply the same word ledger as orbit/polarization changes.
3. Reject every seed with a polynomial inverse before PBW work.
4. Quantize the word factors exactly and saturate the correction support they
   generate.
5. Require a surviving all-relation Weyl endomorphism; a principal-symbol or
   finite-order obstruction is not enough.
6. Only then test nonsurjectivity by an invariant stable under the recorded
   canonical orbit.

The HC4 and DC2 programmes thus share candidate generation, exact canonical
words, and orbit bookkeeping. They remain separate at the two places where
the mathematics is genuinely different: Hessian-pencil recognition for HC4,
and existence plus nonsurjectivity of a rank-two Weyl endomorphism for DC2.

## 5. Reproduction

Verify the common rank-two/rank-three conventions, word inverses, and
symplectic identities with

```bash
.venv/bin/python scripts/verify_canonical_transform_search.py
```

The established finite searches retain their existing commands:

```bash
PYTHONHASHSEED=0 .venv/bin/python \
  scripts/search_hc4_mixed_canonical_pivots.py \
  --output artifacts/generated-results/hc4_mixed_canonical_pivot_search.json

.venv/bin/python scripts/search_dc2_hessian_symbol_candidates.py \
  --shortlist 8 \
  --output artifacts/generated-results/dc2_hessian_symbol_optimization.json
```

The first command tests Hessian and descent gates; the second tests
Hamiltonian-Hessian/Cayley, symplectic, and bounded quantization gates. A
shared canonical primitive does not merge their conclusions.
