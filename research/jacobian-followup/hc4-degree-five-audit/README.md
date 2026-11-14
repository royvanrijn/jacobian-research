# Degree five Hessian theorem reproduction and closure

**The degree-at-most-five HC4 problem is closed here by reproduction of
David Niedbala Giraudin's November preprints.** Every four-variable polynomial
potential of degree at most five with nonzero constant Hessian determinant
has a polynomial gradient inverse. This is an external theorem, not a
novelty claim by this repository.

The sources are the [6 November full statement](https://hexagonmath.org/pdf/2611.00161v1)
and its [5 November ternary-top dependency](https://hexagonmath.org/pdf/2611.00134v1).
Their original programs pass all **276** reported checks in this environment.
The [normalization review](NORMALIZATION_REVIEW.md) accounts for the complete
case split and the imported theorems. An additional independent calculation
excludes the entire 35-coefficient rank-two passive-block class, avoiding
any dependence on leaf-specific diagonal normalizations.

The result is registered as `HC4D5EXT1`, with `HC4D5Z1` recording the extra
block exclusion. It closes `OP-HC4-D5`. Our older projective and denominator
calculations remain preserved evidence. Full independent reimplementation,
formal verification and external review are not claimed.

| Replay | Checks | Elapsed |
| --- | ---: | ---: |
| Binary preliminaries | 19 | 5.3 s |
| Binary branches | 148 | 41.8 s |
| Unary part 1 | 12 | 4.4 s |
| Unary part 2 | 21 | 4.5 s |
| Unary part 3 | 47 | 5.3 s |
| Ternary dependency | 29 | 110.5 s |

These are exact symbolic identities and characteristic-zero ideal
calculations, not a search over sampled coefficient values. The independent
block calculation used 35 coefficient variables, two inverse variables and
167 equations; its Groebner basis is `[1]`.

## Scope for the remaining research

The closure concerns the **degree of the potential**. A degree-five graph
substituted into the Meng–Yang construction can yield a much higher-degree
potential, so `OP-HC4-SCHUR` stays open. Under cotangent doubling, this result
only supplies the already-known plane degree-at-most-four range; it does
not improve the retained JC2 degree frontier. Tang's unrestricted structural form is reproduced separately in the
[structural review](../tang-hc4-audit/STRUCTURAL_REVIEW.md). Its consequences
now cover potential degree 101 and marked-collision graph degree 48, while
unrestricted HC4 remains unresolved. The subsequent
[all-degree graph theorem](../hc4-polynomial-graphs/README.md) closes the
standard polynomial-graph collision route independently of these degree bounds.

## Evidence and replay

[Source receipts](sources.json), [extracted-source hashes](source-files.json),
the six `author-*.json` receipts and their full logs preserve the actual
replays. Author source files are byte-identical to the downloaded archives.
The only compatibility change is an import alias from the author's
passagemath-singular module to installed full Sage.

```sh
# Retained-evidence integrity check; does not rerun the CAS calculations:
.venv/bin/python research/jacobian-followup/hc4-degree-five-audit/verify_evidence.py

# Fresh bounded author and independent-block replays; outputs are preserved
# under research/artifacts/local in a newly allocated directory:
.venv/bin/python research/jacobian-followup/hc4-degree-five-audit/verify_evidence.py --fresh-author --fresh-zquad
```

Runtime: Sage 10.9 / Python 3.14.7, python-flint 0.9.0 and SymPy 1.14.0.
The installed Sage launcher is `/home/royvanrijn/.local/bin/sage`.
To supply the wheel used by the compatibility adapter:

```sh
uv pip install --python /home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python --target /tmp/jacobian-hc4-sage-deps python-flint==0.9.0
```

The independent Z calculation is capped at 60 seconds. An optional attempt
to prepare a smaller unsaturated lift ended without a certificate; it is
recorded as incomplete and is not used in the proof.
