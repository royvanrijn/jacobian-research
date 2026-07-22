# From one collision to marked-root Keller maps

This repository verifies a three-dimensional polynomial Keller map with a
three-point collision, explains it through a tangent-map normal form, and
develops weighted, cancellation, decorated-normalization, Hurwitz, and
rank-two symplectic consequences.  In generic degree `N>=4`, the coarse
decorated normalization already gives an `(N-3)`-dimensional family of stable
classes.  Adding one affine root sheet generically recovers the seed exactly.

## The foundational map

Put `u=1+xy` and define

\[
F(x,y,z)=\left(
u^3z+y^2u(4+3xy),
y+3xu^2z+3xy^2(4+3xy),
2x-3x^2y-x^3z
\right).
\]

On `x!=0`, set `t=y+1/x` and `r=2/x`.  For target coordinates `(a,b,c)`,

\[
a=t^2+rt/2-ct^3,\qquad b=r+4t-3ct^2.
\]

The two coordinate Jacobians give

\[
\det\frac{\partial(a,b,c)}{\partial(t,r,c)}
\det\frac{\partial(t,r,c)}{\partial(x,y,z)}
=\frac r2(-2x)=-2,
\]

so `det DF=-2` polynomially everywhere.  Exact substitution gives

\[
F(0,0,-1/4)=F(1,-3/2,13/2)=F(-1,3/2,13/2)=(-1/4,0,0).
\]

The minimal independent replay is:

```bash
python3 scripts/verify_counterexample_independent.py
```

## Canonical proof path

The common geometric normal form is

\[
(W,s)\longmapsto (s,Ws-H(W)).
\]

The [tangent-map core](verified/TANGENT_MAP_CORE.md) supplies the inverse
pencil, generic degree, critical divisor, discriminant normalization, full
Hessian Fitting divisor, incidence form, weighted suspension, and comparison
with cancellation suspension.  Later papers cite this theorem instead of
rederiving those calculations.

The primary dependency chain is:

1. `F1` — foundational Keller collision.
2. `W1` — tangent-map core and weighted suspension.
3. `S1` — stable normalization functoriality.
4. `C1` — universal cancellation construction.
5. `B1` — complete canonical boundary exhaustion.
6. `P1` — cancellation reconstruction residue and parameter faithfulness.
7. `M1` — finite degreewise stable multiplicity.
8. `D1` — decorated-normalization moduli of dimension `N-3`.
9. `F2` — generic affine-mark faithfulness.
10. `H1` — internal Hurwitz--LL compactification theorem.
11. `R1`, `R2` — all-degree rank-two descent and parameter faithfulness.

The exact scopes, dependencies, checkers, and review states live only in [`MATH_STATUS.json`](MATH_STATUS.json).  [STATUS.md](STATUS.md) is generated
from it.  A checker is reproducibility evidence; it is not external review.

## Main papers

- [Foundational counterexample](papers/core-counterexample/main.pdf)
- [Discriminant pencils](papers/discriminant-pencils/main.pdf)
- [Decorated discriminant normalization](papers/decorated-discriminant-normalization/main.pdf)
- [Marked-root degreewise multiplicity](papers/marked-root-multiplicity/main.pdf)
- [Hurwitz--LL rerooting](papers/hurwitz-ll-rerooting/main.pdf)

The degree-five and degree-six calculations are worked regressions generated
from the all-degree results, not independent theorem authorities.  Quartic
models and external-map classifications are likewise examples or audits.

## Reproduction

Create the pinned Python environment in [REPRODUCE.md](REPRODUCE.md), then run:

```bash
make check verify-minimal verify-master verify-theorems verify-papers
```

`make check` validates links and the status graph.  `verify-minimal` needs only
the Python standard library.  `verify-master` covers the cancellation chain,
including the unconditional all-`(m,r)` thick-contact formula.  The conditional
contact-resultant refinements are separate from the `M1` proof path.
`verify-papers` discovers and compiles every `papers/*/main.tex`; CI uses the
same discovery rule and archives every resulting PDF.

The full command catalogue, heavier symbolic runs, optional Lean certificate,
and independent Macaulay2 comparison are documented in
[REPRODUCE.md](REPRODUCE.md).

## Research status

There are three primary continuation problems:

- `OP-MARK`: valuative affine-sheet extension through every collision stratum.
- `OP-CR`: cancellation contact resultants for `r>=5`.
- `OP-LR`: intrinsic target coercivity and no escape for filtered LR contact.

In particular, `b_m=34m+1` is an exact source-only profile in one
determinant-normalized target gauge.  It is not target-minimal and is not a
stable LR obstruction.  The target torus lowers the first two degrees from
`(35,69)` to `(25,49)`; the continuation `24m+1` remains conjectural.

Other questions—arithmetic Galois theory, controlled-boundary suspension
classification, wider quantization, coefficient-scheme gluing,
quadratic--cubic flexibility, the plane degree frontier, and Gaussian or
Image/Vanishing consequences—are retained as parked side programmes in
[STATUS.md](STATUS.md).

## Provenance

The earliest public item located by the repository audit is
[Levent Alpöge's post](https://x.com/__alpoge__/status/2079028340955197566).
Contemporaneous sources attribute the example to Alpöge and Fable; Dean
Cureton supplied a separate
[Lean certificate](https://github.com/deancureton/jacobian), and Alexis
Gallagher developed the weighted-lift viewpoint in a
[same-day explainer](https://jacobianfun.org/jacobian-explained).

The surviving record does not settle the complete discovery history, original
prompt, model conversation, or exact UTC timestamp.  This repository repeats
the contemporaneous attribution without making a priority claim.  See the
[provenance audit](archive/legacy-notes/PROVENANCE_AUDIT.md) for the detailed
record and the canonical sources for precise mathematical claims.

## Repository policy

Proofs belong in canonical sources, statuses belong in `MATH_STATUS.json`, and
open-problem IDs belong in `STATUS.md`.  Papers and notes may cite those IDs
but do not maintain independent continuation queues.  Superseded derivations
remain available under `archive/` and outside primary navigation.
