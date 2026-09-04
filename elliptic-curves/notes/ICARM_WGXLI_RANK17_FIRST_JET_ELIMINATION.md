# Projective first-jet elimination for the five-fibre `wgxli` target

Status: **exact finite-field necessary-condition elimination in every
normalized projective parameter chart; the literal displayed-sign,
displayed-order interpolation is rejected modulo 17, 53, and 67**.

This is not a characteristic-zero nonexistence theorem. It assumes that the
first seventeen public points on curves 351, 356, 376, 377, and 385 are
literally corresponding sections of a rootless-K3 `(8,12;4,6)` model.

## Outcome

After the weighted gauges

```text
t_351=0,  t_356=1,  t_385=-1,  u_351=1,
```

the eliminator exhausts the ordered pair `(t_376,t_377)` in
`P1(F_p)\{0,1,-1}`. This includes the finite chart and both distinct boundary
charts in which exactly one residual parameter is infinity. The results are

```text
field     projective ordered pairs     geometric solutions     timeouts
GF(17)                         210                       0            0
GF(53)                        2550                       0            0
GF(67)                        4160                       0            0
```

For each fixed pair, `msolve` returns the unit ideal `[1]`. Thus the test is
over the algebraic closure, not a search for only `F_p`-rational solutions.
A literal rational model would have to leave the distinct normalized chart
through bad or colliding parameter/scaling reduction at all three primes.

## Projective elimination

Write a base point as `[T:Z]`, and let the five normalized nodes be
`[a_k:b_k]`. For a proposed family

```text
Y^2 = X^3 + A(T,Z) X + B(T,Z)
```

the binary degrees are 8 and 12, while each section has binary degrees 4 and
6. Put

```text
L(T,Z) = product_k (b_k*T-a_k*Z).
```

Five values determine a binary quartic `x_i`. Every binary sextic ordinate
through its five displayed values is uniquely

```text
y_i = ybar_i + L*ell_i,
```

where `ell_i` is a binary linear form. This isolates, and then eliminates,
the 34 free ordinate coefficients for the seventeen sections.

Differentiate the section equation in the local parameter at each node
(`t=T/Z` at finite nodes and `s=Z/T` at infinity). If

```text
alpha_k = A'(t_k),       beta_k = B'(t_k),
c_ik = 2*y_i(t_k)*L'(t_k),
R_ik = x_i(t_k)*alpha_k + beta_k
       - (2*y_i*ybar_i' - 3*x_i^2*x_i' - A*x_i')(t_k),
```

then

```text
c_ik*ell_i(t_k) = R_ik.
```

The five right sides divided by `c_ik` must be values of one binary linear
form. Using the nodes at 0 and 1 eliminates `ell_i` and gives three
compatibility equations per section, hence 51 equations. The five values and
five local derivatives of the binary octic `A` satisfy the unique left-kernel
relation of their `10 x 9` evaluation matrix. One more equation saturates the
four remaining fibre scales. Each modular chart therefore has 15 variables:
four nonlinear scales, ten derivative auxiliaries, and one saturation
inverse. There are 52 nonzero equations at 17 and 53 at both 53 and 67.

This is the projective form of the earlier affine Hermite calculation. It
also explains the pair counts

```text
(p-2)*(p-3),
```

rather than the affine-only `(p-3)*(p-4)` count.

## Published-R17 positive controls

The compact certified published-R17 model and all seventeen reconstructed
sections replay the eliminated equations at both clean primes requested for
the collision test. At each of 53 and 67 the checker uses a nonsingular
finite control chart and both nonsingular infinity orientations; every chart
has all 53 nonzero equations and an exact witness. Thus the clean-prime empty
target charts are not explained by the mod-17 degeneration or by omitting a
projective boundary.

The mod-17 compact control also gives exact algebraic witnesses, but the fixed
control fibre at `q=1` is singular and some reduced control equations vanish
identically. It is deliberately retained only as a degenerate regression
witness. The 53 and 67 controls are the nonsingular, nondegenerate positive
controls used to rule out an implementation-level explanation of the empty
target charts.

## Interpretation

The literal basis is only the first exact gate. The bounded sign,
permutation, and elementary-mutation analysis is recorded in
[`ICARM_WGXLI_RANK17_BOUNDED_REBASING.md`](ICARM_WGXLI_RANK17_BOUNDED_REBASING.md).
Neither numerical height-Gram correlation nor this necessary-condition test
alone proves or disproves unrestricted common-family membership.

## Reproduction

From the repository root, with Sage and `msolve` available:

```bash
sage -python \
  elliptic-curves/cas/eliminate_icarm_wgxli_rank17_first_jet.sage \
  --prime 17 --jobs 4 --threads 1 --pair-timeout 60 \
  --output \
  artifacts/generated-results/elliptic-curves/icarm_wgxli_rank17_first_jet_mod17_v2.json

sage -python \
  elliptic-curves/cas/eliminate_icarm_wgxli_rank17_first_jet.sage \
  --prime 53 --jobs 14 --threads 1 --pair-timeout 60 \
  --output \
  artifacts/generated-results/elliptic-curves/icarm_wgxli_rank17_first_jet_mod53_v2.json

sage -python \
  elliptic-curves/cas/eliminate_icarm_wgxli_rank17_first_jet.sage \
  --prime 67 --jobs 32 --threads 1 --pair-timeout 180 \
  --reuse-unit-outputs \
  --output \
  artifacts/generated-results/elliptic-curves/icarm_wgxli_rank17_first_jet_mod67_v1.json
```

The 67 run was first exhausted with a 60-second per-chart limit; its timeout
charts were then rerun with the displayed 180-second command. Reuse occurs
only when the old output is exactly the unit ideal and its newly rendered
input hash matches the previous artifact. Missing, partial, nonunit, and
hash-mismatched outputs are solved again. Use the same command with `--check`
to replay an artifact. Runtime fields are ignored by the checker; all
solver-input hashes and mathematical statuses are compared.

The generated mod-17, mod-53, and mod-67 artifacts have SHA-256 hashes

```text
67d642871268d339c5b4c8ea55e601c546aa7da089ceb07f5d21664cec5b8994
ca5547ccd1a246f020c057715e798d951863e9b679955f365037d21ea399fa2a
0cba1fc3c6ff80509ea564de9f945b3a3d1eb3047df0c05de64ea5e31fe45041
```
