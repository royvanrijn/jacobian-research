# Bounded rebasing rejection for the five-fibre `wgxli` target

Status: **formal finite rejection inside the declared smallest signed,
permuted, and one-elementary-mutation rebasing bound**.

This bounded theorem remains valid, but its former negative operational
interpretation is superseded.  The complete norm-twelve atlas reconstructs
the family after a unimodular rebasing with seven simultaneous differences,
outside the one-mutation bound below.  See
[`../../elkies-k3/R17_NORM12_RECORD_LINEAGE_SWEEP_2026-09-04.md`](../../elkies-k3/R17_NORM12_RECORD_LINEAGE_SWEEP_2026-09-04.md).

<!-- status-consumer: EC-K3-R17-NORM12-RECORD-LINEAGE-ATLAS 291a539d07b842b9 -->

The five records 351, 356, 376, 377, and 385 share a strong numerical and
private-pipeline fingerprint. They do not give a rootless-K3 `(8,12;4,6)`
realization under the finite rebasing declared below. This is not an
unrestricted `GL(17,Z)` conclusion.

## Declared bound and outcome

The stages are nested:

1. exhaust all relative diagonal signs, modulo global section inversion;
2. permit only fingerprint-indistinguishable label permutations;
3. if those fail, permit one common unimodular elementary column mutation,
   fixing the stable denominator anchors, with column `l1` norm at most 2,
   changed-height ratio at most 1.25 on every fibre, and at least a 1% joint
   Gram-objective improvement.

The first two stages retain only displayed relative signs and displayed
order. The third stage enumerates 352 transformations and retains exactly

```text
P4 -> P4-P1.
```

The displayed basis and that exact mutated basis both fail the complete
projective first-jet necessary condition at 17 and 53:

```text
basis                         prime     charts     solutions     timeouts
displayed signs/order            17        210             0            0
displayed signs/order            53       2550             0            0
P4 -> P4-P1                      17        210             0            0
P4 -> P4-P1                      53       2550             0            0
```

Every empty chart means a unit ideal over the algebraic closure. No candidate
survives two primes, so parameter/scaling reconstruction and held-out section
or whole-fibre validation are not triggered.

## Sign alignment

For each of the ten pairs of the five canonical-height Gram matrices, the
checker exhausts `2^16=65536` relative diagonal signs after fixing one
pairwise sign gauge. It minimizes the scale-fitted relative Frobenius
residual. The joint five-fibre objective is the sum of the ten squared
pairwise residuals.

All ten pairwise minima occur at displayed relative signs and are
simultaneously attainable, so their sum is a rigorous global lower bound for
this finite sign problem and is attained. The joint minimum is
`0.21213241870077748`; any nontrivial pairwise sign raises the joint lower
bound by at least `0.15876042945482075`. A 10% near-optimal band therefore
contains no second relative-sign assignment. Curve pair 356/385 is the
strongest initial pair, with residual `0.1049802750010892`.

The word "relative" matters: a common diagonal sign on all five fibres is
independent inversion of the corresponding generic sections and is a gauge,
not a distinct alignment.

## Fingerprint-bounded permutations

The exact/numerical fingerprint graph uses:

- stable square roots of point-coordinate denominators, including
  integrality;
- exact identity/nonidentity reduction on the supplied models, after
  verifying local minimality at every component prime;
- scale-free absolute height correlations and normalized height rows.

The stable fixed labels are 2, 3, 5, 15, 16, and 17, with modal denominator
roots `1,5,71,1,7,41`. The declared ambiguity graph has respectively

```text
curve                 351     356     376     377     385
perfect matchings       1       1      21       2       1
```

All matchings are scored exactly as declared. Only the identity on each fibre
lies within 10% of its displayed residual against reference curve 356. The
elliptic inversion/group-law replay then verifies all 85 retained rational
points on their canonical short curves.

## One elementary mutation

Signed permutations already cover `P_i -> +/-P_j`. With all other columns
fixed, a replacement column giving a full unimodular basis must contain the
old `P_i` with coefficient `+/-1`; modulo the preceding signed-permutation
gauge, the remaining one-step words are `P_i -> P_i +/- P_j`. The strict
search fixes the six stable anchors, applies one common transformation across
all five fibres, and enumerates

```text
(17-6)*16*2 = 352
```

such shears. Numerical Gram alignment is only a proposal gate. Exactly one
word clears the height and 1% improvement bounds: `P4-P1`. Its joint
objective improves from about `0.21213242` to `0.20795284`; its maximum
changed-height ratio is `1.17813`.

The candidate is then constructed using the exact rational elliptic group
law on all five curves. Every new point lies on its short curve, the basis
matrix has determinant 1, all five fibres are nonsingular, and all 85 point
coordinates are defined at 17 and 53. Only after those exact checks is it
passed to modular elimination.

## Negative controls and conclusion

The signed-residual components

```text
{363,364,378}        {389,390,391}
```

remain separate from the target component and from each other at the pinned
cutoff. They are controls for the fingerprint classification, not candidate
fibres silently folded into the target.

The bounded conclusion is therefore:

> These five records share a numerical/private-pipeline fingerprint, but not
> a rootless-K3 `(8,12;4,6)` realization under the declared section rebasing.

This excludes neither repeated elementary mutations nor a larger bounded
basis search, and it says nothing about unrestricted `GL(17,Z)`, another
family shape, or parameter reductions that are bad/colliding at every tested
prime. It also makes no common-family claim from Gram correlation.

## Reproduction

From the repository root:

```bash
sage -python \
  elliptic-curves/cas/analyze_icarm_wgxli_rank17_rebasing.sage --check

.venv/bin/python \
  elliptic-curves/cas/construct_icarm_wgxli_rank17_bounded_mutation.py --check

sage -python \
  elliptic-curves/cas/eliminate_icarm_wgxli_rank17_first_jet.sage \
  --input \
  artifacts/generated-results/elliptic-curves/icarm_wgxli_rank17_mutation_p4_minus_p1_v1.json \
  --prime 17 --jobs 4 --threads 1 --pair-timeout 60 \
  --work-dir \
  artifacts/local/elliptic-curves/wgxli-r17-first-jet-mutation-p4-minus-p1 \
  --output \
  artifacts/generated-results/elliptic-curves/icarm_wgxli_rank17_mutation_p4_minus_p1_first_jet_mod17_v1.json \
  --check

sage -python \
  elliptic-curves/cas/eliminate_icarm_wgxli_rank17_first_jet.sage \
  --input \
  artifacts/generated-results/elliptic-curves/icarm_wgxli_rank17_mutation_p4_minus_p1_v1.json \
  --prime 53 --jobs 18 --threads 1 --pair-timeout 60 \
  --work-dir \
  artifacts/local/elliptic-curves/wgxli-r17-first-jet-mutation-p4-minus-p1 \
  --output \
  artifacts/generated-results/elliptic-curves/icarm_wgxli_rank17_mutation_p4_minus_p1_first_jet_mod53_v1.json \
  --check

.venv/bin/python \
  elliptic-curves/cas/certify_icarm_wgxli_rank17_bounded_rejection.py --check
```

The literal projective eliminations and their positive controls are
documented in
[`ICARM_WGXLI_RANK17_FIRST_JET_ELIMINATION.md`](ICARM_WGXLI_RANK17_FIRST_JET_ELIMINATION.md).

The generated sign/permutation, exact mutation, mutation mod-17, mutation
mod-53, and aggregate rejection artifacts have SHA-256 hashes

```text
2da02bed8982f94833d88a0fbd9e154cf693b792e3f4f1c99910a300a257b4f7
44a6b1bb54fd84313e194f02511d95fcec16dee0eb03ba23315dfc1adaefd3aa
b1e308be23a6af087f8ff5ea69ae55e0dc3e55e72bfd458f89a30cab5ebeda62
4ede15dca6c6cab88401a1539fbd99dfb372d4c3fb91317b5a985b87d83dbc54
bb223d4b81138cd20c768a80f60ceb0795a49e4b45745a2f990a9b43c7a7791f
```
