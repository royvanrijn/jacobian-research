# NS0024 MW4 modular recovery frontier — 2026-09-01

<!-- status-consumer: EC-K3-NS0024-DIRECT-QQ-INOSE-OBSTRUCTION e87afc1b3529a07f -->
<!-- status-consumer: EC-K3-NS0024-QQ-MARKING-OBSTRUCTION b7f0cf002c0411fe -->

> **Source decision (2026-09-04).** This semistable reconstruction is parked
> as geometric NS0024 work. A full rational NS0024 marking over `QQ` is
> impossible by the Fricke-quotient obstruction; see
> [`NS0024_QQ_MARKING_OBSTRUCTION_2026-09-04.md`](NS0024_QQ_MARKING_OBSTRUCTION_2026-09-04.md).
> It cannot close the arithmetic MW17 milestone over `QQ(t)`.
> The resolved exporter now also supports `--split-section-opens`.  Giving
> each chart condition its own inverse lowers the representative explicit-
> center, one-hyperplane input from about 86 KB to about 59 KB and removes the
> high-degree grouped chart products.  Exact 55-second probes at `p=11,101`
> and `65521` did not finish; that is a performance boundary, not evidence of
> nonexistence.

## Outcome

The requested marked MW4 source is **not yet recovered**.  No
characteristic-zero reconstruction was attempted, and the q4/orbit1 edge-1
compiler was not launched without an exact marked point or family input.

The remaining gap is nevertheless substantially narrower and the modular
search is cheaper:

1. all 400 minimum-pole bases were re-scored by resolved-component depth;
2. among the 120 bases containing the exact q4/orbit1 projection, a depth-13
   basis improves the previous depth-12 choice;
3. its third polynomial section has profile `(4,2,0)`, so the depth-three I7
   and depth-two I5 congruences determine its entire `X` polynomial;
4. exact tangent saturation reduces a representative fixed-P4 quotient from
   dimension 288 to 16;
5. the declared rational MW3 charts were exhaustively enumerated at p=11 and
   p=13;
6. at p=11 every P4 over `GF(11^2)` above every rational resolved-basis MW3
   marking was exhausted, with no exact MW4 completion;
7. the replacement joint ideal leaves the surface and all MW3 coordinates
   algebraic over `GF(p)`, so its closed points may carry the entire marking
   over `GF(p^d)` rather than extending only P4.

The exact compact ledger is
`artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-modp-census.json`.

## Arithmetic realizability obstruction

The foundry proves a geometric K3/Neron--Severi realization and a geometric
MW4-to-MW17 route. It does **not** prove

```text
NS(X)=NS_Q(X),
rank MW(Q(t))=4 at the source,
or rank MW(Q(t))=17 at the rootless target.
```

The later rational-marking obstruction proves that the first equality cannot
hold for geometric `NS0024` over `QQ`; hence full source rank four and full
rootless arithmetic rank seventeen cannot both descend as the required
rank-19 marking. This distinction remains essential for number-field and
geometric applications. A
geometric Mordell--Weil basis may carry a nontrivial Galois action and need
not contribute its full rank over `QQ(t)`.

Every joint closed-point extraction now carries an early arithmetic gate.
The verifier computes three Frobenius orbit sizes separately: the unmarked
surface coefficients `(A,B)`, the four sections, and the full resolved marking
including the chosen `I7/I5` orientation roots.  Only a Frobenius power fixing
`(A,B)` acts on the MW lattice of one surface.  When such a power exists, the
verifier expresses the four conjugate sections in the certified MW4 basis,
checks the resulting integral matrix against the exact height Gram of
determinant `95/14`, and records `rank MW4^(Frob=1)`.  Any extra field needed
only by the orientation roots is recorded without lowering the section rank.

If the unmarked surface itself has orbit size greater than one, the factor is
recorded as inconclusive for rational descent: it is a closed orbit of
conjugate surfaces, not a nontrivial Frobenius action on one prime-field
surface.  If the surface is prime-field fixed but the fixed MW rank is below
four, that is warning evidence against full rational descent of that
realization.  Repeated full fixed rank across good primes is positive evidence
only; neither pattern proves `NS=NS_Q` without a common characteristic-zero
producer and a descent argument.

The same gate would have been required at the rootless endpoint with fixed
rank seventeen. It is now closed negatively for `NS0024`, rather than merely
unverified.

## Improved exact basis

The old objective maximized labelled-node incidences and then minimized total
pair intersection.  It left four tied optima, all with resolved-component
depth 12.  Re-enumerating every minimum-pole basis containing q4/orbit1 gives
the following depth-13 recommendation:

| section | MW quotient coordinate | `P.O` | `(I7,I5,I4)` |
|---|---:|---:|---:|
| Q1 | `(-1,0,1,0)` | 0 | `(6,0,0)` |
| Q2 = q4/orbit1 | `(-1,0,0,0)` | 0 | `(2,1,1)` |
| Q3 | `(-1,-1,0,0)` | 0 | `(4,2,0)` |
| Q4 | `(0,0,0,1)` | 1 | `(6,4,3)` |

Its section-intersection Gram matrix is

```text
-2  1  1  1
 1 -2  0  2
 1  0 -2  2
 1  2  2 -2
```

For Q3, write `center^2=-A/3`.  The congruences

```text
X_Q3 = center mod t^3,
X_Q3 = center mod (t-1)^2
```

give five linear conditions on a degree-four `X`; hence no free
`X_Q3` coefficient remains.  This is the main mathematical reduction.

The exact certificate and all four old optima plus the new recommendation are
emitted by
`scripts/certify_lattice_foundry_ns0024_mw4_basis.sage`.

## Algorithmic reductions

### Resolve component charts before elimination

For a split `I_n` component `k`, use depth `d=min(k,n-k)` and substitute

```text
X = center mod s^d
```

before forming section-square equations.  In particular, the old all-node
sections of component two have `x1=-a1/6`; carrying `x1` as a free variable
left the singular tangent of the Weierstrass node in the global system.

### Saturate exact component-one tangents

Node incidence alone includes higher-depth and embedded branches.  Inverting
the first exceptional tangent at I7, I5 and I4 changed the test quadratic
fixed-P4 algebra as follows:

| formulation | variables | input size | quotient dimension |
|---|---:|---:|---:|
| fully recursive Y | 9 | about 11 MB | 108 on the earlier slice |
| sparse cubic Y, unsaturated | 14 | about 7 KB | 288 |
| sparse cubic Y, three exact tangent inverses | 17 | about 7 KB | 16 |

The saturated test Gröbner basis took about seven seconds instead of growing
past multiple gigabytes.

### Joint extension-defined MW3 formulation

The rational-MW3 staging was useful for a bounded census but is structurally
too restrictive for the next search.  The replacement resolved-depth13 ideal
solves the surface, Q1, Q2, forced Q3, and Q4 jointly over `GF(p)`.  A closed
point of degree `d` therefore places all four sections over one residue field
`GF(p^d)`.

The exporter also supports `--explicit-formal-centers`.  This retains the five
`I5` and four `I4` formal-center jet coefficients as auxiliary variables with
their exact quadratic recurrences.  It is equivalent on the existing
`r1 != 0`, `ri != 0` chart but avoids recursively expanding high powers of
`r1_inverse` and `ri_inverse`.  On the pinned one-hyperplane `p=11` system it
changes the export from 46 variables, 56 equations and about 1.15 MB to 55
variables, 65 equations and about 86 KB; the largest exact-order saturation
row drops from about 295 KB to under 2 KB.  These are implementation-size
measurements, not mathematical evidence for a point.

The optional `--split-section-opens` encoding separates the three Q1, three
Q2, three Q3, and six Q4 nonvanishing conditions.  With the same explicit
centers, one hyperplane, and fixed RUR anchor it uses 66 variables and 76
equations but only about 59 KB of input, versus 55 variables, 65 equations and
about 86 KB for the grouped formulation.  The extra variables are linear
inverse witnesses; the benefit is that no section-chart saturation row is a
large product.  Bounded `msolve` probes in characteristics 11, 101, and 65521
did not close within 55 seconds, with peak memory between 1.3 and 3.2 GB.
This validates the leaner encoding and records its current compute boundary;
it does not assert emptiness, dimension, or a closed point.

The intended dimension count is two after MW3 and one after MW4.  Hence the
first exact-point attempt uses one generic affine hyperplane in the eight
surface variables on the joint MW4 locus.  Two slices belong only to an
MW3-only zero-dimensional diagnostic; taking them before adjoining Q4 would
generically miss the MW4 curve.

The joint formulation also imposes:

* exact I7/I5/I4 orders through the first unmatched formal-branch coefficient;
* exact exceptional tangents for every nonidentity component;
* explicit X/Y open covers for the identity components of Q1 and Q3;
* the genuine Q4 pole and `c != 0,1`;
* sparse low-degree Y coordinates rather than recursively expanded inverse
  expressions.

On the all-Y identity chart at p=11, one deterministic hyperplane gives 45
variables and 55 equations.  The export is about 1.1 MB.  A capped 120-second
`msolve -g 1` probe did not finish; this is a performance observation, not a
negative result.  The earlier recursive form was abandoned after expression
swell, before elimination.

As an exact positive control, fixing the MW3-only ideal to the known resolved
p=11 marking with `(r1,ri)=(3,9)` gives a nonunit zero-dimensional ideal and a
32-element Gröbner basis.  Its forced Q3 is
`X=(1,10,10,8,7)`, `Y=(0,0,0,5,0,7,10)`.  Thus the new equations and the
all-Y open chart contain an independently enumerated resolved MW3 point; this
does not supply Q4.

The earlier bounded scanner remains:

```text
split-node chart -> exact MW3 census -> fixed-MW3 P4 algebra -> exact marking
```

The scanner now supports the original and resolved-depth bases, first-hit and
all-triples modes, and full rational P4 orientation gates.  Cyclic variable
rotation is an exact fallback for the few finite triangular decompositions
whose default order is slow.

## Exact finite-field results

### Rational MW3 census

| prime / basis | marked triples | surfaces | split-node charts |
|---|---:|---:|---:|
| p=11, original | 4 | 4 | 4 |
| p=11, depth-13 | 9 | 9 | 9 |
| p=13, original | 14 | 14 | 11 |
| p=13, depth-13 | 16 | 16 | 13 |

The depth-13 basis therefore exposes more rational modular points at both
primes.  The p=11 and p=13 full rational-P4 scans found no exact completion.

### Complete p=11 `GF(11^2)` P4 census above rational MW3

There are nine p=11 rational depth-13 MW3 markings.

* For `c in GF(11) \ {0,1}`, 81 fixed-c algebras were solved over
  `GF(11^2)`.
* For `c not in GF(11)`, every one of the 55 monic irreducible quadratic
  minimal polynomials was imposed on each MW3 marking: 495 algebras.
* All 576 algebras were decoded exactly.  Every reconstructed section passed
  its Weierstrass identity before component and intersection classification.
* Exact MW4 markings found: **0**.

Thus the following bounded negative is complete:

> In the declared split chart, no depth-13 marked MW4 point over `GF(11^2)`
> lies above an MW3 marking rational over `GF(11)`.

This does not exclude an MW3 marking that is itself defined only over an
extension field.

## Remaining gap and next attack

The immediate target is one exact marked source point, not a positive-dimensional
family.  Any such point must be replayed through all four section identities,
absolute component labels, and the full Gram matrix and then sent directly to
the q4/orbit1 edge-1 compiler.  Only after that yes/no equation-side test should
the containing family be recovered.

The cheapest source-point search not covered by the census is now a closed
point of the joint resolved-depth13 ideal, rather than another P4 scan above
rational MW3 points.  Two viable exact routes are:

1. eliminate a one-hyperplane joint MW4 slice and decode every irreducible
   eliminant factor as a Frobenius orbit over `GF(p^d)`;
2. first recover the positive-dimensional MW3 ideal, then solve Q4 over its
   function/residue fields without specializing MW3 to a rational point.

The first exact marked point should invoke
`compile_lattice_foundry_ns0024_edge1_modp.sage` immediately.  A successful
isolated-point compilation is still only a modular equation-side test; the
one-dimensional marked component, a second-prime match, and characteristic-zero
reconstruction remain later gates.

The arbitrary-degree extractor follows the prior Q80 common-producer method.
With `--fixed-rur-anchor`, the joint exporter adds one fixed integral linear
form involving the surface and all four sections.  The RUR decoder requires a
squarefree eliminant whose degree equals the quotient dimension, substitutes
every coordinate back into the original joint equations, factors the
eliminant over `GF(p)`, and treats each irreducible factor separately as one
closed point/Frobenius orbit over `GF(p^d)`.  It then invokes the independent
joint source verifier and, on the first exact marking, automatically runs the
adapter and edge-1 compiler.  It never guesses `d`, enumerates a whole
`GF(p^d)`, or places unrelated factors in a common splitting field.

## Reproduction entry points

```bash
sage elkies-k3/scripts/certify_lattice_foundry_ns0024_mw4_basis.sage
c++ -O3 -std=c++17 -o /tmp/ns0024-mw4-seed \
  elkies-k3/scripts/search_ns0024_mw4_marked_seed_modp.cpp
/tmp/ns0024-mw4-seed 11 3 3 resolved-all
/tmp/ns0024-mw4-seed 13 3 3 resolved-all
sage -python elkies-k3/scripts/solve_lattice_foundry_ns0024_p4_from_mw3_modp.sage --help
sage -python elkies-k3/scripts/extract_lattice_foundry_ns0024_p4_gb_solution.sage --help
sage -python elkies-k3/scripts/recover_lattice_foundry_ns0024_mw4_family_resolved_modp.sage \
  --prime 11 --surface-hyperplane 1,2,3,4,5,6,7,8,9 \
  --fixed-rur-anchor --explicit-formal-centers \
  --export-msolve /tmp/ns0024-joint-p11.ms
msolve -t 32 -v 2 -g 1 -m 4096 \
  -f /tmp/ns0024-joint-p11.ms -o /tmp/ns0024-joint-p11.rur
sage -python elkies-k3/scripts/extract_lattice_foundry_ns0024_joint_gb_point.sage --help
sage -python elkies-k3/scripts/extract_lattice_foundry_ns0024_joint_rur_point.sage --help
python3 elkies-k3/scripts/summarize_lattice_foundry_ns0024_modp_census.py --help
```

Once certified joint-point artifacts exist at several primes, aggregate their
Frobenius gates with repeated
`--joint-point LABEL artifacts/generated-results/<point>.json`.  The summary
reports the prime-field fixed-rank histogram and keeps its evidentiary boundary
separate from characteristic-zero descent.

The RUR extractor audits every irreducible eliminant factor as a closed point,
including rejected exact-marking factors and the arithmetic-realizability
summary of every accepted factor.  It selects the first accepted factor only
for the existing single-point edge-1 handoff; that selection does not truncate
the factor census.

The long census commands, directory hashes, exact counts, and proof boundary
are recorded in the compact generated ledger cited above.
