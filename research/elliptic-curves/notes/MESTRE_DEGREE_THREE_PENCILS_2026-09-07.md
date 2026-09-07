# Degree-three pencils and the missing fibre-component check

The completed audit below now bounds all554 retained pencils by generic
rational rank11. The initial three-case construction is preserved here as
the input to that full retained-portfolio check.

Three explicit divisor constructions on the new Mestre parent u=11 define
rational Jacobian fibrations of old degree3. Their exact generic ranks over
Q are **8,10,10**, and over Qbar are **9,11,11**. The preliminary visible-curve
calculation allowed14 in every case. Complete geometric frame roots close
that gap before any new Weierstrass equation compilation or point search.

These are constructions on an existing new parent, not three new parents.
They do not improve its current MW11 presentation. No high-rank fibre,
inventory addition, pairwise inequivalence or point-search result is claimed.
Authority: `EC-MESTRE-DEGREE-THREE-PENCIL-GATE-20260907` in
[`MATH_STATUS.json`](../../MATH_STATUS.json).
The [report](../../artifacts/generated-results/elliptic-curves/mestre_degree_three_gate_v1.json)
binds all stages and the independent proof.

## Why degree three, and why a bisection?

The [complete degree-two gate](MESTRE_DEGREE_TWO_RANK_GATE_2026-09-07.md)
excludes generic MW rank above13 at old degree2 on all six determinant468
parents. Their full rational NS rank is18, with an old split I4 fibre.

For any new nef fibre class D with D.F_old=3, its four nonnegative integral
intersections with the components of that I4 sum to3. At least one component
is therefore vertical for the new fibration. Its square-minus-two class is
nonzero modulo D, so Shioda–Tate gives

```
generic MW_Q <= 18 - 2 - 1 = 15.
```

Thus a hypothetical MW16 fibration on these parents must have old degree
at least4. This says nothing about specialized elliptic ranks.

A fibre D=O+P+Q made of three distinct old sections has an additional
obstruction. For its three components to have zero intersection with D,
their pairwise intersections are all1. They give an A2 root lattice in
that fibre. The old I4 component above is a distinct curve in another new
fibre, so adds a further independent root. Any rational Jacobian fibration
of this form has generic MW rank at most13, independently of section-word
size. The same conclusion holds after translating the three old sections.

The bounded initial roster considered1,562 signed words of support at most3.
It found200 sections with O.P=1 and1,726 pairs admitting a visible rational
section of the resulting pencil. Their visible MW upper bounds range6..12.
These are a bounded roster, not all degree-three classes. The general
three-section obstruction above made compiling its best row unnecessary.

Instead use a rational bisection R with O.R=2. Since R is a smooth rational
curve, R^2=O^2=-2. Then D=O+R has square0 and zero intersection with each
of its components. Every other irreducible curve intersects it nonnegatively,
so D is nef. A rational curve C with D.C=1 proves primitivity and supplies
a section. The standard K3 linear-system theorem then gives an elliptic
pencil over Q. Its reducible fibre has root rank1, allowing a possible
generic rank14 after the old I4 contribution.

## Exact rational bisections and translation classes

Start from the previously compiled
[degree-two chord fibration](MESTRE_EXPLICIT_CHORD_FIBRATIONS_2026-09-07.md).
Its four rational quartic points come from the two old I2 fibres. Using the
first as origin gives two generators K,Q. The retained height-six domain
contains56 old sections of degree1 on this alternate fibration; the first
three in the frozen word ordering are added as further generators.

All50 signed support-one-or-two words in these five points are frozen before
curve inspection. Exact group addition and the retained birational inverse
produce two degree-two rational curves, rows27 and35. Their inverse chord
coordinate is identically the alternate base parameter z. They are sections
of that alternate fibration, hence smooth rational square-minus-two curves
on the same K3. Their degrees and nontrivial deck actions are checked exactly.
They have old O intersections1 and0 respectively.

The zero and22 signed old basis translations give46 further exact curve
parametrizations. None has O intersection2. This calculation nevertheless
determines each base bisection's intersection with every old basis section.
There are only90 possible nonidentity component-intersection tuples: each
old I2 count lies between0 and2, and the three old I4 counts sum to at most2.
In each case exactly one tuple yields an integral NS class of square-2.

The two classes, in the existing rational NS basis, are

```
R27 = [3,2,0,0,-1,-2,-1,1,-1,0,0,0,0,0,0,0,0,0]
R35 = [1,1,0,0,-1,-2,-1,1,0,0,0,1,-1,0,0,0,0,0].
```

This enables class selection instead of further blind coordinate trials.
For R-S, its intersection with O equals R.S, and its intersection with
another old section P equals R.(P+S). Applying these exact identities to
the6,340 retained signed height-at-most-six words produces554 O+R-S pencils.
Three have visible MW upper bound14. They are selected before full frame
enumeration. The retained domain's separate enumeration is selection data;
the independent proof here does not certify its completeness, and no complete
classification of all degree-three pencils is asserted.

## Complete frame calculation

The [full geometric NS Gram matrix](MESTRE_FULL_NS_GRAMS_2026-09-07.md)
has rank19 and determinant468, with an explicit Galois involution whose
fixed rank is18. For each selected D and its visible rational section C,
the span of D,C+D is a unimodular hyperbolic plane. Its integral orthogonal
complement is the full rank17 frame. Every norm2 vector in the positive
negative-frame Gram matrix is a geometric vertical root.

| Selected row | Base bisection | Visible MW upper bound | All roots | Geometric / rational root ranks | Exact MW Q / Qbar |
|---|---:|---:|---:|---|---|
| 1 | 27 | 14 | 20 | 8 / 8 | 8 / 9 |
| 2 | 35 | 14 | 14 | 6 / 6 | 10 / 11 |
| 3 | 35 | 14 | 14 | 6 / 6 | 10 / 11 |

The [standalone checker](../cas/verify_mestre_degree_three_gate.sage) uses
rational LDL decomposition and exact recursive integer intervals to enumerate
every vector of norm at most2. It does not invoke the producer's PARI
enumerator. It also checks the unimodular hyperbolic splitting, complete frame,
Galois action, exact bisection equations and inverse maps, both unique NS
classes, and all three translated divisor classes. The first translated
bisection's full coordinates and all eleven old basis intersections are
checked separately. The three ranks follow from Shioda–Tate and the fixed
root spaces; no conjectural Picard bound or specialization score enters.

The [portable bundle](../../artifacts/generated-results/elliptic-curves/mestre_degree_three_proof_bundle_v1.json)
contains the inputs. The copied checker and bundle pass in a fresh directory
without repository imports. All eight supervised stages total
36.553333365824074 seconds, including1.560442331014201 seconds for that proof.
Two initial launcher API errors occurred before any CAS worker started;
there was no arithmetic rerun. Source review, launcher and report CPU are
outside that timing. Every mathematical stage has a protocol, time/memory
limits and retained output under the corresponding local `mestre-*` directory.

## Consequence for the next search

On these explicit parents with certified full NS data, a visible-component
upper bound is insufficient admission evidence for a desired generic rank.
After proving a proposed divisor is a primitive nef pencil with a rational
section, compute the complete geometric frame roots and their Galois-fixed
span before spending effort on coefficients and section bases for that rank
target. Lower-rank fibrations can still merit a separately motivated height
or visibility experiment; these generic bounds do not exclude their fibres.

The three selected constructions fail the intended MW14 improvement. Their
new Weierstrass equations and MW bases remain unconstructed; no larger
translation, parameter or point-search wave is queued.

```sh
sage -python elliptic-curves/cas/verify_mestre_degree_three_gate.sage --input \
  artifacts/generated-results/elliptic-curves/mestre_degree_three_proof_bundle_v1.json
python3 elliptic-curves/cas/report_mestre_degree_three_gate.py --check
```
## Completed retained portfolio: no generic improvement over MW11

The subsequent [full retained-portfolio audit](../../artifacts/generated-results/elliptic-curves/mestre_parent_coverage_followup_v1.json)
closes the other200 retained candidates with visible ceilings12 or13.
Their exact generic rational ranks are:

| Rank | Count |
|---:|---:|
|6|5|
|7|28|
|8|69|
|9|71|
|10|27|

The corresponding geometric ranks are one larger. Independent rational
LDL root enumeration checks all200 frames on the same certified geometric
NS lattice. The three earlier ceiling14 cases below are reused, not rerun.
A separate standalone audit recomputes all554 translated divisor classes,
their visible vertical roots and rational degree-one witnesses. The
remaining351 have visible upper bounds at most11. Consequently **none of
these554 retained pencils has generic rational rank greater than11**.

All554 fibre classes are distinct in the marked NS lattice. They need not
be inequivalent under surface automorphisms. This result neither classifies
all old-degree-three fibrations nor certifies completeness of the retained
height domain. It gives no specialized-rank upper bound or search-yield
comparison. It closes this retained portfolio before further coefficient
construction; no point boxes, new parents or inventory additions follow.

The initial standalone extension failed because partial intersections no
longer uniquely identified every translated class. The original verifier
and failed replay are preserved. Version2 uses the exact rotation of split
fibre components under translation, checks the original rational bisection
coordinates, and passes on the unchanged frame bundle. The producer was
not rerun. This does not invalidate the original three-case proof.

Authority: `EC-MESTRE-RETAINED-DEGREE-THREE-PORTFOLIO-20260907` in
[`MATH_STATUS.json`](../../MATH_STATUS.json).
The additional stages cost **14.940146973 supervised seconds**, including
the1.412671109-second failed verification attempt. The portable inputs are
the [200-frame bundle](../../artifacts/generated-results/elliptic-curves/mestre_retained_pencil_proof_bundle_v1.json)
and [554-class bundle](../../artifacts/generated-results/elliptic-curves/mestre_retained_visible_proof_bundle_v1.json).

```sh
sage -python verify_mestre_degree_three_gate_v2.sage --input mestre_retained_pencil_proof_bundle_v1.json
sage -python verify_mestre_retained_visible_gate.sage --input mestre_retained_visible_proof_bundle_v1.json
python3 elliptic-curves/cas/record_parent_coverage_followup.py mestre --check
```
