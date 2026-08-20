# Fermigier rank-20 residual 2-Selmer / higher-descent experiment

## Target

The pinned Fermigier–Mestre specialization is

```text
u = 28917/20
E = [1, 1, 1,
     -4437412060110743641525245114305,
     3586842216822165612930264910099076801587288127]
```

The repository already has an exact rank-at-least-20 independence certificate,
an exact 20-dimensional mod-2 certificate for the chosen basis, bounded
small-prime saturation work, and independent PARI/eclib descent runners.

The goal here is to compute the **2-Selmer quotient relative to the known
rank-20 subgroup**, expose the surviving covering spaces, and continue only
arithmetically plausible classes through 4- and optionally 8-descent.

## Why this is the residual computation we want

Magma's

```magma
TwoDescent(E :
    RemoveGens := SequenceToSet(known),
    RemoveTorsion := true
);
```

computes locally soluble 2-coverings after factoring out the subgroup generated
by `known` and the torsion image. Since the repository separately proves that
the 20 supplied points are independent modulo 2, this is the quotient that
matters.

If it returns `n` nonzero coverings then `n + 1 = 2^d`, where `d` is the
residual 2-Selmer dimension.

For this experiment `d = 0` is decisive: the Selmer upper bound and the exact
rank >= 20 certificate meet at rank 20.

## What this version improves

The older repository job had the right first step but did too much downstream
work and did not map a 4-cover hit back to E.

This version adds:

1. full Cassels--Tate pairing on the residual 2-Selmer quotient;
2. `FourDescent(... RemoveGensEC := known)` so the known rank-20 directions
   are removed again at the 4-cover level;
3. pairing every 4-cover against all residual 2-covers before point search;
4. `AssociatedEllipticCurve(F : E := E)` to map `PointsQI` hits back to E;
5. `IsLinearlyIndependent(known cat [P])` to distinguish a rank-21 candidate
   from a saturation/index point;
6. optional 8-descent only on covers that survive the earlier filters.

## First run

From the repository root:

```bash
python3 elliptic-curves/cas/build_fermigier_rank20_residual_selmer.py \
  --output artifacts/local/elliptic-curves/fermigier_rank20_residual_selmer.m \
  --point-bound 100000

mkdir -p artifacts/local/elliptic-curves

caffeinate -i magma -b \
  artifacts/local/elliptic-curves/fermigier_rank20_residual_selmer.m \
  2>&1 | tee artifacts/local/elliptic-curves/fermigier_rank20_residual_selmer.log
```

Do not add `SetClassGroupBounds("GRH")` to the first run. A GRH-assisted repeat
can be useful for exploration, but keep it separate from an unconditional
certificate.

## Interpret `R20SEL|` output

`S0_EXACT_RANK20`: no residual 2-Selmer class remains. Exact rank 20.

`S1_EXACT_RANK20_CT`: residual classes exist but the Cassels--Tate pairing has
zero radical. They are Sha obstruction rather than a missing MW direction.
Exact rank 20.

`S2_SATURATION_POINT`: a rational point is found on a surviving 4-cover, but
its image on E is dependent with the 20 points. Enlarge the rank-20 lattice and
rerun. Do not call it rank 21.

`S3_RANK_AT_LEAST21_CANDIDATE`: a point maps to E and the 21-point set is
independent. Immediately run the repository's exact finite-quotient
independence certificate and replay the minimal model + exact conductor.

`S2_UNRESOLVED_HIGHER_DESCENT`: a 4-cover survives the 4x2 pairing but no point
was found at the current `PointsQI` bound. Continue on that cover, not with a
raw x-search on E.

`S4_INCOMPLETE_HIGHER_DESCENT`: a radical 2-class returned no 4-cover or hit the
configured 4-cover cap. Treat this as computationally incomplete, not as a rank
bound.

## If 4-cover survivors remain

Increase the point bound geometrically:

```bash
python3 elliptic-curves/cas/build_fermigier_rank20_residual_selmer.py \
  --output artifacts/local/elliptic-curves/fermigier_rank20_residual_selmer_h1m.m \
  --point-bound 1000000
```

Then consider `10^7` only if the same 4-cover remains unresolved.

`PointsQI` is searching an intersection of two quadrics. A huge-height point on
E can have a much smaller representative on a reduced 4-cover.

## 8-descent

Only after an all-zero 4x2 Cassels--Tate row and no `PointsQI` hit:

```bash
python3 elliptic-curves/cas/build_fermigier_rank20_residual_selmer.py \
  --output artifacts/local/elliptic-curves/fermigier_rank20_residual_selmer_8.m \
  --point-bound 1000000 \
  --eight-descent
```

The 4-cover/2-cover Cassels--Tate pairing is used first because it can eliminate
covers before the much heavier 8-descent.

## Independent diagnostics

eclib/mwrank:

```bash
PYTHONUNBUFFERED=1 caffeinate -i \
  sage -python elliptic-curves/cas/run_fermigier_rank20_mwrank_descent.py \
  --n-aux 22 \
  2>&1 | tee artifacts/local/elliptic-curves/fermigier_rank20_mwrank_naux22.log
```

PARI:

```bash
PYTHONUNBUFFERED=1 caffeinate -i \
  sage -python elliptic-curves/cas/run_fermigier_rank20_pari_descent.py \
  --efforts 0,1,2,4,8 \
  2>&1 | tee artifacts/local/elliptic-curves/fermigier_rank20_pari_descent_e8.log
```

These are useful independent upper-bound computations, but they do not replace
the explicit residual-cover workflow.

## Research direction if a class survives

For each surviving 2/4-cover:

1. retain its exact quartic / quadrics from the `R20SEL` log;
2. minimise and reduce the covering model;
3. use Cassels--Tate to discard obstructed classes;
4. search only survivors with `PointsQI`;
5. if unresolved, use 8-descent;
6. only then consider custom lattice enumeration on that specific reduced model.

Relevant background: Cremona--Fisher--Stoll on minimisation/reduction of
2-, 3- and 4-coverings; Fisher/Schaefer/Stoll on the Cassels--Tate pairing; and
Magma's higher 2-power descent implementation.

## Rank versus saturation

The existing basis has bounded small-prime saturation work, but was not
globally proved saturated because no exact rank upper bound was known.

Keep these outcomes distinct:

- trivial residual quotient -> exact rank 20;
- zero Cassels--Tate radical -> exact rank 20;
- mapped dependent point -> saturation/lattice improvement;
- mapped independent 21st point -> genuine rank growth.
