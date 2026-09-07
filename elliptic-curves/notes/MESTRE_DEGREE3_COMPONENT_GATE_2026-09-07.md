# Mestre degree-three component bounds

Authority: `EC-MESTRE-DEGREE3-COMPONENT-GATE-20260907`.

For each of the six determinant-468 Mestre surfaces, every Q-Jacobian
fibration of degree three relative to the original fibre has arithmetic
generic MW rank at most **15**. If its fibre is a triangle formed from
three old rational sections, the bound improves to **13**.

These are geometric bounds for the entire stated constructions, not
bounded word-search results. Thus an MW14–15 constructor on this cohort
cannot come from an old-section triangle, and an MW16 constructor must
have old degree at least four. The latter also uses the existing
[complete degree-two MW13 bound](MESTRE_DEGREE_TWO_RANK_GATE_2026-09-07.md).
A degree-one new fibration meets at most one of the four old I4
components, leaving three independent vertical components and giving
the same rank-at-most-13 bound.

## Untouched fibre components

The [full NS proof](MESTRE_FULL_NS_GRAMS_2026-09-07.md) supplies arithmetic
Picard rank 18 and a split I4 fibre on every surface. Its four rational
irreducible components are

`C0=F-I4_1-I4_2-I4_3, C1=I4_1, C2=I4_2, C3=I4_3`.

They sum to the old fibre F and have the affine A3 intersection matrix.
For a nef new fibre D with `D.F=3`, the nonnegative integers `D.Ci`
sum to three. At least one Ci, say C, is therefore vertical for D.
It is not a multiple of D, since `F.C=0` and `F.D=3`. Consequently the
rational fibre-root space has rank at least one, and Shioda–Tate gives

`rank MW <= 18-2-1 = 15`.

Now suppose `D=S0+S1+S2`, where the three old rational sections meet
pairwise once. Since `D.C=0`, each `Si.C=0`. The four curves have Gram
matrix `affine A2(-1) + <-2>`, of rank three and with radical D.
Modulo D they supply `A2(-1)+A1(-1)` in the rational fibre-root space.
Thus

`rank MW <= 18-2-3 = 13`.

This argument applies to every such three-section triangle, with no
choice of a fixed anchor and no bound on section words or heights.
The rational section required for a Jacobian fibration is part of the
claim's hypothesis; no zero section or equation is manufactured here.

The [certificate](../../artifacts/generated-results/elliptic-curves/mestre_degree3_component_gate_v1.json)
checks the four component classes on all six NS matrices, all twenty
possible nonnegative intersection profiles of total degree three, and
the root-rank calculation. Rationality and effectivity of the old
components remain inputs from the existing full NS proof.

```sh
sage -python elliptic-curves/cas/verify_mestre_degree3_component_gate.sage
```

The exact check and saved-certificate round trip pass under a 30-second,
one-worker cap using SageMath 10.9. No parameter search is needed.

## Exploratory work retained separately

The degree-three obstruction motivated a finite trial of nontriangle
isotropic classes. Its first 256 eligible classes have only one positive
old root orthogonal modulo three. Nevertheless their **full abstract
frames**, including all geometric roots and their Galois-invariant span,
allow arithmetic ranks only 6–11 if realized. No nef representative or
elliptic equation is claimed for these classes.

A subsequent bounded three-step two-neighbour trial first minimizes
positive root count and then rational root rank, and checks sixteen
complete geometric frames at each step. None of these 48 candidates
allows arithmetic rank above 11. This is not a bound for all neighbours,
all fibrations on these surfaces, or the broad parent problem.

The exact scripts, frame data and terminal logs are retained under
`artifacts/local/elliptic-curves/`, named
`mestre-degree3-one-old-root-pilot.*`, `mestre_degree3_one_old_root_v1.json`
and `mestre-root-descent-pilot.*`. These exploratory calculations are not
dependencies of the component-bound theorem.

The construction target remains open: an explicit MW14–16 family with
its full generic coordinate basis and exact rational specialization to
302. Nontriangle degree-three MW14–15 fibrations, degree-four-and-higher
fibrations, and other surfaces remain eligible.
