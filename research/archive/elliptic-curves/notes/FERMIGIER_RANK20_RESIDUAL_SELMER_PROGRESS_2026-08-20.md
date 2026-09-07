# Fermigier rank-20 residual 2-Selmer progress — 2026-08-20

Status: **research update / computational diagnostics**, not a new rank claim.

## Goal

For the certified rank-at-least-20 Fermigier–Mestre specialization at

```text
u = 28917/20
```

with global minimal model

```text
[1, 1, 1,
 -4437412060110743641525245114305,
 3586842216822165612930264910099076801587288127]
```

determine whether the 2-Selmer information relative to the known 20-dimensional Mordell–Weil subgroup leaves any genuine residual direction that could yield rank at least 21.

The pinned 20-point basis remains the bounded-saturation candidate with SHA-256

```text
6fbdc4367d52ca92cfdfef8b0cc71347b2943784df3780a7a72646b0caff898e
```

and has an exact mod-2 independence certificate of rank 20.

## General-purpose CAS attempts

### eclib / mwrank

The Sage/eclib 2-descent was rerun with 22 auxiliary primes:

```text
R20MWRANK|stage=two_descent|status=start|selmer_only=true|...|n_aux=22
```

It fails before producing mathematical descent information because its quartic search tries to convert a bound of order `3.16e30` to a native C/C++ `long`:

```text
Attempt to convert -0.3157237469e31 to long fails!
2-descent: lower bound ... on c too large
```

Classification: **implementation overflow**, not a rank/Selmer result. Increasing search limits or `n_aux` does not address this failure mode.

### PARI `ellrankinit` / direct `ellrank`

Both

```text
ellrankinit(E)
```

and direct

```text
ellrank(E, 0, known_points)
```

remain stuck for long periods before returning a rank interval. The bottleneck is the cubic number-field/global class-group/unit arithmetic used by the 2-descent, not point-search effort.

### Denis Simon 2-descent

Sage's bundled Simon GP implementation reaches

```text
E[2] = [[0]]
```

and then enters `bnfinit` on the cubic 2-division field. Inspection of `ellQ.gp` confirms that for trivial rational 2-torsion it explicitly executes

```gp
bnf = bnfinit(eqell,1);
rang = ell2descent_gen(ell,bnf,1,help);
```

so this path hits the same global-number-field bottleneck.

### PARI/Sage/Hecke class-group attempts

The cubic field and maximal order themselves are cheap:

```text
Number-field construction: ~0.0015 s
Hecke maximal_order(K):     ~0.108 s
PARI nfinit:                ~0.20 s
```

but all full/global class-group variants remain slow:

- PARI `bnfinit(f,0)`
- PARI `bnfinit(f,1)`
- Sage `K.class_group(proof=False)`
- Hecke `class_group(OK)`
- Hecke `ray_class_group(...; n_quo=2, GRH=true)`

The evidence therefore points specifically to global relation/class-group/regulator/unit computation as the dominant bottleneck, rather than ordinary arithmetic in the cubic field.

## Cubic 2-division field

Using the integral long Weierstrass model, the cubic used for the Kummer representation is

```text
f(x) = x^3
       - 5750886029903523759416717668139307*x
       + 167347710468055045100164888198438918505621536951206
```

The polynomial discriminant has 332 bits. Its bad rational-prime support consists of only 12 primes:

```text
2
3
5
7
13
17
31
79
1049
71889448247
40200713707633
491007790268548705232623905732119
```

There are 25 prime ideals above these rational primes. Prime decomposition is extremely fast (typically sub-millisecond; even the largest bad prime takes only a few milliseconds).

This strongly motivates a custom residual 2-Selmer implementation that avoids a full BNF.

## Kummer valuation skeleton

For each of the 20 known points `P`, compute

```text
alpha_P = x(P) - theta in K*/K*^2
```

and record parity valuations at the 25 bad prime ideals.

Result:

```text
rows=20
columns=25
rank=5
```

Only 7 of the 25 valuation columns occur at all in the known subgroup.

Conclusion: valuation parity alone is far too lossy. Fifteen of the 20 known mod-2 Mordell–Weil directions live in unit/residue/global square-class information.

## Odd local square classes

For every odd bad prime ideal, augment valuation parity with the quadratic character of the local unit part. PARI is used for residue-field arithmetic because Sage's residue-field path triggered a macOS/OpenBLAS `SIGILL` in `matrix_modn_dense_float`/FFLAS.

Result:

```text
valuation_only=5
odd_local_squareclasses=12
known_mod2_rank=20
```

Incremental rank:

```text
[1,2,3,4,5,6,7,8,9,10,11,12,12,12,12,12,12,12,12,12]
```

Thus odd local square classes recover 12 of the 20 known global Kummer directions.

## 2-adic and archimedean square classes

There are two primes above 2, with local metadata

```text
[(e=2,f=1), (e=1,f=1)]
```

Using PARI `nfislocalpower(...,2)` to construct exact local square-class coordinates for the product of the two 2-adic completions gives only two additional independent directions:

```text
odd=12
two_adic=2
odd_plus_2=14
real_places=3
odd_plus_2_plus_real=14
```

So all Selmer-relevant local places together distinguish only 14 of the 20 known global Kummer classes.

This is an important structural result: the remaining six dimensions are genuinely global square-class information, not omitted bad-prime or archimedean data.

## Auxiliary good-prime fingerprints: 20/20 recovered

Small good rational primes were then used only as **witness characters**. They do not impose extra Selmer conditions; they merely distinguish global square classes that are locally identical at the Selmer-relevant places.

Starting from rank 14, a greedy scan found:

```text
q=11   gain 1   rank 15
q=19   gain 1   rank 16
q=23   gain 2   rank 18
q=29   gain 1   rank 19
q=59   gain 1   rank 20
```

Final result:

```text
R20AUXFP|result=faithful_known_kummer_fingerprint
         |dimensions=20
         |auxiliary_primes=[11,19,23,29,59]
```

This is currently the strongest algorithmic progress from the experiment.

We now have a cheap, explicit and faithful 20-dimensional fingerprint of the known Mordell–Weil Kummer image using:

1. odd bad-prime local square classes;
2. the two 2-adic local square classes;
3. real signs;
4. witness primes `11,19,23,29,59`.

All of these computations take milliseconds to low seconds and require no class group or regulator.

The witness primes are **coordinates only**. They must not be interpreted as additional Selmer conditions.

The same run now writes a BNF-free signature-map artifact with all twenty
actual Kummer representatives in ascending power-basis coordinates and their
packed images in the 51 local and 24 witness coordinates.  It is directly
consumable by
[`residual_selmer_quotient.py`](../../../elliptic-curves/cas/residual_selmer_quotient.py), which
reduces future candidate global squareclasses modulo this faithful known
Mordell--Weil image.  This records the coordinate system; it is not a Selmer
upper bound or a class-group computation.

## Explicit BNF-free squareclass candidates from large-prime closures

The retained generators now make relation combinations directly usable as
global squareclasses. The extractor sparsely eliminates the complement of the
declared Selmer set `S`; every retained combination has even valuation at each
prime outside `S`. Thus it is an explicit bounded candidate in `K(S,2)`, with
an exact ascending power-basis representative.

On the one-run `5689:5096` calibration (30,000 sampled principal generators),
the non-S projection kernel has dimension 305. Its 296-dimensional
square-norm kernel has an explicitly verified global-square basis, so the full
norm-compatible span is trivial. The older 15-row support gate had raw
signature rank seven after quotienting by the known rank-20 Mordell--Weil
Kummer image in the fixed 51-local/24-witness coordinate system; this was a
raw diagnostic only. It is **not** a complete calculation of `K(S,2)`, a
residual 2-Selmer basis, a Selmer upper bound, or evidence for seven
additional rational points.

## Current direction: custom mod-2 ideal relation collection

The same collector now has an explicit ERH-certified factor-base gate rather
than relying on a stabilized relation rank.  The cubic field's Bach threshold
is 262,523 rational-prime norm; materializing its 42,251 prime ideals uses no
class-group or regulator routine.  A ten-special, 300,000-generator bounded
calibration produced 483 exact closures, relation rank 27, and an
ERH-conditional S-class quotient model of dimension 42,207 (down from
42,226 before relations). The 23,034 canonical principal relations `(p)` of
the materialized factor base are now stored explicitly; they lower that
ERH-conditional model to 19,204 without any BNF or class-group computation.
Eliminating the non-S columns supplies 464 explicit
generator products; its 455-dimensional square-norm kernel has an explicitly
verified global-square basis. Thus the entire norm-compatible span of the
bounded candidates is trivial. The former individual-row target had 251
products and raw signature rank 11 modulo the known rank-20 image, but was not
a descent filter. This is a conditional, bounded relation collection—not a
complete class computation, local Selmer calculation, or rank conclusion.

Because every generic class-group implementation remains stuck, the current approach is to collect only the information actually needed modulo 2.

A factor base is built from:

- all 25 prime ideals above the Selmer bad primes;
- prime ideals above small auxiliary rational primes.

For small algebraic integers

```text
alpha = a*w0 + b*w1 + c*w2
```

in the maximal order, factor the integer norm, compute exact prime-ideal valuations above its rational factors, and retain the principal-ideal parity relation whenever every odd-exponent prime ideal belongs to the factor base.

Each accepted row is an **exact relation in the ideal class group modulo 2**. The incompleteness is only whether enough relations / a sufficiently complete factor base have been collected.

The first working run used:

```text
factor-base bound = 500
factor-base columns = 188
S columns = 25
coefficient box = 12
```

and immediately produced its first exact rank-gaining relation:

```text
R20REL2|stage=relation|status=rank_gain
       |label=box:0,0,1
       |sampled=1
       |smooth=1
       |rank=1
       |fb_qdim=187
       |S_proj_rank=1
```

The deterministic/random relation collection is currently still running. No completeness or class-group claim should be made until relation-rank stabilization and a certification strategy are established.

## Interpretation and next research steps

The experiments now support a more focused algorithm than generic `ellrank`/`bnfinit`:

```text
cubic field arithmetic (cheap)
    -> construct global K(S,2) generators modulo squares
    -> encode them in the faithful 20-bit Kummer fingerprint
    -> quotient immediately by the known rank-20 image
    -> impose actual local Selmer conditions only on the residual quotient
```

The main remaining global input is the 2-primary `S`-class information and unit square classes. A full class group/regulator appears unnecessarily expensive for this purpose.

Immediate priorities:

1. finish and improve the mod-2 principal-ideal relation collector;
2. replace full integer norm factorization by a factor-base smoothness sieve if norm factorization dominates;
3. grow the factor base and relation set until the observed quotient stabilizes;
4. investigate a rigorous certification of the 2-primary `S`-class quotient without computing the full class group;
5. add compact unit-square representatives;
6. construct the residual `K(S,2)` basis and quotient it using the established 20-dimensional fingerprint;
7. only then build explicit residual 2-covers / search for a new rational point.

## Claim discipline

Nothing in this note changes the public rank status of the curve.

Still proved:

```text
rank(E) >= 20
```

Not yet proved:

```text
rank(E) = 20
rank(E) >= 21
```

The local-square-class and auxiliary-prime fingerprints are exact computations for the known 20-point subgroup. The factor-base relation collector produces exact individual principal-ideal relations, but its completeness is currently heuristic/in progress.
