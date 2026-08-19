# New-family hidden generic basis and half-lattice experiments

This note records the follow-up work on the six-root quartic family with root set

```text
(-47,-43,-31,30,45,46)
```

It is deliberately separate from the earlier construction/classification notes.  The conclusions below distinguish proved subgroup statements from exploratory specialization searches.

## 1. Starting point

The quartic construction supplies twelve automatic sections.  Taking one as origin gives eleven displayed generic Mordell--Weil directions.  At `T=11`, PARI with the eleven mapped automatic points proves the specialization has rank exactly 11.

The fixed six-root K3 surface has finite-minimal coefficient degrees 8 and 12, Euler number 24, reducible-fiber root rank 3, and hence a Shioda--Tate geometric Mordell--Weil rank ceiling 15.  Thus this architecture cannot have generic rank 17 or 18.

## 2. PARI basis exposes a large index defect

At `T=11`, feeding the exact rank-11 curve to PARI returns a much smaller rank-11 basis.  Combining the automatic and PARI points preserves rank 11 but changes the regulator from

```text
automatic regulator = 1.1277251689020845147228811242464072157619384976349336281114516323e16
union regulator     = 7.4686286774084602112431070820914852338495419928427115035791860166e7
```

so

```text
regulator ratio      = 150994944 = 12288^2
sqrt regulator ratio = 12288
```

and therefore the PARI/union subgroup enlarges the automatic subgroup by index

```text
12288 = 2^12 * 3.
```

Exact relations recover the eleven union generators in automatic coordinates with denominators

```text
12, 4, 6, 4, 12, 2, 12, 12, 6, 2, 6
```

and the rational change-of-basis determinant is exactly

```text
1 / 12288.
```

## 3. Saturation at T=11

The union basis was tested with eclib/mwrank saturation restricted to the small primes relevant to the observed index.

```text
2-SATURATION RESULT = (True, 1, '[ ]')
UP-TO-3 SATURATION RESULT = (True, 1, '[ ]')
```

Hence the rank-11 union subgroup is proved 2-saturated and 3-saturated at `T=11`.  Full unrestricted saturation was not attempted again after earlier memory failures; no claim is made here that this proves the full Mordell--Weil group is saturated at every larger prime.

For the half-lattice experiment the 2-saturation statement is the relevant one: the 2048 classes in `L/2L` are the correct parity classes for this specialization.

## 4. The divided points lift generically

The same eleven divisibility relations were tested at `T=13`; all eleven lift exactly over `Q` with the same denominators.

Rational interpolation on the fixed finite-minimal short Weierstrass family then reconstructs all eleven divided points as rational functions

```text
U_0(T), ..., U_10(T)
```

and every one satisfies the elliptic-curve equation identically over `Q(T)`:

```text
section 0 identity = True
section 1 identity = True
...
section 10 identity = True
VERIFIED GENERIC SECTIONS = 11 / 11
```

Section 0 required a targeted reconstruction: its x-coordinate has numerator degree 20 and denominator degree 16; a low-degree brute-force interpolation had missed it.  The degree pair was located modulo several primes and then reconstructed once over `Q` with 48 held-out exact samples.

Thus the family has an explicit hidden generic rank-11 subgroup above the original automatic subgroup.  The exact generic section relations imply an index-12288 inclusion between the displayed automatic subgroup and this hidden subgroup.  This is a subgroup statement; it does **not** by itself prove that the hidden subgroup is the full generic Mordell--Weil group, nor that the generic rank is exactly 11.

## 5. Half-lattice holes

Using the 2-saturated `T=11` union basis, all 2047 nonzero parity classes were ranked by their distance to the Mordell--Weil lattice.  A fast Babai-style CVP approximation was replayed with independent integer embeddings at 46 and 52 bits.

The ranking was completely stable:

```text
top  10: overlap  10/10
top  20: overlap  20/20
top  50: overlap  50/50
top 100: overlap 100/100
top 250: overlap 250/250
```

The first masks at `T=11` were

```text
1155, 385, 1614, 1195, 230, 1400, 259, 1064, 1715, 1646, ...
```

Exact CVP was too expensive and is unnecessary for the observed stable ranking.

## 6. Specialization-specific height geometry

A key correction was to recompute the canonical height matrix after specialization instead of transporting the `T=11` ranking to huge `T`.

At `T=27472911`, independent 46/52-bit rankings were again perfectly stable, but the overlap with the `T=11` ranking collapsed:

```text
46 vs 52 stability:
  top  10: 10/10
  top  20: 20/20
  top  50: 50/50
  top 100: 100/100
  top 250: 250/250

T=11 vs T=27472911:
  top  10: overlap  1/10
  top  20: overlap  1/20
  top  50: overlap  1/50
  top 100: overlap  6/100
  top 250: overlap 36/250
```

The specialized top mask is `230`, followed by `1137, 1598, 798, 431, ...`.

The depth distribution is non-flat:

```text
max depth    = 62.031005952919934
median depth = 32.50699568963816
max/median   = 1.9082355855079147
top10/median = 1.7129048001728786
top20/median = 1.674659032314074
```

This shows that half-lattice geometry varies substantially with specialization; any future fake-descent search should rank holes on the actual specialized height form.

## 7. Fake 2-descent quartics

For a short Weierstrass model

```text
y^2 = x^3 + A x + B
```

and `P=(x_P,y_P)`, the line-through-`-P` construction gives the quartic

```text
w^2 = m^4 - 6*x_P*m^2 - 8*y_P*m - 3*x_P^2 - 4*A.
```

An initial implementation accidentally applied this short-model formula after `global_minimal_model()`, which can produce nonzero `a1,a2,a3`; those early zero-candidate runs are invalid and must not be used as evidence.

The corrected implementation stays on the fixed finite-minimal short model.  For every tested hole, an exact known-point calibration verifies

```text
QUARTIC SELF-CHECK = True.
```

The practical problem is scale.  A representative large specialization produced a known rational slope with roughly

```text
numerator bits   = 570
denominator bits = 506
```

while the ratpoints search box `N=10^6, D=10^4` covers only about 20/14 bits.

Cremona--Stoll `hyperellred` reduced the quartic coefficient maximum from about 519 to 434 bits, but the same calibration point only moved to roughly

```text
538 numerator bits / 516 denominator bits.
```

Thus plain `hyperellratpoints` on these covers is not an effective search primitive at the tested bounds.

The corrected short-model search was also rerun with specialization-specific top-100 half-lattice holes at several high-Nagao integer specializations.  No rank-gaining candidate was found.  This negative result is sufficient to stop spending compute on naive ratpoints enlargement; it does **not** rule out the half-lattice classes themselves when coupled to a stronger 2-cover minimization/sieve implementation.

## 8. Current status and next search direction

The durable positive result is:

1. eleven explicit hidden generic rational sections;
2. an exact index-12288 enlargement of the original automatic rank-11 subgroup;
3. a 2- and 3-saturated specialization basis at `T=11`;
4. stable specialization-specific half-lattice ranking machinery.

The simple fake-2 + generic-ratpoints route is frozen after corrected top-hole searches produced no rank gain.

The next search direction is rational specialization `T=a/b`.  For every good prime `p`, a candidate is represented projectively by `(a:b)` and its local Nagao contribution depends only on `a/b mod p` (with the denominator-zero case represented by the point at infinity).  Therefore each prime admits a precomputed table of `p+1` local symbols, and billions of coprime rational parameters can be ranked using only table lookups and modular residue updates.

A second track should benchmark the complete search pipeline on a known high-generic-rank family before using it to judge new constructions.
