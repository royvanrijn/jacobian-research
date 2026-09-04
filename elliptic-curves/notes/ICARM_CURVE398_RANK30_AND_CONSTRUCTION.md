# ICARM curve 398: rank 30 lower bound and construction boundary

## Status

Curve 398 was previously present only as a construction-recognition target:
the repository excluded it from the published rootless-MW17 chart and later
from the complete six-class norm-twelve atlas.  It did not have its own point,
independence, torsion, conductor, or construction certificate.

It now has an independent exact certificate for

```text
rank E(Q) >= 30.
```

The public construction mechanism is also identified: Elkies and Klagsbrun
found the curve in September 2025 by searching elliptic fibrations on their
Picard-rank-19 K3 surface `X948` which have one reducible fibre of type `I2`
or `III`.  The associated root lattice is `A1`, so Shioda--Tate gives generic
Mordell--Weil rank

```text
19 - 2 - rank(A1) = 16.
```

Their announcement says that hundreds of these MW16 fibrations were searched
with an improved codebase using Drew Sutherland's `smalljac` library.  The
rank-30 specialization is a jump by at least fourteen over the generic rank.

The exact fibration, rational base parameter, sixteen specialized section
vectors, and search transcript for curve 398 are not public in the curve page
or announcement.  Thus the method is known, but the precise construction is
still `UNKNOWN` and is not yet independently reproducible.

The authors report exact rank 30 conditional on GRH for number fields.  This
repository has not replayed that upper bound, so the canonical theorem remains
the unconditional lower bound only.

Public sources:

- [ICARM curve 398](https://elliptic-rank.icarm.cloud/curve/398) and its
  [machine-readable record](https://elliptic-rank.icarm.cloud/curve/398.json);
- [Elkies--Klagsbrun rank 27--30 announcement](https://listserv.nodak.edu/cgi-bin/wa.exe?A2=NMBRTHRY;d593eaaa.2608&S=);
- [Elkies' explicit paper on the published MW17 fibration](https://arxiv.org/abs/2608.25406),
  which supplies the current public equation for `X948` but not curve 398's
  different MW16 fibration.

## Exact curve and rank certificate

The public integral generalized Weierstrass model is

```text
y^2 + x*y = x^3 + A*x + B,

A = -12892599774455576272301592959047823530919513428112484011550,
B = 560755046348395412977088824999503890617558856687636981223935662848386484980312656296132.
```

All thirty public point pairs satisfy this equation exactly.  Under

```text
X = 36*x + 3,
Y = 108*(2*y + x),
```

they lie on the integral short model

```text
Y^2 = X^3 - 27*c4*X - 54*c6,
```

where

```text
c4 = 618844789173867661070476462034295529484136644549399232554401,
c6 = -484492360045013636812204744800499628677331653669724066470531855995232127989814233888689649.
```

The short 2-division cubic has no root modulo 23.  Exact exhaustive reductions
in `E(F_p)/2E(F_p)` at

```text
11, 13, 17, 43, 47, 53, 61, 67, 83, 103, 107, 109, 113,
131, 139, 149, 157, 173, 179, 193, 199, 211, 227, 229, 239
```

produce 31 binary rows with column rank 30.  Infinite descent therefore proves
that the thirty displayed points are independent in `E(Q)`.  Independently,

```text
#E(F_11) = 18,  #E(F_23) = 31,  gcd(18,31) = 1,
```

so rational torsion is trivial.

PARI's exact `ellisomat` enumeration returns one rational isomorphism class
and the degree matrix `Mat(1)`.  Thus curve 398 has no nontrivial rational
isogeny to another elliptic curve over `Q`.  Combined with a direct
`j`-recognition miss for any candidate family, this also rules out reaching
curve 398 through a rational isogeny from that family.

The numerical height matrix and regulator are retained only as diagnostics;
they are not used in the independence proof.

## Local arithmetic

PARI verifies that the displayed integral model is globally minimal and
semistable.  Its minimal discriminant factors as

```text
2^18 * 3^10 * 5^7 * 7^4 * 19^3 * 29^4 *
31^2 * 37^2 * 41^2 * 59^2 * 73^2 * 79^3 *
101^2 * 151^2 * 197^2 *
20084070565614383 *
277131689105980733414153 *
91586120381369539248864736998169886452607476058975046580346704317125237.
```

Every bad prime is multiplicative.  The Kodaira fibres, in the same order, are

```text
I18, I10, I7, I4, I3, I4,
I2, I2, I2, I2, I2, I3,
I2, I2, I2, I1, I1, I1.
```

Consequently every conductor exponent is one, and their product is the public
conductor

```text
2835647668537242470520670933169702570972603599817629325681756626410520535791722336193301955569324372206950385073317702302276416849610.
```

The Tamagawa product is `46,448,640` and the global root number is `+1`,
consistent with the displayed even-rank subgroup.  These local calculations
are exact except for the separately labelled numerical height diagnostic.

At the public snapshot retrieved on 2026-09-04, ICARM marks curve 398 as the
record among curves of rank at least 30 simultaneously for conductor, naive
height, Faltings height, and absolute discriminant.  This is dated database
metadata, not a timeless optimality theorem.

## Point-cloud construction fingerprint

Twenty-eight of the thirty displayed point pairs are integral.  At the fifteen
small multiplicative primes, the counts of available public points reducing
to the singular node are

```text
p:       2  3  5  7 19 29 31 37 41 59 73 79 101 151 197
node:   29 28 29 24 25 28 22 24 21 23 24 26  15  24  21
```

Point 3 is unavailable modulo 2 and point 15 modulo 3 because of their
denominators.  None of the thirty points reduces to the node at any of the
three large `I1` primes.  This exact clustering is a useful target fingerprint
for recovering the hidden specialization; it is not by itself a family or
rank theorem.

## Exact construction exclusions

The forensics now give four sharp negative boundaries.

1. The exact degree-24 `j`-recognition polynomial for the published R17
   fibration is irreducible modulo 1009.  Curve 398 is not a rational fibre of
   that chart, even after quadratic twisting.
2. The complete 43-chart norm-twelve atlas has six rational `j`-map classes,
   and curve 398 misses all six.  This rules out the currently certified
   norm-twelve shared-zero models, not other fibrations on `X948`.
3. The exact rational-isogeny class of curve 398 is a singleton.  This
   supersedes the earlier degree-`3,5,7,11` modular-polynomial checks: no
   nontrivial rational isogeny from any rational fibre can lead to curve 398.
4. The repository has one equation-explicit `A1`/MW16 fibration from the
   fixed-corridor reverse lift.  For this family, the curve-398 `j`-preimage
   polynomial has degree 24 and no projective root modulo 179.  Hence curve
   398 is not a rational fibre of this particular A1 family, twist-stably.

The fourth test is deliberately narrow.  It does not enumerate or exclude the
hundreds of other A1 fibrations in the authors' search.

## Equation-source inventory

A repository and public-source audit finds no second rational A1 equation that
can presently be tested against curve 398.

- The new public `X948` paper gives the published rootless MW17 fibration and
  its sections.  Curve 398 already misses that `j`-map.
- The fixed-corridor reverse lift supplies the only characteristic-zero
  `QQ` A1/MW16 equation currently available here.  It is the family excluded
  modulo 179 above.
- The alternate Q80 route has an exact generic A1 lattice frame and final
  neighbour transport, but its handoff explicitly records
  `generic_characteristic_zero_equation = null` and
  `equation_status = NOT_YET_COMPILED`.  Its characteristic-zero equation
  corridor is a CM24 specialization over `QQ(sqrt(-3))`, not the generic
  determinant-948 family required for rational `j`-recognition.
- Three retained, bounded Q80 shell windows contain 48 generic A1/MW16
  candidates.  They collapse at CM24 to two nef specialization classes of
  old-fibre degrees 47 and 43.  Those computations are bounded windows and do
  not provide generic `QQ` Weierstrass equations.
- The announcement's hundreds of searched A1 fibrations are not supplied as
  equations, lattice markings, or a parameter ledger in any public source
  found in this audit.

This is a complete inventory of the currently available equation sources, not
an exhaustive theorem about all A1 fibrations on `X948`.  It explains exactly
why further rational `j`-recognition cannot proceed without first compiling a
new generic A1 equation or obtaining the authors' census.

## Remaining exact reconstruction target

The next construction proof should be certificate-driven:

1. obtain or reconstruct the authors' A1-fibration census on `X948`;
2. solve exact projective `j`-preimage equations for curve 398 until a rational
   parameter is found;
3. transport and specialize all sixteen generic sections;
4. express them integrally in the public rank-30 subgroup and certify the
   fourteen-dimensional exceptional quotient;
5. record the specialization and the sieve inputs that selected it.

Until those data exist, the exact fibration/parameter/section map stays
`UNKNOWN`; bounded lattice searches and the current finite exclusions do not
close it.

## Reproduction

From the repository root, with PARI/GP available:

```bash
.venv/bin/python elliptic-curves/cas/verify_icarm_curve398_rank30.py --check

.venv/bin/python -m unittest \
  elliptic-curves/tests/test_icarm_curve398_rank30.py -v
```

The checker recomputes point membership, the finite-reduction independence
matrix, torsion witnesses, discriminant and conductor factorizations, every
local reduction and root number, the complete rational isogeny class, the
node-incidence fingerprint, and the modulo-179 A1 exclusion.  The pinned output is
`artifacts/generated-results/elliptic-curves/icarm_curve398_rank30_and_construction_v1.json`,
with whole-file SHA-256
`1fd4f23ff2167321be0e3a7bf12b693f0a9ebe26d1e2125ce131da30ad05bf60`.
