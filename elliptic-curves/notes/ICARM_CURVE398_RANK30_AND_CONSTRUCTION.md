# ICARM curve 398: hidden A1/MW16 recovery and blind rank-30 rediscovery

## Status

Curve 398 has an independent exact certificate for

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

The construction boundary is now closed independently.  A complete modular
screen of the 63,917 minimum-norm-eight neighbours of the equation-explicit
`norm12-orbit-11952` rootless chart leaves two candidates; exact factorization
shows that both have rational curve-398 parameters and specialize to curves
`Q`-isomorphic to curve 398.  Both candidates are now compiled as
`I2+22I1`/MW16 fibrations and all sixteen sections from each saturated generic
basis are transported into the public rank-30 subgroup.  Unexpectedly, the two
specialized integral MW16 subgroups are equal, not transverse.  An independently
redacted half-lattice search from the first basis recovers the full displayed
rank-30 subgroup.

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

## Recovered hidden A1/MW16 fibration

On the rootless alternate-Q80 chart `norm12-orbit-11952`, priority trace
`16875` has equation-section coordinates

```text
w = (0,0,0,0,0,-1,0,-1,-1,1,1,1,0,-1,0,0,1),
w.M.w = 8.
```

Its isotropic class is

```text
D = (2,2,w) = O + P_w,
D.F = 2,  D.O = 0.
```

The residual-chord quartic and its binary-quartic invariants produce a short
Weierstrass equation with coefficient degrees `(8,12)`.  Its finite
discriminant is squarefree of degree 22, while the orders at infinity are
`(ord(c4),ord(c6),ord(Delta))=(0,0,2)`.  Thus the generic fibre configuration
is exactly

```text
I2 at infinity + 22 I1,
```

and Shioda--Tate gives generic Mordell--Weil rank 16.

The exact curve-398 parameter in this pencil is

```text
lambda = -273478312517509127154149830485048828022673347107308547939067553994727903425458545978043182638015899676311550557441827100822466901248 / 243076210150914055804756105904064536659703543720469425499709810733677965174759784940636972086422417178984090368085211
```

The cross-multiplied `j` equation factors over `Q` into degrees `1+23`, and
the specialization at the displayed linear root is `Q`-isomorphic—not merely
quadratically twist-equivalent—to curve 398.

The fixed parity coset contains exactly 166 old sections meeting `D` once.
Taking

```text
(0,0,-1,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
```

as the new zero, the compiler selects sixteen of the other degree-one curves.
Their exact Shioda height Gram is half-integral of rank 16 and determinant
`474=948/2`; this proves saturation of the MW16 basis.  After specialization,
height-dual recovery followed by exact rational group-law replay expresses
all sixteen in the ordered public 30-point group.  The 30-by-16 coordinate
columns have Smith factors all one and maximum absolute entry 65.

The complete discovery screen used 34 primes.  It excludes 63,915 of the
63,917 committed norm-eight classes and leaves priorities `16875` and `63669`.
Exact factorization gives one rational parameter for each survivor, and both
specializations are `Q`-isomorphic to curve 398.  The compiled fibration above
uses the equation-cheaper first survivor; the second is compiled below.

## Exact two-parent collision

The second survivor has priority `63669`, orbit `0x06119`, and trace

```text
w = (-2,1,3,0,1,0,-4,0,-1,-2,3,0,-2,-1,2,-3,2).
```

Its residual-chord invariants again give degrees `(8,12)`, squarefree finite
discriminant of degree 22, and fibre configuration `I2+22I1`.  Its unique
rational curve-398 parameter is

```text
lambda = -541266381922712529166100960678122326542295329017811351186978386511278040283284966392829974955759690589708833207806994323443840 / 1966455527134683136777607542029510829585376789066249361045523577208160221833556912096256713936098199933472678271.
```

The second fixed parity coset contains exactly 180 degree-one old sections.
With

```text
(0,0,1,-1,0,-1,-1,0,1,-1,1,1,1,0,0,-1,0)
```

as zero, the first sixteen independent enumerated sections have Shioda Gram
of rank 16 and determinant `474=948/2`, so they form a saturated generic MW16
basis.  Exact specialization, height-dual coordinate recovery, and rational
group-law replay give a second `16 x 30` integral embedding in the ordered
public basis.

Let `C1` and `C2` be the two embedding matrices and let `G1,G2` be their row
groups in `M30`.  Exact integer solves produce mutually inverse matrices
`U,V in GL(16,Z)`, both of determinant one, with

```text
C2 = U*C1,   C1 = V*C2,   U*V = V*U = I16.
```

Consequently this is equality of integral subgroups, not merely equality after
tensoring with `Q`:

```text
G1 = G2,
rank(G1 intersection G2) = 16,
rank(G1 + G2) = 16.
```

The Smith diagonal of the stacked `32 x 30` generator matrix is sixteen ones
followed by fourteen zero directions.  Thus

```text
M30 / (G1 + G2) = Z^14,
[M30 : G1 + G2] = infinity.
```

The ambiguity of sign in the `Q`-isomorphism of either fibre does not change
its subgroup.  This pair therefore does not explain any of the apparent
fourteen-rank specialization jump by transverse generic directions.  It is a
strong negative control for the proposed rank-32 mechanism: a useful search
must test the relative specialized subgroups, since distinct MW16 fibrations
with `Q`-isomorphic fibres can collide completely.

## Blind MW16 half-lattice calibration

The search input is a separately pinned redacted fixture containing only the
short equation, the sixteen specialized generic points, and their exact
generic MW16 height Gram.  It contains neither the public 30-point fixture,
the generic-to-public embedding, nor any held-out point coordinate.

The generic half-lattice census is complete over all `2^16=65,536` parity
classes.  Exactly twelve have maximum twice-norm 23.  The frozen ledger's
legacy label `half-lattice depth 23/8` is only a fixed-basis chart-priority
score, not an arithmetic depth or covering invariant.  Searching those twelve
first-priority birational charts through reduced-coordinate height `10^5`
raises the exact discovered subgroup

```text
M16 -> M21.
```

The adaptive quotient rule then searches all

```text
12 * (2^5 - 1) = 372
```

nonzero five-bit lifts and raises it to

```text
M21 -> M30.
```

All 384 charts completed the declared bounded search with zero timeouts.  The
five initial and nine adaptive gains are certified independent by exact finite
reduction and every relation is replayed by rational group law.

The chart order has no Selmer interpretation: these pointed quartics are
birational models of the same elliptic curve, not nontrivial 2-coverings.
Likewise, the five-bit quotient Hamming weight only orders charts in the
recorded quotient basis.  Any lattice enlargement or basis/complement change
invalidates the ordering and its calibration; chart identities,
representatives, scores, and order must be recomputed and state-fingerprinted
before reuse.  A miss would imply neither point absence nor Selmer structure.

Only after the blind run stopped did the verifier load the public fixture.  A
Smith completion extends the primitive generic MW16 coordinate rows to a
unimodular basis of the public rank-30 lattice; its last fourteen rows define
the held-out complement.  All fourteen replay exactly in the blind basis, and
mutual integral embeddings prove

```text
blind discovered subgroup = displayed public rank-30 subgroup.
```

The run stops at its predeclared next-wave limit: fourteen quotient bits would
require 196,596 lifts.  Thus it proves full recovery of the displayed subgroup,
not algorithmic stability, rank exactly 30, or saturation in the unknown full
group `E(Q)`.  This is nevertheless the requested cross-fibration calibration:
the generic MW16 half-lattice plus the adaptive quotient rule rediscovers a
rank jump of fourteen outside the R17/MW17 family.

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

## Construction boundary closed

The earlier equation-source audit remains useful history: the published R17
chart, all six shared-zero norm-twelve classes, and the fixed-corridor A1 model
really do miss curve 398.  The missing step was to recognize that every
minimum-norm-eight class on a rootless chart gives a zero-neutral divisor
`D=(2,2,w)=O+P_w`, hence precisely the first complete A1/MW16 layer.  The new
screen and compiler close that boundary without importing the authors'
unpublished fibration census.

The result promotes curve 398 from a provenance target to the major
cross-fibration calibration control.  A rank-32 campaign in this recovered
family would still need sixteen independent quotient directions and all the
same exact finite-reduction and residual-descent gates; the present bounded
success is not evidence for a rank-32 theorem by itself.

## Reproduction

From the repository root, with PARI/GP available:

```bash
.venv/bin/python elliptic-curves/cas/verify_icarm_curve398_rank30.py --check

.venv/bin/python -m unittest \
  elliptic-curves/tests/test_icarm_curve398_rank30.py \
  elliptic-curves/tests/test_curve398_mw16_adaptive_half_lattice.py \
  elliptic-curves/tests/test_curve398_two_parent_collision.py -v

sage -python elkies-k3/scripts/compile_icarm_curve398_hidden_a1_mw16.sage --check

sage -python elkies-k3/scripts/compile_icarm_curve398_second_a1_mw16_collision.sage --check

sage -python elliptic-curves/cas/verify_icarm_curve398_two_parent_collision.sage

sage -python elliptic-curves/cas/verify_curve398_mw16_adaptive_half_lattice_search.sage --check
```

The checker recomputes point membership, the finite-reduction independence
matrix, torsion witnesses, discriminant and conductor factorizations, every
local reduction and root number, the complete rational isogeny class, the
node-incidence fingerprint, and the modulo-179 A1 exclusion.  The pinned output is
`artifacts/generated-results/elliptic-curves/icarm_curve398_rank30_and_construction_v1.json`,
with whole-file SHA-256
`1fd4f23ff2167321be0e3a7bf12b693f0a9ebe26d1e2125ce131da30ad05bf60`.

The recovered construction, blind search, and post-search verification are
respectively pinned in

- `artifacts/generated-results/elliptic-curves/icarm_curve398_hidden_a1_mw16_v1.json`;
- `artifacts/generated-results/elliptic-curves/icarm_curve398_two_parent_collision_v1.json`;
- `artifacts/generated-results/elliptic-curves/curve398_mw16_adaptive_half_lattice_blind_v1.json`;
- `artifacts/generated-results/elliptic-curves/curve398_mw16_adaptive_half_lattice_verification_v1.json`.
