# The 356/385 meet at the local Kummer level

<!-- status-consumer: EC-K3-R17-074D9-LOCAL-KUMMER-SEPARATION 375cb897b59e077f -->
<!-- status-consumer: EC-K3-R17-074D9-QUOTIENT-ARITHMETIC-BLOCK-OBSTRUCTION af45468d1b7d831a -->

Date: 2026-09-04  
Status: exact coordinate-complement separation and exact obstruction to the proposed quotient block decomposition

## Result

The rigid-invisible exceptional blocks at ICARM curves 356 and 385 are
already different at the local Kummer level.  There is no common fingerprint
to turn into congruences on the `074d9` family parameter, so the CRT/inward
search gate is not opened.

The stronger arithmetic-block formalization does **not** exist on the
ten-dimensional quotient as proposed.  The earlier nonpivot blocks are
deterministic RREF coordinate complements, not quotient Hilbert modules.
For a componentwise Hilbert form `B_w` to descend from the displayed
`F_2^12` to `Q=F_2^12/R`, the two-dimensional visible relation space `R`
must lie in `rad(B_w)`.  Exact replay finds failures for both records.  Thus
the graph, its placewise pairing-vector span, and a partition such as `10`,
`6+4`, or `4+3+3` are not defined independently of a complement basis.

This conclusion is deliberately narrower than saying that no abstract common
rank-jump mechanism can exist.  It rejects a mechanism characterized by the
declared local Kummer invariants: valuation parity, local image dimensions,
component-image orders, and componentwise Hilbert symbols at 2, all bad
primes, and the fixed common-good block.

The preserved coordinate-complement certificate is
[`../artifacts/generated-results/elkies-k3-r17-074d9-local-kummer-meet-v1.json`](../artifacts/generated-results/elkies-k3-r17-074d9-local-kummer-meet-v1.json).
The quotient descent certificate is
[`../artifacts/generated-results/elkies-k3-r17-074d9-quotient-arithmetic-blocks-v1.json`](../artifacts/generated-results/elkies-k3-r17-074d9-quotient-arithmetic-blocks-v1.json).

## Exact inputs and quotient blocks

The replay imports the exact `norm12-orbit-074d9` specializations and the
certified displayed quotients.  At curves 356 and 385 the quotient basis is
literally `P18,...,P29`.  Row reduction of the complete rigid-bisection image
gives the following canonical coordinate complements over `F_2`:

| curve | rigid pivots | ten-dimensional coordinate complement |
|---:|---|---|
| 356 | `P19,P25` | `P18,P20,P21,P22,P23,P24,P26,P27,P28,P29` |
| 385 | `P21,P24` | `P18,P19,P20,P22,P23,P25,P26,P27,P28,P29` |

The controls use their certified displayed quotient bases:

- curve 351: `P18,...,P25`;
- curve 376: `P18,...,P22`;
- curve 377: `P18,...,P23`;
- alternate-Q80 curve 12:
  `P2,P11,P4,P3,P6,P8,P17,P10,P28,P24,P19,P15`.

These are quotients inside the displayed public subgroups.  No full
Mordell--Weil group or rank upper bound is asserted.

## Quotient audit and obstruction

The new replay retains all twelve named generators and presents the quotient
by the two exact visible relations.  The ten RREF nonpivot labels are used
only as deterministic coordinates.  At every bad prime `p` it constructs the
following data that genuinely descend:

- the local Kummer image
  `im(kappa_p)/im(kappa_p(R))`, including its kernel in `Q`;
- the displayed component image modulo the exact visible component image,
  its exponent and annihilator filtration;
- the induced mod-two component module
  `D_p/(image(R)+2D_p)`.

It also stores every componentwise Hilbert matrix on the twelve named
generators.  Descent of such a matrix is checked by the exact condition

```text
R * B_w = 0.
```

This condition fails at the following split bad primes:

| curve | obstruction primes | number of obstructed local factors | span dimension of the obstruction rows in `(F_2^12)^*` |
|---:|---|---:|---:|
| 356 | `13,23,37,139` | 8 | 4 |
| 385 | `5,29,37,41,73,109,127` | 14 | 10 |

At each listed rational prime the cubic algebra splits as three linear
factors.  Two component forms have the same nonzero obstruction and the third
has zero obstruction, so their sum is zero.  This is consistent with the
corestricted local Tate pairing being identically trivial on local point
images, but it does not make the individual component forms descend.

If “local Hilbert pairing” is interpreted only as that genuine corestricted
`Q_p`-pairing, every matrix is zero.  Its graph is edgeless for both records,
so a coordinate presentation can be written as `1+...+1`; this is
basis-independent but vacuous, carries no coupling information, and selects
no canonical one-dimensional summands.  The nontrivial componentwise data
needed for the proposed giant-block test are exactly the data that fail to
descend.

Consequently the requested quotient pairing tensor is undefined.  An edge in
a graph formed from one selected complement can change after replacing a lift
by that lift plus a visible relation.  The certificate therefore deliberately
does not report a `10=...` connected-component partition.

There is a smaller valid diagnostic using only the quotient Kummer and
component data.  These distinguish the ten fixed canonical presentation
directions with five bad places for curve 356 (22 minimum sets; the first is
`{2,3,13,23,751}`) and four for curve 385 (the unique minimum is
`{13,29,47,89}`).  This does not repair the missing Hilbert tensor and is not
an indecomposable-module decomposition.

## Kummer representatives

For a global minimal generalized model put

```text
b2 = a1^2+4a2,  b4 = a1*a3+2a4,  b6 = a3^2+4a6
f2(z) = z^3+b2*z^2+8*b4*z+16*b6.
```

Writing `zeta` for the image of `z`, every displayed generator `P` is stored
as

```text
alpha(P) = 4*x(P)-zeta in Q[zeta]^*/Q[zeta]^{*2}.
```

The factor 4 is a global square, so this is the usual `x(P)-theta` class.
Every norm is checked exactly:

```text
N(alpha(P)) = (4*(2*y(P)+a1*x(P)+a3))^2.
```

For all six curves the cubic is irreducible.  Its discriminant is checked to
be `256*abs(Delta_E)`, and the public bad-prime list is independently checked
as a complete factorization with every factor proved prime.  The certificate
stores each Kummer element's three power-basis coefficients, but never uses
those coordinates for a cross-curve comparison.

If `r_p` is the number of field factors of the cubic algebra over `Q_p`, the
ambient local dimension used below is `r_p-1` for odd `p` and `r_p` for
`p=2`.  This is the usual identity
`dim E(Q_p)/2E(Q_p)=dim E(Q_p)[2]` away from 2, with the additional formal-group
dimension at 2.  The displayed-block dimensions are then computed directly
from their local squareclasses, not inferred from the ambient dimension.

## The two-adic separation

Both ten-direction blocks span the full two-dimensional local Kummer image,
but their remaining invariant data are different:

| curve | Kodaira type at 2 | `c_2` | ambient/block dimension | component-image order multiset on the ten directions |
|---:|---|---:|---|---|
| 356 | `I1*` | 4 | `2/2` | `1^3, 4^7` |
| 385 | `I16` | 16 | `2/2` | `1^2, 2^2, 4^1, 8^3, 16^2` |

Here the component datum is the order of the point's image in the Neron
component group.  It is invariant under relabelling or automorphisms of that
finite group.  The order is computed by exact group law and the least divisor
`m` of `c_2` for which `mP` reduces to the nonsingular component.

The componentwise Hilbert data separate the blocks again.  The cubic algebra
at 2 has local degrees `1+2` for both records.  Among the 45 unordered pairs
of invisible directions:

- curve 356 has 21 pairs with symbol `-1` in both local factors and 24 with
  symbol `+1` in both;
- curve 385 has 45 pairs with symbol `+1` in both local factors.

The product over the two factors is `+1` in every case.  That product is the
local Tate pairing under the standard cubic identification and is
universally trivial on images of local points; only the anonymous
componentwise symbols carry discriminatory information here.

## Fixed common-good block

The predeclared common-good primes are

```text
53, 67, 71, 79, 83, 97, 101, 113.
```

The exact dimensions of `E(Q_p)/2E(Q_p)` are:

| `p` | curve 356 | curve 385 | full ten-block fingerprint equal? |
|---:|---:|---:|---|
| 53 | 2 | 1 | no |
| 67 | 1 | 1 | no |
| 71 | 1 | 1 | no |
| 79 | 2 | 2 | no |
| 83 | 2 | 1 | no |
| 97 | 1 | 0 | no |
| 101 | 2 | 1 | no |
| 113 | 2 | 1 | no |

Each selected block spans the ambient local image at every prime in the
table.  The ambient dimensions agree at only three of eight primes, and even
at those three primes the anonymous unit-squareclass/direction multisets are
different.  Thus none of the eight complete block signatures matches.

As a control, the number of matching ambient dimensions for 356/385 is only
`3/8`; curve 376 agrees with curve 385 at `4/8`, while curves 12 and 351 agree
at `5/8`.  The partial dimension agreement of the record pair is therefore
not a record-specific signal.

## All bad places and conclusion

Every bad prime of all six curves is included.  At each prime the artifact
records:

- the anonymous `(e,f)` types of the local cubic factors;
- the exact ambient local Kummer dimension and the image dimensions of the
  full quotient basis and selected comparison block;
- each generator's valuation-parity support, odd-prime unit squareclass
  multiset, and automorphism-invariant component-image order;
- every self and pairwise componentwise Hilbert-symbol multiset, with the
  corestricted local Tate product checked to be `+1`.

After removing the rational-prime labels, the multisets of bad-place block
signatures for curves 356 and 385 are still unequal.  Together with the
two-adic and fixed-common-good separations, this is a negative result for the
proposed shared local fingerprint.  No congruence conditions are extracted,
no CRT class is manufactured, and no inward search is run.  The operational
next step is a direct, checkpointed family-wide parameter search with its own
independent arithmetic gate.

## Replay

```bash
sage -python \
  elkies-k3/scripts/certify_r17_074d9_local_kummer_meet.sage --check
```

The replay uses Sage 10.9/PARI, proves all discriminant factors prime, certifies
the six cubic number fields, and recomputes every local row and Hilbert symbol.
The preserved coordinate-complement v1 bytes are independently replayed by

```bash
sage -python \
  elkies-k3/scripts/certify_r17_074d9_local_kummer_meet.sage \
  --legacy-coordinate-complement --check
```
