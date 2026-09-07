# Curve 302: local filtration of the recovered quotient

This is a post-search, exact calculation on the sealed 302 calibration.  It
uses the public 31-point group only after all generic-only search waves have
ended.  It constructs no point, selects no chart and gives no rank upper
bound.  Mathematical status remains with
[`MATH_STATUS.json`](../../MATH_STATUS.json), in particular the recovery
calibration entry.

Let `D` be the displayed rank-31 public group, `M17` the literal
specialization of the recovered parent, and `M24` the final recovered group.
The public-word identities are integral, and both `M17` and `M24` are
primitive in `D`.  Thus

\[
 \operatorname{rank}(D/M_{17})=14,\qquad
 \operatorname{rank}(M_{24}/M_{17})=7,\qquad
 \operatorname{rank}(D/M_{24})=7.
\]

The relevant distinction is not an arbitrary integral `7+7` splitting.
Work modulo two, apply the cubic Kummer map, and let `loc_S` collect the
places at infinity and all twenty bad rational primes.  The exact output
gives

\[
\begin{array}{c|ccc}
 & M_{17}/2 & M_{24}/2 & D/2\\
\hline
\dim\operatorname{loc}_S & 17 & 21 & 21\\
\dim\ker(\operatorname{loc}_S) & 0 & 3 & 10.
\end{array}
\]

The product of individual local point-image dimensions is 22; the joint
image on `D` has dimension 21.  The generic subgroup fills each individual
local point image, but its joint image has dimension only 17.  Hence the
observed exceptional quotient has the canonical mod-two filtration

\[
0 \longrightarrow \mathbf F_2^{10}
 \longrightarrow D/M_{17}\pmod 2
 \longrightarrow \mathbf F_2^4
 \longrightarrow0.
\]

`M24/M17` maps onto all four joint local patterns and has a
three-dimensional strict-at-`S` kernel.  Since `loc_S(M24)=loc_S(D)`,

\[
(D/M_{24})\pmod2
 \;\cong\;
 \ker(\operatorname{loc}_S|D)
 /\ker(\operatorname{loc}_S|M_{24})
 \;\cong\;\mathbf F_2^7.
\]

So the recovered seven are **four locally visible directions plus three
strict directions**.  Every remaining mod-two quotient class has a strict
representative.  This is the clean separation sought by the controlled
experiment: the adaptive chart policy has already covered all joint bad-place
patterns, whereas its seven-direction residual is entirely in the strict
local sector.

One short strict witness is the word

\[
P_4+P_{12}+P_{18}+P_{23}+P_{27}.
\]

Its retained Kummer class is trivial at `S`; at the good prime 47 the cubic
has root 29 and the product has residue 23, a nonsquare.  This records a
nonzero global Kummer class, rather than inferring one from a zero local row.

## What this prioritizes

The seven residual directions are exactly where the strict-class/half-ideal
constructor should be tested.  The immediate clean endpoint is to construct
one strict class from the equation and generic MW17 data alone, then use the
sealed public data only for retrospective comparison.  This report does not
yet compute half ideals, Artin characters, 4- or 8-divisibility, canonical
coset minima, or a class-group upper bound, so it does not identify a
particular ideal-class block with any individual residual direction.

In particular, it does **not** prove that the strict classes are new
ideal-class blocks, that they have rational representatives without the known
points, or that the full Mordell--Weil group has rank 31.  Those are separate
questions.

## Exact replay

The standalone checker recomputes every local squareclass signature from the
literal public equation, verifies the 24 integral group identities, checks
both Smith primitivity statements, and independently certifies the full
31-dimensional mod-two public Kummer image using complete-splitting good
primes through 751:

```sh
sage -python elliptic-curves/cas/verify_curve302_recovered_quotient_local_filtration.sage --check
```

Its immutable output is
[`curve302_recovered_quotient_local_filtration_v1.json`](../../artifacts/generated-results/elliptic-curves/curve302_recovered_quotient_local_filtration_v1.json).
The exact recovery identities it consumes are in
[`curve302_recovered_public_span_v1`](../../artifacts/generated-results/elliptic-curves/curve302_recovered_public_span_v1/manifest.json),
and the generic-only local arithmetic is bound by
[`rank_jump_curve302_strict_constructor_arithmetic_v1.json`](../../artifacts/generated-results/elliptic-curves/rank_jump_curve302_strict_constructor_arithmetic_v1.json).
