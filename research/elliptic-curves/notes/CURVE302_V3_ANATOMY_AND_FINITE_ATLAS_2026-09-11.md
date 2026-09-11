# Curve302 V3 anatomy: exact chart geometry, finite-atlas cost, and limits of closure

## Status and provenance

This note stores an externally supplied retrospective analysis and separates
what can be verified in this repository from what is currently attested only by
the supplied hash manifest.

The supplied analysis used excerpts of
[`half_lattice_pointed_sieve.py`](../cas/half_lattice_pointed_sieve.py), an
`x1092-elkies-pair-smoke.zip` parent-data bundle, and a library patch named
`jacobian-research-curve302-rank31.patch`. The latter two payloads and the
measurement programs/results are not present in this checkout. The preserved
[attestation](../data/curve302_v3_anatomy_verification_v1.json) declares 17
passing tests, 12,082 primary point/anchor checks, 1,734 same-class checks, and
byte-for-byte replay. Those counts and numerical measurements are therefore
**external attestation only** until every file named in its hash table is
available and matches.

The local [verification artifact](../../artifacts/generated-results/elliptic-curves/curve302_v3_anatomy_local_verification_v1.json)
independently checks the symbolic chart identities, the exact public Curve302
points, good-reduction orders, the mod-2 generic/complement split, and the
reduced parent parameter already recorded in the repository. It deliberately
reports `UNVERIFIED_MISSING_HASHED_PAYLOADS` for the external measurements.

This is a retrospective coordinate diagnostic. It does not replay the
historical V3 scheduler, reinterpret the current raw trajectory ledger, launch
a point search, or supply a rank-25 cross-fibre control panel. Its parameter
heights use denominator/shift normalization before final GL2 Gauss reduction,
so they are not historical V3 box requirements.

## The exact pointed quartic

Let

\[
E:y^2=x^3+Ax+B
\]

be nonsingular over \(\mathbf Q\), and choose an affine anchor
\(Q=(a,b)\) in the known subgroup \(S\). The chart uses

\[
t(P)=\frac{y(P)+b}{x(P)-a}.
\]

This is the slope of the line through \(-Q\). Eliminating \(x\) from the
cubic intersection and removing the known root \(x=a\) gives

\[
C_Q:\quad z^2=F_Q(t)
=t^4-6at^2-8bt-3a^2-4A,
\]

with the mutually inverse affine formulas

\[
x=\frac{t^2-a+z}{2},\qquad y=t(x-a)-b,
\qquad z=2x+a-t^2.
\]

At \(Q\) and \(O\), the slope has a pole. At \(-Q\), for \(b\ne0\), its
removable value is

\[
-\frac{3a^2+A}{2b}.
\]

These points must be handled on the smooth projective model rather than by
dividing blindly by \(x-a\). The local implementation does this explicitly.

Three maps must be distinguished:

1. \(C_Q\to E\) by the displayed \((x,y)\) formulas. This is birational and
   extends to an isomorphism of smooth projective curves.
2. \(C_Q\to\mathbf P^1\) by \(t\). This is the degree-two map associated to
   the divisor \(O+Q\).
3. \(C_Q\to E\), after the first map, by \(P\mapsto2P-Q\). This is the
   degree-four 2-covering associated to that divisor.

V3 uses the first map to recover a point on \(E\). Calling that birational map
a degree-two descent covering conflates it with the other two maps. The third
map has Kummer class \(\delta(Q)\) and rational image
\(-Q+2E(\mathbf Q)\). Its domain is a soluble genus-one curve, but the
degree-two model as a covering can still represent a nonzero class.

The deck involution of the second map is exactly

\[
(t,z)\longleftrightarrow(t,-z),\qquad P\longleftrightarrow Q-P.
\]

Since \(Q\in S\), this gives the equality of integral, unsaturated extensions

\[
S+\mathbf ZP=S+\mathbf Z(Q-P).
\]

Thus the two roots over one rational parameter supply at most one new rational
rank direction. A multiple-gain batch must contain more than this one deck
pair. Historical tied positive candidates must therefore be compared after
quotienting by equal subgroup extensions, not merely as distinct points or
displayed quotient words. The current closure-transition graph's grouping by
exact quotient line is relevant to that refinement, while its earlier
candidate-level 21/27 tie count is not by itself evidence for 21 independent
alternative channels.

There is a separate linear-algebra constraint. A one-vector extension can
increase the dimension of the intersection with a fixed core by at most one.
At a single-gain stage selected because the actual vector increases that
dimension, the actual vector already realizes the maximum possible increase.
The absence of a strictly better single-vector candidate is forced. Double-gain
batches require a distinct batch comparison.

## The covering coordinate and half-lattice height

Put \(R=2P-Q\). Exact group-law elimination gives

\[
x(R)=\frac{N_Q(t)}{F_Q(t)},
\]

where

\[
\begin{aligned}
N_Q(t)={}&at^4+4bt^3+(6a^2+4A)t^2\\
&+4abt+a^3+4B.
\end{aligned}
\]

The local checker verifies this identity by polynomial reduction using
\(b^2=a^3+Aa+B\) and \(z^2=F_Q(t)\). It also verifies the quartic map and the
deck involution symbolically.

For a fixed normalized degree-two covering model, the standard covering-height
bound specializes to

\[
B_1\le h(t)-\frac14 h_x(R)\le B_2.
\]

With \(h_x=2\widehat h+O_E(1)\), the leading term is
\(\frac12\widehat h(2P-Q)\). This is the precise half-lattice link. The
constants depend on both the elliptic curve and the chosen covering model, so
raw parameter heights cannot be compared across anchors or fibres without
controlling normalization and local distortion.

For fixed primitive integral homogeneous quartics \(N(U,V),D(U,V)\) defining
the covering \(x\)-map, let \((n,d)\) be primitive,
\(M=\max(|n|,|d|)\), and
\(g=\gcd(N(n,d),D(n,d))\). Then

\[
h_x(R)=\log\max(|N(n,d)|,|D(n,d)|)-\log g
\]

and hence

\[
h(t)-\frac14h_x(R)
=-\frac14\log\frac{\max(|N(n,d)|,|D(n,d)|)}{M^4}
+\frac14\log g.
\]

The first term measures real-place distortion and the second finite-place
cancellation. This is a concrete exact quantity for future chart comparison;
sharp local bounds for the current V3 chart normalization have not been
computed here.

## A frozen operational cost

An unrestricted minimum over coordinate systems is meaningless because a
target-dependent Mobius transformation can send the target parameter to zero.
Freeze a target-blind anchor bank \(\mathcal A(S)\) and normalization for each
chart, then define

\[
\kappa_{S,\mathcal A}(P)=
\min_{C\in\mathcal A(S)}
\log_2\max(|n_C(P)|,d_C(P)),
\]

where \(n_C(P)/d_C(P)\) is reduced. A complete record should retain the exact
integer height, parameter, root, anchor coordinates or word, and quartic
coefficient size, together with a stated convention for projective exceptional
points.

This measures finite-atlas coordinate accessibility. It does not measure the
existence of an unknown point or runtime. If the banks are nested, \(\kappa\)
is nonincreasing by construction, so a falling curve alone does not prove a
dynamic amplification mechanism.

## Supplied Curve302 diagnostic

The supplied analysis reports reconstructing the X1092 rank-17 source on the
public Curve302 model at

\[
t_0=-\frac{164518}{143797},\qquad u=20677577209=143797^2,
\]

and checking the exact scaling

\[
A_{302}=u^4A(t_0),\qquad B_{302}=u^6B(t_0),
\]

with the corresponding \((u^2x,u^3y)\) point transport. The full source bundle
needed to replay that scaling is absent. The value of \(t_0\), all 31 public
points, and good-reduction orders

\[
\#E(\mathbf F_{17})=26,\qquad \#E(\mathbf F_{31})=43
\]

are independently checked locally. The coprime orders prove trivial rational
torsion. The repository's exact mod-2 embedding of its rank-17 generic basis
also independently yields the claimed greedy public complement indices

```text
0, 2, 3, 5, 8, 12, 13, 16, 18, 21, 25, 28, 29, 30
```

in zero-based original public-point order.

The externally supplied measurement takes the first seven complement points in
that order and withholds the last seven from every anchor bank. These are not
the named historical local/strict/residual axes. At each state, its target-blind
bank consists of all signed basis points and signed pair sums/differences,
growing from 578 to 1,152 charts. It reports exact inverse/forward checks and
deck-pair verification.

The supplied values are:

| Public target, 1-based | Initial M17, 578 charts | M24, 1,152 charts | Generic-only, 1,152 charts |
|---:|---:|---:|---:|
| 17 | 26.219 | 20.371 | 26.219 |
| 19 | 60.212 | 57.106 | 52.038 |
| 22 | 49.776 | 36.275 | 41.601 |
| 26 | 55.785 | 48.487 | 46.257 |
| 29 | 54.107 | 54.107 | 54.107 |
| 30 | 49.037 | 40.000 | 36.804 |
| 31 | 63.436 | 38.280 | 50.212 |

Entries are reported \(\log_2\) exact parameter heights. The equal-count
generic control uses a fixed lexicographic sequence of signed generic triples;
it is one frozen policy, not an optimum over generic banks.

For target 31, the reported decrease from 63.436 to 38.280 bits corresponds to
a parameter-height factor of about 37.4 million. It is not a measured runtime
speedup. The equal-count generic control does better for targets 19, 26, and
30, so the table supports point-specific coordinate improvements rather than
universal superiority of subgroup enlargement.

With the original 578-chart bank held fixed, the supplied same-class ablation
reports:

| Exact target | Minimum parameter-height bits |
|---|---:|
| \(P\) | 26.219 |
| \(P+2T\) | 266.962 |
| \(P+4T\) | 1191.146 |

Here \(T\) is the first specialized generic section. These points have the
same class modulo \(2S\), yet the fixed finite bank gives radically different
costs. This demonstrates the intended distinction if the measurements replay:
finite-atlas \(\kappa\) is not a function of the Kummer class alone. Transporting
or enlarging the bank can change the result.

## Descent interpretation

The anchor class labels the covering. Replacing \(Q\) by \(Q+2T\) and \(P\)
by \(P+T\) preserves \(2P-Q\) and gives an equivalent degree-two model after
the corresponding coordinate change. Within transformations furnished by the
known subgroup, \(Q\bmod2S\) is therefore the natural chart label. Over the
full curve, soluble covering classes are organized by
\(E(\mathbf Q)/2E(\mathbf Q)\), subject to the usual model-equivalence choices.

For the genuine covering map,

\[
\delta(2P-Q)=\delta(Q).
\]

This labels the covering and does not identify the newly found point class
\(\delta(P)\). Pointed quartics therefore do not directly enumerate the unknown
Selmer classes represented by missing points.

Likewise, the ten-dimensional strict subspace observed inside the known
Curve302 quotient is realized arithmetic data. It is not automatically the
full independent Selmer space, and plotting its realized dimension against
already acquired rank partly displays rank-nullity rather than a prospective
predictor.

## Required cross-fibre control

No scientific rank-17-to-25 panel comparison accompanies the supplied data.
A suitable control should:

1. freeze exact curves, initial generic subgroups, independently known held-out
   targets, family labels, and rank lower bounds;
2. use comparable integral/minimal models and a fixed target-blind normalization
   and anchor budget;
3. measure exact inverse-chart heights, coefficient sizes, covering-class
   coverage, local distortion, and cumulative known-target coverage;
4. include equal-budget generic-only expansion and representation controls;
5. count points modulo the current subgroup and deck involution while retaining
   full integer contents; and
6. avoid interpreting a lower-bound-17 fibre with no known extras as a negative
   inventory of missing points.

Selecting every held-out point from V3 successes creates detection bias.
Independent point sources or deliberately masked known subgroups are required
for generalization claims.

## Consequences for a possible V4

No V4 search is launched by this note. Before such a search, reclassify the
historical ties by exact subgroup-extension orbits and freeze the cross-fibre
benchmark. A defensible architecture is

```text
known subgroup
  -> certified anchor classes
  -> normalized covering models with transport maps
  -> explicit coverage and height bounds
  -> bounded parameter search
  -> whole-cloud replay and exact subgroup admission
```

The first implementation priorities are deck/orbit deduplication, full integer
point words, complete coordinate transforms, and a separation between geometric
eligibility and historical execution coverage. An exact inverse-chart query can
decide whether a known point lies in a declared bounded box without rerunning a
point search. Incomplete historical logs still cannot prove which parts of that
box were executed.

Cover-height bounds can be computed before a missing point is known;
\(\kappa(P)\) itself requires \(P\). A coverage planner may improve point
finding, but it is not automatically a rank-jump predictor or an explanation
for why Curve302 has fourteen exceptional dimensions.

## Reproduction and references

Run the local exact checks with:

```bash
sage -python research/elliptic-curves/cas/verify_curve302_v3_anatomy.py
sage -python research/elliptic-curves/cas/verify_curve302_v3_anatomy.py --check
```

Primary references:

- N. Elkies, [*Three lectures on elliptic surfaces and curves of high rank*](https://arxiv.org/abs/0709.2908), especially the half-lattice and “fake 2-descent” discussion on pp. 12–13.
- J. Cremona, T. Fisher, M. Stoll, [*Minimisation and reduction of 2-, 3- and 4-coverings of elliptic curves*](https://arxiv.org/abs/0908.1741), especially Section 3.2 on rational divisors and model equivalence.
- T. Fisher, G. Sills, [*Local solubility and height bounds for coverings of elliptic curves*](https://arxiv.org/abs/1103.4944), especially equation (1.1) and Section 4.
- Z. Klagsbrun, T. Sherman, J. Weigandt, [*The Elkies Curve has Rank 28 Subject only to GRH*](https://arxiv.org/abs/1606.07178), for the independent class-group route to conditional upper bounds.

