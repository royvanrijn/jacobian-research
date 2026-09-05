# The fixed-field experiment at u=-1: sixteen obstructed dimensions

Status: **the Cassels--Tate pairing on the certified eighteen-dimensional
subspace has rank 16 and radical dimension 2**. There are exactly three
nonzero compatible combinations left in this subspace. Their rational-point
realization remains **UNKNOWN**.

The curve, labelled cubic-field identification, and ordered basis
\(w_1,\ldots,w_{18}\) are those of the `parameter_u = "-1"` row in
[`fixed_cubic_field_fermigier_rank20_local_kummer_u2_v1.json`](../../artifacts/generated-results/elliptic-curves/fixed_cubic_field_fermigier_rank20_local_kummer_u2_v1.json).
The [family note](FIXED_CUBIC_FIELD_VARYING_CURVE_EXPERIMENT_2026-09-04.md)
proves their global independence and every local condition. This computation
uses that subspace; it does not require a class group or a complete Selmer
basis.

## Exact result and interpretation

Write \(W=\langle w_1,\ldots,w_{18}\rangle\), and identify the pairing
values \(0,1/2\) with \(0,1\in\mathbf F_2\). All 153 independent entries
were evaluated. The matrix has 68 nonzero upper-triangular entries and

\[
\operatorname{rank}(M)=16,\qquad
R=\ker M,\qquad \dim R=2.
\]

An exact basis of the radical is

\[
\begin{aligned}
r_1={}&w_1+w_3+w_4+w_6+w_7+w_8+w_{10}+w_{11}+w_{15}+w_{18},\\
r_2={}&w_2+w_3+w_4+w_7+w_{12}+w_{13}+w_{14}+w_{16}+w_{17}.
\end{aligned}
\]

The summary certificate contains the full matrix and an invertible change
of basis taking it to eight hyperbolic blocks plus a two-dimensional zero
block. Its eight displayed pairs span a 16-dimensional complement \(H\)
with \(W=H\oplus R\).

For every \(w\in W\setminus R\), some \(v\in W\) has
\(\langle w,v\rangle_{CT}=1/2\). Rational Kummer classes pair trivially
with all Selmer classes, so every such \(w\) maps to a nonzero element of
\(\Sha(E_{-1})[2]\), and its covering has no rational point. Thus:

- **262,140 classes are provably obstructed:** all \(2^{18}-2^2\) elements
  outside \(R\), including all eighteen original basis elements.
- **Three nonzero classes remain compatible:** \(r_1,r_2,r_1+r_2\).
- \(\dim_{\mathbf F_2}\Sha(E_{-1})[2]\ge16\), and
  \(\dim(W\cap\delta(E_{-1}(\mathbf Q)))\le2\).

This does **not** give a rank upper bound for the whole curve. Also, \(R\)
is the radical against \(W\), not necessarily the intersection of \(W\)
with the radical against the full Selmer group. Its elements may pair
nontrivially with classes outside \(W\), or survive this pairing while
remaining nontrivial higher-divisible Sha classes.

The separately certified point \(Q=(A+1,A-B+1)\) has Kummer class
\(\eta\notin W\), as proved in the family note. It remains a rank-one
lower bound. On \(W\oplus\langle\eta\rangle\), the pairing therefore has
rank 16 and radical \(R\oplus\langle\eta\rangle\) of dimension 3. This
does not assert that the nineteen-dimensional space is the full Selmer
group.

## Arithmetic certificate

The computation implements [Fisher's Theorem 3.1 and Remark 3.3](https://antsmath.org/ANTSXV/papers/ANTS-XV_fisher.pdf).
For each basis class and each pairwise sum, eliminate \(d\) from the
original quadrics to obtain a ternary conic, parametrize it over
\(\mathbf Q\), and minimize and reduce the resulting binary quartic.
Every rational transformation is retained. The verifier substitutes the
resulting quadratic polynomials \(\gamma(x)\) into the **original**
cubic-field identity and checks

\[
[\theta](\beta\gamma^2)+k^2g(x)=0,
\qquad [\theta^2](\beta\gamma^2)-k^2g(x)=0.
\]

The two parameter matrices are invertible. Hence \(y^2=g(x)\), with
\(d=ky\), is the declared class's covering, not merely a curve with the
same Jacobian invariant.

The quartics are rescaled to common exact invariants \(I,J\). In
\(L=\mathbf Q[\phi]/(\phi^3-3I\phi+J)\), the certificate records a
square root \(m^2=z(g_i)z(g_j)z(g_{i+j})\). Fisher's formula then gives
a quadratic \(\gamma_{ij}\) and the pairing as a product of local
Hilbert symbols. The verifier derives this quadratic again from \(m\)
and the three quartics, rather than accepting a supplied pairing form.

For every relevant place the evidence includes a rational \(x_v\), the
exact values \(g_i(x_v)\) and \(\gamma_{ij}(x_v)\), and the Hilbert
symbol. It verifies that the first value is a nonzero local square and the
second is nonzero. Elementary rational Hilbert formulas independently
replay the producer's PARI values.

The quadratic is normalized to primitive integral coefficients. Every
prime dividing the first quartic's discriminant or the second quartic's
leading coefficient is covered, along with 2, 3, 5, 7 and infinity.
Exact factor reconstruction and proved primality check this support.
Fisher's Remark 3.3 proves that all omitted places contribute zero.
Changing the quadratic by its nonzero rational content or sign does not
change the global product, by Hilbert reciprocity.

The final replay checks **178 cover maps, 153 matrix entries, eight
additional symmetry/bilinearity entries, and 65 distinct primes**. Four
additional entries directly pair the two radical generators with \(w_1\)
and \(w_6\), obtaining zero. The published nonzero example on \(571a1\)
also replays, including its nontrivial real contribution.

The initial attempts through Magma's public calculator did not complete
within its 60-second limit. They supply no pairing values. The successful
computation uses SageMath 10.9 and PARI 2.17.3 locally, with a 60-second
entry limit, a 900-second main campaign limit, per-cover/per-entry
checkpoints, and bounded additional controls. Uncomputed or failed entries
never become zero entries; the complete matrix is required before a
radical is emitted. Initial producer exceptions were repaired and the
affected entries recomputed before certification.

## Point-solving queue

The only nonzero inherited candidates are now these three combinations.
The masks below refer to the **twenty-dimensional anchor basis**, not to
the eighteen-dimensional surviving basis.

| Class | Anchor mask | Surviving-basis weight | Largest absolute reduced quartic coefficient |
|---|---:|---:|---:|
| \(r_1+r_2\) | 1047173 | 13 | 90955156224555672488895745651920 |
| \(r_1\) | 596921 | 10 | 309800464478559462126352851414620 |
| \(r_2\) | 450876 | 9 | 367648610965534709612552532137443 |

Use \(g_1=r_1+r_2\) and \(g_2=r_1\) as the next explicit-cover
generators; the third candidate is \(g_1+g_2\). This order uses the actual
recorded reduced coefficient heights as a search scheduling heuristic,
not a theorem about which class is soluble. It also shows why searches
limited to small surviving-basis weight miss every remaining candidate.

All three reduced quartics were searched with PARI `hyperellratpoints`
at height argument 10,000; each returned no affine point. Rational points
at parameter infinity were checked separately and none occur. The cheap
certificate replay reruns these three bounded searches. Every result
remains **UNKNOWN** as a point-or-Sha classification inside the radical.

Further point solving should use these three classes, with changes of
presentation or translation by the already known point \(Q\) as useful.
Earlier searches on arbitrary basis elements or low-weight sums are
historical regressions; their bounded misses are no longer the reason
those nonradical classes are excluded.

## Six globally minimal radical models and a bounded higher-descent attack

The follow-up [model certificate](../../artifacts/generated-results/elliptic-curves/fixed_field_radical_models_v1.json)
works exclusively on the three radical masks and their translations by
\(Q\). It **does not decide any of the three classes**. All remain
**UNKNOWN**, and the certified curve rank lower bound remains one.

For each mask, the original cubic-field quadrics are parametrized anew,
first for \(\beta\) and then for \(\beta\eta\). The six resulting
quartics, conic transformations, and their exact maps to the specified
Kummer classes are retained in the
[compressed evidence](../../artifacts/generated-results/elliptic-curves/fixed_field_radical_models_evidence_v1.json.gz).
For a quartic

\[
y^2=a s^4+b s^3t+c s^2t^2+d st^3+e t^4,
\]

the alternative projective model is

\[
XZ-Y^2=0,\qquad
W^2=aX^2+bXY+cY^2+dYZ+eZ^2.
\]

Magma's degree-four `Minimise` and `Reduce` supply new equations and both
change-of-equations and change-of-variables matrices. These are replayed
over \(\mathbf Q\), independently of Magma. In the convention where the
two matrices \(H_1,H_2\) are Hessians, let \(I,J\) be the binary-quartic
invariants of \(\det(sH_1+tH_2)\). Then
\(c_4=I\), \(c_6=J/2\), and
\(\Delta=(c_4^3-c_6^2)/1728\). All six models have **exactly the same
\(c_4,c_6,\Delta\) as the global minimal Jacobian**. Since they are
integral and locally soluble, their level is zero: this certifies global
genus-one minimality, not merely hyperelliptic-model minimization.
It does not assert an absolute minimum for coefficient height over all
presentations. See the
[genus-one minimization reference](https://docs.magma-maths.org/ArithmeticGeometry/ModelsOfGenusOneCurves/g1minred.html).

Each model was passed to `PointsQI` with height argument \(10^7\), seed
one, and `OnlyOne := true`. All six searches completed with no returned
point. This is Elkies's **p-adic lattice method**, rather than another
increase of `hyperellratpoints` height; the height argument refers to the
new projective coordinates. It is not a proof of global insolubility.
The [PointsQI reference](https://magma.maths.usyd.edu.au/magma/handbook/text/1570)
describes its search semantics. The preparation runner retains a height-one
quartic smoke check and explicitly checks the parameter-infinity chart.
The subsequent [search-geometry audit](#search-geometry-audit-and-revised-method)
proves that every nominal height box used here was already empty by
elementary real inequalities. These runs should not be interpreted as
substantial new point-search coverage.

For every possible hit the verifier applies the inverse projective map,
recovers the quartic point in either projective chart, and replays the
conic/quartic maps and the original cubic-field square identity. A hit on
a translated cover is mapped back by subtracting \(Q\), with the chord
identity supplying an explicit square root for the original class.
Realized masks would immediately give an exact binary independence
certificate; valuation parity above 19 adds the independent \(Q\).
There were no hits to certify in this run.

Four-descent was attempted on all six quartics with the known
discriminant primes supplied through `StoreFactor`, without a GRH setting.
Every attempt was incomplete. The evidence retains exact submitted jobs,
raw responses, and failure classifications; a trailing `DONE` after a CAS
error is explicitly rejected. In particular, preloading the factors
removed an initial 153-digit integer-factorization bottleneck and exposed
the harder quartic number-field Selmer computation. The original
\(r_1+r_2\) presentation reached the field

\[
z^4-z^3-143526494875963836950z^2
-796653480445236814916542885700z
-903069777838011393730646221980998372168,
\]

whose discriminant is
`17207612547621358265560224336784329653572551167050221201938192360`.
The public calculator exhausted memory during its class-group calculation.
A separate local PARI pilot on this field, with a 1.024 GB PARI stack and
a 300-second wall limit, also did not return a class group. Its exact input
and output are retained as `local_field_pilot` in the compressed evidence.
`bnfinit` alone would require further certification even if it returned.
None of these resource failures is an empty Selmer set or a Sha obstruction.

Offline replay and the narrow map/failure regressions are:

```sh
sage -python elliptic-curves/cas/run_fixed_field_radical_covers.py
sage -python -m unittest elliptic-curves/tests/test_fixed_field_radical_covers.py \
  elliptic-curves/tests/test_fixed_field_point_realization.py
```

Prepare the six models and standalone Magma jobs without network access:

```sh
sage -python elliptic-curves/cas/run_fixed_field_radical_covers.py --prepare
```

The `.m` files can be run on a local Magma installation. A bounded public
calculator regeneration is explicitly opt-in:

```sh
sage -python elliptic-curves/cas/run_fixed_field_radical_covers.py \
  --workdir artifacts/local/fixed-field-radical-fresh --prepare \
  --online-qi --online-descent --collect
```

Each job is checkpointed; the calculator imposes 60 seconds per job, and
the HTTP client has an 80-second timeout. Existing responses, including
failures, are retained. `--collect` refuses to publish a completed
four-descent without first examining its arithmetic. The remaining
mathematical task is still a rational point or a further obstruction on
one of the three classes. These computations supply neither additional
rank-transport evidence nor a refutation of transport for the radical.

## Search-geometry audit and revised method

The [exact geometry certificate](../../artifacts/generated-results/elliptic-curves/fixed_field_radical_search_geometry_v1.json)
proves that all six recorded models have no rational point of naive
projective height at most \(B=10^7\). This is a **bounded height exclusion**,
not a global obstruction: all three classes remain **UNKNOWN**.

Write old coordinates as new coordinates times the recorded matrix \(S\).
Every \(S\) is integral with determinant \(\pm1\), and its fourth column is
\((1,0,0,0)^t\). Thus the old ordinate \(W\) is precisely the first new
coordinate. A rational point can be written primitively as

\[
[X:Y:Z:W]=[s^2:st:t^2:y],\qquad \gcd(s,t)=1,\qquad y\in\mathbf Z.
\]

Indeed, an integral conic point initially has a common factor multiplying
\((s^2,st,t^2)\); the quartic equation forces that factor to divide \(W\)
too, so primitivity removes it. Unimodularity preserves primitivity.
Consequently a new coordinate height bound \(B\) implies

\[
|y|\le B,\qquad t^2\le B\|S_{*,3}\|_1.
\]

Here column 3 denotes the old \(Z\) coordinate, counting from one.
For \(t\ge1\), set \(f(x)=F(x,1)\). Necessarily
\(0\le f(s/t)=y^2/t^4\le B^2\).

| Class | Presentation | Bound on \(|t|\) | Exact exclusion |
|---|---|---:|---|
| \(r_1\) | original | 4472 | \(f(x)>10^{30}\) for every real \(x\) |
| \(r_1\) | translated by \(Q\) | 3162 | Four rationally isolated bands, none containing an allowed denominator |
| \(r_2\) | original | 3162 | \(f(x)>10^{33}\) for every real \(x\) |
| \(r_2\) | translated by \(Q\) | 3162 | Four rationally isolated bands, none containing an allowed denominator |
| \(r_1+r_2\) | original | 3162 | Four rationally isolated bands, none containing an allowed denominator |
| \(r_1+r_2\) | translated by \(Q\) | 3162 | \(f(x)>10^{31}\) for every real \(x\) |

For the three uniform bounds, rational Sturm sequences prove that
\(f-10^k\) has no real zero, and its value at zero is positive.
For the other three models, disjoint rational intervals enclose all roots
of \(f\) and \(f-B^2\). Exact root counts and signs prove that
\(0\le f\le B^2\) is contained in these intervals. For every interval
\([l,h]\) and every \(1\le t\le3162\), the verifier checks
\(\lceil lt\rceil>\lfloor ht\rfloor\): 37,944 integer inequalities total,
without square testing. The six leading coefficients exceed \(B^2\)
and are nonsquares, also excluding parameter infinity.

The recorded `PointsQI` calls used the default `ExactBound := false`,
so the implementation can also find points outside its nominal box.
This certificate does not describe or exclude that additional exploration.
It does show why the numerical height argument was a poor measure of
progress. Minimal discriminant did not give useful coordinate scaling.

```sh
sage -python elliptic-curves/cas/audit_fixed_field_radical_search_geometry.py
sage -python -m unittest elliptic-curves/tests/test_fixed_field_radical_search_geometry.py
```

Regeneration adds `--write-certificate`. Algebraic root isolation only
proposes interval endpoints; verification uses rational arithmetic.
The six original model/evidence files remain unchanged.

The revised method is to reuse the two-dimensional structure before
requesting more arithmetic:

1. **Keep Q-translations as optional point-search charts.** In the Kummer
   exact sequence, \(r\) and \(r+\delta(Q)\) have the same image in
   \(\Sha[2]\). Their solubility is equivalent by translation of the
   covering map, and their Cassels--Tate pairings agree. New translations
   cannot provide an additional obstruction. Once level zero is certified,
   further finite-prime minimization also cannot improve the level.
2. **For a further obstruction, target only two pairing entries.** A
   certified auxiliary Selmer class \(s\) with
   \((\langle r_1,s\rangle,\langle r_2,s\rangle)\ne(0,0)\) excludes
   exactly two of the three targets by bilinearity. Such an \(s\) must lie
   outside \(W+\langle\delta(Q)\rangle\). No repeated pairing inside
   that known space helps. Finding this witness remains an unresolved
   construction problem; a full ambient Selmer basis is not necessary
   to certify a witness once found.
3. **For points, construct an actual lift before searching four-covers.**
   The substitution \([s:t:y]\mapsto[s^2:st:t^2:y]\) is an embedding of
   the same curve; it is not a further two-covering of it. The existing
   map to \(E\) still has degree four. An actual four-cover over the
   target instead maps to \(E\) with degree sixteen. Such a lift, with
   its exact map, is the missing input to a higher-descent point search.
4. **Share higher-descent arithmetic across the generators.**
   [Fisher's four-descent addition method, Sections 6 and 8](https://www.dpmms.cam.ac.uk/~taf1000/papers/fourdesc.pdf)
   constructs sums of existing lifts without computing every quartic
   field's class group. Once lifts of two independent target classes
   exist, their sum supplies the third. A lift of \(\delta(Q)\) comes
   from the known rational point, allowing translated lifts to be reused
   as well. General addition still requires algebra computations; this
   is not a claim that it is cheap or already implemented here. The hard
   initial lift is not supplied by the restricted radical calculation.

For positive discovery, an individual exact norm witness in the chosen
quartic algebra can suffice to construct a candidate lift. A full class
group is not a prerequisite for certifying a point subsequently found
on it. A failure to find such witnesses in a partial subgroup, however,
cannot certify insolubility. Likewise, elliptic-logarithm searches using
only multiples of \(Q\) cannot discover these classes: modulo two that
subgroup contains only \(0\) and \(\delta(Q)\). A completeness argument
would need additional Mordell--Weil information that is currently absent.

Every eventual point still requires the retained exact map back to its
original Kummer mask and the corresponding independence certificate.
No new large computation was launched for this audit.

## Reproduction

The compact [summary](../../artifacts/generated-results/elliptic-curves/fixed_cubic_u_minus1_cassels_tate_v1.json)
contains the matrix, radical, symplectic pairs and three search quartics.
The compressed [arithmetic evidence](../../artifacts/generated-results/elliptic-curves/fixed_cubic_u_minus1_cassels_tate_evidence_v1.json.gz)
contains every exact witness. It is not an attestation-only matrix audit.

```sh
sage -python elliptic-curves/cas/verify_fixed_cubic_cassels_tate.sage --check
sage -python -m unittest elliptic-curves/tests/test_fixed_cubic_cassels_tate.py
```

Regeneration, separate from the cheap replay, uses:

```sh
sage -python elliptic-curves/cas/run_fixed_cubic_cassels_tate.sage \
  --workdir artifacts/local/fixed-cubic-ct-fresh \
  --entry-seconds 60 --campaign-seconds 900 --point-height 10000
```

The producer verifies the complete matrix before scheduling the radical
covers. Add `--write-certificate` to also run the eight additional controls
and replace the compact certificates. Different valid reduction witnesses
may change individual quartic presentations; the pairing matrix is the
arithmetic invariant being certified.
