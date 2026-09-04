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
