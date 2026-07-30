# Marked multi-pole peak-reduction experiment

## Status

The bounded experiment found no counterexample.  Its factor traces expose a
uniform local pole-change inequality, proved below and recorded canonically
as Proposition 6.6 of
[`JC2_GLOBAL_COX_PACKET_ATTACK.md`](JC2_GLOBAL_COX_PACKET_ATTACK.md).
Consequently the marked multi-pole peak-reduction conjecture is proved for
reduced alternating Jung words, meaning complete triangular factors of
degree at least two.

The operational statement tested is:

> If an alternating triangular word lowers the final marked pole height,
> then its initial coordinate pair admits a complete triangular polynomial
> shear which lowers the height immediately.

Thus a counterexample is an initial pair which is terminal for the complete
single-shear compiler but is lowered by an alternating word.

## Exact state and ledger

The search works in \(\mathbb Q(t)\).  A state consists of

\[
 (u,v;\ p_1,\ldots,p_s;\ \sim),
\]

where the \(p_i\) are the marked boundary valuations and \(\sim\) is a
fixed conductor-pairing relation.  The two charts use

\[
 \{0,\infty\}
 \quad\hbox{and}\quad
 \{0,1,\infty\}.
\]

At every state and after every factor the checker computes:

1. both exact orders and pole orders at every marked boundary;
2. both initial coefficients in the local parameters \(t-p_i\), and in
   \(t^{-1}\) at infinity;
3. the unchanged conductor pairing; and
4. the peak
   \[
   H(u,v)=\sum_i
   \bigl(\max(0,-\nu_i(u))+\max(0,-\nu_i(v))\bigr).
   \]

The pairing is deliberately retained as a marked decoration.  Polynomial
coordinate automorphisms do not change it, so it does not affect the
peak calculation itself.  A future JC(2) endpoint compiler must combine
this automorphism calculation with the actual finite-normalization pairing
and meridian data.

For a shear \(v\mapsto v+c u^d\), cancellation at \(p_i\) is possible only
when

\[
 \operatorname{pole}_{p_i}(v)
 =d\,\operatorname{pole}_{p_i}(u),
 \qquad
 c=-\frac{\operatorname{in}_{p_i}(v)}
          {\operatorname{in}_{p_i}(u)^d}.
\]

The complete-shear test recursively follows all such forced coefficients in
strictly descending degree.  It therefore includes multi-term factors whose
leading prefix preserves or raises \(H\).

## Results

The main all-word certificate uses seed basis elements of pole order at most
two, at most two nonzero basis terms per seed function, factor degrees
\(1,2,3\), and alternating length two.  It checks both affine triangular
peaks and nonlinear second factors.

| chart | terminal initial pairs | words from terminal pairs | all words | globally lowering | delayed first step | delayed with initial complete reduction |
|---|---:|---:|---:|---:|---:|---:|
| two-pole Laurent | 104 | 7,488 | 15,120 | 348 | 40 | 40 |
| three-pole rational | 480 | 34,560 | 54,432 | 694 | 126 | 126 |

The largest observed delayed peak is one.  The three-pole trace in the JSON
artifact includes the genuine sequence

\[
 4\longrightarrow5\longrightarrow3.
\]

At the initial state the complete compiler already exposes a shear lowering
\(4\) to \(3\).  All 166 delayed paths have the same qualitative behavior.
Their factor-degree rows are \((1,1)\) and \((1,2)\); no delayed word in this
grid has both Jung degrees at least two.

Three complementary terminal searches were run.

| seed/factor grid | two-pole terminal words | three-pole terminal words | counterexamples |
|---|---:|---:|---:|
| extended pole-order-two seeds; reduced degrees \(2,3\); length at most 3; polydegree at most 27 | 16,640 | 76,800 | 0 |
| signed coefficients on the small basis; reduced degrees \(2,3\); length at most 3 | 26,880 | 26,880 | 0 |
| every degree-\(2/3\) factor with lower coefficients in \(\{-1,0,1\}\); length 2 | 5,760 | 5,184 | 0 |

These are bounded computations and are not used to establish the
degree-independent theorem below.  They also do not turn the conductor
decoration into the conductor map of a hypothetical Keller normalization.

## Inductive invariant and proof

Consider two consecutive opposite factors

\[
x'=x+P(y),\qquad y'=y+Q(x'),
\]

of degrees \(d,e\ge2\).  At one marked valuation write \(a,b,a',b'\)
for the four relevant pole orders and set

\[
\alpha=a'-a,\qquad\beta=b'-b.
\]

Then

\[
\boxed{\beta\ge\alpha.}
\]

If \(\beta<0\), cancellation in the second factor forces \(b=ea'\).
Since \(db=dea'>a'\), the first factor can have output pole \(a'\) only
when \(a=db\) and its leading terms cancel.  Hence

\[
a-a'=(de-1)a',
\qquad b-b'\le ea'\le(de-1)a',
\]

which gives \(\beta\ge\alpha\).  If \(\beta\ge0\), only
\(\alpha>0\) needs consideration.  In that case \(a'=db>a\), and then
\(ea'=edb>b\), so \(b'=ea'\).  Therefore

\[
\beta=(ed-1)b\ge db\ge db-a=\alpha.
\]

Let \(\Delta_{j,i}\) be the pole change made by factor \(j\) at marked
valuation \(p_i\).  Applying the inequality to every consecutive pair
gives

\[
\Delta_{1,i}\le\Delta_{2,i}\le\cdots\le\Delta_{r,i}.
\]

After summing over \(i\), the complete-factor height increments are
nondecreasing.  A word with negative total increment must therefore have
negative first increment.  Thus the first complete triangular factor
already lowers the current marked multi-pole peak.

This proof is independent of the number of marked poles, residue fields,
initial coefficients, and conductor pairing.  Residues remain necessary to
construct the actual complete cancellation polynomial; the pairing remains
necessary for the later JC(2) conductor and monodromy filters.  Degree-one
affine factors are not entries of reduced Jung polydegree and are handled
separately by affine basis reduction.

## Reproduction

The practical all-word certificate is:

```bash
.venv/bin/python plane-jc/cas/search_marked_multipole_peak.py \
  --max-degree 3 \
  --max-length 2 \
  --include-linear \
  --extended-seeds \
  --max-seed-terms 2 \
  --scan-all \
  --output artifacts/generated-results/marked_multipole_peak_search.json
```

The reduced Jung polydegree-\(27\) terminal search is:

```bash
.venv/bin/python plane-jc/cas/search_marked_multipole_peak.py \
  --max-degree 3 \
  --max-length 3 \
  --extended-seeds \
  --max-seed-terms 2 \
  --output \
    artifacts/generated-results/marked_multipole_peak_reduced_degree27.json
```

The signed-seed and complete-factor checks are:

```bash
.venv/bin/python plane-jc/cas/search_marked_multipole_peak.py \
  --max-degree 3 --max-length 3 --signed-seeds \
  --output artifacts/generated-results/marked_multipole_peak_signed_degree27.json

.venv/bin/python plane-jc/cas/search_marked_multipole_peak.py \
  --max-degree 3 --max-length 2 --complete-factors \
  --output artifacts/generated-results/marked_multipole_peak_complete_factors.json
```

The generated-orbit falsification checks are:

```bash
.venv/bin/python plane-jc/cas/search_marked_multipole_peak.py \
  --max-degree 3 --max-length 3 \
  --extended-seeds --max-seed-terms 2 --orbit-only \
  --output artifacts/generated-results/marked_multipole_peak_orbit_degree27.json

.venv/bin/python plane-jc/cas/search_marked_multipole_peak.py \
  --max-degree 3 --max-length 2 --include-linear \
  --extended-seeds --max-seed-terms 2 --four-pole --orbit-only \
  --output artifacts/generated-results/marked_multipole_peak_orbit_four_pole.json
```

The first generates 154,560 high-complexity reduced endpoints through
polydegree \(27\).  Every last factor raises the peak, so its inverse is an
immediate lowering factor.  The second includes affine factors and four
marked poles.  Among 160,272 tested forward words, 4,942 endpoints above
their base have a nonincreasing last factor; every one admits a different
complete lowering shear.
