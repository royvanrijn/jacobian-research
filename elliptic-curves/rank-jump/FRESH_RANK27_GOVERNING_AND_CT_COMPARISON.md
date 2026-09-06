# Fresh rank-27 fibres: a masked governing/CT comparison

Follow-up: [the fixed103b2 high fibre must contain at least nine strict rational
directions outside its generic subgroup](MATCHED103B2_JUMP_REQUIRES_NINE_STRICT_DIRECTIONS.md).
The new bound uses an equation-defined reciprocity constraint; it does not
fill the missing independent class basis below.

The 16-fibre panel gives two negative controls, but **does not yet compute
CT on the additional Selmer quotient**. Exceptional points and classes
obtained from them were excluded from the arithmetic inputs.

* An identically selected pair of generic sections has an explicit
  degree-eight governing polynomial with splitting-field degree **192 on
  every fibre**, including ordinary controls. This structure does not
  discriminate jumps in this panel.
* Large inherited strict blocks and large simultaneous −1-twist CT
  switches occur on low-gain fibres. Three controls have strict dimensions
  **9, 10, 5** and CT-switch ranks **8, 10, 4**. Several fresh +10 fibres
  have inherited strict dimension zero. These measurements cannot explain
  their large exceptional quotients.

The independent cubic class-group attempts did not supply the missing
additional classes. Thus the requested comparison on an independently
constructed **exceptional quotient** remains **UNKNOWN**, rather than
being filled with the earlier point-derived matrices.

The [comparison JSON](../../artifacts/generated-results/elliptic-curves/rank_jump_fresh_governing_comparison_v1.json)
and [CSV](../../artifacts/generated-results/elliptic-curves/rank_jump_fresh_governing_comparison_v1.csv)
are the compact outputs. This note changes no search protocol or mathematical
status entry.

## Panel and accounting

Selection uses the completed V18 inventory: all eight rank-at-least-27
curves, five distinct same-family low-gain controls, and historic curves
356, 385 and 398. The
[manifest](../../artifacts/generated-results/elliptic-curves/rank_jump_fresh_governing_panel_manifest_v1.json)
binds the exact source bytes and records matching gaps and exposure.

The rank column is a **certified lower bound**, not an exact full rank.
The increment is the recorded bound minus the marked generic rank m;
it is relative to that specified subgroup. Generic mod-two independence
is checked anew from its sections. The source rank certificates supply
the retained high-rank labels; the arithmetic worker does not inspect
those exceptional-point certificates. The historic first-17 identification
uses the retained exact R17 lineage certificate; curve398 uses its recovered
MW16 parent. We do not treat a convenient arbitrary reference subgroup as
a generic family.

Let G be the marked generic Kummer subgroup, and let U be the norm-square
cubic classes which are locally square at 2, infinity and all bad primes,
and unramified outside those places. The table measures

\[
 k=\dim(G\cap U),\qquad
 b=\operatorname{rank}\bigl(\mathrm{CT}_{E^{(-1)}}|_{G\cap U}\bigr).
\]

An interval is a proved bound from incomplete local coverage. It is not a
measured zero. Every row has a degree-192 governing field for its fixed
first generic pair; that pair need not itself lie in U.

| Case | Family / parameter | Rank ≥ | m | Recorded increment | k | b |
|---|---|---:|---:|---:|---:|---:|
| new-40 | 074d9, 2818/1535 | 27 | 17 | +10 | 0–5 | 0–4 |
| matched low | 074d9, 2824/885 | 17 | 17 | 0 | 9 | 8 |
| new-71 | 103b2, 3726/881 | 27 | 17 | +10 | 0 | 0 |
| matched low | 103b2, −1049/2296 | 17 | 17 | 0 | 0 | 0 |
| new-41 | 11952, −2448/11 | 27 | 17 | +10 | 0 | 0 |
| matched low | 11952, −1171/1683 | 17 | 17 | 0 | 10 | 10 |
| new-188 | 11952, 110314/102227 | 27 | 17 | +10 | 0 | 0 |
| matched low | 11952, 130349/28916 | 17 | 17 | 0 | 5 | 4 |
| new-72 | 11952, 2012/211 | 27 | 17 | +10 | 0–5 | 0–4 |
| new-48 | 11952, 2828/2015 | 27 | 17 | +10 | 0 | 0 |
| new-186 | 11952, 4286/1881 | 27 | 17 | +10 | 0–5 | 0–4 |
| new-90 | a1-fibration-01, −1867/270 | 27 | 16 | +11 | 0–6 | 0–6 |
| matched low | a1-fibration-01, −3187/3697 | 17 | 16 | +1 | 0–7 | 0–6 |
| ICARM356 | published R17 lineage | 29 | 17 | +12 | 0–2 | 0–2 |
| ICARM385 | published R17 lineage | 29 | 17 | +12 | 0 | 0 |
| ICARM398 | recovered curve398-p16875 MW16 | 30 | 16 | +14 | 0 | 0 |

The upper bound on b uses alternation: b≤2 floor(k_upper/2). It does not
claim the uncomputed form attains that bound. The low labels are censored
search outcomes; in particular “0” is not an exact-rank-17 theorem.

Controls are matched within family by coefficient-bit difference, then
parameter-height-bit difference. The four smaller11952 highs reuse one
control; they are not four independent negative controls. MW16 has only
an observed +1 control in the retained low cohort. Historical rows are
calibrations, not same-parameter-scale matches to the fresh families.

There is one recorded implementation deviation: the export selected the
large11952 control pool when parameter-height **bit length exceeded13**,
whereas the prose protocol said height exceeded4096. Consequently
4286/1881 used the compact pool. No case was rematched after observing its
arithmetic. This is a descriptive panel, not a preregistered aggregate test.

## What was computed without exceptional points

The [masked inputs](../../artifacts/generated-results/elliptic-curves/rank_jump_fresh_governing_panel_inputs_v1.json)
contain only case tokens, equations and marked generic sections. Labels,
rank annotations, search outputs, exceptional coordinates and their Kummer
classes are absent. Export uses the source certificates to select and
project cases; the arithmetic processes read only the projection and their
own subsequent checkpoints. Thus selection is retrospective, while the
arithmetic dependency is independent of exceptional points.

The [base computation](../../artifacts/generated-results/elliptic-curves/rank_jump_fresh_governing_panel_v1.json)
reconstructs the cubic, verifies the generic classes, computes local
characters and their joint kernel, and constructs half ideals for each
complete strict basis. It factors no point numerator to construct those
ideals. For a generic point written x=a/d², y=b/d³, it starts with

\[
 \gamma=a-d^2\theta,\quad N(\gamma)=b^2,\quad (b,\gamma),
\]

then corrects the ideal at discriminant primes and verifies the exact
identity \(\mathfrak J_\beta^2=(\beta)\) for each strict product β.

The Artin matrix is

\[
 A_{ij}=\chi_{\beta_i}([\mathfrak J_{\beta_j}]).
\]

The [small completion](../../artifacts/generated-results/elliptic-curves/rank_jump_fresh_governing_completion_v1.json)
repairs nonunit residues at the small primes5 and11 with exact local
Frobenius evaluations, checked separately against generic section
characters. It also proves zero full strict kernel when a subset of bad
places already has zero kernel. Otherwise it retains only an upper bound.

On this inherited rational strict subspace, the original CT form vanishes.
The [previously proved scalar-cup identity](INDEPENDENT_SCALAR_CUP_AND_TWIST_BLOCKS.md)
therefore gives the actual transported form

\[
 \mathrm{CT}_{E^{(-1)}}=A+A^{\mathsf T}.
\]

The low controls decompose as

\[
 074d9:\ H^{\oplus4}\oplus0_1,\qquad
 11952_{\rm small}:\ H^{\oplus5},\qquad
 11952_{\rm large}:\ H^{\oplus2}\oplus0_1.
\]

Here H is an alternating plane. Their restricted radicals have dimensions
1,0,1. Rational classes on the twist must lie in these radicals, but a
restricted radical does not prove rational solubility. These are not
whole-twist rank bounds without a full Selmer group.

The [octic certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_fresh_governing_octics_v1.json)
uses the first two generic sections on every masked model. With c=x_Q−x_P,

\[
 h(T)=T^8-4(y_P+y_Q)T^6+6c^2(x_P+x_Q)T^4
       -4c^3(y_Q-y_P)T^2+c^6.
\]

All16 cubics have S3 Galois group, and each selected pair is independently
nonzero modulo2. The
[cochain and faithful group argument](EXPLICIT_PAIR_GOVERNING_OCTIC.md)
then gives joint class-field degree96 and governing degree192. Each explicit
octic has nonzero discriminant. For every admitted inert prime≤199 the
factor pattern and independently evaluated radical norm give the same
Frobenius bit: (1,1,3,3) for0, (2,6) for1.

These bits are not automatically CT values for another family parameter.
They require the appropriate local twist hypotheses. The polynomial also
depends on section representatives and the chosen model. Its degree-eight
minimal governing encoding is not a minimal rational-point carrier for a
jump block. No native-cover count or four-cover-lift count is inferred from
these inherited class counts.

## Paired findings

1. **103b2, +10 versus observed0.** Coefficient sizes differ by only two
   bits. Both inherited strict kernels are zero and both governing fields
   have degree192. The computed data do not distinguish this pair.
2. **11952 at large parameter height, +10 versus observed0.** Both models
   have304-bit short coefficient scale and the same parameter-height bit
   length. The high's inherited kernel is zero; the low's has dimension5
   and CT rank4. The large inherited obstruction block occurs on the low.
3. **074d9, +10 versus observed0.** Parameter-height bit lengths agree;
   coefficients differ by17 bits. Even the high's incomplete local audit
   bounds k≤5 and b≤4, below the low's exact k=9 and b=8. Completing the
   high's remaining factorization cannot reverse that inequality.
4. **11952 at smaller height.** The +10 curve −2448/11 has k=0; its
   observed-zero control has k=10 and CT rank10. This is a ten-dimensional
   rational generic block on the original low-gain curve which is completely
   excluded from rationality after the −1 twist. Simultaneous solubility
   switches exist without a large observed specialization quotient.
5. **MW16 +11 and historic controls.** The fresh MW16 comparison remains
   incomplete (k≤6 versus k≤7). The historic +12/+14 rows have k≤2,0,0.
   They provide no evidence for a large inherited strict block as the
   source of an extreme quotient.

## Exactly what implication is missing

The desired point-independent object is a supply of classes **beyond G**,
with local admissibility established independently of exceptional points.
Only then can the governing/CT calculation measure their shared obstruction.
In a fixed quadratic-twist control the cubic algebra and selected classes
can be held fixed. In these K3 family comparisons the cubic fields themselves
change with t; one cannot identify their additional class spaces merely
because every generic pair has the same abstract degree192 group.

On every row the independent additional Selmer basis, its CT form, the full
CT radical and the minimal simultaneous-solubility carrier remain UNKNOWN.
Seven attempts at complete discriminant factorization exceeded30 seconds. The other
nine equation-only class-group attempts exhausted the fixed256MiB PARI stack.
No class-group upper bound, Sha dimension or new rational point follows.
The separate small-kernel calculation supplies useful information without
rerunning those failed campaigns.

The missing chain is therefore still

\[
 t\ \Longrightarrow\ \text{independent additional globally admissible classes}
 \ \Longrightarrow\ \text{vanishing higher obstructions}
 \ \Longrightarrow\ \text{enough rational lifts}.
\]

Neither a universally soluble cubic norm cochain nor an inherited zero CT
form closes either of the last two implications. Zero CT can leave
higher-divisible Sha; rational solubility needs a further argument.

## Priorities after this panel

1. **Most useful unresolved mechanism — solubility.** A change in higher
   obstructions on an independently constructed additional class block.
   This remains a candidate mechanism, not a finding established by these
   production data.
2. **Required incidence computation.** Obtain a certified additional
   class basis modulo G, or useful bounds, from the cubic S-class arithmetic.
   A focused next target is the103b2 pair above: its same-family inherited
   baselines agree, so any discriminator must occur elsewhere. Before a
   larger class calculation, establish a feasible relation/character method
   and an explicit resource gate.
3. **Weak or refuted positive proxies.** Pair governing degree192 is
   universal under the checked hypotheses. Large inherited strict dimension
   and large inherited CT-switch rank occur on ordinary controls. A zero
   inherited kernel also occurs on an ordinary control. None is a selector.
4. **Falsifiable next test once that basis exists.** On the fixed103b2 pair,
   compute the same additional-class incidence and CT obstruction data with
   labels masked. If the quotient dimensions and obstruction ranks agree,
   this CT-level mechanism fails to discriminate that pair; higher descent
   or a rational construction must carry the difference. No new parameter
   search is required.
5. **Information for Agent1 now.** No candidate-scoring change is justified.
   The cheap zero-kernel certificate can eliminate unnecessary work on an
   inherited strict block, and the separate masked input format makes a
   later quotient computation auditable. Neither is a visibility feature
   promoted to a rank predictor.

## Verification and bounded replay

The [independent verification](../../artifacts/generated-results/elliptic-curves/rank_jump_fresh_governing_panel_verification_v1.json)
passes16 generic-independence/kernel replays,895 local-power checks,24 exact
ideal-square certificates,206 independent Jacobi evaluations and205
CAS-independent modular octic/radical replays. The octic discriminants use
an independent integer Sylvester determinant. The finite192-element group
and all its action/cocycle identities are replayed separately.

The three local protocols cap workers at30–60 seconds, one worker at a time;
no point search, new parameter or active-search mutation occurs. An initial
Sage/Fraction interoperability error was corrected before the final local
certificates; completed equation-only work and genuine resource-limit
failures were retained, with the initial transcripts in ignored checkpoints.

From the repository root:

```sh
sage -python elliptic-curves/rank-jump/verify_fresh_governing_panel.py check
python3 elliptic-curves/rank-jump/summarize_fresh_governing_panel.py check
```

The capture commands in the three scripts write new outputs exclusively.
They preserve capped failures instead of silently treating them as arithmetic
zeros. Shared navigation, active experiments, search outputs and
MATH_STATUS.json are untouched by this work.
