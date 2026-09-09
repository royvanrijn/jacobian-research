# Local non-incidence versus global branch dependence: an explicit separation

## Results and boundary

**New construction, independently verified.** There is an explicit infinite
sequence of orbit8044 conic points such that:

- every term gives17 inherited points plus one independent rational point;
- at all25 primes in the old dependent-conic footprint, the elliptic curve
  and all18 marked point reductions equal those of the dependent split;
- the conic parameters converge over the reals to that dependent split;
-19 disjoint primes uniformly certify rank18 for every integer term.

Thus the original finite footprint has rank17 on every18-column point
list in this family, although each list really has rank18. This is an exact
construction, not a sampled correlation or a numerical height comparison.

**Verified local obstruction.** In contrast, the same conic is nonsplit
throughout explicit local neighbourhoods of302 and all eight fixed null
addresses. At302 one such neighbourhood in the **original parent
coordinate** is

\[
\boxed{v_{29}(t)\ge2\quad\Longrightarrow\quad
          \text{orbit8044 does not split over }\mathbf Q.}
\tag{1}
\]

**Scope.** These statements concern the specified conic and marked branch.
They do not prove rank17 for the limiting elliptic fibre, exclude other
seed constructions at any control, or produce302's seed. The known302
seed is already evidence that conic non-incidence is not a rank bound.
No running search or factory was changed.

## 1. Exact construction

**Verified inherited inputs.** Reuse the old conic slope coordinate `u`,
not the later factory coordinate. The generic parametrization is

\[
T(u)=t_c+\frac{q'(t_c)-2s_cu}{u^2-q_2},\qquad
S(u)=s_c+u(T(u)-t_c),\qquad S(u)^2=q(T(u)).
\]

Its coefficients, rational point `(t_c,s_c)`, residual quadratic lift and
all18 rational section maps are in the unchanged
[old input packet](../../artifacts/generated-results/elliptic-curves/det1092_conic_seed_progression_v1/input.json).
That packet already has an independent generic map and `u=0` rank18 proof.

**Retrospective negative anchor, not an exceptional point oracle.** The
smaller rational conic preimage of the dependent reduced parameter
`s=-528/3635` is

\[
r=-\frac{52073620095012194896969766919966900215391}{198893}.
\]

The reconstructed branch at `r` satisfies
`Q_r=M10-M14-M15` in the one-based inherited17 ordering. The checker
reconstructs this from the conic formula and verifies the rational group
identity, rather than inferring dependence from finite rank17.

Let

```text
S = {7,13,19,29,31,37,43,47,53,59,67,73,79,83,89,97,
     101,103,109,113,127,131,137,139,149}

Tproof = {181,193,199,223,227,229,239,241,251,257,269,
          271,277,281,283,307,311,313,317}.
```

Define the positive integers

\[
\begin{split}
L&=\left(\prod_{p\in S}p\right)
       7^{24}13^{24}29^{24}31^{60}101^{12},\\
P&=\prod_{p\in T_{\rm proof}}p
  =4554842988330072766085778231028493113035977947.
\end{split}
\]

Write `r=a/b` in lowest terms, with `b>0`. The unique integer
`0<q0<LP` satisfying

\[
q_0\equiv1\pmod L,
\qquad aq_0+bL\equiv0\pmod P
\tag{2}
\]

exists because the checker verifies `gcd(P,ab)=gcd(P,L)=1`.
No coprimality of `ab` with `L` is asserted. The second congruence makes
`q0` a unit modulo `P`; the first ensures it is nonzero modulo `L`.

**New explicit seed sequence.** For every integer `n`, put

\[
\boxed{u_n=r+\frac{L}{q_0+LPn},\qquad t_n=T(u_n).}
\tag{3}
\]

Evaluate the old18 section maps at `u_n`. All are defined on a smooth
rational elliptic fibre and give a certified rank18 subgroup. The denominator
in(3) is never zero, since `0<q0<LP`. The `u_n` are distinct and tend to
`r` as `n` tends to either infinity. Since `T` has degree2, the construction
gives infinitely many distinct original parameters `t_n`, tending over the
reals to the dependent parameter `T(r)`.

The [new construction packet](../../artifacts/generated-results/elliptic-curves/det1092_local_shadow_seeds_v2/construction.json)
exports all constants and integer polynomial maps. A fixed rational scaling
of the original short model makes the limiting equation exactly the old
dependent-control equation; this scaling is explicit and verified.

## 2. Why every term is certified, not merely eventually independent

**New constructive application of congruence specialization.** In homogeneous
coordinates the sequence is

\[
[U_n:V_n]=[a(q_0+LPn)+bL:b(q_0+LPn)].
\]

Equation(2) gives

\[
(U_n,V_n)\equiv(a,b)\pmod L,
\quad U_n\equiv0\pmod P,
\quad V_n\equiv bq_0\not\equiv0\pmod p
\ (p\mid P).
\tag{4}
\]

Every section is exported as a primitive integral polynomial triple, and
each coefficient or base map as an integral polynomial pair. At a control
prime `p`, let `m` be the common valuation in a point triple, or the
denominator valuation for a coefficient pair. The exponent of `p` in `L`
is strictly greater than every relevant `m`. Consequently(4) preserves
the normalized marked reductions modulo `p`, including points reducing
to infinity. The larger exponents at five primes compensate for exact
polynomial-map contents; they are conservative bounds, not optimized costs.

At the proof primes, every relevant denominator at `[0:1]` is a unit
and every point triple has nonzero reduction. Homogeneity and(4) therefore
give precisely the reductions of the old positive `u=0` specialization.
These nonzero reductions also rule out a rational pole, undefined point
triple or singular elliptic fibre for any integer `n`.

**Independent exact rank computation.** The verifier enumerates the complete
finite elliptic groups and their doubled subgroups at the proof primes,
without using the constructor's cached cubic characters. Their quotient
matrix has18 independent columns. At181 the finite group has odd order,
excluding rational2-torsion. Any integral relation would thus have all
coefficients even, and repeated division proves independence.

The same independent finite-group procedure gives rank17 on the18 columns
at the25 control primes, with rank17 on the inherited prefix. The actual
point reductions, not just their ranks, are retained and coincide with
the dependent anchor's reductions for every term of(3).

## 3. A precise obstruction to fixed local predictors

**New deduction.** No test based only on these25 reduced elliptic curves and
marked point configurations can correctly decide whether the conic branch
is in the inherited rational span. It receives identical data on the
dependent anchor and every independent term of(3).

The example extends to any prescribed finite precision at these same primes.
To preserve data modulo `p^h`, replace `L` by

\[
L_h=L\prod_{p\in S}p^{h-1}
\]

and solve the same CRT conditions with `L_h`. The disjoint proof primes
remain usable and certify every term. This is a symbolic extension of the
same congruence proof, not a new computed precision survey.

**Established literature and broader consequence.** The generic18 sections
are independent and the family is nonisotrivial. The
[specialization theorem](https://swc-math.github.io/aws/2024/2024SilvermanLecture4Slides.pdf)
therefore makes all but finitely many rational conic addresses productive.
The completed [density/height note](DET1092_SEED_DENSITY_AND_LIMITING_LATTICE_2026-09-08.md)
already verifies its hypotheses. Every finite adelic neighbourhood of `r`
contains infinitely many rational addresses, so it contains productive ones.
Hence no fixed finite-place, finite-precision test can give a neighbourhood
certificate of *branch dependence* at `r`. The explicit construction above
is stronger for its specified panel: it certifies every term, with no unknown
exceptional bound.

**Important exclusions.** This does not obstruct a positive independence
certificate at finitely many suitable primes—the proof itself supplies one.
It does not obstruct adaptive prime selection, global rational identities,
exact halving, all-place conditions, or information incorporating the exact
rational coefficients rather than only their fixed-precision reductions.
In particular the old halving-or-cycle classifier is not contradicted:
finite rank failure was never its dependence proof.

## 4. Local non-incidence is different and genuinely detectable

**Verified application on the unchanged roster.** For each old address `tau`,
compute `q(tau)`. The first obstructing odd prime at most179 supplies the
following exact square obstruction. Here `v=v_p(q(tau))`, and `c` is the
leading unit modulo `p`; if `v` is odd or `c` is nonsquare, the conic fibre
cannot split over `Q_p`, hence cannot split over `Q`.

| Address | p | v | Leading unit c | v_p(t-tau) sufficient for the same obstruction |
|---|---:|---:|---:|---:|
| scale-0131232 |61|0|6 nonsquare|1|
| scale-0257585 |29|2|12 nonsquare|2|
| scale-0487239 |3|5|1|4|
| scale-0177036 |29|3|28|3|
| scale-0043332 |43|0|19 nonsquare|1|
| scale-0590501 |29|2|12 nonsquare|2|
| scale-0290097 |43|0|42 nonsquare|1|
| scale-0748009 |43|0|39 nonsquare|1|
|302, original t=0 |29|2|21 nonsquare|2|

The uniform ball statement is checked without factorization:

\[
q(\tau+\delta)-q(\tau)=q'(\tau)\delta+q_2\delta^2.
\]

The recorded lower bound on `v_p(delta)` makes both terms have valuation
greater than `v`. Thus the obstructed leading squareclass remains unchanged.
The checker verifies each unit using both direct square-residue enumeration
and Euler's criterion; all nine cases complete within the fixed prime pool.

**Meaning for the original controls.** Orbit8044 fails at302 and the eight
nulls already at local incidence. At the dependent conic address it splits,
but the chosen rational branch is inherited. At the three productive conic
addresses—including the first seed of the rank21 cascade—it splits and
its branch is independent. The nine blinded MW16 recoveries change their
classification when the reference subgroup becomes full MW17; they are
not evidence that local splitting alone identifies an exceptional direction.
No later cascade point is used here.

## 5. Remaining task and reproducibility

**Open endpoint.** The exact distinction is now constructive:

\[
\text{local cover obstruction}\quad\ne\quad
\text{marked branch dependence}\quad\ne\quad
\text{absence of other seeds on the elliptic fibre}.
\]

The family(3) supplies independent seeds from formulas alone, but uses
the same conic and therefore cannot reach302 or any of its eight nonsplit
control addresses. An equation-only choice of another useful302 member,
or an exact global marked-class criterion selecting it, is still missing.
No rank upper bound for any control follows from this work.

**Verified computation.** Constructor0.57s, independent replay0.67s, and
local-incidence check0.46s, each under25s. The constructor reuses old positive
prime exposures; the independent checker recomputes only small complete
finite groups. No new rational specialization was evaluated. A projective
versus affine coordinate-format preflight failure is retained in v1;
v2 changes that comparison only, not the mathematical rule.
The first replay's repeat check also exposed a JSON array-versus-tuple
comparison after all arithmetic checks had passed. Both reports are retained;
`independent-replay-v2.json` normalizes JSON types and passes repeated replay.

```bash
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/construct_det1092_local_shadow_seeds.sage
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/verify_det1092_local_shadow_seeds.sage
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/verify_det1092_conic_local_incidence.sage
```

- [Protocol and source hashes](../../artifacts/generated-results/elliptic-curves/det1092_local_shadow_seeds_v2/protocol.json).
- [Exact construction](../../artifacts/generated-results/elliptic-curves/det1092_local_shadow_seeds_v2/construction.json).
- [Independent uniform rank18 replay](../../artifacts/generated-results/elliptic-curves/det1092_local_shadow_seeds_v2/independent-replay-v2.json).
- [Nine local non-incidence balls](../../artifacts/generated-results/elliptic-curves/det1092_local_shadow_seeds_v2/local-nonsplitting.json).
