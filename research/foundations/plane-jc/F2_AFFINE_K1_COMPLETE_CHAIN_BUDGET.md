# F2 `k=1` complete-chain point budget

> **Status.**  Exact conditional complete-chain identity.  On a completed
> F2 model whose only affine nonproperness component is one rational
> one-puncture `k=1` curve, let `u` be the number of sheets fixed by its
> geometric meridian.  After subtracting the full Cartier logarithmic
> determinant cycle and its conormal kernel degree, the global point budget
> is exactly `u-1`.  It is independent of geometric degree, carrier contact,
> source self-intersections, admissible boundary refinements, and the
> squarefree/double terminal row.  For the E8 cusp, every residue-degree-`f`
> ramified row requires point charge at least `2f`.  All `13` simple-inertia
> orbifold actions violate this budget.  Therefore they are excluded whenever
> all remaining point corrections are effective; otherwise a survivor must
> carry an explicitly quantified negative normalization/`Fitt_1` correction
> at an unresolved global attachment.  Establishing that final sign is still
> required for an unconditional exclusion of `(75,125)`.

The algebra and all orbifold rows are replayed by
[`verify_f2_affine_k1_complete_chain_budget.py`](../scripts/verify_f2_affine_k1_complete_chain_budget.py).

## 1. Complete determinant notation

Let

\[
 f:(X,D_X)\longrightarrow(Y,D_Y)
\]

be a common completed logarithmic model.  Put

\[
 L_X=K_X+D_X,\qquad L_Y=K_Y+D_Y,\qquad
 D_{\log}=L_X-f^*L_Y.                            \tag{1.1}
\]

The divisor `D_log` is the determinant divisor of the logarithmic
differential.  It includes the extraction-root cycle, every affine row, and
every exceptional component introduced by total transform.  It must not be
replaced by a list of strict transforms with their node matching discarded.

Suppose the unique affine nonproperness curve `C` has rational normalization
with one logarithmic puncture.  Let its source rows have transverse and
residue degrees `(e_i,f_i)`, and put

\[
 A=\sum_i e_if_i.                                \tag{1.2}
\]

If the meridian fixes `u` sheets, its nontrivial cycles use all other sheets,
so

\[
 \boxed{A=d-u.}                                  \tag{1.3}
\]

Write

\[
 c=L_Y\mathbin\cdot\bar C.                       \tag{1.4}
\]

The contracted extraction-root cycle has zero pairing with `f^*L_Y`.
Therefore projection along every affine row gives

\[
 f^*L_Y\mathbin\cdot D_{\log}=Ac.               \tag{1.5}
\]

Squaring (1.1) yields the complete, refinement-invariant determinant square

\[
 \boxed{
 D_{\log}^2=L_X^2-dL_Y^2-2Ac.}                  \tag{1.6}
\]

## 2. Conormal kernel and global cancellation

The immersed one-puncture conormal formula gives affine kernel degree

\[
 \sum_i e_if_i
 \left(L_Y\mathbin\cdot\bar C+1\right)
 =A(c+1).                                        \tag{2.1}
\]

The extraction-root kernel degree is zero.  Hence the complete generic
cyclic divisorial contribution is

\[
 \operatorname{ch}_2^{\rm div}
 =A(c+1)+\frac12D_{\log}^2.                     \tag{2.2}
\]

The global logarithmic Chern budget is

\[
 B_f=\frac12\left(L_X^2-dL_Y^2+2(d-1)\right).   \tag{2.3}
\]

Subtract (2.2), use (1.6), and cancel:

\[
\begin{aligned}
 B_f-\operatorname{ch}_2^{\rm div}
 &=\frac12\left(L_X^2-dL_Y^2+2(d-1)\right)
   -A(c+1)\\
 &\qquad-\frac12\left(L_X^2-dL_Y^2-2Ac\right)\\
 &=d-1-A\\
 &=\boxed{u-1}.                                  \tag{2.4}
\end{aligned}
\]

This is the desired degree-independent complete-chain identity.  It does not
use the F2 values of either log square or the carrier intersection.  Those
values merely give the checker a concrete regression.  In particular,
neither more negative strict transforms, further target centers, nor more
common-fan refinements enlarge the point budget.

## 3. Why raw negativity was the wrong variable

Write the complete determinant cycle as

\[
 D_{\log}=D_{\rm root}+D_{\rm aff},\qquad
 D_{\rm root}^2=54,                              \tag{3.1}
\]

and put

\[
 S=-D_{\rm aff}^2,\qquad
 I=D_{\rm root}\mathbin\cdot D_{\rm aff}.       \tag{3.2}
\]

Then

\[
 D_{\log}^2=54-S+2I.                            \tag{3.3}
\]

The old `4N` term retained only strict-component self-intersections.  The
stable replacement is the pair `(S,I)`: `S` contains every thickened affine
node correction, while `I` records root--affine matching.  Boundary blowups
replace strict-transform squares and raw node lengths but preserve (3.3).
Using (1.6), the combination is fixed exactly by

\[
 \boxed{S-2I=54-D_{\log}^2.}                    \tag{3.4}
\]

Thus apparent repair by increasing `N` was a bookkeeping artifact unless a
genuinely new determinant component or point module was also introduced.

## 4. The E8 cusp deficit

At the E8 cusp `P^5-Q^3=0`, the branch multiplicity is three.  The isolated
unibranch Fitting ledger gives at least `2f_i` on a complete
residue-degree-`f_i` row.  Put

\[
 F=\sum_i f_i.
\]

Any effective completion must therefore satisfy

\[
 \boxed{u-1\ge2F.}                               \tag{4.1}
\]

For simple inertia, every `e_i=2`; hence `F=R=(d-u)/2` and (4.1) becomes

\[
 u-1\ge d-u,\qquad\text{or equivalently}\qquad
 2u\ge d+1.                                      \tag{4.2}

The complete orbifold atlas violates this in every row:

\[
\begin{array}{c|c|c|c|c}
d&\#&u&R&\delta=2R-(u-1)=d-2u+1\\ \hline
6&1&2&2&3\\
10&1&2&4&7\\
12&1&4&4&5\\
15&1&3&6&10\\
20&1&4&8&13\\
24&2&4&10&17\\
30&1&4&13&23\\
30&1&2&14&27\\
40&2&4&18&33\\
60&1&4&28&53\\
120&1&4&58&113
\end{array}                                                    \tag{4.3}
\]

Consequently:

\[
 \boxed{
 \text{No simple-inertia E8 row admits an effective cyclic completion.}}
                                                               \tag{4.4}
\]

## 5. Exact survivor dichotomy

Equation (4.4) is unconditional once the complete point filtration is known
to be effective.  Without that sign theorem, it gives an equally concrete
dichotomy.  Write the complete signed point identity as

\[
 u-1=\ell_{\rm cusp}+\eta_{\rm other},\qquad
 \ell_{\rm cusp}\ge2R.                          \tag{5.1}
\]

Then every putative survivor must satisfy

\[
 \boxed{\eta_{\rm other}\le-\delta,}             \tag{5.2}
\]

with `delta` taken from (4.3).  In words, it needs a negative
normalization/`Fitt_1` class of at least the displayed magnitude.  A split
node correction, which is positive relative to the glued cyclic module,
cannot help.

The terminal and carrier common-fan packets are log-etale, the extraction
root is cyclic with exact charge `27`, the smooth endpoint is already
explicit, and the outgoing tail has zero correction.  Therefore any class
in (5.2) is forced onto the still-unresolved global affine attachment or an
uncompiled center connecting it to the determinant support.  This replaces
the former unbounded-negativity problem by one finite local sign-and-length
problem.

The present theorem does not yet prove that such a negative correction is
impossible, nor does it address meridian inertia greater than two, additional
affine nonproperness components, or `k>1`.  Those are the remaining claim
boundaries for the full `(75,125)` exclusion.

The subsequent
[`cyclic-submodule positivity theorem`](LOG_COKERNEL_CYCLIC_SUBMODULE_POSITIVITY.md)
closes this sign gap.  Every isolated `Fitt_1` defect of the actual
generically cyclic rank-two cokernel is a positive quotient of its cyclic
Cartier submodule.  Hence the negative corrections demanded by (5.2) cannot
occur, and every one-component simple-inertia E8 completion is excluded.

<!-- status-consumer: LCSP1 8658eebeb1d65671 -->

The multi-component extension and higher-inertia audit are now explicit in
the [`global affine-ramification budget`](F2_AFFINE_GLOBAL_RAMIFICATION_BUDGET.md).
For components with normalization invariants `(g_j,s_j)` it replaces (2.4)
by `d-1+sum_j A_j(2g_j-2+s_j)`.  At the degree-six floor there is one
cubic-inertia E8 equality row: its image is the natural `A_6`, its budget is
`2=2`, and its passport equals the terminal residue passport.  An exact
local model creates a generically split contracted `diag(t,t)` packet, so
the remaining obstruction is positive-dimensional `Fitt_1`, not an
isolated negative correction.

<!-- status-consumer: PF2GRB1 aa3a0efd2e0ff277 -->

## Reproduction

```bash
.venv/bin/python scripts/verify_f2_affine_k1_complete_chain_budget.py
```

The checker evaluates both global models for all carrier contacts
`0<=b<=8`, several smooth-blowup counts, and every simple-inertia orbifold
row.  It verifies the cancellation (2.5) and every deficit in (4.3).
