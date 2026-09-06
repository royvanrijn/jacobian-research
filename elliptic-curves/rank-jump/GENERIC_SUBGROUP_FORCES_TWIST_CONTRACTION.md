# Generic classes already force the production twist contractions

The exact Selmer contractions 8,6,3 from the
[production twist comparison](PRODUCTION_TWIST_INCIDENCE_AND_SOLUBILITY.md)
do **not** need exceptional points. In all three controls, generic rational
classes already span the entire quotient of the original local conditions
by their intersection with the twist conditions. Local duality then forces
both the subgroup inclusion and the exact dimension difference.

This removes the earlier exceptional-point dependency for the **relative**
R17-high comparison. Its two **individual** Selmer dimensions still need
the retained boundary certificate or additional arithmetic. The result also
weakens a tempting explanation of the high jumps: maximal contraction under
this twist is already forced by the generic subgroup and holds for the
observed-zero control as well.

## A generic saturation criterion

Let E and E′ have a fixed identification of their two-torsion and the same
local Selmer conditions outside a finite set S. Write

\[
 L=\bigoplus_{v\in S}\delta_vE(\mathbf Q_v),\qquad
 L'=\bigoplus_{v\in S}\delta_vE'(\mathbf Q_v),\qquad C=L\cap L'.
\]

Assume these local point products are self-orthogonal for the same
nondegenerate Tate pairing. This holds for the standard labelled
quadratic-twist comparison. Let S₀ and S₁ denote the two full Selmer groups,
and let G⊂S₀ be a known global subgroup; in this application G is the
marked generic rational Kummer subgroup. Put

\[
 d=\dim(L/C),\qquad
 e=\dim\bigl((\operatorname{loc}_S G+C)/C\bigr).
\]

**Criterion.** If e=d, then

\[
 \boxed{S_1\subseteq S_0,\qquad
 0\longrightarrow S_1\longrightarrow S_0
 \stackrel{\operatorname{loc}\bmod C}{\longrightarrow}L/C
 \longrightarrow0.}
\]

Furthermore,

\[
 \boxed{S_0=S_1+G,\qquad
 S_0/G\simeq S_1/(G\cap S_1).}
\]

These conclusions require no full relaxed boundary, no class-group
calculation, no derivative class, and no exceptional rational points.

**Proof.** For a∈S₁ and g∈G, global reciprocity gives
⟨loc a,loc g⟩_S=0. Outside S their local classes lie in the same
self-orthogonal point condition, so all omitted local contributions vanish.
Also loc a∈L′ annihilates C⊂L′. Surjectivity of G onto L/C means
loc G+C=L. Hence loc a annihilates L, and L=L^⊥ implies loc a∈L.
Thus a∈S₀. On S₀ the kernel of localization modulo C is exactly S₁,
and G supplies surjectivity. Subtracting a suitable element of G from
any element of S₀ puts it in S₁. The quotient identity follows.

The required self-duality and global reciprocity are the standard Selmer
structure framework described in
[Morgan–Paterson, §3.1](https://arxiv.org/pdf/2011.04374).
The criterion above is the direct linear-duality argument, applied here
to the marked generic subgroup; no statistical theorem from that paper
is used.

There is also a useful fail-closed bound when the gate is incomplete:

\[
 \boxed{2e-d\ \le\ \dim S_0-\dim S_1\ \le\ d.}
\]

Indeed, let S_C=S₀∩S₁. The first quotient S₀/S_C injects into L/C
and contains the e-dimensional generic image. The second quotient
S₁/S_C injects into L′/C and annihilates that generic image. The
pairing between L/C and L′/C is perfect, so its dimension is at most
d−e. Subtracting the two dimensions proves the lower bound; the upper
bound follows from dim(S₀/S_C)≤d. A failed saturation gate therefore
does not establish small rank or disprove maximal contraction.

## Three masked production comparisons

| Retrospectively selected control | Marked generic rank m | Joint generic local dimension | dim L | Local change d | Generic image e | dim(G∩S₁) |
|---|---:|---:|---:|---:|---:|---:|
| A1/MW16-05, 307/206, observed +9 | 16 | 15 | 16 | 8 | 8 | 8 |
| Published R17, −2300/843, observed +7 | 17 | 15 | 17 | 6 | 6 | 11 |
| Published R17, −1561/3133, observed 0 | 17 | 11 | 12 | 3 | 3 | 14 |

The final column uses the previously certified mod-two independence of
the marked generic subgroups: dim(G∩S₁)=m−d. The new contraction
criterion itself only needs the displayed local spanning certificate.
The observed-zero label remains a bounded-search observation, not an
exact-rank theorem.

For these three fibres, respectively, the whole additional Selmer quotient
is therefore

\[
 S_0/G\simeq S_1/G_1,\qquad\dim G_1=8,11,14.
\]

Every rational class on E, including one not yet found, can be corrected
by a generic rational class so that it becomes Selmer on E′ as well.
The correction remains rational on E. It need not be rational on E′:
the production CT blocks explicitly obstruct that inference.

Thus a large specialization quotient cannot be attributed to exceptional
directions uniquely admitted by these original local conditions. After
generic correction, all such directions lie in a common Selmer space.
Their simultaneous rational solubility on the original curve is still
the issue. This statement concerns the specified −1 twist comparison;
it does not rule out other arithmetic incidence events in the family.

## Bounded experiment and independence of exceptional points

The [frozen protocol](GENERIC_LOCAL_CONTRACTION_PROTOCOL.json) allowed only
three 30-second workers, no new local or global candidate enumeration,
and no parameter sweep. The export includes:

* each integral cubic and its complete discriminant factorization;
* only the marked generic rational points, transported exactly;
* the local places and previously certified local twist x witnesses.

The [masked input](../../artifacts/generated-results/elliptic-curves/rank_jump_generic_local_contraction_inputs_v1.json)
contains no exceptional point, derivative class, CT matrix, strict character
basis or jump label. The worker enforces its field whitelist, rechecks the
equations and discriminant factors, constructs the local cubic algebras,
and checks all local dimensions and twist point witnesses anew.
It does not read the earlier retrospective bundles in worker mode.
The controls themselves were chosen retrospectively; masking proves an
input-dependency statement, not prospective selection performance.

The [local certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_generic_local_contraction_v1.json)
stores an explicit generic correction mask for **every generator of L**.
Adding that correction puts the generator in C. Sage vector spaces
independently verify those corrections and both quotient dimensions.
All three cases pass. The
[consequence certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_generic_contraction_consequences_v1.json)
records the resulting full-Selmer identities.

## Why a saturated generic image is not itself a rare-jump signal

There is a simple dimensional caution. In an artificial model where the
g-dimensional generic image is uniform among subspaces of an
ℓ-dimensional binary L, its probability of surjecting onto L/C of
dimension d is

\[
 \frac{2^{d(\ell-g)}{\ell-d\brack g-d}_2}{{\ell\brack g}_2}.
\]

To count, choose its (g−d)-dimensional intersection with C and then
choose the graph of a map from L/C into the remaining (ℓ−g)-dimensional
quotient of C. For the table's dimensions, these fractions are
256/257, 2859118592/2863245995, and 584/585—each above 99.6%.

This is **not an arithmetic probability estimate**. Generic local images
need not be uniformly distributed, and global constraints matter. The
calculation only shows why saturation can be commonplace when the marked
generic image already occupies nearly all of L. Three retrospective
successes are insufficient to turn it into a high-rank selector.

## Revised priorities

1. **Supported solubility mechanism:** a shared arithmetic class space can
   contain a large rational block on one curve and an obstructed block on
   its twist. The independent cup/CT certificates remain evidence for this.
2. **Useful incidence certificate, weak high-jump explanation:** generic
   saturation predicts the exact relative Selmer contraction before
   exceptional points are supplied. The contraction is already accounted
   for by the generic subgroup, including in the observed-zero control.
3. **Missing implication:** an independently computable specialization
   condition must create enough global classes in S₁/G₁ and make many of
   them rational on the original curve. Local saturation supplies neither.
   Full class-group information and higher solubility obstructions remain
   distinct computations.
4. **Potential use for Agent 1:** a generic local spanning certificate can
   replace an expensive full-boundary calculation when only a relative
   Selmer comparison is needed. It is an incidence diagnostic. It does
   not justify promoting a candidate or changing point-search budgets.

The next useful mechanism must distinguish rational images inside the
common Selmer quotient, or independently establish that quotient's size;
another successful local contraction alone would add little evidence.

Replay from the repository root:

```sh
sage -python elliptic-curves/rank-jump/generic_local_contraction.py check
python3 elliptic-curves/rank-jump/generic_contraction_consequences.py check
python3 -m unittest discover -s elliptic-curves/rank-jump -p test_generic_local_contraction.py
```
