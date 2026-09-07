# Exact full-Selmer comparison without a class group

The six fixed-cubic controls now have **exact relative full 2-Selmer
dimensions**, despite the anchor's absolute Selmer dimension remaining
UNKNOWN. If \(s_0=\dim\operatorname{Sel}_2(E_0/\mathbb Q)\), then
\[
\begin{array}{c|rrrrrr}
u&-3&-2&-1&1&2&3\\ \hline
\dim\operatorname{Sel}_2(E_u/\mathbb Q)-s_0&-3&-7&-1&-7&-6&-4.
\end{array}
\]
The uncomputed excess above each retained Selmer subspace is the same
unknown quotient of the anchor Selmer group. This resolves the preceding
ramification uncertainty: its full image has dimension zero at \(u=-3,-2,1\)
and one at \(u=-1,2,3\).

This is an **incidence** theorem. It does not turn Selmer classes into
rational points, give absolute Selmer dimensions, or prove a numerical rank
upper bound. It concerns the six fixed-field controls and their rank-at-least-20
anchor; an explanation of the actual prospective MW17/MW16 jumps remains open.

## Why a class-group computation was unnecessary for this comparison

The [previous endpoint](FULL_SELMER_RAMIFICATION_BLOCK.md) proposed testing
the unramified Selmer part by a cubic-field \(S\)-unit/class-group computation.
Inspection of the existing complete local images revealed a stronger
alternative. The twenty global anchor classes span the full anchor local
point image at every retained place. The explicit transporter \(\eta_u\)
lies locally in the sum of the anchor and specialized point images.
Together these 21 global classes fill half of the relevant local quotient.
Reciprocity then proves that they fill its entire possible global image.

The [frozen protocol](SELMER_COMPARISON_PROTOCOL.json) verifies this observation
on all six existing controls using binary linear algebra only. No class-group
computation, parameter sweep, local witness search or point search was run.
The result does not bypass the absolute Selmer problem; it proves that only
one common absolute unknown remains.

## A general boundary-completeness lemma

Let two elliptic curves \(E_0,E_1\) over \(\mathbb Q\) have a fixed
Galois-equivariant identification of their 2-torsion. This identification
preserves the Weil pairing: every invertible map on a two-dimensional
\(\mathbb F_2\)-space is symplectic. Write \(H_v=H^1(\mathbb Q_v,E_0[2])\)
and let \(L_{0,v},L_{1,v}\subset H_v\) be their local point images.
Both are their own orthogonal complements for local Tate duality.
Set
\[
C_v=L_{0,v}\cap L_{1,v},\qquad
D_v=L_{0,v}+L_{1,v},\qquad
\mathcal B=\bigoplus_vD_v/C_v.
\]
Only finitely many summands are nonzero. Local duality makes
\(\mathcal B\) nondegenerate, since \(D_v^\perp=C_v\), and
\[
\dim\mathcal B=2d,\qquad
d=\sum_v(\dim L_{0,v}-\dim C_v).
\]

Let \(S_C,S_D\) denote the global cohomology classes satisfying the local
conditions \(C_v,D_v\), respectively. Localization gives an injective map
\[
S_D/S_C\longrightarrow\mathcal B.
\]
Its image is isotropic: two global classes have total local cup-product
invariant zero by reciprocity, and outside the changed places their
localizations lie in the same isotropic point image.

**Lemma.** If a certified global subspace \(T\subset S_D\) has image of
dimension \(d\) in \(\mathcal B\), then
\[
S_D=S_C+T.
\]

**Proof.** The image of \(T\) is an isotropic half-dimensional subspace,
hence equals its orthogonal complement. Every class in \(S_D\) is
orthogonal to it by reciprocity, so has boundary in its image. Subtract a
class in \(T\) with that boundary; the difference lies in \(S_C\). \(\square\)

The local-duality and Kummer self-duality inputs are
[Morgan, *On 2-Selmer groups of twists after quadratic extension*, Theorem 3.1 and Example 3.5](https://londmathsoc.onlinelibrary.wiley.com/doi/full/10.1112/jlms.12533).
The elementary dimension argument above suffices; no Selmer completeness
assumption or conjectural class-group bound enters.

## Applying the lemma to the fixed-cubic pencil

Use the labelled identification given by
\(\alpha_u=\theta+u\theta^2\), with \(K=\mathbb Q(\theta)\), and the
twenty-dimensional space \(W\) of independent anchor point classes.
Take
\[
T_u=W+\langle\eta_u\rangle,\qquad
\eta_u=[D(u)(1-u\theta)].
\]
These are the same classes and local signatures as in
[the affine Selmer calculation](AFFINE_SELMER_AND_CT.md).
Each \(T_u\) has dimension 21: a newly bad prime gives \(\eta_u\) an odd
valuation where all classes in \(W\) have even valuations.

At every retained place, the span of the anchor signatures has the
independently known full dimension of \(E_0(\mathbb Q_v)/2E_0(\mathbb Q_v)\).
Its containment in that point image follows from the exact global anchor
points. The local dimension depends on the local field and the 2-torsion
module, so agrees with the dimension independently certified for \(E_u\).
Thus the span is exactly \(L_{0,v}\), including at 2 and infinity.

The checker forms \(C_v,D_v\), verifies every signature in \(T_u\) belongs
to \(D_v\), and reduces all 21 signatures modulo \(C_v\). Concatenating
these remainders over the places gives their exact boundary matrix.
Outside the retained support both curves have good odd reduction and
the classes are unramified, so both local conditions agree.

Write \(W_u=W\cap\operatorname{Sel}_2(E_u)\) and
\(T_u^{\rm Sel}=T_u\cap\operatorname{Sel}_2(E_u)\). The results are:

| \(u\) | \(d\) | \(\dim W_u\) | Anchor boundary rank | Specialized boundary rank \(b_u\) | \(\dim T_u^{\rm Sel}\) |
|---:|---:|---:|---:|---:|---:|
| -3 | 4 | 17 | 3 | 0 | 17 |
| -2 | 8 | 13 | 7 | 0 | 13 |
| -1 | 3 | 18 | 2 | 1 | 19 |
| 1 | 8 | 13 | 7 | 0 | 13 |
| 2 | 8 | 13 | 7 | 1 | 14 |
| 3 | 6 | 15 | 5 | 1 | 16 |

In each case the full 21-class boundary has rank \(d\), so the lemma
applies. Its intersection with the anchor local conditions is precisely
the boundary of \(W\); the additional class is ramified at a prime good
for \(E_0\). Its intersection with the specialized local conditions is
precisely the boundary of \(T_u^{\rm Sel}\), computed by the prior affine
test. Consequently, for \(S_0=\operatorname{Sel}_2(E_0)\) and
\(S_u=\operatorname{Sel}_2(E_u)\),
\[
S_0=W+S_C,\qquad S_u=T_u^{\rm Sel}+S_C,\qquad
W\cap S_C=T_u^{\rm Sel}\cap S_C=W_u.
\]
In particular the inclusion maps give canonical isomorphisms
\[
\boxed{S_0/W\ \simeq\ S_C/W_u\ \simeq\ S_u/T_u^{\rm Sel}.}
\]
Thus all six uncomputed excesses are the same anchor quotient. Taking
dimensions proves the table of full-Selmer differences at the start.
This is stronger than an equality of observed subspace dimensions.

## Exact ramification and the unresolved \(u=2\) covers

The full new-prime ramification quotient is now completely determined.
Every \(S_C\) class lies in \(S_0\) and is unramified at every newly bad
prime. Within \(T_u^{\rm Sel}\), ramification is zero on \(W_u\) and
all-one on its affine coset when that coset exists. Hence
\[
\dim e(S_u)=b_u,\qquad
\ker(e:S_u\to\mathbb F_2^{S_{\rm new}})=S_C=S_0\cap S_u.
\]
The three zero conclusions improve the UNKNOWN cases in
[the earlier ramification certificate](FULL_SELMER_RAMIFICATION_BLOCK.md).
That certificate remains valid as a historical weaker bound.

In particular the proposed \(u=2\) completeness test has an exact answer
in terms of the anchor:
\[
S_2^0=W_2\quad\Longleftrightarrow\quad S_0=W.
\]
If an unramified Selmer class is missing at \(u=2\), an equally missing
anchor Selmer class exists, with the quotient correspondence above.
The two retained affine CT-radical covers of masks 438453 and 91780
remain of UNKNOWN rational solubility. Completing the global boundary
does not complete the global Selmer group or its CT pairing.

## What this says about the high-jump event

Put \(\epsilon=\dim S_0-20\geq0\), which remains UNKNOWN. The complete
accounting is now
\[
\begin{array}{c|rrrrrr}
u&-3&-2&-1&1&2&3\\ \hline
\dim S_u&17+\epsilon&13+\epsilon&19+\epsilon&
13+\epsilon&14+\epsilon&16+\epsilon\\
\text{certified restricted CT rank}&16&12&16&12&12&14\\
\text{resulting MW-rank bound}&1+\epsilon&1+\epsilon&3+\epsilon&
1+\epsilon&2+\epsilon&2+\epsilon.
\end{array}
\]
The last row is an inequality in an **unknown** \(\epsilon\), not a
numerical rank upper bound. It follows because the CT pairing kills
rational point classes, and its rank on a subspace is a lower bound
for the full pairing rank. The irreducible cubic excludes rational
2-torsion. No finiteness assumption on the full Sha group is needed.

This separates the two arithmetic effects in an actual high/low control:
the anchor has twenty rational independent directions, while \(u=-1\)
has full Selmer dimension only one smaller, yet carries a certified
16-dimensional nondegenerate CT block. A near-preservation of Selmer
incidence therefore coexists with a large loss of rational solubility
inside the transported space. At \(u=2\) both effects occur: the full
Selmer dimension drops by six and the retained pairing has rank twelve.
These are exact comparisons; a bounded point-search miss is not used.

The strongest mechanism template remains
\[
\text{specialization geometry}
\ \Rightarrow\ \text{simultaneously soluble classes}
\ \Rightarrow\ \text{independent rational directions}.
\]
Here global duality closes the **incidence comparison** completely relative
to one anchor unknown. The missing implication is still the simultaneous
rational solubility of the surviving arithmetic classes. The existing
genus-two transport identities provide a constructive model of that
implication, but not a general explanation for the prospective high fibres.

The next useful target is to understand the CT change on the common
Selmer space through the explicit transporter/cover geometry. A single
anchor completeness certificate would settle all six absolute Selmer
dimensions simultaneously; repeating six class-group calculations would
not address six independent unknowns. Any prospective use by Agent 1
would need an equation-derived replacement for the public anchor
classes, and a separate solubility test. No selector change is justified.

## Reproducibility

The [certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_selmer_comparison_v1.json)
retains each local intersection and sum, all boundary rows, the strict
kernel, both point-structure kernels, and source hashes.

```sh
python3 elliptic-curves/rank-jump/selmer_comparison.py check
sage -python elliptic-curves/rank-jump/affine_selmer.py verify
sage -python -m unittest discover -s elliptic-curves/rank-jump -p test_selmer_comparison.py
```

The first command replays the new binary certificate. The second
independently checks the retained local arithmetic used as input.
The third also reconstructs every local quotient and global strict kernel
with Sage's independent finite-field linear algebra.
No active search file, status entry or Agent 1 output was changed.
