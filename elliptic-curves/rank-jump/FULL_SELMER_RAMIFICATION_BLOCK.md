# A full-Selmer ramification block and its remaining solubility gap

The six fixed-cubic controls have a common exact incidence constraint:
**every 2-Selmer class ramifies at all newly bad primes together, or at none.**
This applies to the full Selmer group, including classes outside the inherited
20-dimensional space and the one affine extension previously examined.
It does not bound the unramified part of that group.

At \(u=2\), twelve additional exact Cassels–Tate pairings identify the two
affine classes in the radical of the restricted pairing. Their explicit
2-covers remain of UNKNOWN rational solubility after one frozen small point
test. These results concern the fixed-cubic retrospective controls, not a
proved explanation of the prospective R17/MW17 or A1/MW16 jumps.

## Setup and the full-Selmer statement

Use \(f(X)=X^3+AX+B\), \(K=\mathbb Q(\theta)\), \(f(\theta)=0\), and the
20 independent anchor classes \(W=\langle\beta_1,\ldots,\beta_{20}\rangle\)
from [the affine analysis](AFFINE_SELMER_AND_CT.md).
The fixed-cubic pencil has distinguished root
\(\alpha_u=\theta+u\theta^2\) and discriminant
\(\operatorname{disc}(f)D(u)^2\), where \(D(u)=1+Au^2+Bu^3\).
Write
\[
\eta_u=D(u)(1-u\theta),\qquad
\kappa_u=1+u\theta+u^2(A+\theta^2).
\]
Since \(\eta_u=\kappa_u(1-u\theta)^2\), these represent the same squareclass.

Let \(S_{\rm new}\) be the newly bad primes in the retained complete local
support. In these six controls every such prime is odd, is good for the
anchor curve, satisfies \(v_p(D)=1\), and has \(u\) a unit. The distinguished
linear prime of \(K\) above \(p\) is \(\theta=1/u\bmod p\).
Let \(\ell_p:W\to\mathbb F_2\) be the previously certified root character.
All classes in \(W\) have even valuations at the primes above \(p\).

The complete local Kummer point image at \(p\) has only two possible valuation
parity vectors: zero and the nonzero vector of \(\eta_u\). The latter has
zero at the distinguished prime and one at each other prime above \(p\).
Thus every \(\delta\in\operatorname{Sel}_2(E_u/\mathbb Q)\) has a well-defined
ramification bit \(e_p(\delta)\). The new certificate verifies this assertion
on every generator of each complete local point image, at all 19 new primes.

**Theorem (these six controls).** For \(u=-3,-2,-1,1,2,3\),
\[
e(\operatorname{Sel}_2(E_u/\mathbb Q))
 \subseteq \langle(1,\ldots,1)\rangle
 \subseteq\mathbb F_2^{S_{\rm new}}.
\]
The image has dimension exactly one at \(u=-1,2,3\). At \(u=-3,-2,1\)
its dimension remains UNKNOWN, either zero or one.

### Proof

Let \(V_u\subset W\) be the classes locally in the point image of \(E_u\)
away from \(S_{\rm new}\), including the real place. Its basis is the kernel
of the retained old-place constraints. Outside the retained support both
curves have good odd reduction and the inherited classes are unramified,
so no additional constraint on \(V_u\) is required.

For \(\beta\in V_u\) and \(\delta\in\operatorname{Sel}_2(E_u/\mathbb Q)\),
the local Tate pairing vanishes at every place outside \(S_{\rm new}\):
both classes belong to the local Kummer point image, which is isotropic.
Here the pairing is the sum of the Hilbert-symbol invariants over the
factors of \(K\otimes\mathbb Q_p\).
Indeed, \(E[2]\) identifies with the even-sum submodule of \(\mu_2^3\)
indexed by the three nonzero 2-torsion points; coordinate dot product
restricts to the Weil pairing on this submodule. Cup product and
corestriction therefore give precisely this cubic-algebra formula.
Local Kummer isotropy is standard local duality; see
[Morgan, *2-Selmer parity for hyperelliptic curves in quadratic extensions*](https://londmathsoc.onlinelibrary.wiley.com/doi/10.1112/plms.12565).

At \(p\in S_{\rm new}\), \(\beta\) has even valuations. The odd-prime
Hilbert-symbol formula reduces its pairing with a local point class to
the dot product of the unit-square bits of \(\beta\) and valuation bits
of that class. The square norm of \(\beta\) makes the sum of the unit
bits at the other primes equal to the distinguished root character.
Consequently this pairing, encoded in \(\mathbb F_2\), is
\(e_p(\delta)\ell_p(\beta)\). The certificate independently checks this
identity for all 20 anchor classes against all retained local point-image
generators, including factors of residue degree two.

Global Hilbert reciprocity over \(K\) now gives
\[
\sum_{p\in S_{\rm new}} e_p(\delta)\ell_p(\beta)=0
\quad\text{for every }\beta\in V_u.
\]
Form the matrix \(R_u\) whose columns are these root characters restricted
to a basis of \(V_u\). Exact binary elimination gives the following data.

| \(u\) | New primes | \(\dim V_u\) | \(\operatorname{rank}R_u\) | \(\ker R_u\) |
|---:|---:|---:|---:|---|
| -3 | 3 | 19 | 2 | all-zero or all-one |
| -2 | 5 | 17 | 4 | all-zero or all-one |
| -1 | 2 | 19 | 1 | all-zero or all-one |
| 1 | 3 | 15 | 2 | all-zero or all-one |
| 2 | 4 | 16 | 3 | all-zero or all-one |
| 3 | 2 | 16 | 1 | all-zero or all-one |

This proves the containment for arbitrary Selmer classes, without assuming
they lie in \(W+\langle\eta_u\rangle\). At \(u=-1,2,3\), the previously
certified locally soluble classes \(\eta_u\beta_m\), with respective masks
\(m=0,591872,659456\), realize the all-one vector. At the other three
parameters the failure of the particular coset \(\eta_u+W\) to meet the
Selmer group does not exclude a realizing class outside that coset. \(\square\)

This is an **incidence** block. If \(S_u^0=\ker e\), the result gives
\(\dim\operatorname{Sel}_2(E_u)-\dim S_u^0\leq1\), with equality in the
three realized cases. It does not prove \(S_u^0=W_u\), where
\(W_u=W\cap\operatorname{Sel}_2(E_u)\). In particular, counting newly bad
primes as independent new Selmer directions is invalid here.

## Completing the affine CT column at \(u=2\)

Previously \(W_2\) had dimension 13, restricted CT rank 12, and radical
generated by anchor mask 513585. A new affine Selmer class was
\(\zeta=\eta_2\beta_{591872}\), and its pairing with that old radical was zero.
The [frozen protocol](U2_AFFINE_RADICAL_PROTOCOL.json) requests its pairings
with twelve old basis vectors, omitting the highest supported coordinate
of the known radical. Bilinearity reconstructs the last entry.

In the ordered basis retained in the certificate, the completed column is
\[
(0,0,1,1,1,0,0,0,1,1,0,0,1).
\]
Solving the old CT matrix equation gives correction coordinate mask 8049.
The enlarged 14-dimensional subspace has CT rank 12 and radical dimension
two. Its two affine radical classes are exactly
\[
\eta_2\beta_{438453},\qquad \eta_2\beta_{91780}.
\]
Their difference is \(\beta_{513585}\). The masks encode products of the
fixed 20 anchor classes, not point labels or rational rank counts.

The exact quartics, cover maps, all twelve CT witnesses and their local
Hilbert symbols are retained. For each quartic the sole point test was
PARI `hyperellratpoints` at parameter height 10000, with a ten-second
worker cap, plus an exact rational-point-at-infinity test. Both completed
with no point. Each status is **UNKNOWN_BOUNDED_MISS**. This is a
**visibility** outcome in the selected reduced quartic coordinates.
Membership in the restricted CT radical is only survival of a
**solubility obstruction**; these classes could still pair nontrivially
with missing Selmer classes, or fail a higher descent.

## Mechanisms and the next falsifiable step

1. **Incidence: reciprocity couples ramification.** This is now proved for
   the full Selmer group in six controls. It supplies a genuine arithmetic
   block before rational points on the specialized curve are supplied.
   Computing its matrix still uses public anchor classes retrospectively;
   this is not an authorized prospective selector for Agent 1.
2. **Solubility: the genus-two transport construction remains the positive
   model.** The [earlier exact identities](LINEAR_TWIST_SOLUBLE_BLOCKS.md)
   turn simultaneous square conditions into independent rational
   differences. The present theorem shows that their common affine
   ramification shift has only one degree of freedom. It does not make
   those square conditions soluble or explain the actual MW17/MW16 jumps.
3. **Weak explanation: many discriminant factors supply many directions.**
   The new ramification quotient has dimension at most one here despite
   two to five new primes. Variation in the unramified Selmer part and in
   global solubility is still essential.
4. **Unresolved: restricted CT survival implies rational points.** The two
   explicit \(u=2\) classes pinpoint this missing implication; a larger
   search box would mainly retest visibility.

The next useful experiment is to test completeness of the unramified
part at this single \(u=2\) control, using certified cubic-field
\(S\)-unit/ideal-class data at the retained old support. The falsifiable
target is \(S_2^0=W_2\). If true, the 14-dimensional extended Selmer space
is full, and the two recorded affine radical covers become a precise
higher-descent problem. If false, a missing unramified class provides a
new CT test against both covers. This requires a separately frozen
resource limit before expensive class-group work; no such computation
was launched in this experiment. Neither outcome alone proves rational
solubility.

For Agent 1, the eventual usable information is a certified ramification
quotient and the obstruction rank on a complete candidate Selmer space,
with its **incidence** and **solubility** roles kept separate. No search
policy change is justified by the current fixed-field control result.

## Certificates and replay

- [Full-Selmer ramification certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_ramification_block_v1.json).
- [CT pairing inputs](../../artifacts/generated-results/elliptic-curves/rank_jump_u2_affine_ct_inputs_v1.json).
- [Completed column and radical](../../artifacts/generated-results/elliptic-curves/rank_jump_u2_affine_radical_v1.json).
- [Two exact covers and bounded-test records](../../artifacts/generated-results/elliptic-curves/rank_jump_u2_affine_covers_v1.json).

From the repository root:

```sh
python3 elliptic-curves/rank-jump/ramification_block.py check
python3 elliptic-curves/rank-jump/u2_affine_radical.py check
sage -python elliptic-curves/rank-jump/u2_affine_radical.py verify
```

The last command verifies retained equations, points if present, CT witnesses
and Hilbert symbols; it does not repeat the bounded point searches.
The underlying complete local images are retained and replayable through
`affine_selmer.py verify`. The two frozen protocols and their source hashes
are included in the new certificates. All arithmetic stores are isolated
under `artifacts/local/rank-jump-u2-affine-radical-v1`. No live-search file,
status entry, candidate population or prospective search was changed.
