# Strict Selmer incidence and explicit Artin-dual ideal blocks

The unramified character blocks now have an exact **point-blind ambient
description**:
\[
\boxed{\operatorname{Sel}_2^S(E)
 \simeq \operatorname{Hom}(\operatorname{Cl}(\mathcal O_{K,S_K}),\mathbf F_2).}
\]
Here \(S\) contains 2, infinity and every bad prime of \(E/\mathbb Q\);
the strict Selmer group means the classes whose localization is **zero**
at every place in \(S\). This is stronger than local solubility there.

The existing arithmetic reduces the unknown full 2-Selmer dimensions on
three actual fibres to
\[
25+\epsilon+b,\quad24+\epsilon+b,\quad17+\epsilon+b,
\qquad \epsilon\ge0,\quad b\in\{0,1\},
\]
respectively. Each curve has its own \(\epsilon,b\). The excess \(\epsilon\)
is an uncomputed \(S\)-class-group dimension; the remaining boundary
uncertainty is at most one dimension. No numerical full-Selmer or rank
upper bound follows until \(\epsilon\) is bounded.

A separate frozen experiment also produces concrete independent small
ideal classes dual to the eight R17 high and six R17 control characters.
Its fixed prime bound fails to separate one of the ten MW16-05 characters.
That negative endpoint is retained without enlarging the dictionary.

## The strict Selmer identification

Let \(E/\mathbb Q\) have irreducible two-division cubic, and let
\(K=\mathbb Q(\theta)\) be its labelled cubic field. Use the usual Kummer
identification
\[
H^1(\mathbb Q,E[2])=
\ker\bigl(N:K^\times/K^{\times2}\longrightarrow
                     \mathbb Q^\times/\mathbb Q^{\times2}\bigr).
\]
This comes from the norm sequence for the permutation module on the three
nonzero two-torsion points. Its diagonal section exists because the degree
is odd. For the retained short models, a point contributes \(x-\theta\).

Let \(S_K\) be all finite primes of \(K\) above the finite part of \(S\), and
\[
\mathrm{Cl}_S(K)=\mathrm{Cl}(\mathcal O_K)/
                  \langle[\mathfrak p]:\mathfrak p\in S_K\rangle.
\]
This is the ordinary \(S\)-ideal class group; real places are required to
remain real in its class field. By global class field theory, its quadratic
characters correspond exactly to quadratic extensions of \(K\) unramified
everywhere, including infinity, and split at every prime in \(S_K\).
See [Milne, Chapter V, Example 3.9 and Exercise 3.15](https://www.jmilne.org/math/CourseNotes/CFT.pdf).

For such an extension \(K(\sqrt\beta)/K\), every finite valuation of
\(\beta\) is even and \(\beta\) is positive at each real embedding.
Consequently \(N\beta\) is a positive rational number with even valuation
at every prime, hence a rational square. Thus **no further norm-square
filter is needed** for these ordinary unramified characters over
\(\mathbb Q\).

Above \(S\), the associated \(E[2]\)-class is zero. At any remaining
prime \(p\), the curve has good reduction and \(p\ne2\). Its local
Kummer image is exactly \(H^1_{\rm ur}(\mathbb Q_p,E[2])\): reduction
identifies the point quotient with the finite-field quotient, the formal
kernel is 2-divisible, and the finite-field Kummer sequence gives this
unramified cohomology group. Hence every \(S\)-split unramified character
is an elliptic Selmer class.

Conversely, a strict Selmer class is locally square above \(S\), and
its local classes outside \(S\) are unramified by good reduction.
It therefore gives exactly such an unramified, \(S_K\)-split quadratic
extension. This proves the displayed isomorphism in both directions.

The relation between cubic class groups and 2-Selmer groups is classical;
see also [Barrera Salazar–Pacetti–Tornaría, Section 2](https://www.famaf.unc.edu.ar/~apacetti/papers/2selmer.pdf).
The statement here deliberately imposes zero local classes at **all**
bad places and at 2 and infinity. It does not apply an unrestricted
class-group inclusion while dropping that paper's bad-place hypotheses.

The [previous character-block theorem](BAD_PLACE_CLASS_FIELD_BLOCKS.md)
gave a witnessed subspace of this ambient group. The present identification
closes the **incidence** implication for every independently constructed
\(S\)-class character, but not its rational-solubility implication.

## Separate the two remaining dimensions of uncertainty

Let \(T=\mathrm{Sel}_2(E)\), \(U=\mathrm{Sel}_2^S(E)\), and let
\[
L_S=\prod_{v\in S}\delta_v E(\mathbb Q_v),\qquad
I=\operatorname{loc}_S(T)\subset L_S.
\]
Write \(c_S=\dim\mathrm{Cl}_S(K)/2\). There is an exact sequence
\[
0\longrightarrow U\longrightarrow T\longrightarrow I\longrightarrow0,
\qquad \dim T=c_S+\dim I.
\]

For the generic Kummer subspace \(G\) of dimension \(m\), let
\(g=\dim\operatorname{loc}_S(G)\). Quotienting gives
\[
0\longrightarrow U/(U\cap G)\longrightarrow T/G
 \longrightarrow I/\operatorname{loc}_S(G)\longrightarrow0,
\]
so the exact **point-blind strict incidence contribution** is
\[
c_S-(m-g).
\]
Both \(c_S\) and the marked generic contribution are defined before
exceptional points are supplied. The current lower bounds for \(c_S\)
still use those points and must not be relabelled point-blind measurements.

Now include the independent known point space \(W\supset G\), of dimension
\(n\). Let \(k=\dim(W\cap U)\), \(a=\dim\operatorname{loc}_S(W)\), and
\(\ell=\dim L_S\). Then \(n=k+a\), and
\[
\dim T=n+\underbrace{(c_S-k)}_{\epsilon\ge0}
       +\underbrace{(\dim I-a)}_{0\le b\le\ell-a}.
\]
The [numeric certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_strict_selmer_model_v1.json)
specializes this exact sequence using the complete local calculations:

| Fibre | \(m\) | \(n\) | \(k\) | \(a\) | \(\ell\) | Full \(\dim T/G\) |
|---|---:|---:|---:|---:|---:|---|
| MW16-05 \(307/206\) | 16 | 25 | 10 | 15 | 16 | \(9+\epsilon+b\) |
| R17 \(-2300/843\) | 17 | 24 | 8 | 16 | 17 | \(7+\epsilon+b\) |
| R17 \(-1561/3133\) | 17 | 17 | 6 | 11 | 12 | \(\epsilon+b\) |

In all three cases \(\ell-a=1\). Thus, after the retained evidence,
the full Selmer incidence problem has only an unknown \(S\)-class excess
and one possible further local-image dimension. We have not computed
either unknown, identified the remaining local line, or proved that any
new class would be rational. If an independent class-group upper bound
were to prove \(c_S=k\), it would immediately give
\(\dim T\le n+1\); this is a conditional implication, not an available
upper bound.

For rational solubility, replace \(U\) by
\(U\cap\delta E(\mathbb Q)\). An \(S\)-class character is rational exactly
when its image in \(\Sha(E)[2]\) is zero. The new identification does not
compute that kernel. This is still the missing second implication in
a rank-jump mechanism.

## Frozen small-ideal experiment

The [protocol](STRICT_ARTIN_PROTOCOL.json) fixes the three completely
covered curves and all degree-one prime ideals
\[
\mathfrak p_{p,a}=(p,\theta-a)
\]
with odd \(p\le97\), \(p\) not dividing the integral cubic discriminant,
and \(a\) a root of that cubic modulo \(p\). Ideals are ordered by \(p,a\).
This dictionary uses only the curve's cubic, with no point coordinates
or character values.

After freezing it, the experiment evaluates the **retained** strict
characters on those ideals. It uses generic strict characters first,
then the relative character basis from the previous certificate.
The resulting Artin matrix is therefore a retrospective diagnostic.

For an unramified class \(\beta\), its Artin bit at \(\mathfrak p\)
is zero exactly when \(\beta\) is a square in \(K_{\mathfrak p}\).
An invertible evaluation submatrix certifies independence of the
selected ideal classes in \(\mathrm{Cl}_S(K)/2\). This remains valid even
when the ideals are not known to generate the entire class group.

| Fibre | Ideals in fixed pool | Retained characters | Artin rank | Complete dual basis? |
|---|---:|---:|---:|---|
| MW16-05 high | 11 | 10 | 9 | no |
| R17 high | 16 | 8 | 8 | yes |
| R17 low | 18 | 6 | 6 | yes |

The [Artin certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_strict_artin_v1.json)
retains the full matrices, exact prime/root labels, selected ideals and
ideal-product words dual to every character when the rank is full.

For R17 high, eight independent ideal classes can be chosen with labels
\[
(23,16),(29,18),(31,22),(41,16),
(43,4),(43,40),(47,35),(59,50).
\]
Their dual products give two coordinates for the generic character block
and six for the relative block. For example
\[
\mathfrak p_{31,22}\mathfrak p_{47,35}
\]
pairs to one with the eighth retained character and to zero with the
other seven, including both generic strict characters. This is a concrete
ideal-class witness, not a claim about the half ideal of that character.

For MW16-05, the character-coordinate mask **74**, equivalently point-mask
**2441006** in the retained 25-point basis, is a nonzero unramified class
that evaluates to zero on every ideal in the fixed pool. Its global
nontriviality was already certified. The failure to detect it is not
evidence for dimension nine rather than ten. The protocol's prime bound
was not increased.

## Artin pairing is not Cassels–Tate pairing

The Artin matrix pairs **class-group characters with ideal classes**.
A nonzero entry detects an ideal's nontrivial image in a class-group
quotient. The Cassels–Tate pairing instead detects obstructions to
rational solubility of Selmer classes.

Every character evaluated here came from rational point witnesses, so
its elliptic CT row vanishes. Its Artin row can nevertheless be nonzero
and independent. Full Artin rank therefore supplies concrete incidence
coordinates; it is not evidence for nonzero Sha or a computation of CT.

## Consequences and next mathematical gate

1. **Incidence:** the exact strict Selmer object is now an ordinary
   \(S\)-class-group dual. Independently finding one of its characters
   really does construct a Selmer class; no additional norm or local
   solubility test remains for that strictly split class.
2. **Incidence:** small ideals give concrete independent class coordinates
   on the R17 pair, but a fixed small dictionary can miss a proven class.
   Neither a deficient Artin matrix nor a class-count lower bound is an
   upper bound on the full ambient space.
3. **Solubility:** the unknown kernel of
   \(\mathrm{Hom}(\mathrm{Cl}_S(K),\mathbf F_2)\to\Sha(E)[2]\)
   is the remaining target. This is where a common auxiliary construction
   or simultaneous obstruction cancellation would have to act.
4. **Weak explanations:** raw class-group size without subtracting generic
   characters, or treating Artin nonvanishing as a CT obstruction.

The next useful computation is an independently bounded source or upper
bound for \(c_S\) on this completely covered R17 pair. The displayed small
ideals are a fixed candidate generating set for such a class-relation
calculation, but no completeness or principal-relation certificate exists
yet. Do not launch a large class-group calculation or expand this finite
Artin dictionary on the strength of this note alone.

For Agent 1, the candidate **incidence** feature is
\(c_S-\dim G^0_S\), provided \(c_S\) is computed independently of exceptional
points. Its relation to rational rank still requires a **solubility**
feature. The present retrospective values are not ready for a selector.

## Replay and scope

From the repository root:
```sh
python3 elliptic-curves/rank-jump/strict_artin.py check
python3 elliptic-curves/rank-jump/strict_selmer_model.py check
sage -python elliptic-curves/rank-jump/verify_strict_artin.py --index 0
sage -python elliptic-curves/rank-jump/verify_strict_artin.py --index 4
sage -python elliptic-curves/rank-jump/verify_strict_artin.py --index 5
```
The three independent arithmetic replays verify 45 labelled prime ideals,
their norms and all 346 Frobenius bits with PARI's local-power routine.
All passed within the 40-second per-curve replay cap. The dimension
accounting and ideal-word matrices are also exactly replayed.

No new curve, point, class group, discriminant factorization or parameter
was searched for. The active search files and mathematical-status entries
are unchanged.
