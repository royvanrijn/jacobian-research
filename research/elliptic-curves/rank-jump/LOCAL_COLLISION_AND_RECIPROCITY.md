# Root collisions, reciprocity, and loss of a soluble block

The next experiment explains **all nineteen new-prime local conditions** in
the retained fixed-field comparison, including their dependencies. It also
proves that its six nonzero deformations preserve the 2-division field while
introducing a full three-dimensional block of new quadratic **4-division
field** data. These field directions are not Mordell–Weil directions.

This advances the incidence/solubility analysis beyond the first
[retrospective panel](ANALYSIS.md). It does not explain a new MW17/MW16 high
fibre by a positive construction. The strongest result about simultaneous
rational solubility here is negative: **no nonzero inherited class can be
rationally soluble on two distinct nonzero tested deformations**.

This audit follows the comparison panel and negative carrier experiment with
a precise local hypothesis, its bounded exact test, a reciprocity follow-up,
and two local/field lemmas. No candidate search, parameter expansion, class group,
new local CAS calculation or new CT pairing ran. The goal of explaining and
predicting the original large rank jumps remains open.

## Fixed input and experiment

Use the existing anchor

\[
 f(x)=x^3+Ax+B,\qquad
 A=-5750886029903523759416717668139307,\quad
 B=167347710468055045100164888198438918505621536951206.
\]

The twenty certified anchor classes are \(\beta_i=x(P_i)-\theta\), where
\(f(\theta)=0\). They span \(W\cong\mathbf F_2^{20}\). The curve family is

\[
 E_u:y^2=\operatorname{Norm}(x-\theta-u\theta^2),\qquad
 D(u)=1+Au^2+Bu^3=\operatorname{Norm}(1-u\theta).
\]

The parameters are exactly \(-3,-2,-1,0,1,2,3\). The
[local protocol](LOCAL_COLLISION_PROTOCOL.json) was frozen before comparing
its predicted root characters with the retained local-condition matrices.
All input comes from the committed full-span local certificates and the
previously computed CT matrices, pinned through commit `661246f`. The
[portable input](../../artifacts/generated-results/elliptic-curves/rank_jump_local_collision_inputs_v1.json)
contains the projected local matrices, fixed class witnesses and compatible
ordered bases. No running search output is read.

At every new odd prime with \(v_p(D(u))=1\), predict one obstruction functional
by the anchor's Kummer character at the residue root \(\theta=1/u\). Compare
its span with the independently retained exact local obstruction span. A
zero restriction to W was permitted; a single mismatch would refute the rule.
All **19/19** eligible events agree, and every event in this dataset has a
nonzero restriction. This is an incidence statement for the specified W.

## A local collision lemma

Let p be odd, let A,B be p-integral with \(p\nmid\operatorname{disc}(f)\), and
let u be a p-adic unit with \(v_p(D(u))=1\). Identify the 2-torsion modules
through the labelled roots \(\alpha_i=\theta_i+u\theta_i^2\). Let U be the
unramified norm-square Kummer space in the étale cubic algebra of f over
\(\mathbf Q_p\). Equivalently, U is the local Kummer image of the good-reduction
anchor. The good-reduction identification with unramified cohomology is
standard; see [Mazur–Rubin, proof of Theorem 6.1, Case 3](https://bpb-us-e1.wpmucdn.com/sites.harvard.edu/dist/a/189/files/2023/01/SELMER-COMPANION-CURVES.pdf).
Their §5 also describes multiplicative local Kummer images. The following
root-coordinate argument specializes this local theory to our pencil.

There is a unique simple residue root \(\theta_k\equiv1/u\pmod p\), which
lifts to a rational p-adic root. Write \(\chi_k\) for the square character of
the component of a class at that root. Then

\[
 \boxed{U\cap\delta_p(E_u(\mathbf Q_p))=\ker(\chi_k:U\to\mathbf F_2).}
\]

Consequently, on any fixed subspace of U, the local obstruction is the
restriction of this single linear functional. This is not a statement that
all surviving classes have global rational points.

**Proof.** All roots of f lie in an unramified splitting extension and their
pairwise differences are units. For distinct i,j,k,

\[
 \alpha_i-\alpha_j=(\theta_i-\theta_j)(1-u\theta_k).
\]

Exactly one factor \(1-u\theta_k\) has valuation one; the other two are units.
Thus two transformed roots a1,a2 have difference of valuation one and reduce
to a common value a, while the third root b remains distinct. The reduction
is the nodal cubic \(y^2=(x-a)^2(x-b)\).

Consider a local point with unramified Kummer class. A nonintegral point
reduces to O and gives the trivial residue character. A point reducing to
the node cannot give an unramified class: in the split case at least one of
\(x-a_1,x-a_2\) has valuation one, since their difference has valuation one;
in the nonsplit unramified quadratic factor their valuations agree and both
are one. These valuations cannot be made even by multiplying by squares.
At a colliding rational 2-torsion point the special Kummer component is
\(f'_u(a_i)\), also of valuation one, so this case is excluded as well.

For a point with smooth reduction, the b-component has square residue:
away from b it is \(x-b=(y/(x-a))^2\) on the reduced cubic. At b the special
component is \((b-a)^2\), again a square. Hence unramified classes in the
new local image satisfy \(\chi_k=0\).

For sufficiency there are two possibilities for the remaining étale
quadratic factor. If it is a field, its unramified square character and the
b-component character must agree by the norm-square condition. Requiring
the latter to vanish makes both trivial, giving the identity class. If the
quadratic factor splits, the kernel has two residue character patterns:
\((0,0,0)\) and \((1,1,0)\) on \((a_1,a_2,b)\). The latter is realized by a
smooth reduced point with

\[
 x=b+s^2,\qquad y=s(x-a),\qquad x-a=s^2-(a-b)\text{ nonsquare}.
\]

Such an s exists because for every nonzero d in \(\mathbf F_p\),
\(\sum_s\left(\frac{s^2-d}{p}\right)=-1\). The point is smooth and lifts by
Hensel's lemma. Its unramified characters are \((1,1,0)\). Thus every class
in the kernel is realized locally. QED.

For an anchor point, the character is computed with no p-adic search:
use \(x(P)-1/u\) if its reduction is nonzero; if P reduces to that 2-torsion
root, use \(f'(1/u)\); if P reduces to O, use zero. These special-value
conventions concern squareclasses, not the literal value of a function at
its zero. The replay independently evaluates these characters and checks
the full twenty-coordinate obstruction word at each prime.

## Reciprocity explains the correlated cuts

The first output showed that five new primes at u=-2 impose only four
independent conditions, while the two new primes at u=-1 impose only one.
The [follow-up protocol](RECIPROCITY_PROTOCOL.json) tests an explanation by
rational Hilbert reciprocity rather than treating the prime conditions as
independent Bernoulli events.

Set \(s=1/u\). At an eligible new prime, for each anchor point P,

\[
 (-1)^{\chi_{p,s}(P)}=(f(s),x(P)-s)_p.
\]

If \(x(P)-s\) is a unit, this follows immediately from \(v_p(f(s))=1\).
If P reduces to the root s, then \(v_p(x(P)-s)=1\) and
\(f(x(P))=y(P)^2\) implies that the two leading terms of
\(f(s)+(x(P)-s)f'(s)\) cancel. The odd-prime Hilbert formula then gives
\(\left(\frac{f'(s)}p\right)\), precisely the special Kummer value above.
A point reducing to O contributes trivially because x has even valuation
and square leading unit.

Rational Hilbert reciprocity therefore supplies a relation between the new
prime words and boundary words at inherited primes and infinity. The omitted
place check is explicit: outside 2, the support of f(s), and denominators of
s,A,B, both entries are units unless \(x(P)-s\) has nonzero valuation. A pole
has even valuation; a zero forces f(s) to be a nonzero residue square by
\(f(x(P))=y(P)^2\). In both cases the symbol is trivial. The existing support
contains all required primes, which the checker confirms by exact division.

The [reciprocity artifact](../../artifacts/generated-results/elliptic-curves/rank_jump_reciprocity_v1.json)
replays every local symbol and product on all twenty points. At u=-3 the sum
of new words is the inherited p=13 word; at u=1 it is an inherited p=31 word.
For the other four parameters their sum is zero. In all six cases the boundary
word lies in the actual inherited-prime obstruction space. **Modulo that
space, the complete relation space among the new-prime words is exactly the
one-dimensional span of the all-ones word.** There are no additional relations
in this finite dataset.

| u | Rank of inherited finite cuts | New primes | Rank of new cuts alone | Rank of real cut | Total local codimension | Locally admissible dimension |
|---:|---:|---:|---:|---:|---:|---:|
| -3 | 1 | 3 | 3 | 0 | 3 | 17 |
| -2 | 3 | 5 | 4 | 0 | 7 | 13 |
| -1 | 1 | 2 | 1 | 0 | 2 | 18 |
| 1 | 4 | 3 | 3 | 1 | 7 | 13 |
| 2 | 3 | 4 | 3 | 1 | 7 | 13 |
| 3 | 3 | 2 | 1 | 1 | 5 | 15 |

Every row satisfies the exact **observed** formula

\[
 \operatorname{codim}_W W_u
 =r_{\rm inherited}+\#\{\text{new primes}\}-1+r_\infty.
\]

The minus one is explained by reciprocity; equality and independence of the
real cut are verified for these rows, not asserted universally. Higher
valuations, different old-prime behaviour, additional linear relations and
other subspaces need separate checks. The positive-u real cuts are retained
complete local computations; the global cubic Galois group being unchanged
does not preserve the labelled real Kummer image.

This is an interpretable incidence feature for a specified transported class
space: **rank of root-character cuts modulo inherited conditions**, with
reciprocity dependencies included. Counting new primes without the relations
loses arithmetic information. It remains unavailable as a selector on the
production MW17/MW16 families, where an analogous point-blind exceptional
class space has not been constructed.

## The preserved cubic hides a full change in four-division data

Let \(F_0=\mathbf Q(E_0[4])\), and suppose f is irreducible. If u has even
one odd prime of good anchor reduction with \(v_p(D(u))=1\) and u a p-adic
unit, then

\[
 \boxed{[F_0\,\mathbf Q(E_u[4]):F_0]=8,\qquad
 \operatorname{Gal}(F_0\,\mathbf Q(E_u[4])/F_0)\cong(\mathbf Z/2)^3.}
\]

**Proof.** The standard root-difference description of the 4-division field is
\(\mathbf Q(E[4])=\mathbf Q(E[2],i,\sqrt{e_1-e_2},\sqrt{e_1-e_3},\sqrt{e_2-e_3})\);
see [Yelton, Theorem 1(a)](https://arxiv.org/pdf/1704.06190). Applying the
root-difference identity above shows

\[
 F_0\,\mathbf Q(E_u[4])
 =F_0(\sqrt{1-u\theta_1},\sqrt{1-u\theta_2},\sqrt{1-u\theta_3}).
\]

At a good odd prime, the anchor's 4-torsion is finite étale, so F0 is
unramified. At a prime of F0 above p, exactly one of the three integral
factors \(1-u\theta_i\) has valuation one, and the others have valuation
zero. Irreducibility of f makes the Galois action transitive on the roots.
Conjugating this prime supplies three valuation functionals with matrix
\(I_3\) on these factors. Their squareclasses in \(F_0^\times/F_0^{\times2}\)
are therefore independent. Kummer theory gives degree eight and the stated
group. QED.

This is a proof by valuations, not a numerical field-degree estimate.
The [six hypothesis certificates](../../artifacts/generated-results/elliptic-curves/rank_jump_four_division_separation_v1.json)
use primes 6807347, 61, 19, 23, 233 and 23 for u=-3,-2,-1,1,2,3 respectively.
Each checks primality, good reduction, the simple root, and valuation one.
No degree-48 field or class group is constructed. In particular, none of
these deformations is 4-congruent to the anchor, even though each has its
specified 2-congruence. These elementary degree-four data are already
maximally new over F0, once E[2] and i are fixed.

This strengthens the explanation of the failed fixed-field transport.
Preserving a cubic field preserves only the first torsion layer; it did not
preserve the next layer relevant to halving and higher descent. It does
**not** deduce the observed CT ranks 12–16 from the three field characters,
prove any Sha dimension by itself, or make a 4-division field a rank predictor.
Nor would 4-congruence by itself prove simultaneous global solubility.

## Simultaneous solubility does not persist in a hidden common subspace

All local bases are embedded in the original twenty-class space. Let Ru be
the radical of the complete *restricted* CT matrix on Wu. A rationally
soluble inherited class must lie in Ru; the converse is unknown. The new
[intersection certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_local_collision_v1.json)
proves

\[
 R_u\cap R_v=0\quad
 \text{for every distinct }u,v\in\{-3,-2,-1,1,2,3\}.
\]

The locally admissible intersections \(W_u\cap W_v\) nevertheless have
sizes from dimension 7 to 15. For example:

| Pair | Common local dimension | Soluble dimension cap in common space on first fibre | Cap on second fibre | Joint cap |
|---|---:|---:|---:|---:|
| -3, -1 | 15 | 1 | 1 | 0 |
| -2, 1 | 7 | 0 | 0 | 0 |
| -1, 1 | 12 | 0 | 1 | 0 |
| 2, 3 | 11 | 1 | 0 | 0 |

Thus neither a stable labelled sub-block nor a stable sum of individually
bad basis classes can account for a transported rational block on two of
these fibres. At u=0 all twenty classes are rational; any individual
nonzero deformation still has its own one- or two-dimensional unresolved
restricted radical. Directions outside the inherited W remain unrestricted.
The joint local intersection across all six deformations is also zero,
but this alone would miss the much sharper pairwise solubility statement.

There is a useful experimental-design correction. For common C contained in
Wu, test the rectangular pairing **C × Wu**, not just C × C. Restricting both
arguments discards partners that can obstruct classes of C. All fifteen
nonzero parameter pairs have at least one strictly weaker self-pairing test.
On (-2,1), for example, the seven-dimensional common space has self-radical
dimension three on the first fibre, but pairing against all of W[-2] leaves
no candidate at all. The experiment retains those discarded-partner ranks
rather than mislabelling the smaller-space radicals as soluble candidates.

## What changed, and the next highest-value question

1. **Incidence now has an exact mechanism in this pencil.** A simple root
   collision imposes one specific unramified Kummer character; reciprocity
   coordinates the cuts. This replaces a loose “bad-prime pressure” account.
2. **Higher torsion preservation is ruled out sharply.** The fixed E[2]
   identification conceals a full degree-eight separation in E[4] data on
   every tested nonzero fibre. This is a concrete missing structure, not merely
   a warning that equal cubic fields need not preserve points.
3. **A persistent hidden rational sub-block is excluded on the finite panel.**
   Pairwise CT radical intersections are zero even when many classes survive
   the same local conditions.
4. **The positive implication remains missing.** None of these results supplies
   an independent globally soluble exceptional block on the recent R17/MW16
   discoveries. Their visibility clusters still cannot stand in for that result.

The next highest-value test is therefore about **higher descent information
on locally surviving classes**, not larger point boxes or more arbitrary
fixed-field parameters. A concrete bounded direction is to use the retained
Fisher entries to test whether the CT obstruction changes on common subspaces
can be expressed through the three explicit factors \(1-u\theta_i\) and their
local characters. First keep the already fixed pairs and compare exact
rectangular forms; any proposed factor-based formula must reproduce bilinear
entries, not just matrix ranks, and must pass u=0 plus an unseen retained pair.
Local terms of a chosen Fisher presentation are not automatically intrinsic
local pairings: a gauge/reciprocity check is required before assigning an
obstruction to a particular prime. If no invariant formula survives, the
next arithmetic input must be an independently constructed complement class
on the matched R17 high/low pair. No production scoring change is justified.

## Replay

```sh
python3 elliptic-curves/rank-jump/local_collision.py check
python3 elliptic-curves/rank-jump/reciprocity.py check
python3 elliptic-curves/rank-jump/four_division.py check
python3 -m unittest discover -s elliptic-curves/rank-jump -p 'test_*.py'
sage -python -m unittest discover -s elliptic-curves/rank-jump -p test_local_collision.py
```

All deterministic replays pass. The tests cover binary intersections against
exhaustive vectors, coefficient transport, the external-pairing-partner
failure mode, the nodal character-sum step, Hilbert bilinearity, and reciprocity
for small rationals. The optional Sage test independently compares **1176**
Hilbert symbols with PARI 2.17.3. Standard Python skips only that optional test.
The original expensive local-completeness and CT arithmetic certificates are
inherited, not rerun. New exact matrix algebra, root-character evaluations,
Hilbert products and prime hypotheses are recomputed.
