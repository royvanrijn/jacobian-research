# F2 `k=1` complement-monodromy stratification

> **Status.**  Exact computational-topological theorem.  The affine curve
> complement has fundamental group `Z` on every immersed, distinct-image
> `k=1` collision-partition stratum
> `4A_1`, `A_3+2A_1`, `2A_3`, `A_5+A_1`, and `A_7`, and also on the generic
> `A_2+3A_1` and `2A_2+2A_1` nonimmersion strata.  Consequently none of
> these seven strata can be the only ramified component of the nonproperness
> set of a connected plane Keller normalization.  A survivor needs a second
> ramified affine component and a second new source-boundary divisor, raising
> the stratum-conditional component floors from `28/49` to `29/50`.  The
> first topological escape is the `E_6+A_1` stratum: its noncyclic complement
> has an exact transitive degree-six permutation representation in which
> every geometric meridian fixes two sheets.  Thus complement monodromy
> alone does not exclude that stratum.  In its coarsest one-divisor
> `(e,f)=(2,2)` realization, the existing `ch_2` and unibranch ledgers give
> the conditional gate `14+2n-4b-s_X/2>=0`.  Multiple-component completions,
> severe nonimmersion/image-merger strata, and every chart `k>1` remain.

The exact implicit equations, certified braid words, Zariski--van Kampen
presentations, and degree-six escape representation are replayed by
[`verify_f2_affine_k1_complement_monodromy.py`](../scripts/verify_f2_affine_k1_complement_monodromy.py)
under SageMath with the optional `sirocco` package.

## 1. Equisingular transport

Use the normalized `k=1` parametrization

\[
 p=t^3+at,\qquad q=t^5+bt^4+ct^2+dt.             \tag{1.1}
\]

The collision-quartic theorem identifies a nonempty Zariski-open subset
`U_node` of the four-dimensional parameter space on which the affine image
has exactly four distinct ordinary nodes.  Its projective closure is a
rational quintic with one branch at infinity of type `(2,5)`.

The nonimmersion hypersurface is

\[
 \Delta_{\rm imm}=
 25a^4+48a^3b^2-144a^2bc+90a^2d+108ac^2+81d^2=0. \tag{1.2}
\]

It is irreducible: the marked-critical-point parametrization

\[
 a=-3r^2,\qquad d=-5r^4-4br^3-2cr              \tag{1.3}
\]

is dominant.  The coprime carrier-discriminant theorem gives a nonempty
open subset `U_cusp` of (1.2) on which the affine singularity packet is

\[
 A_2+3A_1.                                      \tag{1.4}
\]

Both `U_node` and `U_cusp` are connected over `C`.  The same will hold for
each collision-partition and bicuspidal stratum constructed below.  On each
one the projective divisor consisting of the quintic and the line at
infinity has constant embedded singularity data and the same `(2,5)` branch
with intersection order five against the line at infinity.  There are no
other singularities by the delta-genus ledger.

The complement fundamental group is consequently locally constant on each
stratum.  One direct justification is to choose disjoint Milnor balls around
the moving singular points.  The local families are topologically trivial,
and outside those balls the proper projective family is a stratified
submersion; Thom's first isotopy lemma glues the local trivializations.
Connectedness then transports one Zariski--van Kamp computation across the
whole stratum.

## 2. Four-node certificate

At

\[
 (a,b,c,d)=(1,0,0,0),\qquad
 (p,q)=(t^3+t,t^5),                              \tag{2.1}
\]

elimination gives

\[
 \boxed{F_{\rm node}=
 -P^5+5P^2Q+5PQ^2+Q^3+Q.}                      \tag{2.2}
\]

The collision theorem already certifies that this curve has four distinct
ordinary nodes.  Projection to `P` has three strands.  Certified homotopy
continuation gives the following six braid words in `B_3`; generator indices
are `1,2` and a negative index denotes the inverse:

\[
\begin{aligned}
 &(2,1,2,1,2,1,2,2,-1,-2,-1,-2,-1,-2),\\
 &(2,1,2,1,2,1,-2,-1,-1,-2),\\
 &(2,1,1,1,2,1,2,-1,-2,-1,-1,-2),\\
 &(2,1,-2),\\
 &(2,1,2,-1,-2),\\
 &(2,1,1,-2).
\end{aligned}                                                   \tag{2.3}
\]

The unsimplified van Kampen presentation contains the relators

\[
 x_2^{-1}x_0,\qquad x_1x_0^{-1},                \tag{2.4}
\]

and its remaining relators vanish after `x_0=x_1=x_2`.  Tietze reduction
therefore gives

\[
 \boxed{\pi_1(\mathbb A^2\setminus V(F_{\rm node}))
       =\langle x\mid\ \rangle\cong\mathbb Z.} \tag{2.5}
\]

There is no hidden finite-order relation: sending every geometric meridian
to `1` gives the standard surjection onto the first homology `Z` of an
irreducible affine plane-curve complement.

## 3. Ordinary-cusp certificate

Use the exact nonimmersion witness

\[
 (a,b,c,d)=(-3,0,1,-7),\qquad
 (p,q)=(t^3-3t,t^5+t^2-7t).                    \tag{3.1}
\]

It has one ordinary cusp and three distinct ordinary nodes.  Elimination
gives

\[
\boxed{\begin{aligned}
F_{\rm cusp}={}&P^5+3P^4+87P^3-66P^2Q+15PQ^2-Q^3\\
 &+79P^2-48PQ+6Q^2-7P+3Q.
\end{aligned}}                                                   \tag{3.2}
\]

The five certified three-strand braid words are

\[
\begin{aligned}
 &(1,2,1,2,1,2,1,2,1,2,-1,-2,-2,-1,-2,-1,-2,-1),\\
 &(1,2,1,2,1,2,2,2,-1,-2,-1,-2,-1),\\
 &(1,2,1,2,1,-2,1,1,2,-1,-2,-1,-2,-1),\\
 &(1,2,1,2,1,-2,-1,-1),\\
 &(1).
\end{aligned}                                                   \tag{3.3}
\]

The last word records the vertical tangency; the other four encode the
three nodes and cusp after transport to the common geometric basis.  Exact
Tietze reduction of the associated van Kampen presentation again gives

\[
 \boxed{\pi_1(\mathbb A^2\setminus V(F_{\rm cusp}))
       \cong\mathbb Z.}                         \tag{3.4}
\]

Thus the generic ordinary-cusp face has no extra monodromy freedom despite
its nontrivial local cusp braid.  The three node relations collapse the
global affine complement group to the single meridian.

### 3.1 All immersed collision partitions

The collision polynomial

\[
 R(u)=u^4+bu^3+au^2+(2ab-c)u-(a^2+d)              \tag{3.5}
\]

runs through **every** monic quartic.  Indeed, if its four lower
coefficients are `(r_3,r_2,r_1,r_0)`, the inverse parameter change is

\[
 b=r_3,\qquad a=r_2,\qquad
 c=2r_2r_3-r_1,\qquad d=-r_2^2-r_0.              \tag{3.6}
\]

Suppose a root `u_0` is off diagonal, both normalization branches are
immersed, and their image is distinct from every other collision value.
Then a root of multiplicity `m` gives two smooth branches of intersection
multiplicity `m`, hence an `A_(2m-1)` point.  This follows directly by
using `p` as a branch coordinate: after imposing the divided difference of
`p`, the divided difference of `q` is `-R(u)`, and `u-u_0` is a local
parameter on the collision branch.

Consequently the five partitions of four give exactly the following dense
equisingular strata.  The table also lists the exact witnesses used by the
checker.

\[
\begin{array}{c|c|c|c}
\text{root partition}&\text{affine packet}&(a,b,c,d)&R(u)\\ \hline
1+1+1+1&4A_1&(1,0,0,0)&u^4+u^2-1\\
2+1+1&A_3+2A_1&(0,0,4,-3)&(u-1)^2(u^2+2u+3)\\
2+2&2A_3&(2,0,0,-5)&(u^2+1)^2\\
3+1&A_5+A_1&(3,-3,-17,-9)&u(u-1)^3\\
4&A_7&(6,-4,-44,-37)&(u-1)^4
\end{array}                                                   \tag{3.7}
\]

For a fixed partition, the monic-quartic locus is irreducible.  Removing
the diagonal, nonimmersion, and equal-image subloci leaves a connected
Zariski-open equisingular stratum.  Exact braid monodromy and Tietze
reduction at the five witnesses give

\[
 \boxed{
 \pi_1(\mathbb A^2\setminus C_\lambda)=\mathbb Z
 \quad\text{for every immersed distinct-image partition }\lambda\vdash4.}
                                                               \tag{3.8}
\]

This closes all deeper tacnodal collision partitions in the immersed part
of the `k=1` chart; they are not residual exceptions to the four-node
calculation.

### 3.2 The generic bicuspidal stratum

The two critical points of `p` are both critical for `q` precisely, away
from the double-critical endpoint `a=0`, when

\[
 c=\frac{2ab}{3},\qquad d=-\frac{5a^2}{9}.       \tag{3.9}
\]

On this irreducible two-parameter locus the collision quartic factors as

\[
 \boxed{
 R(u)=\frac{(3u^2+4a)(3u^2+3bu-a)}9.}            \tag{3.10}
\]

The first factor consists of the two diagonal critical roots.  On a
nonempty open subset they give two ordinary cusps; the second factor gives
two distinct ordinary nodes.  Thus the generic packet is `2A_2+2A_1`.
At

\[
 (a,b,c,d)=(-3,0,0,-5),\qquad
 (p,q)=(t^3-3t,t^5-5t),                         \tag{3.11}
\]

one has `R=(u^2-4)(u^2+1)`.  Its exact van Kampen presentation again
reduces to one generator and no relations.  Connected equisingular
transport proves

\[
 \boxed{
 \pi_1(\mathbb A^2\setminus C_{2A_2+2A_1})=\mathbb Z.}        \tag{3.12}
\]

Thus adding a second ordinary cusp still does not create enough global
monodromy to connect the fixed affine sheets.

### 3.3 The first escape: `E_6+A_1`

The first noncyclic complement appears when the two critical points of `p`
coalesce.  The connected stratum is

\[
 a=c=d=0,\qquad b\ne0.                          \tag{3.13}
\]

Rescaling `t`, `P`, and `Q` identifies every member with the witness

\[
 (a,b,c,d)=(0,1,0,0),\qquad
 (p,q)=(t^3,t^5+t^4).                           \tag{3.14}
\]

Here `R=u^3(u+1)`.  The triple diagonal root has local parametrization
orders `(3,4)` and is an `E_6` cusp; the simple off-diagonal root gives one
node.  Thus the affine packet is `E_6+A_1`, and

\[
 F_{E_6+A_1}=P^5+P^4+3P^3Q-Q^3.                \tag{3.15}
\]

The exact three-meridian van Kampen group is noncyclic.  More decisively,
its pinned raw relators admit the following action on six sheets:

\[
 x_1\mapsto(3\ 4)(5\ 6),\qquad
 x_2\mapsto(3\ 5)(4\ 6),\qquad
 x_3\mapsto(1\ 3)(2\ 6).                      \tag{3.16}
\]

The generated action is transitive and every geometric meridian fixes two
sheets.  More structurally, label the six sheets by the six two-element
subsets of a four-element set.  The three displayed meridians are the edge
actions induced by the adjacent transpositions `(3 4)`, `(1 2)`, and
`(2 3)` in `S_4`.  Hence the image is exactly `S_4` in its transitive
six-edge action.  A vertex transposition fixes its own edge and the
complementary edge while exchanging the remaining four edges in two pairs.
It therefore realizes the exact permutation-theoretic pattern
that the cyclic-complement lemma forbids: boundary degree four plus affine
remainder two at total degree six.  This is not a finite Keller cover, but
it proves that complement monodromy by itself cannot eliminate this
stratum.  The `E_6+A_1` locus must instead be attacked by its boundary
attachment and logarithmic `ch_2` packet.

There is already a useful conditional numerical gate.  The coarsest
one-divisor interpretation of (3.16) has one boundary row `(e,f)=(2,2)`
and affine remainder two.  If this divisor has `E^2=-n`, follows `b` of the
eight smooth carrier centers, and the completed source adds `s_X` further
smooth-boundary blowups, its cyclic divisorial term is

\[
 A_{E_6}=4(b-7)-2n.                              \tag{3.17}
\]

At squarefree degree six the root-subtracted global residual is `-10`, so
an exact root-plus-one-row filtration would have total finite length

\[
 \ell(Z)=18+2n-4b-\frac{s_X}{2}.                 \tag{3.18}
\]

The `E_6` branch multiplicity is three.  The isolated unibranch attachment
theorem therefore forces at least `(3-1)f=4` of (3.18), leaving

\[
 \boxed{
 \ell(Z_{\rm rest})=14+2n-4b-\frac{s_X}{2}\ge0.} \tag{3.19}
\]

For the minimal self-intersection `n=1`, this gives
`8b+s_X<=32`; in particular `b>=5` is impossible in this conditional
one-row realization.  If the two transpositions in (3.16) arise from two
distinct residue-degree-one boundary divisors instead, both divisorial
packets must be booked separately and (3.19) is not the applicable formula.
This dichotomy is now the exact source-side question at the first escape.

## 4. Cyclic-complement obstruction

The following lemma is independent of F2.

### Lemma 4.1 -- a fixed sheet needs another branch component

Let

\[
 \pi:S\longrightarrow\mathbb A^2                 \tag{4.1}
\]

be a connected finite normal cover of degree `d>1`.  Let `C` be an
irreducible component of its branch locus and suppose

\[
 \pi_1(\mathbb A^2\setminus C)=\mathbb Z.         \tag{4.2}
\]

If local monodromy around `C` has a fixed sheet, then the branch locus has
another irreducible component.

#### Proof

If `C` were the entire branch locus, the restriction over
`A^2 minus C` would be a connected finite etale cover.  Its monodromy would
be a transitive permutation representation of `Z`, hence generated by one
permutation `sigma`.  Transitivity of the cyclic group generated by `sigma`
is equivalent to `sigma` being a single `d`-cycle.  Such a permutation has
no fixed point.  This contradicts the fixed-sheet hypothesis.  \(\square\)

The lemma only needs `C` to be the entire **branch** locus.  Additional
unramified nonproperness components do not help: the finite cover extends
etale across them, so they do not add monodromy generators.

## 5. F2 consequence

For a plane Keller map, the canonical finite normalization is connected and
the affine source is etale.  Over the generic point of a nonproperness
component `C`, its finite-flat ledger is

\[
 d=\sum_i e_if_i+\sum_j a_j,
 \qquad \sum_j a_j\ge1.                         \tag{5.1}
\]

The terms `a_j` are affine sheets.  Keller etaleness makes them fixed points
of the local meridian around `C`.  Purity supplies at least one row with
`e_i>1`, so `C` is a branch component.

Apply Lemma 4.1 and (2.5), (3.4), (3.8), or (3.12).  On any of the seven
connected cyclic-complement `k=1` strata,

\[
 \boxed{\text{the purity curve cannot be the only ramified affine
 nonproperness component}.}                    \tag{5.2}
\]

A second branch component is itself an affine nonproperness component and
is dominated by a distinct ramified source-boundary divisor.  None of the
certified `27/48` F2 components can serve: they map to target infinity, to a
point, or are already log-etale.  Hence the stratum-conditional source floors
are

\[
 \boxed{
 N_{\rm source}^{\rm squarefree}\ge29,\qquad
 N_{\rm source}^{\rm double}\ge50.}             \tag{5.3}
\]

This is stronger than the unconditional purity floors `28/49`, but it is not
yet an exclusion.  The second ramified component can connect the fixed
affine sheets to the first component's inertia orbits.  Its target chart,
boundary signature, and logarithmic Chern packet must be added to the common
model before the global budget is evaluated.

## 6. Claim boundary and next attack

The result closes the following `k=1` faces as **single-row** completions:

1. every immersed, distinct-image collision partition of four;
2. the generic ordinary-cusp face `A_2+3A_1`; and
3. the generic bicuspidal face `2A_2+2A_1`.

The tacnode, double-tacnode, `A_5`, and `A_7` strata are therefore no longer
gaps.  What remains inside `k=1` is qualitatively smaller: severe
nonimmersion strata beginning with `E_6+A_1`, collision values that merge
into multibranch points, their intersections, and the concentrated monomial
`E_8` cusp `(p,q)=(t^3,t^5)`.  The last two named cusp strata have noncyclic
complements.  Nor does the argument control the 23 charts `k=2,...,24`.

The next finite computation now has a sharper split.  On the cyclic strata,
a two-row logarithmic `ch_2` compiler must include both new affine
divisorial packets rather than continuing the one-packet residual formula.
On the noncyclic `E_6/E_8` cusp strata, complement monodromy must be combined
with the unibranch boundary lower ledger.  For `E_6`, the branch
multiplicity is three, so every complete isolated attachment fiber already
costs at least `2f`; the source incidence decides whether the exact charge
is larger.  Target-value merger strata still require one complement-group
certificate per connected equisingular locus.

The concentrated `E_8` endpoint is now classified separately in
[`F2_AFFINE_K1_E8_MONODROMY.md`](F2_AFFINE_K1_E8_MONODROMY.md).  Among
transitive degree-six actions whose meridian has the required type
`2+2+1+1`, its torus-knot group has exactly one conjugacy class, with image
`A_5`.  Thus topology permits the endpoint but removes all permutation
freedom.  Its preferred longitude equals the meridian in this action, so
the peripheral orbits separate the two transposition cycles: the endpoint
forces two distinct `(e,f)=(2,1)` source-boundary divisors and excludes a
single `(2,2)` divisor.

The degree-six action is not the only fixed-sheet action of its image.
Exhausting all `A_5` coset actions gives precisely the F2 degrees
`6,10,15,30`, with `2,4,6,14` distinct `(2,1)` ramified rows.  Thus
complement topology does not force the degree-six floor.  This atlas is
exhaustive inside the icosahedral quotient, not among all finite quotients
of the E8 cusp group.

<!-- status-consumer: PF2K1E8M1 bbb282c6bcfa62fc -->

For simple inertia the finite-quotient gap is now closed.  Adding `M^2=1`
to the E8 cusp group gives a universal group of order `240` with center of
order four.  Its `30` subgroup classes give exactly `13` fixed-sheet F2
actions in degrees `6,10,12,15,20,24,30,40,60,120`.  The preferred
longitude detects the central gluing and produces `(2,f)` rows with
`f=1,2,4`.  What remains outside this atlas is inertia greater than two.

<!-- status-consumer: PF2K1E8O1 4251750ed4e43c89 -->

## Reproduction

With SageMath and `sirocco` available, run

```bash
sage -python scripts/verify_f2_affine_k1_complement_monodromy.py
```

The checker verifies seven cyclic resultants, their pinned braid words and
unsimplified van Kampen presentations, and their exact Tietze reductions to
one generator with no relation.  It separately verifies the `E_6+A_1`
presentation and the transitive degree-six fixed-sheet action (3.16).

## References

- SageMath,
  [*Zariski--Van Kampen method implementation*](https://doc.sagemath.org/html/en/reference/curves/sage/schemes/curves/zariski_vankampen.html),
  for the certified braid-continuation and presentation routines.
- Thom's first isotopy lemma, applied to the proper stratified projective
  family after isolating the moving plane-curve singularities.
