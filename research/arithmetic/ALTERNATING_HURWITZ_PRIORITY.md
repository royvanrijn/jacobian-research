# Priority comparison for alternating Laurent covers and collinear nodes

The abstract Laurent passport in the bridge note has established precedent.
Edmonds--Kulkarni--Stong already imply its realizability, and Pakovich's
classification gives a second direct route. In every even degree at least ten,
full alternating monodromy then follows from a classical primitive-group theorem.
The explicit knot tuples, their braid stabilizers, and the exact degree-eight
coefficient certificates remain useful concrete data. This comparison does not
assert priority for those explicit formulas or for a newly produced rational
family.

## The exact passport

Write n=2b, b>=4. The branch data are

\[
 \Pi_b=\{\underbrace{(2,2,1^{n-4}),\ldots,(2,2,1^{n-4})}_{b\text{ entries}},
             (n-1,1)\}.
\]

There are q=b+1>=5 branch values. Each finite partition has defect two and
the last has defect n-2. The total defect is 2b+n-2=2n-2, so a connected
realization has genus zero. The entries of the two-part partition are
unequal and distinguish the two source poles.

## Exact conclusions already implied by prior theorems

Edmonds, Kulkarni and Stong, *Realizability of branched coverings of surfaces*,
Transactions AMS 282 (1984), 773--790,
[Proposition 5.3, pages 784--785](https://courses.mai.liu.se/GU/TA1022/Edmonds-Branched%20coverings.pdf),
classify realizability when one partition is (d-1,1). In their notation the
conditions are total defect at least 2d-2, even total defect, and exclusion
of two explicit exceptional patterns. Here d=n>=8 and q>=5, so neither
exception applies. The total defect is exactly 2n-2. The proposition therefore
already supplies a connected sphere cover with the exact passport Pi_b.

Pakovich, *Solution of the Hurwitz problem for Laurent polynomials*,
Journal of Knot Theory and Its Ramifications 18 (2009), 271--302,
[Theorem 1.1, page 273](https://pakovich.github.io/lau.pdf),
proves realizability of every genus-zero Laurent passport with more than three
branch values. Our passport has q>=5, and the partition (n-1,1) has length two,
so it is an immediate specialization. This is not merely broad similarity:
the partitions and Riemann--Hurwitz equality match exactly. Section 3 contains
the constructive constellation argument; Theorem 3.7 on page 291 completes
the q>3 case. The introduction on page 272 explicitly attributes the
(n-1,1) case to Edmonds--Kulkarni--Stong.

These are complex/topological realizability statements. After Riemann existence
and algebraic choice of branch points one obtains an algebraic-number model,
as in the bridge note. Neither statement supplies a rational coefficient model,
an arithmetic alternating group over Q, or a prescribed generating tuple.

For n>=10, the full geometric alternating group also follows from these
realizability results. Let G be the transitive monodromy group of any such
realization. It contains an (n-1)-cycle fixing one point. A block containing
that fixed point is invariant under the cycle, so is either the singleton or
the whole n-point set; thus G is primitive. It contains a double transposition,
and the classical small-support theorem gives A_n <= G for n>=9. All the
branch permutations are even, so G <= A_n. Consequently G=A_n.

The exact small-support formulation is recorded in Gareth A. Jones,
*Characteristically simple Beauville groups, I: cartesian powers of alternating
groups*, [arXiv:1304.5444v1, Proposition 6.1(3), page 9](https://arxiv.org/pdf/1304.5444v1).
Jones explicitly imports this classical result and points to Cameron's proof
sketch. This source comparison imports the theorem; it does not independently
reprove its small-support classification. The primitivity and parity
application above are elementary and complete.

Degree eight is a genuine exception to that automatic group inference:
AGL_3(2) is primitive on eight points, lies in A8, and has double
transpositions and a seven-cycle. Thus passport realizability alone cannot
certify A8 in degree eight. The knot tuple generation proof and the exact
degree-eight monodromy certificate provide that additional information.

The maximum-collinear-node existence statement is a geometric consequence of
the already-realizable passport. Normalize its poles to obtain
f=(H+a)/z, H(0)=0, a!=0. The finite branch partitions give b pairs of simple
critical points with equal value. The discriminant parametrization
(H'(r),rH'(r)-H(r)) then gives b ordinary nodes on t=a, and the degree-n
Bezout bound is saturated, exactly as proved in the bridge note. This geometry
does not require G=A_n, so Edmonds--Kulkarni--Stong already imply its existence
also in degree eight. The collinear-node formulation is not explicitly stated
in the inspected sources; no claim of priority for that formulation follows
from its absence.

## Comparison with the other closest constructions

| Source and exact locator | What it supplies | Difference from the retained bridge |
|---|---|---|
| Edmonds--Kulkarni--Stong, Proposition 5.3, pp. 784--785 | Realizability with an (n-1,1) partition, subject to the complete listed exceptions | Directly implies this exact passport; does not prescribe full group in degree eight or rational coefficients |
| Pakovich, Theorem 1.1, p. 273 | Every genus-zero Laurent passport with q>3 is realizable | Directly implies this exact passport; gives no required arithmetic field or particular Nielsen representative |
| Magaard--Voelklein, Lemma 3.1(i), arXiv:math/0304130v1, p. 5 | For every n>=6 there are genus-zero full A_n covers with all branch cycles double transpositions | All n-1 branch cycles have index two; there is no specified long-cycle two-pole fiber |
| Wewers--Bouw, Proposition 7 and Theorem 10, published pp. 193--194 and 197 | Alternating covers in characteristic p with prescribed all-3-cycle or all-double-transposition inertia in its stated ranges | Does not impose the (n-1,1) pole fiber; its positive-characteristic degree and characteristic restrictions must be retained |
| Mestre, ICM 1994, Section 2, p. 319 | For odd n and general rational P, a rational pencil P-tQ with discriminant Disc(P) times a square; generic branch polynomials have a triple root | The denominator Q generally has degree up to n-1. Its Wronskian is a square and its finite local permutations are 3-cycles |
| Mueller, Theorem 4.8, published pp. 416--418 | Classification of indecomposable two-pole monodromy, allowing the natural alternating/symmetric case with many ramification possibilities | Broad Laurent alternating monodromy is established context; that classification does not isolate this passport or provide its coefficients |

The all-double-transposition result of Magaard--Voelklein is sometimes hidden
by citations to their positive-genus Theorem 3.3. The directly relevant
genus-zero statement is Lemma 3.1, including its explicit degree-raising
operation. Its n=6 base verification is delegated to GAP or a direct check in
the paper; no new replay of that auxiliary base computation is asserted here.

Mestre's square discriminants arise by a different local mechanism. On page
319 he constructs Q,R satisfying P'Q-PQ'=R^2, with R squarefree, and states that
each branch specialization has one triple root and otherwise simple roots.
Thus a square discriminant does not by itself identify our double-transposition
passport. In the current pencil the critical points are simple; the square
occurs because two distinct critical points share a value.

Primary context sources:

- [Magaard--Voelklein, arXiv:math/0304130v1](https://arxiv.org/pdf/math/0304130v1).
- [Wewers--Bouw, published article](https://msp.org/pjm/2005/222-1/pjm-v222-n1-p07-p.pdf).
- [Mestre, Constructions polynomiales et theorie de Galois, ICM 1994](https://www.imj-prg.fr/wp-content/uploads/2020/prix/mestre1994.pdf).
- [Mueller, published article](https://www.numdam.org/item/ASNSP_2013_5_12_2_369_0.pdf).

## What the bounded audit settles

The passport existence theorem is not a new existence theorem. For b>=5,
full alternating monodromy is also already a consequence of established
theorems, and maximum-collinear-node existence follows geometrically for
every b>=4. The contribution that can safely be described here is an explicit
bridge and certificate tying the knot generators and braid action to these
covers, together with concrete polynomial coefficient realizations and
separately proved descent statements.

The inspected sources do not supply or classify the specific degree-eight
coefficient formulas in the current notes. They do not identify their braid
orbits. They do not decide whether a particular exact rational coefficient
family appeared elsewhere. Those narrower explicit-construction priority
questions remain open after this bounded audit. Search failure is not a
nonexistence or originality proof.

The [source ledger](alternating_hurwitz_priority_ledger.json) pins each retained
version, read locator, exact application, and limitation. This is internal
primary-source checking, without external human review or an exhaustive
historical search.
