# A fitted secant pencil changes the generic subgroup

On all six retained high/low controls, the first two marked generic points
produce an oblique split-cubic pencil of arithmetic rank **0**, whose square
base change has arithmetic rank **2**. This is a genuine simultaneous
solubility construction. Its two directions contribute **zero** to the
exceptional quotient of the original rank-16 or rank-17 family.

Thus the construction in [the preceding note](OBLIQUE_BLOCKS_BEYOND_THE_HORIZONTAL_OBSTRUCTION.md)
can be fitted to both high-gain and observed-zero fibres without using any
exceptional points. Its existence after fitting is not evidence for the
original rank jump. The parent family and generic subgroup must be fixed
before a block can be credited to an exceptional quotient.

## Frozen comparison

The [protocol](RETROSPECTIVE_SECANT_PENCIL_PROTOCOL.json) fixes one pencil
per curve, using exactly the first two marked generic points. There are no
new point searches, parameter searches, or adaptive choices. Each capture
worker has a 30-second cap. Inputs and outputs are hash-bound and retained
independently of the live search.

Here (m) is the original generic rank, (w) the certified independent
specialized subgroup rank, and (q=w-m) its certified quotient rank.
Neither (w) nor (q) asserts a full specialized curve-rank upper bound.
Observed-zero controls remain censored by their retained search evidence.

| Original family | Original parameter | (m) | (w) | (q) | Fitted pencil rank, before / after square cover | (w-2) |
|---|---:|---:|---:|---:|---:|---:|
| A1/MW16-05 | 307/206 | 16 | 25 | 9 | 0 / 2 | 23 |
| A1/MW16-05 | -3158/1291 | 16 | 16 | 0 observed | 0 / 2 | 14 |
| A1/MW16-04 | -1647/91 | 16 | 25 | 9 | 0 / 2 | 23 |
| A1/MW16-04 | -2177/2397 | 16 | 16 | 0 observed | 0 / 2 | 14 |
| published R17 | -2300/843 | 17 | 24 | 7 | 0 / 2 | 22 |
| published R17 | -1561/3133 | 17 | 17 | 0 observed | 0 / 2 | 15 |

The last column illustrates the accounting hazard: changing the parent
creates an apparent known specialization gain of 14 or 15 even on the
observed-zero controls. Precisely,

\[
w-2=(m-2)+q.
\]

The first (m-2) directions were already generic for the original family.
The fitted construction has not explained them, or the remaining (q).

## Exact construction and rank proof

Write a retained short model as (E:y^2=F(x)=x^3+Ax+B). For the selected
independent rational points (P=(a,p)), (Q=(b,q)), set

\[
s=\frac{q-p}{b-a},\qquad h=p-sa,\qquad L(x)=sx+h.
\]

Their secant gives the identity

\[
F(x)-L(x)^2=(x-a)(x-b)(x-c),\qquad c=s^2-a-b.
\]

All six cases have (s\ne0), three distinct roots, and

\[
x_0=-h/s,\qquad C=F(x_0)\ne0.
\]

Define the **new** pencil

\[
\mathcal E_t:\quad y^2=F(x)+(t-1)L(x)^2.
\]

At (t=1) it recovers the original curve. At (t=n^2), it has sections

\[
P_n=(a,np),\quad Q_n=(b,nq),\quad T_n=(c,nL(c)),
\qquad P_n+Q_n+T_n=O.
\]

The first two are independent: a generic integral relation would
specialize at the smooth fibre (n=1) to a relation between the two
certified independent input points. Both are anti-invariant under
(n\mapsto-n).

The program computes the discriminant and (c_4) over \(\mathbb Q[t]\)
exactly. In every case the discriminant is squarefree of degree three,
coprime to (c_4), which has degree two; its value at zero is nonzero.
The finite singular fibres are three (I_1)'s. At infinity the minimal
orders are (v(c_4)=2,v(\Delta)=9), giving (I_3^*). This is a rational
elliptic surface, so Shioda--Tate gives geometric generic rank

\[
10-2-7=1.
\]

The constant geometric section (R=(x_0,\sqrt C)) is non-torsion. It is
disjoint from (O); the finite (I_1)'s contribute zero to its height.
The nonidentity simple components of (D_7) have corrections
(1,7/4,7/4). Thus its height is at least (2-7/4=1/4).
This is a lower bound, not an asserted exact height. Since (R) spans
the rank-one geometric space and its constant-field Galois character is
that of \(\sqrt C\), the arithmetic generic rank is

\[
\operatorname{rank}\mathcal E(\mathbb Q(t))=[C\in\mathbb Q^{\times2}].
\]

After (t=n^2), the exact pullback discriminant has six simple finite
zeros, while infinity has minimal orders (0,6), hence fibre (I_6).
The new surface is again rational, of geometric rank (10-2-5=3).
The invariant (R) and anti-invariant independent (P_n,Q_n) span this
space. Consequently

\[
\operatorname{rank}\mathcal E_{n^2}(\mathbb Q(n))
=2+[C\in\mathbb Q^{\times2}].
\]

The six exact rational values of (C) are all nonsquares. This proves
the reported arithmetic ranks (0\to2), without a numerical height
estimate, point search, or full descent on any production curve.
The elliptic-surface and height formulas used here are the standard
ones in [Schütt--Shioda, *Elliptic surfaces*](https://arxiv.org/pdf/0907.0298),
especially the discussions of Shioda--Tate and the Mordell--Weil lattice.

## What would make such a block relevant?

The fitted square cover is an **incidence and solubility construction for
the new pencil**. It is neither a visibility statistic nor an incidence
test for the original exceptional quotient. Its matched-pair result is
negative: all three pairs pass identically.

A relevant shared-cover mechanism must instead supply the following
data on the original parameter line (u):

1. A common nontrivial squareclass (d(u)), specified without exceptional
   points, and algebraic sections defined over its cover.
2. At least two independent directions modulo the original generic
   subgroup, with dependence and specialization collapse controlled.
3. A specialization condition making the cover rationally soluble,
   together with a proof that the directions survive on that fibre.

This would give the desired chain from a specialization condition to a
multidimensional soluble block. Merely fitting a split residual cubic at
one selected fibre supplies none of these original-family implications.

One narrower candidate remains open. A secant of original generic
sections (P(u),Q(u)) defines (x_0(u)=-h(u)/s(u)) and
(C(u)=F_u(x_0(u))). The section
((x_0(u),\sqrt{C(u)})) lies on the **original** family over that quadratic
cover. If several distinct constructions shared a nontrivial squareclass
and were independent modulo the original generic subgroup, a square
specialization could expose a block. No such shared identity or
independence statement has been established here. The six first-pair
specialized values are nonsquares. No selector is proposed from them.

## Ranked conclusions and remaining work

1. **Strongest proved model mechanism — incidence plus solubility:** a
   shared quadratic cover can add two independent rational directions.
   Its relevance to production exceptional quotients still needs a
   construction over the original family, as specified above.
2. **Strongest production structure — incidence:** the retained strict
   Selmer/ideal-class blocks and exact local contraction theorems explain
   substantial shared arithmetic structure. They do not decide whether
   the remaining global classes are rational or Sha; see the
   [production twist analysis](PRODUCTION_TWIST_INCIDENCE_AND_SOLUBILITY.md).
3. **Weak or excluded explanations:** the fitted secant event does not
   distinguish any paired controls; horizontal shared-value blocks are
   excluded for A1/MW16-05 by the existing norm obstruction; halving-field
   compression and chart recovery alone do not establish incidence.
4. **Missing computations and theorems:** original-family shared-cover
   identities and quotient independence; global solubility of the
   relevant covers rather than just Selmer or Jacobian lifting; complete
   production class-group information where absolute Selmer dimensions
   remain conditional. None is replaced by the present fitted pencils.
5. **Information for Agent 1:** retain the original generic subgroup in
   every proposed block score; require a coefficient/generic-section
   construction before exceptional points are supplied; label a local
   test as incidence and a rational-cover test as solubility. The fitted
   secant predicate has no demonstrated selection value and should not
   become a candidate score.

## Evidence and replay

- [Exact frozen inputs](../../artifacts/generated-results/elliptic-curves/rank_jump_retrospective_secant_pencil_inputs_v1.json)
- [Six pencil certificates](../../artifacts/generated-results/elliptic-curves/rank_jump_retrospective_secant_pencils_v1.json)
- [Independent verification](../../artifacts/generated-results/elliptic-curves/rank_jump_retrospective_secant_verification_v1.json)

The producer verifies the fibre-configuration gates and section identities
using Sage. The independent verifier recomputes the polynomial invariants
with rational list arithmetic, checks the original generic-point source
and rank-two fingerprints, uses Sage's separate square test, and checks
the (D_7) correction values from its inverse Cartan matrix. These are
algebraic certificates plus the displayed elliptic-surface argument, not
a formal proof-assistant derivation.

```sh
sage -python elliptic-curves/rank-jump/retrospective_secant_pencils.py check
sage -python elliptic-curves/rank-jump/verify_retrospective_secant_pencils.py check
python3 -m unittest discover -s elliptic-curves/rank-jump -p test_retrospective_secant_pencils.py
```

All active search files, outputs, policies, and mathematical status entries
are outside this change.
