# Attack program for the restricted exact minima

## 1. One coupled frontier, not four unrelated searches

The cubic rank and index questions are linked by an elementary but useful
constraint.  If a nilpotent matrix has rank `r` and nilpotency index `nu`,
then

\[
 \nu\leq r+1.
\]

Indeed, a Jordan block of length `nu` already contributes `nu-1` to the
rank.  The certified rank-17 witness has index 18 and therefore saturates
this inequality.  No coordinate compression of that witness can lower its
rank without first shortening the chain.

This makes the positive index program the master lower-bound attack.  A
theorem that every cubic-homogeneous Keller map of index at most `d` is
invertible would simultaneously give

\[
 \nu_{\rm cub}\geq d+1,\qquad r_{\rm cub}\geq d.
\]

The first missing case is `d=3` for arbitrary cubic-homogeneous corrections.
The published index-three theorem covers power-linear maps only, so it
cannot be used at this step
([Adamus--Bogdan--Crespo--Hajto](https://arxiv.org/abs/1508.02012)).

## 2. Attack A: the full cubic index-three identity ideal

Write the cubic correction as a symmetric vector-valued tensor

\[
 H(x)=T(x,x,x),\qquad N(x)=JH(x)=3T(x,x,-).
\]

The condition `N(x)^3=0` is a finite homogeneous coefficient ideal in the
entries of `T`.  The inverse of `X+H` is a rooted-tree series whose first
terms are

\[
\begin{aligned}
K_3&=-H,\\
K_5&=N H,\\
K_7&=-N^2H-\tfrac12\,\operatorname{Hess}(H)[H,H].
\end{aligned}
\]

The concrete program is:

1. polarize every coefficient of `N(x)^3=0`, retaining the
   `GL(V)`-module labels rather than expanding in a fixed large dimension;
2. generate inverse-tree covariants degree by degree;
3. reduce those covariants modulo the polarized nilpotency ideal;
4. either prove a uniform truncation bound or output the first surviving
   covariant as the exact obstruction to an index-three theorem.

This attack has an unambiguous stop condition.  Vanishing of all tree
modules through a proven degree bound closes the case.  A surviving module
must be realized by an exact tensor or shown incompatible with the Keller
trace identities.  Random examples and one-point matrix powers are not
evidence here.

The useful intermediate theorem is narrower:

> Classify the additional tensor identity, beyond `N^3=0`, that kills
> `K_7`, and prove invertibility on that stratum.

It splits the full problem into a finite representation-theoretic
calculation instead of asking directly for an inverse formula.

The exact triangular calibration

\[
F=(x_0+x_1^3,\ x_1+x_2^3,\ x_2)
\]

has `(JH)^2!=0`, `(JH)^3=0`, and inverse

\[
F^{-1}=(y_0-(y_1-y_2^3)^3,\ y_1-y_2^3,\ y_2).
\]

Its inverse correction has nonzero homogeneous parts in degrees
`3,5,7,9`.  In particular, the degree-nine term `y_2^9` proves that no
uniform index-three inverse bound can be smaller than nine.  This is checked
by
[`verify_index_three_inverse_model.py`](../scripts/verify_index_three_inverse_model.py)
and recorded in
[`index_three_inverse_degree_model.json`](../artifacts/generated-results/index_three_inverse_degree_model.json).
Thus “all inverse-tree covariants above degree nine vanish” is the first
plausible sharp theorem target; degree-seven truncation is already false.

The first universal tree reduction is now implemented in
[`derive_index_three_tree_obstruction.py`](../scripts/derive_index_three_tree_obstruction.py).
It generates the exact inverse recursion through degree eleven.  At degrees
`3,5,7,9,11` there are respectively `1,1,2,4,8` non-plane cubic derivative
trees.  Euler homogeneity strengthens `N^3=0` to

\[
 N^2H=0,
\]

because `NX=3H`.  Closing this identity, together with the directional
derivative of `N^3v=0`, under differentiation and insertion into cubic
contexts gives seven degree-eleven relations of rank five.  They do **not**
yet kill `K_11`.  In the resulting three-dimensional quotient its normal
form is

\[
\boxed{
\Omega_{11}=
-\frac12 C(NH,H,H)
-\frac12 B(B(H,H),H)
-\frac16 N(C(H,H,H)),
}
\]

where `B=D^2H` and `C=D^3H`.  The complete rational tree ledger is
[`index_three_inverse_tree_obstruction.json`](../artifacts/generated-results/index_three_inverse_tree_obstruction.json).

This replaces the broad index-three question by the first exact fork:

- prove `Omega_11=0` using additional polarizations of the coefficient ideal
  of `N^3=0`, then continue to degree 13;
- or realize `Omega_11!=0` by an exact cubic tensor satisfying `N^3=0`,
  disproving degree-nine truncation and exposing the next inverse term.

The computation does not claim that `Omega_11` survives the full
nilpotency ideal.  It proves that Euler, first polarization, differentiation,
and context closure alone are insufficient.

## 3. Attack B: rank three through the kernel bundle

For `rank(JH)=3`, the generic kernel has codimension three and contains the
image because `JH` is nilpotent.  The known rank-two classification should
not be extended by brute-force coefficients.  The new datum at rank three is
the rational kernel map

\[
 x\longmapsto\ker JH(x)\in\operatorname{Gr}(n-3,n).
\]

There are two strata:

- **constant-kernel stratum:** quotient the common kernel first and classify
  the resulting essential three-variable tensor;
- **moving-kernel stratum:** use the quadratic entries of `JH` and the
  inclusions `im(JH)\subseteq ker(JH)` to bound the Pluecker degree of the
  kernel map, then classify its lowest-degree images.

A practical calculation should enumerate Pluecker degree, not ambient
dimension.  The first milestone is to prove that degree zero is invertible
and to eliminate degree one.  Any surviving degree-one normal form becomes
a small exact ansatz for collision equations.  This provides either a
rank-three counterexample or the first improvement `r_cub>=4`.

One substantial four-variable ansatz has now been eliminated.  Consider

\[
H=(0,\lambda x_1^3,\ x_2\Delta+p(x_1,x_2),\
   x_1\Delta+q(x_1,x_2)),\qquad
\Delta=x_1x_3-x_2x_4,
\]

with arbitrary binary cubics `p,q`, a family arising in the homogeneous
nilpotent-Jacobian normal-form analysis
([de Bondt--Yan](https://arxiv.org/abs/1302.6930)).
Exact coefficient extraction from `(JH)^3` contains both `3*lambda` and
`-3*lambda`.  Hence index three forces `lambda=0`.  On that locus only two
components of `H` remain, so `rank(JH)<=2`; moreover the last-coordinate
block is `z -> (I+B)z+c` with `B^2=0`, giving the explicit inverse
`z=(I-B)(y-c)` of degree at most five.

Thus no member of this whole family lies simultaneously on the rank-three
and index-three strata: the rank-three branch has index at least four.  The
exact certificate is
[`verify_index_three_rank_normal_form.py`](../scripts/verify_index_three_rank_normal_form.py),
with artifact
[`index_three_rank_normal_form_exclusion.json`](../artifacts/generated-results/index_three_rank_normal_form_exclusion.json).
The next kernel-bundle step is to determine whether this family exhausts
the Pluecker-degree-one stratum or to write the missing normal forms.

## 4. Attack C: cotangent kernel excess

For the cotangent potential `p(x,y)=y^T H(x)`, put

\[
 N=JH,\qquad A=\sum_i y_i\operatorname{Hess}(H_i),\qquad
 M=\operatorname{Hess}p=
 \begin{pmatrix}A&N^T\\N&0\end{pmatrix}.
\]

If the columns of `K` span `ker(N)`, exact block elimination gives

\[
 \boxed{\operatorname{rank}M
 =2\operatorname{rank}N+\operatorname{rank}(K^TAK)}.
\]

The last summand is the **cotangent kernel excess**.  The checker
[`analyze_cotangent_kernel_excess.py`](../scripts/analyze_cotangent_kernel_excess.py)
verifies the block identity at twelve independent good-prime
specializations across the four certified sources and writes
[`cotangent_kernel_excess_frontier.json`](../artifacts/generated-results/cotangent_kernel_excess_frontier.json).

The current profiles are:

| cubic source | `rank JH` | kernel excess | cotangent Hessian rank |
|---|---:|---:|---:|
| ambient-21 witness | 18 | 2 | 38 |
| index-18 witness | 18 | 2 | 38 |
| rank-17 witness | 17 | 4 | 38 |
| Hessian-rank witness | 18 | 1 | 37 |

For the last row, the excess is exact rather than merely sampled: the
characteristic-zero generic ranks 18 and 37 are independently certified, so
their difference \(37-2\cdot18=1\) is exact.

The combined nine-atom, width-64 search has 140 terminal maps.  Its sampled
kernel-excess histogram is

| excess | terminals |
|---:|---:|
| 1 | 1 |
| 2 | 8 |
| 3 | 17 |
| 4 | 31 |
| 5 | 41 |
| 6 | 23 |
| 7 | 15 |
| 8 | 4 |

The `qb+x2s` witness is the unique excess-one terminal; no terminal has
excess zero.  This is a bounded diagnostic, not a lower bound, but it shows
that the next atom must be synthesized to kill the restricted form rather
than selected by widening the same beam.

Consequently the next strict Hessian improvement has only two immediate
targets:

- rank 18 and excess 0, giving Hessian rank 36;
- rank 17 and excess at most 2, again giving rank at most 36.

At fixed cubic rank 18, no cotangent construction can improve 37 by more
than one.  This replaces raw Hessian-rank search by the geometric condition
`K^TAK=0` on the moving Jacobian kernel.

The first implemented hybrid forces the zero-excess `xxs` gate together
with both rank-reducing atoms.  It reaches 48 terminals but its best profile
is `(rank,index,Hessian rank,dimension)=(22,19,45,28)`, with excess one.
The pairwise beams also terminate without an improvement:

| forced atoms | terminals | best `(rank,index,Hessian,dimension)` | excess |
|---|---:|---:|---:|
| `xxs+v2r` | 32 | `(21,19,44,27)` | 2 |
| `xxs+y2vb` | 34 | `(19,18,43,26)` | 5 |

Thus the known zero-excess and rank-reducing gates are not naively
composable.  The retained records are
[`restricted_bcw_circuit_search_xxs_rank_hybrid_w24.json`](../artifacts/generated-results/restricted_bcw_circuit_search_xxs_rank_hybrid_w24.json),
[`restricted_bcw_circuit_search_xxs_v2r_w16.json`](../artifacts/generated-results/restricted_bcw_circuit_search_xxs_v2r_w16.json),
and
[`restricted_bcw_circuit_search_xxs_y2vb_w16.json`](../artifacts/generated-results/restricted_bcw_circuit_search_xxs_y2vb_w16.json).
New atoms should now be generated from the equations `K^TAK=0`, rather than
guessed from terminal coordinates.

## 5. Attack D: separate the quartic dimension problem

Every cotangent lift doubles the cubic dimension.  Therefore beating the
current 42-variable quartic by this route is exactly the problem of finding
a noninjective cubic-homogeneous source in dimension at most 20.  Hessian
rank and ambient dimension should no longer share one scalar objective:

- the dimension beam must minimize the fixed-point quotient dimension
  before evaluating the cotangent lift;
- the rank beam must minimize `(rank JH, kernel excess)`;
- a direct non-cotangent quartic search is justified only if it uses the
  symmetric Hessian tensor and nilpotency equations from the start.

For direct quartics, the next rigorous lower-bound step is dimension six.
The tractable formulation is the coefficient ideal generated by
`tr((Hess P)^j)` for `1<=j<=6`, stratified by generic Hessian rank and
linear kernel dimension.  The output must be a classification or a
certificate that the relevant saturation is empty; a random HN scan cannot
raise the lower bound.

## 6. Promotion rules

Each attack has three levels:

1. modular profiles select a candidate or an algebraic stratum;
2. exact characteristic-zero computation proves polynomial rank and power
   identities;
3. an independent sparse replay verifies the collision or the positive
   classification certificate.

Only level three changes an endpoint.  Bounded failures are retained because
they identify incompatible circuit factorizations, but they are never
reported as lower bounds.

The immediate order of work is:

1. synthesize circuit atoms from `K^TAK=0`;
2. build the polarized index-three inverse-tree reducer;
3. begin the rank-three kernel-bundle classification at Pluecker degrees
   zero and one;
4. run the dimension-six direct-HN saturation independently of the
   cotangent program.
