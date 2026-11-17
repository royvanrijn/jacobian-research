# F2 `k=1` `E_8` icosahedral monodromy and peripheral atlas

> **Status.**  Exact degree-six permutation theorem and exhaustive `A_5`
> coset atlas.  For the affine monomial cusp `P^5-Q^3=0`, every transitive
> six-sheet action whose geometric meridian has cycle type `2+2+1+1` is
> conjugate to one explicit action with image `A_5`.  Across every transitive
> coset action of this `A_5` quotient, the fixed-sheet degrees in the F2 range
> are exactly `6,10,15,30`.  In every case the preferred longitude equals the
> meridian, so each ramified transposition is a separate `(e,f)=(2,1)`
> source-boundary row.  This is not a classification of all finite quotients
> of the cusp group and does not construct a Keller cover.  It instead closes
> the complete icosahedral subproblem and gives a uniform logarithmic-Chern
> ledger for its remaining candidates.

The enumeration is replayed without external CAS dependencies by
[`verify_f2_affine_k1_e8_monodromy.py`](../scripts/verify_f2_affine_k1_e8_monodromy.py).

## 1. The cusp group and its meridian

At the concentrated endpoint of the normalized `k=1` chart,

\[
 p=t^3,\qquad q=t^5,
\]

the affine image is the weighted-homogeneous cusp

\[
 C_{E_8}:\quad P^5-Q^3=0.                       \tag{1.1}
\]

Its affine complement has the `(3,5)` torus-knot group

\[
 G_{E_8}=\langle a,b\mid a^3=b^5\rangle.       \tag{1.2}
\]

In the abelianization, `[a]=5` and `[b]=3`.  Hence

\[
 m=a^{-1}b^2                                      \tag{1.3}
\]

has linking number one and is a geometric meridian, up to reversing the
common orientation.  Reversal does not change its permutation cycle type.

## 2. Exact degree-six classification

Enumerate all pairs `(A,B)` in `S_6^2` satisfying

\[
 A^3=B^5,                                       \tag{2.1}
\]

retain the transitive pairs, and impose

\[
 \operatorname{type}(A^{-1}B^2)=2+2+1+1.       \tag{2.2}
\]

There are exactly `720` labeled solutions.  They form one simultaneous
conjugacy class.  One representative is

\[
\boxed{\begin{aligned}
 A&=(1\ 2\ 3)(4\ 5\ 6),\\
 B&=(2\ 4\ 5\ 3\ 6),\\
 M=A^{-1}B^2&=(1\ 3)(2\ 4).
\end{aligned}}                                                   \tag{2.3}
\]

Every solution has

\[
 \operatorname{type}(A)=3+3,\qquad
 \operatorname{type}(B)=5+1,
 \qquad A^3=B^5=1.                              \tag{2.4}
\]

For (2.3), `M^2=A^3=(MA)^5=1`, and the exact generated subgroup has order
`60`.  It is therefore the `(2,3,5)` icosahedral group `A_5`, in its
exceptional transitive action on six points.  The equality between the
`720`-element conjugacy orbit of (2.3) and the complete solution set proves
uniqueness up to sheet relabeling.

### 2.1 The peripheral action separates the two rows

Let

\[
 z=a^3=b^5.
\]

For the `(3,5)` torus knot, the preferred longitude is

\[
 \ell=z,m^{-15}.                                \tag{2.5}
\]

Every solution in (2.3)--(2.4) has `z=1` and `M^2=1`.  Hence

\[
 L=M,                                             \tag{2.6}
\]

and the peripheral subgroup has sheet orbits

\[
 \boxed{\{1,3\},\quad\{2,4\},\quad\{5\},\quad\{6\}.} \tag{2.7}
\]

Connected components over a boundary torus are exactly peripheral orbits.
For a ramified orbit its size is `e*f`, while the meridian-cycle length is
`e`.  Each two-element orbit in (2.7) therefore has `(e,f)=(2,1)`.  In
particular the longitude does not exchange the two transposition cycles, so
they cannot be the two residue embeddings of one `(2,2)` divisor.

## 3. Consequence for `(75,125)`

The two fixed points of `M` are exactly the affine-sheet remainder required
by Keller etaleness at geometric degree six.  The other four sheets are two
transpositions.  By (2.7), every E8 survivor therefore has two distinct new
source-boundary divisors

\[
 (e_1,f_1)=(e_2,f_2)=(2,1).                     \tag{3.1}
\]

This raises the stratum-conditional source-component floors to `29/50`, as
on the cyclic-complement strata, although here both source rows lie over one
noncyclic target component.

Write `E_i^2=-n_i`, put `N=n_1+n_2`, let the common target curve follow `b`
smooth carrier centers, and let the completed source add `s_X` further
smooth-boundary blowups.  Each row contributes
`2*(b-7)-2*n_i` to the cyclic divisorial ledger.  The total squarefree
finite length is therefore

\[
 \ell(Z)=18+2N-4b-\frac{s_X}{2}.                 \tag{3.2}
\]

The multiplicity-three cusp lower bound costs at least two on each
residue-degree-one row.  After subtracting those four units, one obtains

\[
 \boxed{
 \ell(Z_{\rm rest})=14+2N-4b-\frac{s_X}{2}\ge0.} \tag{3.3}
\]

Equivalently,

\[
 \boxed{28+4N-8b-s_X\ge0.}                      \tag{3.4}
\]

Since both divisors are exceptional, the minimal numerical choice is
`N=2`.  With `s_X=0`, (3.4) forces `b<=4`.  The separate carrier-jet theorem
supplies the complementary endpoint
gate: **if** the curve has maximal truncated contact `b=8`, the actual seven
fixed centers must lie on a prime codimension-four complete intersection.
Equation (3.4) then forces

\[
 \boxed{4N\ge36+s_X.}                            \tag{3.5}
\]

Thus maximal contact is impossible for the minimal two-row configuration;
it must leave the carrier by `b<=4`.

There is an exact reason that the `n=1` inequality does not yet contradict
the carrier equations.  Eliminate the free E8 transports successively by
(6.17) of the carrier-jet theorem.  For fixed nonzero leading scales, the
conditions that **some** transported E8 target have contact at least `b`
with the fixed carrier centers are

\[
\begin{array}{c|c|c}
b\text{ threshold}&\text{raw-center conditions}&\text{codimension}\\ \hline
b\ge2&\zeta_1=0&1\\
b\ge3&\zeta_1=0&1\\
b\ge4&\zeta_1=0&1\\
b\ge5&\zeta_1=\widehat E_4=0&2\\
b\ge6&\zeta_1=\widehat E_4=0&2\\
b\ge7&\zeta_1=\widehat E_4=\widehat E_6=0&3\\
b\ge8&\zeta_1=\widehat E_4=\widehat E_6=\widehat E_7=0&4.
\end{array}                                                    \tag{3.6}
\]

The repeated rows occur because `eta`, `mu`, and `nu` recover the matches at
orders `2`, `3`, and `5`.  The first normalization-independent compatibility
equation appears only at `b>=5`.  The minimal two-row Chern gate forces
`b<=4`, exactly before that equation.  Therefore the present carrier and
Chern constraints occupy complementary ranges and do not exclude the
minimal E8 packet.  What remains is to expose the centers only to the order
actually reached, fix the transports by independent global data, or bound
the total negativity `N`.  The row decomposition itself is no longer
ambiguous.  Neither the permutation classification nor (3.3) alone excludes
`(75,125)`.

## 4. The complete `A_5` coset atlas

The degree-six action is only one transitive action of its `A_5` image.
Enumerating all `59` subgroups of `A_5` and quotienting by conjugacy gives
the nine classical subgroup classes.  Their coset actions have the following
meridian data.  Here `r` is the number of ramified transpositions and `u` is
the number of fixed sheets.

\[
\begin{array}{c|c|c|c|c}
H&d=[A_5:H]&\operatorname{type}(M)&u&r\\ \hline
A_5&1&1&1&0\\
A_4&5&2^2 1&1&2\\
D_5&6&2^2 1^2&2&2\\
S_3&10&2^4 1^2&2&4\\
C_5&12&2^6&0&6\\
V_4&15&2^6 1^3&3&6\\
C_3&20&2^{10}&0&10\\
C_2&30&2^{14}1^2&2&14\\
1&60&2^{30}&0&30
\end{array}                                                    \tag{4.1}
\]

The positive affine-sheet remainder requires `u>0`.  Thus, within this
icosahedral quotient and in the F2 degree range, the complete list of
one-component fixed-sheet escapes is

\[
 \boxed{d\in\{6,10,15,30\}.}                  \tag{4.2}
\]

This is already a useful negative result: complement topology cannot force
the degree-six floor.  The zero-fixed-sheet actions in degrees `12,20,60`
cannot by themselves supply the affine remainder, although this statement
does not exclude a completion with further affine branch components.

For every row of (4.1), the central element has trivial image and

\[
 L=zM^{-15}=M.                                  \tag{4.3}
\]

Consequently every ramified peripheral orbit is exactly one meridian
two-cycle.  All `r` ramified rows are distinct and have `(e,f)=(2,1)`; none
can be merged into a residue-degree-greater-than-one row.  Conditional on the
same component filtration used in Section 3, the source-component floors are
therefore `27+r` in the squarefree case and `48+r` in the double-row case.
For the fixed-sheet degrees these are

\[
\begin{array}{c|c|c|c}
d&r&\text{squarefree floor}&\text{double-row floor}\\ \hline
6&2&29&\text{--}\\
10&4&31&\text{--}\\
15&6&33&54\\
30&14&41&62
\end{array}.                                                   \tag{4.4}
\]

There is also a uniform Chern gate.  Let the `r` rows have
`E_i^2=-n_i`, set `N=\sum_i n_i`, let all rows follow the same `b` carrier
centers, and retain the `s_X` additional smooth-boundary blowups.  Subtracting
the multiplicity-three cusp lower bound `2r` from the finite point budget
gives

\[
 \boxed{
  2\ell_{\mathrm{rest}}^{\mathrm{sf}}
   =7d-62+4N-4r(b-6)-s_X\ge0,}                 \tag{4.5}
\]

and, on the double row,

\[
 \boxed{
  2\ell_{\mathrm{rest}}^{\mathrm{dbl}}
   =7d-67+4N-4r(b-6)-s_X\ge0.}                 \tag{4.6}
\]

The parity conditions are correspondingly

\[
 s_X\equiv d\pmod 2\quad\text{(squarefree)},
 \qquad
 s_X\equiv d+1\pmod 2\quad\text{(double row)}.               \tag{4.7}
\]

At the minimally negative choice `N=r`, maximal contact `b=8`, and `s_X=0`,
the doubled squarefree remainders for `d=6,10,15,30` are respectively

\[
 -28,\ -8,\ 19,\ 92,                            \tag{4.8}
\]

while the double-row values at `d=15,30` are `14,87`.  Hence minimal
maximal-contact icosahedral packets are excluded in squarefree degrees six
and ten, but extra negativity repairs the inequality and the higher degrees
survive it outright.  Equations (4.5)--(4.8) do not exclude the atlas; they
identify the missing global input precisely: a geometric upper bound on
`N`, a restriction on the attained carrier contact `b`, or compatibility
with the other compiled boundary packets.

The subsequent
[`simple-inertia orbifold atlas`](F2_AFFINE_K1_E8_ORBIFOLD_ATLAS.md)
removes the assumption that the central image is trivial.  The universal
`M^2=1` quotient has order `240`; its order-four center glues transposition
cycles into `(2,2)` and `(2,4)` rows and enlarges the fixed-sheet spectrum
to ten F2 degrees.  Thus this section is exactly the central-trivial
subatlas of that complete simple-inertia theorem.

<!-- status-consumer: PF2K1E8O1 4251750ed4e43c89 -->

## Reproduction

```bash
.venv/bin/python scripts/verify_f2_affine_k1_e8_monodromy.py
```

The command enumerates all `S_6` images of (1.2), checks transitivity and
the meridian condition, proves that the `720` retained labeled actions are
one conjugacy orbit, and verifies that the generated image has order `60`.
It then enumerates all `59` subgroups of that image, obtains the nine
conjugacy classes and every entry of (4.1), evaluates the preferred longitude,
and checks the numerical specializations of (4.5)--(4.8).
