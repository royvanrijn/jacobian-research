# F2 carrier logarithmic node profile

> **Status.**  This note proves the exact carrier-local logarithmic profiles
> for both live F2 `(75,125)` carrier rows.  The two squarefree simple-root
> spectators each require one blowup.  The double-carrier fivefold point
> requires four carrier-centered blowups and one additional fan-alignment
> blowup.  Every node created there is tame log-etale.  The principal arms
> require six further fan-alignment components each before their terminal
> divisors, and every corresponding node has logarithmic determinant `3`.
> Combining these refinements with `PF2LNP1` strengthens the source-boundary
> lower bounds from `19/31` to `27/48` components.  The result does not cover
> the upstream carrier-extraction chain, the outgoing terminal tail, the
> affine purity row, or any uncompiled global resolution center, and it does
> not exclude `(75,125)`.

The upstream item left outside this theorem is subsequently compiled in
[`F2_UPSTREAM_CARRIER_EXTRACTION_PROFILE.md`](F2_UPSTREAM_CARRIER_EXTRACTION_PROFILE.md).
The outgoing item is subsequently closed as a unimodular log-etale fan map in
[`F2_OUTGOING_TERMINAL_TAIL.md`](F2_OUTGOING_TERMINAL_TAIL.md).
Its carrier-zero ladder is unimodular, while its extraction-root node has the
nonzero length-`54` degree-one branch-matching quotient.
The later
[`affine-purity frontier`](F2_AFFINE_PURITY_FRONTIER.md) proves that a
different component with `e>1` over an affine target curve must be added,
raising the global component floors to `28/49` without changing this note's
carrier-local count.
<!-- status-consumer: PF2UCE1 7f15bc756cc73fff -->

The rational orders, common regular fans, exponent matrices, intersection
forms, and canonical integrality are checked by
[`verify_f2_carrier_log_node_profiles.py`](../scripts/verify_f2_carrier_log_node_profiles.py).
The reusable common-refinement routine is in
[`log_node_profiles.py`](../jcsearch/log_node_profiles.py).

## 1. Carrier target coordinates

Use the certified carrier coordinates

\[
 q=y,\qquad v=xy^5,qquad
 P=q^{-15}c(v)^3+O(q^{-14}),qquad
 -Q=q^{-25}\frac95c(v)^5+O(q^{-24}),              \tag{1.1}
\]

where `c(v)=v(v-1)^2R(v)`.  After the seven removable target shears, put

\[
 \pi=\frac{P^3}{(-Q)^2}=q^5U(v)+O(q^6),qquad
 w=q^{36}H(v)+O(q^{37}),                           \tag{1.2}
\]

with

\[
 U=\frac{25}{81c},qquad
 H=\frac{d(v)N(v)}{c(v)^8}.                        \tag{1.3}
\]

The target ray is `(5,36)`.  Its adjacent regular rays are `(1,7)` and
`(4,29)`, and its residue coordinate is

\[
 \zeta=\frac{w^5}{\pi^{36}}.                       \tag{1.4}
\]

At the node `zeta=0`, regular toric coordinates are

\[
 \Pi_0=\frac{\pi^{29}}{w^4}
      =q\frac{U^{29}}{H^4}+O(q^2),qquad
 \Xi_0=\zeta.                                      \tag{1.5}
\]

At the node `zeta=infinity`, use

\[
 \Pi_\infty=\frac{w}{\pi^7}
      =q\frac{H}{U^7}+O(q^2),qquad
 \Xi_\infty=\zeta^{-1}.                            \tag{1.6}
\]

These formulas supply the transverse column that residue inertia alone does
not see.

## 2. The order compiler

Write

\[
 u=\operatorname{ord}(U),\qquad h=\operatorname{ord}(H).
\]

Equations (1.5)--(1.6) give the exact pairs

\[
 \begin{array}{c|c}
 \text{target node}&
   (\text{transverse coefficient order},\text{residue index})\\ \hline
 \zeta=0&(29u-4h,\;5h-36u)\\
 \zeta=\infty&(h-7u,\;36u-5h).
 \end{array}                                        \tag{2.1}
\]

At a simple root of `c`, one has `(u,h)=(-1,-7)`, hence

\[
 (p,e)=(-1,1)\quad\text{over }\zeta=0.             \tag{2.2}
\]

At a double root of `c`, `(u,h)=(-2,-15)`, hence

\[
 (p,e)=(-1,3)\quad\text{over }\zeta=\infty.        \tag{2.3}
\]

At the extra simple root of `N` in the double row, `(u,h)=(0,1)`, so

\[
 (p,e)=(-4,5)\quad\text{over }\zeta=0.             \tag{2.4}
\]

Finally `U` and `H` have infinity orders `5` and `36`.  The normal
coefficient in (1.5) has order `+1`, while the residue map has contact order
three with its finite limiting value.  Keller boundary support therefore
gives the same smooth-endpoint model as in `PF2LNP1`:

\[
 \mathcal T_f^{\log}\simeq R/(t^3),                \tag{2.5}
\]

with zero normalization defect.

## 3. Squarefree carrier

For

\[
 R(v)=\frac{v^2-3v+3}{25},qquad
 \zeta(v)=1+\frac1{(v-1)^3},                       \tag{3.1}
\]

the marked carrier points are:

| source point | target | `(p,e)` | boundary role |
| --- | --- | --- | --- |
| `v=0` | node `zeta=0` | `(-1,1)` | existing carrier endpoint |
| either root of `R` | node `zeta=0` | `(-1,1)` | new simple spectator |
| `v=1` | node `zeta=infinity` | `(-1,3)` | principal arm |
| `v=infinity` | smooth point `zeta=1` | `(+1,3)` | existing carrier endpoint |

Each simple spectator needs one blowup.  At its carrier node the exponent
matrix is `diag(1,1)`, so the log cokernel is zero.  The fact that the
residue inertia is trivial does not remove the boundary branch: the
transverse pole in (2.2) forces it.

## 4. Double carrier

Let `rho^2-3rho+1=0`, `alpha=3(rho+1)/5`, and

\[
 \zeta(v)=
 \frac{v(v-\alpha)^5}{(v-1)^3(v-\rho)^3}.          \tag{4.1}
\]

The marked points are:

| source point | target | `(p,e)` | boundary role |
| --- | --- | --- | --- |
| `v=0` | node `zeta=0` | `(-1,1)` | existing endpoint |
| `v=alpha` | node `zeta=0` | `(-4,5)` | new carrier attachment |
| `v=1,rho` | node `zeta=infinity` | `(-1,3)` | two principal arms |
| `v=infinity` | smooth third branch value | `(+1,3)` | existing endpoint |

At `v=alpha`, four successive carrier-centered blowups give source rays

\[
 (1,0),(4,1),(3,1),(2,1),(1,1).                   \tag{4.2}
\]

The exponent map is

\[
 M_5=\begin{pmatrix}1&-4\\0&5\end{pmatrix}.       \tag{4.3}
\]

The cone between `(2,1)` and `(1,1)` crosses the next target ray.  Its
primitive inverse image is `(3,2)`, so one additional blowup is necessary.
The final fan is

\[
 (1,0),(4,1),(3,1),(2,1),(3,2),(1,1).             \tag{4.4}
\]

Every adjacent logarithmic exponent matrix has determinant `5`.  Thus all
five new nodes are log-etale in characteristic zero.

## 5. Principal-arm common refinement

At `v=1` and, in the double row, at `v=rho`, the local exponent map is

\[
 M_3=\begin{pmatrix}1&-1\\0&3\end{pmatrix}.       \tag{5.1}
\]

Before map alignment, the source rays from the carrier through the terminal
divisor are

\[
 (1,0),(1,1),(1,2),(3,7),(5,12).                  \tag{5.2}
\]

In local coordinates at the carrier target node, the relevant target fan is

\[
 (1,0),(0,1),(-1,6),(-2,11),\ldots,(-7,36).       \tag{5.3}
\]

Pulling those rays back and minimally regularizing gives

\[
\begin{split}
 &(1,0),(1,1),(1,2),(5,11),(4,9),(7,16),(3,7),\\
 &(11,26),(8,19),(13,31),(5,12).                  \tag{5.4}
\end{split}
\]

Thus six components were missing between the carrier and each terminal
divisor.  Every adjacent source cone is regular, maps into one regular target
cone, and has logarithmic exponent determinant

\[
 \det(M_3)=3.                                      \tag{5.5}
\]

All these nodes are log-etale.  This is a common-fan calculation, not an
assumption that the source valuation fan already resolves the map.

## 6. Strengthened source graphs

Insert the common refinements into the `PF2LNP1` principal graphs.

For the squarefree row, six arm components and two simple spectators change

\[
 (19\text{ components},6\text{ leaves})
 \longmapsto(27,8).                                \tag{6.1}
\]

The carrier has weight `-4` and valency five; the terminal has weight `-9`.
The intersection determinant is `1`, the inertia is `(1,26,0)`, and the
canonical vector is integral.

For the double row, the two arm refinements contribute twelve components
and the `alpha` chain contributes five:

\[
 (31\text{ components},10\text{ leaves})
 \longmapsto(48,11).                               \tag{6.2}
\]

The carrier has weight `-7` and valency five; both terminals have weight
`-9`.  The intersection determinant is `-1`, the inertia is `(1,47,0)`, and
the canonical vector is integral.

These are carrier-complete **lower bounds**, not complete compactifications.
Further global centers may only enlarge them.

## 7. Consequence for the universal defect programme

Neither the terminal Belyi packet nor its certified carrier-local completion
supplies a nonzero normalization-defect class.  Residue indices `1`, `3`,
and `5` become invertible log exponent determinants after the necessary
common fan refinement.

After the subsequent outgoing-tail and affine-purity frontier the remaining
candidate support is restricted to:

1. the actual target nonproperness curve, its pullback factorization, and the
   proximity chain of the purity-forced new component;
2. uncompiled global resolution centers; and
3. the global localized-second-Chern remainder and possible cancellation of
   the extraction-root length-`54` matching class.

<!-- status-consumer: PF2OTT1 af25012e34020e11 -->

<!-- status-consumer: PF2APF1 192055eb737d3140 -->

The subsequent extraction theorem supplies such a class, but in degree-one
local cohomology.  The central logical gap is now global: prove that Keller
geometry cannot absorb or cancel that exact class.  More coefficient layers
on either tame packet cannot change this target.
