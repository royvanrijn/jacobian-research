# Global attachment compiler for F2 `(75,125)`

<!-- status-consumer: PF2LNP1 e4f0f231bf7494d5 -->

<!-- status-consumer: PF2CLP1 41625dd5d3f8f898 -->

<!-- status-consumer: PF2UCE1 7f15bc756cc73fff -->

<!-- status-consumer: LCBBC1 b3eb4679f781c55f -->

<!-- status-consumer: LCCT1 2fc6ecea7a7c8b49 -->

<!-- status-consumer: LKGD1 8a357250b5005186 -->

<!-- status-consumer: LTKT1 32ac27318f16c20c -->

<!-- status-consumer: LCHB1 176bf85520516fa6 -->

## Result and claim boundary

The exact checker
[`cas/verify_f2_75_125_global_attachment.py`](cas/verify_f2_75_125_global_attachment.py)
now carries the terminal F2 row back through the Kummer translation to the
original affine source coordinates.  This resolves the previously missing
principal proximity chain.  The terminal attachment part was subsequently
corrected by the exact [`log-node profile theorem`](F2_LOG_NODE_PROFILE.md):
residue incidence alone missed a transverse pole of order two at each
interior node preimage.

The important coordinate point is that the terminal variable `y` is

\[
y_{\rm tr}=y_{\rm old}-X^{-1},
\]

not the original affine coordinate.  Therefore the Laurent normal
`(12,-17)` does **not** become a monomial original-lattice ray
`(-25,17)`.  Instead it gives the nonmonomial original valuation

\[
\boxed{
\nu(x)=-25,\qquad
\nu(y_{\rm old})=5,\qquad
\nu(xy_{\rm old}^5-c)=12,
} \tag{1}
\]

where `c=1` for the fixed principal chain and `c=rho` for the second
double-root chain.

Equation (1) has an exact two-stage resolution:

1. six blowups extract the monomial carrier ray `(-5,1)`;
2. at the marked carrier point `v=xy_old^5=c`, six further blowups extract
   the local ray `(5,12)` and its terminal divisor.

The three interior node preimages require two successive boundary blowups at
each distinct point of every terminal divisor.  Consequently the corrected
principal source lower bound is:

| F2 case | principal arms | components | leaves | carrier valency |
| --- | ---: | ---: | ---: | ---: |
| squarefree | 1 | 19 | 6 | 3 |
| double root | 2 | 31 | 10 | 4 |

Every terminal divisor has valency five and self-intersection `-7`.  The
compiled source matrices
are unimodular, have Hodge inertia `(1,n-1,0)`, and give integral canonical
coefficient vectors by adjunction.

These `19/31` graphs are the pinned compiler's terminal-resolved principal
skeletons, not the latest carrier-complete lower bounds.  The subsequent
[`carrier log-node profile theorem`](F2_CARRIER_LOG_NODE_PROFILE.md) inserts
six common-fan components on every principal arm, the two simple squarefree
spectator branches, and the five-component double-row attachment.  It gives:

| F2 case | components | leaves | carrier valency | terminal weights |
| --- | ---: | ---: | ---: | --- |
| squarefree | 27 | 8 | 5 | `(-9)` |
| double root | 48 | 11 | 5 | `(-9,-9)` |

Both refined intersection forms remain unimodular with integral canonical
vectors.  All marked new nodes are log-etale.  The subsequent upstream
extraction theorem makes the carrier-zero ladder unimodular and finds the
distinct extraction-root cokernel `R/(W^3U^18)`, with length-`54` branchwise
matching quotient.  The blowup-conservation theorem replaces that raw length
by the stable Cartier charge `D_root^2/2=27`.  The kernel-line theorem further
identifies the exact cyclic contribution as `deg(K_root)+27`; since the packet
is contracted, the Gauss-degree theorem makes this `27-e_root<=27` for a
nonnegative kernel-direction degree.  The tangential-coordinate theorem
computes `e_root=0` from full divisibility of `d(f^*z)`, so the cyclic root
contribution is exactly `27`.  The outgoing terminal-tail theorem maps the
remaining source rays to the existing target fan unimodularly, so that tail
is also closed with zero defect.

The subsequent complete-chain theorem shows that the missing affine-row
self-intersections are not independent repair parameters.  For a sole
rational one-puncture `k=1` affine component, the full Cartier determinant
cycle leaves point budget exactly `u-1`.  Every simple-inertia E8 action has
cusp lower `2R>u-1`, so it requires a negative normalization/`Fitt_1` class
at the still-unresolved global attachment.

<!-- status-consumer: PF2K1CB1 5cc386dba344a867 -->

For the unique degree-six cubic-inertia equality row, the contracted-divisor
Smith theorem now replaces that unspecified attachment by one exact packet.
If `T^2=-n`, `v` is its boundary valency, and
`I=(D_log-2T).T`, its remaining isolated point budget is

\[
 P_{\rm other}=2I-v-3n\ge0.
\]

When `I=3`, only a `(-1)` vertex of valency one, two, or three survives,
leaving respectively two, one, or zero point units.  Thus every future
candidate declaration for the cubic row must include `(T,I,v,n)`; the gate
is independent of carrier contact and lower Laurent coefficients.

<!-- status-consumer: LCDSC1 07dcd994b4faf092 -->

The global logarithmic `ch_2` theorem now supplies an exact ledger on the
same partial model.  Combining the `(5,2)` and `(5,36)` target clusters gives
`L_Y^2=-5`; adjunction on the refined squarefree/double source graphs gives
`L_X^2=-6/-11`.  Hence the global budgets are `(7*d-8)/2` and
`(7*d-13)/2`, and subtracting the root class gives `(7*d-62)/2` and
`(7*d-67)/2`.  Their values `-10` and `17/2` at the degree floors are virtual
remainders, not exclusions: the affine purity row changes the common model
and contributes a new divisorial module, and no effective finite-length
filtration has yet been constructed.

The target side remains completely explicit: four blowups extract `(5,2)`,
with boundary weights

\[
(-2,-2,-1,-3,-2),
\]

intersection determinant `1`, and canonical coefficients

\[
(-3,-6,-9,-4,-2).
\]

This is genuine global-attachment progress.  Target-valuation uniqueness also
proves that both double-root packets land on the same extracted target
divisor, so every double-root realization has geometric degree at least
twelve.  This compiler by itself does not determine the simple spectator
inertia or exclude `(75,125)`.  The later affine-purity theorem forces at
least one additional branch-boundary component, raising the source floors to
`28/49`, but does not construct its target curve or pullback factorization. The
subsequent
[`carrier Wronskian classifier`](F2_75_125_CARRIER_WRONSKIAN_CLASSIFIER.md)
does determine the spectator inertia and reduces both cofactor strata to
three exact parameter points.

The pinned output is
[`../artifacts/generated-results/jc2_f2_75_125_global_attachment.json`](../artifacts/generated-results/jc2_f2_75_125_global_attachment.json).

## 1. Returning to the original source coordinates

The selected Kummer chart is

\[
x=X^5,qquad
t=Xy_{\rm old}-1,qquad
y_{\rm tr}=y_{\rm old}-X^{-1},qquad
z=y_{\rm tr}^{-1}=\frac Xt.
\]

Thus

\[
X=tz,qquad y_{\rm tr}=z^{-1}.
\]

The terminal source normal satisfies

\[
\nu(t)=12,qquad \nu(z)=-17.
\]

It follows that

\[
\nu(X)=-5,qquad
\nu(y_{\rm tr})=17,qquad
\nu(X^{-1})=5.
\]

Since the last two orders differ, no leading cancellation is possible in
`y_old=X^{-1}+y_tr`, and hence

\[
\nu(y_{\rm old})=5,qquad \nu(x)=-25. \tag{2}
\]

Put

\[
v=xy_{\rm old}^5=(Xy_{\rm old})^5=(1+t)^5.
\]

In characteristic zero,

\[
\nu(v-1)=\nu((1+t)^5-1)=\nu(t)=12. \tag{3}
\]

The gcd of `25,5,12` is one, so (1) is already a normalized divisorial
valuation.  This also explains why monomial orders alone were insufficient:
the center polynomial `v-1` supplies the missing value coprime to five.

## 2. Six-blowup carrier

The predecessor monomial valuation is the primitive ray

\[
(-5,1).
\]

In the standard `P^2` cone from the affine divisor `y_old=0`, with ray
`(0,1)`, to the line at infinity, with ray `(-1,-1)`, one has

\[
(-5,1)=6(0,1)+5(-1,-1). \tag{4}
\]

The minimal regular subdivision containing coordinate ray `(6,5)` inserts
six exceptional rays.  In regular coordinate order the cone is

```text
(1,0), (2,1), (3,2), (4,3), (5,4), (6,5), (1,1), (0,1).
```

The first endpoint is nonboundary.  The other seven components form the
carrier boundary chain, with self-intersections

```text
(-2,-2,-2,-2,-1,-6,0).
```

The component corresponding to `(6,5)` is the carrier.  Its residue
coordinate is

\[
v=xy_{\rm old}^5.
\]

The fixed squared factor marks `v=1`.  If `R` is squarefree, its two simple
roots `rho_1,rho_2` mark two distinct spectator points on this same carrier.
If `R` has a nonzero double root `rho`, then `rho!=1` because `R(1)=1/25`,
and `v=rho` is a second principal center.

Thus the source incidence question is no longer “where are the Kummer
orbits?”  Their carrier points are known.  What remains unknown for a simple
root is the branch scale and global inertia attached at that point.

## 3. Six-blowup principal arm

At a principal center `v=c`, use

\[
q=y_{\rm old},qquad r=v-c.
\]

The terminal valuation has local orders

\[
(\nu(q),\nu(r))=(5,12). \tag{5}
\]

The minimal regular subdivision of `(5,12)` is

```text
(1,0), (1,1), (1,2), (3,7), (5,12), (2,5), (1,3), (0,1).
```

The first endpoint is the pre-existing carrier and the final endpoint is a
nonboundary tangential curve.  The six new boundary components have weights

```text
(-2,-4,-2,-1,-3,-2).
```

The terminal component is the fourth new component, corresponding to
`(5,12)`, and initially has self-intersection `-1`.

There is also an exact endpoint orientation.  At `c=1`, expansion of
`y_tr=y_old-X^-1` gives, up to a unit,

\[
s_{\rm term}=X^{17}y_{\rm tr}^5
             \sim \frac{r^5}{q^{12}}. \tag{6}

Its orders on the adjacent rays `(3,7)` and `(2,5)` are `-1` and `+1`.
Therefore

- the `(3,7)` neighbor is the endpoint `s=infinity`;
- the `(2,5)` neighbor is the endpoint `s=0`.

Kummer transfer gives the same oriented arm at the second double-root center.

## 4. Three two-blowup interior attachments per arm

The residue map

\[
h(s)=\frac{125s(s+1)^5}{(9s^2+15s+5)^3}
\]

gives five marked boundary-neighbor slots on every terminal divisor.

| source point | target locus | index | different | source graph role |
| --- | --- | ---: | ---: | --- |
| `s=0` | node `h=0` | 1 | 0 | `(2,5)` arm neighbor |
| `s=infinity` | smooth `h=125/729` | 3 | 2 | `(3,7)` arm neighbor |
| `s=-1` | node `h=0` | 5 | 4 | new interior branch |
| first denominator root | node `h=infinity` | 3 | 2 | new interior branch |
| second denominator root | node `h=infinity` | 3 | 2 | new interior branch |

The residue table alone does not determine the number of point blowups.  In
target-node coordinates, the transverse coefficient has order `-2` at each
of the three interior points.  The local form is

\[
 (\pi,\xi)=(\tau w^{-2}\cdot\text{unit},
             w^e\cdot\text{unit}),
 \qquad e=5,3,3.                                  \tag{7}
\]

One blowup leaves a pole of order one.  A second blowup at the intersection
of the strict terminal divisor and the first exceptional curve gives
`tau=w^2*u`.  Each marked point therefore contributes a two-component chain
with weights `(-1,-2)` from the terminal outwards.  The three chains lower
the terminal self-intersection from `-1` to `-7` and raise its boundary
valency from two to five.  At the terminal node the resolved logarithmic
matrix is `diag(1,e)`, so its cokernel vanishes.

For one arm, the resulting principal source lower-bound graph has 19
components, 6 leaves, determinant `1`, and canonical vector

```text
(-1,-2,-3,-4,-5,-2,-3,-4,-3,-6,-9,-4,-2,
 -16,-8,-16,-8,-16,-8).
```

For two arms at `v=1,rho`, the arms share the seven-component carrier chain.
The resulting graph has 31 components, 10 leaves, determinant `1`, and two
valency-five terminal vertices.

## 5. Complete minimal target extraction

Near the `Q`-dominant target point set

\[
a=(-Q)^{-1},\qquad b=P/(-Q).
\]

Here `a=0` is the target line at infinity and `b=0` is nonboundary.  The
minimal regular fan containing `(5,2)` is

\[
(1,0),(3,1),(5,2),(2,1),(1,1),(0,1). \tag{8}
\]

Four blowups give the boundary weights and canonical coefficients stated
above.  For `eta=a^2/b^5`, its orders on
`(3,1),(5,2),(2,1)` are `(1,0,-1)`.  Since `h=eta^{-1}`:

- `h=infinity` is the node `(3,1)|(5,2)`;
- `h=0` is the node `(5,2)|(2,1)`.

There is only one target divisor for all principal packets.  In the global
target coordinates

\[
\operatorname{ord}(a,b)=(5,2),
\]

so `pi=b^3/a` has order one, while `eta=a^2/b^5` has order zero and
nonconstant pullback `eta^{-1}=h(s)`.  Every principal source divisor is
therefore centered at the generic point of the unique divisor `(5,2)`.  Two
double-root packets cannot land on distinct target divisors.

## 6. The normalization cases

Let `d` be the global geometric degree and let each `rho_i` denote all still
unclassified contributions over one target divisor.

- Squarefree:

  \[
  d=6+\rho_T,\qquad \rho_T\ge0.
  \]

- Double root:

  \[
  d=12+\rho_T,\qquad \rho_T\ge0.
  \]

The formerly retained distinct-target alternative is excluded by the target
valuation argument above.  Thus the live global split has only the squarefree
row and the same-target double row.

The principal rows have `e=1`, so none supplies the purity branch component.
A separate component with `e>1` over an affine nonproperness curve and its
positive affine companion are required.  The affine-purity frontier proves
that this component is new and gives the exact finite generic bounds, but its
target equation and the factorization of its pullback remain unknown.  The
subsequent target-curve theorem reduces the target equation to 24 singular
normalization charts `(3k,5k)`, `1<=k<=24`, each on a nonunit
divided-difference collision/critical locus.
On the first chart, the collision equations reduce further to one quartic;
its generic target conductor consists of four ordinary affine nodes and the
fixed `(2,5)` infinity cusp.
The puncture theorem places this chart transversely on `(5,2)`.  A direct
comparison with a principal terminal packet is numerically compatible only
at `lambda=125/729,e=3`, but the terminal neighborhood is already a resolved
morphism and has no extraction slot for a divisor dominating the affine
curve.  Candidate mode must locate that divisor outside the certified
terminal neighborhood and must not copy the terminal index into it.

## 7. Candidate mode

The pinned audit and its gate regressions run with

```bash
.venv/bin/python plane-jc/cas/verify_f2_75_125_global_attachment.py
.venv/bin/python plane-jc/cas/test_f2_75_125_global_attachment.py
```

A proposed completion can be checked with

```bash
.venv/bin/python plane-jc/cas/verify_f2_75_125_global_attachment.py \
  --candidate proposed_completion.json
```

Candidate mode verifies:

- complete smooth geometric SNC source data;
- the original orders `(-25,5,12)`, carrier center, Kummer identification,
  and proximity chain for every terminal component;
- the five distinct marked neighbors per terminal component;
- connected-tree incidence, unimodularity, Hodge signature, and integral
  adjunction;
- exhaustive finite-flat target ledgers and principal packet placement;
- the purity row and affine companion;
- both squarefree spectator profiles, when applicable; and
- one product-one meridian system per target ledger, including every local
  terminal `A_6` triple and declared connectedness/transitivity.

The possible results are `excluded_candidate`, `incomplete`, or
`passes_compiled_necessary_gates_not_an_existence_proof`.  Passing remains a
necessary-condition audit, not construction of a Keller map.

## 8. Subsequent carrier classification and what remains

The carrier Wronskian classifier extracts the deeper target ray `(5,36)`.
It forces

```text
squarefree: R(v)=(v^2-3v+3)/25, carrier (e,f)=(1,3);
double:     rho^2-3*rho+1=0,   carrier (e,f)=(1,6).
```

In the squarefree row the two simple `R` points are unramified points of the
cyclic cubic carrier map `1+1/(v-1)^3`.  In the double row the carrier map is
the terminal degree-six Belyi map after a linear source change.  The carrier
log-node theorem resolves their marked local boundary points but does not add
the unresolved global rows.  Thus the remaining inputs are:

1. impose the three exact carrier parameter points on the complete lower
   Laurent coefficient system;
2. retain the exact cyclic root contribution `27`, then compile noncyclic
   attachments and all possible cancellation centers;
3. solve or further stratify the 24 affine target collision charts, form the
   implicit equation, and factor its pullback into every boundary and affine
   row;
4. determine the common global field degree and global localized-Chern row;
5. supply the corresponding global branch cycles and run candidate mode.

Until one of those declarations fails a hard compiled gate, `(75,125)`
remains unexcluded.

<!-- status-consumer: PF2OTT1 af25012e34020e11 -->
<!-- status-consumer: PF2APF1 192055eb737d3140 -->
<!-- status-consumer: PF2ATC1 9ab722c45c586b73 -->
<!-- status-consumer: PF2K1C1 358a6ba820e8b2f1 -->
<!-- status-consumer: PF2PPA1 b24c4d80c2f8230e -->
<!-- status-consumer: PF2CW1 a7774b0fa736b64c -->
