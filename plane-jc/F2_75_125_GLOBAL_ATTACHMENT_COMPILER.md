# Global attachment compiler for F2 `(75,125)`

## Result and claim boundary

The exact checker
[`cas/verify_f2_75_125_global_attachment.py`](cas/verify_f2_75_125_global_attachment.py)
now carries the terminal F2 row back through the Kummer translation to the
original affine source coordinates.  This resolves the previously missing
principal proximity chain.

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

The three interior node preimages then require three further boundary
branches at distinct points of each terminal divisor.  Consequently the
minimal principal source boundary has:

| F2 case | principal arms | components | leaves | carrier valency |
| --- | ---: | ---: | ---: | ---: |
| squarefree | 1 | 16 | 6 | 3 |
| double root | 2 | 25 | 10 | 4 |

Every terminal divisor has valency five.  The compiled source matrices
are unimodular, have Hodge inertia `(1,n-1,0)`, and give integral canonical
coefficient vectors by adjunction.

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
inertia, supply the purity-forced affine row, or exclude `(75,125)`.  The
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

## 4. Three interior attachments per arm

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

Blowing up the three distinct interior points lowers the terminal
self-intersection from `-1` to `-4` and raises its boundary valency from two
to five.

For one arm, the resulting minimal principal source graph has 16 components,
6 leaves, determinant `-1`, and canonical vector

```text
(-1,-2,-3,-4,-5,-2,-3,-4,-3,-6,-9,-4,-2,-8,-8,-8).
```

For two arms at `v=1,rho`, the arms share the seven-component carrier chain.
The resulting graph has 25 components, 10 leaves, determinant `1`, and two
valency-five terminal vertices.

## 5. Complete minimal target extraction

Near the `Q`-dominant target point set

\[
a=(-Q)^{-1},\qquad b=P/(-Q).
\]

Here `a=0` is the target line at infinity and `b=0` is nonboundary.  The
minimal regular fan containing `(5,2)` is

\[
(1,0),(3,1),(5,2),(2,1),(1,1),(0,1). \tag{7}
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
A separate row with `e>1` over an affine nonproperness curve and its positive
affine companion are still required.

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
the terminal degree-six Belyi map after a linear source change.  Thus the
remaining inputs are:

1. impose the three exact carrier parameter points on the complete lower
   Laurent coefficient system;
2. construct the purity-forced affine ramification row and affine companion;
3. determine the common global field degree and all remaining pullback rows;
4. supply the corresponding global branch cycles and run candidate mode.

Until one of those declarations fails a hard compiled gate, `(75,125)`
remains unexcluded.
<!-- status-consumer: PF2CW1 a7774b0fa736b64c -->
