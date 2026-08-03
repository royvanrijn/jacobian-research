# F2 `j=1` derivation audit for `(75,125)`

> **Status: forced skeleton, not a normal-form certificate.**  The published
> complete-chain theorem determines two consecutive edges and one Puiseux
> translation.  The companion modified-chart audit additionally derives
> `gamma=2`, the `d=3` chart, and the exact linear image of its polynomial
> projection.  That projection has a top-band unit ideal, but the chain still
> does not determine the negative-tail correction needed for a polynomial
> modified pair.  Consequently this note does not claim to eliminate
> `(75,125)`, and the compiler must continue to reject it as a complete front
> end.  The exact tangent audit also corrects an earlier overclaim: the first
> four zero layers do not force divisibility by the common root.  Their true
> source-band kernels already have dimensions `6,6,7,7`, and layer 35 has
> dimension `10`.  Nonlinear rows then force legitimate root continuation
> through descent `7`; the first surviving local defect is the descent-`8`
> double-prime ratio `27*y^2-9*y+1=0`.

The executable record is
[`cas/f2_75_125_frontend.py`](cas/f2_75_125_frontend.py).  It emits the
machine-readable residual certificate
`plane-jc.f2-75-125-residual.v3`.  The exact finite layer-envelope replay is
[`cas/classify_f2_75_125_layers.py`](cas/classify_f2_75_125_layers.py).

## 1. What the F2 table and chain theorem force

Family `F2` at `j=1` has

\[
 A_0=(5,20),\quad A'_0=(1,0),\quad
 A_1=(7/5,2),\quad (m,n)=(3,5).
\]

Thus `dir(A0-A0')=(5,-1)`.  The complete-chain theorem passes from
`L=K[x,y]` to `L^(5)` by choosing a nonzero root whose multiplicity is `2m`
in the `P` edge polynomial and applying

\[
 y\longmapsto y+\lambda x^{-1/5},\qquad \lambda\ne0.
\]

Put `X=x^(1/5)`.  The bracket changes by the chain rule:

\[
 [P,Q]_{X,y}=5X^4[P,Q]_{x,y},
\]

so scalar normalization makes the transformed bracket exactly `X^4`.
The two consecutive forced edges in the integral `(X,y)` lattice are

| edge | `P` endpoints | `Q` endpoints | weight |
| --- | --- | --- | --- |
| translated type II | `(75,60)` to `(21,6)` | `(125,100)` to `(35,10)` | `(1,-1)` |
| final type I | `(21,6)` to `(4,1)` | `(35,10)` to `(1,0)` | `(5,-17)` |

Every displayed endpoint is an actual vertex, hence its coefficient is
nonzero.  These facts are consequences of the regular-corner endpoints; they
do **not** say that the four displayed points are the complete polygons.

## 2. Complete normalization of the terminal edge

Write `s=X^17 y^5`.  Torus rescaling of `X,y,P,Q`, preserving the normalized
bracket, makes the terminal edge

\[
 P_{\rm I}=X^4y(1+s),\qquad
 Q_{\rm I}=-X\left(1+3s+\frac95s^2\right).
\]

Direct differentiation gives

\[
 [P_{\rm I},Q_{\rm I}]_{X,y}=X^4.
\]

Equivalently, before normalization, if the endpoint coefficients are
`a,b` for `P` and `c,d,e` for `Q`, the bracket equations are

\[
 -ac=1,\qquad 2ad-6bc=0,\qquad 5ae-3bd=0.
\]

All five coefficients are nonzero and the last three are uniquely fixed once
the two `P` vertices are normalized.  This block therefore has obstruction
rank zero: it is the mandatory terminal type-I bracket, not the new family
obstruction.

## 3. The forced common-power band

For

\[
 t=Xy,\qquad z=y^{-1},\qquad [t,z]_{X,y}=-z,
\]

the translated type-II edge is the common-power band

\[
 C_{\rm top}=t^7H(t)z^5,\quad \deg H=18,
\]

where the common root is scaled so that \(H(0)=1\) and the leading
coefficient of `H` is nonzero, and

\[
 P_{\rm top}=t^{21}H(t)^3z^{15},\qquad
 Q_{\rm top}=-\frac95t^{35}H(t)^5z^{25}.
\]

The factor `-9/5` is forced by the normalized coefficient at `(35,10)`.
Its bracket vanishes identically.  More generally,

\[
 [p_i(t)z^i,q_j(t)z^j]_{X,y}
 =\bigl(i\,p_iq'_j-j\,p'_iq_j\bigr)z^{i+j}.
\]

Thus the formal top is layer `15+25=40`, not layer 39.  The normalized
right-hand side is `X^4=t^4z^4`, so the 35 missing **zero** layers are
exactly `39,38,...,5`; layer 4 carries the target.  The former constant
bracket formula was copied from the different chart
`t=xy^2,z=y^-1` used for `(72,108)` and does not apply here.

## 4. Exact B0 classification of the 35 layers

Let `x^i y^j` be a source monomial and put
\(\ell=5i-j\).  After `x=X^5`, the Puiseux translation, and scaling the
nonzero translation root to one, it contributes

\[
 t^\ell(1+t)^jz^\ell,\qquad j=5i-\ell.
\]

The total-degree bounds give

\[
0\le i,\quad 0\le 5i-\ell,\quad 6i-\ell\le D,
\]

with `D=75` for `P` and `D=125` for `Q`.  The terminal supporting
halfspaces are

\[
5a-17b\le3\quad(P),\qquad 5a-17b\le5\quad(Q).
\]

On band \(\ell=a-b\), these inequalities force the initial Taylor jets of
the linear combination of \((1+t)^j\) to vanish.  The jet matrix is
\(\binom{j}{k}\); its falling-factorial Vandermonde minors have full rank.
This proves a finite linear band space without assuming a lower polygon or
a `gamma` branch.

For source basis elements
\[
p=t^\ell(1+t)^r,\qquad q=t^m(1+t)^s,
\]
the complete quadratic contribution compresses to

\[
\ell pq'-mp'q=(\ell s-mr)t^{\ell+m}(1+t)^{r+s-1}. \tag{1}
\]

Equation (1) reconstructs every scalar coefficient equation exactly.  The
checker classifies all 35 layers by band pair, \(t\)-exponent, and Kummer
character.  In the missing window it finds:

| quantity | exact count |
| --- | ---: |
| zero layers | 35 (`39` through `5`) |
| contributing band pairs | 665 |
| linear parameters after the terminal jet equations | 978 |
| affine parameters after the five endpoint normalizations | 973 |
| active compressed binomial generators | 165,980 |
| raw scalar Keller coefficient rows | 5,625 |
| structurally active rows after the terminal jets | 5,348 |

The pinned JSON stores every band, all terminal normalizations, every
band-pair incidence row, the five character counts, and SHA-256 digests of
the exact linear equations and compressed quadratic generators.  This is
the promised B0 classification.  It also shows why an unfiltered Gröbner
calculation is premature: any future coefficient-route exclusion must first
cut the 978-parameter jet over-envelope (973 after normalization) to
exhaustive B1 polygon/`gamma` masks.

### 4.1 Joint classification of the common-power top band

The top two band spaces already give a first exact B1 cut.  Set \(u=1+t\)
and

\[
A(u)=(u-1)^2H(u-1).
\]

The `P` top band belongs to the translated source space precisely when
\(A(u)^3\in k[u^5]\), and the `Q` top band gives
\(A(u)^5\in k[u^5]\).  For every fifth root of unity \(\zeta\), these two
identities give

\[
A(\zeta u)^3=A(u)^3,\qquad A(\zeta u)^5=A(u)^5.
\]

Since \(2\cdot3-5=1\), division in \(k(u)\) gives
\(A(\zeta u)=A(u)\).  Hence \(A(u)=G(u^5)\).  The double zero of \(A\) at
\(u=1\) forces \((v-1)^2\mid G(v)\), so every possible common top is

\[
\boxed{
H(t)=(1+u+u^2+u^3+u^4)^2R(u^5),\quad
u=1+t,
}
\]

where

\[
R(v)=av^2+bv+\left(\frac1{25}-a-b\right),\qquad a\ne0.
\]

Conversely this formula satisfies both source-band conditions, has
\(H(0)=1\), and has degree 18.  Thus the apparent 18-coefficient common root
is exactly a two-parameter quadratic family.  Its natural algebraic root
strata are `disc(R) != 0`, `disc(R)=0`, and `R(0)=0`; identifying any of
these with an exhaustive historical `gamma` branch would require a separate
normal-form argument.  This joint top cut lowers the normalized B0 count
from 973 only to 959, so it is substantial but not yet a feasible global
Gröbner system.

### 4.2 Correction: the exact top tangent is much larger

The previous version of this section inserted

\[
p_{15-\delta}=C_0^2U,\qquad
q_{25-\delta}=-\frac95C_0^4V \tag{2}
\]

before proving either divisibility.  Equation (2) does reduce the restricted
slice to the displayed ODE, but the exact source bands do **not** force (2).
The resulting claim that layers `39,...,36` continue a polynomial source root
was therefore invalid.

The unrestricted calculation is simpler.  For an arbitrary exact source-band
element (p=p_{15-delta}), put

\[
q=-3C_0^2p. \tag{3}
\]

Since

\[
P_0=C_0^3,\qquad Q_0=-\frac95C_0^5,
\]

direct substitution in the layer-(40-delta) bracket gives zero.  Moreover,
(C_0) is the band-`5` part of the genuine degree-`25`, height-`1` source
polynomial

\[
\boxed{C_R=x(xy^5-1)^2R(xy^5).} \tag{4}
\]

Consequently multiplication by (C_0^2) maps every exact P band into the
paired exact Q band: source degree, terminal height, and Laurent band increase
by `50`, `2`, and `10`.  This is also checked directly in the factored band
bases

\[
V_\ell=t^\ell u^{j_0}(u^5-1)^\nu
        \mathbf Q[u^5]_{\le N-\nu}. \tag{5}
\]

Thus every P-band direction survives the homogeneous top equation.  The
q-only kernel satisfies

\[
5C_0q'-(25-\delta)C_0'q=0,
\qquad q^5=cC_0^{25-\delta}. \tag{6}
\]

The fixed factor ((u^5-1)^2) in (C_0) makes (6) a source Laurent
polynomial exactly for
(delta=5,10,15,20,25), with respective generators

\[
C_0^4, C_0^3, C_0^2, C_0, 1. \tag{7}
\]

The exact kernel dimensions at descents `1,...,5` are therefore

\[
\boxed{6, 6, 7, 7, 10}, \tag{8}
\]

not `2,2,3,3,6`.  In particular, the extra layer-`35` direction is the
ordinary commuting (C_0^4) term in Q after a formal-root gauge change.  It
is not the modified-Laurent residue `lambda*C^(-1)`.

The actual formal (C_0^{-1}) resonance occurs only at descent `30`, hence
layer `10`; (C_0^{-2}) occurs at descent `35`, layer `5`.  Neither is an
independent homogeneous source-band kernel because a negative power of
(C_0) is not a Laurent polynomial.  Such a coefficient can arise only
through nonlinear cancellation with the forced `F` tail.

For the nonlinear descent, let

\[
T_\delta(q)=5C_0q'-(25-\delta)C_0'q
\]

on the exact Q band.  The new P columns lie in
\(\operatorname{im}T_\delta\) by (3), so a known lower-layer forcing
(f_\delta) is soluble precisely when

\[
\operatorname{rank}[T_\delta\mid f_\delta]
=\operatorname{rank}T_\delta, \tag{9}
\]

equivalently when the corresponding maximal-minor/Fitting class vanishes.
This is the correct triangular handoff.  It retains all 35 nonlinear zero
layers; there is no justified five-layer common-root reduction.

### 4.3 Nonlinear recovery and the first exact residual branch

The larger tangent kernel does not remain arbitrary.  On the primitive
centralizing branch, expand the particular solution

\[
Q=-\frac95P^{5/3}.
\]

If (p_i) is the P coefficient at relative descent (i), its first rows are

\[
q_1=-3C_0^2p_1,
\]

\[
q_2=-3C_0^2p_2-\frac{p_1^2}{C_0}, \tag{9a}
\]

\[
q_3=-3C_0^2p_3-\frac{2p_1p_2}{C_0}
       +\frac{p_1^3}{9C_0^4}. \tag{9b}
\]

Thus source polynomiality first requires

\[
C_0\mid p_1^2,
\qquad
C_0^4\mid p_1(p_1^2-18C_0^3p_2). \tag{9c}
\]

At a multiplicity-one prime of (C_0), (9c) forces
(v(p_1)\ge2).  At a multiplicity-two prime, rows two and three leave only
the exceptional valuation

\[
v(C_0)=2,\qquad v(p_1)=3,\qquad v(p_2)=0. \tag{9d}
\]

Scale the descent parameter so that the residue of P in (9d) is
(1+as+bs^2), with (ab\ne0).  The fourth and fifth coefficients of its
(5/3)-power have primitive numerators

\[
E_4=a^4-9a^2b+27b^2,
\]

\[
E_5=7a^4-60a^2b+135b^2.
\]

Their exact resultant is

\[
\operatorname{Res}_b(E_4,E_5)=1701a^8. \tag{9e}
\]

Hence the exceptional double-prime valuation is impossible whenever both
rows are zero.  Nonnegative centralizer terms cannot repair this leading
residue: a term (C_0^{5-j}) has valuation gap
(j(15/k-2)>0) when the first defect has descent (k\le7).

It follows that a first non-root P defect at any descent `1,...,7` is forced
back into the (C_0^2) root slice.  Equivalently, after legitimate recursive
root changes, P continues as a source cube through band `8`.  This recovers
more common-root continuation than the old divisible tangent argument, but
now from the nonlinear rows and with a proof.

The threshold is sharp for the present ledger.  At first-defect descent `8`,
rows `16,24,32` precede the target descent `36`, while row `40` lies after it.
The first surviving local branch is supported only at multiplicity-two primes
and has

\[
v(C_0)=2,\qquad v(p_8)=3,\qquad v(p_{16})=0,
\]

with normalized ratio (y=b/a^2) satisfying

\[
\boxed{27y^2-9y+1=0},
\qquad \operatorname{disc}=-27. \tag{9f}
\]

Thus the earliest residual lives over (mathbf Q(\sqrt{-3})).  Equation
(9f) is an exact local candidate, not yet a global source-band solution or a
plane counterexample: the target row intervenes before the missing fifth
residue and must now be coupled to all Kummer-conjugate double primes.

The first boundary pivot tested the alternative.  The quadratic `R` produces
four exact cover-level contact partitions, but a contact multiplicity alone
does not determine the first normal order, a toroidal branch scale, or a
finite-normalization ramification row.  Even the strongest naive
contact-to-row promotion survives the coarse finite-flat packet budget.
That negative result is retained in
[`F2_BOUNDARY_HANDOFF.md`](F2_BOUNDARY_HANDOFF.md), but it is no longer the
current stopping point.  The subsequent
[`Kummer-orbit transfer`](F2_KUMMER_ORBIT_TRANSFER.md) groups the nonzero
centers into rigid orbits, excludes the zero-root strata, and reduces the
row to one principal chain or two copies.  The
[`terminal residue calculation`](F2_TERMINAL_RESIDUE_COVER.md) then supplies
the exact degree-six target row.  Thus sequential descent through the lower
thirty layers remains stopped, while the F2 programme is reopened at the
global source/target gluing stage.
<!-- status-consumer: PF2GC1 33dbc5ff48b5d064 -->

## 5. Why the older `(50,75)` calculation does not fill the gap

Section 5 of the 2014 polynomial-system paper treats the `j=0`, `(m,n)=(2,3)`
member.  It states two preliminary cases `gamma=2,3`, but explicitly says
that no proofs are supplied for that first reduction.  Its later coefficient
systems use `P=C^2`, `Q=C^3+...`; replacing those exponents by `3,5` is not a
consequence of the printed argument and would not establish exhaustive lower
supports.  For the F2 chain itself, the generated-corner identity now fixes
`gamma=2` and hence `d=3`; see
[`F2_MODIFIED_CHART_BRIDGE.md`](F2_MODIFIED_CHART_BRIDGE.md).  This removes
the old gamma ambiguity but does not justify the Laurent cut.

## 6. Residual coefficient route and current obstruction

If the direct coefficient route is resumed, the JSON certificate leaves
these proof obligations:

1. extend the exact binomial-jet image from its nonnegative-`xi` projection to
   a complete positive/negative Laurent-band ledger;
2. derive the polynomial tail correction forced by positive-negative bracket
   cross terms (naive truncation is ruled out by the projected top-band unit
   ideal); and
3. retain the full exact P-band kernel and compile each nonlinear forcing in
   the cokernel of `T_delta`; the actual formal `lambda*C0^(-1)` location is
   descent `30` (layer `10`), where it is not an independent source mode.

This remains a valid route to a coefficient-level F2 exclusion, but it is no
longer the next programme-wide obstruction.  The Kummer and terminal-residue
theorems bypassed it for the selected principal chain and produced the exact
row `(e,f)=(1,6)`.  The current obstruction is global source/target gluing,
including the same-target versus distinct-target double-packet split,
spectator incidence, the three forced interior attachment points over the
target toric nodes, the endpoint-over-smooth branch point, and the
purity-forced affine ramification row.

The exact all-`r` formal compiler, its 14/22-function `r=3` windows, the
direct Artinian residue Fitting ideals, and the proof that the exact `r=3`
polynomial projection has unit top-band ideal are recorded in
[`F2_MODIFIED_LAURENT_FAMILY.md`](F2_MODIFIED_LAURENT_FAMILY.md).  That
calculation eliminates every branch of the naive projection but does not fill
the tail-correction obligations 1--3 above.
<!-- status-consumer: PF2MCB1 6ff13314e0090f52 -->

## 7. Kummer-character bridge gate

The terminal bracket \(X^4\) suggests the Kummer coordinate \(u=X^5/5\).
A monomial-Jacobian block descends to a constant-Jacobian map in \(k[u,y]\)
only when every \(X\)-exponent is divisible by five.  Here the terminal
characters are

\[
P:\{1,4\},\qquad Q:\{0,1,3\}\pmod 5.
\]

Thus the five-term terminal block does not descend to the certified
support-six Keller theorem.  The missing lower bands must be compiled with
their characters before that **coefficient/support** route can be used; they
are not a prerequisite for globally gluing the separately certified terminal
target row.  The exact gate and the obstruction to controlling affine support
from coarse Newton data are proved in
[`AFFINE_SUPPORT_NEWTON_BRIDGE.md`](AFFINE_SUPPORT_NEWTON_BRIDGE.md).

## 8. Reproduction

```bash
.venv/bin/python plane-jc/cas/classify_f2_75_125_layers.py
.venv/bin/python plane-jc/cas/test_f2_75_125_frontend.py
.venv/bin/python plane-jc/cas/f2_75_125_frontend.py
```

Intentional regeneration of
[`../artifacts/generated-results/jc2_f2_75_125_character_layers.json`](../artifacts/generated-results/jc2_f2_75_125_character_layers.json)
uses `--refresh`.  The corrected `v2` artifact has whole-file SHA-256
`3b30a686a27c1adebc37e354cf5d0a2a60a07ef811c96bfbc02a426886889cf7`.
The final command prints the machine-readable partial front-end certificate.
