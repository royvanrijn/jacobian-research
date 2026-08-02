# F2 `j=1` derivation audit for `(75,125)`

> **Status: forced skeleton, not a normal-form certificate.**  The published
> complete-chain theorem determines two consecutive edges and one Puiseux
> translation.  It does not determine the rest of either Laurent polygon.
> Consequently this note does not claim to eliminate `(75,125)`, and the
> compiler must continue to reject it as a complete front end.

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
calculation is premature: the next theorem must cut the 978-parameter jet
over-envelope (973 after normalization) to the exhaustive B1
polygon/`gamma` masks.

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

### 4.2 Exact upper descent and the first surviving mode

The layer equations are nevertheless triangular at the top.  Write
\(C_0=t^7H(t)\).  At descent \(r\), after subtracting the contributions
already generated by higher common-root bands, put

\[
p_{15-r}=C_0^2U,\qquad
q_{25-r}=-\frac95C_0^4V,\qquad W=5U-3V.
\]

The layer \(40-r\) equation reduces exactly to

\[
(5-r)C'_0W-5C_0W'=0. \tag{2}
\]

Since \(\operatorname{ord}_{t=0}C_0=7\), a nonzero homogeneous solution
would have

\[
\operatorname{ord}_{t=0}W=\frac{7(5-r)}5.
\]

For \(r=1,2,3,4\) this is nonintegral.  Therefore \(W=0\), and layers
39 through 36 force continuation of a single common root.  The successive
new root bands have dimensions `2,2,3,3`; their exact forms are recorded in
the JSON.

At \(r=5\), equation (2) becomes \(W'=0\).  A genuine one-dimensional mode
survives:

\[
\lambda C_0^2z^{10}
\]

in `P`, alongside a five-parameter common-root correction.  It commutes
with the top `Q` term because

\[
[C_0^2z^{10},C_0^5z^{25}]=0.
\]

Thus the first four zero layers do not obstruct F2, and layer 35 does not
produce a unit ideal: it opens the expected lower-power \(C_0^2\) branch.
Only layers 34 through 5 remain unclassified at this stronger level.  This
is the decisive diagnostic for the present approach: closure is plausible
only if propagation of this \(C_0^2\) mode creates an early nonzero de Rham
class; otherwise the 30-layer continuation is too broad for a final direct
elimination.

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

## 5. Why the older `(50,75)` calculation does not fill the gap

Section 5 of the 2014 polynomial-system paper treats the `j=0`, `(m,n)=(2,3)`
member.  It states two preliminary cases `gamma=2,3`, but explicitly says
that no proofs are supplied for that first reduction.  Its later coefficient
systems use `P=C^2`, `Q=C^3+...`; replacing those exponents by `3,5` is not a
consequence of the printed argument and would not establish exhaustive lower
supports.

## 6. Residual system and next obstruction

The JSON certificate contains three proof obligations:

1. prove support control after `y -> y+lambda*x^(-1/5)`;
2. determine, prove exhaustive, and normalize every `gamma` branch for
   `(m,n)=(3,5)` (the old values `2,3` may not simply be assumed);
3. refine the exact B0 band envelope to exhaustive B1 polygon masks and
   impose the common-power coefficient relations before compiling the first
   feasible Gröbner/de Rham block.

This is the next genuinely new obstruction.  The terminal edge is already
solved and contributes no de Rham obstruction; the missing object is a
theorem-level exhaustive lower-boundary classification, not a larger
Gröbner-basis calculation.

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
their characters before the support theorem or the log-boundary ledger can
be used.  The exact gate and the obstruction to controlling affine support
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
uses `--refresh`.  The final command prints the machine-readable partial
front-end certificate.
