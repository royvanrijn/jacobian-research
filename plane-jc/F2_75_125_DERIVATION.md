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
> double-prime ratio `27*y^2-9*y+1=0`.  Its four fixed Kummer-orbit
> supports are now excluded by the normalized `Q_1` endpoint.  The only
> remaining support is a nonzero double root of `R`; that branch passes an
> exact normalized local target-jet test.  The corner inequalities have also
> been carried through the entire lower Laurent tail: there are `240` zero
> layers from `40` through `-200` (apart from target layer `4`) and `2,418`
> jet-reduced B0 parameters.  Thus the fifth fractional-power numerator at
> descent `40` is not an equation by itself: it lies on a genuine lower-tail
> Fitting row.  The exact first Kummer return has now been compiled.  Its
> omitted fifth-binomial term cancels the apparent quartic-terminal top
> conflict, and the remaining packet is a surjective `8 x 10` map with unit
> minor `3^5*5^16*e^13` over the movable-double-root base.  Both roots of
> `27*y^2-9*y+1` survive, and the same source-return/Bezout mechanism absorbs
> every first-appearance edge forcing through `v^10`.  Every one of the `53`
> pivot columns in those unit minors has now been lifted to an original
> polynomial source combination strictly inside both certified supporting
> edges.  The determinants and lifts generalize to `18*r-1` strict-interior
> pivots for every `r>=2`, proving uniform return surjectivity through order
> ten.  Thus this freedom is genuine within the corner envelope, not an
> artifact of independent Laurent coefficients.  Exact three-point CRT then
> removes the controllable `w=0` block and leaves a rank-`24` global Hermite
> module at `w=1,w0` over the rank-two candidate algebra.  The complete
> fixed-endpoint block can now be eliminated: its leading target coordinate
> is the terminal normalization, and the other ten coordinates have constant
> determinant `75000`.  Exact source lifts preserve the controlled `w=0`
> block, so only `13` global Hermite coordinates remain after substitution.
> The substitution is now carried exactly: it produces degree-at-most-eight
> circuits with `1,061` active source coordinates, not a thirteen-variable
> ring.  Unit minors eliminate `134` endpoint-disjoint new-Q coordinates on
> layers `39..29`, leaving `927` active coordinates.  The live algebraic gap
> begins precisely at descent `12` (layer `28`), where the P3/Q13 power row
> couples to the endpoint solution; its Schur/Fitting system must be carried
> through descent `37` and the thirteen residual functionals.

The carrier-specialization continuation in
[`F2_75_125_CARRIER_SPECIALIZATIONS.md`](F2_75_125_CARRIER_SPECIALIZATIONS.md)
now fixes every exposed Schur and Hermite kernel/cokernel over
`Q[rho]/(rho^2-3*rho+1)`.  The `53` successive coordinates are the normalized
factor-quotient block; the full coupled Laurent cokernel has `347`
coordinates, including `294` divisibility/local-jet conditions.  The
downstream
[`nonlinear forcing compiler`](F2_75_125_NONLINEAR_FORCING.md) now substitutes
the ten endpoint circuits and upper tangent solutions into all of them,
appends the terminal `7+6` block and the descent-eight incidence rows, and
routes the rational squarefree carrier through later spacings `9..90`.  The
localized circuit-ideal test and the squarefree target-and-tail compilers
remain open.
<!-- status-consumer: PF2CS1 666da98d2d24669e -->
<!-- status-consumer: PF2NF1 cfd1da5136c0b6d0 -->

The executable record is
[`cas/f2_75_125_frontend.py`](cas/f2_75_125_frontend.py).  It emits the
machine-readable residual certificate
`plane-jc.f2-75-125-residual.v3`.  The exact finite layer-envelope replay is
[`cas/classify_f2_75_125_layers.py`](cas/classify_f2_75_125_layers.py), and
the carried endpoint/upper-power reduction is
[`cas/reduce_f2_75_125_endpoint_system.py`](cas/reduce_f2_75_125_endpoint_system.py).

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
checker classifies all 35 upper zero layers by band pair, \(t\)-exponent,
and Kummer character.  In that window it finds:

| quantity | exact count |
| --- | ---: |
| zero layers | 35 (`39` through `5`) |
| contributing band pairs | 665 |
| linear parameters after the terminal jet equations | 978 |
| affine parameters after the five endpoint normalizations | 973 |
| active compressed binomial generators | 165,980 |
| raw scalar Keller coefficient rows | 5,625 |
| post-jet support-row upper bound | 5,344 |

The `5,344` count is deliberately labelled an upper bound: independent
coefficient supports in two jet-reduced bands can still cancel as a bilinear
common-power pair.  Exact reconstruction instead uses the stored jet
equations together with all `165,980` compressed source-basis generators.

The pinned JSON stores every upper band, all terminal normalizations, every
band-pair incidence row, the five character counts, and SHA-256 digests of
the exact linear equations and compressed quadratic generators.  It now also
stores the complete lower B0 envelope described in Section 4.4.  This is a
necessary coefficient system, not an assertion that every allowed monomial
occurs on an actual lower Newton polygon.

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

There is now an exact orbit filter.  At descent (3\cdot8=24), the
exceptional valuation makes the absolute Q-band-one coefficient vanish at
its supporting prime.  The complete source band is

\[
q_1=t u^4S(u^5),\qquad \deg S\le20,
\]

and the terminal normalization is (S(1)=[t]q_1=-1).  Therefore (q_1) is
nonzero at each of the four nontrivial roots of (u^5=1).  None of the fixed
multiplicity-two primes can support (9f).  The only remaining possibility is

\[
\boxed{R(w)=\alpha(w-w_0)^2,\quad w_0\ne0,1,\quad S(w_0)=0.} \tag{9g}
\]

Normalization fixes

\[
R(w)=\frac{(w-w_0)^2}{25(1-w_0)^2}.
\]

Hence the tiny exact core of this earliest branch is the finite étale
rank-two algebra

\[
\mathcal A_8=
B[y]/(27y^2-9y+1),\qquad
B=\mathbf Q[w_0,w_0^{-1},(w_0-1)^{-1}]. \tag{9g'}
\]

In the basis `(1,y)`, multiplication by `y` has matrix

\[
\begin{pmatrix}0&-1/27\\1&1/3\end{pmatrix},
\]

with norm `1/27` and discriminant `-1/27`.  Thus (9g') is reduced and splits
into two components after adjoining (\sqrt{-3}).

The first target jet does not remove (9g).  Indeed, the exact bands

\[
q_1=t u^4S(w),\qquad
p_3=t^3u^2(w-1)A(w),\qquad w=u^5,
\]

admit the degree-one interpolants

\[
S(w)=\frac{w-w_0}{w_0-1},
\]

\[
A(w)=\frac15+\frac{w-1}{w_0-1}
\left(\frac1{15w_0^2}-\frac15\right). \tag{9h}
\]

They satisfy the endpoint conditions (S(1)=-1), (A(1)=1/5), and at
(u_0=1+t_0), (w_0=u_0^5),

\[
q_1(t_0)=0,
\qquad
3p_3(t_0)q_1'(t_0)=t_0^4. \tag{9i}
\]

Equation (9i) verifies compatibility of the terminal (P_3/Q_1) summand at
the movable double root.  It is not a solution of every other summand in the
global layer-four polynomial identity.

### 4.3a The complete target cokernel

The full target row can now be written without expanding its 34 scalar
coefficients blindly.  Its new bands are

\[
p_{-21}=t^{-21}u^{21}A(w),\quad \deg A\le9,
\qquad
q_{-11}=t^{-11}u^{11}S(w),\quad \deg S\le19.
\]

After removing the P follower, their combined variable is

\[
S+3w^2(w-1)^4R(w)^2A.
\]

On the double-root stratum put

\[
E(w)=(w-1)(w-w_0).
\]

Up to a unit in
`B=QQ[w0,w0^-1,(w0-1)^-1]`, the new-band target operator is

\[
\mathcal L_{36}(S)=15w^2E^5
\left(5wE S'+(11E+22wE')S\right). \tag{9j}
\]

Thus twelve target conditions are the divisibility by

\[
w^2(w-1)^5(w-w_0)^5.
\]

After this factor is removed, the remaining map sends polynomials of degree
at most `19` to polynomials of degree at most `21`.  In monomial bases its
column `j` has entries

\[
(5j+11)w_0,qquad -(5j+33)(1+w_0),qquad 5(j+11)
\]

in rows `j,j+1,j+2`.  The first twenty rows form a full-rank triangular
minor.  Solving them leaves exactly two residual equations,

\[
\rho_{20}=f_{20}+128(1+w_0)s_{19}-145s_{18},
\qquad
\rho_{21}=f_{21}-150s_{19}. \tag{9k}
\]

Consequently the complete target cokernel has rank `14`: twelve local jets
and the two rows (9k).  Equation (9i) checks only one of those fourteen
coordinates.

The complete old B0 band span does fill this cokernel before the earlier
Keller rows are imposed.  Two explicit sets of fourteen old source-basis
generators give `34 x 34` minors whose gcd is `w0^12`; this is a unit over
`B`.  This is a span statement, not simultaneous nonlinear solvability.

There is a particularly transparent missing coordinate.  Both the
`P_3/Q_1` contribution and (9j) vanish at `w=0`, whereas the target has value
one.  Across all old bands only `(P_-1,Q_5)` and `(P_5,Q_-1)` can supply the
constant term, so every globalization must satisfy

\[
5\left(A_{P,5}(0)B_{Q,-1}(0)
-A_{P,-1}(0)B_{Q,5}(0)\right)=1. \tag{9l}
\]

In fact all lowest-`u` target coefficients combine into one edge equation.
Writing `x=z^5`,

\[
P=A(x)+uz^{-1}B(x)+O(u^2),\qquad
Q=C(x)+uz^{-1}D(x)+O(u^2),
\]

gives

\[
\boxed{A'D-BC'=1/5}. \tag{9m}
\]

It forces the complete edge to break from a nonconstant common-power pair:
if `A=G^3` and `C=(-9/5)G^5`, the left side is divisible by `G^2`.  But the
edge equation itself survives.  With the unit `e=-R(0)`, an exact witness is

\[
A=x+e^3x^3,quad B=e/5,quad
C=-\frac95e^5x^5,quad D=\frac15-\frac35e^3x^2. \tag{9n}
\]

Thus the local target jet does not globalize using only the new target
bands, but the necessary off-grid edge correction has an exact solution.

The sparse correction has an exact all-family explanation.  Put

\[
m=2r-1,\qquad
\beta_r=\frac{2r^2}{(r-1)(2r-1)}.
\]

For every `r>=2`, the four polynomials

\[
\begin{aligned}
A&=x+e^r x^r,& B&=\frac{r-1}{10}e,\\
C&=-\beta_r e^{2r-1}x^{2r-1},&
D&=\frac15(1-r e^r x^{r-1})
\end{aligned} \tag{9o}
\]

satisfy (A'D-BC'=1/5).  The witness (9n) is exactly the `r=3`
specialization.  This does not produce a four-term Keller pair.  In the
auxiliary transverse coordinate (v=uz^{-1}),

\[
\det\frac{\partial(A+Bv,C+Dv)}{\partial(x,v)}
=\frac15+
\frac{r(r-1)^2}{50}e^{r+1}x^{r-2}v. \tag{9p}
\]

If (P=A+Bv) is kept linear in (v), cancellation of (9p) forces

\[
[v^2]Q=\frac{BD'}{2A'}
=-\frac{r(r-1)^2e^{r+1}x^{r-2}}
{100(1+r e^r x^{r-1})}. \tag{9q}
\]

The denominator has degree `r-1` and the numerator degree `r-2`, so (9q)
is not a polynomial for any `r>=2`.  Thus the literal sparse polynomial
escape fails at the very next transverse coefficient.

It does, however, complete exactly as a formal series.  Set

\[
\kappa=(5B)^{-1},\qquad
H(A(x))=C(x)+\kappa x,
\]

and define

\[
P=A(x)+Bv,\qquad Q=-\kappa x+H(P). \tag{9r}
\]

Then (partial(P,Q)/\partial(x,v)=B\kappa=1/5), while differentiation
of the defining composition gives (Q(x,0)=C(x)) and
(Q_v(x,0)=D(x)).  For `r=3`, the beginning of the forced series is

\[
H(s)=e^{-1}s-e^2s^3+\frac65e^5s^5-3e^8s^7
+10e^{11}s^9+\cdots.
\]

This tail cannot terminate.  If (H) were a polynomial of degree (q),
then (\deg H(A)=rq), but (\deg(C+\kappa x)=2r-1); the equation
(rq=2r-1) has no integral solution.  The edge escape is therefore
formally integrable but not a polynomial counterexample in its sparse form.
A more general completion must turn on a (v^2) coefficient in `P` and
then further off-grid bands.  Excluding those corrections is precisely the
simultaneous target-and-tail Fitting problem, so (9q) is a genuine pruning
theorem rather than an F2 exclusion.

That next correction can also be classified exactly.  Write

\[
P=A+Bv+P_2v^2,\qquad Q=C+Dv+Q_2v^2.
\]

The coefficient of (v) gives the Bezout equation

\[
A'Q_2-C'P_2=\frac12BD'. \tag{9s}
\]

A uniform particular solution is

\[
\begin{aligned}
P_2^{(0)}&=-\frac{r(r-1)^3}{200}
 e^{r+2}x^{r-2},\\
Q_2^{(0)}&=-\frac{r(r-1)^2}{100}e^{r+1}x^{r-2}
+\frac{r^2(r-1)^2}{100}e^{2r+1}x^{2r-3}.
\end{aligned} \tag{9t}
\]

Since (gcd(A',C')=1), every polynomial solution within the exact edge
degree bounds is

\[
P_2=P_2^{(0)}+A'T_2,qquad
Q_2=Q_2^{(0)}+C'T_2,qquad
T_2=\alpha_2+\beta_2x. \tag{9u}
\]

Thus the forced second-order repair is an affine two-parameter family, not
an obstruction.  It still cannot terminate at order two.  For `r>=3`, the
coefficient of (v^2) successively gives

\[
[x^0]=\frac15\beta_2,qquad
[x^{r-2}]_{\beta_2=0}
=\frac{3r(r-1)}5e^r\alpha_2,
\]

and then

\[
[x^{2r-4}]_{\alpha_2=\beta_2=0}
=-\frac{3r^2(r-1)^4}{1000}e^{2r+2}\ne0. \tag{9v}
\]

For `r=2`, the analogous pivots occur at (x^2,x,1) and end in
(-3e^6/250\ne0).  Hence no member of the all-`r` sparse section has a
quadratic transverse polynomial completion.

For `r=3`, the exact recursion has now been continued through the last order
before a Kummer return.  At each order (i=3,4), a particular solution is
adjoined and the full freedom is again

\[
(P_i,Q_i)\longmapsto(P_i+A'T_i,Q_i+C'T_i),qquad
T_i=\alpha_i+\beta_i x.
\]

If the expansion stops cubically, five coefficients of the (v^3) row force

\[
\beta_2=\alpha_2=\beta_3=0,qquad
\alpha_3=\frac6{125}e^6,
\]

and then leave the nonzero residue

\[
\frac{108}{125}e^{12}. \tag{9w}
\]

Thus the cubic coefficient ideal is the unit ideal.  The quartic calculation
is only slightly larger.  The (v^4) row has a nine-element exact Groebner
basis.  Reducing three coefficients of the (v^5) row modulo it gives, up
to invertible scalars,

\[
9\alpha_2e^5+25\alpha_4,qquad
249\alpha_2e^5+275\alpha_4,qquad
1875\beta_4+608e^{10}. \tag{9x}
\]

The first two have determinant (-3750), hence
(alpha_2=alpha_4=0).  A Groebner-basis row then gives
(4375\beta_4+2763e^{10}=0); substituting this into the last expression in
(9x) leaves

\[
-\frac{4033}{7}e^{10}\ne0. \tag{9y}
\]

Therefore the quadratic, cubic, and quartic transverse **terminations** are
all empty.  An exact source escape must reach

\[
\boxed{v^5},
\]

which is exactly the first order at which the invariant return
(w=xv^5) appears and the terminal binomial-jet relations couple coefficients
that were independent at orders below five.  The coupled return can now be
computed exactly, and it survives.

Write

\[
A(x)=\sum_d a_dx^d,\qquad C(x)=\sum_jc_jx^j,
\]

and let

\[
J_4=[v^4]\det\frac{\partial(P,Q)}{\partial(x,v)}
\]

after adjoining (P_5v^5,Q_5v^5).  For a residue-zero source band
(\ell=5d), the coefficients below the return are the expansion of
(a_d(1-u)^{5d}).  Reconstructing all such bands rather than treating
(P_5,Q_5) as independent gives the exact fifth-binomial correction

\[
\boxed{\mathcal R_5=5J_4+\Omega_5(A,C)}, \tag{9z}
\]

where (\mathcal R_5) packages the (u^4) coefficients on Laurent layers
(0,5,\ldots,40) and

\[
\Omega_5(A,C)=25\sum_{d,j}a_dc_j
\left(d\binom{5j}{5}-j\binom{5d}{5}\right)x^{d+j}. \tag{9aa}
\]

For (9n), this is

\[
\Omega_5=-28125e^5x^6(85+231e^3x^2). \tag{9ab}
\]

There is a sign worth making explicit.  Since (t=u-1=-1) at the edge,
the edge parameter in (9n) is

\[
e=-R(0)=-\frac{w_0^2}{25(1-w_0)^2}. \tag{9ac}
\]

The common-power top bands fix the two leading return coefficients:

\[
[x^4]P_5=-\frac{3e^3(1003w_0+2)}{w_0},\qquad
[x^6]Q_5=\frac{18e^5(5314w_0+1)}{w_0}. \tag{9ad}
\]

Substitution into (9z) makes its (x^8), or layer-40, coefficient vanish
identically, as it must from the common-power bracket.  The remaining
coefficients (x^0,\ldots,x^7) are eight linear equations in the ten free
return coefficients

\[
[x^{0..3}]P_5,\qquad [x^{0..5}]Q_5.
\]

Their `8 x 10` matrix has a maximal minor

\[
\boxed{3^5 5^{16}e^{13}}. \tag{9ae}
\]

This is a unit in

\[
B=\mathbf Q[w_0,w_0^{-1},(w_0-1)^{-1}],
\]

and remains a unit after the rank-two base change

\[
B\longrightarrow B[y]/(27y^2-9y+1).
\]

Thus the first-return cokernel is zero and its solution space is affine of
relative dimension two.  In particular, both conjugate descent-eight
branches survive the complete first-appearance (v^5) packet.  The nonzero
quartic-terminal residue (9y) was an artifact of setting the return bands to
zero; (9aa) is the source term that absorbs it.

Nor does the obstruction merely move to the next coefficient.  Conditional
on the already exposed source jets, orders (v^6,\ldots,v^9) each supply
(5+7) fresh P/Q coefficients.  At (v^{10}) the fixed common-power top
returns, but the lower coefficients again supply (5+7) freedoms.  The common
Bezout map

\[
(P_*,Q_*)\longmapsto A'Q_*-C'P_*,\qquad
\deg P_*\le4,\quad\deg Q_*\le6,
\]

has a `9 x 12` matrix with unit minor

\[
81e^{10}. \tag{9af}
\]

Consequently every first-appearance forcing through (v^{10}) is absorbed.
This identifies the failure of the standalone edge strategy: deeper
transverse termination tests do not approach an exclusion.  The information
not seen at (u=0) is the simultaneous target and layer-zero Hermite data at
(w=1) and (w=w_0), together with the rest of the lower Laurent tail.  Those
global rows, rather than (v^{11}) truncation, are the next genuine gap.

### 4.3b Source-lift decision and the exact global quotient

There is a necessary audit before accepting that conclusion: the pivot
coefficients used in (9ae) and (9af) might conceivably be freedoms of the
Laurent coefficient box which do not lift to polynomials satisfying the two
corner inequalities.  In fact they all lift.

For a source band on Laurent layer (\ell), write its exact factored form as

\[
t^\ell u^{j_0}(u^5-1)^\nu K(u^5),\qquad \deg K\le N. \tag{9ag}
\]

If a transverse coefficient of order (n) and edge degree (d) is selected,
then (\ell=5d-n) and

\[
k=\frac{n-j_0}{5}\in\{0,\ldots,N\}.
\]

Its triangular source lift is obtained by varying (K) in the direction
((u^5)^k).  In the original monomial basis this is the finite combination

\[
t^\ell\sum_{s=0}^{\nu}(-1)^{\nu-s}\binom{\nu}{s}
u^{j_0+5(k+s)},                                      \tag{9ah}
\]

whose term with index (s) comes from

\[
(i,j)=\bigl(i_{\min}+k+s,\;5(i_{\min}+k+s)-\ell\bigr).
\]

The `v^5` minor uses the two P columns of degrees `1,2` and the six Q
columns of degrees `0,...,5`.  At each of orders `6,...,10`, the common
Bezout minor uses the two P columns of degrees `3,4` and the seven Q columns
of degrees `0,...,6`.  Thus (9ah) gives exactly

\[
8+5\cdot9=53
\]

source lifts.  Two selected bands have a fixed terminal normalization:
`Q_7[x^4]` on layer `13` and `Q_9[x^2]` on layer `1`.  There one replaces
((u^5)^k) by the unit-equivalent direction
((u^5-1)(u^5)^k); a follower coefficient exists in both bands, so this
preserves the normalized value at (u^5=1) without changing the pivot rank.

Expanding all `53` combinations gives nonnegative source exponents.  On P
their largest total degree is `47<75`; on Q it is `73<125`.  Moreover, if
(\nu_{\rm eff}) is (\nu), or (\nu+1) in the two normalized bands, their
first surviving translated exponent is (a_{\min}=\ell+\nu_{\rm eff}) and

\[
h-\bigl(17\ell-12a_{\min}\bigr)>0.                 \tag{9ai}
\]

The minimum slack in (9ai) is `1`.  Hence every pivot direction satisfies
the terminal jet relations strictly below the terminal supporting line as
well as the total-degree line.  This decides the backtracking fork:

The calculation is uniform in the family parameter.  On the all-`r` sparse
edge section put (m=2r-1) and

\[
A=x+e^rx^r,\qquad
C=-\frac{2r^2}{(r-1)(2r-1)}e^{2r-1}x^{2r-1}.      \tag{9ai1}
\]

Set (c=2r^2/(r-1)).  Then

\[
A'=1+re^rx^{r-1},\qquad -C'=ce^{2r-1}x^{2r-2}.
\]

After moving the inhomogeneous fifth-binomial correction to the forcing
side, the first-return linear map is (25(A'Q_5-C'P_5)).  Select P degrees
`1,...,r-1` and Q degrees `0,...,2r-1`.  Ordered by monomial degree, block
elimination gives the `(3r-1) x (3r-1)` determinant

\[
(-1)^{r-1}25^{3r-1}r
\left(\frac{2r^2}{r-1}\right)^{r-1}
e^{\,2r(r-1)+1}.                                  \tag{9ai2}
\]

Indeed, the low monomials of the Q columns form an identity block.  The P
columns of degrees `2,...,r-1` then give diagonal high monomials.  The
degree-`1` P column overlaps the last low Q pivot; eliminating it creates
the last high entry (-r*c^{r-1}*e^{2r(r-1)+1}).  The column permutation and
the common factor (25^{3r-1}) give precisely the sign and scalar in (9ai2).

For each of orders `6,...,10`, select P degrees `3,...,r+1` and Q degrees
`0,...,2r`.  The common Bezout map (A'Q-C'P) has determinant

\[
(-1)^{r-1}
\left(\frac{2r^2}{r-1}\right)^{r-1}
e^{(r-1)(2r-1)}.                                  \tag{9ai3}
\]

Here the Q low monomials again form an identity block, while all selected P
monomials lie in the complementary high rows and have diagonal coefficient
(c*e^{2r-1}); moving the P block past the Q block contributes
((-1)^{r-1}).  This proves (9ai3) directly for symbolic (r), while the
checker also reconstructs the exact matrices for `r=2,...,8`.

Both are units over (\mathbf Q[e,e^{-1}]) for every (r\ge2).  The source
lift argument (9ag)--(9ai) applies to these columns using the family terminal
halfspace

\[
(7r-4)\ell-(5r-3)a\le h,
\]

where (h=r) on P and (h=2r-1) on Q.  There are

\[
(3r-1)+5(3r)=18r-1
\]

lifts through order `10`.  Their maximum source total degrees are
`13r+8` on P and `26r-5` on Q, leaving gaps `12r-8` and `24r-20` from the
degree bounds `25r` and `25(2r-1)`.  Their minimum terminal slack is again
one.  The only normalized selected columns are the order-`7` Q direction of
degree (r+1) on layer (5r-2), and the order-`9` Q direction of degree `2` on
layer `1`; multiplying each by (w-1) gives the required follower lift.

> **Source-lift conclusion.**  The surjectivity through `v^10` is supplied
> by genuine polynomial source directions inside the two-edge corner
> envelope.  Uniformly for every `r>=2`, the sparse-edge staircase is
> surjective through the second Kummer return.  No refinement using only
> those two certified edges can delete these directions.  A new lower Newton
> edge, not presently implied by the corner chain, could still cut some of
> them.

It is therefore useful to quotient the controlled local block exactly.  Put

\[
B=\mathbf Q[w_0,w_0^{-1},(w_0-1)^{-1}].
\]

Confluent evaluation at (0,1,w_0), with target multiplicities `(2,5,5)`,
has determinant

\[
2^{10}3^4w_0^{10}(w_0-1)^{25},                    \tag{9aj}
\]

and with layer-zero multiplicities `(3,6,6)` it has determinant

\[
2^{17}3^6 5^2w_0^{18}(w_0-1)^{36}.                \tag{9ak}
\]

Both are units in (B), so these are exact Chinese-remainder decompositions,
not dimension counts.  On target layer four, after removing (t^4), write the
local character-zero remainder as (G(w)).  Its coordinates (G(0)) and
(G'(0)) occur at (u)-orders `0` and `5`, hence at source transverse orders
`1` and `6`.  The edge Bezout identity controls the first, while the
`9 x 12` unit block (9af) controls the second.  Notice that the second target
coordinate is the order-`6` target return; the order-`5` packet itself lies
on Laurent layers divisible by five.

For layer zero, (J_0=dH/du) and, for
(H=h_0+h_1w+h_2w^2+\cdots),

\[
[u^4]J_0=5h_1,\qquad [u^9]J_0=10h_2.              \tag{9al}
\]

These are controlled by the surjective packets at source orders `5` and
`10`.  Triangularity means that all four removals are simultaneous with the
earlier edge equations.  Removing them changes each rank-`14` obstruction
to rank `12`.  Define the rank-`24` coordinate module over (B) by

\[
\begin{aligned}
\mathcal M_B={}&
 B[\epsilon_1]/(\epsilon_1^5)
 \oplus B[\epsilon_{w_0}]/(\epsilon_{w_0}^5)
 \oplus B^2 \\
&\oplus B[\eta_1]/(\eta_1^6)
 \oplus B[\eta_{w_0}]/(\eta_{w_0}^6),             \tag{9am}
\end{aligned}
\]

Over the descent-eight candidate algebra

\[
\mathcal A_8=B[y]/(27y^2-9y+1),
\]

the smallest presently justified global residual module is
(\mathcal M_{\rm glob}=\mathcal M_B\otimes_B\mathcal A_8).  It has rank
`24` over (\mathcal A_8), equivalently underlying rank `48` over (B).  The
two copies of (B) in (9am) are the terminal target residues (\rho_{20}) and
(\rho_{21}).  This is not yet an exclusion: the earlier triangular Laurent
solutions must still be imposed.  The fixed endpoint admits one further
exact elimination.

#### The complete fixed-endpoint elimination

Put (\delta=w-1).  The corner-derived band factors show that only `22`
band pairs can reach target orders (\delta^0,\ldots,\delta^4), and only `25`
nonzero-weight pairs can reach the six layer-zero coordinates through
(\delta^5).  Their exact local dependency cone uses `81` Taylor coordinates
on each side.  Three values are already normalized, so it contains `159`
free coordinates.  No lower band is silently assumed: closing the Laurent
dependencies requires precisely layers `40` down to `3`, with

\[
-22\le \ell_P\le15,
\qquad
-12\le \ell_Q\le25.                              \tag{9an}
\]

The new endpoints in (9an) enter together on layer `3`: (Q_{-12}) first
pairs with (P_{15}), while (P_{-22}) pairs with (Q_{25}).

For the elimination itself use the four Taylor coefficients of (K_{P_3})
in orders `1..4` and the six coefficients of (K_{P_{-1}}) in orders `0..5`.
They occur in the endpoint rows only through the target pairs

\[
(P_3,Q_1),\qquad(P_{-1},Q_5),
\]

and the layer-zero pairs

\[
(P_{-1},Q_1),\qquad(P_3,Q_{-3}).
\]

The normalized values (K_{P_3}(1)=1/5) and (K_{Q_1}(1)=-1) make the target
(\delta^0) coordinate identically `1`.  On the other four target rows and
six layer-zero rows, the coefficient matrix of the ten selected Taylor
variables is square and has

\[
\boxed{\det M_{w=1}=75000=2^3\cdot3\cdot5^5.}     \tag{9ao}
\]

This determinant is independent of every higher Taylor coefficient of
(Q_1,Q_5,Q_{-3}) and of all forcing from the other band pairs.  Thus (9ao)
is an affine-linear ideal elimination, not a generic-rank assertion.

It also respects the already removed (w=0) coordinates.  For a selected
endpoint order (j), use the degree-seven follower

\[
L_j(\delta)=\delta^j+(j-7)(-1)^j\delta^6
                    +(j-6)(-1)^j\delta^7.          \tag{9ap}
\]

The polynomial (L_j) has the prescribed Taylor coefficients through order
five and is divisible by (w^2=(1+\delta)^2).  All ten followers lie inside
the exact degree bounds `11` for (K_{P_3}) and `12` for (K_{P_{-1}}).
Consequently they preserve both controlled jets at (w=0).

All eleven (w=1) coordinates are therefore discharged: one is an identity
and ten solve for source variables.  The residual global coordinate module
has rank

\[
5+2+6=13.                                           \tag{9aq}
\]

over (\mathcal A_8): five target jets at (w=w_0), the two triangular target
residues, and six layer-zero coordinates at (w=w_0).  The solved variables
have now been substituted into the earlier Keller equations on layers
`40..3` and into these thirteen rows, as described next.

#### Carrying the substitution into the Laurent/Fitting system

The substitution is most compact in a straight-line presentation; expanding
it does not produce a polynomial ring on thirteen variables.  Keeping all
`22` target and `25` layer-zero forcing pairs, and writing (H(0)) as its exact
quadratic source-band circuit, the complete forcing vector (b) gives

\[
 x=-M_{w=1}^{-1}b.                                \tag{9aq'}
\]

Direct polynomial arithmetic verifies (M_{w=1}x+b=0).  The ten entries of
(x) have `1,489` straight-line terms.  If (H(0)) is counted as one node their
degrees are

\[
(2,3,4,5,1,2,3,4,5,6);
\]

after its quadratic source definition is expanded, their degrees are bounded
by

\[
(3,4,5,6,2,3,4,5,6,7).
\]

Consequently the carried bracket/Fitting circuits have total degree at most
eight.  This is an exact substitution identity, not a degree-only surrogate.

After removing the common monomial (t^L u^{\chi(L)}), every retained Laurent
row is one polynomial in (w).  The top row `40` vanishes identically on the
common-power top.  Layers `39..5` and `3` have `1,172` structural coefficient
slots: `353` on layers `39..29` are unchanged by (9aq'), while `819` on
layers `28..5` and `3` acquire the pivot circuits.  The target row is not
duplicated among those zero rows.  Instead, for

\[
G_{\rm red}=J_{4,\rm red}-1,
\qquad
D=w^2(w-1)^5(w-w_0)^5,
\]

the five Taylor coefficients at (w_0) are retained, and the quotient
(f=\operatorname{quo}_w(G_{\rm red},D)) has the two triangular residues
(\rho_{20},\rho_{21}) from the operator (N) above.  The six layer-zero rows
are

\[
H_{\rm red}(w_0)-H_{\rm red}(0),\quad
H_{\rm red}^{(j)}(w_0)\ (1\le j\le5).             \tag{9aq''}
\]

There are `1,061` active source coordinates before further triangular power
elimination.  In particular, the `30` coordinates in (P_{-21},Q_{-11})
cannot be deleted merely because their target image is killed: they re-enter
layer `3` against (Q_{24}) and (P_{14}).  The descent-eight component itself
is imposed by

\[
\operatorname{ord}_{w_0}K_{P_7}\ge3,qquad K_{Q_1}(w_0)=0,
\]

and, if (a=[(w-w_0)^3]K_{P_7}), by the exact globally scaled ratio

\[
K_{P_{-1}}(w_0)
=25^3w_0(w_0-1)^6y a^2.                           \tag{9aq'''}
\]

Here the factors in (9aq''') follow from
(C_0\sim t_0^5(w-w_0)^2/25),
(p_8\sim t_0^7u_0^3(w_0-1)^3a(w-w_0)^3), and
(u_0^5=w_0); they are not an untracked local rescaling.

The endpoint-disjoint part of the power elimination can also be completed.
For descents (1\le\delta\le11), put (\ell=25-\delta) and write the new
Q-band as
(t^\ell u^{j_0}(w-1)^\nu S(w)).  After a forced common factor is removed,
the new-Q operator is tridiagonal:

\[
N_\delta(w^k)=A_k w^k+B_k w^{k+1}+C_k w^{k+2},
\]

\[
\begin{aligned}
A_k&=w_0(5k+j_0),\\
B_k&=(w_0+1)(2\ell-5k-j_0)-5\nu w_0,\\
C_k&=5k+j_0+5\nu-4\ell.
\end{aligned}                                      \tag{9aq''''}
\]

Off resonance, rows and columns (0,\ldots,d) have determinant

\[
w_0^{d+1}\prod_{k=0}^{d}(5k+j_0),
\]

and at (\delta=5,10), deleting the one centralizer direction gives the unit
minor (w_0^d5^d d!).  Thus `134` of the `136` new-Q coordinates on layers
`39..29` eliminate over (B), two source centralizers remain, and the block
leaves `219` Fitting-coordinate slots.  The full carried presentation then
has `927` active source coordinates and `1,056` displayed generator slots
(neither count is asserted minimal).

This identifies the real coupling boundary.  At descent `12`, layer `28`,
the new pair is (P_3,Q_{13}); from there down, substituting (9aq') changes the
new-Q operator itself.  Subtracting all later uncoupled operator ranks would
therefore be invalid.  The next exact object is the coupled Schur/Fitting
circuit on descents `12..37`, followed by the thirteen functionals
(9aq'')--not another independent endpoint rank count.  No unit ideal and no
counterexample is obtained at this stage.

The terminal toroidal coordinates clarify the repeated numeral five.  If

\[
s=t^{17}z^{12},\qquad q=t^{10}z^7,
\]

then (t=s^{-7}q^{12}) and

\[
w-1=(1+t)^5-1=5s^{-7}q^{12}+O(q^{24}).             \tag{9ar}
\]

Thus (w=1) is transverse to the whole terminal divisor.  The terminal
residue five-cycle instead acts in the independent (s)-cover at (s=-1).
Their exact order-five normalizers and arithmetic characters are compared in
[`F2_A6_SIMPLE_SPECTATOR_GLUING.md`](F2_A6_SIMPLE_SPECTATOR_GLUING.md); only
an additional toroidal node theorem can identify the two local inertia
groups.

### 4.4 The complete lower B0 tail and why (E_5=0) is invalid

The same corner inequalities used above do not stop at layer `4`.  Without
assuming a lower polygon, they give the complete necessary band intervals

\[
-75\le\ell_P\le15,
\qquad
-125\le\ell_Q\le25. \tag{10}
\]

Applying the same full-rank binomial jets on every band gives:

| side | bands | source coefficients before jets | jet equations | dimension |
| --- | ---: | ---: | ---: | ---: |
| P | 91 | 706 | 53 | 653 |
| Q | 151 | 1,901 | 136 | 1,765 |

Thus the complete B0 envelope has `2,418` linear parameters, or `2,413`
after the five endpoint normalizations.  Its bracket layers run from `40`
down to `-200`.  Apart from target layer `4`, all `240` layers are zero;
`204` of them lie below the target.  The exact compressed record contains
`13,741` band-pair incidences, `1,327,026` nonzero source-basis generators,
and `41,685` raw scalar coefficient rows.  The inexpensive post-jet
support-row upper bound is `41,388`; the jet equations plus compressed
generators, rather than that upper bound, are the exact system.

This changes the interpretation of the descent-eight resultant.  Its fifth
multiple is descent `40`, bracket layer `0`.  At that row the corner envelope
contains new bands

\[
p_{-25}=t^{-25}u^{25}A(u^5),\quad \deg A\le8
\quad(\dim=9),
\]

\[
q_{-15}=t^{-15}u^{15}S(u^5),\quad \deg S\le18
\quad(\dim=19). \tag{11}
\]

The target-tail cross pair (p_{11},q_{-11}) also contributes to layer zero.
Accordingly the relevant linear operator is

\[
T_{40}(q)=5C_0q'+15C_0'q
=5C_0^{-2}(C_0^3q)'. \tag{12}
\]

The primitive numerator (E_5) was valid only when the fifth multiple lay
strictly before the target and no noncentral tail had entered.  For descent
eight, imposing (E_5=0) would silently delete the genuine bands (11).  The
correct next condition is the layer-zero maximal-minor/Fitting condition for
(12), with the target row and every lower source band retained.

That Fitting condition now has an exact small form.  On layer zero every pair
is `(P_ell,Q_-ell)`, and hence

\[
J_0=\frac d{dt}H_{40},\qquad
H_{40}=\sum_\ell \ell P_\ell Q_{-\ell}. \tag{12a}
\]

Every product in (12a) is a polynomial in `w` of degree at most `33`.  On
the double-root stratum set

\[
K=w^3(w-1)^6(w-w_0)^6. \tag{12b}
\]

Up to a unit, the new `P_-25/Q_-15` columns contribute exactly

\[
15K\left(S_{-15}+3w^2(w-1)^4R^2A_{-25}\right),
\qquad \deg(\cdots)\le18.
\]

Since (12a) only requires `H_40` to be constant, the complete first
post-target obstruction is

\[
H_{<40}\in B+K B[w]_{\le18}. \tag{12c}
\]

Equivalently, it is the residue of `H_<40` in

\[
\mathcal A_{40}/B,
\qquad
\mathcal A_{40}=B[w]/(K).
\]

The Artinian algebra has length `15`, with local lengths `(3,6,6)` at
`0,1,w0`; quotienting by constants leaves a rank-`14` obstruction module.
An explicit Hermite form of its fourteen equations is

- `H'(0)=H''(0)=0`;
- `H(1)=H(0)` and `H^(j)(1)=0` for `1<=j<=5`; and
- `H(w0)=H(0)` and `H^(j)(w0)=0` for `1<=j<=5`.

The old numerator `E5` is recovered precisely as the first movable-root
coordinate of the pure eight-step staircase:

\[
[w-w_0]H_{<40}=-\frac{a}{27}E_5. \tag{12d}
\]

Off-grid pairs can change this coordinate.  Indeed, before imposing the
earlier triangular rows, the `ell=-1` and `ell=1` source-basis products span
all fourteen nonconstant residue classes: their factors are respectively
`w(w-w0)` and `w(w-1)`, and their difference is `(1-w0)w`.  Therefore layer
zero alone cannot exclude the branch.  Equations (9z)--(9af) perform the
first `u=0` return packet and show that it has zero cokernel.  After (9ao),
the remaining exact calculation is to substitute its ten pivot solutions
together with the other solutions on layers `40..3` into the five target
jets and six layer-zero coordinates at (w=w_0), and into the two target
residues.  The lower tail is retained; the fixed (w=1) block no longer needs
to be carried as an independent quotient.

There is an important remaining branch ledger.  Descent `8` is the
**earliest** survivor, not the only possible position of the first non-root
defect.  The full P envelope permits relative descents through `90`.  Counting
only nonlinear multiples strictly before target descent `36` gives:

| first-defect spacing | primitive zero multiples before target | status |
| --- | --- | --- |
| `1..7` | through at least the fifth | excluded by (9e) |
| `8` | `2,3,4` | earliest residual, reduced to (9f)--(9g) |
| `9..11` | `2,3` | open target/tail branches |
| `12..17` | `2` | open target/tail branches |
| `18..90` | none | no nonlinear multiple precedes the target |

The machine record enumerates every spacing separately.  Therefore killing
(9g) would remove the first escape but would still not prove F2 exclusion;
the later regimes must either be eliminated by the full Fitting descent or
shown to be absorbed by an exact polynomial common root.

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
the exact degree-six target row.  The full lower coefficient ledger is now
compiled, but its nonlinear triangular elimination is stopped at the
earliest movable-double-`R` Fitting branch and the later-spacing regimes in
the table above.  The earliest coefficient escape lies on the boundary
two-principal-packet stratum, which gives a precise first meeting point for
the coefficient and gluing routes.
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

1. first specialize the complete `2,413`-variable normalized B0 system to the
   nonzero double-root stratum (9g) and one of the two roots of (9f), then
   process first-defect spacings `9..90` if that earliest stratum is empty;
2. triangularly compile the upper zero-row solutions into the explicit
   fourteen-coordinate target cokernel (9j)--(9k), then reduce the result in
   the length-15 layer-zero algebra (12b), while retaining all bands through
   layer `-200`; and
3. either obtain a unit Fitting ideal, or reconstruct a global source pair
   and then test its nonvanishing vertices and exact Newton row.  Consistency
   of the B0 over-envelope alone would not be a counterexample because it may
   use supports absent from the actual lower polygon.

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
calculation eliminates every branch of the naive projection.  The full
source-band compiler here explains exactly why its compact fifth residue
cannot be promoted to an exclusion: at descent `40` the target correction
and the genuine lower Laurent bands have already entered.
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
.venv/bin/python plane-jc/cas/reduce_f2_75_125_endpoint_system.py
.venv/bin/python plane-jc/cas/test_f2_75_125_frontend.py
.venv/bin/python plane-jc/cas/f2_75_125_frontend.py
```

Intentional regeneration of
[`../artifacts/generated-results/jc2_f2_75_125_character_layers.json`](../artifacts/generated-results/jc2_f2_75_125_character_layers.json)
uses `--refresh`.  The complete `v9` artifact has whole-file SHA-256
`96e4fd2ff853fcba9d41a72973ea55a1acb2cd41eefad21e69efd6ae73df8b8b`.
The carried endpoint artifact
[`../artifacts/generated-results/jc2_f2_75_125_endpoint_reduction.json`](../artifacts/generated-results/jc2_f2_75_125_endpoint_reduction.json)
has whole-file SHA-256
`9834ed2ba4e64a2b034a83cd0604140206f1a8192a561509c208d02e0a0ca189`.
The final command prints the machine-readable partial front-end certificate.
