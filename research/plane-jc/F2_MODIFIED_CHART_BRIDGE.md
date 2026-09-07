# The F2 modified-chart bridge and projected top-band obstruction

> **Status.**  The F2 corner chain fixes `gamma=2`, hence the historical
> `d=3` chart, and derives all finite support envelopes used below.  Retaining
> the translated binomial-jet relations gives exact source-image ranks, not a
> box of independent coefficients.  In the literal nonnegative-`xi`
> projection, the complete `r=3` top-band ideal is the unit ideal; therefore
> every projected polynomial branch is excluded.  This does **not** yet
> exclude `(75,125)`: the full Laurent powers have negative `xi` tails, and
> the calculation proves that simply discarding those tails cannot preserve
> the top bracket.  Any valid modified-Laurent theorem must derive the exact
> tail correction rather than assume a bracket-preserving truncation.  The
> exact tangent audit further corrects the old divisible-slice descent:
> layers 39 through 36 do not force a polynomial common-root continuation.

The exact replay is
[`cas/verify_f2_modified_chart_bridge.py`](cas/verify_f2_modified_chart_bridge.py).
Its pinned output is
[`../artifacts/generated-results/jc2_f2_modified_chart_bridge.json`](../artifacts/generated-results/jc2_f2_modified_chart_bridge.json),
with SHA-256
`ac7dbc170cafbcf028079b9ccdb41afd78c333e3850abc0761f64cc056e7d7b8`.

## 1. The chain fixes `gamma=2`

Put

\[
 m=2r-1.
\]

The complete F2 chain starts with

\[
 A_0=(5,20),\qquad A'_0=(1,0),\qquad A_1=(7/5,2).
\]

For the selected nonzero root, the generated-corner identity is

\[
 A_1=A'_0+\gamma(1/5,1).
\]

Comparing either coordinate gives

\[
 \boxed{\gamma=2}. \tag{1}
\]

Thus the formal `d=2` system is a useful historical regression, but it is
not a second F2 branch.  In the historical complementary notation,

\[
 d=5-\gamma=3. \tag{2}
\]

## 2. Direct monomial bridge

Let `(X,Y)` be the certified coordinates after `x=X^5` and the selected
Puiseux translation.  Put

\[
 X=\xi y^2,\qquad Y=\xi^{-2}y^{-7}. \tag{3}
\]

The exponent map and coordinate Jacobian are

\[
 (a,b)\longmapsto(a-2b,2a-7b), \tag{4}
\]

\[
 \frac{\partial(X,Y)}{\partial(\xi,y)}=-3\xi^{-2}y^{-6}. \tag{5}
\]

The certified bracket `X^4` therefore becomes

\[
 \boxed{[P,Q]_{\xi,y}=-3\xi^2y^2}. \tag{6}
\]

The forced terminal vertices map as follows:

\[
\begin{array}{c|cc}
 &\text{certified }(X,Y)&\text{modified }(\xi,y)\\ \hline
P&(7r,2r)&(3r,0)\\
P&(4,1)&(2,1)\\
Q&(7m,2m)&(3m,0)\\
Q&(7r-3,m)&(3r-1,1)\\
Q&(1,0)&(1,2).
\end{array} \tag{7}
\]

Hence the candidate polynomial degrees are consequences of the chain:

\[
 \deg_\xi P=3r,\qquad \deg_\xi Q=3m. \tag{8}
\]

With `t=XY` and `z=Y^(-1)`, the full common root transforms by

\[
 t=\xi^{-1}y^{-5},\qquad z=\xi^2y^7,
\]

\[
 \boxed{t^7H(t)z^5=\xi^3H(\xi^{-1}y^{-5})}. \tag{9}
\]

In particular, its negative `xi` tail is intrinsic: for F2 the source-band
argument gives `deg H=18`, so (9) ranges from `xi^3` through `xi^-15`.

## 3. Corner-derived supports and their exact linear image

An original source monomial on band `ell=5i-j` contributes

\[
 t^\ell(1+t)^jz^\ell.
\]

If `k` is a surviving `t` exponent, (3) sends it to

\[
 (\xi\text{-exponent},y\text{-exponent})
 =(2\ell-k,7\ell-5k). \tag{10}
\]

The total-degree and terminal-halfspace inequalities give

\[
 k\ge
 \max\left(
 \ell,
 \left\lceil
 \frac{(7r-4)\ell-h}{5r-3}
 \right\rceil
 \right), \tag{11}
\]

where `h=r` on P and `h=m` on Q.  Equations (10)--(11) determine all
possible support positions.

They do not make those positions independent.  On a fixed band the source
coefficients form

\[
 t^\ell\sum_j a_j(1+t)^j. \tag{12}
\]

The rows below the bound (11) vanish, and the surviving rows are the image
of the kernel of that binomial-jet matrix.  The checker constructs this
kernel and its projection over `QQ` exactly.

For `r=3` the literal nonnegative-`xi` images are

\[
\begin{array}{c|ccc}
 &\text{possible positions}&\text{image rank}&\text{linear relations}\\ \hline
P&83&74&9\\
Q&215&196&19.
\end{array} \tag{13}
\]

The leading-minus-one support is `{-2,-5}`, uniformly in the checked family
samples, so the normalized translation is forced to have

\[
 \xi=x-G(y),\qquad G=g_{-2}y^{-2}+g_{-5}y^{-5}. \tag{14}
\]

After (14), the familiar `r=3` P mask has 80 possible positions and the Q
mask has 212.  These are support-box sizes, not scalar dimensions of the
source image.

The distinction is already visible on the P top diagonal.  Write

\[
 c_i=[\xi^{9-i}y^{-5i}]P_+\qquad(0\le i\le9), \tag{15}
\]

where `P_+` denotes literal retention of nonnegative `xi` powers.  Every
source image satisfies

\[
\begin{aligned}
0={}&3470256c_0-1118621c_1+299936c_2-65511c_3\\
&+11250c_4-1430c_5+120c_6-5c_7. \tag{16}
\end{aligned}
\]

The formal terminal point has `c_0=1` and all other top-diagonal
coefficients zero, so (16) evaluates to `3470256`.  Thus the old terminal
resonance is not merely missing an endpoint: it is outside the exact source
image.

## 4. The formal terminal resonance

For every `r>=2`, the ambient support box contains

\[
\begin{aligned}
 P_T&=x^{3r}+a x^2y,\\
 Q_T&=x^{3m}+\frac mr a x^{3r-1}y
 +\frac{m(r-1)}{2r^2}a^2xy^2.
\end{aligned} \tag{17}
\]

Direct differentiation gives

\[
 [P_T,Q_T]
 =\frac{3(r-1)m}{2r^2}a^3x^2y^2. \tag{18}
\]

This exactly explains the terminal binomial section of the formal
modified-series compiler.  It is not a plane candidate.  For `r=3`, (16)
excludes it from the literal source projection before any nonlinear tangent
analysis.

The former cubic point with bracket `(35/9)y^6x^2` is excluded even earlier:
the forced F2 weight in (6), after (14), is proportional to
`y^2(x-G)^2`, not `y^6x^2`.

## 5. Ambient tangent calculation, correctly interpreted

For comparison with the older residue calculation, put

\[
 z=y^{-3},\qquad x=yX,\qquad
 P=y^9p(X,z),\qquad Q=y^{15}q(X,z).
\]

The bracket equation becomes

\[
3\bigl(5p_Xq-3pq_X-z(p_Xq_z-p_zq_X)\bigr)
=\mu z^6(X-g_2z-g_5z^2)^2. \tag{19}
\]

Linearizing all 80 support-box coordinates and the three target parameters
at the terminal point gives a matrix of shape `89 x 83`, rank `77`, and a
six-dimensional kernel:

\[
a_{4,2},\quad a_{7,3},\quad 6a_{2,1}+a_{9,3},\quad
\tfrac15a_{7,2}+\mu,\quad -2a_{8,3}+g_2,\quad a_{8,4}+g_5. \tag{20}
\]

Pulling these modes back through (14) and imposing the exact source-image
relations leaves only

\[
a_{4,2},\qquad a_{7,3},\qquad \tfrac15a_{7,2}+\mu. \tag{21}
\]

This does not define a three-dimensional source tangent component: the base
point itself violates (16).  Equations (20)--(21) are an audit of the old
support-box over-approximation, not a surviving branch classification.

## 6. The top-diagonal unit ideal

The preceding single relation kills the terminal point.  The complete top
diagonal kills every branch in the bracket-preserving polynomial projection.

For an exponent `n` (equal to `r` on P and `m` on Q), write

\[
 P_{+,\mathrm{top}}=\xi^{3n}B_n(s),
 \qquad s=\xi^{-1}y^{-5},
 \qquad \deg B_n\le3n. \tag{22}
\]

Put

\[
 \Phi(t)=\frac{(1+t)^5-1}{t}.
\]

The exact binomial-jet image in (12) is equivalent to

\[
 B_n(t)\equiv
 \Phi(t)^{2n}R_n((1+t)^5)\pmod {t^{3n+1}},
 \qquad \deg R_n\le2n. \tag{23}
\]

It has dimension `2n+1` inside the `3n+1` top-diagonal positions, hence
exactly `n` linear relations.

Now make the conditional cutting hypothesis that the polynomial projection
has the target bracket.  Its highest diagonal bracket must vanish.  If B and
D denote the P and Q top polynomials, this gives

\[
 mB'D-rBD'=0. \tag{24}
\]

Because `gcd(r,m)=1` and `B(0)=D(0)=1`, (24) implies

\[
 B=H^r,\qquad D=H^m,
 \qquad H=1+c_1t+c_2t^2+c_3t^3. \tag{25}
\]

Set

\[
 v=(1+t)^5-1,qquad
 \tau=(1+v)^{1/5}-1,
\]

\[
 L(v)=H(\tau)\left(\frac{\tau}{v}\right)^2. \tag{26}
\]

Equation (23) for the exponent `n` is exactly the gap system

\[
 [v^q]L(v)^n=0,
 \qquad 2n+1\le q\le3n. \tag{27}
\]

For `n=3`, the checker compares these three equations directly with the
three left-kernel relations of the `10 x 7` source matrix.  Their coefficient
spans agree, with change-of-basis determinant `5^(-42)`.

Let `I_r` contain (27) for `n=r` and `n=2r-1`.  Exact rational Gröbner
calculations give:

- `r=2`: the five equations have Gröbner basis `(1)`;
- `r=3`: the three P equations define an Artinian algebra of length `27`,
  with leading ideal `(c1,c2,c3^27)`;
- in that length-27 algebra, the first Q equation `[v^11]L^5` is already a
  unit.  Its residue has degree `26`, the P eliminant has degree `27`, and
  their gcd is one;
- their nonzero resultant, equivalently the determinant/norm of
  multiplication by that Q residue, has numerator/denominator lengths
  `505/700` digits and SHA-256
  `cd610d23ca92bc10c7788b72a603d3f71353d36f107bb87695669f4decb27314`;
- `r=4`: the four P equations alone have Gröbner basis `(1)`.
- `r=5`: the five P equations alone also have Gröbner basis `(1)` in an
  optional exact `QQ` Singular replay.

Therefore

\[
 \boxed{I_3=(1)}. \tag{28}
\]

There is no tiny residual `r=3` plane candidate in this projection: all of
its top-diagonal branches are already gone.  The Artinian length `27` is the
intermediate P-only scheme; the first Q gap has invertible Fitting
determinant and removes all 27 points at once.

## 7. What remains after the elimination

Equation (28) is also a warning about the desired theorem.  The full common
root (9) has powers whose complete Laurent top brackets cancel.  Their
positive parts cannot have cancelling top bracket, by (28).  Thus a valid
family-level modified-Laurent result cannot be the statement “discard all
negative `xi` powers and retain the bracket.”

The exact full-Laurent seed is nevertheless tiny.  Put

\[
 v=(1+s)^5-1,\qquad
 \Phi(s)=v/s.
\]

The source invariance and `gcd(r,2r-1)=1` argument give

\[
 \boxed{
 H(s)=\Phi(s)^2\left(\frac1{25}+\rho_1v+\rho_2v^2\right)
 },\qquad \rho_2\ne0. \tag{29}
\]

Indeed, `s^2H(s)=v^2(1/25+rho1*v+rho2*v^2)` is a polynomial in
`(1+s)^5`.  The nonzero endpoint is exactly

\[
 \rho_2\xi^{-15}y^{-90}. \tag{30}
\]

The complete seed is the monomial-chart image of the honest degree-25 source
polynomial

\[
x(xy^5-1)^2
\left(\frac1{25}+\rho_1(xy^5-1)+\rho_2(xy^5-1)^2\right). \tag{30a}
\]

This explains why its negative tail cannot simply be deleted.  It also exposes
a correction to the former upper-descent claim.  Exact source bands do not
force the divisibilities used in that calculation.  At `r=3`, the first-five
homogeneous kernel dimensions are `6,6,7,7,10`: every P-band direction has the
Q follower `q=-3*C0^2*p`, and the extra descent-five mode is the ordinary
commuting `C0^4` term.  The formal `lambda*C0^(-1)` resonance occurs only at
descent `30`, layer `10`, and is not an independent homogeneous source-band
mode.  It can appear only through nonlinear cancellation with the `F` tail.
The exact rank/Fitting statement is recorded in
[`F2_75_125_DERIVATION.md`](F2_75_125_DERIVATION.md).

The first omitted positive/negative cross term is explicit.  If

\[
B=\operatorname{trunc}_{\le3r}H^r,
\qquad
D=\operatorname{trunc}_{\le3m}H^m,
\]

then comparison with the zero Wronskian of the full powers gives

\[
[s^{3r}](mB'D-rBD')
=-m(3r+1)[s^{3r+1}]H^r. \tag{31}
\]

The checker verifies (31) in the `r=3` two-parameter source seed; its right
side is a nonzero polynomial in `rho1,rho2`.  Equation (31) is the first
piece of the required tail correction.  Vanishing of this one polynomial is
not enough—the unit-ideal calculation uses the whole projected Wronskian.

The remaining exact problem is narrower and different:

> **Tail-correction lemma.**  Derive, from the full Laurent pair and the
> corner chain, the correction contributed by positive-negative Laurent
> cross terms (and by `lambda*C^(-1)+F`) to the polynomial modified pair.
> Prove its complete support and show whether it can satisfy the lower-band
> polynomiality and endpoint conditions.

The companion Laurent-layer classifier now supplies the complete necessary
source-band envelope through bracket layer `-200`, so “complete support” here
no longer means guessing missing source bands.  What remains unproved is the
nonlinear pushforward of those exact bands into this modified `xi` chart and
the resulting specialized Fitting elimination.  In particular, the earliest
descent-eight defect survives its first local target jet; its fifth multiple
already contains genuine lower-tail variables and cannot be tested by setting
the primitive fractional-power numerator `E5` to zero.

The companion classifier now makes the first two Fitting rows explicit in
the original Laurent chart.  The full target has a rank-fourteen cokernel,
and layer zero is the rank-fourteen quotient of the length-fifteen algebra
`B[w]/(w^3*(w-1)^6*(w-w0)^6)` by constants.  Old B0 source-basis generators
span both cokernels before their earlier Keller equations are imposed.  Thus
the remaining bridge is specifically the simultaneous triangular
specialization of those generators, not the construction of either
cokernel.  The lowest-`u` edge escape has also been pushed one order farther:
its uniform all-`r` four-term witness is the first jet of an exact formal
shear, but the shear is necessarily infinite, and keeping P linear in the
transverse coordinate forces a nonpolynomial Q coefficient at order two.
The complete polynomial order-two repair is now classified and cannot
terminate quadratically.  At `r=3`, cubic and quartic termination ideals are
also `(1)`.  The exact tail correction at `v^5` is now derived in
[`F2_75_125_DERIVATION.md`](F2_75_125_DERIVATION.md).  Its fifth-binomial
wrap makes the first-return matrix surjective over both descent-eight
branches, and the edge recursion remains surjective through `v^10`.
Therefore (28) remains a conditional projection exclusion, not an exclusion
of `(75,125)` or of the infinite F2 family.  The exact ten-variable fixed-
`w=1` substitution and the endpoint-disjoint power rows through layer `29`
have now been carried.  The missing information is localized to the coupled
Schur/Fitting rows from layer `28` through layer `3` and the thirteen
remaining `w=w0`/residue functionals.  The lower Laurent tail remains
retained.
The exact `r=4,5` samples suggest a stable family obstruction, but no
induction from these samples to all `r` has been proved.

## Reproduction

```bash
.venv/bin/python plane-jc/cas/verify_f2_modified_chart_bridge.py
# Optional exact family sample (requires Singular):
.venv/bin/python plane-jc/cas/verify_f2_modified_chart_bridge.py --extended-r5
# Intentional regeneration after review:
.venv/bin/python plane-jc/cas/verify_f2_modified_chart_bridge.py --refresh
```

Expected markers include:

```text
F2_MODIFIED_EXACT_SOURCE_PROJECTION_RANK_PASS
F2_MODIFIED_TERMINAL_OUTSIDE_SOURCE_IMAGE_PASS
F2_R3_PROJECTED_TOP_GAP_ARTINIAN_LENGTH=27
F2_R3_PROJECTED_TOP_GAP_RESULTANT_NONZERO_PASS
F2_R3_PROJECTED_TOP_GAP_UNIT_IDEAL_PASS
```
