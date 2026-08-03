# The F2 modified-Laurent family: exact residue and terminal theorem

> **Status.**  This note proves an all-parameter theorem for the two forced
> corner-chain edges and their terminal residue cover.  It also gives an
> exact all-parameter coefficient system **conditional** on a modified
> common-power chart.  The companion modified-chart bridge now proves that
> the F2 chain fixes `gamma=2`, hence the candidate value `d=3`, derives its
> `y^2(x-G)^2` bracket and finite Laurent support envelopes, and removes the
> `d=2` system as a non-F2 branch.  Retaining the exact binomial-jet links
> gives source-image ranks `74/83` for P and `196/215` for Q at `r=3`.
> Moreover, the literal nonnegative-`xi` polynomial projection has a
> top-diagonal unit ideal: all of its `r=3` branches are excluded.  This also
> proves that naive deletion of the negative Laurent tail cannot preserve the
> bracket.  The missing theorem must derive the tail correction, not assume
> a bracket-preserving truncation.  The corrected exact tangent calculation
> also shows that the formerly claimed first-four-layer common-root descent
> had assumed unproved divisibility: its true `r=3` kernel dimensions are
> `6,6,7,7,10`, and the formal `lambda*C^(-1)` resonance is at layer `10`,
> not layer `35`.  The nonlinear rows recover a source root through descent
> `7` and leave the exact descent-`8` local quadratic
> `27*y^2-9*y+1=0`.  Consequently this work does not yet
> exclude `(75,125)` or the infinite F2 family.  The formal residue
> calculation can now be pushed in the opposite direction: every all-`r`
> compact Fitting residue has an explicit endpoint-binomial point, and both
> historical-tail `r=3` residues have smooth positive-dimensional solutions
> with nonzero `lambda` and nonzero first visible `F` coefficient.  Thus the
> finite power equations alone cannot give the desired exclusion.  Two exact
> bracket gates now cut this geometry: the certified `X^4` weight excludes
> every `B(t) in K[t^d]` congruence section for `d=2,3`, and on the
> historical-tail `d=3` cubic section a monomial leading bracket coefficient
> leaves only a sharp `lambda=0`, exponent-`7*ell` ray.  The remaining gap is
> the chart/coordinate bridge, not a guessed support mask for `F_(-6)`.

The executable calculations are
[`cas/generate_f2_modified_system.py`](cas/generate_f2_modified_system.py),
[`cas/verify_f2_modified_chart_bridge.py`](cas/verify_f2_modified_chart_bridge.py),
[`cas/verify_f2_terminal_residue_cover.py`](cas/verify_f2_terminal_residue_cover.py),
and
[`cas/verify_f2_a6_simple_spectator_gluing.py`](cas/verify_f2_a6_simple_spectator_gluing.py).
The canonical generated systems are pinned in
[`../artifacts/generated-results/jc2_f2_modified_laurent_family.json`](../artifacts/generated-results/jc2_f2_modified_laurent_family.json).
Its SHA-256 is
`bca206498c153e41a2f31344015df2ce63890f8b12228b6d0ba2c0970eb87c85`;
the recorded software assumptions are `.python-version`, `requirements.txt`,
and exact characteristic-zero SymPy arithmetic.

## 1. Outcome and claim boundary

Put

\[
r=j+2,\qquad m=2r-1.
\]

The family row is the `F_2` row of the published table reproduced in
[`external/arxiv-1708.07936/Some_algorithms.tex`](external/arxiv-1708.07936/Some_algorithms.tex).

The F2 table and complete-chain theorem give the degree family

\[
\boxed{(\deg P,\deg Q)=(25r,25m)}. \tag{1}
\]

They also force two consecutive edges and a terminal degree-`2r` residue
cover with geometric monodromy `A_(2r)`.  These statements are unconditional
within the certified F2 corner-chain setup.

The separate polynomial ansatz

\[
P=C^r,\qquad Q=C^m+\lambda C^{-1}+F \tag{2}
\]

has an exact finite formal compiler.  The corner-derived bridge selects its
`d=3,h=2` F2 candidate but does not yet prove that the polynomial cut satisfies
(2).  Under (2), the compiler:

1. reproduces both published `r=2` systems;
2. constructs the `r=3,d=2` window with `14` coefficient functions;
3. constructs the `r=3,d=3` window with `22` coefficient functions;
4. eliminates all power equations triangularly;
5. presents the formal-series obstruction by an Artinian Fitting ideal;
6. gives closed formulas in `r,d`;
7. proves that every compact residue ideal is proper by an all-parameter
   endpoint-binomial section;
8. identifies a rational fourfold for `r=3,d=2` and a smooth formal
   coefficient-torus branch for `r=3,d=3`; and
9. proves a uniform congruence-support bracket gate, then removes the
   single-monomial `F_(-6)` assumption from the cubic-invariant gate; and
10. proves that the exact `r=3` projected top-band ideal is `(1)`, via a
    length-27 P algebra and a nonzero Q multiplication determinant.

The `d=2` system is retained as a published `r=2` regression and formal
comparison; the generated-corner formula shows that it is not an F2 branch
for `r=3` or for the F2 family.  The numbers `14` and `22` count coefficient functions of `y`, not scalar
variables.  The parameter `lambda` and the coefficients of `F` are additional
unknowns.  Turning these formal systems into finite scalar systems requires
proved Laurent-`y` masks and nonvanishing conditions.

In particular, the new residue geometry is a negative result for a tempting
strategy: neither the determinant in the `d=2` chart nor the Fitting ideal in
the `d=3` chart is an obstruction on the unrestricted coefficient space.

## 2. What the corner chain derives for every `r`

Let `X=x^(1/5)` after the forced Puiseux translation.  The chain forces the
actual vertices

\[
\begin{aligned}
P &: (25r,20r)\longrightarrow(7r,2r)\longrightarrow(4,1),\\
Q &: (25m,20m)\longrightarrow(7m,2m)\longrightarrow(1,0).
\end{aligned} \tag{3}
\]

The change of variables contributes

\[
[P,Q]_{X,y}=5X^4[P,Q]_{x,y}, \tag{4}
\]

so scalar normalization makes the bracket `X^4`.  All vertices in (3) have
nonzero coefficients.  Equation (3), however, is not a claim that the
displayed vertices are the complete Newton polygons.

### 2.1 Exact source-monomial envelope

For an original monomial `x^i y^j`, put

\[
\ell=5i-j.
\]

After the substitution and translation it contributes

\[
t^\ell(1+t)^jz^\ell,
\qquad t=Xy,\quad z=y^{-1}. \tag{5}
\]

The original total-degree bound `i+j<=D` is therefore exactly

\[
0\le i,\qquad 0\le5i-\ell,\qquad6i-\ell\le D, \tag{6}
\]

where `D=25r` for `P` and `D=25m` for `Q`.  The terminal supporting
halfspaces are

\[
ma-(7r-4)b\le r\quad(P), \tag{7}
\]

\[
ma-(7r-4)b\le m\quad(Q). \tag{8}
\]

Equations (6)--(8) are the exact corner-derived B0 jet envelope.  They do not
prove that every allowed monomial occurs, or exclude cancellation among
translated monomials.  Thus they are not exhaustive Laurent support masks.

### 2.2 Uniform common-power top band

The translated type-II edge is uniformly

\[
C_{\rm top}=t^7H(t)z^5,\qquad \deg H=18,
\]

\[
P_{\rm top}=C_{\rm top}^{,r},\qquad
Q_{\rm top}=-B_r C_{\rm top}^{,m},
\qquad
B_r=\frac{2r^2}{(r-1)(2r-1)}. \tag{9}
\]

Put `u=1+t` and

\[
A_0(u)=(u-1)^2H(u-1).
\]

Source-band membership gives

\[
A_0^r,A_0^m\in k[u^5]. \tag{10}
\]

Since `gcd(r,m)=1`, applying `u -> zeta*u` for a fifth root of unity shows
that `A_0(zeta*u)/A_0(u)` has both `r`-th and `m`-th power one.  It is therefore
one, and `A_0` itself belongs to `k[u^5]`.  Divisibility at `u=1` then yields

\[
H(t)=(1+u+u^2+u^3+u^4)^2R(u^5), \tag{11}
\]

where `R` is quadratic and `R(1)=1/25`.

There is an important correction to the upper bracket descent.  The formal
substitution

\[
p=C_0^{r-1}U,\qquad
q=-B_rC_0^{m-1}V \tag{12}
\]

does reduce a **divisible slice** to the ODE printed in the earlier version of
this note, but exact source bands do not force either divisibility in (12).
For arbitrary (p) in the P band (5r-\delta), the genuine unrestricted
tangent follower is

\[
\boxed{q=-B_r\frac mr C_0^{m-r}p.} \tag{12a}
\]

It lies in the exact Q band (5m-\delta): since (m-r=r-1), multiplication
by (C_0^{r-1}) adds exactly source degree (25(r-1)), terminal height
(r-1), and band (5(r-1)).  Substitution verifies that (12a) kills the
linearized bracket identically.  Hence every exact P-band direction survives;
the first four descents do not force continuation of a polynomial source
root.

The q-only kernel has the exact equation

\[
5C_0q'-(5m-\delta)C_0'q=0,
\qquad q^5=cC_0^{5m-\delta}. \tag{12b}
\]

Because (C_0=t^5((1+t)^5-1)^2R((1+t)^5)) contains a fixed factor of
multiplicity two, (12b) has a Laurent-polynomial source solution precisely at

\[
\delta=5j,\qquad 1\le j\le m,
\]

where it is generated by (C_0^{m-j}).  The later formal resonances

\[
\delta=5(m+k),\quad 1\le k\le r-1,
\qquad C_0^{-k}, \tag{12c}
\]

are not homogeneous source-band modes.  In particular, the modified-Laurent
term `lambda*C^(-1)` is at descent

\[
5(m+1)=10r,
\]

namely bracket layer (5(r-1)), not at descent five.  It can occur only via
nonlinear cancellation with the `F` tail.  For `r=3`, the corrected first-five
kernel dimensions are `6,6,7,7,10`; the extra layer-35 term is the ordinary
commuting (C_0^4), while `lambda*C_0^(-1)` is located at layer `10`.

Thus (9)--(11) remain a family-level top-band theorem, but no triangular
common-root continuation has been proved below it.  The exact nonlinear
handoff is instead a cokernel/Fitting condition for the operator in (12b).

For the first member `r=3`, the next nonlinear rows can now be evaluated.
They force any first non-root P defect at descent at most `7` back into the
source-root slice.  The first surviving local branch is at descent `8`, is
supported at multiplicity-two primes of `C_0`, and reduces to

\[
27y^2-9y+1=0,
\qquad \operatorname{disc}=-27. \tag{12d}
\]

This is the earliest exact tail-correction candidate over
`QQ(sqrt(-3))`; the target descent `36` intervenes before its fifth residue,
so (12d) is not yet a global source pair.

### 2.3 The corner-derived `d=3` candidate chart

The generated-corner identity

\[
A_1=A'_0+\gamma(1/5,1)
\]

and `A1=(7/5,2)` force `gamma=2`.  The direct monomial substitution

\[
X=\xi y^2,\qquad Y=\xi^{-2}y^{-7}
\]

sends the certified bracket `X^4` to `-3*xi^2*y^2`, the terminal P/Q
vertices to leading `xi` degrees `3r,3(2r-1)`, and the common root to

\[
\xi^3H(\xi^{-1}y^{-5}).
\]

The exact source-band and terminal-halfspace inequalities give the
leading-minus-one coefficient support `{-2,-5}`.  Hence the normalized
translation is

\[
\xi=x-G,\qquad G=g_{-2}y^{-2}+g_{-5}y^{-5}.
\]

This derives `d=3`, `h=2`, and the finite coefficient **envelopes** from the
corner chain.  It also shows that the formal `d=2` comparison is not an F2
branch.  The envelopes are not independent coefficient spaces: the exact
`r=3` nonnegative projections have `9` P relations and `19` Q relations.
Their top-diagonal equations prove that cutting off the negative `xi` tail
cannot produce a polynomial pair with the same bracket.  The complete
derivation, all `r=3` masks, the unit-ideal certificate, and the remaining
tail-correction lemma are in
[`F2_MODIFIED_CHART_BRIDGE.md`](F2_MODIFIED_CHART_BRIDGE.md).
<!-- status-consumer: PF2MCB1 6ff13314e0090f52 -->

## 3. Conditional modified-Laurent family theorem

Assume now the chart (2).  Put

\[
n=rd,\qquad t=x^{-1},\qquad
C=x^dA(t),\qquad
A(t)=1+\sum_{q\ge2}a_qt^q. \tag{13}
\]

The series `A` is infinite.  Only the finite window

\[
L=dm+n-1=d(3r-1)-1 \tag{14}
\]

is needed for the rows considered here.  In particular,
`d-L=1-d(3r-2)` is the last inspected `x`-exponent of `C`, not a lower
support bound.

Polynomiality of `P=C^r` supplies

\[
[t^{n+k}]A^r=0,
\qquad1\le k\le dm-1. \tag{15}
\]

The `k`-th equation has the new term `r*a_(n+k)`.  Thus (15) triangularly
determines

\[
a_{n+1},\ldots,a_L. \tag{16}
\]

There are

\[
dm-1=d(2r-1)-1 \tag{17}
\]

power rows and `n-1=rd-1` remaining `Q` rows.  The total number of inspected
coefficient functions is

\[
L-1=d(3r-1)-2. \tag{18}
\]

The four regression cases are therefore:

| `r,d` | inspected `a_q` | power rows | `Q` rows |
| --- | ---: | ---: | ---: |
| `2,2` | 8 | 5 | 3 |
| `2,3` | 13 | 8 | 5 |
| `3,2` | 14 | 9 | 5 |
| `3,3` | 22 | 14 | 8 |

### 3.1 Exact `r=2` regression

The comparison source is Section 5 of
[`external/arxiv-1406.0886/Polynomial_system.tex`](external/arxiv-1406.0886/Polynomial_system.tex).

For `d=2`, write

\[
C=x^2+u+a x^{-1}+b x^{-2}+\cdots.
\]

After triangular elimination the three published left-hand sides are

\[
R_1=3ab, \tag{19}
\]

\[
R_2=\lambda-\frac32ua^2+\frac32b^2, \tag{20}
\]

\[
R_3=-\frac12a(a^2+6ub). \tag{21}
\]

They satisfy `R_i=-F_(-i)` in the historical `gamma=3` chart.

For `d=3`, write

\[
C=x^3+s x+u+a x^{-1}+b x^{-2}+c x^{-3}+\cdots.
\]

The five rows are

\[
R_1=\frac32(sa^2+2ac+b^2), \tag{22}
\]

\[
R_2=-\frac32(ua^2-2bc), \tag{23}
\]

\[
R_3=\lambda-3uab-\frac32sb^2-\frac12a^3+\frac32c^2, \tag{24}
\]

\[
R_4=-\frac32(2uac+ub^2+2sbc+a^2b), \tag{25}
\]

\[
R_5=\frac12\left(
-2\lambda s+3u^2a^2+6usab-6ubc+3s^2b^2-3sc^2-3a^2c-3ab^2
\right). \tag{26}
\]

The published assignments are `R_1=R_2=0`, `R_3=-y^3`, and
`R_i=-F_(-i)` for `i=4,5`.  The compiler reproduces these polynomials, but
does not supply the omitted preliminary proof that the two historical
`gamma` charts are exhaustive.

## 4. Compact triangular elimination in polynomial coordinates

Define

\[
B(t)=\operatorname{trunc}_{\le n}(A(t)^r)
=1+\sum_{i=2}^{n}b_it^i,
\qquad
S=B^{-1/r},
\qquad
M=B^2S=B^{2-1/r}. \tag{27}
\]

The change from `a_2,...,a_n` to `b_2,...,b_n` is triangular with diagonal
`r`.  Equations (15) give

\[
A=B^{1/r}\pmod {t^{L+1}}. \tag{28}
\]

Since `m=2r-1`, the remaining core becomes

\[
C^m+\lambda C^{-1}
=x^{dm}\left(M+\lambda t^{2n}S\right). \tag{29}
\]

Consequently row `k`, for `1<=k<n`, is

\[
R_k=M_{dm+k}+\lambda S_{k-d}, \tag{30}
\]

where `S_j=0` for `j<0`.  Lambda first occurs in row `k=d`, with coefficient
one.  Exact all-`r` recurrences are

\[
rjS_j=
\sum_{i=1}^{\min(j,n)}((r-1)i-rj)b_iS_{j-i}, \tag{31}
\]

\[
rjM_j=
\sum_{i=1}^{\min(j,n)}((3r-1)i-rj)b_iM_{j-i}, \tag{32}
\]

with `b_0=1`, `b_1=0`, and `b_i=0` for `i>n`.

### 4.1 The bracket derives the visible `F` rows

Suppose the chosen modified chart has bracket `x`-degree `h`, and let
`f(y)x^ell` be the first term of `F` that contributes nontrivially to the
bracket.  Comparing its bracket with the leading degree `n` part of `P`
gives

\[
n+\ell-1=h,
\qquad
\ell=h+1-n. \tag{33}
\]

Thus the first visible `F` row is

\[
k_0=n-h-1. \tag{34}
\]

Within rows `1,...,n-1`, the zero rows and freely adjustable visible rows are

\[
1,\ldots,n-h-2, \tag{35}
\]

and

\[
n-h-1,\ldots,n-1, \tag{36}
\]

respectively.  There are `h+1` rows in (36).  They are not the whole support
of `F`; its formal tail continues below this finite window.

Lambda is an independent zero-row pivot exactly when

\[
d<k_0
\quad\Longleftrightarrow\quad
d(r-1)>h+1. \tag{37}
\]

Only in this case does the zero-row system itself imply

\[
\lambda=-M_{2n}. \tag{38}
\]

The historical `r=2` chart has `h=2`, hence three visible rows.  The
corner-chain chart certified for `r=3` has bracket `X^4`, hence `h=4` and
five visible rows.  No proved coordinate bridge currently carries the
certified `X^4` chart to an `h=2` modified chart.

### 4.2 A uniform congruence-support bracket gate

There is one support statement which needs no Laurent-`y` mask.  Suppose

\[
B(t)\in K[y,y^{-1}][t^d].
\]

Then `M=B^(2-1/r)` and `S=B^(-1/r)` also have only `t`-degrees divisible
by `d`.  Assume `0<=h<=dr-2`, so the first bracket-visible row lies in the
compiled window.  It is

\[
k_0=dr-h-1.
\]

The two series indices in that row are

\[
d(2r-1)+k_0=d(3r-1)-(h+1),
\]

and

\[
k_0-d=d(r-1)-(h+1).
\]

Both are congruent to `-(h+1)` modulo `d`.  Therefore

\[
\boxed{
d\nmid(h+1)\quad\Longrightarrow\quad F_{-k_0}=0
\text{ on }B(t)\in K[y,y^{-1}][t^d].
} \tag{38-cong}
\]

But the leading bracket term is `dr*F_(-k0)'*x^h`.  Hence a nonzero bracket
of degree `h` excludes the whole congruence section whenever
`d` does not divide `h+1`.  This theorem is uniform in `r`.

The checker substitutes the congruence section into all six compiled cases.
For the corner-derived chart, `h=4` and `h+1=5`.  Whenever `dr>=6`, the only
integer `d>=2` whose congruence section can escape this test is therefore
`d=5`.  In particular, neither `d=2` nor `d=3` divides `5`, so both proposed
`r=3` congruence sections are incompatible with the certified nonzero `X^4`
bracket.  This is a genuine test of the certified bracket weight inside the
conditional common-power charts; it does not eliminate the full `d=2` or
`d=3` residue, and it does not prove that either chart is forced by the
corner chain.

For the historical-tail choice `h=2`, the same gate excludes the `d=2`
congruence section, while `d=3` is the unique one of these two cases not
removed because `3` divides `h+1`.  Section 5.4 analyzes exactly that
surviving cubic section.

### 4.3 A universal endpoint-binomial section

The compact residue is never the unit ideal.  Indeed, for every `r>=2` and
`d>=2`, let `n=dr`, choose `a!=0`, and set

\[
B(t)=1+a t^n,\qquad
\lambda=-\frac{(2r-1)(r-1)}{2r^2}a^2. \tag{38a}
\]

Both `M=B^(2-1/r)` and `S=B^(-1/r)` are supported in degrees divisible by
`n`.  In the row window `1<=k<n`,

\[
dm+k=2n-d+k
\]

is divisible by `n` only for `k=d`; likewise `S_(k-d)` can be nonzero in
this window only for `k=d`.  At that pivot,

\[
M_{2n}=\binom{2-1/r}{2}a^2
=\frac{(2r-1)(r-1)}{2r^2}a^2, \tag{38b}
\]

so (39) gives

\[
\boxed{R_k=0\quad(1\le k<n)}. \tag{38c}
\]

All visible `F` coefficients may be set to zero at this point.  Hence every
compiled Artinian/Fitting residue ideal is proper, uniformly in `r,d,h`.
This is stronger than checking the four generated small cases, although the
checker also substitutes (39) into all six compiled `h=2,4` records.

The lacunarity locates the first lower coefficient missed by the compact
window.  Rows `n,...,n+d-1` still vanish, while row `k=n+d`, at `x`-exponent
`-(n+d)`, equals

\[
\binom{2-1/r}{3}a^3
-\lambda\frac{a}{r}
=\boxed{\frac{(2r-1)(r-1)}{3r^3}a^3}. \tag{38d}
\]

Thus any realization of this section must acquire a matching lower `F`
coefficient exactly at that band.

There is an important claim boundary.  Equation (38c) also sets the first
bracket-visible `F` coefficient to zero unless its row coincides with the
lambda pivot.  A nonzero bracket can therefore exclude this special section.
Its role is to prove that the residue ideal itself can never eliminate the
family and to identify the first lower band which a support mask must test;
it is not by itself a plane-map candidate.

## 5. The `r=3` residual systems

Write

\[
\mu_j=[t^j]B^{5/3},\qquad s_j=[t^j]B^{-1/3}. \tag{39}
\]

### 5.1 Conditional historical tail `h=2`

For `d=2`, `n=6`, the zero-row residue is

\[
\boxed{\mu_{11}=0},\qquad \lambda=-\mu_{12}. \tag{40}
\]

The visible coefficients are then solved by

\[
F_{-3}=-\mu_{13},
\]

\[
F_{-4}=-(\mu_{14}-\mu_{12}s_2),
\]

\[
F_{-5}=-(\mu_{15}-\mu_{12}s_3). \tag{41}
\]

For `d=3`, `n=9`, lambda is `-mu_18` and the exact residual ideal is

\[
\boxed{
(\mu_{16},\mu_{17},\mu_{19},\mu_{20}-\mu_{18}s_2)
}. \tag{42}
\]

The visible rows solve

\[
F_{-6}=-(\mu_{21}-\mu_{18}s_3),
\]

\[
F_{-7}=-(\mu_{22}-\mu_{18}s_4),
\]

\[
F_{-8}=-(\mu_{23}-\mu_{18}s_5). \tag{43}
\]

Without Laurent-`y` restrictions, (40) is one equation in five `b_i`, while
(42) is four equations in eight `b_i`.  Neither residual algebra is
Artinian.  These are tiny exact formal candidates, not plane-counterexample
certificates.

### 5.2 The actually certified bracket degree `h=4`

For `r=3,d=2`, all five `Q` rows are visible `F` rows, so the generic formal
residue is empty.  For `r=3,d=3`, only rows `1,2` remain as zero-row
conditions.  The certified bracket degree therefore weakens rather than
strengthens the formal obstruction.

### 5.3 Exact geometry of the `d=2` residue

The determinant (40) is not close to an exclusion.  In polynomial
coordinates its generator is a nonzero scalar times

\[
\begin{aligned}
E={}&7b_2^4b_3-12b_2^3b_5-36b_2^2b_3b_4-12b_2b_3^3
 +54b_2b_4b_5+27b_3^2b_5+27b_3b_4^2\\
&+54(b_2b_3-3b_5)b_6.
\end{aligned} \tag{43a}
\]

Consequently, on `b_2*b_3-3*b_5!=0`, it solves rationally as

\[
b_6=-\frac{
7b_2^4b_3-12b_2^3b_5-36b_2^2b_3b_4-12b_2b_3^3
+54b_2b_4b_5+27b_3^2b_5+27b_3b_4^2}
{54(b_2b_3-3b_5)}. \tag{43b}
\]

Thus this open part of the residue is a rational fourfold.  An exact smooth
coefficient-torus point is

\[
(b_2,b_3,b_4,b_5,b_6)
=\left(1,1,1,1,\frac{55}{108}\right). \tag{43c}
\]

At (43c), the residue Jacobian has rank one and

\[
\lambda=\frac{7525}{34992},\qquad
(F_{-3},F_{-4},F_{-5})
=\left(\frac{25}{243},-\frac{85}{2916},-\frac{701}{8748}\right). \tag{43d}
\]

Every displayed value is nonzero.  In particular, even the necessary
nonvanishing of the first bracket-visible coefficient `F_(-3)` does not
remove this component.  Any `d=2` exclusion must use the actual
Laurent-`y` masks or further bracket equations, not the determinant alone.

### 5.4 Cubic-invariant reduction of the `d=3` residue

There is an exact three-parameter congruence stratum

\[
B(t)=1+u t^3+v t^6+w t^9. \tag{43e}
\]

All four generators in (42) vanish identically on (43e).  Put

\[
A=u^2-3v,\qquad D=2u^3-9uv+27w. \tag{43f}
\]

The lambda pivot and visible tail reduce to

\[
\boxed{
\lambda=-\frac5{6561}(3A^3+D^2),\qquad
F_{-6}=\frac5{6561}A^2D,\qquad
F_{-7}=F_{-8}=0.
} \tag{43g}
\]

These are the two translation invariants of the cubic

\[
\bar P=X^3+uX^2+vX+w.
\]

Let `bar Q` be the polynomial part at infinity of `bar P^(5/3)`.  Treating
`u,v,w` as functions of `y`, exact differentiation gives

\[
[\bar P,\bar Q]_{X,y}
=(3X+u)\lambda'(y)+3F_{-6}'(y). \tag{43h}
\]

Since `X=x^3`, the corresponding full bracket is three times `x^2` times
(43h).  The ansatz requires `lambda` to be a scalar, so on this stratum the
bracket collapses to

\[
[P,Q]_{x,y}=9x^2F_{-6}'(y). \tag{43i}
\]

This converts a four-generator Fitting condition into one exact invariant
support gate.  In fact, no single-monomial support assumption on `F_(-6)`
is needed.  Suppose only that the leading coefficient of the required
nonzero bracket is a Laurent monomial

\[
[P,Q]_{x,y}=\mu y^s x^2+\text{lower powers of }x,
\qquad \mu\ne0. \tag{43i'}
\]

Comparing the `x^2` coefficient in (43i) gives

\[
9F_{-6}'=\mu y^s. \tag{43i''}
\]

If `s=-1`, this has no Laurent-polynomial antiderivative.  Otherwise, with
`N=s+1!=0`,

\[
A^2D=\alpha y^N+\beta,
\qquad \alpha\ne0. \tag{43i'''}
\]

The integration constant `beta` cannot be nonzero.  Indeed, for
`f=alpha*y^N+beta`,

\[
Nf-yf'=N\beta. \tag{43i''''}
\]

Thus `f` is squarefree in the Laurent UFD, while `A^2` divides `f`; hence
`A=a*y^p` is a Laurent unit.  If `beta!=0`, then

\[
D=b y^{N-2p}+c y^{-2p},\qquad bc\ne0.
\]

The three terms of `D^2` have exponents

\[
2N-4p,\qquad N-4p,\qquad -4p. \tag{43i'''''}
\]

They are pairwise distinct because their differences are `N,N,2N`.
The single term `3A^3` can cancel at most one of them.  Constancy of
`3A^3+D^2` would therefore force two of (43i''''') to be zero, but every
pair of those equations forces `N=0`, a contradiction.  Hence `beta=0`.

Now `A^2D` is a Laurent unit, so

\[
A=a y^p,\qquad D=b y^q.
\]

Constancy of `3A^3+D^2` forces either `p=q=0`, or

\[
3p=2q\ne0,\qquad 3a^3+b^2=0.
\]

In the nonconstant case there is an integer `ell` with

\[
p=2\ell,\qquad q=3\ell,\qquad N=2p+q=7\ell,\qquad \lambda=0. \tag{43j}
\]

Consequently the congruence stratum is excluded unless

\[
s=7\ell-1\quad(\ell\ne0),\qquad \lambda=0. \tag{43j'}
\]

This is sharp at the finite compiled-system level: `u=0`,
`v=y^(2*ell)`, and `w=y^(3*ell)/3` give

\[
A=-3y^{2\ell},\qquad D=9y^{3\ell},\qquad
\lambda=0,\qquad F_{-6}=\frac5{81}y^{7\ell}. \tag{43j''}
\]

Together with `F_(-7)=F_(-8)=0`, these values solve all eight fiber rows in
the `r=3,d=3,h=2` window.  The checker verifies the first member explicitly:

\[
B(t)=1+y^2t^6+\frac13y^3t^9,\qquad
\operatorname{coeff}_{x^2}[P,Q]_{x,y}=\frac{35}{9}y^6.
\tag{43j'''}
\]

The chart bridge now settles the weight comparison: F2 has `gamma=2`, so
its derived modified bracket is `mu*y^2*(x-G)^2`.  Therefore
(43j'')--(43j''') is **not** an F2 candidate.  It remains a sharp witness
for the unrestricted cubic invariant gate, but its `y^6` leading weight
excludes it from the forced F2 chart.

The correct weight-two resonance is instead uniform in `r`:

\[
B(t)=1+a y t^{3r-2},
\]

\[
P_T=x^{3r}+a x^2y,
\]

\[
Q_T=x^{3(2r-1)}+\frac{2r-1}{r}a x^{3r-1}y
+\frac{(2r-1)(r-1)}{2r^2}a^2xy^2,
\]

and

\[
[P_T,Q_T]=
\frac{3(r-1)(2r-1)}{2r^2}a^3x^2y^2. \tag{43j''''}
\]

This is exactly the formal terminal edge with its upper common-root data
deleted.  For `r=3` it is already outside the exact source image.  On the P
top diagonal, the primitive source relation (16) of the companion chart note
evaluates to `3470256`, rather than zero.  Thus (43j'''') is not a plane
counterexample and is not the base point of a source-compatible projected
branch.

For example, the historical `d=3` bracket weight `s=2` would give `N=3`
and therefore excludes this whole stratum; a constant leading bracket
coefficient gives `N=1` and also excludes it.  The bracket weight and
coefficient envelopes are now corner-derived; the remaining conditional
step is the Laurent-to-polynomial cut.  The argument no longer depends on an
unproved support mask for `F_(-6)`.  Equations (43j'')--(43j''') remain a
sharp non-F2 test point for the invariant calculation, while (43j'''') is
only an ambient support-box point excluded by the exact source relations.

The full `d=3` residue is nevertheless much larger than (43e).  At
`u=v=w=1`, equivalently `b_3=b_6=b_9=1` with the other `b_i` zero, the four
residue equations have Jacobian rank four; the pivot minor in
`(b_2,b_4,b_5,b_7)` is

\[
\frac{6988960000}{847288609443}\ne0. \tag{43k}
\]

At this point

\[
\lambda=-\frac{1880}{6561},\qquad F_{-6}=\frac{400}{6561}. \tag{43l}
\]

Fix `b_3=b_6=b_9=1` and `b_8=epsilon`.  The formal implicit-function theorem
and exact coefficient comparison give the unique two-jet

\[
\begin{aligned}
b_2&=3\epsilon+O(\epsilon^3),&
b_4&=3\epsilon^2+O(\epsilon^3),\\
b_5&=2\epsilon+O(\epsilon^3),&
b_7&=\epsilon^2+O(\epsilon^3).
\end{aligned} \tag{43m}
\]

Over `Q((epsilon))`, every `b_2,...,b_9` is therefore nonzero, while lambda
and `F_(-6)` remain units.  The `d=3` Fitting scheme has a smooth
four-dimensional component meeting the full coefficient torus.  As in
`d=2`, only the unproved Laurent-`y` masks and the complete bracket/support
ledger can turn this surviving formal geometry into either a contradiction
or a plane candidate.

## 6. Direct Artinian and Fitting presentation

Let

\[
R=\mathbf Q[1/r,b_2,\ldots,b_n],
\qquad
J=R[t]/(t^{L+1}),
\qquad
Y=B^{1/r}, \tag{44}
\]

and put

\[
U=\langle1,t,\ldots,t^{dm}\rangle,
\qquad
T_h=\langle t^{L-h},\ldots,t^L\rangle. \tag{45}
\]

Multiplication by the unit `Y` preserves `T_h`.  Multiplying (29) by `Y`
shows that the residual condition is exactly

\[
\boxed{B^2\in YU+T_h+R\,t^{2n}\subset J}. \tag{46}
\]

Let

\[
\epsilon=
\begin{cases}
1,&d<n-h-1,\\
0,&d\ge n-h-1.
\end{cases} \tag{47}
\]

The base coefficient matrix has `dm+1` columns from `YU`, `h+1` terminal
columns, and the lambda column precisely when it is independent.  Its rank is

\[
dm+h+2+\epsilon. \tag{48}
\]

Since `rank_R(J)=dm+n`, the quotient rank is

\[
q=n-h-2-\epsilon. \tag{49}
\]

Append the coefficient vector of `B^2`.  When `q>0`, condition (46) is the
vanishing of the maximal minors of size

\[
dm+h+3+\epsilon. \tag{50}
\]

Equivalently it is

\[
\operatorname{Fitt}_{q-1}
\operatorname{coker}[YU,T_h,t^{2n},B^2]. \tag{51}
\]

The base columns have distinct unit leading terms, so row reduction identifies
the ideal (51) with the zero-row equations in Section 5.  In particular:

- `r=3,d=2,h=2`: the appended matrix is `16 x 16`; its determinant is a
  unit times `mu_11`.
- `r=3,d=3,h=2`: the appended matrix is `24 x 21`; its `21 x 21` minors
  generate the four-element ideal (42), namely `Fitt_3`.
- `r=3,d=2,h=4`: `q=0`, so there is no generic residue condition.
- `r=3,d=3,h=4`: `q=2`, giving `Fitt_1` and two zero-row generators.

There is also a useful but different Toeplitz determinant.  Differentiating
the fixed-`(lambda,F)` coefficient map gives multiplier

\[
G(t)=\frac{2r-1}{r}BS
-\frac{\lambda}{r}t^{2n}\frac{S}{B}, \tag{52}
\]

and matrix

\[
T_{kq}=[t^{dm+k-q}]G(t),
\quad1\le k<n,\quad2\le q\le n. \tag{53}
\]

Its determinant generates `Fitt_0` of the corresponding square coefficient-map
Kahler module.  It proves generic nondegeneracy of that map; it is not the
residual obstruction ideal (51).

## 7. Family-level terminal residue theorem

The last edge of (3), unlike the full modified chart, has a complete uniform
normalization.  Put

\[
s=X^{7r-4}y^{2r-1},
\]

\[
P_I=X^4y(1+s),
\]

\[
Q_I=-X(1+A_rs+B_rs^2), \tag{54}
\]

where

\[
A_r=\frac{2r}{r-1},
\qquad
B_r=\frac{2r^2}{(r-1)(2r-1)}. \tag{55}
\]

For undetermined coefficients `A,B`, direct differentiation gives

\[
X^{-4}[P_I,Q_I]
=1+(2r-(r-1)A)s+(rA-(2r-1)B)s^2. \tag{56}
\]

Thus (55) is the unique normalized solution and

\[
\boxed{[P_I,Q_I]=X^4}. \tag{57}
\]

In `t=Xy,z=y^(-1)` coordinates, the primitive edge direction and normal are

\[
(7r-4,5r-3),
\qquad
\nu=(5r-3,4-7r). \tag{58}
\]

The Bezout identity

\[
5(7r-4)-7(5r-3)=1
\]

proves primitivity.  The pole orders of `(P_I,Q_I)` along `nu` are
`(r,2r-1)`.  At target infinity put `a=(-Q_I)^(-1)` and `b=P_I/(-Q_I)`.
Their pullback orders are

\[
(2r-1,r-1). \tag{59}
\]

Since

\[
(2r-1)-2(r-1)=1,
\]

`a/b^2` is a transverse uniformizer.  The transverse index is therefore
one, while the residue coordinate is

\[
h_r(s)=\frac{b^{2r-1}}{a^{r-1}}
=\frac{P_I^{,2r-1}}{(-Q_I)^r}
=\frac{s(1+s)^{2r-1}}{(1+A_rs+B_rs^2)^r}. \tag{60}
\]

The numerator and denominator are coprime of degree `2r`, and

\[
\frac{h_r'}{h_r}
=\frac1{s(1+s)(1+A_rs+B_rs^2)}. \tag{61}
\]

It follows that the complete passport is

\[
\boxed{
(2r-1,1)\mid(r,r)\mid(3,1^{,2r-3})
}. \tag{62}
\]

The three branch values are `0`, infinity, and `B_r^(-r)`.

### 7.1 Geometric and arithmetic monodromy

On `2r` letters, set

\[
\alpha=(1\ 2\ \cdots\ 2r-1),
\qquad
\gamma=(1\ r\ 2r), \tag{63}
\]

\[
\beta=(1\ 2r\ r-1\ r-2\ \cdots\ 2)
(r\ 2r-1\ 2r-2\ \cdots\ r+1). \tag{64}
\]

With composition from right to left,

\[
\alpha\beta\gamma=1. \tag{65}
\]

These permutations have the three types in (62), and all are even.  The
action is transitive.  The `(2r-1)`-cycle fixes one point and is transitive on
the others, so the group is two-transitive and hence primitive.  For `r>=3`,
the literal `3`-cycle and Jordan's theorem imply that the group contains
`A_(2r)`.  Since all generators are even,

\[
\boxed{G_{\rm geom}=A_{2r}}. \tag{66}
\]

The `r=2` case is checked directly and gives `A_4`.

For

\[
f_u(s)=s(1+s)^{2r-1}
-u(1+A_rs+B_rs^2)^r,
\]

exact resultant calculation gives

\[
\operatorname{disc}_s(f_u)
=\frac{(-1)^{r+1}}{(2r-1)^{2(r-1)^2-1}}
u^{2r-2}(1-B_r^ru)^2. \tag{67}
\]

For completeness, put `g=1+A_r*s+B_r*s^2`,
`a=1-B_r^r*u`, and `G=g(-1)=-1/(2r-1)`.  At every root `alpha` of `f_u`,
equation (61) gives

\[
f_u'(\alpha)=\frac{(1+\alpha)^{2r-2}}{g(\alpha)}. \tag{68}
\]

The product formulas obtained from `f_u(-1)=-uG^r` and the two roots of
`g` are

\[
\prod_\alpha(1+\alpha)=-\frac{uG^r}{a},
\qquad
\prod_\alpha g(\alpha)=\frac{G^{2r-1}}{a^2}. \tag{69}
\]

Hence

\[
\operatorname{Res}(f_u,f_u')
=a^3u^{2r-2}G^{2r^2-4r+1}. \tag{70}
\]

Using
`disc(f_u)=(-1)^r*a^(-1)*Res(f_u,f_u')` gives (67) directly.

Its squareclass over `Q(u)` is

\[
\boxed{(-1)^{r+1}(2r-1)}. \tag{71}
\]

Thus the arithmetic group is `S_(2r)` unless (71) is a rational square, in
which case it is `A_(2r)`.  Equivalently, the exceptional arithmetic-`A`
case occurs when `r` is odd and `2r-1` is a square.

The regular geometric closure has genus

\[
g_{\rm reg}=1+\frac{|A_{2r}|}{2}
\left(1-\frac1{2r-1}-\frac1r-\frac13\right). \tag{72}
\]

For `r=3`, equations (54)--(72) specialize to the certified degree-six
`A_6` residue cover and its genus-`25` regular closure.

## 8. The extracted `A_6` permutation-gluing audit

The companion
[`F2_A6_SIMPLE_SPECTATOR_GLUING.md`](F2_A6_SIMPLE_SPECTATOR_GLUING.md)
tests the suggested spectator model.  Its bridge assumption is explicit:
each of the two simple squarefree-`R` Kummer orbits supplies one transposition
on the same target `P^1`, the terminal `A_6` triple fixes every outside sheet,
and there are no other branch cycles.

Under that assumption, the meridian product makes the two transpositions
equal.  Transitivity and Riemann--Hurwitz force degree `7` and genus `0`.
There are exactly `30` normalized tuples, `6` simultaneous-conjugacy classes,
and `30,240` fully sheet-labelled tuples.  Every tuple generates `S_7`.
The certified endpoint/interior labels prune this to three classes of
signature `(5,3,1)` under the strongest naive requirement that the connector
anchor avoid both toric source endpoints.  Thus even the marked
simple-spectator model survives; it does not give a contradiction.

More generally, repeated star connectors construct connected genus-zero
witnesses in every degree `6+k`.  This proves that product one,
transitivity, Riemann--Hurwitz, and remaining degree alone cannot exclude the
row once additional simple branch cycles are allowed.  The actual Kummer
orbits have not yet been proved to be transpositions or placed on a common
target component, so the enumeration remains conditional.

## 9. Exact conclusion

The requested family-level exclusion is not presently justified.  What is
proved is sharper than a failed search:

1. the corner chain gives the exact degree family, two forced edges, source
   envelope, and common top band; its root has the explicit polynomial source
   lift `x*(x*y^5-1)^2*R(x*y^5)`;
2. the assumed modified-Laurent family has a closed triangular compiler;
3. its `r=3,h=2` residue is the determinant (40) or the four-generator
   Fitting ideal (42);
4. every compact residue ideal is proper, the `d=2` residue contains a
   rational fourfold, and the `d=3` residue contains a smooth fourfold meeting
   the full coefficient torus;
5. the uniform gate (38-cong) excludes every `d=2,3` congruence section
   under the certified `h=4` bracket weight, and also excludes `d=2` for
   `h=2`;
6. the surviving `d=3,h=2` congruence stratum reduces to the two cubic
   invariants (43g); its integration constant is impossible, and a monomial
   leading bracket leaves exactly the `lambda=0`, exponent-`7*ell` ray
   (43j') with the tiny finite candidate (43j'');
7. the full certified-`h=4` charts still have zero or two generic residue
   equations, so the congruence cut does not eliminate their unrestricted
   residues;
8. the terminal edge gives an unconditional degree-`2r`, geometric-`A_(2r)`
   residue theorem; and
9. the coarse and marked `A_6` spectator gluing constraints admit exact
   solutions; and
10. in the literal polynomial projection, the `r=3` P top-gap scheme has
    Artinian length `27`, the first Q gap is a unit in that algebra, and the
    combined top-band ideal is `(1)`; `r=4` and the optional exact Singular
    `r=5` sample are already killed by P alone; and
11. the corrected all-`r` top tangent retains every P-band direction.  Its
    nonnegative centralizer modes are `C^(m-1),...,1`; the negative modes
    `C^(-1),...,C^(-(r-1))` are formal resonances but not independent source
    kernels.  For `r=3`, `lambda*C^(-1)` first occurs at layer `10`; and
12. the `r=3` nonlinear rows recover a source root through descent `7` and
    isolate the first exact local defect at descent `8` by the quadratic
    `27*y^2-9*y+1`, over `QQ(sqrt(-3))`.

An exclusion now requires a tail-correction theorem deriving the actual
polynomial modified pair from the full positive and negative Laurent bands,
or stronger global incidence data that rules out the surviving boundary
gluings.  Equations (40) and (42) remain the smallest exact **conditional
formal-series residue equations**, but their smooth loci are neither source
projections nor verified planar candidates.  There is no residual `r=3`
candidate in the literal projection: its exact top-band ideal is already the
unit ideal.  The coefficient-side target is therefore the negative-tail
cross term omitted by that projection, not another search inside the old
14- or 22-function support boxes.

## 10. Reproduction

Run:

```text
.venv/bin/python plane-jc/cas/generate_f2_modified_system.py --include-equations
.venv/bin/python plane-jc/cas/generate_f2_modified_system.py --include-equations --output artifacts/generated-results/jc2_f2_modified_laurent_family.json
.venv/bin/python plane-jc/cas/verify_f2_modified_chart_bridge.py
# Optional; requires Singular:
.venv/bin/python plane-jc/cas/verify_f2_modified_chart_bridge.py --extended-r5
.venv/bin/python plane-jc/cas/verify_f2_terminal_residue_cover.py
.venv/bin/python plane-jc/cas/verify_f2_a6_simple_spectator_gluing.py
```

The expected final markers include:

```text
F2_50_75_PUBLISHED_REGRESSION_PASS
F2_75_125_D2_SYSTEM_PASS
F2_75_125_D3_SYSTEM_PASS
F2_MODIFIED_ARTINIAN_RESIDUE_FITTING_PASS
F2_R3_D2_RESIDUE_RATIONAL_FOURFOLD_PASS
F2_R3_D3_LAURENT_BRACKET_GATE_PASS
F2_R3_D3_RESIDUE_SMOOTH_TORUS_BRANCH_PASS
F2_R3_D3_CUBIC_INVARIANT_REDUCTION_PASS
F2_MODIFIED_ENDPOINT_BINOMIAL_BRANCH_PASS
F2_MODIFIED_CONGRUENCE_SUPPORT_BRACKET_PASS
F2_MODIFIED_EXACT_SOURCE_PROJECTION_RANK_PASS
F2_MODIFIED_TERMINAL_OUTSIDE_SOURCE_IMAGE_PASS
F2_R3_PROJECTED_TOP_GAP_ARTINIAN_LENGTH=27
F2_R3_PROJECTED_TOP_GAP_RESULTANT_NONZERO_PASS
F2_R3_PROJECTED_TOP_GAP_UNIT_IDEAL_PASS
F2_TERMINAL_FAMILY_GEOMETRIC_MONODROMY=A_(2*r)
F2_TERMINAL_RESIDUE_COVER_PASS
F2_A6_MARKED_INTERIOR_FILTER_CONJUGACY_CLASSES=3
F2_A6_SIMPLE_SPECTATOR_GLUING_CONDITIONAL_PASS
```
