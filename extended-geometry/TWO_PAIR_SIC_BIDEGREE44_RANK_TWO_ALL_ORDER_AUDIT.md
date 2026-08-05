# Prerequisite audit for the bidegree-\((4,4)\) rank-two problem

## 1. Outcome

The proposed all-order calculation cannot yet be specialized to a known
rank-two survivor.  The current exact input is weaker:

1. the first thirteen moments have a semistable common zero on the
   rank-at-most-two determinantal variety;
2. the Hilbert-series proof is existential and records neither coordinates
   nor a residue field for its semistable point.

Consequently there is no exact coefficient point at which to derive and
solve a scalar recurrence.  Recurrence derivation and evaluation of
\(\mu_{14}\) are therefore parked.  Neither an all-order rank-two witness
nor an exact tail obstruction is proved here.

The split-symbol theorem has removed the former rank-one prerequisite
from the SIC search: the complete rank-one Segre cone is safe in every
degree.  The remaining squarefree rank-one Rabinowitsch membership is
relevant only to the stronger finite-prefix claim that the existential
thirteen-moment point itself must have exact rank two.  It is not a gate
for searching directly on exact-rank-two factor charts.

This corrects a possible misreading of equation (6.1) in the
[rank frontier](TWO_PAIR_SIC_BIDEGREE44_RANK_FRONTIER.md).  That displayed
rank-two matrix is a transversality point used to prove Jacobian rank
thirteen.  It is not a moment zero:
\[
\begin{aligned}
 \mu _1&=7414,\\
 \mu _2&=3675739680,\\
 \mu _3&=12167497410877440,\\
 \mu _4&=148010006143680629760000.
\end{aligned} \tag{1.1}
\]

## 2. Parked holonomic infrastructure

The determinantal parameterization can nevertheless be put into the form
needed by creative telescoping.  Write
\[
 C=UW,\qquad
 U=(u_{iq})\in\operatorname {Mat}_{5\times2},\qquad
 W=(w_{qj})\in\operatorname {Mat}_{2\times5}. \tag{2.1}
\]
With
\[
 \Phi_C(x,y)=\sum_{i,j=0}^4c_{ij}x^iy^j,
\]
define the Laurent polynomial
\[
\begin{aligned}
 A_q(u)&=\sum_{i=0}^4u_{iq}u^{4-i},\\
 B_q(u,t)&=\sum_{j=0}^4
   w_{qj}u^{j-4}t^j(1-t)^{4-j},\\
 P_{U,W}(u,t)&=A_1(u)B_1(u,t)+A_2(u)B_2(u,t).
\end{aligned} \tag{2.2}
\]
Direct multiplication gives
\[
 P_{U,W}(u,t)
 =\Phi_C\left(1,u,t,\frac{1-t}{u}\right). \tag{2.3}
\]
The formal beta identity from the full-rank witness therefore gives
\[
 \boxed{\frac{\mu_m(C)}{(4m+1)!}
 =\operatorname {CT}_u\int_0^1P_{U,W}(u,t)^m\,dt.} \tag{2.4}
\]
Equivalently, the factorial-normalized ordinary generating function is
\[
 \boxed{
 {\cal G}_{U,W}(s)
 =\sum_{m\geq0}\frac{\mu_m(C)}{(4m+1)!}s^m
 =\operatorname {CT}_u\int_0^1
   \frac{dt}{1-sP_{U,W}(u,t)}.} \tag{2.5}
\]

Formula (2.5) is valid holonomic infrastructure, but it is not an active
recurrence calculation.  Closure of holonomic functions under coefficient
extraction and definite integration implies that every exact algebraic
specialization of \(U,W\) has a P-finite moment sequence.  No such
specialization is currently available from the semistable fiber.  Once an
explicit exact-rank-two point or component is known, an all-order result
will still require an explicit scalar operator, telescoping certificates
including the endpoint terms, its singular-step audit, and enough exact
initial values.

### 2.1 Quotient the internal gauge first

The factor variables in (2.1) should not be submitted to a moment solver
with their four-dimensional internal gauge still present.  Let \(I\) be a
pair of rows with \(\det U_I\ne0\).  Applying (2.1) with
\(G=U_I^{-1}\) gives the unique representative

\[
 U'=
 \begin{pmatrix}
  1&0\\0&1\\a_{20}&a_{21}\\a_{30}&a_{31}\\a_{40}&a_{41}
 \end{pmatrix},\qquad
 B=(U_IV^{\mathsf T}),\qquad C=U'B,                       \tag{2.6}
\]

after relabelling the pivot rows as \(0,1\).  Thus this quotient chart has
six \(a\)-coordinates and ten \(b\)-coordinates, exactly

\[
 6+10=16=\dim X_{4,2}.                                   \tag{2.7}
\]

There is no residual internal \(\operatorname {GL}_2\)-gauge.  Exact rank
two is imposed by localizing at one column minor

\[
 \Delta_{pq}=b_{0p}b_{1q}-b_{0q}b_{1p}\ne0.              \tag{2.8}
\]

The ten row-pair charts and ten column-minor opens cover the exact-rank-two
stratum.  No ambient \(3\times3\) determinantal ideal is needed.  On the
displayed chart the first moment already has the constant pivot

\[
\begin{aligned}
 \mu _1={}&24b_{00}+6b_{11}
 +4(a_{20}b_{02}+a_{21}b_{12})\\
 &+6(a_{30}b_{03}+a_{31}b_{13})
 +24(a_{40}b_{04}+a_{41}b_{14}),                         \tag{2.9}
\end{aligned}
\]

so \(b_{00}\) should be eliminated before any higher moment is formed.

The same beta transform treats mixed sequences without changing the
denominator.  For

\[
 M_{e,a,b}=\xi _1^a\xi _2^{e-a}z_1^bz_2^{e-b}
 \qquad(0\le a,b\le e),
\]

direct coefficient extraction gives

\[
 \boxed{
 \frac{\mathcal E_2(M_{e,a,b}F^m)}{(4m+e+1)!}
 =\operatorname {CT}_u\int_0^1
 u^{b-a}t^b(1-t)^{e-b}P_{U',B}(u,t)^m\,dt.}              \tag{2.10}
\]

Hence a low-degree multiplier search changes only the numerator of the
same relative period.  The checker verifies (2.10) exactly for every
monomial multiplier of degrees \(e=0,1,2\) and orders \(0\le m\le4\).

### 2.2 An exact-rank-two all-order control family

There is a useful fixed-flag control inside another quotient chart.  Put

\[
 U=\begin{pmatrix}
 0&0\\0&0\\1&0\\0&1\\a_{40}&a_{41}
 \end{pmatrix},\qquad
 B=\begin{pmatrix}
 b_{20}&b_{21}&0&0&0\\
 b_{30}&b_{31}&b_{32}&0&0
 \end{pmatrix},                                         \tag{2.11}
\]

and localize at \(b_{21}b_{32}\ne0\).  The coefficient matrix \(C=UB\)
then has exact rank two, and every nonzero entry satisfies \(i>j\).  In the
relative period, the corresponding monomial has \(u\)-degree \(j-i\le-1\).
Over the function field

\[
 K=\mathbb Q(a_{40},a_{41},b_{20},b_{21},b_{30},b_{31},b_{32})
\]

with \(b_{21}b_{32}\) inverted, this proves

\[
 \nu_m:=\frac{\mu_m}{(4m+1)!},\qquad
 \operatorname {CT}_u P^m=0\quad(m\ge1),\qquad
 \boxed{\nu_{m+1}=0\quad(m\ge0).}                        \tag{2.12}
\]

Thus the scalar recurrence has forward coefficient one and no exceptional
integer step; its initial condition is \(\mu _1=0\).  The valuation argument
is stronger here than running a general creative-telescoping reduction.

For the multiplier in (2.10), the numerator has \(u\)-degree at most \(e\).
Therefore

\[
 \boxed{
 \mathcal E_2(M_{e,a,b}F^m)=0\qquad(m>e).}               \tag{2.13}
\]

The exact searches at \(e=1,2\) do find nonzero low-order mixed values,
but every sequence terminates at the bound (2.13).  This seven-parameter
exact-rank-two family is consequently SIC-safe.  It lies in the known
one-sided nullcone and is a control family, not a semistable component of
the unresolved rank-two moment fiber.

### 2.3 Exact component exhaustion on two separated rows

A larger two-channel chart can be closed without passing to a recurrence.
Fix the row pivot \(U=(e_0,e_4)\) and take the dense coefficient support

\[
 \{0\}\mathbin{\times}\{1,2,3,4\}\ \cup
 \{4\}\mathbin{\times}\{0,1,2,3\}.                       \tag{2.14}
\]

The internal \(\operatorname {GL}_2\)-gauge is already gone.  Overall
scaling and the contraction-preserving diagonal torus have coefficient
weights \((1,i-j)\).  The anchor weights \((1,-1)\) and \((1,1)\) are
independent up to a finite isogeny, so over the algebraic closure the
dense orbit can be normalized by

\[
 c_{01}=c_{43}=1.                                       \tag{2.15}
\]

Write the other six coefficients as \(z_0,\ldots,z_5\) in row-major
order and saturate by \(z_0\cdots z_5\).  In particular

\[
 \det C_{\{0,4\},\{0,4\}}=-c_{04}c_{40}=-z_2z_3\ne0,    \tag{2.16}
\]

so this is an exact-rank-two coefficient torus, not a calculation in the
ambient determinantal ideal.

The exact characteristic-zero moment scheme has the following sharp
finite-prefix profile:

\[
 \begin{array}{c|c}
  \text{equations}&\text{localized scheme over }\mathbb Q\\ \hline
  \mu_1,\ldots,\mu_7&\text{zero-dimensional of degree }604,\\
  \mu_1,\ldots,\mu_8&\text{empty (unit ideal).}
 \end{array}                                             \tag{2.17}
\]

The first row is a full rational-univariate calculation, not a collection
of modular specializations.  The second row proves that \(\mu_8\) removes
every exact component of the seven-moment scheme.  Consequently:

> **Proposition 2.1.** The dense two-separated-row chart (2.14) contains
> no all-order pure-moment point and hence no rank-two SIC
> counterexample.

There is therefore no surviving component on which to derive a relative
period recurrence or search for an infinite mixed tail.

The coordinate boundary is now closed as well.  The eight positions in
(2.14) have distinct diagonal-torus weights
\[
 \{-4,-3,-2,-1,1,2,3,4\}.                              \tag{2.18}
\]
Transpose and simultaneous reversal both negate the weights.  There are
135 proper support orbits under this involution, including the empty
support.  Their matrix strata are
\[
 \begin{array}{c|rrrrr}
  \text{stratum}&0&\text{one-sided rank one}&
  \text{mixed rank one}&\text{rank-one/rank-two}&
  \text{exact rank two}\\ \hline
  \text{orbit count}&1&15&2&3&114 .
 \end{array}                                             \tag{2.19}
\]
For each of the three coefficient tori containing both ranks, every
nonzero \(2\times2\) minor chart and the closed all-minors-zero locus are
localized separately over \(\mathbb Q\).

All 119 mixed proper support orbits have a sharp exact moment cutoff at
most seven:
\[
 \begin{array}{c|rrrrrr}
  \text{sharp cutoff}&2&3&4&5&6&7\\ \hline
  \text{orbit count}&56&17&15&8&18&5 .
 \end{array}                                             \tag{2.20}
\]
Sharp means that the preceding localized moment ideal is nonunit and
adjoining the displayed moment gives the unit ideal over \(\mathbb Q\).
The remaining proper supports are the zero form or one-sided rank-one
forms, hence SIC-safe.  Combining (2.17)--(2.20) gives:

> **Corollary 2.2.** Every form supported on the complete
> two-separated-row coordinate subspace (2.14) is SIC-safe.  Its only
> all-order pure-moment points are the zero form and one-sided rank-one
> forms.

Other row-pivot charts are not covered by this corollary.

## 3. Rank two does not supply a small generic cutoff

For a generic rank-two factor point, the exponent support of (2.2) has
Newton polygon
\[
 \operatorname {conv}\{(-4,0),(0,0),(4,4),(-4,4)\}. \tag{3.1}
\]
Its Euclidean area is \(24\), hence its normalized two-dimensional volume
is \(48\).  Thus the determinantal rank condition does not by itself
collapse the Laurent support to the six-point, normalized-volume-eight
support of the known full-rank witness.

The number \(48\) is not asserted to be the scalar recurrence order.
The period realization is relative to the endpoints \(t=0,1\), and the
coefficient scaling curve passes through the zero polynomial at \(s=0\).
The ordinary-point and scalar cyclic-vector gates from
[the holonomic algorithm note](HOLONOMIC_HYPERGEOMETRIC_ALGORITHMS.md)
therefore remain necessary.  In particular, the thirteen zero moments do
not propagate merely from holonomicity or from (3.1).

## 4. Mandatory order of work

The required order is:

1. parameterize the exact-rank-two open by requiring both factors in
   (2.1) to have rank two, equivalently by inverting suitable
   \(2\)-by-\(2\) minors;
2. impose a dimension-sized pure-moment system and compute an explicit
   closed point, or an explicit positive-dimensional component and its
   function field, inside that exact-rank-two open;
3. specialize (2.5), produce a checkable
   creative-telescoping recurrence, and evaluate \(\mu_{14}\) and any
   later bridge values it requires;
4. only after that characteristic-zero certificate exists, store primitive
   integral period/raw-moment companion forms, their singular-step ledgers,
   and selected good-prime \(p\)-curvatures.  Prime-power conclusions must
   use the uncancelled step ledger rather than only the generic curvature
   characteristic polynomial.

If the recurrence propagates zero, a fixed mixed multiplier must still be
tested before the point becomes an SIC witness.  If a tail moment is
nonzero, that exact value is the desired obstruction for that closed
point; excluding rank two globally would require treating every
semistable component.

Step 4 is calibrated by
[`TWO_PAIR_SIC_FROBENIUS_CURVATURE_BRIDGE.md`](TWO_PAIR_SIC_FROBENIUS_CURVATURE_BRIDGE.md).
It does not unpark the present scalar recurrence: the missing exact closed
point or component in step 2 remains the prior obstruction.

Independently, certifying the remaining squarefree rank-one target-only
membership
\[
\lambda^4(\lambda-1)^4
\bigl(p(8c-3d^2)\bigr)^5
\in(f_3,f_4,f_5,f_6)
\]
would sharpen the finite-prefix Hilbert result by forcing its existential
semistable point off rank one.  The exponent \(5\) is the common least
exponent modulo \(101,103,107\), but the rational lift remains open.  This
classification problem no longer precedes steps 1--3.

## 5. Reproduction

Run

```bash
python3 scripts/verify_two_pair_sic_bidegree44_rank_two_all_order_audit.py
python3 scripts/verify_two_pair_sic_bidegree44_rank_two_direct_chart.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree44_two_row_channel.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree44_two_row_boundaries.py
```

The dependency-free checker verifies the factor identity (2.3), the
beta/constant-term identity (2.4) through order four at the displayed
exact rank-two chart point, the four nonzero values in (1.1), and the
Newton polygon and normalized volume in (3.1).  The finite replay audits
why recurrence work is parked; it is not an all-order recurrence
certificate.  The second checker verifies the gauge quotient (2.6)--(2.9),
the mixed relative period (2.10), and the all-order pure and mixed
valuation certificates (2.12)--(2.13).  The last checker constructs the
degree-\(604\) exact scheme in (2.17) and certifies the eighth-moment unit
ideal over \(\mathbb Q\).  The boundary checker proves
(2.18)--(2.20), including the separate rank-one and exact-rank-two
localizations on the three mixed-rank coefficient tori.
