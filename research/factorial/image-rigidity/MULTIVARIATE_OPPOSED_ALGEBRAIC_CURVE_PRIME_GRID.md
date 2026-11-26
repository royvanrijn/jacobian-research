# Prime restrictions on algebraic opposed height curves

Fix an ordinary two-monomial core. Append \(r\) copies of
\((L,L+d)\) and \(s\) copies of \((Y+e,Y)\), where \(d,e\ge1\)
are fixed integers and the signed gap \(h=rd-se\) is fixed.
Assume the ordinary heights lie on a prescribed nonzero
\(F(L,Y)=0\) with \(F\in\mathbb Z[L,Y]\). Consider adjacent
orders \(m,m'\ge2\), either of which may be marked.

Mark the \(Y\) row after reversing both columns. For an integer
\(q\) with \(1\le q\le\min(m,m')\), choose any nonempty subset
\[
 T\subset\{t:1\le t\le qe,\ m\nmid t\}.
\]
Put \(x_t=mY+t\) and \(v_u=mL+u\), the latter for
\(1\le u\le md\). Define the admissible index set
\[
 U(T)=\{u:1\le u\le md,\ F(-u/m,-t/m)=0
                                      \ \text{for some }t\in T\}.
\]
Let \(b=\deg F\), let \(W\) be the least common multiple of the
nonzero integers \(\left|m^b F(-u/m,-t/m)\right|\) for those
\(u,t\), and take \(W=1\) if the list is empty. Use the marked
cutoff \(D_0=\operatorname{lcm}(1,\ldots,mm'e)\) and the fixed
integer \(B_C\) from the
[opposed mark-product theorem](MULTIVARIATE_OPPOSED_MARK_PRODUCT_HEIGHT_CONSTRAINT.md).
Let \(D_1=\operatorname{lcm}(D_0,W)\).

**Theorem.** At a common factor of degree at least \(q\), the
ordinary heights on \(F=0\) satisfy
\[
 \left(\prod_{t\in T}\frac{x_t}{\gcd(x_t,D_1)}\right)^s
 \ \bigm|\
 \left(\prod_{u\in U(T)}v_u\right)^r
 W^{rmd}D_0^{(m+1)s}B_C^{m+1}.                                 \tag{1}
\]
The list \(U(T)\) and all constants on the right can be computed
from the prescribed algebraic equation and the fixed data.
In particular, if one eligible mark has \(U(\{t\})=\varnothing\)
at \(q=2\), it excludes this entire algebraic height curve for
large \(Y\) at the marked pair, uniformly in the admissible
counts with that signed gap.

**Consequence for boundary branches.** Suppose \(r,s\to\infty\)
and \(Y\asymp L^\beta\) for \(\beta>0\). At any accumulating
sequence of ordinary collisions on \(F=0\), the first two
blocks at each marked order must satisfy
\[
                  |U(T)|\ \ge\ \frac{d}{e}\,\beta\,|T|
                              \quad\text{for every }T.          \tag{2}
\]
Thus a prescribed rational polynomial height trajectory
\(Y=P(L)\) of degree \(\beta\ge2\) cannot give a positive
diffuse limit of collisions when the gap of the taller row
satisfies \(e<d\). With fixed \(h\) this yields a computable
cutoff uniform over all sufficiently large counts for that
polynomial trajectory. At equal gaps, a necessary condition
is that every polynomial \(P(L)+t/m\), for an eligible
first-two-block offset, split into distinct linear factors
with roots of the form \(-u/m\) in the active factorial strip.

## Proof

For each \(t,u\), substitute \(L=(v_u-u)/m\) and
\(Y=(x_t-t)/m\) into \(m^bF\). This gives an integer
polynomial in \(v_u,x_t\) with constant
\(m^b F(-u/m,-t/m)\). At a prime \(p\) for which
\(v_p(x_t)>v_p(D_1)\), every unallowed index \(u\), meaning
\(F(-u/m,-t/m)\ne0\), has
\[
                              v_p(v_u)\le v_p(W).
\]
Otherwise both arguments of the integer polynomial
would be divisible by \(p^{v_p(W)+1}\). The equation would
force its nonzero constant to be divisible by the same
power, contrary to the definition of \(W\).

Fix a prime with \(v_p(x_t)>v_p(D_1)\). Such a prime is charged
for at most one mark \(t\in T\), since
\(p^{v_p(D_0)+1}>mm'e\) and distinct mark differences have
absolute value less than \(mm'e\).
The prime-by-prime estimate in the
[proof of the mark-product theorem](MULTIVARIATE_OPPOSED_MARK_PRODUCT_HEIGHT_CONSTRAINT.md#proof-of-the-product-condition),
applied to the reversed columns, gives at a degree-at-least-\(q\)
common factor
\[
 s\big(v_p(x_t)-v_p(D_0)\big)
 \le r\sum_{u=1}^{md}v_p(v_u)
          +(m+1)\big(sv_p(D_0)+v_p(B_C)\big).
\]
The left side bounds from above the valuation of
\((x_t/\gcd(x_t,D_1))^s\). On the right, the contributions
of unallowed indices total at most \(rmd\,v_p(W)\).
All the other indices lie in \(U(T)\). Summing the disjoint
prime exponents proves (1).

For \(q=2\), an ordinary common root implies a common factor
of degree at least two by even-order real positivity.
If \(U(\{t\})=\varnothing\), take logarithms of (1) and
divide by \(s\). Its left side is at least
\(\log(mY)-\log D_1\). Its right side is bounded independently
of the heights since \(s\ge1\) and \(r/s\le(e+|h|)/d\).
This proves the asserted exact exclusion.

More generally the same logarithmic inequality, before
division by \(\log L\), is
\[
 |T|\big(\log(mY)-\log D_1\big)
 \le \frac{r}{s}|U(T)|\log(m(L+d))
       +\frac{r}{s}md\log W
       +(m+1)\log D_0+\frac{m+1}{s}\log B_C .
\]
Using \(r/s\to e/d\) proves (2). If \(Y=P(L)\) is a rational
polynomial of degree \(\beta\), any eligible mark has at
most \(\beta\) admissible active offsets because
\(P(-u/m)+t/m=0\) has that degree. At \(e<d\) this
contradicts (2). For an explicit cutoff, first choose a
count threshold so that \(r/s\le (d+e)/(2d)<1\), and then
use \(P(L)\ge c L^\beta\) for an explicit positive rational
\(c\) beyond a computable \(L\). The last displayed
inequality gives a computable upper bound on \(\log L\)
at a collision. Since an absolute cutoff leaves a finite
set of adjacent pairs, taking their maximum gives an
all-degree threshold. At \(e=d\), condition (2) for a
singleton needs all \(\beta\) roots among the indicated
active offsets, proving the last assertion.
