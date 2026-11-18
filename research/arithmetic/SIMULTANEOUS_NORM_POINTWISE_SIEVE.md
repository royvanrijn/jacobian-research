# Simultaneous norms: pointwise sieve upper bounds and a fixed-prime criterion

The forbidden-prime sieve has a pointwise upper bound at the natural
candidate-count scale. Its large-prime tail can therefore be made small by
choosing one fixed prime cutoff, provided a matching global lower bound is
available. The remaining distribution input then concerns each fixed prime,
rather than primes in a growing range.

The upper bounds below are unconditional consequences of Henriot's corrected
upper-bound theorem. The lower bound and fixed-prime marginal estimates in the
conditional criterion remain open. This note does not prove the every-target
binary conjecture or the original additive \(o(N)\) formulation of H2.

## Definitions

Write

\[
 G(x,y)=x^2+y^2,\qquad E(u,v)=u^2-uv+v^2,
\]

and let \(r_G(n),r_E(n)\) count all ordered integer representations by these
forms. Put \(h(n)=r_G(n)r_E(n)/24\) for positive \(n\), and define
\[
C(N)=\sum_{a+b=N}h(a)h(b),\qquad R(N)=576C(N).
\]
All component values \(a,b\) are positive. For forbidden \(p\ge5\), let \(B_p(N)\) count
coordinate witnesses with \(p\mid ab\). Define
\[
\mathcal T(N)=
 \prod_{\substack{q\mid N\\q\equiv1\ (12)}}(1+9/q)
 \prod_{\substack{q\mid N\\q\ge5,\ q\not\equiv1\ (12)}}(1+1/q).
\tag{S}
\]
Every factor is at least one. The target is always \(N\equiv2\pmod{12}\),
and all implied constants are independent of that target unless otherwise
specified. The function \(h\) is multiplicative, with

\[
\begin{array}{c|c}
 p\text{ class}&h(p^e)\\\hline
 p\equiv1\ (12)&(e+1)^2\\
 p\equiv5,7\ (12)&\mathbf1_{2\mid e}(e+1)\\
 p\equiv11\ (12),\quad p=2,3&\mathbf1_{2\mid e}.
\end{array}
\tag{Eul}
\]

For forbidden primes, at least one of the quadratic norm forms is inert, so a
positive component value divisible by \(p\) is divisible by \(p^2\).
For the admissible target class, both component values are automatically
\(1\pmod{12}\); the argument and small-prime local counts are in the
[local-density note](SIMULTANEOUS_NORM_LOCAL_DENSITIES.md).
After all forbidden primes have been excluded, the surviving values lie in

\[
 C_{12}=\{c\ge1:p\mid c\Longrightarrow p\equiv1\pmod{12}\}.
\]

**Upper-bound theorem.** There is a constant \(C>0\) such that
\[
R(N)\le C N\mathcal T(N),\qquad
B_p(N)\le C\frac{N}{p^2}\mathcal T(N)
       \quad(p\le N^{1/8},\ p\text{ forbidden}).
\tag{U}
\]
Consequently, for every \(Y\ge2\),
\[
\sum_{p>Y}B_p(N)
\ll \frac{N\mathcal T(N)}Y+N^{15/16}.
\tag{T}
\]
The same bounds hold after imposing any balance condition, since its bad
counts are subsets of these counts. The first bound on \(R\) is unrestricted.

## Exact rescaling of a bad component

For \(p\equiv11\pmod{12}\), \(h(p^2m)=h(m)\). For \(p\equiv5,7\),
\(h(p^2m)\le3h(m)\). These follow term by term from the Euler factors:
local values at exponent \(2j\) are respectively \(1\) and \(2j+1\),
and odd exponents vanish. Therefore, when \(p\nmid N\),
\[
B_p(N)\le 3456\sum_{1\le m<N/p^2}h(m)h(N-p^2m).
\tag{D}
\]
The constant is \(576\times2\times3\): either component can be bad.
When \(v_p(N)=1\), the count is exactly zero. When \(p^2\mid N\), both
values are divisible by \(p^2\), and
\[
B_p(N)\le9 R(N/p^2).
\tag{D'}
\]
For \(p\equiv11\), this last identity is exact with constant one.
All these statements include singular coordinate representations.

## The sharp Henriot calculation directly from corrected Theorem 5

Use Henriot's Theorem 5 in
[*Nair--Tenenbaum bounds uniform with respect to the discriminant*,
Math. Proc. Camb. Phil. Soc. 152 (2012), 405--424](https://doi.org/10.1017/S0305004111000752),
as replaced by **New Theorem 5**, p.377, in the
[2014 erratum, Math. Proc. Camb. Phil. Soc. 157, 375--377](https://doi.org/10.1017/S0305004114000280).
The original definitions and hypotheses are also available in
[arXiv:1102.1643v1](https://arxiv.org/html/1102.1643v1#S3).
We use the corrected \(\breve\rho\), not the original \(\hat\rho\), and
do not appeal to a sharpness or lower-bound theorem.

For \(k=1\), or \(k=p^2\le N^{1/4}\) with forbidden \(p\nmid N\), take
\[
Q_1(t)=t,\qquad Q_2(t)=N-kt,
\qquad F(u,v)=h(u)h(v).
\]
The two factors and their product are primitive, their total degree is two,
and their product coefficient norm is \(\ll N\).
The multiplicative function \(F\) is its own minimal majorant. At every prime,
\(h(q)=c(q)=4\mathbf1_{q\equiv1(12)}\), \(h(q^e)\le(e+1)^2\le4^e\),
and \(h(n)\le\tau(n)^2\ll_\epsilon n^\epsilon\). Thus
\(F\in\mathcal M_2(4,B_\epsilon,\epsilon)\) for every fixed sufficiently
small positive \(\epsilon\).

Corrected Theorem 5 bounds an interval sum by a constant times
\[
y\prod_{2<q\le x}(1-\rho_Q(q)/q)
\sum_{n_1n_2\le x}h(n_1)h(n_2)
\frac{\breve\rho(n_1,n_2)}{[n_1\operatorname{rad}n_1,n_2\operatorname{rad}n_2]}.
\]
All terms are nonnegative. Every prime divisor of \(n_1n_2\le x\) is at
most \(x\). Drop the product bound while retaining that prime support and
factor the resulting larger sum, using the multiplicativity proved in the
erratum. The resulting Euler factor is
\[
E_q(N,k)=\sum_{a,b\ge0}h(q^a)h(q^b)P_q(a,b).
\]
Here \(P_q(0,0)=1\). If an exponent is positive, its linear value has
that exact valuation; if one exponent is zero and the other positive, the
zero-exponent value is required to be a unit. These are precisely the
cross-coprimality conditions defining \(\breve\rho\). The normalized
congruence count is Haar probability modulo \(q^{\max(a,b)+1}\).
These terms are not a partition of the unit mass because the zero vector
has no conditions.

For \(q\nmid Nk\), the roots of the two linear values are distinct,
so only vectors with exactly one positive exponent occur. The probability
of exact valuation \(r\ge1\) for either factor is
\((1-1/q)q^{-r}\). Consequently
\[
E_q(N,k)=1+2c(q)/q+O(q^{-2}).
\]
For \(q\mid N\), \(q\nmid k\), the corrected cross-coprimality kills
vectors \((a,0),(0,b)\) with \(a,b>0\). The vector \((1,1)\) has probability
\((q-\nu_q)/q^2\), with \(\nu_q=2\) if \(v_q(N)=1\), and \(1\) otherwise.
All vectors with \(r=\max(a,b)\ge2\) contribute in total
\[
\ll\sum_{r\ge2}(2r+1)(r+1)^4q^{-r}\ll q^{-2},
\]
uniformly in \(v_q(N)\), since at least one unit-slope linear equation
specifies a residue modulo \(q^r\). Thus
\[
E_q(N,k)=1+c(q)^2/q+O(q^{-2})\qquad(q\mid N,\ q\nmid k).
\tag{E}
\]
At the variable leading-coefficient prime \(p\nmid N\), the second value
is always a unit modulo \(p\). Only vectors \((a,0)\), \(a>0\), occur,
and \(h(p)=0\), so \(E_p=1+O(p^{-2})\).
The factor at 2 is bounded uniformly by the same convergent valuation sum;
it is outside the first product. The finite factor at 3 is also harmless.

The universal nonexceptional baseline at \(q>2\) is
\[
A_q=(1-2/q)(1+c(q)/q)^2.
\]
Its partial products are bounded above and below by positive constants:
\(\log A_q=2(c(q)-1)/q+O(q^{-2})\), and fixed-progression Mertens gives
\(\sum_{q\le x}(c(q)-1)/q=O(1)\). Comparing (E) with this baseline at
\(q\mid N\), \(q\le x\), gives
\[
\frac{(1-1/q)E_q(N,k)}{A_q}
=\begin{cases}
1+9/q+O(q^{-2}),&q\equiv1\pmod{12},\\
1+1/q+O(q^{-2}),&q\text{ forbidden}.
\end{cases}
\tag{F}
\]
These ratios have product \(\ll\mathcal T(N)\). The bounds are uniform
in the target valuation: the higher-valuation sums above are uniform, and
finitely many small primes can be absorbed by their uniformly bounded
Euler factors. At \(p\mid k\), the root count is one, and the ratio to
the bad-prime baseline is bounded uniformly for \(p\ge5\).
Only primes at most \(x\) occur in this direct Euler majorant. Primes of
\(N\) above \(x\) cause no loss because their factors in \(\mathcal T\)
are at least one.

**Why the erratum matters.** With the original \(\hat\rho\), a zero
exponent does not impose the cross-coprimality above. At \(q\mid N\),
vectors \((1,0),(0,1)\) would then contribute spurious extra first-order
terms. The coefficient \(9/q\) in (F) uses the corrected theorem.

## Uniform interval parameters, including the endpoints

We prove
\[
\sum_{1\le m<N/k}h(m)h(N-km)\ll (N/k)\mathcal T(N)
\tag{C}
\]
for the \(k\)'s above. Put \(M=N/k\). First remove the two regions
\(km\le N^{1/2}\) or \(N-km\le N^{1/2}\). They contain
\(O(N^{1/2}/k+1)\) integers. Taking a fixed small \(\epsilon<1/100\),
their total is
\[
\ll (N^{1/2}/k+1)N^{2\epsilon}\ll M,
\]
since \(k\le N^{1/4}\).

Cover the remaining left half by the bands
\(km\in(N/2^{j+2},N/2^{j+1}]\), stopping after the first band whose
lower endpoint is at most \(N^{1/2}\). Cover the right half by the
same bands for \(N-km\). The last bands extend only to
\(N^{1/2}/2\), so all values remain positive and no short truncated
interval is introduced. Each resulting interval \((x,x+y]\) has
\(y\le x\), \(x\gg N^{1/4}\), and \(y\gg N^{1/4}\). Their total
length is \(O(M)\). As \(x\le N\), these inequalities imply
\(y>x^{1/8}\) and \(x\ge C_0\|Q\|^{1/8}\) for large \(N\).
Thus we may use Henriot's parameters \(\alpha=\delta=1/8\), then choose
the growth \(\epsilon\) small enough to meet the theorem's numerical
restriction

\[
 0<\epsilon<\frac{\alpha}{50g(g+\delta^{-1})}=
 \frac1{8000};
\]

for example \(\epsilon=1/20000\) works with \(g=2\).
Rounding intervals changes \(O(\log N)\) points and costs
\(O(N^{2\epsilon}\log N)=o(M)\), or can be avoided with half-open real
endpoints. The corrected upper bound and (E)--(F) give
\(\ll y\mathcal T(N)\) on each interval. Summing proves (C).

For \(k=1\), (C) proves the bound on \(R\). For \(p\nmid N\), apply
(C) with \(k=p^2\) in (D). For \(p^2\mid N\), use (D') and the bound on
\(R(N/p^2)\), observing \(\mathcal T(N/p^2)\le\mathcal T(N)\).
This proves (U). All arguments are pointwise in the target.

For \(Y<p\le N^{1/8}\), sum (U) and use \(\sum_{m>Y}m^{-2}\ll1/Y\).
For larger primes use the elementary sharp tail, which we recall here.
Each bad prime satisfies \(p^2\mid a\) or \(p^2\mid b\); hence there are at most
\(2\lfloor N/p^2\rfloor\) positive additive pairs. By (Eul), the coordinate
multiplicity \(576h(a)h(b)\ll_\varepsilon N^\varepsilon\). Therefore

\[
 \sum_{p>P}B_p(N)\ll_\varepsilon
 N^\varepsilon\sum_{p>P}\frac N{p^2}
 \ll_\varepsilon \frac{N^{1+\varepsilon}}P.
\tag{ET}
\]

This includes primes dividing the target and all singular coordinate tuples.
Taking \(P=N^{1/8}\), \(\varepsilon=1/16\) gives
\(\sum_{p>N^{1/8}}B_p(N)\ll N^{15/16}\). This proves (T), including
\(Y>N^{1/8}\) by the same existing tail.

## Consequence: fixed-prime marginal estimates can suffice

Fix the unbalanced positive interval or a fixed balance interval \(I\).
Suppose the remaining pointwise inputs are
\[
R_I(N)\ge cN\mathcal T(N),\qquad
\frac{B_{p,I}(N)}{R_I(N)}-\beta_p(N)\longrightarrow0
\quad\text{for every fixed forbidden }p,
\tag{H*}
\]
uniformly along all admissible targets. The marginal statement includes
varying \(v_p(N)\); it is not just a statement on unit targets.

Then (T), divided by the lower bound in (H*), gives a tail of size
\(O(1/Y)+o(1)\) relative to \(R_I(N)\). For each fixed \(Y\), sum the
finitely many marginal errors and then let \(Y\to\infty\). Therefore
\[
\sum_pB_{p,I}(N)
\le\left(\sum_p\beta_p(N)+o(1)\right)R_I(N),
\]
and the existing local bound yields
\[
R_{I,\mathrm{good}}(N)
\ge(10-\pi^2-o(1))R_I(N)>0.
\tag{G}
\]
No error estimate uniform in a growing prime range is assumed in (H*).
This is an alternative sufficient route, not a proof of (H*) and not the
original additive \(o(N)\) version of H2. Bare \(R_I\gg N\) does not
cancel the unbounded factor \(\mathcal T(N)\), even on reduced targets.

## The normalization is the natural local one

The next calculation concerns the raw coordinate densities, not the unweighted
indicator of \(C_{12}\). It shows that the target factor in (H*) has the
expected local size for this particular candidate count.

Fix \(q\ge5\), write \(x=q^{-1}\), and put \(f_q(e)=h(q^e)\). The
four-coordinate density for a fixed nonzero common value \(t\in\mathbb Z_q\)
is

\[
 w_q(t)=\kappa_q f_q(v_q(t)),\qquad
 \kappa_q=
 \begin{cases}
 (1-x)^2,&q\equiv1\ (12),\\
 1-x^2,&q\equiv5,7\ (12),\\
 (1+x)^2,&q\equiv11\ (12).
 \end{cases}
 \tag{L1}
\]

Indeed, a split binary norm represents a value of valuation \(e<k\) modulo
\(q^k\) in \((e+1)(q-1)q^{k-1}\) ways. An inert binary norm has no
representations for odd \(e\), and \((q+1)q^{k-1}\) for even \(e\).
Multiplying the two counts and normalizing by \(q^{2k}\) proves (L1).

Normalize Haar measure on \(\mathbb Z_q\) to have mass one. For
\(r=v_q(N)\), the raw eight-coordinate, three-equation density is

\[
 \alpha_q(N)=\kappa_q^2 I_{q,r},\qquad
 I_{q,r}=\int_{\mathbb Z_q}
 f_q(v_q(a))f_q(v_q(N-a))\,da.
 \tag{L2}
\]

This formula includes singular lifts. To see that the normalized finite
counts converge to this integral, away from the residue balls of \(0,N\)
the split and inert formulas just used agree with the integrand for
\(k>r\). At the zero residue a split norm has
\([q+k(q-1)]q^{k-1}\) representations and an inert norm has
\(q^{2\lfloor k/2\rfloor}\). The normalized common-value block is thus
\(O_q((k+1)^2)\). The other block has fixed valuation \(r\), so the
two exceptional balls contribute \(O_{q,r}(k^2q^{-k})\), tending to zero.
The integral tails converge by \(f_q(e)\le(e+1)^2\) and geometric decay.

Partitioning by the two valuations gives the exact formula

\[
\begin{split}
 I_{q,r}={}&(1-x)\sum_{0\le e<r}f_q(e)^2x^e
       +(1-2x)f_q(r)^2x^r\\
       &+2(1-x)f_q(r)\sum_{e>r}f_q(e)x^e.
\end{split}
\tag{L3}
\]

For the first sum the two valuations both equal \(e<r\). For the middle
term both equal \(r\), excluding the two leading unit residues that
would make either component have higher valuation. The final term covers
the two possibilities in which one valuation exceeds \(r\).
This also proves the formula for odd target valuations, where forbidden-prime
coefficients at odd exponents vanish.

For a unit target, (L3) becomes

\[
 I_{q,0}=1-2x+2(1-x)\sum_{e\ge1}f_q(e)x^e.
 \tag{L4}
\]

The three generating functions are

\[
 \sum_{e\ge0}f_q(e)x^e=
 \begin{cases}
 (1+x)/(1-x)^3,&q\equiv1\ (12),\\
 (1+x^2)/(1-x^2)^2,&q\equiv5,7\ (12),\\
 (1-x^2)^{-1},&q\equiv11\ (12).
 \end{cases}
\]

Consequently

\[
 \alpha_q(1)=
 \begin{cases}
 (1-x)^2(1+4x-x^2),&q\equiv1\ (12),\\
 (1-x)(1-x+3x^2+x^3),&q\equiv5,7\ (12),\\
 (1-x)(1+x)^3,&q\equiv11\ (12).
 \end{cases}
 \tag{L5}
\]

Here the notation \(\alpha_q(1)\) denotes the common density for every
\(q\)-adic unit target; no unit class has been suppressed.

Put \(c=f_q(1)\in\{0,4\}\). Uniformly for every \(r\ge1\), (L3) gives

\[
 I_{q,r}=1+(c^2-1)x+O(x^2),\qquad
 I_{q,0}=1+(2c-2)x+O(x^2).
 \tag{L6}
\]

For \(r=1\) the middle term supplies \(c^2x\); for \(r\ge2\) the
\(e=1\) term does so. All other terms have exponent at least two, with
absolute total bounded by a constant times
\(\sum_{e\ge2}(e+1)^4x^e\ll x^2\). In the final term, one may use
\(f_q(r)f_q(e)\le(e+1)^4\) for \(e>r\). These bounds are independent
of \(r\) and \(q\ge5\). Dividing the two expansions proves

\[
 \frac{\alpha_q(N)}{\alpha_q(1)}
 =1+\frac{(c-1)^2}{q}+O(q^{-2})
 =\begin{cases}
 1+9/q+O(q^{-2}),&q\equiv1\ (12),\\
 1+1/q+O(q^{-2}),&q\equiv5,7,11\ (12),
 \end{cases}
 \tag{L7}
\]

uniformly for every positive target valuation.

To control the finitely many small primes uniformly as well, the unit--unit
stratum gives \(I_{q,r}\ge1-x\) for \(r\ge1\), while
\(I_{q,0}\ge1-2x\). Cauchy--Schwarz gives

\[
 I_{q,r}\le(1-x)\sum_{e\ge0}f_q(e)^2x^e
 \le(1-x)\sum_{e\ge0}(e+1)^4x^e.
 \tag{L8}
\]

These are positive finite bounds independent of the target valuation.

Let \(\chi_{12}\) be the primitive quadratic character of discriminant
12. Formula (L5) gives
\(\alpha_q(1)=1+2\chi_{12}(q)/q+O(q^{-2})\).
Define the positive baseline product by the convergent normalization

\[
 A_0=L(1,\chi_{12})^2
       \prod_{q\ge5}\alpha_q(1)(1-\chi_{12}(q)/q)^2.
 \tag{L9}
\]

The product after the \(L\)-value is absolutely convergent and strictly
positive, since its factors are \(1+O(q^{-2})\).
The original unit-target product has this value, by fixed-character
prime-number/Mertens theory. At 2 and 3 the admissible target densities are
\(9/4\) and \(64/27\), respectively, as proved in the
[local-density note](SIMULTANEOUS_NORM_LOCAL_DENSITIES.md#why-these-are-the-relevant-prime-restrictions).
Thus the full predicted raw singular series is

\[
 \mathfrak S_h(N)=\frac94\frac{64}{27}\, A_0
          \prod_{\substack{q\mid N\\q\ge5}}
                  \frac{\alpha_q(N)}{\alpha_q(1)}.
 \tag{L10}
\]

Combining (L7)--(L8), absorbing finitely many small primes by their uniform
positive bounds and summing \(O(q^{-2})\) over the remaining primes, gives

\[
 c_1\mathcal T(N)\le\mathfrak S_h(N)\le c_2\mathcal T(N)
 \tag{L11}
\]

with absolute \(0<c_1\le c_2<\infty\) for all admissible targets.
This proves the asserted local normalization. It does not transfer the
local densities to any global integral-point count.

## Uniform truncation of deep valuations

First, an exact rescaling controls a common deep valuation. For a forbidden
\(p\) and integer \(L\ge1\), let \(D_{p,L}(N)\) count coordinate witnesses
for which \(p^{2L}\mid a,b\). This count is zero unless
\(p^{2L}\mid N\). When the target is divisible by that power, (Eul) gives

\[
 h(p^{2L}m)\le(2L+1)h(m),\qquad
 D_{p,L}(N)\le(2L+1)^2R(N/p^{2L}).
 \tag{V1}
\]

For \(p\equiv11\pmod{12}\), one has the exact identity
\(h(p^{2L}m)=h(m)\). The global upper bound in (U), with its constant
enlarged to cover bounded targets, therefore gives

\[
 D_{p,L}(N)\ll(2L+1)^2\frac N{p^{2L}}\mathcal T(N).
 \tag{V2}
\]

An individual component can also have a deep valuation while the other does
not. The following bound controls this additional possibility uniformly.
Let \(\widetilde D_{p,L}(N)\) count witnesses with
\(p^{2L}\mid a\) **or** \(p^{2L}\mid b\).

**Deep-component bound.** With an absolute implied constant, independent of
\(p,L,N\),

\[
 \widetilde D_{p,L}(N)
 \ll(2L+1)^2\frac N{p^{2L}}\mathcal T(N)
 \qquad(N\ge p^{8L}).
 \tag{V3}
\]

This includes every target valuation, including \(v_p(N)\) tending to
infinity. It is an upper bound, not an equidistribution assertion.

To prove it, first extend the primitive dilated correlation (C) to every
\(k=p^{2j}\), not just \(p^2\): if \(p\nmid M\),
\(M\equiv2\pmod{12}\), \(j\ge1\), and \(p^{2j}\le M^{1/4}\), then

\[
 \sum_{1\le m<M/p^{2j}}h(m)h(M-p^{2j}m)
 \ll \frac M{p^{2j}}\mathcal T(M).
 \tag{V4}
\]

Every step of the proof of (C) is uniform in this exponent. The polynomial
\(t(M-p^{2j}t)\) remains primitive with coefficient norm \(\ll M\).
At its exceptional leading-coefficient prime, the second factor is always
a unit and \(h(p)=0\); the calculation \(E_p=1+O(p^{-2})\) uses only
\(p\mid k\), not the exponent of \(p\) in \(k\). The endpoint bounds
and interval parameters use only \(k\le M^{1/4}\). Thus the same constant
works for all \(p,j,M\). Bounded targets can again be absorbed by increasing
that constant.

Put \(K=2L\) and \(r=v_p(N)\). If \(r\ge K\), a component divisible
by \(p^K\) forces the other to be divisible by the same power, so (V2)
applies. If \(r<K\) is odd, no such witness exists: a deep component
would force the other component to have the forbidden odd valuation \(r\).
For the remaining case \(r<K\) even, write

\[
 a=p^K m,\qquad b=p^r(M-p^{K-r}m),\qquad M=N/p^r,
 \qquad p\nmid M.
 \tag{V5}
\]

The second parenthesis is a unit. For \(p\equiv5,7\pmod{12}\), its
weight is therefore exactly
\(h(b)=(r+1)h(M-p^{K-r}m)\), while
\(h(a)\le(K+1)h(m)\). For \(p\equiv11\pmod{12}\), both factors
are one. The hypothesis \(N\ge p^{4K}\) ensures

\[
 p^{K-r}\le (N/p^r)^{1/4}=M^{1/4},
\]

because \(N\ge p^{4K-3r}\) is sufficient. Since \(K-r\) is even,
(V4) bounds the count in (V5) by

\[
 \ll(K+1)(r+1)\frac{M}{p^{K-r}}\mathcal T(M)
 \le(K+1)^2\frac N{p^K}\mathcal T(N).
\]

Reversing the two components gives the same bound for the other event.
This proves (V3), with constants uniform in all parameters. No extra factor
of \(L\) is introduced by the target valuation.

### A precise finite-stratum sufficient condition

For a fixed forbidden \(p\), one can now remove the unbounded valuation
from the global marginal problem, at a uniformly small cost. The remaining
finite-stratum global correlations are still open.

Fix \(K=2L\), and let \(R_I^{<K}(N)\) and \(B_{p,I}^{<K}(N)\) be
\(R_I(N),B_{p,I}(N)\) restricted to \(v_p(a),v_p(b)<K\). With the
same restriction on the local integral, put

\[
 f_{p,K}(e)=f_p(e)\mathbf1_{e<K},\qquad
 I_{p,r}^{<K}=\int_{\mathbb Z_p}
 f_{p,K}(v_p(a))f_{p,K}(v_p(N-a))\,da,
\]

\[
 \beta_p^{<K}(N)=1-\frac{g_r}{I_{p,r}^{<K}},\qquad
 g_r=\begin{cases}1-2/p,&r=0,\\1-1/p,&r\ge1.\end{cases}
 \tag{F1}
\]

The function \(f_p\) here is the forbidden-prime coefficient in (Eul):
it is \((e+1)\mathbf1_{2\mid e}\) for \(p\equiv5,7\pmod{12}\),
and \(\mathbf1_{2\mid e}\) for \(p\equiv11\pmod{12}\).
Both local denominators contain the unit--unit contribution \(g_r\),
which is at least \(3/5\).

Formula (L3), with \(f_p\) replaced by \(f_{p,K}\), shows that (F1)
depends only on the finite list of target strata

\[
 r=0,1,\ldots,K-1,\qquad r\ge K.
 \tag{F2}
\]

In the last stratum its denominator is simply
\((1-p^{-1})\sum_{e<K} f_p(e)^2p^{-e}\); in particular it does not
require an estimate for each individual deeper target valuation.

The corresponding local truncation error is uniformly bounded by

\[
 0\le\beta_p(N)-\beta_p^{<K}(N)
 \ll(K+1)^2p^{-K}.
 \tag{F3}
\]

For completeness, if \(r<K\) is odd, the removed local contribution is
zero. If \(r<K\) is even, exactly one component can have valuation at
least \(K\), and the removed contribution to \(I_{p,r}\) is

\[
 2f_p(r)(1-p^{-1})\sum_{e\ge K}f_p(e)p^{-e}
 \ll(K+1)^2p^{-K}.
\]

If \(r\ge K\), the two deep events coincide. Cauchy--Schwarz on the
ball \(p^K\mathbb Z_p\), preserved by \(a\mapsto N-a\), bounds their
contribution by

\[
 (1-p^{-1})\sum_{e\ge K}f_p(e)^2p^{-e}
 \ll(K+1)^2p^{-K}.
\]

The constants are absolute for \(p\ge5\). The unit--unit contribution
is unchanged, so substituting these bounds in (F1), with denominators at
least \(3/5\), proves (F3).

Suppose the natural lower bound \(R_I(N)\ge cN\mathcal T(N)\) in
(H*) holds. Here \(I\) is unrestricted or a fixed macroscopic balance
interval, and \(c>0\) is fixed. For all sufficiently large \(L\), (V3)
then implies that \(R_I^{<K}(N)>0\) for sufficiently large \(N\) and

\[
 \left|\frac{B_{p,I}(N)}{R_I(N)}-
       \frac{B_{p,I}^{<K}(N)}{R_I^{<K}(N)}\right|
 \ll_c(K+1)^2p^{-K}.
 \tag{F4}
\]

Indeed, let \(A\) be the number removed from \(R_I\). Every removed
point is also bad at \(p\), so the two truncated counts are \(R_I-A\)
and \(B_{p,I}-A\). Their ratio difference is nonnegative and at most
\(A/R_I\), since

\[
 \frac{B_{p,I}}{R_I}-\frac{B_{p,I}-A}{R_I-A}
 =\frac{A(R_I-B_{p,I})}{R_I(R_I-A)}\le\frac A{R_I}.
\]

The unrestricted bound (V3) dominates \(A\), proving (F4) uniformly in
the target valuation once \(N\ge p^{4K}\).

Consequently the following is a sufficient, finite-stratum replacement for
the fixed-prime marginal input in (H*): for each fixed sufficiently large
\(L\), require

\[
 \frac{B_{p,I}^{<2L}(N)}{R_I^{<2L}(N)}
       -\beta_p^{<2L}(N)\longrightarrow0
 \tag{F5}
\]

as \(N\to\infty\) in **each** of the finitely many strata (F2).
The final stratum includes every \(N\) divisible by \(p^{2L}\), with
no upper bound on \(v_p(N)\). Limits along all integers in each such
stratum, rather than just along a selected sequence, are required. A finite
number of these limits is uniform over their union. Combining (F3)--(F5)
and then letting \(L\to\infty\) proves the original fixed-prime marginal
limit in (H*). It suffices to establish (F5) for an unbounded sequence of
fixed \(L\). The order of limits is first \(N\to\infty\), then
\(L\to\infty\).

This reduction leaves a definite global problem. Define

\[
 C_{i,j,I}(N)=
 \sum_{\substack{p^{2i}u+p^{2j}v=N\\p\nmid uv\\
                 (p^{2i}u,p^{2j}v)\text{ satisfies }I}}
 h(u)h(v),\qquad u,v>0.
\]

Then

\[
 R_I^{<2L}(N)=576\sum_{0\le i,j<L}
 f_p(2i)f_p(2j)C_{i,j,I}(N),
 \tag{F6}
\]

and \(B_{p,I}^{<2L}\) is the same sum with the term \(i=j=0\)
omitted. Thus each fixed truncation is a finite family of binary weighted
correlations with fixed coefficients and unit restrictions. The pointwise
relative estimates needed for these correlations are not supplied by
Henriot's upper bound, local densities or the deep-valuation estimate.
In particular, (F5) remains an open global hypothesis; no fixed-prime
marginal is proved here beyond the already exact zero at \(v_p(N)=1\).

### What remains at the first fixed primes

At \(p=11\), and more generally \(p\equiv11\pmod{12}\), the local
coefficient is one on even valuations. Summing (L3) gives, for every
\(r\ge0\),

\[
 I_{p,r}=\frac{1-p^{-r-1}}{1+p^{-1}}.
\]

For \(r\ge1\), substituting the unit--unit contribution into (F1) without
truncation yields

\[
 \beta_p(N)=\frac{p^{-2}-p^{-r-1}}{1-p^{-r-1}}.
 \tag{F7}
\]

Since the exact unrestricted bad count is \(B_p(N)=R(N/p^2)\) whenever
\(r\ge2\), the desired marginal limit requires the pointwise relation

\[
 \frac{R(N/p^2)}{R(N)}
 -\frac{p^{-2}-p^{-r-1}}{1-p^{-r-1}}\longrightarrow0.
 \tag{F8}
\]

In particular, along targets with \(v_p(N)\to\infty\), it requires
\(R(N/p^2)/R(N)\to p^{-2}\). Neither the local formula nor the matching
upper and conjectural lower bounds establish this dilation regularity.

At \(p=5,7\), and more generally \(p\equiv5,7\pmod{12}\), define

\[
 u_p(m)=h(m/p^e)\mathbf1_{2\mid e},\qquad e=v_p(m).
\]

Then \(h(p^2m)=h(m)+2u_p(m)\) exactly. For positive \(M\), write
\(C_{f,g}(M)=\sum_{a+b=M}f(a)g(b)\), with positive summands and no
balance restriction. Swapping the summands gives
\(C_{h,u_p}=C_{u_p,h}\). When \(p^2\mid N\), exact rescaling therefore
gives

\[
 B_p(N)=576\left\{
 C_{h,h}(N/p^2)+4C_{h,u_p}(N/p^2)
                    +4C_{u_p,u_p}(N/p^2)\right\}.
 \tag{F9}
\]

These identities expose the unresolved issue already at the first fixed
primes: a pointwise dilation law at \(11\), and additional pointwise binary
correlations at \(5,7\). They supply no asymptotic for those correlations.

## Centered power windows: the same criterion for theta greater than one half

Fix \(1/2<\theta<1\), put \(H=N^\theta\), and restrict
\(|a-b|\le H\). Write \(R_H,B_{p,H}\) for the corresponding coordinate
counts. There is \(\eta_\theta>0\) such that, uniformly for \(Y\ge2\),
\[
\sum_{p>Y}B_{p,H}(N)
\ll_\theta \frac{H\mathcal T(N)}Y+H N^{-\eta_\theta}.
\tag{PT}
\]
This is a pointwise upper bound, not a lower count or a representation theorem.

Take \(\kappa=\theta/8\). For \(p\le N^\kappa\), the relevant interval
for the scaled variable \(m\) is centered at \(N/(2p^2)\), with length
\(\asymp H/p^2\). Its left endpoint \(x\asymp N/p^2\) and length
\(y\asymp H/p^2\) satisfy Henriot's hypotheses with fixed parameters,
for example short-interval exponent \(\alpha=\theta/4\) and
coefficient-height exponent \(\delta=1/4\):
\[
\log_N(y/x^\alpha)
\ge \theta-\alpha-2\kappa(1-\alpha)>0,
\qquad x\gg N^{1-2\kappa}\gg\|Q\|^{1/4}.
\]
Rounding endpoints is harmless as the interval has a fixed positive power
length. The corrected local-factor proof above gives
\[
R_H(N)\ll_\theta H\mathcal T(N),\qquad
B_{p,H}(N)\ll_\theta \frac H{p^2}\mathcal T(N)
\quad(p\le N^\kappa).
\tag{PU}
\]
For primes dividing \(N\), use the exact zero at valuation one, or scale
both component values by \(p^2\) and apply the same short-interval estimate
at target \(N/p^2\), width \(H/p^2\). The parameter inequalities are
uniform in \(p\le N^\kappa\).

For larger primes, a separate elementary window tail is useful. Each of
\(p^2\mid a\) and \(p^2\mid b\) selects at most \(H/p^2+1\) possible
additive pairs. Each pair has \(O_\varepsilon(N^\varepsilon)\) coordinate
multiplicity. No prime greater than \(\sqrt N\) occurs. Therefore
\[
\sum_{p>P}B_{p,H}(N)
\ll_\varepsilon N^\varepsilon(H/P+\sqrt N).
\tag{PW}
\]
This bound includes singular points and prime divisors of the target.
Set \(P=N^\kappa\), and choose
\(0<\varepsilon<\min(\kappa,\theta-1/2)\). Since \(\mathcal T\ge1\),
(PU)--(PW) prove (PT), with any positive
\(\eta_\theta<\min(\kappa-\varepsilon,\theta-1/2-\varepsilon)\).

Consequently the fixed-prime criterion (H*) remains sufficient in these
windows after replacing its lower bound by
\(R_H(N)\ge c_\theta H\mathcal T(N)\), and requiring the marginal limits
inside the same centered windows. These global pointwise inputs remain open.
No threshold improvement for the known semigroup density-one theorem follows
from this upper-bound argument alone.

## What is proved and what remains

The unconditional results are (U), the fixed-cutoff tail (T), the local
comparison (L11), the common- and individual-valuation bounds (V2)--(V3),
and the centered-window upper and tail bounds (PU)--(PT). The conditional implication (H*) to (G)
is proved as well, as is the conditional finite-stratum reduction
(F1)--(F6). No pointwise lower bound for \(R_I\) or \(R_H\), and no
fixed-prime marginal limit, is supplied here. Those are indispensable open
global inputs.

The imported analytic estimate was checked against the 2014 erratum's
actual New Theorem 5, including its cross-coprimality conditions and growth
restriction. Two separate internal derivation reviews examined that
specialization, the local Euler majorant, interval parameters and the
fixed-prime criterion. They do not constitute external mathematical review
or formal verification. The argument uses no new finite coverage computation;
the pre-existing local formula checker has its own scope in the local-density
note. No novelty claim is attached to these consequences of Henriot's bound.
