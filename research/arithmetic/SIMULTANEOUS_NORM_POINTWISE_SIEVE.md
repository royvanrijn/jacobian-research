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

## A uniform bound on a common deep valuation

There is a further exact truncation bound. For a fixed forbidden \(p\)
and integer \(L\ge1\), let \(D_{p,L}(N)\) count coordinate witnesses
for which \(p^{2L}\mid a,b\). This count is zero unless
\(p^{2L}\mid N\). When the target is divisible by that power, (Eul) gives

\[
 h(p^{2L}m)\le(2L+1)h(m),\qquad
 D_{p,L}(N)\le(2L+1)^2R(N/p^{2L}).
 \tag{V1}
\]

For \(p\equiv11\pmod{12}\), one has the exact identity
\(h(p^{2L}m)=h(m)\).
The global upper bound in (U), with its constant enlarged to cover bounded
targets, therefore gives

\[
 D_{p,L}(N)\ll(2L+1)^2\frac N{p^{2L}}\mathcal T(N).
 \tag{V2}
\]

Under the natural lower bound in (H*), this is a uniformly vanishing
relative contribution as \(L\to\infty\). This statement concerns a
common deep valuation of the two component values. It is not a proof of
fixed-prime equidistribution or of an unrestricted finite-level reduction.

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
comparison (L11), the common-valuation bound (V2), and the centered-window
upper and tail bounds (PU)--(PT). The conditional implication (H*) to (G)
is proved as well. No pointwise lower bound for \(R_I\) or \(R_H\), and no
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
