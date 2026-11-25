# An obstruction to a fixed set of Newton polygon primes

Let
\[
P_{m,a,d}(t)=\sum_{k=0}^m c_{m,k}t^k,\qquad
c_{m,k}=\binom mk\frac{(am+dk)!}{(am)!}.
\]
For every finite set of primes \(S\), there are arbitrarily large \(n\)
for which the Newton polygons of \(P_{n,6,14}\) and \(P_{n+1,6,14}\)
at every \(p\in S\) admit degree two in the common-factor degree test.

In that test, a shared slope \(u/v\) in lowest terms contributes any
multiple of \(v\) up to the smaller of the two edge widths. Both polygon
multiplicativity and the integrality of the edge heights are required
for this necessary condition on a rational factor degree. Degree two
also obeys [real-root parity](BINOMIAL_RESIDUAL_LIFTING.tex) and, when
\(n>2\), the [strict degree bound](BINOMIAL_IRREDUCIBILITY_REDUCTION.tex). The example satisfies
the support and order lower bounds in the [residual region](ORDINARY_BINOMIAL_RESIDUAL_REGION.md)
at every \(N=n+1\ge25\). The construction supplies no upper bound
relative to the explicit cutoff \(N_*\). A degree admitted by this test is
only a possibility based on slope and width.

## Stabilizing a supporting line

Fix \(a,d>0\), write \(B=a+d\), and take \(v_p(p)=1\).
The extreme indices attaining
\[
L_{m,p}(\lambda)=\min_{0\le k\le m}\{v_p(c_{m,k})-\lambda k\}
\]
give the endpoints of the edge at slope \(\lambda\) when their difference
is positive. Write \(\mu_p=d/(p-1)\).

Suppose \(m_R\equiv m_0\pmod{p^R}\), where \(m_0\ge2\) is fixed,
\(m_R\to\infty\), and \(\log_p(Bm_R)=O(R)\). If \(\lambda<\mu_p\), the
edge at \(\lambda\) for \(m_0\) persists with the same width at \(m_R\)
for all sufficiently large \(R\). In fact the factorial interval of
length \(dk\) gives
\[
v_p(c_{m_R,k})\ge v_p((dk)!)
\]
by counting multiples of each power of \(p\) and using the integrality
of \(\binom{m_R}{k}\). Legendre's formula gives
\(v_p((dk)!)=\mu_p k+O(\log(k+1))\), so only a fixed bounded set of
indices can attain a value at most \(L_{m_0,p}(\lambda)\). For these
indices \(k\le m_0\), the valuations of the nonzero factors in the
binomial numerator and the factorial interval agree with their
values at \(m_0\) once \(R\) is large. At the remaining fixed indices
\(k>m_0\) the factor \(m_R-m_0\) in the binomial numerator forces the
coefficient valuation to tend to infinity. The extreme minimizing
indices therefore agree.

For \(\lambda>\mu_p\), use the coefficients in reverse order and divide
by the top coefficient. Put
\[
r_{m,j}=\frac{c_{m,m-j}}{c_{m,m}}
=\frac{\binom mj}{\prod_{i=0}^{dj-1}(Bm-i)},\qquad
M=\min_{0\le j\le m_0}\{v_p(r_{m_0,j})+\lambda j\}.
\]
The edge persists if the inequalities
\[
T(j)+\lambda j>M\qquad(j>m_0,\ dj>Bm_0)
\tag{1}
\]
hold, where
\[
T(j)=v_p(m_0!)+v_p((j-1-m_0)!)-v_p(j!)
-v_p(B)-v_p((Bm_0)!)-v_p((dj-1-Bm_0)!).
\tag{2}
\]
To prove sufficiency, the denominator interval in \(r_{m_R,j}\) contains
at most \(dj/p^\ell+1\) multiples of \(p^\ell\), and none for
\(p^\ell>Bm_R\). Hence
\[
v_p(r_{m_R,j})+\lambda j
\ge (\lambda-\mu_p)j-\lfloor\log_p(Bm_R)\rfloor.
\]
For large \(R\), indices \(j>CR\) for a sufficiently large fixed \(C\)
therefore have weighted valuation greater than \(M\). For all
\(j\le CR\), every nonzero linear factor in either product at \(m_0\)
has absolute value \(O(R)<p^R\) for large \(R\), so its valuation at
\(m_R\) is the same. The fixed denominator \(j!\) of the binomial
coefficient is retained exactly. For \(j\le m_0\) both products are
nonzero at \(m_0\). For \(m_0<j\le Bm_0/d\), only the numerator has a
zero at \(m_0\) and its valuation tends to infinity on this finite
interval. For \(j>m_0\) and \(dj>Bm_0\), the numerator factor
\(m_R-m_0\) cancels the denominator factor \(B(m_R-m_0)\) up to \(1/B\).
The other factors give (2). Consequently (1) leaves exactly the
same extreme minimizing reverse indices, with the same difference.

The test (1) is finite. Let \(\delta=\lambda-\mu_p>0\). Since
\[
\frac{m_0!(j-1-m_0)!}{j!}=\frac1{j\binom{j-1}{m_0}},
\]
the first three terms of (2) are at least
\(-2\lfloor\log_p j\rfloor\). Legendre's formula gives a contribution
of at most one for each power of \(p\) in the binomial valuation, and
\(v_p(j)\le\lfloor\log_p j\rfloor\). The last two factorials in
(2) have arguments whose sum is \(dj-1\), so their valuations sum
to at most \(dj/(p-1)\). Thus the left side of (1) is at least
\[
\delta j-2\lfloor\log_p j\rfloor-v_p(B).
\]
Choose \(e\ge0\) so that
\[
\delta p^e-2e-v_p(B)>M,\qquad \delta(p-1)p^e\ge2.
\]
For \(j\ge p^e\), each band \(p^h\le j<p^{h+1}\) has a lower-bound
minimum at least that of its predecessor, starting with \(h=e+1\).
It suffices to verify (1) for the remaining indices \(j<p^e\).

## Local certificates and the Chinese remainder theorem

For \((a,d)=(6,14)\), the rows give a lower order \(n_0\) and choices
of slope with width available at both orders \(n_0,n_0+1\). Each
selection is different from \(\mu_p\). Each width is a multiple of
its reduced slope denominator and the selected widths total two.

| \(p\) | \(n_0\) | \((\lambda,\text{width})\) |
|---:|---:|:---|
| 2 | 6 | \((27/2,2)\) |
| 3 | 15 | \((13/2,2)\) |
| 5 | 7 | \((4,2)\) |
| 7 | 2 | \((2,2)\) |
| 11 | 3 | \((3/2,2)\) |
| 13 | 5 | \((1,2)\) |
| 17 | 2 | \((1,2)\) |
| 19 | 2 | \((1,2)\) |
| 23 | 4 | \((0,1), (1,1)\) |
| 29 | 5 | \((0,1), (1,1)\) |
| 31 | 4 | \((1/2,2)\) |
| 37 | 12 | \((1/2,2)\) |
| 41 | 7 | \((0,2)\) |
| 43 | 8 | \((0,2)\) |

The [checker](verify_binomial_newton_prime_sets.py) constructs the actual
integer coefficients for both orders and checks the extreme minimizing
indices using rational arithmetic. For every selected slope above the
mean, it computes \(M\), finds \(p^e\) by the two displayed inequalities,
and checks every indispensable index in (1) below \(p^e\) by Legendre's
formula. Run the checker from the repository root.

```
python3 research/factorial/image-rigidity/verify_binomial_newton_prime_sets.py
```

The Python 3.12.3 replay checks 14 primes, 32 edges, and 2166 reverse values.
For every prime \(p>46\), take \(n_0=2,\lambda=0\) and width two. The
coefficients of indices 0 through 2 are units because the largest needed
factorial argument is \(3a+2d=46\) and the binomial coefficients are units.
All coefficients have nonnegative valuations, so this horizontal edge
has the required width and lies strictly below \(\mu_p\).

For a nonempty finite \(S\), let \(Q_R=\prod_{p\in S}p^R\).
Choose \(Q_R\le n_R<2Q_R\) with \(n_R\equiv n_0(p)\pmod{p^R}\)
for each \(p\in S\). Both \(n_R\) and \(n_R+1\) have logarithmic size
\(O(R)\) and satisfy the needed congruences for the two base orders.
The stabilization result applies simultaneously to the finitely many
selected slopes and retains the degree-two choices for every
sufficiently large \(R\). The empty set of primes imposes no degree
constraint.

## Exact residual gcds at all sufficiently large primes

For the fixed support \((a,d)=(9,10)\), every finite set of primes
\(p\ge73\) admits arbitrarily large common lower orders \(n\) such that
at each selected prime the edges of slope \(1/2\) have residual gcd
of degree one. The reduced slope denominator is two, so the test
using residual gcds allows degree two at every selected prime.

Fix such a prime and write \(B=19\). Let \(\delta\) be the representative
of \(p\bmod 10\) in \([-5,5]\), so \(0<|\delta|\le3\).
Choose \(q\in\{1,11\}\) with \(qp\equiv\delta\pmod{20}\), and put
\(\ell=(qp-\delta)/10\), an even positive integer. Choose an integer
\(i\in[1,9]\) with \(r=10+\delta-i\in[1,9]\). There are
\(9-|\delta|\ge6\) choices. Put
\[
 D=qp-i=10(\ell-1)+r,\qquad K=\ell/2-(q-1).
\]
In particular \(2\le K<p\). Given a unit \(u\bmod p\), choose its
positive representative \(u\le Bp\) subject to
\(p^K u+D\equiv0\pmod B\), and set \(Bm=p^Ku+D\).
Write \(u_*=p^Ku\), so \(v_p(u_*)=K\).

For each \(i\) write \(m\bmod p=(kp-i)/B\) with
\(1\le k\le B-1\); reduction modulo \(B\) determines \(k\) from
\(i\). Discard \(k=1,2\), and if \(q=11\) discard \(k=11\).
The low digits of \(\ell\) and \(\ell+2\) in base \(p\) are
\((p-\delta)/10\) and \((p-\delta)/10+2\). The low digit of \(m\)
exceeds \((p-\delta)/10+3\), since \(k\ge3\) and \(p\ge73\), and is below \(p-1\).
When \(q=11\) the digit of \(m\) in position \(p\) is nonzero.
Indeed \(Bm\equiv D\pmod{p^2}\) as \(K\ge2\), and if
\(m\bmod p^2=m\bmod p\) then \(B(m\bmod p)=D\), since both sides
are smaller than \(p^2\). This would force \(k=11\). It follows from
the base-\(p\) coefficient rule, obtained by expanding
\((1+z)^p=1+z^p\) in \(\mathbb F_p[z]\), that
\[
 v_p\binom m\ell=v_p\binom{m+1}{\ell+2}
 =v_p\binom{m+1}2=0 .
\]

At order \(s\in\{m,m+1\}\), the reverse edge is determined by the maxima of
\[
 -v_p(r_{s,j})-j/2
 =v_p\left(\prod_{v=Bs-10j+1}^{Bs} v\right)-v_p\binom sj-j/2
 \quad (0\le j\le s).
\]
For \(j<p^2/20\), the factorial interval can contain only one
multiple of \(p^2\), namely \(u_*\), because it lies within distance
less than \(p^2\) of \(u_*\) at both ends. The first multiple of
\(p\) enters the lower interval at \(1+\lfloor(p-i)/10\rfloor>2\)
and the upper interval at index \(2\), since \(10\le B-i<20\).
Successive multiples of \(p\) enter at indices separated by at least
\(\lfloor p/10\rfloor>2\). The deep multiple \(u_*\) enters at
indices \(\ell,\ell+2\), respectively, giving total factorial
valuations \(K+q-1=\ell/2\) and \(K+q=(\ell+2)/2\), respectively.
The corresponding binomial coefficients are units. Between the
jumps the factorial valuation is constant. The spacing of the
unit jumps, the initial values and the values at the deep jump
show that the maxima at indices \(j<p^2/20\) are zero, attained
exactly at \(0,\ell\) and \(0,2,\ell+2\), respectively.
For \(j\ge p^2/20\), counting multiples of each power of \(p\)
gives a factorial valuation at most
\(10j/(p-1)+\lfloor\log_p(Bs)\rfloor\). Here \(Bs<p^{K+2}\)
and \(K<p\), so this bound is strictly less than \(j/2\).
These are therefore exactly the reverse edge indices.

At an even edge index, reduce the weighted coefficient \(p^{j/2}r_{s,j}\)
and use the monomial \(Y^{j/2}\). Equivalently, use \(Y=z^2\) after
scaling the reciprocal variable by \(x=\sqrt p\,z\). The residuals are
\[
 1+(L/u)Y^{\ell/2},\qquad
 1+AY+(U/u)Y^{\ell/2+1}
 \quad\text{in }\mathbb F_p[Y],
\]
where \(L,U,A\) are nonzero and independent of \(u\bmod p\).
The high multiple contributes the factor \(u^{-1}\) at both long
indices. All other multiples of \(p\) in the corresponding
intervals are \(u_*+j p\), with unit parts \(j\bmod p\), and
the binomial coefficients here depend only on \(m\bmod p^2\).

For completeness, put \(W=Bm\). The factorial interval at upper
index \(\ell+2\) extends the one at lower index \(\ell\) by one
factor \(u_*+r-10\) at the bottom and the block \((W+1)_{19}\)
at the top. The interval at upper index \(2\) is precisely this
top block with the extra factor \(W\) at the bottom. The top
block contains one ordinary multiple of \(p\), and its normalized
unit cancels when comparing \(A\) with \(U/L\). The ratio of the
two binomial coefficients at the long indices is
\[
 \frac{\binom{m+1}{\ell+2}}{\binom m\ell}
   =\frac{(m+1)(m-\ell)}{(\ell+1)(\ell+2)} .
\]
Since \(W\equiv-i\) and \(u_*+r-10\equiv\delta-i\), equality
\(A=U/L\) would imply in \(\mathbb F_p\)
\[
 i(f-2)=\delta(f-2B/10),\qquad
 f=(1-\delta/10)(2-\delta/10).
\]
Here \(f-2=\delta(\delta-30)/100\ne0\), so at most one choice
of \(i\) is excluded. The other exclusions also each remove at
most one choice. Fix a remaining \(i\), put \(y=(U/L-A)^{-1}\),
and choose \(u\equiv -Ly^{\ell/2}\pmod p\). The lower residual
vanishes at \(y\). Subtracting \( (U/L)Y\) times the lower residual
from the upper yields \(1+(A-U/L)Y\), so the gcd is exactly linear.

The reverse stabilization above applies with \(\lambda=1/2\)
and with zero minimum. Its tail bound
\((1/2-10/(p-1))j-2\lfloor\log_p j\rfloor-v_p(B)\) is positive
for every \(j\ge p\), since \(\lfloor\log_p j\rfloor\le j/p\)
there and \(p\ge73\). Both base orders exceed \(p\).
At each fixed edge index the normalized \(p\)-unit of \(r_{s,j}\)
also stabilizes. Indeed all its nonzero linear factors at the base
order retain their valuations and their normalized units under
sufficiently deep \(p\)-adic congruences, and its factorial denominator
\(j!\) is fixed. The same CRT argument therefore stabilizes the
actual residuals at every selected prime.

The construction leaves primes below \(73\) unspecified and
imposes no upper bound on \(n\).

### First ramified lifting condition

The residual gcd need not lift to a shared root. For the same choices
of \(p,\delta,i,\ell\) above, put
\[
 m_0=-i/19,\qquad l_0=-\delta/10,\qquad
 \tau=\frac{38(m_0-l_0)}{(l_0+1)(l_0+2)(\delta-i)}
\]
and define the rational numbers
\[
 \begin{split}
 b_0&=\frac{m_0}{(-i-9)_{10}},&
 b_1&=\frac{m_0+1}{(10-i)_{10}},\\
 w_0&=\frac{m_0-l_0}{(l_0+1)(r-19)_{10}},&
 w_1&=\frac{m_0-1}{3(-i-10)_{10}},&
 w_2&=\frac{m_0-l_0-1}{(l_0+3)(r-20)_{10}} .
 \end{split}
\]
All denominators are units at the primes in the construction.
A common root of valuation \(1/2\) requires the reduction to vanish of
\[
 \Delta_{\delta,i}
   =-l_0\bigl((\tau-1)b_1+w_1-\tau w_2\bigr)
    -\bigl(2-(l_0+2)\tau\bigr)(b_0-w_0).
\]

To see this, write \(\pi=\sqrt p\) and \(x=\pi z\). In the two
scaled reverse polynomials, the zeroth order residues in \(\pi\) are
the displayed polynomials in \(Y=z^2\). The only reverse indices
of excess valuation \(1/2\) are \(1,\ell+1\) in the lower order
and \(1,3,\ell+3\) in the upper order. This follows by the same
jump spacing and tail bounds, with the respective extra indices
\(\ell+1,\ell+3\) still having unit binomial coefficients.
Taking the ratio of adjacent reverse coefficients shows that
the reduced coefficients at these indices, after dividing by
\(\pi\), are
\[
 (b_0,(L/u)w_0),\qquad (b_1,Aw_1,(U/u)w_2),
\]
respectively. These formulas use \(m\equiv m_0\), \(\ell\equiv l_0\),
and \(U/(LA)=\tau\).

A putative common root reduces to \(z_0\) with \(z_0^2=y\).
Choose a lift of \(z_0\) in the unramified extension so that the
zeroth order polynomials evaluated at it vanish modulo \(p\).
The derivative of the lower residue is a unit, since
\(\ell\not\equiv0\pmod p\) and \(z_0\ne0\). Taylor expansion first
forces the root to lie within \(\pi\) of the lift. Reducing the
two Taylor equations for its first correction modulo the maximal
ideal gives a common scalar correction times each zeroth order
derivative, plus the respective first order value.
Eliminate the correction and substitute
\((L/u)y^{\ell/2}=-1\) and \(y=(U/L-A)^{-1}\).
After division by \(Az_0^2\), this compatibility determinant is
\(\Delta_{\delta,i}\). The same conclusion holds at the stabilized
CRT orders. The reverse tail bounds used for stabilization are
strict by more than \(1/2\) outside the fixed finite set of
indices, so these first order coefficients also stabilize.

For the 28 choices \(|\delta|\in\{1,3\}\),
\(\max(1,\delta+1)\le i\le\min(9,9+\delta)\), direct rational
substitution gives a nonzero numerator of absolute value at most
\(418112777\). For a finite-order example, take
\[
 n=146603499986155091679122137738,
\]
which has 97 bits and lies above the current cutoff in (T5). Its residues modulo \(79^5,83^5,97^6\) are
\(625251304,3566867498,579419484714\), respectively. The selected
values of \(i\) are \(8,9,5\). At each prime the proof of the
local edges applies directly, with \(19(n+1)<p^{18}\) in the
large-index bound. The exact residual coefficients and their
first ramified lifting values are as follows.

| \(p\) | Lower residual | Upper residual | Shared \(Y\) | \(\Delta_{\delta,i}\bmod p\) |
|---:|:---|:---|---:|---:|
| 79 | \(1+78Y^4\) | \(1+58Y+20Y^5\) | 1 | 0 |
| 83 | \(1+76Y^4\) | \(1+67Y+Y^5\) | 21 | 17 |
| 97 | \(1+77Y^5\) | \(1+31Y+51Y^6\) | 51 | 0 |

The residual test allows degree two simultaneously at the three
primes. The first ramified lifting condition at \(83\) excludes
the shared slope there.
