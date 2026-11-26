# A signed endpoint test for diffuse factorial tails

A fully shared tail attached to an ordinary mirrored core can be
certified in every degree by comparing two quantities from its
original integer rows. A single sign condition on a polynomial
of the signed endpoints controls the exact adjacent radial scale.
Signed lower endpoints grouped by gap control the nonsymmetric
part of the factorial coefficients. The test permits arbitrarily
many rows, different gaps, and opposite orientation orders in
different gap groups.

## Criterion

Remove equal rows. Fix any ordinary mirrored core, including
private mirrored pairs. The added rows are fully shared. Opposite
copies of an exact reflected pair can be moved into the core. Write
their lower and upper exponents as \(a_i,b_i=a_i+d_i\), with
positive integers \(a_i,d_i\). Let \(\epsilon_i=1\) if the row
increases from the first monomial to the second and \(-1\)
otherwise. Assume the total signed gap is zero. Define
\[
 \begin{split}
 P_d(X)&=\sum_{d_i=d}\epsilon_i X^{a_i},&
 T(X)&=\sum_d(1-X^d)P_d(X),\\
 A&=\sum_d d^3\int_0^\infty u\,|P_d(e^{-u})|\,du,&
 D&=\left|\sum_i\epsilon_i
               \log\frac{b_i+1/2}{a_i+1/2}\right|.
 \end{split}
\]
Take \(J,B,\eta_1\) from the
[uniform exact skew test](MULTIVARIATE_ORDERED_OPPOSITE_PAIR_RADIAL_CLOSURE.md#an-exact-skew-test-for-any-diffuse-tail-of-a-mirrored-core).
In particular \(J\ge4\) is a computable absolute upper cutoff
for the remaining terminals, and \(B,\eta_1>0\) are uniform in
the core and the tail. These constants may be unevaluated.

**Theorem.** Suppose \(T\) is nonzero and has one sign on \((0,1)\).
If the tail satisfies
\[
             \max_i\frac{d_i}{a_i+d_i/2}\le\eta_1,\qquad
                        A\le\frac{D}{B J^3},
                                                               \tag{1}
\]
then its extension of the core satisfies ordinary SFC(2) in
every degree over every field of characteristic zero. Both
conditions concern the physical coefficient ratio, with no
independent rescaling of that ratio between the two orders.

The sign of \(T\) can be checked by univariate real root
isolation. The number \(D\) is the absolute logarithm of a
positive rational number. For \(A\), split each integral at the
positive \(u\) corresponding to roots of \(P_d\) in \((0,1)\)
and integrate each monomial \(u e^{-a_i u}\) on the resulting
intervals. Its antiderivative is
\(-e^{-a_i u}(u/a_i+1/a_i^2)\). In particular \(A\) is rational
if each \(P_d\) has one sign, even if different \(P_d\) have
opposite signs. Strict instances of the numerical inequality
in (1) can be certified with rational enclosures. If the
tail is replicated with any positive integer multiplicity,
both \(A\) and \(D\) acquire that same factor.

## The physical radial scale

Write \(g(x)=\log\Gamma(1+x)\) and \(q(u)=u/(e^u-1)\), with
\(q(0)=1\). For adjacent \(n,N=n+1\) put
\[
 \Phi(x)=g'(Nx)-g'(nx)-\log(N/n),\qquad
 K(x)=N g''(Nx)-n g''(nx).
\]
The trigamma series gives
\[
 K(x)=\int_0^\infty e^{-xu}
             [q(u/N)-q(u/n)]\,du .
\]
Indeed \(g''(x)=\sum_{h\ge1}(x+h)^{-2}\), and integrating
\(u e^{-(x+h)u}\) and summing the geometric progression gives
the Laplace form. All arguments are positive.
[Trigamma series, DLMF 5.15.1](https://dlmf.nist.gov/5.15.E1)

The difference defining \(\Phi\) tends to zero at infinity.
To see this directly, use
\(1/(1+y)\le g''(y)\le1/y\) for \(y>0\), by integral
comparison in the series, and integrate \(xg''(t x)\) over
\(n\le t\le N\). Since \(\Phi'=K\), it follows that
\[
          \Phi(x)=-\int_0^\infty e^{-xu}
                  [q(u/N)-q(u/n)]\,\frac{du}{u}.
\]
For an increasing row the difference of its two log endpoint
scales is \(\int_{a_i}^{b_i} (g'(Nx)-g'(nx))\,dx\).
Reversal negates this expression. Upon summing, the
\(\log(N/n)\) terms cancel by the zero signed gap, and the
mirrored core contributes zero. Thus the exact scale difference
\(C=\log(G_N/G_n)\) from that test is
\[
                C=-\int_0^\infty [q(u/N)-q(u/n)]
                                     T(e^{-u})\,\frac{du}{u^2}.
                                                               \tag{2}
\]
The integrals converge, since the bracket is \(O(u)\) near zero
and \(T(1)=0\), and all \(a_i\ge1\) at infinity.

For real \(v\ge0\), the elementary bound
\(-q'(v)\ge e^{-v}/2\) follows for \(v>0\) by multiplying out.
It becomes \(2v-3+4e^{-v}-e^{-2v}\ge0\), whose left side
vanishes at zero and has derivative \(2(1-e^{-v})^2\).
Continuity gives the bound at zero. Therefore
\[
             q(u/N)-q(u/n)\ge
                         \frac{u}{2nN}\,e^{-u/n}.
\]
Under the one-sign condition on \(T\), equation (2) and \(n\ge2\)
give
\[
 |C|\ge\frac1{2nN}\int_0^\infty e^{-u/2}|T(e^{-u})|\,\frac{du}{u}
       =\frac{D}{2nN}>0.                                         \tag{3}
\]
For the equality, integrate
\((e^{-(a_i+1/2)u}-e^{-(b_i+1/2)u})/u
=\int_{a_i+1/2}^{b_i+1/2} e^{-x u}\,dx\) in \(u\).
The constant sign makes the integral of the absolute value
equal the absolute value of the summed integral. Its positivity
comes from the nonzero polynomial \(T\).

## The coefficient asymmetry

At order \(j\), continue the endpoint-normalized, binomial-stripped
log coefficient of each increasing row as a function of \(v=k/j\).
The row contributes
\[
 g(j(a_i+d_i v))-(1-v)g(ja_i)-v g(jb_i).
\]
Let \(J_j^o(v)\) be the odd part of the sum of these row logarithms,
with the reversed rows using \(1-v\). Oddness here means taking
half the difference under \(v\mapsto1-v\). The mirrored core
does not enter. In the Laplace integral for the second derivative
of \(J_j^o\), put \(y=v-1/2\) and use the variable \(u\) scaled
by \(j\). It gives exactly
\[
 (J_j^o)''(v)=-j\sum_d d^2\int_0^\infty
    e^{-du/2}\sinh(dyu)\,q(u/j)P_d(e^{-u})\,du.
\]
All derivatives can pass inside the integral because the row
arguments are at least \(j\min_i a_i\) on the coefficient interval.
Since \(0<q\le1\) and \(|y|\le1/2\),
\[
        |e^{-du/2}\sinh(dyu)|\le du/2,\qquad
                            \|(J_j^o)''\|_\infty\le jA/2.
\]
The first bound uses \(\sinh(t)\le t e^t\) for \(t\ge0\).
The odd chord has zero endpoint values. The interval Dirichlet
Green kernel for the second derivative has absolute integral
\(v(1-v)/2\le1/8\), so the skew norm in the cited radial test
satisfies
\[
                                \delta_j\le jA/16.              \tag{4}
\]

At a possible collision with \(3\le N<J\), that test requires
\(|C|\le2B(\delta_n+\delta_N)\). The curvature and relative-gap
reductions in the test apply under (1). But (1) and (4) give
\[
               2B(\delta_n+\delta_N)\le B J A/4
                         \le D/(4J^2)< |C|,
\]
by (3). The same cited reduction covers all \(N\ge J\), and
terminal two is separated by positivity of the even moment at
every real coefficient ratio. This proves the theorem over
\(\mathbb C\). The coefficients are rational, so the rational
Bézout identity gives the assertion over every
characteristic-zero field.

## An opposed-order ordinary block

Take three copies of each of \((H,H+1)\) and \((H+2,H+1)\),
and one copy each of \((H+3,H+1)\) and \((H+2,H+4)\), with
\(H\) a positive integer. Then
\[
 P_1=3X^H(1-X),\quad P_2=-X^{H+1}(1-X),\quad
 T=X^H(1-X)^2(3-X-X^2)>0\quad(0<X<1).
\]
Both gap groups balance. An ordered matching of the gap units
with \(d_+\ge d_-\) in every match would have to match equal gaps,
since the positive and negative gap-unit measures have equal
first moments in \(d\). The same holds for \(d_+\le d_-\).
In the gap-one group the increasing midpoint is smaller,
whereas in the gap-two group it is larger. Hence the block
has no matching with a single midpoint ordering as in the
[ordered gap-mass criterion](MULTIVARIATE_ORDERED_GAP_MASS_RADIAL_CLOSURE.md).

Integration gives
\[
 \begin{split}
 A&=3\left(\frac1{H^2}-\frac1{(H+1)^2}\right)
     +8\left(\frac1{(H+1)^2}-\frac1{(H+2)^2}\right)
                      \le \frac{22}{H^3},\\
 D&\ge\int_0^\infty e^{-(H+1/2)u}(1-e^{-u})^2\,\frac{du}{u}
                      \ge \frac1{(H+5/2)^2}.
 \end{split}
\]
The first inequality uses the derivative of \(x^{-2}\).
The last one uses \(1-e^{-u}\ge u e^{-u}\). This gives an
effective cutoff on \(H\) for (1), valid for any multiplicity
and any mirrored core.

## Higher signed midpoint cancellations

There is a sufficient comparison in terms of the signed lower
heights within each fixed gap. For each nonzero group \(P_d\)
let \(\nu_d=\sum_{d_i=d}\epsilon_i\delta_{a_i}\), and let \(a_{0,d}\)
be its smallest supported lower height after moving exact opposite
copies into the core. Choose an integer
\(p_d\ge1\), which may vary by gap, and suppose \(\nu_d\)
annihilates all polynomials of degree less than \(p_d\).
Assume also that, simultaneously for all groups, the signed
functions
\[
 H_d(t)=(-1)^{p_d}\int
             \frac{(a-t)_+^{p_d-1}}{(p_d-1)!}\,d\nu_d(a)
\]
are nonnegative almost everywhere. Reversing every sign is
also allowed. The truncated zeroth power is the indicator of
\(a>t\).

Taylor's formula with integral remainder, with basepoint
\(a_{0,d}\), gives
\[
 P_d(e^{-u})=u^{p_d}\int_{a_{0,d}}^\infty e^{-tu}H_d(t)\,dt.
\]
The compact support of \(H_d\) follows from the moment
conditions and the finite support of \(\nu_d\).
Every nonzero group is therefore positive on \(0<X<1\).
The group contributions \(A_d,D_d\) to \(A,D\) are
nonnegative, without intergroup cancellation. Using
\(1-e^{-du}\ge d u e^{-du}\), integrate the nonnegative
kernels in \(u\) first to obtain
\[
 \begin{split}
 A_d&=d^3(p_d+1)!\int H_d(t)t^{-p_d-2}\,dt,\\
 D_d&\ge d\,p_d!\int H_d(t)(t+d+1/2)^{-p_d-1}\,dt .
 \end{split}
\]
Consequently,
\[
        \frac{A_d}{D_d}
          \le \frac{d^2(p_d+1)}{a_{0,d}}
            \left(1+\frac{d+1/2}{a_{0,d}}\right)^{p_d+1}.
                                                               \tag{5}
\]
The bound is the ratio of the positive integrands at the smallest
height, since that ratio decreases with \(t>0\). The maximum of
(5) bounds \(A/D\). In particular all such profiles satisfy (1)
if
\[
       p_d+1\le\varepsilon a_{0,d}/d^2\quad\hbox{for every }d,
 \qquad
       0<\varepsilon\le\min\{\eta_1,\,1/2,\,1/(3BJ^3)\}.
\]
Indeed the exponential factor in (5) is at most
\(e^{3\varepsilon/2}<3\), and \(d/(a_{0,d}+d/2)\le\varepsilon\).

For explicit integer rows at gap \(d\), take the signed
multiplicities \((-1)^k\binom pk\) at lower heights
\(a+kh\), \(0\le k\le p\), with \(a,h\ge1\) integers. The
corresponding transform is
\[
       P_d(e^{-u})=e^{-a u}(1-e^{-h u})^p
        =u^p\int_{[0,h]^p}e^{-(a+v_1+\cdots+v_p)u}\,dv .
\]
This expression gives the same comparison leading to (5)
directly, with \(t\ge a\). The threshold therefore allows the
number \(p\) of canceled moments to grow linearly with \(a/d^2\).
These cancellations are in heights at a fixed gap, rather than
in gaps at a fixed normalized midpoint.
