# Stationary central roots and the compactification of weak-tail jets

Fix ordinary two-monomial core moments \(Q_n,Q_N\) with
\(N=n+1\), \(n\ge2\). Append a nonempty finite shared tail
with signed gaps \(d_i\ne0\) and midpoints \(M_i\).
Use the weak curvature \(W=\sum_i d_i^2/M_i\) and put
\[
 w_i=\frac{d_i^2}{M_i W},\qquad u_i=\frac1{d_i},\qquad
 x=\sum_i w_i u_i,\quad \mu=\sum_i w_i^2 u_i,\quad
 \nu=\sum_i w_i^2u_i^2,\quad \xi=\sum_i w_i^3u_i^2 .
                                                               \tag{1}
\]
The signed gap \(h=\sum_i d_i\) can be fixed at any integer.
The results allow arbitrary row counts and independently
varying heights and gaps.

For a core root \(\alpha\) of odd multiplicity \(r=2l+1\)
at order \(m\), use the expansion
\[
 Q_m(\alpha e^y)=b y^r(1+c_1y+c_2y^2+c_3y^3+O(y^4)),
                      \qquad b\ne0 .
\]
The coefficients are exact algebraic numbers for ordinary
core moments. Normalize the tail endpoints by a positive
variable dilation and scalar at each order.

**Theorem.** The central root in this cluster, defined by
the [weak-tail Hermite expansion](MULTIVARIATE_WEAK_TAIL_MULTIPLE_ROOT_RESPONSE.md),
has angular displacement
\[
 -\frac{W}{m}\operatorname{Im}c_1+
 W^2\mathcal A_{m,r}(\mu,\nu,\xi)+O_{Q_m,m}(W^3),
                                                               \tag{2}
\]
uniformly in the tail, where
\[
\begin{aligned}
 \mathcal A_{m,r}(\mu,\nu,\xi)
   =\operatorname{Im}\Bigg[&
     -\frac{r+2}{m^2}\left(c_3-c_1c_2+\frac{c_1^3}{3}\right)
     -\frac{\mu}{3m^2}\left(lc_1^2-(r+2)c_2\right)\\
     &+\left(\frac{\nu}{2m^2}-\frac{\mu}{2m}
                        +\frac{l(\mu^2-2\xi)}{3m^2}\right)c_1
                           \Bigg].
\end{aligned}                                                    \tag{3}
\]
Thus a pair with equal first angular velocities is separated
by every sufficiently weak tail if the difference of (3) is
nonzero throughout the compact set of its limiting statistics
described below. This extends the simple-root second-order
test to any odd multiplicities.

At a physical integer-twist contact one only needs this
nonvanishing on the intersection of that set with
\[
 2[N c_{1,n}-n c_{1,N}]+x+
                   \frac{\mu}{3}[n(r_N-1)-N(r_n-1)]=0.           \tag{4}
\]
Equation (4) is the exact first physical condition of the
[multiple-root comparison](MULTIVARIATE_WEAK_TAIL_MULTIPLE_ROOT_RESPONSE.md#an-additional-physical-test-at-a-stationary-odd-pair).
All the other branches are already separated there.
The compact test is effective whenever its minimum absolute
margin is positive. Below the absolute cutoff it can be
combined with the earlier even-multiplicity and first-order
tests to certify all degrees for the original weak extensions
of a fixed core.

## The compact set for weak-tail statistics of any fixed jet order

Write \(\mathcal U=\{0\}\cup\{\pm1/d:d\in\mathbb Z_{>0}\}\).
For any fixed finite collection of integer pairs \(p\ge2,q\ge0\), set
\[
                   S_{p,q}=\sum_i w_i^p u_i^q.
\]
The exact limiting set of \((x,(S_{p,q}))\) as \(W\downarrow0\)
at fixed signed gap \(h\) is the following. Take any finite
or countable list of \(w_i>0\), \(u_i\in\mathcal U\) with
\(\sum_i w_i\le1\), and one real number \(s\) with
\[
             |s|\le 1-\sum_i w_i,\qquad
 x=\sum_i w_i u_i+s,\qquad S_{p,q}=\sum_i w_i^p u_i^q .           \tag{5}
\]
Use \(u^0=1\). The number \(s\) describes the contribution of
weights tending individually to zero. It is the only residual
mass term needed for the factorial jets below.

To prove necessity, order each tail by decreasing weights
and pad it by zeros. Take a diagonal subsequence on which
each weight and label converges and on which \(x\) converges.
The weight sum of every finite prefix is at most one. The
tail after \(K\) entries contributes at most
\((K+1)^{1-p}\) in absolute value to \(S_{p,q}\), since
\(|u_i|\le1\) and \(w_i\le1/(K+1)\) for \(i>K\).
This gives the formula for \(S_{p,q}\) in the limit.
The remaining contribution to \(x\) has absolute value at
most the remaining mass. Letting \(K\) increase gives (5).

Conversely, retain a finite initial part of the weight list.
Split the unused mass into small weights with labels \(+1,-1\)
to approximate the desired total contribution to \(x\), with
their contribution to every \(S_{p,q}\) tending to zero.
Replace occurrences of label zero among the retained weights
by reciprocals of arbitrarily large integer gaps. For each
finite list the resulting signed gap is an integer, so correct
it to \(h\) by adding unit-gap rows of arbitrarily small total
weight and renormalizing. For any finite selection of positive
weights and nonzero reciprocal labels, choose the corresponding
signed gaps and set the midpoint of row \(i\) to a nearest
allowed midpoint to \(d_i^2/(t w_i)\). As \(t\downarrow0\)
the actual normalized weights converge to the chosen weights
and \(W\sim t\). This proves the claimed realizability by
a diagonal sequence of ordinary shared tails.

The set (5) admits explicit finite algebraic outer
approximations. Retain at most \(K\) weights and allow only
labels \(0,\pm1,\ldots,\pm1/D\). Replacing the omitted small
reciprocal labels by zero changes \(x\) by at most \(1/D\)
and changes every \(S_{p,q}\), \(q\ge1\), by at most \(1/D\).
Allow additionally \((K+1)^{1-p}\) error in \(S_{p,q}\) for
the weights after the first \(K\). For the first statistic
retain the interval for \(s\) in (5) for the entire remaining
mass. The center of each such error box is itself an
attained limit, because finite weighted lists are a case
of (5). The outer sets converge to the exact limiting
set with the explicit stated coordinate errors. In such
a finite test the weights and \(s\) are real variables on
a simplex and an interval. The labels have finitely
many choices. Real algebraic feasibility can therefore
certify disjointness of equations such as (3) and (4).
If those equations have no simultaneous solution on (5),
a sufficiently fine outer test succeeds by compactness.

For the second jet, \(S_{2,1}=\mu\), \(S_{2,2}=\nu\), and
\(S_{3,2}=\xi\). There is also a useful inexpensive outer
test using \(\lambda=2\xi-\mu^2\):
\[
                \mu^2\le\lambda,\qquad \nu^2\le\lambda
                         \le\nu\le1 .                           \tag{6}
\]
Indeed Cauchy–Schwarz gives \(\mu^2\le\xi\) and
\(\nu^2\le\xi\). For the upper bound let \(y_i=w_i u_i\).
Then \(\nu-\lambda=\sum_i(1-2w_i)y_i^2+(\sum_iw_i y_i)^2\).
This is
\[
 \sum_{i<j}(w_j y_i+w_i y_j)^2
       +2\sum_i y_i^2\!\!\sum_{\substack{j<k\\j,k\ne i}} w_jw_k
\]
when \(\sum_iw_i=1\), and is nonnegative. The bounds persist
in (5) by taking limits. Equations (3) are affine in
\((\mu,\nu,\lambda)\). For two stationary branches the terms
\(-\mu\,\operatorname{Im}(c_1)/(2m)\) cancel between them.
Equations (4), (6), and the first-order joint constraints
on \((x,\mu)\) can therefore give a simple exact certificate
before resorting to finer outer approximations.

## Factorial expansion and the second central coefficient

Let \(G_m(k)=\sum_{j\ge1}g_j k^j\) interpolate the log
endpoint-normalized tail multiplier. For this proof take
the polynomial on a fixed integer grid containing
\(0,\ldots,m\) and at least the indices \(0,\ldots,6\).
Extend the log-Gamma chords to that grid before interpolating.
For sufficiently small \(W\) all factorial arguments in
that fixed real range are positive. This allows the same
Taylor coefficients at low orders even when \(m<4\).
On \(Q_m\) only the values at \(0,\ldots,m\) are used.

The required coefficients are
\[
 \begin{aligned}
    g_2&=\frac{W}{2m}
              +W^2\left(\frac{\mu}{4m}-\frac{\nu}{4m^2}\right)
                                                      +O_m(W^3),\\
    g_3&=-W^2\frac{\mu}{6m^2}+O_m(W^3),\\
    g_4&= W^3\frac{\xi}{12m^3}+O_m(W^4),\qquad
                          g_j=O_m(W^{j-1})\quad(j\ge2).
 \end{aligned}                                                   \tag{7}
\]
Here is the general expansion underlying (7).
Write \(s=k/m-1/2\), and use the Bernoulli numbers with
\(B_1=-1/2\). The Taylor coefficients of \(q(t)=t/(e^t-1)\)
in the integral for the trigamma give, through any fixed
order,
\[
 \frac{d^2}{dk^2}\log T_{m,k}
 = \sum_{\substack{a,b\ge0\\a+b<P}}
     \frac{(-1)^a B_b\binom{a+b}{a}}{m^{b+1}}\,
          W^{a+b+1} s^a S_{a+b+1,a+2b}
                                      +O_{m,P}(W^{P+1}),        \tag{8}
\]
where \(S_{1,0}=1\). Indeed \(M_i^{-1}=W w_i u_i^2\)
and \(d_i/M_i=W w_i u_i\). For the argument
\(x_i=m M_i(1+Ww_i u_i s)\), use the expansions of
\[
 d_i^2\psi'(1+x_i)
     =d_i^2\int_0^\infty e^{-x_i t}q(t)\,dt .
\]
The Taylor error on \(t\le1\) and the exponentially
small integral remainder on \(t\ge1\) give the uniform
trigamma expansion with computable constants. Expansion
of \((1+Ww_i u_i s)^{-b-1}\) then gives (8).
The per-row errors are bounded by a constant times
\((Ww_i)^{P+1}\), whose sum is at most a constant times
\(W^{P+1}\). Integrating twice with zero endpoint values
gives the log multiplier jets. The coefficient of \(k^2\)
through order two, of \(k^3\) through order two, and of
\(k^4\) through order three gives (7) after interpolation
of the displayed polynomials.

Similarly, the actual endpoint log scale correction after
subtracting \(h\log(N/n)\) has finite jets in these statistics.
It comes from integrating \(x\psi'(1+jx)-1/j\) over \(n\le j\le N\)
and over each oriented row interval
\(M_i+d_i[-1/2,1/2]\).
The same expansion now starts at \(b\ge1\), and its
\(a,b\) term uses \(W^{a+b}S_{a+b,a+2b-1}\) times a
fixed integral in \(j\) and the normalized row position.
The only nonconstant statistic of weight exponent one
is \(S_{1,1}=x\). This is why a single signed residual
mass parameter suffices for the full finite jets.

For the angular central calculation the real drift
\(g_1\) is an exact translation in logarithmic coordinate.
Remove it. Set \(t=\sqrt W\), and treat
\[
           a=g_2/W,\quad M=g_3/W^2,\quad
                      T=g_4/W^3,\quad U=g_5/W^4
\]
as fixed real parameters at a particular tail. They lie
in a common compact set and \(a\to1/(2m)\).
The scaled operator on \(u=y/t\) starts with
\(a\partial_u^2+t M\partial_u^3+t^2T\partial_u^4+
t^3U\partial_u^5\).
Let \(H_j(u)=\exp(a\partial_u^2)u^j\).
For the odd central branch the scaled solution is
\(u=t A+t^3D+O(t^5)\), with uniform error.
To justify the parity and uniformity, freeze all real
parameters \(g_j/W^{j-1}\), replace \(t\) by a free
complex variable, and divide the scaled expression by
\(t^r\). Its singularity at zero is removable by the
same order count as in the Hermite expansion. Replacing
\(t,u\) by \(-t,-u\) multiplies the expression by \((-1)^r\).
The central zero is uniformly simple at \(t=0\), so its
analytic root is odd in \(t\), with a common local Taylor
bound over the compact parameter range. Evaluating this
analytic expansion at the actual \(t\) is valid even
when the tail parameters vary with \(W\).

For completeness the coefficient of \(t\) gives
\[
                        A=-2a c_1-Ml/a .
\]
The coefficient of \(t^3\), divided by \(H_r'(0)\),
is \(D\) plus the sum of the following contributions,
where every derivative is evaluated at zero:
\[
\begin{aligned}
 &\frac{A^3}{6}\frac{H_r'''}{H_r'}
 +\frac{A^2}{2}\frac{c_1 H_{r+1}''+M H_r^{(5)}}{H_r'}\\
 &+\frac{A}{H_r'}\left(c_2 H_{r+2}'+c_1 M H_{r+1}^{(4)}
                        +T H_r^{(5)}+\frac{M^2}{2}H_r^{(7)}\right)\\
 &+\frac1{H_r'}\left(c_3 H_{r+3}+c_2 M H_{r+2}'''
          +c_1 T H_{r+1}^{(4)}+\frac{c_1 M^2}{2}H_{r+1}^{(6)}
          + U H_r^{(5)}+ MT H_r^{(7)}+\frac{M^3}{6}H_r^{(9)}\right).
\end{aligned}
\]
Use \(H_j^{(k)}(0)=j!a^{(j-k)/2}/((j-k)/2)!\)
when \(j-k\) is even and nonnegative, and zero otherwise.
Substitution of \(A\) and collection of the nonreal terms
give
\[
 \operatorname{Im}D=
 \operatorname{Im}\left[
 -4a^2(r+2)(c_3-c_1c_2+c_1^3/3)
 +2M(lc_1^2-(r+2)c_2)+l(3M^2/a^2-4T/a)c_1\right].
\]
The angular log displacement is
\(W\operatorname{Im}A+W^2\operatorname{Im}D+O(W^3)\).
In its first term retain the order-\(W\) correction to
\(a\) in (7); in its second use the leading values of
\(a,M,T\) there. This proves (2) and (3).

The same uniform analytic expansion and (8) produce any
fixed finite order of the root jet, with coefficients
polynomial in finitely many of the statistics (5).
When using these higher coefficients it is necessary
to keep track of which lower coefficients vanish
identically: a lower coefficient depending on the tail
can itself approach zero at a compensating rate. The tests
(2) and (4) avoid this issue because stationarity makes
the first angular coefficient identically zero in the
difference of the two branches.

## A stationary control that passes the first physical test

Take \(H(z)=z^2+2z+17/16\), \(\alpha=-1+i/4\) and
\[
 V(z)=z^2+\frac{66}{31}z+\frac{119}{496},\qquad
 Q_{12}=H^3(1+A H^3),\qquad Q_{13}=H^3(V+z H^3).
\]
All coefficients are positive for \(A>0\), and \(Q_{12}\)
is positive on the real line. The common conjugate pair
has multiplicity three in each parent. Its two normalized
central velocities satisfy \(c_{1,12}/12=c_{1,13}/13\),
because \(\alpha V'(\alpha)/V(\alpha)=c_{1,12}/12\).
Thus both first-order tests allow cancellation, for
example at \(x=\mu=0\) in (4).
Increasing the positive rational number \(A\) changes
neither \(c_{1,12}\) nor \(c_{2,12}\), while it changes
\(c_{3,12}\) by \(A[\alpha(\alpha-\overline\alpha)]^3\).
The imaginary part of the latter factor is nonzero.
Since the other terms in (3) are uniformly bounded
on (6), the second angular difference has one strict
sign for all large enough \(A\). For example \(A=1000\)
gives the exact upper-minus-lower score
\[
 \frac{220696375}{63857664}+\frac{337}{55432}\mu
                     -\frac{\nu}{624}+\frac{\lambda}{936},
\]
which is larger than three on (6). This controls the
criterion at exact first-order stationarity. The two
control polynomials serve as an algebraic fixture.
