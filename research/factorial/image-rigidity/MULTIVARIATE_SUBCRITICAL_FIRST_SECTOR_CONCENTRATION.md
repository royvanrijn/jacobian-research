# First-sector bounds at and below shared critical curvature

For ordinary exponent columns \(a,b\), let
\[
 P_n(z)=\sum_{k=0}^n p_{n,k}z^k,\qquad
 p_{n,k}=\binom nk\prod_i((n-k)a_i+kb_i)! .
\]
Rows with equal exponents can be omitted. For a fully shared profile
write \(d_i=b_i-a_i\) and define
\[
 W(t)=\sum_i\frac{d_i^2}{a_i+d_i t},\qquad
 \rho(t)=t(1-t)W(t),\qquad
 \gamma=\max_{0\le t\le1}\rho(t).
\]

**Theorem.** Suppose the profile is fully shared and \(\gamma\le1\).
For any \(n\ge1\) and \(r>0\), give \(K\in\{0,\ldots,n\}\) mass
proportional to \(p_{n,k}r^k\). Uniformly in \(r\) one has
\[
 \operatorname{Var}(K)\le34n^{3/2},\qquad
 \operatorname{Var}(K)\le\frac{16n}{1-\gamma}\quad(\gamma<1).
                                                               \tag{1}
\]
Consequently \(P_n\) has no zero with
\[
 |\arg z|\le\frac1{\sqrt{34}\,n^{3/4}},
 \qquad\text{or, if }\gamma<1,\qquad
 |\arg z|\le\frac{\sqrt{1-\gamma}}{4\sqrt n}.
\]
For every \(n\ge2048\) it has no zero in the entire first sector
\(|\arg z|\le2\pi/(n+1)\). These statements permit balanced and
mixed orientations at the critical value \(\gamma=1\).

There is also an exact smaller-order criterion. For fully shared
profiles let \(V=\sum_i d_i^2/\min(a_i,b_i)\).
If \(n\ge2\) and \(V\le 4n/(n+1)\), the order-\(n\) polynomial
has no root in the closed first sector. In particular,
\(V\le 26/7\) excludes consecutive common roots in the first sector
at every order, using the existing ordinary terminals through
thirteen for \(n\le12\).
If there is exactly one private unit row, on the \(a\) side, and
all other rows are shared, then the sufficient condition
\(\sum_{\rm shared}d_i^2/b_i\le1\) excludes the first sector
for every \(n\ge2\) by the same log-concavity argument.

## Critical curvature and a continuous coefficient bound

For an unequal row, the function
\(\rho_i(t)=d_i^2t(1-t)/(a_i+d_i t)\) has maximum
\((\sqrt{a_i}-\sqrt{b_i})^2\), attained at
\(\sqrt{a_i}/(\sqrt{a_i}+\sqrt{b_i})\), and is strictly concave.
Thus \(\gamma\le1\) and positive integer endpoints imply
\(1/4\le b_i/a_i\le4\). A maximizer \(t_*\) of \(\rho\)
lies in \([1/3,2/3]\), since every individual maximum lies
in that interval. If every gap is zero, take \(t_*=1/2\).

The exact tangent remainder for each rational row term is
\[
 \rho_i(t_*)+\rho_i'(t_*)(t-t_*)-\rho_i(t)
  =\frac{d_i^2 a_i b_i (t-t_*)^2}
           {(a_i+d_i t_*)^2(a_i+d_i t)} .
\]
For the indicated ranges,
\(a_i b_i/(a_i+d_i t_*)^2\ge4/9\) and
\(a_i+d_i t\le2(a_i+d_i t_*)\). For example, the square
of \((a_i+d_i t_*)/\sqrt{a_i b_i}\) is at most the square
of \((1/3)/2+(2/3)2=3/2\). The second comparison follows
since each endpoint is at most twice a convex combination
with both weights at least \(1/3\) and endpoint ratio at most
four. Hence \(W(t)\ge W(t_*)/2\ge2\gamma\).
The tangent terms cancel in the sum, giving
\[
          1-\rho(t)\ge (1-\gamma)+\tfrac89\gamma(t-t_*)^2
                  \ge\tfrac89(t-t_*)^2 .                      \tag{2}
\]

Extend the logarithm of \(p_{n,k}r^k\) to \(F(x)\) on \([0,n]\)
by Gamma functions, and put \(t=x/n\). The trigamma series
gives
\[
 \psi_1(s+1)<\frac1{s+1/2},\quad
 \psi_1(s+1)\ge\frac1{s+1},\quad
 \psi_1(s+1)\ge\frac1s-\frac1{2s^2}\quad(s>0).
\]
The first two comparisons use midpoint and right-endpoint
integral sums of \(u^{-2}\). For the last, the convex trapezoidal
sum gives \(1/(s+1)+1/[2(s+1)^2]\), which is at least the
displayed bound. It follows that, for \(n\ge64\),
\[
 F''(x)\le \frac{33}{n^2}-\frac{32}{9}\frac{(x-nt_*)^2}{n^3},
 \qquad
 F''(x)\le \frac{33}{n^2}-\frac{4(1-\gamma)}{n}.                \tag{3}
\]
To verify this first on \(1/8\le t\le7/8\), the row curvature
is at most \(W(t)/n\); the negative Gamma curvature from the
binomial is at least \(1/[nt(1-t)]-33/n^2\), since
\(\tfrac12[t^{-2}+(1-t)^{-2}]\le32+32/49<33\).
Use (2) and \(t(1-t)\le1/4\).

For \(t\le1/8\), the ratio of the row denominators at \(1/3\)
and \(t\) is at most \(2/(1+3t)\). This follows by maximizing
\((1+s/3)/(1+s t)\) over \(-3/4\le s\le3\). Therefore
\(W(t)\le9\gamma/(1+3t)\). The binomial curvature times \(n\)
is at least \(1/(t+1/64)+64/65\). The difference
\[
              \frac1{t+1/64}-\frac9{1+3t}
\]
decreases on \([0,1/8]\), as
\((1+3t)/(t+1/64)\ge88/9>3\sqrt3\). On \([0,1/10]\), its
value plus \(64/65\) is at least
\(320/37-90/13+64/65>128/81\), a bound for
\((32/9)(t-t_*)^2\). On \([1/10,1/8]\), its value plus
\(64/65\) is at least \(56/99+64/65>(32/9)(17/30)^2\),
again bounding \((32/9)(t-t_*)^2\).
This proves the first part of (3) outside the central interval.
For the second part, the same difference with \(9\gamma\) in
place of \(9\), including \(64/65\), is at least \(4(1-\gamma)\)
by linearity in \(\gamma\). Check \(\gamma=1\) from the
positive margin just proved and \(\gamma=0\) from
\(64/9+64/65>4\). Reverse \(a,b\) for \(t\ge7/8\).

## Compare two gaps with the same midpoint

Let \(K,L\) be independent with the tilted law and put
\(c=(K+L)/2\), \(d=|K-L|/2\). Conditional on \(c\), the allowable
\(d\)'s lie on the integer or half-integer grid of spacing one.
Their weights are proportional to
\(\exp(S_c(d))\), with a factor of two for \(d>0\), where
\(S_c(y)=F(c+y)+F(c-y)\). For every allowable \(d\ge1\) let
\(d'\le d/2\) be the largest point on the same grid at most
\(d/2\). This is allowable, and each \(d'\) has at most two
preimages. From the first inequality of (3),
\[
 S_c'(y)=\int_{c-y}^{c+y}F''(x)\,dx
        \le \frac{66y}{n^2}-\frac{64y^3}{27n^3}.
\]
Here the integral of \((x-nt_*)^2\) is at least \(2y^3/3\).
Integration from \(d'\) to \(d\) gives
\[
          S_c(d)-S_c(d')\le
                         \frac{33d^2}{n^2}-\frac{5d^4}{9n^3}.
                                                               \tag{4}
\]
A uniform penalty \(-E\) in (4) over a tail of \(d\) bounds
its probability by \(4e^{-E}\), since multiplicity and the
orientation factors each cost at most two.

For \(n\ge64\), put \(s=d/n^{3/4}\). At \(s\ge4\) the right
side of (4) is at most \(-5s^4/18\), since
\(33/\sqrt n\le33/8\). Integration of the tail gives
\[
 \frac{\mathbb E[d^2]}{n^{3/2}}
       \le16+8\int_4^\infty s e^{-5s^4/18}\,ds<17.
\]
For the last estimate, \(s^4\ge16s^2\) on the interval, so
the integral times eight is at most
\((9/10)e^{-640/9}<1\). Since
\(\operatorname{Var}(K)=2\mathbb E[d^2]\), this proves the
critical bound in (1). For \(n<64\) the trivial variance
bound \(n^2/4\) suffices.

If \(\epsilon=1-\gamma>0\), \(n\ge64\) and
\(n\ge33/(2\epsilon)\), the second inequality of (3) gives
\(F''\le-2\epsilon/n\), so
\[
                    S_c(d)-S_c(d')
                          \le-\frac{3\epsilon d^2}{2n}.
\]
The same conditional comparison bounds
\(\Pr(d\ge x)\le\min(1,4e^{-3\epsilon x^2/(2n)})\).
For \(x<1\) the right side is one, so the grid comparison
is needed only for \(d\ge1\). Integrating gives
\(\operatorname{Var}(K)\le4(1+\log4)n/(3\epsilon)\).
If \(n\) fails one of the two lower bounds, then
\(n^2/4\le16n/\epsilon\). This proves the second part of (1).

The angular bounds in the theorem follow from
\[
 \Re\mathbb E e^{i\theta(K-\mathbb E K)}
                  \ge1-\frac{\theta^2}{2}\operatorname{Var}(K),
 \qquad
 \frac{P_n(re^{i\theta})}{P_n(r)}=\mathbb E e^{i\theta K}.
\]

## A sharper first-sector threshold

The differential inequality for \(S_c'\) also shows that
\[
           S_c(d)-\frac{33d^2}{n^2}+\frac{16d^4}{27n^3}
\]
is nonincreasing in \(d\). Thus the conditional gap law is
dominated in likelihood-ratio order by weights
\(\exp(33d^2/n^2-16d^4/(27n^3))\) on its full half-grid,
with the same orientation factors as before. Beyond the finite
physical support the ratio is zero. In particular the conditional
second moment is at most that of the reference law, by the
nonpositive covariance of its increasing square function and
the decreasing density ratio. The full reference on signed gaps
is even, on either \(\mathbb Z\) or \(\mathbb Z+1/2\) according
to the midpoint parity. The even second moment on this
full grid agrees with the reference on the corresponding
half-grid with orientation factors.

For \(n\ge2048\), scale this full grid by \(h=n^{-3/4}\le1/300\).
The reference weight at a scaled point \(s\) is
\(f_0(s)=\exp(Cs^2-Ds^4)\), where \(C=33/\sqrt n\le3/4\)
and \(D=16/27\). Write \(I_j=\int_{\mathbb R}s^j f_0(s)\,ds\).
The integral of \((s f_0(s))'\) is zero, so
\(4D I_4=I_0+2C I_2\). Cauchy–Schwarz for \(I_4\) gives
\(I_2/I_0<21/20\), since
\(4D(21/20)^2-2(3/4)(21/20)-1>0\) and the quadratic
has only one positive root. Further \(I_0\ge2e^{-D}>1\).

For either shifted grid, the error between an integral and
its grid sum times \(h\) is bounded by \(h/2\) times total
variation. Indeed on each centered grid cell integrate the
difference from the center value and bound it by \(h/2\)
times the integral of the absolute derivative on that cell.
The function \(f_0\) has maximum \(\exp(C^2/(4D))<4/3\), so
its variation on the line is less than \(16/3\). The function
\(f_2(s)=s^2 f_0(s)\) has just one maximum on each open
half-line. Its maximum is less than \(8/3\), by
\(Cs^2\le Ds^4/2+C^2/(2D)\) and
\(\sup_s s^2 e^{-Ds^4/2}\le1/\sqrt D<4/3\).
Here \(e^{C^2/(2D)}<e^{1/2}<2\).
Its variation is therefore less than \(32/3\). The full-grid
reference second moment after scaling is bounded by
\[
 \frac{(21/20)I_0+16h/3}{I_0-8h/3}
    \le\frac{21/20+16/900}{1-8/900}<\frac{27}{25}.
\]
This estimate holds for every midpoint and both parities.
Hence \(\operatorname{Var}(K)=2\mathbb E[d^2]
\le (54/25)n^{3/2}\) for \(n\ge2048\). At
\(|\theta|\le2\pi/(n+1)\) the characteristic function bound
above has loss less than \(20(54/25)/\sqrt n<1\), completing
the claimed first-sector exclusion.

## Exact log concavity at smaller curvature

For any positive unimodal coefficient list \(c_0,\ldots,c_n\)
at a fixed radius, split the list by its distinct superlevel
heights. The resulting positive sum consists of interval
indicators \([a,b]\) all containing a common mode. For
\(0<\theta\le2\pi/(n+1)\), the geometric sum on \([a,b]\)
has argument \((a+b)\theta/2\) and a strictly positive
multiplying factor, except that the full interval can
vanish at the boundary. Because the intervals have a
common mode, their midpoint range has length at most
\(n/2\), so the nonzero sums lie in a cone of opening
at most \(n\theta/2<\pi\). Hence the total cannot vanish
unless the list is constant and the point is on the
boundary. This is the sector form of the
[unimodal root-of-unity criterion, Proposition 1.7](https://arxiv.org/pdf/1009.6022).
Strictly log-concave coefficients remain strictly log-concave
after any positive radial dilation, hence unimodal and nonconstant
for \(n\ge2\).

At any interior index \(1\le k<n\), write
\(u_i=(n-k)a_i+kb_i\), \(e_i=|d_i|\). The row's second
logarithmic difference is
\[
 e_i^2\int_{-1}^1(1-|v|)
                \psi_1(u_i+1+e_i v)\,dv
                   < \frac{e_i^2}{u_i-e_i+1/2}\quad(e_i>0).
\]
This follows from integrating the centered second difference
against the triangular kernel and using the midpoint trigamma
bound. The binomial contributes
\(-\log(1+1/k)-\log(1+1/(n-k))\), whose magnitude is at least
\((k+1/2)^{-1}+(n-k+1/2)^{-1}\) by Jensen. Thus an exact
rational sufficient test for strict log concavity is
\[
             \sum_i\frac{d_i^2}{u_i-|d_i|+1/2}
                    \le \frac1{k+1/2}+\frac1{n-k+1/2}
                      \qquad (1\le k<n).
\]
For a fully shared profile, \(u_i-|d_i|\ge n\min(a_i,b_i)\)
at every interior index. The left side is at most \(V/n\),
and the right side is at least \(4/(n+1)\).
This proves the fully shared claim, since \(4n/(n+1)\ge26/7\)
for \(n\ge13\). The smaller physical orders are covered by
the [terminal thirteen theorem](MULTIVARIATE_TERMINAL_THIRTEEN_ALL_PROFILES.md)
and its predecessors.

A single private unit row on the \(a\) side cancels the
\((n-k)!\) denominator of the binomial. The remaining
factor \(1/k!\) gives second log difference
\(-\log(1+1/k)\). Suppose \(\sum_{\rm shared}d_i^2/b_i\le1\),
so \(d_i^2\le b_i\) for each shared row.

For an increasing gap \(e_i=d_i>0\), integrality gives
\(a_i=b_i-e_i\ge e_i\), including \(e_i=1\) since the row is
shared. Hence \(a_i\ge b_i/2\), \(e_i\le b_i/2\), and the
centered denominator \(u_i+1/2\) is at least \(b_i(k+1/2)\).
In the triangular integral the perturbation \(e_i v\) is
the sum of two independent uniforms on
\([-e_i/2,e_i/2]\). Equivalently it is their average after
both are rescaled to \([-b_i/2,b_i/2]\), followed by
a contraction of factor \(2e_i/b_i\le1\).
For any convex function, Jensen bounds the expectation
under averaging and contraction by its expectation under
a single such uniform, because the uniform has mean zero.
Apply this to the reciprocal with center \(b_i(k+1/2)\),
after bounding the trigamma by the centered reciprocal.
The row difference is strictly less than
\[
 \frac{e_i^2}{b_i}\int_{-b_i/2}^{b_i/2}
             \frac{dv}{b_i(k+1/2)+v}
          =\frac{e_i^2}{b_i}\log(1+1/k).
\]
For a decreasing gap, \(a_i=b_i+e_i\) gives
\(u_i-e_i\ge b_i(k+1)\). The direct reciprocal bound is then
strictly smaller than \(e_i^2/[b_i(k+1)]\), which is at most
\((e_i^2/b_i)\log(1+1/k)\).
Summing either row bound and using \(\sum d_i^2/b_i\le1\)
proves strict log concavity, hence the first-sector exclusion.
