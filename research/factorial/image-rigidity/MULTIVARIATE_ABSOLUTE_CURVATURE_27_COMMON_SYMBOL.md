# Absolute-curvature angular separation for factorial-symbol blends

For distinct fully shared ordinary exponent columns \(a,b\),
let \(L_j\) be their endpoint-normalized log-factorial
symbol per degree. Put
\[
               S=\sum_i\frac{(b_i-a_i)^2}{\max(a_i,b_i)}
\]
over the unequal rows. At \(n\ge1,\ N=n+1\) use
independently chosen \(t_n,t_N\in[0,1]\) and set
\[
       \widehat L_m=(1-t_m)L_n+t_m L_N,\qquad
       P_m(z)=\sum_{k=0}^m\binom mk e^{m\widehat L_m(k/m)}z^k,
                       \qquad m=n,N .
\]

**Theorem.** If \(S\ge27\) the polynomials have disjoint
root arguments at every consecutive degree, under
independent positive rescalings. In particular this holds
with either common symbol on the two grids and for every
interpolation between the symbols. No reflection of rows
or unit-circle root hypothesis is needed.

The proof retains the numerical margins of the
[absolute-curvature comparison](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md).
It supplies the same discrete curvature inequality on
both grids for either symbol, then uses the matched
Gamma increment to bound the change of each of the
two analytic halves.

## Curvature on both grids

Write \(\delta_i=b_i-a_i\),
\(x_i(q)=a_i+q\delta_i\), and
\(V=\sum_i \delta_i^2/\min(a_i,b_i)\).
For \(z\ge0\) the polygamma series and the convex
trapezoid bound give
\[
 \psi'(1+z)\ge \frac1{1+z}+\frac1{2(1+z)^2}
                                 \ge \frac1{z+2/3}.
\]
The first inequality follows by comparing the
integral of \(t\mapsto(1+z+t)^{-2}\) on each
unit interval to the average of its endpoint
values. The second follows by clearing positive
denominators and using \(1+z\ge1\).
Therefore the row's contribution to \(L_j''(q)\)
is at least \(\delta_i^2/(x_i(q)+2/(3j))\).

For either \(m,j\in\{n,N\}\) and
\(1\le k<m\), express the second difference
of \(mL_j(k/m)\) as \(1/m\) times the average
of \(L_j''\) against the symmetric triangular
density of radius \(1/m\) centered at \(k/m\).
Convexity of the reciprocal gives a bound
\[
    \Delta_k^2\big(m L_j(k/m)\big)\ \ge\
     \frac1m\sum_i\frac{\delta_i^2}
                      {x_i(k/m)+2/(3j)}
                                 \ \ge\ \frac S m .
\]
Indeed \(x_i(k/m)\le\max(a_i,b_i)-|\delta_i|/m\)
and \(2/(3j)\le1/m\) for \(n\ge2\).
Either convex blend has the same lower bound.
Using the two endpoint values now gives
\[
                c_{m,k}:=\binom mk e^{m\widehat L_m(k/m)}
                 \le\binom mk e^{-S k(m-k)/(2m)}.            \tag{1}
\]
For \(n=1\) the first polynomial is \(1+z\) and
the second has two endpoints one and a middle
coefficient strictly less than two by strict
convexity of both symbols. It has no root on the
negative real axis. Thus only \(n\ge2\)
is needed below.

## Matched logarithms of the two halves

Use the same split into halves as in the
absolute-curvature proof, allowing the whole
even middle coefficient in one half. On all
coefficient regions of that proof, (1) gives
the same root annulus and same strict bound
on the magnitude of each half's angular log
derivative. In particular the small-half
bound on \(|z|\le e^{1/m}\) gives this derivative
less than \(1/4\). In the
[large-slope region](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md#an-endpoint-coefficient-estimate)
the coefficient bound is \(8^{-k}\) in each
half. Here \(n\ge255\) and \(e^{1/n}/8<1/7\),
so on \(|z|\le e^{1/m}\) each half differs
from one by at most \(1/6\) and its angular
log derivative has magnitude at most
\(7/30<1/4\), by geometric summation. In
these cases the root annulus follows by the
same dominance at \(|z|=e^{1/m}\) and the
reciprocal circle. Integrating the angular
derivatives of the four half logs gives
strictly less than \(\theta\) error in the
adjacent phase equations at any common
upper argument \(0<\theta\le\pi\).

In the remaining cluster region use
\(D=N e^{-S/2}\), \(K=\lceil32D\rceil+12\)
and the existing bounds
\[
    D>1/6,\quad D/n<1/720000,\quad n>59000,\quad
    (1+V)K/n<297/2500,\quad S/n<1/2100.
\]
The
[high-curvature symbol interpolation estimate](MULTIVARIATE_UNIFORM_EVENTUAL_MIRRORED_COMMON_SYMBOL_INTERLACING.md#unrestricted-curvature-of-ordinary-integer-profiles)
uses the kernel
\[
       0\le\partial_\nu\frac1{\nu(e^{u/\nu}-1)}
                         \le\frac1{2\nu^2}.
\]
Its endpoint-log coefficient estimates apply to
either order, either orientation of the halves
and all blend parameters as follows.

For the left half at fixed \(m\), use the
continuous coefficient logarithm
\[
 f_t(x)=\log\Gamma(m+1)-\log\Gamma(m+1-x)-x\log m
        +m\big((1-t)L_n(x/m)+tL_N(x/m)
                 -x[(1-t)L_n(1/m)+tL_N(1/m)]\big).
\]
It gives \(\log(c_{m,k}k!/B_t^k)\) at \(x=k\), where
\(B_t=c_{m,1}\). For the right half reverse the
rows. With \(L=\max(1,\max_i|\delta_i|/a_i)\) in the
chosen orientation and \(\alpha=2L/n\), the same
bounds as there hold for \(x\in[0,K]\), \(2\le r\le K\):
\[
 \begin{split}
 H_r&=(r-2)!\alpha^{r-2}4(1+V)/n,\\
 |\partial_x^r f_t|&\le H_r,\qquad
 |\partial_t\partial_x^r f_t|
             \le H_r(r-1)/(2n^2).
 \end{split}
\]
The kernel estimates use \(x_i(x/m)\ge a_i/2\),
\(m-x\ge m/2\), \(m\ge n\), and
\(\sum_i |\delta_i|^r/a_i^r\le V L^{r-2}\).
The relative bound on \(\partial_t\log B_t\) is
\(V/(4nN)<1/n\), by integrating the curvature
increment on the full interval and comparing
to its endpoint chord. Also (1) at \(k=1\)
gives \(B_t\le D e^{S/(2n)}<101D/100\).
These estimates do not depend on whether a
blend increases or decreases the coefficient
relative to the actual symbol at that order.

Consequently the forward differences, exponential
edge weights, and their blend derivatives have
the same majorants as the
[cluster connected-hypergraph proof](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md#sharper-connected-sum).
The truncated log changes by at most \(13D/n\)
on \(|z|\le6/5\) over any blend interval of length
at most one. The analytic log reconstruction
costs less than \(1/120\) at either endpoint.
The polynomial tail estimate used there follows
again from (1), including for the degree-\(N\)
blend. Thus each half log differs from its
actual-order half log by less than \(1/55\).
Its individual norm is at most \(M=3(D+1)\).
The absolute curvature proof bounds the
actual adjacent log difference by \(1/55\).
Including both blends yields
\[
         \|\log H_m\|_{6/5}\le M,\qquad
         \|\log H_N-\log H_n\|_{6/5}\le\delta:=3/55          \tag{2}
\]
for either orientation.

## Off-circle comparison

By Rouché the root log-radii \(q_m\) are
confined to \(|q_m|\le7(D+1)/m\), as in the
[absolute-curvature radial comparison](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md#paired-radius-and-phase-comparison).
For a proposed common upper argument \(\theta\),
let \(J_m(q,\theta)\) be the difference of the
two half logs at \(e^{q+i\theta}\) and
\(e^{-q-i\theta}\). They are evaluated at
modulus at most \(21/20\). Cauchy's estimates
applied to (2) give
\[
 \begin{aligned}
 |J_N-J_n|&\le2\delta,&
 |\partial_\theta(J_N-J_n)|&\le14\delta,\\
 |\partial_qJ_m|&\le14M,&
 |\partial_q\partial_\theta J_m|&\le210M .
 \end{aligned}
\]
The two real root equations imply
\[
    (n-14M)|q_N-q_n|\le 7(D+1)/n+2\delta ,
                       \qquad n|q_N-q_n|<1/9.
\]
For the last bound use \(M/n\le21D/n\),
\(7(D+1)/n\le49D/n<1/10000\) and
\(14M/n<1/2000\).
Since the logs are real on positive radii,
integrating the two angular derivative
estimates gives
\[
 \big|\Im(J_N(q_N,\theta)-J_n(q_n,\theta))\big|
 \le (14\delta+210M|q_N-q_n|)\theta<4\theta/5 .
\]
The subtracted phase equation would place an
integer multiple of \(2\pi\) between zero
and \(2\theta\le2\pi\), a contradiction.
Real coefficient conjugation and positivity
exclude the remaining angle choices.
