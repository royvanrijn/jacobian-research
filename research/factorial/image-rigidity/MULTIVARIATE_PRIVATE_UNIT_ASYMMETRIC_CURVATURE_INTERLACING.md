# All-degree interlacing with private variables on both sides and residual curvature five

Consider an ordinary profile with a private variable
in each monomial. Define the residual curvature by
\[
                     S_*=\sum_{\max(a_i,b_i)>0}
                           \frac{(b_i-a_i)^2}{\max(a_i,b_i)}-2 .
\]
If there is a private unit pair, \(S_*\) is the ordinary
curvature of all rows outside that pair. No symmetry is
assumed.

**Theorem.** If \(S_*\ge5\), every consecutive pair
has strictly interlacing root arguments, listed from the
positive real axis. In particular the ordinary profile
satisfies SFC(2) in every degree over every
characteristic-zero field.

The same angular statement applies whenever the
endpoint-normalized coefficients of the geometric
moment polynomial at each degree satisfy the envelope
in (1) below. Only their positive moment representation
and the envelope are used.

## A common radial bound

At degree \(m\), reserve one private row on each side
and extract a factorial from each, using the
[positive Gamma factorization](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md#the-geometric-resolvent-without-reflection-symmetry).
The two extracted factorials cancel the binomial
coefficient. Normalize the endpoints to one by a
positive scaling of variable and polynomial, and write
\[
                        P_m(z)=1+\sum_{k=1}^{m-1}c_{m,k}z^k+z^m.
\]
The usual [factorial endpoint curvature estimate](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md#an-endpoint-coefficient-estimate)
applies on each unreserved row. For a reserved row of
exponent \(d\), the remaining factor is \((dj)!/j!\),
with \(j=k\) or \(j=m-k\). Each rising factorial of
shape \(h/d\), \(1\le h<d\), in that factor has log second
difference at an interior index \(j\) at least
\[
                  \log\frac{j+h/d}{j-1+h/d}
                                  \ge\frac1{j+h/d}\ge\frac1m .
\]
Thus the residual contribution of this row is at
least \((d-1)/m\), exactly the bound used after
subtracting one from its full curvature \(d\).
After both reservations the endpoint chord estimate
therefore gives
\[
                0<c_{m,k}\le \exp\left(-\frac{5k(m-k)}{2m}\right).
                                                               \tag{1}
\]

At radius \(R=e^{1/m}\), group indices \(k,m-k\) for
\(1\le k\le m/2\). Their envelope is at most
\(e^{-5k/4}\), and \(R^k+R^{m-k}\le 1+e\).
Doubling the even middle term only increases the bound.
Thus
\[
 \left|\sum_{k=1}^{m-1}c_{m,k}z^k\right|
     \le (1+e)\frac{e^{-5/4}}{1-e^{-5/4}}
     < \frac{15}{4}\frac{29}{71}<\frac85<e-1.
\]
Rouché's theorem, with \(1+z^m\), and the reciprocal
argument therefore give
\[
                            e^{-1/m}<|z|<e^{1/m}
                                                               \tag{2}
\]
at every root.

## Controlling both endpoint phases

By the [asymmetric angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md#the-geometric-resolvent-without-reflection-symmetry),
at each upper-half-plane root \(z=r e^{i\theta}\),
\[
               m\theta=(2k-1)\pi+p_\theta(a)+p_\theta(b),
             \qquad p_\theta(u)=\arg(1-u e^{-i\theta}).
                                                               \tag{3}
\]
Here \(a,b\) are the two endpoint ratios of the positive
resolvent. We show that their phases are each below
\(\theta/2\) whenever
\[
                       m\ge4,\qquad \frac{\pi}{m-1}
                                        \le\theta\le\frac{\pi}{2}.
                                                               \tag{4}
\]

Either ratio is \(B/A\), with
\[
         B=\mathbb E[s/D],\qquad A=\mathbb E[1/D],\qquad
                           D=1-2cs+s^2,\quad c=\cos\theta.
\]
On the left this uses \(s=rT\) under the degree-\(m\)
left endpoint measure. On the right it uses \(s=1/(rT)\)
under the right endpoint measure, which is the left
measure reweighted by \(T^m\) and normalized. The two
corresponding moments \(u_j=\mathbb E s^j\) are
\(r^j c_{m,j}\) and \(r^{-j}c_{m,m-j}\), respectively.
By (1) and (2), both satisfy
\[
  u_j\le \exp\left(-\frac{5j(m-j)}{2m}+\frac{j}{m}\right),
   \qquad
  u_1\le e^{-5/2+7/(2m)},\qquad u_2\le e^{-5+12/m}.             \tag{5}
\]
In (4), set \(C_\theta=1/[2(1-c)]\le(m-1)^2/9\).
The latter inequality follows from sine concavity on
\([0,\pi/6]\), since
\(\sin(\theta/2)\ge\sin(\pi/(2(m-1)))\ge3/(2(m-1))\).
We have
\[
           A\ge1+2c u_1-u_2,\qquad
           s/D\le C_\theta,\qquad
           s/D\le s+C_\theta s^2 .
\]
For the first inequality use \(1/D\ge2-D\).
The last two follow from
\(D-2(1-c)s=(s-1)^2\) and
\(D-2(1-c)(2c-s)=(s-(2c-1))^2\).

It suffices to show \((1+2c)B<A\) for both endpoints.
For \(4\le m\le13\),
\[
      (1+2c)B-A\le u_1+\bigl(1+(m-1)^2/3\bigr)u_2-1 .
\]
Set \(g(x)=(1+(x-1)^2/3)e^{-5+12/x}\).
For \(x\ge4\), the sign of its derivative is the sign of
\(x^3-7x^2+12x-24\), a strictly increasing function on
that interval. Hence \(g\) has no interior maximum on
a closed subinterval of \([4,\infty)\).
On \([4,7]\), both endpoint values of \(g\) are below
\(11/20\), and \(u_1<1/5\).
On \([8,13]\), both endpoint values of \(g\) are below
\(17/20\), and \(u_1<13/100\).
This proves the strict bound up to thirteen.

For \(m\ge14\), on \(s\le2/3\) we have
\[
          s/D\le \frac{s}{(1-s)^2}\le s+12s^2.
\]
On the complementary event use \(s/D\le C_\theta\)
and apply Markov's inequality with \(k=\lfloor m/2\rfloor\).
With \(k\ge7\) and \(m-1\le2k\), (5) gives
\[
 \begin{split}
      B&\le u_1+12u_2+T_m,\qquad
          T_m=C_\theta(3/2)^k u_k,\\
      u_1&<53/500,\qquad u_2<2/125,\\
      T_m&\le \frac{(2k)^2}{9}\,e^{1/2}
                           \bigl((3/2)e^{-5/4}\bigr)^k<1/10 .
 \end{split}
\]
For the tail, replace \(e^{1/2}\) by \(5/3\) and
\((3/2)e^{-5/4}\) by \(43/100\). The resulting expression
decreases for \(k\ge7\), and its value at seven is below
\(1/10\). Consequently
\[
    (1+2c)B-A\le u_1+37u_2+3T_m-1
                 <\frac{53}{500}+37\frac{2}{125}+\frac3{10}-1<0 .
\]
This proves \((1+2c)B<A\) throughout (4).
Comparing the tangents
\((B/A)\sin\theta/(1-(B/A)c)\) and \(\tan(\theta/2)\)
then gives \(p_\theta(B/A)<\theta/2\).

The exponential comparisons can be checked by the positive
Taylor sum through degree twelve for each lower bound on
\(e^x\). The upper bounds \(e<11/4\) and \(e^{1/2}<5/3\)
follow from the Taylor sum through degree six with the
positive geometric tail estimate.

## Strict adjacent angular order

Consider an upper-order root with label \(k\) also
present in the lower order \(m-1\). It suffices to put
its angle below the lower box's left endpoint,
\((2k-1)\pi/(m-1)\). If its angle exceeds \(\pi/2\),
its own upper box implies \(4k>m+1\), which already
puts that upper box below the lower left endpoint.
Otherwise, if the desired strict bound has yet to
hold, the angle is in (4).

Both phases are then below \(\theta/2\), so (3) implies
\((m-1)\theta<(2k-1)\pi\), as required.
The other side of the interlacing, between a lower
\(k\)-th root and an upper \((k+1)\)-st root, follows
from the boxes. The pairs with upper order two and
three also follow from the boxes, reflection and the
unique odd negative real root. Thus every pair
strictly alternates. The ordinary polynomials have
rational coefficients before positive normalization,
so the resulting rational Bézout identity gives
SFC(2) over every characteristic-zero field.
