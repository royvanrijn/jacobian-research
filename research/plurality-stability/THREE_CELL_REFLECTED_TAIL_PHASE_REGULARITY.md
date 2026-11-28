# Smooth tail profile of a reflected three-cell maximizer

Fix \(0<r<1\) and let \(h\) be the increasing graph height of
a balanced three-cell Gaussian maximizer invariant under a central
reflection that exchanges two labels. Write \(\beta=U_3-U_1\)
for the constant pair score difference on the right graph.
The [reflection reduction](THREE_CELL_REFLECTION_REDUCTION.md)
gives positive upper and lower bounds on \(h'\) and a discretely
homogeneous limit
\[
    H(x)=\lim_{n\to\infty}r^n h(r^{-n}x),\qquad
    H(rx)=rH(x),\qquad |h(x)-H(x)|=O(1).
\]
Then \(H\) is smooth on \(x>0\). For every integer \(k\ge1\),
\[
 h^{(k)}(x)=O(x^{1-k}),\qquad
 h^{(k)}(x)-H^{(k)}(x)=O(x^{-k}) .
\]
Moreover, with \(\sigma=\sqrt{1-r^2}\) and
\(z_\beta=\Phi^{-1}((1+\beta)/2)\),
\[
 h(x)-H(x)
  = \frac{\sigma z_\beta}{r-1}\sqrt{1+H'(x)^2}+O(x^{-1}).
                                                               \tag{TP1}
\]
These bounds give a smooth slope as a function of logarithmic
phase. They allow that slope to vary over a period.

## Derivatives of the pressure equation

The [local affine limit argument](THREE_CELL_REFLECTED_SHAPE_RIGIDITY.md#local-affine-limits-at-infinity)
gives bounded derivatives of every fixed order, \(h''(x)\to0\),
and \(|h(rx)-r h(x)|\le C\) on the tail. Extend \(h\) smoothly
to the whole line with bounded derivatives of each positive order,
using its analytic extension near zero and a smooth cutoff on the
negative half-line.

Let \(Z\) be standard normal, \(S=rx+\sigma Z\), and
\(\psi(t)=2\Phi(t/\sigma)-1\). On the right graph the pressure
equation is
\[
              \mathbb E\psi(h(S)-r h(x))=-\beta+e(x).
                                                               \tag{TP2}
\]
For each fixed order \(j\ge0\), \(e^{(j)}(x)=O(e^{-c_j x^2})\)
with \(c_j>0\). Indeed the actual pair contrast is its binary
expression on \(S>0\); only \(S<0\) uses the reflected branch
and the extending graph. The correction there is a convolution
over a fixed negative half-line with horizontal density
\(p_\sigma(s-rx)\). Its vertical derivatives are uniformly bounded
by Gaussian smoothing. Differentiating at \(y=h(x)\) uses only
bounded derivatives of \(h\); horizontal derivatives of the
Gaussian insert polynomial factors and preserve a Gaussian tail.

Put
\[
 U=h(S)-r h(x),\qquad B(x)=\mathbb E\psi'(U),\qquad
 P_x(ds)=B(x)^{-1}\psi'(h(s)-r h(x))p_\sigma(s-rx)\,ds .
\]
The bound on \(h(rx)-r h(x)\) and the slope bound make \(B(x)\)
bounded below for large \(x\), by restricting to \(|S-rx|\le1\).
The kernel \(P_x\) is a probability density dominated by
\(C p_\sigma(s-rx)\). For each \(j\ge1\) define the random
displacement derivative
\[
                    \Delta_j=r^j h^{(j)}(S)-r h^{(j)}(x).
\]
Differentiation under the Gaussian expectation is justified by the
uniform bounds on derivatives of \(h\) and \(\psi\). For \(k\ge1\)
the \(k\)-th derivative of (TP2) gives
\[
 \mathbb E[\psi'(U)\Delta_k]=
   -\mathbb E
     \sum_{\substack{n_1,\ldots,n_{k-1}\ge0\\
                      \sum j n_j=k,\ N=\sum n_j\ge2}}
       c_{\mathbf n}\,\psi^{(N)}(U)\prod_j\Delta_j^{n_j}
      +O(e^{-c_kx^2}),                                    \tag{TP3}
\]
where \(c_{\mathbf n}=k!/\prod_j (j!)^{n_j}n_j!\) and the sum
is empty for \(k=1\). Thus the leading term of the equation
for \(h^{(k)}(x)\) is \(r^{k-1}\int h^{(k)}\,dP_x\).
Every other term is a product of at least two lower-order
displacement derivatives, with orders adding to \(k\).

For \(j\ge1\) let \(E_j(X)=\sup_{x\ge X}|h^{(j)}(x)|\).
If \(0<a<r\) is fixed, the mean-value estimate on \(S\ge a x\)
is
\[
 |\Delta_j-(r^j h^{(j)}(rx)-r h^{(j)}(x))|
       \le C_j E_{j+1}(a x)|S-rx| .                        \tag{TP4}
\]
All finite Gaussian moments of \(|S-rx|\) are uniform, both
under the original Gaussian and under \(P_x\). The omitted
event \(S<a x\) has Gaussian-small probability, and the
displacements are bounded there. In particular (TP4) controls
every fixed finite \(L^p\) norm of the displacements once their
values at \(S=rx\) are controlled.

## Curvature decay

For \(k=1\), (TP3) gives
\[
                  h'(x)=\int h'(s)\,P_x(ds)+O(e^{-c x^2}).
\]
Taking \(j=1\) in (TP4) shows
\(|h'(x)-h'(rx)|\le C E_2(a x)+O(e^{-c x^2})\).
Also \(\mathbb E|\Delta_1|^2\le C E_2(a x)^2+O(e^{-c x^2})\).
The square is the only nonlinear term of (TP3) for \(k=2\).
Consequently
\[
 E_2(X)\le r E_2(aX)+C E_2(aX)^2+O(e^{-cX^2}).
                                                               \tag{TP5}
\]
For any \(1/2<p<1\) choose \(a<r\) sufficiently close to \(r\)
so that \(r<a^p\). Since \(E_2(X)\to0\), the coefficient
\(r+C E_2(aX)\) in (TP5) is eventually at most
\(\theta a^p\) with \(\theta<1\). Multiplying by \(X^p\)
and iterating on \(X,aX,a^2X,\ldots\) down to a fixed compact
interval gives \(E_2(X)=O(X^{-p})\).

For \(D>0\) fixed and \(x\) large, put
\(t(x)=rx-D\sqrt{\log x}\). This is increasing and lies in
\([rx/2,rx]\) for all sufficiently large \(x\). In the leading
expectation of (TP3), keep \(S\ge t(x)\). The excluded probability
under \(P_x\) is \(O(x^{-K})\) for any prescribed \(K\) if \(D\)
is chosen sufficiently large. With the preliminary bound just
proved, the square in the second derivative equation gives
\[
 E_2(X)\le rE_2(t(X))+C X^{-2p}.
                                                               \tag{TP6}
\]
This gives \(E_2(X)=O(X^{-1})\) by the following iteration. Whenever
\[
 E(X)\le r^{k-1}E(t(X))+C X^{-q},\qquad q>k-1,
                                                               \tag{TP7}
\]
multiply by \(X^{k-1}\). The coefficient of
\(t(X)^{k-1}E(t(X))\) becomes
\((rX/t(X))^{k-1}=1+O(\sqrt{\log X}/X)\), and the error becomes
\(O(X^{-(q-k+1)})\). Along the iteration, each new \(X\) is at
most \(r\) times the last, and the terminal \(X\) lies in a fixed
compact interval bounded away from zero. In reverse order the
earlier values grow at least geometrically. Both the sum of the
coefficient excesses and the sum of the errors are bounded
uniformly, so \(X^{k-1}E(X)\) is bounded. Apply this with
\(k=2,q=2p>1\) in (TP6). The first derivative equation now gives
\(h'(x)-h'(rx)=O(x^{-1})\).

## Higher derivatives and convergence of the phase

We prove inductively for \(k\ge3\) that \(E_k(X)=O(X^{1-k})\)
and
\[
       r^{k-1}h^{(k-1)}(rx)-r h^{(k-1)}(x)=O(x^{-(k-1)}).
                                                               \tag{TP8}
\]
The previous section gives the base of the induction:
\(E_2(X)=O(X^{-1})\) and the corresponding order-one
displacement at \(S=rx\) is \(O(x^{-1})\).
Assume the bounds on \(E_j\) through \(j=k-1\) and on the
displacements at \(rx\) through \(j=k-2\). By (TP4),
\(\Delta_j=O(x^{-j})\) in every fixed finite Gaussian \(L^p\)
norm for \(j\le k-2\). For \(\Delta_{k-1}\) the crude bound
in the same norms is \(O(x^{2-k})\), using the bound on \(E_{k-1}\)
on \(S\ge a x\). Every product in the \(k\)-th equation of
(TP3) either contains only \(j\le k-2\), hence has expected
absolute value \(O(x^{-k})\) by Hölder's inequality, or is
\(\Delta_{k-1}\Delta_1\), with expected absolute value
\(O(x^{1-k})\). For \(k=3\) these are \(\Delta_1^3\) and
\(\Delta_2\Delta_1\); for \(k=4\) the square \(\Delta_2^2\)
is already in the first class.

Choosing any fixed \(a<r\) gives the preliminary recurrence
\[
                E_k(X)\le r^{k-1}E_k(aX)+C X^{1-k}.
\]
Take \(k-2<p<k-1\) and \(a\) sufficiently close to \(r\) that
\(r^{k-1}<a^p\). Multiplication by \(X^p\) and iteration with
that strict coefficient inequality yields \(E_k(X)=O(X^{-p})\).
The \((k-1)\)-st equation of (TP3) uses only products of
displacements with order at most \(k-2\), so
\(\mathbb E[\psi'(U)\Delta_{k-1}]=O(x^{-(k-1)})\).
Since \(B(x)\) is bounded below, (TP4) and the preliminary bound
on \(E_k\) now give
\[
 |r^{k-1}h^{(k-1)}(rx)-r h^{(k-1)}(x)|=O(x^{-p}).
\]
Using (TP4) again, \(\Delta_{k-1}\) has the same improved bound
in every fixed finite Gaussian \(L^p\) norm, with the moment
exponent allowed to differ from this decay exponent \(p\).
Thus the only borderline product in the \(k\)-th equation is
now \(O(x^{-p-1})\) in expected absolute value, and
\(p+1>k-1\). Restricting its leading expectation to
\(S\ge t(x)\) gives (TP7) for \(E_k\) with \(q=p+1\).
We obtain \(E_k(X)=O(X^{1-k})\). Substituting this bound into
the preceding use of (TP4) proves (TP8) and completes the
induction.

For any integer \(j\ge1\), (TP8) with \(k=j+1\) yields
\[
              h^{(j)}(x)-r^{j-1}h^{(j)}(rx)=O(x^{-j}).
                                                               \tag{TP9}
\]
Thus the \(j\)-th derivatives of \(r^n h(r^{-n}x)\) and
\(r^{n-1}h(r^{-(n-1)}x)\) differ by \(O(r^n x^{-j})\),
uniformly on each compact interval of positive \(x\) for
large \(n\). Summing these differences proves smooth convergence
to \(H\) on compact intervals and
\(h^{(j)}(x)-H^{(j)}(x)=O(x^{-j})\) for large \(x\).

If \(H\) is nonlinear, \(H'\) is a nonconstant periodic
function of \(\log x\). Its logarithmic derivative has
both signs. Since \(x h''(x)=x H''(x)+O(1/x)\), the
curvature of the true graph therefore takes both signs
arbitrarily far out, with size at least a positive
multiple of \(1/x\) on geometric sequences. Eventual
convexity or concavity forces a linear phase.

## The pressure-dependent height correction

Write \(b(x)=h(x)-H(x)\) and \(m(x)=H'(x)\). By the preceding
bounds, \(b=O(1)\), \(b'=O(x^{-1})\), and \(H''=O(x^{-1})\).
Discrete homogeneity gives \(H(rx)=rH(x)\) and \(m(rx)=m(x)\).
For \(S=rx+\sigma Z\) with \(S\ge a x\), Taylor approximation
about \(rx\) gives
\[
 h(S)=rH(x)+b(rx)+\sigma m(x)Z
                 +O\bigl((|Z|+Z^2)/x\bigr).
\]
The omitted event is Gaussian-small, and the first derivative
of \(\psi\) is bounded. Plugging the approximation in (TP2)
and integrating the linear Gaussian expression gives
\[
 \beta=2\Phi\!\left(
      \frac{r b(x)-b(rx)}
           {\sigma\sqrt{1+m(x)^2}}\right)-1+O(x^{-1}).
\]
Its argument is bounded, so inversion yields
\(r b(x)-b(rx)=\sigma z_\beta\sqrt{1+m(x)^2}+O(x^{-1})\).
Define \(b_*(x)=\sigma z_\beta\sqrt{1+m(x)^2}/(r-1)\),
which has \(b_*(rx)=b_*(x)\).
For \(v=b-b_*\), this says \(r v(x)-v(rx)=O(x^{-1})\).
Evaluating at \(x/r\) and iterating outward gives
\(v(x)=r^n v(x/r^n)+O(x^{-1}\sum_{j=1}^n r^j)\).
Since \(v\) is bounded, \(v(x)=O(x^{-1})\). This proves (TP1).
