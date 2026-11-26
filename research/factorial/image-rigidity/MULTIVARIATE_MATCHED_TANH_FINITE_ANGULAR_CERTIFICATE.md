# An all-degree angular certificate from one pair of matched tanh moments

Let \(a,b\) be distinct nonnegative integer exponent columns, omitting
equal rows. For a row of smaller exponent \(\ell\) and positive gap
\(d=|b_i-a_i|\), give each of its \(d\) Gauss factors the row's
orientation \(\epsilon\in\{+,-\}\). For factor \(h\in\{1,\ldots,d\}\)
put \(v=d/(2\ell+d)\), \(c=2h/d-1\), and
\[
 T_j=\frac{(j-1)v}{j+cv},\qquad
 W=\sum_{\rm rows}2dv .
\]
These are the [exact matched Gauss coordinates](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid).
Private rows have \(\ell=0\). Fix a rational \(K>0\) with \(W\le K\).
Take the observations at the single reference pair \(t=T_2,\ s=T_3\).
The observations at that pair are the positive mixed moments
\[
 M_{p,q}^{\epsilon}=\sum_{\text{factors in orientation }\epsilon}t^p s^q,
                      \qquad p\ge1,\quad q\ge0 .
\]

**Theorem.** Let \(J\ge4\) be any integer at least the
[effective absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
For every rational precision \(\eta>0\), a computable finite set of
these moments gives the endpoint-normalized log coefficients in
every degree \(2\le j<J\) to absolute error at most \(\eta\).
The approximations are explicit rational linear combinations of
the moments. The error is uniform in the dimension and the
exponents under the sole bound \(W\le K\).

The resulting rational coefficient boxes give one finite certificate
for disjoint consecutive root arguments in every degree. On a fixed
ordinary profile the certificate succeeds at sufficiently high
precision if and only if its consecutive root arguments are disjoint
in every degree. On any family with a uniform positive angular gap
through \(J-1\), one finite moment index set and a positive input
error tolerance suffice for the whole family. The construction
does not assume a separation margin for the general ordinary class.

## Uniform transfer from the reference pair

For all \(j\ge2\), the row bound \(2d^2/(2\ell+d)\le K\) implies
\[
 0<T_j=\frac{(j-1)d}{(j-1)d+2j\ell+2h}
 \le\frac{K(j-1)}{K(j-1)+4}.
 \tag{1}
\]
Indeed \(K(2j\ell+2h)\ge4d\) follows from
\(d(j\ell+h)\ge2\ell+d\), using \(j\ge2\) and \(d,h\ge1\).
The grid also gives \(t\le v\le3t\) and \(T_j\le v\) for
all \(j\ge2\). Thus \(\sum t\le K/2\) and \(T_j\le3t\).
In particular \(2\sum t\le W\le6\sum t\), so the required
cap itself can be certified from the unsigned first moment.

Define
\[
 \delta=\frac{3s}{4t}-1,\qquad
 \kappa_j=\frac{2(j-1)}{j},\qquad \beta_j=\frac{2(3-j)}j .
\]
The matched map gives
\[
       \delta=\frac{ct}{3-ct}=\frac{cv}{6+2cv},\qquad
       T_j=\kappa_j t\,\frac{1+\delta}{1+\beta_j\delta}.
       \tag{2}
\]
For the second equality, the exact grid gives
\((j-1)/T_j=(3-j)/t+2(j-2)/s\).
Insert \(s=(4/3)t(1+\delta)\) and simplify.
In particular \(|\delta|\le t/(3-t)\le t/2\).
Also \(|cv|\le1\) gives \(|\delta|\le1/4\), and \(|\beta_j|\le2\).
Let \(\mathcal C_L\) denote the Chebyshev polynomial with
\(\mathcal C_L(\cos \phi)=\cos(L\phi)\). For \(L\ge1\), the expression
\[
 R_L(t)=\frac{1-\mathcal C_L(1-2t)}{2L^2t}
\]
is a rational polynomial, with \(R_L(0)=1\) by continuity.
In fact for \(1-2t=\cos\phi\), it is
\((\sin(L\phi/2)/(L\sin(\phi/2)))^2\).
The bound \(|\sin(Lx)|\le L|\sin x|\) follows from the geometric sum
for the ratio. Hence \(0\le R_L(t)\le1\) and \(tR_L(t)\le1/L^2\)
on \([0,1]\).

Put \(D=(3s-4t)/4=t\delta\). For integers \(L\ge1\) and \(H\ge0\) let
\[
 \delta_L=D\,\frac{1-R_L(t)}{t},\qquad
 Y_j(t,s)=\kappa_j t(1+\delta_L)\sum_{r=0}^H(-\beta_j\delta_L)^r .
\]
The quotient defining \(\delta_L\) is a rational polynomial
because \(R_L(0)=1\), and \(Y_j\) is divisible by \(t\).
On the Gauss grid \(\delta_L=\delta(1-R_L(t))\), so
\[
 |\delta-\delta_L|\le tR_L(t)/2\le1/(2L^2).
\]
Both \(|\beta_j\delta|\) and \(|\beta_j\delta_L|\) are at most
\(1/2\). The derivative of \((1+x)/(1+\beta_j x)\) has modulus
at most \(12\) between \(\delta\) and \(\delta_L\).
The geometric remainder at \(\delta_L\), including the prefactor
\(\kappa_j(1+\delta_L)\), is bounded by \(5\,2^{-H-1}\).
Thus
\[
 |T_j-Y_j(t,s)|\le t\,\omega,\qquad
 \omega=\frac{12}{L^2}+5\,2^{-H-1}.                            \tag{3}
\]
The error is weighted by \(t\) even for factors arbitrarily
close to the origin. This makes it summable uniformly in the
number of rows.

## Approximating the exact log coefficients

For \(0\le h'<j\) put \(b_{j,h'}=2h'/(j-1)-1\).
An increasing Gauss factor contributes to the normalized log
coefficient at \(k\) the exact expression
\[
 \ell_{j,k}(y)=\sum_{h'=0}^{k-1}\log(1+b_{j,h'}y)
       -\frac{k}{j}\sum_{h'=0}^{j-1}\log(1+b_{j,h'}y),
             \qquad y=T_j .
 \tag{4}
\]
In fact, if the order-\(j\) Gauss base is \(B\), then
\(T_j=(j-1)/(2B+j-1)\) and \((B+h')/(B+(j-1)/2)=1+b_{j,h'}T_j\).
The normalization cancels the scale \(B+(j-1)/2\).
The decreasing contribution is \(\ell_{j,j-k}(y)\).
For both orientations the log endpoint coefficients vanish.
The binomial coefficient is included after exponentiation.

Let \(\ell_{j,k}^{[P]}\) be (4) with each logarithm replaced
by its degree-\(P\) Taylor polynomial at zero. Put
\(\rho=K(J-1)/(K(J-1)+4)\) and \(\sigma=(1+\rho)/2\).
Take \(\omega\le(1-\rho)/2\), so \(|Y_j|\le\sigma\).
The log series remainder for \(|b|\le1\) at \(T_j\) is at most
\(T_j\rho^P/(1-\rho)\).
The derivative of the Taylor polynomial for \(\log(1+b y)\)
has modulus at most \(1/(1-\sigma)\) for \(|y|\le\sigma\).
The sums of absolute weights in (4) are at most \(2j\).
By (3) and \(T_j\le3t\), the error on summing over the
whole profile is therefore at most
\[
    JK\left(\frac{3\rho^P}{1-\rho}+
                         \frac{\omega}{1-\sigma}\right).       \tag{5}
\]
Choose \(L,H,P\) so that (5) is at most \(\eta\).
The polynomials \(\ell_{j,k}^{[P]}(Y_j(t,s))\) and
\(\ell_{j,j-k}^{[P]}(Y_j(t,s))\) are divisible by \(t\),
as their outer polynomials have zero constant term.
Their signed or separate orientation sums therefore use
only the prescribed positive mixed moments. If the moments
are measured with error, the sum of the absolute rational
coefficients of these two polynomials gives the error multiplier.
One can reserve half of \(\eta\) for this contribution.

## Certifying every degree

For \(2\le j<J\), construct rational boxes for the normalized
coefficients using (5), rigorous enclosures of the exponential,
and the bounds \(0\le [z^k]P_j\le\binom jk\) from log convexity.
Keep both endpoints equal to one. For each \(2\le n<N=n+1<J\)
ask whether polynomials in the two boxes have roots on a common
ray. Writing \(z=x+iy\) and the upper root as \(Rz\) with real
\(R>0\) gives four polynomial equations in real variables,
with polynomial inequalities for the boxes. All such roots
have moduli in \([2^{-j},2^j]\) in degree \(j\), by the monic
root bound and reversal. One may therefore restrict \(z\)
to its compact annulus and \(R\) to \([2^{-n-N},2^{n+N}]\).
This is a finite real algebraic feasibility test with
rational data. Infeasibility at all these pairs proves the
desired separation. At \(N=2\) the lower root is real and the
quadratic factorial moment is strictly positive on the real axis.
At \(N\ge J\) the absolute angular theorem applies.

Here is a quantitative stopping condition. If every trial
coefficient differs from its exact coefficient by at most
\(\binom jk \varepsilon\) with \(0<\varepsilon<1\), then every trial
root \(\widehat z\) in degree \(j\) is within
\(\varepsilon^{1/j}(1+|\widehat z|)\) of some exact root \(z\).
Indeed the exact monic polynomial has modulus at most
\(\varepsilon(1+|\widehat z|)^j\) there. Factor it and take a nearest root.
The inequality
\[
 \left|\frac{\widehat z}{|\widehat z|}-\frac{z}{|z|}\right|
 \le2\frac{|\widehat z-z|}{|\widehat z|}
\]
and the annulus then bound unit-projection error by
\(2(1+2^J)\varepsilon^{1/J}\). Thus a minimum exact chord gap
\(\gamma>0\) between unit projections across the pairs below \(J\)
makes every test infeasible as soon as
\(4(1+2^J)\varepsilon^{1/J}<\gamma\).
Such \(\varepsilon\) can be achieved uniformly by decreasing (5)
and the enclosure width of each exponential. If an exact pair
has a common ray, it is feasible in every box containing its
coefficients.

The [absolute curvature sixteen theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_16.md)
permits using \(K=32\) for the entire remaining curvature range.
If \(S=\sum d_i^2/\max(a_i,b_i)\ge16\) it already gives separation.
Otherwise \(W\le2S<32\), which meets the uniform cap.
The certificate therefore requires one fixed reference pair
even for the full bounded-curvature obligation. The remaining
problem is to prove that the relevant ordinary moment data
give a positive angular gap throughout that obligation.
