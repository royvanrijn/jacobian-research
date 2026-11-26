# Effective all-degree stability of comparable shared cores under mixed tails

Start with an ordinary two-monomial support with a private
variable on each side. Orient its columns so the shared
rows of the core have \(a_i\le b_i\).
Append arbitrary ordinary shared rows with \(a_i,b_i\ge1\),
possibly of mixed orientations, and put
\[
        \Delta=\sum_{\rm tail}(b_i-a_i),\qquad
        \tau=\sum_{\rm tail}\frac{(b_i-a_i)^2}{\max(a_i,b_i)}.
\]
Suppose \(\Delta\ge0\). There is an explicit positive bound
on \(\tau\), in terms of an absolute ordinary terminal cutoff
alone, that guarantees SFC(2) for the augmented support
in every degree. In particular this gives an effective
strengthening of the
[balanced-tail stability theorem](MULTIVARIATE_COMPARABLE_CORE_BALANCED_TAIL_STABILITY.md)
and allows tails with nonzero total degree difference.

More precisely, take any \(J\ge4\) at least the
[computable absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Set
\[
 \begin{gathered}
 \sigma=1/128,\quad M=2^\sigma\exp(20000\sigma(J+\sigma)),
       \quad p=2^{1-J}/J,\quad \delta=(p/(12M))^{1/\sigma},\\
 \rho=\frac{p}{2(1+\delta^{-1})},\quad L=3J\,3^J,\qquad
       \tau_0=\min\left\{\frac1{10},
                          \frac{\rho}{96LJ^2\,3^J}\right\}.
 \end{gathered}
\]
The sufficient condition is \(\tau<\tau_0\). No bound
on the number of appended rows or the integer \(\Delta\)
is required. The constant is small, and \(J\) is unevaluated.

## One physical comparison point

At \(\tau=0\) the comparable-core theorem applies directly.
For \(\tau>0\) the absolute eventual theorem handles terminals \(N\ge J\),
and the [absolute curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md)
handles total ordinary curvature at least \(27\).
The first and second moments cannot share a root, because
the first root is real and the second moment is strictly
positive on the real axis. Thus consider \(N=n+1<J\),
\(n\ge2\), with total curvature below \(27\).

Reserve the two private factorials as in the
[both-private geometric representation](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md#the-geometric-resolvent-without-reflection-symmetry).
Let \(U,V\) be the residual positive core bases,
\(T=V/U\), and use \(\mathbb E_j\) for the normalized
tilt by \(U^j\). Put
\[
  g_j=(\mathbb E_j T^j)^{1/j},\qquad
  \Pi_j(w)=\mathbb E_j G_j(wT/g_j),\qquad
                         G_j(w)=1+w+\cdots+w^j .
\]
For the tail let \(D_j\) be its positive endpoint scale.
The normalization of the full polynomial at physical
argument \(z\) is \(w_j=g_jD_jz\).
Its normalized coefficients are those of \(\Pi_j\)
multiplied by the normalized tail coefficients
\(t_{j,k}\) at the corresponding indices. All these
coefficients are at most one by Hölder.

The
[tail curvature and center bounds](MULTIVARIATE_COMPARABLE_CORE_BALANCED_TAIL_STABILITY.md#the-balanced-tail-disappears-at-both-orders)
also give, before setting the degree difference to zero,
\[
 1-t_{j,k}\le\tau j/4,\qquad
       \left|\log(D_N/D_n)-\Delta\log(N/n)\right|
                               \le\frac{\tau}{nN(1-\tau)}
\]
when \(\tau\le1/10\). The curvature estimate uses the
second derivative bound \(\tau/[j(1-\tau)]\) and the
endpoint chord. Put \(\kappa=D_N/D_n\) and
\(\widehat\kappa=\max(1,\kappa)\).
Since \(\Delta\ge0\), we have
\[
                          0\le\log(\widehat\kappa/\kappa)\le\tau .
                                                               \tag{1}
\]

Suppose \(z\) is a common root, and put
\(w=w_n=r e^{i\theta}\) and \(C=g_N/g_n\).
The full normalized leading and constant coefficients
are one, and the other coefficients are at most one.
Comparison with the leading term at radius two and
reversal give \(1/2<|w_j|<2\).
Since \(w_N=C\kappa w\), (1) implies
\(rC\widehat\kappa<3\) and \(C\widehat\kappa<6\).
The coefficient error on either full root is at most
\(\tau J2^J\) in comparison with the normalized core.
To evaluate the upper core polynomial at \(C\widehat\kappa w\)
instead of \(w_N\), integrate its derivative along the
positive rescaling; the log-radius interval is at
most \(\tau\) and \(|w\Pi_N'(w)|<2J\,3^J\) on the
disk of radius three. Thus
\[
       |\Pi_n(w)|<\varepsilon,\qquad
       |\Pi_N(C\widehat\kappa w)|<\varepsilon,\qquad
                                \varepsilon=L\tau .             \tag{2}
\]

## The retained likelihood-ratio comparison

Let \(\mu\) be the law of \(X=T/g_n\) under \(U^n\),
and \(\nu\) its law under \(U^nV\).
The
[conditional Efron comparison](MULTIVARIATE_BOTH_PRIVATE_COMPARABLE_SHARED_INTERLACING.md#why-comparable-shared-ordinary-rows-satisfy-the-test)
makes \(d\nu/d\mu\) nondecreasing. Also the law of
\(\log X\) under \(\mu\) is log-concave, by convolution
of the independent log-Gamma factors. Thus the law
obtained by multiplying \(X\) by \(\widehat\kappa\ge1\)
under \(\mu\) is higher in likelihood-ratio order
than \(\mu\). Indeed for its log-density \(f\), the
ratio \(f(y-a)/f(y)\) is increasing for \(a\ge0\),
by log concavity. Dilation preserves the comparison
between \(\mu\) and \(\nu\); transitivity shows that
the law \(\widehat\nu\) of \(\widehat\kappa X\) under
\(U^nV\) is higher in likelihood-ratio order than
\(\mu\). If \(X\) is constant this comparison is just
the order of two point masses.

The same \(\widehat\nu\) is the size-biased law of
\(\widehat\kappa X\) under \(U^N\).
For \(b=\mathbb E_N(\widehat\kappa X)>0\), the identity
\(G_N(s)=1+sG_n(s)\) gives
\[
               \Pi_N(C\widehat\kappa w)
                           =1+w b\,\mathbb E_{\widehat\nu}G_n(wX),
                                                               \tag{3}
\]
with \(X\) on the right denoting the variable of
integration under \(\widehat\nu\).
For later use, if \(I(s)=\operatorname{Im}G_n(se^{i\theta})\),
\[
 b\,\mathbb E_{\widehat\nu} I_+(rX)
        \le \sum_{k=1}^n r^k\mathbb E_N[(\widehat\kappa X)^{k+1}]
        \le C\widehat\kappa\sum_{k=1}^n (rC\widehat\kappa)^k
        <6J\,3^J .                                           \tag{4}
\]
The middle step uses the normalized upper core
coefficient bounds through degree \(N\).

## A quantitative reserve at the imaginary sign change

Conjugate the alleged common root into the upper
half-plane. The
[both-private angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md)
exclude common real roots and place this angle in the
same positive label \(k\) at both orders, with
\[
        (2k-1)\pi<n\theta<N\theta<(N+1)\theta<2k\pi,
                              \qquad 0<\theta<\pi/2.
\]
Set \(q=\sin\theta,\ c=\cos\theta,\
A=-\sin(N\theta),\ B=-\sin(n\theta)\).
We have \(q\ge2/n\ge2/J\). Also
\(N\theta-(2k-1)\pi\in(\theta,\pi-\theta)\), so \(A\ge q\).
The imaginary function \(I(s)\) has exactly one
positive-to-negative crossing, by its formula
\[
           I(s)=\frac{s[q+s^n(A-Bs)]}{1-2cs+s^2}.
\]
If \(H(s)=\operatorname{Re}G_n(se^{i\theta})\),
the exact geometric identity
\((1-se^{i\theta})G_n(se^{i\theta})=1-s^Ne^{iN\theta}\)
also gives
\[
                      (1/s-c)I(s)=qH(s)+As^n.                 \tag{5}
\]
For \(s=rX\) under \(\mu\), we have
\(A\mathbb E_\mu s^n=Ar^n\ge p\), by the choice
of normalization and \(r>1/2\).
The
[uniform Gamma moment bound](MULTIVARIATE_COMPARABLE_CORE_BALANCED_TAIL_STABILITY.md#uniform-normalized-endpoint-data)
at \(s=-\sigma\) gives
\(\mathbb E_\mu[(rX)^{-\sigma}]\le M\).
On \(0<s\le\delta<1/2\) we have
\(|(1/s-c)I(s)|\le(1+s)/(1-s)\le3\), by expanding the
geometric sum. Its expectation on this event is
at most \(3M\delta^\sigma=p/4\).
By (2), \(|\mathbb E_\mu H(rX)|<\varepsilon\).
Taking expectations in (5), and bounding the
factor \(1/s-c\) on \(s\ge\delta\), gives
\[
              \mathbb E_\mu |I(rX)|\ge\rho
                           \qquad\text{if }\varepsilon\le p/4.
\]
If also \(\varepsilon\le\rho/2\), then
\(a_+=\mathbb E_\mu I_+(rX)\ge\rho/4\), because
\(|\mathbb E_\mu I(rX)|<\varepsilon\).

Finally use the likelihood-ratio comparison on the
ordered product region \(x>y\) with the nonnegative
test \(I_-(rx)I_+(ry)\), where \(I_\pm\) are the
positive and negative parts of \(I\). Its support
lies in that region by the single crossing. This gives
\[
        a_+\,\mathbb E_{\widehat\nu} I_-(rX)
             \ \ge\ \mathbb E_\mu I_-(rX)\,
                       \mathbb E_{\widehat\nu} I_+(rX).
\]
Consequently, by (2) and (4),
\[
                 b\,\mathbb E_{\widehat\nu} I(rX)
                         \le \frac{24J\,3^J\,\varepsilon}{\rho}.
\]
On the other hand, the imaginary part of (3) after
division by \(w\), together with (2), requires
\[
                 b\,\mathbb E_{\widehat\nu} I(rX)
                                \ge(q-\varepsilon)/r .
\]
The sufficient bound on \(\tau\) makes
\(\varepsilon<\rho/(96J^2\,3^J)\), which implies
all preceding smallness conditions and \(\varepsilon<q/2\).
The upper estimate is below \(1/(4J)\), whereas the
lower estimate is above \(1/(2J)\). This contradiction
excludes a common physical root.

The argument uses the absolute eventual theorem only
to fix \(J\). The lower orders are handled by the
uniform inequalities above. Rational Bézout gives
the assertion over every characteristic-zero field.
