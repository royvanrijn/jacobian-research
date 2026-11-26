# Matched tanh angular separation near the private-unit core and all integer twists for one shared unit gap

Take two ordinary exponent columns. Suppose their only private rows
are a unit row on each side. Shared rows can have either
orientation. At adjacent orders \(n,N=n+1\), put
\[
                 A=\sum_{\rm shared\ rows}\frac{(b_i-a_i)^2}
                                                 {\min(a_i,b_i)}.
\]
If \(n\ge2\) and \(A\le1/(nN)\), the two factorial moments have
disjoint root arguments under independent positive rescalings.
This condition can be certified from a positive first matched
tanh moment. With an absolute terminal cutoff it gives an
all-degree angular criterion.

For the core with precisely one unequal shared row \((a,a+1)\),
\(a\ge1\), there is a stronger all-degree coprimality statement.
Write its ordinary moments as \(Q_m\). Then for every integer
\(h\) and every \(n\ge1\) the polynomials
\[
                         Q_n(z),\qquad Q_N((N/n)^h z)
\]
are coprime. Consequently the
[small-tail stability theorem](MULTIVARIATE_ARBITRARY_CORE_TWISTED_TAIL_STABILITY.md)
applies to this core for every integer signed gap of the appended
shared tail. The proof combines a geometric root comparison with
a width-one Newton polygon at a prime dividing one of the degrees.

## The geometric angular estimate

The two private factorials cancel the binomial coefficient, up to
a scalar \(m!\) in degree \(m\). After endpoint normalization the
coefficients \(c_k\), \(0\le k\le m\), therefore satisfy
\(c_0=c_m=1\), \(0<c_k\le1\), and \(\log c_k\) is convex in \(k\).
More precisely its continuous extension \(\ell(x)\), \(0\le x\le m\),
satisfies \(\ell''(x)\le A/m\). To see this, for a shared row use
\[
 (\log\Gamma(ma_i+x(b_i-a_i)+1))''
       =(b_i-a_i)^2\psi'(ma_i+x(b_i-a_i)+1)
       \le\frac{(b_i-a_i)^2}{m\min(a_i,b_i)}.
\]
The last inequality follows from the integral comparison of
the positive series for \(\psi'\). Normalizing the endpoints
subtracts an affine function. The endpoint chord bound gives
\(-\ell(k)\le A k(m-k)/(2m)\le Am/8\), so
\[
          \sum_{k=1}^m |c_k-c_{k-1}|
                        =2(1-\min_k c_k)\le Am/4 .              \tag{1}
\]
The equality uses the decrease and then increase of \(c_k\),
with a possible interval of minima.

Let \(F_m(z)=\sum_{k=0}^m c_k z^k\) and \(r=1/(nN)\). Then
\[
      (1-z)F_m(z)=1-z^{m+1}
                          +\sum_{k=1}^m(c_k-c_{k-1})z^k.
\]
If \(\zeta^{m+1}=1\), on the boundary of \(|z-\zeta|\le r\)
put \(\alpha=(m+1)r\). For \(m=n,N\) we have
\(\alpha\le2/3\). Expanding \(z^{m+1}=(1+(z/\zeta-1))^{m+1}\)
and using the triangle inequality gives
\[
 |1-z^{m+1}|\ \ge\ \alpha-
       \sum_{k\ge2}\frac{\alpha^k}{k!}
 \ \ge\ \alpha-\frac{\alpha^2}{2(1-\alpha/3)}
 \ >\ \frac{\alpha}{2}.
\]
The last sum bound uses \(k!\ge2\cdot3^{k-2}\) for \(k\ge2\).
By (1) the perturbation is bounded by
\((Am/4)(1+r)^m\le Am/2\le mr/2\), since \(mr\le1/n\) and
\(e^{1/2}<2\). Rouché's theorem gives one root per disk.
The disks are disjoint within each degree because consecutive
\((m+1)\)-th roots of unity are at distance
\(2\sin(\pi/(m+1))\ge4/(m+1)>2r\).
The disk at \(1\) contains the extra root supplied by \(1-z\).
All other disks each contain one root of \(F_m\).

A disk of radius \(r<1\) centered on the unit circle changes the
argument by at most \(\arcsin r<\pi r/2\). The reference arguments
of the nontrivial \((n+1)\)-th and \((n+2)\)-th roots of unity
have circular distance at least \(2\pi/((n+1)(n+2))\).
The two angular errors total less than \(\pi/(nN)\), which is
at most this reference distance for \(n\ge2\). This proves the
claimed angle separation.

## A first matched moment criterion for every degree

For each shared Gauss factor recover the midpoint rate from its
matched tanh values at any fixed reference lower order \(j\ge2\):
\[
             v_i=\left(\frac{j}{T_i'}-
                               \frac{j-1}{T_i}\right)^{-1}.
\]
On a row of smaller exponent \(\ell\ge1\) and gap \(d\) the
[exact Gauss grid](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid)
has \(d\) factors, each with \(v_i=d/(2\ell+d)\). Thus
\[
                  A=\sum_{\rm shared\ factors}
                                  \frac{2v_i}{1-v_i}.
\]
Only the shared factors are summed. The private units have
already been reserved. If \(M_1=\sum_{\rm shared\ factors}T_i\)
and \(c=(j+1)/(j-1)\), the same grid gives \(v_i\le cT_i\).
When \(cM_1<1\), therefore, \(A\le2cM_1/(1-cM_1)\).
Take \(J\ge4\) at least the
[computable absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
The single condition
\[
                      M_1\le
              \frac{1}{c\,[2(J-2)(J-1)+1]}                     \tag{2}
\]
ensures \(A\le1/((J-2)(J-1))\). The geometric estimate then
covers every \(2\le n<N<J\). The first and second moments
have disjoint root arguments since the former has one negative
root and the latter is strictly positive on the real line.
The absolute theorem covers \(N\ge J\). Inequality (2)
places no restriction on shared orientations or row count.
The cutoff is unevaluated. This argument uses the matched
first moment to bound the change of the entire geometric
polynomial.

## The Newton polygon for the single shared unit-gap row

For this core, divide \(Q_m\) by its constant coefficient and
take the reciprocal polynomial:
\[
         R_m(x)=x^m Q_m(1/x)/Q_m(0)
                    =\sum_{k=0}^m (ma+1)_k\,x^{m-k}.
\]
This polynomial is monic integral. Twists \(h\ge0\) are covered
by the
[comparable-core twisted comparison](MULTIVARIATE_ARBITRARY_CORE_TWISTED_TAIL_STABILITY.md#the-exceptional-twists-for-comparable-both-private-cores).
For a twist \(h=-b\), \(b\ge1\), a common root \(z\) would make
both \(1/z\) and \((N/n)^b/z\) algebraic integers. Since
\(\gcd(n^b,N^b)=1\), this makes their common scaling
\(y=1/(n^b z)\) integral as well. In particular the corresponding
reciprocal roots at orders \(m=n,N\) are \(m^b y\).
They are nonreal, since the even moment has no real root.

The geometric estimate covers \(a\ge nN\), since here \(A=1/a\).
Suppose \(a<nN\). Choose the largest prime-power divisor
\(q=p^s\) of \(nN\). It divides exactly one of \(n,N\), call it
\(m\), and write \(m=q u\), \(p\nmid u\). Put \(d_j=[x^j]R_m\)
and \(f(j)=v_p(d_j)+s j\). Directly from the coefficient ratios,
\[
        f(j+1)-f(j)=s-v_p(m(a+1)-j),\qquad 0\le j<m .
\]
At \(j=0\) this is nonpositive. If \(q\nmid j\), then
\(v_p(m(a+1)-j)=v_p(j)<s\), so \(f(j)\) strictly increases.
Besides \(0\), the only possible minimum indices are \(tq+1\)
for \(0\le t<u\). For consecutive such indices,
\[
 f(tq+1)-f((t-1)q+1)
      =C-v_p(u(a+1)-t),\qquad C=sq-\frac{q-1}{p-1}.           \tag{3}
\]
Indeed a full block contributes \(s q\), minus
\(\sum_{j=1}^{q-1}v_p(j)=(q-1)/(p-1)-s\) and the valuation
\(s+v_p(u(a+1)-t)\) of its final factor. Hence if
\(u(a+1)<p^C\), every increment in (3) is positive.
The last minimizer of \(f\) is then \(1\). The Newton polygon of
\(R_m(p^s x)\) can have at most one root of nonnegative valuation,
counted with multiplicity. Equivalently \(R_m(x)\) has at most one
root of valuation at least \(s\). But every conjugate of \(m^b y\)
has valuation at least \(bs\ge s\), and a nonreal algebraic integer
has at least two conjugates. This gives the contradiction.

## The elementary bound for the largest prime power

Put \(L(q)=\operatorname{lcm}(1,\ldots,q)\). We have
\(nN\mid L(q)\) by the choice of \(q\), and for \(n\ge2\),
\(N\le\sqrt{3 L(q)/2}\). Since \(a<nN\),
\[
            u(a+1)\le u nN
                      \le \frac{\sqrt{3L(q)^3/2}}{q}.
\]
We claim \(q p^C>\sqrt{3L(q)^3/2}\) for every prime power \(q=p^s\).
For \(q\ge18\), first note that
\(\log L(q)\le 2q\log 2+\log q\). For an integer \(t\ge2\),
each prime power in \((\lceil t/2\rceil,t]\) contributes its prime
as a factor of \(\binom{t}{\lfloor t/2\rfloor}\), by the factorial
valuation formula. Each such prime occurs at most once in this
interval. Hence
\(\log L(t)-\log L(\lceil t/2\rceil)\le t\log 2\).
Iterate with \(t_i=\lceil q/2^i\rceil\) until \(i=\lceil\log_2 q\rceil\).
The sum of the \(t_i\) before the final value \(1\) is at most
\(2(q-1)+\lceil\log_2 q\rceil\le2q+\log_2 q\).
This proves the bound. Also
\(\log p/(p-1)\le\log2\), so
\[
  \log\!\left(\frac{q p^C}{\sqrt{3L(q)^3/2}}\right)
       \ge q\log(q/16)+\log2-\tfrac12\log q-\tfrac12\log(3/2)>0.
\]
The last expression is increasing for \(q\ge18\).
At \(q=18\), use \(\log(9/8)\ge1/9\) and \(e^4>27/4\).
For prime powers below \(18\), the integers in the last two columns
verify \(2(qp^C)^2>3L(q)^3\) by squaring and cubing.

| \(q\) | \(L(q)\) | \(qp^C\) |
|---:|---:|---:|
| 2 | 2 | 4 |
| 3 | 6 | 27 |
| 4 | 12 | 128 |
| 5 | 60 | 3125 |
| 7 | 420 | 823543 |
| 8 | 840 | 1048576 |
| 9 | 2520 | 43046721 |
| 11 | 27720 | 285311670611 |
| 13 | 360360 | 302875106592253 |
| 16 | 720720 | 9007199254740992 |
| 17 | 12252240 | 827240261886336764177 |

This gives the strict bound needed in (3). The first pair
\(n=1,N=2\) is separate because real and nonreal roots cannot
match under a positive scaling. All conclusions on coprimality
descend to every characteristic-zero field, since the twists and
the moment coefficients are rational.
