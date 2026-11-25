# Exclusion of the high tail at a quadruple leading root

Work over a field of characteristic zero. On \(2xy+t^2=1\), use
the normalized sphere functional \(\Lambda\), which takes the
angular constant term and integrates in \(t\) over \([-1,1]\) with
measure \(dt/2\).

## Proposition (High-tail exclusion)

Put \(U=1+x\), \(v=Ut\), and fix scalars \(d,a,b,c,e\). Consider
the regular winding-one phase
\[
 f=\frac{U}{x}\left(P_0(v)+x^5(av+b)+cx^6+ex^7\right),
 \qquad P_0(v)=(v^2-1)(v+d).
\]
If \(\Lambda(f^m)=0\) at each of \(m=5,6,7,8,9\), then
\(a=b=c=e=0\). In particular, a phase of least sphere degree
seven in this family is never eventually pure.

### A linear moment block

Use the formal change \(X=x/U,\ T=Ut\) at \(x=0\). It preserves
\(dx/x\wedge dt=dX/X\wedge dT\). The original endpoints become
\(T=\pm(1-X)^{-1}\). Each differs from \(\pm1\) by a series
divisible by \(X\). Along \(T=\pm1+Xv_1\) the transformed phase
is regular in \(X\), because \(P_0(\pm1)=0\). Thus integrating
any power of the phase over either endpoint correction introduces
a factor \(X\) from \(dT=X\,dv_1\) and has angular constant term
zero. This permits use of the original moment functional on
\[
 \frac{P_0(T)}{X}
 +\frac{X^4(aT+b)}{(1-X)^5}
 +\frac{cX^5}{(1-X)^6}
 +\frac{eX^6}{(1-X)^7}.
\]
The first term alone has zero moments at all positive orders.
For \(5\le m\le9\), a summand with at least two factors from
the remaining terms has order in \(X\) at least
\(8-(m-2)=10-m>0\). Only the summands linear in the tail
contribute. For \(n=m-1\), define
\[
 I_n=\frac12\int_{-1}^1P_0(T)^n\,dT,\qquad
 J_n=\frac12\int_{-1}^1T P_0(T)^n\,dT .
\]
Then \(\Lambda(f^m)/m\) is the product of
\[
 \left(\binom n4 J_n,\ \binom n4 I_n,\ \binom n5 I_n,\
       \binom n6 I_n\right)
 \tag{1}
\]
with the column \((a,b,c,e)^t\). We use \(\binom nk=0\) for
\(k>n\). For example the \(c\) entry selects degree \(m-6\)
in \((1-X)^{-6}\), giving \(\binom n5\).

### The rank certificate

Form the matrix in (1) with rows \(n=4,\ldots,8\), in that
order. Let \(A(d)\) be the determinant with the \(n=7\) row
removed and \(B(d)\) the determinant with the \(n=5\) row removed.
To compute them use
\[
 \begin{aligned}
 I_n &= \sum_{\substack{0\le j\le n\\ j\ {\rm even}}}
    \sum_{h=0}^n
    \frac{(-1)^{n-h}\binom nj\binom nh}{2h+j+1}d^{n-j},\\
 J_n &= \sum_{\substack{0\le j\le n\\ j\ {\rm odd}}}
    \sum_{h=0}^n
    \frac{(-1)^{n-h}\binom nj\binom nh}{2h+j+2}d^{n-j}.
 \end{aligned}
 \tag{2}
\]
The denominators in (2) are at most \(25\) in these rows.
For \(z=d^2\), reduction modulo \(29\) gives the following
coefficient lists in increasing order of powers of \(z\):
\[
 \begin{aligned}
 A&:\quad (-10,14,-12,-8,12,-7,5,-4,-5,5,11,-4),\\
 B&:\quad (-12,10,-4,7,-3,-8,-8,8,10,-4,-1,-7,1).
 \end{aligned}
\]
Start the Euclidean algorithm on the monic versions of \(B,A\)
in \(\mathbb F_{29}[z]\), making each nonzero remainder monic
before the next division. The following table records the
remainders before that normalization. The last remainder is the
nonzero constant \(5\).

| Degree | Leading coefficient | Constant coefficient |
| ---: | ---: | ---: |
| 10 | 23 | 24 |
| 9 | 4 | 5 |
| 8 | 24 | 27 |
| 7 | 25 | 23 |
| 6 | 16 | 17 |
| 5 | 12 | 10 |
| 4 | 19 | 6 |
| 3 | 26 | 10 |
| 2 | 14 | 1 |
| 1 | 23 | 25 |
| 0 | 5 | 5 |

Coprimality persists under substituting \(z=d^2\) by substituting
in a Bézout relation. Since \(\deg I_n\le n\) and
\(\deg J_n\le n-1\), the degrees of \(A,B\) in \(d\) are at
most \(22,24\), respectively. The reductions attain those bounds.
Thus the resultant of \(A,B\) in \(d\) is nonzero modulo \(29\),
so it is nonzero over \(\mathbb Q\). The two minors cannot vanish
together at any \(d\) over a characteristic-zero field. Therefore
the five moment equations force \(a=b=c=e=0\).

If all four coefficients vanish, the term of negative angular
exponent is \(x^{-1}P_0(t)\). It reduces to \(-2y(t+d)\)
on the sphere, and every other term has total degree at most six.
So the phase has sphere degree at most six. In least degree seven,
eventual purity would imply the
vanishing of every positive moment by the
[quadruple-root boundary corollary](WINDING_ONE_CUBIC_LOGARITHMIC_BOUNDARY.md#consequences-in-degree-seven).
That would force all four coefficients to vanish, which contradicts
the degree-seven hypothesis.
