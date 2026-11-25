# An annulus uniform in the allocation of bounded private degree

Let \(a,b\) be distinct ordinary nonnegative integer exponent
columns in any dimension, with equal rows removed. Set
\[
 P_{\rm tot}=\sum_{a_i b_i=0}\max(a_i,b_i),\qquad
 V_{\rm sh}=\sum_{a_i b_i>0}\frac{(b_i-a_i)^2}{\min(a_i,b_i)} .
\]
For any positive integer bound \(P\ge P_{\rm tot}\), put
\(R_P=\max(100,P\max(4,P))\).
At every order \(m\ge1\), all roots of the endpoint-normalized
moment polynomial lie in
\[
        R_P^{-P_{\rm tot}}e^{-V_{\rm sh}/2}
               \ \le\ |z|\ \le\
                R_P^{P_{\rm tot}}e^{V_{\rm sh}/2}.         \tag{1}
\]
For a sequence with total curvature
\(S=\sum_i(b_i-a_i)^2/\max(a_i,b_i)\) bounded, one may take
\(P=\max(1,\lceil\sup S\rceil)\) and use
\(V_{\rm sh}\le S(S+2)\). Thus (1) in particular gives a
fixed normalized annulus for every bounded-curvature
one-private-side family, with arbitrary allocation of its
private degree and arbitrary dimension.

## A single rising factorial of fixed positive base

For \(1/P\le T\le1\), let
\(Q_m(t)=\sum_{k=0}^m\binom mk(T)_k t^k\), with
\((T)_k=T(T+1)\cdots(T+k-1)\) and \((T)_0=1\).
The reciprocal variable \(w=1/t\) gives the polynomial
\[
 B_m(w)=\frac{w^m Q_m(1/w)}{(T)_m}
        =\sum_{k=0}^m \frac{w^k}{k!}\,
          \frac{(-m)_k}{(1-m-T)_k}.                        \tag{2}
\]
This coefficient formula is evaluated only for \(k\le m\),
so it also applies at \(T=1\). The consecutive ratios of
its positive coefficients, writing the coefficient of
\(w^k\) first, are
\[
                 (k+1)\frac{m+T-k-1}{m-k},
                     \qquad 0\le k<m .
\]
They lie in \([T,m]\) because \(0<T\le1\).
Eneström–Kakeya therefore puts all roots in
\(T\le |w|\le m\). We next improve the lower bound to
\(m/100\) for \(m\ge\max(4,P)\).

For \(0<T<1\) the exact exponential transform is
\[
 e^{-w}B_m(w)=\sum_{k\ge0}
        \frac{(-w)^k}{k!}\frac{(1-T)_k}{(1-m-T)_k}.       \tag{3}
\]
To see this directly, for all nonnegative integers \(k\),
the identity
\[
 \sum_{\ell=0}^k (-1)^\ell\binom k\ell
          \frac{(a)_\ell}{(b)_\ell}
           =\frac{(b-a)_k}{(b)_k}
\]
follows for \(0<a<b\) by integrating the binomial expansion
of \((1-x)^k\) under the beta probability density with
parameters \(a,b-a\). The two sides are rational in \(a,b\),
so the identity holds wherever their denominators are
nonzero. Apply it to the Taylor coefficients of the left
side of (3) with \(a=-m,b=1-m-T\). The sum terminates
at \(\ell=m\) on that side, giving (3).

Let \(c=1/100,\ q=\lfloor m/2\rfloor\) and \(|w|\le cm\).
For \(1\le k\le q\), \((1-T)_k\le k!\), and each of the
\(k\) denominator factors in (3) has modulus at least \(m/2\).
The sum of these terms in absolute value is at most
\(2c/(1-2c)=1/49\). For \(q<k\le m\), the first \(q\)
denominator factors have the same bound and the remaining
\(k-q\) have product at least \(T(k-q-1)!\). Hence the
absolute sum in this range is at most
\[
                 T^{-1}(2c)^q c m e^{c m}\le1/100.         \tag{4}
\]
Indeed \(T^{-1}\le m\), \(q\ge m/4\) for \(m\ge4\),
and \(e^{1/100}(1/50)^{1/4}<1/2\), so the expression is
bounded by \(m^2/(100\,2^m)\le1/100\).

For \(k>m\), the first \(m\) denominator factors in (3)
have absolute product \((T)_m\ge T(m-1)!\).
The remaining factor \((1-T)_{k-m}\) cancels in the
numerator and leaves a product of \(m\) factors bounded
above by \(k^m\). Thus the absolute terms in this range
are bounded by
\[
                   v_k=\frac{k^m(cm)^k}{T(m-1)!\,k!}.
\]
For \(k\ge m\), \(v_{k+1}/v_k\le e c<3/100\).
Also \(v_m\le m T^{-1}(e^2 c)^m\) by \(m!\ge(m/e)^m\).
Since \(T^{-1}\le m\), \(e^2c<1/8\), and \(m^2\le2^m\),
the sum for \(k>m\) is at most \(2m^2 8^{-m}\le1/128\).
Together with (4), this gives
\[
       |e^{-w}B_m(w)-1|\le 1/49+1/100+1/128<1
                     \qquad (|w|\le m/100).               \tag{5}
\]
The estimates used only \(T^{-1}\le m\), so they extend to
\(T=1\) by continuity of the finite polynomial \(B_m\) in (2). Thus \(B_m\) has no zeros in the
disk in (5).

For \(m\ge\max(4,P)\), the center
\(G=((T)_m)^{1/m}\) satisfies \(m/e^2\le G\le m\).
For the lower bound use
\((T)_m\ge T(m-1)!\ge m!/(Pm)\),
\(m!\ge(m/e)^m\) and \((Pm)^{1/m}\le e\).
The upper bound uses \(T+k\le m\) for \(k<m\).
Since the root in the normalized variable is \(z=G/w\),
we get \(e^{-2}\le |z|\le100\).
For \(m<\max(4,P)\) the bounds \(T\le G\le m\) and
\(T\le|w|\le m\) give \(1/(Pm)\le |z|\le Pm\). In either case the single
normalized factor and its reciprocal have all root moduli
in \([R_P^{-1},R_P]\).

## Combining shared and private rows

A private row \((0,d)\) has
\[
                (dk)!=d^{dk}\prod_{\ell=1}^{d}(\ell/d)_k .
\]
Indeed both products group the integers \(1,\ldots,dk\)
by their residues modulo \(d\). The geometric factors
\(d^{dk}\) cancel after normalization by the two endpoints.
Thus this row is a binomial-basis coefficient product of
\(d\) normalized rising factorial factors with bases
\(\ell/d\in[1/P,1]\). A private row \((d,0)\) instead uses
the reciprocals of these normalized factors.

A shared row has all normalized root moduli between
\(e^{-V_i/2}\) and \(e^{V_i/2}\), where
\(V_i=(b_i-a_i)^2/\min(a_i,b_i)\), by the
[shared row estimate](MULTIVARIATE_UNIFORM_ROOT_ANNULUS.md#shared-rows-and-large-degree-imbalance).
The full normalized polynomial is the binomial-basis
coefficient product of all these factors. The
[closed-disc composition lemma](GAP_TWO_GRACE_GLOBAL_LOCALIZATION.md#1-a-self-contained-closed-disc-convolution-lemma)
multiplies their upper bounds, and application to all their
reciprocals multiplies their lower bounds. This proves (1)
with exactly \(P_{\rm tot}\) private factors.
Finally, in a shared row write the exponents as \(l,l+d\)
with \(l\ge1\). The row curvature satisfies
\(d^2/(l+d)\ge d/l-1\), hence \((l+d)/l\le S+2\).
Summing the shared-row comparison proves
\(V_{\rm sh}\le S(S+2)\).
