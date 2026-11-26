# First-sector concentration above shared or private-endpoint criticality

Let \(P_n(z)=\sum_k p_{n,k}z^k\) be the ordinary factorial moment
polynomial of two exponent columns, and for \(r>0\) let \(K\) have
weights proportional to \(p_{n,k}r^k\). Equal rows may be deleted.
The estimates below hold uniformly over \(r\) and do not use any
restriction on the signed total degree.

**Fully shared profiles.** Put
\[
 \rho(t)=t(1-t)\sum_i\frac{(b_i-a_i)^2}{(1-t)a_i+tb_i},
 \qquad \epsilon=\max_{[0,1]}\rho-1 .
\]
If \(0\le\epsilon\le1/100\) and \(n\ge1024\), then
\[
 \operatorname{Var}(K)\le
           \frac{16}{3}\epsilon n^2+\frac83n^{3/2}+33n .        \tag{1}
\]
In particular the closed sector \(|\arg z|\le2\pi/(n+1)\) is
zero-free if \(\epsilon\le1/200\) and \(n\ge16384\), or if
\(\epsilon\le1/125\) and \(n\ge262144\).
The [subcritical theorem](MULTIVARIATE_SUBCRITICAL_FIRST_SECTOR_CONCENTRATION.md)
covers \(\max\rho\le1\) with a sharper first-sector bound.

**One private unit.** Suppose the only private row is \((1,0)\).
Put \(\gamma=\sum_{\rm shared}(b_i-a_i)^2/b_i\).
If \(\gamma\le17/16\), set \(\epsilon=\max(0,\gamma-1)\)
and \(c_0=1/2-\epsilon\). For \(n\ge64\) we have
\[
 \operatorname{Var}(K)
 \le 2 \max\left\{
      \frac{(8+\epsilon n+n^{2/3}/2)^2}{c_0^2},\
      \frac98 n^{4/3}\right\}.                                \tag{2}
\]
When \(\gamma\le1\) this gives \(\operatorname{Var}(K)\le8n^{4/3}\)
at every \(n\ge1\), using the trivial variance bound at \(n<64\).
Consequently the sector \(|\arg z|\le1/(\sqrt{8}\,n^{2/3})\)
is zero-free. The full closed first sector is zero-free if
\(\epsilon\le 1/25\) and \(n\ge8192\), or throughout
\(\gamma\le17/16\) if \(n\ge524288\).
For \(\gamma\le1\), the earlier log-concavity test already excludes
the closed first sector at every \(n\ge2\).

In both cases the zero-free conclusions follow by applying
\[
 \Re\mathbb E e^{i\theta(K-\mathbb E K)}
       \ge 1-\theta^2\operatorname{Var}(K)/2,\qquad
 P_n(re^{i\theta})/P_n(r)=\mathbb E e^{i\theta K}.
\]
For the first sector it suffices that \(\operatorname{Var}(K)<n^2/20\),
since \(\pi^2<10\).

## Conditional gap comparison

Let \(F(x)\) be the continuous Gamma extension of
\(\log(p_{n,x}r^x)\) on \([0,n]\). With \(K,L\) independent
and identically distributed, condition on \(q=(K+L)/2\) and set
\(d=|K-L|/2\). The logarithmic gap weights, apart from the
orientation factor of two for \(d>0\), are \(S_q(d)=F(q+d)+F(q-d)\).
As in the [critical comparison](MULTIVARIATE_SUBCRITICAL_FIRST_SECTOR_CONCENTRATION.md#a-sharper-first-sector-threshold),
an upper bound for \(S_q'\) gives likelihood-ratio domination by
integrating that bound and extending the dominating weights over
the full integer or half-integer grid. The orientation factor
is the same in both laws. The decreasing likelihood ratio
implies an upper bound on the conditional second moment. At
the end \(\operatorname{Var}(K)=2\mathbb E[d^2]\).
We also use the centered-cell quadrature estimate from that proof,
with error at most half the grid spacing times total variation.

## The fully shared curvature bound

Write \(u_i(t)=a_i+d_i t\), \(d_i=b_i-a_i\), and
\(W(t)=\sum_i d_i^2/u_i(t)\). Let \(w\) maximize \(\rho\).
A single row has maximum \((\sqrt{a_i}-\sqrt{b_i})^2\), so
\(\max\rho\le101/100\) forces each endpoint ratio to be in
\([4/17,17/4]\). Indeed the smaller endpoint is an integer
at least one and \((1+\sqrt{101/100})^2<17/4\).
The row maxima lie at \(\sqrt{a_i}/(\sqrt{a_i}+\sqrt{b_i})\),
so \(w\in[13/40,27/40]\).

For \(p_\lambda(t)=t(1-t)(1+(\lambda-1)t)\) the row bound
\[
                  p_\lambda(t)\le \sqrt{\lambda/8}
            \qquad(4/17\le\lambda\le17/4)
\]
follows as follows. Reversal reduces to \(\lambda\ge1\).
For each \(t\), \(p_\lambda(t)/\sqrt\lambda\) is convex
as a function of \(\sqrt\lambda\), so test the two endpoints.
For \(\lambda=17/4\) the polynomial increases up to \(1/2\)
and is concave thereafter. Its tangent at \(5/8\) gives
\(p_{17/4}(t)\le1458/2048<5/7\) on \([1/2,1]\), and the
increasing interval has no larger values. Since
\((5/7)^2<17/32\), the endpoint passes.

The [exact tangent remainder](MULTIVARIATE_SUBCRITICAL_FIRST_SECTOR_CONCENTRATION.md#critical-curvature-and-a-continuous-coefficient-bound)
of row \(i\), divided by \(t(1-t)(t-w)^2\) and by its
contribution \(\rho_i(w)=d_i^2w(1-w)/u_i(w)\), is
\[
           \frac{a_i b_i}{w(1-w)u_i(w)t(1-t)u_i(t)}\ge8 .
\]
The linear terms cancel at \(w\). Thus
\[
      \frac{1-\rho(t)}{t(1-t)}
                     \ge -\frac{\epsilon}{t(1-t)}+8(t-w)^2 .
\]
On \(1/8\le t=x/n\le7/8\) the trigamma bounds of the critical
proof give \(nF''(x)\le-(1-\rho(t))/(t(1-t))+33/n\).
For \(1/4\le t\le3/4\) the excess term is at most
\(16\epsilon/3\). For \(t\le1/4\) in the remaining interval
its excess over \(16\epsilon/3\) is absorbed by \(2(t-w)^2\).
Indeed on each interval \([u,v]\) with successive endpoints
\(1/8,1/6,7/32,1/4\) one has the rational comparisons
\[
       \frac1{100}\left(\frac1{u(1-u)}-\frac{16}{3}\right)
                    \le 2(13/40-v)^2.
\]
The right-hand intervals follow by reversal.

For \(t\le1/8\), comparison to \(W(1/3)\) gives
\(W(t)\le (75/8)(1+\epsilon)/(1+13t/4)<(19/2)/(1+13t/4)\).
The binomial curvature times \(n\) is at least
\(1/(t+1/1024)+1024/1025\).
The difference \(1/(t+1/1024)-(19/2)/(1+13t/4)\)
is decreasing there, since the ratio between its
denominators in reverse order is at least
\(1440/129>\sqrt{247/8}\). Testing at \(t=1/10\) and \(1/8\)
gives on the respective intervals \([0,1/10]\) and \([1/10,1/8]\)
a margin larger than \(6(27/40)^2\) and \(6(23/40)^2\).
These bound \(6(t-w)^2\) on those intervals. Other endpoints
again follow by reversal. We conclude
\[
          F''(x)\le \frac{33}{n^2}+\frac{16\epsilon}{3n}
                                -\frac{6(x-nw)^2}{n^3}.        \tag{3}
\]

The integral of \((x-nw)^2\) over \([q-y,q+y]\) is at
least \(2y^3/3\). Integrating the resulting bound on
\(S_q'\), the conditional law is dominated by the full-grid
weight \(e^{-n(t^2-\mu)^2}\), where
\(t=d/n\) and \(2\mu=33/n+16\epsilon/3\).
The continuum weight \(f_0(t)=e^{-n(t^2-\mu)^2}\) has
mean \(t^2\) at most \(\mu+1/(2\sqrt n)\). To see this,
integrate \((t f_0(t))'\) and use Cauchy–Schwarz to get
\((\mathbb E t^2)^2-\mu\mathbb E t^2\le1/(4n)\).
Also \(\mu<1/10\) and the continuum mass satisfies
\[
 I_0\ge 2e^{-1}\bigl(\sqrt{\mu+1/\sqrt n}-\sqrt\mu\bigr)
                         >2e^{-1}/\sqrt n .
\]
The variation of \(f_0\) is at most four. The function
\(t^2 f_0(t)\) has one maximum on each positive or negative
half-line and its maximum is at most \(\mu+1/\sqrt n\).
The latter follows by writing \(u=t^2\) and bounding
\((u-\mu)e^{-n(u-\mu)^2}\) above by \(1/\sqrt n\)
when \(u\ge\mu\). Its total variation is at most four times
that upper bound. The quadrature estimate with spacing \(1/n\)
therefore gives for the discrete mean
\[
   \mathbb E t^2\le \mu+
   \frac{1/2+6\mu+3/\sqrt n}{(1-3/\sqrt n)\sqrt n}
                      \le\mu+\frac4{3\sqrt n},
\]
using \(e<3\), \(n\ge1024\), and \(\mu\le1/10\).
This proves (1) after averaging over \(q\). Its two stated
first-sector ranges follow by substitution.

## The one-private endpoint bound

The private row cancels the \((n-k)!\) denominator.
Only \(1/k!\) contributes negative binomial log curvature.
Put \(W(t)=\sum_{\rm shared}d_i^2/u_i(t)\). The bound on
\(\gamma\) implies \(d_i^2\le(17/16)b_i\), and integrality
gives \(1/2\le b_i/a_i\le2\). For an increasing gap of size
\(e\), this follows from \(b_i\ge2e\), which holds directly
for \(e=1\), by rounding for \(e=2\), and from the square
bound for \(e\ge3\). For a decreasing gap the square bound
and integrality give \(e\le b_i\).
Consequently
\[
          \gamma-tW(t)=(1-t)\sum_{\rm shared}
                 \frac{d_i^2 a_i}{b_i u_i(t)}
                            \ge \frac{\gamma}{2}(1-t).
\]
Thus \((1-tW(t))/t\ge-\epsilon+c_0(1-t)\) for \(0<t\le1\).
For \(t=x/n\ge1/4\), apply the trigamma bounds to \(F\) to get
\[
                    F''(x)\le\frac8{n^2}+\frac{\epsilon}{n}
                                        -\frac{c_0(n-x)}{n^2}.
                                                               \tag{4}
\]
For \(t\le1/4\) this follows directly because
\(W(t)\le2\gamma\le17/8\), whereas the negative \(1/k!\)
curvature times \(n\) is at least \(64/17\) for \(n\ge64\).

In the conditional comparison write \(m=n-q\), so \(d\le m\).
When \(c_0m\le 8+\epsilon n+n^{2/3}/2\), this bound on \(d\)
gives the first term of (2). Otherwise (4) gives
\(S_q'(y)\le-y n^{-4/3}\). The dominating full-grid weight is
then \(e^{-t^2/2}\) on the parity grid scaled by
\(t=d/n^{2/3}\). At grid spacing \(h\le1/16\), its centered
quadrature errors for mass and second moment are at most
\(h\) and \(2h\), since the total variations are two and
\(8/e<4\), respectively. The continuum mass and second
moment coincide and exceed two. The discrete second moment
is at most \((2+2h)/(2-h)\le34/31<9/8\).
This proves (2), and its first-sector ranges follow by
substitution, using \(8192^{1/3}>20\) and
\(524288^{1/3}>80\).
