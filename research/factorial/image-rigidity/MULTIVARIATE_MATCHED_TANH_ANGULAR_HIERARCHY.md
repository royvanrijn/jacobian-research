# Matched curvature bounds at the angular comparison scale

Let \(n\ge2\), \(N=n+1\), and let \(a,b\) be distinct ordinary
nonnegative exponent columns. Omit equal rows and write
\(d_i=b_i-a_i\), \(m_i=(a_i+b_i)/2\), and
\(\delta=\max_i|d_i|/m_i\). Let \(L_j(s)\) be the endpoint-normalized
log-factorial symbol per degree, as in the
[matched midpoint computation](MULTIVARIATE_MATCHED_TANH_ORDINARY_ANGLE_LIFT.md#exact-curvature-of-the-assembled-symbol).
The exact factorial polynomial, up to a scalar and positive change
of variable, has coefficient at \(k\) equal to
\(\binom jk e^{jL_j(k/j)}\).

Put
\[
 K(s)=L_n''(s),\qquad G(s)=L_N''(s)-L_n''(s),\qquad
 C=K(1/2),\quad J=G(1/2).
\]
Both \(C\) and \(J\) are strictly positive. For every integer \(p\ge0\),
\[
 \frac{|K^{(p)}(1/2)|}{p!}\le C\delta^p,\qquad
 \frac{|G^{(p)}(1/2)|}{p!}\le (p+1)J\delta^p.                 \tag{1}
\]
The constants \(C\) and \(J\) are the exact midpoint curvatures of
the lower symbol and of the matched increment. The second bound
controls the shape of the increment relative to its own size.

For \(q=\delta |w|<1\), subtracting the constant and linear terms
from these two curvatures gives
\[
 \begin{aligned}
 |K(1/2+w)-C-wK'(1/2)|&\le \frac{C q^2}{1-q},\\
 |G(1/2+w)-J-wG'(1/2)|&\le \frac{J q^2(3-2q)}{(1-q)^2}.
 \end{aligned}                                                \tag{2}
\]
If \(K^+\) and \(G^+\) denote the averages under reflection about
\(1/2\), summing only even powers gives
\[
 |K^+(1/2+w)-C|\le\frac{C q^2}{1-q^2},\qquad
 |G^+(1/2+w)-J|\le\frac{J q^2(3-q^2)}{(1-q^2)^2}.              \tag{3}
\]

## Proof of the curvature hierarchy

The [exact matched kernel](MULTIVARIATE_HOLOMORPHIC_MATCHED_LOG_TRANSPORT.md#the-exact-curvature-kernel)
represents the row contributions to \(K\) and \(G\) by
\[
 \begin{aligned}
 K_i(1/2+w)&=d_i^2\int_0^\infty e^{-(m_i+d_iw)u}h(u/n)\,du,\\
 G_i(1/2+w)&=d_i^2\int_0^\infty e^{-(m_i+d_iw)u}q_n(u)\,du,\\
 h(t)&=\frac{t}{e^t-1},\qquad q_n(u)=h(u/N)-h(u/n).
 \end{aligned}
\]
The integrals converge and can be differentiated on compact subsets
of \(|w|\delta<1\). The function \(h\) is positive and strictly
decreasing for \(t>0\). It is strictly convex because
\[
 h''(t)=\frac{e^t}{(e^t-1)^2}\bigl(t\coth(t/2)-2\bigr)>0,
\]
where \(\tanh(t/2)<t/2\). Consequently \(q_n(u)>0\) and
\[
 \frac{q_n(u)}u=\int_{1/N}^{1/n}-h'(vu)\,dv
\]
is decreasing in \(u>0\).

If a probability law on \(u\ge0\) is tilted by a positive decreasing
weight \(b(u)\), its expectation of \(u^p\) cannot increase.
For independent original samples \(U,V\), the sign follows by
integrating
\((U^p-V^p)(b(U)-b(V))\le0\).
For \(K_i\), tilt the exponential law with rate \(m_i\) by \(h(u/n)\).
For \(G_i\), tilt the Gamma law with shape two and rate \(m_i\) by
\(q_n(u)/u\). The original \(p\)-th moments are \(p!/m_i^p\) and
\((p+1)!/m_i^p\). Differentiating in \(w\) multiplies the integrands
by \((-d_i u)^p\). Divide by the respective row masses, apply the
moment bounds, and sum to obtain (1). Summing the corresponding
geometric series gives (2) and (3).

For a quantitative comparison put \(S_*=\sum_i d_i^2/m_i\).
The inequality \(h(t)\ge e^{-t}\) gives
\[
 C\ \ge\ \sum_i\frac{d_i^2}{m_i+1/n}\ \ge\ S_*/2,\qquad
 \delta\le S_* .
\]
Here \(m_i\ge1/2\), \(n\ge2\), and \(|d_i|\ge1\). For fully
shared ordinary rows \(m_i\ge3/2\), so \(S_*\le4C/3\). Also
\(-h'(t)\le1/2\) follows from convexity and \(h'(0)=-1/2\). Thus
\[
 J\le\frac1{2nN}\sum_i\frac{d_i^2}{m_i^2}
    \le\frac{S_*\delta}{2nN}.
\]
In particular, for weak curvature at a fixed \(w\), the reflected
shape corrections in (3) are \(O(C^3)\) for the lower symbol and
\(O(C^4/(nN))\) for the increment. The estimates are uniform in the
number of rows and include private rows. Both leading matched
midpoint clocks lie in the physical range
\[
       \frac{n}{N}<
          \frac{(C+J)/(2N)}{C/(2n)}<1 .
\]
The lower inequality uses \(J>0\). The upper follows from the strict
decrease of \(\psi'(1+jm_i)\) in \(j\), since
\(L_j''(1/2)/j=\sum_i d_i^2\psi'(1+jm_i)\).

## Using the midpoint hierarchy in a paired phase lift

Take the common Gaussian strength \(S=C\) and put
\(E_j(s)=j[L_j(s)+C s(1-s)/2]\).
On \(|s|\le T\), set \(r=T+1/2\) and \(q=\delta r\). When \(q<1\),
the hypotheses on the correction in a
[paired phase lift](MULTIVARIATE_DIFFUSE_ABSTRACT_PAIRED_HALF_LIFT.md)
hold with the curvature tolerance
\[
 \max\{\|E_n''\|/n,\ \|E_N''\|/N,\ \|E_N''-E_n''\|\}
 \ \le\
 r\,|L_n'''(1/2)|+\frac{C q^2}{1-q}+\frac{N J}{(1-q)^2}.       \tag{4}
\]
Indeed \(E_N''-E_n''=K-C+NG\), and (1) bounds \(\|G\|\)
by \(J/(1-q)^2\). For a mirrored profile the first term vanishes
and \(C q^2/(1-q^2)\) suffices in place of \(C q^2/(1-q)\).
This criterion uses the exact skew after orientation cancellation
and the exact matched curvature increment.

The rational matched tanh test also extends to the explicit
Gaussian strip \(16/3\le S\le6\). With the lower-order parameters
\(r_i=nT_i/(n-1)\), let \(A=\sum r_i\),
\(P=\sum\epsilon_i r_i^2\), and \(V_p=\sum r_i^p\), with indices
now over Gauss factors. For \(n\ge10^7\), sufficient conditions for
disjoint arguments are
\[
 \frac{16}3\le 2A\le 6,\qquad
        262154|P|+46000000000 V_3+\frac{4V_2}{n}\le10^{-18}.    \tag{5}
\]
Use the [explicit paired lift](MULTIVARIATE_DIFFUSE_ABSTRACT_PAIRED_HALF_LIFT.md#explicit-application-down-to-curvature-163)
at \(T=65538\) with Gaussian strength \(S=2A\). Its analytic
tolerance is \(10^{-18}\).
On \(|s|\le T\) one has \(|2s-1|\le131077\).
The bound on \(V_3\) in (5) implies each \(r_i<10^{-9}\), and
the matched upper parameters are smaller than \(2r_i\).
Thus \(\Re(1+\epsilon_i r(2s-1))\ge3/4\) along every matched
parameter segment. Write \(d=3/4\). Apply the
[midpoint reciprocal expansion and grid drift estimates](MULTIVARIATE_MATCHED_TANH_ORDINARY_ANGLE_LIFT.md#size-and-adjacent-change-of-the-log-coefficient-correction)
with this disk size. The signed term has coefficient \(2\cdot131077\).
The unsigned third-moment coefficient is at most
\[
 \frac{2\cdot131077^2}{d}
       +\frac{2}{3d^3}\left(\frac14+\frac5{12}\right)
       <46000000000 .
\]
The adjacent parameter drift has coefficient at most \(4/n\)
on \(V_2\). These bounds give the hypotheses of the explicit lift.

For example take two reflected rows with gap \(D=10^{16}\) and
smaller exponent \(\ell=9D^2/25\). Then
\(D/(2\ell+2D)\le r_i\le D/(2\ell)\) for each of the \(2D\) factors.
They give \(P=0\), \(2D^2/(\ell+D)\le2A\le50/9\),
\(V_2\le625/(162D)\), and \(V_3\le15625/(2916D^2)\).
These inequalities verify (5) for every \(n\ge10^7\).

## An exact collision for an unrestricted common Laplace symbol

Positive Laplace curvature alone allows an exact common root at
consecutive degrees, even with the same analytic symbol at both
orders and simple unit-circle roots at each order. For \(c>0\) set
\[
 h_0=4\,\operatorname{arcosh}2,\quad
 D(s)=\cosh(h_0(s-1/2))-\cosh(h_0/2),\quad
 L(s)=-c s(1-s)+\alpha c^2D(s).
\]
This is an endpoint-normalized real symmetric symbol with positive
Laplace curvature when \(\alpha>0\). Use the samples
\(\binom jk e^{jL(k/j)}\) at \(j=3,4\).
In the variable \(x=z+z^{-1}\), the cubic after removal of its
factor \(1+z\) gives \(x=x_3=1-3e^{3L(1/3)}=-2+2c+O(c^2)\).
For the quartic, division by \(z^2\) gives
\(Q(x)=x^2+4u x+6v-2\). Since \(D(1/4)=-5\) and \(D(1/2)=-6\),
\[
 u=e^{-3c/4-20\alpha c^2},\qquad
 v=e^{-c-24\alpha c^2},\qquad
 4u^2-6v+2=(3/2-16\alpha)c^2+O(c^3).
\]
Using \(Q(x_3)=(x_3+2u)^2-(4u^2-6v+2)\) shows that
\(c^{-2}Q(x_3)\) extends real analytically at \(c=0\) with value
\(16\alpha-5/4\). Its \(\alpha\) derivative there is \(16\).
The implicit function theorem gives \(\alpha(c)=5/64+O(c)>0\)
for all sufficiently small positive \(c\), with \(Q(x_3)=0\).
The two quartic roots in \(x\) are
\(-2+c+O(c^2)\) and \(-2+2c+O(c^2)\). They and \(x_3\)
lie in \((-2,2)\), so the two polynomials have simple unit-circle
roots and share the conjugate pair determined by \(x_3\).

This example also identifies a factorial constraint visible in just
the quartic samples. For any ordinary \(L_j\), let \(L_j^+\) denote
its reflection average and let \(C=L_j''(1/2)\) at that order.
The proof of (1) for \(K\) applies with \(n=j\).
On \([0,1]\) the reflected curvature is at least \(C\) by the
positive Laplace integral with the factor \(\cosh(d_iwu)\).
Therefore \(-L_j(1/2)\ge C/8\).
From (3) and the zero endpoint values it follows that
\[
 L_j^+(1/4)-\tfrac34 L_j(1/2)=O(C^3)
 \qquad (C\longrightarrow0).
\]
Indeed, after subtracting the endpoint-normalized quadratic of
curvature \(C\), the absolute second derivative on \([0,1]\) is
\(O(C^3)\) since \(\delta\le2C\), and integrating the endpoint
Green kernel gives a bound of half this norm times \(s(1-s)\).
For the colliding symbol the same sample difference is
\(-\alpha(c)c^2/2\), whereas the midpoint sample would force
\(C\le-8 L(1/2)=2c+O(c^2)\).
Thus the colliding samples violate a quantitative ordinary
factorial constraint for sufficiently small \(c\). This obstruction
applies to the relaxed common-symbol class.

The estimates above give analytic control of the exact linked symbol
and an angular criterion in a fixed curvature strip. A general
all-degree comparison needs a gap estimate between consecutive
degrees in the physical clock range and a transfer that preserves
those gaps for exact factorial symbols.
