# Binary variance for prescribed-prime semigroups in power-length windows

Fix \(q\ge3\). Let \(C_q\) consist of \(1\) and the positive integers whose
prime factors are all \(1\bmod q\). Write \(f=\mathbf1_{C_q}\),
\(L=\operatorname{lcm}(2,q)\), \(\alpha=1/\varphi(L)\), \(Z=\log X\) and
\(\mu_X=c_q Z^{\alpha-1}\), where
\[
c_q=\frac{1}{\Gamma(\alpha)}
\lim_{s\downarrow1}(s-1)^\alpha\prod_{p\equiv1\bmod q}(1-p^{-s})^{-1}>0.
\]
All supported integers are \(1\bmod L\). Write
\(\mathcal B=\{p:p\nmid L,\ p\not\equiv1\bmod q\}\) and define
\[
S_q(N)=L\mathbf1_{N\equiv2\bmod L}\prod_{p\in\mathcal B}
\begin{cases}
1+1/(p-1),&p\mid N,\\
1-1/(p-1)^2,&p\nmid N.
\end{cases}
\]
The primes in \(\mathcal B\) are odd. This product is bounded below by a
positive \(q\)-dependent constant on the compatible progression.

**Theorem VqPower.** Fix \(5/6<\theta<1\) and \(H=X^\theta\). For ordered
positive summands, put
\[
r_{q,H}(N)=\sum_{\substack{a+b=N\\ |a-b|\le H}}f(a)f(b).
\]
Then
\[
\sum_{X<N\le2X}\big|r_{q,H}(N)-H\mu_X^2S_q(N)\big|^2
 \ll_{q,\theta} XH^2\mu_X^4 Z^{-1/4}. \tag{P1}
\]
Thus \(r_{q,H}(N)/(H\mu_X^2)=S_q(N)+O(Z^{-1/16})\) outside
\(O_{q,\theta}(X/Z^{1/8})\) targets. For any fixed \(\varepsilon>0\), at most
\(O_{q,\theta}(\varepsilon^{-2}X/Z^{1/4})\) targets have normalized error
at least \(\varepsilon\). At \(q=12\) the scale is
\(c_{12}^2 H/Z^{3/2}\). In particular the theorem applies at \(H=X^{0.99}\).

The proof first estimates the **one-variable** exponential sum over each
interval of length comparable to \(H\). Its major arcs require arithmetic
progression asymptotics on shorter intervals. For the other frequencies a
factorization treats every supported integer. Parseval then estimates the
difference between two local binary convolutions. Averaging over interval
starts gives the sharp balanced-pair cutoff of (P1).

The companion audit derives a conditional number-field prime-ideal estimate
and identifies the additional hypotheses needed for the other proposed
transfers. [Eisenstein Hecke transfer audit](EISENSTEIN_HECKE_TRANSFER_AUDIT.md).

## A short-interval estimate for character twists

Throughout the proof choose fixed real numbers satisfying
\[
2/3<\gamma<2\theta-1,\qquad 1-\gamma<\tau<1/3,\qquad
1-\tau<\gamma_0<\gamma .
\]
For every fixed \(A,B>0\), \(t^{\gamma_0}\le y\le t\), \(d\le(\log t)^B\)
and character \(\chi\bmod d\) whose primitive character has conductor
not dividing \(L\), we claim
\[
 \sum_{t<n\le t+y}f(n)\chi(n)\ll_{q,A,B,\theta}y(\log t)^{-A}. \tag{P2}
\]
For the untwisted sequence, the sum equals
\(\int_t^{t+y}D(v)\,dv+O_{q,A,\theta}(y(\log t)^{-A})\), with a smooth
density
\[
 D(v)=c_q(\log v)^{\alpha-1}\{1+O_q(1/\log v)\}. \tag{P3}
\]

To see this, put \(Q_0=\operatorname{lcm}(d,L)\). On \(\Re s>1\) the
Dirichlet series factors as
\[
 F_\chi(s)=\sum_{n\ge1}\frac{f(n)\chi(n)}{n^s}
       =G_\chi(s)\prod_{\psi\bmod L}L(s,\chi\psi)^\alpha .
\]
Use the imprimitive product characters modulo \(Q_0\). Orthogonality makes
the linear terms of the Euler logarithms agree at each prime. Both linear
terms vanish at primes dividing \(Q_0\). The remainder \(\log G_\chi\) is an
absolutely convergent series on \(\Re s>1/2\), uniformly bounded for
\(\Re s\ge0.59\) by \(O_q(\sum_p p^{-1.18})\). In (P2) all the product
characters are nonprincipal. In (P3) there is one principal factor. The
exponents are real and positive with sum one. On every branch the magnitude
of the product is therefore the geometric mean of the corresponding
\(L\)-function magnitudes.

Take \(T\in[t^\tau,2t^\tau]\) and \(\sigma_0\in[0.59,0.60]\), avoiding
the finitely many zero ordinates and abscissae on the relevant boundaries.
Truncated Perron for the difference of the two summatory functions, at
\(b=1+1/\log t\), has error \(O(t\log^2t/T+1)\). Indeed the coefficients
are bounded by one. Shift endpoints by at most one to half-integers in the
truncation estimate and sum the usual kernel bound near each endpoint
harmonically. This changes a count and the integral of (P3) by \(O(1)\).
On the shifted interval, still denoted \((t,t+y]\), its Perron kernel
\[
 K(s)=\{(t+y)^s-t^s\}/s
\]
is \(O(yt^{\sigma-1})\) throughout the rectangle. On the two horizontal
edges it is also \(O(t^\sigma/T)\).

For \(Q_0\ll_q(\log t)^B\), the number of product zeros with
\(\Re\rho\ge v,\ |\Im\rho|\le T\), counted with multiplicity, satisfies
\[
 N(v,T)\ll_q (Q_0^2T)^{(5/2)(1-v)}\log^{13}(2Q_0T),\qquad 1/2\le v\le1 .
\]
For \(v\le4/5\) use Montgomery's exponent \(3(1-v)/(2-v)\), since
\(3/(2-v)\le5/2\). For \(v\ge4/5\) use Jutila's exponent
\((2+\varepsilon)(1-v)\) with \(\varepsilon=1/4\). Both sum over primitive
characters with conductor at most \(Q_0\), including the zeta factor via
the conductor-one character. Reflection gives the bound at negative heights.
Our list repeats each inducing primitive character at most \(\varphi(L)\)
times, and imprimitive Euler factors have zeros only on \(\Re s=0\).
[Montgomery, Theorem 1](https://doi.org/10.1007/BF01404638),
[Jutila, Theorem 1, (1.8)](https://journals.msp.org/mscand/article/download/1943/1942/1974),
and [Inoue, Lemma 1](https://www.numdam.org/item/10.5802/jtnb.1162.pdf).

All zeros in the rectangle have \(\Re\rho\le1-\delta\), with
\(\delta\gg_{q,B}(\log t)^{-3/4}\) for large \(t\). For heights at least
10 in absolute value, the Vinogradov–Korobov zero-free denominator is
\(O(\log(2Q_0)+(\log T)^{2/3}(\log\log T)^{1/3})\). This also covers
\(\zeta\) via a fixed imprimitive principal character. At bounded height
apply the classical zero-free region and \(1-\beta\gg_\xi Q_0^{-\xi}\)
for a possible real exceptional zero, taking \(\xi>0\) sufficiently small
in terms of \(B\). This is where ineffectivity enters.
[Khale, Theorem 1.1](https://doi.org/10.1093/qmath/haae010)
and [Siegel's theorem](https://doi.org/10.4064/aa-1-1-83-86).

On the rectangle, away from a fixed small disk about \(1\) in the untwisted
case, every branch satisfies
\[
 |F_\chi(\sigma+iu)|\ll_q(\log t)^{O(1)}
 (Q_0(T+2))^{\max(0,1-\sigma)/2}. \tag{P4}
\]
Here convexity needs only logarithmic losses. For a primitive factor of
conductor \(Q\), fix the output height \(u\) and set
\(\rho=1/\log(10Q(3+|u|))\). On the line \(1+\rho\) the Euler series costs
\(O(\rho^{-1})\). On the line \(-\rho\) the functional equation and
Stirling give \(O(\rho^{-1}(Q(3+|v|))^{1/2+\rho})\) at height \(v\).
For \(\zeta\) apply these bounds to \((s-1)\zeta(s)/(s+1)\) to remove the
pole. Multiply by \(\exp((s-iu)^2)\) and apply the three-lines theorem.
The Gaussian bounds the left-edge height factor by a constant times
\((3+|u|)^{1/2+\rho}\), using \(3+|v|\le(3+|u|)(1+|v-u|)\).
For \(\sigma\le1+\rho\) the resulting exponent is
\((1+\rho-\sigma)/2\). The \(\rho/2\) costs an absolute constant.
For \(\sigma>1+\rho\) the Euler-series bound suffices. On the rectangle
the finite Euler factors cost at most \(2^{\omega(Q_0)}\le Q_0\), a log
power. The principal pole multiplier is bounded away from zero outside
the disk. Taking the geometric mean gives (P4) on every branch.

Cut horizontally leftward from each distinct zero to the edge \(\sigma_0\),
taking the union of overlapping slits. In the untwisted case cut also from
the pole \(1\) to that edge. Continue the Euler-product branches from the right
into the slit rectangle. Small detours around zeros vanish by positivity of
the exponents. At the pole the detour vanishes since \(\alpha<1\). The combined
integral on the two sides of all zero slits is bounded by
\[
 \ll y(\log t)^{O(1)}
 \int_{\sigma_0}^{1-\delta}
 t^{-(1-v)}(Q_0(T+2))^{(1-v)/2}(Q_0^2T)^{(5/2)(1-v)}\,dv
 \ll_{q,A,B,\theta} y(\log t)^{-A}. \tag{P5}
\]
At each \(v\) the number of banks is at most twice the number of zeros
with real part at least \(v\). Absorbing all \(Q_0\) factors into the log
loss leaves the power \(t^{-(1-3\tau)(1-v)}\). The saving is positive and
\(\delta\log t\gg(\log t)^{1/4}\).
For the fixed untwisted product choose \(r>0\) so small that all its factors
are zero-free in \(|s-1|\le4r\). No zero slit can meet \(|s-1|<r\):
a zero to its right with \(|\Im\rho|<r\) would lie in the larger disk.
The part of the pole slit to the left of \(1-r\) costs \(O_q(yt^{-c})\)
for some \(c>0\) by (P4). The left vertical edge costs at most
\(y t^{-2/5+(6/5)\tau}(\log t)^{O(1)}\), using
\(\sigma_0\le3/5\) in \(-(1-\sigma_0)+\tau(1+(1-\sigma_0)/2)\).
The horizontal edges cost \(O(t^{1-\tau}(\log t)^{O(1)})\) by the second
kernel bound and \(\tau<1/3\). The same bound covers Perron truncation.
Since \(\gamma_0>1-\tau\), these edge and truncation errors are
\(O(y(\log t)^{-A})\) for every fixed \(A\).

For (P3), \(B_0(s)=(s-1)^\alpha F(s)\) extends holomorphically across
\(|s-1|\le2r\) with \(B_0(1)>0\). The remaining two sides of the real slit
\([1-r,1]\) give \(\int_t^{t+y}D(v)\,dv\) with
\[
 D(v)=\frac{\sin(\pi\alpha)}{\pi}
      \int_0^r u^{-\alpha} B_0(1-u)v^{-u}\,du .
\]
Expanding \(B_0(1-u)=B_0(1)+O_q(u)\) and using the gamma reflection formula
gives (P3) with \(c_q=B_0(1)/\Gamma(\alpha)\). This proves both short
estimates.

## Rational frequencies and arithmetic coefficients

Write \(e(v)=\exp(2\pi i v)\). Fix any \(B>0\), and let \(a/d\) be reduced,
\(d\le Z^B\). For a real interval \(J=(t,t+y]\) with \(t\asymp X\),
\(y\asymp X^\gamma\), split its integers by \(b\bmod d\), then write
\(g=(b,d)\), \(D_0=d/g\) and \(b'=b/g\bmod D_0\). In this class put
\(n=gn'\). Complete multiplicativity gives \(f(n)=f(g)f(n')\).
Among the characters modulo \(D_0\) those with conductor dividing \(L\)
are the annihilator of the subgroup attained by residues \(1\bmod L\).
Their character sum at the unit \(b'\) is
\[
 E_{D_0}(b')=\varphi((D_0,L))\mathbf1_{b'\equiv1\bmod(D_0,L)}
\]
by the Chinese remainder theorem. These twists on \(f\)-support all equal
\(f(n')\mathbf1_{(n',D_0)=1}\). Remove this coprimality by inclusion-exclusion
over squarefree \(s\mid D_0\) made up of allowed primes. The nonexceptional
twists use (P2) on \(J/g\), and inclusion-exclusion uses (P3) on \(J/(gs)\).
Here \(gs\le d^2\) is polylogarithmic, so each interval satisfies the
\(\gamma_0\) length condition, with a harmless increase of \(B\) for the
rescaled variables. Use an arbitrarily large \(A\) before summing the errors.

The density after rescaling is \(D(v/(gs))/(gs)\). Uniformly for \(v\asymp X\),
\[
 D(v/(gs))=\mu_X\{1+O_{q,B}(\log\log X/\log X)\}.
\]
The sum over all \(b\) of absolute coefficients multiplying these derivative
terms is at most
\[
 \sum_{g\mid d}g^{-1}\prod_{p\mid d}(1+1/p)\ll(\log\log(3d))^2 ,
\]
since \(\sum_{b'\ {\rm unit}}E_{D_0}(b')/\varphi(D_0)=1\).
Thus the main class mass relative to \(\mu_X\) equals
\[
 \nu_d(b)=\frac{f(g)}{g}\frac{E_{D_0}(b')}{\varphi(D_0)}
            \prod_{\substack{p\mid D_0\\ f(p)=1}}(1-1/p).
\]
Put \(b_{a/d}=\sum_{b\bmod d}\nu_d(b)e(ab/d)\), including \(b_{0/1}=1\).
We have proved, uniformly in the stated ranges,
\[
 \sum_{n\in J}f(n)e(an/d)
    =y\mu_X\{b_{a/d}+O(Z^{-1+o(1)})\}. \tag{P6}
\]

The local content of \(\nu_d\) will identify the Fourier model. At allowed
\(p\), the class modulo \(p^k\parallel d\) is uniform on all residues,
including zero. At \(p\in\mathcal B\) it is uniform on units. At
\(p^j\parallel L\) it is uniform on residues \(1\bmod p^{\min(j,k)}\).
These assertions follow from the displayed formula, using \(g\equiv1\bmod L\)
whenever \(f(g)=1\). CRT gives their product. In particular a nonzero
\(b_{a/d}\) requires \(d=d_0d_1\) with \(d_0\mid L\) and \(d_1\) squarefree
on \(\mathcal B\). Its magnitude is \(1/\varphi(d_1)\), so
\[
 |b_{a/d}|\le \varphi(L)/\varphi(d). \tag{P7}
\]

For the rest of the proof take \(I=(x,x+h]\), \(X/4\le x\le X\),
\(H\le h\le2H\), and put
\[
 F_I(\beta)=\sum_{n\in I} f(n)e(\beta n),\qquad
 J_I(v)=\int_I e(vt)\,dt,\qquad T_D=\lfloor h^2/(XZ^B)\rfloor .
\]
Choose a Dirichlet approximant \(a/d\) with \(d\le T_D\),
\(|v|=|\beta-a/d|\le1/(dT_D)\), using circular distance. When \(d<Z^B\),
partition \(I\) into equal pieces of length asymptotic to \(X^\gamma\).
The phase \(e(vn)\) varies on each piece by
\(O(X^{\gamma-(2\theta-1)}Z^B)\), which is smaller than any inverse log power.
Apply (P6) to each piece and restore the continuous phase. The result is
\[
 F_I(\beta)=\mu_X b_{a/d}J_I(v)+O(h\mu_X Z^{-1+o(1)}). \tag{P8}
\]
The same formula holds for any prescribed \(d\le Z^{1/4}\) and
\(|v|\le Z^{1/4}/h\), by the same partition.


## Bilinear sums on the remaining frequencies

For each fixed \(J_0>0\), choose \(B\) sufficiently large in terms of \(J_0\).
We prove
\[
 |F_I(\beta)|\ll_{q,\theta,J_0}h Z^{-J_0}
 \quad\text{if}\quad Z^B\le d\le T_D,\quad |\beta-a/d|\le d^{-2}. \tag{P9}
\]
A coefficient bounded by \(Z^C d_2(n)^C\), with fixed \(C\) and the usual
divisor-counting function \(d_2\), has square sum \(O(UZ^{C'})\) over
\(U<n\le2U\), \(n\ll X\), for fixed \(C'\). This follows by dominating the
fixed divisor moment by a higher divisor function and summing its convolution.

First consider a sum \(B_{II}\) of \(a_m b_k e(\beta mk)\) over
\(M<m\le2M,\ K<k\le2K,\ mk\in I\), with both coefficients bounded in this
sense and \(M,K\gg X^{1/6}\). Interchange the variables to arrange
\(M\le K\), so \(M\ll\sqrt X\). Apply Cauchy–Schwarz in \(k\).
Two values \(m,m+j\) can then contribute only for \(|j|\ll hM/X\).
For each such pair the remaining \(k\) range is an interval of length
\(O(h/M)\). Its exponential sum at nonzero \(j\) is
\(O(\min(h/M,\|j\beta\|^{-1}))\), with \(\|\cdot\|\) distance to the nearest
integer. For every \(j\), \(\sum_m |a_m a_{m+j}|\ll MZ^{C'}\).
Grouping differences in blocks of length at most \(d/2\) gives
\[
 \sum_{1\le j\le D}\min(V,\|j\beta\|^{-1})
 \ll (D/d+1)(V+d\log(2d)),
\]
since points within one block are \(1/(2d)\)-separated. With
\(D\asymp hM/X,\ V\asymp h/M\), including the diagonal, this gives
\[
 \frac{|B_{II}|^2}{h^2}
 \ll Z^{O(1)}
 \left\{\frac{X}{hM}+\frac1d+\frac{M\log X}{h}
                         +\frac{Xd\log X}{h^2}\right\}. \tag{P10}
\]
Here \(X/(hM)\ll X^{5/6-\theta}\). Each term admits the saving needed in
(P9) when \(B\) is sufficiently large.

For a linear inner sum, suppose the small index \(\ell\ll X^{1/2}\) has a
divisor-bounded coefficient. The other variable ranges over \(I/\ell\)
with weight \(1\), a fixed character modulo \(L\), or such a character times
\(\log k\). Remove the logarithm by partial summation. Splitting into classes
modulo \(L\) gives geometric sums with phase \(\ell L\beta\). The reduced
approximant to \(L\beta\) has denominator \(d'=d/(d,L)\asymp_L d\) and
error \(O_L((d')^{-2})\), so spacing holds on blocks of length \(c_L d'\) for
a sufficiently small \(c_L>0\). For \(U<\ell\le2U\) with \(U\ge Z^{B/2}\),
Cauchy–Schwarz and
\[
 \sum_{U<\ell\le2U}\min(V,\|\ell L\beta\|^{-1})^2
       \ll_L (U/d'+1)(V^2+d'V),\qquad V\asymp h/U ,
\]
give a bound relative to \(h\) of
\(Z^{O(1)}(1/d'+1/U+U/h+d'/h)^{1/2}\). The squared spacing estimate
follows by summing \(\min(V,C_Ld'/j)^2\) within each block, apart from its
one largest term. If \(U<Z^{B/2}\), then \(\ell<c_Ld'\) for large \(X\)
and \(\|\ell L\beta\|\gg_L1/d'\). The relative bound
\(Z^{O(1)} U d'/h\) is power-small since \(d'/h\le h/(XZ^B)\).
These estimates hold for every prefix of the inner interval, as required
by partial summation.

It remains to produce precisely these sums for the full function \(f\).
Write \(P^+(n)\) and \(P^-(n)\) for the largest and smallest prime
factors of \(n>1\).

1. **When \(P^+(n)\le X^{5/6}\)**, list the prime factors of \(n\) in
   decreasing order, with multiplicity. Let \(m\) be the shortest initial
   product exceeding \(X^{1/6}\), and \(k=n/m\). If the product uses at
   least two factors, both the preceding product and the last factor are
   at most \(X^{1/6}\), so \(m\le X^{1/3}\) in that case.
   Thus \(X^{1/6}<m\le X^{5/6}\) and
   \(k\gg X^{1/6}\). Conversely the conditions
   \(m>X^{1/6}\), \(m/P^-(m)\le X^{1/6}\),
   \(P^+(m)\le X^{5/6}\) and \(P^+(k)\le P^-(m)\) select that unique
   factorization, including at repeated primes. The coefficients are
   \(f(m)f(k)\). Separate the last inequality by expanding the indicator
   of \([0,\lceil2X\rceil]\) on a cyclic group of order greater than
   \(8X\) and at most \(9X\), evaluated at \(P^-(m)-P^+(k)\).
   Its normalized Fourier coefficients have total absolute value
   \(O(\log X)\) by the geometric series. Each mode introduces bounded
   separate factors on \(m\) and \(k\). Apply (P10) on \(O(\log^2 X)\)
   dyadic blocks.

2. **When \(P^+(n)>X^{5/6}\)** this prime is unique and has exponent one.
   Write \(n=up\), \(f(u)=1,\ p\equiv1\bmod L\). Keep \(u\le x/X^{5/6}\).
   On this range \(up\in I\) already enforces \(p>X^{5/6}\). The omitted
   band \(x/X^{5/6}<u\le(x+h)/X^{5/6}\) contributes at most
   \(O((h/X^{5/6}+1)(h/X^{1/6}+1))\). Replace the prime indicator by
   \(\Lambda(p)/\log(x/u)\), keeping its congruence. On primes,
   \(\log p=\log(x/u)+O(h/X)\), and
   \(\sum_{u\ll X^{1/6}}(h/u+1)\ll h\log X\), so this replacement
   costs \(O(h^2(\log X)^{O(1)}/X)\). The extra prime powers cost
   \(O(X^{7/12}(\log X)^{O(1)})\), by summing
   \(O((X/u)^{1/2}\log X)\) over \(u\ll X^{1/6}\). All three errors
   are \(O(hZ^{-J_0})\) for fixed \(\theta<1\).

   Detect the progression by characters modulo \(L\). Let
   \(\mu_{\rm M}\) denote the Möbius function, \(*\) Dirichlet convolution
   in this paragraph, and take \(U_0=V_0=X^{1/6}\). At arguments
   exceeding \(V_0\), Vaughan's identity takes the form
   \[
   \Lambda=(\mu_{\rm M})_{\le U_0}*\log
          -(\mu_{\rm M})_{\le U_0}*\Lambda_{\le V_0}*1
          +((\mu_{\rm M})_{>U_0}*1)*\Lambda_{>V_0}.
   \]
   This follows from \(\mu_{\rm M}*1=\mathbf1_{\{1\}}\) and
   \(1*\Lambda=\log\). Each fixed character factors over the convolution
   variables. In the first two terms combine \(u\) with the bounded small
   product to obtain \(\ell\ll X^{1/6}U_0V_0=X^{1/2}\).
   The remaining variable has the linear weight treated above. In the last
   term both displayed convolution indices exceed \(X^{1/6}\). Combine
   \(u\) with the first and apply (P10). The cutoff on \(u\) and the
   factor \(1/\log(x/u)\) belong to its coefficient. All resulting
   coefficients are bounded by fixed divisor and logarithm powers.

These finitely many types of sums, with their dyadic and Fourier losses,
prove (P9).

## Comparing the local binary convolutions

Choose \(J_0\) sufficiently large, then its \(B\) from (P9). Assign each
\(\beta\in\mathbb R/\mathbb Z\) one of the Dirichlet approximants used above,
breaking ties by a fixed order. Let \(P(\beta)=\mu_X b_{a/d}J_I(\beta-a/d)\)
when its assigned \(d<Z^B\), and zero otherwise. Equations (P8)–(P9) give
\[
 \|F_I-P\|_\infty\ll h\mu_X Z^{-1+o(1)} .
\]
Parseval gives \(\int |F_I|^2=\sum_I f\ll h\mu_X\) by (P6) at \(d=1\).
Hereafter unmarked Fourier integrals are over the unit circle with normalized
Lebesgue measure. Bound \(|J_I(v)|\ll\min(h,1/|v|)\) for \(|v|\le1/2\).
Since each \(\beta\) is assigned a single approximant, (P7) gives
\[
 \int |P|^2\ll h\mu_X^2 \sum_{d\le Z^B}\frac{\varphi(L)^2}{\varphi(d)}
       \ll_{q,B}h\mu_X^2\log\log X .
\]
The last estimate follows by expanding
\(d/\varphi(d)=\sum_{s\mid d}\mu_{\rm M}(s)^2/\varphi(s)\) and summing
harmonically. As \(\alpha<1\), \(\mu_X\log\log X=O(1)\). Consequently
\[
 \int |F_I^2-P^2|^2
 \le 2\|F_I-P\|_\infty^2\left(\int |F_I|^2+\int |P|^2\right)
 \ll h^3\mu_X^3 Z^{-2+o(1)}. \tag{P11}
\]

Put \(z=Q=Z^{1/4}\), and take \(X\) large enough that \(z\) exceeds all
prime factors of \(L\). Define the periodic density
\[
 W=L\prod_{\substack{p\in\mathcal B\\ p\le z}}p,\qquad
 w_z(t)=L\mathbf1_{t\equiv1\bmod L}
       \prod_{\substack{p\in\mathcal B\\ p\le z}}
            \frac{\mathbf1_{p\nmid t}}{1-1/p}.
\]
Its period \(W\le L\exp(O(z\log(2z)))=X^{o(1)}\) and its mean is one.
The CRT computation after (P6) identifies its Fourier coefficient for
testing against \(e(at/d)\) as \(b_{a/d}\) whenever \(d\le Q\), and on
each of its other modes as the same \(b_{a/d}\) for that mode.
Let \(M_I(\beta)=\mu_X\sum_{n\in I}w_z(n)e(\beta n)\).

The major arcs have circular radius \(Q/h\) around reduced \(a/d\) with
\(d\le Q\). They are disjoint for large \(X\) and have measure \(O(Q^3/h)\).
On each, (P8) holds at its center. Every other Fourier mode of \(w_z\) is at
distance \(\gg1/(QW)\) from the entire arc. All coefficients have modulus
at most one and there are at most \(W\) modes, so their total contribution
to \(M_I\) is \(O(\mu_X QW^2)\). For the resonant mode the discrete
interval kernel differs from \(J_I(v)\) by \(O(1+Q)\) by summation over
unit intervals. This proves \(|F_I-M_I|\ll h\mu_X Z^{-1+o(1)}\) there,
with both transforms \(O(h\mu_X)\). Dividing the integral of
\(|F_I^2-M_I^2|^2\) over these arcs by \(h^3\mu_X^4\) gives
\(O(Q^3 Z^{-2+o(1)})\).

On the minor arcs bound \(\int |P|^4\) by the sum of its single-kernel
integrals over the assigned regions. At centers of denominator \(\le Q\)
the remaining kernel has \(|v|>Q/h\), integral
\(\int |J_I(v)|^4\,dv\ll h^3/Q^3\) there. Other centers use the full
\(O(h^3)\) bound. By (P7),
\[
 \sum_{a/d\ {\rm reduced}} |b_{a/d}|^4\ll_q\sum_d\varphi(d)^{-3}<\infty,
 \qquad
 \sum_{\substack{a/d\ {\rm reduced}\\d>Q}} |b_{a/d}|^4\ll_q Q^{-2+o(1)} ,
\]
where the fractions are taken modulo one. Equation (P11) now gives
\[
 \frac{\int_{\rm minor}|F_I|^4}{h^3\mu_X^4}
   \ll Z^{-2+o(1)}/\mu_X+Q^{-3}+Q^{-2+o(1)}
   \ll Z^{-1/2+o(1)}.
\]
For \(M_I/\mu_X\) take disjoint radius \(1/(3W)\) neighborhoods of its modes.
On each the sum of nonresonant kernels is \(O(W^2)\), and the full transform
therefore has fourth power at most eight times the resonant fourth power plus
\(O(W^8)\). Off these neighborhoods the whole transform is \(O(W^2)\).
Integrate the resonant geometric kernels as above, using
\(\left|\sum_{n\in I}e(vn)\right|\ll\min(h,1/\|v\|)\).
The normalized minor fourth moment of \(M_I\) is
\(O(Q^{-3}+Q^{-2+o(1)}+W^8/h^3)\).

Let \(\star\) denote ordinary additive convolution on the integers and
\[
 E_I(N)=((f\mathbf1_I)\star(f\mathbf1_I))(N)
       -\mu_X^2((w_z\mathbf1_I)\star(w_z\mathbf1_I))(N).
\]
Adding the major and minor bounds and applying Parseval proves, uniformly in
the allowed \(x,h\),
\[
 \sum_N|E_I(N)|^2\ll h^3\mu_X^4 Z^{-1/2+o(1)}. \tag{P12}
\]

## Recovering the balanced pairs and the four-form energy

Let \(\mathbf E_X u=X^{-1}\sum_{X<N\le2X}u(N)\) and
\(R(N)=r_{q,H}(N)/(H\mu_X^2)\). Define
\(S_z(N)=W^{-1}\sum_{t\bmod W}w_z(t)w_z(N-t)\).
CRT gives the same formula as \(S_q\) with the product restricted to \(p\le z\).
The local condition excludes one possible residue for \(t\) when \(p\mid N\)
and two otherwise. Its periodic second moment is
\[
 K_z=L\prod_{\substack{p\in\mathcal B\\p\le z}}(1+(p-1)^{-3})\le
 K_q=L\prod_{p\in\mathcal B}(1+(p-1)^{-3})<\infty .
\]
Indeed the local factor has mean one and variance \(1/(p-1)^3\) over \(N\bmod p\).
Since \(W=X^{o(1)}\), \(\mathbf E_X S_z^2=K_z+O(WK_z/X)\).
For an interval of length \(u\), the sum of \(w_z(t)w_z(N-t)\) is
\(uS_z(N)+O(W^3)\), by completing periods and using \(0\le w_z\le W\).

Take \(\eta=Z^{-1/8}\). For the upper comparison use \(h=(1+\eta)H\) and
all starts \(x\) satisfying
\[
 (N+H)/2-h<x<(N-H)/2 .
\]
They have total length \(\eta H\) and lie in \([X/4,X]\) for large \(X\).
Both summands of every balanced pair at \(N\) belong to each such \(I\).
In the model convolution the possible first summands form an interval of
length at most \(h\). Thus the excess of \(R\) over \((1+\eta)S_z\) is
bounded above by \((|E_I(N)|/\mu_X^2+O(W^3))/H\).
For the lower comparison take \(h=H\) and
\(|2x+H-N|<\eta H\), again a set of starts of length \(\eta H\) inside
\([X/4,X]\). Every pair in the interval is balanced, and the possible
first summands form an interval of length
\(H-|2x+H-N|\ge(1-\eta)H\). The excess of \((1-\eta)S_z\) over \(R\)
satisfies the same bound.

Average each squared positive excess over its starts and interchange \(x,N\).
Apply (P12) for every \(x\), enlarging the sum to all \(N\) on the right.
The resulting bounds are
\[
 \begin{split}
 \mathbf E_X (R-(1+\eta)S_z)_+^2+
 \mathbf E_X ((1-\eta)S_z-R)_+^2
  &\ll Z^{-1/2+o(1)}/\eta+W^6/H^2,\\
 \mathbf E_X |R-S_z|^2
  &\ll \eta^2+Z^{-1/2+o(1)}/\eta+W^6/H^2\ll Z^{-1/4}.
 \end{split} \tag{P13}
\]
The second line uses the bounded second moment of \(S_z\).

For the passage to actual integers, write the product in \(S_q\) as
\(A B(N)\) where \(A=\prod_{p\in\mathcal B}(1-(p-1)^{-2})>0\) and
\(B(N)=\prod_{p\mid N,\ p\in\mathcal B}(1+1/(p-2))\).
For the corresponding truncations, \(A-A_z=O(1/z)\).
On squarefree products of primes in \(\mathcal B\) put
\(u(p)=2/(p-2)+(p-2)^{-2}\) and extend multiplicatively, with zero
on all other integers except \(u(1)=1\). Then
\(\sum_d u(d)/d<\infty\), and
\[
 \mathbf E_X (B-B_z)^2\le\mathbf E_X(B^2-B_z^2)
 \le 2\sum_{d:\,P^+(d)>z}\frac{u(d)}d
 \ll_q\sum_{p>z}p^{-2}\ll 1/z .
\]
The middle bound uses \(\#\{N\le2X:d\mid N\}\le2X/d\). Extract a prime
larger than \(z\) from the squarefree divisor for the last estimate.
The same expansion gives \(\mathbf E_X B^2=O_q(1)\), so
\(\mathbf E_X|S_q-S_z|^2\ll_q1/z\). With \(z=Z^{1/4}\), (P13) now proves
(P1). Chebyshev's inequality gives the stated exception counts.

Finally \(K_q-K_z\ll_q z^{-2}\) and the normalized \(L^2\) distance from
\(R\) to \(S_z\) is \(O(Z^{-1/8})\). Hence the four-form energy has its
specified arithmetic coefficient,
\[
 \sum_{X<N\le2X}r_{q,H}(N)^2
  =\{K_q+O_{q,\theta}(Z^{-1/8})\}XH^2\mu_X^4 . \tag{P14}
\]
The sum in (P14) is exactly the sum of
\(f(a)f(N-a)f(c)f(N-c)\) over the two balanced windows.
This proves the upper bound (E) in the stated range of fixed exponents.
[Fixed-width theorem and general power-window reduction](PRESCRIBED_PRIME_BINARY_VARIANCE.md).

The proof yields a logarithmic saving in the exceptional proportion and
uses ineffective constants from Siegel's theorem. The separate
[D12Power lower-count theorem](PRESCRIBED_PRIME_POWER_SAVING_EXCEPTIONS.md) gives
a power-saving exceptional set at \(q=12,H=X^{99/100}\) through a weighted
rough-kernel subset. A power-saving variance for the full count, the
every-compatible-target assertion and an effective finite cutoff remain open here. This is an internally reviewed
derivation. External mathematical review and the priority audit remain
outstanding.
