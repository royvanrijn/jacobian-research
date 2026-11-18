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

**Theorem VqPower.** Fix \(9/10<\theta<1\) and \(H=X^\theta\). For ordered
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

## A short-interval estimate for character twists

Throughout the proof choose fixed real numbers satisfying
\[
4/5<\gamma<2\theta-1,\qquad 1-\gamma<\tau<1/5,\qquad
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
 N(v,T)\ll_q (Q_0^2T)^{3(1-v)/(2-v)}\log^{13}(2Q_0T),\qquad 1/2\le v\le1 .
\]
Montgomery's estimate sums over all primitive characters with conductor
at most \(Q_0\). Our list repeats each at most \(\varphi(L)\) times, and
imprimitive Euler factors have zeros only on \(\Re s=0\).
[Montgomery, *Zeros of L-functions*, Theorem 1](https://doi.org/10.1007/BF01404638)
and [Inoue, Lemma 1, p. 282](https://www.numdam.org/item/10.5802/jtnb.1162.pdf).

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
 (Q_0(T+2))^{2\max(0,1-\sigma)}. \tag{P4}
\]
For nonprincipal \(L\)-functions truncate the defining sum at
\(M=\lceil(Q_0(T+2))^2\rceil\). The tail is \(O(Q_0(1+|s|)M^{-\sigma})\)
by partial summation and periodic character cancellation. The head is
\(O(\log(2M)M^{\max(0,1-\sigma)})\). Euler summation for the fixed principal
factor adds \(M^{1-s}/(s-1)\), covered away from the disk.

Cut horizontally leftward from each distinct zero to the edge \(\sigma_0\),
taking the union of overlapping slits. In the untwisted case cut also from
the pole \(1\) to that edge. Continue the Euler-product branches from the right
into the slit rectangle. Small detours around zeros vanish by positivity of
the exponents. At the pole the detour vanishes since \(\alpha<1\). The combined
integral on the two sides of all zero slits is bounded by
\[
 \ll y(\log t)^{O(1)}
 \int_{\sigma_0}^{1-\delta}
 t^{-(1-v)}(Q_0(T+2))^{2(1-v)}(Q_0^2T)^{3(1-v)}\,dv
 \ll_{q,A,B,\theta} y(\log t)^{-A}. \tag{P5}
\]
At each \(v\) the number of banks is at most twice the number of zeros
with real part at least \(v\). The net power of \(t\) in the integrand is
at most \(-(1-5\tau-o(1))(1-v)\), and \(\delta\log t\gg(\log t)^{1/4}\).
For the fixed untwisted product choose \(r>0\) so small that all its factors
are zero-free in \(|s-1|\le4r\). No zero slit can meet \(|s-1|<r\):
a zero to its right with \(|\Im\rho|<r\) would lie in the larger disk.
The part of the pole slit to the left of \(1-r\) costs \(O_q(yt^{-c})\)
for some \(c>0\) by (P4). The left vertical edge costs at most
\(y t^{-0.40+1.82\tau+o(1)}\). The horizontal edges and Perron truncation
cost \(O(t^{1-\tau}(\log t)^{O(1)})\). All are absorbed by (P5).

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
