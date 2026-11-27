# Erdős #365: a subquarter bound for consecutive powerful numbers

Let \(N(X)=\#\{n\le X:n,n+1\text{ are powerful}\}\), where all
positive prime valuations of a powerful integer are at least two. We prove
\[
 N(X)\ll_\varepsilon X^{12127/49920+\varepsilon}
 \qquad (X\ge2,\ \varepsilon>0).                         \tag{1}
\]
The exponent is approximately \(0.242929\). The proof combines
Plymesser's published preprint estimate for reduced determinant windows
with a new use of a fixed quadratic norm fiber and the first moment of
quadratic class group \(3\)-torsion. The polylogarithmic bound in
[Erdős #365](https://www.erdosproblems.com/365) remains open. An exact
Pell-period count gives an [extended limit for \(N(X)/\log X\)](#logarithmic-limit-and-the-remaining-count).
The [small-unit asymptotic](#a-degree-one-projection-for-larger-pell-units) includes \(\epsilon_D\le D^{125/106}\).

[Reuss, Theorem 4](sources/reuss-consecutive-1212.3150v2.pdf) gives
the exponent \(29/100\).
[Plymesser, Theorem 1.1 and Remark 9.1](sources/plymesser-powerful-products-hexagon2611.00602-v1.0.1.pdf)
give \(1/4\) and identify the binding corner of that proof. There one
cube kernel is bounded, the other has size \(X^{1/4}\), and its
square variable has size \(X^{1/8}\). The fiber estimate below handles
both orientations of the corner. The optimization then uses the same
estimate throughout the box domain.

## Exact count by Pell periods

**Proposition (Exact Pell-period count).** For each squarefree \(D>1\),
write \(\epsilon_D=T_D+U_D\sqrt D\) for the least positive solution of
\(T^2-DU^2=1\), viewed as a unit greater than one in
\(\mathbb Z[\sqrt D]\). Define
\[
 R_D=\frac{D}{\gcd(D,U_D)},\qquad
 L_D=\begin{cases}R_D&4\mid U_D,\\ 2R_D&4\nmid U_D,\end{cases}
 \qquad \tau_D=L_D\log\epsilon_D.
\]
For \(Y=\operatorname{arcosh}(2X+1)\), put
\(A(Y)=\#\{D>1:D\text{ squarefree},\ \tau_D\le Y\}\). Then
\[
 N(X)=\sum_{\substack{D>1\\D\ {\rm squarefree}}}
       \left\lfloor \frac{Y}{\tau_D}\right\rfloor
     =\sum_{k\ge1} A(Y/k).                                            \tag{2}
\]

To see this, use the squarefree part \(D\) of \(n(n+1)\). It exceeds one
because two positive consecutive integers cannot both be squares.
Since the two integers are coprime, both are powerful exactly when
their product is powerful. With \(n(n+1)=D V^2\), the pair corresponds
to \(T=2n+1,\ U=2V\) in \(T^2-DU^2=1\). The product is powerful exactly
when \(D\mid V\). Also \(V\) is even for a powerful pair, since the even
member has valuation at least two. Equivalently
\(\operatorname{lcm}(4,D)\mid U\). Conversely this last divisibility
forces \(T\) odd by the Pell equation modulo four, and gives
\(n(n+1)=D(U/2)^2\) with \(D\mid U/2\), using \(4\mid U\) if \(D\) is
even. This recovers the unique pair from \(D,T,U\).

All positive Pell solutions are
\(\epsilon_D^j=T_{D,j}+U_{D,j}\sqrt D\) with \(j\ge1\). For any \(p\mid D\),
the binomial expansion gives
\(U_{D,j}\equiv jT_D^{j-1}U_D\bmod p\), and \(T_D\) is a unit modulo
\(p\). Thus \(D\mid U_{D,j}\) exactly when \(R_D\mid j\). If \(D\) is even,
\(T_D\) is odd and \(U_D\) is even by the equation modulo eight, since
an odd \(U_D\) would force \(8\mid D\). Hence \(R_D\) is always odd.

If \(4\mid U_D\), it divides the second coefficient of every power. If
\(U_D\equiv2\bmod4\), then \(T_D\) is odd and the binomial expansion
modulo four gives \(U_{D,j}\equiv jT_D^{j-1}U_D\), so \(4\mid U_{D,j}\)
exactly for even \(j\). If \(U_D\) is odd, then \(D\) is odd and \(T_D\)
is even. Modulo two the second coefficient is odd in every odd power,
and the second coefficient of the square is \(2T_DU_D\), divisible by
four. This proves that the admissible exponents are precisely the
multiples of \(L_D\). In particular each \(D\) has a first pair. For
\(j=kL_D\) its smaller member is \((\cosh(k\tau_D)-1)/2\), giving (2).
The sums are finite since \(\tau_D\ge\log\epsilon_D>\tfrac12\log D\).

[Walker, Theorems 2.2, 2.5, 3.2 and 3.5](https://www.fq.math.ca/Scanned/14-2/walker.pdf)
gives the divisibility and power conditions for separate Pell families.
[Aktaş–Murty, §1 and Theorem 1.1](https://doi.org/10.1142/S1793042117500142)
uses the map from a powerful pair to a square-upper pair and connects
the stronger exponent \(C<2\) to prime-radicand unit exceptions.
[Sayim, Part II, Proposition II.2.3 and §5.3](https://raw.githubusercontent.com/berkay-yueksel-sayim/erdos-364-365-powerful-numbers-pell/142612d9d464761ddf0af18bc53076e1aa124343/papers/part-II_pairs-of-consecutive-and-near-consecutive-powerful-numbers.pdf)
gives separate family rates and asks whether their sum converges.
Equation (2) uses the single kernel of the coprime product, with the
condition modulo four included in its period.

[PI1](PELL_PERIOD_RING_CLASS_INERTIA.md#the-full-relative-order-and-the-powerful-pair-period)
identifies the inverse period multiplier with a relative narrow ring class
kernel divided by the modulus. PI3 shows that its odd inertia can lie
arbitrarily deep in the ring class group even when it is generated by a
power of a single actual prime class. The ring class identification gives
an exact weight and leaves the remaining summation open.

## Logarithmic limit and the remaining count

**Corollary (Logarithmic limit).** In the extended positive reals,
\[
 \lim_{X\to\infty}\frac{N(X)}{\log X}
  =\sum_{\substack{D>1\\D\ {\rm squarefree}}}\frac1{\tau_D}.           \tag{3}
\]
Indeed, \(Y/\log X\to1\), and \(\lfloor Y/\tau_D\rfloor/Y\) tends to
\(1/\tau_D\) and is bounded by \(1/\tau_D\). Use dominated convergence
if the series is finite. If it diverges, each finite partial sum bounds
the limit inferior from below.

For comparison with Sayim's separate family rates, \(n+1\) is a square
exactly when \(j\) is even. When \(j=2h\) the Pell square gives
\(n+1=T_{D,h}^2\). Conversely if \(n+1=a^2\) then \(n=Db^2\), and the
unit for \(T=2n+1\) is \((a+b\sqrt D)^2\). Thus for odd \(L_D\) the
square-upper and remaining pairs each have rate \(1/(2\tau_D)\).
The remaining pairs belong to one ordered-kernel family: writing
\(n_1=d a^2,\ n_1+1=e b^2\) for the first pair, the odd powers of
\(b\sqrt e+a\sqrt d\) still have integral coefficients of
\(\sqrt e,\sqrt d\), and their squares are the corresponding Pell
powers with kernel \(D=de\). For even \(L_D\) all the pairs have a
square upper member. This identifies the sum in (3) with the rates in
Sayim's convergence question, using a cutoff on the lower member
instead of the upper. The uniform bound on each floor term in (2)
shows that convergence alone suffices for the logarithmic asymptotic.

For each fixed \(C>1\), (2) gives \(N(X)=O((\log X)^C)\) if and only if
\(A(Y)=O(Y^C)\). In one direction use \(A(Y)\le N((\cosh Y-1)/2)\).
In the other sum \(\sum_{k\ge1} O((Y/k)^C)\), since \(A\) vanishes
below a fixed positive threshold. In particular the polylogarithmic
question reduces to the number of kernels whose first pair occurs
by the given height. There is also an exact truncation estimate. For
any fixed integer \(q\ge2\),
\[
 N(X)-\sum_{1\le k<q} A(Y/k)
   \le q N\big((\cosh(Y/q)-1)/2\big)
   \ll_{q,\varepsilon} X^{12127/(49920q)+\varepsilon}.
\]
For the first inequality with the factor \(q\), group indices
\(qh\le k<q(h+1)\) and use monotonicity of \(A\). The second inequality
uses (1) and \((\cosh(Y/q)-1)/2=O_q(X^{1/q})\).

The forward problem is therefore to control the small periods
\(L_D\log\epsilon_D\). For any fixed \(\delta>0\), kernels with
\(R_D\ge D^\delta\) or \(\log\epsilon_D\ge D^\delta\) give a
polylogarithmic contribution by (2), since these conditions imply
\(\tau_D\ge D^\delta\) and there are at most \(Y^{1/\delta}\) such
kernels below period \(Y\). The following result removes the kernels
with small common divisor from the convergence question in (3).

## Logarithmic asymptotic for kernels with small common divisor

**Proposition (Restricted logarithmic asymptotic).**

Fix \(B>0\) and \(0\le b<1/2\). Let \(\mathcal D_{B,b}\) consist of the
squarefree \(D>1\) with \(\gcd(D,U_D)\le B(\log D)^b\), and let
\(N_{B,b}(X)\) count pairs with \(D\) in this set. Then
\[
 \begin{split}
 S_{B,b}&=\sum_{D\in\mathcal D_{B,b}}\frac1{\tau_D}<\infty,\\
 N_{B,b}(X)&=S_{B,b}Y
    +O_{B,b}\big(Y(\log Y)^{b-1/2}\sqrt{\log\log Y}\big)
       \qquad (Y=\operatorname{arcosh}(2X+1),\ Y\to\infty).           \tag{4}
 \end{split}
\]

[Diao, Theorem 1.3](https://arxiv.org/html/2610.28001v1#S1.SS3)
bounds the number of nonsquare \(D<Z\) with
\(\epsilon_D\le D^\theta\) by \(O(Z^{1-c/\theta^2})\), uniformly for
\(3\le\theta\le c\sqrt{\log Z/\log\log Z}\) with some absolute \(c>0\).
The unit in that theorem is the least positive solution of \(T^2-DU^2=1\)
as in (2). This uses the stated uniform growing range, beyond the
qualitative almost-all statement of Corollary 1.4.

On \(W\le D<2W\), take
\(\theta=\kappa\sqrt{\log W/\log\log W}\) for an absolute
\(\kappa>0\) small enough that the theorem applies with \(Z=2W\) and
\(c/\kappa^2\ge3\). For large \(W\), the kernels with
\(\log\epsilon_D\le\theta\log W\) therefore number
\(O(W/(\log W)^3)\). A kernel in \(\mathcal D_{B,b}\) has
\(1/\tau_D\le B(\log D)^b/(D\log\epsilon_D)\). For the small-unit
kernels use \(\log\epsilon_D>\tfrac12\log D\), giving total reciprocal
weight \(O_B((\log W)^{b-4})\) on the interval. For the other kernels,
\(\log\epsilon_D>\theta\log W\) gives total weight
\(O_B((\log W)^{b-3/2}\sqrt{\log\log W})\). Summing these bounds over
doublings from \(W=Z\) to infinity shows that the rate tail beyond
\(Z\) is \(O_{B,b}((\log Z)^{b-1/2}\sqrt{\log\log Z})\).
In (2), keep the floors for \(D\le\sqrt Y\), where their difference
from the linear terms sums to \(O(\sqrt Y)\). The rest of the
difference is bounded by \(Y\) times the rate tail, proving (4).

In particular all kernels with \(\gcd(D,U_D)=1\) contribute a finite
rate, with (4) at \(B=1,\ b=0\). For every fixed \(b<1/2\), any
divergence of the full rate must come from kernels with
\(\gcd(D,U_D)>B(\log D)^b\). The full series and the polynomial bound
for \(A(Y)\) remain open.

## Logarithmic asymptotic for small fundamental units

**Proposition (Small-unit logarithmic asymptotic).**

Fix \(1/2\le\alpha<64/57\), and count only pairs whose product kernel \(D\)
has \(\epsilon_D\le D^\alpha\). Their count satisfies
\[
 N_\alpha(X)= S_\alpha Y+O_\alpha(Y^{1-\delta_\alpha}),\qquad
 S_\alpha=\sum_{\substack{D>1\ {\rm squarefree}\\\epsilon_D\le D^\alpha}}
                 \frac{1}{\tau_D}<\infty,\qquad
 Y=\operatorname{arcosh}(2X+1),                                      \tag{5}
\]
for some \(\delta_\alpha>0\). This restriction permits an unbounded
\(\gcd(D,U_D)\).

Write \(H_D=\gcd(D,U_D)\) and let \(h\) be the odd part of \(H_D\). On
\(W\le D<2W\), also restrict the Pell trace to \(P\le T_D<2P\).
Then \(U_D\ll P/\sqrt W\), and \(H_D\le 2h\).
Since \(\tau_D\ge D\log\epsilon_D/H_D\), it suffices to bound the
sum of \(h\) on this box. We obtain, for every \(\varepsilon>0\),
\[
 \sum_{D\ {\rm in\ the\ box}} h
 \ll_{\alpha,\varepsilon} W^\varepsilon
     \left(\sqrt W+\frac{P^{4/3}}{\sqrt W}+P^{57/64}\right).          \tag{6}
\]

For \(h\le4P^{1/3}\), fix \(U=U_D\) and a divisor \(h\mid U\).
The condition \(h\mid D\) implies \(T_D^2\equiv1\bmod U^2h\).
This congruence has at most \(4\cdot 2^{\omega(U)}\) roots, by the
prime-power bounds for square roots of one and the Chinese remainder
theorem. The traces in \(W\le D<2W\) range over an interval of length
\(O(U\sqrt W)\), since \(T_D^2=1+DU^2\). Multiplying the number of
traces by \(h\) gives
\(O_\varepsilon(U^\varepsilon(\sqrt W/U+h))\) for this choice.
Sum over \(h\mid U\), \(h\le4P^{1/3}\), and \(U\ll P/\sqrt W\).
The divisor bound gives the first two terms in (6), with logarithms
absorbed into \(W^\varepsilon\). For fixed \(T_D,U_D\), the value of
\(D\) is unique.

For \(h>4P^{1/3}\), write \(h=de\) by assigning each prime in \(h\)
to the one of \(T_D-1,T_D+1\) it divides. Its cube divides that factor
because \(T_D^2-1=DU_D^2\). If \(T_D\) is odd, divide both factors by
two before writing them as \(d^3u,e^3v\), which have difference one.
For \(T_D\) even, use the original odd factors, with difference two.
In either case \(d^3u,e^3v\asymp P\), \(d,e\ll P^{1/3}\), and the two
products are coprime. Each trace determines its squarefree kernel
\(D\) from \(T_D^2-1\), so this bounds the kernel count by a count of
solutions of \(e^3v-d^3u=s\), where \(s=\pm1,\pm2\) allows an exchange
of the sides.

Use [Reuss, Theorem 1](sources/reuss-consecutive-1212.3150v2.pdf) with
\(k=3,\ l=1\). To fit its strict windows, select \(D_0,E_0,x\) by
rounding \(7d/10,7e/10,(d^3u)/4\) down to the grid \(2^{j/16}\) with
\(j\in\mathbb Z\), after orienting the two scales so \(D_0\le E_0\).
There are \(O((\log P)^2)\) choices since \(x\asymp P\). Both cube
kernels exceed two by \(de>4P^{1/3}\). Hence \(D_0,E_0>1\), and
\(1<d/D_0,e/E_0<2\). Also
\[
 \frac{u}{x/D_0^3}=(D_0/d)^3(d^3u/x)
\]
lies strictly between one and two using the grid bounds
\((7/10)2^{-1/16}<D_0/d\le7/10\) and
\(4\le d^3u/x<4\cdot2^{1/16}\). The same holds for \(v/(x/E_0^3)\)
as \(e^3v/(d^3u)=1+O(P^{-1})\).
The hypothesis \(D_0E_0\gg x^{1/3}\) holds by the cutoff. Further,
\(de=h\mid U_D\) and \(U_D\ll P/\sqrt W\), so at
\(P\ll W^\alpha\) we have
\(D_0E_0\ll x^{1-1/(2\alpha)}\). Thus both \(\log(D_0E_0)\) and
\(\log((x/D_0^3)(x/E_0^3))\) are comparable to \(\log x\).
The comparability reflects \(\alpha<3/2\).

For the exponent calculation put \(a=\log D_0/\log x\),
\(b=\log E_0/\log x\), \(\psi=a+b\), and
\(m=(9/8)\psi(2-3\psi)\). The inequalities \(0\le a\le b\le1/3\) and
\(\psi\ge1/3\) may be used up to \(O(1/\log x)\), which changes the
resulting exponents by only a bounded multiplicative factor in the
count. Reuss bounds the unweighted count in the box by \(x^\varepsilon\)
times \(x\) to the smaller of the exponents
\[
 F_1=\max\big((\psi+m)/2,b\big),\qquad
 F_2=\max\big((2-3\psi+m)/2,1-3a\big).
\]
There is also a congruence estimate with exponent
\[
 F_3=a+1-3b+\max(0,b-3a).
\]
Indeed, fix \(d,v\). Since \(d\) is odd and the shift has absolute
value one or two, \(v\) is a unit modulo \(d^3\) on any occupied fiber.
The congruence \(v e^3\equiv s\bmod d^3\) has at most \(3^{\omega(d)}\)
unit residue classes. In each class the number of \(e\) is
\(O(1+E_0/D_0^3)\), and then \(u\) is fixed. There are
\(O(D_0 x/E_0^3)\) choices of \(d,v\).

Weighting by \(h\asymp x^\psi\) costs one factor \(x^\psi\).
For \(1/3\le\psi\le1/2\) use \(F_1\). The linear term gives
\(\psi+b\le5/6\), and \((3\psi+m)/2\) increases to \(57/64\) at
\(\psi=1/2\). For \(\psi\ge1/2\), we have \(b\le3a\), since
\(b\le1/3\) and \(a\ge1/6\). Therefore
\[
 \begin{split}
 \psi+\min(F_2,F_3)
 &\le \max\left(\psi+\frac{2-3\psi+m}{2},\
                   \psi+\min(1-3a,1-3\psi+4a)\right)\\
 &\le \max\left(\psi+\frac{2-3\psi+m}{2},\
                    1-\frac{2\psi}{7}\right)\le\frac{57}{64}.
 \end{split}
\]
For the second line the two affine functions of \(a\) cross at
\(a=3\psi/7\), where their minimum is largest. The first expression
decreases from \(57/64\) at \(\psi=1/2\), and \(1-2\psi/7\le6/7\).
Summing the boxes gives the last term of (6).

We have \(P\ll W^\alpha\). On division by \(W\), (6) and summation of
the \(O(\log W)\) trace ranges give a reciprocal rate
\(O_\alpha(W^{-\eta})\) over \(W\le D<2W\) for any
\[
 0<\eta<\min\left(\frac12,\frac32-\frac43\alpha,1-\frac{57}{64}\alpha\right).
\]
Here \(64/57<9/8\), so all three terms in the minimum are positive
in the stated range. Summing the dyadic tail gives
\(\sum_{D>Z,\ \epsilon_D\le D^\alpha}1/\tau_D=O_\alpha(Z^{-\eta})\)
for a fixed such \(\eta\). In (2), sum the floor errors for
\(D\le Z=Y^{1/(1+\eta)}\) trivially and use \(Y\) times the tail for
larger \(D\). This proves (5) with \(\delta_\alpha=\eta/(1+\eta)\).

Thus for every fixed \(\alpha<64/57\) in the stated range, the unresolved
rate lies among kernels with \(\epsilon_D>D^\alpha\) and with
\(\gcd(D,U_D)>B(\log D)^b\) for any fixed \(B>0,\ 0\le b<1/2\) from
(4). The full rate and the polynomial count in (2) remain open.

## Extending the small-unit range by fixing square factors

**Theorem (Small-unit asymptotic with the square factors fixed).**
The logarithmic asymptotic (5) holds for every fixed
\(1/2\le\alpha\le500/437\), with an error of order
\(O_\alpha(Y^{1-\delta})\) for some absolute \(\delta>0\).

Use the trace and kernel boxes \(T_D\asymp P\) and \(D\asymp W\)
from (6), and keep the notation \(h\) for the odd part of
\(\gcd(D,U_D)\). For \(h\le P^{1/4}\) the argument for small
\(h\) in (6), with the smaller cutoff, bounds the weighted count by
\[
 W^\varepsilon\left(\sqrt W+\frac{P^{5/4}}{\sqrt W}\right).
\]
For \(P^{1/4}<h\le4P^{1/3}\) use \(T_D^2\equiv1\bmod h^3\).
A dyadic block \(h\asymp K\) contributes
\(O_\varepsilon(P^\varepsilon(P/K+K^2))\), by summing the at
most \(2^{\omega(h)}\) residue classes for each odd squarefree \(h\).
The total in this range is \(O_\varepsilon(P^{3/4+\varepsilon})\).
A trace gives a unique squarefree part of \(T_D^2-1\), so this
count needs no multiplicity factor for \(D\).

For \(h>4P^{1/3}\), split \(h=de\) across the two coprime factors
as in (6). Keep their remaining square factors explicit, writing
\[
                 e^3 v^2 g-d^3 u^2 f=s,\qquad s\in\{\pm1,\pm2\}.
\]
Here \(D=defg\) with \(f,g\) squarefree and \(uv\asymp U_D/h\).
For odd \(T_D\) this factorization uses \((T_D-1)/2,(T_D+1)/2\).
For even \(T_D\) it uses the odd integers \(T_D-1,T_D+1\).
The two products on the left are coprime in both cases.
Write \(d,e,u,v,f,g\) on dyadic scales
\(P^a,P^b,P^c,P^{d_*},P^i,P^j\), respectively. Orient the cube
factors so \(a\le b\) and put
\[
 \lambda=\frac{\log W}{\log P},\qquad
 \psi=a+b,\qquad
 \sigma=c+d_*=1-\lambda/2-\psi,\qquad
 i=1-3a-2c,\quad j=1-3b-2d_*.
\]
Relations between the exponents from the dyadic sizes hold
up to \(O(1/\log P)\). All exponent inequalities below have
fixed slack. In particular \(1/3\le\psi\), \(b\le1/3\), and
\(\lambda\ge437/500\), up to errors of that size.

### Counting with the square factors fixed

Fix the integers \(u,v\). Put
\(D_0=P^a,\ E_0=P^b,\ F_0=P^i,\ G_0=P^j\).
In the range used below \(\log(D_0E_0)\) and \(\log(F_0G_0)\)
are comparable to \(\log P\). The number of tuples for these
fixed square factors is bounded by
\[
 \ll_\varepsilon P^\varepsilon
       \big((F_0G_0M_0)^{1/2}+F_0+G_0\big),\qquad
 \log M_0=\frac98\,\frac{\log(D_0E_0)\log(F_0G_0)}{\log P}.
 \tag{9}
\]
This is Reuss's determinant and lattice estimate with the
denominators of \(g/f\) in place of those of \(v^2g/(u^2f)\).
The following checks give uniformity in the two fixed integers.
[Reuss, §§2.1–2.3](sources/reuss-consecutive-1212.3150v2.pdf)

On the scales \(D_0/E_0\) and \(G_0/F_0\), respectively,
the ratios \(d/e\) and \(g/f\) satisfy \(t=C z^3+O(P^{-1})\)
in the normalized coordinates \(z,t\), with \(C\asymp1\).
Indeed, the original ratio error is
\(s/(v^2 f e^3)\), of order \((G_0/F_0)/P\).
Divide the first normalized coordinate into intervals of length
\(M^{-1}\), with \(M=M_0P^{\varepsilon'}\) for a sufficiently small
positive \(\varepsilon'\). Reuss's Lemma 7 calculation applies to
these two scales. In a rectangle of monomials of bidegrees at
most \(A,B\), the denominator-clearing cost after normalizing
the columns is bounded, for \(n=(A+1)(B+1)\) rows, by a constant
times \((D_0E_0)^{nA/2}(F_0G_0)^{nB/2}\). The monomials in the
interval coordinate and the vertical error have weights
\(i_0\log M+j_0\log P\). The first \(n\) such weights sum to
\(\tfrac{2\sqrt2}{3}n^{3/2}\sqrt{\log M\log P}+O(n\log P)\),
by the elementary triangle count in §2.1 of Reuss. Taking
\(A/B\) close to \(\log(F_0G_0)/\log(D_0E_0)\) and taking
both degrees large in terms of \(\varepsilon'\) makes the
cleared determinant smaller than one. Its integral entries
and the cofactor construction give a nonzero auxiliary form
of bounded bidegree and coefficient height \(P^{O_{\varepsilon'}(1)}\).

On each interval change basis in the lattice for \((g,f)\) as
in Reuss, §2.2. Write \(T_1\) for the larger coordinate bound
in the reduced lattice box. The primitive pair \((g,f)\)
stays primitive in the lattice basis, and \((d,e)\) is primitive.
An absolutely irreducible component of the auxiliary form
of degree at least two in \((d,e)\) has
\(O_\varepsilon(P^\varepsilon T_1)\) points, using Heath-Brown's
Lemma 2 with the two pairs exchanged. The degree and coefficient
bounds are uniform. An irreducible rational component that is
geometrically reducible contributes only its rational points in
intersections of conjugate components, bounded by Bézout.
[Heath-Brown, Lemma 2](https://www.impan.pl/shop/en/publication/transaction/download/product/82579)

For degree one in \((d,e)\), write
\((d,e)=\pm(A_1(g,f),A_2(g,f))/q\) with coprime homogeneous
binary forms \(A_1,A_2\) and \(q\) their positive value gcd.
The homogeneous resultant identities in the two projective
charts and primitivity of \((g,f)\) force \(q\) to divide a
fixed nonzero integer of size \(P^{O(1)}\) for this component.
There are \(O_\varepsilon(P^\varepsilon)\) choices for \(q\).
For each choice and sign, substitution gives an equation of
bounded positive homogeneous degree in the lattice coordinates
of \((g,f)\), with nonzero right side \(\pm s q^3\).
If the homogeneous form vanishes identically there are no
solutions. Otherwise the equation has \(O(T_1)\) integer
points in the reduced box, by fixing one coordinate and
counting polynomial roots, including the bounded number of
exceptional coordinate values.

A component independent of \((d,e)\) fixes \(O(1)\) primitive
pairs \((g,f)\). For each such pair the possible fractions
\(d/e\) lie in an interval of length \(O((D_0/E_0)/P)\).
Their minimum spacing is \(\gg E_0^{-2}\), giving \(O(1)\)
choices because \(D_0E_0\ll P^{2/3}\). A component independent
of \((g,f)\) similarly fixes \(O(1)\) pairs \((d,e)\).
For each such pair the linear equation in \(f,g\) has
\(O(1+P/((D_0E_0)^3u^2v^2))=O(1)\) solutions.
These facts also show why the original coprime equation has
bounded multiplicity over each of the fixed ratios in the
degree-zero cases.

Finally the normalized centers for \(g/f\) arise within
\(O(M^{-1})\) of grid points of spacing \(M^{-1}\), each
used \(O(1)\) times because the derivative of the normalized
cubic is bounded away from zero. Reuss's §2.3 sum therefore
bounds \(\sum T_1\) by
\(O_\varepsilon(P^\varepsilon((F_0G_0M)^{1/2}+F_0+G_0))\).
The bound includes the two axial shortest-vector cases.
The factors of the auxiliary form have fixed degrees and
polynomial height, and the lattice basis change has
polynomially bounded entries. Absorbing \(\varepsilon'\)
into \(\varepsilon\) proves (9).

### Weighted count and summation

For \(1/3\le\psi\le241/500\), keep the first Reuss estimate
from (6). The weighted exponent is at most
\[
 \max\left(\psi+b,\frac{3\psi+(9/8)\psi(2-3\psi)}2\right)
 \le 0.87320325<437/500.
\]
The second entry is increasing on this interval.

For \(\psi\ge\psi_0=241/500\), sum (9) over the
\(O(P^\sigma)\) choices for the fixed squares and multiply
by the weight \(h\asymp P^\psi\). The resulting exponent is
\(1-\lambda/2\) plus the largest of
\[
 F_{\rm lat}=\frac{(\lambda-\psi)(1+9\psi/8)}2,\qquad i,\qquad j .
\]
Here \(0<\lambda-\psi\asymp1\), and
\(\psi\le1-\lambda/2\le563/1000\).
Set \(T_*=3\lambda/2-1\) and \(\eta=1/10000\).
The entry \(j\le1-3\psi_0/2=0.277\) is below \(T_*-\eta\).
The lattice entry also lies below \(T_*-\eta\).
The difference \(T_*-F_{\rm lat}\) increases with \(\lambda\)
on this range, and at \(\lambda=437/500\) the lattice entry
is decreasing in \(\psi\) and equals \(0.302281\) at \(\psi_0\).
At that endpoint \(T_*=0.311\).

If \(i>T_*-\eta\), use a congruence instead. For fixed
\(u,v,d,g\), the equation determines \(e\) in
\(O_\varepsilon(P^\varepsilon)\) classes modulo \(d^3u^2\),
since \(g v^2\) and \(s\) are units there and a unit has at
most \(3^{\omega(du)}\) cube roots. In the \(s=\pm2\) case
all the variables are odd, so two causes no exception.
Here \(a\ge\psi_0-1/3\), so \(b<3a\), and each class yields
\(O(1)\) choices. The exponent before the weights and fixed
square choices is therefore \(a+j\). From
\(i=1-3a-2c>T_*-\eta\) and \(c+d_*=\sigma\), we find
\[
 \begin{split}
  a+j &< -b+1-\lambda/2+\eta\\
      &< 5/3-\lambda-\psi_0+4\eta/3
       \le T_*-\eta .
 \end{split}
\]
For the second inequality use
\(a<2/3-\lambda/2+\eta/3\). For the last, use
\(\psi_0-(8/3-5(437/500)/2)=1/3000>7\eta/3\).
Thus in all cases the weighted exponent is at most
\(\lambda-\eta\) in this range, before arbitrarily small
losses from the determinant estimate and dyadic summation.

Combining the ranges gives a power-saving bound for the
weighted sum over \(D\asymp W\) when \(P\ll W^{500/437}\).
Indeed the terms \(\sqrt W\), \(P^{5/4}/\sqrt W\), and
\(P^{3/4}\) are each \(O(W^{1-\eta})\).
The final range contributes at most \(W P^{-\eta}\)
before the small losses, and \(P\gg\sqrt W\).
Since \(\tau_D\ge D\log\epsilon_D/\gcd(D,U_D)\),
some \(\eta_0>0\) gives a reciprocal rate tail
\(O_\alpha(Z^{-\eta_0})\) for \(D>Z\) in this unit range.
The floor truncation in (5) now proves the asserted
asymptotic. The rate of the remaining large-unit kernels
and the full polylogarithmic count remain open.

## Pell solutions on a fixed algebraic curve

**Proposition (Logarithmic contribution of fixed curves).**

Let \(C\) be a fixed affine algebraic curve over \(\mathbb Q\) in the
Pell surface \(T^2-DU^2=1\). Write \(\mathcal D_C\) for the squarefree
\(D>1\) for which \(C\) contains any positive integral Pell solution.
Then
\[
 A_C(Y)=\#\{D\in\mathcal D_C:\tau_D\le Y\}\ll_C\sqrt Y,\qquad
 \sum_{D\in\mathcal D_C}\frac1{\tau_D}=S_C<\infty.
                                                                    \tag{7}
\]
The number of pairs in (2) contributed by all these kernels is
\(S_CY+O_C(\sqrt Y)\). A finite union of fixed curves has the same
bounds.

It suffices to work with a geometrically irreducible curve with
nonconstant \(D\), infinitely many admissible integral points, and
\(U\) not identically zero. Constant \(D\) and finite sets of points
give finite sets of kernels. For a \(\mathbb Q\)-irreducible curve
that is geometrically reducible, rational points lie on intersections
of conjugate components and are finite. Discarding the finite singular
points does not affect the bounds.

By [Alvanos–Bilu–Poulakis, Theorem 1.2](https://arxiv.org/pdf/0907.2097v1),
the smooth projective model has genus zero and either one pole of the
coordinate functions or two poles conjugate over a real quadratic
field. The nonsingular rational point guarantees that the model is
\(\mathbb P^1\) over \(\mathbb Q\). We use the orders of the zeros
of \(T^2-1=DU^2\) to bound \(\gcd(D,U)\) in both cases. At any point
away from the poles where \(U\) has a zero of order \(v\ge1\) and \(D\)
has order \(d_0\ge0\), one has \(T=\pm1\) there. The derivative of \(T\)
in a local coordinate then has a zero of order \(d_0+2v-1\), at least
\(v+\min(d_0,v)\). This is the same critical-point calculation used
for polynomial Pell–Abel equations, applied here to the common
divisor of the two coordinate functions. [Bogatyrev–Gendron, §1(A)](https://arxiv.org/html/2306.00884v2)

**One pole.** Put the rational pole at infinity with parameter \(x\).
Now \(D,T,U\in\mathbb Q[x]\). Let \(d=\deg D>0,\ u=\deg U\), and let
\(G\) be the polynomial gcd of \(D,U\). The degree identity
\(\deg T=u+d/2\) makes \(d\) even and at least two. The local
calculation gives \(UG\mid T'\), so \(\deg G\le d/2-1\).
The denominators of the rational parameters of the integral points
are bounded. In fact if a valuation of \(x\) is sufficiently negative,
the leading term of the nonconstant polynomial \(D(x)\) has strictly
smallest valuation, which is negative. Outside the finite set of
primes of the coefficients and the leading coefficient this holds
as soon as \(x\) has a denominator at that prime. Thus all parameters
lie in \(M^{-1}\mathbb Z\) for a fixed integer \(M>0\).

Evaluate a polynomial Bézout identity \(G=F_1D+F_2U\) at a parameter
of an integral point. A fixed denominator clears the values of
\(F_1,F_2\), so \(\gcd(D(x),U(x))\ll_C |G(x)|\). Here \(G(x)\ne0\)
because \(D(x)>1\). For large \(|x|\),
\(D(x)\asymp_C |x|^d\) on the positive branch and
\(\gcd(D(x),U(x))\ll_C|x|^{d/2-1}\).
Every positive Pell solution is a power of the least one, and \(U_D\)
divides its second coefficient. Hence \(R_D=D/\gcd(D,U_D)\) is
\(\gg_C |x|^{d/2+1}\). As \(\tau_D\gg R_D\), the number of parameters
with \(\tau_D\le Y\) is \(O_C(Y^{1/(d/2+1)})\), which is
\(O_C(\sqrt Y)\).

**Two conjugate poles.** With \(K\) their quadratic field, use a
coordinate \(z\) taking the poles to zero and infinity, with conjugation
sending \(z\) to \(z^{-1}\) at rational points. For example if the poles
in a rational coordinate \(x\) are \(\omega,\overline\omega\), take
\(z=(x-\omega)/(x-\overline\omega)\). The functions \(D,T,U\) are
Laurent polynomials in \(z\) over \(K\). Their pole orders at zero and
infinity are equal for each function by conjugation. Denote the orders
at infinity by \(d,m,u\) respectively, so \(d\ge1,\ u\ge0\) and
\(m=u+d/2\). Let \(G\in K[z]\) be the monic gcd in the Laurent ring
normalized to have nonzero constant term. Write \(g=\deg G\). Both
poles of \(z\,dT/dz\) have order \(m\), giving \(2m\) zeros away from
the poles. The local calculation shows that these zeros include the
\(2u\) zeros of \(U\) and all \(g\) zeros of \(G\), with multiplicities
added. Consequently \(g\le2m-2u=d\).

The finite-place valuations of the possible \(z\) are bounded and
vanish outside a fixed finite set of primes of \(K\). At each finite
place use a Laurent polynomial having a pole at zero to rule out
sufficiently positive \(v(z)\), by dominance of its lowest term, and
one having a pole at infinity to rule out sufficiently negative
\(v(z)\). Thus the principal fractional ideals \((z)\) range over a
finite set, and the values of \(z\) lie in finitely many cosets of the
unit group of \(K\). If \(r=\max(|z|,|z|^{-1})\) at a real embedding,
each coset has \(O_C(\log(2Z))\) values with \(r\le Z\), by the rank-one
unit group of a real quadratic field. Also \(D(z)\asymp_C r^d\) for
large \(r\) on the positive branch.

The Laurent Bézout identity \(G=F_1D+F_2U\) can be evaluated after
clearing the denominators of \(F_1(z),F_2(z)\) with one fixed rational
integer \(M\). This follows from the finite-place bounds just proved.
The algebraic integer \(MG(z)\) is divisible in \(\mathcal O_K\) by
\(H=\gcd(D(z),U(z))\). Its norm is nonzero and divisible by \(H^2\).
One argument of \(G\) and its conjugate is bounded in absolute value
by one and the other has absolute value \(r\). Therefore
\[
 H^2\le \big|N_{K/\mathbb Q}(MG(z))\big|\ll_C r^g,\qquad
 R_D\ge D(z)/H\gg_C r^{d/2}.
\]
This gives \(A_C(Y)\ll_C\log(2Y)\) in the two-pole case.

In both cases, partial summation of \(A_C(Y)\ll_C\sqrt Y\) gives
\(\sum_{\tau_D>Y,\ D\in\mathcal D_C}1/\tau_D=O_C(Y^{-1/2})\).
In (2) the floor errors for \(\tau_D\le Y\) are at most \(A_C(Y)\),
while the omitted linear terms are bounded by \(Y\) times the rate
tail. This proves (7) and its error estimate. Over a finite union use
the union of the distinct kernels in defining the rate.

**Corollary (Zariski density of first-pair Pell triples).**

For each squarefree \(D>1\), take the Pell triple with exponent \(L_D\)
in (2), which gives the first powerful pair for that kernel. These
triples are Zariski dense in \(T^2-DU^2=1\). In fact on each fixed curve
with nonconstant \(D\), the preceding proof gives
\(\gcd(D,U)\ll_C\sqrt D\) on its positive integral points. Since the
first-pair triples satisfy \(D\mid U\), only finitely many kernels can
occur on that curve. The same conclusion holds for curves with
constant \(D\). A proper algebraic subset of the irreducible Pell
surface is a finite union of curves and points, so it cannot contain
the first-pair triples for infinitely many \(D\). The same applies
over \(\mathbb C\), since expanding any vanishing polynomial in a
\(\mathbb Q\)-basis of its coefficients yields rational vanishing
polynomials, at least one of which is nonzero on the surface.

A fixed finite union of algebraic curves therefore contributes at
most a logarithmic count. A growing family requires a bound uniform
in the curves and their number.

## Fixed fibers from a trace equation

**Proposition (Uniform fixed fiber).**

Let \(d,v\ge1\), with \(d\) squarefree, and let
\(\sigma=\pm1\). For \(K=\mathbb Q(\sqrt{-\sigma d})\) quadratic write
\(h_3(K)=|\operatorname{Cl}(\mathcal O_K)[3]|\). Use \(h_3=1\) in
the split case \(d=1,\sigma=-1\). The number of positive integer
solutions \((u,e)\) of
\[
 d^3u^2+\sigma=v^2e^3
\]
is \(O(2^{\omega(v)}h_3(K))\) with an absolute constant. In particular it
is \(O_\varepsilon(d^{1/2+\varepsilon}v^\varepsilon)\).

In the quadratic case put \(\alpha=1+du\sqrt{-\sigma d}\). Its trace
is \(2\) and its norm is \(\sigma v^2e^3\). Factor its nonzero principal
ideal as \((\alpha)=\mathfrak D\mathfrak A^3\) with cube free integral
ideal \(\mathfrak D\) and integral ideal \(\mathfrak A\). A prime above
an odd rational prime cannot divide both \(\alpha\) and its conjugate,
by the trace identity. Hence an odd prime dividing \((\alpha)\) lies
above a split prime \(p\), and at most one of the two primes above
\(p\) occurs. Its exponent is \(2v_p(v)+3v_p(e)\) by the norm identity.
The residue modulo \(3\) is fixed by \(v\), with at most two possible
choices of prime when this residue is nonzero. There are at most nine
choices of residual exponents at the primes above \(2\). Thus there
are at most \(9\cdot2^{\omega(v)}\) possibilities for \(\mathfrak D\).

For each \(\mathfrak D\), the class of \(\mathfrak A\) solves
\([\mathfrak A]^3=[\mathfrak D]^{-1}\). Each nonempty fiber of this map
has \(h_3(K)\) elements. Fix an integral ideal \(\mathfrak J\) for
one such class and a generator \(\gamma\) of
\(\mathfrak D\mathfrak J^3\). Then \(\mathfrak A=\beta\mathfrak J\) for
some \(\beta\in\mathfrak J^{-1}\), and
\(\alpha=\eta\gamma\beta^3\) for a unit \(\eta\).
An arbitrary cube in the unit can be absorbed into \(\beta\), so
Dirichlet's unit theorem leaves at most nine choices for \(\eta\)
(one free generator at most and a cyclic group of roots of unity).

Write \(L(r,s)=r w_1+s w_2\) in a \(\mathbb Z\)-basis of
\(\mathfrak J^{-1}\). The trace equation becomes
\[
 F(r,s)=2,\qquad F=\operatorname{Tr}_{K/\mathbb Q}(\eta\gamma L^3).
\]
This binary cubic has integral coefficients because
\(\gamma(\mathfrak J^{-1})^3=\mathfrak D\subseteq\mathcal O_K\).
The two linear forms \(L,\overline L\) are independent over \(K\):
otherwise \(w_1/w_2\) would be rational. Over \(K\) the cubic is an
invertible linear change of a diagonal binary cubic with two nonzero
coefficients. It has distinct projective roots in characteristic zero.
Its content divides \(2\) if there is a solution, and every solution
has \(\gcd(r,s)=1\), since its cube divides \(2\).
[Stewart, the consequence of Theorem 1 on p. 796](https://uwaterloo.ca/pure-mathematics/sites/default/files/uploads/documents/s0894-0347-1991-1119199-x_0.pdf)
therefore bounds the number of pairs by an absolute constant after
dividing out the content. The theorem allows reducible forms with
nonzero discriminant. Each pair determines \(\alpha\), hence \(u\)
and at most one positive \(e\).

In the split case \(u>1\) and the two nonzero integers \(1+u,1-u\)
have sum \(2\) and product \(-v^2e^3\). Factor each as \(a_j b_j^3\)
with \(a_j>0\) cube free and \(b_j\in\mathbb Z\), retaining the sign
in \(b_j\). The same odd prime argument and the at most nine choices
at \(2\) give \(O(2^{\omega(v)})\) coefficient pairs. The equation
\(a_1b_1^3+a_2b_2^3=2\) again meets Stewart's hypotheses. This covers
both signs and the exceptional square radicand. Finally
\(h_3(K)\le h_K\ll_\varepsilon d^{1/2+\varepsilon}\) by Minkowski's bound
and counting ideals of norm up to \(O(\sqrt d)\) with the divisor
bound. Also \(2^{\omega(v)}\ll_\varepsilon v^\varepsilon\).

A related general input is
[Evertse–Silverman, Theorem 1(a), pp. 237–238](sources/evertse-silverman-superelliptic-counts-1986.pdf).
Take the base field \(\mathbb Q\), \(f(T)=(d^3T^2+\sigma)/v^2\),
and \(S\) consisting of infinity and primes dividing \(2dv\).
Its discriminant \(-4\sigma d^3/v^4\) is an \(S\)-unit and both
roots lie in \(K\) (in \(\mathbb Q\) in the split case).
That theorem gives \(O_\varepsilon((dv)^\varepsilon h_3(K))\),
which also suffices for (1). The trace proof keeps the explicit
\(2^{\omega(v)}\) dependence and avoids including primes dividing \(d\)
in the branch budget.

## Averaging over the cube kernel

For each sign the fundamental discriminant of \(K\) has absolute
value at most \(4d\), and the squarefree radicands parametrize distinct
quadratic fields. The first moment in
[Bhargava–Shankar–Tsimerman, Theorem 2](sources/bst-2013-davenport-heilbronn.pdf)
therefore gives \(\sum_{d\le D}h_3(K)\ll D\) over the squarefree
kernels, including the split convention. For any \(D,V\ge1\), the
number of solutions with \(d\le D,\ v\le V\) is consequently
\(O_\varepsilon(DV^{1+\varepsilon})\) for either sign. For fixed
\(D\) and \(V=X^{1/8}\) this gives the bounded-kernel
\(O_{D,\varepsilon}(X^{1/8+\varepsilon})\) bound.

## Global bound from four estimates

Every positive powerful number has a unique representation \(d^3u^2\)
with \(d\) squarefree and \(u>0\). Use \(e^3v^2-d^3u^2=1\) and divide
\(n=d^3u^2\) into dyadic ranges of size comparable to \(x\). Group the
kernels at dyadic scales \(D,E\), and apply Plymesser's Lemma 2.3 to use
\(O(1)\) exact windows per group, with comparable \(D,E,x\) and
\(U=x^{1/2}D^{-3/2},\ V=x^{1/2}E^{-3/2}\).
To use open endpoints in all four variables, round \(3d/4,3e/4,n/4\)
down to the grid \(\{2^{j/8}:j\in\mathbb Z\}\) to obtain \(D,E,x\).
Then \(d/D,e/E\) are strictly between \(1\) and \(2\), and
\(u/U=(D/d)^{3/2}(n/x)^{1/2}\) lies between \(1\) and \(2\)
because \((3/4)2^{-1/8}<D/d\le3/4\) and \(4\le n/x<4\cdot2^{1/8}\).
The same bounds with the extra factor \(((n+1)/n)^{1/2}\) give
\(1<v/V<2\) for \(n\ge2\). This grid still uses \(O(1)\)
choices per dyadic group.
Write the kernel exponents as \(a,b\), and set
\[
 s=a+b,\quad t=1-3s/2,\quad m=9st/8,\qquad
 \bar u=(1-3a)/2,\quad \bar v=(1-3b)/2.
\]
A nonempty window has exponents within \(O(1/\log x)\) of the square
\([0,1/3]^2\). Clipping to the square changes every expression below
by \(O(1/\log x)\), which costs only a constant factor. Assume
\(a\le b\) by exchanging the pairs and using the corresponding sign
in the fiber. The determinant input can always be applied to the
original orientation with shift \(1\) since all of its exponents
are symmetric. Each counted solution has \(\gcd(d,e)=\gcd(u,v)=1\),
by the shift equation.

The available box count exponents, up to \(\varepsilon\), are

1. \(s\) and \(t\) by Plymesser's Lemma 4.1 with shift \(1\). The
   two estimates fix the kernel pair and the square variable pair
   respectively. The first is the uniform quadratic norm count in a
   polynomial height range. The second uses Stewart on the binary cubic
   with right-hand side \(1\), whose coefficients are coprime since
   \(\gcd(u,v)=1\).
2. \(a+\bar v\) by the averaged fixed fiber, with exponent \(O(\varepsilon)\)
   absorbed by applying that estimate with smaller \(\varepsilon\).
3. \(\max(m,B)\) by Plymesser's Theorem 6.3 when the actual
   \(s,t\ge 1/10\), where
   \[
     B=m+\tfrac3{10}\max((s-m)/2,a-m,b-m)
          +\tfrac15\max((t-m)/2,\bar u-m,\bar v-m).
   \]

The third input uses Reuss's **Lemma 7**, with \(k=3,l=2,h=1\),
and \(M=\lceil x^{(1+\delta)m}\rceil\) for a fixed sufficiently small
\(\delta>0\). In the region where it is needed below,
\(s,t\ge 6/25-O(1/\log x)\), so \(\log M\asymp\log x\),
\(M\le x^{1/2}\) and \(x\ge 3M\) for large \(x\).
The construction in Sections 5–6 uses the two reduced lattices and
primitive pairs with gcd set \(\{1\}\). The fixed ratios use Lemma 4.1.
The bidegree-one cases use Proposition 5.3 with bounded bidegree and
polynomial coefficient height. The remaining components have bidegree
\((\alpha,\beta)\) with \(\alpha,\beta\ge2\). Lemma 5.5 embeds their
actual integral coordinates in dimension 12 on an affine curve of
degree \(2\alpha+3\beta\ge10\). This is the input for
[Castryck–Cluckers–Dittmann–Nguyen, Theorem 3 and Remark 4.1.5](https://arxiv.org/pdf/1904.13109v1).
The coefficient height bound uses Gelfond's inequality as in
Plymesser's Theorem 3.5. The first minimum tails of Lemmas 6.1–6.2
are at this same \(M\), and the two moment exponents after
Cauchy–Schwarz are \(3/5,2/5<2\). This application uses the
reduced count of Theorem 6.3 and Lemma 4.1. It needs neither the
separate Reuss box bound in Theorem 3.1 nor the anisotropic estimate
of Theorem 3.4 or the optimization in Proposition 7.1(ii).

### Exact optimization

Let \(\theta=12127/49920\). If any of \(s,t,a+\bar v\) is at most
\(\theta\), use that estimate. Otherwise set \(\Delta=b-a\ge0\).
The fiber inequality gives
\[
 a+\bar v=\tfrac12-\tfrac14s-\tfrac54\Delta>\theta,\qquad
 \Delta<(2-s-4\theta)/5 .
\]
Here \(6/25<\theta\) and \(s,t>\theta\), so
\(6/25\le s\le38/75\). On this interval \(m\ge(26/25-s)/5\):
the difference is concave, and it is positive at both endpoints
by substitution. Thus \(\Delta\le m\).
The first maximum in \(B\) is \((s-m)/2\) and the second is
\((t-m)/2+\tfrac12\max(0,3\Delta/2-m)\). Since \(3s+2t=2\),
\[
 B=\max\left(\tfrac1{10}+\tfrac34m,\;
              \tfrac1{10}+\tfrac{13}{20}m+\tfrac3{20}\Delta\right).
\]
The first term is at most \(77/320<\theta\), as \(m\le3/16\).
The second is at most
\[
 \tfrac4{25}-\tfrac3{25}\theta+\tfrac{561}{800}s-\tfrac{351}{320}s^2
 \ \le\ \tfrac4{25}-\tfrac3{25}\theta+
             \frac{(561/800)^2}{4(351/320)}=\theta .
\]
The last equality defines the displayed rational exponent. Also
\(m\le3/16<\theta\). There are \(O((\log X)^3)\) windows over
all dyadic sizes and kernel scales. Their count gives (1) by
taking a smaller exponent loss in each window and absorbing the
logarithms.

## Internal review of the subquarter argument

At an odd prime, a self-conjugate prime ideal dividing \((\alpha)\)
would divide its conjugate and hence \(2\). This rules out both the
inert and ramified cases. For a split prime with nonzero cube residue
the choice is just which of its two primes occurs. The trace-form
change of variables has nonzero determinant because the ideal basis
has irrational ratio. Thus the form has nonzero discriminant before
primitivization, and its content divides 2 on an occupied fiber.
Stewart's primitive-form bound therefore applies with right-hand side
1 or 2, including when the form is reducible.

In the cached Plymesser source, §2.3 uses \(D,E\ge1/2\) and \(x>1\)
with \(U,V\) determined by those scales. The grid above has \(D,E>1/2\)
for positive kernels and \(x>1\) once the finitely many small \(n\) are
discarded. Theorem 6.3 permits each fixed \(h,c,\mathcal G\). The
application uses \(h=1\), \(c=1/10\) and \(\mathcal G=\{1\}\), justified
before applying the determinant estimate. Clipping does not change that
eligibility in the remaining region. The fixed-trace average handles
both assignments of the smaller kernel while the determinant window
keeps its original orientation. In the exact optimization the
quadratic's stationary point \(s=187/585\) lies in the allowed interval;
substitution gives the stated exponent, with \(1/4-\theta=353/49920\).
See the [primary source audit](PLYMESSER_2026_POWERFUL_COUNTS_SOURCE_SCOPE.md)
for the imported geometric and counting estimates. These are internal
checks. Independent review of the trace-fiber argument and the complete
original-width window theorem remains open.

## A degree-one projection for larger Pell units

**Theorem (Larger small-unit range for the Pell rate).**
Set
\[
 \lambda_*=\frac{172+\sqrt{1153}}{243},\qquad \alpha_*=\lambda_*^{-1}.
\]
For every fixed \(1/2\le\alpha<\alpha_*\), there is \(\eta_\alpha>0\)
such that
\[
 \sum_{\substack{W\le D<2W,\ D\ \mathrm{squarefree}\\
                  \epsilon_D\le D^\alpha}}\gcd(D,U_D)
        \ll_\alpha W^{1-\eta_\alpha}.                          \tag{365.L1}
\]
In particular this holds at \(\alpha=125/106\). The kernels in
this unit range have a finite summed rate \(S_\alpha\) and contribute
\(S_\alpha Y+O_\alpha(Y^{1-\delta_\alpha})\) consecutive powerful
pairs in (5), for some \(\delta_\alpha>0\).

Fix the two square factors of a Pell trace. On a short interval the
ratio of the cube factors satisfies an auxiliary equation that is
linear in that ratio. Substituting its rational expression into
the original shifted equation gives a binary form. The cubic powers
in that form prevent it from being a power of a single linear factor.
This gives a subpower count on each such interval.

### The determinant with six columns

Use \(T_D\asymp P\), \(D\asymp W\), and \(h\) the odd part of
\(\gcd(D,U_D)\), as in the preceding estimates. The contribution with
\(h\le4P^{1/3}\) to the sum of \(h\) is bounded there by
\[
 W^\varepsilon\left(\sqrt W+P^{5/4}/\sqrt W\right)
        +O_\varepsilon(P^{3/4+\varepsilon}).
\]
For \(h>4P^{1/3}\), use the same factorization
\[
 e^3v^2g-d^3u^2f=s,\qquad s\in\{\pm1,\pm2\},\qquad D=defg,\qquad h=de .
\]
The two terms on the left are coprime and of order \(P\).
Write \(D_0,E_0,F_0,G_0\) for the dyadic sizes of \(d,e,f,g\).
Put \(\lambda=\log W/\log P\) and
\(\psi=\log(D_0E_0)/\log P\), so that
\(F_0G_0\asymp P^{\lambda-\psi}\).
In a box the total number of pairs \(u,v\) is
\(O(P^{1-\lambda/2-\psi})\), since its dyadic product has
this order.

Fix \(u,v\). In the normalized coordinates
\(z=(d/e)/(D_0/E_0)\), \(t=(g/f)/(G_0/F_0)\), one has
\(t=Cz^3+O(P^{-1})\) with \(C\asymp1\) in an occupied box.
Both \(t\) and \(z\) range over compact positive intervals. In
particular \(z=C^{-1/3}t^{1/3}+O(P^{-1})\). In an interval
\(|t-t_0|\le\Delta\), take six rows of the evaluation matrix
for
\[
                         1,\ t,\ t^2,\ z,\ zt,\ zt^2 .
\]
Write \(x=t-t_0\), and replace the \(t^j\) and \(zt^j\) columns
by \(x^j\) and \(zx^j\), respectively. In the three columns
containing \(z\), expand each as a smooth column plus an
error column with entries \(O(P^{-1}|x|^j)\), for \(j=0,1,2\).
A term with \(k\) error columns has order at most
\[
          P^{-k}\Delta^{(6-k)(5-k)/2+k(k-1)/2}.
\]
Indeed, expanding in the error columns leaves a smooth
minor of size \(6-k\), whose determinant has the factor
\(\Delta^{(6-k)(5-k)/2}\) by divided differences.
The distinct error-column indices have sum at least
\(k(k-1)/2\). This argument includes repeated abscissae,
for which a smooth minor with two identical abscissae vanishes.

Before normalization the rows, multiplied by \(ef^2\),
have integral entries. Restoring the column scales therefore
bounds the cleared determinant by a constant times
\[
 P^{6\lambda-3\psi}
       \left(\Delta^{15}+P^{-1}\Delta^{10}
                       +P^{-2}\Delta^7+P^{-3}\Delta^6\right).
                                                                  \tag{365.L2}
\]
Here the clearing factor is \(O((D_0E_0)^3(F_0G_0)^6)\).
If this bound is less than one, the integer matrix has
rank at most five. A nonzero coefficient vector of
polynomial height in \(P\) follows by selecting a maximal
nonzero minor and using cofactors. The interval therefore
has a nonzero auxiliary form of bidegree at most \((1,2)\)
in \((d:e),(g:f)\).

### Counting on the auxiliary equation

A fixed primitive ratio \(g/f\) leaves \(d/e\) in an interval
of length \(O((D_0/E_0)P^{-1})\), by the original equation.
Distinct reduced fractions in the box are separated by
\(\gg E_0^{-2}\). Thus there are \(O(1+D_0E_0/P)=O(1)\)
possibilities, using \(d,e\ll P^{1/3}\).
For fixed \(d/e\) the equation is linear in \(f,g\), with
coprime coefficients by the coprimality of the two shifted
terms. Its solution count is at most
\(O(1+P/((D_0E_0)^3u^2v^2))=O(1)\), using
\(D_0E_0\asymp h>4P^{1/3}\). These observations handle the
auxiliary forms depending on only one ratio and any
factors in the second ratio alone.

After removing the common homogeneous factors of the two
coefficients of \((d,e)\), write
\[
                         (d,e)=\pm(A_1(g,f),A_2(g,f))/q ,
\]
where \(A_1,A_2\) are coprime homogeneous integer forms of
equal degree \(k\le2\), and \(q=\gcd(A_1(g,f),A_2(g,f))>0\).
If either form is zero, the original auxiliary equation
fixes finitely many ratios \(g/f\) or admits no positive pairs.
Degree \(k=0\) fixes \(d/e\). The forms have polynomial
coefficient height in \(P\), since cancellation,
denominator clearing and factoring take place in bounded
degree. The identities from the resultants on the two
affine charts show that \(q\) divides a fixed nonzero
integer of polynomial size in \(P\): one identity gives
\(q\mid R_1 f^{j_1}\), the other \(q\mid R_2 g^{j_2}\), with
\(R_1,R_2\) nonzero, and \(\gcd(f,g)=1\).
There are \(O_\varepsilon(P^\varepsilon)\) choices of \(q\).
For \(k\ge1\), substitution gives
\[
 B(g,f)=\pm s q^3,\qquad
 B(X,Y)=v^2 X A_2(X,Y)^3-u^2 Y A_1(X,Y)^3.                       \tag{365.L3}
\]

To check the form, cancel the common homogeneous factors in
the ratio \(v^2 X A_2^3/(u^2 Y A_1^3)\).
Since \(A_1,A_2\) have no common projective zero, any
cancellation is at \(X=0\) or \(Y=0\) with multiplicity one.
If there are \(c\) such cancellations, the resulting map of
projective lines has degree \(3k+1-c>0\), and its fibers
over \(0,\infty\) together have at most \(2k+2-c\) distinct
points. In particular \(B\ne0\). If \(B\) were a power of
one linear form over \(\mathbb C\), the fiber over \(1\)
would have just one point. But for a degree-\(d\) map of
projective lines, the total number of distinct points over
\(0,1,\infty\) is at least \(d+2\), by Riemann–Hurwitz. Here
\(2k+3-c<3k+3-c\), a contradiction.

The number of primitive solutions of (365.L3) in the dyadic
box is \(O_\varepsilon(P^\varepsilon)\) for each \(q\) and
sign. If \(B\) has two distinct rational irreducible factors,
fix the nonzero value of each of those two primitive factors.
Each value divides the right side, and the two resulting
affine curves intersect in a bounded number of points by
Bézout. They have no common component because their top
homogeneous forms are coprime. If \(B\) is a constant times
a power of a single irreducible form, fix the value of that
form, which again divides the right side when the factors
are taken primitive. For degree at least three the form has
nonzero discriminant, so Stewart's primitive Thue bound
gives \(O_\varepsilon(P^\varepsilon)\) solutions after
dividing out content. Degree one was excluded above.
For a quadratic form \(aX^2+bXY+c'Y^2\) of discriminant
\(\Delta_0=b^2-4ac'\) nonsquare, a solution of value \(m\ne0\)
gives the quadratic integer \(2aX+(b+\sqrt{\Delta_0})Y\)
of norm \(4am\). There are at most \(\tau(|4am|)\) integral
ideals of that norm. In a real quadratic field the two
embeddings have modulus bounded by a power of \(P\), and the
norm gives an inverse polynomial lower bound. The logarithms
of the absolute values of real units are spaced by a uniform
positive constant, so at most \(O(\log P)\) generators of
each ideal occur. In an imaginary quadratic field at most six
units occur. The transformation is injective, giving the
same subpower estimate. All size bounds are uniform because
the coefficients, variables and right sides have polynomial
height in \(P\).
[Stewart, consequence of Theorem 1 on p. 796](https://uwaterloo.ca/pure-mathematics/sites/default/files/uploads/documents/s0894-0347-1991-1119199-x_0.pdf)

### Cutoffs and the rate

For an explicit cutoff take \(\lambda\ge106/125\) and split
at \(\psi_0=457/1000\), allowing \(O(1/\log P)\) changes in the
dyadic exponents. At \(\psi\le\psi_0\), use Reuss's first
estimate from (6), with \(m=(9/8)\psi(2-3\psi)\).
The weighted exponent is at most
\[
       \max\!\left(\psi_0+\frac13,\,
                     \frac{21}{8}\psi_0-\frac{27}{16}\psi_0^2\right)
                  < \frac{106}{125}-\frac1{10000}.
\]
The quadratic is increasing on this interval of \(\psi\).
At \(\psi\ge\psi_0\), set
\(\mu=3\lambda/2-1-\eta\), where \(\eta=1/10000\),
and partition the \(t\)-range into \(O(P^\mu)\) intervals
of length \(O(P^{-\mu})\). At
\((\lambda,\psi)=(106/125,457/1000)\), the four decay
exponents \(15\mu,1+10\mu,2+7\mu,3+6\mu\) in (365.L2)
exceed \(6\lambda-3\psi\) by at least \(1/500\).
Each excess is increasing in both \(\lambda\) and \(\psi\).
Hence the determinant vanishes for large \(P\).
Summing the per-interval count over \(u,v\) and inserting
the weight \(h\) costs \(P^{1-\lambda/2}\), giving a total
\(O_\varepsilon(P^{\lambda-\eta+\varepsilon})\).

For the stated open endpoint choose
\(\lambda_0\in(\lambda_*,106/125]\) with \(\lambda_0<1/\alpha\)
and take \(\psi_0\) just above \(3(1-\lambda_0)\).
The quadratic bound at \(3(1-\lambda_0)\) is strictly less
than \(\lambda_0\); equality defines \(\lambda_*\) in this
range. The linear bound \(\psi_0+1/3<\lambda_0\) also holds.
At \(\psi=3(1-\lambda)\) and \(\eta=0\), the one-error
exponent in (365.L2) equals \(6\lambda-3\psi\).
The other three exponents exceed it by
\(15\lambda/2-6\), \(4-9\lambda/2\) and \(6-6\lambda\),
respectively, all positive at \(\lambda_0\).
Thus the slight increase of \(\psi_0\) allows a sufficiently
small fixed \(\eta>0\) with slack in every estimate. Monotonicity
then handles \(\lambda\ge\lambda_0\).
The bound for \(h\le4P^{1/3}\) has a power saving over \(W\)
since \(\lambda_0>5/6\), and the unit cutoff implies
\(P\ll W^\alpha\). Summing the dyadic boxes proves (365.L1),
using \(P\gg\sqrt W\) to write the saving as a power of \(W\).

By \(1/\tau_D\le \gcd(D,U_D)/(D\log\epsilon_D)\), (365.L1)
gives a reciprocal-rate tail \(O_\alpha(Z^{-\eta_\alpha})\)
for \(D>Z\) after adjusting \(\eta_\alpha>0\).
Truncating the floor formula (2) as in (5) gives the
power-saving error. The large-unit tail of the unrestricted
series in (3) remains open.

## Near-minimal units with a smaller period multiplier

**Theorem (A summable range with subpower period multipliers).**
Write \(H_D=\gcd(D,U_D)\), \(Q_D=U_D/H_D\), and recall
\(R_D=D/H_D\). Define
\[
 \mathcal L(x)=\log\log(e^e x),\qquad
 \mathcal F(x)=\log x\left(\frac{\log(2\mathcal L(x))}{\mathcal L(x)}
                    \right)^{1/3}.
\]
For each fixed \(0\le\theta<1/10\), there are \(C_\theta,c_\theta>0\)
such that the squarefree kernels in
\[
 \mathcal D_\theta
   =\{D>1:\ Q_D\le R_D^\theta,\ \log R_D\ge C_\theta\mathcal F(D)\}
                                                                  \tag{365.M1}
\]
satisfy
\[
      \sum_{\substack{D>Z\\D\in\mathcal D_\theta}}\frac1{\tau_D}
          \ll_\theta \exp(-c_\theta\mathcal F(Z))\qquad(Z\ge2).
                                                                  \tag{365.M2}
\]
Their consecutive powerful pairs have count
\(S_\theta Y+O_\theta(Y\exp(-c'_\theta\mathcal F(Y)))\),
where \(S_\theta=\sum_{D\in\mathcal D_\theta}1/\tau_D\),
\(Y=\operatorname{arcosh}(2X+1)\), and \(c'_\theta>0\).
Here \(Q_D\) measures the extra size of the unit beyond
\(U_D\ge H_D\). The lower cutoff for \(R_D\) in (365.M1) is
subpower in \(D\). At that cutoff, the permitted upper bound
for \(\epsilon_D\) is \(D^{3/2-o(1)}\).

Use the trace boxes \(T_D\asymp P,\ D\asymp W\) above, and
let \(h\) be the odd part of \(H_D\). The hypotheses give
\[
 P\asymp\sqrt W\,H_D Q_D,\qquad
 \sqrt W\ll P\ll W^{3/2}.
\]
If \(h\le4P^{1/3}\), then \(H_D\le8P^{1/3}\) and
\(P\ll W^{1/2+\theta}P^{(1-\theta)/3}\). This implies
\(P\ll W^{(3/2+3\theta)/(2+\theta)}\ll W^{6/7}\);
the earlier small-unit theorem already gives a power-saving
rate tail on this range.

For \(h>4P^{1/3}\) use the coprime shifted equation
\[
 e^3v^2g-d^3u^2f=s,\qquad s\in\{\pm1,\pm2\}.
\]
Again \(D=defg\) and \(de=h\). In addition \(fg\asymp R_D\)
and \(uv\asymp Q_D\), since \(H_D/h\in\{1,2\}\) and
\(uv=U_D/h\) or \(U_D/(2h)\). Subdivide into dyadic boxes
and let \(D_0,E_0,F_0,G_0,U_0,V_0\) be the respective sizes
of \(d,e,f,g,u,v\). Put
\[
 \psi=\frac{\log(D_0E_0)}{\log P},\qquad
 \beta=\frac{\log(F_0G_0)}{\log P},\qquad
 \sigma=\frac{\log(U_0V_0)}{\log P}.
\]
Thus \(\psi\le2/3+O(1/\log P)\), \(0<\beta\le1+O(1/\log P)\),
and \(\sigma\le\theta\beta+O(1/\log P)\), the second bound
using \(P^2\asymp (D_0E_0)^3F_0G_0(U_0V_0)^2\) and
\(D_0E_0\gg P^{1/3}\). The cutoff supplies
\(\beta\log P\ge C_\theta\mathcal F(W)-O(1)\).

### A determinant with growing degree

Choose \(\theta<\theta_1<\theta_2<1/10\) and then a fixed
\(b_0<1/3\) sufficiently close to \(1/3\) that
\(b_0^2/\theta_2>1/(1-\theta_2)\).
For sufficiently large \(P\) one has
\(\sigma\le\theta_1\beta\) and \(1-\psi\ge b_0\).
Set
\[
 \mu=(1-\theta_2)\beta,\qquad
 K=\text{nearest integer to }\frac12+\frac{b_0}{2\theta_2\beta},
\]
breaking ties arbitrarily. In this range \(K\ge2\) and
\(K\ll_\theta 1/\beta\).

Fix \(u,v\) and use the normalized ratios \(z,t\) from the
previous proof, so \(z=C^{-1/3}t^{1/3}+O(P^{-1})\), with
\(C\asymp1\). Use the \(2K\) monomials \(t^j,zt^j\)
for \(0\le j<K\) on intervals of length at most
\(\Delta=P^{-\mu}\) in \(t\). The clearing factor for
the \(2K\)-row determinant is
\(P^{K\psi+K(K-1)\beta}\exp(O(K^2))\):
restore the two ratio scales in the columns and multiply
each resulting row by \(ef^{K-1}\). Expand in \(j\) error
columns after centering \(t\) as in (365.L2). The decay
exponent of each term is at least
\[
 j+\mu\left(\frac{(2K-j)(2K-j-1)}2+\frac{j(j-1)}2\right).
\]
The constants outside these powers are \(\exp(O(K^2\log(2K)))\).
In detail, each smooth minor has the divided-difference
factor \(\Delta^{(2K-j)(2K-j-1)/2}\) and a determinant
of normalized derivatives of order at most \(2K-1\).
The functions \((t-t_0)^i\) and \(t^{1/3}(t-t_0)^i\)
are analytic on a fixed-radius complex neighborhood of
the real ratio range, with modulus at most \(\exp(O(K))\)
there. Cauchy's estimate bounds each derivative divided
by its factorial by \(\exp(O(K))\). Expanding the error
minor bounds it by \(j!\) times the product of the
\(P^{-1}\Delta^i\) factors for its distinct indices.
The column choices and Laplace expansion add at most
\(\exp(O(K\log(2K)))\) terms.

On putting \(\ell=K-j\), the excess of the decay over
the clearing exponent is bounded below by
\[
        K b_0-\theta_2\beta K(K-1)-\frac1{4(1-\theta_2)\beta}
 \ \ge\  \frac{b_0}2+
          \frac1{4\beta}
           \left(\frac{b_0^2}{\theta_2}-\frac1{1-\theta_2}\right).
\]
The first inequality on the excess uses
\(\mu\ell^2-\ell\ge-1/(4\mu)\).
The displayed second inequality follows by rounding
the vertex of the quadratic in \(K\) to the nearest integer.
The surplus dominates the prefactor, since the lower bound
on \(\beta\) makes \(K\ll_\theta(\log\log P)^{1/3}\).
Each cleared determinant is therefore zero for sufficiently
large \(P\). Every interval admits an actual nonzero
auxiliary form of bidegree at most \((1,K-1)\) in the
two ratios. Before taking cofactors the integral entries
are \(O(P^{O(1)})\), since \(d,e\ll P^{1/3}\),
\(f,g\le 2P^\beta\) and \(K\beta=O_\theta(1)\).
The coefficient logarithmic height after taking cofactors
is consequently \(O_\theta(K\log P)\).

### Dependence of the point count on the degree

Repeat the degree-one projection argument leading to (365.L3).
The fixed-ratio cases remain \(O(1)\) per ratio and there are
\(O(K)\) exceptional ratios. For the other ratios, the
homogeneous forms \(A_1,A_2\) have degree \(1\le k<K\)
and logarithmic height \(O_\theta(K\log P)\). For this height
bound one can take primitive factors: each coefficient of a
degree-\(d\) factor is bounded by \(2^d\) times its Mahler measure,
and the measure of an integral factor is at most the measure of
the polynomial being factored. The latter is bounded by its
coefficient sum. These elementary estimates follow from the
root product and Jensen's formula. The two chart resultants
therefore have logarithmic size \(O_\theta(K^2\log P)\), by
the Sylvester determinant, and \(q\) divides their product
as before. In (365.L3) the degree is \(3k+1\le3K\) and
the right-hand side has logarithmic size \(O_\theta(K^2\log P)\).
The same ramification argument excludes a power of a single
linear form for every \(k\ge1\).

For \(n\ge1\), both \(\log\tau(n)\) and \(\omega(n)\) are
\(O(\log(3n)/\log\log(3n))\). To see this for \(\log\tau(n)\),
split the primes at \(\sqrt{\log(3n)}\). Below the cutoff
there are at most \(\sqrt{\log(3n)}\) factors and each contributes
\(O(\log\log(3n))\). Above the cutoff the sum of exponents is
at most \(2\log n/\log\log(3n)\), and each summand
\(\log(a+1)\) is at most \(a\log 2\).
Stewart's bound for the irreducible-factor case of degree at
least three is \(2800r^{1+\omega(|m|)}\) at degree \(r\) and
nonzero right side \(m\), after dividing out content.
The reducible-factor case counts pairs of nonzero divisor
values with at most \(O(K^2)\) intersections. The quadratic
case counts ideals of norm with logarithmic size
\(O_\theta(K^2\log P)\) and at most \(O_\theta(K^2\log P)\)
unit generators per ideal in the given height range.
Consequently, including the choice of \(q\), the total solutions
on one interval for the fixed \(u,v\) are at most
\[
        \exp\left(O_\theta\left(
              \frac{K^2\log(2K)\log P}{\log\log P}\right)\right).
                                                                  \tag{365.M3}
\]
All the degree dependence here is explicit enough to let
\(K\) vary with \(P\).
[Stewart, consequence of Theorem 1 on p. 796](https://uwaterloo.ca/pure-mathematics/sites/default/files/uploads/documents/s0894-0347-1991-1119199-x_0.pdf)

Since \(K\ll_\theta1/\beta\) and \(P\) lies between bounded
powers of \(W\), choosing \(C_\theta\) large in (365.M1)
makes the exponent in (365.M3) at most
\((\theta_2-\theta_1)\beta\log P/2\). Indeed,
\(K^2\log(2K)/(\beta\log\log P)\) is then
\(O_\theta(C_\theta^{-3})\).
There are \(O(P^\mu)\) intervals and \(O(P^\sigma)\) pairs
\(u,v\) in the box. The number of traces in it is at most
\[
                O_\theta\!\left(P^{(1-\kappa_\theta)\beta}\right),
       \qquad \kappa_\theta=\frac{\theta_2-\theta_1}{2}>0 .
\]
Each kernel has reciprocal period at most
\(1/(R_D\log\epsilon_D)\ll P^{-\beta}/\log P\).
Thus the rate in the box has the saving
\(\exp(-\kappa_\theta\beta\log P)\). The dyadic box count
is polynomial in \(\log W\), so it is absorbed into this
saving. Summing over \(W\ge Z\) preserves the shape of the
exponential bound with a smaller constant. For instance,
on the dyadic indices \(j=\lfloor\log_2 W\rfloor\) the
decay exponent dominates \(c\sqrt{j}\) for large \(j\).
One can retain half of the exponent at \(Z\) using
eventual monotonicity of \(\mathcal F\) and sum the other
half over the indices. This proves (365.M2). In (2) keep
the floors for \(D\le\sqrt Y\), for an error \(O(\sqrt Y)\).
The rest of the error is at most \(Y\) times the rate tail
at \(\sqrt Y\), which gives the claimed logarithmic law.

## A larger square factor range for polynomial period multipliers

**Theorem (Weighted estimate with larger square factors).**
For squarefree \(D>1\) retain \(H_D=\gcd(D,U_D)\),
\(R_D=D/H_D\), and \(Q_D=U_D/H_D\). For \(W\ge2\),
\[
 \sum_{\substack{W\le D<2W,\ D\ {\rm squarefree}\\
                  R_D\ge D^{5/16},\ Q_D\le R_D^{1/5}}}
       \frac{H_D}{D\log\epsilon_D}\ll W^{-1/5000}.                 \tag{365.P1}
\]
The hypotheses permit \(\epsilon_D\) up to order \(D^{5/4}\), since
\(\epsilon_D\asymp \sqrt D\,H_D Q_D\le D^{3/2}R_D^{-4/5}\).
The rate in this range is finite, and (2) gives a restricted
logarithmic asymptotic with a power saving in \(Y\).

Use the trace boxes \(T_D\asymp P\), \(D\asymp W\) and let \(h\)
be the odd part of \(H_D\). If \(h\le4P^{1/3}\), then
\[
 P\asymp \sqrt W\,H_D Q_D
 \ll W^{7/10}H_D^{4/5}\ll W^{7/10}P^{4/15}.
\]
Thus \(P\ll W^{21/22}\). For sufficiently large \(W\) these
kernels have \(\epsilon_D\le D\), and the estimate from (6) at
\(\alpha=1\) gives (365.P1) on this part.

For \(h>4P^{1/3}\) use the factorization in the proof of
(365.M2) and the same dyadic scales. In an occupied box,
\[
 \sigma\le\frac{\beta}{5}+O(1/\log P),\qquad
 2=3\psi+\beta+2\sigma+O(1/\log P).
\]
As \(\log D=(\psi+\beta)\log P+O(1)\) and
\(\log R_D=\beta\log P+O(1)\), the lower bound on \(R_D\)
gives \(38\beta+10\sigma\ge10-O(1/\log P)\).
Consequently \(1/4-O(1/\log P)\le\beta\le1+O(1/\log P)\).
The upper bound uses \(\psi\ge1/3-O(1/\log P)\) and \(\sigma\ge0\).

For fixed \(u,v\) subdivide the normalized \(t\)-range into intervals
of length at most \(P^{-\mu}\), with
\(\mu=\beta-\sigma-1/1000\). Use the determinant with \(2K\)
columns from the preceding proof, taking \(K=5\) for
\(\beta\le13/50\) and \(K=4\) for \(\beta>13/50\). Put
\[
 A_j=\frac{(2K-j)(2K-j-1)+j(j-1)}2\qquad(0\le j\le K).
\]
The excess of the decay exponent for the term with \(j\) error
columns over the clearing exponent is
\[
 j+\mu A_j-K\psi-K(K-1)\beta
 \ \ge\
 j-\frac{2K}{3}
 +\beta\!\left(\frac45 A_j+\frac7{15}K-K(K-1)\right)
 -\frac{A_j}{1000}-O(1/\log P).
\]
Indeed the coefficient of \(\sigma\) before bounding it is
\(2K/3-A_j<0\). The minimum of the displayed affine expression
before the last two subtractions can be checked at the endpoints
for each \(j\).

| \(K\) | Interval for \(\beta\), up to \(O(1/\log P)\) | Minimum before subtractions | Largest \(A_j\) |
| --- | --- | --- | --- |
| \(5\) | \([1/4,13/50]\) | \(1/20\) | \(45\) |
| \(4\) | \([13/50,1]\) | \(2/75\) | \(28\) |

For \(K=4\) the minimum occurs at \(j=2\), where \(A_2=16\).
For \(j\ne2\) in that row the minimum is at least \(1/15\).
Thus all the excesses are positive by a fixed margin for large \(P\).
Every interval has a nonzero auxiliary equation of bidegree at most
\((1,K-1)\), of polynomial coefficient height in \(P\).

The count on an auxiliary equation in the preceding proof is
\(O_\varepsilon(P^\varepsilon)\) for these two fixed degrees.
For the irreducible-form case this uses the primitive Thue bound of
[Stewart, p. 796](https://uwaterloo.ca/pure-mathematics/sites/default/files/uploads/documents/s0894-0347-1991-1119199-x_0.pdf),
with both the degree and the polynomial size of the right side as
in (365.M3). The \(O(P^\sigma)\) choices of \(u,v\) and
\(O(P^\mu)\) intervals give a count
\(O_\varepsilon(P^{\beta-1/1000+\varepsilon})\) per box.
The weight in (365.P1) is \(1/(R_D\log\epsilon_D)\),
which is \(\ll P^{-\beta}/\log P\). There are
only polynomially many boxes in \(\log W\) since
\(\sqrt W\ll P\ll W^{5/4}\) and all factor sizes are bounded
by a fixed power of \(W\). Using a sufficiently small fixed
\(\varepsilon>0\) proves (365.P1).
Summing the dyadic bounds gives a rate tail \(O(Z^{-1/5000})\).
Keeping the floors in (2) for \(D\le\sqrt Y\) gives the claimed
power saving in \(Y\).

## Very small period multipliers and logarithmic square factors

**Theorem (Weighted tail for very small period multipliers).**
Fix \(a,b\ge0\) with \(2a+b<1\). Then
\[
 \sum_{\substack{D>Z,\ D\ {\rm squarefree}\\
              \log R_D\le(\log D)^a,\ Q_D\le(\log D)^b}}
       \frac{H_D}{D\log\epsilon_D}
 \ll_{a,b} (\log Z)^{2a+b-1}\log(2\log Z)\qquad(Z\ge3).            \tag{365.S1}
\]
Here the units satisfy \(\epsilon_D=D^{3/2+o(1)}\) along the
specified range as \(D\to\infty\), since
\(\epsilon_D\asymp D^{3/2}Q_D/R_D\).

For each \(D\) assign the odd primes in \(H_D\) to the factor
of \(T_D^2-1\) in which they occur. As in the proof of (365.M2),
this gives a tuple of positive integers with
\[
 e^3v^2g-d^3u^2f=s,\qquad
 s\in\{\pm1,\pm2\},\qquad D=defg,\qquad
 fg\asymp R_D,\quad uv\asymp Q_D
\]
and coprime shifted terms. This factorization uses no lower bound
on \(de\). Fix \(f,g,u,v,s\). The form \(v^2g X^3-u^2f Y^3\)
has nonzero discriminant and content one in an occupied tuple,
and \(\gcd(d,e)=1\). The primitive-solution bound of Stewart
on p. 796 gives \(O(1)\) pairs \((e,d)\), uniformly in the
coefficients, since the degree is three and the right side has
absolute value at most two. Each pair fixes \(D\).

On \(\exp(L)\le D<\exp(2L)\), with \(L\ge1\), the weight is
\(\ll1/(fg L)\) since \(\log\epsilon_D>\tfrac12\log D\).
The hypotheses give \(fg\ll\exp((2L)^a)\) and \(uv\ll (2L)^b\).
Summing after the uniform primitive count bounds the block by
\[
 \frac1L
 \sum_{fg\ll \exp((2L)^a)}\frac1{fg}
 \sum_{uv\ll(2L)^b}1
 \ \ll_{a,b}\ L^{2a+b-1}\log(2L).
\]
The first sum is bounded by the product of two harmonic sums,
and the second by summing the allowed lengths in \(v\) over
\(u\). Now double \(L\) starting at \(\log Z\); the block bound
is summable because \(2a+b<1\). This proves (365.S1).

## A joint count in the period ratio and unit quotient

The same primitive Thue estimate also gives a count without an upper
bound on the fundamental unit or on the common divisor. Write
\(d(n)\) for the divisor function in this section.

**Proposition (Uniform joint fibers).**
For all positive integers \(R,Q\),
\[
 \#\{D>1:D\text{ squarefree},\ R_D=R,\ Q_D=Q\}
                    \ll d(R)d(Q).                              \tag{365.J1}
\]
The constant is absolute. In fact the bound also holds when all positive
integral Pell solutions are counted, without requiring fundamentality.

Here is the normalization, including the prime two. Put
\(H=\gcd(D,U)\), let \(h\) be its odd part, and take the coprime factors
\(A=T-1,B=T+1\) if \(T\) is even and
\(A=(T-1)/2,B=(T+1)/2\) if \(T\) is odd. For each odd prime in \(h\),
its cube divides exactly one factor. Assigning these primes gives
\[
 A=d^3fu^2,\quad B=e^3gv^2,\quad de=h,\quad D=defg,
 \qquad g v^2e^3-f u^2d^3=s.                         \tag{365.J2}
\]
All variables are positive, \(\gcd(d,e)=1\), and the coefficients
\(fu^2,gv^2\) are coprime. The exact possibilities are

| Parity | \(fg\) | \(uv\) | \(s\) |
|---|---|---|---|
| \(T\) even | \(R\) | \(Q\) | 2 |
| \(T\) odd, \(D\) odd | \(R\) | \(Q/2\) | 1 |
| \(T\) odd, \(D\) even | \(2R\) | \(Q\) | 1 |

The second row requires even \(Q\). For example, in the last row
\(H=2h\) and \(AB=D(U/2)^2=h^3(2R)Q^2\), which explains both
factors of two. The Pell congruences exclude \(T\) even with \(D\) even.

For each row there are \(O(d(R)d(Q))\) ordered choices of \(f,g,u,v\),
since \(d(2n)\le2d(n)\). Fix those choices. The binary cubic
\(gv^2X^3-fu^2Y^3\) has content one and nonzero discriminant.
[Stewart, consequence of Theorem 1, printed p. 796](https://uwaterloo.ca/pure-mathematics/sites/default/files/uploads/documents/s0894-0347-1991-1119199-x_0.pdf)
bounds its primitive solutions at either fixed right side \(s=1,2\)
by an absolute constant, including reducible cubics. Each solution
\((e,d)\) recovers \(D\), \(T\), and then positive \(U\) uniquely.
This proves (365.J1). This is an application of that classical uniform
Thue bound, not a new uniform Thue theorem.

### The critical logarithmic boundary

**Corollary (A logarithmic saving on the multiplier boundary).**
Fix \(0\le a<1/2\) and \(\eta>0\). Put \(b=1-2a\) and
\(\ell(D)=\log\log(e^eD)\). Then
\[
 \sum_{\substack{D>Z,\ D\text{ squarefree}\\
       \log R_D\le(\log D)^a,\quad
       Q_D\le(\log D)^b/\ell(D)^{2+\eta}}}
       \frac{H_D}{D\log\epsilon_D}
                \ll_{a,\eta}\ell(Z)^{-\eta}\qquad(Z\ge3).     \tag{365.J3}
\]
Thus the strict inequality \(2a+b<1\) in (365.S1) can be replaced
by equality with the displayed extra saving. This remains a restricted
sector; it supplies no bound for unrestricted large \(Q_D\).

To prove it, use the blocks \(e^L\le D<e^{2L}\). The weight is
\(\ll1/(R_DL)\). In (365.J2) the coefficient products satisfy
\(fg\ll\exp((2L)^a)\) and
\(uv\ll L^b/(\log L)^{2+\eta}\), for sufficiently large \(L\).
The uniform primitive count just proved bounds the block by
\[
 \frac1L\sum_{fg\ll\exp((2L)^a)}\frac1{fg}
       \sum_{uv\ll L^b/(\log L)^{2+\eta}}1
 \ll_{a,\eta}\frac{L^{2a+b-1}}{(\log L)^{1+\eta}}
 =\frac1{(\log L)^{1+\eta}}.
\]
Here \(\sum_{uv\le V}1\ll V\log(2V)\), and \(b>0\).
Summing over \(L=2^j\log Z\) gives (365.J3); bounded initial
blocks are absorbed in the constant. The exact count (2) and dominated
convergence give a logarithmic asymptotic on this sector.

### Polynomial period counts with a period-dependent quotient cutoff

**Corollary (Polynomial quotient sector).**
Uniformly for \(Y\ge2,V\ge1\),
\[
 \#\{D>1:D\text{ squarefree},\ \tau_D\le Y,\ Q_D\le V\}
                         \ll YV\log(2V).                       \tag{365.J4}
\]
Indeed \(\tau_D\ge R_D\log\epsilon_D>\tfrac12R_D\log D\)
and \(D\ge R_D\). Consequently
\(R_D\ll Y/\log(2Y)\); the finitely small values of \(R_D\)
cause no exception. Sum (365.J1), using
\(\sum_{n\le B}d(n)\ll B\log(2B)\).

For any fixed \(B>0\), the fixed sector
\(\{D:Q_D\le\tau_D^B\}\) therefore has period count
\(O_B(Y^{1+B}\log(2Y))\), and its powerful-pair count is
\(O_B((\log X)^{1+B}\log\log X)\), by (2) and
\(\sum_{k\ge1}k^{-1-B}<\infty\). To finish the weaker global
polynomial-period target by this route, one still has to bound the
complementary kernels with \(Q_D>\tau_D^B\) for some fixed \(B\).

The full rate in (3) remains open. A sufficient estimate on the
complement of the proved summable ranges is a bound
\(O((\log W)^{-1-\eta})\) for some \(\eta>0\) on its weighted sum
\(\sum_{W\le D<2W} H_D/(D\log\epsilon_D)\). In particular the
estimates above still leave large \(Q_D\) together with large
\(H_D\) outside their unit and multiplier cutoffs.

The [joint endpoint control script](verify_arithmetic_endpoint_joint_controls.py)
was run internally with Python 3.12.3 on 27 November 2026 from the
repository root. Its exact prime-two normalization checks passed for
all 3041 squarefree \(2\le D\le5000\): 770 even-trace, 1256
odd-trace/odd-kernel, and 1015 odd-trace/even-kernel cases. It also
checks the finite witness-cover controls specified in the
[squarefree return note](ERDOS_675_SQUAREFREE_RETURNS.md#joint-square-witnesses-and-their-least-crt-representatives).
These controls check vulnerable normalizations and examples; the
unrestricted statements (365.J1)--(365.J4) use the written proofs above.
