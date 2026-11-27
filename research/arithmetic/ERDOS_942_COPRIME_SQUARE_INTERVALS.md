# Erdős #942 and coprime occupants of square intervals

Write \(h(q)=\#\{t:q^2\le t<(q+1)^2,\ t\text{ powerful}\}\).
The left endpoint \(q^2\) is powerful, so the count in the
open interval is \(h(q)-1\).
[Erdős #942](https://www.erdosproblems.com/942), from
[Erdős, printed p. 35](https://www.renyi.hu/~p_erdos/1976-39.pdf#page=11),
asks for matching logarithmic power upper and infinitely often lower
bounds. The logarithmic upper bound remains open. The uniform power
upper bound with exponent \(48/211\) is proved below.

The [short power window theorem CPW1](ABC_COPRIME_SHORT_POWER_WINDOWS.md#a-positive-density-primitive-construction)
gives positive density of intervals with any fixed number of pairwise
coprime nonsquare powerful occupants. Its constants depend on that number.
The [coprime occupancy density](#limiting-density-of-large-coprime-occupancy)
below quantifies how rapidly the proportion falls with multiplicity.
[Narumi–Tachiya, Theorems 1–2 and Lemmas 1–2](https://arxiv.org/html/2611.20656v1)
use simultaneous approximation to obtain ordinary occupancy of order
at least \(\log q/\log\log q\) infinitely often. Their approximation
conditions do not impose pairwise coprimality.

## Uniform upper bound from local conics

**Theorem (Uniform square-interval upper bound).**
For every \(\varepsilon>0\),
\[
                          h(q)\ll_\varepsilon q^{109/474+\varepsilon}. \tag{8}
\]
The bound also holds for a window of bounded length in the square-root
coordinate. Thus an interval \([x,x+C\sqrt{x}]\) with fixed \(C>0\)
contains \(O_{C,\varepsilon}(x^{109/948+\varepsilon})\) powerful integers.
The earlier \(6/25\) square-root exponent used Filaseta–Trifonov in
both ranges. The proof below retains that estimate in the outer
blocks, uses local conics in the balanced blocks, and counts a
second integer quantity attached to pairs of points between these
ranges. Rational approximation to a local slope gives a thinner
integer projection where that condition has a wider error.
[Hughes, §§8–9](https://raw.githubusercontent.com/scottdhughes/erdos942/37feb687a3b2a75e5cae00297b73532da4b088e9/paper/main.pdf)

Write an occupant of the open square interval uniquely as \(a^2b^3\)
with \(b\ge2\) squarefree. At least one of \(a,b\) is at most
\((q+1)^{2/5}\). Use that variable for the abscissa \(u\) and sum
dyadically over \(N<u\le2N\) with \(N\ll q^{2/5}\). Ignoring
squarefreeness only enlarges the count. In the \(b\)-aspect use
\(f(u)=q u^{-3/2}\), and in the \(a\)-aspect use
\(f(u)=q^{2/3}u^{-2/3}\). Write \(f(u)=X u^{-s}\) and \(T=XN^{-s}\).
In both cases the vertical error is \(\Delta\ll T/q\). For example
in the \(a\)-aspect, the mean value theorem gives
\(((q+1)/u)^{2/3}-(q/u)^{2/3}\ll q^{-1/3}N^{-2/3}\).
In a window \([x,x+C\sqrt{x}]\) take \(q=\sqrt{x}\), so the upper
endpoint has square root \(q+O_C(1)\). The same two bands suffice with
\(C\)-dependent constants in \(\ll\). Bounded \(N\) contributes
\(O_C(1)\). For the other blocks there is at most one ordinate
per integer abscissa in the band, since \(\Delta<1/2\).

### Unbalanced blocks

Use Theorem 7 of
[Filaseta–Trifonov](sources/filaseta-trifonov-fractional-parts-1996.pdf)
with order parameter \(k=2\), \(s=3/2\) or \(2/3\). Its hypotheses
\(N^s\le X\) and \(\Delta\le c_s N^{-1}\) hold in these ranges for
sufficiently large \(N\) or \(q\), respectively. The two terms are
\[
             X^{1/5}N^{(2-s)/5}+\Delta X^{1/15}N^{(27-s)/15}.
\]
The second is \(O_C(q^{11/75})\) in both aspects using
\(N\ll q^{2/5}\). Put
\[
 p=\frac{109}{474},\qquad
 \beta_0=10p-2=\frac{71}{237},\quad
 \beta_*=1-3p=\frac{49}{158},\quad
 \beta_1=\frac32-5p=\frac{83}{237}.
\]
For \(N\le q^{\beta_0}\) in the \(b\)-aspect the first term is
\(q^{1/5}N^{1/10}\le q^p\). For \(N\le q^{1/3}\) in the
\(a\)-aspect it is \(q^{2/15}N^{4/15}\le q^{2/9}\).

### Intermediate cube-kernel blocks

Let \(\beta=\log N/\log q\), and suppose
\(\beta_0\le\beta\le\beta_1\) in the \(b\)-aspect.
By Lemma 9 of Filaseta–Trifonov, take their subset \(S'\) of the
abscissae with original cardinality \(O(\#S'+1)\). For an integer
\(a\ge1\), use \(t(a)\) for the number of consecutive pairs in \(S'\)
whose gap is \(a\), as in their proof. Put \(A=Nq^{-p}\), so gaps
\(a>A\) number \(O(q^p)\). For small gaps use the cutoff
\[
 A_0=\begin{cases}
 q^{(p-1)/4}N^{7/8}, & \beta\le\beta_*,\\
 C_0 N^{1/4}, & \beta\ge\beta_*,
 \end{cases}
\]
using either one at the common endpoint and a large enough fixed
\(C_0\) in the second case.
The individual pair estimate in the proof of Filaseta–Trifonov's
Theorem 7, using Lemma 10, gives before their final summation
\[
 \sum_{a\le A_0}t(a)
 \ll 1+\Delta q^{1/3}N^{1/2} A_0^{4/3}+qN^{-7/2}A_0^4.
\]
Its inputs are the separation in Lemma 9, \(q>N^2\),
\(A_0\le N^{2/3}\), and \(\Delta\ll_C N^{-3/2}\).
For the first cutoff the powers of \(q\) in the last two terms are
\(p/3+\beta/6\) and \(p\). For the second they are
\(1/3-2\beta/3\) and \(1-5\beta/2\), both less than \(p\) for
\(\beta\ge\beta_*\). Hence this count is \(O_C(q^p)\).

For \(A_0<a\le A\), write the two points as \((u,v),(u+a,w)\).
Use the integer \(k=(4u+5a)w-(4u-a)v\), whose on-curve value is
\(\Phi_a(u)=(4u+5a)f(u+a)-(4u-a)f(u)\). This is the first-order
Padé expression of Filaseta–Trifonov, with sign and scale chosen
for \(s=3/2\). Inverting a pair invariant to recover the abscissa
also occurs in [Pandey, §§2.4–2.6](https://arxiv.org/html/2401.13981v3)
for the reciprocal-square phase.
Put
\[
 K=q a^3N^{-7/2},\qquad V=q aN^{-5/2},\qquad E=N^{-1/2},
 \qquad L=NE/K,\qquad B_a(u)=f(u)-f(u+a).
\]
Then \(k=\Phi_a(u)+O_C(E)\). On \(u=Nx\), \(1/2\le x\le3\),
\[
 \frac{\Phi_a(Nx)}{K}=\frac58 x^{-7/2}+O(a/N),\qquad
 \frac{B_a(Nx)}{V}=\frac32 x^{-5/2}+O(a/N),
\]
with the same errors in derivatives through order four in \(x\).
The first expansion follows by canceling the terms of orders zero,
one and two in \(a/u\). In both gap ranges \(K\gg_C q^{1/30}\)
and \(a/N=o(1)\). Thus \(k\asymp K\) is positive. On the expanded
interval the inverse \(U_a(k)=\Phi_a^{-1}(k)\) exists for all pairs,
and \(|u-U_a(k)|\ll_C L\).
Let \(g_a(k)=B_a(U_a(k))\). In the variable \(\xi=k/K\),
\[
 \frac{U_a(k)}N=c_1\xi^{-2/7}+O(a/N),\qquad
 \frac{g_a(k)}V=c_2\xi^{5/7}+O(a/N)
                                                                  \tag{10}
\]
for fixed \(c_1,c_2>0\), with the same errors in derivatives through
order three. This follows by differentiating the inverse of the
first expansion, whose derivative is bounded away from zero
on the compact \(x\)-range. In particular, \(g_a'''\asymp V/K^3
=:\Lambda\). In addition to \(k\in\mathbb Z\), a pair requires
the two integers \(u\) and \(b=v-w\), the latter satisfying
\(b=B_a(u)+O_C(\Delta)\).

#### Counting the ordinate difference for \(\beta_*\le\beta\le\beta_1\)

As \(|B_a'|\ll V/N\), we have
\[
              |b-g_a(k)|\ll_C VE/K+\Delta\ll_C N^{1/2}/a^2=O_C(\omega),
              \qquad \omega=N^{1/2}/a^2.
\]
The tolerance is less than \(1/4\) by the choice of \(C_0\), so \(b\)
is unique for each \(a,k\).
[Robert–Sargos, Corollary 2](sources/robert-sargos-third-derivative-2003.pdf)
bounds the number of possible \(k\) by
\[
 \ll_\varepsilon K\omega + K^{1+\varepsilon}\Lambda^{3/19}
                         +K^{3/4+\varepsilon}+\Lambda^{-1/3}.       \tag{11}
\]
The derivative hypothesis holds by (10) on an interval of length
comparable to \(K\). If \(\Lambda\) is bounded below, the trivial
\(O(K)\) count suffices instead and is covered by (11).

For a given \(a,k,b\) the congruence \(4ub\equiv-k\bmod a\)
follows from \(k=6aw-(4u-a)b\). If soluble, it gives a class
with spacing \(a/\gcd(a,4b)\), and \(\gcd(a,4b)\mid k\). The
number of \(u\) is therefore \(O_C(1+L\gcd(a,k)/a)\). The
additional multiplicities cost \(O_{C,\varepsilon}(N^{1/2}q^\varepsilon)\)
over all \(a\), since \(\sum_{1\le k\ll K}\gcd(a,k)\ll K\tau(a)\)
by summing \(d\lfloor O(K)/d\rfloor\) over \(d\mid a\).

The four powers of \(q\) obtained by summing (11), evaluating at
\(a=A=q^{\beta-p}\) and bounding the number of gaps by \(A\), are
\[
 1-\beta-2p,\qquad
 \frac{13}{19}+\frac{\beta}{2}-\frac{52p}{19},\qquad
 \frac34+\frac{5\beta}{8}-\frac{13p}{4},\qquad
 \frac23+\beta-\frac{11p}{3}.
\]
Each term increases as a function of \(a\) before evaluation.
The four maxima on the present range are
\(p,\ p,\ 35/158,\ 247/1422\), respectively.

#### Rational-slope projection for \(\beta_0\le\beta\le\beta_*\)

For each \(a\) in this range let \(Q=\lfloor q^{1/64}\rfloor\)
and subdivide the inverse interval in \(k\) into
\(J\ll 1+(V/N)Q^2\) pieces. On each piece the slope
\(B_a'(U_a(k))\) varies by \(O(Q^{-2})\), since its derivative
is \(O(V/(NK))\). Dirichlet approximation at one point of the
piece supplies integers \(1\le\ell\le Q\) and \(r\) with
\(\left|\ell B_a'(U_a(k))-r\right|\ll Q^{-1}\) throughout
the piece. This is the local rational-slope choice applied to
the two recovered coordinates. [Huxley–Ivić, §6](https://emis.de/ft/41261)

The integer \(\ell b-r u\) must now approximate
\[
         H_a(k)=\ell g_a(k)-r U_a(k)
\]
within \(O_C(Q\Delta+L/Q+QVL^2/N^2)=:\omega'\). This follows
by Taylor expansion of \(B_a(u)\) about \(U_a(k)\). To check the
new third derivative, first subtract
\(B_a'(U_a(k))U_a'''(k)\) from \(g_a'''(k)\).
The leading term from (10) is
\[
            -\frac{630}{343} c_2\,\frac{V}{K^3}\xi^{-16/7}.
\]
Indeed \(B_a'(U_a(k))=g_a'(k)/U_a'(k)\), with leading term
\(-\frac52(c_2/c_1)(V/N)\xi\).
The derivative error from the rational slope, after dividing
by \(\ell\), is \(O(N/(\ell QK^3))=o(V/K^3)\), using the
bounds below. Hence \(-H_a'''\asymp\ell\Lambda\) on the piece.

For clarity, the following rational bounds suffice. Writing
\(\alpha=\log a/\log q\), we have
\[
 \frac{p-1+7\beta/2}{4}\le\alpha\le\beta-p.
\]
At the four corners of this trapezoid, and hence throughout it,
\[
 a\le q^{81/1000},\quad K\le q^{161/1000},\quad
 q^{-7/1000}\le V/N\le q^{22/1000},\quad L\le q^{4/1000},\quad
 q^{-161/1000}\le\Lambda\le q^{-149/1000}.
\]
These use the exponents \(\alpha,\ 1-7\beta/2+3\alpha,\
1-7\beta/2+\alpha,\ 4\beta-1-3\alpha,\ -2+8\beta-8\alpha\)
in the displayed order. In particular \(J\ll q^{54/1000}\),
\(J\ll K\), and \(\omega'\ll_C q^{4/1000-1/64}\), since
\(\Delta\ll_C N^{-3/2}\) and \(V/N^2\ll q^{-27/100}\).
Also \(\sum_{a\le A}K\ll q^{236/1000}\), by the maximum
\(1+\beta_*/2-4p=223/948<236/1000\) for its exponent.

Apply the same Robert–Sargos count to \(-H_a\) in each piece
of length \(M_i\), rounding endpoints and counting pieces with
\(M_i<1\) directly. The first two error and main terms use
\(\sum M_i\ll K\), and
\(\sum M_i^{3/4}\le J^{1/4}(\sum M_i)^{3/4}\).
For every candidate \(a,k\) there are \(O_C(1+L)\) choices of
the integer \(u\). Summing over the pieces and then over \(a\)
gives, up to \(q^\varepsilon\), the following exponent bounds:
\[
 \begin{array}{c|l}
 \text{term}&\text{exponent bound}\\ \hline
 M_i\omega' & 236/1000+8/1000-1/64\\
 M_i(\ell\Lambda)^{3/19}
       &236/1000+4/1000+(3/19)(1/64-149/1000)\\
 M_i^{3/4}&81/1000+ (3/4)(161/1000)+(1/4)(54/1000)+4/1000\\
 (\ell\Lambda)^{-1/3}&81/1000+161/3000+54/1000+4/1000 .
 \end{array}
\]
They are all strictly below \(p\). The count of rounded or
short pieces is bounded by the last row. Combining these two gap
ranges with the discarded small and large gaps proves
\(O_{C,\varepsilon}(q^{p+\varepsilon})\) throughout
\(\beta_0\le\beta\le\beta_1\).

### Balanced blocks and local conics

Use \(q^{\beta_1}\le N\ll q^{2/5}\) in the \(b\)-aspect and
\(q^{1/3}\le N\ll q^{2/5}\) in the \(a\)-aspect. Put
\(\sigma=p\) and \(\sigma=11/50\) respectively. Partition the block into \(O(q^{\sigma+\eta})\)
intervals of length \(H\le Nq^{-\sigma-\eta}\) for any small fixed
\(\eta>0\). The normalized points \((x,y)=(u/N,v/T)\) satisfy
\[
                     y=x^{-s}+e,\qquad |e|\ll_C q^{-1}.
\]
For six points in one interval take the integer determinant of
the columns \(1,u,v,u^2,uv,v^2\). After removing the factor \(N^4T^4\),
the normalized row is a sum of three rows with respective factors
\(1,e,e^2\). The middle row has nonzero entries only in the
\(y,xy,y^2\) columns, and the last only in the \(y^2\) column.
In a nonzero term let \(m_0,m_1,m_2\) denote the counts of these
three kinds of row. Thus \(m_0+m_1+m_2=6,\ m_1\le3,\ m_2\le1\).
Within each kind, remove the scalar \(e^j\) and take divided
differences among the rows to extract their Vandermonde product.
The remaining entries are bounded by fixed derivative bounds.
The term before restoring \(N^4T^4\) is therefore
\[
   \ll_C q^{-m_1-2m_2-(\sigma+\eta)
                           \left(\binom{m_0}{2}+\binom{m_1}{2}\right)} .
\]
For either \(\sigma\) the smallest decay exponent at \(\eta=0\)
is \(1+10\sigma\), attained at \((5,1,0)\). In the \(b\)-aspect
\(N^4T^4=q^4/N^2\le q^{4-2\beta_1}=q^{1+10\sigma}\). In the \(a\)-aspect
it is \(q^{8/3}N^{4/3}\ll_C q^{16/5}=q^{1+10\sigma}\).
The positive \(\eta\) makes the restored determinant \(o(1)\).
Integrality makes it zero for large \(q\). Hence all the band
points in an interval lie on a single polynomial of degree at most
two. A relation from minors of a maximal independent set has
integer coefficients of size \(q^{O(1)}\), and we divide out their
common gcd. Blocks with fewer than six points can be counted directly.

Every line has at most two band points. To see this, take
\(u_1<u_2<u_3\). The determinant in \(1,u,v\) differs from the
on-curve determinant by \(O(\Delta(u_3-u_1))\). The on-curve
magnitude is comparable to
\[
             (T/N^2)(u_2-u_1)(u_3-u_2)(u_3-u_1)
\]
by the second divided difference of \(f\).
Since \(\Delta/(T/N^2)\ll N^2/q=o(1)\) and the gaps are integral,
the band determinant is nonzero. In particular the span of every
triple is \(\gg R\), where \(R=(N^2/T)^{1/3}\). These balanced
ranges have \(R\to\infty\), \(\Delta N=o(1)\), \(T\gg_C N\), and
\(H\ll\sqrt N\).

A quadratic of nonzero discriminant contributes
\(O_\varepsilon(q^\varepsilon)\) band points. For completeness
write it as \(Au^2+Buv+Cv^2+Du+Ev+F\) with
\(\mathfrak d=B^2-4AC\ne0\). A bounded integral change of variables
achieves \(A\ne0\). For \(J=2BD-4AE\), \(K=D^2-4AF\), completing
squares gives
\[
 (2\mathfrak d v+J)^2-\mathfrak d(4Au+2Bv+2D)^2=J^2-4\mathfrak d K .
\]
The substitution is injective and has polynomially bounded height.
A nonzero right-hand side gives either integer factorization or a
quadratic-field norm equation whose principal ideals are counted by
the divisor bound. Imaginary quadratic units are finite of order at
most six. In a real quadratic field, the two embeddings of a
generator have absolute values \(\ll q^{O(1)}\), and the nonzero
norm gives the lower absolute bound \(q^{-O(1)}\) in either
embedding. The positive fundamental unit is at least the golden
ratio, so there are \(O(\log q)\) allowable generators per ideal.
The zero-norm case gives at most one point or a union of lines.

### Integral points on a parabolic conic

The argument for an arbitrarily oriented parabola uses a quadratic
congruence and the curvature of the band. The quadratic major arcs
treated in [Trifonov, §3.2](sources/trifonov-lattice-curve-squarefull-II-2609.17715v1.pdf)
are graphs of quadratic polynomials in the horizontal coordinate.
The local conic determinant here also allows tilted axes.

The rank-one quadratic part has the form \(A(ru+tv)^2\) with
\(A\in\mathbb Z\setminus\{0\}\) and \(\gcd(r,t)=1\). Indeed a
rational rank-one quadratic is proportional to the square of one
primitive integer linear form, and the gcd of \(r^2,2rt,t^2\) is
one. Complete \(z=ru+tv\) to an integral unimodular coordinate pair
\((z,w)\) with polynomially bounded coefficients. In these
coordinates the conic has equation
\[
                       A z^2+B z+Cw+D=0 .
\]
If \(C=0\) it is a union of at most two lines. Otherwise \(w\) is
determined by \(z\) and \(\gcd(A,B,C,D)=1\).
The condition \(Az^2+Bz+D\equiv0\bmod |C|\) places \(z\) in
\(O_\varepsilon(q^\varepsilon)\) residue classes modulo a common
integer \(Q\ge\sqrt{|C|}\). In fact, for \(p^e\parallel |C|\)
the solutions occupy at most four classes modulo
\(p^{\lceil e/2\rceil}\). If \(p\mid A\), either \(B\) is a unit
and the root lifts uniquely, or \(D\) is a unit and there are no
roots. If \(p\nmid A\) complete the square for odd \(p\).
For \(p=2\) with \(A\) odd, the odd \(B\) case has simple roots
modulo two, and the even \(B\) case permits completing the square
as \((Az+B/2)^2\). Roots of a square congruence with right-hand
side zero modulo \(p^e\) are divisible by \(p^{\lceil e/2\rceil}\).
Otherwise the right-hand side must have even valuation \(2j<e\),
and the unit square congruence modulo \(p^{e-2j}\) has at most four
roots. These lift to at most four classes modulo
\(p^{\lceil e/2\rceil}\) in the original variable. Apply the
Chinese remainder theorem and the divisor bound.

Retain every other point of the conic in order of \(u\). By the
three-point span bound, the resulting gaps are \(\gg R\).
For any triple of retained points write \(d_{ijk}=v[u_i,u_j,u_k]\)
for its second divided difference. It is comparable to \(T/N^2\).
For any four retained points of span \(W\), all ratios of their
second divided differences are \(1+O(W/N)\). Indeed on-curve
second divided differences have this property by the third
derivative bound, and the error \(O(\Delta/R^2)\) is smaller:
\[
       \frac{\Delta/R^2}{(T/N^2)(R/N)}=\Delta N=o(1).
\]
If \(t=0\), the parabola is the graph of a quadratic polynomial in
\(u\). Four retained points cannot lie on it. In fact
\(f[u_2,u_3,u_4]-f[u_1,u_2,u_3]\) has magnitude
\(\asymp (T/N^3)W\) by the third divided difference, while the
combined ordinate errors are \(O(\Delta/R^2)=o((T/N^3)W)\).

Suppose \(t\ne0\). The \(z\)-ordinates have all their second
divided differences of sign \(t\). They are distinct since \(w\)
is determined by \(z\). Thus the points split into at most two
subsequences on which \(z\) is monotone. On one of these
subsequences consider any four retained points, and set
\(t_{ij}=(z_j-z_i)/(u_j-u_i)\). The unimodular map to \((z,w)\)
and the quadratic expression for \(w\) give the exact identity
\[
                 d_{ijk}=-\frac{A}{C}t_{ij}t_{ik}t_{jk}.
\]
In particular
\[
              \frac{d_{234}}{d_{123}}=
                        \frac{t_{24}}{t_{12}}\frac{t_{34}}{t_{13}}.
\]
The two ratios on the right are either both at least one or both
at most one. This follows because the slopes of the \(v\)-ordinates
increase when both endpoints move to the right, and the \(z\)-slopes
have a common sign. Since their product is \(1+O(W/N)\), each ratio
is \(1+O(W/N)\). On the other hand
\(t_{34}-t_{13}=t\,W d_{134}\) by the definition of second divided
differences. Consequently \(P=|t_{13}|\gg |t|T/N\).

All other secant slopes of \(z\) in this subsequence have magnitude
comparable to \(P\). The \(v\)-slopes vary by
\(O(TH/N^2+\Delta/R)\) over the interval, so the \(z\)-slopes
vary by at most \(O(P(H/N+N\Delta/(TR)))=o(P)\).
Thus the \(z\)-range of the subsequence has length \(O(PH)\).
The identity for \(d_{ijk}\) also gives
\[
          \sqrt{|C|}\ \ge\ \sqrt{|C/A|}
                     \ \asymp\ (N/\sqrt T)P^{3/2}
                     \ \gg\ P\sqrt N .
\]
Together with \(H\ll\sqrt N\) and the quadratic residue classes,
this bounds the size of the subsequence by \(O_\varepsilon(q^\varepsilon)\).
Subsequences of length at most three already meet that bound.
Accounting for the discarded points changes the bound only by a
constant.

Every local conic now contributes \(O_\varepsilon(q^\varepsilon)\)
points. The balanced block therefore has
\(O_\varepsilon(q^{\sigma+\varepsilon})\) points, and adding the
unbalanced blocks and the square endpoint proves (8). The bound
\(7/30\) from using the same conics for \(N\ge q^{1/3}\) and
Filaseta–Trifonov below that cutoff is superseded numerically by (8).
For that wider conic range take \(\sigma=7/30\) in the \(b\)-aspect,
since \(N^4T^4\le q^{10/3}\). Taking the unprojected ordinate
condition alone with \(p=3/13\) gives the intermediate bound for
\(q^{4/13}\le N\le q^{9/26}\) and the earlier \(3/13\) exponent.
The uniform logarithmic-power question remains open.

### Contact calculation for the narrower balanced range

The earlier parabolic cluster estimate used the five rows of
derivatives of orders zero through four of
\(1,x,x^{-s},x^2,x^{1-s},x^{-2s}\). At \(x=1\) these are the falling
factorial evaluations at the six distinct exponents
\(\lambda=(0,1,-s,2,1-s,-2s)\). The null vector has entries
\(c_i=1/\prod_{j\ne i}(\lambda_i-\lambda_j)\) by Lagrange
interpolation. Its discriminant \(c_5^2-4c_4c_6\) is \(8/175\)
for \(s=3/2\) and \(59049/22400\) for \(s=2/3\). Separate
dilations of the two coordinates give nonvanishing at every
\(x\in[1,2]\). For five abscissae separated by \(N^{3/8+\eta}\),
the normalized divided-difference errors are
\(O(N^{-5/2}(N/N^{3/8+\eta})^4)=O(N^{-4\eta})\).
Thus sufficiently localized points with those separations force a
nonparabolic conic. In a local interval the parabolic points are
covered by four intervals of length \(N^{3/8+\eta}\), with a
three-point span lower bound \(R\) in each. This gives the earlier
\(N^{599/1000+\varepsilon}\) balanced-block estimate when
\(T\le N^{1+1/100}\), by taking local interval length
\(N^{1-11/20-2/(5\cdot100)-\eta}\).

## Uniform upper bound from two local cubics

**Theorem (Uniform upper exponent \(48/211\)).** For every
\(\varepsilon>0\),
\[
                    h(q)\ll_\varepsilon q^{48/211+\varepsilon}.
\]
For each fixed \(C>0\), the number of powerful integers in
\([x,x+C\sqrt x]\) is \(O_{C,\varepsilon}(x^{24/211+\varepsilon})\).

Use the same dyadic decomposition, with \(p=48/211\) in the
\(b\)-aspect and the three thresholds
\[
 \beta_0=10p-2=\frac{58}{211},\qquad
 \beta_*=1-3p=\frac{67}{211},\qquad
 \beta_c=\frac{2(71p-13)}{19}=\frac{70}{211}.
\]
The unbalanced estimate applies for \(\beta\le\beta_0\).
The \(a\)-aspect is unchanged, with exponent at most \(2/9\).
All statements below hold with the same window tolerance
\(\Delta\ll_C T/q\), so the proof applies to both counts.

### Extending the intermediate counts

For \(\beta_0\le\beta\le\beta_c\), use \(A=Nq^{-p}\) and the
two small-gap cutoffs and estimates above with the new \(p\).
In the upper part \(\beta_*\le\beta\le\beta_c\), the first two
exponents from (11) reach \(p\) at \(\beta_*,\beta_c\),
respectively. The other two are bounded by their values at
\(\beta_c\), which are \(46/211,104/633\). The congruence multiplicity still costs
\(O_{C,\varepsilon}(N^{1/2}q^\varepsilon)\). The size
\(K\gg q^{1/30}\) follows from \(K\gg qN^{-11/4}\) when
\(a>C_0N^{1/4}\).

In the lower part \(\beta_0\le\beta\le\beta_*\), use the same
rational-slope projection with \(Q=\lfloor q^\theta\rfloor\).
The gap domain
\[
 \frac{p-1+7\beta/2}{4}\le\alpha=\log_q a\le\beta-p
\]
is a triangle. Take \(\theta\) affine in \((\beta,\alpha)\), with
values \(1/200,\ 33/1000,\ 1/40\) at the vertices
\[
 (\beta_0,10/211),\qquad
 (\beta_*,143/1688),\qquad
 (\beta_*,19/211),
\]
respectively. The derivative check can use the numerator \(r\) in
the slope approximation. If \(r=0\), then \(H_a'''=\ell g_a'''\).
If \(r\ne0\), then \(\ell V/N\gg1\) since
\(|\ell B_a'(U_a(k))-r|\ll Q^{-1}\) and \(|B_a'|\ll V/N\).
In that case the derivative error after dividing by \(\ell\) is
\(O(N/(\ell QK^3))=O(\Lambda/Q)\). In either case the
absolute third derivative is comparable to \(\ell\Lambda\) and
has constant sign on the piece.

For the count, let \(k_0,s_0,l_0,\lambda_0\) be the powers of \(q\)
in \(K,V/N,L,\Lambda\), respectively, as given above. Set
\[
 j_0=\max(0,s_0+2\theta),\qquad m_0=\max(0,l_0),\qquad
 w_0=\max(\theta-3\beta/2,l_0-\theta,\theta+s_0-\beta+2l_0).
\]
Thus \(J\ll q^{j_0}\), \(1+L\ll q^{m_0}\), and
\(\omega'\ll_C q^{w_0}\). At the vertices,
\(k_0>\max(j_0,1/30)\), \(w_0<0\), and \(\theta+\lambda_0<0\),
and these inequalities hold throughout the triangle. Grouping
gaps dyadically costs the factor \(q^\alpha\) at each scale.
In the order used in (11), after summing over pieces, the four
exponents including this factor and the abscissa multiplicity are
\[
 \begin{split}
 &\alpha+k_0+m_0+w_0,\qquad
 \alpha+k_0+m_0+\tfrac3{19}(\theta+\lambda_0),\\
 &\alpha+\tfrac34 k_0+\tfrac14 j_0+m_0,\qquad
 \alpha+j_0-\tfrac13\lambda_0+m_0.
 \end{split}
\]
These are convex piecewise affine functions of the two coordinates.
Their maxima over the three vertices, in the same order, are
\[
 \frac{47787}{211000},\qquad \frac{455007}{2004500},\qquad
 \frac{45919}{211000},\qquad \frac{12833}{63300},
\]
each strictly below \(p\). The last term also accounts for the
short pieces. This proves the required intermediate bound.

### Canceling the single-error term of the cubic determinants

Suppose \(\beta_c\le\beta\le 2/5\) in the \(b\)-aspect, with
the endpoint enlarged by \(O(1/\log q)\) if needed for
\(N\ll q^{2/5}\). Put \(T=qN^{-3/2}\) and use intervals of
normalized length \(\rho\le q^{-p-\eta}\), as above. Write
\(x=u/N,\ y=v/T=x^{-3/2}+e\) with \(|e|\ll_C q^{-1}\).
Use the six conic columns \(1,u,v,u^2,uv,v^2\), and
adjoin either \(u^3\) or \(u^2v\). For any seven points
let \(D_0,D_1\) be the corresponding determinants after
normalization of both coordinates. The factors restoring the
integer determinants are \(N^7T^4=q^{4+\beta}\) and
\(N^6T^5=q^{5-3\beta/2}\), respectively.

Expand every row in powers of its error. For a term with \(m_i\)
rows of degree \(i\), \(i=0,1,2\), the same grouped divided
differences as for conics give the normalized bound
\[
 \ll_C q^{-m_1-2m_2}\rho^{\binom{m_0}{2}+\binom{m_1}{2}}.
\]
Here \(m_0+m_1+m_2=7\) and \(m_2\le1\).
For \(D_0\) we also have \(m_1\le3\), since there are only
three columns involving \(y\). Over those terms the least
decay exponent at \(\eta=0\) is \(1+15p>4+2/5\).
It follows that \(D_0=0\) by integrality of the restored determinant.

For \(D_1\), the first-order terms with \((m_1,m_2)=(1,0)\)
need an additional factor \(\rho\). To obtain it, take
\(x_0\) in the short interval and, for the column monomials
\(x^iy^j\), form the seven-by-seven matrix whose first six
rows are the derivatives of orders \(0,\ldots,5\) of
\(x^{i-3j/2}\) at \(x_0\) and whose last row is the values
\(j x_0^{i-3(j-1)/2}\). Write its determinant as \(F_0(x_0)\)
or \(F_1(x_0)\) for the two choices of columns. At \(x_0=1\),
the derivative entries are \(\prod_{d=0}^{r-1}(i-3j/2-d)\)
in row \(r\), and the last row is \(j\). These matrices give
\[
 F_0(1)=10276875/32,\qquad F_1(1)=1063125/64.
\]
Scaling the columns and rows at another \(x_0\in[1,2]\)
shows \(F_0(x_0)\ne0\) and \(F_1(x_0)/F_0(x_0)=O(1)\).

In a term with only the error of point \(i\), put the error
row last and apply divided differences to the other six.
After removing their Vandermonde product, the leading
determinant is \(F_j(x_0)/\prod_{r=0}^{5}r!\) for \(D_j\),
with remainder \(O(\rho)\). Indeed each divided difference
of order \(r\) differs from the corresponding derivative at
\(x_0\) divided by \(r!\) by \(O(\rho)\), and the error row
also differs from its value at \(x_0\) by \(O(\rho)\).
The row sign and Vandermonde are shared by the two determinants.
Consequently each single-error term in
\(D_1-(F_1(x_0)/F_0(x_0))D_0\) is \(O_C(q^{-1}\rho^{16})\).
The correction is a real multiple of the already vanishing
\(D_0\), so this combination equals \(D_1\).
For the other terms use \(m_1\le4\) in the grouped bound.
After the single-error cancellation every decay exponent at
\(\eta=0\) is at least \(2+11p=5-3\beta_c/2\).
This beats or equals \(5-3\beta/2\), and the positive
\(\eta\) gives strict savings. The restored integer determinant
is again zero.

If the common six conic columns have rank less than six
on the entire cell, there is already a conic relation.
Otherwise, the vanishing of both sets of seven-by-seven
minors gives two polynomials through all cell points whose
degree-three terms are nonzero multiples of \(u^3\) and
\(u^2v\), respectively. Their greatest common divisor over
\(\mathbb Q\) has degree at most two. Outside its zeros the
two remaining factors have at most nine common points by
Bézout. The points on the common divisor satisfy a conic
relation of polynomially bounded height, obtained using
minors of their six conic columns. The earlier conic counts
apply over this enlarged range: \(N^2/q=o(1)\),
\(\Delta N=o(1)\), \(T\gg N\), and \(N\rho\ll\sqrt N\).
Each cell therefore contains \(O_\varepsilon(q^\varepsilon)\)
points. There are \(O(q^{p+\eta})\) cells, which proves the
balanced-block bound. Summing dyadic blocks and counting
the left square endpoint completes the proof.

## Coprime occupancy with quadratic logarithmic height

**Theorem (Coprime occupancy with exponent one half).** For every sufficiently
large integer \(m\), there is an open interval \((q^2,(q+1)^2)\) containing
at least \(3m/8\) pairwise coprime nonsquare powerful integers, with
\(q\ge m\) and \(\log q=O(m^2\log m)\). In particular the coprime occupancy
is \(\gg(\log q/\log\log q)^{1/2}\) along an unbounded sequence of \(q\).

Use the same Laurent cancellation as below, with prime shifts of size
\(O(m\log m)\). Interpolation over the small prime completions allows
these shifts to stay small. Once the argument \(B^2\) is fixed, we use
the rank of the evaluation forms in the remaining coefficients to
bound the number of pairs that share a larger prime.

### Cancellation and small-prime interpolation

Choose \(P_0=\lceil C m\log m\rceil>4m+5\) with a sufficiently large
absolute constant \(C\). The elementary prime count below supplies \(m\)
distinct odd primes \(d_1,\ldots,d_m\le P_0\) from one residue class
modulo \(8\). Put \(n=2m\), \(c_j=\binom{-3/2}{j}\) and
\[
 R_0(T)=\prod_i(T-d_i)=\sum_{j=0}^m a_jT^j,\qquad
 S(T)=\sum_{\ell=0}^m s_\ell T^\ell,\qquad
 R_0(T)S(T)=\sum_{k=0}^n r_k T^k .
\]
The coefficients of \(S\) will be rational with power-of-two denominators.
Set
\[
 \lambda_k(Z,d)=\frac{1}{c_{k+1}}\sum_{j=0}^k c_j d^j Z^{k-j},
 \qquad
 P(Z)=\sum_{k=0}^n\frac{r_k}{c_{k+1}} Z^k,\qquad
 L_i(Z)=\sum_{k=0}^n r_k\lambda_k(Z,d_i).
\]
The coefficient of \(Z^{-1}\) in \(P(Z)(1+d_i/Z)^{-3/2}\) is
\(d_iR_0(d_i)S(d_i)=0\), and \(L_i(Z)\) is the polynomial part.
For \(E=2T\frac{d}{dT}+3\) and \(\Phi=E^{-1}(R_0S)\), we also have
\(L_i(0)=\Phi(d_i)\). Indeed \(E\) multiplies \(T^k\) by \(2k+3\),
and \(c_k/c_{k+1}=-1+1/(2k+3)\).

For an odd prime \(p\le P_0\), construct
\(\Phi_p\in\mathbb Z_p[T]\) of degree less than \(2m\) such that
\(E\Phi_p(d_i)=0\) for all \(i\) and \(\Phi_p(d_i)\) is a unit
whenever \(d_i\) is. For the shifts in a nonzero class \(a\bmod p\)
use the series \((1+(T-a)/a)^{-3/2}\) in the quotient by the product
of \((T-d_i)^2\) over that class, taking \(a\) to be an integer
representative. The binomial coefficients are \(p\)-integral. If
there are \(r\) shifts in the class, \((T-a)^{2r}\) is divisible by
\(p\) in the quotient, so this series converges there. Evaluation
modulo each \((T-d_i)^2\) gives its value and first derivative at \(d_i\).
The value is \(1\bmod p\) and the series satisfies \(Ef=0\).
For the zero class use the zero polynomial. The defining products for
distinct classes are comaximal over \(\mathbb Z_p\), since they have
unit resultants. The polynomial Chinese remainder theorem gives
\(\Phi_p\).

At \(p=2\), let \(a\) be an odd representative of the common class
modulo \(8\), and write \(U=(T-a)/8\). In the finite free
\(\mathbb Z_2\)-algebra obtained by quotienting by
\(\prod_i(U-(d_i-a)/8)^2\), use \((1+8U/a)^{-3/2}\).
Its \(j\)-th coefficient has valuation at least \(j\), so the series
converges. Its polynomial representative of degree less than \(2m\)
gives a polynomial \(\Phi_2\) in \(T\) with coefficient valuations at
least \(-6m\). Evaluation modulo the squared linear factors gives
unit values at every \(d_i\) and the same differential equation
\(E\Phi_2(d_i)=0\), by the change of variable. The quotient maps
used in both constructions are continuous on the finite free
coefficient modules, so taking the limits preserves the value and
derivative identities.

For each \(p\le P_0\), \(S_p=E\Phi_p/R_0\) is a polynomial of degree
less than \(m\). Division by the monic \(R_0\) preserves the coefficient
valuation bounds. Put
\[
 e_p=1+\lfloor\log_p(4m+3)\rfloor,\qquad D_0=2^{6m},\qquad
 K=D_0\prod_{p\le P_0}p^{e_p},\qquad X=m^{12}.
\]
Choose \(A_\ell\in[0,K)\) by the integer Chinese remainder theorem
so that
\[
 D_0s_\ell=A_\ell+Kh_\ell,\qquad X\le h_\ell<2X,\qquad h_\ell\in\mathbb Z,
\]
gives \(v_p(S-S_p)\ge e_p\) coefficientwise at every \(p\le P_0\).
For this, approximate the \(p\)-integral coefficient of \(D_0S_p\)
to precision \(p^{e_p+v_p(D_0)}\). Applying \(E^{-1}\) loses at
most \(e_p-1\) powers of \(p\). Thus every \(L_i(0)\) is integral at
these primes, with residue \(\Phi_p(d_i)\).

Require \(B\) to be divisible by \(K_B=\prod_{p\le P_0}p^{b_p}\),
where \(b_2=10m+10\) and \(b_p=1+2\lfloor\log_p(4m+3)\rfloor\)
at odd \(p\). Then \(L_i(B^2)\equiv L_i(0)\bmod p\) at all these
primes. To see this, each term of positive \(Z\)-degree in
\(L_i(Z)\) has valuation at least \(1\) upon substituting \(Z=B^2\).
At odd \(p\) the \(r_k,c_j\) are integral and
\(v_p(c_{k+1})\le 2\lfloor\log_p(4m+3)\rfloor\), by the binomial
formula and the factorial valuation formula. At \(2\) we have
\(v_2(r_k)\ge-6m\), \(v_2(c_j)\ge-2j\) and \(v_2(c_{k+1})\le0\),
also by the binomial formula. Beyond \(P_0\), the denominators in
\(L_i\) are units. Consequently \(L_i(B^2)\) is an integer. Also,
every prime \(p\le P_0\) can divide at most one of
\((B^2+d_i)L_i(B^2)\), namely the one with \(d_i=p\) if such a
shift is present.

### Choosing the square argument

The primorial bound below gives \(\log K,\log K_B=O(m\log m)\).
Indeed the extra powers beyond the primorial contribute at most
\(O(m\log m)\) in total, since only \(p\le4m+3\) contribute them.
The coefficient heights of \(P\) and all \(L_i\) are bounded by a
uniform \(H=\exp(O(m\log m))\) for all the allowed \(h_\ell\).
This follows from \(\sum |a_j|\le\prod_i(1+d_i)\), the bounds
on \(s_\ell\) and \(d_i\), and \(1\le |c_j|\le 2j+1\).
The common leading coefficient has modulus
\(|s_m/c_{n+1}|\ge 1/(2n+3)\), since \(s_m\ge1\).
Hence the Laurent remainder estimate below supplies a threshold
of logarithm \(O(m\log m)\) above which \(B^2+d_i\) is strictly
between consecutive squares, \(P(B^2)\) and all \(L_i(B^2)\)
are nonzero, and
\[
 \big|(B^2+d_i)^{3/2}|L_i(B^2)|-|B^3P(B^2)|\big|<\tfrac12 .
\]
The same threshold can ensure \(|B^3P(B^2)|\ge m+2\), using
leading-term dominance uniformly in \(S\).

Choose \(B\) beyond this threshold, with \(K_B\mid B\) and with
no prime \(P_0<p\le Y=m^3\) dividing any \(B^2+d_i\).
Write \(B=K_B b\). For each of these primes at most \(2m\) classes
of \(b\) are forbidden, so their densities \(v_p\) satisfy
\(v_p<1/2\). The bound \(\pi(x)\ll x/\log x\) and partial summation
give \(\sum_{P_0<p\le Y}1/p=O(1)\). Thus
\(V=\prod_p(1-v_p)\ge \exp(-2\sum_p v_p)=\exp(-O(m))\).
As in the progression estimate below, truncate inclusion-exclusion
at a sufficiently large odd order \(O(m)\). With \(S_0=\sum_p v_p\),
the absolute main-term tail is bounded by
\(\sum_{j>k}S_0^j/j!\) at order \(k\), which is at most \(V/2\)
for \(k\) a sufficiently large odd integer of size \(O(S_0+1)\).
The integer counting error is at most the sum of the forbidden
class counts for the subsets used, bounded by
\(\sum_{j\le k}(2mY)^j=\exp(O(m\log m))\).
An interval of \(b\) of length \(\exp(O(m\log m))\), starting beyond
the prescribed threshold, then contains an admissible \(b\). This
gives \(\log B=O(m\log m)\). The numbers \(B^2+d_i\) are pairwise
coprime. At \(p\le P_0\) only \(d_i=p\) can contribute, and at
\(p>P_0\) the differences of the shifts are nonzero modulo \(p\).

### Rank and large-prime collisions

Fix \(Z=B^2\). The \(m\) forms \(L_i(Z)\) in the coefficients of \(S\)
have rank \(m\) over every \(\mathbb F_p\) with \(p>P_0\).
For a relation with weights \(w_i\), put
\(\mu_k=\sum_i w_i\lambda_k(Z,d_i)\). Vanishing of the relation
on coefficients of \(S\) says \(\sum_{j=0}^m a_j\mu_{j+\ell}=0\)
for \(0\le\ell\le m\). Since \(R_0\) is monic and the shifts
are distinct modulo \(p\), the Vandermonde matrix on the first
\(m\) terms and this recurrence yield
\(\mu_k=\sum_i C_i d_i^k\) for \(0\le k\le2m\).
The identity
\[
 (2k+3)\lambda_k(Z,d_i)+2(k+1)Z\lambda_{k-1}(Z,d_i)
   =-2(k+1)d_i^k,\qquad 1\le k\le2m,
\]
therefore implies
\[
 \sum_i\bigl((2k+3)C_i+2(k+1)(ZC_i/d_i+w_i)\bigr)d_i^k=0
 \qquad(1\le k\le2m).
\]
The sequences \(d_i^k\) and \(k d_i^k\) are independent on these
indices. Indeed, relabeling \(k=j+1\) identifies relations
between them with relations between evaluations and first
derivatives on polynomials of degree less than \(2m\), since the
\(d_i\) are nonzero. These evaluation and derivative functionals
are independent by the Chinese remainder theorem modulo
\(\prod_i(T-d_i)^2\). For each \(i\), subtracting the coefficient
of \(k\) in the last displayed relation from its constant
coefficient gives \(C_i=0\), and then \(w_i=0\).

After the substitution \(D_0s_\ell=A_\ell+Kh_\ell\), clear
denominators by a common factor supported at primes at most \(P_0\),
for instance \(D_0\,4^{2n+1}(2n+3)!\). Denote the resulting
integer affine forms by \(G_i\). Their linear parts still have
rank \(m\) over every \(p>P_0\), and hence over \(\mathbb Q\).
Their coefficient logarithmic heights are \(O(m^2\log m)\)
because \(Z=B^2\) and their degree before evaluation at \(Z\)
was at most \(2m\). Primes \(p>P_0\) divide \(L_i(Z)\) exactly
when they divide \(G_i\).

Choose the \(h_\ell\) independently and uniformly in their
integer ranges. For each pair \(i,j\) and each prime
\(P_0<p\le X\), use a nonsingular two-column minor of the two
linear parts and fix all other coordinates. The probability
that \(p\) divides both \(G_i,G_j\) is at most
\((1/p+1/X)^2\le 4/p^2\). The expected number of pairs with a
collision of this kind is \(O(m^2/P_0)=o(m)\).

For \(p>X\), any collision has some coordinate in which at least
one of the two linear coefficients is nonzero modulo \(p\).
For each such candidate coordinate \(\ell\), fix the other
coordinates, and write the forms as \(u h_\ell+v\) and
\(u' h_\ell+v'\). We can skip coordinates with \(u=u'=0\).
The eliminant \(u'v-uv'\) is a nonzero affine form in the other
coordinates, since the two linear parts are independent over
\(\mathbb Q\). Its probability of taking the value zero is at
most \(1/X\). Otherwise its logarithmic size is \(O(m^2\log m)\),
so it has \(O(m^2)\) prime divisors larger than \(X\). For each
such prime at which at least one of \(u,u'\) is a unit, at most one
choice of \(h_\ell\) gives the collision. Summing over
coordinates gives \(O(m^3/X)\) per pair.

Finally consider a prime \(p>P_0\) dividing \(Z+d_i\) and \(G_j\)
for \(j\ne i\). It must be greater than \(Y\). For \(p\le X\),
rank one of the single linear part gives probability at most
\(2/p\). For each \(i\) the sum of \(1/p\) over such divisors of
\(Z+d_i\) is bounded by \(\log(Z+d_i)/(Y\log Y)\).
For \(p>X\) the probability is at most \(1/X\) and the number
of these divisors is at most \(\log(Z+d_i)/\log X\).
Summing over ordered pairs gives \(O(m^3/Y+m^3/X)\).
Together with pairwise coprimality of the bases and the
small-prime conditions, the expected number of pairs sharing a
prime in \((Z+d_i)L_i(Z)\) is \(o(m)\).

Take a choice with at most \(m/4\) such pairs, and remove at most
one index per pair to retain at least \(3m/4\) pairwise coprime
powerful integers \((Z+d_i)^3L_i(Z)^2\). Each is nonsquare
because \(Z+d_i\) is nonsquare. Their square roots lie within
\(1/2\) of the same real number, so they meet at most two open
unit intervals with integer endpoints. At least \(3m/8\) lie
in one such interval \((q,q+1)\). The threshold ensured \(q\ge m\),
and the coefficient and argument bounds give
\(\log q=O(n\log B+\log H)=O(m^2\log m)\). This proves the
claimed lower bound since \(x/\log x\) is increasing for \(x>e\).
Taking a sufficiently large constant multiple of a prescribed
multiplicity parameter gives the hypotheses and conclusion of the
cubic logarithmic height theorem below.

## Earlier cubic logarithmic height bound

**Theorem (Coprime occupancy with exponent one third).** For each integer
\(m\ge2\), there is an interval \((q^2,(q+1)^2)\) containing at least
\(\lceil m/2\rceil\) pairwise coprime nonsquare powerful integers, where
\[
 m\le q,\qquad \log q\ll m^3\log(m+1)
\]
with an absolute constant. In particular the coprime occupancy is
\(\gg \bigl(\log q/\log\log q\bigr)^{1/3}\) along an unbounded
sequence of \(q\).

In the construction below, small primes are handled by scaling the shifts.
For the intermediate primes, extra coefficients in the cancellation
polynomial allow all the values to be units at a prescribed argument.
This reduces the size of the shifts, which determines how far out the
square-root approximations become accurate.

Let \(n=2m+1\), \(d=2m(n+1)\), \(C=8(m+1)\), and put
\[
 A=4^{2n+1}(2n+3)!,\qquad M=A\prod_{p\le C}p.
\]
Use \(m\) shifts \(d_i=M t_i\), where the \(t_i\) will be primes
larger than \(2d\) and of size \(m^{O(1)}\). Let
\[
 \begin{split}
 S(T)&=T^{m+1}+\sum_{\ell=0}^m s_\ell T^\ell,\qquad
 R(T)=T\prod_i(T-d_i)S(T)=\sum_{k=0}^n r_{k+1}T^{k+1},\\
 P(Z)&=\sum_{k=0}^n p_k Z^k,\qquad
 p_k=c_{n+1}r_{k+1}/c_{k+1},\qquad c_j=\binom{-3/2}{j},\\
 L_i(Z)&=\sum_{k=0}^n p_k\sum_{j=0}^k c_j d_i^j Z^{k-j},
 \qquad F_i(Z)=(Z+d_i)L_i(Z).
 \end{split}                                                     \tag{942.1}
\]
The three polynomials in \(Z\) are monic. Since \(R(d_i)=0\),
the same Laurent coefficient calculation as in the construction below
shows that the \(Z^{-1}\) term of \(P(Z)(1+d_i/Z)^{-3/2}\) vanishes.
Take all \(s_\ell\equiv0\bmod M\). Then all coefficients of \(R\)
below the leading one are multiples of \(M\). The denominators of the
rational factors \(c_{n+1}c_j/c_{k+1}\) in (942.1) divide \(A\),
by the binomial formula for \(c_j\). Each nonleading term in a
coefficient of \(P\) or \(L_i\) has a factor \(r_{k+1}\) with \(k<n\)
or a positive power of \(d_i\). The polynomials are therefore integral
and all their nonleading coefficients are divisible by \(M/A\).
This holds also for \(F_i\). Thus \(B\equiv1\bmod p\) keeps all
\(F_i(B^2)\) nonzero modulo every \(p\le C\).

### Choosing the shifts and the intermediate prime conditions

First work with \(S(T)=T^{m+1}\). For each pair \(i\ne j\) the
resultant of \(F_i,F_j\) is a nonzero polynomial in the shifts.
Indeed, specialize only \(d_i=a\ne0\) and set the other shifts to zero.
With \(\lambda=c_{n+1}/c_n\), this gives
\(P(Z)=Z^{n-1}(Z-\lambda a)\) and \(F_j(Z)=Z^n(Z-\lambda a)\).
The polynomial part formula gives
\[
 L_i(0)=a^n(c_n-\lambda c_{n-1}),\qquad
 L_i(\lambda a)=c_n a^n .
\]
Both values are nonzero since \(c_j/c_{j-1}\) varies with \(j\).
Also \(Z+a\) is nonzero at both arguments since \(\lambda\ne-1\).
The product of the pairwise resultants has total degree \(O(m^4)\)
in the shifts, using a Sylvester matrix of size \(2(n+1)\) and the
bound \(n+1\) on each coefficient degree. Substituting \(d_i=M t_i\)
preserves nonvanishing as a polynomial in the \(t_i\).
The grid lemma below therefore allows a choice of the \(t_i\)
among any \(O(m^4)\) distinct primes above \(2d\), with the implicit
constant chosen larger than the degree bound. There are sufficiently
many such primes of size \(m^{O(1)}\). Indeed
\(\binom{2u}{u}\ge4^u/(2u+1)\), and the exponent of each prime \(p\)
in this binomial coefficient is at most \(\lfloor\log_p(2u)\rfloor\):
each power \(p^j\le2u\) contributes zero or one by the factorial
valuation formula. Comparing logarithms gives
\(\pi(2u)\gg u/\log(2u)\). Choose \(u\) a sufficiently large
fixed power of \(m\) times an absolute constant. Fix the resulting
vector of shifts, for which the pairwise resultants are nonzero at
\(s_0=\cdots=s_m=0\).

Fix now a prime \(C<p\le2d\). The \(d_i\) and the \(c_j\),
\(0\le j\le n+1\), are units modulo \(p\). We will take
\(B\equiv0\bmod p\), so it suffices to make each \(L_i(0)\) nonzero.
Write
\[
 T\prod_j(T-d_j)=\sum_{k=0}^m a_k T^{k+1},\qquad a_m=1 .
\]
For \(0\le\ell\le m\), the coefficient of \(s_\ell\) in \(L_i(0)\)
is
\[
 c_{n+1}d_i^\ell
 \sum_{k=0}^m a_k d_i^k\frac{c_{k+\ell}}{c_{k+\ell+1}}
 =
 c_{n+1}d_i^\ell
 \sum_{k=0}^m\frac{a_k d_i^k}{2(k+\ell)+3}.                 \tag{942.2}
\]
The constant \(-1\) in the quotient
\(c_j/c_{j+1}=-1+1/(2j+3)\) disappears because
\(\sum_k a_k d_i^k=0\).
If all \(m+1\) coefficients (942.2) vanished, the rational function
\[
 X\longmapsto\sum_{k=0}^m\frac{a_k d_i^k}{X+2k+3}
\]
would vanish at \(X=0,2,\ldots,2m\) in \(\mathbb F_p\).
None of these points is a pole because \(p>C\). Its numerator over
the common denominator has degree at most \(m\), so it would be the
zero polynomial. The poles \(-2k-3\) are distinct, and the numerator
evaluated at \(-2m-3\) is \(d_i^m\prod_{k<m}(2k-2m)\ne0\).
This contradiction shows that \(L_i(0)=0\) excludes an affine
hyperplane of coefficient vectors modulo \(p\). The \(m\) hyperplanes
together have at most \(m p^m<p^{m+1}\) points. Choose a vector
outside their union for each intermediate prime.

The prescribed residues of \(\mathbf s\), including \(0\bmod M\),
give a Cartesian grid of progressions with common step
\(K=M\prod_{C<p\le2d}p\), with \(\log K=O(m^2)\) by the
primorial bound below. For our fixed shifts the resultant product is
a nonzero polynomial in \(\mathbf s\), as seen at the zero vector.
Its total degree is \(O(m^3)\) since the coefficients of each
\(F_i\) are affine in \(\mathbf s\). Take the first sufficiently many
nonnegative integers in each progression and apply the grid lemma
again. We obtain \(\log(2+\max_\ell|s_\ell|)=O(m^2)\), preserving all
intermediate prime conditions and all nonzero resultants.

### Choosing the square interval

The chosen shifts satisfy \(\log(2+\max_i d_i)=O(m\log(m+1))\).
It follows from (942.1) that the logarithms of all coefficient
heights of \(P,L_i,F_i\) are \(O(m^2\log(m+1))\). To see this,
the product \(T\prod_i(T-d_i)\) has coefficient sum bounded by
\(\prod_i(1+d_i)\), its product with \(S\) costs at most the
coefficient sum of \(S\), and the polynomial part formula uses
at most \(O(m^2)\) terms with one power \(d_i^j\), \(j\le n\).
The bounds \(1\le|c_j|\le2j+1\) control the other factors.

Let \(R_1\) be the product of the absolute values of the pairwise
resultants of the chosen integral \(F_i\). The Sylvester determinants
give \(\log(2R_1)=O(m^5\log(m+1))\). Enforce the chosen residues
\(B\equiv1\bmod p\) for \(p\le C\) and \(B\equiv0\bmod p\) for
\(C<p\le2d\). Their common step is \(\prod_{p\le2d}p\), of logarithm
\(O(m^2)\). A prime \(p>2d\) outside \(R_1\) cannot divide two
different \(F_i(B^2)\). For the remaining primes \(p>2d\) dividing
\(R_1\), exclude the roots of the monic polynomial
\(\prod_i F_i(B^2)\) of degree \(d\). There are at most \(d\)
forbidden residues per prime. The progression estimate below
uses only the step, the bound \(d\) on the number of excluded classes,
and \(J=1+\lceil\log_2 R_1\rceil\) as an upper bound for the number
of primes. In particular it works for the present prescribed residue
just as for \(B\equiv1\) on the small primes. In its notation \(Y=100(d+1)J=m^{O(1)}\), and only primes
\(2d<p\le Y\) contribute to the small-prime sum \(S\).
The primorial bound gives \(\pi(x)\ll x/\log x\) by counting primes
above \(\sqrt x\) with weight at least \((\log x)/2\).
Partial summation therefore gives
\(\sum_{2d<p\le Y}1/p=O(1)\) because \(\log Y/\log(2d)=O(1)\).
Thus \(S=O(d)\), and the error factor in the progression estimate
has \(\log E=O(d\log(m+1))\). With this bound the estimate supplies
an admissible \(B\) within \(\exp(O(m^2\log(m+1)))\) beyond any
positive threshold.

Use the Laurent remainder estimate below with degree \(n\), shift bound
\(D=\max_i d_i\) and coefficient bound \(H\) for \(P,L_i\).
A threshold \(B\) greater than a sufficiently large constant times
\((n+2)^2 H D^{n+2}\) ensures the common square-root error is less than
one, and its logarithm is \(O(m^2\log(m+1))\).
It also ensures that every \(B^2+d_i\) is strictly between \(B^2\)
and \((B+1)^2\) and that all polynomial values in (942.1) are
positive. We can therefore choose \(B\) with
\(\log B=O(m^2\log(m+1))\).

The integers \((B^2+d_i)^3 L_i(B^2)^2\) are pairwise coprime
since their prime supports are those of the \(F_i(B^2)\). They
are powerful and nonsquare, and their square roots differ from
the common integer \(Q=B^3P(B^2)\) by less than one.
At least half lie in one of the two open intervals of length one
with endpoint \(Q\). In either case its lower endpoint \(q\) is
at least \(m\) by the choice of threshold, and
\(\log q=O(n\log B+\log H)=O(m^3\log(m+1))\).
This proves the theorem. It gives a lower bound for coprime occupancy;
the logarithmic uniform upper bound in Erdős #942 remains open.

## Growing coprime multiplicity

**Theorem (Growing coprime occupancy).**
For each integer \(m\ge2\) there is an interval \((q^2,(q+1)^2)\)
containing at least \(\lceil m/2\rceil\) pairwise coprime nonsquare
powerful integers, with
\[
 m\le q,\qquad \log q\ll m^4
\]
with an absolute constant. Consequently along an unbounded sequence
of \(q\) the coprime occupancy is
\(\gg(\log q)^{1/4}\).
In particular the lower bound \((\log q)^{1/5}\) holds for all
sufficiently large \(q\) in that sequence.

We seek integers of the form \((B^2+d_i)^3L_i(B^2)^2\), with one
common integer square root approximation \(B^3P(B^2)\).
Truncating a binomial expansion alone leaves an error of order
\(B\) in the square root. We choose \(P\) so the first discarded
coefficient vanishes for all \(i\). This reduces the error to
\(O(1/B)\) with the construction data fixed. Resultants then let
us choose \(B\) with all the powerful values pairwise coprime.

### Cancellation at all selected shifts

Let \(c_j=\binom{-3/2}{j}\), so \(c_0=1\) and
\[
 c_j=(-1)^j(2j+1)\binom{2j}{j}/4^j,\qquad
 c_{j+1}/c_j=-(2j+3)/(2j+2).
\]
For a vector of shifts \(\mathbf d=(d_1,\ldots,d_m)\), define
\[
 R(T)=T\prod_i(T-d_i)=\sum_{j=1}^{m+1}r_jT^j,\qquad
 P(Z)=\sum_{k=0}^m p_k Z^k,\qquad p_k=c_{m+1}r_{k+1}/c_{k+1}.
\]
In particular \(P\) is monic and \(p_k\) is homogeneous of degree
\(m-k\) in the shifts. Define the monic polynomial
\[
 L_i(Z)=\sum_{k=0}^m p_k\sum_{j=0}^k c_j d_i^j Z^{k-j},\qquad
 F_i(Z)=(Z+d_i)L_i(Z).
\]
These polynomials depend on \(\mathbf d\). In the Laurent expansion
of \(P(Z)(1+d_i/Z)^{-3/2}\) at infinity, \(L_i\) is the polynomial
part. The coefficient of \(Z^{-1}\) is
\(\sum_k p_k c_{k+1}d_i^{k+1}=c_{m+1}R(d_i)=0\).

For example, with \(m=2\) and \((d_1,d_2)=(24,48)\), this yields
\[
 P=Z^2+84Z+1680,\quad L_1=Z^2+48Z-264,\quad L_2=Z^2+12Z-48 .
\]
At \(B=90725\) the two resulting powerful integers are coprime
and both exceed \((B^3P(B^2)-1)^2\) while lying below
\((B^3P(B^2))^2\). These claims can be checked by evaluating the
three integer polynomials and taking a gcd. The proof for arbitrary
\(m\) follows.

### Making the resultants nonzero

For any pair \(i\ne j\), the resultant in \(Z\) of \(F_i,F_j\)
is a nonzero polynomial in the shifts. To prove this, specialize
\(d_i=a\ne0\) and all other shifts to zero. Put
\(\lambda=c_{m+1}/c_m\). Then
\[
 P(Z)=Z^{m-1}(Z-\lambda a),\qquad
 F_j(Z)=Z^m(Z-\lambda a),
\]
so it suffices to check \(F_i\) at \(0\) and at \(\lambda a\).
The polynomial part formula gives
\[
 L_i(0)=a^m(c_m-\lambda c_{m-1}),\qquad
 L_i(\lambda a)=c_m a^m .
\]
For the second equality the two truncated sums from the two
terms of \(P\) cancel except for the term of index \(m\).
The first quantity is nonzero because \(c_m/c_{m-1}\ne\lambda\).
The factors \(Z+a\) are nonzero at both points because
\(a\ne0\) and \(\lambda\ne-1\). This proves the assertion about
the resultant.

The coefficient of \(Z^k\) in \(F_i\) has degree at most
\(m+1\) in the shifts. Its Sylvester determinant with \(F_j\)
thus has degree at most \(2(m+1)^2\). The product over pairs
has degree at most \(m(m-1)(m+1)^2\). A nonzero polynomial of
total degree \(D\) cannot vanish on a Cartesian grid with more
than \(D\) entries in each coordinate. This follows by expanding
in one coordinate, using the univariate root bound, and inducting
on the number of coordinates. Thus we can choose positive integers
\(t_i\le 1+(m+1)^4\) for which all the resultants in the vector
\(\mathbf t\) are nonzero. In particular the \(t_i\) are distinct.

### Choosing the integer parameter

Put
\[
 A=4^{2m+1}(2m+3)!,\qquad d=2m(m+1),\qquad
 M=A\prod_{p\le 2d}p,\qquad d_i=Mt_i.
\]
For the vector \(\mathbf t\), all coefficients of \(P,L_i,F_i\)
have denominator dividing \(A\). Indeed
\((2k+3)\binom{2k+2}{k+1}\) divides \((2m+3)!\) for
\(0\le k\le m\), and the powers of \(4\) from
\(c_{m+1}/c_{k+1}\) and an additional \(c_j\) divide
\(4^{2m+1}\). This also clears the denominators in the sums.
Scaling \(\mathbf t\) by \(M\) multiplies the coefficient of
\(Z^k\) by \(M^{m-k}\) in \(P,L_i\), and \(M^{m+1-k}\) in
\(F_i\). Thus the polynomials for \(\mathbf d\) have integer
coefficients, and all the coefficients below the leading one are
divisible by every prime at most \(2d\).

Use \(F_i(Z;\mathbf t)\) to denote the unscaled polynomials, and put
\[
 R_0=\prod_{i<j}
 \left|\operatorname{Res}_Z(A F_i(Z;\mathbf t), A F_j(Z;\mathbf t))\right|.
\]
This is a positive integer. For a prime \(p\) not dividing \(M R_0\),
no two of the integers \(F_i(B^2;\mathbf d)\) can be zero modulo
\(p\). Otherwise \(B^2/M\) modulo \(p\) would be a common root of
the corresponding unscaled polynomials, by homogeneity and
\(p\nmid M\). Since \(A\) is also a unit, the prime would divide
\(R_0\).

Every prime dividing \(M\) is at most \(2d\). Impose
\(B\equiv1\bmod \operatorname{rad}(M)\). Then all
\(F_i(B^2;\mathbf d)\) are one modulo the primes dividing \(M\).
For every other prime dividing \(R_0\), exclude the zeros
of the monic polynomial \(\prod_i F_i(B^2;\mathbf d)\) in \(B\).
It has degree \(d\), so at most \(d\) classes are forbidden
modulo any such prime \(p>2d\).

For the unscaled vector, each coefficient of \(R(T)\) has modulus
at most \(\prod_i(1+t_i)\), and \(1\le |c_j|\le2j+1\).
The formulas for \(P,L_i,F_i\) give
\(\log\max(2,H(A F_i(Z;\mathbf t)))=O(m\log(m+1))\),
where \(H\) is maximum absolute coefficient.
The Sylvester determinants have size \(2(m+1)\), so
\(\log R_0=O(m^4\log(m+1))\). This estimate uses the unscaled
polynomials, which keeps \(M\) out of the product of resultants.
Also \(\log M=O(m^2)\). To see the primorial bound, the product of
primes in \((l,2l]\) divides \(\binom{2l}{l}\le4^l\).
Sum the logarithmic bound over intervals between powers of two.

Write \(B=1+\operatorname{rad}(M)b\) and let \(b\) run through
\(T\) consecutive integers. Each \(p\mid R_0\), \(p\nmid M\)
excludes \(b_p\le d\) residue classes of \(b\). Put
\[
 J=1+\lceil\log_2 R_0\rceil,\qquad Y=100(d+1)J .
\]
There are at most \(J\) primes to consider. For those at most \(Y\)
put \(g_p=b_p/p\), \(S=\sum g_p\) and \(V=\prod(1-g_p)\).
We have \(g_p<1/2\), \(S=O(d\log Y)\) by the harmonic sum,
and \(V\ge e^{-2S}\). By the Chinese remainder theorem any
\(j\) prescribed prime exclusions intersect in a fraction
\(\prod g_p\) of classes modulo their product. The same
fraction holds in a progression whose step is coprime to all
those primes. In \(L\) successive terms of such a progression
the count differs from \(L\prod g_p\) by at most \(\prod b_p\).

Use inclusion–exclusion truncated at consecutive even and
odd orders of size at least \(10(S+1)\) and at most \(10(S+1)+3\).
The main term tails are at most
\(\sum_{j\ge10(S+1)}S^j/j!\le V/8\), by \(j!\ge(j/e)^j\)
and \(V\ge e^{-2S}\). The total counting error is at most
\(E=(1+dY)^{O(S+1)}\), with an absolute constant in the exponent,
by summing \(\prod b_p\) up to those orders. Thus the count
passing the primes at most \(Y\) is between \((7/8)VL-E\)
and \((9/8)VL+E\), uniformly over those progressions.
Take \(E\ge1\). Its logarithm is \(O(m^2(\log(m+1))^2)\).

For \(Y<p\le\sqrt T\), apply the upper small prime count to
each forbidden step \(p\) progression with \(L\le T/p+1\).
Summing the main terms costs at most \(2VTdJ/Y+2VdJ\),
and the errors cost at most \(EdJ\). The primes exceeding
\(\sqrt T\) cost at most \(dJ(\sqrt T+1)\), even without
a small prime restriction. Take
\[
 T=\left\lceil\big(100(d+1)(J+1)(E+1)/V\big)^2\right\rceil .
\]
The total exclusions are smaller than the lower small prime
count at \(L=T\), since \(dJ/Y\le1/100\) and
\(V\sqrt T\ge100(d+1)(J+1)(E+1)\). An admissible \(B\)
therefore occurs within \(O(\operatorname{rad}(M)T)\) beyond
any prescribed positive threshold. The logarithm of this
length is \(O(m^2(\log(m+1))^2)\).

### Passing to one open interval

Set \(D=\max_i d_i\) and let \(H\ge1\) be the maximum coefficient
modulus over \(P,L_i\) for \(\mathbf d\). The same coefficient
bounds give \(\log(2D)=O(m^2)\) and
\(\log(2H)=O(m^3)\). For \(Z=B^2\ge2D\) the binomial
series converges absolutely. The contribution to the remainder
from the \(k\)-th term after omitting its polynomial part
and the \(Z^{-1}\) term is at most
\[
 |p_k| Z^k\sum_{j\ge k+2}(2j+1)(d_i/Z)^j
 \ll (k+2)|p_k|D^{k+2}Z^{-2}.
\]
The \(Z^{-1}\) terms sum to zero by the cancellation identity.
Consequently for positive \(B\), as soon as \(L_i(B^2)>0\),
\[
 \left| B^3P(B^2)-(B^2+d_i)^{3/2}L_i(B^2)\right|
 \ll (m+2)^2 H D^{m+2}/B .
\]
Choose a common threshold for \(B\) larger than an absolute
constant times \((m+2)^2 H D^{m+2}\). This also ensures
\(d_i<2B+1\), \(B>m+2\), and positivity of \(P(B^2)\) and
all \(L_i(B^2)\), by monicity. The logarithm of the threshold
is \(O(m^3)\). Apply the preceding progression count beyond this
threshold. Its length bound gives \(\log B=O(m^3)\).

Each \(T_i=(B^2+d_i)^3L_i(B^2)^2\) is powerful. The base \(B^2+d_i\)
lies strictly between \(B^2\) and \((B+1)^2\), so \(T_i\) is
nonsquare. Its prime support equals that of \(F_i(B^2;\mathbf d)\).
Thus the \(T_i>1\) are pairwise coprime and distinct. Their
square roots lie at distance strictly less than one from the
integer \(Q=B^3P(B^2)\), and none equals \(Q\). At least half
lie on the same side of \(Q\). The two possible intervals
have lower square roots \(q=Q-1\) and \(q=Q\).
Since \(P\) is monic and the threshold makes its leading term
dominate, \(q\ge m\). Also \(\log Q=O(m\log B+\log H)=O(m^4)\).
This proves the theorem.

The argument is an elementary polynomial construction with an elementary
coprimality sieve. It gives no logarithmic upper bound on \(h(q)\).

## Limiting density of large coprime occupancy

**Theorem (Coprime occupancy density).** Let \(d_m^{\mathrm{cop}}\) be
the natural density of integers \(q\ge1\) such that \((q^2,(q+1)^2)\)
contains \(m\) or more pairwise coprime powerful integers. Define
\(d_m\) by the same count with no coprimality test. Both densities
exist and are positive for every \(m\ge1\). Put
\[
 s=\frac32,\qquad A=\int_0^\infty\frac{du}{1+u^{s}}
     =\frac{4\pi}{3\sqrt3},\qquad \delta=\frac{6}{\pi^2}.
\]
Then, as \(m\to\infty\),
\[
\begin{aligned}
 \log d_m^{\mathrm{cop}}
   &=-s m\log m-(s+1)m\log\log m
      +\bigl(s(1+\log A)-\gamma\bigr)m+o(m),\\
 \log d_m
   &=-s m\log m+s(1+\log(\delta A))m+o(m),
\end{aligned}                                                    \tag{942.D1}
\]
where \(\gamma\) is Euler's constant. The density limit is taken at
fixed \(m\) before this asymptotic in \(m\). The ordinary limiting
distribution is due to Shiu and was generalized to successive
\(k\)-th powers by Xiong and Zaharescu.
[Shiu, Mathematika 27 (1980), 171–178](https://doi.org/10.1112/S0025579300010056)
[Xiong–Zaharescu, Indag. Math. 22 (2011), 77–86](https://doi.org/10.1016/j.indag.2011.08.003)

### Fixed kernels and the coprimality sieve

Index intervals by their upper square \(n^2\), with \(n=q+1\).
For each squarefree \(b\ge2\), let \(w_b=b^{-s}\) and
\(a_b(n)=\lfloor nw_b\rfloor\). The unique possible occupant
with cube kernel \(b\) is \(a_b(n)^2b^3\). It occurs exactly
when \(0<\{nw_b\}<w_b\). The probability in the density limit of
any specified presence and absence pattern on a finite set of kernels
is the product of \(w_b\) over the present kernels and \(1-w_b\) over
the absent kernels. Indeed, \(1\) and the \(w_b=\sqrt b/b^2\) are
linearly independent over \(\mathbb Q\). In the multiquadratic field
generated by the square roots of the primes dividing the kernels,
the \(w_b\) have distinct nontrivial sign characters. Weyl's criterion
gives equidistribution on the finite torus of fractional parts.
Applying the criterion to \(\{nw_b/D\}\) for any fixed integer \(D\)
also makes each \(a_b(n)\bmod D\) uniform independently of the
presence tests and of the other residues. In fact the integer and
fractional parts of \(D\{nw_b/D\}\) yield the required residue and
fractional part, respectively.

Let \(S\) be a fixed set of \(m\) present kernels. If two share a
prime, their occupants cannot be pairwise coprime. If the kernels
are pairwise coprime, the additional density factor for coprimality is
\[
 C_m(S)=C_m\prod_{p\mid\prod_{b\in S}b}
                   \bigl(1+(m-1)/p\bigr)^{-1},\qquad
 C_m=\prod_p (1-1/p)^{m-1}\bigl(1+(m-1)/p\bigr).                \tag{942.D2}
\]
For \(p\) outside the kernel product, at most one of the \(m\)
floors may vanish modulo \(p\). If \(p\) divides one kernel, the
other \(m-1\) floors must be nonzero modulo \(p\). The factors follow
by applying the fixed-modulus equidistribution with the product of
the primes used in the test.

To pass to all primes, fix two distinct kernels \(b,c\) and put
\(\beta=w_c/w_b\), a quadratic irrational. Writing the two positive
floors as \(du,dv\) when they have a common divisor \(d\), one has
\(u\ll_b X/d\) and \(|\beta u-v|\ll_{b,c}1/d\) for \(n\le X\).
The inequality \(\|\beta j\|\gg_\beta 1/j\) follows by clearing
denominators in the quadratic equation of \(\beta\) and evaluating it
at a nearest rational with denominator \(j\). It gives \(d\ll_{b,c}
\sqrt X\). It also separates the points \(\{\beta u\}\) for
\(u\ll_b X/d\) by \(\gg_{b,c}d/X\) on the circle. Thus there are
\(O_{b,c}(X/d^2+1)\) choices of \(u\) at each such \(d\), with
\(O_b(1)\) integers \(n\) for each \(d,u\) since the first floor
is fixed. Summing over \(d>P\) bounds the excluded count by
\(O_{b,c}(X/P+\sqrt X)\). This argument also applies if one specifies
the absence of any other finitely many kernels, since ignoring
such tests only increases the exclusion bound. Taking \(P\) past all
prime divisors of the kernels and then to infinity proves (942.D2).

For a kernel cutoff \(B\), the occupancy events are finite unions
of the finite tests just analyzed. Their natural densities exist
by first using the modular tests and then bounding the discrepancy
by the finite sum of large-divisor errors. Among \(n\le X\) there
are \(O(X B^{-1/2})\) centres with any kernel \(b>B\), since each
positive integer \(a\) with \(a b^s<X\) places its root in just one
open unit interval. This bound permits passage from the finite to
the infinite events. It also shows that the event for which exactly
\(S\) is present and its occupants are pairwise coprime has density
\[
 C_m(S)\prod_{b\in S}w_b\prod_{b\notin S}(1-w_b),
\]
where all products over \(b\) use the squarefree integers \(b\ge2\).
The last product lies between a fixed positive constant
\(c_*=\prod_b(1-w_b)\) and \(1\). These exact-presence events are
disjoint as \(S\) varies. Bounding the event of at least \(m\)
coprime occupants above by the union of all \(m\)-kernel tests,
and below by the disjoint exact-presence events, gives
\[
 c_*\sum_{\substack{|S|=m\\S\ {\rm pairwise\ coprime}}}
            C_m(S)\prod_{b\in S} w_b
 \ \le\ d_m^{\mathrm{cop}}\
 \le\ \sum_{\substack{|S|=m\\S\ {\rm pairwise\ coprime}}}
            C_m(S)\prod_{b\in S} w_b.                          \tag{942.D3}
\]
Both union bounds can first be applied to finite sets of kernels.
The large-kernel count gives the upper bound in the limit. The
lower bound follows by summing over any finite number of the disjoint
events. For ordinary occupancy the same argument gives
\(c_* e_m((w_b)_b)\le d_m\le e_m((w_b)_b)\) with \(e_m\) denoting the
\(m\)-th elementary symmetric sum.

### The contribution of prime kernels

On truncating the product for \(C_m\) at \(p\le m\), the logarithm
of the remaining factors is \(O(m^2\sum_{p>m}p^{-2})=o(m)\). Indeed
each factor is the probability that at most one among \(m\) uniform
residues is zero, with failure probability at most
\(\binom m2/p^2<1/2\) for \(p>m\). The prime number theorem bounds
the tail sum by \(O(1/(m\log m))\). On \(p\le m\), the sum of
\(\log(1+(m-1)/p)\) is \(o(m)\): the primes \(p\le\sqrt m\)
contribute \(O(\sqrt m\log m)\), and the intervals
\(m/2^{j+1}<p\le m/2^j\) above \(\sqrt m\) contribute in total
\(O((m/\log m)\sum_{j\ge0}(j+1)2^{-j})\). Mertens' product
theorem now gives
\[
                  \log C_m=-m\log\log m-\gamma m+o(m).
                                                                    \tag{942.D4}
\]
[Mertens, J. Reine Angew. Math. 78 (1874), 46–62](https://www.degruyterbrill.com/document/doi/10.1515/crll.1874.78.46/html)

Put \(E_m=e_m((p^{-s})_p)\). For the upper bound in (942.D3), choose
one designated prime from each kernel in a coprime \(m\)-set,
for example the largest. These primes are distinct. The selected
primes determine labels for the kernels, and any additional prime
can occur in at most one labelled kernel. Summing over the selected
primes and dropping the ordering restrictions on the remaining
primes yields the upper bound
\[
       C_m E_m \prod_p (1+m p^{-s})= C_m E_m\,\exp(o(m)).
\]
The logarithm of the product is \(O(m^{1/s}\log m)=o(m)\), by
bounding primes by integers below \(m^{1/s}\) and summing
\(m k^{-s}\) over integers above it. For the lower bound, take
only sets \(S\) of distinct prime kernels. The logarithmic cost
of the product correction to \(C_m\) in (942.D2) is \(o(m)\)
uniformly over such sets. The primes \(p\le m\) already have the
preceding \(o(m)\) bound. Those with \(m<p\le m\log m\) have
total cost
\[
 \ll m\sum_{m<p\le m\log m} p^{-1}
   \ll \frac{m\log\log m}{\log m}
\]
by dyadic subdivision and the prime count. At most \(m\) of the
remaining primes can occur, each with cost at most \(1/\log m\).
It follows from (942.D3) that
\[
             \log d_m^{\mathrm{cop}}=\log C_m+\log E_m+o(m).
                                                                    \tag{942.D5}
\]

### Elementary symmetric sums at large order

For either the sequence \(t_1<t_2<\cdots\) of primes or that of
squarefree integers at least two, let \(\rho_j=t_j^{-s}\) and
\(N(Y)=\#\{j:t_j\le Y\}\). Then \(N(Y)\sim Y/\log Y\) or
\(N(Y)\sim\delta Y\), respectively. The second relation follows
from \(\mu^2(t)=\sum_{d^2\mid t}\mu(d)\), summing through \(Y\).
Choose \(Y=Y(m)\) by
\[
 m=\sum_j(1+(t_j/Y)^s)^{-1},\qquad z=Y^s.
\]
The sum is continuous and strictly increasing in \(Y\), tending
to zero at zero and to infinity. As \(Y\to\infty\),
\[
 m\sim A N(Y),\qquad
 \log\prod_j(1+z\rho_j)\sim s A N(Y).
\]
Indeed \(N(Yu)/N(Y)\to u\) for every fixed \(u>0\), which
gives the integral limits for \(1/(1+u^s)\) and
\(\log(1+u^{-s})\) on compact subintervals of \((0,\infty)\).
Above \(u=R\), dyadic subdivision bounds the normalized tail
by \(O(R^{1-s})\) using \(N(Yu)\ll u N(Y)\) for \(u\ge1\).
Below \(u=\eta\), the portion \(t_j\le\sqrt Y\) contributes
\(O(\sqrt Y\log Y)=o(N(Y))\) before normalization. For
\(\sqrt Y<t_j\le\eta Y\), the counting asymptotics give
\(N(Yu)\ll u N(Y)\) for \(u\ge Y^{-1/2}\), so dyadic
subdivision bounds the normalized tail by
\(O(\eta(1+\log(1/\eta)))\). Integration by parts gives
\(\int_0^\infty\log(1+u^{-s})\,du=s A\).

The tilted elementary symmetric coefficients
\[
  \frac{z^k e_k((\rho_j)_j)}{\prod_j(1+z\rho_j)}
\]
are the probabilities for the sum of independent Bernoulli variables
with success probabilities \(z\rho_j/(1+z\rho_j)\). The sum has
mean \(m\) and variance at most \(m\), so one coefficient at
distance \(O(\sqrt m)\) from \(m\) has probability
\(\gg m^{-1/2}\). For these indices the successive probability ratios
lie between
\[
  \frac{z}{k+1}\sum_{j>k}\rho_j
       \quad\text{and}\quad
  \frac{z}{k+1}\sum_{j\ge1}\rho_j .
\]
This follows by adding one unused index to each \(k\)-set, with every
\((k+1)\)-set counted \(k+1\) times. The lower bound is
\(\gg1\), using \(2m\le j\le3m\) where \(t_j\asymp Y\).
The upper bound is polynomial in \(m\), since
\(Y=O(m\log m)\). Moving at most \(O(\sqrt m)\) indices
therefore costs at most \(\exp(O(\sqrt m\log m))\). The probability
at \(m\) lies between \(\exp(-o(m))\) and \(1\), proving
\[
                    \log e_m((\rho_j)_j)=s m-s m\log Y+o(m).
\]
For primes \(Y\sim m\log m/A\), while for squarefree integers
\(Y\sim m/(\delta A)\). Substituting these two values and (942.D4)
into (942.D5) and the ordinary occupancy bounds proves (942.D1).

The argument gives an asymptotic for densities and leaves the
location of the first intervals with growing coprime occupancy,
and the uniform logarithmic upper bound, open.
