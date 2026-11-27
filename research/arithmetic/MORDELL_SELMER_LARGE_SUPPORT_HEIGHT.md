# Rational-point height and large cubic ray support

Use the family and rational isogeny \(\psi:E\to E'\) from
[the inert-prime Mordell Selmer calculation](MORDELL_INERT_PRIME_SELMER_MOMENTS.md).
Thus \(R\) is squarefree with all prime factors inert in
\(\mathbf Q(\sqrt{-2})\) and \(2\bmod3\), and
\(E':y^2=x^3+54R^2\). The tame support of a class in the
forward Selmer group is the product of primes \(p\mid R\)
where its cubic ray character is relatively ramified.
For a finite rational point on \(E'\), write
\(x=a/c^2,\ y=b/c^3\) with \(c>0\) and
\(\gcd(a,c)=\gcd(b,c)=1\).
Use the multiplicative height \(H_x=\max\{|a|,c^2\}\)
in this specified model.

**Theorem MH1.** For a rational point on \(E'\) the tame support
of its forward Kummer class is \(m=\gcd(|a|,R)\). For all
\(X,H,U\ge1\), the number of pairs \((R,P)\), with
\(X\le R<2X\), \(P\in E'(\mathbf Q)\) a finite point,
\(H_x(P)\le H\), and \(R/m\le U\), is
\[
 \ll_\epsilon (XHU)^\epsilon
       \left( X^{2/3}
            + \frac{H^{5/3} U^{22/9}}{X^{4/3}}
            + \frac{H^{11/8} U^2}{X} \right).
 \tag{1}
\]
In particular for fixed \(\eta>0\), apart from
\(O_\eta(X^{1-\sigma})\) curves in the dyadic interval for
some \(\sigma>0\), a fully supported Selmer line has no
rational representative of \(x\)-height at most \(X^{7/5-\eta}\).
The same holds for support \(m\ge R/X^\kappa\) and
\(H\le X^{7/5-22\kappa/15-\eta}\) for fixed \(\kappa,\eta\ge0\)
with \(\eta>0\).

For the support claim, the dual Kummer description for
the isogeny from \(E'\) back to \(E\) evaluates
\(y-3R\sqrt6\) on \(E'\) for the forward class into
\(\operatorname{Sel}_{\psi}(E/\mathbf Q)\).
[Bhargava–Elkies–Shnidman, Remark 25](https://math.huji.ac.il/~shnidman/j%3D0fixedtypos.pdf)
At \(p\mid R\) the quadratic field \(\mathbf Q(\sqrt6)\)
splits, and its two local evaluations are represented by
\(b\pm3R c^3\sqrt6\), since the denominator is a cube.
If \(p\nmid a\), the point equation shows \(b\) is a
unit at \(p\), so both evaluations are units. Every
unit of \(\mathbf Q_p\) is a cube for \(p\equiv2\bmod3\).
If \(p\mid a\), then \(c\) is a unit, \(v_p(b)=1\), and
the difference of the two evaluations has valuation one.
Their product is \(a^3\), so their valuations are one and
\(3v_p(a)-1\), in some order. This is a nonzero ramified
local class. In the character identification of Theorem
MS1 this is exactly a nonzero ray coordinate at \(p\).
There are no finite points with \(a=0\) since \(54R^2\)
is nonsquare.

To count the points, put \(S=R/m\), \(a=mt\), and \(b=mu\).
The last division is valid from the point equation.
It yields
\[
                    u^2-54 S^2 c^6=m t^3,\qquad
 \gcd(c,t)=\gcd(S,t)=1,\quad 0<|t|\le HS/X,\quad c\le\sqrt H.
 \tag{2}
\]
The relation \(\gcd(S,t)=1\) follows from the squarefree
support definition. Fix \(S\) and \(T=|t|\), and temporarily
discard the prime conditions once \(X/S\le m<2X/S\) is imposed.
The sign of \(t\) and the sign of \(u\) cause only bounded
multiplicity. For each residue \(r\bmod T^3\) with
\(r^2\equiv54S^2\bmod T^3\) one has \(u\equiv r c^3\bmod T^3\).
There are \(O_\epsilon(T^\epsilon)\) choices for \(r\):
at primes away from \(6S\) there are at most two roots
per prime power, and \(\gcd(S,T)=1\). Prime \(2\) cannot
divide \(T\) because the valuation of \(54S^2\) at \(2\)
would be one. The exponent of \(3\) in \(T\) is at most
one since the valuation of \(54S^2\) at \(3\) is three.
For \(v_3(T)=1\), there are three roots modulo \(27\).

Write \(\alpha=\sqrt{54}\) and take \(u\ge0\). Equation (2)
bounds \(|u-\alpha S c^3|\) by a constant times each of
\((X/S)^{1/2}T^{3/2}\) and \(X T^3/(S^2 c^3)\).
For a dyadic \(c\)-interval of length at most \(M\) with
\(M\le c<2M\), it places \((u-r c^3)/T^3\) among integers
within distance
\[
 \Delta_M\ll\min\{(X/S)^{1/2}T^{-3/2},\ X/(S^2 M^3)\}
 \tag{3}
\]
of \((\alpha S-r)c^3/T^3\). If \(\Delta_M\ge1/4\), this
costs \(O(M\Delta_M)\). Otherwise it is bounded by the
number of fractional parts within \(\Delta_M\) of zero.
The dyadic sum of the main terms \(M\Delta_M\) is
\(O(X^{2/3}/(ST))\), by splitting at
\(M=X^{1/6}S^{-1/2}T^{1/2}\).

To bound the remaining discrepancy, for integral \(j\ge1\)
and all rational \(a/q\) the elementary quadratic-irrational
inequality gives
\[
  \left|j(\alpha S-r)/T^3-a/q\right|
                         \gg (j S T^3 q^2)^{-1}.
\]
Indeed apply \(|\alpha-b/v|\gg v^{-2}\) with
\(b/v=(r+a T^3/(jq))/S\); the bound for \(\alpha\)
follows by multiplying by \(\alpha+b/v\) when the difference
is small. Dirichlet approximation with cutoff
\(Q\asymp M^{3/2}(j S T^3)^{1/2}\) gives a coprime \(a,q\)
with \(q\le Q\) and \(q\gg Q/(jS T^3)\). Weyl's inequality
for the cubic exponential sum over the \(c\)-interval
then gives
\[
 \left|\sum_c e\!\left(j(\alpha S-r)c^3/T^3\right)\right|
 \ll_\epsilon M^{1+\epsilon}
         \left((jS T^3/M^3)^{1/8}+M^{-1/4}\right).
\]
The shifted interval is allowed because translating the
summation variable preserves the leading coefficient.
If its length is below \(M\), the same bound with \(M\)
follows from Weyl's inequality with the actual length,
since each of the three terms multiplied by the length
to the first power is increasing in that length.
[Kerr, equation (1) for the classical Weyl inequality](https://arxiv.org/pdf/2101.08058)
Use the
[Erdős–Turán inequality in the form of the opening discrepancy bound](https://londmathsoc.onlinelibrary.wiley.com/doi/10.1112/jlms.70519)
with cutoff \(\max(1,\lfloor(M/(S^{1/3}T))^{1/3}\rfloor)\)
when \(M\ge S^{1/3}T\). It gives discrepancy
\[
              \ll_\epsilon M^\epsilon
                   (M^{2/3}S^{1/9}T^{1/3}+M^{3/4}).
\]
The same bound holds trivially for \(M<S^{1/3}T\).
This estimate is uniform in the root \(r\).

After the dyadic sum \(M\le\sqrt H\), the cost per \(S,T\),
up to the divisor and height subpowers, is
\[
           \frac{X^{2/3}}{S T}
           +H^{1/3}S^{1/9}T^{1/3}+H^{3/8}.
\]
Each \(S,t,c,u\) determines \(m\) by (2) and therefore
determines \(R\) and the point. Summation over
\(T\le HS/X\), then \(S\le U\), gives (1).
The principal range \(U=1\) has error exponents below one
as soon as \(H\le X^{7/5-\eta}\). For
\(U=X^\kappa\) the same comparison gives the asserted
general range. The dyadic family has order
\(X(\log X)^{-3/4}\) by the unweighted squarefree prime
support estimate in Corollary MS3, before enforcing any
condition modulo \(9\). In particular the exceptional
curves have a power-smaller count. Nonzero Selmer lines
with the full tame support exist for every \(R>1\), by
Corollary MS2. The calculation bounds low-height
representatives and leaves open rational solubility above
that height.

An exact support control is \(R=115\), with the integral
points \((-5,845),(-69,621),(115,1495)\) on \(E'\).
Their support products are \(5,23,115\), respectively.
At \(R=5\), the point \((-5,35)\) gives full support at
height \(5\), consistent with the exceptional curves in
the asymptotic bound.

## Distribution of the support among all Selmer lines

**Corollary MH2.** Give equal weight to each pair \((R,L)\)
with \(N(E)\le Y\) and \(L\) a nonzero
\(3\)-Selmer line of \(E\). If \(m\) is its tame support,
then as \(Y\to\infty\) the variable
\(\log(R/m)/\log R\) converges in distribution to
\(\operatorname{Beta}(1/4,1/2)\); the finite case \(R=1\)
may be omitted. In particular, for any \(0<\theta<7/10\)
and \(0<\kappa<(21-30\theta)/22\), at least a limiting
lower proportion
\[
 \frac{\displaystyle\int_0^\kappa v^{-3/4}(1-v)^{-1/2}\,dv}
      {\displaystyle\int_0^1 v^{-3/4}(1-v)^{-1/2}\,dv}
 \tag{4}
\]
of the lines have no rational representative on \(E'\)
of \(x\)-height at most \(N(E)^\theta\).
This is a statement about the simultaneous conductor
sum, with the lines weighted equally across all curves.

For a divisor \(m\mid R\) with \(m>1\), the number of
nonzero lines with tame support exactly \(m\) is
\(2^{\omega(m)}3^{\epsilon_R}/2\) by Corollary MS2.
The lines with \(m=1\) contribute at most one per curve,
which is negligible compared to the total number of lines
by Corollary MS3 at \(\nu=0,1\). We therefore analyze the
factorization \(R=mS\) with the weight \(2^{\omega(m)}\),
keeping the congruence classes of \(R\bmod9\) up to sign.
This is a weighted divisor distribution. The classical
comparison of Bareikis and Manstavičius treats a
conditional divisor probability on each integer. Here
each integer instead contributes its actual number of
Selmer lines, and only the prescribed prime set is allowed.
[Bareikis–Manstavičius, Theorem 1](https://web.vu.lt/mif/e.manstavicius/wp-content/uploads/2016/08/EM-5.pdf)

Write \(f_a(n)\) for \(a^{\omega(n)}\) on the
squarefree integers supported on \(\mathcal P\), and
zero elsewhere, for \(a=1,2\). For a character \(\chi\)
of \((\mathbf Z/9\mathbf Z)^\times/\{\pm1\}\) the relevant
two-variable Euler factors are
\[
 1+2\chi(p)p^{-s}+\chi(p)p^{-t}
  =(1+2\chi(p)p^{-s})(1+\chi(p)p^{-t})
   \left(1-\frac{2\chi(p)^2p^{-s-t}}
                    {(1+2\chi(p)p^{-s})(1+\chi(p)p^{-t})}\right).
 \tag{5}
\]
The correction product has coefficients \(c_\chi(d,e)\)
with \(\sum_{d,e}|c_\chi(d,e)|(de)^{-3/4}<\infty\):
at \(p\in\mathcal P\) the geometric series in both
denominators converge absolutely and the local norm of
the correction is \(O(p^{-3/2})\).
Thus the weighted sum for \(mS\le X,\ S\le X^\kappa\)
is the sum, against \(c_\chi(d,e)\), of
\[
 \sum_{\substack{uv\le X/(de)\\v\le X^\kappa/e}}
                 f_2(u)\chi(u) f_1(v)\chi(v).
 \tag{6}
\]
It suffices to restrict \(d,e\le(\log X)^K\) for sufficiently
large fixed \(K\). Indeed (6) has absolute bound
\(O(X(de)^{-1}(1+\log X)^{1/4})\) by the partial sum bound
\(O(x)\) for \(f_2\) and the harmonic sum bound
\(O((1+\log x)^{1/4})\) for \(f_1\). Absolute convergence
with exponent \(3/4\) then gives the asserted tail control.

For the principal character let
\[
 C_a=G_a/\Gamma(a/4)\qquad(a=1,2),
\]
using the constants in Corollary MS3 with \(z=a\).
The one-variable Selberg–Delange estimate there gives
\(\sum_{n\le x} f_a(n)\sim C_a x(\log x)^{a/4-1}\).
For \(d,e\) in the logarithmic range, sum over \(u\)
in (6) first; its upper endpoint is at least
\(X^{1-\kappa}/d\), so the one-variable asymptotic is
uniform. The remaining sum uses the harmonic measure
of \(f_1\), for which partial summation gives
\[
 (\log X)^{-1/4}
       \sum_{v\le X^r}\frac{f_1(v)}v\longrightarrow 4C_1 r^{1/4}
       \qquad(0\le r\le1).
\]
The limit is continuous. Since \(\kappa<1\), the
factor \((1-r)^{-1/2}\) from the \(u\)-sum is bounded
on a neighborhood of \([0,\kappa]\).
It follows that (6), uniformly in the logarithmic
range, is asymptotic to
\[
 \frac{C_1 C_2 X}{de}(\log X)^{-1/4}
                          \int_0^\kappa r^{-3/4}(1-r)^{-1/2}\,dr .
\]
The constant from the \(c_\chi\)-sum is the principal
correction product in (5) at \(s=t=1\).
For a nonprincipal \(\chi\), the prime mean of
\(f_2\chi\) is zero, so the \(u\)-sum has arbitrary
logarithmic saving by the same theorem used in MS3.
With the harmonic \(v\)-sum and the correction tail,
these terms contribute \(o(X(\log X)^{-1/4})\).
Character orthogonality gives a factor \(1/3\) for
each class modulo \(9\) up to sign.

The total unconstrained tame weight on \(R\) is
\(3^{\omega(R)}\). Its one-variable Euler product is
the product of (5) at \(s=t\); its constant from
Selberg–Delange gives the same normalization as
replacing the integral endpoint by one, by the beta
integral and \(G_3=G_1G_2\sum_{d,e}c_1(d,e)/(de)\).
This proves the CDF ratio in each of the three
congruence classes. The contribution of \(R\le
X^{1-\delta}\) for fixed \(\delta>0\) is negligible
by that total asymptotic. Consequently using
\(\log R\) in the denominator has the same continuous
limiting distribution as using \(\log X\).
Finally use \(X=(Y/N_j)^{1/2}\) separately in the
three classes from (4) of Corollary MS3 and multiply
by \(3^{\epsilon_R}/2\). The limiting ratio is the
same in every class, so remains the same after
adding their contributions.

For (4), apply Theorem MH1 with \(U\ll X^\kappa\)
and \(H\ll X^{2\theta}\) in dyadic intervals
\(X\le R<2X\). The number of representatives in
these ranges is power-smaller than \(X\) by the
strict inequality on \(\kappa\). Each rational
representative accounts for at most one Selmer
line. Summing over the dyadic intervals up to
the conductor bound therefore removes a
negligible fraction of the pairs counted by the
limit law. This establishes (4) without inferring
global insolubility of any of these lines.
