# A power-saving exceptional set for balanced sums in \(C_{12}\)

Let \(\mathcal P=\{p:p\equiv1\bmod 12\}\), let \(C_{12}\) be the positive
integers whose prime factors all belong to \(\mathcal P\), and include \(1\)
in \(C_{12}\). Write \(f=\mathbf1_{C_{12}}\) and count ordered pairs by
\[
 r_{12,H}(N)=\sum_{\substack{a+b=N\\ |a-b|\le H}}f(a)f(b),
 \qquad a,b\ge1 .
\]

**Theorem D12Power.** There is an absolute \(c>0\) such that, for
\(H=X^{0.99}\) and sufficiently large \(X\),
\[
 \#\left\{X<N\le2X:\ N\equiv2\bmod12,\quad
    r_{12,H}(N)<c\,\frac{H}{(\log X)^{3/2}}\right\}
      \ll X^{1-10^{-6}}. \tag{D12P}
\]
The congruence on \(N\) is necessary. The conclusion is a lower bound for
the **unweighted** count at its semigroup scale. The full-count variance
and its arithmetic main term are treated separately.
[Power-window variance theorem](PRESCRIBED_PRIME_POWER_WINDOW_VARIANCE.md).

The proof counts \(a=mk,\ b=m'l\), with small \(m,m'\in C_{12}\) and
\(k,l\in C_{12}\) free of primes below a fixed power of \(X\). On those
remaining factors the auxiliary weight \(w=f\tau_4\) is bounded. Its
four-character factorization supplies a power saving uniformly in the
divisibility conditions used by the sieve. Each small cofactor is uniquely
recoverable from its summand.

| Parameter | Value or choice | Use |
|---|---|---|
| \(H\) | \(X^{0.99}\) | Length of the summand interval |
| \(\rho\) | \(10^{-4}\) | Rational denominators up to \(X^\rho\) in the main term |
| \(\kappa\) | \(\rho/10^4\) | Cofactor and divisibility parameters at most \(X^\kappa\) |
| \(\delta\) | \(\rho/100=10^{-6}\) | Exception exponent retained in (D12P) |
| \(s\) | Fixed sufficiently large for the sieve bound in Section 4 | Level exponent |
| \(r,z\) | \(r=\kappa/(2s),\quad z=X^r\) | Prime cutoff, with \(z^s\le X^\kappa\) |

The proof dependencies are
\[
 \text{weighted progression means}\ \longrightarrow\
 \text{weighted binary correlations}\ \longrightarrow\
 \text{lower-bound sieve on }kl\ \longrightarrow\
 \text{unweighted count}.
\]
In what follows \(X^{o(1)}\) losses are uniform in the bounded parameters.
Fixed powers of logarithms and of divisor functions at arguments
\(O(X^{O(1)})\) are absorbed into these losses.

## 1. Weighted progression means and their local coefficients

Define \(\tau_4(n)\) as the number of ordered four-factor decompositions of
\(n\), \(w(n)=f(n)\tau_4(n)\), and \(w_B(n)=w(Bn)\). Write \(\mathcal X\)
for the four characters modulo \(12\), taking every character to be zero
at integers not coprime to \(12\). For \(B\in C_{12}\),
\[
 \sum_{n\ge1}\frac{w_B(n)}{w(B)n^t}
     =E_B(t)\prod_{\psi\in\mathcal X}L(t,\psi),\qquad \Re t>1 .
 \tag{1}
\]
At \(p\equiv5,7,11\bmod12\), the four character values are two \(+1\)'s
and two \(-1\)'s. Thus the correction factor in \(E_B\) is
\((1-p^{-2t})^2\). The factors at \(2,3\) and allowed \(p\nmid B\) are
one. For \(p^v\parallel B\), the correction with \(u=p^{-t}\) is
\[
 (1-u)^4\sum_{i\ge0}
       \frac{\binom{v+i+3}{3}}{\binom{v+3}{3}}u^i .
\]
The coefficient ratio is a cubic polynomial in \(i\). Multiplication by
\((1-u)^4\) leaves a polynomial of degree at most three with constant one
and uniformly bounded coefficients, as seen from \(i=0,1,2,3\).
Consequently \(E_B(t)=\sum_j e_B(j)j^{-t}\) satisfies
\[
 \sum_j |e_B(j)|j^{-0.6}\ll_\varepsilon B^\varepsilon
 \quad(\varepsilon>0).
\]
Indeed the product away from \(B\) converges absolutely, and the remaining
factors are \(1+O(p^{-0.6})\) in this norm. The numbers \(|e_B(j)|\) are
bounded by a fixed divisor-function power uniformly in \(B\), and
\(w(B)\ll_\varepsilon B^\varepsilon\). For \(B\notin C_{12}\) one has
\(w_B=0\).

Let \(\chi\bmod D\) have primitive conductor not dividing \(12\). Each
\(\chi\psi\) is nonprincipal with period dividing \(12D\), so its sum on
any interval is \(O(D)\). In the cumulative four-character convolution to
\(Y\), split the four positive indices into dyadic intervals. On each box,
sum last over the index with largest dyadic lower endpoint. Its upper bound
depends on the other three indices, whose dyadic lengths have product at
most \(Y^{3/4}\). Indeed the product of all four lower endpoints is at
most \(Y\) and the largest is at least its fourth root. Summing over
\(O(\log^4(2Y))\) boxes and convolving with \(e_B(j)\chi(j)\) gives
\[
 \left|\sum_{n\le Y} w_B(n)\chi(n)\right|
 \ll_\varepsilon D (2BY)^\varepsilon Y^{3/4}. \tag{2}
\]

For the untwisted mean, absorb the two Euler factors of the principal
character at \(2,3\) into \(E_B\), pull out \(\zeta\), and call the remaining
coefficient sequence \(a_B\). The three nonprincipal characters give
\[
 A_B(v):=\sum_{n\le v}a_B(n)
 \ll_\varepsilon B^\varepsilon v^{2/3}\log^3(2v)
\]
by the three-variable version of the same box argument.
Also \(\sum_{n\le v}|a_B(n)|\ll_\varepsilon B^\varepsilon v\log^3(2v)\) by
convolution with \(\tau_3\). Split \(a_B*1\), a Dirichlet convolution, at
the \(a_B\)-index \(V=Y^{3/4}\). For indices up to \(V\), replacing the
floor by its argument costs \(O_\varepsilon(B^\varepsilon V\log^3(2Y))\).
The tail of \(Y\sum a_B(v)/v\) is
\(O_\varepsilon(B^\varepsilon YV^{-1/3}\log^3(2Y))\) by partial summation.
The original terms with index above \(V\) sum to
\(\sum_{j\le Y/V}(A_B(Y/j)-A_B(V))\), with the same bound.
Restoring \(w(B)\) and evaluating the residue at one yields
\[
 \begin{split}
 \sum_{n\le Y}w_B(n)&=C_B Y+O_\varepsilon((2BY)^\varepsilon Y^{3/4}),\\
 C_B&=C\prod_{p^v\parallel B}
 \frac{\sum_{i\ge0} w(p^{v+i})p^{-i}}
      {\sum_{i\ge0} w(p^i)p^{-i}},\qquad
 C=\mathop{\rm Res}_{t=1}\sum_{n\ge1} w(n)n^{-t}>0 .
 \end{split} \tag{3}
\]
For forbidden \(B\) use \(C_B=0\). Positivity of \(C\) also follows directly
from (1) at \(B=1\), using the positive residue of the principal character
and the positive values at one of the three fixed real nonprincipal
\(L\)-functions.

For allowed \(m,b\le X^\kappa\) let
\[
 A_{m,b}(n)=\mathbf1_{mb\mid n}\,w(n/m)
\]
with zero value when the divisibility fails. Thus \(n=mk\) and \(b\mid k\).
Write \(e(t)=\exp(2\pi i t)\). For an interval \(J\) of length
\(\asymp X^{0.97}\) with left endpoint \(\asymp X\), and reduced \(a/D\)
with \(D\le X^\rho\),
\[
 \sum_{n\in J} A_{m,b}(n)e(an/D)
       =|J| B_{m,b}(a/D)+O(X^{0.80}). \tag{4}
\]
To see this and identify the coefficient, put \(n=mk\) and split \(k\)
into classes modulo \(M=\operatorname{lcm}(b,D)\le X^{\rho+\kappa}\).
Retain the classes \(c\) with \(b\mid c\). In a class \(c\bmod M\), put
\(t=(c,M)\), \(M'=M/t\) and \(c'=c/t\). Character orthogonality for the
unit class of \(k/t\bmod M'\) uses (2) unless the primitive conductor divides
\(12\). The latter characters sum at \(c'\) to
\(\varphi((M',12))\mathbf1_{c'\equiv1\bmod(M',12)}\), by CRT and
orthogonality on the image of the units \(1\bmod12\). On the support of
\(w(tn')\), their twists are \(w(tn')\mathbf1_{(n',M')=1}\) whenever
\(t\in C_{12}\). Inclusion-exclusion over \(j\mid\operatorname{rad}(M')\)
gives (3) with \(B=tj\), with zero mass for forbidden \(tj\).
The resulting class mass per unit length in \(k\) is
\[
 \frac{\varphi((M',12))\mathbf1_{c'\equiv1\bmod(M',12)}}
      {t\varphi(M')}\sum_{j\mid\operatorname{rad}(M')}
             \frac{\mu_{\rm M}(j) C_{tj}}j ,
\]
where \(\mu_{\rm M}\) is the Möbius function. Each main term is exactly
linear in interval length. Even bounding the number of classes, characters
and divisors together by \(O(M^3)\) and each cumulative error by
\(O(MX^{3/4+o(1)})\) leaves \(O(X^{0.80})\). Scaling \(J\) by \(1/m\)
accounts for the length in (4).

The local version of this computation uses normalized Haar measure on
\(\mathbb Z_p\). For \(p\nmid12\) and \(p^j\parallel12\), respectively, set
\[
 W_p(k)=\frac{w(p^{v_p(k)})}
             {(1-1/p)\sum_{i\ge0}w(p^i)p^{-i}},\qquad
 W_p(k)=p^j\mathbf1_{k\equiv1\bmod p^j}.
\]
Each has integral one. Values at the single point \(k=0\) are immaterial.
Factor the additive character \(k\mapsto e(amk/D)\) into its local
characters \(e_p(amk/D)\) by CRT, taking the trivial character at primes
away from \(D\). Then
\[
 B_{m,b}(a/D)=\frac Cm\prod_p\int_{\mathbb Z_p}
   \mathbf1_{v_p(k)\ge v_p(b)} W_p(k)e_p(amk/D)\,dk . \tag{5}
\]
Indeed, by (3), requiring \(p^v\mid k\) multiplies the normalized mass by
\(p^{-v}\sum_i w(p^{v+i})p^{-i}/\sum_i w(p^i)p^{-i}\), the integral of
\(W_p\) on that ball for \(p\nmid12\). Subtract successive tails to obtain
an exact valuation. The surviving unit classes are equidistributed by the
character calculation above, with the sole congruence restriction at the
primes of \(12\). Every surviving valuation divisor there scales by a unit
\(1\bmod12\). Multiplying the local masses gives the displayed class formula,
and testing it against \(e(amk/D)\) proves (5).

At \(p\nmid12\) these integrals depend only on the local conductor, since
both weights are radial. Put \(u=v_p(m)\), \(v=v_p(b)\), and \(t=v_p(D)\).
For \(t\le u+v\) their magnitude is at most one. For \(t>u+v\) it is
\(O((1+t)^3/p^{t-u})\). In fact at conductor \(p^{t-u}\) the shell of
valuation \(t-u-1\) contributes \(-p^{-(t-u)}W_p(p^{t-u-1})\), the lower
shells contribute zero, and the tail from valuation \(t-u\) has the asserted
bound by summing a cubic against the geometric Haar masses. At \(p\mid12\)
only conductors dividing \(12\) occur, with a fixed congruence phase.
The untested divisibility factors at primes away from \(D\) have mass at
most one. Multiplying and absorbing fixed divisor-function powers gives
\[
 |B_{m,b}(a/D)|\ll (mb)D^{-1+o(1)}. \tag{6}
\]
The elementary bound \(|B_{m,b}|\le C/m\) from (5) will also be used.


## 2. Fourier estimates and uniform weighted correlations

Fix \(I=(x,x+H]\) with \(X/3\le x\le X\), and define
\[
 F_{m,b}(\beta)=\sum_{n\in I}A_{m,b}(n)e(\beta n),\qquad
 J_I(v)=\int_I e(vt)\,dt,\qquad T_D=\lfloor H^2/X^{1+\rho}\rfloor .
\]
For each \(\beta\) on the unit circle choose a reduced Dirichlet approximant
\(a/D\) with \(D\le T_D\), \(|v|=|\beta-a/D|\le1/(DT_D)\), breaking ties
by a fixed order and using the signed circular difference. When \(D<X^\rho\),
partition \(I\) into equal intervals of length asymptotic to \(X^{0.97}\).
On each piece \(e(vn)\) varies by \(O(X^{-0.01+\rho})\). Freezing that
factor, applying (4) and restoring the continuous phase therefore gives
\[
 F_{m,b}(\beta)=B_{m,b}(a/D)J_I(v)+O(HX^{-0.008}). \tag{7}
\]
The phase error is \(O(HX^{-0.01+\rho+o(1)})\) since
\(A_{m,b}\le X^{o(1)}\). The cumulative error from (4) is
\(O(HX^{-0.17})\).

For \(X^\rho\le D\le T_D\) we prove
\[
 |F_{m,b}(\beta)|\ll H X^{-\rho/8+o(1)}. \tag{8}
\]
Expand \(w(bn')\) using (1), so the original index is
\(n=mbj n_1n_2n_3n_4\) with coefficient
\(w(b)e_b(j)\prod_{i=1}^4\psi_i(n_i)\). Put \(K=X^{0.01}\).
For a fixed \(j\) there are at most \(H/(mbj)+1\) possible residual
products, each with at most \(X^{o(1)}\) decompositions. Hence the tail is
\[
 \begin{split}
 \left|\sum_{j>K}\cdots\right|
 &\ll X^{o(1)}
 \left(\frac{H}{mb}\sum_{j>K}\frac{|e_b(j)|}j+
                         \sum_{j\le3X}|e_b(j)|\right)\\
 &\ll X^{o(1)}(HK^{-0.4}+X^{0.6})
 \end{split}
\]
by the \(0.6\)-norm of \(e_b\). This bound holds on each individual short
interval, including its endpoint error.

For \(j\le K\), split the five variables into \(O(\log^5 X)\) dyadic boxes.
In each contributing box the product of the four character-variable lower
endpoints is \(\gg X/(32mbK)\), so its largest endpoint is \(\gg X^{0.24}\).
If it exceeds \(X^{0.7}\), combine the other indices, including \(mbj\),
into \(\ell\ll X^{0.30}\). The longest variable carries a single fixed
character. All combined coefficients are \(\ll X^{o(1)}\), by the uniform
divisor bound on \(e_b\) and the number of divisors. Split the longest
index by its class modulo \(12\). Its allowable values in the dyadic box
and the product strip form an interval in that class. The geometric sum
has frequency \(12\ell\beta\). The reduced approximant to \(12\beta\) has
denominator \(D'=D/(12,D)\asymp D\) and error \(O(1/(D')^2)\).

On \(U\le\ell<2U\) with \(U<X^{\rho/2}\) we have
\(\|12\ell\beta\|\gg1/D'\), where \(\|\cdot\|\) denotes distance to an integer.
The contribution relative to \(H\) is
\(O(X^{o(1)}UD'/H)\ll X^{-0.01-\rho/2+o(1)}\).
For \(U\ge X^{\rho/2}\), Cauchy–Schwarz and spacing give a relative bound
\[
 X^{o(1)}(1/D'+1/U+U/H+D'/H)^{1/2}
       \ll X^{-\rho/4+o(1)} .
\]
Indeed a block of \(cD'\) consecutive indices, for a sufficiently small
absolute \(c>0\), has spacing \(\gg1/D'\) modulo one. Summing the squared
geometric bound on these blocks yields
\[
 \sum_{U\le\ell<2U}\min(V,\|12\ell\beta\|^{-1})^2
       \ll(U/D'+1)(V^2+D'V),\qquad V\asymp H/U .
\]
The \(V^2\) covers the nearest point of each block, and the other terms
are bounded by summing \(\min(V,O(D'/i))^2\) over \(i\ge1\).

If the largest character-variable endpoint is at most \(X^{0.7}\), isolate
that variable and combine all the others. Both resulting indices have
size \(\gg X^{0.23}\) and \(\ll X^{0.77}\). Subdivide the combined index
dyadically. Write the sum in the larger index as \(\sum_k b_k T_k\) and
use \(|\sum_k b_k T_k|^2\le(\sum_k|b_k|^2)(\sum_k|T_k|^2)\).
This removes the larger coefficient even when it is the isolated character.
Equivalently, if its unit-class indicator is retained in the second factor,
split \(k\) into classes modulo \(12\). The geometric phases then use
\(12\beta\) with denominator \(D'=D/(12,D)\asymp D\) and spacing on blocks
of length \(cD'\), as in the linear bound above. Both calculations give
the same estimates up to constants.
For smaller dyadic size \(M\), expansion leaves differences
\(|j'|\ll HM/X\) of two smaller indices. For each such pair the common range of the larger index is an
interval of length \(O(H/M)\). For each \(j'\) the absolute sum of the
coefficient products on the smaller side is \(O(MX^{o(1)})\).
For \(j'\ne0\), grouping differences into blocks of length at most \(D/2\)
gives
\[
 \sum_{1\le j'\le Y}\min(V,\|j'\beta\|^{-1})
       \ll (Y/D+1)(V+D\log(2D)).
\]
Here \(|\beta-a/D|\le D^{-2}\) gives spacing at least \(1/(2D)\) in a block.
Using \(Y\asymp HM/X\), \(V\asymp H/M\), and adding the diagonal gives
\[
 \frac{|\text{box sum}|^2}{H^2}
 \ll X^{o(1)}
 \left(\frac X{HM}+\frac1D+\frac{M\log X}{H}
                         +\frac{XD\log X}{H^2}\right)
 \ll X^{-\rho+o(1)} .
\]
The original interval \(I\) supplies the only joint product restriction,
so the geometric inner sums retain their interval ranges. The box totals
and the tail bound prove (8).

Define \(P_{m,b}(\beta)=B_{m,b}(a/D)J_I(v)\) if the assigned \(D<X^\rho\),
and zero otherwise. Equations (7)–(8) show
\[
 \|F_{m,b}-P_{m,b}\|_\infty\ll HX^{-\rho/8+o(1)}.
\]
Parseval gives \(\int|F_{m,b}|^2\ll HX^{o(1)}\). For the other transform,
each frequency is assigned to exactly one approximant. By (6) and
\(\int_\mathbb R|J_I(v)|^2\,dv=H\),
\[
 \int |P_{m,b}|^2
 \ll H (mb)^2\sum_{D\le X^\rho}D^{-1+o(1)}
 \ll H X^{4\kappa+o(1)} ,
\]
where integrals without a specified domain are over the unit circle.
Thus for two parameter choices, denoted by subscripts \(1,2\),
\[
 \int|F_1F_2-P_1P_2|^2
 \ll H^3 X^{-\rho/4+4\kappa+o(1)} . \tag{9}
\]

For each reduced \(a/D\) with \(D<X^\rho\), its assigned region contains
all \(|v|<1/(2DT_D)\). A different approximant of denominator \(D'\le T_D\)
would imply \(1\le D'/(2T_D)+D/T_D<1\). Hence extending the real-kernel
integral for this region to all \(v\in\mathbb R\) costs \(O(DT_D)\) by
\(|J_I(v)|\ll\min(H,1/|v|)\). The at most \(X^{2\rho}\) modes together
cost \(O(HX^{-0.01+2\rho})\) since \(X^\rho T_D/H\le H/X\) and (5)
bounds both coefficients. Real Fourier inversion gives
\[
 \int_\mathbb R J_I(v)^2 e(-Nv)\,dv
    =H_N,\qquad H_N=(H-|N-(2x+H)|)_+ .
\]
The coefficient at \(N\) of the model product, using the inverse transform,
is therefore
\[
 H_N\sum_{D<X^\rho}\sum_{\substack{a\bmod D\\(a,D)=1}}
 B_{m,b}(a/D)B_{m',b'}(a/D)e(-aN/D)+O(HX^{-0.01+2\rho}). \tag{10}
\]

Write \(\mathfrak S_{m,b;m',b'}(N)\) for the denominator-grouped series in
(10) extended to all \(D\). Nonzero modes have \(D=D_0D_1\) with
\(D_0\mid12\), \((D_1,12)=1\). The local factors at \(D_1\) in (5)
do not depend on the unit numerator, while the \(D_0\) phase depends only
on the CRT component there. Summing primitive numerators therefore gives
the Ramanujan sum \(c_{D_1}(N)\) on that component and at most twelve
terms on the other. Since \(|c_{D_1}(N)|\le(D_1,N)\), the absolute value
of each denominator group is
\[
 \ll mbm'b'\,D^{-2+o(1)}(D_1,N).
\]
For \(D>R\) extract \(R^{-1/2}\), use
\((D_1,N)\le\sum_{v\mid(D_1,N)}v\), and sum \(D^{-3/2+\varepsilon}\)
over multiples of \(v\). For fixed \(\varepsilon<1/2\) the resulting
divisor sum is \(O(\tau(N))\). In particular, for \(R=X^\rho,\ N\asymp X\)
the tail is \(O(X^{4\kappa-\rho/2+o(1)})\), uniformly in all parameters.

Let \(\star\) denote additive convolution on integers. There are \(O(H)\)
targets in the central interval \(|N-(2x+H)|\le H/2\). Equation (9),
Parseval, (10) and this tail yield
\[
 \sum_{\substack{X<N\le2X\\|N-(2x+H)|\le H/2}}
 \left|(A_{m,b}\mathbf1_I\star A_{m',b'}\mathbf1_I)(N)
                -H_N\mathfrak S_{m,b;m',b'}(N)\right|^2
         \ll H^3 X^{-10\delta}. \tag{11}
\]
Here \(\rho/4-4\kappa>10\delta\) and \(\rho-8\kappa>10\delta\).
The estimates are uniform for \(m,b,m',b'\le X^\kappa\) in \(C_{12}\).
If \(b\) or \(b'\) has a forbidden factor, the weighted count and its
local mass are zero.

## 3. The singular series on simultaneous good targets

Take grid starts \(x\) spaced \(H/2\). Their half-open central target
blocks \([2x+H-H/2,\,2x+H+H/2)\) tile the line. For large \(X\), blocks
meeting \((X,2X]\) have \(X/3\le x\le X\). For each fixed parameter
quadruple, (11) bounds the number of errors exceeding \(HX^{-\delta}\)
per block by \(O(HX^{-8\delta})\). There are at most \(X^{4\kappa}\)
quadruples and \(O(X/H)\) blocks. Outside \(O(X^{1-7\delta})\) targets,
therefore, all those correlations on the assigned block satisfy
\[
 (A_{m,b}\mathbf1_I\star A_{m',b'}\mathbf1_I)(N)
   =H_N\mathfrak S_{m,b;m',b'}(N)+O(HX^{-\delta}). \tag{12}
\]
Here \(H_N\ge H/2\), and any two summands in \(I\) satisfy the requested
balance. The simultaneous assertion permits sieve and cofactor choices
that depend on \(N\).

Fix a compatible target satisfying (12). Retain cofactor pairs
\[
 m=du,\qquad m'=dv,\qquad d\mid N,\quad d\ \text{squarefree},\qquad
 d,u,v\in C_{12},\quad d,u,v\le\sqrt z,\quad (u,v)=1,\quad(uv,N)=1 .
 \tag{13}
\]
Then \(d=(m,m')\), so each pair determines its parameters. Put \(N'=N/d\);
it is \(2\bmod12\). Solutions to \(uk+vl=N'\) form a progression
\[
 k(t)=k_0+vt,\qquad l(t)=l_0-ut .
\]
The two conditions \(k,l\equiv1\bmod12\) select the same single class of
\(t\bmod12\). At \(p\nmid12\), passing to this class preserves normalized
\(p\)-adic Haar measure.

The local measures in (5) are pushed forward by \(k\mapsto mk\) and
\(l\mapsto m'l\). For example the first has density in \(a\) equal to
\[
 p^{v_p(m)}\mathbf1_{p^{v_p(m)}\mid a}\,
 W_p(a/m)\,\mathbf1_{v_p(a/m)\ge v_p(b)}
\]
on the divisibility ball, and zero off it. Each density belongs to
\(L^2(\mathbb Z_p)\), since the weight at valuation \(i\) grows at most
cubically in \(i\). Its convolution with the second density, evaluated at
\(N\), equals the local integral of their product at \(a,N-a\).
The local Fourier series is absolutely convergent by Cauchy–Schwarz on
Fourier coefficients, and the convolution is continuous by \(L^2\)
continuity of translations. This identifies the Fourier expression at
every \(N\). CRT on finite prime sets, followed by the absolute
**denominator-group** convergence after (10), identifies the global
series with the product of these local integrals times \(C^2/(mm')\).
In taking finite-prime exhaustions one can include all primes of
\(12mm'bb'\) from the start.

On the common divisibility ball write \(a=a_0+[m,m']t\), where
\([m,m']\) denotes the least common multiple. This ball exists locally
because \(d\mid N\). Its Haar scale \(p^{-v_p([m,m'])}\), combined with
the two density multipliers, leaves \(p^{v_p(d)}\). The resulting quotient
forms are \(k(t),l(t)\) above up to Haar-preserving changes of parameter.
At \(p^j\parallel12\) the two congruence weights integrate to \(p^j\).
We obtain
\[
 \begin{split}
 \mathfrak S_{m,b;m',b'}(N)
 &=\frac{12C^2}{duv}\prod_{p\nmid12} I_p(b,b'),\\
 I_p(b,b')&=\int_{\mathbb Z_p}W_p(k(t))W_p(l(t))
    \mathbf1_{p^{v_p(b)}\mid k(t)}
    \mathbf1_{p^{v_p(b')}\mid l(t)}\,dt .
 \end{split} \tag{14}
\]
Write \(I_p=I_p(1,1)\). If \(p\nmid uv\), the proportion of \(t\)
giving two units is at least \(1-2/p\). If \(p\mid uv\), one form is
always a unit since \((uv,N')=1\), and the other has unit slope.
The values of \(W_p\) on units are uniformly bounded above and below.
Thus \(I_p\gg1\). Also \(I_p\ll1\) by the bounded \(L^2\) norms of
\(W_p\) and the same slope alternatives.

At \(p\nmid12uvN'\), a forbidden prime contributes
\((1-2/p)/(1-1/p)^2=1+O(p^{-2})\). At an allowed prime set
\(v_0=(1-1/p)^3=W_p(\text{unit})\). The functions \(W_p(k)-v_0\) and
\(W_p(l)-v_0\) have disjoint supports because their roots modulo \(p\)
are distinct. Both have integral \(1-v_0\), giving
\[
 I_p=2v_0-v_0^2=1-(1-v_0)^2=1+O(p^{-2}). \tag{15}
\]
This also verifies convergence and positivity of the unrestricted base
product in (14).


## 4. Lower sieve on the two rough factors

For a fixed pair (13), give \(k,l\) in the summand interval the nonnegative
weight \(w(k)w(l)\). At \(p\nmid12\) let \(g(p)\) be the \(I_p\)-normalized
local mass of \(p\mid k(t)l(t)\), and take \(g(p)=0\) at \(p\mid12\).
A forbidden prime has \(g(p)=0\) because its weight is supported on units.
For squarefree \(j\), expand \(j\mid kl\) using
\[
 \mathbf1_{p\mid kl}
   =\mathbf1_{p\mid k}+\mathbf1_{p\mid l}
                        -\mathbf1_{p\mid k}\mathbf1_{p\mid l}
\]
at each \(p\mid j\). The at most \(3^{\omega(j)}\) terms impose
\(b\mid k,\ b'\mid l\) with \(b,b'\mid j\). By (12) and (14), the base
weighted count with \(j\mid kl\), for \(j\le X^\kappa\), equals
\[
 H_N\mathfrak S_{m,1;m',1}(N)\prod_{p\mid j}g(p)+O(HX^{-\delta+o(1)}).
 \tag{16}
\]
This uses the local expansion inside each factor of (14). In particular the
main terms multiply after the local inclusion-exclusion.

We need a sieve dimension bound **uniform in \(N,u,v\)**. For an allowed
prime, the \(I_p\)-mass surviving both unit restrictions is at least
\((1-2/p)(1-1/p)^6\), so \(1-g(p)\gg1\). On a root residue modulo \(p\),
a linear form with unit slope which vanishes there has conditional second
moment of \(W_p\) bounded by
\(O(\sum_{i\ge0}(i+2)^6p^{-i})=O(1)\). A form nonzero there has constant
unit weight. If one slope is divisible by \(p\), its form is a unit on
the entire progression. Applying Cauchy–Schwarz on each of at most two
root residues gives \(I_p g(p)\ll1/p\), hence \(g(p)\ll1/p\).

At \(p\mid N'\) with \(p\nmid uv\) the two roots coincide. For example,
when \(p^2\mid N'\), both forms have valuation one on \(1-O(1/p)\) of the
common root residue. The two factors \(w(p)=4\) then both enter its
conditional mass. The bound below accommodates this possible first-order
coefficient \(16\) in \(p g(p)\). Choose an absolute \(K_0\) with
\(-\log(1-g(p))\le K_0/p\) at all primes, increasing it to cover the
finitely many small ones. Mertens' estimate on every interval gives an
absolute \(C_0\) such that
\[
 \prod_{y_1<p\le y_2}(1-g(p))^{-1}
 \le \left(1+\frac{C_0}{\log y_1}\right)
       \left(\frac{\log y_2}{\log y_1}\right)^{K_0},
       \qquad 3/2\le y_1\le y_2. \tag{17}
\]
Indeed \(\sum_{y_1<p\le y_2}1/p=
\log(\log y_2/\log y_1)+O(1/\log y_1)\), with a bounded adjustment for
\(y_1<2\). The large uniform \(K_0\) is sufficient.

Choose \(s>1+2/(\exp(0.53/K_0)-1)\) fixed and sufficiently large,
\(r=\kappa/(2s)\), \(z=X^r\), and \(P(z)=\prod_{p\le z}p\).
The fundamental lemma gives lower weights \(|\lambda_j^-|\le1\)
supported on \(j\mid P(z),\ j\le z^s\) such that
\[
 \sum_{j\mid n}\lambda_j^-\le\mathbf1_{(n,P(z))=1},\qquad
 \sum_j\lambda_j^-\prod_{p\mid j}g(p)
       \ge 0.011\prod_{p\le z}(1-g(p)).
\]
Apply the cited theorem with prime density parameter \(\nu(p)=p g(p)\).
Its hypotheses are \(0\le g(p)<1\), (17), and sufficiently large
\(z^s\) for the fixed constants. Its lower factor is \(11/10^3\).
[Koukoulopoulos, *The Distribution of Prime Numbers*, Chapter 18,
Axiom 2, and Theorem 19.1(a)–(c), pp. 195–196](https://dms.umontreal.ca/~koukoulo/documents/publications/primes.pdf#page=206).

For \(j\le z^s\) all its divisors \(b,b'\) lie below \(X^\kappa\).
Use the first sieve inequality on the nonnegative sequence in (16) and
the second on the resulting main term. The sifted count for (13) is at
least
\[
 0.011 H_N\mathfrak S_{m,1;m',1}(N)
     \prod_{p\le z}(1-g(p))-O(H X^{-\delta+\kappa+o(1)}). \tag{18}
\]
The error allows all \(j\le z^s\) and all their signed assignments in
(16).

To evaluate the survival factors, let \(\nu_p^\ast\) be the number of root
residues for \(k(t)l(t)\bmod p\) for \(p\nmid12\). There is one root if
\(p\mid uv\) or \(p\mid N'\), and two otherwise. This follows from unit
slopes when \(p\nmid uv\) and the always-unit form when \(p\mid uv\).
At a forbidden prime,
\[
 I_p=\frac{1-\nu_p^\ast/p}{(1-1/p)^2}.
\]
At an allowed prime \(p\le z\), retaining the unit residues gives the
**exact** factor
\[
 I_p(1-g(p))
     =(1-\nu_p^\ast/p)(1-1/p)^6
     =(1-1/p)^8\frac{1-\nu_p^\ast/p}{(1-1/p)^2}. \tag{19}
\]
By Mertens' estimate for the fixed class \(1\bmod12\),
\(\prod_{p\le z,\,p\in\mathcal P}(1-1/p)^8\asymp(\log z)^{-2}\).
The exponent is \(8/\varphi(12)=2\).
[Loughran–Matthiesen, Lemma 2.4(2)](https://arxiv.org/html/1904.12845v2#S2).

For every \(p\nmid12\) the normalized two-root factor is
\(1-1/(p-1)^2\), and the product of these lower factors is positive.
A single root increases it by \((p-1)/(p-2)\ge1+1/p\). In particular
this gain applies at every allowed \(p\le z\) dividing \(N'\).
For allowed \(p>z\) outside \(uvN'\), (15) gives an absolute positive
lower bound for their product for large \(z\). At allowed \(p>z\) dividing
\(uvN'\), the two-unit mass yields \(I_p\ge1-O(1/p)\). There are
\(O(\log X/\log z)\) such primes, so their product is at least \(1/2\)
for large \(X\), uniformly in the target and cofactors. Together with
\(H_N\ge H/2\), (14) and (19) bound the **main term** in (18) from below by
\[
 \gg \frac{H}{duv(\log z)^2}
       \prod_{\substack{p\le z\\p\in\mathcal P,\ p\mid N'}}(1+1/p).
 \tag{20}
\]

## 5. Recovering the cofactor mass and the unweighted pairs

There are at most \(z^2\le X^{2\kappa}\) distinct pairs \(m,m'\) retained
in (13). The total subtracted error in (18) is
\(O(HX^{-\delta+3\kappa+o(1)})=o(H/(\log X)^2)\).
It remains to sum (20) uniformly, including targets having many allowed
prime factors.

Choose a sufficiently large absolute \(C_1\), put \(y=z^{1/C_1}\), and
first allow coprime \(u,v\) formed from primes \(p\in\mathcal P,\ p\le y\),
\(p\nmid N\). Before imposing their size cutoffs the harmonic mass is
\[
 \sum_{u,v}\frac1{uv}
  =\prod_{\substack{p\le y\\p\in\mathcal P,\ p\nmid N}}
       (1-p^{-2})(1-p^{-1})^{-2}
  \gg (\log X)^{1/2}
       \prod_{\substack{p\le y\\p\in\mathcal P,\ p\mid N}}(1-1/p)^2 .
 \tag{21}
\]
Here \(\log y\) is a fixed positive multiple of \(\log X\). Under the
probability measure obtained by normalizing this Euler product, a given
prime divides a specified coordinate with probability \(1/(p+1)\).
Conditional on that event its expected exponent is \(p/(p-1)\).
The expected logarithm of each coordinate is consequently at most
\(\sum_{p\le y}\log p/(p-1)\le2\log y+O(1)\). Markov's inequality and
the union bound give \(u,v\le\sqrt z\) simultaneously on at least half
the mass for \(C_1\) sufficiently large.

Let \(T_N=\{p\in\mathcal P:p\mid N,\ p\le\log X\}\).
All prime divisors of \(N\) above \(\log X\) have reciprocal sum \(O(1)\),
since there are at most \(\log(2X)/\log 2\) of them. Both \(y\) and \(z\)
exceed \(\log X\) for large \(X\). The product on the right of (21) is
therefore bounded below by a fixed multiple of
\(\prod_{p\in T_N}(1-1/p)^2\).

Restrict \(d\) to squarefree products of primes in \(T_N\). Each \(p\in T_N\)
not dividing \(d\) divides \(N'=N/d\), so the enhancement in (20) is at
least
\[
 \prod_{p\in T_N}(1+1/p)\prod_{p\mid d}(1+1/p)^{-1}.
\]
Summing \(d^{-1}\prod_{p\mid d}(1+1/p)^{-1}\) unrestrictedly gives
\(\prod_{p\in T_N}(1+1/(p+1))\). Under these normalized weights,
\(\mathbb E\log d=\sum_{p\in T_N}\log p/(p+2)=O(\log\log X)\).
Thus \(d\le\sqrt z\) retains at least half this mass, uniformly in \(N\).
The three factors at each \(p\in T_N\) multiply to
\[
 (1-1/p)^2(1+1/p)(1+1/(p+1))=1-3/p^2+2/p^3 .
\]
Their product over any \(T_N\) has an absolute positive lower bound.
This restores the potential loss from the coprimality condition in (21).
Summing (20) gives \(\gg H(\log X)^{1/2}/(\log z)^2
\gg H/(\log X)^{3/2}\), and the summed error in (18) is smaller.
The cofactor restoration also appears in the logarithmic-exception
argument for more general windows.
[Semigroup density-one theorem](PRESCRIBED_PRIME_SEMIGROUP_DENSITY_ONE.md).

Finally, \(m\le z\), and each retained \(k\) has every prime factor
greater than \(z=X^r\). Since \(k\le2X\) for large \(X\),
\(\Omega(k)\le 2/r\), counting multiplicity. Therefore
\[
 w(k)\le\tau_4(k)\le4^{\Omega(k)}\le4^{2/r}.
\]
In \(a=mk\), the integer \(m\) is exactly the part of \(a\) composed of
prime powers with prime at most \(z\). This proves uniqueness of the
cofactor for each summand, including multiplicities. The same holds for
\(m'l\). Each retained ordered pair is supported on \(C_{12}\) and receives
weight at most \(4^{4/r}\). Dividing by that fixed constant proves the
unweighted lower bound. The \(O(X^{1-7\delta})\) exceptional targets from
(12) give (D12P) with the stated conservative exponent \(10^{-6}\).

## Scope and review

The weighted mean estimates, the short-strip Fourier bounds, the
denominator-grouped singular series, the lower-sieve hypotheses and
the restoration at primes dividing \(N\) received independent internal
mathematical audits. The argument treats \(q=12\) and \(H=X^{0.99}\).
A variance estimate for the full unweighted count with a power-saving
error would require a suitably refined main term, accounting for logarithmic
secondary terms and the changing local density across the dyadic interval.
That problem and an every-target lower bound remain open here. No practical cutoff has been computed.
External mathematical review and historical priority remain outstanding.
