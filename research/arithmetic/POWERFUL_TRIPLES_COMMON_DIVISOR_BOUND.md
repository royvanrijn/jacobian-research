# Powerful triples of large order with arbitrary common divisor

An integer is \(s\)-full if each of its nonzero prime valuations
is at least \(s\).

## Theorem. The sharp exponent with arbitrary common divisor

For all sufficiently large fixed \(s\), the
number of ordered positive \(s\)-full triples \(a+b=c\le X\),
with arbitrary common divisor, is \(O_{s,\varepsilon}(X^{1/s+\varepsilon})\).
Its exponent is sharp, since
\[
(a,b,c)=(2^s t^s,2^s t^s,2^{s+1}t^s)
\]
gives \(\gg_s X^{1/s}\) distinct triples for positive integers \(t\).
Any triple of orders \(p,q,r\) with minimum \(s\) belongs to
this count.

## Lemma. Keeping the deficient primes after division

Write \(g=\gcd(a,b,c)\) and \(u=a/g,v=b/g,w=c/g\), so \(u+v=w\)
is primitive. If a prime \(\ell\) divides \(g\), its valuation in
\(g\) is the minimum of three positive valuations at least \(s\).
Thus \(g\) is \(s\)-full. Let \(D\) be the product of primes whose
positive valuation on one of \(u,v,w\) is less than \(s\). Such a
prime must divide \(g\), as otherwise its valuation stayed
unchanged. We also have
\[
 \operatorname{rad}(uvw)\le D(uvw)^{1/s}<D w^{3/s}.
 \tag{1}
\]

For squarefree \(D\), the number of \(s\)-full \(g\le G\) divisible
by \(D\) is
\[
 \ll_{s,\alpha}G^{1/s}D^{-1+\alpha}\qquad(\alpha>0).
 \tag{2}
\]
Write \(g=\prod_{e=s}^{2s-1}z_e^e\) by the standard powerful
factorization. For a prime with valuation \(qs+j\), \(0\le j<s\),
use exponent \(q\) in \(z_s\) when \(j=0\), and otherwise use
exponent \(q-1\) in \(z_s\) together with one factor in \(z_{s+j}\).
In any representation each prime of \(D\) divides
some root \(z_e\). Assign it to one such root. Write \(D_e\) for
the product of the primes assigned to \(e\), and put \(z_e=D_e y_e\).
For a fixed assignment the number of tuples is at most
\[
 \left(\frac{G}{\prod_e D_e^e}\right)^{1/s}
      \prod_{e=s+1}^{2s-1}\left(\sum_{n\ge1} n^{-e/s}\right)
 \ll_s G^{1/s}/D .
\]
The series converge, and \(s^{\omega(D)}\ll_{s,\alpha}D^\alpha\)
bounds the number of assignments. This proves (2).
[Canonical powerful factorization](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#fourier-counting-in-a-box)

## Excluding the remaining limiting profile

Suppose the stated upper bound fails for an unbounded sequence of
orders \(s\). For each fixed such \(s\), some \(\varepsilon_s>0\)
gives failures at unbounded \(X\). Dyadically split into bins
\(g\asymp G, w\asymp Y, D\asymp H\). Let \(M_s(H,Y)\) be the count
of primitive quotients in the last two bins. Equation (2) bounds
the number of triples in the bin by
\(\ll_{s,\alpha}G^{1/s}H^{-1}X^\alpha M_s(H,Y)\), since \(D\le g\le X\).
Choosing \(\alpha\) small relative to \(\varepsilon_s\) and absorbing
the \(O((1+\log X)^3)\) bins gives bins with
\[
 M_s(H,Y)\gg_s Y^{1/s}H X^{\varepsilon_s/2},
\]
after reducing \(\varepsilon_s\) if necessary. Here \(GY\ll X\).
In particular \(Y\) tends to infinity on a subsequence. At fixed
\(s\) further select a subsequence with \(h_s=\lim \log_Y H\), which
exists since \(1\le H\ll Y^3\). Along this subsequence the quotient
count has limiting exponent at least \(1/s+h_s\).

For every fixed \(\rho>0\), (1) puts these quotients in the \(abc\)
count with height \(2Y\) and budget \(3/s+h_s+\rho\) for sufficiently
large \(Y\).
Using the \(2\lambda/3\) bound and then taking \(\rho\) down to zero
at fixed \(s\) gives \(h_s\le3/s\). Next, for each fixed \(\eta>0\)
and all sufficiently large \(s\), the limiting one-third count
applies at budgets just above \(3/s+h_s\). It gives
\[
 1/s+h_s \le (1/3+\eta)(3/s+h_s).
\]
For \(0<\eta<2/3\) this implies
\(s h_s\le3\eta/(2/3-\eta)\), and hence \(s h_s\to0\).
[The elementary and limiting small-budget bounds](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#exceptional-set-estimates-and-dyadic-coordinates)

Now form anatomic boxes for the primitive quotients, with truncation
loss \(o(1/s)\). The compression index of the anatomic factorization
can be chosen at least \(s\). Thus the roots at indices \(e<s\)
group exactly the primes in \(D\), by their actual valuations in
the quotients. Their total logarithmic root mass in the box is
at most \(h_s+o(1/s)\). The total root mass \(L_s\) of a populated
box is at most \(3/s+h_s+o(1/s)\) by (1) and the anatomic radical
comparison. After paying the dyadic and truncation losses at
fixed \(s\) and then sending \(s\to\infty\), a box attaining the
quotient count has exponent at least \(1/s+h_s-o(1/s)\).
The limiting upper bound \(D_{\mathrm{box}}/L_s\le1/3+o(1)\) forces
\(D_{\mathrm{box}}/L_s\to1/3\) and \(sL_s\to3\).

Consequently the total normalized mass below degree \(s\) tends to
zero. At every other coordinate the scaled degree \(eL_s\) is at
least \(sL_s\), so all three weak-limit measures are supported
on \([3,\infty)\). Balance fixes each mass at \(1/3\), and the
numerical height bounds each first moment by one. Thus all three
are point masses at three, contradicting the endpoint lower bound
\(1/2\) on their summed nonatomic mass.
[Anatomic factorization, Lemma 2.2](sources/bernert-browning-lichtman-teravainen-2410.12234v2.pdf) · [The endpoint diffusion theorem](ABC_WIDE_COFACTOR_CONGRUENCES.md#necessary-conditions-at-normalized-exponent-13)

## Asymptotic expansion for comparable orders

Let the fixed orders be \(s_A,s_B,s_C\), with
\(m=\min s_i\), \(M=\max s_i\). There is an absolute \(\theta>0\) such that for
\(M/m\le1+\theta\) and large \(m\), the count with arbitrary
common divisor has an asymptotic expansion
\[
 A_{\boldsymbol s}(X)=
 \sum_{j=0}^J C_{\boldsymbol s,j}X^{1/(M+j)}
 +O_{\boldsymbol s,J,\varepsilon}\bigl(X^{1/(M+J+1)+\varepsilon}\bigr),
 \qquad (-1)^j C_{\boldsymbol s,j}>0 .
 \tag{PA}
\]
Here \(A_{\boldsymbol s}(X)\) counts ordered positive triples
with the three prescribed fullness orders and \(c\le X\).
Every fixed \(J\) is allowed for sufficiently large orders.
In fact, with the constant \(\kappa\) below, it suffices that
\(M\ge16(J+1)^2\) and \(M\ge100(J+1)/\kappa\), in addition
to an absolute lower cutoff on \(m\). In particular the number
of terms allowed is \(\gg\sqrt M\). The implied constant need
not be uniform in the orders or \(J\).

This makes the maximum order \(M\) govern the leading exponent
even when the orders are unequal. All three equal orders give
the count in the opening theorem. The single-integer
asymptotic underlying the calculation is classical. The proof
uses a coprimality-uniform form and a summable estimate for
primitive additive quotients. [Cohen–Davis, §§6–7](https://acta.bibl.u-szeged.hu/14284/1/math_031_fasc_003_004_363-371.pdf)

### Uniform bound for primitive quotients

For a primitive quotient \(u+v=w\), let \(P=\operatorname{rad}(uvw)\)
and let \(D\) be the squarefree product of primes whose valuation
on their quotient letter is below its required order. For
\(w\asymp Y,D\asymp H\), there are absolute
\(\kappa,\theta>0\) such that the number of these quotients,
at sufficiently large comparable fixed orders, is
\[
 \ll_{\boldsymbol s,\epsilon}
 Y^{(1-\kappa)/M+\epsilon}H^{m/M}.
 \tag{PQ}
\]
Indeed otherwise choose counterexample orders with
\(m\to\infty,\ M/m\to1,\ \kappa\to0\). At each fixed
triple of orders, failure for some fixed positive exponent
allowance gives arbitrarily large \(Y\). Pass to a subsequence
on which \(h=\lim\log_Y H\) exists, since \(1\le H\ll Y^3\). The quotient exponent is at least
\((1-\kappa)/M+(m/M)h\). Write \(\lambda_0=\sum_i1/s_i\).
The radical bound is \(\operatorname{rad}(uvw)\le D w^{\lambda_0}\).
The \(2\lambda/3\) bound, at budgets just above \(h+\lambda_0\),
first gives \(h=O(1/M)\). Then the limiting
\((1/3+\eta)\lambda\) bound gives
\[
 (m/M-\tfrac13-\eta)Mh
 \le(\tfrac13+\eta)M\lambda_0-(1-\kappa).
\]
Let the orders tend to infinity and then \(\eta\downarrow0\),
to obtain \(Mh\to0\).

Use anatomic boxes with truncation loss \(o(1/M)\) and
compression index above all the required orders. Every
coordinate on a letter below its order consists only
of its deficient primes, so their total logarithmic
root mass is at most \(h+o(1/M)\). The total root mass
is at most \(\lambda_0+h+o(1/M)\), while the chosen box
exponent is at least \(1/M-o(1/M)\). The limiting
one-third upper coefficient forces \(ML\to3\) and
normalized exponent \(d\to1/3\). The absence of
normalized low root mass then forces all three weak
measures to be supported in \([3,\infty)\).
They each have mass \(1/3\) and first moment at most
one, hence all are point masses at three. This
contradicts endpoint diffusion. All height limits and
the choice of \(h\) precede passage to varying orders,
so the contradiction gives fixed positive \(\kappa,\theta\).

### Counting the coprime part of the common divisor

Write the common divisor as \(g=f b\) with all prime factors
of \(f\) dividing \(P\) and \(\gcd(b,P)=1\). The latter
factor \(b\) is \(M\)-full. At each prime \(\ell\mid P\), let
\(K_\ell=\max_i(s_i-v_\ell(u_i))\), with
\((u_A,u_B,u_C)=(u,v,w)\). A positive exponent of \(f\)
must be at least \(K_\ell\), and exponent zero is allowed
exactly when \(\ell\nmid D\). Primitivity gives \(K_\ell\ge m\).
Thus, for any
fixed \(\sigma>0\), its allowed series satisfies
\[
 \sum_f f^{-\sigma}
 \le D^{-m\sigma}\prod_{p\mid P}(1-p^{-\sigma})^{-1}
 \ll_{\sigma,\epsilon}D^{-m\sigma}w^\epsilon.
 \tag{SF}
\]
This follows by comparing each local factor outside a
finite set with \(p^{\epsilon/3}\), and using \(P\le w^3\).

Let \(B_M(Y;P)\) count \(M\)-full integers at most \(Y\)
coprime to \(P\), and put \(a_i=M+i\) for \(0\le i\le J\).
Write \(\zeta_P(t)=\zeta(t)\prod_{p\mid P}(1-p^{-t})\).
The Dirichlet series of its indicator is
\[
 \prod_{p\nmid P}\left(1+\sum_{\ell\ge M}p^{-\ell s}\right)
   =G_M(s;P)\prod_{i=0}^J\zeta_P(a_i s),
\]
where \(G_M\) has absolutely summable Dirichlet coefficients
on each real \(\sigma>1/(M+J+1)\), uniformly in \(P\).
Indeed the local product with \(\prod_i(1-p^{-a_i s})\)
cancels all positive exponents below \(M+J+1\), for
\(M>J\). The local factors and hence \(G_M\) are positive
at each of the real points \(1/a_i\).

For completeness, if \(I\subset\{0,\ldots,J\}\) is nonempty,
put \(d=|I|\), \(A=\sum_{i\in I}a_i\). The count of
\(\prod_{i\in I}x_i^{a_i}\le Y\) with all \(x_i\) coprime
to \(P\) is
\[
 \rho_P\sum_{i\in I}Y^{1/a_i}
       \prod_{j\in I\setminus\{i\}}\zeta_P(a_j/a_i)
 +O_{M,J}\left(2^{d\omega(P)}
                   Y^{e_I}\right),
 \tag{DH}
\]
where \(\rho_P=\prod_{\ell\mid P}(1-1/\ell)\) and the
error exponent in (DH) can be taken to be
\(e_I=(d-1)(1+J/M)/A\).
For \(d=1\) the elementary coprime count gives error
\(O(2^{\omega(P)})\). For \(d\ge2\), put \(U=Y^{1/A}\).
At least one coordinate is at most \(U\), so apply
inclusion-exclusion over the nonempty sets \(K\) of these
coordinates. For proper \(K\), apply induction to the
complement \(B=I\setminus K\). The resulting error sums
converge with the required exponent since \(a_j e_B<1\)
for each \(j\in K\). Indeed, for \(b=|B|\),
\[
 a_j e_B\le(1-1/b)(1+J/M)^2<1
\]
under \(M\ge16(J+1)^2\). Summing these errors over
the small coordinates gives \(U^{|K|+e_B\sum_{i\in B}a_i}\),
at most \(U^{(d-1)(1+J/M)}\), with the stated
\(2^{d\omega(P)}\) factor.

For each main term use, at real \(t>0,\ t\ne1\),
\[
 \sum_{\substack{x\le U\\(x,P)=1}}x^{-t}
   =\zeta_P(t)+\rho_P\frac{U^{1-t}}{1-t}
       +O_t(2^{\omega(P)}U^{-t}).
\]
The unrestricted weighted sum follows by partial summation
with the zeta constant continued to \(t<1\).
Inclusion-exclusion over the divisors of \(P\) gives the
display uniformly. Extending the divisor sums from
\(\ell\le U\) to all divisors costs
\(O_t(2^{\omega(P)}U^{-t})\) in both leading terms.

For a fixed pole index \(i\), put
\(Z_j=\zeta_P(a_j/a_i)\) and
\(R_j=\rho_P U^{1-a_j/a_i}/(1-a_j/a_i)\) for \(j\ne i\).
The alternating sum over nonempty \(K\subseteq I\setminus\{i\}\)
of the products of \(Z_j+R_j\) on \(K\) and \(Z_j\) elsewhere
is \(\prod_{j\ne i}Z_j+(-1)^d\prod_{j\ne i}R_j\).
The first term gives the coefficient in (DH).
Upon multiplication by \(Y^{1/a_i}\), the second has
coefficient
\((-1)^d\rho_P^d a_i^{d-1}/\prod_{j\ne i}(a_i-a_j)\)
on \(U^d\). Lagrange interpolation gives
\(\sum_i a_i^{d-1}/\prod_{j\ne i}(a_i-a_j)=1\).
Thus these terms cancel the all-small box, whose leading
term is \((-1)^{d+1}\rho_P^d U^d\).
In any product using a remainder in the \(j\)-th slot,
the factor \(U^{-a_j/a_i}\) cancels that slot from the
initial \(U^{A/a_i}\). Each of the other \(d-1\) slots
then contributes an exponent at most \(1+J/M\).
This gives the same \(U^{(d-1)(1+J/M)}\) bound for
the remainder in the main-term calculation. This proves (DH).
Finally \(e_I<1/(M+J+1)\) for the full \(I\) under
\(M\ge16(J+1)^2\), using
\((1-1/d)(1+(J+1)/M)^2<1\).

Convolve (DH) with the absolutely summable coefficients
of \(G_M\), including main term tails beyond \(Y\).
For \(1/(M+J+1)<\sigma<1/(M+J)\), the result is
\[
 B_M(Y;P)=\sum_{j=0}^J b_j(P)Y^{1/a_j}
                 +O_{M,J,\sigma}(2^{(J+1)\omega(P)}Y^\sigma),
\]
where
\[
 b_j(P)=\rho_P G_M(1/a_j;P)
                \prod_{i\ne j}\zeta_P(a_i/a_j).
\]
Exactly \(j\) of the zeta arguments lie in \((0,1)\);
their values are negative and the other factors are
positive, so \((-1)^j b_j(P)>0\).
The coefficients are uniformly bounded in \(P\),
and the estimate extends to \(0<Y<1\), where the
actual count is zero and \(Y^{1/a_j}\le Y^\sigma\).

### Summing the quotient and common factors

The sought count is \(\sum_{u+v=w\ \mathrm{primitive}}
\sum_{f\ \mathrm{allowed}} B_M(X/(wf);P)\).
To justify summing its error and main terms, choose
\(\sigma>1/(M+J+1)\) sufficiently close to that
endpoint, and put \(h=\log_Y H\). By (SF), the quotient
summability condition for the error is bounded by
the dyadic sums of \(Y^{-\sigma+\epsilon}H^{-m\sigma}\)
times the quotient count. For \(h\le20/M\), (PQ) supplies
exponent at most
\[
 -\kappa/M+(1/M-\sigma)(1+mh)+O(\epsilon)
 \le-\kappa/M+\frac{21(J+1)}{M^2}+O(\epsilon)<0 .
\]
For \(h>20/M\), use the elementary two-smallest-radicals
count \(O_\epsilon(Y^{2\lambda_0/3+\epsilon}H^{2/3})\).
This is uniform in \(H\): the product of the two smallest
radicals is at most \((D w^{\lambda_0})^{2/3}\), and fixing
their complete values determines the third. In that branch
the dyadic exponent is at most
\[
 (2/3)\lambda_0-\sigma-(m\sigma-2/3)h+O(\epsilon)<0
\]
with a fixed margin of order \(1/M\), since
\(\lambda_0\le3.3/M\) and \(m\sigma-2/3\ge1/6\).
The latter follows from \(M/m\le1.1\),
\(M\ge16(J+1)^2\), and \(\sigma>1/(M+J+1)\).
At each fixed order take \(\epsilon\) sufficiently small.
The logarithmically many \(H\)-bins do not affect
summability. The same bound implies absolute convergence
of all main coefficients, since \(1/a_j>\sigma\).
Explicitly,
\[
 C_{\boldsymbol s,j}
 =\sum_{u+v=w\ {\rm primitive}} b_j(\operatorname{rad}(uvw))
                       \,w^{-1/a_j}\sum_{f\ {\rm allowed}} f^{-1/a_j}.
 \tag{PC}
\]
The terms have a common strict sign since every quotient
allows positive choices of \(f\), and the quotient \(1+1=2\)
already gives a nonzero term. Choosing \(\sigma\) as close
to \(1/(M+J+1)\) as needed proves (PA).

In particular, after dividing by \(X^{1/M}\), the count
coming from any fixed primitive quotient converges to its
term in \(C_{\boldsymbol s,0}\). The limiting distribution
on primitive quotients has total mass one after normalization
by that constant. Its height tail is a power saving by the
same dyadic summability argument.

## An arbitrary-order upper bound with common divisors

Let \(s_A,s_B,s_C\) be any sufficiently large fixed fullness
orders, let \(M=\max_i s_i\), and put \(\lambda_0=\sum_i1/s_i\).
There is an absolute \(0<\delta<1/3\) such that the full-triple
count with arbitrary common divisor satisfies
\[
 A_{\boldsymbol s}(X)
 \ll_{\boldsymbol s,\varepsilon}
 X^{\max\{\,1/M,\ (1/3-\delta)\lambda_0\,\}+\varepsilon}.
 \tag{GU}
\]
No comparability of the orders is needed. In the strongly
asymmetric range this materially improves the \(1/\min_i s_i\)
bound obtained by forgetting the larger orders.

For a primitive quotient \(u_A+u_B=u_C=W\), let \(\mathcal D\) be the
set of its deficient primes, and define
\[
 K_\ell=\max_i(s_i-\nu_\ell(u_i)),\qquad
 Z=\prod_{\ell\in\mathcal D}\ell^{K_\ell}.
\]
Here a deficient prime has nonzero valuation \(e<s_i\)
on its unique quotient letter \(i\). If that letter has
maximum order then \(K_\ell/M\ge1-e/s_i\). Otherwise
a maximum-order letter has valuation zero and
\(K_\ell=M\). Consequently, prime by prime,
\[
 \operatorname{rad}(u_Au_Bu_C)
 \le Z^{1/M}\prod_i u_i^{1/s_i}
 \le Z^{1/M} W^{\lambda_0}.
 \tag{GU1}
\]
For a prime of degree \(e<(1-\vartheta)s_i\) the same comparison
gives \(K_\ell/M\ge\vartheta\). Consequently the sum of the
logarithms of those primes on any letter is at most
\(\log Z/(M\vartheta)\), for every \(0<\vartheta<1\).

In dyadic bins \(W\asymp Y,\ Z\asymp H\), the quotient
count satisfies the uniform estimate
\[
 Q(Y,H)\ll_{\boldsymbol s,\varepsilon}
      Y^{(1/3-\delta)\lambda_0+\varepsilon} H^{1/M}.
 \tag{GU2}
\]
If not, choose contrary sequences with all three orders
tending to infinity and \(\delta\to0\). For each fixed
order triple take an unbounded sequence of heights witnessing
failure with a fixed positive exponent allowance.
Let \(z=\lim \log_Y H/M\) on a subsequence of these
dyadic bins. By (GU1) the radical count is at budget
\(\lambda_0+z+o(1)\), while the claimed exponents
would reach \((1/3-\delta)\lambda_0+z\).
The elementary \(2\lambda/3\) bound first gives
\(z\le(1+3\delta)\lambda_0\), by applying it at any fixed
budget above \(\lambda_0+z\) and decreasing the excess.
For any fixed \(\eta>0\), the limiting one-third bound
gives \((2/3-\eta)z\le(\eta+\delta)\lambda_0\) for
sufficiently large orders. Thus \(z=o(\lambda_0)\).

Extract anatomic boxes with truncation loss
\(o(\lambda_0)\) and compression above all prescribed
orders. Total root mass satisfies \(L\le(1+o(1))\lambda_0\)
and the attained exponent is at least
\((1/3-o(1))\lambda_0\). The one-third box bound gives
\(L/\lambda_0\to1\) and \(D_{\rm box}/L\to1/3\).
By the weighted low-degree bound following (GU1),
the normalized root mass at unscaled indices below
\((1-\vartheta)s_i\) is \(o(1)\) on each letter.
The prelimit first moment is at least
\((1-\vartheta)s_iL(1/3-o(1))\) by endpoint balance,
and at most \(1+o(1)\). Hence
\(\limsup s_i\lambda_0\le3/(1-\vartheta)\).
Send \(\vartheta\downarrow0\). The reciprocals of
the three \(s_i\lambda_0\) sum to one, so all tend
to three. Every weak measure is therefore supported
on \([3,\infty)\) and, by mass and moment, equals
\((1/3)\delta_3\). This contradicts endpoint
diffusion and proves (GU2).

Finally at fixed quotient, split the common divisor
into its part \(f\) supported on
\(P=\operatorname{rad}(u_Au_Bu_C)\) and its part
coprime to \(P\). The latter is \(M\)-full.
At every deficient prime the part \(f\) must incorporate
the full exponent \(K_\ell\). The same Euler product argument
as in (SF) therefore gives, for any fixed \(\sigma>0\),
\[
 \sum_{f\ {\rm allowed}} f^{-\sigma}
 \le Z^{-\sigma}\prod_{\ell\mid P}(1-\ell^{-\sigma})^{-1}
 \ll_{\sigma,\varepsilon}Z^{-\sigma}W^\varepsilon .
 \tag{GU3}
\]
The \(M\)-full integer upper bound gives at most
\(O((X/W)^{1/M}Z^{-1/M}W^\varepsilon)\)
common divisors for this quotient. Multiply by (GU2)
and sum over \(Y\le X\) and \(1\le H\ll Y^{3M}\).
There are only logarithmically many dyadic scales at
fixed orders. The remaining height exponent is
\((1/3-\delta)\lambda_0-1/M\), giving the maximum
in (GU). Relative to the older \(1/\min_i s_i\)
bound, this can give close to a factor-three
improvement when one order is much smaller than
both others.

## An asymptotic with a weighted primitive height constraint

For a primitive quotient \(u_A+u_B=u_C=W\) write
\(h_i=\log u_i/\log W\). Fix the orders \(s_i\ge2\),
put \(M=\max_i s_i\), and count triples whose primitive
quotient satisfies
\[
 \sum_i \frac{M h_i}{s_i}\le3+\tau .
 \tag{WH}
\]
There are absolute \(\tau,\kappa>0\) such that this count
has the expansion (PA) for sufficiently large \(M\), with
the coefficients in (PC) summed over the quotients in
(WH). Its error is
\(O_{\boldsymbol s,J,\varepsilon}(X^{1/(M+J+1)+\varepsilon})\)
under \(M\ge16(J+1)^2\) and \(M\ge100(J+1)/\kappa\), in
addition to an absolute lower cutoff on \(M\).
The sector is empty if \(M/s_C>3+\tau\). Otherwise it
contains the quotient \(1+1=2\), so all coefficients have
the same strict alternating signs as in (PA).

This allows substantially separated orders. When
\(s_B=s_C=M=2s_A\), for instance, every quotient with
\(u_A\le W^{1/2}\) satisfies (WH). For near-equal orders
with \(M\sum_i1/s_i\le3+\tau\), the restriction holds
for all primitive quotients.

### A quotient bound in the weighted height sector

The first inequality in (GU1) gives
\[
 \operatorname{rad}(u_Au_Bu_C)
 \le W^{(3+\tau)/M}Z^{1/M}.
\]
For quotients in (WH) with \(W\asymp Y\), \(Z\asymp H\),
there are absolute \(\tau,\kappa>0\) such that, for all
sufficiently large \(M\),
\[
 Q_{\mathrm{wh}}(Y,H)
 \ll_{\boldsymbol s,\varepsilon}
       Y^{(1-\kappa)/M+\varepsilon} H^{1/M}.
 \tag{WH1}
\]
Suppose otherwise and take a sequence with \(M\to\infty\) and
\(\tau,\kappa\to0\). At each fixed order triple take the height
limit, passing to a subsequence with \(z=\lim\log_Y H\).
The elementary bound gives
\((1-\kappa)+z\le 2(3+\tau+z)/3\), so \(z\) is bounded.
The limiting one-third count, at budget approaching
\((3+\tau+z)/M\), gives \(z\to0\) by letting the fixed
allowance decrease after the order limit.

Take anatomic boxes with truncation error \(o(1/M)\) and
compression index above all the orders. The bound on
radicals gives \(L\le(3+o(1))/M\), and the count exponent
is at least \(1/M-o(1/M)\), so \(ML\to3\) and \(d\to1/3\).
At this endpoint [the short-height bound](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-strict-saving-when-the-short-term-is-small)
keeps the limiting heights of all three letters above
an absolute positive constant. By (WH), every \(s_i/M\)
is then bounded below. Extract limits \(a_i=\lim s_i/M>0\)
and \(h_i^*=\lim h_i\).
The weighted low-degree inequality after (GU1) shows
that the normalized root mass below \((1-\vartheta)s_i\)
tends to zero for every fixed \(\vartheta>0\). Hence
each weak measure is supported on \([3a_i,\infty)\).
It has mass \(1/3\) and first moment at most \(h_i^*\),
so \(h_i^*\ge a_i\). But \(\sum_i h_i^*/a_i\le3\) by
(WH), forcing equality in all three. The measures
are point masses, contrary to endpoint diffusion.
This proves (WH1) uniformly without an order-ratio
hypothesis.

### Summing the expansion in this sector

Use the same coprime-to-\(P\) \(M\)-full expansion as in
(PA). By (GU3), its error and the terms extending past
the height cutoff are bounded by \(X^\sigma\) times the
sum of \(W^{-\sigma+\varepsilon}Z^{-\sigma}\) over the
quotients in (WH), at any fixed
\(1/(M+J+1)<\sigma<1/(M+J)\).
For \(z=\log_Y H\le20\), (WH1) bounds the dyadic
exponent by
\[
 -\kappa/M+(1/M-\sigma)(1+z)+O(\varepsilon)<0.
\]
For \(z>20\), the elementary two-smallest-radicals
estimate bounds it by
\[
 \frac{2(3+\tau)}{3M}-\sigma
              -\left(\sigma-\frac{2}{3M}\right)z
              +O(\varepsilon)<0.
\]
The margins follow from the stated conditions on
\(M,J\), taking \(\tau\) and then \(\varepsilon\) small.
The dyadic series converges, giving exactly the (PA)
expansion over the restricted primitive quotients.

## A fixed-target estimate with a repaired target

Let \(r_{p,q}(N)\) count ordered positive \(p\)-full and
\(q\)-full summands of \(N\), with no coprimality condition.
Choose an auxiliary order \(r\) for the target and write
\(M=\max(p,q,r)\), \(\lambda_0=1/p+1/q+1/r\),
\(\beta=(1/3-\delta)\lambda_0\) with \(\delta\) as in
(GU), and \(B=\max\{p,q,1+\min(p,q)\}\).
Define the repair multiplier by
\[
 \mathcal R(N)=\prod_{\substack{\ell\mid N\\v_\ell(N)<r}}
                \ell^{\max(p,q,r-v_\ell(N))}.
\]
If the orders are sufficiently large for (GU2) and
\(M\beta<1\), then
\[
 r_{p,q}(N)\ll_{p,q,r,\varepsilon} N^\varepsilon
 \min\!\left\{N^{1/B},
    (N\mathcal R(N))^{1/(M+(1-M\beta)B)}\right\}.
 \tag{FT}
\]
For an \(r\)-full target the repair multiplier is \(1\).
In particular for \(p=q=r=s\) and sufficiently large \(s\)
the exponent at \(s\)-full targets is at most
\(1/(s+3\delta(s+1))\).

To prove (FT) fix \(g\mid N\) and put \(W=N/g\).
In the primitive quotient, \(Z\mid g\mathcal R(N)\).
Indeed \(K_\ell\) already fits in \(v_\ell(g)\)
when \(\ell\mid g\) and the target meets its required
valuation. When \(\ell\mid g\) but the target is
deficient there, the only unsatisfied entry is the
target's, with shortfall \(r-v_\ell(N)\). When
\(\ell
mid g\), a quotient deficiency must be on
the target because both original summands are full.
The excess is then at most
\(\max(p,q,r-v_\ell(N))\).
Summing (GU2) over the dyadic sizes of \(Z\) bounds
the count for this \(g\) by
\[
 (N\mathcal R(N))^{1/M}W^{\beta-1/M}N^\varepsilon.
\]

Write \(k=\min(p,q)\). On the primitive summands extract
the product of prime powers of each positive valuation
below \(k\). Each extracted factor divides \(g\):
a prime of deficient quotient valuation divides \(g\),
and \(v_\ell(g)\ge k\) there. The two remaining parts
are \(k\)-full, so the
[coefficient-uniform fixed-target graph bound](POWERFUL_SUM_GRAPHS.md#the-exponents-k-and-k1-suffice)
counts the choices for them at \(W^{1/(k+1)+\varepsilon}\)
for fixed extracted factors. The divisor bound sums
their choices. One can also extract the deficient
prime powers on the stronger summand below its own
order \(\max(p,q)\). These primes all divide \(g\),
and their possible exponents are bounded at fixed
orders, giving \(N^{o(1)}\) possible extracted
factors. Counting the stronger full part alone gives
\(W^{1/\max(p,q)}N^\varepsilon\).
Consequently the first bound for each \(g\) is
\(W^{1/B}N^\varepsilon\). Interpolate its minimum with
the decreasing bound above over \(1\le W\le N\), and
sum over the divisor-many \(g\). This gives (FT).

## An explicit pointwise saving for arbitrary full orders

For every sufficiently large fixed \(k\), let \(N\) be
a prescribed \(k\)-full target. The number of ordered
pairs of positive \(k\)-full summands, with arbitrary
common divisor, is
\[
 O_{k,\varepsilon}\bigl(N^{(1-1/307)/k+\varepsilon}\bigr).
 \tag{PT}
\]
This applies unchanged to higher summand fullness orders,
without any comparability assumption. More generally,
let \(D_k(N)\) be the product of the complete prime powers
of \(N\) whose valuation is below \(k\). If for a fixed
\(d_0<1/306\) one has \(D_k(N)\le N^{d_0}\), the same proof
gives the saving \((1-306d_0)/(307k)\) below \(1/k\).

Use the notation of the preceding fixed-\(g\) argument.
In the primitive target \(w=N/g\), let \(n_0\) be the
\(k\)-full part formed by taking all complete prime
powers of valuation at least \(k\). The remaining
factor divides \(gD_k(N)\): if the original target
valuation was already at least \(k\), a reduced
valuation below \(k\) forces a positive gcd valuation,
which is at least \(k\) from original summand fullness.
For the actual target defect primes the reduced part
divides \(D_k(N)\).

Put \(T_0=(1-306 d_0)/307\). If \(g>N^{T_0}\), the
coefficient-uniform pointwise graph gives
\[
 w^{1/(k+1)+\varepsilon}
 \le N^{(1-T_0)/k+\varepsilon}.
\]
Otherwise write \(t=\log_N g\le T_0\), so
\(n_0\ge N^{1-2t-d_0}\). It follows that
\[
 w\le n_0^{306/305},
\]
because \(307t+306d_0\le1\).
Extract the two quotient coefficients \(\kappa_i\) of
valuations below \(k\), as before fixed at divisor
cost. Use the canonical compression on the remaining
full factors to write the two quotient summands as
\[
       \kappa_1 R_1^k W_1,\qquad \kappa_2 R_2^k W_2,
\]
with \(R_i\) the product of the canonical roots.
If some \(\kappa_i W_i>n_0^{1/305}\), the corresponding
root product has upper bound \((w n_0^{-1/305})^{1/k}\).
Enumerate that side by its root-product tuples, with
divisor-many allocations per product and fixed
coefficient; the other quotient summand is determined.
This costs exponent at most
\[
 \frac1k\left(1-t-\frac{1-2t-d_0}{305}\right)
 \le \frac1k\left(1-\frac{1-d_0}{305}\right)
 \le\frac{1-T_0}{k}.
\]

Otherwise set \(T=R_1/R_2\) and
\(V=(\kappa_1W_1)/(\kappa_2W_2)\).
These are unit rationals modulo \(n_0\), with
\(VT^k\equiv-1\bmod n_0\),
\[
 H(T)\le n_0^{306/(305k)},\qquad H(V)\le n_0^{1/305}.
\]
Use the canonical adaptive rational bound with these
parameters, at exponent \((1-1/153)/k+\varepsilon\)
in base \(n_0\). This is at most the desired exponent
in base \(N\). Every exact ratio determines the two
primitive quotient values by their fixed sum \(w\).
Summing over the divisor-many \(g\) and coefficients
proves the pointwise bound. For a full target take
\(d_0=0\), obtaining (PT).

The [adaptive rational estimate](../../unpublished/powerful-number-representations/main.tex)
is the second range of the adaptive rational lemma in Section 9.4
(`lem:adaptive-rational` in the source). Its modulus is the
\(k\)-full core \(n_0\).

## A joint bound for height and valuation deficiency

Prescribe fullness orders \(s_i\ge2\) on a set \(I\) of at least
two of the three letters. Put \(M=\max_{i\in I}s_i\).
For the primitive quotient \(u_A+u_B=u_C=W\), write
\(h_i=\log u_i/\log W\). The sector with fixed parameters \(b,\rho\)
has
\[
 \operatorname{rad}\!\left(\prod_{i\notin I}u_i\right)\le W^b,
 \qquad
 b+\sum_{i\in I}h_i/s_i\le\rho .
 \tag{JD1}
\]
Take \(b=0\) if all three letters have prescribed orders. On the
prescribed letters define
\[
 \mathcal D=\{\ell:\ 0<v_\ell(u_i)<s_i
                          \text{ for some }i\in I\},\qquad
 D=\prod_{\ell\in\mathcal D}\ell^{1-v_\ell(u_i)/s_i},\qquad
 Z=\prod_{\ell\in\mathcal D}\ell^{\max_{j\in I}(s_j-v_\ell(u_j))}.
\]
In the product defining \(D\), \(i\) is the unique letter belonging to
\(\ell\). The real-valued \(D\) measures the positive fractional
valuation deficit. The integer \(Z\) is the least multiplier that
repairs every deficiency on the prescribed letters. The inequality
at the deficient primes in (GU1) gives \(D\le Z^{1/M}\).
Also \(\operatorname{rad}(u_Au_Bu_C)\le W^\rho D\) directly from
(JD1), writing \(\ell=\ell^{v_\ell(u_i)/s_i}\ell^{1-v_\ell(u_i)/s_i}\)
at each deficient prime.

There are absolute \(\rho_*,\delta>0\) such that, for fixed
\(0<\rho\le\rho_*\), the number of primitive quotients in (JD1) with
\(W\asymp Y\) and \(D\asymp H\) is
\[
 O_{\mathbf s,b,\rho,\varepsilon}\!
       \left(Y^{(1/3-\delta)\rho+\varepsilon} H^{3/4}\right).
 \tag{JD2}
\]
The unprescribed letter can be the largest term or either summand.
The constants are uniform in the dyadic bins.

Suppose that (JD2) fails for a sequence with \(\rho,\delta\to0\).
At each fixed parameter choice first pass to an unbounded sequence
of heights, taking a subsequence on which \(\log_Y H\to z\).
The elementary \(2\lambda/3\) bound at budgets tending down to
\(\rho+z\) gives
\[
 (1/3-\delta)\rho+\tfrac34 z\le\tfrac23(\rho+z),
 \qquad z\le(4+12\delta)\rho .
\]
Thus the radical budget tends to zero in the outer limit. For any
fixed \(\eta>0\), the limiting one-third count then gives
\((3/4-1/3-\eta)z\le(\eta+\delta)\rho\). Take \(\eta\downarrow0\)
after the parameter limit to get \(z=o(\rho)\).

Form anatomic boxes with losses \(o(\rho)\) and compression above
the prescribed orders. Also subdivide the complete values dyadically
so that each numerical height has a limit \(h_i^*\) on a further
subsequence. The root budget is \(L\le(1+o(1))\rho\), and the count
exponent is at least \((1/3-o(1))\rho\). The limiting one-third box
bound forces \(L/\rho\to1\) and \(d\to1/3\). The endpoint balance
gives limiting root mass \(1/3\) on every letter. The short-height
bound keeps each \(h_i^*\) positive.

For any \(0<\theta<1\) and \(i\in I\), the logarithms of the primes
on \(u_i\) with positive valuation below \((1-\theta)s_i\) sum to at
most \((\log D)/\theta\). Those low coordinates therefore have
normalized root mass \(o(1)\). Put \(t_i=s_i\rho\). By (JD1),
\(t_i\ge h_i\). The moment bound on each letter, applied first
outside the low coordinates, gives
\((1-\theta)t_i(L/\rho)(1/3-o(1))\le 1+o(1)\).
Consequently the \(t_i\) have finite positive subsequential limits
\(t_i^*\), and the weak measure on letter \(i\) is supported in
\([t_i^*,\infty)\). The first moment is at most \(h_i^*\), so
\(h_i^*/t_i^*\ge1/3\).

If a letter is unprescribed, its radical bound in (JD1) and the
letterwise radical comparison of the anatomic reduction give its
root mass at most \(b+o(\rho)\). Therefore \(\liminf b/\rho\ge1/3\)
when \(|I|=2\). With all three letters prescribed we have \(b=0\).
In both cases (JD1) yields
\[
 \sum_{i\in I} h_i^*/t_i^*\le |I|/3 .
\]
Every summand equals \(1/3\). Each prescribed measure has mass
\(1/3\), first moment \(t_i^*/3\), and support above \(t_i^*\).
Hence it is \((1/3)\delta_{t_i^*}\). Two such measures contradict
pairwise diffusion, proving (JD2).

### A global bound with a tail for the mandatory multiplier

For \(T\ge1\) restrict the original triples to the sector (JD1) with
\(Z\ge T\), allowing arbitrary common divisor and imposing fullness
only on \(I\). They are counted by at most
\[
 O_{\mathbf s,b,\rho,\varepsilon}\!
 \left(
     X^{\max\{1/M,\,(1/3-\delta)\rho\}+\varepsilon}
     T^{-1/(4M)}
 \right).
 \tag{JD3}
\]
At fixed quotient split the common divisor into its part \(f\)
supported on the quotient radical and its coprime part. The
coprime part is \(M\)-full. The Euler estimate (GU3) applies as
written, with the new \(Z\), because precisely the primes in
\(\mathcal D\) have mandatory positive minima. Thus the number
of multipliers is
\(O((X/W)^{1/M}Z^{-1/M}W^\varepsilon)\).
For \(W\asymp Y\) and \(Z\asymp K\), sum (JD2) over all
\(D\)-bins \(H\ll K^{1/M}\). The resulting count is bounded by
\[
 O\!\left(X^{1/M}Y^{(1/3-\delta)\rho-1/M+\varepsilon}
                      K^{-1/(4M)}\right),
\]
using the convergent geometric sum of the powers of \(H\).
The sum over \(K\gg T\) yields the displayed tail. Sum over
\(Y\le X\) and absorb its logarithmic factor into \(X^\varepsilon\).
For three prescribed orders \(\rho\) may be set to
\(\sum_i1/s_i\) when it is sufficiently small, giving (GU) with
a power tail in the required repair. It can be smaller when the
smaller fullness order is on a short primitive summand.

For example prescribe orders \(p,q\) on the summands and let the
primitive target have radical at most \(W^b\). If
\(b+\log(u_A)/(p\log W)+\log(u_B)/(q\log W)\le\rho\),
(JD3) uses \(\max(p,q)\) as its common-divisor order. The
unprescribed target's prime valuations enter only through \(b\).

### Asymptotic terms in the same sector

Write \(\beta=(1/3-\delta)\rho\). For any fixed \(J\ge0\) with
\(M\ge16(J+1)^2\) and \(\beta\le 1/(M+J+1)\), the count in
(JD1) has the expansion (PA) with its coefficients summed over
the sector and error \(O_{\mathbf s,b,\rho,J,\varepsilon}
(X^{1/(M+J+1)+\varepsilon})\).
Take
\[
 \max\{\beta,3/(4M),\,1/(M+J+1)\}<\sigma<1/(M+J)
\]
arbitrarily close to \(1/(M+J+1)\). The condition on \(M\)
puts \(3/(4M)\) strictly below this endpoint. The summability
estimate (GU3) reduces both the scalar expansion error and its
tail extension to \(X^\sigma\sum W^{-\sigma+\varepsilon}Z^{-\sigma}\).
Using (JD2) and summing \(D\ll K^{1/M}\) in a \(Z\)-bin
gives dyadic exponents \(\beta-\sigma+O(\varepsilon)<0\)
in \(Y\) and \(3/(4M)-\sigma<0\) in \(K\). This proves
absolute convergence and the claimed error. Whenever the
primitive sector is nonempty, the coefficients have strict
alternating signs by (PC).

For a concrete sector with only two fullness conditions,
take \(p=q=M\) on the summands and \(b=1/M,\rho=3/M\).
Every primitive quotient with \(\operatorname{rad}(W)\le W^{1/M}\)
satisfies (JD1). For fixed \(J\) and sufficiently large \(M\), its
arbitrary-common-divisor count has the stated expansion. The
sector contains \(1+(W-1)=W\) with \(W=2\cdot3^e\) for any
integer \(e\) large enough that \(W\ge6^M\). Its target has
valuation one at two and meets the radical bound.

### Separate tails for the prescribed letters

For each \(i\in I\), let \(D_i,Z_i\) be the products in \(D,Z\)
over the deficient primes on that letter. Fix
\(1/3<\gamma_i<1\) with \(\sum_{i\in I}(1-\gamma_i)<1\).
There are \(\delta,\rho_*>0\), depending only on these exponents,
for which the quotient bound (JD2) in independent bins \(D_i\asymp H_i\)
becomes
\[
 O_{\mathbf s,b,\rho,\boldsymbol\gamma,\varepsilon}
 \left(Y^{(1/3-\delta)\rho+\varepsilon}\prod_{i\in I}H_i^{\gamma_i}\right).
 \tag{JD4}
\]
Indeed fix the two remaining complete terms after omitting any \(j\in I\).
Their radical product is at most a constant times \(Y^\rho\prod_{i\ne j}H_i\).
Choose the radical product, then allocate its primes to the two letters,
then choose each term on its fixed radical support. The last two choices
cost \(Y^\varepsilon\) by the divisor bound and (CR). The omitted term
is determined by addition. In a failing sequence put \(z_i=\lim\log_Y H_i\)
and omit an index with largest \(z_i\). The resulting elementary count gives
\[
 \max_i z_i-\sum_{i\in I}(1-\gamma_i)z_i
           \le (2/3+\delta)\rho .
\]
The coefficient \(1-\sum_i(1-\gamma_i)>0\) gives
\(\max_i z_i=O(\rho)\). The limiting one-third count then gives
\(\sum_i(\gamma_i-1/3-\eta)z_i\le(\eta+\delta)\rho\) for fixed \(\eta>0\).
Take \(\eta\downarrow0\) after the parameter limit to obtain
\(\sum_i z_i=o(\rho)\). The valuation and moment argument proving (JD2)
now applies.

Put \(M_{-i}=\max_{j\in I\setminus\{i\}}s_j\), with value zero for
an empty index set, and put
\[
 \theta_i=\frac{s_i-1}{s_i\max(s_i-1,M_{-i})}.
\]
For a deficient prime on \(i\), the ratio of its exponent in \(D_i\)
to its exponent in \(Z_i\) is
\((s_i-e)/(s_i\max(s_i-e,M_{-i}))\), which is largest at \(e=1\).
Thus \(D_i\le Z_i^{\theta_i}\), with \(\theta_i\le1/M\).
The argument for (JD3) now gives the simultaneous tail
\[
 O_{\mathbf s,b,\rho,\boldsymbol\gamma,\varepsilon}
 \left(
 X^{\max\{1/M,(1/3-\delta)\rho\}+\varepsilon}
 \prod_{i\in I} T_i^{-(1/M-\gamma_i\theta_i)}
 \right)
 \qquad (Z_i\ge T_i\ge1).
 \tag{JD5}
\]
The expansion proof likewise uses \(\max_i\gamma_i\theta_i<1/(M+J+1)\)
in place of \(3/(4M)<1/(M+J+1)\). Taking one \(\gamma_i\) close to
\(1/3\) and the others close to one gives a stronger tail in the
chosen direction. With all \(\gamma_i\) equal to \(\gamma\), the
admissible range is \(\gamma>2/3\) for three prescribed letters
and \(\gamma>1/2\) for two.

The equal slopes \(2/3\) and \(1/2\), respectively, cannot replace
these strict inequalities for all deficiency bins with a positive
saving below \(\rho/3\). For three prescribed orders \(M\) and
\(\rho=3/M\), take all three quotient terms squarefree and comparable
to \(Y\). There are \(\gg Y^2\) primitive pairs of summands of this
kind. A squarefree sieve verifies this by choosing admissible
residues for finitely many primes, then bounding the larger-prime
loss by \(O(Y^2/K+Y^{3/2})\) for primes above \(K\), including
the coprimality condition. The local density products are bounded
below because the bad density at a prime is \(O(p^{-2})\).
In this family \(D\asymp Y^{3(1-1/M)}\), and the exponent
at \(\gamma=2/3\) would be \(2-(1+3\delta)/M+\varepsilon\).
For two prescribed summand orders \(M\), take \(W=2^j\) and
odd squarefree summands comparable to \(W\). There are \(\gg W\)
such ordered pairs by the same squarefree sieve, now excluding
at most two residues per odd \(p^2\). They are automatically
primitive. Take \(b=M^{-2}\) and \(\rho=2/M+b\), so (JD1) holds
for large \(j\), and \(D\asymp W^{2(1-1/M)}\). The exponent
at \(\gamma=1/2\) would be \(1-(1/3+2\delta)/M+O(M^{-2})+\varepsilon\).
Each argument may use finitely many dyadic bins because the
height ratios stay bounded above and below by fixed constants.
[Fixed-radical count](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

## One prescribed fullness order and captured scaled degrees

Prescribe one letter \(A\) of an ordered positive triple to be
\(M\)-full. Write \(u_A,u_U,u_V\) for the three letters of its
primitive quotient, \(W\) for the largest quotient letter, and
\(h_i=\log u_i/\log W\). Impose the combined sector conditions
\[
 \operatorname{rad}(u_Uu_V)\le W^b,\qquad b+h_A/M\le\rho .
 \tag{JW1}
\]
Fix \(a>0\) and a positive integer \(K\). Require, on one of the two
unprescribed letters \(T\), that there exist at most \(K\) intervals
in \([0,\infty)\), of total length at most \(\eta\), whose union \(J_T\)
satisfies
\[
 \sum_{\substack{\ell\mid u_T\\ \rho v_\ell(u_T)\in J_T}}\log\ell
      \ge a\rho\log W .
 \tag{JW2}
\]
The intervals and the choice of \(T\) may depend on the quotient.

There is an absolute \(\eta_0>0\) such that for any fixed
\(0<\eta\le\eta_0 a\), there are \(\delta,\rho_*>0\), depending
at most on \(a,K\), for which (JD2) holds in this sector with one
prescribed letter, for every fixed \(0<\rho\le\rho_*\).
Here
\[
  D=\prod_{0<v_\ell(u_A)<M}\ell^{1-v_\ell(u_A)/M},
  \qquad Z=\prod_{0<v_\ell(u_A)<M}\ell^{M-v_\ell(u_A)}=D^M .
\]
Consequently (JD3) holds for this sector, including its \(Z\)-tail,
and the same sector expansion (PA) holds under its conditions
\(M\ge16(J+1)^2\) and
\((1/3-\delta)\rho\le1/(M+J+1)\). The sign statement for nonempty
sectors is unchanged.

### The weighted quotient bound and the interval capture

Prove (JD2) by the same order of limits as before. At each
fixed parameter choice first take an unbounded sequence of
heights exhibiting a failure with positive exponent allowance,
and extract \(z=\lim\log_Y H\) for the \(D\)-bin. Then send
\(\rho,\delta\to0\) across the failing choices.
The inequalities of (JD2) use only
\(\operatorname{rad}(u_Au_Uu_V)\le W^\rho D\) and the
coefficient \(3/4>2/3\), so give \(z=o(\rho)\) and a
limiting \(d=1/3\) anatomic profile with \(L/\rho\to1\).
Take truncation loss \(o(\rho)\) and the compression index \(k\)
with \(k\ge M\) and \(k\rho\to\infty\).
Dyadic numerical splitting supplies height limits \(h_i^*\).

From (JW1) and the letterwise radical comparison, the combined
root mass on \(U,V\) is at most \(b+o(\rho)\). Balance forces
\(\liminf b/\rho\ge2/3\). Put \(t_A=M\rho\). The first-moment
argument in (JD2) gives a subsequence with \(0<t_A^*<\infty\)
and the scaled measure on \(A\) supported on
\([t_A^*,\infty)\), and \(h_A^*/t_A^*\ge1/3\).
Equation (JW1) gives the reverse inequality since
\(b/\rho+h_A/t_A\le1\). Thus
\(\nu_A=(1/3)\delta_{t_A^*}=(1/3)\delta_{3h_A^*}\).
The [height-adjusted density estimate](POWERFUL_TRIPLES_ABC_ENDPOINT_TRANSFER.md#uniform-diffuse-density-at-a-height-adjusted-atom)
applies to \(U,V\).

After pigeonholing which unprescribed letter satisfies (JW2),
discard captured primes of actual valuation at least \(k\)
from its captured record. Their total logarithm is at most
\(\log W/k=o(\rho\log W)\). All retained primes occur at their
exact squarefree anatomic positions. For each fixed budget and
compression, there are finitely many subsets of these positions.
Pigeonhole the retained set of captured positions before taking
the height limit. This costs a fixed factor, even when the
intervals in (JW2) depend on the quotient. Choose a representing
interval union for this set. After changing from scales
\(\rho e\) to \(Le\), its total length is at most
\(\eta L/\rho=\eta+o(1)\), and the captured normalized root
mass is at least \(a-o(1)\). Intersect the intervals with
\([0,R+1]\), losing vanishing mass by weak convergence and the
absolute compact support of the limiting measure. Pad with
empty intervals and pass to a subsequence on their finitely
many endpoints. The limiting closed intervals have total
length at most \(\eta\). The captured weak submeasure is
supported on their union, since any compact set at positive
distance from it eventually misses the selected coordinate
degrees. It is dominated by \(\nu_T\) and has
mass at least \(a\). This contradicts the density bound if
\(\eta\le a/(2C_0)\). This proves (JD2), and the calculations
for (JD3) and the sector expansion then apply verbatim.

For example, prescribe only fullness order \(M\) on the
target and take \(b=2/M\), \(\rho=3/M\). The primitive quotient
with summands \(2\cdot3^{2M}\) and \(1\) and their sum \(W\)
qualifies at sufficiently large \(M\). The combined summand
radical is \(6\le W^{2/M}\), and the exponent of \(3\) in
the first summand has scaled value \(6\) and supplies
\(\log3>(1/7)\rho\log W\). Thus use a single interval of
arbitrarily small positive width containing \(6\).
Taking \(g=Z\) repairs the target, and \(g\) is coprime
to the first primitive summand. The repaired first summand
still has valuation one at \(2\).


### Smaller deficiency exponents

In the one-prescribed-letter sector, fixing the two unprescribed
complete primitive terms fixes the third. Their radicals have
combined product at most \(Y^{b+o(1)}\) in a dyadic height bin.
Choosing a product radical has that many possibilities, allocating
its primes to the two terms costs \(Y^{o(1)}\), and there are
\(Y^{o(1)}\) complete terms on fixed radical supports by the fixed
positive exponent Rankin estimate. Thus
\(Q(Y,H)\ll_\varepsilon Y^{b+\varepsilon}\) uniformly in the
deficiency bin \(H\). Replace \(3/4\) in (JD2) by any fixed
\(\gamma>1/3\). The elementary bound first gives
\(z=O(\rho)\) because \(b\le\rho\). The limiting one-third
estimate then gives \(z=o(\rho)\) because \(\gamma>1/3\).
The remainder of the proof is unchanged.

For the narrower sector \(b\le(2/3+\tau)\rho\), a fixed
\(\gamma<1/3\) is also possible for sufficiently small fixed
\(\tau>0\). Otherwise choose failing parameters with
\(\gamma\uparrow1/3\) from below and \(\tau,\rho,\delta\to0\).
As before \(z=O(\rho)\). The failed lower estimate is now
\((\rho+z)/3-o(\rho)\). Hence \(L/(\rho+z)\to1\) and
the normalized count exponent tends to \(1/3\). The two
unprescribed letters have combined normalized mass \(2/3\)
by balance, and total unnormalized mass at most \(b+o(\rho)\).
The constraint on \(b\) thus forces \(z=o(\rho)\).
The height-adjusted one-atom and window-capture contradiction
then finishes the proof.

In either version the resulting repair tail in (JD3) is
\(T^{-(1-\gamma)/M}\), and the scalar expansion only needs
\(\gamma/M<1/(M+J+1)\) in place of \(3/(4M)<1/(M+J+1)\).

## Targets of small radical and two prescribed full summands

Let \(p,q\) be sufficiently large fixed fullness orders. Put
\[
 S=\max(p,q),\quad k=\min(p,q),\quad
 B=\max(S,k+1),\quad \lambda=1/p+1/q .
\]
Use the two-prescribed-letter case of (JD4) with
\(\gamma_A=\gamma_B=\gamma=3/5\). Shrink its absolute saving if
necessary so \(\sigma=1/3-\delta\ge19/60\). Write
\(\beta=\sigma\lambda\), \(a_0=\gamma/S\), and
\(\chi=\beta-a_0>0\). In particular \(\chi\ge1/(30S)\).
All constants below may depend on the fixed orders and the
positive exponent allowance. None depends on the prescribed
target or the gcd cutoff.

### Fixed targets, short summands and large common divisors

There is an absolute \(\nu_*>0\) such that for every \(N\ge2\)
with \(\operatorname{rad}(N)\le N^{\nu_*}\), the number of ordered
\(p\)-full and \(q\)-full summands \(A,B'\) of \(N\) satisfying
\(A\le N^\alpha\) and \(\gcd(A,B')\ge G\), where \(0\le\alpha\le1\)
and \(G\ge1\), is
\[
 \ll_{p,q,\varepsilon} N^\varepsilon
  \min\left\{N^{1/B},\
   N^{\sigma(\alpha/p+1/q)}\operatorname{rad}(N)^\sigma
                        G^{-\chi}\right\}.
 \tag{LR1}
\]
The lower cutoff on the orders is absolute. The estimate imposes
no fullness on the target, and the summand orders may have any
ratio. In particular a prescribed \(r\)-full target with
\(1/r\le\nu_*\) and no short-summand restriction has bound
\[
 N^{\min\{1/B,\;\sigma(1/p+1/q+1/r)\}+\varepsilon}.
\]
When \(p=q=s\) and \(r/s\) is large, this approaches
\((2/3-2\delta)/s\) in the exponent. This improves the graph
and the \(k\)-full-target explicit bounds in that regime.

Fix \(g\mid N\), let \(W=N/g\), and set
\(u=A/g,\ v=B'/g\). They form a primitive quotient.
For the two prescribed letters, let \(D_i,Z_i\) be their
fractional deficiency and mandatory repair factors from (JD4).
Since the original summands meet their prescribed orders,
\(Z_1Z_2\mid g\). This follows prime by prime: a deficient quotient
prime on one letter requires at least its missing valuation there
and the entire prescribed order on the other summand. Also
\(D_i\le Z_i^{1/S}\), so after summing the deficiency bins the
quotient penalty is at most \(g^{a_0}\) times a subpower of \(N\).
The part of \(g\) that can appear in \(Z_1Z_2\) is more precisely
\[
 g_{\rm rem}=\prod_{\ell^e\parallel N,\ v_\ell(g)=e}\ell^e.
\]
Indeed a prime in \(W\) is coprime to both primitive summands.
Thus \(g^{a_0}\) can be replaced by \(g_{\rm rem}^{a_0}\) in
the following estimate at large \(W\). The exact factors
\(D_i\le Z_i^{\theta_i}\) can also be used.

Take a small absolute \(b_*>0\) with
\(b_*+\lambda\) strictly below the budget cutoff in (JD4),
after choosing the orders sufficiently large, and take
\(4\nu_*\le b_*\). If \(W>N^{1/4}\), the target radical exponent
is \(b=\log_W\operatorname{rad}(W)\le b_*\). Factoring out the gcd
gives \(h_1=\log_W u\le\log_W(N^\alpha/g)\) and \(h_2\le1\).
Apply (JD4) with two prescriptions and sector budget arbitrarily
close from above to \(b+h_1/p+1/q\). The constants are uniform
as these exponents vary with \(N,g\): round \(b\) and the
allowable upper bound on \(h_1\) upwards on a finite grid whose
spacing is sufficiently small compared to the given exponent
allowance. The parameter choices are then finite, and there is
fixed slack to the absolute budget cutoff. Summing the \(D_i\)
bins, and using that the fixed \(W\) is part of its dyadic bin,
bounds the count at this gcd by
\[
 W^{\sigma(h_1^{\max}/p+1/q)}
      \operatorname{rad}(W)^\sigma g^{a_0}N^\varepsilon
 \le N^{\sigma(\alpha/p+1/q)}
      \operatorname{rad}(W)^\sigma g^{-\chi}N^\varepsilon .
\]
Here \(h_1^{\max}=\min(1,\log_W(N^\alpha/g))\).
If the latter is negative the gcd gives no solutions.

For \(W\le N^{1/4}\), the earlier coefficient-uniform fixed-gcd
bound gives \(N^\varepsilon W^{1/B}\le N^{1/(4B)+\varepsilon}\).
If there are solutions then \(g\le N^\alpha\), so
\[
 N^{\sigma(\alpha/p+1/q)}g^{-\chi}
 \ge N^{\sigma(1-\alpha)/q+\alpha a_0}
 \ge N^{\sigma/S}
 \ge N^{1/(4B)} .
\]
This proves the same fixed-gcd bound, even retaining
\(\operatorname{rad}(W)^\sigma\ge1\). Finally sum over the
divisor-many \(g\mid N\), use \(g\ge G\) and
\(\operatorname{rad}(W)\le\operatorname{rad}(N)\), and retain
the independent coefficient-uniform bound \(N^{1/B+\varepsilon}\).

### An aggregate count under an original-target radical cap

For any fixed \(0\le\nu\le1/S\), count ordered positive triples
with \(p\)-full and \(q\)-full summands as above, largest term
\(N\le X\), and \(\operatorname{rad}(N)\le X^\nu\).
Their number is
\[
 \ll_{p,q,\nu,\varepsilon}
   X^{\beta+\nu\max(\sigma,1-S\beta)+\varepsilon}.
 \tag{LR2}
\]
This includes the stronger per-target condition
\(\operatorname{rad}(N)\le N^\nu\).
For equal summand orders \(p=q=s\), its exponent is
\(2\sigma/s+(1-2\sigma)\nu\), strictly below \(1/s\) for
\(\nu<1/s\). It is also stronger than summing (LR1)
over all these targets using the radical count for their number.

Put \(\tau=\max(\sigma,1-S\beta)\) and
\(a=(1-\tau)/S\); thus \(0<a\le\beta\), with \(0<\tau<1\).
Fix a primitive quotient \(u+v=W\), write \(R=\operatorname{rad}(W)\),
\(P=\operatorname{rad}(uvW)\), and \(Z=Z_1Z_2\).
Split the common divisor \(g=f g_0\) uniquely into its \(P\)-supported
part \(f\) and the coprime part \(g_0\).
The latter is \(S\)-full, and
\(\operatorname{rad}(g_0)\operatorname{rad}(Wf)=\operatorname{rad}(N)\).
The elementary height count and fixed-radical count therefore bound
the number of such \(g_0\) by the smaller of
\((X/(Wf))^{1/S}\) and \(X^\nu/\operatorname{rad}(Wf)\),
up to a factor \(X^\varepsilon\). If either range is below one,
there are no choices. Interpolating their minimum gives
\[
 X^{a+\tau\nu+\varepsilon}W^{-a} f^{-a}
                          \operatorname{rad}(Wf)^{-\tau}.
\]
The deficient primes on summands do not divide \(W\), are
mandatory in \(f\), and appear in \(Z\) with exponent at most
\(S\). Hence
\(\operatorname{rad}(Wf)\ge R\operatorname{rad}(Z)\ge R Z^{1/S}\).
The Euler sum over allowed \(P\)-supported \(f\) is
\(\ll Z^{-a}W^\varepsilon\), by (GU3).
Since \(a+\tau/S=1/S\), the resulting multiplier count is
\[
 X^{a+\tau\nu+\varepsilon}W^{-a+\varepsilon}
                              R^{-\tau} Z^{-1/S}.
 \tag{LR3}
\]

Use dyadic height and radical bins \(W\asymp Y\), \(R\asymp H_R\),
and separate \(Z_i\) bins. In the low-radical branch
\(H_R\le Y^{b_*}\), apply (JD4) with the two summands prescribed,
using the finite-grid step above for \(b\), and with \(h_i\le1\).
The constants are uniform in radical bins. Summing independent
\(D_i\) bins up to the respective \(Z_i^{1/S}\) bounds gives
quotient exponent bound \(Y^{\beta+\varepsilon}H_R^\sigma
Z^{\gamma/S}\). Multiplying by (LR3), the \(Z_i\) bins sum
geometrically, their exponents being \(-(1-\gamma)/S\).
The \(R\) factor is harmless since \(\sigma-\tau\le0\).
Summing over \(Y\le X\) yields
\(X^{a+\tau\nu+\max(0,\beta-a)+\varepsilon}
 = X^{\beta+\tau\nu+\varepsilon}\).

In the high-radical branch \(H_R>Y^{b_*}\), use only the radical
bound for \(g_0\). Counting all \(P\)-supported integers \(f\)
up to \(X\) by the fixed-support Rankin estimate costs
\(X^\varepsilon W^\varepsilon\). This gives multiplier bound
\(X^{\nu+\varepsilon}W^\varepsilon R^{-1}Z^{-1/S}\).
For the quotient, fix the two terms of smallest radical. Their
radical product is at most \((W^\lambda RD)^{2/3}\). The same
fixed-radical estimate gives quotient bound
\(Y^{2\lambda/3+\varepsilon}H_R^{2/3}Z^{2/(3S)}\) after summing
the \(D_i\) bins. The \(Z_i\) exponents are again strictly
negative after multiplication. The remaining radical factor
\(H_R^{-1/3}\) and height factor \(Y^{2\lambda/3+\varepsilon}\)
give a summable height exponent, since one may fix \(b_*\)
and require \(4\lambda\le b_*\). Thus this branch contributes
at most \(X^{\nu+\varepsilon}\). It is covered by (LR2), because
\(\nu\le1/S\) and \(1-\tau\le S\beta\).

## A quarter coefficient at very small target radical

Prescribe orders \(p,q\) on two summands \(u,v\) of a primitive
quotient of height \(W\), and leave its target \(W\) unprescribed.
Write \(h_1=\log u/\log W,\ h_2=\log v/\log W\). Assume
\[
 R=\operatorname{rad}(W)\le W^b,\qquad
 b+h_1/p+h_2/q\le\rho .
\]
Let \(D=D_1D_2\) be the fractional deficiency factor defined
as in (JD2). For \(W\asymp Y,\ D\asymp H\), the elementary
bound
\[
 Q(Y,H)\ll_\varepsilon Y^{(\rho+b)/2+\varepsilon}H^{1/2}
 \tag{LQ2}
\]
counts the target and whichever prescribed summand has smaller
radical. The product of the two prescribed radicals is at most
\(W^{\rho-b}D\), and their chosen smaller radical is at most
its square root. Choose the two radicals and then choose
numbers of the fixed supports, at subpower cost.

For any fixed \(\gamma>1/2\) and \(\eta>0\), uniformly at
sufficiently small fixed \(\rho\), one has
\[
 Q(Y,H)\ll Y^{\rho/4+b/4+\eta\rho+\varepsilon}H^\gamma.
 \tag{LQ3}
\]
Otherwise take a sequence of failing fixed parameter choices
with \(\rho\to0\), taking the height limit first at each fixed
choice, and pass to \(z=\lim\log_Y H\). Formula (LQ2) first
forces \(z=O(\rho)\). Take anatomic boxes with compression
above both orders and losses \(o(\rho)\). Then
\(L\le\rho+z+o(\rho)\), \(L_C\le b+o(\rho)\), and the
count exponent is \(\gg\rho\), giving \(L\asymp\rho\).
The [profile inequality (LQ1)](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-quarter-inequality-for-a-letter-of-small-radical-mass)
gives
\[
 D_{\rm box}\le(\rho+z+b)/4+o(\rho),
\]
contrary to the failed lower bound in (LQ3).
The root cap \(L_C\) uses the letterwise radical comparison.

There is a strict improvement when \(b\) is sufficiently
small relative to \(\rho\). Fix \(\gamma>1/2\). There are
\(\tau,\delta,\rho_*>0\), with \(\delta<1/16\) if desired,
such that in the sector \(b\le\tau\rho,\ \rho\le\rho_*\),
\[
 Q(Y,H)\ll Y^{(1/4-\delta)\rho+\varepsilon}H^\gamma .
 \tag{LQ4}
\]
Otherwise choose failures with \(\tau,\delta,\rho\to0\).
Repeat the foregoing argument. Comparing with (LQ1) gives
\[
 (\gamma-\tfrac14)z\le \delta\rho+b/4+o(\rho),
\]
so \(z=o(\rho)\), \(L/\rho\to1\), \(d\to1/4\), and \(l_C=0\).
The cheap-letter argument for (LQ1) excludes a supported
summand because \(l_C<d\). Since the support threshold at \(d=1/4\)
is \(1/2\), this forces \(l_1=l_2=1/2\).

Write \(s_1=p,s_2=q\). At every fixed \(\theta>0\), the
positive valuations below \((1-\theta)s_i\) carry root
logarithm at most \((\log D)/\theta=o(\rho\log W)\).
They remain exact degree coordinates by the choice of
compression. The anatomic first moment on letter \(i\) is
at most its numerical height \(h_i\), since the initial
coefficient is a positive integer. This comparison can be
made with no additive outer error by using lower dyadic
root endpoints and upper complete endpoints in the inner
height limit. Hence, with \(t_i=s_i\rho\),
\[
 \frac{h_i}{t_i}\ge(1-\theta)(L/\rho)(l_i-o(1)).
\]
This bound remains valid if \(t_i\to0\). The sum of its two
left sides is at most \(1\) by the sector condition. Therefore
both tend to \(1/2\). On each positive-height summand, \(t_i\)
tends to \(2h_i^*>0\), its weak degree measure is supported
on \([2h_i^*,\infty)\), and its mass \(1/2\) and available
first moment \(h_i^*\) force the measure
\((1/2)\delta_{2h_i^*}\) with no escaped moment. It is thus
finite-scale. A possible zero-height summand is itself a
wide cofactor. The positive-height target \(C\) has modulus
cost at most \(l_C=0\). Under it (GW) releases both
summands, giving the contradiction to \(d=1/4\).
In the box comparison the order of limits is the same as in
(JD2): heights first at each fixed exponent allowance and
compression, and only then the small-budget subsequence.

### Fixed targets with a short-summand cap and gcd cutoff

Put \(\lambda=1/p+1/q\), \(S=\max(p,q)\),
\(k=\min(p,q)\), and \(B=\max(S,k+1)\).
For any fixed \(\gamma>1/2\), some \(\delta>0\) gives the
following bound, uniformly in \(N\ge2\), \(0\le\alpha\le1\),
and \(1\le G\le N\). The count is of ordered \(p\)-full
and \(q\)-full summands \(A,B'\) of \(N\) with the indicated
shortness and common-divisor restrictions.
\[
 r_{p,q}(N; A\le N^\alpha,\ \gcd(A,B')\ge G)
 \ll_{p,q,\gamma,\varepsilon} N^\varepsilon
 \min\left\{
   (N/G)^{1/B},
   N^{\sigma_*(\alpha/p+1/q)}
        \operatorname{rad}(N)^\gamma G^{-\sigma_*\lambda}
 \right\},
 \qquad \sigma_*=\tfrac14-\delta .
 \tag{LQ5}
\]
The orders are sufficiently large in terms of \(\gamma\).

Fix a possible gcd \(g\) and put \(W=N/g\), \(R=\operatorname{rad}(W)\),
\[
 A_0=\min\bigl(1,\log_W(N^\alpha/g)\bigr),\qquad
 \Lambda=A_0/p+1/q.
\]
If \(g>N^\alpha\), this gcd contributes zero. Otherwise
\(h_1\le A_0,\ h_2\le1\). Choose an internal deficiency exponent
\(1/2<\gamma_1<\gamma\), small \(c_0>0\) within the (LQ4)
sector width, and small absolute \(b_*>0\) below the budget
cutoffs. When \(b=\log_W R\le c_0\Lambda\), apply (LQ4)
with \(\rho=b+\Lambda\). After summing the \(D\)-bins this gives
\(W^{\sigma_*\Lambda+\varepsilon} R^{\sigma_*}
R_{\rm rem}^{\gamma_1}\), where \(R_{\rm rem}=\operatorname{rad}(N)/R\). Every deficient
prime lies in \(g\) by original fullness and is coprime to \(W\)
by primitivity, so the bound \(D\le R_{\rm rem}\) follows from
the fractional exponents in \(D\).

If \(c_0\Lambda<b\le b_*\), use (LQ3) instead, giving
\(W^{(1/4+\eta)\Lambda+\varepsilon}
R^{1/2+2\eta}R_{\rm rem}^{\gamma_1}\).
Choose \(\delta,\eta>0\) so small that
\(\delta+\eta<c_0(\gamma-1/2-2\eta)\) and shrink the
(LQ4) saving to that \(\delta\). The lower bound
\(R>W^{c_0\Lambda}\) absorbs the excess height exponent.
Both branches are therefore bounded by
\(W^{\sigma_*\Lambda+\varepsilon}\operatorname{rad}(N)^\gamma\).
Choose grid spacings for \(b,A_0\) sufficiently fine in the
fixed orders and exponent allowance, so the constants are
uniform, preserving the strict sector slack in the low
branch and the strict inequality in the medium branch.

Finally, if \(b>b_*\), the fixed-gcd graph bound
\(W^{1/B+\varepsilon}\) already fits this estimate, since
\(W^{1/B}\le R^\gamma\) for large enough orders. Hence
every branch satisfies the same bound. Since
\(W^{A_0}\le N^\alpha/g\),
\[
 W^{\sigma_*\Lambda}
 \le N^{\sigma_*(\alpha/p+1/q)} g^{-\sigma_*\lambda}.
\]
Sum over the divisor-many admissible \(g\), using \(g\ge G\).

### Aggregate counterpart at very small original target radical

Fix \(\gamma=3/5\) and let \(\tau,\delta\) be as in (LQ4).
Set \(c_0=\min(1/8,\tau/4)\) and decrease
\(\delta\) to at most \(c_0/16\). Put
\[
  \sigma_*=\tfrac14-\delta,\qquad
  \beta_*=\sigma_*\lambda,\qquad
  t_*=\max(\sigma_*,1-S\beta_*).
\]
For sufficiently large fixed orders and
\(0\le\nu\le\min(1/S,c_0\lambda/16)\), the aggregate count over
\(N\le X\) with \(\operatorname{rad}(N)\le X^\nu\) is
\[
 \ll_{p,q,\nu,\varepsilon} X^{\beta_*+t_*\nu+\varepsilon}.
 \tag{LQ6}
\]

Use the multiplier count (LR3) from the low-radical proof,
with \(a=(1-t_*)/S\). It gives, at any one quotient,
\(X^{a+t_*\nu+\varepsilon}W^{-a+\varepsilon}
R^{-t_*}Z^{-1/S}\), where \(R=\operatorname{rad}(W)\)
and \(Z\) is the mandatory summand repair. The same proof
with \(t=1\) gives the pure radical multiplier count
\(X^{\nu+\varepsilon}W^\varepsilon R^{-1}Z^{-1/S}\).
Indeed in that case count the \(P\)-supported \(f\) at
subpower cost instead of using its weighted Euler sum.

For dyadic \(W\asymp Y\), if \(R\le Y^{c_0\lambda}\)
up to the dyadic bin convention, use (LQ4) with
\(b=\log_Y R+o(1)\) rounded upwards on a finite grid and
\(\rho=\lambda+b\). The defining inequality \(b\le\tau\rho\)
has fixed slack. After summing the \(D\)-bins, the quotient
bound is \(Y^{\beta_*+\varepsilon}R^{\sigma_*} Z^{\gamma/S}\).
The weighted multiplier count gives (LQ6) in this branch
because \(a\le\beta_*,\ t_*\ge\sigma_*\), and the sum in the
\(Z\)-bins is geometric.

For \(Y^{c_0\lambda}<R\le Y^{b_*}\), where \(b_*>0\) is a
small absolute threshold, use (LQ3) with sufficiently small
fixed \(\eta\le c_0/16\) and the pure radical multiplier count.
After summing the \(Z\)-bins, the contribution has height
exponent at most
\[
 \nu+\left(\tfrac14+\eta
              -c_0(\tfrac12-2\eta)\right)\lambda .
\]
Here \(1/2+2\eta<\gamma\) can be ensured by reducing \(\eta\).
This exponent is at most \(\beta_*+t_*\nu\), since
\((1-t_*)\nu\le c_0\lambda/16\), \(\delta,\eta\le c_0/16\),
and \(c_0(1/2-2\eta)>3c_0/16\).
For \(R>Y^{b_*}\), use the elementary two-smallest-radicals
quotient estimate and the pure radical multiplier bound as
in the high-radical branch of (LR2). This gives at most
\(X^{\nu+\varepsilon}\), which is covered by (LQ6) because
\(\nu\le1/S\) and \(1-t_*\le S\beta_*\).
Dyadic cutoff factors and finite-grid approximations
are absorbed into the desired exponent allowance.

In particular if \(p=q=s\), the exponent in (LQ6) is
\((1/2-2\delta)/s+(1/2+2\delta)\nu\). For example it is
strictly below \(1/(2s)\) when \(0<\nu\le\delta/s\),
which satisfies the imposed range for the chosen \(\delta\).


With the additional condition \(A\le N^\alpha\), the same
proof of (LR2) and (LQ6) applies with \(\lambda\) replaced by
\(\alpha/p+1/q\) in their exponents and stated ranges.
After division by \(g\), one has \(u\le W^\alpha\), so all
three quotient branches use this smaller height budget.
The multiplier bounds are unchanged.

## One full summand and a small target radical

The quarter saving extends to one prescribed full summand when the
other summand has a radical bound. Write \(\sigma=1/4-\delta\).
There are absolute positive \(\delta,c_*\) such that for sufficiently
large fixed \(p\), the number of representations \(A+B=N\) with
\(A\) \(p\)-full, \(\operatorname{rad}(B)\le N^{1/p}\), and
\(\operatorname{rad}(N)\le N^{c_*/p}\), is
\(O_{p,\varepsilon}(N^{(1/2-\delta)/p+\varepsilon})\) with
arbitrary common divisor. The aggregate count with \(N\le X\),
\(\operatorname{rad}(B)\le X^{1/p}\), and
\(\operatorname{rad}(N)\le X^\nu\) is
\(O_{p,\nu,\varepsilon}(X^{2\sigma/p+(1-2\sigma)\nu+\varepsilon})\)
for \(0\le\nu\le c_*/p\).

### The primitive quotient with a fractional repair cost

For \(u+v=W\) a primitive quotient, prescribe fullness
order \(p\) only on \(u\). Allow the other summand \(v\)
arbitrary valuations. Put \(h=\log u/\log W\), and
assume
\[
 \operatorname{rad}(v)\le W^{b_B},\qquad
 \operatorname{rad}(W)\le W^{b_C},\qquad
   h/p+b_B+b_C\le\rho .
\]
Let \(D=\prod_{0<v_\ell(u)<p}\ell^{1-v_\ell(u)/p}\),
and let \(Z=D^p\) be the mandatory integer multiplier.
Fix any \(\gamma>1/4\). There are \(\tau,\delta,\rho_*>0\)
such that, for fixed \(0<\rho\le\rho_*\) and \(b_C\le\tau\rho\),
the quotient count in \(W\asymp Y,D\asymp H\) is
\[
 \ll_{p,b_B,b_C,\rho,\varepsilon}
      Y^{(1/4-\delta)\rho+\varepsilon}H^\gamma.
 \tag{QF2}
\]
The implied constants do not depend on the bins.

In a failing sequence with \(\rho,\delta,\tau\to0\), fix
each parameter choice, take its unbounded height sequence
and pass to \(z=\lim\log_Y H\) before changing parameters.
Counting the two unprescribed complete terms by their
radicals gives \(Q\ll Y^{b_B+b_C+\varepsilon}\), so
\(z=O(\rho)\). Form anatomic boxes with losses \(o(\rho)\)
and compression greater than \(p\), and split numerical
heights before passing to the limiting profile.
The complete radical bound is
\(\operatorname{rad}(uvW)\le D u^{1/p}\operatorname{rad}(vW)\).
Thus \(L\le\rho+z+o(\rho)\), \(L_C\le b_C+o(\rho)\),
and (LQ1) gives \(D_{\rm box}\le(\rho+z+b_C)/4+o(\rho)\).
The failed estimate and \(\gamma>1/4\) imply \(z=o(\rho)\),
\(L/\rho\to1\), \(d\to1/4\), and \(l_C=0\). Thus each
summand has normalized root mass \(1/2\).

The unprescribed summand root bound gives
\(\liminf b_B/\rho\ge1/2\). On \(u\) the logarithmic mass
of the primes of degree below \((1-\theta)p\) is at most
\((\log D)/\theta\), for fixed \(\theta>0\).
These are exact original degree positions. The root
moment at each fixed parameter choice is at most its
actual numerical height, using lower dyadic root
endpoints in the height limit. Hence
\[
 h/(p\rho)\ \ge\ (1-\theta)(L/\rho)(1/2-o(1)).
\]
Since \(b_B/\rho+h/(p\rho)\le1\), both ratios tend to \(1/2\).
If the limiting numerical height \(h^*>0\), then
\(p\rho\to2h^*\), the limiting degree measure is supported
on \([2h^*,\infty)\), and its half mass and moment at most
\(h^*\) make it \((1/2)\delta_{2h^*}\) with no escaped
moment. Thus \(u\) is finite-scale. If \(h^*=0\), then
\(u\) is numerically short and is releasable anyway.
Either outcome contradicts (LQ8) with \(l_C=0\). This proves (QF2).
[At \(l_C=0\), both opposite letters must be diffuse](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-global-inequality-for-a-low-mass-letter-and-diffuse-degrees)

Without the condition \(b_C\le\tau\rho\), the same
elementary argument, (LQ1), and the height limit give, for
every fixed \(\eta>0\) at sufficiently small fixed \(\rho\),
\[
 Q(Y,H)\ll
       Y^{\rho/4+b_C/4+\eta\rho+\varepsilon}H^\gamma .
 \tag{QF3}
\]
Indeed \(z=O(\rho)\) as above, and the box upper exponent
already contradicts any failure since \(\gamma>1/4\).

### Fixed-target count with arbitrary common divisor

For sufficiently large fixed \(p\), let \(N\) be a target
with \(\operatorname{rad}(N)\le N^{c/p}\), with \(c>0\)
small and fixed. Count \(A+B=N\) with \(A\) \(p\)-full,
\(\operatorname{rad}(B)\le N^{1/p}\), and arbitrary gcd.
Then
\[
 r(N)\ \ll_{p,\varepsilon}
 N^{2\sigma/p+\varepsilon}\operatorname{rad}(N)^\gamma,
 \qquad \sigma=1/4-\delta .
 \tag{QF4}
\]
In particular decreasing \(c\) so that \(\gamma c\le\delta\)
gives fixed-target exponent at most \((1/2-\delta)/p\).

For \(g=\gcd(A,B)\mid N\), put \(W=N/g\). If
\(W>N^{1/8}\), the quotient has
\(b_B\le\log_W N^{1/p}\le8/p\) and
\(b_C=\log_W\operatorname{rad}(W)\le8c/p\).
Use \(\rho=1/p+b_B+b_C\), rounding the upper bounds upward
on an arbitrarily fine finite grid so that the implied
constants are uniform in \(N,g\). Choose \(c<\tau/16\)
to preserve the sector margin. Every deficient prime of
\(A/g\) divides \(g\), since \(A\) is full, and it is
coprime to \(W\). Thus \(D\le R_{\rm rem}\), where
\(R_{\rm rem}=\operatorname{rad}(N)/\operatorname{rad}(W)\).
Summing the deficiency bins in (QF2) gives
\[
 W^{\sigma/p}N^{\sigma/p}
       \operatorname{rad}(W)^\sigma R_{\rm rem}^\gamma N^\varepsilon
 \le N^{2\sigma/p+\varepsilon}\operatorname{rad}(N)^\gamma .
\]
For \(W\le N^{1/8}\), remove the deficient prime powers of
\(A/g\). Their supports divide \(g\) and the allowable
exponents are bounded by \(p\), so the coefficient choices
cost at most \(N^\varepsilon\) by the fixed-support bound. Count the remaining full part in
\(W^{1/p+\varepsilon}\) ways. This lower-height branch fits
(QF4), and summing over the divisors \(g\) proves the bound.

### Aggregate count with both original radical caps

For sufficiently large fixed \(p\), count \(A+B=N\le X\)
with only \(A\) prescribed \(p\)-full, and
\(\operatorname{rad}(B)\le X^{1/p}\),
\(\operatorname{rad}(N)\le X^\nu\).
For some sufficiently small fixed \(c>0\), and
\(0\le\nu\le c/p\), the count is
\[
 \ll_{p,\nu,\varepsilon}X^{2\sigma/p+(1-2\sigma)\nu+\varepsilon}.
 \tag{QF5}
\]
For \(\nu\le\delta/p\) as well, this exponent is strictly
below \(1/(2p)\).

To verify the multiplier bound, at fixed quotient \(u,v,W\)
put \(R_B=\operatorname{rad}(v)\), \(R_C=\operatorname{rad}(W)\)
and \(P=\operatorname{rad}(uvW)\). Split \(g=f g_0\), where
\(f\) is \(P\)-supported and \(g_0\) is coprime to \(P\).
Then \(g_0\) is \(p\)-full. The deficient primes in \(Z\)
all divide \(f\), and none divides \(vW\), so they occur in
both original radicals. The height count and fixed-radical
count bound the number of possible \(g_0\) by the minimum
of
\[
 (X/(Wf))^{1/p},\qquad
 \frac{X^{1/p}}{R_B\operatorname{rad}(Z)},\qquad
 \frac{X^\nu}{R_C\operatorname{rad}(Z)}
\]
up to \(X^\varepsilon\). The count is zero if any range
is below one. For each fixed radical of \(g_0\), the integer
multiplicity is \(O_\varepsilon(X^\varepsilon)\).
The allowed support factors of \(f\) have total \(f^{-a}\) sum
\(O_{a,\varepsilon}(X^\varepsilon)\) for each fixed \(a>0\), by
the product \(\prod_{\ell\mid P}(1-\ell^{-a})^{-1}\).
Use weights \(\sigma,\sigma,1-2\sigma\) on the three multiplier bounds.
Summing the allowed \(f\) by the Euler bound
at exponent \(a=\sigma/p\), and using
\(\operatorname{rad}(Z)\ge Z^{1/p}\), gives the multiplier
count
\[
 X^{2\sigma/p+(1-2\sigma)\nu+\varepsilon}
 W^{-\sigma/p+\varepsilon}
       R_B^{-\sigma}R_C^{-(1-2\sigma)} Z^{-1/p}.
\]
Using only the last radical count and counting \(f\) at
subpower cost also gives
\(X^{\nu+\varepsilon}R_C^{-1}Z^{-1/p}\).
The caps here are on the original values. In particular, the
primes of \(Z\) enter both radicals independently of \(g_0\).

Take dyadic bins \(W\asymp Y,R_B,R_C,Z\). As usual, set
\(b_i=\log_Y R_i\) to within finite-grid exponent loss.
Fix a small absolute \(b_*>0\) in the allowed budget range.
When both \(b_i\le b_*\) and \(b_C\le c_0(1/p+b_B)\), choose
\(c_0<\tau/4\) and use (QF2) with \(\rho=1/p+b_B+b_C\).
Its quotient estimate
\(Y^{\sigma/p}R_B^\sigma R_C^\sigma Z^{\gamma/p}Y^\varepsilon\)
cancels the height and \(B\)-radical terms. The remaining
target-radical exponent is negative since \(\sigma<1/4\),
and the repair bins sum geometrically if \(1/4<\gamma<1\).
The other bins contribute logarithms.

When both \(b_i\le b_*\) but \(b_C>c_0(1/p+b_B)\), use
(QF3) and the pure target-radical multiplier. The height
and radical part of their product, beyond \(X^\nu\), is
at most
\[
 (Y^{1/p}R_B)^{\,1/4+\eta}
 R_C^{-(1/2-2\eta)}
 \le (Y^{1/p}R_B)^{\,1/4+\eta-c_0(1/2-2\eta)}.
\]
The repair bins still sum geometrically. With \(\delta,\eta\)
small compared with \(c_0\), the displayed exponent is a
positive number strictly below \(\sigma\) by a fixed margin.
Since \(Y\le X\), \(R_B\le X^{1/p}\), and \(\nu\le c/p\)
with \(c\) sufficiently small, this branch is covered by
(QF5). One may take \(c_0\le1/8\) and
\(\delta,\eta\le c_0/16\), then \(c\le c_0/16\).

Finally if either \(b_B\) or \(b_C\) exceeds \(b_*\),
then \(Y\le X^{1/(pb_*)}\) since both original caps are at
most \(X^{1/p}\). The elementary quotient estimate fixing
the target and the prescribed summand gives
\(Y^{1/p+\varepsilon}R_C D\). The pure target-radical
multiplier cancels \(R_C D\) since \(D=Z^{1/p}\).
The branch contributes at most
\(X^{\nu+1/(p^2b_*)+\varepsilon}\), smaller than the exponent in (QF5)
for large \(p\) and sufficiently small fixed \(c\).
Every use of a variable \(b_i\) is by a finite upper grid
with fixed slack; dyadic factors and the logarithmic
number of bins are absorbed into the exponent allowance.

### Scalar asymptotics in the one-full primitive quotient sector

Impose only \(p\)-fullness on the original first summand. On its
primitive quotient \(u+v=W\), impose the conditions of (QF2),
namely
\[
 \operatorname{rad}(v)\le W^{b_B},\qquad
 \operatorname{rad}(W)\le W^{b_C},\qquad
 (\log_W u)/p+b_B+b_C\le\rho,\qquad b_C\le\tau\rho .
\]
Use fixed \(1/4<\gamma<1\), and let \(\sigma_0=1/4-\delta\) be the
constant of (QF2) for this \(\gamma\). For fixed \(p,J\) with
\[
 p\ge16(J+1)^2,\qquad
 \max\{\sigma_0\rho,\gamma/p\}<1/(p+J+1),
\]
and \(\rho\) in the allowed range, the original count with height
at most \(X\) and arbitrary common divisor is
\[
 \sum_{j=0}^{J} C_j X^{1/(p+j)}
 +O_{p,b_B,b_C,\rho,J,\varepsilon}
       (X^{1/(p+J+1)+\varepsilon}),
 \qquad (-1)^j C_j>0
 \tag{QF6}
\]
if the primitive sector is nonempty. The coefficients are sums
of the scalar coefficients from (PC) over the indicated quotients
and their allowed support factors. The radical caps here are on
the primitive quotient. The original target and second summand
may have large radicals.

For a fixed quotient write \(P=\operatorname{rad}(uvW)\) and
\(g=f g_0\) as in the preceding aggregate count. Thus \(f\) is
\(P\)-supported, \(g_0\) is coprime to \(P\), and \(g_0\) is
\(p\)-full. The allowable \(f\) repair \(u\), so they are divisible
by \(Z=D^p\). Conversely such a supported \(f\) with \(p\)-full
\(fu\), together with any coprime \(p\)-full \(g_0\), gives an
original triple from this quotient. Each original triple determines
its quotient and both factors of its gcd uniquely. For fixed
\(a>0\), (GU3) gives
\[
 \sum_{f\text{ allowed}} f^{-a}\ll_{a,\varepsilon} W^\varepsilon Z^{-a}.
\]
Indeed it suffices to sum over \(f=Z f'\), with \(f'\) supported on
\(P\), and apply the fixed-support Euler product.

Use (PC) with \(M=p\) to count \(g_0\) up to \(X/(Wf)\). Choose
\(a\) just above \(1/(p+J+1)\) and below \(1/(p+J)\). The hypotheses
ensure that
\(a>\max\{\sigma_0\rho,\gamma/p,3/(4p)\}\).
The scalar expansion error, and the tails from extending each
coefficient sum past \(Wf=X\), are bounded by \(X^a\) times
\(\sum_{\text{quotients}} W^{-a+\varepsilon} Z^{-a}\).
Count a dyadic box \(W\asymp Y,Z\asymp H\) by (QF2), at cost
\(O(Y^{\sigma_0\rho+\varepsilon}H^{\gamma/p})\). Both dyadic
sums converge, since \(a>\sigma_0\rho\) and \(a>\gamma/p\).
Taking \(a\) sufficiently close to the endpoint gives the error
in (QF6). The same convergent majorant proves absolute convergence
of each coefficient sum, because each \(1/(p+j)>a\).
All terms of each coefficient sum have the same strict sign by
(PC), giving the alternating signs.

For example, take \(b_B=1/p\), \(b_C=c_*/p\),
\(\rho=(2+c_*)/p\) and \(\gamma=1/2\), with \(c_*>0\)
sufficiently small. This works for every fixed \(J\) and all
sufficiently large \(p\). Set \(W=4^t\), \(v=1\), \(u=W-1\)
with \(3\nmid t\) and \(t\) large enough that \(2\le W^{c_*/p}\).
Since \(v_3(4^t-1)=1\), the mandatory repair \(g=Z\) has
\(v_3(g)=p-1\). Hence the two original letters \(gv\) and \(gW\)
both have valuation \(p-1\) at three, while \(gu\) is \(p\)-full.

## An exact-power summand at a full fixed target

For every \(c_0,\eta>0\) there is an absolute-in-the-orders cutoff
\(\lambda_*(c_0,\eta)>0\) with the following property. Fix positive
integers \(p,q,r\), put \(\lambda=1/p+1/q\), and assume
\(\lambda\le\lambda_*\) and \(r\lambda\ge c_0\).
For each prescribed \(r\)-full \(N\), count ordered positive
solutions
\[
 A+B=N,\qquad A=x^p,\qquad B\text{ is }q\text{-full},
\]
with arbitrary common divisor. For every \(\varepsilon>0\)
this count is
\[
 O_{c_0,\eta,p,q,r,\varepsilon}
     \left(N^{\eta\lambda+\varepsilon}\right).
 \tag{EP}
\]
This permits unequal \(p,q\), with the same absolute normalized
cutoff, and an \(r\) at least a fixed multiple of their harmonic scale. For a nonsymmetric admissible representation, take \(g=6^{pJ}\) with \(pJ\ge\max(q,r)\),
and set \(A=g,\ B=2g,\ N=3g\).

### Fixing the gcd and the complete gcd-supported parts

Take \(0<\eta\le1\), set \(\tau=\eta/8\), and choose
\(0<\kappa\le c_0\tau/8\). These are fixed before the
small-order limit. Fix \(g=\gcd(A,B)\mid N\) and write \(W=N/g\).
Freeze the complete \(\operatorname{rad}(g)\)-supported parts
\(A_s,B_s\) of both original summands. There are \(N^\varepsilon\)
possibilities, uniformly in \(g\mid N\), by counting integers
at most \(N\) on a prescribed radical support. For example
the fixed-support Rankin product with a sufficiently small
positive exponent proves this bound. Both \(A_s,B_s\) are
multiples of \(g\). With fixed positive integers
\(A_0=A_s/g,\ B_0=B_s/g\), the primitive quotient equation
becomes
\[
 A_0y^p+B_0 V=W,
 \qquad V\text{ is }q\text{-full}.
\]
Neither \(y\) nor \(V\) has a prime on the support of \(g\).
In particular \(A_0y^p\) and \(B_0V\) are coprime, and both
are coprime to \(W\). The coefficients \(A_0,B_0\) enter the counts through unit
residues and the fixed equation.

If \(W\le N^\tau\), count either \(y\le W^{1/p}\) or the
\(q\)-full integer \(V\le W\), since either determines
the solution for fixed coefficients. Their minimum is at
most \(W^{\lambda+\varepsilon}\). Thus this branch and the
divisor-many outer choices fit (EP).

Suppose \(W>N^\tau\). Put \(k=\lfloor\kappa/\lambda\rfloor\)
and let \(n\mid W\) consist of all complete prime powers
in \(W\) of valuation at least \(k\). If a prime power in
\(W/n\) has original target exponent \(e_N\) and reduced
exponent \(0<e_W<k\), then its exponent in \(g\) is at least
\(r-k\). Consequently
\[
 \log(W/n)\le\frac{k}{r-k}\log g,\qquad
 \frac{k}{r-k}\le\frac{\kappa}{c_0-\kappa}<\tau/2.
\]
Thus \(n>N^{\tau/2}\) is a fixed \(k\)-full modulus determined by
\(N\) and \(g\).

### Boxing the second summand

At fixed \(p,q,r\) take \(K=\lceil\max(q,\lambda^{-2})\rceil\).
Write \(V\) in its exact degrees \(q,\ldots,K^2\), with one
additional degree-\(K\) root combining all primes of actual
valuation \(e>K^2\). For such a prime write
\(e=K\lfloor e/K\rfloor+(e\bmod K)\). The product of the
remainders has logarithm at most
\[
 K\sum_{v_\ell(V)>K^2}\log\ell\le\log N/K .
\]
Every degree-root base is at most \(N^{1/q}\), and the sum
of all logarithms of these roots is at most \((\log N)/q\),
because all represented positive degrees are at least \(q\).
The squarefree exact degree slots and the extra packed
slot may be kept separate even if the packed slot degree
also occurs among the exact ones. The number of slots is
finite at each fixed order choice. Dyadically box these roots
at a power-of-\(\log N\) cost depending only on that choice.
All coefficient, divisor, and boxing losses can be made
arbitrarily small compared with \(\lambda\) in the
\(N\)-exponent at each fixed order choice.

We justify the remaining uniform assertion by contradiction.
If (EP) failed for arbitrarily small order budgets, choose
fixed order triples \((p_j,q_j,r_j)\) with \(\lambda_j\to0\),
and for each such choice an unbounded sequence of targets
witnessing failure with a positive exponent allowance. Take
the target-height limit first for each order triple,
selecting over the divisor-many \(g\), the two supported
coefficients, and the fixed-order root boxes. Thus a selected
box can be assumed to contain solutions of count exponent
at least \((\eta/2)\lambda_j\) in \(N\). At each fixed \(j\)
pass to a further subsequence so that the dyadic root height
exponents have limits. All normalized root masses and
numerical weights below refer to these limits before \(j\)
tends to infinity.

Set \(\xi=\eta/32\). There are at most \(1/\xi\) root slots
of \(V\) of normalized mass at least \(\xi\), since the
sum of the normalized masses is at most \(1/(q\lambda)\le1\).
The complement of those slots has numerical weight tending
to zero. Otherwise along a subsequence its weight tends to
some \(w>0\). Keep just its slots of degree at least
\(w/(4\lambda)\); the discarded part has weight at most
\(w/4+o(1)\). Greedy grouping of the remaining slots into
blocks of mass at most \(2\xi+o(1)\) gives boundedly many
blocks, since each completed block except the last costs
at least \(\xi-o(1)\). One block has positive lower
numerical weight \(w'>0\) on a subsequence. Freeze the
values of the roots in that block at exponent cost at most
\((2\xi+o(1))\lambda\). Their full product \(m\) divides
\(V\), is a modulus of size at least \(N^{w'/2}\), and
\(A_0,W\) are units modulo \(m\) by quotient primitivity.
Then \(y^p\equiv W/A_0\bmod m\). The \(p\)-th power map on
units has kernel of size \(N^{o(1)}\) at each fixed \(p\),
by CRT and the unit groups of prime powers. Since
\(y\le N^{1/p}\le N^\lambda<m\) for small enough \(\lambda\),
this congruence gives at most \(N^\varepsilon\) values of
\(y\), each determining the solution. This contradicts
the selected exponent \((\eta/2)\lambda\), since
\(2\xi=\eta/16\).

### Counting under the full target core

It follows that the selected boxes express \(V\) by boundedly
many roots of height \(N^{O(\lambda)}\) and a cofactor
of height \(N^{o(1)}\); the latter includes the packing
remainder. Apply the amplified bounded-root modulus
estimate (GW) under the fixed \(n\):
\[
 \frac{y^p}{V}\equiv-\frac{B_0}{A_0}\pmod n .
\]
All variables and coefficients are units by primitivity.
One exact ratio \(y^p/V\) determines the positive primitive
pair \(A_0y^p,B_0V\), or equivalently determines it by the
fixed sum \(W\). Select one representation per exact ratio.
The number of short roots is bounded by \(m_*=2+\lceil
1/\xi\rceil\), independently of \(j\). Pad with units to keep the dimension fixed. Choose the amplification order \(\ell\) sufficiently
large in terms of \(\kappa,m_*,\eta\), before letting
\(\lambda_j\to0\), so that
\[
 2(m_*^2+m_*)/(\kappa\ell)<\eta/4.
\]
The relative narrow height constants in (GW) are bounded
in terms of \(\ell,\kappa,\tau,m_*\), using
\(\log_N n\ge\tau/2\). The amplified wider relative
height tends to zero along the chosen order sequence.
The exponents of the short variables are fixed in each
inner box and the required lower cutoff on \(k\) is
independent of them. Thus (GW), after amplification,
gives normalized exponent at most
\(2(m_*^2+m_*)/(\kappa\ell)\) up to arbitrarily small
inner height allowances. This contradicts the selected
\(\eta/2\) lower exponent and proves (EP).

For this amplification, colour primes
by their sign in the exact rational output, retain a
subpower fraction so no cancellation can occur in
products of outputs, and choose one actual coordinatewise
tuple product for each exact output product. The divisor
bound gives at most a subpower of factorizations of any
one such product, since its numerator and denominator
each split into fixed \(\ell\) parts. Apply the
unamplified (GW) bound to the selected amplified tuples
with exact-fibre bound one and target raised to the
\(\ell\)-th power. This is valid here since all exact
outputs have polynomial height and one output
determines the solution.

### Full fixed targets with an approximate valuation basis

[The bounded-basis valuation condition (AB1)](ABC_WIDE_COFACTOR_CONGRUENCES.md#adaptive-finite-basis-approximation-of-prime-valuations)
also works when one releases the first summand modulo an opposite full block.


Fix \(c_0,\eta>0\) and \(m,H\). For sufficiently small
fixed \(\lambda>0\), let \(N\) be \(r\)-full and
\(r\lambda\ge c_0,\ q\lambda\ge1\), with \(q,r\) fixed
before the target height limit. Count \(A+B=N\) with \(B\)
\(q\)-full, arbitrary gcd, and (AB1) on \(A\) with \(X_0=N\).
There are \(\delta_*,\lambda_*>0\) such that at
\(\delta\le\delta_*,\lambda\le\lambda_*\) the count is

\[
 O_{c_0,\eta,m,H,\lambda,q,r,\varepsilon}
       \bigl(N^{\eta\lambda+\varepsilon}\bigr).
 \tag{EP1}
\]

Use the gcd, coefficient and full quotient-core construction
from the exact-power theorem. Freeze the \(g\)-supported
parts in full; the approximate valuation condition restricts
to the primes outside \(g\). For \(W\le N^{\eta/8}\) count
only the genuinely \(q\)-full part of \(B\), at cost
\(W^{1/q+\varepsilon}\). Otherwise the fixed \(W\)
has a full modulus at numerical weight bounded below by
a positive constant. In a failure let \(\lambda,\delta\to0\),
always taking each fixed parameter height limit first and
pigeonholing the finite exponent lists in (AB2). The
remaining original parts of \(A\) and \(B\) are written as
fixed coefficients times their narrow/wide factorizations.
The part of \(A\) without its coefficient is releasable
by (AB2).

In the high-block stopping proof for \(B\), if positive
numerical weight persists outside finitely many heavy
roots, a bounded-cost block of \(B\) has weight bounded
below and fullness comparable to \(1/\lambda\). Freeze
it and release \(A\) in the resulting unit congruence
by (GW), using the exact-output amplification above.
Choose the amplification order sufficiently large after
selecting the limiting block weight and fullness constant,
and before taking the order limit. The freeze costs less than the selected normalized exponent,
and the release can have arbitrarily small additional cost. Therefore \(B\) is releasable. Apply (GW) to the
complete ratio of the unfixed parts under the fixed full
target core and absorb the fixed supported coefficients
into the unit residue. A fixed exact ratio recovers the
pair by its fixed sum and coefficients. This completes
the generalized fixed-target estimate.


## The exponent of the minimally repaired primitive count

Prescribe arbitrary orders \(s_A,s_B,s_C\ge1\), with
\(M=\max s_i\). Allow a fixed nonempty sector \(\mathcal S\)
of primitive quotients \(u_A+u_B=u_C=W\). It can include
radical or valuation constraints on the quotients. For a
deficient prime \(p\) on its unique quotient letter \(i\),
meaning \(0<v_p(u_i)<s_i\), set
\(k_p=\max_j(s_j-v_p(u_j))\). Put
\[
             Z(\boldsymbol u)=\prod_{\text{deficient }p}p^{k_p},
       \qquad H(\boldsymbol u)=W Z(\boldsymbol u).
\]
As in (GU1), \(Z\) is the least multiplier that repairs the
prescribed orders. Let \(Q_{\mathcal S}(Y)\) count the primitive
quotients in \(\mathcal S\) with \(H\le Y\), and let
\(A_{\mathcal S}(Y)\) count all triples having the prescribed
orders and quotient in \(\mathcal S\), up to height \(Y\).
Then, for fixed orders,
\[
 \limsup_{Y\to\infty}\frac{\log A_{\mathcal S}(Y)}{\log Y}
 =
 \max\left\{\frac1M,\
       \limsup_{Y\to\infty}\frac{\log Q_{\mathcal S}(Y)}{\log Y}\right\}.
 \tag{MPA}
\]
The relevant positive Dirichlet series for summing over
primitive quotients is \(\sum_{\boldsymbol u\in\mathcal S}
H(\boldsymbol u)^{-s}\). Its abscissa of convergence is the
second exponent when \(\mathcal S\) is infinite.

For the two lower bounds, use \(Z\) as common divisor, or
choose one quotient and multiply its minimally repaired triple
by arbitrary \(M\)-th powers. For the upper bound, for each
quotient let \(P=\operatorname{rad}(u_Au_Bu_C)\). Write any
allowed multiplier \(g\) as \(Z f b\), with \(f\) supported
on \(P\) and \(b\) coprime to \(P\). Then \(b\) is \(M\)-full.
There are \(Y^{o(1)}\) choices of \(f\le Y\), uniformly in
\(P\le Y^3\), by the fixed-support Rankin bound. The number of
possible \(b\) is \(O_M((Y/H)^{1/M})\) after dropping \(f\)
from the height constraint. It follows that
\[
 A_{\mathcal S}(Y)\ll_{\boldsymbol s,\varepsilon}
 Y^{1/M+\varepsilon}
 \sum_{\substack{\boldsymbol u\in\mathcal S\\H(\boldsymbol u)\le Y}}
                      H(\boldsymbol u)^{-1/M}.
\]
Dyadic summation proves (MPA), and the same summation gives the
Dirichlet-series assertion. The term \(1/M\) in the unrestricted
count persists even when the primitive quotient series is
summable strictly below \(1/M\).

## A reverse transfer by repairing only the target

Take orders \((1,1,M)\). For \(R=\operatorname{rad}(c)\), the
least multiplier and the resulting target are exactly
\[
  Z_M(c)=\prod_{p\mid c}p^{(M-v_p(c))_+},\qquad
  H_M(c)=c Z_M(c)=\prod_{p\mid c}p^{\max(v_p(c),M)}.
 \tag{MPR}
\]
In particular \(R^M\le H_M(c)\le R^M c\). The map
\((a,b,c)\mapsto (Z_Ma,Z_Mb,H_M(c))\) is injective because the
primitive quotient recovers the input. If all target valuations
are at most \(M\), the new target is exactly \(R^M\), regardless
of their individual sizes. With scaled spectra of mass \(l_C\),
limit measure \(\nu_C\), numerical height \(h_C\), and escaped
moment \(e_C=h_C-\int r\,d\nu_C(r)\), \(ML\to q>0\) gives
\[
 \lim \log_X H_M
     = q l_C+\int(r-q)_+\,d\nu_C(r)+e_C .
\]
This follows by writing the excess integral as the full moment
minus the integral of \(\min(r,q)\), so large valuations that
escape the degree scale are included.

Let \(0<\xi<1/3\). For fixed \(\lambda,M\) define
\(\mathcal S_{\lambda,\xi}\) by
\[
 c^{(1/3-\xi)\lambda}\le\operatorname{rad}(c)
                  \le c^{(1/3+\xi)\lambda},\qquad
 \operatorname{rad}(ab)\le c^{(2/3+\xi)\lambda}.
\]
Write \(Q_{\lambda,\xi,M}(Y)\) for its minimally repaired count.
The following two existence statements are equivalent.

- The unrestricted exceptional-set count has an absolute
  coefficient gap below \(1/3\) at all sufficiently small fixed
  budgets.
- There are constants \(K,\xi,\sigma>0\), with \(0<\xi<1/3\),
  \(0<\sigma<1\) and
  \((1-\sigma)(1/3+\xi+1/K)<1/3\), such that for every
  sufficiently small fixed \(\lambda\), and \(M=\lceil K/\lambda\rceil\),
  \[
      Q_{\lambda,\xi,M}(Y)\ll_{\lambda,M,\varepsilon}
                          Y^{(1-\sigma)/M+\varepsilon}.
   \tag{MRG}
  \]

For the reverse transfer to abc, the sector covers every
strict exception except those with a letter whose radical is
at most \(c^{(1/3-\xi/2)\lambda}\). Indeed a target above the
upper cap pushes one summand below that threshold by
primitivity, and an excessive \(\operatorname{rad}(ab)\)
pushes the target below the lower cap. These omitted families
already have a fixed saving by the [one-letter quarter bound](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-quarter-inequality-for-a-letter-of-small-radical-mass).
On the remaining family, \(c\le X\) gives
\(H_M(c)\le X^{1+M(1/3+\xi)\lambda}\).
Using (MRG) and \(1/M\le\lambda/K\) gives the desired power
saving by the displayed condition on the constants.
In the other direction, the sector has original total radical
at most \(c^{(1+2\xi)\lambda}\), and \(H_M(c)\le Y\) implies
\(c\le Y^{1/(M(1/3-\xi)\lambda)}\) by the lower radical cap.
Apply the hypothesized abc gap \(\delta_0>0\) at budget
\((1+3\xi)\lambda\), absorbing the bounded \(c\) cases. The
resulting exponent is at most
\((1/3-\delta_0)(1+3\xi)/(M(1/3-\xi))+\varepsilon\).
First take \(\xi\) small in terms of \(\delta_0\), then choose
\(\sigma>0\) leaving a strict margin, and finally take \(K\)
large enough for the required inequality.

This places the remaining abc problem in a specific
one-prescribed-order tail of the powerful-triple count.
The [one-order weighted sector](#one-prescribed-fullness-order-and-captured-scaled-degrees)
has an extra restriction on the unprescribed degree spectra.
In the transfer above, (MPR) forces
\[
 (1/3-\xi)\lambda-1/M
 \ \le\ \frac{\log Z_M(c)}{M\log c}
 \ \le\ (1/3+\xi)\lambda ,
\]
so these quotients generally incur a positive normalized
repair cost. A bound for small repairs, or a power bound on
the unrestricted padded count with its unavoidable
\(Y^{1/M}\) term, does not give (MRG).

Under one-target repair, homogeneous determinants inherit common
row content. For the uncentred Taylor determinant the new fullness
can be entirely absorbed by this content, and the original
primitive valuation is sharp on an explicit family.
[Primitive content after repair](ABC_MINIMAL_REPAIR_PROJECTIVE_CONTENT.md)

A different count applies when the primitive quotient itself
factors into independent cyclotomic norm values. Opposite
primitive \(q\)-th powers give limiting coefficient
\(1/(\tau(q)+2)\) before repair. For same-side weighted
powers, two nonlinear rational factors of the primitive
target can give a repaired exponent \(1/(2M)\) when the
coefficient radical cost is negligible. With an irreducible
quadratic norm, bounded squarefree kernels alone still allow
the full repaired exponent, even on a diffuse target profile.
[Primitive cyclotomic split and its norm counterexample](ABC_MINIMAL_REPAIR_CYCLOTOMIC_SPLIT.md)


In the primitive weighted norm case, a normal-closure
determinant cancels the reference ideal denominator and
bounds the number of rays. For prime-power binomial fields,
the matching class-group torsion is subpower even as the
coefficients vary. For quadratic norms with square norm
defect, a separate factorization through an actual reference
solution gives rational factors and cancels its common ideal
without any short-height requirement.
[Primitive norm ray packing](ABC_PRIMITIVE_NORM_RAY_PACKING.md)
· [Splitting after cancelling the reference ideal](ABC_MINIMAL_REPAIR_IDEAL_CONTENT_SPLIT.md)

The minimally repaired exponent is also the exponent of
occurring balanced three-support unions. A finite number
of random external-prime tests approximates the count
with a power-small error while retaining the exact height
and all three radical conditions. On the target marginal
alone, the corresponding exact reduction is a weighted
sum of fixed-target counts.
[True support occurrence and modular averaging](ABC_MINIMAL_REPAIR_SUPPORT_OCCURRENCE.md)
