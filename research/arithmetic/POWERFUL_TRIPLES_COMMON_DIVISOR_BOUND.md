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
