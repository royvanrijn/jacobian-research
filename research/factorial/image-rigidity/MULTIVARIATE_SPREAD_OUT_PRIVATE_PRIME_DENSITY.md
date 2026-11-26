# All-degree separation for most spread-out ordinary tails

Fix \(J\ge4\) at least the [absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md). Its value is unevaluated. Fix positive constants \(c,C,b\) with \(c<C\). For each \(R\), append \(R\) ordinary unequal rows with independently chosen smaller exponents \(a_i\), each uniform on an integer interval \(I_i\) satisfying
\[
 I_i\subset[cR,CR],\qquad T_i=|I_i|\ge bR/\sqrt{\log R}.
\]
The positive integer gaps \(d_i\) and orientations are fixed in advance and may vary by row and by \(R\). These random choices give the ordinary counting measure on ordered tuples of native integer exponents.

**Theorem.** For any fixed ordinary core, the proportion of tuples whose ordinary binomial moments are coprime at every adjacent pair is at least \(1-\exp(-R^{1-o(1)})\). If even one tuple in the sampling box has total curvature below twelve, the proof gives a stronger conclusion. Outside a set of the same size, each terminal \(2\le m<J\) has at least
\[
 \frac{c_1 R}{\log R}\exp(-c_2\sqrt{\log R})
 \tag{1}
\]
distinct primes at which the two coefficient Newton polygons have disjoint slopes. The positive constants may depend on \(J,c,C,b\) and the core. All slopes are computed in the original common argument, so rational coprimality also gives coprimality over every characteristic-zero field.

The theorem allows a varying fixed part of the profile as well. For each pair, strip the smaller-end factorial from each of its rows. Form the two full stripped ramps of every fixed row. It is enough that the number of distinct primes \(p\ge cR\) dividing any of their factors at that pair is \(o(R/(\log R)^{3/2})\), uniformly over the finitely many pairs. There is no multiplicity charge for repeated fixed ramps.

## An exact prime certificate

At order \(j\), a row with smaller exponent \(a\) and gap \(d\) contributes the first \(kd\) entries of
\[
 ja+1,\ldots,ja+jd
\]
at coefficient index \(k\) in the increasing orientation, and the first \((j-k)d\) entries in the decreasing orientation. The removed factorial \((ja)!\) only changes the scalar of that moment.

Fix a terminal \(m<J\) and write \(n=m-1\). Suppose \(p=ma_i+1\) is prime, \(p>J^2d_i+J\), and \(p\) divides no full stripped factor of any other row at order \(m\) or \(n\). The other entries of row \(i\) at order \(m\) differ from \(p\) by positive integers less than \(p\). At order \(n\),
\[
                    m(na_i+h)-np=mh-n,\qquad 1\le h\le nd_i .
\]
The right side is nonzero since \(m\nmid n\), and its magnitude is less than \(p\). Hence \(p\) divides no entry at order \(n\). The binomial coefficients are \(p\)-units since \(p>J\).

The order-\(m\) ordinates are therefore one everywhere except at the endpoint where the first entry of the marked row is absent. There the ordinate is zero. Its lower polygon consists of one edge of slope \(1/m\) or \(-1/m\) according to orientation. The order-\(n\) polygon is horizontal. The Newton polygon product rule rules out a nonconstant rational common factor. The rule is stated for integer polynomials in reversed coefficient order in the cited source. A common rational factor can be cleared of denominators, and reversing the order reverses slopes for both parents. [Filaseta and Rocks, Theorem 2.1](https://www.impan.pl/shop/publication/transaction/download/product/91656)

## Many private primes for every tested pair

Write \(S=\sum_i d_i^2/M_i\) over all unequal rows, where \(M_i\) is the larger exponent. If every tuple has \(S\ge12\), the [absolute curvature twelve theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md) proves the conclusion. Otherwise fix any tuple with \(S<12\), and put \(d_{\max}=\max_i d_i\) over the sampled rows. All sampled larger exponents are at most \(CR+d_{\max}\). This gives
\[
       d_{\max}^2<12(CR+d_{\max}),\qquad
       \sum_i d_i^2<12(CR+d_{\max}).
\]
In particular \(d_{\max}=O(\sqrt R)\), \(\sum_i d_i^2=O(R)\), and \(\sum_i d_i=O(R)\) by Cauchy--Schwarz. The bounds are properties of the fixed gaps, so they hold across the sampling box. The candidate primes \(p=ma_i+1\) satisfy \(p>J^2d_{\max}+J\) for all sufficiently large \(R\).

Fix \(m\) and call a prime \(p\ge cR\) forbidden when it divides a fixed ramp at either order. For a fixed core there are no such primes for all sufficiently large \(R\). For the stated varying fixed part their number is \(o(R/(\log R)^{3/2})\).

The progression \(ma_i+1\), with \(a_i\in I_i\), contains \(\gg T_i/\log R\) primes uniformly over the indicated windows. In fact, subtract the fixed-modulus estimates
\[
             \pi(x;m,1)=\frac{\operatorname{Li}(x)}{\varphi(m)}
                              +O_m\bigl(x/(\log x)^2\bigr)
\]
at the endpoints \(x=m(\min I_i-1)+1\), \(y=m\max I_i+1\). Both endpoints are comparable to \(R\), and the main term is at least \(mT_i/(\varphi(m)\log y)\). Here \(T_i\ge bR/\sqrt{\log R}\) makes the endpoint error \(O(R/(\log R)^2)\) smaller in order. For \(m=2\) the progression counts all odd primes and the same estimate follows from the ordinary prime number theorem with this error. [Bennett, Martin, O'Bryant and Rechnitzer, Theorem 1.3 and (1.15)](https://personal.math.ubc.ca/~gerg/papers/downloads/EBPAP.pdf)

Deleting the forbidden primes leaves a uniform probability at least \(c_3/\log R\) that row \(i\) offers a prime \(p=ma_i+1\). For fixed \(p\) and any row \(h\ne i\), each possible factor \(ja_h+u\), at order \(j=n,m\), imposes a single residue class modulo \(p\). This uses \(p>J\). The probability of that class in the interval is at most \(1/p+1/T_h\). Thus the chance that this row hits \(p\) is at most
\[
      q_h(p)=(m+n)d_h(1/p+1/T_h).
\]
The gap bounds give \(\max_h q_h(p)=O(\sqrt{\log R/R})\) and \(\sum_h q_h(p)\le K\sqrt{\log R}\) for a fixed \(K\). Conditioning on the candidate prime of row \(i\) leaves all other choices independent. Their probability of avoiding \(p\) is at least
\[
        \prod_{h\ne i}(1-q_h(p))\ \ge\ \exp(-2K\sqrt{\log R})
\]
for sufficiently large \(R\), by \(\log(1-t)\ge-2t\) on \(0\le t\le1/2\).

Let \(Z_m\) count distinct primes satisfying the certificate. Each prime satisfying it has just one offering row. Summing the preceding conditional lower bound over rows therefore gives
\[
         \mu_m=\mathbb E Z_m
            \ \ge\ c_3 R\,(\log R)^{-1}\exp(-2K\sqrt{\log R}).
\]
Changing one \(a_i\) can change the status only of a prime \(p\ge cR\) dividing one of this row's old or new full stripped factors. There are at most \(2(m+n)d_i\le4Jd_i\) such factors. Each factor is at most \(J(CR+d_{\max})\), which is less than \((cR)^2\) for sufficiently large \(R\). It thus has at most one prime divisor at least \(cR\). Consequently \(Z_m\) changes by at most \(4Jd_i\) under a change of row \(i\)'s base.

For completeness, expose the independent bases in sequence and consider the conditional expectations of \(Z_m\) after each exposure. At any given prefix, when row \(i\) is exposed its conditional expectation ranges over an interval of width at most \(4Jd_i\), by averaging the one-row change bound over the same future choices. Its difference from the preceding conditional expectation has absolute value at most \(4Jd_i\). For a centered random variable \(Y\) conditional on the preceding choices, with \(|Y|\le L\), interpolation of \(e^{tY}\) between \(Y=-L,L\) gives conditional expectation at most \(\cosh(tL)\le e^{t^2L^2/2}\). The last inequality follows termwise from the power series and \((2h)!\ge 2^h h!\). Apply this at each step to the negative conditional-expectation difference, with \(L=4Jd_i\). Iteration and Markov's inequality give
\[
 \Pr[ Z_m<\mu_m/2]
 \le \inf_{t>0}\exp\left(-t\mu_m/2+
                        \frac{t^2(4J)^2}{2}\sum_i d_i^2\right)
 = \exp\left(-\frac{\mu_m^2}{8(4J)^2\sum_i d_i^2}\right).
\]
The exponent in magnitude is at least \(c_4R(\log R)^{-2}\exp(-c_5\sqrt{\log R})\) for positive constants \(c_4,c_5\).
A finite union over \(2\le m<J\) gives (1), and the absolute angular theorem supplies all remaining orders. A tuple with (1) retains its certificates at all orders after changing any subset of bases whose corresponding gaps sum to less than the quantity in (1) divided by \(4J\), provided their new smaller exponents remain in \([cR,CR]\). Conditioning on any subset of the tuple box of relative size at least \(\exp(-R^\delta)\), for a fixed \(0<\delta<1\), preserves the exceptional proportion bound \(\exp(-R^{1-o(1)})\). Indeed, divide the unconditioned exceptional proportion by the size of the conditioning event. The resulting added term \(R^\delta\) in the exponent is smaller in order than \(c_4R(\log R)^{-2}\exp(-c_5\sqrt{\log R})\).

## A uniformly positive diffuse limit

Given any fixed core and any target tail curvature \(\tau>0\), let \(R\) be even and use \(R/2\) rows of each orientation with \(d_i=1\). Take every window to start at \(A=\lfloor R/\tau\rfloor\) and have length \(\lceil R/\sqrt{\log R}\rceil\). Every tuple has tail signed degree zero, largest relative gap \(O(R^{-1})\), and
\[
           \sum_{i=1}^R\frac1{a_i+1}
                         =\tau+O((\log R)^{-1/2})
\]
uniformly, by bounding the sum using the two ends of the window. This covers approach to any prescribed positive diffuse strength, including strengths exceptional for a fixed core's Gaussian comparison. The proof uses exact factorial coefficients throughout. The exceptional set of tuples may still contain supports at those strengths. Uniform exclusion of correlated or adversarial height choices remains open.

## Exact checks

The [finite checker](verify_multivariate_spread_out_private_prime_density.py) enumerates 10,290 ordered tuples and order pairs. It replays 5,896 private prime certificates with both orientations, removes 2,680 offers blocked by fixed rows, verifies 14 rational expectation identities by independent residue counting, and compares 12 marked arrays with direct factorial products. These checks concern the finite certificate and counting identity. The asymptotic estimate is proved above.
