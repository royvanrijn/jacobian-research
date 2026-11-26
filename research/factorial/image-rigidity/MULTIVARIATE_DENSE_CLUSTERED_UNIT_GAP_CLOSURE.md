# All-degree closure for dense clusters of unit-gap rows

Fix an ordinary two-monomial core and append unequal rows of gap one. Let the appended smaller exponents lie in \([M,M+w]\), and let \(s\) be the number of occupied heights. At each height \(a\), write \(r_+(a)\) and \(r_-(a)\) for the two orientation multiplicities. The multiplicities are arbitrary. Divide the occupied heights into two classes according to whether \(r_+(a)=r_-(a)\), and let \(L\in\{1,2\}\) be the number of nonempty classes.

**Theorem.** Fix \(0<\theta\le1\) and \(0<\alpha<\theta/L\). Along any such family with
\[
                         w\longrightarrow\infty,\qquad
         s\ge\theta w,\qquad w\le M^\alpha,
\]
every sufficiently large profile has coprime consecutive factorial moments in every degree over every characteristic-zero field. One may use \(L=2\) uniformly, even when the number of classes varies along the family. If all occupied heights have equal orientation counts, or if every occupied height has unequal counts, one may use \(L=1\).

The fixed core may be replaced by arbitrary retained rows. For each adjacent pair with terminal \(4\le N<J\), where \(J\ge4\) is at least the [absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md), discard the smaller-end factorial at each retained row. Define \(B_N\) to be the number of distinct primes greater than \(2J^2(w+1)\) in the two full retained ramps. The same conclusion holds if \(B_N=o(w)\) uniformly over the finitely many terminals. In particular repetitions of a retained row have the same charge as a single copy.

The proof marks one factor at each occupied height and excludes every prime shared with a different height by an integer determinant. The aggregate size of the marks guarantees a useful large prime. This allows an unbounded number of occupied heights and needs no congruences between them. The single-height and fixed affine-template results have different hypotheses. [Common-height closure](MULTIVARIATE_COMMON_HEIGHT_BOUNDED_GAP_TAIL_CLOSURE.md) [Affine count-simplex certificate](MULTIVARIATE_AFFINE_TAIL_COUNT_POLYTOPE.md) [Aligned unit gaps](MULTIVARIATE_ALIGNED_UNIT_GAP_CENTRAL_NEWTON.md)

## The two marks and their exact polygons

Consider a terminal \(4\le N<J\). Let \(m\) be the even order in the pair and \(n\) its odd neighbor, with \(m\ge4\). Put
\[
        t_a=\begin{cases}
                 1,&r_+(a)\ne r_-(a),\\
                 m/2+1,&r_+(a)=r_-(a),
              \end{cases}
       \qquad x_a=ma+t_a,\qquad X=2J^2(w+1).
\]
In both cases \(1\le t_a<m\). Stripping the order-\(j\) factorial \((jb)!\) at an appended row of smaller exponent \(b\) leaves factors \(jb+h\), \(1\le h\le j\). At coefficient index \(k\) the increasing row takes the first \(k\) factors, and the decreasing row takes the first \(j-k\) factors.

Suppose \(p>X\) divides \(x_a\). In the even order any different factor \(mb+h\) has nonzero difference from \(x_a\) of magnitude less than \(p\). Literal equality is possible only at \(b=a,h=t_a\), since \(1\le h,t_a\le m\). For the odd order, the cross difference is
\[
           m(nb+h)-n x_a=mn(b-a)+mh-nt_a.
\]
It is nonzero because \(\gcd(m,n)=1\) and \(m\nmid t_a\), and has magnitude at most \(J^2w+2J^2\le X\). The prime therefore sees only the one literal factor \(x_a\) in appended rows. If \(p\) is absent from the retained ramps, the odd polygon is horizontal. The binomial coefficients are units because \(p>J\).

Write \(E=v_p(x_a)>0\). For the mark \(t_a=1\), the even ordinates at the endpoints are \(r_-(a)E\) and \(r_+(a)E\); every interior ordinate is their sum. The polygon is a single edge of slope \((r_+(a)-r_-(a))E/m\), which is nonzero. For \(t_a=m/2+1\), the ordinates are \(r_-(a)E\) to the left of \(m/2\), zero at \(m/2\), and \(r_+(a)E\) to the right. Both multiplicities are positive, and the two edges have nonzero slopes of opposite signs. In either case the slope lists are disjoint. The Newton polygon product rule now excludes a nonconstant common rational factor. [Filaseta and Rocks, Theorem 2.1](https://www.impan.pl/shop/publication/transaction/download/product/91656)

## A product bound for the chosen marks

Any prime \(p>X\) can divide at most one of the \(x_a\), since these integers are distinct and their span is at most \(mw+m<X\). Thus at most \(B_N\) marks are hit by a prime \(p>X\) from a retained ramp. If a common factor exists, every remaining mark must be \(X\)-smooth, meaning that all its prime divisors are at most \(X\).

In each class \(t_a=t\), normalize by the exact integer \(g_t=\gcd(m,t)\). The resulting mark
\[
              y_a=(ma+t)/g_t=(m/g_t)a+t/g_t
\]
is an integer at least \(M\) and at most \(Y=J(M+w+1)\). Its affine coefficients are coprime. Both normalized classes are still \(X\)-smooth after the deletion above. In a fixed class, for any prime power \(p^e\), the possible heights with \(p^e\mid y_a\) lie in a single residue class modulo \(p^e\), or there are none. Indeed \(p\mid m/g_t\) implies \(p\nmid t/g_t\). Otherwise the leading coefficient is a unit modulo \(p^e\). Thus there are at most \(w/p^e+1\) marked heights per class with this divisibility, and no power exceeding \(Y\) occurs. Taking the log of the product of the remaining \(y_a\) gives the necessary inequality
\[
 (s-B_N)_+ \log M
      \le L w\sum_{p\le X}\frac{\log p}{p-1}
                              +L\pi(X)\log Y .                 \tag{1}
\]
This argument charges only distinct obstructing primes of the retained rows. For example,
\[
 B_N\le \frac{2J}{\log X}
           \sum_{\rho\ {\rm distinct}}d_\rho\log(JM_\rho),
\]
where \(d_\rho\) is the nonzero gap and \(M_\rho\) the larger exponent of a retained type. Each stripped factor of size at most \(JM_\rho\) has at most \(\log(JM_\rho)/\log X\) prime divisors greater than \(X\). A fixed core has \(B_N=0\) for large \(w\).

## Elementary prime estimates

Put \(\vartheta(x)=\sum_{p\le x}\log p\) and \(\psi(x)=\sum_{p^e\le x}\log p\). The primes in \((u,2u]\) divide \(\binom{2u}{u}\) for integer \(u\ge1\), and this coefficient is at most \(4^u\). Summing over dyadic intervals gives \(\vartheta(x)=O(x)\), and summing over powers gives
\(\psi(x)=O(x+\sqrt{x}\log x)=O(x)\). The estimate
\(\pi(X)=O(X/\log X)\) follows by splitting the primes at \(\sqrt X\).

For \(U=\lfloor X\rfloor\), the valuation formula for \(U!\) gives
\[
 \sum_{p^e\le U}\frac{\log p}{p^e}
            \le\frac{\log(U!)+\psi(U)}{U}\le\log U+O(1).
\]
The first omitted term of the \(p\)-series is less than \(1/U\), and its entire omitted tail is at most \(2/U\). Summing these tails against \(\log p\) over \(p\le U\) costs \(O(1)\). Consequently \(\sum_{p\le X}\log p/(p-1)\le\log X+O(1)\).

Since \(\alpha<\theta/L\le1\), we have \(w\le M\) for the indicated families. The right side of (1), divided by \(w\log M\), is at most
\[
                L\,\frac{\log w}{\log M}
                    +O_J\left(\frac1{\log M}+\frac1{\log w}\right)
                      \le L\alpha+o(1).
\]
The left side divided by the same quantity is at least \(\theta-o(1)\). This contradicts \(L\alpha<\theta\) for all sufficiently large \(w\). The first two adjacent pairs follow from [ordinary low-order coprimality](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md), and \(N\ge J\) follows from the absolute angular cutoff. The coprimality conclusion transfers to characteristic-zero fields by a rational Bezout identity.

## Positive diffuse strength with arbitrary cluster width

For any \(\tau>0\), take all \(w+1\) heights in the band and put
\(\lfloor\tau M/(2(w+1))\rfloor\) rows of each orientation at each height. Choose any fixed \(0<\alpha<1\) and let \(w=\lfloor M^\alpha\rfloor\). The integer multiplicity is positive for large \(M\). The signed tail degree vanishes, the largest relative gap tends to zero, and the curvature tends to \(\tau\), because \(a+1=M(1+o(1))\) uniformly and \(w=o(M)\). This is a deterministic all-degree result at any prescribed positive diffuse strength, with arbitrarily many genuinely different heights. The theorem also permits adversarially varying multiplicities and the indicated additional retained rows. It leaves wider or sparsely occupied clusters outside its hypotheses.

## Exact checks

The [finite checker](verify_multivariate_dense_clustered_unit_gap_closure.py) tests 576 order and cluster combinations. It checks 898 large-prime certificates, including both neighboring odd orders and both mark branches, plus four certificates of depth two or three. It also tests 300 normalizations with nontrivial gcd, 4,568 residue-count bounds and 490 retained-overlap examples. These are finite replays of the arithmetic. The product estimate and the uniform height conclusion are proved above.
