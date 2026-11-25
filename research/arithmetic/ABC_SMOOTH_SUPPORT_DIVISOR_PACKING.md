# Divisor packing with a polylogarithmic prime pool

A count by support unions allows every prime of a triple to contribute to the enumeration cost. The full powers on only part of those primes can already determine the primitive triple. For a fixed prime pool of polylogarithmic cardinality, the two effects combine into a discounted bound for that part of the radical. This applies to all smooth supports below the logarithmic power two, and combines with the balanced endpoint dispersion to treat a larger prime pool and allow some primes outside it.

The [effective radical bound](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
gives a stronger small-budget pool theorem. The divisor
packing here also works at fixed budgets away from zero
and tracks the angular multiplicity of a partial record.

## The count from three divisors

Let \(X\ge 3\), fix \(\kappa>0\), and let \(\mathcal P\) be any set of at most \(C_{\mathcal P}(\log X)^\kappa\) primes. Its elements and the set itself may vary with \(X\). Put
\[
             b_\kappa=\max(0,1-1/\kappa).
\]
Write \(R_{\mathcal P}(D)=\prod_{p\mid D,\ p\in\mathcal P}p\) and \(R_{\mathrm{out}}(D)=\prod_{p\mid D,\ p\notin\mathcal P}p\) for any positive integer \(D\). Let \(T_1,T_2,T_3\le X\) be positive pairwise coprime integers with \(T_3=\sigma T_1+\tau T_2\), for fixed signs. Suppose \(T_1\le X^{h_1}\), \(T_2\le X^{h_2}\) for \(h_1,h_2\in[0,1]\). Restrict to triples which admit divisors \(F_i\mid T_i\) with
\[
 \begin{split}
             F_1F_2F_3&>2X^{h_1+h_2},\\
             R_{\mathcal P}(F_1F_2F_3)&\le X^s,\qquad
             R_{\mathrm{out}}(F_1F_2F_3)\le X^u .
 \end{split}                                                     \tag{SP1}
\]
Here \(s,u\in[0,3]\). The number of these triples is
\[
             O_{\kappa,C_{\mathcal P},\varepsilon}
                    (X^{b_\kappa s+u+\varepsilon}).
                                                                    \tag{SP2}
\]
The constant is uniform in the caps and the prime pool. The same bound holds if a single extra full prime power with prime in \(\mathcal P\) may be omitted from the counted radical caps in (SP1), as long as it is included in the product on the first line.

Indeed, an ordered divisor record \((F_1,F_2,F_3)\) determines at most one triple. If two triples have this record, write their first two cofactor pairs as \((x,y)\) and \((x',y')\). Both satisfy \(\sigma F_1 x+\tau F_2 y\equiv0\bmod F_3\). Since \(\gcd(F_1,F_3)=1\), it follows that \(F_3\mid xy'-x'y\). The height caps and the first line of (SP1) give
\[
    |xy'-x'y|\le \frac{2X^{h_1+h_2}}{F_1F_2}<F_3 .
\]
Thus the determinant vanishes. Both cofactor vectors are positive and primitive, so they agree, and the additive relation fixes \(T_3\).

For any small fixed \(a>b_\kappa\) with \(a>0\), the number of squarefree integers at most \(X^s\) supported on \(\mathcal P\) is at most
\[
      X^{as}\prod_{p\in\mathcal P}(1+p^{-a})
               =X^{as+o(1)} .
\]
In fact, the sum of \(p^{-a}\) over a set of \(B\) primes is at most \(\sum_{j=1}^{B}(j+1)^{-a}\), and \(B\le C_{\mathcal P}(\log X)^\kappa\). For \(a<1\) the bound is \(O_a(B^{1-a})=o(\log X)\), and the bound for \(a\ge1\) is \(O_a(1+\log(1+B))\). This is the squarefree version of Rankin's elementary smooth-number argument. The outside radical has at most \(X^u\) choices by its size. The combined radical \(R\) has \(O(\log X/\log\log X)\) prime factors since \(R\le X^3\). Assigning its primes to the three letters costs \(3^{\omega(R)}=X^{o(1)}\). The number of positive exponent vectors on those primes for a product at most \(X^3\) is at most
\[
             \binom{\lfloor 3\log X/\log 2\rfloor+\omega(R)}
                        {\omega(R)}
                       = X^{o(1)} .
\]
This follows by bounding the sum of the exponents and using \(\omega(R)/\log X\to0\) in the binomial estimate. A single additional prime in \(\mathcal P\), its exponent and its letter have only polylogarithmically many choices, independently of its size. Taking \(a\) sufficiently close to \(b_\kappa\) and adding the exponent allowances gives (SP2). The Rankin bound and exponent enumeration are uniform in the selected support.
[Rankin's smooth-number argument](https://doi.org/10.1112/jlms/s1-13.4.242)

## A mixed condition on the powers supported in the pool

Let \(P_i=\prod_{p\mid T_i,\ p\in\mathcal P}p^{v_p(T_i)}\), the complete pool-supported part of letter \(i\). Suppose fixed parameters \(H_0>h_1+h_2\), \(S\ge0\) satisfy
\[
              \log_X(P_1P_2P_3)\ge H_0,\qquad
              \log_X\operatorname{rad}(P_1P_2P_3)\le S .
\]
For each fixed \(J\in(h_1+h_2,H_0]\), the number of triples is
\[
     O_{\kappa,C_{\mathcal P},J-h_1-h_2,\varepsilon}
         \left(X^{b_\kappa S J/H_0+\varepsilon}\right).
                                                                    \tag{SP3}
\]
Primes outside \(\mathcal P\) can carry arbitrary additional weight and have no cost in this estimate.

To check this, view a prime power \(p^e\) on its own letter as contributing numerical weight \(e\log_X p\) at radical cost \(\log_X p\). Take these powers from the \(P_i\) in descending order of \(e\). Selecting a fraction \(J/H_0\) of each power has weight at least \(J\) and cost at most \(S J/H_0\). The greedy fractional selection attains weight \(J\) at no greater cost. To see optimality, let \(e_*\) be the boundary order. The greedy vector maximizes numerical weight minus \(e_*\) times cost coordinatewise, since it takes all orders above \(e_*\) and none below it. This rules out reaching the same weight at a lower cost. The last power may be fractional. Restore it in full, pay its prime and exponent as the extra pool power in (SP2), and use \(X^J>2X^{h_1+h_2}\) for large \(X\). This rounding works even when \(\mathcal P\) contains polynomial-sized primes.

More generally suppose one may choose three divisors \(E_i\mid T_i\) on primes outside \(\mathcal P\), with total numerical weight at least \(W\) and total radical cost at most \(U\). For any fixed \(J\in[0,H_0]\) with \(J+W>h_1+h_2\), the bound becomes
\[
     O_{\kappa,C_{\mathcal P},J+W-h_1-h_2,\varepsilon}
         \left(X^{U+b_\kappa S J/H_0+\varepsilon}\right).
                                                                    \tag{SP3'}
\]
Use the same selection for pool weight \(J\), taking no pool power if \(J=0\), and multiply the selected part on each letter by \(E_i\). Disjoint support gives the displayed cost by (SP2).

In particular, for any fixed \(\lambda>0\), count the shell \(X/2<c\le X\) of primitive positive \(a+b=c\) with \(\operatorname{rad}(abc)\le c^\lambda\) and all its primes in \(\mathcal P\). Dyadic bins for the three complete letters cost only \(X^{o(1)}\). Use the two smallest height caps from a bin as \(h_1,h_2\). If \(H\) is the total numerical weight of the three letters, then \(H\ge1-o(1)\) and \(h_1+h_2\le 2H/3+o(1)\). Formula (SP3), or its greedy argument with these actual weights, gives
\[
                 O_{\lambda,\kappa,C_{\mathcal P},\varepsilon}
                    (X^{(2/3)b_\kappa\lambda+\varepsilon}).
                                                                    \tag{SP4}
\]
One can take \(J=h_1+h_2+\eta\) with a fixed \(\eta>0\) as small as needed relative to the final exponent allowance, and use the lower bin value for \(H_0\). If \(b_\kappa=0\), the bound is subpower. If \(1<\kappa<2\), the coefficient is below \(1/3\), with no small-budget restriction in this statement. The usual smooth-prime choice \(\mathcal P=\{p:p\le(\log X)^\kappa\}\) meets the cardinality bound. The same exponent applies to \(c\le X\) with the cutoff \((\log c)^\kappa\), by summing dyadic shells.

## Uniform saving beyond the logarithmic power two

There exist absolute \(\kappa_*>2\), \(\chi_*>0\), \(\delta_*>0\) and \(\lambda_*>0\) with the following property. For all fixed \(0<\lambda<\lambda_*\), and any prime pool \(\mathcal P\) of size at most \(C_{\mathcal P}(\log X)^{\kappa_*}\), the shell count of primitive positive \(a+b=c\) with
\[
             X/2<c\le X,\qquad
             \operatorname{rad}(abc)<c^\lambda,\qquad
             R_{\mathrm{out}}(abc)\le c^{\chi_*\lambda}
\]
is
\[
       O_{\lambda,C_{\mathcal P},\varepsilon}
                     (X^{(1/3-\delta_*)\lambda+\varepsilon}).
                                                                    \tag{SP5}
\]
The constants are uniform in the choice of the prime pool. In particular one can take all primes at most \((\log X)^{\kappa_*}\) and then sum shells with the cutoff \((\log c)^{\kappa_*}\).

Suppose no positive \(\delta_*,\lambda_*\) work for the pool parameters chosen below. The existing global limiting bound produces a sequence of actual balanced anatomic boxes at exponent \(d\to1/3\) and with total root exponent \(L/\lambda\to1\). Here \(L\) is the sum of the base-\(X\) root-height exponents, and the degree measures are normalized by \(L\). At each fixed budget use first the inner height limit and retain only triples satisfying the pool restriction. The unrestricted count of the same boxes also has exponent tending to \(1/3\), so the uniform pooled degree dispersion (PQ2) applies.
[Balanced endpoint and pooled degree dispersion](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-pooled-degree-inequality-and-primitive-shift-fibres)

Let \(\eta_{\rm PQ}>0\) be the absolute lower bound in (PQ2), and fix \(0<\eta\le\min(\eta_{\rm PQ}/8,1/16)\). In any sufficiently late box, select the top half of the pooled anatomic root mass on two letters. The selected numerical weight exceeds half the combined numerical weight by at least \(2\eta\). Select the top half of the root mass on the third letter as well, supplying at least half its numerical weight.

This fractional selection transfers to actual prime powers of each tuple in the box. At fixed truncation, the dyadic root coordinates have actual logarithmic masses equal to the recorded ones up to inner \(o(1)\). In the decomposition \(T_i=T_{i,0}\prod_e T_{i,e}^e\), the coefficient \(T_{i,0}\) contributes \(o(1)\) numerical weight in the outer limit, and every contribution of \(p\) at an artificial degree \(e\) satisfies \(e\le v_p(T_i)\). If all selected fractions of those occurrences of \(p\) account for numerical weight \(W_p\), take the fraction \(W_p/(v_p(T_i)\log_X p)\) of the actual power. This is at most one, because the sum of the original exponents contributed by the occurrences is at most \(v_p(T_i)\). Its radical cost \(W_p/v_p(T_i)\) is no greater than the corresponding selected total anatomic root cost. This remains valid if an artificial root contributes several copies of a prime. Thus for each actual tuple in the box there is a fractional selection of its own prime powers with cost at most \(L/2+o(\lambda)\) and weight at least \(H/2+\eta\), where \(H=\sum_i\log_X T_i\). Coefficient and box losses were absorbed by the chosen margin.
[Anatomic decomposition and comparison with the actual radical](https://arxiv.org/html/2410.12234v2#S2)

Put \(d_0=\eta/12\) and \(t=2/3-d_0\). Combine the half selection and the all-prime selection with coefficients \(2(1-t)\) and \(2t-1\), respectively. The result is a fractional selection of cost at most \(t\lambda+o(\lambda)\) and numerical weight at least
\[
              tH+2(1-t)\eta>\frac23 H+\frac{\eta}{3}
\]
since \(H\le3\). The sum of the two smaller complete-letter height caps is at most \(2H/3+o(1)\) with the inner limit at fixed box. Sort all prime powers of the actual tuple in descending order of valuation and greedily select them to exceed the sum of those two caps by, say, \(\eta/4\). Before rounding the boundary power, the cost is at most \(t\lambda+o(\lambda)\), by the fractional optimality already used in (SP3).

If the boundary prime is in \(\mathcal P\), put it in the extra pool slot and leave it out of the paid caps. If it is outside \(\mathcal P\), pay for it with the rest of the selected outside radical. In either case the paid pool-supported radical cost is at most \(t\lambda+o(\lambda)\). The entire outside radical cost is at most \(\chi_*\lambda\) by hypothesis, including any last outside prime whose full rounding might have a positive cost. The rounded numerical weight still gives the product inequality (SP1). This argument needs neither a size bound on each prime in the pool nor smallness of every outside valuation.

Apply (SP2) in each of the selected boxes. Write \(\kappa_*=2+d_0\) and \(\chi_*=d_0/12\). The coefficient of \(\lambda\), before arbitrary exponent allowances, is at most
\[
                  b_{\kappa_*}(2/3-d_0)+\chi_*+o(1)
                                      \le 1/3-d_0/4+o(1).
\]
Indeed \(1/3-b_{\kappa_*}(2/3-d_0)=(2d_0/3+d_0^2)/(2+d_0)\ge d_0/3\). This contradicts the selected count exponent. The strict inequalities use the inner height limits at the fixed budgets before taking the outer sequence. They prove existence of the thresholds in (SP5), with an explicit choice of \(\kappa_*,\chi_*\) from the dispersion margin. This argument does not compute a numerical threshold for \(\delta_*\).

The support restriction in (SP5) is on a fixed pool for each shell. It also allows a fixed collection of \(X^{o(1)}\) candidate pools by summation. A pool can be the union of the prime supports of at most \((\log X)^{\kappa_*-1}\) fixed polynomial-height triples, since each uses \(O(\log X/\log\log X)\) primes. The result then bounds the subfamily whose radical is almost entirely supported on that union. Choosing pools separately for the members of the family would require charging the number of different choices.

## Consequence for minimally repaired primitive quotients

Use the central repair sector with
\[
 c^{(1/3-\xi)\lambda}\le R_c\le c^{(1/3+\xi)\lambda},
 \qquad
 R_aR_b\le c^{(2/3+\xi)\lambda},
 \qquad
 M=\lceil K/\lambda\rceil .
\]
For sufficiently small absolute \(\xi>0\), additionally require outside radical at most \(c^{\chi_*\lambda}\) for the choice \(\mathcal P=\{p:p\le(\log c)^{\kappa_*}\}\), or for any fixed pool of at most \(C_{\mathcal P}(\log Y)^{\kappa_*}\) primes. Then the number of distinct least-repaired primitive quotients with repair height at most \(Y\) is
\[
       O_{\lambda,M,C_{\mathcal P},\varepsilon}(Y^{(1-\sigma)/M+\varepsilon})
\]
for some absolute \(\sigma>0\).

Indeed \(H_M(c)=\prod_{p\mid c}p^{\max(M,v_p(c))}\ge R_c^M\) forces \(c\le Y^{1/[M(1/3-\xi)\lambda]}\). The sector satisfies the radical bound at budget \((1+3\xi)\lambda\), and its outside cap is no larger than the cap in (SP5) for this slightly enlarged budget. Apply the shell-summed (SP5) at the upper bound for \(c\). The exponent before the arbitrary allowance is
\[
    \frac{(1/3-\delta_*)(1+3\xi)}{M(1/3-\xi)} < \frac{1-\sigma}{M}
\]
by taking \(\xi\) sufficiently small in terms of \(\delta_*\). Each original primitive triple determines one quotient record, so the count uses the least repair height.

For the fixed pool at scale \(Y\), first count the triples with \(c\le Y^{1/(4M\lambda')}\), where \(\lambda'=(1+3\xi)\lambda\). The count up to any height \(Z\) at budget \(\lambda'\) is \(O_{\lambda',\varepsilon}(Z^{\lambda'+\varepsilon})\) by enumerating the total radical, assigning its primes to the letters, and counting exponent vectors as in (SP2). This gives \(O(Y^{1/(4M)+\varepsilon})\) for the small-\(c\) part. On the rest, the size of the fixed pool is \(O_{\lambda,M,C_{\mathcal P}}((\log c)^{\kappa_*})\), so the same pool can be used throughout the dyadic shells of (SP5).
[Minimal target repair and the original radical sector](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md)

## Comparison with counts of smooth solutions

Drappeau's Theorem 1.5 gives a GRH-conditional asymptotic for the full smooth \(a+b=c\) count when the lower smoothness parameter is \((\log X)^{8+\varepsilon}\). Jain's Theorem 1 gives a representation of sufficiently large integers by one power-smooth and one polylogarithmically smooth summand. Both statements concern bounds on the largest prime factor. Formulas (SP3) and (SP5) charge the product of distinct primes in selected divisors. The Euler product bounds the number of these radicals, and the fixed divisor record gives a determinant obstruction to a second primitive triple.
[Drappeau, Theorem 1.5](https://arxiv.org/pdf/1203.1742#page=4) · [Jain, Theorem 1](https://arxiv.org/html/2611.21136v1#S1)


## Charged radicals and angular packing

Let the height bound be \(X\). Assign every prime \(p\le X\) a number
\(q_p\in(0,1]\), and put
\[
 s(p)=-\log_X q_p,\qquad C(F)=\sum_{p\mid F}s(p),\qquad
 B=\sum_{p\le X}q_p .
\]
A squarefree support \(R\) with \(C(R)\le b\) contributes at least
\(X^{-b}\) to the product \(\prod_{p\le X}(1+q_p)\). Hence the number
of such supports is at most \(X^b e^B\). With product height at most
\(X^3\), allocation to three letters and recovery of exact
exponents cost \(X^{o(1)}\), uniformly as in (SP2).

For each triple in a fixed complete-letter height bin, allow a choice
of divisors \(F_i\mid T_i\), and write \(w=\log_X(F_1F_2F_3)\).
If
\[
 C(F_1F_2F_3)+(h_1+h_2-w)_+\le b ,
 \tag{SP6}
\]
the number of triples is at most \(X^{b+o(1)}e^B\). This also counts
records whose product height is below the uniqueness threshold in
(SP1). Indeed, for a fixed record put \(Y_i=X^{h_i}/F_i\) for
\(i=1,2\). The cofactor pair \((x,y)=(T_1/F_1,T_2/F_2)\) is
primitive. Any two different pairs satisfy
\(F_3\mid xy'-x'y\) by the congruence used in (SP1). Their directions
after scaling by \((Y_1^{-1},Y_2^{-1})\) have angular distance
\(\gg F_3/(Y_1Y_2)\), since the scaled vectors have lengths at most
\(\sqrt{2}\) and the scaled determinant is at least
\(F_3/(Y_1Y_2)\) in absolute value. The angle interval has bounded
length. Therefore the record admits
\(O(1+X^{h_1+h_2-w})\) points. Binning \(w\) on a mesh of size
\(1/\log X\) and using (SP6) in each bin proves the count. The same
proof applies in an angular interval of length \(O(X^{-v})\) for
the original pairs scaled by their complete-letter caps, with
\((h_1+h_2-w-v)_+\) in (SP6). The scaled directions in the record
proof agree with those original scaled directions.

The charge can be specified by the incidence frequencies of the
family itself. For any nonempty family \(\mathcal S\) of primitive
triples with height at most \(X\), let \(f_p\) be the proportion of
rows in \(\mathcal S\) containing \(p\) in their radical. Take
\[
 q_p=\max(f_p,p^{-1}),
 \qquad s(p)=\min\{\log_X p,-\log_X f_p\}
\]
with the second entry interpreted as \(+\infty\) when \(f_p=0\).
Then \(B=o(\log X)\) uniformly in \(\mathcal S\), since
\[
 \sum_p f_p
   =\frac{1}{|\mathcal S|}\sum_{(a,b,c)\in\mathcal S}\omega(abc)
   =O\!\left(\frac{\log X}{\log\log X}\right)
 \qquad\text{and}\qquad
 \sum_{p\le X}p^{-1}=O(\log\log X).
\]
For the first bound compare the product of \(j\) distinct primes
with \((j+1)!\). The second follows by bounding
\(\sum_{p\le X}p^{-1-1/\log X}\) by
\(\log\zeta(1+1/\log X)\). Thus empirical frequencies can be used
without enumerating families or prime pools. For example, if
\(|\mathcal S|\ge X^d\), then, outside \(X^{d-\tau+o(1)}\) rows,
every admissible choice of the three divisors satisfies
\[
 C(F_1F_2F_3)+(h_1+h_2-w)_+>d-\tau
\]
for any fixed \(\tau>0\) in the chosen bin. Complete-letter dyadic
height bins contribute only a power of \(\log X\).

In particular every pair of complete letters has combined cost at
least \(d-\tau\) outside such an exception, since the two values
determine the third. Hence the total score is at least
\(3(d-\tau)/2\) there. If \(\mathcal S\) is an exceptional family at
budget \(\lambda\) with \(d\ge (1/3-\eta)\lambda\), then for any
\(0\le\alpha<1/2\) the mass of primes with \(f_p<p^{-\alpha}\)
is at least
\[
 \sum_{\substack{p\mid abc\\ f_p<p^{-\alpha}}}\log_X p
   \ \ge\ \frac{1/2-\alpha-3(\eta+\tau/\lambda)/2}{1-\alpha}\lambda
\]
on the same rows. Indeed the other primes have charge at most
\(\alpha\log_X p\), and the total actual radical weight is at
most \(\lambda\). The packing estimate also tests the numerical
location of this mass through (SP6).

