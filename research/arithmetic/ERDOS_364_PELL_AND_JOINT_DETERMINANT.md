# Pell descent and a joint surface cover for Erdős #364

Let \(C(X)\) count integers \(2\le n\le X\) for which \(n-1,n,n+1\) are all powerful. An integer is powerful if every positive prime valuation is at least two. The Pell reduction separates the candidates whose centers are fundamental Pell traces. The [averaged pair argument](ERDOS_364_AVERAGED_PAIR_COMPONENT_COUNT.md) gives exponent \(6/25-10^{-7}\) for \(C(X)\) using the stated fixed-shift estimates, improving the \(293/1220\) contact bound recorded below. The joint determinant gives a separate surface cover in every dyadic kernel box.

## Pell descent and minimal repair

A candidate has \(4\mid n\), since every other residue of \(n\) modulo 4 puts a term in \(2\bmod 4\). The three terms are then pairwise coprime. Write \(n+i=a_i^2b_i^3\) with \(b_i\) squarefree and \(a_i,b_i\) positive, allowing common primes in \(a_i,b_i\). The \(b_i\) are pairwise coprime. The total kernel \(D=b_{-1}b_0b_1\) gives \(n^3-n=D^3(a_{-1}a_0a_1)^2\). Thus in the model \(DY^2=n^3-n\) the extra divisibility is \(D\mid Y\), and on \(z^2=x^3-D^2x\) the associated point has \(x=Dn\), \(z=D^2Y\).

For the Pell reduction, put \(q=b_{-1}b_1\) and \(Z=q a_{-1}a_1\). Then \(n^2-qZ^2=1\) and \(q\mid Z\). Conversely, if \(n\ge2\) is even and powerful, \(q\) is squarefree and \(n^2-qZ^2=1\) with \(0<Z\) and \(q\mid Z\), then the coprime positive factors of \(n^2-1=q^3(Z/q)^2\) are powerful. Here \(q\) is uniquely determined as the squarefree part of \(n^2-1\) and is greater than 1, since \((n-1)^2<n^2-1<n^2\). Reduction modulo 8 gives \(q\equiv7\pmod8\), since \(4\mid n\) and \(Z\) is odd. At fixed \(D\) there are at most \(3^{\omega(D)}\) assignments of its primes to the three kernels, each of which fixes \(q\). Each \(q\) gives \(O(\log(2X))\) traces up to \(X\) by the Pell parametrization below. Since \(D<X\) for a candidate, the count with fixed \(D\) is \(O_\varepsilon(X^\varepsilon)\).

Write \(\epsilon=T_1+U_1\sqrt q\) for the least positive Pell unit greater than 1 and \(\epsilon^k=T_k+U_k\sqrt q\) for \(k\ge1\). Every positive Pell solution is a power of \(\epsilon\): division by a suitable power puts the unit in \([1,\epsilon)\), and a norm-one unit in that interval has nonnegative integer coordinates, so it equals 1 by minimality. Modulo 2, a Pell unit is either 1 or \(\sqrt q\), the latter of order two. Thus a candidate \(T_k\) even requires \(T_1\) even, \(U_1\) odd and \(k\) odd. At each \(p\mid q\), expansion gives \(U_k\equiv kT_1^{k-1}U_1\pmod p\), with \(T_1^2\equiv1\pmod p\). Consequently \(q\mid U_k\) if and only if \(R=q/\gcd(q,U_1)\) divides \(k\).

For \(k>1\), we have \(k\ge3\) and \(\epsilon>2\sqrt q\,U_1\ge 2q^{3/2}/R\). If \(R=1\), use \(\epsilon^k\ge\epsilon^3\). Otherwise \(3\le R\le q\), and \(x\log(2q^{3/2}/x)\) is increasing on \([3,q]\), since its derivative is at least \(\log(2\sqrt q)-1>0\). Using \(k\ge R\) gives in both cases

\[
T_k>\tfrac12\epsilon^k>\frac4{27}q^{9/2}.
\]

Let \(F(X)\) count candidate centers with \(k=1\). These correspond one-to-one to squarefree \(q\equiv7\pmod8\) with \(q\mid U_1\), \(U_1\) odd and \(T_1\le X\) powerful. In particular a dyadic box with \(b_-\asymp B_-\), \(b_+\asymp B_+\) contains \(O(B_-B_+)\) fundamental candidates. The Pell characterization is classical. Reinhart's Definition 1.4 calls the \(k=1\) condition a strong counterexample and Example 3.1 gives a prime \(q\equiv7\pmod8\) with \(q\mid U_1\) and \(U_1\) odd. Sayim's Part I, Lemma I.1.2 and Theorem I.2.1, gives the same index divisibility and uses the same lower bound on \(U_1\) to separate heights. [Reinhart, *A counterexample to the Pellian equation conjecture of Mordell*](https://arxiv.org/pdf/2402.09827) · [Sayim, *Consecutive Powerful Numbers and Pell Equations*, Part I](https://doi.org/10.5281/zenodo.23127702)

For the repair, set \(t=T_1\). If \(p\parallel t-1\) or \(p\parallel t+1\), then \(p\) is odd and \(p\parallel t^2-1\), so \(p\mid q\), \(p\nmid U_1\), and \(p\mid k\) by the Pell congruence. For the center, the Chebyshev polynomial \(\mathcal T_k\) given by \(\mathcal T_k(\cos v)=\cos kv\) has integer coefficients by its recurrence, is odd for odd \(k\), and has derivative \((-1)^{(k-1)/2}k\) at zero. It satisfies \(T_k=\mathcal T_k(t)\) by the same recurrence for the Pell traces. Thus

\[
T_k=(-1)^{(k-1)/2}kt+t^3P_k(t),\qquad P_k\in\mathbb Z[t].
\]

If \(p\parallel t\) and \(p\nmid k\), it follows that \(p\parallel T_k\), impossible. This identity also gives \(v_2(t)=v_2(T_k)\) for the even \(t\) at issue, so \(t\ge4\). The exact products \(f_i=\prod_{p\parallel t+i}p\) are pairwise coprime and satisfy \(f_{-1}f_0f_1\mid\operatorname{rad}(k)\), with each \((t+i)/f_i\) powerful.

There is no extra kernel factor in counting \((t,k)\), since \(q\) is the squarefree part of \(t^2-1\). We have \(t<\epsilon<(2X)^{1/k}\) and \(k=O(\log(2X))\) since \(\epsilon>2\sqrt7\). The elementary estimate for powerful integers up to \(Y\) is \(O(Y^{1/2})\), by summing \((Y/b^3)^{1/2}\) over \(b\). Counting \(t=f_0w\) over \(f_0\mid\operatorname{rad}(k)\) therefore gives \(O(k(2X)^{1/(2k)})\) for each \(k\). The term \(k=3\) and the sum over \(k\ge5\) give the unconditional bound \(C(X)-F(X)\ll X^{1/6}\).

## Power-of-two centers and forced non-Wieferich primes

An even center \(b^k\), with \(b\ge2\) even and \(k\ge2\), has coprime neighbours \(b^k-1,b^k+1\). Both are powerful exactly when \(b^{2k}-1\) is powerful. Granville applied this observation to \(b=2\), showing that a finite set of base-two non-Wieferich primes would give a powerful triple. Keeping the order at each prime gives a more selective criterion. [Granville, pp. 215–217](https://dms.umontreal.ca/~andrew/PDF/PowerfulFLT.pdf)

**Proposition (Exact order criterion for a power center).**

For each odd prime \(p\nmid b\), put \(r_p=\operatorname{ord}_p(b)\) and \(d_p=r_p/\gcd(r_p,2)\). Call \(p\) base-\(b\) Wieferich if \(b^{p-1}\equiv1\bmod p^2\). The neighbours of \(b^k\) are both powerful if and only if \(p\mid k\) for every base-\(b\) non-Wieferich prime \(p\nmid b\) with \(d_p\mid k\). For a prime cutoff \(B\), start with \(M=1\) and process such primes in increasing order through \(B\), multiplying \(M\) by \(p\) when \(d_p\mid M\). The resulting \(M_B\) is the least positive integer satisfying the order criterion at all primes through \(B\), and divides every other such integer. For \(b=2\), the sequence \(M_B\) stabilizes if and only if there is a consecutive powerful triple with center a power of two.

Indeed \(p\mid b^{2k}-1\) if and only if \(d_p\mid k\). In that case, binomial lifting at the odd prime gives
\[
 v_p(b^{2k}-1)=v_p(b^{r_p}-1)+v_p(k)
             =v_p(b^{p-1}-1)+v_p(k).
\]
All primes dividing the product are odd and coprime to \(b\). This proves the criterion. Since \(d_p<p\), later primes cannot change whether \(d_p\) divides the product. At its step, if \(d_p\) divides the current product, every solution to the truncated criterion must contain \(p\). Otherwise the current product itself passes that step. Induction gives both minimality and feasibility for \(M_B\). If an exponent works at every prime, then every \(M_B\) divides it and the sequence stabilizes. If the sequence stabilizes, its limit passes the criterion at every prime. At \(b=2\) it is divisible by \(3\), forced at the first step, so the center itself is powerful.

**Theorem (Wieferich burden of a power-of-two counterexample).**

If \(2^k-1,2^k,2^k+1\) are all powerful, then \(k\) is divisible by an explicit product of \(795\) distinct odd primes that exceeds \(10^{6283}\). Moreover there are at least \(2^{795}\) distinct base-two Wieferich primes, all exceeding the greatest prime divisor of \(k\), and at least one base-two Wieferich prime exceeds \(2k\).

The [ordered sweep](verify_erdos_364_power_centers.cpp) checks all \(50{,}847{,}533\) odd primes through \(10^9\) using an odd least-prime-factor sieve. It marks \(789\) forced primes for \(b=2\), with largest \(999{,}375{,}343\). The current \(M\) is squarefree and odd. The test \(d_p\mid M\) is equivalent to \(b^{\gcd(2M,p-1)}\equiv1\bmod p\). The exponent for this test is computed as twice the product of the marked odd primes dividing \(p-1\). The modular test at \(p^2\) then recognizes the non-Wieferich primes among these candidates.

The [additional certificate](verify_erdos_364_power_center_certificate.py) uses six simple prime divisors of \(2^{2j}-1\) beyond that cutoff. The indices are \(j=43,127,129\) and all divide \(M_{10^9}\), since \(3,43,127\) appear in the sweep. Each additional prime exceeds \(j\), so \(2^{2j}-1\) having valuation one at that prime forces it to divide \(k\) by lifting. The checker gives the full factorization of \(p-1\) and an integer witness \(a\) for each such prime. It verifies \(a^{p-1}\equiv1\bmod p\) and \(\gcd(a^{(p-1)/\ell}-1,p)=1\) at each prime \(\ell\mid p-1\). For any prime divisor of \(p\) these tests force the order of \(a\) to be \(p-1\), proving primality. The script also verifies the two congruences at \(p,p^2\) and the bound on the product.

To obtain the lower bound on the number of Wieferich primes, let \(q\) be the greatest prime divisor of \(k\), so \(q\ge7\), and write \(w\) for the number of distinct odd prime factors of \(k\). For every divisor \(j\) of the odd part of \(k\) that is divisible by \(q\), use primitive primes of orders \(j\) and \(2j\). There are at least \(2^{w-1}\) such divisors, by taking \(q\) times a product of any subset of the other odd prime factors. Each order exceeds or equals \(q\), so its prime is greater than \(q\) and cannot divide \(k\). It nevertheless divides \(2^{2k}-1\). The order criterion therefore requires it to be base-two Wieferich. The orders \(j,2j\) are distinct as \(j\) varies over odd divisors. This supplies at least \(2^w\ge2^{795}\) such primes. A primitive prime of order \(2k\) supplies the last assertion.

Here the required primitive primes admit a short direct proof. For an index \(n\) whose greatest prime divisor \(q\) is at least \(7\), the cyclotomic value \(\Phi_n(2)\) has an odd prime divisor of order \(n\). To see this, for any odd prime \(\ell\) write \(f=\operatorname{ord}_\ell(2)\). If \(f\nmid n\), the cyclotomic valuation is zero. If \(n=fm\), then Möbius inversion of \(2^n-1=\prod_{d\mid n}\Phi_d(2)\) gives the valuation
\[
 \sum_{j\mid m}\mu(m/j)\big(v_\ell(2^f-1)+v_\ell(j)\big).
\]
For \(m=1\) the prime has order \(n\). For \(m>1\) the sum is one if \(m\) is a power of \(\ell\) and zero otherwise. In the former case \(\ell=q\), since \(f\mid\ell-1\). Thus without a prime of order \(n\), the odd number \(\Phi_n(2)\) would divide \(q\). On the other hand, the cyclotomic product gives
\[
 \Phi_n(2)\ \ge\ 2^{\varphi(n)}\prod_{d\ge1}(1-2^{-d})
           \ >\ 2^{\varphi(n)}/4\ \ge\ 2^{q-1}/4\ >\ q .
\]
The infinite product exceeds \(1/4\) by retaining its first two factors and using the bound \(1-\sum_{d\ge3}2^{-d}\) on the remaining product. This applies to all the indicated orders \(j,2j,2k\).

Replay both checks from the repository root with GCC 13.3.0 and Python 3.12.3.

~~~sh
g++ -std=c++20 -O3 research/arithmetic/verify_erdos_364_power_centers.cpp -o /tmp/erdos364-centers
/tmp/erdos364-centers 2 1000000000 > /tmp/erdos364-base2.csv
python3 research/arithmetic/verify_erdos_364_power_center_certificate.py /tmp/erdos364-base2.csv
~~~

The order criterion does not establish non-stabilization. A full exclusion of power-of-two centers, and the existence question for general centers in Erdős #364, remain open.

### Odd valuations at least three

**Theorem (Cubic Wieferich burden for power-of-two centers).**
Suppose \(2^k-1,2^k,2^k+1\) are powerful. Put
\(m=\operatorname{rad}(k/2^{v_2(k)})\), \(w=\omega(m)\), and let
\(q\) be the largest prime divisor of \(m\). There are at least
\(2^{w+1}-2w-1\) distinct primes \(\ell\nmid k\) with
\(v_\ell(2^{\ell-1}-1)\) odd and at least three. At least \(2^w-1\)
lie above \(q\), with at least \(2^{w-1}-1\) in each class
\(3,7\bmod8\), and one exceeds \(2m\). In particular the total
is at least \(2^{796}-1591\), with at least \(2^{795}-1\) above
\(2^{127}-1\). At the least admissible exponent \(k\), the
prime exceeding \(2m\) exceeds \(2k\).

For each divisor \(j>1\) of \(m\), consider the two squarefree indices
\(j,2j\). The Möbius product gives
\[
 \Phi_j(2)=\prod_{d\mid j}(2^d-1)^{\mu(j/d)}\equiv7\pmod8,\qquad
 \Phi_{2j}(2)=\prod_{d\mid j}(2^d+1)^{\mu(j/d)}\equiv3\pmod8 .
\]
For the first congruence \(d=1\) gives \(1\), every other \(d\)
gives \(-1\bmod8\), and \(\sum_{d\mid j}\mu(j/d)=0\) while
\(\mu(j)=\pm1\). For the second congruence use
\(\Phi_{2j}(X)=\Phi_j(-X)\). The \(d=1\) factor gives \(3\bmod8\)
and every other factor gives \(1\bmod8\).

Write \(f_\ell=\operatorname{ord}_\ell(2)\) for an odd prime \(\ell\).
The cyclotomic valuation formula above shows that at a squarefree
index \(t\), a prime \(\ell\) dividing \(\Phi_t(2)\) has either
\(t=f_\ell\) or \(t=\ell f_\ell\). In the latter case the prime
divides \(t\) and has valuation one. Discard every index in the
chosen set for which some prime of \(k\) divides \(\Phi_t(2)\).
Each odd prime of \(k\) discards at most two indices. The prime
\(3\) divides \(k\) by the order criterion and discards at most
one, since \(f_3=2\) and the index \(2\) is excluded. Of the
\(2^{w+1}-2\) chosen indices, at least \(2^{w+1}-2w-1\) remain.

At any retained index \(t\), all primes in \(\Phi_t(2)\) have
order \(t\) and lie outside \(k\). For each such prime, lifting
gives
\[
 v_\ell(\Phi_t(2))=v_\ell(2^t-1)
     =v_\ell(2^{2k}-1)=v_\ell(2^{\ell-1}-1),
\]
since \(t\) divides \(2k\) and \(\ell-1\), and \(\ell\nmid k\).
The product \(2^{2k}-1\) is powerful, so all of these exponents
are at least two. The cyclotomic value is a nonsquare modulo
\(8\), hence it has a prime with odd exponent at least three.
Primes obtained at distinct indices are distinct because
their orders are distinct. The condition \(v_\ell(2^{\ell-1}-1)\ge3\)
is the order-two Wieferich condition in Skula's terminology.
[Skula, Proposition 2.2](https://math.colgate.edu/~integers/t19/t19.pdf)

When \(\ell\) has order \(j\) with \(j\) odd, \(2\) is a square
modulo \(\ell\). The supplementary law for \(2\) gives
\(\ell\equiv1,7\bmod8\), so the factor \(\Phi_j(2)\equiv7\bmod8\)
requires a prime in the class \(7\) to odd exponent. At order
\(2j\), the order divides \((\ell-1)/2\) exactly when
\(\ell\equiv1\bmod4\). The supplementary law gives
\(\ell\equiv1,3\bmod8\) in the two cases, respectively. Thus
\(\Phi_{2j}(2)\equiv3\bmod8\) requires an odd exponent in the
class \(3\).

For the \(2^{w-1}\) choices of \(j\) divisible by \(q\), a prime
of order \(j\) or \(2j\) exceeds \(q\). The only possible prime
of \(k\) in those cyclotomic values is therefore the one at an
index \(t=q f_q\) with an imprimitive divisor \(q\). At most one
index is lost across both parities. In particular at least one
of \(m,2m\) is retained. An odd prime of either order satisfies
\(2m\mid \ell-1\). The previously certified exponent divisors
give \(w\ge795\) and \(q\ge2^{127}-1\). If any exponent works,
the least one is the odd squarefree limit of the order sweep,
so there \(m=k\).

The power-of-two existence question still requires proving
non-stabilization or finding an admissible exponent.

## Transfers from the preprint counting statements

Plymesser's [powerful-products paper](sources/plymesser-powerful-products-hexagon2611.00602-v1.0.1.pdf), Lemma 8.3, states that for every \(\varepsilon>0\) there is \(\eta>0\) such that fixed positive integer coefficients \(f_i\le Y^\eta\) give \(O_\varepsilon(Y^{\theta+\varepsilon})\) starts, where \(\theta=77/320\) and the three terms belong to the respective \(f_i\) times the powerful integers. For fixed \(k\) the possible coefficients from the repair form a finite set. Apply the lemma with start \(t-2\) and \(Y=(2X)^{1/k}\), increasing the threshold for \(X\) to fit \(Y^\eta\). Using it at \(k=3,5\), and the elementary middle-term bound for \(k\ge7\), gives

\[
C(X)=F(X)+O_\varepsilon(X^{77/960+\varepsilon}).
\]

The [quadratic companion](sources/plymesser-powerful-quadratics-hexagon2611.00699-v1.0.1.pdf), Theorem 1.1, states the exponent \(293/1220+\varepsilon\) for powerful values of a fixed positive-leading integer quadratic with nonsquare discriminant. For \(k=3\), the powerful integer \(T_3=t(4t^2-3)\) has factors with gcd dividing 3. If \(3\nmid t\), then \(4t^2-3\) is powerful. If \(t=3u\), then \(4t^2-3=3(12u^2-1)\), and \(12u^2-1\) is coprime to \(3t\), hence powerful. Both quadratics have discriminant 48. The bound at \(k=3\) becomes \(O_\varepsilon(X^{293/3660+\varepsilon})\), and the \(k=5\) use of Lemma 8.3 and the elementary \(k\ge7\) bound are smaller. This yields the stronger transfer

\[
C(X)=F(X)+O_\varepsilon(X^{293/3660+\varepsilon}).
\]

The [uniform small-coefficient estimate](ERDOS_137_UNIFORM_SMALL_COEFFICIENT_COUNTS.md#pell-repair-corollary) gives the sharper remainder \(O_\varepsilon(X^{(6/25-10^{-7})/3+\varepsilon})\) by applying the same repair at indices 3 and 5. It uses the fixed-shift pair inputs and the leading-form affine estimate.

The two exponents in this section are conditional on their respective Plymesser inputs. The coefficient and discriminant hypotheses are checked here. The [source audit](PLYMESSER_2026_POWERFUL_COUNTS_SOURCE_SCOPE.md) records which proof passages and formal hypotheses have been inspected.

## The only maximizer of the pairwise bound

Use \(x_1,x_2,x_3\in[0,1/3]\) for the dyadic exponents of the three squarefree kernels, measured at center height \(X\). For a pair \((a,b)\) of exponents put \(s=a+b\), \(h=1-3s/2\), \(d=|a-b|\), \(w=s-1/3\), \(m=(9/8)sh=3/16-(27/16)w^2\). Define \(\Phi(a,b)\) as the minimum of \(s,h\) and, on the gate \(s,h\ge1/10\), \(a,b\ge1/100\), also \(B_{\rm red}\). This is the \(\Phi_3\) of [Plymesser, §§6–7](sources/plymesser-powerful-products-hexagon2611.00602-v1.0.1.pdf). In this notation the exact formula for its third term is

\[
B_{\rm red}=\theta-\frac{81}{64}w^2
   +\frac3{20}(d-m)_+ +\frac1{10}(3d/2-m)_+,
\]

where \(v_+=\max(v,0)\). Indeed the two reduced-coordinate maxima are \((s-m)/2+(d-m)_+/2\) and \((h-m)/2+(3d/2-m)_+/2\), with weights \(3/10\) and \(1/5\).

For \(0\le\eta\le1/1500\), suppose every pair has \(\Phi\ge\theta-\eta\). Sort the three exponents as \(a\le b\le c\). If \(a<1/100\), then \(b>0.229\) and \(b+c<0.51\) by the bounds on pair sums from \(s,h\ge\theta-\eta\). Thus \((b,c)\) lies in the gate, with \(s>0.458\), \(h>0.239\), \(m>0.12\), and \(d<0.052\). Both positive parts vanish, contradicting \(B_{\rm red}\ge\theta-\eta\) since \(s>0.458\). Therefore all pairs lie in the gate.

Put \(\kappa=(8/9)\sqrt\eta\) and \(e=20\eta/3\). When \(|w|>\kappa\), at least one positive part must be active. Separating \(d\le m\) and \(d\ge m\) gives, respectively,

\[
d\ge\frac18+\frac{117}{16}w^2-e,\qquad
d\ge\frac5{32}+\frac{45}{16}w^2-\frac{10}3\eta.
\]

Both imply \(d\ge1/8+(45/16)w^2-e\). Write \(w_1=a+b-1/3\), \(w_2=b+c-1/3\) and \(L=c-a=w_2-w_1\). If both \(|w_j|>\kappa\), sum the two bounds to get \(L\ge1/4+(45/32)L^2-2e\). This is impossible since \(1/4+(45/32)L^2-L\ge1/4-8/45>2e\). If exactly one exceeds \(\kappa\), the corresponding difference \(d\) is at most \(L\le |w|+\kappa\). This is impossible because \(1/8+(45/16)w^2-|w|\ge1/8-4/45>e+\kappa\). Hence both adjacent sums are within \(\kappa\) of \(1/3\). The inequalities \((a+b)/2\le b\le(b+c)/2\) show that all three exponents are within \(3\kappa/2=(4/3)\sqrt\eta\) of \(1/6\).

At \(\eta=0\) this proves uniqueness, and at \(a=b=c=1/6\) the value is exactly \(\theta\). In particular, for any \(0<\rho\le1/100\), the minimum of the pairwise exponents is at most \(\theta-\rho^2/4\) outside the box \(|x_i-1/6|\le\rho\). If one uses the pairwise counting estimates of [Plymesser, Lemmas 4.1 and 6.4 and Theorem 6.3](sources/plymesser-powerful-products-hexagon2611.00602-v1.0.1.pdf), this also bounds the count in the exterior by \(O_{\rho,\varepsilon}(X^{\theta-\rho^2/4+\varepsilon})\). The shifts of the three pairs are 1 or 2, their two coordinate gcds are 1, and the dyadic windows in §8.1 use only an \(O(1/\log X)\) perturbation of the exponents. Each solution to a selected pair determines the center, so the pair count bounds the triple count. Summation over dyadic boxes of kernels and centers preserves the exponent with the displayed \(\varepsilon\).

## Transferring the contact determinant to the pair count

The fixed-shift equation for a primitive pair is \(e^3v^2-d^3u^2=h\), with \(h\) fixed and nonzero. In the windows of the powerful-products paper, \(d,e,u,v\) have sizes \(D,E,U,V\) with \(D^3U^2=E^3V^2=X\). Put \(c=D/E\), \(r=(d/e)/c\) and \(z=(v/u)/c^{3/2}\). Then \(r\in[1/2,2]\) and

\[
z=r^{3/2}\bigl(1+h/(d^3u^2)\bigr)^{1/2}=r^{3/2}+O_h(X^{-1}).
\]

Fix \(57/305<\mu<1/5\). Partition into the multiplicative slope intervals of §5 of the powerful-products paper with \(M=\lceil X^\mu\rceil\). In normalized coordinates each interval has width \(O(X^{-\mu})\). At its left endpoint \(r_0\), put \(w=r-r_0\) and \(v'=z-(r_0+w)^{3/2}\). Their Taylor weights are \(\mu,1\). The Taylor coefficients and remainders for the latter graph are uniform over \(r_0\in[1/2,2]\).

[Plymesser's quadratic paper, Proposition 5.7](sources/plymesser-powerful-quadratics-hexagon2611.00699-v1.0.1.pdf) is a statement about any real polynomial graph and any basis of a rectangular polynomial space. For sufficiently large fixed \(J\), depending on \(\mu\), it gives a total weight strictly greater than \(JN\) for any nonsingular Taylor coefficient minor on the rectangle of bidegree \((3J,2J)\), where \(N=(3J+1)(2J+1)\). Apply it to the Taylor polynomial of \((r_0+w)^{3/2}\) and the basis \((r_0+w)^i z^j\). As in §5.3 of that paper, truncate at a fixed order for which the error is of weight greater than \(JN\). Expansion of the evaluation determinant then gives \(O(X^{-JN-\gamma})\) for some fixed \(\gamma>0\). Only minors with nonsingular Taylor coefficient matrices contribute to the finite main sum.

When clearing denominators, the column indexed by \((i,j)\) is multiplied by \(c^{i+3j/2}\), whose product over all columns is \(c^{3JN}\). Row clearing costs at most \(O((E^{3J}U^{2J})^N)\). The exact aggregate scale is therefore \(O((D^{3J}U^{2J})^N)=O(X^{JN})\), even for unequal \(D,E\). The cleared evaluation is an integer determinant of size \(O(X^{-\gamma})\), so it vanishes. Each original integer entry is bounded by \(\max(d,e)^{3J}\max(u,v)^{2J}\ll_J X^{2J}\) from the two products of size \(X\). The cofactors of a maximal-rank submatrix therefore produce a bounded-bidegree integer auxiliary form of coefficient size \(X^{O_\mu(1)}\).

The per-interval component proof of [the powerful-products paper, Proposition 5.6](sources/plymesser-powerful-products-hexagon2611.00602-v1.0.1.pdf) uses the determinant only for the existence of such a form. Its lattice strips require \(X\ge3|h|M\), which holds here for fixed \(h\). Its low-degree component bound allows coefficient height \(X^{O_\mu(1)}\), and its lifted-curve bound is uniform in the coefficients. The tail estimates of Lemmas 6.1–6.2 then use the actual \(M\) to sum the same per-interval estimate. Consequently the box exponent from Theorem 6.3 is now

\[
\mu+\frac3{10}\max\left(\frac{s-\mu}{2},\,a-\mu,\,b-\mu\right)
+\frac15\max\left(\frac{h_0-\mu}{2},\,\frac{1-3a}{2}-\mu,\,\frac{1-3b}{2}-\mu\right),
\qquad h_0=1-3s/2.
\]

Here \(h_0\) is the sum of the two square-part exponents. If \(|a-b|\le2\mu/3\), both first terms in the maxima are active and the exponent is \(1/10+3\mu/4\). In the context of the triple, boxes with old pairwise exponent at most \(\theta-1/2000\) already have an exponent below \(293/1220\). In every other box the stability proof with \(\eta=1/2000\) gives \(\max_i x_i-\min_i x_i\le2\kappa<0.04\). Thus any pair there has \(|a-b|\le2\mu/3\), even after the window perturbation for large \(X\). Taking \(\mu\) sufficiently close to \(57/305\) for a given \(\varepsilon>0\), and summing the windows, gives

\[
F(X)\le C(X)\ll_\varepsilon X^{293/1220+\varepsilon},\qquad
\frac{77}{320}-\frac{293}{1220}=\frac9{19520}.
\]

This transfer uses the fixed-shift pair estimates and the component and lattice bounds of the first paper, and the contact-weight Proposition 5.7 of the second. It does not apply the theorem on powerful quadratic values to \(n^2-1\), whose discriminant is a square. The deduction is written here with the source estimates as inputs. The relevant written proof passages were inspected, with no replay of either formal development.

## A joint determinant in all dyadic boxes

Work at center height \(X\le n\le 2X\), with each kernel in a dyadic interval \([B_i,2B_i)\) and \(B_i\ge1\). Put \(A_i=X^{1/2}/B_i^{3/2}\) and normalize the slopes by

\[
r_\pm=\frac{b_\pm/B_\pm}{b_0/B_0},\qquad
s_\pm=\frac{a_\pm/A_\pm}{a_0/A_0},\qquad t=1/n.
\]

The two identities are \(r_\pm^3s_\pm^2=1\pm t\), so \(s_\pm=r_\pm^{-3/2}(1\pm t)^{1/2}\). Partition \((r_-,r_+)\in[1/2,2]^2\) into rectangles of side lengths \(X^{-\mu_-},X^{-\mu_+}\), with fixed positive \(\mu_-,\mu_+\). In each rectangle expand at its two left endpoints and \(t=0\). One expansion variable records the single displacement, with weight 1. The two slope variables have weights \(\mu_-\) and \(\mu_+\).

Write \(V\subset\mathbb P_b^2\times\mathbb P_a^2\) for the hypersurface with equation

\[
G_0=b_-^3a_-^2-2b_0^3a_0^2+b_+^3a_+^2=0.
\]

It is irreducible. In the chart \(b_0a_0\ne0\) with unscaled ratios, solving for \((a_+/a_0)^2\) gives \((2-r_-^3s_-^2)/r_+^3\) in the rational function field of the other ratios. This has odd valuation along \(r_+=0\) and is not a square. The quadratic is primitive over the polynomial ring in the other ratios. No homogeneous factor can be supported on the excluded coordinate hyperplanes, since neither \(b_0\) nor \(a_0\) divides \(G_0\).

Take bidegrees \((3L,2L)\) for a positive integer \(L\). Let \(N_L=\binom{3L+2}{2}\binom{2L+2}{2}\). Choose the quotient basis by omitting multiples of the leading monomial of \(G_0\). Its size is

\[
H=N_L-N_{L-1}=36L^3+O(L^2).
\]

Factor out from each row \((b_0/B_0)^{3L}(a_0/A_0)^{2L}\), bounded uniformly at fixed \(L\). The remaining scale of the column for a monomial \(m=b^{\mathbf j}a^{\mathbf k}\) is \(\sigma(m)=\prod_i B_i^{j_i}A_i^{k_i}\). Over all monomials of bidegree \((3L,2L)\) the aggregate exponents of each \(B_i,A_i\) are \(LN_L, (2/3)LN_L\) by symmetry. Since \(B_i^3A_i^2=X\), their scale product is \(X^{LN_L}\). The omitted monomials are the leading monomial of \(G_0\), of scale \(X\), times the entire rectangle of bidegree \((3L-3,2L-2)\). The omitted scale product is \(X^{LN_{L-1}}\). Thus over the quotient basis

\[
\prod_m \sigma(m)=X^{LH}.
\]

In any nonempty box \(B_i^3\le2X+1\), so every \(A_i\) is bounded below by an absolute positive constant. Together with \(B_i\ge1\) this gives \(\sigma(m)\gg_L1\) and bounds the scale product of any column subset by \(O_L(X^{LH})\). After removing row and column factors, the entries are analytic functions of the two scaled slopes and \(t\) on a fixed neighborhood of the expansion centers, including when one of the kernel intervals is bounded. Any \(h\)-row, \(h\)-column minor in the original integer matrix is therefore \(O(X^{LH-W(h)})\), where \(W(h)\) sums the \(h\) smallest values of \(\mu_-i+\mu_+j+k\) for \(i,j,k\ge0\), with multiplicity. Indeed, Taylor expand each normalized column to a fixed total degree high enough that every remainder is \(O(X^{-W(h)})\). A determinant term with repeated Taylor monomials vanishes by multilinearity, and distinct monomials have total weight at least \(W(h)\). The constants are uniform over all boxes and cells at fixed degrees and widths.

The number of weights at most \(v\) is \(v^3/(6\mu_-\mu_+)+O_{\mu_-,\mu_+}(v^2+1)\), by counting unit cubes against the weighted simplex. Inverting this estimate for the ordered weights and summing gives

\[
\frac{W(h)}h=\frac34(6h\mu_-\mu_+)^{1/3}+O_{\mu_-,\mu_+}(1).
\]

Since \(H=36L^3+O(L^2)\), the condition

\[
\mu_-\mu_+>\frac8{729}
\]

gives \(W(H)>LH\) for sufficiently large fixed \(L\). It also gives \(W(\lceil\lambda H\rceil)>LH\) for a fixed \(\lambda<1\) sufficiently close to 1, since the ratio \(W(\lceil\lambda H\rceil)/(LH)\) tends to \((9/2)(\mu_-\mu_+)^{1/3}\lambda^{4/3}\). All these minors vanish by integrality for large \(X\), so the evaluation rank is less than \(\lceil\lambda H\rceil\). Its kernel gives relations linearly independent modulo \(G_0\), with integer coefficient vectors of size \(X^{O_L(1)}\) obtained from minors of a maximal-rank integer submatrix.

In each nonempty cell, one such relation intersects the irreducible \(V\) properly, in dimension at most two. The components have bounded degree and number for fixed \(L\), by Bézout in the Segre embedding after multiplying the bihomogeneous equations by powers of \(a_0,b_0\) to give equal bidegrees. Those coordinate hyperplanes are outside the positive chart. The number of rectangles is \(O(X^{\mu_-+\mu_+})\). Taking equal widths and summing over the dyadic boxes gives a surface cover with \(O_\varepsilon(X^{4\sqrt2/27+\varepsilon})\) surfaces of bounded degree, uniformly for all kernel sizes.

Using just the products of the two slopes gives \((r_-r_+)^3(s_-s_+)^2=1-t^2\). For that two-variable determinant the Taylor weights are \(\mu,2\), but the two product denominators have twice the original exponents. In the symmetric box, bidegree \((B,A)\) in these product coordinates has \(H\sim BA\) and cost \(2(B/6+A/4)\). The leading average Taylor weight is \((4/3)\sqrt{\mu H}\). Since \(BA\le6(B/6+A/4)^2\), the rectangular monomial bound then requires \(\mu>3/8\) for a leading-order strict margin. The joint cover retains both slopes with the shared displacement.

## Second determinant for surfaces with a curve projection

Consider a geometrically integral surface in \(V\) among the joint auxiliary components, with a curve image of degree \(e\) in one of the blocks \(b\), \(a\), or \(c_i=b_i a_i\). Work on the positive chart. If a geometric component is not defined over \(\mathbb Q\), its rational points lie on its intersection with a distinct Galois conjugate, a bounded-degree locus of dimension at most one. We may therefore consider the surface and its image over \(\mathbb Q\).

The surface is the torus part of the pullback of the image curve on \(V\), with its projective closure understood. Indeed the generic fiber over a torus point is a smooth conic, a smooth plane cubic, or an open in a line, respectively. For the third statement substitute \(a_i=c_i/b_i\) in \(G_0\), giving a linear equation in the \(b_i\) with coefficients proportional to \(c_i^2\). The fibers over boundary points of the \(b\) and \(a\) projections have dimension at most one, since the nonzero diagonal quadratic or cubic never vanishes on a whole fiber plane. In the third projection, the indeterminacy locus has disjoint supports of \(a,b\) and dimension at most one. A fiber with exactly two nonzero \(c_i\) fixes the ratios of both participating \(b_i\) and \(a_i\) and leaves at most one parameter. Exactly one nonzero product is incompatible with \(G_0=0\). The image curve meets the coordinate boundary in finitely many points, so no additional surface is supported there. The pullback is generically reduced along the surface, by the smooth generic fiber over the reduced image curve.

In the ambient product, let \(H_b,H_a\) be the hyperplane classes and \(L_0=3H_b+2H_a\). Each surface is consequently the unique dimension-two component of the complete intersection of \(V\) with the projection equation of class \(eH_b,eH_a\), or \(e(H_b+H_a)\). Its Hilbert function in bidegree \((3L,2L)\) has leading coefficient \(L_0^3\) times that projection class, divided by two. The three leading terms are \(18eL^2,27eL^2,45eL^2\), with \(O_e(L)\) error. This follows also from the bigraded complete intersection exact sequence; the ambient dimensions for the four terms are products of binomial coefficients. The coordinate boundary has smaller dimension, so it does not affect the leading term.

In a cell of width \(X^{-\mu}\), use the two normalized ratios of the projected block and \(t=1/n\) as analytic coordinates. For the \(a\) block invert by \(r_\pm=(1\pm t)^{1/3}s_\pm^{-2/3}\). For \(c\), write \(v_\pm=r_\pm s_\pm=r_\pm^{-1/2}(1\pm t)^{1/2}\), so \(r_\pm=(1\pm t)v_\pm^{-2}\) and \(s_\pm=v_\pm^3/(1\pm t)\). Each inverted chart has uniformly bounded derivatives near the positive ranges, and the two projected ratios vary by \(O(X^{-\mu})\). The Taylor weights are \((\mu,\mu,1)\).

For each power \(t^j\), the restrictions of polynomials of total degree at most \(k\) in the two centered plane ratios span at most \(ek+O_e(1)\) dimensions, by reduction modulo the image equation. This bound uses no coefficient size for that equation. Thus the cumulative rank of Taylor monomial evaluations through weight \(w\) is at most \(e w^2/(2\mu)+O_{e,\mu}(w+1)\), by summing over \(j\). A nonzero term of the evaluation determinant must select independent evaluations. For a basis of \(H\) restricted sections its required weight sum is at least \((2/3)\sqrt{2\mu/e}\,H^{3/2}-O_{e,\mu}(H)\). The infimum for the leading mean weight, divided by \(L\), is consequently bounded below by \(4\sqrt\mu,\sqrt{24\mu}\), or \(\sqrt{40\mu}\), according to the block.

Select the basis of restricted sections from the original monomials. In the near-binding region the row cost per evaluation is at most \(O(X^{17L/15})\): the exponent divided by \(L\) is \(3\max\beta_i+2\max\alpha_i=1+3(\max\beta_i-\min\beta_i)\le17/15\). At \(\mu=11/105\) all three Taylor bounds strictly exceed this cost. Taylor truncation at a sufficiently high fixed order and integer determinant vanishing produce a bounded-degree form that restricts nontrivially to the surface. Its zeros cut that surface into boundedly many curves in the same cell by Bézout. Thus the shared surfaces still requiring a surface count are those dominating all three blocks. This cut leaves the arithmetic count on the resulting curves as a separate step.

After removing constant-block curves by pair-ratio separation, a positive-genus joint curve whose \(\mathcal O(3,2)\)-degree is below 12 must have pullback degrees \((2,2)\). Each image is a line, because a degree-one map to a conic would give genus zero. In the product of those lines the \((2,2)\) curve must be a factor of the restricted \((3,2)\) progression identity. The restriction is nonzero, since at a torus \(b\)-point a smooth diagonal conic cannot contain the whole \(a\)-line. The remaining \((1,0)\) factor gives a fixed \(b\)-point where the diagonal conic contains that line. Rank one forces a coordinate line; rank two forces a constant ratio between two participating \(a_i\), and rank three admits no line. The rank-two ratio must be rational to contain any positive rational point and is covered by pair-ratio separation. Hence among the curves left by these joint cuts, only positive-genus lift degrees 12–21 can still require a new curve estimate.

## Point count required on the surfaces

Linearly independent auxiliary forms modulo \(G_0\) may have a common divisor on \(V\). The rank bound leaves open whether the base locus is a union of curves, and a curve cover would still require a uniform arithmetic point count including exceptional components. The coefficient bound in the determinant proof is polynomial in \(X\) with a degree-dependent exponent, so a polynomial dependence on the coefficients in a point-counting theorem cannot be absorbed into \(X^\varepsilon\).

Let \(\Psi=\min_{i<j}\Phi(\beta_i,\beta_j)\), using the clipped kernel exponents \(\beta_i\) in \([0,1/3]\) for the old pairwise estimates. Clipping and changes of window affect the actual exponents by \(O(1/\log X)\). With \(\eta=1/1600\), the stability proof gives the following exact region:

\[
\Psi\ge \frac6{25}\quad\Longleftrightarrow\quad
 \frac{14}{45}\le\beta_i+\beta_j\le\frac{16}{45}\quad(i<j).
\]

In the forward direction the two adjacent sums in sorted order lie within \(\kappa=1/45\) of \(1/3\), hence so does the third. Conversely, the pair-sum conditions give span at most \(2/45\), every exponent within \(1/30\) of \(1/6\), and \(m\ge14/75\) for all pairs. The gate holds, both positive parts vanish, and \(B_{\rm red}\ge\theta-81/(64\cdot45^2)=6/25\). The other two terms \(s,h\) in \(\Phi\) are strictly greater than \(6/25\) throughout this region. Consequently \(\Psi>6/25\) corresponds exactly to all three pair sums in the open interval \((14/45,16/45)\). The old pairwise estimates already handle equality, with the window perturbations absorbed into \(X^\varepsilon\). For these boxes \(q\gg X^{14/45}\), so the Pell height inequality makes every candidate fundamental for sufficiently large \(X\). Also, at fixed kernel vector, the residues \(n\equiv-i\pmod{b_i^3}\) pin down \(n\) in the height band. Indeed their combined modulus is \((b_-b_0b_+)^3\gg X^{7/5}\) by summing the three pair-sum inequalities and using pairwise coprimality.

Take \(\mu_-=\mu_+=\mu=11/105\), which satisfies \(\mu^2-8/729=1/893025>0\). The strict margin allows a sufficiently large fixed auxiliary degree. There are \(K=O(X^{22/105})\) cells per box. In the boxes with \(\Psi>6/25\), remove the union of the candidate points already counted by the component estimates below, applying each estimate to a bounded number of its components per cell. Their total is \(O_\varepsilon(X^{6/25+\varepsilon})\). Let \(N_R\) count the remaining candidates in cell \(R\). A bound of \(O_\varepsilon(X^{16/525+\varepsilon})\) per cell, or for the mean over the cells of the box, suffices for a joint-cover proof of the \(6/25\) count. The independent averaged pair proof above already attains that exponent.

An alternative sufficient bound counts the ordered distinct collisions between those remaining points:

\[
P(X)=\sum_R N_R(N_R-1)\ll_\varepsilon X^{142/525+\varepsilon}.
\]

Indeed \((\sum_R N_R)^2\le K(\sum_R N_R+P(X))\), so \(\sum_R N_R\ll K+\sqrt{KP(X)}\). The exponent of the square-root term is \((22/105+142/525)/2=6/25\). For a collision of two candidates \(n,n'\), define \(u_i=b_i b'_0\), \(u'_i=b'_i b_0\), \(v_i=a_i a'_0\), \(v'_i=a'_i a_0\) for \(i=\pm1\). Then

\[
u_i^3v_i^2-u_i'^3v_i'^2=i(n'-n),\qquad
|u_i-u'_i|\ll B_iB_0X^{-\mu_i},\qquad |v_i-v'_i|\ll A_iA_0X^{-\mu_i}.
\]

The difference bounds follow from the cell intervals on \(r_i\) and their induced intervals on \(s_i\), since \(t=O(X^{-1})\). Either sign's two products determine the ordered centers: their difference gives \(n'-n\), and their product formula \(u_i^3v_i^2=(n+i)n'\) is strictly increasing in \(n\) for fixed \(n'-n\) on the positive domain. Neither difference in the size bounds vanishes for distinct candidates, by the pair-ratio separation below. The open collision count concerns the actual packed products with their common middle factors and the original squarefreeness and primitivity conditions.

A shared geometrically irreducible surface whose two images under \((r_-,s_-)\) and \((r_+,s_+)\) both have dimension at most one has constant \(t\), so it contributes at most one center. If \(t\) were nonconstant, the two fields generated by the images, each containing \(t\), would both be algebraic over the one-variable field generated by \(t\) over the algebraic closure of the constants. Their compositum is the field of the surface in this chart, of transcendence degree two, a contradiction. Thus the surfaces that require a new count dominate at least one of these pair projections.

### Candidate application of fixed-quadratic square-product counts

For a fixed separable \(P\in\mathbb Z[T]\) of degree two, de la
Bretèche, Wang and Xu bound the ordered distinct-input pairs with
nonzero square product by \(O_P(N^{1/2})\) in Proposition 5.1.
Proposition 6.1 bounds the pairwise distinct four-input count by
\(O_{P,\varepsilon}(N^{3/2+\varepsilon})\). Lemma 6.2 tracks a
specified divisor of the squarefree kernel of values of the same
fixed quadratic. [arXiv:2607.06398v1](https://arxiv.org/html/2607.06398v1)

Set \(h=n'-n\). For \(i=\pm1\) the two packed values in the collision
are \(P_{h-i}(n+i)\) and \(P_{h+i}(n)\), where \(P_a(T)=T(T+a)\).
In one dyadic shift band and one cell, test whether the Pell
and squarefree-kernel divisor arguments can handle these values.
A direct use of the pair or quadruple count would also require
encoding the relevant configurations as nonzero square products of
values of a single separable quadratic, with controlled fibres.
Powerfulness alone supplies no squarefree-kernel match. Treat
\(h=\pm1\) separately where the polynomial has a double root.
The polynomials and discriminants vary with \(h\), whereas the
source's constants fix \(P\). Track that dependence and keep the
opposite shifts, common middle factors and cell widths during
summation. The required total over all ordered survivor collisions
in a box is \(O_\varepsilon(X^{142/525+\varepsilon})\).

## A divisibility bound on shared surfaces

Fix one index \(i\) and let \(j,k\) be the other two. Use the integer coordinates
\[
x=n+i,\qquad w_1=b_i^2b_j a_i a_j,\qquad w_2=b_i^2b_k a_i a_k.
\]
They are evaluations of sections of bidegree \((3,2)\) on the affine normalization by the shift 1. Suppose the component's points obey a nonzero integer polynomial \(P(x,w_1,w_2)\) of fixed bounded degree and coefficient height \(X^{O(1)}\). If one monomial of its support uniquely minimizes both weights \((2,1,1)\) and \((3,2,2)\), with corresponding coefficient \(c\), the component contributes \(O_\varepsilon(X^\varepsilon)\) centers. This includes the case of nonzero constant term.

To see this, take any prime \(p\mid a_i b_i\), with \(e=v_p(a_i)\) and \(f=v_p(b_i)\). Coprimality across the three terms gives the coordinate valuations \(e(2,1,1)+f(3,2,2)\). The monomial singled out has valuation, before its coefficient, at least \(e+f\) below every other monomial. The zero polynomial value forces \(p^{e+f}\mid c\). Hence \(a_i b_i\mid c\), leaving \(O_\varepsilon(X^\varepsilon)\) divisor pairs and at most one center for each. The condition is on the actual support of an equation satisfied by the component. The general determinant does not ensure such an equation.

## Pair-ratio separation and binomial components

Fix distinct indices \(i,j\). For two centers \(n,n'\) in a height band \([X,2X]\), put \(u=b_i b'_j\), \(u'=b'_i b_j\), \(v=a_i a'_j\), \(v'=a'_i a_j\). The same product identity gives

\[
u^3v^2-u'^3v'^2=(i-j)(n'-n),\qquad
u^3v^2,u'^3v'^2\ge (X-1)^2.
\]

If \(u=u'\) and \(v\ne v'\), the left side has absolute value at least \(u^3(v+v')\ge2u^{3/2}(X-1)\). Equal kernel ratios give the same pair of kernels by primitivity, so \(u=b_i b_j\ge2\) since distinct positive squares have gap at least three. This lower bound exceeds \(2X\) for \(X\ge2\). If \(v=v'\) and \(u\ne u'\), the left side has absolute value at least \(3v^2\min(u,u')^2\ge3v^{2/3}(X-1)^{4/3}>2X\) for large \(X\). If both products agree the centers coincide. Thus a fixed rational ratio in either coordinate block admits at most one candidate center per sufficiently large height band, in any box. This gives the fixed-pair count without a separate Pell or Thue bound.

Write \(x=(x_-,x_0,x_+)\) for either the kernel vector or the square-part vector. In a nontrivial homogeneous binomial constraint \(\prod_i x_i^{m_i}=c\) with fixed bounded integer exponents summing to zero, at least two coordinates are active. If a point exists, \(c=A/B\) is positive rational with \(A,B\) coprime and of polynomial height in \(X\), by evaluation at that point. Pairwise coprimality forces every prime of an active coordinate to divide \(AB\). Each such prime can occur at at most one active coordinate, and its valuation is fixed once that index is specified by \(m_i v_p(x_i)=v_p(c)\). There are at most \(3^{\omega(AB)}\ll_\varepsilon X^\varepsilon\) choices for the active coordinates, each fixing at most one center by pair-ratio separation.

More generally, in either block suppose an integer equation of bounded degree and polynomial coefficient height has a unique monomial of least exponent in coordinate \(x_i\), with coefficient \(c_i\). For each \(p^e\parallel x_i\) every other coordinate is a unit modulo \(p\). After factoring out the least monomial valuation, all other terms have an extra factor \(p^e\). Thus \(x_i\mid c_i\). If this holds for two distinct coordinates, whose minimizing monomials may differ, divisor counting and the fixed-pair result again give \(O_\varepsilon(X^\varepsilon)\) centers on the component.

## Averaging the reduced lattice heights

For equal width exponent \(\mu\), let \(M=\lceil X^\mu\rceil\). For each cell with lower endpoints \(r_-^0,r_+^0\), form the rank-three lattice given by mapping \(b\in\mathbb Z^3\) to

\[
\left(\frac{b_0}{B_0},\ M\left(\frac{b_-}{B_-}-r_-^0\frac{b_0}{B_0}\right),\ M\left(\frac{b_+}{B_+}-r_+^0\frac{b_0}{B_0}\right)\right).
\]

Use the analogous lattice for \(a\), with \(A_i\) and slopes \((r_i^0)^{-3/2}\). The normalized lattice images of genuine candidates are in a fixed cube, since \(1/n\ll1/M\). Let \(L_b,L_a\) be the reciprocals of the first sup-norm minima, capped below by 1. For \(T>1\) and \(T\ge C\max B_i/M\) with a sufficiently large absolute \(C\), the number of cells with \(L_b\ge T\) is \(O(B_0B_-B_+/T^3)\). A short vector has \(q_0=|b_0|\) between 1 and \(B_0/T\): if \(b_0=0\), each of the other integer coordinates must have absolute value less than 1. Reverse the sign so that \(b_0=q_0\). The other two coordinates are positive and comparable to \(q_0 B_i/B_0\), giving \(O(B_-B_+q_0^2/B_0^2)\) choices. Each such vector constrains each interval endpoint to a range of length \(O(B_0/(MTq_0))\), hence meets \(O((B_0/(Tq_0))^2)\) cells. Summing over \(q_0\) proves the tail bound. The same argument applies on the \(a\) side because the endpoint map \(r\mapsto r^{-3/2}\) preserves spacing up to constant factors on \([1/2,2]\).

For \(S=B,A\), put \(L_{*,S}=\max(1,\max S_i/M,(S_0S_-S_+/M^2)^{1/3})\). Summing the tails in dyadic ranges gives \(\sum_R L_b^p\ll_p M^2 L_{*,B}^{p}\) for \(0<p<3\), and the analogous result for \(a\). Hölder yields \(\sum_R L_b^p L_a^q\ll_{p,q}M^2 L_{*,B}^p L_{*,A}^q\) for \(p,q\ge0\), \(p+q<3\), including an exponent zero by the separate moment bound. In the localized region write \(\alpha_i=(1-3\beta_i)/2\) and let bars denote coordinate averages. The spans imply \(\max\beta_i-\bar\beta\le4/135\) and \(\max\alpha_i-\bar\alpha\le2/45\). For \(4/45\le\mu\le1/5\), the reference heights satisfy \(L_{*,B}\ll X^{\bar\beta-2\mu/3}\) and \(L_{*,A}\ll X^{\bar\alpha-2\mu/3+(2/45-\mu/3)_+}\), up to fixed factors from the window conventions. The exponents are positive, and \(3\bar\beta+2\bar\alpha=1\) gives

\[
L_{*,B}^3L_{*,A}^2\ll X^{\max\{1-10\mu/3,\ 49/45-4\mu\}}.
\]

For the arithmetic lift, choose three independent vectors at successive minima in each lattice. Their index is bounded by the upper bound in Minkowski's second theorem and Hadamard's determinant inequality. Cramer's rule and the same bound on the product of the minima give coordinates \(O(L_b)\) and \(O(L_a)\) for vectors in the fixed cube. Multiplication by the respective index gives integer coordinates with the same bounds. This is an invertible rational linear change of the original integer coordinates on each side. Take its 60 monomials of bidegree \((3,2)\). They are integers of height \(O(L_b^3L_a^2)\), and the two fixed differences become two affine linear equations. On \(V\) these equations give the affine chart of its Segre–Veronese image where \(b_0^3a_0^2-b_-^3a_-^2\ne0\), scaled so this linear form is 1. Each lifted point determines its center by the linear form for \(b_0^3a_0^2\).

If a cell contains a bounded number of integral curve components of bounded degrees, let \(d=3d_b+2d_a\) denote the degree in this embedding. Here \(d_b,d_a\) are the degrees of the pullbacks of the hyperplane classes to the curve, so they include the degree of a nonconstant projection onto its image. The affine chart preserves \(d\), and [Castryck–Cluckers–Dittmann–Nguyen, Theorem 3](sources/castryck-cluckers-dittmann-nguyen-dimension-growth-2020.pdf) bounds each component by \(O_{d,\varepsilon}(X^\varepsilon(L_b^3L_a^2)^{1/d})\). The lifts have polynomial height even in the exceptional cells. For \(d\ge2\) the relevant mixed moment has total order \(5/d<3\). Thus summing over cells gives exponent

\[
E_d(\mu)=2\mu+\frac1d\max\left\{1-\frac{10\mu}{3},\ \frac{49}{45}-4\mu\right\}.
\]

The unequal-window term comes from the square-part reference height. The triple \((2/15,8/45,8/45)\), approached from the strict region, attains \(\max\alpha_i-\bar\alpha=2/45\). At \(\mu=11/105\), the height exponent is \(211/315\) and

\[
E_{22}(11/105)=\frac{1663}{6930}=\frac6{25}-\frac1{34650}.
\]

This handles degree at least 22 in the sum over cells. For degree 21 the unequal-window term alone would require \(\mu\le889/8550<2\sqrt{2}/27\) to reach \(6/25\). The same source's Theorem 2 separately handles surface components whose projection to \(\mathbb P_b^2\) is an integral curve of bounded image degree at least 8. The CRT gives at most one center per projected point, the lattice coordinates on that side give projective height \(O(L_b)\), and the summed exponent is at most \(22/105+(2/8)(8/45-22/315)=149/630<6/25\). These deductions apply to the indicated components of bounded multiplicity per cell.

## Rational curves with controlled parametrizations

Suppose a rational curve in \(V\) is given by two triples \(P_i,Q_i\) of integral homogeneous binary forms, with one fixed degree per triple, no common zero within either triple, and coefficient heights at most \(X^K\) for fixed \(K\). Count only the candidate points with rational preimages under the supplied map. A birational normalization covers all smooth rational points, and singular points can be treated separately. Put \(M_i=P_i^3Q_i^2\) and \(G=M_0-M_-=M_+-M_0\). If the curve carries a genuine parametrized point, \(G\not\equiv0\). The function \(f=M_0/G\) is the center at those points, so a constant \(f\) contributes at most one center.

If \(f\) has degree \(d>0\), let \(Z\) be the geometric projective roots of \(G\), let \(p\) count the poles of \(f\), and let \(e\) count points of \(Z\) over \(-1,0,1\). All poles are in \(Z\), and \(p+e\le|Z|\). The multiplicity of every point outside \(Z\) in the three finite fibers is at least two: \(f+i=M_i/G\) and every zero of \(M_i\) has multiplicity at least two. These fibers together contribute at least \(3d/2-e/2\) to the ramification, and the poles contribute \(d-p\). Riemann–Hurwitz for \(f:\mathbb P^1\to\mathbb P^1\) gives \(d\le 2p+e-4\), hence \(|Z|\ge3\).

At a primitive integral parameter \((s,t)\), the gcd \(g_P\) of the three values \(P_i(s,t)\) divides a fixed nonzero integer \(R_P\) of size \(X^{O(1)}\), and likewise for \(Q\). For a positive-degree triple use the homogeneous resultant of one nonzero coordinate form with a bounded integer linear combination of all three that avoids its zeros. Such a combination exists by testing coefficients \((1,j,j^2)\) for more than twice the degree many integers \(j\). The resultant identities in the two parameter charts and \(\gcd(s,t)=1\) give the gcd divisibility. For degree zero use the gcd of the constants. A genuine point has primitive coordinate vectors, so the shift 1 gives \(|G(s,t)|=g_P^3g_Q^2\). The primitive squarefree part of \(G\) is an integral binary form of degree at least three, divides \(G\), and takes a nonzero signed divisor value of \(R_P^3R_Q^2\). [Stewart's Theorem 1, consequence on p. 796](https://uwaterloo.ca/pure-mathematics/sites/default/files/uploads/documents/s0894-0347-1991-1119199-x_0.pdf) applies to content one and nonzero discriminant, including reducible forms. Summing its uniform bound over the divisor values gives \(O_\varepsilon(X^\varepsilon)\) genuine parametrized centers at bounded degree and polynomial coefficient height.

## Polynomial-height normalizations from candidate points

Let a geometrically integral genus-zero curve over \(\mathbb Q\) in \(V\) have bounded \(\mathcal O(3,2)\)-degree \(D\). Its singular points contribute \(O_D(1)\) centers. Indeed a finite birational projection to a plane curve of degree \(D\) is an isomorphism over the smooth locus of the plane curve. The plane singular locus has at most \(D(D-1)\) points by intersecting with a nonzero partial derivative, and each projection fiber has at most \(D\) points by the linear-section degree bound. This uses the finite birational projection of [Castryck–Cluckers–Dittmann–Nguyen, Proposition 4.3.2](sources/castryck-cluckers-dittmann-nguyen-dimension-growth-2020.pdf).

If there are at most \(D^2\) smooth candidates, the count is already bounded in terms of \(D\). Otherwise select \(D^2+1\) distinct smooth candidate points. The normalization has a rational point and hence is \(\mathbb P^1_{\mathbb Q}\). Fix the first three preimages at \([1:0],[0:1],[1:1]\) by a rational change of parameter. In a polynomial system for the coefficients of the two binary-form triples, prescribe their degrees \(d_b,d_a\) from this normalization, with \(3d_b+2d_a=D\). In each triple fix to 1 a coefficient that is nonzero in the actual map, and require a nonzero resultant of two fixed bounded integer combinations by adjoining its reciprocal as a variable. Constant triples just need the nonzero fixed coefficient. Prescribe the selected projective images by cross multiplication at their preimages. For all marks after the three fixed ones use unknown affine parameters in the chart omitting \([1:0]\), which none of them equals. This gives a system with bounded dimension, number of equations and degrees in terms of \(D\), and log coefficient height \(O_D(\log X)\).

Every complex solution gives a morphism with image degree at most \(D\) containing the \(D^2+1\) distinct selected points. By the projective Bézout degree bound the image is the original curve. Its pullback of the embedding line bundle has degree \(D\), so the map has degree one onto that curve. The three specified preimages on its smooth locus fix the parameter change on the normalization. The coefficient scales, other preimages, and resultant reciprocals are then unique. Thus the affine system has a single point, which is rational because the original normalized map and preimages are rational. [Krick–Pardo–Sombra, Corollary 2.11](sources/krick-pardo-sombra-nullstellensatz-math9911094.pdf) bounds the log height of this point by \(O_D(\log X)\) from the height of the equations. Clearing denominators gives coefficient heights \(X^{O_D(1)}\). The preceding rational-curve count applies to all smooth rational points via this normalization, giving a uniform \(O_{D,\varepsilon}(X^\varepsilon)\) count on the full genus-zero curve including its bounded singular locus.

The residual arithmetic count concerns shared surfaces outside the proved component classes and any remaining curves of degree at most 21. A bounded-degree geometrically integral curve in this residual can only require a new curve count in positive genus. The sufficient collision bound retains both shifts, squarefreeness and primitivity on the surviving points.
