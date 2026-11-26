# Valuation spans and multiplicative cosets in small-budget abc

A set of prime valuations that determines a complete letter or a primitive pair ratio gives a strict saving below the one-third coefficient if its prime partition sum is sublogarithmic at any fixed multiplier \(a<1\). In particular this applies to arbitrary rational multiplicative cosets of rank \(O((\log X)^\kappa)\) for every fixed \(\kappa\). The same valuation projection bounds fractional packings of separating prime sets in balanced saturating families, without a rank restriction.

## Counting low-rank cosets by separating valuations

For a positive reduced rational \(u=a/b\), write
\(\mathcal R(u)=\operatorname{rad}(ab)\), with \(\mathcal R(1)=1\).
For a fixed-length tuple \(\boldsymbol u=(u_1,\ldots,u_d)\), use \(\mathcal R(\boldsymbol u)\) for the radical of the product of the reduced numerators and denominators. Let \(H_0>0\), \(\rho\ge0\) and \(\kappa>0\) be fixed. Set
\[
                       \beta_\kappa=(1-1/\kappa)_+ .
\]
For every \(\varepsilon>0\), every multiplicative subgroup \(\Gamma\le(\mathbb Q_{>0}^{\times})^d\) of rank at most \((\log X)^\kappa\), and every coset \(g\Gamma\), one has
\[
 \#\{\boldsymbol u\in g\Gamma:\ \text{all reduced numerators and
 denominators are at most }X^{H_0},\
                   \mathcal R(\boldsymbol u)\le X^\rho\}
       \ll_{d,H_0,\rho,\kappa,\varepsilon}
               X^{\beta_\kappa\rho+\varepsilon}.
 \tag{VS1}
\]
The rank here is the dimension of the rational span of the valuation vectors. There is no bound on \(g\) or on generators of \(\Gamma\).

In fact the valuation vectors have finite support in the coordinates \((j,p)\) for coordinate index \(j\) and rational prime \(p\). Their rank \(r\) means that \(r\) coordinate functionals are linearly independent on their span and hence separate its vectors. Use the same functionals on the coset, where differences of valuation vectors belong to the span. The radical of the nonzero recorded valuations is supported on at most \(r\) distinct primes. For any such squarefree radical \(Q\le X^\rho\), the at most \(d\,\omega(Q)\) recorded nonzero positions, signs and absolute valuations have \(X^{o(1)}\) possibilities, uniformly in \(Q\). There are at most \(2^{O_d(\omega(Q))}\) position and sign choices and the sum of the absolute valuations is \(O_{d,H_0}(\log X)\). The estimate follows from stars and bars and \(\omega(Q)=O_\rho(\log X/\log\log X)\).

For any fixed \(0<\sigma<1\), Rankin's inequality bounds the number of possible radicals by
\[
       X^{\rho\sigma}\prod_{p\in S}(1+p^{-\sigma})
           \le X^{\rho\sigma}
                   \exp\big(O_\sigma(r^{1-\sigma}+1)\big),
\]
where \(S\) is the set of the selected primes and \(p_j\ge j+1\) for the \(j\)-th prime. Choose \(\sigma>\beta_\kappa\) sufficiently close to \(\beta_\kappa\) that \(\rho(\sigma-\beta_\kappa)\) is smaller than the prescribed exponent allowance. Then \(\kappa(1-\sigma)<1\). This proves (VS1), including \(\rho=0\), where only one radical is allowed. The same proof is uniform for \(\rho\) in any fixed compact interval.

For \(\kappa>1\) and \(0<\rho\le\min(1,H_0)\), the exponent \(\beta_\kappa\rho\) in this general count is attained up to \(o(1)\). Generate a subgroup by all primes up to \(y=(\log X)^\kappa\). Products of \(k=\lfloor\rho\log X/\log y\rfloor\) distinct such primes are at most \(X^\rho\). The prime number theorem and \(\binom Mk\ge(M/k)^k\) give \(X^{\rho(1-1/\kappa)+o(1)}\) such products. This tests the projection count for integers. It does not give any such family of primitive additive rows.

## Strict saving from a separating prime pool

Fix \(0\le a<1\) and a function \(g(X)=o(\log X)\). Count positive primitive \(a_1+a_2=a_3\le X\) with \(\operatorname{rad}(a_1a_2a_3)<a_3^\lambda\). Let \(Z\) denote either a prescribed complete letter or a prescribed ratio of two different letters. Suppose all its possible values in the counted family are separated by a prescribed set \(S_X\) of prime valuations, and
\[
                       \sum_{p\in S_X}p^{-a}\le g(X).
 \tag{VS2}
\]
Separation means that distinct possible values of \(Z\) have distinct vectors of valuations on \(S_X\). For every fixed \(\zeta>0\), all sufficiently small fixed \(\lambda>0\), and every \(\varepsilon>0\), the count is respectively
\[
 \begin{split}
 &O_{a,g,\zeta,\lambda,\varepsilon}
       \left(X^{(\gamma_{\rm letter}(a)+\zeta)\lambda+\varepsilon}\right),
       &&\gamma_{\rm letter}(a)=\frac{1+a}{2(2+a)},\\
 &O_{a,g,\zeta,\lambda,\varepsilon}
       \left(X^{(\gamma_{\rm ratio}(a)+\zeta)\lambda+\varepsilon}\right),
       &&\gamma_{\rm ratio}(a)=\frac{a}{1+2a}.
 \end{split}
 \tag{VS3}
\]
These coefficients have deficits from \(1/3\) equal to
\((1-a)/(6(2+a))\) and \((1-a)/(3(1+2a))\), respectively. The constants are uniform in the separating set and in the family. Separating only inside each of at most \(X^{\xi\lambda}\) prescribed record fibres, with the same partition-sum bound, costs an additional \(\xi\lambda\) in the exponent.

Work first on a height shell \(X/2<a_3\le X\). Write \(\mathcal R_S(Z)=\prod_{p\in S,\ v_p(Z)\ne0}p\). Call a row heavy if
\(\log\mathcal R_{S_X}(Z)\ge s\lambda\log X\), where \(s>0\) is fixed. Its full triple radical contains the selected primes dividing \(Z\). In the pair-ratio case this follows from primitivity, so there is no cancellation between those two letters. By placing price \(a\log_Xp\) on \(S_X\) and natural price \(\log_Xp\) elsewhere, (EB5) bounds the heavy count by
\[
                 X^{(1/3-(1-a)s/3+\eta)\lambda+\varepsilon}
\]
up to a fixed multiplier, for any fixed \(\eta>0\) at sufficiently small fixed budgets. This step uses the full three-letter radical. [Discounted prime pool count (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

On the light branch, (VS2) and Rankin at price \(a\) give at most \(X^{as\lambda+o_X(1)}\) distinct projected radicals. For \(a=0\) the squarefree support count is simply \(2^{|S_X|}=X^{o(1)}\). Fixing the signs and powers on a projected support costs subpower, as in (VS1). Since the valuations separate \(Z\), the light dictionary has exponent at most \(as\). For a ratio this counts all rows, since the reduced ratio recovers the corresponding positive coprime letters and hence the third. Equating \(as\) with the heavy coefficient gives \(s=1/(1+2a)\) and \(\gamma_{\rm ratio}\).

For a light complete-letter dictionary, bin the actual radical caps of the three letters in units of \(\lambda\log X\) with fixed mesh \(\tau>0\). Write \(t,u,v\) for the caps, with \(t\) the designated letter. They sum to at most \(1+3\tau\). Let \(d=\log|\mathcal F|/(\lambda\log X)\) for the retained light family in this bin. The fibre bound (FL1) gives
\[
                    d\le as+(u+v)/4+\eta+o_X(1),
\]
since \(\min(u,v)\le(u+v)/2\). The two-letter bound (GD3) bounds \(2d\) separately by \(t+u\) and \(t+v\), with arbitrary fixed positive allowances at sufficiently small budgets. These bounds apply to exact radical bins by anatomic truncation, taking its coefficient and truncation losses smaller than the chosen normalized allowance. Adding the two orientations gives
\(t\ge4d-1-3\tau-O(\eta)-o_X(1)\). Thus \((u+v)/4\le 1/2-d+O(\tau+\eta)+o_X(1)\), and
\[
                    d\le 1/4+as/2+O(\tau+\eta)+o_X(1).
\]
Equating the last main coefficient with the heavy coefficient gives \(s=1/[2(2+a)]\) and \(\gamma_{\rm letter}\). If the individual caps are already arbitrarily close to \(1/3\), the simpler light bound \(1/6+as\) gives the stronger balanced deficit \((1-a)/(6(1+2a))\). In both cases the mesh and the fixed allowances are taken sufficiently small in terms of \(\zeta\), then the inner height limit is taken at each fixed budget. [Uniform fixed-letter bound (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools) · [Two-letter bound (GD3)](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#group-dictionary-count) · [Exact valuations in the height-first convention](ABC_SECANT_ENDPOINT_TRANSFER.md#the-height-first-setting)

For a cumulative bound, a dyadic shell of height \(Y\ge X^\theta\) with fixed \(\theta>0\) inherits a condition of the same form, uniformly with \(\widetilde g_\theta(Y)=\sup_{Y\le x\le Y^{1/\theta}}g(x)=o(\log Y)\), after ignoring primes larger than \(Y\). Summing shells costs only a subpower factor. The remaining rows of height less than \(X^\theta\) are covered by the existing unrestricted limiting bound, taking \(\theta\) small enough for the desired exponent, including when \(\gamma_{\rm ratio}=0\). These choices give the cumulative (VS3). [Unrestricted limiting count](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)

By selecting a valuation coordinate basis as in (VS1), every coset of rank \(O((\log X)^\kappa)\) satisfies the hypothesis for some fixed \(a<1\), for any fixed \(\kappa>0\). One may use \(a>\beta_\kappa\), arbitrarily close to \(\beta_\kappa\), and absorb the difference in \(\zeta\). For \(\kappa\ge1\) the resulting deficits are \(1/[6(3\kappa-1)]\) for a letter and \(1/[3(3\kappa-2)]\) for a ratio, before arbitrary allowances. A family supported on a prescribed list of cosets pays the list size in the exponent. The claim also applies to a coset condition on the actual tuple \((a_1,a_2,a_3)\) by projecting to a ratio, since that multiplicative homomorphism has image rank at most the original rank. More generally a separating basis of size at most \(X^\mu\) whose primes are all at least \(X^\nu\) satisfies (VS2) whenever \(\mu<\nu\) and \(a\ge\mu/\nu\). All these bounds concern primitive rows and thus apply to the primitive quotients of least repairs.

## Fractional packing of separating prime sets

Let \(Z\) be a random variable supported on any finite family of positive rationals of height at most \(X^{H_0}\). Its law can be the marginal law of a uniform family of primitive rows. A separating set \(S\) distinguishes the distinct possible values of \(Z\).
Uniformly in \(S\) and the law, one has
\[
                   H(Z)\le \mathbb E\log\mathcal R_S(Z)+o_{H_0}(\log X).
 \tag{VS4}
\]
Indeed every particular value of \(\mathcal R_S\) allows only \(X^{o_{H_0}(1)}\) signed valuation records by the fixed-support count in (VS1). The weight \(\mathcal R_S(z)^{-1}\) summed over distinct \(z\) is at most \(X^{o_{H_0}(1)}\sum_{Q\le X^{2H_0}}1/Q=X^{o_{H_0}(1)}\). The nonnegativity of relative entropy with respect to the probability law proportional to these weights proves (VS4).

Give the separating sets nonnegative weights \(w_S\) with \(\sum_{S\ni p} w_S\le1\) at each prime. For any finite such packing with weight \(P=\sum_S w_S\), (VS4) gives
\[
                   P\big(H(Z)-o_{H_0}(\log X)\big)
                             \le \mathbb E\log\mathcal R(Z).
 \tag{VS5}
\]
The error before multiplication by \(P\) is uniform. A prime-coordinate basis of the valuation difference span of the family separates it. Thus (VS5) applies in particular to fractional matroid base packings, with no bound on the dimension. The argument also works conditionally on any recorded datum. If each conditional support has a fractional packing of weight at least a fixed \(P\), use one scaled to weight \(P\) in each fibre and average (VS5), to obtain the same estimate for \(H(Z\mid\text{record})\) in place of \(H(Z)\).

Now use the uniform law on balanced families of primitive rows in height shells with the normalized logarithmic size tending to \(1/3\) and all three radical caps tending to \(1/3\) in units of \(\lambda\log X\). The inner height limit at fixed \(\lambda\) precedes the small-budget limit. For a ratio of two letters \(H(Z)/(\lambda\log X)\to1/3\), since the ratio determines the row, and the expected normalized ratio radical is at most \(2/3+o(1)\). For a complete letter \(T\), the fixed-letter bound (FL1) implies \(H(T)/(\lambda\log X)\ge1/6-o(1)\), because every fibre has normalized logarithmic size at most \(1/6+o(1)\); its expected normalized radical is at most \(1/3+o(1)\). Consequently in both cases
\[
                              P\le 2+o(1).
 \tag{VS6}
\]
In particular there cannot be three disjoint separating prime sets along a balanced saturating sequence. The statement covers arbitrary different choices of separating sets at each shell.

For the valuation matroid on the finite prime set \(E\) of the family, let \(r(A)\) denote the rank of the valuation coordinate functionals from \(A\subseteq E\) on its difference span. At positive rank, the maximum fractional base packing equals
\[
              \min_{A\subseteq E:\,r(A)<r(E)}
                        \frac{|E\setminus A|}{r(E)-r(A)}.
\]
Every base uses at least \(r(E)-r(A)\) coordinates outside \(A\), which proves the upper bound. For the reverse, replace each element by \(k\) parallel copies and use integral matroid base packing. A subset of copies projecting onto \(A\) has rank \(r(A)\) and its complement has size at least \(k|E\setminus A|\). Thus one can pack at least the integer part of \(k\) times the displayed minimum. Project these bases back, giving each weight \(1/k\), and let \(k\) grow. It follows from (VS6) that a balanced saturating family has a set \(A\) with
\[
           r(A)<r(E), \qquad
           |E\setminus A|\le (2+o(1))(r(E)-r(A)).
 \tag{VS7}
\]
In particular \(|E\setminus A|<3(r(E)-r(A))\) eventually. This constraint relates rank loss to the number of supporting primes. The radical weight of this set is a separate quantity. [Edmonds and Fulkerson, Theorem 2c](https://nvlpubs.nist.gov/nistpubs/jres/69B/jresv69Bn3p147_A1b.pdf#page=5)

## Counting with several separating sets

Fix \(P\ge2\). Consider all primitive triples at total radical
budget \(\lambda\) in a dyadic height shell. Suppose there is a
fractional packing of weight at least \(P\) of separating sets
for a prescribed complete letter, or for a prescribed pair
ratio. Then for every \(\eta>0\) and all sufficiently small
fixed \(\lambda>0\), the count is \(O(X^{(\gamma+\eta)\lambda+\varepsilon})\)
for every \(\varepsilon>0\), uniformly in the separating sets.
The respective coefficients are
\[
 \gamma_{\rm letter}(P)=\max\{1/4,\,2/(P+4)\},\qquad
 \gamma_{\rm ratio}(P)=\min\{1/P,\,2/(P+4)\}.                \tag{VS8}
\]
In particular three disjoint separating sets for any single
letter or pair ratio give coefficient \(2/7\). These
claims impose no separate upper bound on the multiplicative rank
or on the size of a prime pool.

Indeed the original packing remains a packing after restricting
to a subfamily. Bin the three individual radicals with a
small fixed mesh in units of \(\lambda\log X\), and absorb
the mesh error into \(\eta\). Write \(d\) for the normalized
logarithmic size of the binned family and \(t\) for the
designated letter's normalized radical cap. The letter's
normalized entropy is at most \(t/P+o_X(1)\) by (VS5).
By (FL1), conditional on that letter the fibre entropy is at
most \((1-t)/4+O(\eta)+o_X(1)\), using the sum of the two
remaining caps. The two-letter counting bound (GD3), with
the two free letters as the selected pair, separately gives
\(d\le (1-t)/2+O(\eta)+o_X(1)\). Maximize the minimum of the
two bounds as \(0\le t\le1\). When \(2\le P\le4\) their
leading terms meet at \(t=P/(P+4)\); when \(P\ge4\) the
first one is at most \(1/4\) throughout.

For the ratio let \(s\) be the sum of its two normalized
radical caps. Its entropy equals the entropy of the entire
row, so (VS5) gives \(d\le s/P+o_X(1)\). Applying (GD3)
to each pair involving the omitted letter gives
\(d\le 1/2-s/4+O(\eta)+o_X(1)\). The two leading terms
meet at \(s=2P/(P+4)\) when \(2\le P\le4\); for \(P\ge4\)
the first is at most \(1/P\). The radical boxes and the
entropy loss are taken at fixed budget before any outer limit.
The shell bounds also sum when the packing is specified at
each shell or on the whole family up to \(X\).

If each of the three complete letters has a packing of
weight at least \(P\), the coefficient improves to
\(2/(3P)\). In fact each pair of letters recovers the row,
so \(d\) is at most the sum of their two normalized
entropies. Add the three inequalities and bound the sum of
the letter entropies by \(1/P\) using (VS5) and the total
radical cap. The same bound follows if each of the three
pair ratios has such a packing: each ratio already carries
the full row entropy, and adding their three radical caps
costs at most twice the total budget. These observations allow
the individual packings to use overlapping prime sets on
different letters or ratios.

[Pointwise fixed-letter estimate (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)
· [Two-letter count (GD3)](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#group-dictionary-count)

## Literature scope and the remaining extraction problem

Agrawal, Bloom and Petridis, v2, Theorem 8, give lower bounds of exponent \(12/7-o(1)\) for the larger of \(|AA|\) and either \(|A+A|\) or \(|A-A|\) in the few-prime regime \(k=o(\log |A|/\log\log |A|)\). Proposition 1 uses small mixed product sets to extract low-rank coset covers. Their Theorems 6 and 7 turn efficient covers into additive-energy or sumset bounds under the rank condition \(r\le\delta\log M\), in addition to the stated lower bound on \(\delta\log M\). In the present height regime, (VS1) shows that rank \(O(\log X)\) cosets contain only \(X^{o(1)}\) distinct actual values, so an \(M\)-coset cover already needs \(M\ge |A|X^{-o(1)}\) on the exponent scale. [Agrawal–Bloom–Petridis, v2, Proposition 1 and Theorems 6–8](https://arxiv.org/html/2512.04931v2)

The exact product-set lower bound is automatic here, even for height-bounded rationals. Give each prime an independent numerator or denominator color and retain a rational if all its prime positions agree with their color. Each reduced rational of height at most \(X^{H_0}\) uses \(O_{H_0}(\log X/\log\log X)\) primes, so some coloring keeps \(|A|X^{-o(1)}\) members of a rational set \(A\). Products of retained rationals have no cross-cancellation. Their reduced numerator and denominator fix the corresponding uncancelled integer products, each with only \(X^{o(1)}\) ordered factorizations for any fixed number \(m\) of factors. Thus \(|A^{(m)}|\ge|A|^mX^{-o(1)}\), in particular for the integer letter and primitive pair ratio dictionaries. For dictionaries of positive normalized exponent at fixed budget, the large-product branch in Theorems 8 and 10 is already met at their exponent scales. Their separate mixed estimate for \(|A+AA|\) in Theorem 4 and Corollary 1 concerns the image of all triples in the set and does not bound this sparse equation graph. Formulas (VS3) and (VS6) give separate estimates for the actual arithmetic values.

A dictionary can evade every cheap separating set while having prescribed diffuse degree data. Fix \(n=X^\lambda\), positive slot exponents \(\alpha_j\), and exact multipliers \(e_j\), with the finite list fixed at each budget. Choose disjoint prime pools \(P_j\) in fixed-factor ranges about \(n^{\alpha_j}\), each of size \(n^{\alpha_j+o_X(1)}\), and take all \(T=\prod_j p_j^{e_j}\) with \(p_j\in P_j\). The radical has exponent \(\sum_j\alpha_j\) in base \(n\), the height has exponent \(\sum_j e_j\alpha_j\), and the dictionary has its full size exponent \(\sum_j\alpha_j\). Each \(p\in P_j\) occurs on exactly the fraction \(1/|P_j|\) of the dictionary. If a set omits two primes of one \(P_j\), switching between them while holding the other slots fixed gives two indistinguishable members. Therefore every separating set contains all but one prime of every slot and its sum in (VS2) is at least \(n^{(1-a)\alpha_j+o_X(1)}\) for each \(j\). The example is compatible with choosing increasingly many slots and diffuse limits after the height limit. Its exact multiplicative structure and marginal counts provide no additive rows.

An unrestricted saving still requires one of these configurations to occur on a sufficiently large part of each putative saturating family, or a different arithmetic mechanism. At a putative endpoint each such multiplicative rank must exceed every fixed power of the logarithm in the iterated sense of (VS3). Any separating set must also evade (VS2) with a uniform partition-sum bound, while its valuation span still satisfies (VS6). Exact product dictionary sizes alone give no contradiction to those simultaneous restrictions.
