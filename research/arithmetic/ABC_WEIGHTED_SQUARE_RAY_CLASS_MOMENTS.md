# Ray class moments for two weighted square letters

A full divisor of one square root gives a ray condition on a norm-one ideal built from the third letter. In the real quadratic case a short interval for the two logarithmic embeddings makes the ray condition effective even when the fundamental unit is large. This gives a uniform count on the balanced endpoint sector where two selected letters have small combined odd-valuation radical.

## A uniform count with a divisor on one square root

Fix \(f,\mu,\ell\ge0\), \(\nu>0\), and \(\eta>0\). Let \(\lambda\) be sufficiently small and fixed before \(X\to\infty\). Choose two roles in a positive primitive triple \(a+b=c\le X\), and write their complete values as
\[
 U=A x^2,\qquad V=B y^2,
 \qquad D=AB\le X^{f\lambda},
\]
where \(A,B\) are positive squarefree. Denote the third complete letter by \(T\) and suppose \(\operatorname{rad}(T)\le X^{\ell\lambda}\). If the selected letters are opposite, take \(V>U\) and assume \(D>1\). Either relation \(T=U+V\) or \(T=V-U\) is permitted. Require at least one integer \(M\) satisfying
\[
 M\mid x\ \text{or}\ M\mid y,\qquad
 \gcd(M,2D)=1,\qquad M\ge X^\nu,\qquad
 \operatorname{rad}(M)\le X^{\mu\lambda}.
\]
Uniformly in these data the number of triples is
\[
 \ll_{f,\mu,\ell,\nu,\eta,\lambda,\varepsilon}
 X^{(\mu+\ell/2+5f/4+\eta)\lambda+\varepsilon}.
 \tag{WR1}
\]
The coefficients can vary with the triples. No degree-grid condition is imposed on \(T\). The same bound permits the cap on \(T\) to count only primes outside a prescribed integer \(P\le X^B\), for fixed \(B\). The prescribed support can depend on the fixed coefficients and divisor record.

For the algebra, let \(K=\mathbb Q(\sqrt D)\) in the opposite case and \(K=\mathbb Q(\sqrt{-D})\) in the same-side case. Write \(s=\sqrt D\) or \(\sqrt{-D}\), respectively, and \(\alpha=B y+x s\). In both cases \(N_{K/\mathbb Q}\alpha=B T>0\). Also \(x,y\) are coprime, and \(T\) is coprime to \(Dxy\). If a prime ideal above a rational prime dividing \(T\) and outside \(2D\) divides both \(\alpha\) and \(\bar\alpha\), it divides \(2B y\) and \(2s x\), contradicting \(\gcd(x,y)=1\). Thus such a rational prime is split and its entire norm valuation lies on one prime ideal. There are no inert contributions outside \(2D\). Since \(M\) is prime to \(2D\), both \(\alpha\) and \(\bar\alpha\) are units at its primes. Their quotient is \(1\) modulo \(M\mathcal O_K\) when \(M\mid x\) and \(-1\) when \(M\mid y\).

Fix \(D,A,B\) and a choice of \(M\) and its role. The coefficient cost is \(X^{f\lambda+o_X(1)}\), since \(D\) determines at most \(2^{\omega(D)}\) allocations of its primes. The divisor cost is \(X^{\mu\lambda+o_X(1)}\) by enumerating its radical and then the integers at most \(X\) supported there. For each split prime outside \(2D\), choose one of the two primes of \(K\) above it as the allowed orientation. Independent fair choices retain an \(X^{-o_X(1)}\) proportion of the rows, since each row specifies at most \(\omega(T)=O(\log X/\log\log X)\) choices. Fix an orientation assignment with at least this many rows. This ensures that every product of allowed prime ideals is disjoint from its conjugate outside \(2D\).

Choose a fixed cutoff \(J\) at the fixed budget, with \(1/(J+1)<\eta\lambda/4\). The full part of \(T\) on valuations larger than \(J\) has radical at most \(X^{1/(J+1)}\), so it can be fixed at cost \(X^{1/(J+1)+o_X(1)}\). Contributions of primes over \(2D\) to \((\alpha/\bar\alpha)\) can also be fixed at subpower cost. At odd primes dividing \(D\) conjugation fixes the ramified prime ideal, so its exponent cancels. At \(2\) there are at most two exponents of size \(O(\log X)\). For the version with a prescribed \(P\), first record the full part of \(T\) supported on \(P\) at subpower cost. The choices of its prime ideals above each rational prime also have subpower cost, since there are at most two primes of \(K\) over it and \(\omega(P)=o(\log X)\). Put these fixed ideal factors into the same record. In both versions the records fix a single norm-one fractional ideal \(\mathfrak d\), prime to \(M\), such that
\[
   (\alpha/\bar\alpha)=
     \mathfrak d\prod_{1\le e\le J}
              (\mathfrak a_e/\bar{\mathfrak a}_e)^e,            \tag{WR2}
\]
where \(\mathfrak a_e\) is the squarefree ideal above the remaining primes outside \(2D\) of exact valuation \(e\) in \(T\). The ideals \(\mathfrak a_e\) are integral and prime to \(M\).

## The compact group and its finite torsion

Write \(I^1(M)\) for the multiplicative group of fractional ideals of absolute norm one and prime to \(M\mathcal O_K\). For \(\gamma\in K^\times\), \(\gamma\equiv1\bmod M\) means that \(\gamma-1\in M\mathcal O_{K,\mathfrak p}\) at every prime \(\mathfrak p\mid M\). No real sign condition is imposed. In the real field write \(r(\gamma)=\log|\sigma_1(\gamma)/\sigma_2(\gamma)|\) for a fixed order of the two embeddings. Define the additive compact group
\[
 G=\bigl(I^1(M)\times\mathbb R\bigr)/
       \{((\gamma),-r(\gamma)):\ \gamma\equiv1\bmod M,\
                                      |N\gamma|=1\}.
                                                               \tag{WR3}
\]
In the imaginary field omit the real coordinate and quotient \(I^1(M)\) by the indicated principal ideals. Write \(g(\mathfrak i)\) for the image of \((\mathfrak i,0)\), or of \(\mathfrak i\) in the imaginary case. Put \(h(t)=[(1),t]\) in the real case. In the imaginary case use \(h(t)=0\).

The subgroup in (WR3) is closed. In each fixed ideal component its real coordinates, when nonempty, form a translate of the logarithms of the units congruent to one. These logarithms form a discrete lattice of rank one in the real case, since that congruence cuts out a finite-index subgroup of the units. The image \(H=h(\mathbb R)\) is therefore a circle in the real case and trivial in the imaginary case. The quotient \(G/H\) is a subgroup of the ordinary ray class group \(\mathrm{Cl}_M(K)\) with finite modulus \(M\mathcal O_K\), so \(G\) is compact. Its \(e\)-torsion has size at most \(e\,|\mathrm{Cl}_M(K)[e]|\). [Degree-zero Arakelov class groups, Schoof](https://library.slmath.org/books/Book44/files/14schoof.pdf)

The ray class exact sequence has kernel on projection to \(\mathrm{Cl}(K)\) given by \((\mathcal O_K/M\mathcal O_K)^\times\) modulo the image of the units. In a finite abelian group a quotient has at most as many elements killed by \(e\) as the original group, since the order of the \(e\)-torsion equals the order of the quotient by \(e\)-multiples. It follows that
\[
   |G[e]|\le e\,|\mathrm{Cl}(K)[e]|\,
                         |(\mathcal O_K/M\mathcal O_K)^\times[e]|
          \le |\mathrm{Cl}(K)[e]|X^{o_X(1)} .                 \tag{WR4}
\]
For the last inequality, \(e\le J\) is fixed in the height limit. The Chinese remainder theorem reduces the local factor to primes \(\ell\mid M\), all odd and unramified in \(K\). In a local factor of residue degree \(d_{\mathfrak p}\) at modulus \(\ell^k\) there are at most \(e\) residue roots. For a principal unit \(u\ne1\), odd unramified lifting of the exponent gives \(v_{\mathfrak p}(u^e-1)=v_{\mathfrak p}(u-1)+v_\ell(e)\). This follows by applying the binomial expansion first to the part of \(e\) prime to \(\ell\) and then successively to multiplication of the exponent by \(\ell\). Hence the principal-unit kernel has at most \(\ell^{d_{\mathfrak p}v_\ell(e)}\) classes, independently of \(k\). The total local contribution at \(\ell\) is at most \(e^4\), since the sum of the local degrees is two. Multiplication over \(\ell\mid M\) gives \(X^{o_X(1)}\) uniformly for \(M\le X\), including nonsquarefree moduli. [Ray class exact sequence, Stevenhagen](https://websites.math.leidenuniv.nl/algebra/Stevenhagen-CFT.pdf)

## The short logarithmic condition and its moments

In (WR2) the generator \(\gamma=\alpha/\bar\alpha\) has norm one. If it is \(-1\) modulo \(M\), multiply it by \(-1\) to obtain a generator congruent to one. Its ideal and \(r\) are unchanged. Hence the left side of (WR2) maps to \(h(r(\gamma))\) in the real field and to \(0\) in the imaginary field.

In the real field, each embedding of \(\alpha\) has absolute value at most \(2\sqrt{BX}\), and the nonzero integral norm is at least one. Thus \(|r(\gamma)|=O(\log X)\) uniformly. Put \(W=\nu\log X/10\) and partition this range into \(O_\nu(1)\) intervals of width at most \(W\). Fix one interval. Let \(U_0=h([-W,W])\subset H\) and let \(K_0\) be the autocorrelation of its indicator for probability Haar measure on \(G\), divided by \(\operatorname{meas}(U_0)\). This function is supported in \(h([-2W,2W])\), is at most one and is at least \(1/2\) on \(h([-W/2,W/2])\). If the interval wraps around the circle, the claimed lower bound still follows by projecting the overlap of intervals on the line. If \(U_0=H\) the value is one on \(H\). Its Fourier coefficients are nonnegative and sum to one by Parseval. In the imaginary field take \(K_0=\mathbf 1_{\{0\}}\), with the same coefficient property. Translating \(K_0\) to test (WR2) changes Fourier coefficients by factors of modulus one.

Bin the ideal norms with \(N\mathfrak a_e\le X^{\lambda\theta_e}\), choosing a finite cap mesh at the fixed budget with \(\sum_e\theta_e\le\ell+\eta/4\). For each slot let \(\mathcal A_e\) be its dictionary of allowed ideals and let
\[
       F_e(\chi)=\sum_{\mathfrak a\in\mathcal A_e}
                         \chi(e\,g(\mathfrak a/\bar{\mathfrak a})).
\]
The global orientation choice makes \(N\mathfrak a\) determine \(\mathfrak a\), so \(\#\mathcal A_e\le X^{\lambda\theta_e}\). Slots of zero cap can be omitted. Set \(H_0=\sum_e\theta_e\) and \(m_e=\lceil H_0/\theta_e\rceil\). These give \(m_e\theta_e\le2H_0\) and \(\sum_e 1/m_e\le1\). If every slot is omitted, there is just one tuple of ideals.

For a moment of order \(2m_e\), expand \(|F_e|^{2m_e}\) against the probability weights \(\widehat K_0(\chi)\). Fix the first product \(\mathfrak A\) of \(m_e\) ideals. Any second product \(\mathfrak B\) that contributes satisfies
\[
  e\,g\bigl((\mathfrak A/\bar{\mathfrak A})
                        (\bar{\mathfrak B}/\mathfrak B)\bigr)
          \in h([-2W,2W]).
\]
The preimage of this arc is contained in at most \(|G[e]|\) translates of \(h([-2W/e,2W/e])\). Indeed if \(e\,g=h(t)\), then \(g-h(t/e)\) is killed by \(e\). In the imaginary case read these sets as points. Two second products \(\mathfrak B_1,\mathfrak B_2\) in one translate differ by an ideal quotient
\[
   (\mathfrak B_1/\bar{\mathfrak B}_1)
                   (\bar{\mathfrak B}_2/\mathfrak B_2)=(\zeta)
\]
where \(\zeta\equiv1\bmod M\), \(|N\zeta|=1\), and \(|r(\zeta)|\le4W/e\) in the real case. This follows directly from the kernel in (WR3). In the imaginary case \(|\sigma(\zeta)|=1\). Put \(D_1=N(\mathfrak B_1\mathfrak B_2)\). Its rational integer ideal clears the denominator of \((\zeta)\), so \(D_1(\zeta-1)\) is integral and lies in \(M\mathcal O_K\). Its absolute values at both embeddings are \(O(X^{4H_0\lambda+\nu/5})\). For sufficiently small fixed \(\lambda\) this is \(o(X^\nu)\), uniformly since \(H_0\le\ell+\eta/4\). A nonzero element of \(M\mathcal O_K\) has absolute norm at least \(M^2\), so \(\zeta=1\). The orientation condition now gives \(\mathfrak B_1=\mathfrak B_2\). A fixed product of a fixed number of ideals of polynomial norm has subpower many ordered factorizations by the divisor bound. Therefore
\[
 \left(\sum_\chi \widehat K_0(\chi)|F_e(\chi)|^{2m_e}\right)^{1/(2m_e)}
    \le X^{\lambda\theta_e/2+o_X(1)}
                         |\mathrm{Cl}(K)[e]|^{1/(2m_e)} .
                                                               \tag{WR5}
\]

Fourier inversion, absolute values and Hölder bound the number of tuples obeying (WR2) by a constant times the product of the norms in (WR5). With fixed ideal data, the possible generators \(\gamma=\alpha/\bar\alpha\) differ by units. In the real field their logarithms lie in a range of length \(O(\log X)\). Distinct nonzero unit logarithms in bounded degree have uniform positive separation, since a unit with both absolute values within a fixed bounded box belongs to a finite list by the integer-coefficient bounds on its minimal polynomial. There are therefore \(O(\log X)\) choices, with at most two signs for a fixed log. The whole unit group is bounded in the imaginary case. Each \(\gamma=(By+xs)/(By-xs)\) determines \(x/y\), and \(\gcd(x,y)=1\) recovers the pair.

The factor from nontrivial class groups is explicitly bounded by
\[
                   \prod_{e=1}^J |\mathrm{Cl}(K)[e]|^{1/(2m_e)}.
\]
The uniform class-number bound \( |\mathrm{Cl}(K)|\ll_\varepsilon D^{1/2+\varepsilon}\) follows by counting integral ideal representatives of bounded norm given by Minkowski's theorem. Because \(\sum_e1/(2m_e)\le1/2\), the factor costs at most \(X^{(f/4)\lambda+o_X(1)}\) in the exponent sense. Add the coefficient, modulus and tail costs to prove (WR1). All choices of \(J\), the cap mesh and the number of slots precede the height limit at fixed \(\lambda\).

For a refinement, \(\mathrm{Cl}(K)[2^{v_2(e)}]\) has subpower size by quadratic genus theory and the fixed cutoff \(J\). With
\[
        \Theta=\max_{p\ {\rm odd\ prime}}
                       \sum_{\substack{e:p\mid e}}\frac1{m_e}
\]
and maximum zero when all retained orders are powers of two, the torsion cost in this bin is at most \(\Theta f/4\). Indeed each odd-primary factor of the class number appears in the weighted product with total exponent at most \(\Theta/2\). This uses exact valuation indices.

## Exclusion at the balanced endpoint

In the original budget \(\operatorname{rad}(abc)<c^\lambda\), suppose that on some two letters the product of the primes of odd exact valuation is at most \(c^{f\lambda}\), where \(f<2/15\) is fixed. There is a constant \(\delta(f)>0\) such that this sector is
\[
            O_{f,\lambda,\varepsilon}
                   \left(X^{(1/3-\delta(f))\lambda+\varepsilon}\right)
\]
for all sufficiently small fixed \(\lambda\). The selected two letters can vary from row to row by a finite union over the pairs.

Otherwise take a saturating sequence with inner height limits at fixed budgets and with \(\lambda\to0\). The one-letter quarter inequality gives balanced radical caps \(1/3+o(1)\) on each letter. The pair bound for the nonatomic degree masses gives total nonatomic mass at least \(1/3\) on the selected two letters. On those two letters the roots \(x,y\) retain the even-valuation powers away from their odd-kernel supports. The removed odd mass is at most \(f\), and removing the prime \(2\) costs vanishing normalized mass. After splitting the degree measures by parity and taking a subsequence, the retained even powers therefore have positive nonatomic mass. At one of the two letters select a compact positive-degree band of arbitrarily small positive normalized mass. Its full product in the square root has positive numerical weight, bounded below along the subsequence, and supplies \(M\) with an arbitrarily small fixed \(\mu>0\) and a fixed \(\nu>0\). The square root halves the numerical weight and preserves the radical of this even band. The choice may be made after the limiting profile is fixed, then held fixed in the inner height limits.

Apply (WR1) with the third cap approaching \(1/3\). Its resulting coefficient can be made as close as desired to \(1/6+5f/4<1/3\), a contradiction. When the selected opposite kernels give \(D=1\), the original values are exact opposite squares and the cyclotomic argument supplies a strict gap. This exclusion uses no hypothesis on the numerical height balance of the chosen pair. Its hypothesis concerns the actual primitive odd-valuation supports and is preserved when counting the primitive quotients in the minimal-repair height band. [Pair nonatomic bound and one-letter quarter inequality](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees) · [Opposite exact squares](ABC_MINIMAL_REPAIR_CYCLOTOMIC_SPLIT.md#additive-identities-for-opposite-even-powers)

In particular any balanced saturating sequence has combined odd-valuation radical mass at least \(2/15\) on each pair, and hence total odd-valuation mass at least \(1/5\). In the minimally repaired sector with total primitive budget at most \((1+2\xi)\lambda\) and lower target cap \(c^{(1/3-\xi)\lambda}\), apply the preceding exclusion at budget \((1+2\xi+\zeta)\lambda\) for any fixed \(\zeta>0\). The exponent for repair height \(Y\) is at most
\[
 \frac{(1/3-\delta(f))(1+2\xi+\zeta)}{M(1/3-\xi)}
       +\varepsilon.
\]
Here \(f\) still bounds the primitive kernel cost relative to the enlarged budget, and \(c\le Y^{1/(M(1/3-\xi)\lambda)}\). For \(\xi,\zeta\) sufficiently small in terms of \(f\), this is strictly below \(1/M\). It counts only the forced repair of each primitive quotient. [Minimal-target-repair equivalence](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)
