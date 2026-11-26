# Recovering pairs from least-repaired target differences

Two binary equations recover a pair of least-repaired primitive triples. The first equation recovers the repaired targets from a recorded difference or a nonlinear factor of that difference. The second recovers the summands from a recorded power class of their ratio. Common factors at both stages are recovered by divisor bounds, even when the repaired triples have large and different common divisors.

## Pair records and their count

Fix an integer \(M\ge3\). For each of two distinct ordered positive primitive triples \(a_i+b_i=c_i\), let
\[
 \begin{split}
 R_i&=\operatorname{rad}(c_i),\qquad W_i=\prod_{p\mid c_i}p^{(v_p(c_i)-M)_+},\\
 H_i&=H_M(c_i)=R_i^M W_i,\qquad Z_i=H_i/c_i,\qquad U_i=Z_i a_i,\qquad V_i=Z_i b_i .
 \end{split}
\]
Thus \(U_i+V_i=H_i\). Put \(S=\gcd(U_1,U_2)\) and \(Q_i=U_i/S\). For an integer \(k\ge2\), put \(D=\gcd(V_1,V_2)\) and write uniquely
\[
                    V_i=D J_i y_i^k,
\]
where \(J_i\) is the \(k\)-th-power-free part of \(V_i/D\). The pairs \((Q_1,Q_2)\), \((J_1,J_2)\) and \((y_1,y_2)\) are coprime. In particular \(J_1/J_2\) is the representative of the power class obtained by taking the exponent remainders on the reduced numerator and denominator separately. Its radical measures exactly where the two repaired valuations on this summand differ modulo \(k\).

Write \(F=Q_2H_1-Q_1H_2\). Always \(F\ne0\) for distinct ordered primitive rows. If \(F=0\), then \(U_1/H_1=U_2/H_2\), hence \(a_1/c_1=a_2/c_2\), which recovers the same primitive rows. Define the record
\[
                (k,Q_1,Q_2,J_1,J_2,W_1,W_2,T),                 \tag{WD1}
\]
where \(T\) can be \(F\). More generally, it can be a nonzero value \(G(R_1,R_2)\) of an integer homogeneous form \(G\) prescribed by the other entries of the record. The form may also depend on the height bound. Its degree is bounded, its coefficients have polynomial height, and it is either a nondegenerate quadratic or has at least three distinct complex projective linear factors. In particular one can take \(T=H_1-H_2\) when this shift is nonzero, or take one nonlinear irreducible factor of
\[
                    Q_2W_1 X_1^M-Q_1W_2 X_2^M.             \tag{WD2}
\]
Evaluate the factor on the original, possibly noncoprime, pair \((R_1,R_2)\). The form (WD2) has no repeated roots and has a nonlinear factor since it has a nonreal root for \(M\ge3\). An irreducible nonlinear factor never vanishes on the positive rational pairs here. Its primitive integral coefficients are of polynomial height. One way to see the last assertion is to use that multiplication has a positive minimum sup norm on pairs of coefficient vectors of sup norm one in bounded degrees. The integer cofactor has norm at least one.

For \(Y\ge2\) restrict to \(H_i\le Y\). Each fixed record has \(O_{M,\varepsilon}(Y^\varepsilon)\) completions, uniformly in \(k\) and in its numerical entries. In the form version the constant also depends on the degree and the polynomial height bound of \(G\). If \(P\) is any fixed prime pool of product at most \(Y^B\), define \(\operatorname{rad}_{\overline P}\) to retain only primes outside it. For any \(A\ge1\), the number of pairs admitting some \(k\ge2\) with
\[
      \operatorname{rad}_{\overline P}(Q_1Q_2J_1J_2W_1W_2|T|)
                       \le A
\]
is
\[
                         O_{M,B,\varepsilon}(A Y^\varepsilon).\tag{WD3}
\]
The form version requires a prescribed rule for the forms. If this rule depends on the order, it uses \(\min(k,\lceil\log_2Y\rceil+1)\). A fixed finite list of rules costs only its size. One may interchange the two summand roles before applying the statement.

### Proof of recovery and enumeration

First recover \((R_1,R_2)\) from its equation \(G(R_1,R_2)=T\). Here \(G\) is (WD2) for the default record. If its degree is \(d\), the gcd \(g=\gcd(R_1,R_2)\) has \(g^d\mid T\). After choosing \(g\) from the divisor list of \(|T|\), the remaining primitive pair satisfies \(G(r_1,r_2)=T/g^d\). For a form with at least three distinct projective factors, Evertse's coefficient-uniform Thue–Mahler bound gives at most \(7^{d^3(2\omega(|T|)+3)}\) primitive pairs, using the primes of \(|T|\). Signs and reducible forms are allowed. Since \(|T|\le Y^{O(1)}\), this and the divisor list cost \(Y^\varepsilon\) with any prescribed allowance. For a nondegenerate quadratic \(G(U,V)=aU^2+bUV+cV^2\) with \(a\ne0\), complete the square to obtain
\[
               (2aU+bV)^2+(4ac-b^2)V^2=4aT .
\]
The coefficients, nonzero target, and variable bounds are polynomial in \(Y\). Apply the uniform binary quadratic count of Vaughan–Wooley, in the form of Browning–Verzobio's Lemma 2.3. If \(a=0\), use \(V(bU+cV)=T\) with \(b\ne0\) instead. This also bounds the quadratic completions by a subpower. [Evertse, Theorem 6](https://ir.cwi.nl/pub/2463/2463D.pdf#page=14) · [Browning–Verzobio, Lemma 2.3](https://arxiv.org/html/2608.24512v1#S2.SS1)

Knowing the target roots and the two \(W_i\) gives the actual \(H_i\) and hence \(F\), even if another form value was recorded. The repaired additive equations give
\[
      F=Q_2V_1-Q_1V_2=D(Q_2J_1y_1^k-Q_1J_2y_2^k).           \tag{WD4}
\]
In particular \(D\mid F\) and has only divisor-many possibilities. For each such choice use the resulting binary equation of nonzero level \(F/D\) to recover \((y_1,y_2)\). At fixed \(k\ge3\) this uses the same Thue–Mahler bound, since the binomial has \(k\) distinct complex projective roots and the pair is coprime. At \(k=2\) the binary quadratic bound applies with both coefficients nonzero. This is the Pell step including the split quadratic case. To make the bound uniform in growing \(k\), first fix a large integer cutoff in terms of the desired allowance. Above it, enumerate \(y_1\le Y^{1/k}\), after which the positive \(y_2\) is determined. Below it there are only finitely many degrees. Finally \(V_i\) is known and \(U_i=H_i-V_i\). The gcd of the actual repaired summands recovers \(Z_i\) because the input triple was primitive. Thus the reconstruction has no common multiplier index.

For (WD3), orders only up to \(\lceil\log_2 Y\rceil+1\) are needed. At any larger order both roots in the displayed decomposition are one, and all the residue parts agree with those for \(\lceil\log_2Y\rceil+1\). One can always truncate before applying the form rule. Equivalently, to allow rules depending on an unbounded order, require the form at such an order to be chosen by its truncated value. The default and difference-factor versions already have this property. Enumerate the joint outside radical, which has at most \(A\) choices. Each integer entry is then supported on this radical and the fixed pool. All nonzero entries have magnitudes bounded by a fixed power of \(Y\), so a uniform fixed-support count gives \(Y^\varepsilon\) possible tuples, with any prescribed allowance. For example apply the Rankin bound at a small fixed power \(t>0\), using \(\prod_{p\mid Q}(1-p^{-t})^{-1}=Y^{o(1)}\) for any support product \(Q\le Y^{O(1)}\). The signs and the order choices cost a subpower. Combining the two counts proves (WD3). [Uniform fixed-support enumeration](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

## Several nonlinear factors at one cost

Fix the coefficient entries \(Q_i,J_i,W_i\), and fix the shared target radical \(g=\gcd(R_1,R_2)\). Write \(R_i=gr_i\). Let \(G_1,\dots,G_t\) be the distinct primitive irreducible nonlinear factors of (WD2). Enlarge the fixed pool by the primes of the nonzero pairwise homogeneous resultants. It still has polynomial-height product for fixed \(M\). Outside this pool the values \(G_j(r_1,r_2)\) have pairwise disjoint prime supports. Indeed, primitiveness of the input pair makes a common vanishing modulo a prime a common projective root, which is detected by the resultant. If their joint outside radical is at most \(A\), some factor has outside radical at most \(A^{1/t}\). Fixing that factor value recovers the primitive roots \(r_i\) at subpower cost, then the known \(g\) recovers \(R_i\). Thus the bound for these coefficient and gcd records is \(O_{M,B,\varepsilon}(A^{1/t}Y^\varepsilon)\). The gcd may vary over a prescribed dictionary, or be paid for by its outside radical along with the other coefficient records. This disjointness argument uses reduced roots. A prime on the shared gcd of the original roots could divide every unreduced factor value.

For instance, if \(Q_2W_1/(Q_1W_2)\) is a rational \(d\)-th power for some \(d\mid M\), write the two coefficients of (WD2) as \(E u^d,E v^d\) with coprime positive \(u,v\). The displayed form is \(E\) times the difference of two \(d\)-th powers of \(uX_1^{M/d}\) and \(vX_2^{M/d}\). Every divisor \(e\ge3\) of \(d\) supplies a nonlinear homogeneous cyclotomic factor at those two arguments. Its roots are simple and distinct from the roots for another divisor, so these evaluated forms can also be used as the \(G_j\) in the argument. For equal coefficients one can take \(d=M\). This buys a factor in the target-record cost when a low joint radical for these nonlinear factors is available.

## Cost in the small-budget formulation and the remaining incidence

If \(c_i\le X\), then \(\operatorname{rad}(W_1W_2)\le X^{2/(M+1)}\), since every prime in a tail has original valuation at least \(M+1\). With \(c_i\le X\) and \(\operatorname{rad}(c_i)\le X^{(1/3+\xi)\lambda}\), all repaired heights are at most \(X^{1+M(1/3+\xi)\lambda}\). Take \(M=\lceil K/\lambda\rceil\) with \(K\) fixed before the height limit. For any prescribed outside pool of polynomial height, the number of distinct original rows lying on at least one pair satisfying
\[
          \operatorname{rad}_{\overline P}(Q_1Q_2J_1J_2|T|)
                             \le X^{s\lambda}
\]
is \(O(X^{(s+2/K)\lambda+\varepsilon})\). Each pair supplies at most two rows. If a prescribed dictionary of \(O(X^{\rho\lambda})\) pools is allowed, add \(\rho\lambda\) to the exponent. Thus \(s+\rho+2/K<1/3\) is a fixed small-budget saving on this whole incidence sector. If both repaired heights are at most \(Y\) and both original targets also satisfy \(R_i\ge c_i^{(1/3-\xi)\lambda}\), one can take \(X=Y^{1/[M(1/3-\xi)\lambda]}\). The resulting exponent in repaired height is at most \((s+\rho+2/K)/[M(1/3-\xi)]+\varepsilon\), strictly below \(1/M\) when \(s+\rho+2/K<1/3-\xi\). The tail charge in original height is independent of its distribution among the roots. [Least target repair](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

The needed input for a global application is a pair relation reaching enough near-saturating rows with inexpensive joint records. Private primes give a specific constraint on such an extraction. For \(T=F\), let a prime outside the prescribed pool be supported on exactly one letter of the six original letters and have valuation \(e>0\) there. Its contribution to the joint radical in (WD3) is as follows.

| Primitive letter | When the prime contributes | Reason |
| --- | --- | --- |
| \(a_i\) | Always | It occurs in \(Q_i\). |
| \(b_i\) | \(k\nmid e\) | Its repaired valuation is unchanged. It occurs in \(J_i\) exactly at this defect and does not divide \(F=Q_2V_1-Q_1V_2\). |
| \(c_i\) | \(e\ne M\) | For \(e<M\) it occurs in \(Q_i\) and \(F\), and for \(e>M\) it occurs in \(W_i\). For \(e=M\) it cannot divide either summand of its row, and the weighted difference is a unit. |

For example, for a private prime on \(c_1\) with \(e<M\), the competing terms in \(Q_2H_1-Q_1H_2\) have valuations \(M\) and \(M-e\). The prime enters \(J_1\) if \(k\nmid M-e\), with no extra cost when it is already charged to \(Q_1\). For all target datum forms, the radical from private primes of the two \(a_i\), and from private primes of the \(c_i\) of valuation below \(M\), is already charged to \(Q_1Q_2\). An arbitrary pair with nearly disjoint balanced supports therefore has cost well above the desired threshold, even before charging \(J_i\) or \(T\). In the extreme case of disjoint supports and no prescribed pool this lower cost is \((4/3-o(1)-2/K)\lambda\log X\) when both original pairs of caps are attained near \(X\). An affordable pair must have repeated support or a pool that already captures much of this contribution. Control of an affordable target difference or nonlinear factor is a separate requirement.

The power-class input is essential for a subpower completion estimate with just the left ratio and the exact targets. Fix distinct primes \(R_1,R_2\) near a growing \(N\), and take \(c_i=R_i^M\), so no repair multiplier is needed. For each positive \(S<\min(c_1,c_2)\) coprime to \(R_1R_2\), use \(S+(c_i-S)=c_i\). These give \(\gg_M N^M\) distinct actual pairs with \(Q_1=Q_2=1\), \(W_1=W_2=1\) and the same nonzero \(F=c_1-c_2\). Both targets have a single prime in their supports. This tests completion without the right power-class record. It imposes no small radical cap on the two original summands.

The fixed-target #943 argument first finds a reduced kernel triple and then recovers a common scale by a Pell or Thue equation at a known target. Here least repair makes two moving targets recoverable from one binary datum even when their common radical and the common divisor of the right summands are unknown. The binary estimates are classical. The pair count and the incidence requirement concern the primitive least repairs. [Fixed-level mixed-power lifting](ERDOS_943_MIXED_POWER_KERNEL_GCD_WEIGHTS.md#2-an-exact-fixed-level-lifting-problem) · [Pell recovery](ERDOS_943_PELL_RECOVERED_KERNEL_GCD_WEIGHT.md#2-the-actual-cubic-and-norm-equations)
