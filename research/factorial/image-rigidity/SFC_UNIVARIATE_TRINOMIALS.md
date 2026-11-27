# Three consecutive factorial moments of univariate trinomials

Let \(\mathcal L(x^m)=m!\), extended linearly over a field of characteristic zero. The target is to exclude simultaneous zeros of \(\mathcal L(f^n),\mathcal L(f^{n+1}),\mathcal L(f^{n+2})\) for every \(n\ge1\) and \(f=Ax^a+Bx^b+Cx^c\) with \(0\le a<b<c\) and \(ABC\ne0\). This is the three-term univariate case of [Edo–van den Essen, Conjecture 2.8](https://arxiv.org/pdf/1304.3956v2). Their two-term results in §3.2 and their special families do not give this support-three statement. The [univariate quadratic proof](DENSE_UNIVARIATE_QUADRATIC_STRONG_FACTORIAL.md) covers \((0,1,2)\). Neither factoring out \(x^a\) nor substituting a power of \(x\) preserves the order-\(j\) moments in general: the weights become \((aj+m)!\) or \((dm)!\), respectively.

## A fixed window and a distant third exponent

**Proposition.** Fix \(a<b\) and \(n\ge1\). Suppose the two binary forms \(P_j(A,B)=\mathcal L((Ax^a+Bx^b)^j)\), for \(j=n,n+1\), have no common projective zero over \(\mathbb C\). Then there is \(c_0(a,b,n)\) such that for every integer \(c\ge c_0\) the three moments of the trinomial cannot all vanish. The same cutoff works over every characteristic-zero field and for all coefficient triples other than the zero triple.

Write \(M_j(A,B,C)=\mathcal L(f^j)\) and \(t_j=((jc)!)^{-1/j}\). In the expansion of \(M_j(A,B,t_j z)\), the coefficient of \(z^j\) is 1 and the constant coefficient is \(P_j\). For \(0<k<j\) the coefficient of \(z^k\), uniformly on \(\max(|A|,|B|)\le1\), has absolute value at most a constant depending on \(j\) times
\[
\frac{(kc+(j-k)b)!}{((jc)!)^{k/j}}=\exp\!\big(kc\log(k/j)+O_{j,b}(\log c)\big).
\]
The bound follows by choosing the \(k\) factors from the highest term and then expanding the other \(j-k\) factors. The equality follows by Stirling's formula at the two positive arguments, with \(j,k,b\) fixed. Thus the intermediate coefficients tend to zero. Stirling also gives \(\log t_j=-c\log c-c\log j+c+O_j(\log c)\). In particular \(t_{n+2}/t_j\to0\) for \(j=n,n+1\).

Suppose simultaneous zeros exist for unbounded integers \(c\). A zero with \(A=B=0\) would have \(C=0\), by the highest-term moment. Otherwise homogeneity lets us normalize \(\max(|A|,|B|)=1\). The monic polynomial \(M_{n+2}(A,B,t_{n+2}z)\) has uniformly bounded coefficients, so the Cauchy root bound gives \(|C|=O(t_{n+2})\). For \(j=n,n+1\), we have \(C/t_j\to0\), and the displayed coefficient bounds give \(M_j(A,B,C)-P_j(A,B)\to0\). Passing to a convergent subsequence of the normalized pairs gives a common projective zero of the two binary forms, a contradiction.

For each fixed \(c\) the equations and projective charts are defined over \(\mathbb Q\). A solution in a characteristic-zero field gives a proper ideal in the corresponding affine chart over \(\mathbb Q\), hence a point over \(\overline{\mathbb Q}\) by the weak Nullstellensatz and then a complex solution. This proves the field-uniform statement.

The hypothesis at \(n=1\) also has a direct proof. Indeed \(P_1=0\) forces the coefficients to be complex-proportional to a nonzero real pair. After scaling to that pair, \(P_2=\int_0^\infty (Ax^a+Bx^b)^2 e^{-x}\,dx>0\). For arbitrary \(a,b,n\), the [retained all-order binomial theorem](ORDINARY_BINOMIAL_RESIDUAL_REGION.md) supplies the hypothesis. Consequently **for every fixed \(a,b,n\)** the trinomial window is nonzero for all sufficiently large \(c\), uniformly in the coefficients and characteristic-zero field. A bound uniform in the starting order, and the windows left below the resulting cutoffs, remain open.

## Computing a cutoff

The separation cutoff can be computed directly from the two lower moment polynomials. In each affine chart \(A=1\) or \(B=1\), use the Euclidean algorithm over \(\mathbb Q\) to find a Bézout identity \(UP_n+VP_{n+1}=1\). Choose a rational upper bound \(L\ge1\) for the sum of the absolute coefficients of \(U,V\) in either chart. On the respective unit disks at least one lower moment has modulus \(\ge1/L\).

Write \(R_{jk}(c)=(kc+(j-k)b)!/((jc)!)^{k/j}\) for \(0<k<j\). Its successive quotient is at most \(\big((kc+(j-k)b+k)/(jc+1)\big)^k\), which is \(<1\) for all \(c>b+k/(j-k)\). Thus we can find an integer \(c_1\) in this range where all \(R_{jk}(c)\le1\), and the inequalities persist for \(c\ge c_1\). In the monic polynomial at order \(n+2\), the \(k\)-th coefficient is bounded by \(\binom{n+2}{k}2^{n+2-k}R_{n+2,k}(c)\) for intermediate \(k\). Its constant coefficient is bounded by \(2^{n+2}((n+2)b)!\). These bounds and the Cauchy bound give a rational bound \(K\) on all its root moduli for \(c\ge c_1\).

For \(j=n,n+1\) the successive quotient of \(t_{n+2}/t_j\) is at most \(j(c+1)/((n+2)c+1)\), also \(<1\) for all sufficiently large \(c\) with a computable starting point. At a root of \(M_{n+2}\) in the normalized charts we have \(|C/t_j|\le Kt_{n+2}/t_j\). Increase \(c_1\) until this bound is \(<1\) and, for both \(j\), smaller than \(1/(L3^j)\). Then the contribution of all positive powers of \(C/t_j\) in \(M_j\) has modulus \(<1/L\), since the sum of the bounds on their coefficients is at most \(3^j\). The same inequalities hold for every larger \(c\). All comparisons in this search use positive rational powers of factorials and can be made exactly by clearing powers. The limiting estimates guarantee termination.
