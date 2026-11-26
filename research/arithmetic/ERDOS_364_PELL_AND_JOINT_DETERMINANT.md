# Pell descent and a joint surface cover for Erdős #364

Let \(C(X)\) count integers \(2\le n\le X\) for which \(n-1,n,n+1\) are all powerful. An integer is powerful if every positive prime valuation is at least two. The Pell reduction separates the candidates whose centers are fundamental Pell traces. Combining estimates from two papers by Plymesser lowers the exponent for \(C(X)\) from \(77/320\) to \(293/1220\), conditional on the inputs stated below. A joint determinant gives a surface cover in every dyadic kernel box. A point count on the surfaces in the near-binding boxes would lower the exponent to \(6/25\).

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

## Transfers from the preprint counting statements

Plymesser's [powerful-products paper](sources/plymesser-powerful-products-hexagon2611.00602-v1.0.1.pdf), Lemma 8.3, states that for every \(\varepsilon>0\) there is \(\eta>0\) such that fixed positive integer coefficients \(f_i\le Y^\eta\) give \(O_\varepsilon(Y^{\theta+\varepsilon})\) starts, where \(\theta=77/320\) and the three terms belong to the respective \(f_i\) times the powerful integers. For fixed \(k\) the possible coefficients from the repair form a finite set. Apply the lemma with start \(t-2\) and \(Y=(2X)^{1/k}\), increasing the threshold for \(X\) to fit \(Y^\eta\). Using it at \(k=3,5\), and the elementary middle-term bound for \(k\ge7\), gives

\[
C(X)=F(X)+O_\varepsilon(X^{77/960+\varepsilon}).
\]

The [quadratic companion](sources/plymesser-powerful-quadratics-hexagon2611.00699-v1.0.1.pdf), Theorem 1.1, states the exponent \(293/1220+\varepsilon\) for powerful values of a fixed positive-leading integer quadratic with nonsquare discriminant. For \(k=3\), the powerful integer \(T_3=t(4t^2-3)\) has factors with gcd dividing 3. If \(3\nmid t\), then \(4t^2-3\) is powerful. If \(t=3u\), then \(4t^2-3=3(12u^2-1)\), and \(12u^2-1\) is coprime to \(3t\), hence powerful. Both quadratics have discriminant 48. The bound at \(k=3\) becomes \(O_\varepsilon(X^{293/3660+\varepsilon})\), and the \(k=5\) use of Lemma 8.3 and the elementary \(k\ge7\) bound are smaller. This yields the stronger transfer

\[
C(X)=F(X)+O_\varepsilon(X^{293/3660+\varepsilon}).
\]

The two exponents in this section are conditional on their respective Plymesser inputs. The coefficient and discriminant hypotheses are checked here. The [source audit](PLYMESSER_2026_POWERFUL_COUNTS_SOURCE_SCOPE.md) records which proof passages and formal hypotheses have been inspected.

## The only maximizer of the pairwise bound

Use \(x_1,x_2,x_3\in[0,1/3]\) for the dyadic exponents of the three squarefree kernels, measured at center height \(X\). For a pair \((a,b)\) of exponents put \(s=a+b\), \(h=1-3s/2\), \(d=|a-b|\), \(w=s-1/3\), \(m=(9/8)sh=3/16-(27/16)w^2\). Define \(\Phi(a,b)\) as the minimum of \(s,h\) and, on the gate \(s,h\ge1/10\), \(a,b\ge1/100\), also \(B_{\rm red}\). This is the \(\Phi_3\) of [Plymesser, §§6–7](sources/plymesser-powerful-products-hexagon2611.00602-v1.0.1.pdf). In this notation the exact formula for its third term is

\[
B_{\rm red}=\theta-\frac{81}{64}w^2
   +\frac3{20}(d-m)_+ +\frac1{10}(3d/2-m)_+,
\]

where \(v_+=\max(v,0)\). Indeed the two reduced-coordinate maxima are \((s-m)/2+(d-m)_+/2\) and \((h-m)/2+(3d/2-m)_+/2\), with weights \(3/10\) and \(1/5\).

For \(0\le\eta\le1/1600\), suppose every pair has \(\Phi\ge\theta-\eta\). Sort the three exponents as \(a\le b\le c\). If \(a<1/100\), then \(b>0.23\) and \(b+c<0.51\) by the bounds on pair sums from \(s,h\ge\theta-\eta\). Thus \((b,c)\) lies in the gate, with \(s>0.46\), \(h\ge0.24\), \(m>0.12\), and \(d<0.05\). Both positive parts vanish, contradicting \(B_{\rm red}\ge\theta-\eta\) since \(s>0.46\). Therefore all pairs lie in the gate.

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

## Point count required on the surfaces

Linearly independent auxiliary forms modulo \(G_0\) may have a common divisor on \(V\). The rank bound leaves open whether the base locus is a union of curves, and a curve cover would still require a uniform arithmetic point count including exceptional components. The coefficient bound in the determinant proof is polynomial in \(X\) with a degree-dependent exponent, so a polynomial dependence on the coefficients in a point-counting theorem cannot be absorbed into \(X^\varepsilon\).

For an explicit sufficient estimate take \(\rho=1/1000\), \(\mu_- =\mu_+=53/500\). These values satisfy the strict determinant inequality, and every candidate in the box is fundamental for large \(X\), since \(q\ge X^{1/3-2\rho}\) and nonfundamental solutions have \(q\ll X^{2/9}\). A uniform \(O_\varepsilon(X^{7/250+\varepsilon})\) bound per cell on the points of its auxiliary surface cover would give \(O_\varepsilon(X^{6/25+\varepsilon})\) in the box. The per-cell estimate remains open. To improve the transferred \(293/1220\) bound globally along this route also requires controlling the near-binding boxes outside the chosen joint cover. Neither counting route settles existence of fundamental or nonfundamental candidates.
