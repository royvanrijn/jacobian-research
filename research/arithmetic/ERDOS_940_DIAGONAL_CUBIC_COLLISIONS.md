# Erdős #940 and the collision fibers of \(x^3+2y^3+4z^3\)

Let \(F(\mathbf x)=x_1^3+2x_2^3+4x_3^3\), \(w=(1,2,4)\), and let \(r_P(n)\) count its representations in \([P,2P]^3\). The required scale is \(P^3\). Indeed \(\sum_n r_P(n)\asymp P^3\), so an upper bound \(O(P^3)\) for \(\sum_n r_P(n)^2\) implies \(\gg P^3\) distinct values by Cauchy–Schwarz. Scaling all roots by \(2\) takes every value \(n\) to \(8n\), with summands \(8x_1^3,16x_2^3,32x_3^3\), all cube-full. This gives the implication to the density-zero assertion in Erdős #940 without assuming that \(2\) and \(4\) themselves are cube-full. The density deficit in one progression established in Beyer de Ryke's preprint concerns the complement and is compatible with positive upper density of this fixed form. [Erdős #940](https://www.erdosproblems.com/940) · [Beyer de Ryke, Theorem 1.1 and Corollary 1.2](https://arxiv.org/abs/2609.35772)

## The primitive difference and the exceptional fibers

For a collision \(F(\mathbf u)=F(\mathbf v)\) with \(\mathbf u\ne\mathbf v\), write \(\mathbf u-\mathbf v=g\mathbf d\) with \(g>0\) and \(\mathbf d\) primitive, and \(\mathbf s=\mathbf u+\mathbf v\). Then
\[
 Q_{\mathbf d}(\mathbf s)=m,\qquad
 Q_{\mathbf d}(\mathbf s)=\sum_i w_id_i s_i^2,\qquad
 m=-g^2 F(\mathbf d)/3.                                           \tag{1}
\]
The parity is \(s_i\equiv gd_i\bmod2\), and the exact inequalities are \(2P+g|d_i|\le s_i\le 4P-g|d_i|\). The identity follows from \(4(u_i^3-v_i^3)=(u_i-v_i)(3s_i^2+(u_i-v_i)^2)\). Conversely those conditions recover exactly one ordered pair from \((g,\mathbf d,\mathbf s)\). If \(H=\|\mathbf d\|_\infty\), then \(gH\le P\). We focus on \(g=1\) and a dyadic shell \(H\asymp P\).

**Theorem.** For \(d_1d_2d_3\ne0\), the fibers (1) have nonzero determinant and nonzero target. If a fiber has a rational affine line, then
\[
 S(\mathbf d)=6d_1d_2d_3F(\mathbf d)
\]
is a square. For each fixed signed squarefree class \(D\ne0\), there are \(O_\varepsilon(H^{38/21+\varepsilon})\) directions in \(\|\mathbf d\|_\infty\le H\) with \(S(\mathbf d)\) of class \(D\), uniformly for \(|D|\ll H^6\), even without primitivity. Each fiber at \(gH\le P\) has \(O_\varepsilon(P^{1+\varepsilon})\) points uniformly. The square-discriminant fibers contribute \(O_\varepsilon(P^{59/21+\varepsilon})\) collisions in total over all \(g,\mathbf d\). Collisions with at least one equal coordinate, excluding the exact diagonal, contribute \(o(P^3)\) in total.

The terms \(d_1^3,2d_2^3,4d_3^3\) have \(2\)-adic valuations in distinct classes modulo \(3\) when nonzero. Hence a nonzero \(\mathbf d\) has \(F(\mathbf d)\ne0\). Homogenization gives a rank-four quadratic form of determinant \(-m\det Q_{\mathbf d}=(2g/3)^2 S(\mathbf d)\). A rational affine line gives a rational totally isotropic plane for that form, which makes its determinant a square. This implication does not require any assertion about isotropy given the square condition.

For the square-class counts, fix a signed squarefree integer \(D\ne0\) with \(|D|\ll H^6\) and consider primitive \(\mathbf d\) for which \(S(\mathbf d)=D z^2\) for an integer \(z\ne0\). Write
\[
 |d_i|=e_i x_i t_{ij}t_{ik}u_i^2,\qquad \{i,j,k\}=\{1,2,3\},
\]
where \(e_i\) is the part of the squarefree kernel supported on primes dividing \(6D\), and \(t_{ij}=t_{ji}\) is the common part of the squarefree kernels of \(|d_i|,|d_j|\) away from those primes. The \(x_i\) and the three \(t_{ij}\) are pairwise coprime and squarefree. For a prime \(p\nmid 6D\) dividing exactly two coordinates, primitivity makes the third coordinate and \(F(\mathbf d)\) units at \(p\). The square-class condition forces the two valuations to have the same parity. Consequently a prime dividing \(x_i\) divides only \(d_i\) among the three coordinates and occurs there to odd valuation. It must divide \(F(\mathbf d)\). The entire squarefree kernel of \(|F(\mathbf d)|\) away from \(6D\) is exactly \(x_1x_2x_3\), by parity at each such prime. Thus
\[
 |F(\mathbf d)|=e_F x_1x_2x_3 h^2,\qquad
 |d_1d_2d_3|=e_1e_2e_3 x_1x_2x_3
              (t_{12}t_{13}t_{23}u_1u_2u_3)^2,
\]
where \(e_F\) is squarefree and supported on primes dividing \(6D\).

Fix the signs and the four \(e\)'s, with \(O_\varepsilon(H^\varepsilon)\) possible choices uniformly in \(D\), by the divisor bound for \(6|D|\). Use dyadic sizes \(X_i\le x_i<2X_i\) and \(T_{ij}\le t_{ij}<2T_{ij}\), and let \(M=X_1X_2X_3\). Select \(j\) for which \(X_j\) is largest, with remaining indices \(i,k\). For fixed \(t_{12},t_{13},t_{23}\), there are \(O(\sqrt{HX_\ell/(T_{\ell m}T_{\ell n})})\) possible \((x_\ell,u_\ell)\) for either \(\ell=i,k\), where \(m,n\) are the other two indices. This follows by taking \(O(X_\ell)\) choices of \(x_\ell\) and using \(|d_\ell|\le H\) to bound \(u_\ell\). Once \(d_i,d_k\) are fixed, \(x_j\) divides the nonzero integer \(w_i d_i^3+w_k d_k^3\). Its two terms have distinct \(2\)-adic valuations. The divisor bound gives \(O_\varepsilon(H^\varepsilon)\) choices of \(x_j\), and each allows at most \(O(\sqrt{H/(X_jT_{ij}T_{jk})})\) choices of \(u_j\). Summing over \(O(T_{12}T_{13}T_{23})\) choices of the \(t\)'s cancels the product of the three \(T\)'s in the denominator. The number of directions in the dyadic block is at most
\[
 \ll_\varepsilon H^{3/2+\varepsilon}\sqrt{X_iX_k/X_j}
 \le H^{3/2+\varepsilon} M^{1/6}
\]
up to a constant depending on \(\varepsilon\).

A second estimate counts directions on the projective cubics \(F(\mathbf d)=\lambda d_1d_2d_3\), where \(\lambda=F(\mathbf d)/(d_1d_2d_3)\). In this block both \(h\) and \(t_{12}t_{13}t_{23}u_1u_2u_3\) are \(O(\sqrt{H^3/M})\). Their squared ratio determines \(\lambda\) after the signs and four \(e\)'s are fixed. This gives \(O(H^3/M)\) possible \(\lambda\)'s. The cubic for each rational \(\lambda\) is smooth unless \(\lambda=6\). Indeed, a common zero of the three partial derivatives in projective space must have all coordinates nonzero, and multiplying the derivative equations gives \(\lambda^3=216\). At \(\lambda=6\) the cubic is the norm of \(d_1+d_2\theta+d_3\theta^2\) in \(\mathbb Q(\theta)\), where \(\theta^3=2\). The norm vanishes only at the zero vector. For every other \(\lambda\), the smooth plane cubic is absolutely irreducible. Clearing denominators and applying Heath-Brown's coefficient-uniform bound gives \(O_\varepsilon(H^{2/3+\varepsilon})\) primitive directions per curve, including both signs. The second block bound is therefore \(O_\varepsilon(H^{11/3+\varepsilon}/M)\). The minimum of the two bounds is \(O_\varepsilon(H^{38/21+\varepsilon})\), by splitting at \(M=H^{13/7}\). Summing over the dyadic sizes, signs and \(e\)'s absorbs their total cost in \(H^\varepsilon\). The same bound holds without primitivity, since \(S(k\mathbf d)=k^6 S(\mathbf d)\) and \(\sum_{k\ge1} k^{-38/21}\) converges. [Heath-Brown, Theorem 3 and (1.9)](sources/heath-brown-curves-math0405392v1.pdf)

To bound a nonempty fiber, fix the midpoint variable outside a same-sign pair of coefficients, and change the sign of the whole equation if needed. The resulting binary equation has positive coefficients \(A,B\) and target \(N\). If \(N>0\), the expression \(A s_i+s_j\sqrt{-AB}\) is an integer of norm \(AN\) in an imaginary quadratic field. Each principal ideal of that norm yields at most six generators, and the number of ideals of norm \(AN\) is at most \(\tau(AN)\), uniformly in the field. When \(N=0\) at most one pair occurs. Here \(AN\ll P^4\), giving the divisor bound after summing the fixed coordinate. For nonsquare \(S\), Corollary 1.3 of Browning–Verzobio gives the independent uniform bound
\[
 O_\varepsilon\!\left(P^{7/6+\varepsilon}H^{-1/6}+P^{1/2+\varepsilon}\right).
\]
Its coefficient size is \(\max_i|w_i d_i|\asymp H\) after removing a common factor at most \(4\) by primitivity. At \(H\asymp P\) the first term is still \(P^{1+\varepsilon}\). [Browning–Verzobio, equation (8) and Corollary 1.3](https://arxiv.org/pdf/2408.11453)

For a dyadic shell \(H/2<\|\mathbf d\|_\infty\le H\), the square-class bound for \(D=1\), the uniform fibre count, and the \(O(P/H)\) choices of \(g\) give \(O_\varepsilon(P^{2+\varepsilon}H^{17/21+\varepsilon})\) ordered collisions with square \(S\). Summing the shells gives \(O_\varepsilon(P^{59/21+\varepsilon})\). More generally, allowing \(O(H^\gamma)\) prescribed classes per shell gives \(O_\varepsilon(P^{59/21+\gamma+\varepsilon})\) for \(\gamma\ge0\). This is \(o(P^3)\) for \(\gamma<4/21\) and sufficiently small \(\varepsilon\).

The divisor bound also controls small primitive differences. At \(g=1\), the nondegenerate directions with \(\|\mathbf d\|_\infty\le P^{2/3-\delta}\) contribute \(O_\varepsilon(P^{3-3\delta+\varepsilon})\). Without fixing \(g\), a shell of primitive size \(H\) has \(O(P/H)\) possible multipliers and contributes \(O_\varepsilon(P^{2+\varepsilon}H^2)\). Summing for \(H\le P^{1/2-\delta}\) gives \(O_\varepsilon(P^{3-2\delta+\varepsilon})\).

If one coordinate agrees, the other coordinates solve \(A(x^3-x'^3)+B(y^3-y'^3)=0\), with \(A/B\) a noncube rational power of 2. For an off-diagonal collision no two signed terms cancel, since cancelling the two terms of the same coefficient forces the other coordinates also to agree. Cross-cancellation would require \(A/B\) to be a cube. Salberger's Corollary 0.7 counts the primitive quadruples by \(O(P^{\sqrt3}(\log(2P))^4)\). Including dilations has the same order because \(\sqrt3>1\). The equal coordinate contributes \(O(P)\), giving \(O(P^{1+\sqrt3}(\log(2P))^4)=o(P^3)\). The diagonal contributes \(|[P,2P]\cap\mathbb Z|^3\). [Salberger, Corollary 0.7](https://londmathsoc.onlinelibrary.wiley.com/doi/full/10.1112/plms.12508)

### Small and unbalanced primitive differences

For fixed \(g\), the number of collisions with three nonzero primitive difference coordinates and \(|d_i|\le H_i\) for \(1\le H_1\le H_2\le H_3\le P/g\) is
\[
 \ll_\varepsilon P^{3/2+\varepsilon}(H_1H_2)^{3/4}H_3^{1/2}.       \tag{5}
\]
The bound holds after permuting the coordinate labels and also if the vectors are allowed to be nonprimitive while \(g\) remains fixed. For a balanced box and fixed \(g\) this is \(O_\varepsilon(P^{3/2+\varepsilon}H^2)\). Hence \(H\le P^{3/4-\delta}\) contributes \(o(P^3)\) for each fixed \(g\). Across all \(g\), the collisions with \(0<\min_i|d_i|\le K\) and three nonzero differences contribute \(O_\varepsilon(P^{11/4+\varepsilon}K^{3/4})\). In particular one can discard \(K\le P^{1/3-\delta}\).

To prove (5), let \(e(t)=\exp(2\pi i t)\), \(I_h=\{x\in\mathbb Z:x,x+gh\in[P,2P]\}\), and
\[
 G_{g,H}(\alpha)=
 \sum_{0<|h|\le H}\ \sum_{x\in I_h}
 e\bigl(\alpha((x+gh)^3-x^3)/g\bigr).
\]
The Fourier coefficient at a nonzero integer \(n\) counts choices of a signed divisor \(h\mid n\) and at most one \(x\in I_h\) for that \(h\), by monotonicity of the difference of cubes on \(I_h\). There is no zero coefficient. The divisor bound and the \(O(PH)\) terms give \(\int_0^1 |G_{g,H}|^2\ll_\varepsilon P^{1+\varepsilon}H\). Cauchy–Schwarz in \(h\) gives \(|G_{g,H}|^2\le 2H(D+U(\alpha))\), where \(D=\sum_h|I_h|\ll PH\) and \(U\) is the sum over \(h\) and \(x\ne y\) in \(I_h\) with phase \(3h(x-y)(x+y+gh)\alpha\). The three nonzero factors in each phase determine \(h,x,y\), and \(x+y+gh>0\). All Fourier coefficients of \(U\) are \(O_\varepsilon(P^\varepsilon)\) by the ternary divisor bound, since the phase integers have size \(O(P^3)\). Therefore orthogonality gives
\[
 \int_0^1 |G_{g,H}|^4
 \ll H\left(D\int_0^1 |G_{g,H}|^2+
               \int_0^1 |G_{g,H}|^2U\right)
 \ll_\varepsilon P^{2+\varepsilon}H^3.
\]
In the last step the integral of \(|G_{g,H}|^2U\) is at most \(O_\varepsilon(P^\varepsilon)(\sum_n\widehat G_{g,H}(n))^2\), using the nonnegative Fourier coefficients. Finally orthogonality counts collisions by the integral of \(\prod_i G_{g,H_i}(w_i\alpha)\), before imposing primitivity. Hölder with exponents \(4,4,2\) gives (5), since multiplication by each integer \(w_i\) preserves the circle integral. Taking the two larger widths to be \(P/g\) gives the bound for the minimum coordinate after summing \(\sum_g g^{-5/4}\).

## Local densities and shared divisors

Now take \(g=1\), \(3\mid F(\mathbf d)\) and \(d_1d_2d_3\ne0\). Other primitive directions either give no fiber or were handled above. In base variables \(\mathbf v\) use the integral polynomial
\[
 T_{\mathbf d}(\mathbf v)=\sum_i w_i d_i (v_i^2+d_i v_i)+F(\mathbf d)/3.
\]
Then \(4T_{\mathbf d}(\mathbf v)=Q_{\mathbf d}(2\mathbf v+\mathbf d)-m\). Let \(\beta_p\) be the limit of \(p^{-2k}\) times the number of zeros of \(T_{\mathbf d}\bmod p^k\). For odd \(p\) this is the midpoint density without any parity restriction. At \(p=2\), some \(d_i\) is odd. For that variable the derivative of \(d_i v_i(v_i+d_i)\) is odd, and each equation after dividing the coefficient \(w_i\) has at most two solutions modulo each \(2^k\). It follows that \(\beta_2\le8\). At \(3\), primitive \(\mathbf d\) satisfying \(3\mid F\) has at least two unit coordinates. Cubes of units are \(\pm1\bmod9\); checking the three cases of a single zero coordinate modulo \(3\), and the case of three units, gives \(v_3(F)=1\). Therefore \(m\) is a unit and every zero modulo \(3\) is nonsingular in some midpoint variable. In particular \(\beta_3\le3\).

At any prime \(p\ge5\) the following bounds are uniform:

| Coefficients divisible by \(p\) | Condition on \(m\) | Density or bound |
| --- | --- | --- |
| None | \(p\nmid m\) | \(1+\left(\frac{S(\mathbf d)}p\right)/p\) |
| None | \(p\mid m\) | At most \(1+1/p\) |
| Exactly one, say \(w_i d_i\) | \(p\nmid m\) | \(1-\left(\frac{-w_j w_k d_j d_k}p\right)/p\) |
| Exactly one, say \(w_i d_i\) | \(p\mid m\) | At most \(2+\min(v_p(d_i),v_p(m))\) |
| Exactly two | \(m\) is a unit | \(1+\left(\frac{-3}p\right)\) |

The unit-target formulas follow from the finite-field counts for a quadratic form of rank three, two or one, followed by lifting. For rank three at a target divisible by \(p\), only the origin is singular. It contributes zero unless \(p^2\mid m\), when its contribution is \(p^{-1}\) times the density with target \(m/p^2\). The nonsingular contribution is \(1-p^{-2}\), yielding the bound by induction. In rank two at target zero modulo \(p\), the nonsingular contribution is \((1+\chi)(1-1/p)\) for \(\chi=\left(\frac{-w_j w_k d_j d_k}p\right)\). In this family \(\chi=1\) when \(p\mid d_i,F(\mathbf d)\) and \(d_jd_k\) is a unit. If the two unit variables are multiples of \(p\) and the third coefficient has valuation at least two, dividing by \(p^2\) lowers that valuation and the target valuation by two with density multiplier \(1\). At valuation one the remaining forms alternate, on their singular parts, between valuation patterns \((0,0,1)\) and \((0,1,1)\), with respective multipliers \(1/p\) and \(1\) and with a drop of one in the target valuation. Their densities are at most two by induction. This gives the stated bound for each step of size two. In rank one \(m/(w_i d_i)\equiv-d_i^2/3\) at the single unit coefficient, giving the last formula. These recursions are finite because \(m\ne0\).

**Theorem.** On the primitive nonsquare fibers with \(3\mid F(\mathbf d)\) and \(d_1d_2d_3\ne0\), the singular series has bounded fixed moments on the scale of the difference box:
\[
 \sum_{\|\mathbf d\|_\infty\le H}\mathfrak S(\mathbf d)^j\ll_j H^3
 \qquad(j\ge1\text{ fixed}).                                    \tag{2}
\]
Here \(\mathfrak S(\mathbf d)=\prod_p\beta_p\), interpreted by factoring out the primitive quadratic \(L(1,\chi)\) attached to the square class of \(S(\mathbf d)\). For \(H\asymp P\) the associated real density is \(O(P/H)\). Hence the total real-density times singular-series expression on those fibers is \(O(P^3)\).

To prove (2), let \(G_i\) be the part of \(\gcd(|d_i|,|F(\mathbf d)|)\) supported on primes \(p\ge5\) not dividing the other two coordinates. The large factors in the table are bounded by
\[
 W(\mathbf d)=\prod_i\tau(G_i)^2\prod_{i<k}\tau(\gcd(|d_i|,|d_k|)).
\]
Each factor here has bounded fixed moments on the scale \(H^3\). For pair gcds, the condition \(q\mid d_i,d_k\) gives \(O(H^3/q^2)\) directions, using nonzero coordinates and \(\lfloor H/q\rfloor\). For \(q\mid G_i\), there are \(O(H/q)\) choices for \(d_i\) and \(O(H)\) choices for one of the other coordinates, a unit modulo \(q\). The remaining unit coordinate satisfies a cube-root congruence with at most \(3^{\omega(q)}\) roots. Since \(q\le H\), this gives \(O(3^{\omega(q)}H^3/q^2)\) directions. Expanding each fixed positive power of the divisor function as a sum over divisors gives coefficients of size polynomial in the valuation at each prime. The resulting Euler series with \(q^{-2}\) converge. Hölder across the six factors yields bounded moments of \(W\).

With \(D_0\) the signed squarefree kernel of \(S(\mathbf d)\), write \(\chi\) for the primitive character of the corresponding fundamental discriminant. Apart from the bounded factors at \(2,3\), the table gives
\[
 \mathfrak S(\mathbf d)
 \ll L(1,\chi)\, W(\mathbf d)\,
       \prod_{p\mid S(\mathbf d)}(1+1/p)^2.                       \tag{3}
\]
Indeed at \(p\nmid 6S\) the local factor times \(1-\chi(p)/p\) equals \(1-p^{-2}\). At \(p\mid S\) the table gives the factors on the right. The finite product in (3) has bounded fixed moments. For \(p>c\log H\), the bound \(|S|\ll H^6\) gives \(\sum_{p\mid S,\ p>c\log H}1/p=O(1/\log\log H)\). For small \(p\), the mod-\(p\) proportion with \(p\mid S\) is \(O(1/p)\), by the elementary zero count for a nonzero polynomial of fixed degree. The Chinese remainder theorem and a small fixed \(c>0\), so that \(\prod_{p\le c\log H}p\le H^{1/2}\), give a bound by the product of \(1+O_j(p^{-2})\), up to an absolute constant from counting each residue class in the box.

It remains to bound the \(L\)-factor moments in this thin family. Apply Granville–Soundararajan, Proposition 2.2, with conductor bound \(Q\asymp H^6\), parameter \(A=100\), and \(y=(\log Q)^{8A}\). All but \(Q^{2/A}\) primitive characters have \(L(1,\chi)\ll\prod_{p\le y}(1-\chi(p)/p)^{-1}\). For each exceptional real character the fixed-square-class bound above is \(O_\varepsilon(H^{38/21+\varepsilon})\), and \(L(1,\chi)\ll\log(2H)\) by summing the Dirichlet series through one period and partial summation in the tail. Their total fixed moment is \(o(H^3)\). The same exceptional classes contribute \(o(P^3)\) to the actual collision count, even summed over all \(g\) and dyadic \(H\), since \(Q^{2/A}\ll H^{3/25}\) and \(3/25<4/21\). For all other characters, Mertens bounds the product over \(c\log H<p\le y\) by an absolute constant. At smaller primes majorize the factor by \((1-\left(\frac{S}p\right)/p)^{-1}\) if \(p\nmid S\), and \((1-1/p)^{-1}\) if \(p\mid S\).

For \(p\ge5\), the mean of \(\left(\frac{S(\mathbf d)}p\right)\) over \(\mathbb F_p^3\) is \(O(1/p)\). One elementary proof changes variables \(u_i=w_i d_i^3\ne0\) and expands the multiplicity of a cube using characters \(\rho_i\) of order dividing three. Since \(\chi_p(d_i)=\chi_p(u_i/w_i)\) for the Legendre symbol \(\chi_p\), this produces at most \(27\) complete sums of the form
\[
 \sum_{u_i\ne0}(\chi_p\rho_1)(u_1)(\chi_p\rho_2)(u_2)
             (\chi_p\rho_3)(u_3)\chi_p(u_1+u_2+u_3).
\]
Expanding the last character using its Gauss sum gives three nonprincipal Gauss sums of modulus \(\sqrt p\), divided by one of the same modulus, and an outer sum of modulus at most \(p-1\). This gives \(O(p^2)\) for each. Thus the \(j\)-th moment of the small-prime majorant at each prime is \(1+O_j(p^{-2})\), also accounting for the \(O(1/p)\) frequency of \(p\mid S\). The same small-modulus CRT bound applies. This proves bounded \(L\)-moments, and Hölder in (3) proves (2). [Granville–Soundararajan, Proposition 2.2](https://arxiv.org/pdf/math/0206031v1)

For the real density in (2), integrate in a coordinate with \(|d_i|=H\) on the intersection where both triples are in \([P,2P]^3\). The derivative of \(T_{\mathbf d}\) in that coordinate has magnitude \(|w_i d_i(2v_i+d_i)|\ge2HP\) and constant sign. The other two coordinates range over a region of area at most \(P^2\). This gives \(O(P/H)\) by the coarea formula. The result controls the arithmetic and real expression even at the critical scale, with no assertion that it uniformly dominates the actual count.

## The lattice discrepancy at the critical scale

**Theorem.** For \(g=1\) and \(H\asymp P\), a positive-density family of balanced primitive directions admits a congruence relaxation modulo \(R\asymp H^3\), with zero-frequency contribution \(O(P^3)\) after summing over the family. The relaxation has an explicit nonzero-frequency remainder, for which the separate square-root estimates of the local conic sums do not yield \(O(P^3)\).

Choose \(|d_i|\in[H/2,H]\) and a sign pattern that includes both signs. Use \(\mathbf d\equiv(-2,3,-1)\bmod72\), a difference realized by the collision \((3,4,1),(5,1,2)\). Require the parts \(q_i\) of \(|d_i|\) away from \(2,3\) to be squarefree, pairwise coprime and coprime to \(F(\mathbf d)\). At each \(p\ge5\), requiring \(p^2\nmid d_i\), no pair both divisible by \(p\), and \(p\nmid(d_i,F)\) forbids \(O(p^{-2})\) of the residue classes modulo \(p^2\). There are permitted classes at every prime, such as those with all coordinates units. The total of the conditions omitted above any fixed cutoff \(K\) removes \(O(H^3/K)\) vectors by the same divisibility count used for (2), including the bound \(O(H^3/p^2)\) whenever \(p^2\mid d_i\). Apply the Chinese remainder theorem first for the fixed primes below \(K\) and then let \(K\) grow. The infinite product of the local proportions is positive. The congruence modulo \(72\) makes the \(2,3\)-parts bounded and makes \(\mathbf d\) primitive at those primes, so each \(q_i\asymp H\).

Write \(R=q_1q_2q_3\). Let \(V\) be a fixed smooth nonnegative compactly supported weight on base variables scaled by \(P\), equal to at least \(1\) on \([1,2]^3\). Every collision in the box is counted by the weighted congruence \(T_{\mathbf d}(\mathbf v)\equiv0\bmod R\). Its exact value is
\[
 \frac{P^3}{R^3}\sum_{\mathbf k\in\mathbb Z^3}
                 \widehat V(P\mathbf k/R)\, C_{\mathbf d}(\mathbf k),
 \qquad
 C_{\mathbf d}(\mathbf k)=\sum_{\substack{\mathbf v\bmod R\\T_{\mathbf d}(\mathbf v)=0\bmod R}}
                                      e(\mathbf k\cdot\mathbf v/R).          \tag{4}
\]
The zero-frequency term is \(\widehat V(0) P^3/R\) times the product of the rank-two unit-target densities in the table. It has bounded average because \(\prod_{p\mid R}(1+1/p)\) does, by \(\sum_{q\mid d_i}\mu^2(q)/q\) and the divisor count in the box. At a prime dividing \(q_i\), \(v_i\) is free. Thus \(C_{\mathbf d}(\mathbf k)\) vanishes unless \(q_i\mid k_i\) for all \(i\).

At the effective frequency height \(R/P\asymp P^2\), there are order \(P^3\) lattice frequencies satisfying those three divisibilities. For a local frequency nonzero on the two conic coordinates, completion by quadratic Gauss sums gives \(p\) times a Kloosterman sum, of magnitude at most \(2p^{3/2}\), including the free coordinate. This follows by the additive-character delta expansion for the diagonal conic with unit target. When the local frequency vanishes on the conic, the factor is of order \(p^2\). Even on the frequencies nonzero on the conic at every prime of \(R\), the direct triangle bound from the former estimates is \(O_\varepsilon(P^{3/2+\varepsilon})\) per fiber on the effective range. For the density goal, the required upper bound for the sum of the signed nonzero-frequency contributions in (4) is \(O(P^3)\) across the family. Equation (4) gives a congruence upper bound for exact solutions, not an asymptotic for those solutions. The discrepancy is at a modulus of order \(P^3\) with all three physical windows of length \(P\). [Weil, On Some Exponential Sums](https://doi.org/10.1073/pnas.34.5.204)

### Kloosterman expansion by divisors

On the same balanced family, the nonzero-frequency terms from reduced denominators \(r\le P^{2-\delta}\) total \(O_{\delta,\varepsilon}(P^{3-\delta/2+\varepsilon})\) for \(0<\delta<1\). Terms with exactly one nonzero frequency coordinate total \(O(P^3)\) in absolute value over all denominators and directions. With a product weight, terms with exactly two nonzero frequency coordinates total \(O(P^3\log^2(2P))\) after summing the two frequencies with their signs and retaining the full transform at the inactive modulus.

Use midpoint variables and take \(W=\prod_i W_i\), where the fixed smooth nonnegative factors are compactly supported and at least \(1\) on \([2,4]\). The count weighted by \(W(\mathbf s/P)\) for \(Q_{\mathbf d}(\mathbf s)\equiv m\bmod R\) is an upper bound for actual collisions. The prime divisors of \(R\) are all odd and coprime to \(m\), and each divides exactly one coefficient \(A_i=w_i d_i\). For \(r\mid R\) put \(r_i=(r,q_i)\), \(L_i=r/(Pr_i)\), and define
\[
 B_r(\boldsymbol\ell)=\sum_i r_i\,\overline{A_i/r_i}\,\ell_i^2\bmod r,
 \qquad
 \kappa_{\mathbf d}(r)=\prod_{p\mid r}\left(\frac{-A_jA_k}{p}\right),
\]
where the inverse in summand \(i\) is taken modulo \(r/r_i\), and \(j,k\) are the two unit coefficients at \(p\). With \(\mathrm{Kl}_r(a,b)=\sum_{t\bmod r}^{*}e((at+b\bar t)/r)\), the weighted congruence count is exactly
\[
 \frac{P^3}{R}\sum_{r\mid R}\frac{\kappa_{\mathbf d}(r)}{r}
    \sum_{\boldsymbol\ell\in\mathbb Z^3}
    \widehat W(\ell_1/L_1,\ell_2/L_2,\ell_3/L_3)
       \mathrm{Kl}_r(-m,-\overline 4 B_r(\boldsymbol\ell)).          \tag{6}
\]
Here the modulus-one Kloosterman sum is \(1\). To see (6), insert the additive delta for the congruence in the Poisson expression and partition the delta parameter by its reduced denominator \(r\), writing the parameter as \(t/r\). The complete sum in coordinate \(i\), initially modulo \(R\), vanishes unless its frequency is \((R/r)r_i\ell_i\). In that case it is \((R/r)r_i\) times the Gauss sum modulo \(r/r_i\) with coefficients \(tA_i/r_i,\ell_i\). Completing the square in each of the three sums gives the phase \(-\bar t\,\overline4 B_r(\boldsymbol\ell)/r\). At every \(p\mid r\), the two quadratic Gauss factors contribute \(p(\frac{-A_jA_k}{p})\) with no dependence on \(t\). Thus the three coordinate sums have total prefactor \(R^3\kappa_{\mathbf d}(r)/r\). Summing \(t\) gives (6). The terms with \(\boldsymbol\ell=0\) together give the zero-frequency term for the midpoint weight, which has total \(O(P^3)\) by the same local product bound as in (4).

For a nonempty support \(J=\{i:\ell_i\ne0\}\), the Weil bound and the rapid decay of the factors \(\widehat W_i\) give an absolute bound for the inner sum divided by \(r\) of
\[
 \ll 2^{\omega(r)}r^{-1/2}\prod_{i\in J} L_i.
\]
The last factor \(r^{-1/2}\prod_{i\in J}L_i\) is \(O(\sqrt r/P)\) for \(|J|\le2\), using \(r_i\le P\), and is \(r^{3/2}/P^3\) for \(|J|=3\). Since \(P^3/R\asymp1\) and \(R\ll P^3\), summing over \(r\le P^{2-\delta}\) and then over \(O(P^3)\) directions gives the stated power saving.

For the support of size one, use a second moment for square Kloosterman arguments. For \(p\) odd and \(p\nmid ab\), the Gauss identity
\[
 \mathrm{Kl}_p(a,b u^2)=\sum_{v\bmod p}\chi_p(v^2-4ab)e(vu/p)
\]
follows by opening the sum in \(t\) and evaluating the quadratic Gauss sums in \(v\) and \(t\). Thus the complete mean of \(|\mathrm{Kl}_p(a,bu^2)|^2/p\) is at most \(1\) by Parseval. Its nonzero normalized Fourier coefficients have size \(O(p^{-1/2})\). Parseval expresses each as \(p^{-1}\) times a quadratic character sum of \((v^2-4ab)((v+h)^2-4ab)\) for shift \(h\ne0\). The two pairs of roots have different sums, so the polynomial is nonsquare and Weil's character bound applies. The Chinese remainder theorem and Poisson summation now give, for odd squarefree \(s\) with \((s,ab)=1\) and a fixed nonnegative smooth rapidly decaying weight \(\psi\),
\[
 \sum_{u\in\mathbb Z}\psi(u/L)
       \frac{|\mathrm{Kl}_s(a,b u^2)|^2}{s}
       \ll_{\psi,\varepsilon} L+s^{1/2+\varepsilon}.               \tag{7}
\]
Indeed the complete mean is at most \(1\), and a Fourier mode of reduced denominator \(d\mid s\) has coefficient \(O(C^{\omega(d)}d^{-1/2})\) for an absolute \(C\). The sum of the nonzero modes of that denominator is \(O_\psi(C^{\omega(d)}\sqrt d)\), by \(L\sum_{j\ne0}|\widehat\psi(jL/d)|\ll_\psi d\). Summing \(d\mid s\) proves (7). [Weil, On Some Exponential Sums](https://www.pnas.org/doi/pdf/10.1073/pnas.34.5.204)

Suppose only \(\ell_i\ne0\), and let \(j,k\) be the other indices. The part of the Kloosterman sum modulo \(r_i\) is the Ramanujan sum \(\mu(r_i)\). Modulo \(r_j,r_k\) its factors have unit square arguments \(b\ell_i^2\) with \(b\) a unit. If either \(r_j\) or \(r_k\) is below \(P^{9/10}\), the pointwise Weil bound gives \(O_\varepsilon(P^{-1/20+\varepsilon}/r_i)\) for the sector, including the factor \(P^3/R\). This costs \(o(P^3)\) after summing all such divisors and directions. Otherwise \(L_i=r_jr_k/P\ge P^{4/5}\). Apply Cauchy–Schwarz and (7) at both moduli to bound the sector by \(O(\sqrt{r_jr_k}/(Pr_i))\). The sum over all divisor choices is bounded by a constant times
\[
 \sum_{r_i\mid q_i}\frac1{r_i}\ \prod_{\ell=j,k}\sum_{f\mid q_\ell} f^{-1/2}.
\]
This expression has bounded average over the direction box, since \(q_\ell\mid |d_\ell|\) and \(\sum_{1\le n\le P}\sum_{f\mid n} f^{-1/2}\ll P\). This proves the bound for one nonzero coordinate.

For two nonzero coordinates \(i,j\), keep the full conic transform at the inactive modulus \(q=q_k\). In the midpoint version of (4), write the frequencies as \(k_\ell=q_\ell n_\ell\) with \(n_i,n_j\ne0\) and \(n_k=0\). The three free-coordinate sums give a common factor \(R\). At \(q_j\), the transform on the conic depends only on \(n_i\), with additive phase \(e(\bar q n_i z/q_j)\) in coordinate \(z=s_i\). The projection has multiplicity
\[
 D_j(z)=\prod_{p\mid q_j}\bigl(1+\chi_p(a_j+b_jz^2)\bigr)
       =\sum_{s\mid q_j}\chi_s(a_j+b_jz^2),
\]
where \(a_j=m/A_k\), \(b_j=-A_i/A_k\) are units modulo \(q_j\), and \(\chi_s=\prod_{p\mid s}\chi_p\). Expand both active-coordinate projections, using \(s\mid q_j\) and \(t\mid q_i\). For the \(s\)-term the sum is zero unless \(n_i=(q_j/s)\ell_i\), when it equals \((q_j/\sqrt{s}) f_s(\ell_i)\), with
\[
 f_s(\ell)=s^{-1/2}\sum_{z\bmod s}\chi_s(a_j+b_jz^2)e(\bar q\ell z/s).
\]
The unitary finite Fourier transform of \(f_s\) has sup norm at most \(1\). The \(t\)-term for \(n_j=(q_i/t)\ell_j\) has the same form. At \(q\) the remaining complete sum is the full conic transform \(C_q(\bar s\ell_i,\bar t\ell_j)\), with \(C_q(x,y)=\sum_{A_i\xi^2+A_j\eta^2=m\bmod q}e((x\xi+y\eta)/q)\).

Write \(L_i=sq/P\), \(L_j=tq/P\) for the frequency widths in these variables. Aggregate \(\widehat W_i(\ell_i/L_i)f_s(\ell_i)\) with \(\ell_i\ne0\) by the residue \(\bar s\ell_i\bmod q\) to form a vector \(U\), and form \(V\) using \(t\) for the \(j\)-coordinate. The unitary Fourier transform of \(U\) has sup norm \(O(\sqrt{s/q})\). In fact, expand \(f_s\) by its finite transform, then apply Poisson summation in \(\ell_i\). Including \(\ell_i=0\) gives \(L_i/\sqrt{sq}\) times a sum of bounded Fourier coefficients multiplied by values of \(W_i\) on a translate of a lattice of spacing \(L_i/s\). Compact support leaves \(O(1+s/L_i)\) terms, yielding \(O(\sqrt{s/q})\) since \(q\le P\). Removing \(\ell_i=0\) changes the bound by at most \(O(|f_s(0)|/\sqrt{q})=O(\sqrt{s/q})\). The bound for \(V\) is \(O(\sqrt{t/q})\).

Fourier inversion in \(C_q\) now bounds the \(U,V\) pairing by \(q\) times the product of the two Fourier sup norms times the number of conic points modulo \(q\). This number is at most \(q\prod_{p\mid q}(1+1/p)\) by the rank-two unit-target formula. The prefactor for this pair \(s,t\) in Poisson summation is \(P^3/(Rq\sqrt{st})\), giving \(O((P^3/R)\prod_{p\mid q}(1+1/p))\). Sum over \(s\mid q_j,t\mid q_i\), then over directions. The two divisor counts have averages \(O(\log(2P))\) and the product over \(p\mid q_k\) has bounded average, giving \(O(P^3\log^2(2P))\). This sums the two nonzero frequencies with their signs and keeps the full transform at the inactive modulus, including its zero-delta terms.

The binary norm description isolates the same issue in a different variable. With \(A,B\asymp H\) the same-sign coefficients, each fixed third midpoint gives a norm \(A(m-C s_3^2)\), with the sign chosen to keep \(A,B>0\). The divisor bound counts all ideals with that norm, whereas the physical equation specifies a principal generator with prescribed coordinates. The quaternionic sieve available for the norm sumsets starts from a positive definite quaternary norm shell in a prescribed order and uses a Hecke bound for its residue statistics. The present fiber is an affine indefinite ternary equation with changing discriminant and target, with the third-coordinate window and the norm shape coupled. The shell bound provides no corresponding estimate for the generator condition or for (4). [Quaternionic sieve and its shell hypotheses](PRIMITIVE_NORM_HEEGNER_CORE_SUMSETS.md#quaternionic-sieve)

The six-variable collision cubic \(F(\mathbf u)-F(\mathbf v)=0\) is smooth in projective space, since its six partial derivatives vanish together only at the origin. In midpoint coordinates it is a diagonal Perazzo form with the additional term \(F(\mathbf d)\). Removing that term makes the whole midpoint plane singular and changes (1) to a family of cones. The Perazzo cancellation theorem is proved for a singular cubic of that type and does not give the lattice estimate here. [Browning–Munshi–Wang, Theorem 1.1 and §2](https://arxiv.org/pdf/2604.19045)

## Smooth inputs, marking, and bounded multiplicity

**Theorem.** For fixed smoothness exponent \(\eta>0\), requiring a fixed number of disjoint prime groups of harmonic mass bounded above and below, with the smallest marking prime tending to infinity and the largest at most \(P^{o(1)}\), can retain \(\gg_\eta P\) roots that are \(P^\eta\)-smooth. If instead the number of disjoint groups tends to infinity with order \(\log\log P\), their primes are at most \(\exp((\log P)^\delta)\) for a fixed \(\delta<1\), and each has bounded harmonic mass, requiring at least one factor in every group retains only \(o(P)\) roots, even without smoothness.

For the fixed groups each fixed product \(d\) of marked primes is \(P^{o(1)}\). Dickman's estimate gives \((\rho(1/\eta)+o(1))P/d\) smooth multiples, uniformly in such products. The first moment of the product of the falling factorials of the marked prime counts is \(\gg P\), because the distinct-prime reciprocal masses in each group are bounded below. Its second moment is \(O(P)\), using \(\lfloor P/[d,e]\rfloor\) and expanding the pair of lists by the number of shared primes in each of the finitely many groups. Cauchy–Schwarz gives the claimed positive proportion. [Hildebrand–Tenenbaum, equations (1.1) and (1.5)](https://tenenb.perso.math.cnrs.fr/PPP/Psi%2B.pdf)

For growing numbers of groups, count how many groups have no divisor of a uniform \(n\le P\). The probability of missing a group is its independent Euler product plus \(o(1)\), uniformly in the group. The same holds for any two groups jointly, with the product of their Euler products. To justify this, truncate inclusion–exclusion at adjacent odd and even integers of order \(\log\log P\). Since the total reciprocal mass for two groups is bounded, the remainders in the products are bounded by a factorial tail tending to zero. The number of squarefree divisors in the truncations is at most \(\exp(O((\log P)^\delta\log\log P))=P^{o(1)}\), so interval errors after normalization are \(o(1)\). The expected number of missing groups is proportional to the number of groups and its variance is smaller than the square of that number. Chebyshev gives the zero-missing probability \(o(1)\). The actual weighted-dilation theorem demands a positive number of distinct marks in every group. Its full support therefore cannot be imposed as a positive-density restriction on the roots. Marking a fixed number of groups leaves the discrepancy in (4). [Weighted-dilation preprint, Definition 3.1 and Theorem 5.1](sources/weighted-dilation-graphs-sep24-2026.pdf)

**Corollary.** The fixed form has \(\gg X^{1-0.24871567/3}\) distinct positive values up to \(X\). Take all three roots \(P^\eta\)-smooth with \(\eta>0\) small and let \(f(\alpha)\) denote their cubic Weyl sum. Orthogonality gives the collision energy as \(\int_0^1 |f(\alpha)f(2\alpha)f(4\alpha)|^2\,d\alpha\). Hölder bounds it by \(\int_0^1|f(\alpha)|^6\,d\alpha\). By orthogonality this last integral is at most Wooley's mixed moment with one pair unrestricted, which is \(O(P^{3+0.24871567})\). The number of smooth input triples is \(\asymp_\eta P^3\), so Cauchy–Schwarz and \(F\le7P^3\) give the stated exponent. [Wooley, Theorem 1.2](https://arxiv.org/pdf/1502.01944)

For a family with \(\asymp P^3\) input triples, \(\sum_n r(n)^{1+\delta}\ll P^3\) for any fixed \(\delta>0\) would give a positive proportion of input triples with bounded multiplicity. Indeed the mass of \(r(n)>K\) is at most \(K^{-\delta}\sum r(n)^{1+\delta}\). Conversely, \(\gg P^3\) occupied values already give \(\gg P^3\) triples with bounded multiplicity, since there are at most \((\sum r(n))/K\) values with \(r(n)>K\). Interpolating the known second moment with the first for \(0<\delta<1\) only gives \(O(P^{3+0.24871567\delta})\). A bound for an actual high-multiplicity tail is needed to remove that loss.

A weaker sufficient target concerns only nearby representations. For \(0<\delta<1\), let \(B(P,\delta)\) count input triples in \([P,2P]^3\) that have another representation of the same value within sup-norm distance \(\delta P\). For two triples at the same level, factoring the cube differences gives \(|z-z'|\le |x-x'|+2|y-y'|\) in this box. Partition the \((x,y)\) box into \(M_\delta=O(\delta^{-2})\) cells of side at most \(\delta P/3\). Two representations in the same cell are \(\delta P\)-close. Each level with \(r_P(n)>2M_\delta\) has at most \(M_\delta\) singleton cells, so at least half its representations are counted by \(B(P,\delta)\). Hence \(\sum_{r_P(n)>2M_\delta}r_P(n)\le 2B(P,\delta)\). A small enough upper bound for \(P^{-3}B(P,\delta)\) at a fixed \(\delta\) would prove positive density by retaining the triples in the bounded-multiplicity levels. The estimates in (5) for sublinear differences do not establish such a bound for fixed \(\delta\).

**Open question.** Prove positive upper density for the fixed \(F\). The immediate estimate missing in this approach is a constant-order average upper bound on actual nonsquare fibers at \(g=1,\ H\asymp P\). On the squarefree-coefficient subfamily, (6) leaves the signed discrepancy at large reduced denominators with three active frequencies, and a logarithmic loss for two active frequencies. The square-discriminant fibers and the sparse exceptional character classes contribute \(o(P^3)\) over all dyadic \(H\) and \(g\). The remaining nonsquare fibres still require a lattice count, including the ranges with \(g>1\), or a direct bound on the high-multiplicity tail. For a complete second-moment route, an \(O(PH^2)\) bound per dyadic \(H\) and per \(g\le P/H\) would sum to \(O(P^3)\), since \(\sum_H H\ll P\).

A finite enumeration on \([1,P]^3\) checked (1) for every collision. At \(P=160\) the distinct values numbered \(0.866706\,P^3\), the energy was \(1.30083\,P^3\), and \(0.999266\) of triples had multiplicity at most four. Among ordered collisions with \(g=1\), three nonzero differences and \(\|\mathbf d\|_\infty\ge P/4\), four had square \(S\). These finite figures only check the reduction and give examples in the exceptional locus.
