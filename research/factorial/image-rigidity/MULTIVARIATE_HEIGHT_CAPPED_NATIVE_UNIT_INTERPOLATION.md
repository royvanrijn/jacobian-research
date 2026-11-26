# Arbitrary unit residues at a native factorial mark

A prime dividing an actual factorial ramp can force a quadratic reduction of one parent. Even when the prime is the entire primitive anchor modulus, the remaining coefficient units of both parents are unrestricted under the addition of boundedly many sufficiently low unit-gap rows. This holds with the complete normalized valuation arrays and the upper exponent envelopes fixed. The linked tests therefore need information about the simultaneous distribution across native moduli.

## Unit interpolation in prescribed height intervals

Fix consecutive orders \(n,N=n+1\), with \(n\ge2\), and let
\[
 \mathcal R=\{k/j:j\in\{n,N\},\ 1\le k\le j\},\qquad
 s=|\mathcal R|=n+N-1 .
\]
Let \(1/2<\theta\le1\), \(c>0\), and fix an integer \(K\) with \(K(\theta-1/2)>s\). For every sufficiently large prime \(p\), prescribe two arithmetic progressions \(I_+,I_-\) of positive integer bases, each of length between \(cp^\theta\) and \(p\) and with step prime to \(p\). For any original ordinary row profile and any collection of unit multipliers \(Z_{j,k}\in\mathbb F_p^\times\), \(j=n,N\), \(1\le k\le j\), one may add \(K\) rows \((A,A+1)\) with \(A\in I_+\) and \(K\) rows \((A+1,A)\) with \(A\in I_-\) so that

- Every added factorial ramp factor in the two orders is a unit at \(p\). After stripping the smaller-end factorials, the valuation of every coefficient of both original parents is unchanged.
- The unit part of the \(k\)-th successive coefficient ratio is multiplied by \(Z_{j,k}\) modulo \(p\).
- The sum of the oriented gaps of the added rows is zero.

The unit part of a nonzero rational number \(u\) is \(p^{-v_p(u)}u\). The theorem thus applies to coefficients and ratios with arbitrary original valuations. The \(2K\) new rows may all have bases below any prescribed height cap containing the two progressions. If both progressions have bases tending to infinity, the new rows have total \(\sum d_i^2/\max(a_i,b_i)=o(1)\) and maximum relative gap \(o(1)\).

**Proof.** An increasing unit row changes the successive ratio by \(jA+k\), and a decreasing unit row changes it by \(1/(jA+j-k+1)\). Equal counts cancel the \(j\)-factors when these expressions are written using fractions in \(\mathcal R\). Introduce independent left and right copies of \(\mathcal R\), with an edge from \(k/j\) on the left to \((j-k+1)/j\) on the right for each desired ratio. The only fraction appearing in both orders in a given copy is \(1\). The two edges at the left \(1\) end at separate vertices of right degree one, and the analogous assertion holds at the right \(1\). These two stars are disjoint because \(j\ge2\). All other edges are isolated. Thus one can assign units \(L_r,R_r\) to the vertices so that \(L_{k/j}/R_{(j-k+1)/j}=Z_{j,k}\) on every edge.

It suffices to show that for either progression \(I\), any vector \(L\in G=(\mathbb F_p^\times)^{\mathcal R}\) is the product of \(K\) vectors
\[
                       f(A)=(A+r)_{r\in\mathcal R},
                       \qquad A\in I,\quad f(A)\in G .           \tag{1}
\]
For \(p>nN\) the distinct fractions remain distinct in \(\mathbb F_p\). Remove their at most \(s\) forbidden roots from \(I\), obtaining a set \(D\) of distinct residues. For a nontrivial character \(\chi=(\chi_r)_r\) of \(G\), the complete additive twists of \(A\mapsto\prod_r\chi_r(A+r)\) restricted to all allowed residues have magnitude at most \(s\sqrt p+s\). In detail, write all the \(\chi_r\) as powers \(e_r\), \(0\le e_r<p-1\), of one primitive character. Put \(g=\gcd(p-1,(e_r)_r)\), use its \(g\)-th power of order \((p-1)/g>1\), and use the polynomial \(\prod_r(X+r)^{e_r/g}\). The gcd of its nonzero root multiplicities together with \((p-1)/g\) is one. The mixed Weil bound for this polynomial and a linear additive twist depends on at most \(s\) distinct roots. Omitting a root belonging to a trivial component costs at most one. The irreducibility hypotheses of Lemmas 1 and 2 of Mauduit and Sárközy are precisely the multiplicity gcd condition just verified and the irreducibility of \(Y^p-Y-X\), which is linear in \(X\). [Mauduit and Sárközy, *On finite pseudorandom binary sequences I*, Lemmas 1 and 2](https://matwbn.icm.edu.pl/ksiazki/aa/aa82/aa8244.pdf)

The sum of the absolute additive Fourier coefficients of a progression of length at most \(p\), with coefficients normalized by \(p^{-1}\), is \(O(\log p)\). The geometric-series bound at nonzero frequency is \(O(p/\min(t,p-t))\) before normalization for the interval of step one, and a step prime to \(p\) permutes the frequencies. Completing the sum over \(I\) therefore gives
\[
                \left|\sum_{A\in D}\chi(f(A))\right|
                   =O_s(\sqrt p\log p).
\]
By multiplicative character orthogonality the number of labelled \(K\)-tuples representing \(L\) in (1) is
\[
                  \frac{|D|^K}{(p-1)^s}
                       +O_s\big((C_s\sqrt p\log p)^K\big).
\]
The constant \(C_s\) can be chosen independent of the progressions, \(\chi\), and \(L\). The first term is \(\gg_{s,K,c}p^{K\theta-s}\), and strictly dominates the error by the choice of \(K\). This proves representation independently for \(L\) and \(R\). Equation (1) also excludes every \(p\)-divisible factor \(jA+k\) in the new ramps. Multiplying the successive ratios yields the stated coefficient claim after normalization by the constant coefficient. All statements apply equally when the common variable is reversed, by interchanging the two row orientations.

## An entire odd-corner anchor modulus can pass its residue tests

Use the notation and orientation convention of [linked quadratic specializations](MULTIVARIATE_ODD_CORNER_LINKED_QUADRATIC_ANCHORS.md). Fix \(m=2q\ge6\), \(n=m+1\), a fixed finite ordinary core, and \(h\in[1,q]\) with \(\gcd(m,h)=1\). There are arbitrarily large integers \(B\) for which \(p=mB+h\) is prime. Start with the core, an increasing row \((B-1,B)\), and a decreasing row \((B-1+q,B-1)\). Add any desired fixed-gap background below \(B/2\), provided its largest upper exponent is below \(p/n\). Take the largest background gap \(o(B)\), so for sufficiently large \(B\) the corner inequalities hold. These hypotheses include backgrounds of unbounded row count.

After stripping the smaller factorials and using the reversed index, the even parent has first three coefficients \(p\)-units and all subsequent coefficients divisible by \(p\). In the decreasing top row its first hit is \(mB+h\), reached at index three because \(m+h\le3q\). The increasing top row never reaches \(p\). All core and background row arguments are below \(p\). At the odd order the top two baselines exceed \(p\), and their complete ramps lie strictly between \(p\) and \(2p\) for sufficiently large \(B\). Thus every coefficient of the stripped odd polynomial is a unit. The prime \(p\) is the entire primitive modulus \((mB+h)/\gcd(m,h)\) at this anchor.

Apply unit interpolation with both base intervals inside \([B/4,B/3]\), \(\theta=1\) and \(K=2s+1\). All added bases and upper exponents are below the top rows. Their count is bounded in terms of \(m\), their absolute curvature is \(O_m(1/B)\), and their net oriented gap is zero. Normalize each parent by its own constant coefficient, a \(p\)-unit in the stripped form. Put \(H(z)=(1+z)(1+2z)\). The arbitrary ratio multipliers allow both normalizations to have prescribed reductions
\[
 \bar P=H,\qquad
 \begin{cases}
 \bar Q=H(1+z)^{m-1}, & \text{zero remainder},\\
 \bar Q=(1+3z)^{m+1}, & \text{coprime reductions}.
 \end{cases}                                                   \tag{2}
\]
Every coefficient in each displayed target of its stated degree is a unit for \(p\) sufficiently large. For the even parent one chooses the first two ratio multipliers to assign its surviving coefficients and may prescribe arbitrary units for the other ratios. Each construction preserves the full stripped valuation arrays of the same base profile, including at indices not surviving reduction, and preserves the upper envelopes. In this primary construction every added factorial argument is below \(p\), so the valuations of the raw coefficients are unchanged as well. The marked and unmarked leading coefficients are units in both (2), and the roots of \(H\) are distinct at \(p>3\). Hence both scalar remainders of the linked marked pair vanish at the entire modulus in the first construction. In the second construction the two polynomials have no common root at the same prime. The native rational specialization at \(X=-h/m\) gives exactly the physical coefficients modulo \(p\), as its denominator is a unit and \(mB+h=0\bmod p\).

The background can realize any prescribed limiting curvature \(\tau>0\) and signed bias \(|\lambda|\le\tau\), with any prescribed total signed gap for all rows added to the core. Here the signed bias is \(\sum_i(b_i-a_i)/\max(a_i,b_i)\). To see this, take \(A_1=\lfloor\sqrt B\rfloor\). Use \(\lfloor|\lambda|(A_1+1)\rfloor\) unit rows in the sign of \(\lambda\), and \(\lfloor(\tau-|\lambda|)(A_1+1)/2\rfloor\) pairs of opposite unit rows, all with smaller exponent \(A_1\). Correct the signed gap of the top and filler rows to its prescribed value at a base \(A_2=\lfloor B/3\rfloor\) using \(O(\sqrt B)\) unit rows. The correction, top pair, and unit calibration have total curvature and signed bias tending to zero, and all added rows have maximum relative gap tending to zero. This verifies (2) in arbitrary weak-tail limits with these parameters, including positive prescribed curvature.

## Simultaneous selected native factors

One can also prescribe the reductions at finitely many other linked anchors at one selected large prime per other anchor. Fix a primary anchor \((m,h)\) as above and finitely many other distinct anchors \((j_i,h_i)\) of the same pair of orders. Choose distinct sufficiently large primes \(\ell_i\), and impose \(j_iB+h_i\equiv0\pmod{\ell_i}\). They may be taken larger than the bounds of the fixed core, \(m(m+1)\), and all the cross determinants \(|j_i h-mh_i|\). Dirichlet's theorem then gives arbitrarily large primes \(p=mB+h\) under all these congruences, since \(p\bmod\ell_i\) is nonzero by the cross-determinant condition, and \(h\) is coprime to \(m\).

At \(\ell_i\) the reduction is marked exactly through degree two and unmarked in full degree for the core and top rows. This follows by the rational-specialization factor test in the linked-anchor note. The nonzero top-row factors are units once \(\ell_i\) exceeds the finitely many nonzero numerators there. Filler and repair bases may be chosen congruent to zero modulo \(\prod_i\ell_i\) at distance \(O(\prod_i\ell_i)\) from the stated heights, making their normalized ramps units at the selected primes. Thus unit interpolation at each fixed \(\ell_i\) supplies prescribed bases modulo \(\ell_i\) for each of the \(2K\) labelled calibration rows, with any of the compatible or incompatible residual targets in (2) after using its own marked order. For example, for an odd marked anchor use the same \(H\) and an unmarked degree-\(m\) target \(H(1+z)^{m-2}\) or \((1+3z)^m\). Take \(\ell_i\) beyond any prime divisor of the finitely many target coefficients. Combine the prescribed base classes by CRT. At the main prime \(p\), each resulting class intersects \([B/4,B/3]\) in an arithmetic progression of length \(\gg_{m,(\ell_i)_i}p\). The proof of (1) works separately for each labelled row with its own progression. Its nontrivial term is bounded by the product of the same \(K\) individual estimates, and its main term is the product of the corresponding lengths divided by \((p-1)^s\). This finishes the simultaneous interpolation. Taking a diagonal sequence allows the auxiliary primes to grow too.

These examples control one entire native modulus and selected prime factors of other moduli. For an unrestricted original profile a decisive incompatibility may occur at a different prime or at higher precision. The constructions make no claim of a common factor of the two rational polynomials.

## Unit interpolation at prime powers with a uniform row count

Fix adjacent orders \(n,N=n+1\), \(n\ge2\), and use \(\mathcal R\) and \(s=n+N-1\) from the [prime interpolation theorem](#unit-interpolation-in-prescribed-height-intervals). Fix \(0<c<1\) and an integer \(R>1/c\). For every sufficiently large prime \(p\), the following holds at every \(e\ge1\). Take two positive-base arithmetic progressions of step prime to \(p\) and of length between \(cp^e\) and \(p^e\). Any units modulo \(p^e\) prescribed as multipliers of the unit parts of all successive coefficient ratios can be realized by adding \(2s+1+Rs\) unit-gap rows of each orientation, within the respective progressions. All the new normalized ramps consist of \(p\)-units and the added rows have net oriented gap zero.

**Proof.** The bipartite forest of sample fractions is unchanged. Assign left and right labels in \((\mathbb Z/p^e)^\times\) by propagating units along the edges, as in the prime case. For either side, it remains to construct a product of sample vectors \((A+r)_{r\in\mathcal R}\) equal to an arbitrary unit vector \(L\) over \(\mathbb Z/p^e\). Denominators of all sample fractions are invertible.

Choose \(s\) distinct allowed residues \(\alpha_1,\dots,\alpha_s\bmod p\) and reserve \(R\) rows for each. Fix the \(2s+1\) other bases in the given progression so that their sample product times \(\prod_i(\alpha_i+r)^R\), coordinate by coordinate, represents \(L\bmod p\). This is the prime interpolation theorem. If the progression is longer than \(p\), use its first \(p\) entries when finding these fixed bases. For \(e=1\), the \(R\) reserved bases for each group can be chosen equal, and the result follows.

For \(e\ge2\), divide out the product of the fixed bases. Write the remaining target as \(U\). The \(s\times s\) matrix
\[
                 C_{r,i}=\frac{1}{\alpha_i+r}\quad(\bmod p)
\]
is invertible. Indeed, if a linear combination \(\sum_i c_i/(x+\alpha_i)\) vanishes at all \(s\) distinct values \(x=r\), its numerator of degree at most \(s-1\) vanishes identically. The poles are distinct, so each \(c_i=0\). The relative Jacobian of \(G_r(A_1,\dots,A_s)=\prod_i(A_i+r)^R\) at the residue tuple is \(RC\), also invertible when \(p>R\). At step \(t\ge1\), replacing \(A_i\) by \(A_i+p^t\delta_i\) multiplies the \(r\)-th component by \(1+p^t R\sum_i C_{r,i}\delta_i\) modulo \(p^{t+1}\). Recursive correction therefore produces integers \(\beta_i\equiv\alpha_i\bmod p\) satisfying \(G(\boldsymbol\beta)=U\bmod p^{e-1}\).

For the last digit, take the \(R\) actual bases of group \(i\) separately from \(A_{i,\nu}=\beta_i+p^{e-1} t_{i,\nu}\) in the given progression. The allowed \(t_{i,\nu}\bmod p\) form an arithmetic progression of length at least \(\lfloor cp\rfloor\), with step a unit, because the original step is prime to \(p\) and the original progression has length at least \(cp^e\). The sum of \(R\) choices fills \(\mathbb F_p\) for sufficiently large \(p\). In fact the sum before taking residues is an arithmetic progression of the same step and length at least \(R(\lfloor cp\rfloor-1)+1\ge p\).

Collecting the last digits multiplies \(G_r(\boldsymbol\beta)\) by
\[
           1+p^{e-1}\sum_i C_{r,i}\sum_{\nu=1}^R t_{i,\nu}
                         \pmod {p^e}.
\]
The matrix \(C\) gives the required sums uniquely, and the preceding progression argument realizes them. Thus the full unit vector is realizable modulo \(p^e\). All chosen bases avoided the forbidden roots modulo \(p\). Using the chosen number of rows on both sides cancels the common order factors in the unit ratios and completes the proof.

### Corner moduli of unbounded depth

Fix even \(m=2q\ge6\), odd \(m+1\), and any fixed finite ordinary core. Choose a sufficiently large fixed prime \(p\equiv1\pmod m\). For \(e\to\infty\) write \(B=(p^e-1)/m\), and include the two top rows \((B-1,B)\) and \((B-1+q,B-1)\) as in the [linked quadratic anchor theorem](MULTIVARIATE_ODD_CORNER_LINKED_QUADRATIC_ANCHORS.md). At the even anchor \(h=1\), the entire primitive modulus is \(mB+1=p^e\).

The even decreasing top row ramp has the literal factor \(p^e\) first at reversed index three. It is the only multiple of \(p\) in its full ramp if \(p>m q\). The even increasing top row ramp lies just below this multiple and its length is less than \(p\). All even coefficients at indices \(0,1,2\) are units and all higher even coefficients vanish modulo \(p^e\). The odd top row ramps have no multiples of \(p\) when \(p\) also exceeds the finitely many nonzero numerators of their linear factors at the specialization \(X=-1/m\). The fixed core can be made a unit contributor by taking \(p\) above all its factorial arguments at the two orders. For additional filler or signed-gap repair rows, choose unit-gap bases congruent to zero modulo \(p\). This guarantees that all their stripped factors at the two orders are units. Moving the filler bases at heights \(\sqrt B\) and \(B/3\) by at most \(p\) preserves the estimates of the [prescribed diffuse-strength construction](#an-entire-odd-corner-anchor-modulus-can-pass-its-residue-tests).

Use the lifting theorem with bases in \([B/4,B/3]\), a fixed positive proportion of \(p^e\) for all large \(e\). After each parent's constant coefficient is normalized to one, prescribe the even reduction \(H=(1+z)(1+2z)\) and either \(H(1+z)^{m-1}\) or \((1+3z)^{m+1}\) at the odd order, now modulo \(p^e\). The prescribed nonzero coefficients are units at \(p\) when \(p\) exceeds their fixed prime divisors. The same full stripped valuation arrays are retained in both constructions. The added calibration rows contribute \(O_{m,c}(1/B)\) curvature, and all upper factorial envelopes remain fixed.

In the first case both scalar remainders vanish modulo the entire \(p^e\). In the second case the reductions are already coprime modulo \(p\). Thus deeper precision at a single complete native anchor supplies no profile-uniform incompatibility, even with its prime fixed and depth arbitrarily large. The other linked anchors impose further arithmetic conditions on the same original profile.

## Exact checks

The [finite checker](verify_multivariate_height_capped_native_unit_interpolation.py) constructs both outcomes at the full native moduli \(2437,3617,4201,2437^2,2437^3\), for even orders \(6,8,10,6,6\), respectively. It uses discrete logarithms and CRT for the first residues and the displayed lift at higher precision, then recomputes the ordinary stripped ramp valuations and residues directly. The checks cover both full valuation arrays, each target coefficient, and evaluation at both simple quadratic roots. The uniform bound on the number of added rows uses the character-sum argument for the first residues.
