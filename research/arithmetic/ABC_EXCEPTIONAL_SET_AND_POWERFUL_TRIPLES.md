# The abc exceptional set and powerful triples

For \(\lambda>0\), let \(N_\lambda(X)\) count ordered positive primitive triples \(a+b=c\le X\) with \(\operatorname{rad}(abc)<c^\lambda\). The \(abc\) conjecture is the assertion that \(N_\lambda(X)\) is bounded in \(X\) for every \(0<\lambda<1\). Let \(T_{p,q,r}(X)\) count ordered positive primitive triples \(a+b=c\le X\) for which the three terms are respectively \(p\)-full, \(q\)-full and \(r\)-full. An integer is \(s\)-full when every positive prime valuation is at least \(s\), so \(1\) is \(s\)-full for every \(s\). [Bernert–Browning–Lichtman–Teräväinen, §1](https://arxiv.org/html/2410.12234v2#S1) · [Browning–Verzobio, §1](https://arxiv.org/html/2608.24512v1#S1)

## Radical and powerfulness scopes

For an \(s\)-full integer \(n\), prime by prime, \(\operatorname{rad}(n)\le n^{1/s}\). Consequently, if \(\lambda_0=1/p+1/q+1/r\), then
\[
\operatorname{rad}(abc)\le a^{1/p}b^{1/q}c^{1/r}<c^{\lambda_0}.
\]
Thus \(T_{p,q,r}(X)\le N_{\lambda_0}(X)\). The \(abc\) conjecture would give finiteness whenever \(\lambda_0<1\). Square-full triples have \(\lambda_0=3/2\). For example, \(9+16=25\) is primitive with radical \(30\), so it is outside every \(N_\lambda\) with \(\lambda\le 1\). Conversely, \(2+3^{10}\cdot109=23^5\) is primitive with radical \(15042<23^4\). Its first two terms have primes of valuation one, even with the fixed strict cutoff \(\lambda=4/5\).

## Exceptional-set estimates and dyadic coordinates

The de Bruijn radical count and fixing two letters gives \(N_\lambda(X)\ll_\varepsilon X^{2\lambda/3+\varepsilon}\) for \(\lambda>0\). The three-author paper states \(X^{33/50}\) for fixed \(0<\lambda<1.001\). Runbo Li's optimization of that version states \(X^{56/85+\varepsilon}\) for \(0<\lambda<1+\delta(\varepsilon)\) with \(\delta(\varepsilon)>0\). Bernert's simpler argument gave \(X^{(23\lambda+3)/40+\varepsilon}\) for \(0<\lambda<2\). The four-author revision proves that last bound for \(0<\lambda\le2\) and gives \(X^{3/5+\varepsilon}\) for fixed \(0<\lambda<1\). Kane determines the power order \(X^{\lambda-1+o(1)}\) in the separate range \(2\le\lambda\le3\). [BLT v1, Theorem 1.2](https://arxiv.org/html/2410.12234v1) · [Runbo Li, Theorem 1.3](https://arxiv.org/abs/2507.02885) · [Bernert, Theorem 1](https://arxiv.org/html/2506.13364v1) · [BBLT v2, Proposition 1.1 and Theorems 1.2–1.3](https://arxiv.org/html/2410.12234v2#S1) · [Kane, Theorems 1–2](https://cseweb.ucsd.edu/~dakane/ABC2.pdf)

BBLT Proposition 2.1 groups primes by valuation and moves large valuations into bounded-index factors and small coefficients. Its boxes have terms \(s_0\prod_{i=1}^{d}s_i^i\), with \(d\) fixed after specifying the exponent loss. Their three coordinate products encode the radical budget up to that loss. Proposition 2.3 fixes two complete letters, Proposition 3.1 uses additive energy of grouped powers, and Proposition 4.1 freezes coordinates and counts primitive solutions to a linear equation. Proposition 6.1 combines selected instances of these bounds on the logarithmic box scales. The [omitted-coordinate graph bound](ABC_OMITTED_COORDINATE_GRAPH.md) supplies another inequality for boxes with exponent-one factors. [BBLT v2, §§2–4 and 6](https://arxiv.org/html/2410.12234v2)

November preprints use collision counts modulo parts of the third letter. Hongmin Li's original archived collision paper stated \(59969/100033\) at \(\lambda=1\), followed by \(54/91\) in the radical-budget paper and \(589903/1000008\) in the squares and twisted-energy paper. The first paper's aiXiv record now points to a withdrawal. These are version-specific claims. [Li, first archive](https://zenodo.org/records/24054892) · [aiXiv withdrawal record](https://aixiv.online/abs/2611.02367) · [Li, radical budget, v2](https://rstation.ai/problems/math/papers/17/v2) · [Li, squares and twisted energy, v3](https://rstation.ai/problems/math/papers/28/v3)

Chunping Zeng's Theorem 1.1 claims \(U=11680007/20000012<0.584\) for \(0<\lambda\le1\), using grouped power differences, 32 third-letter modulus choices and an exact rational tree. The program witness at \(D=7/12\) and the strict separation in §8 compare exponent programs. The archived rational tree and the arithmetic argument have separate review obligations. The 20 November archive proposes the stricter \(5845/10009\) and adds a pair-deviation split. The [certificate audit](ABC_PAIR_DEVIATION_CERTIFICATE_AUDIT.md) checks that split against the exact ancestor rows. [Zeng, Theorem 1.1 and §§3–8](https://hexagonmath.org/pdf/2611.00621v1) · [19 November certificate archive](https://doi.org/10.5281/zenodo.24354214) · [20 November update](https://doi.org/10.5281/zenodo.24407593)

Pointwise bounds for the size of \(c\) are a separate route to finiteness. With \(R=\operatorname{rad}(abc)\), Stewart–Yu give \(\log c\ll R^{1/3}(\log R)^3\). Pasten's Theorem 1.4 gives \(\log c\le \eta^{-1}\exp(C\sqrt{\log R\log\log R})\) when \(a\le c^{1-\eta}\), and the same expression with \(\min\{P^+(a),P^+(b),P^+(c)\}\) in place of \(\eta^{-1}\) in the second part. Soulé's November preprint claims \(c\ll_\varepsilon R^{24+\varepsilon}\) in Theorem B. If its proof is accepted, that bound gives finiteness for fixed \(\lambda<1/24\) by choosing \(\varepsilon\) so that \(\lambda(24+\varepsilon)<1\). [Stewart–Yu, Theorem 1](https://uwaterloo.ca/pure-mathematics/sites/default/files/uploads/documents/s0012-7094-01-10815-6.pdf) · [Pasten, Theorem 1.4](https://arxiv.org/html/2312.03566v1) · [Soulé, Theorem B](https://arxiv.org/html/2611.05217v1)

At the endpoint \(\lambda=1\), Bright gives infinitely many triples with \(c/R\) growing at least as \(\exp(C\sqrt{\log c}/\log\log c)\) for an explicit positive constant \(C\). For any fixed \(0<\lambda<1\) this lower bound on \(c/R\) has a subpower scale in \(c\). [Bright, Theorem 3.1](https://cs.curtisbright.com/reports/abc-lower.pdf)

Carella's v3 preprint claims infinitely many strict exceptions, but the identity used in equation (4.35) already fails at a prime square. For the weight \(z^{\omega(n)}\), the Euler factor at \(p\le y\) is \(1+zp^{-s}/(1-p^{-s})\), whose coefficient at \(p^{-2s}\) is \(z\). The displayed factor \((1-p^{-s})^{-z}\) has coefficient \(z(z+1)/2\). They differ for \(z\ne0,1\), including the varying values near \(1\) used in the Perron argument. The claimed short-interval estimate in (4.40) requires a proof for the correct generating function. [Carella v3, §4.3 and equation (4.35)](https://arxiv.org/pdf/2608.16764v3)

## Fourier bound for powerful triples

For an integer \(s\ge2\) define
\[
\Phi(s)=\max_{1\le k\le s}\frac{4k-1}{3k(2s+k-1)}.
\]
For all fixed \(p,q,r\ge2\) and \(\varepsilon>0\),
\[
 T_{p,q,r}(X)\ll_{p,q,r,\varepsilon} X^{\Phi(p)+\Phi(q)+\Phi(r)+\varepsilon}.
 \tag{1}
\]
The same upper bound counts all positive triples with these powerfulness conditions, with or without primitivity. The proof combines the Fourier moments in BBLT Proposition 3.1 with the ordered-coordinate optimization in Browning–Verzobio Lemma 4.3. [BBLT v2, Proposition 3.1](https://arxiv.org/html/2410.12234v2#S3) · [Browning–Verzobio, §4](https://arxiv.org/html/2608.24512v1#S4)

### Fourier counting in a box

For each \(s\)-full term, write
\[
 n=\prod_{i=s}^{2s-1} x_i^i,
\]
with all \(x_i\) for \(i>s\) squarefree and pairwise coprime. Indeed, a valuation \(v\ge s\) divisible by \(s\) is absorbed into \(x_s\). For \(v=qs+j\), \(0<j<s\), put the prime into \(x_{s+j}\) and its \((q-1)\)-st power into \(x_s\). This specifies a unique tuple. Separate the coordinates into dyadic intervals \([P_i,2P_i)\) with \(P_i\) integral powers of two starting at \(1\). In a box that contains an admissible term, \(\prod_i P_i^i\le X\). Dropping the squarefreeness, coprimality and height restrictions within the box enlarges the count. The monomials in the enlarged box are \(O_s(X)\).

For one letter set \(P=\prod_i P_i\), \(M=\max_i P_i\), and \(S(t)=\sum_{\mathbf x} e(t\prod_i x_i^i)\) over the enlarged box. The box form of BBLT's moment argument gives
\[
 \int_0^1|S(t)|^2\,dt\ll X^\varepsilon P,\qquad
 \int_0^1|S(t)|^4\,dt\ll X^\varepsilon P^3/M.
 \tag{2}
\]
Here each integer of size \(O_s(X)\) has at most \(O_{s,\varepsilon}(X^\varepsilon)\) tuples giving that value by the divisor bound, which proves the second moment. For the fourth moment, choose \(j\) with \(P_j=M\) and write the monomial as \(D x_j^j\). Cauchy–Schwarz over the other coordinates, followed by integration against \(|S(t)|^2\), bounds the moment by \(P/M\) times the number of solutions of \(D(x_j^j-y_j^j)=V-W\), where \(V,W\) range over full coordinate tuples and both \(x_j,y_j\) range over the \(j\)-th interval. If \(V=W\) the power difference is zero, giving \(O_{s,\varepsilon}(P^2X^\varepsilon)\) solutions. For distinct \(V,W\), each choice of the two full tuples allows at most divisor-many choices of \(D\) and its coordinates. A nonzero difference of two \(j\)-th powers, \(j\ge2\), fixes their positive gap up to a divisor and determines the smaller entry by strict monotonicity. This gives the same bound for the nonzero part and proves (2). The signs of the difference cause no change, since the two positive entries can be swapped.

Cauchy–Schwarz gives \(\int |S|^3\ll X^\varepsilon P^2/\sqrt{M}\). Apply Hölder with exponent three to the three letters and use additive orthogonality for \(a+b=c\). With \(\alpha_i=\log_X P_i\) on each letter, its contribution to the logarithmic exponent is
\[
 \frac23\sum_i\alpha_i-\frac16\max_i\alpha_i,
 \qquad \sum_i i\alpha_i\le1 .
 \tag{3}
\]
Constants and divisor losses in (2) can be taken with a smaller exponent and absorbed into any prescribed error in (1).

### Exact optimization and fixed offsets

For \(t\ge0\), let \(k(t)=\#\{i:\alpha_i>t\}\). The expression in (3) is the integral of \(2k(t)/3-1/6\) on the set \(k(t)>0\). The weight at these indices is at least \(k(t)(2s+k(t)-1)/2\), the sum of the \(k(t)\) least allowed indices. Bounding the ratio of the two integrands and integrating proves the upper bound \(\Phi(s)\). Equality for the continuous box-scale optimization is realized by putting equal positive mass at the first \(k\) indices, normalized to make the weighted sum one. This is also Browning–Verzobio's Lemma 4.3 with \(\theta=3/4\), multiplied by \(2/3\). Summing the box bounds costs \(O_{p,q,r}((1+\log X)^{p+q+r})\), absorbed into \(X^\varepsilon\). [Browning–Verzobio, Lemma 4.3](https://arxiv.org/html/2608.24512v1#S4)

For real \(k>0\) the expression defining \(\Phi(s)\) increases up to \(k=(1+\sqrt{8s-3})/4\) and decreases afterwards. The maximizing integer therefore satisfies \(k=\sqrt{s/2}+O(1)\) as \(s\) grows. Substituting this into (3)'s exact optimum and expanding gives
\[
 \Phi(s)=\frac{2}{3s}-\frac{\sqrt 2}{3s^{3/2}}+\frac{1}{2s^2}+O(s^{-5/2}).
\]
In particular for fixed \(u\ge v\ge0\), \(p=r+u,\ q=r+v\), (1) gives
\[
 T_{p,q,r}(X)\ll_{u,v,r,\varepsilon}
 X^{\,1/p+1/q-\sqrt{2}\,r^{-3/2}+(3/2+(u+v)/3)r^{-2}
       +O_{u,v}(r^{-5/2})+\varepsilon}.
\]
The saving from \(1/p+1/q\) is asymptotically \(\sqrt{2}\,r^{-3/2}\). It improves the leading \(r^{-2}\) saving in Browning–Verzobio's Theorem 1.1 and the leading \(2r^{-2}\) saving in the [powerful-number graph manuscript](../../unpublished/powerful-number-representations/main.tex). The count here is over all targets \(c\le X\). The fixed-target graph theorem controls a different quantity. [Browning–Verzobio, Theorem 1.1](https://arxiv.org/html/2608.24512v1#S1)

## Stronger saving from the grouped moments

The full divisibility grouping improves (1) at large orders. For fixed \(u\ge v\ge0\) and \(p=r+u,\ q=r+v\), it gives the exponent
\[
 \frac1p+\frac1q-\left(\frac1{\sqrt{2}}-o_{u,v}(1)\right)
                  \frac{\sqrt{\log r}}{r^{3/2}}+\varepsilon
 \tag{4}
\]
as \(r\to\infty\), with the same height count and fixed-order convention as in (1).

For a letter of order \(s\), group every coordinate with index divisible by a chosen \(g\ge2\). The argument for (2) applies to their combined \(g\)-th power. The number of coordinate tuples giving the same combined root is bounded by the divisor function with a fixed number of factors. Thus we may replace the maximum in (3) by \(M=\max_{g\ge2}\sum_{g\mid j}\alpha_j\). Put
\[
 S=\sum_{j=s}^{2s-1}\alpha_j,\quad
 V=\sum_{j=s}^{2s-1}(j-s)\alpha_j,\quad
 \Delta=\frac{2}{3s}-\frac23 S+\frac16 M,\quad
 L_s=\frac{\sqrt{\log s}}{s^{3/2}}.
\]
The height constraint gives \(sS+V\le1\), hence \(\Delta\ge2V/(3s)+M/6\). It suffices to consider \(\Delta\le L_s\). Then \(M\le6L_s\) and \(1/s-(3/2)L_s\le S\le1/s\). Every composite \(j\le 2s-1\) has a prime factor at most \(\sqrt{2s}\). The composite indices therefore have total mass at most \(\pi(\sqrt{2s})M=O(1/(s\sqrt{\log s}))\). Thus the prime indices have total mass \(P=(1+o(1))/s\), and each has mass at most \(M\). [BBLT v2, Proposition 3.1](https://arxiv.org/html/2410.12234v2#S3)

Take \(H=\sqrt{s\log s/8}\). For the prime indices, compare the offset \(j-s\) with \(H\) and integrate the count of primes before each offset. This gives
\[
 V\ge HP-M\int_0^H \#\{j\text{ prime}:s\le j\le s+t\}\,dt.
\]
Brun–Titchmarsh bounds the integrand by \(2t/\log t+1\) for \(t>1\). The integral up to \(H/\log H\) is \(o(H^2/\log H)\) by the trivial bound \(1+t\). Above that threshold the logarithm in the denominator is \((1+o(1))\log H\), so the integral is at most \((1+o(1))H^2/\log H=(1+o(1))s/4\). Therefore
\[
 \Delta \ge \frac{2HP}{3s}-o(1)M
 \ge \left(\frac{1}{3\sqrt 2}-o(1)\right)L_s.
\]
Adding the deficits at \(s=p,q,r\) and using \(2(1/p+1/q+1/r)/3=1/p+1/q+O_{u,v}(r^{-2})\) proves (4). [Montgomery–Vaughan, equation (1.12)](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf)

## Replay of the 19 November modulus-locking certificate

The author archive SHA-256 is f34c459e5dccbfa71f2e9729afefd2c206f3202107b14f11842747de584c6538. In the extracted archive, the manifest and all three original checks passed with the commands below. [Zenodo archive](https://doi.org/10.5281/zenodo.24354214)

```sh
sha256sum -c SHA256SUMS
python3 -B check_extended.py cert_fulloc_5840.json.gz --faults
python3 -B check_extended.py cert_old_5845.json.gz --faults
python3 -B check_obstructions.py
```

The main tree returned 65,576 leaves and exact upper \(11680007/20000012\), and the old control returned 49,033 leaves and upper \(11697/20012\). Both rejected the three published fault injections. The obstruction checker returned \(403/690\) on the old witness, \(53/92\) at the same point under FullOC, and \(7/12\) on the symmetric FullOC witness. The replay checks the finite rational programs. Validation of every arithmetic input and its application to \(N_\lambda(X)\) is a separate obligation. [Zeng, §§6–8 and Appendix A](https://hexagonmath.org/pdf/2611.00621v1)

## From the exact FullOC optimum to the exceptional-set count

The [repaired FullOC certificate](ABC_PAIR_DEVIATION_CERTIFICATE_AUDIT.md) bounds its program by \(D\le7/12\). The symmetric feasible point in Zeng's equation (8.7) attains \(7/12\), so this is the optimum of the program. The upper certificate uses only the coprime exponent pairs. The lower witness satisfies the full family that also includes noncoprime pairs. [Zeng, §§5–8](https://hexagonmath.org/pdf/2611.00621v1)

In particular, for \(0<\lambda\le1\) and each \(\varepsilon>0\),
\[
 N_\lambda(X)\ll_\varepsilon X^{7/12+\varepsilon}.
 \tag{5}
\]
It suffices to work with \(Y/2\le c\le Y\) and \(\operatorname{rad}(abc)<Y\). Apply BBLT Lemma 2.2 to each term with parameter \(\eta>0\). The three exceptional coefficients cost at most \(Y^{3\eta}\) choices and the dyadic coordinate scales cost a fixed power of \(\log Y\) at each \(\eta\). By pairwise coprimality and the lemma's radical comparison, the product of the coordinate scales is at most \(Y^{1+3\eta}\). The weighted scale of each letter is at most \(Y\), and that of the largest letter is at least a constant times \(Y^{1-\eta}\) for fixed \(\eta\), since \(c\ge Y/2\). This keeps the radical budget for the limiting program and the height denominator in the geometry estimate. [BBLT v2, Lemma 2.2](https://arxiv.org/html/2410.12234v2#S2)

If (5) failed, there would be a positive gap above \(7/12\) in the count exponent of arbitrarily large height boxes. Choose a sequence \(\eta\to0\) and heights beyond all constants at each fixed \(\eta\), retaining a positive gap after the choice of coefficients and dyadic boxes. For each box retain the first six logarithmic coordinates per letter, their total over all indices, and the box count exponent \(D\). The three height bounds and the pair counts bound these parameters independently of the truncation length. Take a convergent subsequence and sort the letters by their degree-four mass plus twice their degree-five mass, with signs in the equations carried through. The pair, Fourier and geometry estimates and Zeng's signed collision inequalities give the limiting program rows. The weighted height tail is bounded from below by seven times the tail mass and the geometry denominator comes from the largest complete value, regardless of its position in the sorting. The row errors vanish in the limit, and the alternatives in each row are closed because they are finite maxima, minima and positive parts. This yields a feasible program point with \(D>7/12\), a contradiction. Summing dyadic heights with an arbitrarily small exponent loss proves (5). [Zeng, §§2–7](https://hexagonmath.org/pdf/2611.00621v1)

## The limit of the six-index projection

For every integer \(r\ge7\), at the budget \(\lambda=3/r\), put all of each letter's \(1/r\) mass in degree \(r\). In the BBLT six-index projection this gives side totals \(1/r\), individual coordinates \(1,\ldots,6\) all zero, and the point \(D=2/r\). The three height inequalities use \(7/r\le1\), and the pair bounds hold with equality. Every projected geometry saving and every displayed Fourier maximum is zero. Thus Proposition 6.1's projected inequalities alone allow the full trivial exponent \(2/r\) at this budget. The Fourier estimate with the degree \(r\) retained gives a strictly smaller exponent. An optimization intended to recover this powerful-triple saving must retain a bound on the higher coordinates. [BBLT v2, inequalities (6.1)–(6.10)](https://arxiv.org/html/2410.12234v2#S6)

## Exceptional-set improvement from grouped roots

The [grouped-root certificate](ABC_PAIR_DEVIATION_CERTIFICATE_AUDIT.md#grouped-root-certificate)
gives the stronger bound
\[
 N_\lambda(X)\ll_\varepsilon X^{1747/3000+\varepsilon}
 \qquad (0<\lambda\le1).
 \tag{6}
\]
The sparse-root sets in the collision estimate separate the number of root
values from their numerical heights. For instance, in the six-index
symmetric profile, the grouped square root \(a_2a_4^2\) has size exponent
\(90/396\) and cardinality exponent at most \(54/396\). The separate bounds
give an exponent \(37/66\) on that profile. The finite certificate uses these
inequalities on the full remaining set of exponent profiles. The compactness
argument above applies unchanged to the enlarged finite list of arithmetic
bounds. The tail enters the new inequalities through the original lower
bound for the frozen modulus, and all additional maximum terms are
continuous. Thus a height-averaged count above the exponent in (6) would give
a feasible limiting profile above the certified threshold.

For a fixed \(\lambda\), one can replace \(1747/3000\) in (6) by
\(U_{\mathrm{cert}}(\lambda)\). The radical product constraint now gives
\(L\le\lambda\) in the limiting program, so the same arithmetic inequalities
admit the budget-dependent leaf bounds. [Exact envelope](ABC_PAIR_DEVIATION_CERTIFICATE_AUDIT.md#dependence-on-the-radical-budget)

### Using all higher coordinates of the largest letter

Let \(C\) denote the actual largest term, with \(C\asymp X\). Write \(r_T\)
for the sum of coordinate masses above degree six on letter \(T\). Freeze
all coordinates of \(C\) above degree six and an additional set
\(O_C\subseteq\{1,\ldots,6\}\). The frozen mass is
\(f_C=r_C+\sum_{i\in O_C}\alpha_{C,i}\). The residual product has weighted
logarithmic height at most \(\sum_{i\le6,\ i\notin O_C}i\alpha_{C,i}\),
up to the dyadic error. Since \(M_Cw=C\), the modulus exponent in the
grouped-root collision inequality can be taken to be
\[
 w_C=1-\sum_{i\le6,\ i\notin O_C}i\alpha_{C,i}
\]
in the limit. This argument uses the complete largest value, so the lower
modulus exponent is uniform in how the tail mass is distributed. On \(B\)
one may also release every coordinate above degree six along with a
chosen set of early coordinates. Its free mass then includes \(r_B\).
The proof of the collision inequality permits these choices even when
the tail indices have no common divisor. The residue group whose root is
counted can consist of early coordinates alone.
[Grouped-root collision proof](ABC_OMITTED_COORDINATE_GRAPH.md#from-energy-to-the-box-exponent)

At the symmetric \(4/7\) profile of the strengthened six-index program,
the first six masses are \((168,0,300,200,25,0)/2100\), the total per
letter is \(700/2100\), and each tail has mass \(7/2100\). Release the
first six coordinates of \(A\) with root groups \(\{3\},\{4\}\).
Release degrees \(1,3,5\) and the full tail of \(B\). Lock degree \(5\)
and the full tail on \(C\), and use the single cubic coordinate of \(B\)
as the residue group. In units of \(1/2100\), the frozen masses are
\[
 f_A=7,\quad f_B=200,\quad f_C=32,\quad w_B=800,\quad w_C=232.
\]
Here \(232=2100-(168+900+800)\) uses the largest height. Both root size
penalties vanish, so \(K=2(693-300-200)+300+200=886\).
The cubic residue cap is \(500-300+(300-232)=268\). The diagonal and
off-diagonal exponents are \(907/2100\) and \(1150/2100\), respectively.
Thus every lift of this profile to an actual dyadic box obeys
\(D\le23/42\) in the limiting exponent estimate. The full-tail choice
also excludes a neighborhood of the profile by continuity.
