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

## Analytic savings at every radical budget

Let \(0<\lambda\le1\) be fixed. For \(t\ge3\) put
\[
S(t)=\sum_{\substack{q\le t\\q\ {\rm prime}}}(t-q),
\qquad
A_\lambda(t)=\frac{(t-1)(3-\lambda)}{2}.
\]
Then
\[
N_\lambda(X)\ll_{\lambda,t,\varepsilon}
X^{B_\lambda(t)+\varepsilon},
\qquad
B_\lambda(t)=
\frac{3+\bigl(4S(t)+4A_\lambda(t)-1\bigr)\lambda}
     {6S(t)+7A_\lambda(t)}.
\tag{7}
\]
Taking \(t=6/\lambda\) gives \(B_\lambda(t)<2\lambda/3\) for every such
\(\lambda\). The published affine estimate \((23\lambda+3)/40\) beats
\(2\lambda/3\) when \(\lambda>9/11\). Formula (7) uses the Fourier and
geometry inputs of that proof with a degree cutoff that can grow as
\(\lambda\) decreases. Hongmin Li's parameter-dependent collision bounds give stronger
estimates at some larger budgets. The radical-budget manuscript
states its envelope on \([0.70,1.10]\), and the twisted-energy
manuscript recomputes budgets \(0.80,0.90,1,1.05,1.10\).
[Li, radical budget, v2](https://rstation.ai/problems/math/papers/17/v2) ·
[Li, squares and twisted energy, v3](https://rstation.ai/problems/math/papers/28/v3) [BBLT v2, Theorem 1.2 and Propositions 3.1 and 4.1](https://arxiv.org/html/2410.12234v2#S1)

There is also an absolute \(c>0\) such that for all sufficiently small
fixed \(\lambda>0\),
\[
N_\lambda(X)\ll_{\lambda,\varepsilon}
 X^{2\lambda/3-c\lambda^{3/2}\sqrt{\log(1/\lambda)}+\varepsilon}.
\tag{8}
\]
One may take any fixed \(c<1/(18\sqrt{14})\) for sufficiently small
\(\lambda\).

The extra saving comes from concentration of the occupied degrees near
one degree in any profile whose exponent approaches \(2\lambda/3\).
Brun–Titchmarsh limits the prime degrees available in that interval.
The order of the saving is the same as in the grouped Fourier bound for
powerful triples after setting \(\lambda\) of order \(1/r\).
[Montgomery–Vaughan, equation (1.12)](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf)

### Fourier and geometry constraints

In an anatomic box, let \(x_{\ell,i}\) be the logarithmic mass of the
coordinate of degree \(i\) on letter \(\ell\), and put
\[
p_i=\sum_\ell x_{\ell,i},\quad L=\sum_i p_i,\quad
F=\sum_\ell\max_{k\ge2}\sum_{k\mid i}x_{\ell,i},\quad
W=\sum_i(i-1)p_i .
\]
The total weighted height gives \(W\le3-L\). In the proof of
Proposition 3.1, before replacing the product of the three coordinate
products by the radical budget, the Fourier bound gives
\(D\le 2L/3-F/6\). The proof of Proposition 4.1 uses that same product
rather than the budget for its first term and the height of the largest
letter for its second. Thus for any selection of coordinates of total
mass \(s\) and total \((i-1)\)-weighted mass \(w\), it gives
\[
D\le\max\{L-s,\;L-1+w\}.
\tag{9}
\]
Here the constraints are written without the vanishing exponent errors.
In the finite box estimate, \(L\le\lambda+O(\eta)\) and both bounds have
\(O(\eta)\) errors after taking \(X\) large for fixed decomposition
parameter \(\eta\). For fixed \(\lambda\), one can divide the coordinate
masses by \(1+O_\lambda(\eta)\) to restore \(L\le\lambda\), and decrease
the rescaled \(D\) by \(O_\lambda(\eta)\) to enforce the displayed bounds.
The change in the constant \(1\) in (9) is also \(O_\lambda(\eta)\).
The estimates are uniform over selections at each fixed truncation.
Then let \(\eta\) tend to zero, with the box and dyadic decomposition
losses absorbed in \(\varepsilon\). [BBLT v2, Proposition 2.1 and proofs of Propositions 3.1 and 4.1](https://arxiv.org/html/2410.12234v2)

### Selecting the released coordinates

Let \(R(u)\) be the largest mass of a fractional selection of the
coordinates with total \((i-1)\)-weighted mass at most \(u\). Selecting
the coordinates in increasing order of \(i\) maximizes this mass.
Include all degree-one coordinates, whose weight is zero. At most one
positive-weight coordinate is split, and its mass is at most \(F\).
The integral prefix therefore loses at most \(F\). At every
\(0\le u<1-L+D\), the second term of (9) is strictly below \(D\), and
the integral prefix gives \(R(u)\le L-D+F\). For \(D>0\), the value
\(u_*=1-L+D\) is positive. Continuity of the fractional optimum gives
\[
R(u_*)\le L-D+F. \tag{10}
\]

Write \(P=L-p_1\) for the mass of the positive-weight coordinates. The
uniform fractional selection on these coordinates gives
\(R(u)\ge p_1+\min(1,u/W)P\) for \(W>0\), with \(R(u)=L\) if \(W=0\).
Together with (10) and \(W\le3-L\) this gives
\[
 D-F\le P\,\frac{2-D}{3-L},
 \qquad
 P\ge\frac{3-\lambda}{2}(D-F).
\tag{11}
\]
For the second inequality, use \(D\ge0\) when \(D\ge F\) and \(P\ge0\)
otherwise. The case \(u_*\ge W\) already gives \(D\le F\).
If \(D=0\) the bounds on the count exponent follow directly.

Assign each \(i\ge2\) to its least prime divisor \(q\). For each \(q\) the
assigned mass is at most \(F\) by the divisibility constraints. Since
\((t-i)_+\le(t-q)_+\), this gives
\[
 (t-1)P\le W+\sum_{i\ge2}(t-i)_+p_i
              \le 3-L+FS(t).
\]
Using \(F\le4L-6D\) and (11) yields
\[
 (6S(t)+7A_\lambda(t))D
 \le3+\bigl(4S(t)+4A_\lambda(t)-1\bigr)L.
\]
The coefficient of \(L\) is positive for \(t\ge3\), so \(L\le\lambda\)
proves (7). Its gap from \(2\lambda/3\) is
\[
 \frac{\bigl(1+2A_\lambda(t)/3\bigr)\lambda-3}
      {6S(t)+7A_\lambda(t)}.
\]
At \(t=6/\lambda\) the numerator is \(3-2\lambda+\lambda^2/3>0\).

Keeping the factor \(2-D\) in (11) gives a sharper implicit bound.
If \(D\ge4\lambda/7\), then \(D-F\ge7D-4L\ge0\), and
\[
 (t-1)(7D-4\lambda)
 \le(2-D)\left(1+\frac{(4\lambda-6D)S(t)}{3-\lambda}\right).
\tag{12}
\]
Here \((4L-6D)/(3-L)\) increases with \(L\), so one can use \(L\le\lambda\)
after substituting \(F\le4L-6D\) into (11) and the prime-divisor estimate.
On \([4\lambda/7,2\lambda/3]\) the left-hand side increases and the
right-hand side decreases. Thus their intersection bounds \(D\)
whenever it lies in the interval. At \(\lambda=1/2\) and \(t=17\),
(12) gives \(D<32311/100000\).

### Concentration near a single degree

Write \(\delta=2\lambda/3-D\) and assume
\(\delta\le\lambda^{3/2}\sqrt{\log(1/\lambda)}\), since a larger \(\delta\)
already gives (8). Then \(L=\lambda+O(\delta)\), \(F\le6\delta\), and
\(D>F\) for sufficiently small \(\lambda\). By (10), \(W>u_*\).
Set \(\theta=u_*/W\). The uniform fractional selection and (10) give
\(1/4\le\theta\le1/2\) for sufficiently small \(\lambda\).
Moreover
\[
 0\le R(u_*)-\theta L
 \le L-D+F-\frac{u_*L}{3-L}
 = F+\frac{2L-3D}{3-L}
 \ll\delta.
\]
Choose a positive threshold \(v\) equal to the weight per unit mass of
the marginal coordinate in the greedy selection. Its fractional
optimum is
\[
 R(u_*)=u_*/v+\sum_i p_i(1-(i-1)/v)_+ .
\]
Let \(A_+\) and \(A_-\) be the sums of \(p_i(1-(i-1)/v)_+\) and
\(p_i((i-1)/v-1)_+\), respectively. Their difference is \(L-W/v\).
The exact identity
\(R(u_*)-\theta L=(1-\theta)A_++\theta A_-\)
therefore gives
\[
 \sum_i p_i\lvert1-(i-1)/v\rvert\ll\delta,\qquad v\asymp 1/\lambda .
\]
For the second assertion use \(\lvert L-W/v\rvert\ll\delta\) and
\(u_*\le W\le3-L\).

Consequently an interval of half-width
\(H=C(1+\delta/\lambda^2)\) around \(v+1\), with an absolute \(C\) chosen
large enough, contains at least \(3\lambda/4\) mass by the preceding
weighted deviation bound. Its indices are comparable to \(1/\lambda\).
Every composite there has a prime divisor \(q\ll\lambda^{-1/2}\),
so the entire composite mass there is at most
\[
 F\,\pi(O(\lambda^{-1/2}))
 \ll \frac{\delta}{\sqrt{\lambda}\log(1/\lambda)}
 =o(\lambda).
\]
Here the prime counting estimate follows already from Brun–Titchmarsh
with modulus one. Each individual prime degree has mass at most \(F\).
Hence the interval contains \(\gg\lambda/\delta\) primes. Counting all
integer indices first gives \(\lambda\ll\delta H\) and
\(\delta\gg\lambda^{3/2}\), since \(\delta=o(\lambda)\). In particular
\(H\gg\lambda^{-1/2}\). Brun–Titchmarsh on the same interval now gives
\[
 \lambda\ll \frac{\delta H}{\log H}
         \ll \frac{\delta^2}{\lambda^2\log(1/\lambda)},
\]
which proves (8). [Montgomery–Vaughan, equation (1.12)](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf)

### Leading constant for the small-budget saving

Write \(\ell=\log(1/\lambda)\), \(\delta=2\lambda/3-D\), and
\(\delta_L=2L/3-D\). It suffices to consider
\(\delta\le\lambda^{3/2}\sqrt\ell\). Then \(L=\lambda+O(\delta)\)
and \(F\le6\delta_L\le6\delta\). Let \(u_*=1-L+D\) and
\(\theta=u_*/W\), as above. By \(W\le3-L\) and (10),
\[
0\le R(u_*)-\theta L
\le F+L-D-\frac{Lu_*}{3-L}
= F+\frac{3\delta_L}{3-L}.
\]
The same bounds give \(\theta=1/3+o(1)\) because
\(u_*\to1\) and \(R(u_*)\le L/3+O(\delta)\). Hence \(W\sim3\).
For the fractional threshold \(v\), the absolute deviation bound above
gives \(W/v=L+O(\delta)\), so \(v\sim3/\lambda\).

The preceding localization gives an interval of length
\(H\asymp1+\delta/\lambda^2\) near \(v+1\) containing a fixed
positive fraction of the total mass. Every degree in it has mass at
most \(F\). Counting integer degrees first gives
\(\lambda\ll FH\ll\delta(1+\delta/\lambda^2)\), so
\(H\gg\lambda^{-1/2}\). The composite degrees there have total mass
\(O(F/(\sqrt\lambda\,\ell))=o(\lambda)\). Brun–Titchmarsh on the
remaining prime degrees gives \(\lambda\ll FH/\log H\). Since
\(H\ll\lambda^{-1/2}\sqrt\ell\), we obtain
\(F\asymp\lambda^{3/2}\sqrt\ell\).

For each degree define the nonnegative cost
\[
c_i=(1-\theta)(v-(i-1))_++\theta((i-1)-v)_+,
\qquad
M(h)=\sum_{c_i\le h}p_i .
\]
The fractional threshold identity gives
\(\sum_i p_i c_i=v(R(u_*)-\theta L)\). Let
\(h_0=\lambda\ell/(18F)\asymp\lambda^{-1/2}\sqrt\ell\).
The interval \(c_i\le h\) has length \((9/2+o(1))h\).
For \(h_0/\ell\le h\le h_0\) its log-length is
\((1/2+o(1))\ell\), uniformly in \(h\). Brun–Titchmarsh and
the per-degree prime mass cap therefore give prime mass at most
\((18+o(1))Fh/\ell\). The composite mass is \(o(\lambda)\)
throughout these intervals because they lie at degrees
\(O(1/\lambda)\). For \(0\le h<h_0/\ell\), monotonicity gives
\(M(h)\le M(h_0/\ell)=o(\lambda)\). Consequently, uniformly
on \([0,h_0]\),
\[
M(h)\le o(\lambda)+(18+o(1))Fh/\ell.
\]
Integrating over that interval and using \(v\sim3/\lambda\) yields
\[
v(R(u_*)-\theta L)=\int_0^\infty(L-M(h))\,dh
\ge (1-o(1))\frac{\lambda^2\ell}{36F},
\qquad
F(R(u_*)-\theta L)\ge(1-o(1))\frac{\lambda^3\ell}{108}.
\]
The initial upper bound and \(F\le6\delta_L\) give
\(F(R(u_*)-\theta L)\le(42+o(1))\delta_L^2\).
Thus the coefficient of \(\lambda^{3/2}\sqrt\ell\) in the saving
can be taken to be \(1/(18\sqrt{14})-o(1)\) as \(\lambda\to0\).

### A proportional saving from powerful-modulus compression

For every fixed \(0<\lambda\le10^{-6}\) and \(\varepsilon>0\),
\[
 N_\lambda(X)\ll_{\lambda,\varepsilon}
 X^{(2/3-10^{-10})\lambda+\varepsilon}.
 \tag{16}
\]
Use the divided-difference congruence bound in Section 9.3 of
[Additive representations by powerful integers](../../unpublished/powerful-number-representations/main.tex).
For \(k\ge12288\), a \(k\)-full modulus \(n\), and a unit \(K\bmod n\),
it gives \(O_{k,\rho}(n^{1/(2k)+\rho})\) positive unit rational pairs
satisfying \(VT^k\equiv K\pmod n\) when
\[
 H(T)\le n^{\tau/k},\qquad H(V)\le n^{\tau 2^{-16}},
 \qquad \tau=65/64 .
\]
The dimension \(16(\lfloor k/4096\rfloor+1)\le k/192\) for this range,
and the remaining numerical inequalities of that lemma hold there.
Its constants may depend on \(k\).

Write \(\eta=10^{-10}\), \(\zeta=3\cdot10^{-6}\), and
\(D=(2/3-\eta)\lambda\). Work in \(Y/2<c\le Y\).
Choose \(0<\rho\le\min(\eta\lambda/100,\varepsilon/10)\).
Factor each letter as \(\ell=\ell_0\prod_{i=1}^d\ell_i^i\) by
BBLT Lemma 2.2, so that \(\ell_0\le Y^\rho\).
Fix the three coefficients and the dyadic intervals for
the \(\ell_i\). Count the triples admitting coordinates in
that box. Measure every mass \(x_{\ell,i}\) by the logarithm
to base \(Y\) of the lower endpoint. As above let \(p_i=\sum_\ell x_{\ell,i}\),
\(L=\sum_i p_i\), and \(S_\ell=\sum_i x_{\ell,i}\).
The radical comparison gives \(L\le\lambda+3\rho\), and the
height upper bounds give \(W=\sum_i(i-1)p_i\le3-L\).
There are \(O_\rho(Y^{4\rho})\) choices of coefficients and intervals
for large \(Y\). [BBLT v2, Lemma 2.2](https://arxiv.org/html/2410.12234v2#S2)

Retain the actual product of the coordinate scales in the moment
and geometry estimates. Fixing two letters, grouped Fourier
counting, and the geometry estimate give the following box
exponents, each with loss at most \(\rho\):
\[
 S_{\ell_1}+S_{\ell_2},\qquad
 2L/3-F/6,\qquad
 \max\{L-s,L-1+w\}.
\]
The last expression holds for every selection of mass \(s\)
and weighted mass \(w\) as in (9). The geometry denominator
is \(\gg_d Y\): the complete value of \(c\) is at least \(Y/2\)
and the lower and upper weighted scales differ by a bounded
factor at the fixed cutoff \(d\).
The divisor losses can be taken as small as \(\rho\) for each fixed \(d\).
[BBLT v2, proofs of Propositions 3.1 and 4.1](https://arxiv.org/html/2410.12234v2)

If any displayed exponent is at most \(D\), the box has the
required bound. In the remaining boxes the inequalities
\(D\le S_{\ell_1}+S_{\ell_2}\), \(D\le2L/3-F/6\), and (9)
hold with this fixed \(D\). In particular,
\[
 (1-3\eta/2)\lambda\le L\le\lambda+3\rho,\qquad
 S_\ell\le L-D,\qquad J:=2L/3-D\le2\eta\lambda,\qquad
 F\le6J .
\]

Apply (10) at \(u_*=1-L+D\). Since \(D>F\), we have
\(W>u_*\). Set \(\theta=u_*/W\). The weighted height bound
and \(R(u_*)\le L-D+F\) give
\[
 \tfrac14\le\theta\le\tfrac13+15\eta,\qquad
 0\le R(u_*)-\theta L\le F+\frac{3J}{3-L}\le8J .
\]
For the fractional threshold \(v\), the identity in the concentration
argument gives \(E:=\sum_i p_i|1-(i-1)/v|\le 64\eta\lambda\).
We have \(W\ge(1-2\lambda)/(1/3+15\eta)\) and
\(\lvert W/v-L\rvert\le E\). Consequently the center
\(j_*=v+1\) and its weighted absolute deviation satisfy
\[
 \frac{3-1/1024}{\lambda}\le j_*\le\frac{3+1/256}{\lambda},
 \qquad A:=\sum_i p_i|i-j_*|=vE\le 200\eta .
\]
Indeed \(W/(L+E)\le v\le3/(L-E)\);
use \(L+E\le (1+65\eta)\lambda\),
\(L-E\ge(1-131\eta/2)\lambda\), and \(\lambda\le10^{-6}\).

Set \(k=\lfloor(1-\zeta)j_*\rfloor\).
Then \((3-1/256)/\lambda\le k\le(3+1/256)/\lambda\), in
particular \(k\ge12288\). Fix the entire tuple \((c_i)\)
in its box and use the \(k\)-full divisor
\[
 n=\prod_{i\ge k}c_i^i .
\]
The discarded weighted mass is at most \(A/\zeta\), since
\(i<k\) implies \(i\le j_*\) and \(|i-j_*|\ge\zeta j_*\).
For large \(Y\), the loss from \(c_0\) and the dyadic endpoints
gives \(\sum_i i x_{c,i}\ge1-2\rho\).
It follows that \(n\ge Y^{1-2\rho-A/\zeta}\ge Y^{99/100}\).

Compress the other letters by
\[
 T=\prod_i\frac{a_i}{b_i},\qquad
 V=\frac{a_0}{b_0}\prod_i\left(\frac{a_i}{b_i}\right)^{i-k}.
\]
They are unit fractions at \(n\) by primitivity, and
\(VT^k=a/b\equiv-1\pmod n\). This ratio and the fixed \(c\)
recover \(a,b\). We have
\(\log_Y H(T)\le L-D+\rho\le(1/3+1/1024)\lambda\), where \(\rho\)
absorbs the dyadic factors for large \(Y\).
For \(V\) the bound is
\[
 \frac{\log H(V)}{\log Y}
 \le 2\rho+\sum_i|i-k|p_i
 \le 2\rho+A+(\zeta j_*+1)L
 \le 2\rho+200\eta+4\zeta+2\lambda
 <3/200000 .
\]
One \(\rho\) pays for the exceptional coefficient and the other
for the bounded dyadic factors at fixed \(d\).
The exact comparisons
\[
 (3+1/256)(1/3+1/1024)<(99/100)\tau,
 \qquad 3/200000 < (99/100)\tau 2^{-16}
\]
put both heights within the congruence lemma. The fixed tuple
therefore admits \(O_{k,\rho}(Y^{128\lambda/767+\rho})\)
pairs. The number of \(c\)-tuples is \(O_d(Y^{S_c})\), where
\(S_c\le(1/3+2\eta)\lambda\). Since
\(1/3+2\eta+128/767<501/1000\), each such box contributes
\(O_{\lambda,\rho}(Y^{501\lambda/1000+\rho})\), below \(Y^{D+\rho}\). The range of \(k\) is finite for fixed
\(\lambda\), so the modulus-count constants can be maximized
over that range. Summing the boxes and then the dyadic ranges
proves (16).

### A matching order for the Fourier and two-group estimates

There are budgets tending to zero for which the Fourier, geometry,
height and two-letter counting inequalities permit
\(D\ge 2\lambda/3-O(\lambda^{3/2}\sqrt{\log(1/\lambda)})\).
Indeed, partition \([M,2M]\) into intervals of length at most
\(H=\lceil\sqrt{M\log M}\rceil\). By the prime number theorem, some
interval contains \(J\gg H/\log M\) primes. Let their minimum and mean
be \(q_{\min}\) and \(\bar q\), respectively. Set
\(\lambda=3/\bar q\), and give each of these prime degrees the same mass
\(\lambda/(3J)\) on each letter. Every letter has weighted height
exactly one, the total mass is \(\lambda\), and the grouped Fourier
quantity is \(F=\lambda/J\). For any selection of released coordinates,
\(w\ge(q_{\min}-1)s\), so the right-hand side of (9) is at least
\(L-1/q_{\min}\). Thus the profile satisfies all the stated inequalities
at
\[
D=\min\{2\lambda/3-\lambda/(6J),\;\lambda-1/q_{\min}\}.
\]
Both deficits from \(2\lambda/3\) are \(O(H/M^2)\), by the lower bound
for \(J\) and \(0\le\bar q-q_{\min}\le H\). The same order of deficit also suffices for the grouped-root,
direct residue-energy and common-shift collision constraints.

To check these further inequalities, write \(u=\lambda/3\) and
\(E=uH/M+u/J\). Each letter now has mass \(u\) and weighted height
one. Any subset of coordinates with mass \(m\) and weighted height
\(h\) satisfies \(m=uh+O(E)\), uniformly in the subset. Every group
with common divisor \(g>1\) has at most one occupied degree, so has
mass at most \(u/J\). For a group with \(g=1\), its root height is
\(h\), and the inequality
\[
 -m+(h-\tau)_+\ge-u\tau-O(E)
\]
holds for every \(\tau\ge0\). This follows from \(m=uh+O(E)\) and
\(0<u<1\).

Use the exact frozen weights \(w_T\) of this profile, so
\(f_T=uw_T+O(E)\). For the grouped-root collision on \(A,B\),
each energy group saves at most \(uw_B/2+O(E)\) in the energy term.
Thus \(K\ge2(u-f_A)-uw_B-O(E)\). For the direct residue-energy
bound each entire released letter has \(K_T\ge2u-uw_C-O(E)\), so
\(f_C+(K_A+K_B)/2\ge2u-O(E)\).

Every available difference cap on \(T\) is at least \(s_T-O(E)\)
because a divisor \(e>1\) can remove at most one occupied degree
from the complement. Every residue cap is at least
\(s_T-uw_C-O(E)\), by applying the displayed inequality for \(g=1\)
and the one-degree bound for \(g>1\). The target-count cap
\(u-f_C\) is also at least \(s_T-uw_C-O(E)\), since \(s_T\le u\)
and \(f_C=uw_C+O(E)\). The off-diagonal upper expression in the
grouped-root inequality is therefore at least
\[
 f_A+f_B+f_C+\tfrac12(u-f_C+K+\delta)
 \ge 2u+\tfrac12(f_B-uw_B+f_C-uw_C)-O(E)
 \ge 2u-O(E),
\]
using \(s_B=u-f_B\).

For the common-shift inequality,
\(h_T=1-w_T\) and hence \(\sigma=(1-w_A-w_B)_+\).
Its paired cap satisfies
\(\rho\ge s_A+s_B-uw_C-O(E)\).
Put \(w=w_A+w_B\). The off-diagonal expression in this case
is at least
\[
 \frac{3u}{2}+\frac{uw+(1-w)_+}{2}-O(E)\ge2u-O(E).
\]
If \(w>1\), absence of off-diagonal collisions leaves the diagonal
expression \(f_A+f_B+u\ge2u-O(E)\).
The early-plus-tail geometry loosenings obey the same comparison
with \(L-1/q_{\min}\). An exact tagged height-one bound for a group
containing the entire tail is its actual weight on this profile,
while the untagged lower frozen-modulus bounds can only weaken the
collision estimate. These arguments hold for all released and
frozen subsets, including sets that use the entire tail. The family
therefore gives a matching order for these exponent inequalities
even with the new height information.
[Grouped-root, common-shift and direct congruence inequalities](ABC_OMITTED_COORDINATE_GRAPH.md#counting-collisions-by-their-common-shift)

## Exceptional-set bounds from the two tall terms

The [exact certificates with two tall terms](ABC_PAIR_DEVIATION_CERTIFICATE_AUDIT.md#exact-certificates-at-fixed-radical-budgets)
prove that for every \(\varepsilon>0\),
\[
 N_{7/10}(X)\ll_\varepsilon X^{2/5+\varepsilon},
 \qquad N_{4/5}(X)\ll_\varepsilon X^{9/20+\varepsilon}.
\]
These bounds also apply to smaller budgets by inclusion of the sets.
Reweighting the tree leaves gives \(N_{3/4}(X)\ll_\varepsilon
X^{1773/4012+\varepsilon}\).


At budget \(7/10\), Li's radical-budget preprint reports \(109/250\).
His later twisted-energy preprint reports \(24811/50480\) at \(4/5\)
and \(677/1240\) at \(9/10\). The new certificates lower the first
two of those exponents. Reweighting also lowers the
radical-budget preprint's reported \(35/76\) at \(3/4\). [Li, radical budget, v2, Theorem 1.2](https://rstation.ai/problems/math/papers/17/v2) ·
[Li, squares and twisted energy, v3, Theorem 1.2](https://rstation.ai/problems/math/papers/28/v3)


To transfer the certificates, apply the same dyadic and anatomic
reduction used for (5). Keep the shorter addend as the short label
and order the two tall labels as in the certificate. Both tall
monomials have weighted logarithmic height tending to one because
their complete values are at least \(Y/4\) and the exceptional
coefficients cost \(Y^{O(\eta)}\). In a convergent subsequence of the
projected masses, the complete-tail upper heights and the tall
whole-tail lower heights therefore hold. The grouped collision,
geometry, higher residue moments and the reduced-fraction count give their finite
disjunctions in the limit. A count exponent above a verified bound
would produce a feasible exponent point above that bound. The
argument also transfers the reweighted certificates since the
limiting budget row is \(L\le\lambda\).


At budget one, the [height refinements of the earlier tree](ABC_PAIR_DEVIATION_CERTIFICATE_AUDIT.md#height-refinements-of-the-budget-one-tree)
give
\[
 N_\lambda(X)\ll_\varepsilon X^{29/50+\varepsilon}
 \qquad (0<\lambda\le1).
\]
The transfer uses the same two tall terms, with the three choices of
the shorter label made inside the original ordering.

For the full range \(0<\lambda\le1\), the higher residue moments and
the reduced-fraction count give the stronger uniform bound
\(N_\lambda(X)\ll_\varepsilon X^{287/500+\varepsilon}\).
The analytic proof classifies boxes by their largest coordinate
masses and applies over all degrees in the anatomic factorization.
[Uniform analytic bound](ABC_OMITTED_COORDINATE_GRAPH.md#a-uniform-analytic-exponent-of-287500)

## A small-budget bound from higher moments and two tall terms

For all sufficiently small fixed radical budgets \(\lambda>0\), the
anatomic box exponent satisfies
\[
D\le \lambda/2+O(\lambda^{3/2}).
\]
Consequently \(N_\lambda(X)\ll_{\lambda,\varepsilon}
X^{\lambda/2+C\lambda^{3/2}+\varepsilon}\) for an absolute constant
\(C\).

Work in the full anatomic coordinates, with mass \(\alpha_{T,i}\)
in degree \(i\), letter totals \(L_T\), and \(L=\sum_T L_T\le\lambda\).
If some frozen coordinates of letter \(C\) have weight at least
\(2\lambda\) and total mass \(f_C\), the uniform higher-moment bound
gives residue energy exponents at most \(L_A,L_B\) for the two
free letters. The direct residue bound therefore gives
\[
 D\le (L-L_C)/2+f_C. \tag{A}
\]
This applies whichever letter is locked. [Uniform energy bound and its box transfer](ABC_OMITTED_COORDINATE_GRAPH.md#uniform-energy-bounds-from-interpolation)

### One possible high-degree coordinate per letter

Put \(K=\lceil\lambda^{-1/2}\rceil\). It suffices to consider
\(D>L/2+3\lambda/K\). For each letter \(T\), set
\[
\theta_T=D-(L-L_T)/2-3\lambda/K>0.
\]
The coordinates of degree \(i\ge K\) with mass at most \(\theta_T\)
have total weight less than \(2\lambda\). Otherwise freeze them
until their weight first reaches \(2\lambda\). Before the last
coordinate the frozen mass is less than \(2\lambda/K\), and the
last costs at most \(\theta_T\). This contradicts (A). On each
letter at most one coordinate can have mass greater than
\(\theta_T\), since \(2\theta_T>L_T\). Call the exceptional
high-degree coordinate dominant when it exists. The total
weight outside the dominant coordinate on any letter is
\(O(\sqrt\lambda)\), since low degrees contribute at most
\(KL_T\) and the other high degrees less than \(2\lambda\).

Order the positive triple so that \(a\le b\) and the largest
term \(c\) is on the dyadic scale \(X\). The two tall letters
\(b,c\) have height one in the limiting normalized box, and thus
each has a dominant coordinate of degree at least \(1/(2\lambda)\).
Write their masses as \(\beta,\gamma\).

If all three letters have dominant coordinates, freeze every other
coordinate. The remaining equation is a generalized Fermat surface
in three primitive roots, with arbitrary nonzero coefficients. Use
Browning–Verzobio's coefficient-uniform Theorem 1.2, distinguishing
one tall degree as \(r\) and the other as \(p\). The root box
sides are at most \(X^\lambda\). Because \(p,r\ge1/(2\lambda)\),
each of its three terms has exponent \(O(\lambda^{3/2})\). Hence
\[
D\le L-\sum_T\beta_T+O(\lambda^{3/2}),
\]
where \(\beta_T\) are the three dominant masses. As
\(\sum_T\beta_T\ge\sum_T\theta_T=3D-L-9\lambda/K\), this gives
\(D\le L/2+O(\lambda^{3/2})\). [Browning–Verzobio, Theorem 1.2](https://arxiv.org/html/2608.24512v1#S1)

### Two-tall strip bound

In the other case the short letter has no dominant coordinate and
\(a\le X^{h}\) with \(h=O(\sqrt\lambda)\). Freeze the non-dominant
coordinates on the tall pair. Their two remaining positive roots
\(u,v\) lie in dyadic ranges of scales \(U,V\), with \(UV\le X^\lambda\)
and \(\gcd(u,v)=1\). The tall powers obey
\[
0<Cv^q-Bu^p\le X^h,\qquad Bu^p,Cv^q\asymp X,
\qquad p,q\ge1/(2\lambda).
\tag{B}
\]
The number of such pairs is \(O_{p,q,\lambda}(X^{O(\lambda^{3/2})})\),
uniformly in the positive frozen coefficients \(B,C\).

First let \(g=\gcd(p,q)\ge3\). The relative difference between
the ratios \(u^{p/g}/v^{q/g}\) of any two pairs is
\(O(X^{-1+h})\) by (B). A nonzero relative difference is at least
\(\gg X^{-2/g}\), since its reduced rational denominator is at
most the product of two root powers, each \(O(X^{1/g})\).
For small \(\lambda\), all ratios are therefore equal. Coprimality
of the two roots makes each ratio a reduced fraction, so the pair
is unique.

Suppose \(g\le2\). Normalize by \(t=u/U\), \(s=v/V\), and put
\(s_0(t)=\kappa t^{p/q}\) for the exact positive power graph.
On solutions of (B), \(s-s_0(t)=O(X^{-1+h})\), and \(\kappa\)
is bounded above and away from zero by constants depending on the
fixed powers. The same bound gives at most one integer \(v\)
for each \(u\). Take
\(m=\lfloor(10\lambda)^{-1/2}\rfloor-1\) and \(k=(m+1)^2\).
The \(k\) monomials \(t^is^j\) with \(0\le i,j\le m\) restrict
to distinct powers on the exact graph, because
\(p/g,q/g\gg1/\lambda\) and \(m=O(\lambda^{-1/2})\).
On \(k\) distinct normalized first coordinates \(t_0<\cdots<t_{k-1}\),
apply the triangular divided-difference row transform to their
monomial matrix. Its determinant divides the original determinant
by the ordinary Vandermonde in the \(t_j\).

For the exact graph this transformed determinant is bounded away
from zero. Indeed, with \(x_j=t_j^{1/(q/g)}\), the original
determinant is a positive column factor times a generalized
Vandermonde in the \(x_j\). It factors as the ordinary Vandermonde
in the \(x_j\) times a Schur polynomial with nonnegative
coefficients. The polynomial is at least one on \(x_j\ge1\) and
the quotient of the two ordinary Vandermondes is bounded away
from zero on the dyadic range. [Macdonald, Chapter I, Sections 3 and 5](https://math.berkeley.edu/~corteel/MATH249/macdonald.pdf)

The transformed exact entries
are bounded by derivatives of the fixed monomials up to order
\(k-1\). For the perturbed graph the entry in row \(j\) differs
by \(O(X^{-1+h}U^j)\) from the exact one, using the explicit
Lagrange formula for a divided difference and the spacing
\(t_{i+1}-t_i\ge1/U\). This error tends uniformly to zero since
\(\lambda k\le1/10\). The transformed perturbed determinant is
therefore nonzero and bounded above.

If all \(k\) values of \(t\) lie in an interval of length
\(\Delta=X^{-2\lambda/(m+2)}\), undo the row transform and then
multiply the columns by \(U^iV^j\). The result is a nonzero
integer determinant of absolute value
\[
\ll (UV)^{mk/2}\Delta^{k(k-1)/2}
\le X^{-mk\lambda/2}.
\]
The inequality uses \(k-1=m(m+2)\), and is impossible for
large \(X\). Each normalized interval of length \(\Delta\)
has at most \(k-1\) points, giving the stated strip bound. [Huxley, introduction on determinant methods for points near curves](https://www.numdam.org/item/ASNSP_1994_4_21_3_357_0/)

Restoring the frozen coordinates on the tall pair, which determines
the whole short term by subtraction, now gives
\[
D\le L_B+L_C-\beta-\gamma+O(\lambda^{3/2}).
\]
Release all coordinates other than the two dominant coordinates
in the geometry row as well. The total released
weight is \(O(\sqrt\lambda)\), so its second term is negative and
its first term yields \(D\le\beta+\gamma\). Averaging these bounds
proves \(D\le L/2+O(\lambda^{3/2})\). For the anatomic reduction, run the
argument with the budget \(\lambda+O(\eta)\) at each fixed truncation
and subtract \(O(\eta)\) from the box count exponent to absorb the
arithmetic losses. The tall heights are then at least \(1-O(\eta)\),
which still forces both tall dominant degrees to be \(\gg1/\lambda\).
The determinant estimates retain their strict margins, and the
weighted height of the short term remains \(O(\sqrt\lambda)\) in the
two-tall case. Taking \(\eta\) sufficiently small in terms of
\(\lambda,\varepsilon\), and summing the height shells, gives the
asserted count.
[BBLT v2, Proposition 4.1](https://arxiv.org/html/2410.12234v2#S4)

### The leading coefficient for whole-coordinate bounds

For integers \(p\) tending to infinity, take \(\beta=1/(p+1)\) and
\(\lambda=6\beta\). Give each letter mass \(\beta\) in degree one and
mass \(\beta\) in degree \(p\). All heights are one. These profiles
permit \(D=3\beta=\lambda/2\) under the whole-coordinate bounds used
above.

The pair and Fourier upper exponents are \(4\beta\) and \(7\beta/2\).
A geometry release including degree \(p\) on any letter contributes
weighted excess at least \(1-2\beta\), and a release with only
degree-one coordinates has mass at most \(3\beta\). For the direct
residue row, a nonempty lock costs at least \(\beta\), and the ideal
energy on each released letter is at least its diagonal mass
\(2\beta\). The empty lock gives no congruence saving. In a grouped
collision, releasing at most \(\beta\) mass on the second letter gives
diagonal exponent at least \(3\beta\). Releasing it entirely sets
the energy modulus to one, so the off-diagonal exponent is at least
\(3\beta\). For a paired shift, the diagonal exponent suffices
unless both free letters are entirely released. Then the shift
exponent is one and the off-diagonal bound is larger than \(3\beta\).

For a partial-modulus rational row with common power \(g>1\), each
free group releases at most \(\beta\) mass. For \(g=1\), a total
release of \(3\beta\) has weighted height at least \(1+\beta\), and
releasing \(4\beta\) has height two. With lock height at most one
and nonempty lock cost at least \(\beta\), these cases also have
right-hand side at least \(3\beta\). With an empty lock its height
penalty cancels at least the released mass. Freezing the three
degree-one coordinates for the generalized Fermat estimate costs
\(3\beta\). This gives the stated obstruction within these
exponent bounds.
