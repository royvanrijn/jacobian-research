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
\begin{aligned}
 N_{7/10}(X)&\ll_\varepsilon X^{2/5+\varepsilon},&
 N_{4/5}(X)&\ll_\varepsilon X^{9/20+\varepsilon},\\
 N_{9/10}(X)&\ll_\varepsilon X^{13/25+\varepsilon},&
 N_1(X)&\ll_\varepsilon X^{11/20+\varepsilon}.
\end{aligned}
\]
These bounds also apply to smaller budgets by inclusion of the sets.
Reweighting the tree leaves gives \(N_{3/4}(X)\ll_\varepsilon
X^{1773/4012+\varepsilon}\).


At budget \(7/10\), Li's radical-budget preprint reports \(109/250\).
His later twisted-energy preprint reports \(24811/50480\) at \(4/5\)
and \(677/1240\) at \(9/10\). The new certificates lower all
three exponents. Reweighting also lowers the
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
the reduced-fraction count also give a uniform analytic bound
\(N_\lambda(X)\ll_\varepsilon X^{4/7+\varepsilon}\).
The analytic proof classifies boxes by their largest coordinate
masses and applies over all degrees in the anatomic factorization.
[Uniform analytic bound](ABC_OMITTED_COORDINATE_GRAPH.md#a-uniform-analytic-exponent-of-47)

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

The energy estimate needed in (A) has a direct integer-moment proof.
For a free letter of total mass \(y>0\), choose
\(h_i=\lceil y/\alpha_{T,i}\rceil\) on its positive-mass coordinates.
Then \(\sum_i 1/h_i\le1\) and \(h_i\alpha_{T,i}\le2y\).
If the modulus has height exponent \(w\ge2\lambda\ge2y\),
all extra moment terms in (6) vanish.
For uniform constants at fixed truncation, freeze the coordinates
of sufficiently small mass in both energy copies. Their total
cost can be made arbitrarily small and the remaining moment orders
are bounded in terms of the truncation and the allowed loss.
Congruence energies of the free letters bound the number of pairs
by Cauchy--Schwarz, and the chosen coordinates of the locked letter
specify the modulus. Once the free pair is fixed, the locked letter
is determined by the equation.

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

## Scalable amplification and stability near the half bound

### Residue counts at shrinking residual heights

Let \(k,\ell\ge1\) be integers with \(k\ge1024\ell^2\). For a \(k\)-full
modulus \(n\), a unit \(K\bmod n\), and positive rational unit pairs,
\[
 VT^k\equiv K\pmod n,\qquad
 H(T)\le n^{2/k},\qquad H(V)\le n^{1/(16384\ell^3)}
\]
has \(O_{k,\ell,\varepsilon}(n^{1/(\ell k)+\varepsilon})\) solutions.
Use the amplification, residue partition and Newton comparison in
Section 9.4 of
[Additive representations by powerful integers](../../unpublished/powerful-number-representations/main.tex),
with
\[
 \alpha=1,\qquad m=16\ell-1,\qquad
 h=\lfloor k/(256\ell^2)\rfloor,\qquad D=16\ell(h+1).
\]
Here \(\alpha\) is the partition depth, \(\ell\) is the amplification
order, and the rectangle uses monomials \(V^iT^j\) with
\(0\le i\le m,\ 0\le j\le h\).
Write \(q=k\log Q_1/\log n\), so \(1/2\le q\le1\).
Distinct amplified pairs have distinct rational \(T\), since \(T\)
determines \(V\bmod n\) and \(2\ell/(16384\ell^3)<1\).
The comparison weights \(mk-ki+j\) are distinct since \(h<k\).

The determinant expansion in comparison and error rows gives the
full Vandermonde divisibility at the base depth, since \(D-1\le k\).
The per-row exponents satisfy the exact bounds
\[
 \ell\left(2h/k+\frac{m}{16384\ell^3}\right)
 \le \frac{9}{1024\ell},\qquad
 \frac{q(D-1)}{2k}\ge\frac1{64\ell}-\frac1{4k}.
\]
The gap is at least \(27/(4096\ell)\), absorbing the factorial of
the dimension for large \(n\). For any \(2D\) amplified rows in
one class, the weighted prime average and minimum-energy selection
give a \(D\)-row subset at one prime with maximum normalized
row valuation at most
\[
 \frac{D-1}{k}+
 \frac{(2D-1)(4\ell-q)}{k}\frac{\log n}{\log n_g},
\]
where \(n_g\) is the full-power part over eligible primes in the
Newton comparison. In the limit as \(n\) grows for fixed \(k,\ell\),
this is at most \(8\ell D/k\le5/8\), using
\(D/k\le1/(16\ell)+16\ell/k\). The product of the excluded prime
powers is bounded in \(k,\ell\) by a fixed power of the diagonal
Vandermonde constant. The deeper residue partition multiplies
\(Q_1\le n^{1/k}\) by a bounded factor as well.
The Newton comparison therefore forces nonvanishing, while the
size bound forces vanishing. There are fewer than \(2D\) amplified
pairs per class, and undoing amplification gives the stated count.
All bounded ranges of \(n\) are absorbed in the constant.

For integers \(k\ge2\) and \(L\ge1\), let \(N_{k,L}(X)\) count
primitive triples \(a+b=c\le X\) with every nonzero prime valuation
in \([k,k+L]\). Then
\[
 N_{k,L}(X)\ll_{k,L,\varepsilon}
 X^{1/k+O(L^{1/3}k^{-4/3})+\varepsilon}.
 \tag{17}
\]
For each letter \(m\) put \(R_m=\operatorname{rad}(m)\) and
\(W_m=m/R_m^k\). At fixed \(c\), use \(n=c\), \(T=R_a/R_b\) and
\(V=W_a/W_b\). Then \(VT^k=a/b\equiv-1\bmod c\),
\(H(T)\le c^{1/k}\) and \(H(V)\le c^{L/k}\).
For \(k\ge16384L\), take
\(\ell=\lfloor(k/(16384L))^{1/3}\rfloor\).
The two conditions on \(k\) and the residual height hold, and
each compressed pair determines the triple at fixed \(c\).
Summing over \(k\)-full \(c\le X\) gives exponent \(1/k+1/(\ell k)\)
before the \(X^\varepsilon\) loss. The trivial count covers
\(k<16384L\) after increasing the absolute constant in (17).
For fixed \(L\) and \(\lambda=3/k\), these triples have radical
budget \(\lambda\), and (17) has exponent \(\lambda/3+o(\lambda)\).

The power in the congruence can be an independent positive
integer \(q>h=\lfloor k/(256\ell^2)\rfloor\).
For \(VT^q\equiv K\bmod n\), the comparison weights become
\(mq-qi+j\), which are distinct by \(q>h\).
The denominator bounds, partition at depth \(\alpha=1\),
row-energy estimate, and Newton precision are unchanged.
The Vandermonde constant is now fixed by \(k,\ell,q\), so the
bound is \(O_{k,\ell,q,\varepsilon}(n^{1/(\ell k)+\varepsilon})\),
uniformly in the unit \(K\).

### The stability supplied by the half-bound proof

Use limiting boxes for the moment, with \(L\le\lambda\).
If \(D=\lambda/2+o(\lambda)\), the half-bound proof gives
\(L=\lambda+o(\lambda)\). Put
\(\Delta=(L/2-D)_+\) and \(K_0=\lceil\lambda^{-1/2}\rceil\).
On a letter with \(L_T\gg\lambda\), call a coordinate heavy if its
degree is at least \(K_0\) and its mass exceeds
\[
 L_T/2-\Delta-3\lambda/K_0.
\]
The stopping argument for (A) applies with this threshold.
The other high degrees together have height less than \(2\lambda\),
and the low degrees contribute height at most \(K_0 L_T\).
Thus the height outside the heavy coordinates is \(O(\sqrt\lambda)\).
There are at most two heavy coordinates, each of mass at least
\(L_T/2-o(\lambda)\). If there are two, each has mass
\(L_T/2+o(\lambda)\) and the remaining mass is \(o(\lambda)\).

If all three letters have mass \(\gg\lambda\), the short letter must
also have a heavy coordinate. Otherwise it has height \(O(\sqrt\lambda)\).
Choose a heavy coordinate carrying weight \(\gg1\) on each tall
letter, which is possible since at most two cover its height
up to \(O(\sqrt\lambda)\). The strip estimate with just these
two coordinates free gives
\(D\le (L_B+L_C)/2+o(\lambda)\), contradicting \(L_A\gg\lambda\).
When every letter has a heavy coordinate, free one per letter,
choosing those of large height on the tall pair. The generalized
Fermat surface bound gives \(D\le L-\sum_T\beta_T+o(\lambda)\),
where \(\beta_T\) are their masses.
Together with the lower bounds for the \(\beta_T\), this forces
each chosen mass to be \(L_T/2+o(\lambda)\). The argument controls
the number and mass of the heavy coordinates. Their degrees are
allowed to vary across the letters.

### Separated degrees and a limit of single-base compression

Take \(N\) tending to infinity and let the three high degrees be
\(d_T=\lfloor\gamma_T N\rfloor\), with
\((\gamma_A,\gamma_B,\gamma_C)=(1,\sqrt2,\sqrt3)\).
On each letter place mass \(\beta_T=1/(d_T+1)\) at degree one
and mass \(\beta_T\) at degree \(d_T\).
The heights are one, \(L_T=2\beta_T\), and set
\(\lambda=L=2\sum_T\beta_T\) and \(D=\sum_T\beta_T\).
For large \(N\), \(2\max_T\beta_T\le\sum_T\beta_T\), so the
two-letter counts allow this value of \(D\). The Fourier row has
\(F=\sum_T\beta_T\), and allowing a high coordinate in a geometry
release makes the second term at least \(L-2\beta_T\ge D\).
Releasing only degree-one terms gives the bound \(L-\sum_T\beta_T=D\).
Every nonempty whole-coordinate residue lock costs at least
\(\beta_T=L_T/2\), so even the diagonal energy exponents on the
free letters give a right side at least \(D\).
Freezing the degree-one coordinates for the surface count costs
exactly \(D\). These profiles satisfy the counting inequalities
used in the half-bound proof with substantial mass at separated
degrees of order \(1/\lambda\).

On such profiles, consider fixing a complete letter, taking any
nontrivial \(k\)-full divisor \(n\) of it with height exponent
\(w=\log_X n\), and compressing the other two letters by a single
base \(k\). Allow arbitrary integer multiplicities \(a,b\) in \(T\)
for the high-degree coordinates of the free letters, of degrees
\(p\le q\), respectively. Their residual exponents are \(p-ka\)
and \(q-kb\). In the coprime anatomic coordinates the height caps
of the scalable lemma force
\[
 |p-ka|\le(p+1)w\nu,\qquad |q-kb|\le(q+1)w\nu,\qquad
 \frac{b}{q+1}\le\frac{2w}{k},\qquad
 \nu=\frac1{16384\ell^3}.
\]
The high coordinates have disjoint prime supports, so each of
these separate contributions is bounded by the corresponding
rational height. In particular \(w\ge q/((q+1)(2+\nu))\ge12/25\)
for large \(N\).

Fixing the letter costs \(L_T=2\beta_T\). To improve \(D\) with
the congruence bound requires \(w/(\ell k)\le D-2\beta_T\).
But \((D-2\beta_T)N\le6/5\) for every choice of the fixed letter
when \(N\) is large, and \(q+(q+1)\nu\le2N\).
Thus \(\ell k/N\ge2/5\) and \(1\le b\le5\ell\).
The full divisor has \(k\le d_T\le2N\), so
\(k\ge1024\ell^2\) implies \(\ell^2\le N/512\).

Let \(\gamma\) be the limiting ratio \(p/q\).
It has \(\gamma^2=P/Q\) with \(P,Q\in\{1,2,3\}\) and \(P<Q\).
The residual inequalities give
\[
 |a/b-p/q|\le3\nu,\qquad
 |p/q-\gamma|\le1/(N-1).
\]
The first follows by writing \(a/b=(p+O((p+1)\nu))/
(q+O((q+1)\nu))\) and using \(p\le q\) and large \(q\);
the second follows from the degree floors.
In particular \(a/b\le2\). Since \(P b^2-Q a^2\) is a nonzero
integer,
\[
 |\gamma-a/b|
 \ge\frac{1}{9b^2}\ge\frac{1}{225\ell^2}.
\]
The two upper bounds instead give
\[
 |\gamma-a/b|\le
 \frac{3}{16384\ell^3}+\frac{1}{N-1}
 \le\left(\frac3{16384}+\frac1{500}\right)\frac1{\ell^2},
\]
a contradiction for sufficiently large \(N\).
This is a limit of counting by fixing one complete letter and
using the scalable congruence bound with single-base monomial
compression. When \(k\asymp1/\lambda\), the short-height
condition on \(T\) holds, and \(H(V)\le n^{o(1)}\),
one can let \(\ell\) tend to infinity slowly enough to satisfy
both \(1024\ell^2\le k\) and the residual cap. The congruence
cost is then \(o(\lambda)\) beyond the cost of fixing the target.
For example, put the high roots of the common-degree profile
above in \(T\), the degree-one factors in \(V\), and use
the high-power part of \(c\) as the modulus. This gives exponent
\(\lambda/3+o(\lambda)\) in that profile.


## A strict small-budget saving for two tall powers and a short letter

Assume \(0<\lambda\) is sufficiently small. In a fixed anatomic box,
suppose \(a\le b\), the short letter \(a\) has weighted height
\(O(\sqrt\lambda)\), and each tall letter \(B,C\) has a chosen
coordinate outside which its weighted height is \(O(\sqrt\lambda)\).
Write \(p,q\) for the degrees of those two coordinates and
\(\beta_B,\beta_C\) for their masses. Both degrees are
\(\gg1/\lambda\) because the two complete tall terms have height
one and \(\beta_T\le\lambda\). In this part of the count one has
\[
 D\le \frac49\lambda+O(\lambda^{3/2}).
\]
The constants may depend on the fixed constants in the assumed
short and complementary height bounds.

### Forcing large prime factors

Write \(S=\beta_B+\beta_C\), \(f_T=L_T-\beta_T\) for \(T=B,C\).
For a sufficiently large fixed \(C_0\) take
\(d_*=4\lambda/9+C_0\lambda^{3/2}\). Apply the pair bounds, the
geometry release of everything except the two chosen coordinates,
and the two-tall strip bound above. If any already gives the
target, use that estimate. In the remaining boxes
\[
 L_T\le L-d_*,\quad S\ge d_*,\quad
 f_B+f_C\ge d_*-E,\qquad E=O(\lambda^{3/2}).
\]
For either tall letter \(T\), let
\(\theta_T=d_*-(L-L_T)/2\). Locking the whole chosen coordinate
gives the target if \(\beta_T\le\theta_T\), since its weight is
at least \(2\lambda\). Thus assume \(\beta_T\ge\theta_T\) on
both tall letters. For \(U\) the other tall letter,
\[
 f_T\ge d_*-L_U+\beta_U-E
 \ge 2d_*-\tfrac12(L+L_U)-E
 \ge \tfrac52d_*-L-E,
 \qquad f_T\le L-\tfrac32d_* .
\]
In particular \(f_B\le3f_C\) and \(f_C\le3f_B\) for large enough
\(C_0\). Also
\[
 \theta_T-\beta_T/2\ge \tfrac94d_*-L-E/2
 \gg \lambda^{3/2}.
\]

For the degree \(e_T=p\) or \(q\) on \(T\), a divisor of its chosen
coordinate with mass in \([2\lambda/e_T,\theta_T]\) gives a lock
of weight at least \(2\lambda\) and cost at most \(\theta_T\).
Bound all tuples with such a divisor by (A), summing over the
divisor's dyadic scales. In the remaining tuples the chosen
coordinate equals \(r_T s_T\), where \(r_T\) is prime and
\(s_T\le X^{2\lambda/e_T}\) up to dyadic constants. Indeed, a
prime factor in the divisor interval itself suffices. If all
prime factors are less than its lower endpoint, the first
partial product reaching that endpoint is still in the interval.
A prime exceeding the upper endpoint leaves a cofactor below
the upper endpoint, since \(\theta_T>\beta_T/2\). That cofactor
must be below the lower endpoint. Fix the two cofactors
\(s_B,s_C\) at total cost \(O(\lambda^2)\) in the exponent.
Their contribution to the weighted heights is \(O(\lambda)\).

### A graph for the remaining prime factors

The tall terms now have the form \(P r^p,Q s^q\), with positive
coefficients \(P,Q\) obtained from their nonchosen coordinates
and the fixed cofactors. Their numbers of possible values have
exponents \(f_B,f_C\), respectively. Their numerical height
exponents are \(O(\sqrt\lambda)\). Every solution satisfies
\[
 0<Qs^q-Pr^p\le X^{h},\qquad h=O(\sqrt\lambda),\qquad
 Pr^p,Qs^q\asymp X .
\]
Graph the possible \(P,Q\), using a selected witness \((r,s)\)
for each edge. The two-tall strip estimate bounds its witness
multiplicity by \(X^{O(\lambda^{3/2})}\), and \(a\) is determined
by the difference of the two tall values.

To bound common neighbours of \(P_1,P_2\), put
\(U=r_1/r_2,\ V=s_1/s_2\) using the two witnesses.
The relation \(U^p/V^q=P_2/P_1\) holds up to relative error
\(O(X^{-1+h})\). If \(g=\gcd(p,q)\ge5\), the rational
\(U^{p/g}/V^{q/g}\) has numerator and denominator
\(O(X^{2/g})\) and lies within relative error \(O(X^{-1+h})\)
of a fixed value. Two distinct such rationals would be at
distance \(\gg X^{-4/g}\) from each other, so at most one occurs.
Unless a prime repeats, factor its numerator and denominator
to recover the four primes up to finitely many assignments.
Within either witness the two primes are distinct by
primitivity. If one repeats across the two witnesses, the
fraction still recovers the two primes on the other side.

For \(g\le4\) the rational pairs \(U,V\) themselves have count
\(X^{O(\lambda^{3/2})}\). To see this, use the divided-difference
determinant from the strip estimate, now on \(U,V\) in fixed
positive compact intervals. Both rationals have height
\(O(X^\lambda)\). Distinct first coordinates have spacing
\(\gg X^{-2\lambda}\), and for a fixed \(U\) the second
coordinate is unique at the error \(O(X^{-1+h})\). Take a
square array of monomials \(U^iV^j\), \(0\le i,j\le m\), with
\(m=\lfloor(20\lambda)^{-1/2}\rfloor-1\). Their powers on
the exact graph are distinct, since \(p/g,q/g\gg1/\lambda\).
Divide the rows by the ordinary Vandermonde using divided
differences. The transformed determinant on the exact graph
is bounded above and away from zero, as in the strip proof.
The errors tend to zero since \(2\lambda(m+1)^2\le1/10\).
After restoring denominators, the determinant is an integer.
With \(n=(m+1)^2\), clearing the denominators costs at most
\(O(1)X^{2\lambda mn}\). Choose first-coordinate intervals of
length \(X^{-6\lambda/(m+2)}\). The Vandermonde saving makes
the integer determinant less than one, so each interval
has fewer than \(n\) points. Factoring the reduced ratios of
primes now gives the same recovery, including the possibility
of one repeated prime.

If neither prime repeats, one recovered witness determines
\(Q\) to an interval of length \(O(X^{h-q\log_X s_1})=o(1)\).
If both repeat, it forces \(P_1=P_2\) by the same height
comparison. With exactly one repeat, the remaining pair
consists of a prime \(z\ll X^\lambda\) and a coefficient
\(Q\ll X^{O(\sqrt\lambda)}\) near a fixed positive monomial
in \(z\) of degree \(p\) or \(-q\). The normalized error is
again \(O(X^{-1+h})\). For this pair take monomials of
degree at most \(m_1\) in \(z\) and \(m_2\) in \(Q\), with
\(m_1\asymp\lambda^{-3/4}\), \(m_2\asymp\lambda^{-1/4}\)
and \((m_1+1)(m_2+1)\lambda\le1/10\).
The powers along the graph are distinct because \(m_1<p,q\).
In the negative-power case multiply the exact evaluation
rows by a positive power of the normalized first coordinate
to apply the same Schur quotient. In either case the quotient
is bounded away from zero. Divided differences amplify the
error by at most \(X^{(m_1+1)(m_2+1)\lambda}\).
The logarithmic scale per row to restore the monomials is
\(O(m_1\lambda+m_2\sqrt\lambda)=O(\lambda^{1/4})\).
Intervals for the normalized first coordinate of length
\(X^{-C\lambda^{5/4}}\), for a sufficiently large fixed \(C\),
therefore make the nonzero integer determinant less than one.
This bounds the repeated-prime cases by \(X^{O(\lambda^{5/4})}\).

The resulting common-neighbour bound holds in either
direction. For \(f_B\le f_C\), Cauchy–Schwarz on the degrees
of the \(Q\) vertices gives an edge exponent at most
\(\max(f_C,f_B+f_C/2)+O(\lambda^{5/4})\).
The same statement holds after exchanging the two letters.
Since the larger frozen mass is at most three times the
smaller, this is at most \(3(f_B+f_C)/4+O(\lambda^{5/4})\).
But \(f_B+f_C=L_B+L_C-S\le\lambda-d_*\). Restoring the
cofactors and multiplicity therefore gives an exponent
strictly below \(d_*\) for small enough \(\lambda\). The
argument uses all original degree indices, with the same
fixed-truncation and shell losses as in the half-budget bound.

## A strict small-budget saving for a common high degree

For a primitive anatomic box at radical budget \(\lambda\), suppose
each letter has a distinguished coordinate of the same degree \(k\)
and its other coordinates have total weighted height \(O(\sqrt\lambda)\).
The two tall letters have limiting height one. Then the box exponent
satisfies \(D\le(5/11+o(1))\lambda\) as \(\lambda\to0\).
This applies to the stated common-degree boxes, with no height-one
assumption on the short letter.
[BBLT, Propositions 2.1 and 4.1](https://arxiv.org/html/2410.12234v2)
[Browning–Verzobio, Theorem 1.2](https://arxiv.org/html/2608.24512v1#S1)

### Masses and divisor locks

Let \(L\le\lambda\) be total mass and \(L_T=b_T+\ell_T\), with
\(b_T\) the distinguished mass. Denote the short letter by \(A\),
write \(s=b_A\), and put \(t=1/k\). The tall masses are
\(t+O(\lambda^{3/2})\), \(s\le t\), \(k\gg1/\lambda\), and
\(S=\sum_T b_T=2t+s+O(\lambda^{3/2})\). We may restrict to
\(L\asymp\lambda\) by the pair bounds. Fix \(0<\delta<1/22\),
put \(v=1/2-\delta\), and suppose \(D\ge vL\).
Geometry releasing all non-distinguished coordinates gives \(D\le S\).
Freezing those coordinates and applying the coefficient-uniform
generalized Fermat bound gives \(D\le L-S+O(\lambda^{3/2})\).
In particular
\[
vL\le S\le(1/2+\delta)L+o(L).
\]
The higher-moment residue lock, freezing the whole distinguished
coordinate of \(T\), gives \(D\le(L-L_T)/2+b_T\) when \(k b_T\ge2\lambda\).
If that fails, only \(A\) can fail it. Its height is then
\(O(\sqrt\lambda)\), and the half-budget strip argument, paired with
the geometry release of all coordinates except the two tall roots,
gives \(D\le(L-L_A)/2+o(L)\). Hence in all cases
\[
\ell_T\le b_T+2\delta L+o(L).
\]

Fix a sufficiently small \(\eta>0\), depending on \(\delta\), and put
\(\tau=2\lambda/k=O(\lambda^2)\) and
\(\theta_T=(\ell_T+b_T)/2-(\delta+\eta)L\).
If the distinguished root has a divisor of mass in
\([\tau,\theta_T]\), freeze that divisor at degree \(k\).
Its modulus weight is at least \(2\lambda\), so the higher-moment
bound gives \((L-L_T)/2+\theta_T=(v-\eta)L\) for this subfamily.
Summing over dyadic divisor ranges costs \(X^{o(1)}\).
Remove these subfamilies on each letter.

Where \(\theta_T>2\tau\), the remaining primes of its root have
mass below \(\tau\) or above \(\theta_T\). The product of the
small ones, counted with multiplicity, has mass below \(\tau\)
by greedily stopping at the first crossing of \(\tau\).
Freeze these tiny remainders at cost \(O(\lambda^2)\).
If \(\ell_T>(2\delta+3\eta)L\), then
\(\theta_T>b_T/2+\eta L/2\). The remaining high root is therefore
prime whenever \(b_T\gg L\). Call these letters prime. A short
letter satisfying this test has \(b_A\gg L\) by the bound on
\(\ell_A\). A tall letter always has \(\theta_T>2\tau\),
even when it fails the prime test, because
\(t\ge S/3+o(L)\ge vL/3+o(L)\) and \(v/6>\delta\).
On a composite letter, a nonconstant reduced ratio of two remaining
high roots leaves a common scalar of mass at most
\[
g_T\le (b_T-\ell_T)/2+(\delta+\eta)L+o(L).
\]
Indeed the numerator or denominator contains a prime above
\(\theta_T\). On a prime letter take \(g_T=0\).

### The residue fibres and the graph

Let a selected prime modulus root \(P\) have mass \(\rho\gg L\).
Its modulus \(P^k\) has weight \(k\rho\gg1\).
The whole low cofactor \(x\) has numeric height \(X^{a_0}\) with
\(a_0=O(\sqrt\lambda)\), even though it has at most \(X^{\ell_T+o(1)}\)
possible values. For a unit \(\gamma\) and root range \(z\ll X^b\),
both congruences
\[
xz^k\equiv\gamma\pmod{P^k},
\qquad x\equiv\gamma z^k\pmod{P^k}
\]
have at most \(X^{(b-\rho)_++O(\lambda^{3/2})+\varepsilon}\) solutions.

Partition the \(z\)-range into blocks of length \(X^h\), with
\(h=\rho-C\sqrt\lambda/k>0\) for a sufficiently large constant \(C\).
For the inverse congruence use \(x^i z^j\) with \(0\le i\le R\),
\(0\le j\le ki\), and group by \(j-ki\). In each group of length
\(m\), the unit \(z\)-factor leaves consecutive powers of \(xz^k\).
Expand those in powers of \(xz^k-\gamma\). The resulting integer
column operations gain \(P^{k\binom m2}\), with
\(\sum\binom m2=kR^3/6+O_k(R^2)\).
For the archimedean bound translate \(z\) to the block origin,
using that the support is downward closed in \(j\). Its leading
base-\(X\) exponent is \((ka_0/3+k^2h/6)R^3\). Hence
\(2a_0+kh<k\rho\) forces vanishing for a sufficiently large fixed \(R\).
For the positive congruence take \(j\le k(R-i)\), group by
\(j+ki\), and expand powers of \(x/z^k-\gamma\).
The archimedean leading exponent is \((ka_0/6+k^2h/6)R^3\),
so \(a_0+kh<k\rho\) suffices.

The points in a block lie on a bounded-degree auxiliary curve.
Primitive factors of degree less than \(k+1\) in the inverse case
or \(k\) in the positive case have only boundedly many intersections
modulo \(P\) with the absolutely irreducible graph. Each \(z\)-residue
has at most one lift in the block, and then \(x\) has at most one
lift modulo \(P^k\). For larger absolutely irreducible factors,
the coefficient-uniform Bombieri–Pila bound in a box of side
\(X^{\max(a_0,h)}\) gives \(X^{O(\lambda^{3/2})+\varepsilon}\) points.
Non-absolutely irreducible rational factors contribute only boundedly
many points by intersecting their conjugates. The number of blocks
is \(O(1+X^{b-h})\).
[Bombieri–Pila, Theorem 5](https://people.maths.ox.ac.uk/pila/Ovals.pdf)

After freezing tiny remainders, fix the selected prime root on \(C\).
The graph joins the two remaining low cofactors. Color the high
prime factors on the free sides in disjoint sets at constant loss,
since their number is bounded in terms of the fixed thresholds.
For each edge the high-root ratio has a prescribed \(k\)-th power
modulo \(P^k\). Its numerator and denominator are \(X^{O(\lambda)}\),
so cross multiplication gives at most one exact ratio for each
modular root. It determines the coprime high pair. For a codegree,
eliminate the shared cofactor and apply the same lift to the
ratio of two such ratios. Prime colors factor it separately
on each side. At least one of the two ratios is nontrivial.
Freeze its common scalar at cost \(g_T\), then apply a residue
fibre to the other scalar and the unknown cofactor. Thus either
codegree exponent is at most
\[
e_C=(\max_{T\ne C} b_T-b_C)_++\max_{T\ne C}g_T+o(L).
\]
The [omitted-coordinate graph estimate](ABC_OMITTED_COORDINATE_GRAPH.md#box-bound) gives
\[
D\le b_C+\max\{\ell_A,\ell_B,
\min(\ell_A+\ell_B/2,\ell_B+\ell_A/2)+e_C/2\}+o(L),
\]
with \(A,B\) here denoting the two free letters. Fixed signs are
absorbed in the residue constants, so the modulus letter can
also be the original short addend.

### Checking all mass configurations

In this case analysis normalize all masses by \(L\), and omit
\(O(\eta)+o(1)\) where inequalities have fixed strict margins.
Write \(x=\ell_A\) for the short low mass and \(y,z\) for
the tall low masses. Thus \(S=2t+s\in[v,1/2+\delta]\),
\(s\le t\), and \(\ell_T\le b_T+2\delta\). If a tall prime
is available and \(s\) is sufficiently small, freeze the
entire short term and the tall prime root. The other tall
term is counted by a residue fibre, at total cost
\(t+L_A\le1/4+5\delta/2+3s/2<v\). Hence any short prime
modulus needed below has macroscopic mass.

* If neither tall root is prime, then \(y,z\le2\delta\)
  and \(x\ge v-4\delta\). The short root is prime and
  \(s\ge1/2-7\delta\). Use the short modulus. With \(y\le z\)
  its nonlinear score is at most
  \(s/2+3t/4+3y/4+z/2+\delta/2\le5S/12+3\delta<v\).
  The linear scores are at most \(s+2\delta<v\).

* If the short root is not prime and a tall one is, freeze
  the whole short term and that tall prime root, then use a
  residue fibre for the other tall term. The cost is at most
  \(t+s+x\le2S/3+2\delta\le1/3+8\delta/3<v\).

* Suppose the short root and at least one tall root are prime.
  For a tall modulus choose the other tall low mass as small
  as available. The linear term \(t+x\) is at most
  \(2S/3+2\delta<v\). When both tall roots are prime the
  other linear term is at most \((1-s-x)/2<v\).
  When that term comes from a composite tall root its low
  mass is at most \(2\delta\), so it also falls below \(v\).
  For the short modulus the linear terms are at most
  \(s+t+2\delta<v\).

  If all three roots are prime, put \(A_0=1-S\).
  The two nonlinear scores are bounded respectively by
  \[
  t+\min((3x+A_0)/4,A_0/2),
  \qquad (s+t)/2+3(A_0-x)/4 .
  \]
  Their minimum is at most \(1/2-(t+s)/4\le1/2-v/8<v\),
  by balancing the affine parts.

  If one tall root is composite, let \(y\le2\delta\) be
  its low mass and set \(e=(t-y)/2+\delta\). Then \(y<x\).
  The tall- and short-modulus nonlinear bounds are
  \(t+y+x/2+e/2\) and \((1-t+y-x)/2+e/2\).
  Their average is
  \((1+2t+2y)/4+\delta/2\le3/8+7\delta/4<v\).

All inequalities retain a positive margin for \(\eta\) sufficiently
small depending on \(\delta<1/22\). Let \(\delta\) tend to \(1/22\).

### Extension when only the two tall degrees agree

The \(5/11\) leading bound also holds for a box whose two tall
distinguished degrees are \(k\) and whose short distinguished
degree is \(j\ge k\), under the same \(O(\sqrt\lambda)\) bound
on all three nondistinguished heights. Use the notation of the
common-degree proof, with \(t=1/k\), short distinguished mass \(s\),
and low masses \(x,y,z\) on the short and the two tall letters
respectively. The short mass satisfies \(s\le1/j\le t\).

The geometry bound and Browning–Verzobio give
\(D\le S+o(L)\) and \(D\le L-S+o(L)\) with \(S=2t+s+o(L)\).
For Browning–Verzobio distinguish a tall exponent as \(r\) and the
other as \(p\), so its contribution is \(O(\lambda^{3/2})\)
despite the unequal third exponent. The whole-root residue lock
gives \(\ell_T\le b_T+2\delta L+o(L)\) on each letter.
If the short coordinate has weight \(js<2\lambda\), instead
use its \(O(\sqrt\lambda)\) total height and the two-tall strip
estimate to obtain the same inequality. In the cheap-divisor
split use \(\tau_T=2\lambda/k_T\) with \(k_T=k\) on the tall
letters and \(k_T=j\) on the short one. Its mass is \(O(\lambda^2)\)
and its amplified weight is \(2\lambda\), so the frozen remainders
add only \(O(\lambda)\) to the numeric low-cofactor heights.

Only the **short prime root** is used as a graph modulus in this
extension. Thus both free exponents in the graph are exactly
\(k\), and all lifting, ratio and codegree arguments are those
of the common-degree proof. The modulus is the short prime
base to the power \(j\). Its base mass is \(\rho=s+O(\lambda^2)\)
and its modulus weight is \(j\rho\ge k\rho\). The residue
determinant still uses the free degree \(k\), with blocks
\(h=\rho-C\sqrt\lambda/k\). Its two required inequalities
\(2a_0+k h<j\rho\) and \(a_0+k h<j\rho\) both follow as
before. No power associated to a frozen remainder is included
in the modulus. Those powers are fixed unit coefficients in the
free congruence. The fibre block overhead is
\((t-s)_++O(\lambda^{3/2})\).

The other option, whenever a tall prime is available, is to
freeze the entire short term and that tall prime root. The other
tall term is constrained by one inverse fibre. Its bound,
in units of \(L\), is
\[
F=t+s+x.
\]
This also handles \(s\) below a small fixed constant, since
\(x\le s+2\delta\) and \(t\le(1/2+\delta-s)/2\).
A short modulus used below thus has macroscopic base mass.

Take \(v=1/2-\delta\) with \(0<\delta<1/22\), and suppose
\(D\ge vL\). Normalize by \(L\) and suppress the small
threshold margins and \(o(1)\) corrections as in the previous proof.
If the short root is not prime, some tall root must be prime,
since otherwise all three low masses are at most \(2\delta\),
incompatible with \(1-S\ge v\). Then
\[
F\le 2S/3+2\delta\le 1/3+8\delta/3<v.
\]
If neither tall root is prime, the short graph gives exactly
the preceding bound \(5S/12+3\delta<v\), with all linear
scores below \(v\).

It remains that the short root and a tall root are prime.
All linear scores from the short-modulus graph are at most
\(s+t+2\delta\le1/3+8\delta/3<v\). Put \(A_0=1-S\).
When both tall roots are prime, its nonlinear score is at most
\[
G=(s+t)/2+3(A_0-x)/4.
\]
Taking the convex combination of the two available bounds
\[
\min(F,G)\le(3F+4G)/7
=\frac{5(s+t)+3A_0}{7}
\le\frac{3+S/3}{7}<v
\]
uses \(s\le t\), \(S\le1/2+\delta\), and \(\delta<1/22\).

If one tall letter is composite, let \(y\le 2\delta\) be its
low mass. Its common-scalar cost in a nontrivial ratio is at
most \((t-y)/2+\delta\). The codegree overhead also includes
\(t-s\). Since the other tall low mass is bigger than \(y\),
the resulting nonlinear score is at most
\[
G=(1-t+y-x)/2+(t-y)/4+\delta/2 .
\]
Now
\[
\min(F,G)\le(F+2G)/3
=1/3+t/6+s/3+y/6+\delta/3
\le5/12+5\delta/6<v .
\]
The last inequality again uses \(s\le t\) and \(\delta<1/22\).
All cases allow \(\delta\) to tend to \(1/22\).

### Refining the common-degree exponent to \(4/9\)

Keep the hypotheses and analytic bounds of the common-degree proof.
One can take any \(0<\delta<1/18\), with \(v=1/2-\delta\), to obtain
\(D\le(4/9+o(1))\lambda\) for these boxes. The fibre proofs,
the divisor-lock removal and the coloring argument are unchanged.
In particular every tall letter still has \(\theta_T>2\tau\):
indeed \(v/6>\delta\). Work in units of \(L\), with
\(S=2t+s\in[v,1/2+\delta]\), short low mass \(x\), and tall
low masses \(y,z\). Thus \(x\le s+2\delta\), \(y,z\le t+2\delta\).
A prime letter has low mass above \(2\delta\), with a small
fixed positive margin as in the original proof. A composite
letter has low mass at most \(2\delta\), up to the margins.
The cases below all have strict inequalities for \(\delta<1/18\).

If neither tall letter is prime, the original short-modulus
bound \(5S/12+3\delta\) still gives a strict saving, since
\(\delta<1/18\). The short root is prime and has macroscopic
mass because \(x\ge1/2-5\delta\) and \(s\ge1/2-7\delta\).

Suppose the short letter is not prime. Then \(x\le2\delta\),
and some tall letter is prime because \(1-S\ge v>6\delta\).
First try to freeze the short term and the tall prime.
Its cost is \(F=t+s+x\), and if \(F<v\) this already suffices.
Otherwise
\[
s\ge1/2-3\delta-2x,\qquad t\le2\delta+x.
\]
These inequalities also imply \((s+x)/2-\delta\) is positive
by a fixed margin, because \(s+x\ge v-t\) and
\(t\le1/4+\delta/2\). Thus the short composite root is subject
to the large-prime factor decomposition in the divisor argument.
Its common-scalar penalty is at most \((s-x)/2+\delta\).
Use a tall prime modulus for the graph.

If both tall roots are prime, choose the smaller remaining
tall low mass \(y\le(1-S-x)/2\). The nonlinear graph bound is
\[
t+x+y/2+(s-x)/4+\delta/2
\le 1/4+t/2+x/2+\delta/2
\le1/4+7\delta/2<v.
\]
The linear term \(t+x\) is below \(v\) by
\(t+x\le1/4+5\delta/2\). The other is
\(t+y\le(1-s-x)/2<v\), since
\(s+x\ge v-t\ge1/4-3\delta/2>2\delta\).

If the other tall root is composite, write \(y\le2\delta\) for
its low mass and use the larger of the two scalar penalties
\[
g_A=(s-x)/2+\delta,\qquad g_B=(t-y)/2+\delta.
\]
If \(g_B\) is larger, the nonlinear bound is at most
\(5t/4+3y/4+x/2+\delta/2\) when \(x\ge y\), and
\(5t/4+x+y/4+\delta/2\) when \(x\le y\).
Both are at most \(8\delta<v\), using \(t\le2\delta+x\)
and \(x,y\le2\delta\).
If \(g_A\) is larger, the nonlinear bound is at most
\[
t+s/4-x/4+\delta/2+\min(x+y/2,y+x/2)
\le1/8+21\delta/4<v,
\]
where \(t+s/4=S/2-s/4\) and
\(s\ge1/2-3\delta-2x\). Both linear terms are at most
\(1/4+5\delta/2<v\).

It remains that the short root and at least one tall root are
prime. The short root may be used as a modulus. If \(s\) is
sufficiently small, freeze the whole short term with a tall prime
as in the original proof, so \(s\) is macroscopic whenever
needed below. Take the taller low masses in the order
\(y\le z\), and for the tall modulus freeze the letter with
low mass \(z\). If there is a composite tall letter it is the
one of low mass \(y\le2\delta\), and the other is prime.
In either case \(x>2\delta\). The linear terms that can obstruct
the two modulus choices are \(t+x\) for the tall modulus and
\(s+z\) for the short modulus. Its other tall linear term
\(t+y\) satisfies \(t+y\le(1-s-x)/2<v\),
and the other short term is \(s+y\le s+z\).
The two possible obstructions cannot coexist: they would
give \(x+z\ge2v-t-s\), incompatible with \(x+z\le1-S\)
as \(t\ge v/3>2\delta\).

When all roots are prime and both modulus choices have
linear terms below \(v\), the preceding nonlinear comparison
gives \(1/2-v/8<v\). If \(t+x\ge v\), use the short modulus.
Its nonlinear score is at most
\[
3/4-t-s/4-3x/4\le3/4-7v/8<v.
\]
If \(s+z\ge v\), use the tall modulus, whose nonlinear score is
at most
\[
t+y+x/2\le1-v-t-x/2\le1/2-v/3<v.
\]

When the tall letter with low mass \(y\) is composite, set
\(e=(t-y)/2+\delta\). The tall and short nonlinear scores are
bounded by
\[
G_t=t+y+x/2+e/2,\qquad G_s=(1-t+y-x)/2+e/2.
\]
If both linear terms are safe, their average is at most
\(1/4+t/2+3\delta/2\), and \(G_s\le1/2-t/4\) since
\(x>2\delta\) and \(y\le2\delta\). Therefore their minimum
is at most \(5/12+\delta/2<v\), by balancing these two bounds.
If \(t+x\ge v\), substitute \(x\ge v-t\) into \(G_s\).
If \(s+z\ge v\), substitute \(x+y\le1-v-2t\) into \(G_t\).
In both cases the available modulus has nonlinear score at
most \(5/16+13\delta/8<v\). This exhausts the cases and proves
the refinement.

## A \(3/7\) small-budget bound with a short letter

Suppose the short letter has weighted height \(O(\sqrt\lambda)\).
On each of the two tall letters, keep one distinguished coordinate,
of degree \(k_B,k_C\), and suppose all its other coordinates
together have weighted height \(O(\sqrt\lambda)\). The degrees may
be unrelated. The exponent in such a box is at most
\((3/7+o(1))\lambda\) as \(\lambda\to0\).

Let \(L\le\lambda\) be total mass, \(L_A,L_B,L_C\) the letter
masses, and \(S=b_B+b_C\) the mass of the two tall roots.
The tall heights imply \(k_T b_T=1+O(\sqrt\lambda)\).
As in the half-budget proof, geometry gives \(D\le S\)
and the two-tall strip count gives \(D\le L-L_A-S+o(\lambda)\).

We also have \(D\le\min(L_B,L_C)+O(\lambda^{5/4})\). To see this, fix
one complete tall term \(N\asymp X\). The other is \(Q z^d\),
with root mass \(b=O(\lambda)\), whole-cofactor numeric height
at most \(X^{a_0}\) for \(a_0=O(\sqrt\lambda)\), and
\[
|Qz^d-N|\le X^{O(\sqrt\lambda)}.
\]
It is an \(X^{-1+O(\sqrt\lambda)}\) relative approximation
to the real monomial \(Q=N z^{-d}\). For a dyadic root scale
\(Z\asymp X^b\), normalize the
coordinates to \(u=z/Z\) and \(v=QZ^d/N\). The normalized
points satisfy \(v=u^{-d}(1+O(X^{-1+O(\sqrt\lambda)}))\)
on a fixed compact positive \(u\)-interval. Each integer
\(z\) determines at most one \(Q\) for sufficiently small
\(\lambda\), since consecutive \(Qz^d\) values are separated
by order \(X^{1-O(\sqrt\lambda)}\). The divided-difference
determinant used for the narrow strip counts
\(X^{O(\lambda^{5/4})+\varepsilon}\) such pairs. In detail, put
\(n\asymp b^{-1}\) with \(nb<1/4\) and use the \(n\) smallest
weights \(i b+j a_0\) of \(z^iQ^j\), taking a positive upper
bound \(a_0\asymp\sqrt\lambda\) if the actual cofactor is smaller.
The largest first monomial exponent is \(O(\sqrt{a_0}/b)<d\),
so the powers \(i-dj\) on the exact graph are distinct.
A common shift makes the powers nonnegative for the Schur
factorization. The determinant divided by the Vandermonde
in \(u\) stays nonzero since the perturbation in a
row of order \(r\) is \(O(X^{-1+O(\sqrt\lambda)+rb})\).
The sum of selected weights is \(O(\sqrt{b a_0}\,n^{3/2})\).
On each normalized root interval of length \(X^{-C b\sqrt{a_0}}\),
this weight sum is smaller than the saving from the
Vandermonde for large fixed \(C\). Restoring the monomial
scales therefore makes the nonzero determinant an integer
of absolute value less than one. Each interval contains
at most \(n-1\) points.
This proves the asserted fibre count. Restoring whichever
tall term was fixed gives the claimed minimum bound.

Assume \(D\ge D_0=(3/7+\epsilon_0)L\) for fixed
\(\epsilon_0>0\), and restrict to \(L\asymp\lambda\), as is
allowed by the pair bounds. Write \(C\) for the tall letter
with the smaller root mass, \(b_C\le b_B\). Set
\[
\tau=\frac{2\lambda}{k_C},\qquad
\theta=D_0-\frac{L-L_C}{2}-\eta L
\]
for a sufficiently small fixed \(\eta>0\).
The minimum bound above gives \(L_C\ge D_0-o(L)\), so
\(\theta\gg L\) and \(\tau=O(\lambda^2)\). The whole-root lock
already gives the target if \(b_C\le\theta\), so the root has
mass above \(\tau\) in the remaining argument.
Discard the subfamily in which this root has a divisor of
mass in \([\tau,\theta]\), by freezing that divisor as a
degree-\(k_C\) modulus and using the higher-moment energy.
Its exponent is at most \(D_0-\eta L\).
On the remaining subfamily there is a prime factor \(P\)
of mass \(\rho>\theta\). Otherwise every prime factor is
smaller than \(X^\tau\), and a product first crossing
\(X^\tau\) would be a forbidden divisor.

Freeze the complete short term and \(P\), at cost
\(L_A+\rho\). The free tall cofactor and root satisfy one
inverse residue fibre modulo \(P^{k_C}\). Its block exponent
can be taken as \(h=\rho-K\sqrt\lambda/k_B\), giving
\(b_B-h=b_B-\rho+O(\lambda^{3/2})\). Indeed
\(k_C/k_B\ge1-O(\sqrt\lambda)\) by the tall weights and
\(b_C\le b_B\), so \(2a_0+k_B h<k_C\rho\) for sufficiently
large \(K\). Also \(h<\rho\), as required to lift the
low-degree components uniquely. The short and free tall
terms determine the remaining tall term. Thus the cost of
freezing \(P\) cancels against the root-block count,
giving
\[
D\le L_A+\max(b_B,b_C)+o(L).
\]

Finally the whole-root locks on both tall letters give
\(b_T\ge D-(L-L_T)/2\). Since \(L_T\ge D-o(L)\),
\(\min(b_B,b_C)\ge(3D-L)/2-o(L)\).
Together with the strip bound
\(S\le L-L_A-D+o(L)\), this gives
\[
D\le L_A+S-\min(b_B,b_C)+o(L)
 \le(3L-5D)/2+o(L).
\]
Hence \(D\le3L/7+o(L)\), as claimed.

## Recovering a primitive ratio from a dominant divisor

In an anatomic box at small radical budget \(\lambda\), suppose
each letter has one chosen coordinate of degree at least two
and complementary weighted height \(O(\sqrt\lambda)\).
Write \(\beta_a,\beta_b,\beta_c\) for the masses of the
chosen coordinates. Order the two addends so that \(a\le b\).
Then the count exponent satisfies
\[
 D\le\beta_a+\min(\beta_b,\beta_c)+O(\lambda^{7/6}).
 \tag{18}
\]
In particular these boxes have
\(D\le(6/13)\lambda+O(\lambda^{7/6})\).
If the three complete logarithmic heights are bounded below
by a fixed positive constant, the bound improves to
\(D\le(2/5)\lambda+O(\lambda^{7/6})\).

### Fixing two roots

Use labels \(F,J,M\) for a fixed root, a free root and a modulus
root, and write the terms as \(Ax^p,By^q,Cz^r\), respectively.
The equation is
\(\epsilon_F Ax^p+\epsilon_J By^q+\epsilon_M Cz^r=0\)
with each \(\epsilon_T\) equal to \(1\) or \(-1\).
Fix \(x,z\) and take \(n=z^r\). The remaining coefficients
and the complete modulus letter may vary. By primitivity,
\(x,y,A,B\) are units modulo \(n\). The equation gives
\[
 VT^q\equiv -(\epsilon_F/\epsilon_J)x^p\pmod n,\qquad
 T=y,\qquad V=B/A .
\]
Each pair \((T,V)\) determines the exact reduced ratio of the
\(F,J\) letters, \(x^p/(VT^q)\). Since the two complete values
are coprime, its reduced numerator and denominator are those
values. The signed equation then determines the third value.
This argument counts primitive triples represented by the box,
irrespective of the number of their anatomic factorizations.

Suppose \(J,M\) have complete logarithmic heights at least a
fixed \(b>0\) in a dyadic shell \(Y/2<c\le Y\).
For sufficiently small \(\lambda\), and sufficiently large \(Y\)
at a fixed anatomic cutoff, the complementary-height assumption
gives \(z^r\ge Y^{b/2}\) and \(q,r\ge b/(4\lambda)\).
Indeed, the chosen coordinate carries height at least \(b/2\)
and its mass is at most \(2\lambda\), allowing for the
decomposition losses. With \(k=\lfloor b/(4\lambda)\rfloor\)
we have \(k\le r,\ q\ge k\) and \(H(T)\le Y^{2\lambda}\le n^{2/k}\).
The \(B/A\) ratio has height at most \(Y^{C_1\sqrt\lambda}\)
for a fixed \(C_1\). This includes the small exceptional
coefficients and bounded factors from the dyadic endpoints.
Take \(\ell=\lfloor c_1\lambda^{-1/6}\rfloor\), with \(c_1>0\)
small enough in terms of \(b,C_1\).
Then \(1024\ell^2\le k\) and
\(H(V)\le n^{1/(16384\ell^3)}\).
Also \(q>h=\lfloor k/(256\ell^2)\rfloor\).

The independent-exponent version of the scalable congruence
bound gives \(O_{k,\ell,q,\varepsilon}(Y^{1/(\ell k)+\varepsilon})\)
possibilities for each fixed \(x,z\). Since
\(1/(\ell k)=O(\lambda^{7/6})\), the box exponent is at most
\(\beta_F+\beta_M+O(\lambda^{7/6})\).
In a finite box the degrees are bounded at each fixed
decomposition loss, so all congruence constants can be
maximized over the finite set of degrees. Coefficient
choices, dyadic powers of \(\log Y\), and approximations to the
weighted heights are absorbed by taking the decomposition loss
sufficiently small in terms of \(\lambda,\varepsilon\).
The constants in the exponent may depend on \(b\) and the
fixed complementary-height bound.

Take \(F=a\) and use either tall letter as \(M\), leaving the
other as \(J\), to obtain (18). The surface bound after
fixing the complementary coordinates gives
\[
 D\le L-\beta_a-\beta_b-\beta_c+O(\lambda^{3/2}).
\]
The residue locks (A) on the two tall letters give
\[
 \beta_b+\beta_c
 \ge 2D-L+(L_b+L_c)/2\ge 5D/2-L,
\]
where the last inequality fixes the two tall letters.
Adding (18) to the surface bound yields
\(2D\le L-(\beta_b+\beta_c)/2+O(\lambda^{7/6})\),
which gives the \(6/13\) bound.
When all heights are bounded below, the same congruence
counts give \(D\le\beta_i+\beta_j+O(\lambda^{7/6})\)
for every pair of distinct letters \(i,j\).
Adding the three inequalities to twice the surface bound
gives \(5D\le2L+O(\lambda^{7/6})\).
[BBLT, Propositions 2.1 and 4.1](https://arxiv.org/html/2410.12234v2)
[Browning--Verzobio, Theorem 1.2](https://arxiv.org/html/2608.24512v1#S1)

The separated-degree profiles above have height one on each
letter and are cut by the \(2/5\) bound. On a comparable-height
profile with two single-heavy letters and one double-heavy
letter, use a single-heavy letter as \(J\), the other as \(F\),
and lock the larger-height heavy coordinate on \(M\).
At the near-half masses the cost is
\((L-L_J)/2+o(\lambda)\), which is below half the budget
for \(L_J\asymp\lambda\).
With just one single-heavy letter available as \(J\), fixing
both heavy coordinates on \(F\) can cost \(L_F+o(\lambda)\).
In particular the mixed case with \(L_F,L_J,L_M\sim\lambda/3\)
can still have cost \(\lambda/2+o(\lambda)\).

## Three-variable congruences at a prime power

In a near-half box with two heavy coordinates on each letter,
the complementary masses are \(o(\lambda)\). After fixing the
low cofactors, one root on the first letter and one on the
modulus letter, the three remaining heavy variables give a
congruence \(U^pV^qW^{-r}\equiv K\) with three large exponents.
The sign and the fixed factors enter the unit \(K\).
The following tests treat positive rational units over
\(n=P^e\), with \(e\ge k\), for a prime \(P\).

### An exact curve with many points in one class

Choose \(p=k\), \(q=k+2\), \(r=k+1\) and \(K=1\).
Since \(p+q=2r\), the points
\[
 (U,V,W)=(t,t,t^2)
\]
satisfy the equation exactly. Take \(P>k+2\) and \(e=ka\)
with \(a\ge1\). For \(P^a/2< b<P^a\) and \(P\nmid b\), set
\(t=(2b-P^a)/b\). These are distinct reduced positive fractions
of height at most \(P^a=n^{1/k}\), all congruent to \(2\bmod P^a\).
Their triples have all heights at most \(n^{2/k}\), are units,
and lie in the single class \((2,2,4)\bmod P^a\).
There are \(\gg P^a=n^{1/k}\) such points.
Thus a uniform \(n^{o(1/k)+\varepsilon}\) count for the unit
surface would fail even at coprime, distinct exponents and
on one residue class. The exact curve supplies the obstruction
to nonvanishing of arbitrary large evaluation matrices.
The anatomic application has additional coprimality constraints
on its heavy roots. For pairwise coprime positive integer roots,
fixing the exact rational value instead gives divisor-many
solutions: its denominator determines \(W^r\), and the primes
of its numerator can be assigned to \(U\) or \(V\).

### Auxiliary surfaces and a transverse curve reduction

Fix a height constant \(C>0\). Let all three heights be at
most \(n^{C/k}\) and assume \(P\nmid r\).
For all sufficiently large \(k\) in terms of \(C\), and large
\(n\) for fixed parameters, partition by the three residues
modulo \(P^a\), \(a=\lfloor e/k\rfloor\).
There are at most \(rP^{2a}\) admissible classes: the two
residues \(U,V\) determine at most \(r\) lifts of \(W\)
on this modulus, by the unit derivative of \(K W^r\).
In each class the rational points lie on an auxiliary
polynomial of degree bounded in terms of \(C\).
On a factor of that polynomial whose reduction meets the
modular surface transversely at the class, they lie on
algebraic curves of degree bounded in terms of \(C\).

For the first assertion, use the expression
\(f=U^pV^q-K W^r\) with a unit lift of \(K\) in \(\mathbb Z_P\).
Fix integral lifts \(u_0,v_0\) of the two residues and lift
the \(W\)-residue by Hensel's lemma to an exact root of
\(f(u_0,v_0,W)\). The root is an integral power series in
\(U-u_0,V-v_0\) on the residue class: the coefficients can
be found recursively by dividing by the unit derivative.
The value at each pair \(U,V\) agrees with the given \(W\)
modulo \(P^e\). Thus for each monomial of total degree at
most \(M\) in the three variables its evaluation column
agrees modulo \(P^e\) with an integral power series in
two variables, both divisible by \(P^a\).

Put \(D=\binom{M+3}{3}\). Let \(d_1,\ldots,d_D\) be the
least total degrees in the two-variable monomial list,
counting multiplicities, and write \(S_2(D)=\sum_j d_j\).
The determinant of \(D\) evaluations has valuation at least
\(aS_2(D)\) if \(ad_D\le e\).
In fact, Taylor expansion and determinant multilinearity
give weights \(a(i+j)\) for distinct monomial columns.
Using \(s\) error columns costs \(se\), enough to replace
the \(s\) largest among the first \(D\) weights.
On clearing rational denominators the integer determinant
is at most \(D!n^{3CMD/k}\), and its divisibility is at least
\(n^{S_2(D)/(2k)}\). Since \(S_2(D)\asymp D^{3/2}\),
choose \(M\) in terms of \(C\) so that
\(S_2(D)/2>3CMD\), then take \(k\ge d_D\).
All such determinants vanish for large \(n\).
The rational evaluation matrix therefore has a nontrivial
kernel, giving the auxiliary polynomial.

For the second assertion, fix an irreducible factor \(G\)
over \(\mathbb Q\), with primitive integer coefficients.
Suppose the gradients of \(G,f\) are independent modulo
\(P\) at the class. Choose a unit two-by-two Jacobian minor
and lift the other two coordinates as power series in
the third. The multivariate Hensel recursion gives integral
coefficients and uniqueness modulo \(P^e\), since the
evaluations satisfy \(G=0\) exactly and \(f\equiv0\bmod P^e\).
Choose two variables omitting a variable on which \(G\)
depends, and use all \(D'=\binom{M'+2}{2}\) monomials of
total degree at most \(M'\) in those two variables.
They are now power series in a single variable divisible
by \(P^a\). For \(k\ge D'-1\), the same argument gives
valuation at least \(a\binom{D'}2\).
The cleared determinant is bounded by
\(D'!n^{2CM'D'/k}\). Choose \(M'\) with
\((D'-1)/4>2CM'\) to force vanishing for large \(n\).
This gives a nonzero polynomial independent of the omitted
variable, so it is coprime to \(G\).
Their common rational points lie on curves or points of
bounded degree. Rank-deficient classes and the arithmetic
counts on the resulting curves remain to be controlled
for the two-heavy boxes.
