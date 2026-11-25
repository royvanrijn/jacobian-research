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

For the small-budget limit, the graph argument gives
\[
N_\lambda(X)\ll_{\lambda,\eta,\varepsilon}X^{(1/3+\eta)\lambda+\varepsilon}
\]
for each fixed \(\eta>0\) and sufficiently small fixed \(\lambda>0\). It combines counts of signed monomials modulo full integers with counts in real short strips. Both estimates allow a bounded number of anatomic factors with unbounded degrees in the limit. When the normalized count exponent exceeds \(1/3\), the mass split forces a saving from one of the two counts. [Small-budget proof](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)

A radical bottleneck on one letter gives a stronger count at small
budgets. If that letter's radical is at most \(c^{\vartheta\lambda}\),
the limiting coefficient is at most
\(\min(1/3,(1+\vartheta)/4)\). The box inequality underlying this
bound is \(4d\le1+\min_i l_i\).
[The quarter inequality and its count consequence](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-quarter-inequality-for-a-letter-of-small-radical-mass)

If \(l_C<d\), the diffuse mass on each opposite letter is at least
\(d-c_C\), where \(c_C\) is the nonatomic mass on \(C\). At equality
in the quarter bound this excludes an opposite letter whose numerical
height lies on finitely many scaled degrees. Applied to powerful
sums, it gives an exponent below \(1/(2p)\) at small-radical targets
when one summand is \(p\)-full and the other has radical at most
\(N^{1/p}\), including representations with common divisors.
[The diffuse mass inequality](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-global-inequality-for-a-low-mass-letter-and-diffuse-degrees)
· [The one-full transfer](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#one-full-summand-and-a-small-target-radical)

There are absolute \(h_*,\sigma>0\) such that, for sufficiently small fixed \(\lambda\), the subcount with \(\min(a,b)\le c^{h_*}\) is \(O_{\lambda,\varepsilon}(X^{(1/3-\sigma)\lambda+\varepsilon})\). The proof isolates a movable interval of scaled valuation degrees and uses a full modulus on the other tall term. [Short-height sector bound](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-strict-saving-when-the-short-term-is-small)

At normalized exponent \(1/3\), if the total nonatomic degree mass
is \(1/2\), all three numerical heights equal \(1\) and all three
limiting degree measures have bounded scaled degrees and first
moment \(1\). [Variable-mesh proof](ABC_VARIABLE_MESH_DEGREE_GROUPING.md#consequences-at-the-minimum-total-nonatomic-mass)

At the balanced endpoint, residual upper quantiles constrain the
mass that can be locked on the three letters. They include escaped
numerical weight. Together with a selective congruence modulus,
they force positive scaled-degree mass above an explicit cutoff
when the largest opposite-coordinate masses are small.
[Residual upper quantiles](ABC_RESIDUAL_QUANTILE_SPECTRUM.md)

For every actual balanced endpoint, an affordable pair of full
divisors reduces the count to distinct targets over a fixed record.
In each sufficiently short numerical interval such a fibre contributes
at arbitrarily small normalized exponent. This forces almost the
entire conditional target entropy to lie in the choice of numerical
interval, for any fixed positive interval exponent small enough for
the chosen loss. The affine rank bound extends fixed-record sparsity to every
balanced profile in short windows centred at linear functions of
boundedly many small integer parameters. This includes ordinary
and polynomial progressions with arbitrary coefficients.
Counting separated targets, or covering them by sufficiently
few such systems, remains open. A separate digit construction
shows that progression sparsity alone cannot give a power
saving on additive energy.
[Fixed-record interval entropy](ABC_WIDE_COFACTOR_CONGRUENCES.md#numerical-interval-entropy)
· [Affine rank and the additive obstruction](ABC_AFFORDABLE_FIBRE_AFFINE_RANK.md)

There is a reverse powerful-number transfer. Repair only \(c\)
to order \(M\), by its least possible common multiplier.
The new height is \(\prod_{p\mid c}p^{\max(M,v_p(c))}\).
For \(M\) a sufficiently large fixed multiple of \(1/\lambda\),
a uniform exponent saving below \(1/M\) in the corresponding
minimal-repair count on the balanced radical sector is equivalent
to the desired uniform small-budget saving. Unrestricted
multipliers always supply a leading \(1/M\) exponent from a
fixed quotient, so the transfer requires control of the
minimally repaired quotients.
[The minimally repaired quotient exponent](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#the-exponent-of-the-minimally-repaired-primitive-count)
· [Reverse transfer](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

Across overlapping affordable records, primitive complete
heights and the two free radical caps alone allow simultaneous
saturation of every fixed slot-subset count, even with
minimal additive energy in the record fibres. The
small-budget construction has unrestricted target radical.
[Overlapping record counterexamples](ABC_OVERLAPPING_FULL_RECORD_COUNTEREXAMPLES.md)

Pointwise bounds for the size of \(c\) are a separate route to finiteness. With \(R=\operatorname{rad}(abc)\), Stewart–Yu give \(\log c\ll R^{1/3}(\log R)^3\). Pasten's Theorem 1.4 gives \(\log c\le \eta^{-1}\exp(C\sqrt{\log R\log\log R})\) when \(a\le c^{1-\eta}\), and the same expression with \(\min\{P^+(a),P^+(b),P^+(c)\}\) in place of \(\eta^{-1}\) in the second part. Soulé's November preprint claims \(c\ll_\varepsilon R^{24+\varepsilon}\) in Theorem B, which would give finiteness for fixed \(\lambda<1/24\). In its genus-six fibration the degree used to prove this is \(-13k/2\) in place of the claimed lower bound \(-3k\). The corrected coefficient leaves the proof of the pointwise bound open. [Torsion calculation](ABC_SOULE_TORSION_AUDIT.md#torsion-lengths-and-the-degree-calculation) [Stewart–Yu, Theorem 1](https://uwaterloo.ca/pure-mathematics/sites/default/files/uploads/documents/s0012-7094-01-10815-6.pdf) · [Pasten, Theorem 1.4](https://arxiv.org/html/2312.03566v1) · [Soulé, Theorem B](https://arxiv.org/html/2611.05217v1)

At the endpoint \(\lambda=1\), Bright gives infinitely many triples with \(c/R\) growing at least as \(\exp(C\sqrt{\log c}/\log\log c)\) for an explicit positive constant \(C\). For any fixed \(0<\lambda<1\) this lower bound on \(c/R\) has a subpower scale in \(c\). [Bright, Theorem 3.1](https://cs.curtisbright.com/reports/abc-lower.pdf)

Carella's v3 preprint claims infinitely many strict exceptions, but the identity used in equation (4.35) already fails at a prime square. For the weight \(z^{\omega(n)}\), the Euler factor at \(p\le y\) is \(1+zp^{-s}/(1-p^{-s})\), whose coefficient at \(p^{-2s}\) is \(z\). The displayed factor \((1-p^{-s})^{-z}\) has coefficient \(z(z+1)/2\). They differ for \(z\ne0,1\), including the varying values near \(1\) used in the Perron argument. The claimed short-interval estimate in (4.40) requires a proof for the correct generating function. [Carella v3, §4.3 and equation (4.35)](https://arxiv.org/pdf/2608.16764v3)

### A radical-count contradiction to the claimed smooth-number lower bound

For any \(X,Y\ge1\) and fixed \(\varepsilon>0\),
\[
 \#\{n\le X:\operatorname{rad}(n)\le Y\}
     \ll_\varepsilon Y X^\varepsilon .
 \tag{CR}
\]
Indeed, fix the radical \(r\), write \(n=rd\), and use
\[
 \#\{d\le X/r:p\mid d\Rightarrow p\mid r\}
 \le (X/r)^\varepsilon\prod_{p\mid r}(1-p^{-\varepsilon})^{-1}
 \ll_\varepsilon X^\varepsilon .
\]
The large-prime factors are bounded by \(p^\varepsilon\)
whenever \(p^\varepsilon\ge2\). The remaining primes contribute
a constant, and summing over \(r\le Y\) proves (CR).

In Carella v3 Corollary 4.1, set \(L_0=\log\log x\). Each
\(x^{1/(7L_0)}\)-smooth integer with at most \(L_0+3\sqrt{L_0}\)
distinct prime factors has radical at most \(x^{1/7+3/(7\sqrt{L_0})}\).
The count in the entire interval \((x,2x]\) is therefore
\(O_\varepsilon(x^{1/7+\varepsilon+o(1)})\), already too small
for the asserted \(x^{3/5}\rho_D(7L_0)\) lower bound in (4.42)
on the shorter interval \((x,x+x^{3/5}]\). Here \(\rho_D\)
is the Dickman function. Its delay equation and monotonicity give
\(u\rho_D(u)=\int_{u-1}^u\rho_D(t)\,dt\ge
\tfrac12\rho_D(u-\tfrac12)\) for \(u\ge1\). Iterating gives
\(\rho_D(u)\ge(2u)^{-\lceil2u\rceil}\), so
\(\rho_D(7L_0)=x^{-o(1)}\). This directly contradicts (4.42),
which the strict-exception argument uses.
[Carella v3, Corollary 4.1 and §5](https://arxiv.org/pdf/2608.16764v3)

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
 N_{9/10}(X)&\ll_\varepsilon X^{1/2+\varepsilon},&
 N_1(X)&\ll_\varepsilon X^{1070000267/2000000686+\varepsilon}.
\end{aligned}
\]
These bounds also apply to smaller budgets by inclusion of the sets.
The [combined overlay](ABC_PAIR_DEVIATION_CERTIFICATE_AUDIT.md#a-strict-improvement-through-unit-radical-budget)
gives the displayed unit exponent and improves the earlier unit tree
on the entire interval from \(991/1000\) through one.
Reweighting the tree leaves gives \(N_{3/4}(X)\ll_\varepsilon
X^{1773/4012+\varepsilon}\).

The [sparse retargeted certificates](ABC_PAIR_DEVIATION_CERTIFICATE_AUDIT.md#sparse-retargeting-above-budget-910)
give these additional exponents. Each exponent applies through the listed budget.


| Radical budget | Count exponent |
| --- | --- |
| \(91/100\) | \(51/100\) |
| \(23/25\) | \(257500283/500000550\) |
| \(37/40\) | \(34666673/66666680\) |
| \(93/100\) | \(21/40\) |
| \(187/200\) | \(211000210979/400000400800\) |
| \(47/50\) | \(53/100\) |
| \(189/200\) | \(213/400\) |
| \(19/20\) | \(854400813351/1600001523200\) |
| \(191/200\) | \(1069/2000\) |
| \(24/25\) | \(106980070711/200000132200\) |
| \(193/200\) | \(52194191487/97561060000\) |
| \(482503/500000\) | \(133749764986241/250000000750000\) |
| \(387/400\) | \(38861/72800\) |
| \(97/100\) | \(267/500\) |
| \(970691/1000000\) | \(2606510113/4872000000\) |
| \(389/400\) | \(106950853/200000400\) |
| \(39/40\) | \(7132667863/13333337200\) |
| \(60941/62500\) | \(835936493193/1562500009375\) |
| \(391/400\) | \(855199842919/1600000001600\) |
| \(48877/50000\) | \(53499982624411/100000000300000\) |
| \(979/1000\) | \(2137999607311/4000000004000\) |
| \(49/50\) | \(1069/2000\) |
| \(980011/1000000\) | \(534529057150099/1000000005000000\) |
| \(1961/2000\) | \(3566601579/6666670000\) |
| \(981/1000\) | \(1069/2000\) |
| \(544721207/554858000\) | \(1069/2000\) |
| \(393/400\) | \(213885997/400001200\) |
| \(196501/200000\) | \(171793213/321280000\) |
| \(983/1000\) | \(534899976623/1000000003000\) |
| \(983003/1000000\) | \(5802594141/10847000000\) |
| \(123/125\) | \(8916501783/16666670000\) |
| \(12300001/12500000\) | \(835921978968881/1562500009375000\) |
| \(197/200\) | \(4279923132521/8000005856000\) |
| \(985001/1000000\) | \(534991037592609/1000000003000000\) |
| \(493/500\) | \(106998000419/200000000800\) |
| \(986001/1000000\) | \(267495625400507/500000002000000\) |
| \(987/1000\) | \(534990173281/1000000324000\) |
| \(4935001/5000000\) | \(2431773951451/4545455000000\) |
| \(247/250\) | \(5349904295969/10000008030000\) |
| \(9880003/10000000\) | \(53499028958498837/100000000500000000\) |
| \(989/1000\) | \(21399616776161/40000031360000\) |
| \(4945001/5000000\) | \(441830515859/825857500000\) |
| \(99/100\) | \(26749499961979/50000000001500\) |
| \(9900001/10000000\) | \(213996134798720789/400000197200000000\) |
| \(991/1000\) | \(267494705491/500000001000\) |
| \(9910001/10000000\) | \(891649018497667/1666666670000000\) |
| \(57804903499/58300100000\) | \(53499/100000\) |
| \(199/200\) | \(17833135833/33333333400\) |
| \(9999/10000\) | \(2674998941999/5000000010000\) |


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
bounded degree. The next argument treats classes of arbitrary rank and counts
the coprime integer roots on the auxiliary curves.

### Auxiliary curves in every residue class

Fix \(C>0\) and positive integer exponents \(p,q,r\). For all
\(k\ge k_0(C)\) and \(k\)-full \(n\), consider positive pairwise
coprime integer roots \(U,V,W\le n^{C/k}\), all units modulo \(n\),
with
\[
 U^pV^qW^{-r}\equiv K\pmod n .
\]
The count is \(O_{C,k,p,q,r,\varepsilon}(n^{2/k+\varepsilon})\),
uniformly in the unit \(K\). The proof below gives the same bound
over all \(k\)-full moduli, including powers of primes of bad
reduction of the auxiliary curves. For any fixed amplification
order \(\ell\), if \(k\ge k_0(\ell C)\), this yields
\(O_{C,k,p,q,r,\ell,\varepsilon}(n^{2/(\ell k)+\varepsilon})\).

First divide the exponent vector by \(g=\gcd(p,q,r)\). The
primitive monomial must lie in one of \(O_g(1)^{\omega(n)}\)
unit residue classes, since the kernel of the \(g\)-th power map
on each prime-power unit group has size at most \(2g\).
It suffices to work in one such class, with primitive exponent
vector \(v=(p,q,-r)\) and target \(K\). We allow constants depending
on the fixed exponents, but the threshold for \(k\) will depend
only on \(C\).

Choose the first auxiliary degree \(M=M(C)\) as in the
two-parameter determinant above. If \(\|v\|_1\le2M\), the
monomial has rational height at most \(n^{2CM/k}\).
Two values in the same residue cannot differ when \(k>4CM\),
by clearing their reduced denominators. Its exact denominator
determines \(W^r\). On its numerator, each prime is assigned
to \(U\) or \(V\) and its exponent then determines the root
valuation. This gives \(n^\varepsilon\) possibilities.
We may thus assume \(\|v\|_1>2M\).

At each prime \(P\mid n\) permute the variables so the exponent
of the last variable is a \(P\)-adic unit. The two-variable
description of the modular surface then holds on every unit
residue class by Hensel lifting. For each of the three choices
of last variable, select a direction \((1,B)\) with \(B\in\mathbb Z\)
so that the restrictions of all \(S=\binom{M+3}{3}\) monomials
of degree at most \(M\) have distinct exponents \(\nu_j\).
In the original order, with last exponent \(-r\), the exponent
for \(U^iV^jW^l\) is \(i+jl_0+l(p+ql_0)/r\) with \(l_0=B\).
The projected linear weights can coincide identically only if
two exponent triples differ by an integer multiple of \(v\).
Their difference has norm at most \(2M\), which excludes
this case. Each pair excludes at most one \(B\).

Write \(e=e_P\) and \(a=\lfloor e/k\rfloor\). In a class at
depth at least \(a\), let \(x,y\) be the relative deviations
of the two free coordinates from their unit centers. The
lifted third coordinate has the form
\(w_0(1+x)^{p/r}(1+y)^{q/r}\) in the original order.
The binomial series is integral because the two ratios
belong to \(\mathbb Z_P\). The exact lift agrees with
the third coordinate modulo \(P^e\).
The Taylor coefficient matrix of orders \(0,\ldots,S-1\)
on \(x=t,\ y=(1+t)^B-1\) is the binomial Vandermonde
\[
 \det\bigl(\tbinom{\nu_j}{i}\bigr)_{0\le i<S,\,1\le j\le S}.
\]
It is a nonzero rational of nonnegative \(P\)-adic valuation
\(b_P\), with \(b_P=0\) outside a finite set determined by
the exponents and \(M\). Center factors multiply its columns
by units. For any primitive integer polynomial \(G\) of
degree at most \(M\), its restriction therefore has a
coefficient along this direction of order at most \(S-1\)
and valuation at most \(b_P\), by the adjugate identity.

If \(b_P=0\), partition at depth \(a\). In coordinates
\(x=t,\ 1+y=(1+s)(1+t)^B\) the restriction is regular in
\(t\) of some order \(m\le S-1\) modulo \((P,s)\).
If \(b_P>0\), put \(L_P=b_P+1\), partition at depth \(a+L_P\)
and scale \(x,y\) by \(P^{L_P}\). Some original coefficient
of total degree at most \(S-1\) has valuation at most \(b_P\),
since substitution of the direction is integral and preserves
the order filtration. After scaling, the minimum coefficient
valuation \(b'_P\) is at most \(b_P+(S-1)L_P\). All terms of
degree at least \(S\) have strictly larger valuation.
Divide by \(P^{b'_P}\) and choose a direction over a finite
unramified extension on which the lowest nonzero homogeneous
part modulo \(P\) does not vanish. An invertible integral
linear change with this direction as the \(t\)-axis gives
regularity in \(t\), again of order at most \(S-1\). The
parameters \(s,t\) lie in \(P^a\) in the extension and the
valuation of rational determinants is unchanged. Set \(b'_P=0\)
in the first case. We may drop prime powers at the finite
exceptional primes for which \(e<2(b_P+(S-1)L_P)\).
The dropped product and all refinements of the partitions
cost a constant in the fixed parameters. Write \(n_1\) for
the retained part, so \(\log n_1=\log n+O(1)\).
The first and second determinants both use these partitions.

There are \(O_{p,q,r,M}(1)^{\omega(n)} n^{2/k}\) classes, since
each pair of free residues determines at most
\(\max(p,q,r)\) third residues at each prime. The first
determinant gives an auxiliary \(G\) for each class as above:
sum the valuation lower bounds with weights \(\log P\) and
use \(a\ge e/(2k)\) and \(n_1\gg n\). We can choose an
integer polynomial with coefficients \(\ll n^{O_C(1/k)}\).
Indeed, the rows with cleared denominators have entries
bounded by \(n^{3CM/k}\); a nonzero minor of maximal rank
and the adjacent cofactors give a nonzero null vector.
Pass to primitive irreducible factors of \(G\), whose
coefficient bounds have the same form. For example, a
Kronecker substitution followed by a fixed-degree univariate
factor height bound suffices.

Fix such a factor and work at each retained prime as above.
Its restriction to the lifted surface vanishes at the
projected points modulo \(P^{e-b'_P}\) after normalization.
Order zero has no such points. For positive order \(m\),
formal Weierstrass preparation and division over the
complete local ring \(\mathcal O[[s]]\) express any evaluation
column, modulo \(P^{e-b'_P}\), as an integral series in
\[
 s^i t^j,\qquad i\ge0,\quad 0\le j<m.
\]
In detail, the restricted equation is a unit times a monic
Weierstrass polynomial of degree \(m\) in \(t\). At each
point the unit stays a unit. The remainder of division has
\(t\)-degree below \(m\), and the multiple of the polynomial
vanishes to the stated precision. All series converge since
\(s,t\in P^a\). This uses regularity of a bounded Taylor order
and makes no first-derivative assumption.

Use \(J=\binom{M'+2}{2}\) monomials in two ambient variables,
omitting a variable on which the factor depends. The sum
of the \(J\) least total degrees of the remainder monomials
is at least
\[
 \sum_{j=0}^{J-1}\lfloor j/m\rfloor\ \ge\ \frac{J(J-m)}{2m}.
\]
For \(k\ge 2J\), each of the \(J\) least weights is at most
\(a(J-1)\le e/2\), so each error column costing \(e-b'_P\)
can replace a column of any such weight. Multilinearity
therefore gives valuation at least \(aJ(J-m)/(2m)\)
at every retained prime. Take \(M'\) large in terms of \(C,S\)
so \((J-S)/(4S)>2CM'\). Then the resulting divisibility
exceeds the cleared real bound \(J!n^{2CM'J/k}\) for large
\(n\). The resulting auxiliary \(H\) has coefficients
\(\ll n^{O_C(1/k)}\) by the same null-vector argument.
It is independent of the omitted variable and hence coprime
to the irreducible factor \(G\). Their common zero set
consists of curves and points of bounded degree.

### Counting the coprime roots without good reduction

Work on the intersection of a factor \(G\) and its second
auxiliary polynomial \(H\). Choose an integer linear change
\(x=U+jV+j^2W,\ y=V,\ z=W\) with \(j=O_C(1)\) for which \(x\)
is nonconstant on every geometric curve component. A component
excludes at most two \(j\), as follows by comparing two distinct
points on it. There are only \(O_C(1)\) components by Bezout.
This change has determinant one.

Form the two resultants of the transformed polynomials by
eliminating \(z\) and \(y\), respectively. They are nonzero by
coprimality in the corresponding rational function polynomial
rings. If both polynomials are independent of the eliminated
variable, their plane intersection is finite, so every curve
component has constant \(x\), a case already excluded.
If there are no curves, Bezout already counts the intersection.
The resultant degrees are \(O_C(1)\) and their coefficient heights
are \(\ll n^{O_C(1/k)}\) by their Sylvester determinants. Remove
factors only in \(x\) and take primitive squarefree parts
\(A(x,y)\), \(B(x,z)\). The removed factors account for
\(O_C(1)\) points on the original intersection, because each
fiber of \(x\) on each curve is finite of bounded degree.
Both partial derivatives \(A_y,B_z\) are nonzero generically
along the projection of every curve component. Indeed,
\(A\) is squarefree over \(\mathbb Q(x)\), and the same
holds for \(B\). A component projects dominantly in \(x\)
and therefore cannot have its projection inside the common
zeros of the squarefree equation and its derivative.
Their exceptional points likewise cost \(O_C(1)\).

Regard \(A,B\) as integer polynomials in the original variables.
In the torus put
\[
 \omega=p\,dU/U+q\,dV/V-r\,dW/W,\qquad
 J_0=\det\begin{pmatrix}
 \nabla A\\ \nabla B\\ pVW& qUW& -rUV
 \end{pmatrix}.
\]
On each curve at a point where \(A_yB_z\ne0\), the two equations
\(A=B=0\) have a one-dimensional tangent space and the
logarithmic differential vanishes on that tangent space
precisely when \(J_0=0\). If \(J_0\) vanishes identically
on a geometric component that meets the torus, then
\(\phi=U^pV^qW^{-r}\) is constant there, by characteristic
zero. Its exact value, if attained at an integer solution,
is rational and gives divisor-many coprime roots as above.
On all other components the zeros of \(J_0\) cost \(O_C(1)\)
by Bezout.

At a remaining solution \(u=(U,V,W)\), take any integer
representative of \(K\bmod n\) and put \(f=U^pV^q-KW^r\).
Modulo each \(P^{e_P}\), the Jacobian of \((A,B,f)\) at \(u\)
equals \(U^pV^q/(UVW)\) times \(J_0(u)\), in determinant.
The prefactor is a unit at all primes dividing \(n\).
Because \(J_0\) has bounded degree and small coefficients,
its nonzero integer value satisfies
\[
 \sum_{P\mid n} v_P(J_0(u))\log P
 \le \log|J_0(u)|
 \le (E(C)/k)\log n+O_{C,p,q,r}(1).
\]
Thus the bad reduction of the curve costs only this much
precision, regardless of the height of \(K\). In particular,
the estimate concerns the bounded-degree logarithmic Jacobian,
which contains no coefficient \(K\), rather than the
potentially very large value of the original Jacobian.

At each prime with \(\eta_P=v_P(J_0(u))\le e_P/3\), the
actual Jacobian determinant has valuation \(\eta_P\).
Strong multivariate Hensel lifting produces an exact solution
of \(A=B=f=0\) in \(\mathbb Z_P^3\) congruent to \(u\) modulo
\(P^{e_P-\eta_P}\), hence modulo \(P^{\lceil2e_P/3\rceil}\).
For clarity, multiplying the original residual vector by the
inverse Jacobian gives a correction in \(P^{e_P-\eta_P}\).
The Taylor error has valuation at least \(2(e_P-\eta_P)>e_P\),
and the Jacobian valuation stays \(\eta_P\); iteration
gives the lift. Its Jacobian remains nonzero, so it is
isolated among the algebraic solutions. There are at most
\(O_{C,p,q,r}(1)\) such isolated exact points over
\(\mathbb Q_P\), by Bezout. This bound can depend on the
original exponents. No uniform local contact-order bound
in those exponents is required.

For a given \(u\), record either the isolated lift at each
prime or that the prime was discarded. The discarded prime
powers together have logarithm at most
\(3E(C)\log n/k+O_{C,p,q,r}(1)\). Thus the lifted congruences
have a common modulus with logarithm at least
\((2/3)(1-3E(C)/k)\log n-O_{C,p,q,r}(1)\).
For large \(k\) in terms of \(C\), this exceeds the threshold
for coordinatewise reconstruction of roots at height
\(n^{C/k}\). Any two solutions with the same record coincide.
There are at most \(O_{C,p,q,r}(1)^{\omega(n)}\) records, hence
\(n^\varepsilon\) roots in a shallow class. Including the
constant components, the finite exceptional points and
the shallow classes proves the asserted \(2/k\) exponent.
The finitely bounded ranges of \(n\) at each fixed choice
of parameters enter the implied constant.

To amplify, color primes with three colors, one for each
coordinate. Keep a triple if all the primes of each of
its roots have the matching color. Pairwise coprimality
makes the expected retention \(3^{-\omega(UVW)}\), which
is \(n^{-o(1)}\) uniformly in the height box. Some coloring
keeps this many triples. Coordinatewise products of \(\ell\)
retained triples remain pairwise coprime and solve the
congruence with \(K^\ell\). The number of input tuples
giving any one product triple is at most the product of
three \(\ell\)-fold divisor functions, or \(n^\varepsilon\).
Apply the unamplified bound with height constant \(\ell C\)
and take the \(\ell\)-th root to obtain \(2/(\ell k)\).

### Transfer to boxes with two heavy roots on each letter

For sufficiently small fixed \(\lambda\), consider boxes
in which every letter has two selected coordinates, and
the total mass outside them is \(o(\lambda)\) as
\(\lambda\to0\). If the box approaches the half-budget
bound, their masses are \(L_T/2+o(\lambda)\) by the
previous stability statement. Choose a tall modulus
letter \(M\), and lock whichever of its heavy coordinates
carries larger weighted height. It has weighted height
at least \(1/2-o(1)\) when the complementary weighted
height is \(O(\sqrt\lambda)\). Freeze this root \(z\),
and freeze one heavy root \(x\) of a second letter \(F\).
Freeze the complementary coordinates of the two free
letters. One remaining root of \(F\), and the two
remaining roots of the other letter \(J\), give a
three-variable congruence with modulus \(n=z^d\).
Choose the sign of the monomial so the two roots of \(J\)
are in the numerator. The root of \(F\) is in the
denominator. All three roots are pairwise coprime by
the anatomic factorization, as well as units modulo \(n\).
Their values determine the two free letters, hence also
the modulus letter.

Take the fullness threshold \(k=d\) itself. If the locked
root mass \(\beta_M\gg\lambda\), then all the free roots
have height \(Y^{O(\lambda)}\le z^C=n^{C/d}\) for a fixed
\(C\). This does not require a positive lower bound on the
weighted height \(d\beta_M\). The heavy degree \(d\) tends
to infinity with \(1/\lambda\), so a sufficiently small
\(\lambda\) permits any fixed amplification order \(\ell\).
For that order the modular loss is \(n^{2/(\ell d)}=z^{2/\ell}\),
or \(O(\lambda/\ell)\) in the \(Y\)-exponent. At a fixed
anatomic cutoff the degrees run through a finite set, and
the small finite-box height errors can be absorbed into
\(Y^\varepsilon\). The resulting bound is
\[
 D\le (L_F+L_M)/2+o(\lambda)+O(\lambda/\ell).
\]
For a near-half sequence with all three letters double-heavy,
choose \(J\) of largest mass and choose \(M\) of positive
fractional mass among the other two. The latter is possible
by the pair counts, which give
\(\max_i L_i\le\lambda/2+o(\lambda)\), while
\(L=\lambda+o(\lambda)\). Taking \(\ell\) fixed and large
before letting \(\lambda\to0\) cuts this configuration.

If there are two double-heavy letters and one single-heavy
letter of positive fractional mass, take the single-heavy
letter as \(M\), even if its logarithmic height tends to
zero. The cofactors on \(M\) are left unfixed. If the third
letter instead has mass \(o(\lambda)\), freeze it in full,
lock a heavy root on one double-heavy letter, and let the
two roots on the other double-heavy letter vary. In the
modular lemma set \(W=1\). Both choices give a strict saving
below half for these near-half configurations.
The profiles with fewer than two double-heavy letters use a
second degree cutoff to shrink the residual coefficient height.

### A fixed-proportion saving below the half-budget exponent

There are absolute \(\eta>0,\lambda_0>0\) such that, for
\(0<\lambda<\lambda_0\),
\[
 N_\lambda(X)\ll_{\lambda,\varepsilon}
 X^{(1/2-\eta)\lambda+\varepsilon}.
\]
This proof gives existence of the gap. It does not optimize
or evaluate \(\eta\) or \(\lambda_0\).

Suppose no such gap holds and use near-half boxes along a
sequence \(\lambda\to0\), with \(D/\lambda\to1/2\).
As in the half-budget proof, \(L/\lambda\to1\). Pass to a
subsequence on which the three mass fractions \(L_T/\lambda\)
converge. By pair counting their limits are at most \(1/2\),
so at least two are positive.

Apply the heavy-coordinate stopping argument at both
\(K_0=\lceil\lambda^{-1/2}\rceil\) and
\(K_1=\lceil\lambda^{-1/4}\rceil\). For any positive-mass
letter the same threshold
\(L_T/2-(L/2-D)_+-3\lambda/K_i\) is used at cutoff \(K_i\).
Every root heavy at \(K_0\) is heavy at \(K_1\).
There are at most two heavy roots per positive-mass letter,
each of mass at least \(L_T/2-o(\lambda)\).
The weighted height outside them is \(O(K_i\lambda+\lambda)\),
by the stopping argument for (A). In particular a letter
with just one heavy root at \(K_1\) has complementary height
\(O(\lambda^{3/4})\). If it has a heavy root at \(K_0\), that is
the same root.

The following one-free-root count will apply in both cases
below. Lock a root of degree \(d\ge K_0\) and mass
\(\beta_M\asymp\lambda\), taking \(n=z^d\) as modulus.
The free letter \(J\) has one heavy root of degree \(q\),
and the other free letter is either fixed in full or
has one heavy root whose value is fixed. In the latter
case the unknown coefficient is the ratio of the two
complements, as in (18). When these complements have
weighted height \(O(\lambda^{3/4})\), the coefficient
has height \(Y^{O(\lambda^{3/4})}\le n^{O(\lambda^{1/4})}\).
Assume \(q\gg d\), with a constant along the subsequence.
Choose \(k=\lfloor c d\rfloor\), with \(c>0\) small enough
in terms of the fixed lower mass fraction and the degree
ratio. Then \(k\le d\), \(H(T)\le Y^{2\lambda}\le n^{2/k}\),
and \(q>h=\lfloor k/(256\ell^2)\rfloor\).
For any fixed \(\ell\), sufficiently small \(\lambda\)
now satisfies the scalable lemma, giving cost
\(O(\lambda/\ell)\) in the box exponent.
Here \(n^{1/(\ell k)}=z^{d/(\ell k)}\), so this remains true
if the modulus has vanishing logarithmic height.
The fixed-root ratio argument recovers the triple. If one
complete letter was fixed, the free root and its coefficient
already recover the other free letter and hence the triple.

First suppose all three limiting mass fractions are positive.
By the near-half stability argument, all three letters have
a heavy root at \(K_0\), including the short letter.
If two or more letters are double-heavy at \(K_1\), use the
three-variable transfer. If a letter is single-heavy, it
can be the modulus letter. It still has a \(K_0\) root.
In the other cases at least two letters are single-heavy
at \(K_1\). Use two such letters as \(F,J\), choosing a
tall one as \(J\). Use the remaining letter as \(M\).
Choose the tall root of \(M\) when \(M\) is tall, meaning
one of its \(K_0\) roots carrying a fixed positive weighted
height. In the short letter use any \(K_0\) heavy root.
The free degree \(q\gg1/\lambda\), while \(d=O(1/\lambda)\),
so the one-free-root count applies. The chosen heavy masses
are \(L_i/2+o(\lambda)\): the lower bounds come from stopping
and the upper bounds from releasing one chosen root per
letter in the surface estimate, with the two tall
degrees of order \(1/\lambda\). For a double-heavy letter
the upper bound also follows from the mass budget directly.
Therefore
\[
 D\le (L-L_J)/2+o(\lambda)+O(\lambda/\ell).
\]
The positive limiting fraction for \(J\) and a sufficiently
large fixed \(\ell\) contradict \(D/\lambda\to1/2\).

Now one letter has mass \(o(\lambda)\). Fix it in full at
cost \(Y^{o(\lambda)}\). Each other letter has mass
\(\lambda/2+o(\lambda)\). If one of them has no \(K_0\) heavy
root, its weighted height is \(O(\sqrt\lambda)\), so it is
the short letter. Freeze the complementary coordinates of
one heavy root of the other positive-mass letter, which is
tall. Its whole value changes in steps of size
\(Y^{1-O(\lambda)}\) as the root increases by one, while
the allowed difference from the fixed tall letter is
\(Y^{O(\sqrt\lambda)}\). There is at most one choice of
that root. The resulting cost is at most
\(\lambda/4+o(\lambda)\).

We may assume these two letters both have \(K_0\) roots.
If either is double-heavy at \(K_1\), lock a root on the
other and set \(W=1\) in the three-variable lemma, as
above. Otherwise they are both single-heavy at \(K_1\).
Lock the root with the smaller degree, so \(q\ge d\) and
the one-free-root count gives \(D\le\beta_M+O(\lambda/\ell)
+o(\lambda)\). A second count fixes the negligible letter
and the complements of the two heavy roots, leaving
\(\epsilon_1 A x^d+\epsilon_2 B y^q=N\) with \(N\ne0\).
This affine curve is absolutely irreducible by Eisenstein
in \(\mathbb C[y][x]\) at a simple root of
\(\epsilon_2 B y^q-N\). The coefficient-uniform
Bombieri--Pila bound in a square of side \(Y^{2\lambda}\)
gives exponent \(O(\lambda^{3/2})\) because \(d,q\ge K_0\).
The second count therefore gives
\(D\le L-\beta_M-\beta_J+o(\lambda)\).
Adding the two bounds and using
\(\beta_J\ge\lambda/4-o(\lambda)\) gives
\(2D\le3\lambda/4+O(\lambda/\ell)+o(\lambda)\),
again a contradiction.
[Bombieri--Pila, Theorem 5](https://people.maths.ox.ac.uk/pila/Ovals.pdf)

All these bounds transfer at a fixed anatomic truncation.
After choosing a sufficiently large fixed \(\ell\) along
each hypothetical sequence, choose the decomposition losses
small relative to \(\lambda\) and the desired exponent error.
The degrees are then finite at each fixed budget and
truncation. The modular constants can depend on them, and
every bound is uniform in the actual frozen coefficients and
roots. The pair bounds, stopping inequalities and weighted
heights tolerate the resulting \(o(\lambda)\) losses.
If the stated uniform gap failed, taking successively
smaller target gaps and budgets in this finite-box reduction
would yield precisely the near-half sequence ruled out above.

## Height concentration and prime roots near the half bound

Consider a sequence of boxes as \(\lambda\to0\), with
\(D=\lambda/2+o(\lambda)\) and \(L_T\gg\lambda\) on all three
letters. The two tall letters have height one. Then the short
height \(h_A\) satisfies \(1-h_A=O(\sqrt\lambda)\). If both
tall letters have two heavy coordinates, the conclusion strengthens
to \(1-h_A=O(\lambda)\).

Use the [near-half stability statement](#the-stability-supplied-by-the-half-bound-proof).
On each letter there are one or two heavy coordinates. Each
has mass \(L_T/2+o(\lambda)\). To see this for any choice of
one heavy coordinate per letter, freeze the others and use
Browning–Verzobio with two exponents at least
\(\lceil\lambda^{-1/2}\rceil\). The surviving equation is
primitive and all three root ranges have mass at most \(\lambda\).
Theorem 1.2 gives an exponent \(O(\lambda^{5/4})\) for the free
roots. Thus the sum of the three chosen masses is at most
\(L-D+o(\lambda)=L/2+o(\lambda)\). Each has the opposite
bound \(L_T/2-o(\lambda)\) by stability. If two roots are heavy
on one letter, their complementary mass is \(o(\lambda)\).
The complementary weighted height is \(O(\sqrt\lambda)\)
whether the letter has one or two heavy roots.

### Counting a tall product in a short real interval

Fix one complete tall term \(N\asymp X\) and put \(\tau=1-h_A\).
The other tall value lies within \(O(X^{1-\tau})\) of \(N\).
For a two-heavy tall letter, freeze its complementary coordinates
at cost \(o(\lambda)\). Its remaining roots satisfy
\[
 Q x^p y^q=N(1+O(X^{-\tau}))
\]
with positive \(Q\) fixed. At least one of \(p,q\) is
\(\gg1/\lambda\), since together their terms carry weighted
height \(1-O(\sqrt\lambda)\). Suppose \(\tau/\lambda\to\infty\).
Let \(m\asymp\sqrt{\tau/\lambda}\) be a positive integer chosen
so that \(m^2\lambda\le\tau/16\).

Write \(g=\gcd(p,q)\). If \(p/g,q/g<m\), then
\(1/g=O(m\lambda)=o(\tau)\). The integer
\(x^{p/g}y^{q/g}\) is within \(O(X^{1/g-\tau})\) of
\((N/Q)^{1/g}\). Consequently there is at most one possible
value for this integer, and the divisor bound gives
\(O_\varepsilon(X^\varepsilon)\) pairs with that value.

Otherwise use the \(n=m^2\) monomials \(x^i y^j\) with
\(0\le i,j<m\). They have distinct exponents
\((q/g)i-(p/g)j\) on the exact monomial graph. For a dyadic
scale \(Z\) of \(x\), put \(u=(x/Z)^{g/q}\) and normalize \(y\)
by \(Y_0=(N/(QZ^p))^{1/q}\). On a fixed compact positive
\(u\)-interval, the normalized monomial evaluations differ
from \(u^{(q/g)i-(p/g)j}\) by \(O(X^{-\tau})\).
Each \(x\) gives at most one \(y\), because successive
\(y\) values change the tall product by relative order
at least \(X^{-\lambda}\). Distinct \(u\) values are separated
by \(\gg X^{-\lambda}\), with a constant allowed to depend on
the fixed degrees.

Divide an \(n\)-point evaluation determinant by the Vandermonde
in \(u\) using divided-difference rows. On the exact graph
the quotient is bounded away from zero. Indeed a common
shift of the distinct integer exponents gives a Schur
polynomial times a monomial, and the Schur polynomial is
positive on the compact positive interval. The exact quotient
is also bounded above. In the near-graph determinant the
error in a divided difference of order \(r\) is
\(O(X^{-\tau+r\lambda})\), so both conclusions persist for
large \(X\). Return to the integer evaluations. Restoring the
monomial scales costs at most \(X^{O(\lambda m n)}\) in the
determinant. On a subinterval in \(u\) of length
\(X^{-C\lambda/m}\) with a sufficiently large fixed \(C\),
the Vandermonde more than cancels this cost. Any \(n\) points
in such an interval would give a nonzero integer determinant
of absolute value less than one. Summing the intervals gives
\[
 O_{\varepsilon}(X^{O(\lambda\sqrt{\lambda/\tau})+\varepsilon})
 \tag{19}
\]
pairs, with degree and box parameters fixed before \(X\to\infty\).

For a one-heavy tall letter write the free term as \(Qx^p\).
The low cofactor \(Q\) is now variable and has numeric height
at most \(X^{a_0}\), with \(a_0=C_0\sqrt\lambda\) for fixed
\(C_0>0\). Assume \(\tau/a_0\to\infty\).
Work on \(u=x/Z\) and normalize \(Q\) by \(N/Z^p\).
The exact graph has powers \(i-pj\) on the monomials
\(x^iQ^j\). Use a rectangle of side counts
\[
 m_1\asymp\sqrt{\tau a_0}/\lambda,\qquad
 m_2\asymp\sqrt{\tau/a_0},
 \qquad n=m_1m_2\le\tau/(16\lambda).
\]
Both side counts grow and \(m_1<p\), since \(p\gg1/\lambda\).
The powers \(i-pj\) are therefore distinct. Each \(x\) gives
at most one \(Q\), since consecutive cofactor values change
the product by relative order at least \(X^{-a_0}\).
The same divided-difference bounds apply, with monomial scales
costing at most \(X^{O(n\sqrt{\tau a_0})}\). Intervals
of length \(X^{-C\lambda\sqrt{a_0/\tau}}\) now suffice.
This gives \(X^{O(\lambda\sqrt{a_0/\tau})+\varepsilon}\)
pairs for the free tall term.

If \(\tau/\sqrt\lambda\to\infty\), either tall term can thus
be counted at exponent \(o(\lambda)\) after fixing the other
term and paying any two-heavy complementary mass. The same
holds when \(\tau/\lambda\to\infty\) and both tall terms are
two-heavy. In either case
\[
 D\le\min(L_B,L_C)+o(\lambda)
 \le (L-L_A)/2+o(\lambda),
\]
which contradicts \(L_A\gg\lambda\). Applying this to a
subsequence where the respective height-gap ratio is unbounded
proves the stated concentration.

### Two further restrictions on the heavy coordinates

Any two heavy degrees on one letter must be coprime in a
near-half sequence. If \(d,e\) on one letter have
\(g=\gcd(d,e)\ge2\), combine their factors into the \(g\)-th
power of \(x^{d/g}y^{e/g}\). Free this combined root along
with one heavy root on each other letter. The latter two
degrees are at least \(\lceil\lambda^{-1/2}\rceil\), so Theorem
1.2 again costs \(O(\lambda^{5/4})\) in free roots. The
combined root has mass \(L_T+o(\lambda)\) or less and
determines its factors up to divisor-many choices.
Counting the complementary coordinates costs at most
\((L-L_T)/2+o(\lambda)\) in all. This contradicts the assumed
near-half exponent.

Also, up to subfamilies of exponent strictly below the
putative \(D\), each heavy root is a prime times an integer
of mass \(o(\lambda)\). To see this, take a slack
\(\xi=\xi(\lambda)\to0\) large enough that \(\xi\lambda\)
dominates the deviations of \(D\) from \(L/2\) and of every
heavy mass \(\beta_T\) from \(L_T/2\). A divisor of mass in
\([2\lambda/d,\beta_T-\xi\lambda]\) in degree \(d\) supplies
the lock (A), with exponent at most
\((L-L_T)/2+\beta_T-\xi\lambda<D\).
On the remaining subfamily the root must have a prime factor
of mass above \(\beta_T-\xi\lambda\). Otherwise each prime
has mass less than \(2\lambda/d\). Multiply the prime factors,
with multiplicity, until their product first has mass at least
\(2\lambda/d\). Its mass is less than \(4\lambda/d\), inside
the forbidden interval, since \(d\ge\lceil\lambda^{-1/2}\rceil\)
and \(\beta_T\gg\lambda\).
The product of all remaining primes has mass at most
\(\xi\lambda\) up to vanishing box errors. Its count can be
paid at exponent \(o(\lambda)\) when fixing that cofactor.

[Browning–Verzobio, Theorem 1.2](https://arxiv.org/html/2608.24512v1#S1)

## A partial-power modulus for clustered high degrees

Consider the small-budget anatomic boxes above. Denote the
coordinate masses on each letter \(T\) by \(\alpha_{T,j}\)
and suppose their total is at most \(\lambda\). Choose a set
\(I\) of coordinates of one letter \(M\), and put
\[
 b=\sum_{j\in I}\alpha_{M,j},\qquad
 w=\sum_{j\in I}j\alpha_{M,j}.
\]
Assume \(w\ge w_0>0\) and \(j\ge c_0/\lambda\) for every
\(j\in I\). On the two free letters \(F,J\), choose a common
integer \(q\ge c_0/\lambda\) and bounded nonnegative integers
\(m_{T,j}\). Suppose
\[
 E=\sum_{T\in\{F,J\}}\sum_j
       |j-qm_{T,j}|\alpha_{T,j}=o(1)
       \quad\text{as }\lambda\to0.
\]
The box exponent satisfies
\[
 D\le b+O\bigl(\lambda E^{1/3}+\lambda^{3/2}\bigr).
 \tag{20}
\]
The constant can depend on \(w_0,c_0\) and the multiplier bound.
As usual the point count has an additional \(X^\varepsilon\)
factor, and the finite truncation losses can be chosen in
terms of the fixed parameters.

Freeze the coordinates of \(M\) indexed by \(I\), at cost
\(X^{b+\varepsilon}\). Their product with the full multiplicities
is a modulus \(n=\prod_{j\in I}x_{M,j}^{j}\asymp X^w\).
For the free letters write their complete values as
\[
 N_F=T_F^q V_F,\qquad N_J=T_J^q V_J,\qquad
 T_H=\prod_j x_{H,j}^{m_{H,j}}\quad(H=F,J).
\]
Here \(V_H\) may be rational, including when a rounded
multiplicity exceeds the true one. Small exceptional
coefficients in the finite anatomic decomposition go
into \(V_H\). Put \(T=T_F/T_J\) and \(V=V_F/V_J\).
Primitivity makes \(T,V\) units modulo \(n\), and the
signed equation implies \(VT^q\equiv K\bmod n\) for the
fixed sign \(K\in\{1,-1\}\).
Furthermore \(H(T)\ll X^{C_0\lambda}\), with \(C_0\) a bound
on the \(m_{T,j}\), and \(H(V)\ll X^E\), up to the
arbitrarily small anatomic loss. Each pair \((T,V)\)
determines the positive rational ratio \(VT^q=N_F/N_J\).
Its reduced numerator and denominator recover the two
coprime complete terms, and their signed sum determines
the remaining term.

Choose an integer \(k\asymp1/\lambda\) at most \(q\) and
each modulus exponent, with \(k\le w/((C_0+1)\lambda)\).
Then \(n\) is \(k\)-full and \(H(T)\le n^{2/k}\) for
sufficiently large \(X\). Choose
\[
 \ell=\left\lfloor c_1
           \min(E^{-1/3},\lambda^{-1/2})\right\rfloor,
\]
omitting the \(E^{-1/3}\) constraint when \(E=0\).
For sufficiently small fixed \(c_1>0\), we have
\(k\ge1024\ell^2\) and \(H(V)\le n^{1/(16384\ell^3)}\).
The strict margins absorb the exceptional coefficients
and dyadic constants. We also have
\(q>\lfloor k/(256\ell^2)\rfloor\).
The scalable congruence bound with an independent
power \(q\) now counts \(O(X^{w/(\ell k)+\varepsilon})\)
pairs for each modulus. Since \(w\le1\) in limiting boxes
and \(1/\ell=O(E^{1/3}+\sqrt\lambda)\), this proves (20).

For instance, consider a near-half configuration with two
heavy coordinates on each letter, of masses \(L_T/2+o(\lambda)\),
and each complete logarithmic height bounded below. Suppose
the four heavy degrees on the two free letters differ by
\(o(1/\lambda)\) from a common \(q\asymp1/\lambda\).
Use multiplier one on these four coordinates and zero on
their complements. The complementary weighted heights are
\(O(\sqrt\lambda)\), so \(E=o(1)\). Choose the larger-height
heavy coordinate of \(M\). Its height is bounded below
and its mass is \(L_M/2+o(\lambda)\); its degree therefore
exceeds \(c_0/\lambda\) for some \(c_0>0\).
The count exponent is at most \(L_M/2+o(\lambda)\).
If every heavy degree belongs to the same cluster, take
\(M\) of minimal mass to obtain \(D\le\lambda/6+o(\lambda)\).
For heavy degree deviations \(O(\lambda^{-1/2})\), (20)
gives error \(O(\lambda^{7/6})\) in terms of the actual
selected coordinate mass.

The same argument also covers bounded multiples of a common
base. For example, with two heavy coordinates on each letter,
the two degrees on one free letter may be near \(q,2q\) and
on the other near \(2q,3q\), where all the differences are
\(o(1/\lambda)\). Their product bases use multiplicities
\(1,2\) and \(2,3\), respectively. This gives the same bound
on fixing one high coordinate of the modulus letter. The
exponent on the chosen modulus coordinate can be unrelated
to the four free exponents.

## An affine surface after fixing a complete letter

Let \(p,q\ge1\), \(r\ge5\), and let \(A_1,A_2,A_3\) be
nonzero integers of absolute value at most \(X\). For
\(2\le H\le X\), the number of positive integer solutions
up to \(H\) in each coordinate of
\[
 f(x,y,z)=A_1 x^p y^q+A_2 z^r+A_3=0
\]
is
\[
 O_{p,q,r,\varepsilon}\bigl(X^\varepsilon H^{3/\sqrt r}\bigr).
 \tag{21}
\]

The polynomial is absolutely irreducible. Over
\(\overline{\mathbb Q}(y)[z]\), its leading coefficient
as a polynomial in \(x\) is a unit. Its constant term
has a simple zero at any root of \(A_2 z^r+A_3\).
Eisenstein and Gauss's lemma apply. The polynomial is
also primitive as a polynomial in \(x,z\) over
\(\overline{\mathbb Q}[y]\).

Apply Browning–Verzobio's affine-surface covering theorem
with modulus one and side polynomial \(g=x_2+1\), which imposes no congruence condition.
For equal sides \(H\), choose lexicographic order with the
\(z\) exponent first. The leading monomial is \(z^r\),
so their parameter is \(H^{1/\sqrt r}\).
For any fixed \(\epsilon>0\) the theorem covers the
integer points by \(O_{p,q,r,\epsilon}(H^{1/\sqrt r+\epsilon})\)
auxiliary hypersurfaces of bounded degree \(O_{p,q,r,\epsilon}(1)\),
each coprime to \(f\), and
\(O_{p,q,r,\epsilon}(H^{2/\sqrt r+2\epsilon})\) exceptional
points. Bézout gives \(O_{p,q,r,\epsilon}(1)\) curve
components per auxiliary hypersurface, with the same
upper bound on their degrees.
[Browning–Verzobio, Theorem 1.2 and (6)](https://arxiv.org/html/2408.11453v2#S1)

On an absolutely irreducible curve component of projective
degree \(d\), work on its smooth projective normalization.
A coordinate vanishing identically gives no positive points.
If \(z\) is constant, the equation fixes the nonzero value of
\(x^p y^q\) for every positive point on that component.
The value, when integral, is at most \(2XH^r\), so the
divisor bound gives \(O_{r,\varepsilon}(X^\varepsilon)\)
positive points.

Otherwise \(F=(A_2/A_3)z^r\) is nonconstant. Its degree as
a map to the projective line is at least \(r\). Its zeros
and poles and the zeros of \(1+F=-(A_1/A_3)x^p y^q\) lie
among the zeros and poles of \(x,y,z\), a set \(S\) of
size at most \(6d\). Each coordinate is a ratio of two
linear forms on the projective curve, so its zeros and
poles have total degree at most \(d\) each. The genus
\(g\) of the normalization is at most \((d-1)(d-2)/2\)
by a birational plane projection. Apply the three-term
function-field \(abc\) inequality
\(\deg F\le 2g-2+|S|\). Equivalently, apply Riemann–Hurwitz
to \(F\), counting the ramification over \(0,-1,\infty\)
alone. This gives
\[
 r\le d^2+3d\le4d^2.
\]
Since \(r\ge5\), the degree satisfies both \(d>1\) and
\(1/d\le2/\sqrt r\).
The coefficient-uniform bound for integer points on
affine curves gives \(O(d^2H^{1/d}(\log H)^\kappa)\)
points on each such component defined over \(\mathbb Q\).
Components exchanged by Galois conjugation instead
contribute only \(O(d^2)\) rational intersection points
per component. Absorb the logarithmic factors and choose
the intermediate \(\epsilon\) small in terms of
\(\varepsilon\). This proves (21).
[Binyamini–Cluckers–Kato, Corollary 8](https://arxiv.org/html/2401.03982v2#S3)

In an anatomic box, fix one complete letter \(M\)
with \(X^{L_M+\varepsilon}\) choices. On a second
letter \(F\), leave two coordinates free. On the remaining
letter \(J\), leave one free coordinate of degree
\(r\gg1/\lambda\). Freeze the complementary coordinates,
of total mass \(f_F+f_J\), including any exceptional
coefficients at arbitrarily small additional cost.
The remaining equation has the shape of (21).
The three root heights are at most \(X^\lambda\)
in a limiting box and the coefficient heights are
at most \(X\), apart from bounded dyadic factors
which can be absorbed in the ambient scale. Thus
\[
 D\le L_M+f_F+f_J+O(\lambda^{3/2}).
 \tag{22}
\]
The count uses only a high degree on the singly released
letter. The original ordering of the three terms affects
only the signs of the coefficients.

In particular, consider a near-half box with one letter
of mass \(o(\lambda)\) and two letters of logarithmic
heights bounded below. If both remaining letters have
two heavy coordinates, fix the negligible-mass letter
and leave both heavy coordinates of one letter free.
On the second letter leave a heavy coordinate of
macroscopic weighted height free and fix the other.
The other root costs \(L_J/2+o(\lambda)\), while the
complementary masses outside those heavy coordinates
are \(o(\lambda)\). The pair bounds force each of the
two macroscopic masses to be \(\lambda/2+o(\lambda)\),
so (22) gives \(D\le\lambda/4+o(\lambda)\).
If one remaining letter is single-heavy and the other
double-heavy, leave the three heavy coordinates free.
The single-heavy coordinate has degree \(\gg1/\lambda\),
and its complementary mass is at most
\(L_J/2+o(\lambda)=\lambda/4+o(\lambda)\).
This gives the same exponent by (22).

### Radical shares in near-half boxes

Every near-half sequence has \(L_T\gg\lambda\) on all three
letters. Suppose instead, on a subsequence, that one letter
has mass \(o(\lambda)\). The pair bounds force the other two
masses to be \(\lambda/2+o(\lambda)\). If the negligible-mass
letter is tall and the short height gap is larger than
\(\sqrt\lambda\) by an unbounded factor, fix that tall term.
The real determinant estimate (19), or its one-heavy version,
counts the other tall value at exponent \(o(\lambda)\).
The negligible-mass term and this tall value determine the
short term. Such a subsequence cannot attain the half
exponent. Thus both macroscopic-mass terms have height
\(1-O(\sqrt\lambda)\), and there is at least one heavy root
on each by the stability argument.

The affine-surface count treats the cases with two heavy
roots on at least one of these terms. If both terms have
only one heavy root, denote their masses by \(\beta_B,\beta_C\)
and their degrees by \(k_B,k_C\). Fixing the negligible-mass
term and the complementary coordinates gives a binary
equation \(A x^{k_B}+B y^{k_C}=C\) with \(ABC\ne0\).
Eisenstein at a simple root of \(C-B y^{k_C}\) proves
absolute irreducibility over \(\overline{\mathbb Q}\).
Both degrees are \(\gg1/\lambda\), so the coefficient-uniform
Bombieri–Pila bound counts the free roots with exponent
\(O(\lambda^2)\). It follows that
\(\beta_B+\beta_C\le L-D+o(\lambda)=\lambda/2+o(\lambda)\).
The lower stability bounds on the two masses then give
\(\beta_B=\beta_C=\lambda/4+o(\lambda)\).

Apply the divisor-lock argument for prime factors to the
root with smaller mass, say the root on \(C\). Its large
prime factor \(P\) has mass \(\rho=\beta_C+o(\lambda)\)
after discarding strata of exponent below \(D\).
Freeze \(P\) and the negligible-mass term. The remaining
free value \(Qz^{k_B}\) satisfies an inverse residue
fibre with unit right side modulo \(P^{k_C}\).
The complementary numeric height \(a_0\) is \(O(\sqrt\lambda)\),
and \(k_C/k_B\ge1-O(\sqrt\lambda)\) because the two
heavy roots each carry height \(1-O(\sqrt\lambda)\) and
\(\beta_C\le\beta_B\).
The [translated inverse determinant](#a-37-small-budget-bound-with-a-short-letter)
uses blocks of root exponent
\(h=\rho-K\sqrt\lambda/k_B\) with a sufficiently large
fixed \(K\). It applies since \(2a_0+k_Bh<k_C\rho\) and
\(h<\rho\). The block count contributes
\(\beta_B-\rho+O(\lambda^{3/2})\), and the curve count
per block has exponent \(O(\lambda^{3/2})\).
Restoring the frozen prime costs \(\rho\), giving
\(D\le\lambda/4+o(\lambda)\) in this last case as well.
Thus a negligible mass on any letter contradicts a
near-half exponent.

For three comparable masses use the height concentration
above, so all heights are \(1-O(\sqrt\lambda)\).
If all three letters have two heavy roots, order the masses
\(L_1\le L_2\le L_3\). In the affine estimate (22), fix
the first complete letter, retain one heavy coordinate
of macroscopic weighted height on the second, and retain
both heavy roots on the third. This costs
\[
 D\le L_1+L_2/2+o(\lambda)
   =L/2-(L_3-L_1)/2+o(\lambda).
\]
Thus surviving all-double profiles have
\(L_T=\lambda/3+o(\lambda)\) on every letter.
If exactly one letter has one heavy root, call its mass
\(L_S\) and write the double-heavy masses as \(L_U,L_V\).
Fixing the cheaper complete double-heavy letter gives
\(D\le\min(L_U,L_V)+L_S/2+o(\lambda)\) by (22);
hence \(L_U=L_V+o(\lambda)\). Fixing the complete
single-heavy letter instead gives
\(D\le L_S+\min(L_U,L_V)/2+o(\lambda)\), so
\(L_S\ge L_U-o(\lambda)\).

More generally two letters cannot each have one coordinate
carrying all but \(o(1)\) of its weighted height, even
if a second heavy coordinate on one of those letters has
nondominant height. To see this, label these two letters
\(F,J\) and the third \(M\). Fix the dominant root on \(F\)
and a heavy root carrying height \(\gg1\) on \(M\), at
total cost \((L_F+L_M)/2+o(\lambda)\). The
[primitive ratio congruence](#recovering-a-primitive-ratio-from-a-dominant-divisor)
has an independent high exponent \(q\gg1/\lambda\) on the
free root of \(J\), and residual height \(X^{o(1)}\) because
the entire cofactor on each of \(F,J\) has height \(o(1)\).
Choose \(k\asymp1/\lambda\) small enough to meet the
root height test, and then take the amplification order
\(\ell\to\infty\) slowly enough to meet the residual
height test. The congruence count costs only \(o(\lambda)\).
The count is bounded by \((L-L_J)/2+o(\lambda)\), which
contradicts \(L_J\gg\lambda\).

## Releasing four roots on three letters

Let \(p,q,r,s\ge1\) be fixed and put \(R=\min(r,s)\).
For nonzero integers \(A,B,C\), count positive integer
solutions with \(x,y,z,w\le H\) and \(\gcd(z,w)=1\) of
\[
 A x^p y^q+B z^r+C w^s=0.
 \tag{23}
\]
Their number is
\[
 O_{p,q,r,s,\varepsilon}(H^{c R^{-1/6}+\varepsilon})
 \tag{24}
\]
for an absolute constant \(c\), uniformly in the coefficients.

### Low-degree subvarieties

Let \(Z\) be a geometrically integral curve or surface
of projective degree \(d\) in the projective closure of
(23). Keep only subvarieties meeting the affine open
set where all four coordinates are nonzero. Suppose
\(d<(R/12)^{1/3}\). Either \(z,w\) are both constant on
\(Z\), or there is a multiplicative identity
\[
 z^u w^v=\gamma,\qquad uv\ne0,\qquad |u|,|v|\le d,
 \tag{25}
\]
on \(Z\).

To prove this, take divisors of the rational functions
\(z,w\) on the projective normalization of \(Z\) over
\(\overline{\mathbb Q}\). For a surface we use a general
hyperplane section of \(Z\), integral of degree \(d\) by
Bertini. We may require that the section meets any
finitely prescribed divisors of the normalization
transversely at points where the normalization and
the prescribed divisors are smooth and where no other
coordinate zero or pole divisor passes. Indeed the
normalization is finite and regular at all generic
points of divisors. The map of a selected divisor to
its image curve is generically separable, and a general
hyperplane meets that image transversely outside the
finitely many exceptional points. The pulled-back
hyperplane then has a nonzero derivative along the
selected divisor. Its curve passes transversely to
the divisor at that point and has the same generic
function field as the original hyperplane section.
Consequently the orders of \(z,w\) at any such divisor
are realized on the normalized section. In particular
they have absolute value at most \(d\), because each
coordinate is the ratio of two linear forms.

On the resulting smooth projective curve, consider
\(F=(B/C)z^r/w^s\) whenever it is nonconstant. The fibers
\(F=0,-1,\infty\) lie among the zeros and poles of
\(x,y,z,w\), using \(1+F=-(A/C)x^p y^q/w^s\).
There are at most \(8d\) such places. Riemann–Hurwitz,
counting ramification over these three fibers, bounds
the degree of \(F\) by \(2g-2+8d\), with \(g\) the genus.
A birational plane projection gives
\(g\le(d-1)(d-2)/2\). At every place we obtain
\[
 |r\,\operatorname{ord}(z)-s\,\operatorname{ord}(w)|
 \le \deg F \le d^2+5d \le 6d^2.
 \tag{26}
\]

If the two divisors on \(Z\) are linearly independent,
two rows \((a_i,b_i)\) of their orders have a nonzero
integer determinant. On a curve these are already
places, and on a surface choose the section to retain
both rows. Thus \(F\) on that curve is nonconstant.
Multiply \(r a_i-s b_i\) by the opposite \(b\)-entry
and subtract. The determinant has absolute value
at least one, so \(r\le12d^3\). The same argument
with the \(a\)-entries gives \(s\le12d^3\), contrary
to the degree hypothesis. If exactly one of \(z,w\)
is constant, retain a nonzero order of the other.
Then (25) gives \(R\le6d^2\), also a contradiction.
When both coordinate divisors are nonzero and
proportional, one nonzero row supplies an integer
relation with both coefficients nonzero and bounded
by \(d\). The corresponding rational function has
zero divisor and is constant on the normal projective
variety. This proves (24).

The number of relevant points on such a subvariety
is \(O_{p,q,r,s,d,\varepsilon}(H^\varepsilon)\).
When (24) holds and both exponents have the same sign,
invert if necessary and apply the divisor bound to
the fixed integer \(z^{|u|}w^{|v|}\le H^{2d}\).
If the signs differ, reading the reduced numerator
and denominator determines \(z,w\) by coprimality.
If both coordinates are constant, they already have
one possible value each. For fixed \(z,w\), (23) fixes
\(x^p y^q\), whose value at an admissible point is a
positive integer at most \(H^{p+q}\). Apply the
divisor bound again.

### Covering the hypersurface

The affine polynomial in (23) is absolutely irreducible.
Over \(\overline{\mathbb Q}(y,w)[z]\), regard it as a
polynomial in \(x\). The leading coefficient is a
unit and the constant term \(Bz^r+Cw^s\) is squarefree,
so Eisenstein applies at any of its irreducible factors.
The polynomial is primitive over
\(\overline{\mathbb Q}[y,w,z]\), since its leading and
constant coefficients are coprime there. Hence its
projective closure is geometrically integral, of degree
\(\Delta=\max(p+q,r,s)\).

Use Theorem 0.10 in Salberger's *Counting rational
points on projective varieties*. For a fixed
\(\eta>0\), the relevant height-\(H\) points are covered
by \(O_{\Delta,\eta}(H^{4/\Delta^{1/3}+\eta})\)
proper auxiliary intersections, whose surface
components have bounded degrees depending only
on \(\Delta,\eta\). A geometric component that is
not defined over \(\mathbb Q\) contains rational
points only in its intersection with a distinct
Galois conjugate. Such intersections supply a bounded
number of bounded-degree curves or isolated points.
The same observation disposes of nonrational curve
components by finite intersections.
[Salberger, Theorems 0.9(b) and 0.10](https://londmathsoc.onlinelibrary.wiley.com/doi/full/10.1112/plms.12508)

Write \(d_0=(R/12)^{1/3}\), and assume \(d_0>1\).
All surface or curve components below degree \(d_0\)
have the divisor estimate above. On a geometrically
integral surface defined over \(\mathbb Q\) of degree
\(d\ge d_0\), Theorem 0.9(b) gives
\(O_{d,\eta}(H^{3/(2\sqrt d)+\eta})\) bounded-degree
curves and \(O_{d,\eta}(H^{3/\sqrt d+\eta})\) exceptional
points. Every curve of degree \(e\ge d_0\) has at most
\(O_e(H^{1/e}(\log H)^\kappa)\) integral points in
the affine chart, by the coefficient-uniform curve
bound. Consequently the total exponent is at most
\[
 \frac4{\Delta^{1/3}}+
 \max\left(\frac3{\sqrt{d_0}},
           \frac{3}{2\sqrt{d_0}}+\frac1{d_0}\right)
 +\varepsilon.
\]
This is \(O(R^{-1/6})+\varepsilon\).
For bounded \(R\), a larger absolute constant in
(24) makes the trivial estimate sufficient.
[Binyamini–Cluckers–Kato, Corollary 8](https://arxiv.org/html/2401.03982v2#S3)

### Application to anatomic boxes

Leave two chosen roots free on one letter \(F\), and
one chosen root free on each of the other two letters.
Assume the degree of each singly released root is
\(\gg1/\lambda\). Freeze all complementary coordinates
and small exceptional coefficients. The remaining
equation is (23), with \(z,w\) coprime since they
come from different primitive letters. The root
heights are at most \(X^\lambda\) in the limiting
box, so (24) gives
\[
 D\le L-\beta_{F,1}-\beta_{F,2}-\beta_J-\beta_M
       +O(\lambda^{7/6}).
 \tag{27}
\]
Here \(L\) is total box mass, \(\beta_{F,1},\beta_{F,2}\)
are the two released masses and \(\beta_J,\beta_M\)
are the singly released masses. The implied constants
in the exponent are uniform when the lower bounds
on the two degrees are uniform multiples of \(1/\lambda\).

## A uniform small-budget exponent below one half

There is an absolute \(C\) such that for sufficiently small
fixed \(\lambda>0\) the exceptional-set count satisfies
\[
 N_\lambda(X)\ll_{\lambda,\varepsilon}
 X^{(18/37)\lambda+C\lambda^{7/6}+\varepsilon}. \tag{28}
\]
The four-root bound makes it possible to release both heavy
coordinates of one letter while keeping a power free on
each of the other two.

Put \(v=18/37\), \(t=5/37\), and suppose a box has exponent
\(D\ge v\lambda\). Write \(K_0=\lceil\lambda^{-1/2}\rceil\).
On a letter of mass \(L_T\), call a coordinate heavy if
its degree is at least \(K_0\) and its mass is greater than
\[
 \theta_T=D-(L-L_T)/2-3\lambda/K_0.
\]
At the masses where we apply this definition \(\theta_T>0\).
The stopping argument for (A) bounds the complementary
weighted height by \(O(\sqrt\lambda)\). Indeed the high
degrees of mass at most \(\theta_T\) have total weight
less than \(2\lambda\), and the low degrees have weight
at most \(K_0 L_T\). Also
\(\theta_T\ge L_T/2-\lambda/74-O(\lambda^{3/2})\).
For \(L_T\ge t\lambda\) there are at most two heavy
coordinates, since \(3\theta_T>L_T\) for sufficiently
small \(\lambda\). If there are two, their complementary
mass is at most \(L-2D+O(\lambda^{3/2})\le\lambda/37+
O(\lambda^{3/2})\). A letter of height \(1-O(\lambda^{1/4})\)
in this mass range has a heavy coordinate carrying
height bounded below, and its degree is \(\gg1/\lambda\).

Suppose first that every \(L_T\ge t\lambda\).
If the short height gap \(\tau=1-h_A\) is at least
\(\lambda^{1/4}\), fix one complete tall term and use
the real determinant estimate (19) or its one-heavy
version on the other. The cost is at most
\[
 \min(L_B,L_C)+(L-2D)+o(\lambda)
 \le \left(\frac{1-t}{2}+1-2v\right)\lambda+o(\lambda)
 =\frac{17}{37}\lambda+o(\lambda).
\]
In this application the two-heavy complementary
mass costs at most \(L-2D+O(\lambda^{3/2})\)
and the one-heavy cofactor is left free. Thus we
may take all three heights to be \(1-O(\lambda^{1/4})\).
If each letter has just one heavy coordinate, the
[dominant-divisor count](#recovering-a-primitive-ratio-from-a-dominant-divisor)
gives \(D\le(2/5)\lambda+O(\lambda^{7/6})\).

Otherwise take a letter \(F\) with two heavy coordinates,
and leave both free along with one heavy coordinate of
macroscopic weighted height on each of the other two.
The sum of the four free masses is at least
\[
 2\theta_F+\sum_{T\ne F}\theta_T
 =4D-\frac{3L-L_F}{2}-O(\lambda^{3/2}).
\]
Their root box has side lengths at most \(X^\lambda\).
The two singly released powers each have degree
\(\gg1/\lambda\), so the four-root count (24) costs
\(O(\lambda^{7/6})\) in count exponent. The roots of the
two singly released letters are coprime by primitivity.
Adding the cost of freezing the complement and moving
the \(D\)-terms to one side gives
\[
 D\le L/2-L_F/10+O(\lambda^{7/6})
 \le (18/37)\lambda+O(\lambda^{7/6}).
\]

Now suppose one letter has mass \(s<t\lambda\).
The other two masses belong to
\([D-s,L-D]\subset[13\lambda/37,19\lambda/37]\) by
the pair bounds. Thus the other two letters have
at most two heavy coordinates. Fix the small-mass
letter. If it is tall and the short height gap is
at least \(\lambda^{1/4}\), the real determinant
estimate on the other tall letter bounds \(D\) by
\(s+L-2D+o(\lambda)\le(6/37)\lambda+o(\lambda)\).
In the surviving boxes both macroscopic-mass terms
have heights \(1-O(\lambda^{1/4})\). The small-mass
letter is already fixed at cost \(s\).

If either other letter has two heavy roots, leave
them free and keep one heavy coordinate of
macroscopic weighted height free on the last letter.
The affine-surface count (22) applies. The first
complementary mass is at most \(\lambda/37+
O(\lambda^{3/2})\), and the second is at most
\(19\lambda/74+\lambda/74+O(\lambda^{3/2})\).
Adding \(s\) gives \(D\le(16/37)\lambda+o(\lambda)\).

The two remaining letters may therefore be taken to
have one heavy coordinate each, of masses
\(\beta_B,\beta_C\) and degrees \(k_B,k_C\).
Fixing their low cofactors and the small-mass term
gives a binary equation \(Ax^{k_B}+By^{k_C}=C\),
with \(ABC\ne0\). It is absolutely irreducible by
Eisenstein at a simple root of \(C-By^{k_C}\).
Its degree is \(\gg1/\lambda\), so Bombieri–Pila
counts the two free roots at exponent \(O(\lambda^2)\).
Restoring the complementary coordinates yields
\[
 \beta_B+\beta_C\le L-D+O(\lambda^2),\qquad
 \beta_T\ge (3D-L-s)/2-O(\lambda^{3/2}).
\]
The second inequality follows from \(\beta_T>\theta_T\)
and \(L_T\ge D-s\).

Choose \(C\) to have the smaller root mass.
The divisor lock (A) saves strictly below \(D\)
if this root has a divisor of mass
\([2\lambda/k_C,\theta_C]\). In the remaining
subfamily it has a prime factor \(P\) of mass
\(\rho>\theta_C\gg\lambda\). Indeed, otherwise
all prime factors have mass below \(2\lambda/k_C\),
and a product first reaching \(2\lambda/k_C\)
gives the excluded divisor. Fix \(P\) and the
complete small-mass term. The free tall product
\(Qz^{k_B}\) lies in a unit residue class modulo
\(P^{k_C}\). Here \(k_C/k_B\ge1-O(\lambda^{1/4})\)
by the two weighted heights and
\(\beta_C\le\beta_B\). The translated inverse
determinant uses blocks of length \(X^h\) with
\(h=\rho-C_0\lambda^{1/4}/k_B\). For fixed
\(C_0\) sufficiently large it satisfies
\(2a_0+k_Bh<k_C\rho\) and \(0<h<\rho\), where
\(a_0=O(\sqrt\lambda)\) bounds the cofactor height.
The cost of the root blocks is
\(\beta_B-\rho+O(\lambda^{5/4})\)
and that within each block is \(O(\lambda^{3/2})\).
Restoring the frozen prime therefore gives
\[
 D\le s+\beta_B+O(\lambda^{5/4}).
\]
Combine this with the bounds on \(\beta_B+\beta_C\)
and \(\beta_C\) to get
\(7D\le3(L+s)+O(\lambda^{5/4})\).
As \(s<t\lambda\), this is at most
\((126/37)\lambda+O(\lambda^{5/4})\).
The two mass ranges thus give the same leading
coefficient \(18/37\).

All bounds are in limiting box exponents. At fixed
\(\lambda\), take the anatomic truncation loss
arbitrarily small before summing dyadic boxes.
The degree-dependent constants in the determinant
bounds are finite at each such truncation, and are
absorbed in the implied constant and \(X^\varepsilon\).
As with the half-budget bound, the continued
intersection of the inequalities holds by taking
limits of the normalized coordinate scales and
subtracting the vanishing decomposition loss.
This proves (28).

## Releasing a bounded number of high-degree coordinates

Fix a positive integer \(m\). Partition at most \(m\)
variables into three disjoint sets \(I_1,I_2,I_3\),
which may include empty sets. Let \(a_1,a_2,a_3\)
be nonzero integer coefficients. Consider
\[
 \sum_{T=1}^{3}a_T\prod_{j\in I_T}x_j^{e_j}=0,
 \qquad e_j\ge E.
 \tag{29}
\]
For fixed exponents \(e_j\), uniformly in the
coefficients, the number of positive pairwise
coprime integer tuples with \(x_j\le H\) is
\[
 O_{m,(e_j),\varepsilon}
   \bigl(H^{C_m E^{-1/(m^2-1)}+\varepsilon}\bigr)
 \tag{30}
\]
when \(m\ge2\). The constant \(C_m\) depends only on
\(m\). With at most one variable the statement holds
with \(H^\varepsilon\).

If at most one of the three monomials involves
variables, the equation fixes its value, and the
divisor bound proves the claim (or there are no
positive points). If at least two monomials involve
variables, the affine polynomial is absolutely
irreducible. Regard it as a polynomial in one
variable from the first nonempty group. Over the
rational function field in all other variables
except one from the second nonempty group, its
constant coefficient is a binomial with a nonzero
constant term, hence squarefree in that remaining
variable. Eisenstein applies at any irreducible
factor of the binomial, and Gauss's lemma uses
disjointness of the sets of variables in the
monomials. The projective closure is then
geometrically integral.

Let \(Z\) be any geometrically integral positive-
dimensional subvariety of that closure, meeting the
open set where every affine coordinate is nonzero.
Let \(d\) be its degree. If every coordinate divisor
on the projective normalization is linearly
independent over \(\mathbb Q\) from the others,
select an invertible valuation matrix \(M\), with
one column per coordinate and one row per selected
prime divisor. Restrict to a smooth curve
normalizing a general linear section of \(Z\).
All selected valuation rows can be realized on
the section. Indeed the finite normalization maps
each selected divisor to a subvariety of the same
dimension. At a point where the divisor, the
normalization and the image are smooth, its
differential along the divisor is injective on a
dense open subset by generic separability in
characteristic zero. A general complementary
linear section meets the image transversely
in that open subset, avoiding the images of its intersections with the
other coordinate zero or pole divisors. Its pullback meets the selected divisor
transversely by the differential condition.
Bertini gives an integral image curve of degree
\(d\), and the pullback is birational to it because
the normalization is birational and the generic
section avoids its exceptional locus. Thus the
restriction retains all selected orders. Each
entry of \(M\) has absolute value at most \(d\),
as each coordinate is a ratio of linear forms on
a curve of degree \(d\).

Choose a ratio \(F\) between two of the three
signed monomials in (29), with at least one of them
involving variables. Then \(F\) is nonconstant on
this curve since the valuation matrix is invertible.
Together \(F\) and \(1+F\) have zeros and poles only
at coordinate zeros and poles, of which there are
at most \(2md\). Riemann–Hurwitz for \(F\) over
\(0,-1,\infty\) gives \(\deg F\le 2g-2+2md\), with
\(g\le(d-1)(d-2)/2\). If \(\mathbf{a}\) is the vector
of monomial exponents in \(F\), including signs
and zeros, the entries of \(M\mathbf{a}\) are
therefore \(O_m(d^2)\). The integer determinant
is nonzero, and the adjugate formula with
\(|M_{ij}|\le d\) yields
\[
 \max_j |a_j|\le C'_m d^{m+1}.
\]
At least one of the nonzero \(|a_j|\) is at least
\(E\).

Consequently, for \(d<d_0=(E/C'_m)^{1/(m+1)}\)
the coordinate divisors have an integer relation
\(\sum_j u_j\operatorname{div}(x_j)=0\) with
\(\mathbf u\ne0\). The coefficients can be chosen
bounded in terms of \(m,d\) by taking minors of a
matrix of independent valuation rows. To bound
the valuation entries by \(d\) even in this case,
use the same section realizing any finite set of
prescribed rows. The rational function
\(\prod_j x_j^{u_j}\) has zero divisor and hence
is constant on the normal projective variety.
At a positive integral point, separate the
positive and negative powers. The resulting
numerator and denominator are already coprime,
so they are determined by the constant. Each
\(x_j\) with \(u_j\ne0\) divides the appropriate
fixed integer, which is at most
\(H^{\sum_j |u_j|}\).
The divisor bound therefore fixes all those
variables at a cost \(O_{m,d,\eta}(H^\eta)\),
for arbitrary \(\eta>0\). Each fixing leaves
fewer variables in the same equation (29), with
modified nonzero coefficients and the same
lower bound on the remaining exponents.

For a subvariety defined over \(\mathbb Q\) of
dimension \(k\ge1\) and degree \(d\ge d_0\), use
a finite birational linear projection to a
hypersurface in \(\mathbb P^{k+1}\) of the same
degree. The image of a point of height \(H\)
has height \(O_{m,d}(H)\), uniformly in the
subvariety. Apply the fixed-degree hypersurface
covering theorem there and pull back the
auxiliary hypersurfaces. The subvariety is
covered, on its rational points in question, by
\(O_{m,d,\eta}(H^{(k+1)/d^{1/k}+\eta})\)
collections of bounded-degree components of
smaller dimension. The pullback is proper
because the auxiliary polynomial does not
vanish on the image. Bézout bounds all new
component degrees in terms of the old degree
and \(\eta\). A geometric component that is not
defined over \(\mathbb Q\) instead has all its
rational points in its intersection with a
distinct Galois conjugate, which reduces the
dimension at no exponent cost. Zero-dimensional
components have boundedly many points.
[Castryck–Cluckers–Dittmann–Nguyen, Proposition 4.3.2 and its proof](https://msp.org/ant/2020/14-8/ant-v14-n8-p09-s.pdf)
[Salberger, Theorem 0.10](https://londmathsoc.onlinelibrary.wiley.com/doi/full/10.1112/plms.12508)

Iterate those two reductions. One either reduces
the dimension, or fixes a nonempty set of
variables and restarts with a strictly smaller
number of variables. There are at most \(m^2\)
steps along any branch, and all degrees that
occur are bounded in terms of the original
exponents and \(\eta\). Since \(k\le m-1\) and
\(d\ge d_0\), each high-degree step has exponent
at most \(m/d_0^{1/(m-1)}+\eta\).
Each low-degree step costs only \(\eta\).
Choose \(\eta\) sufficiently small in terms of
\(\varepsilon,m\), giving (30).
For bounded \(E\) the trivial estimate is absorbed
by increasing \(C_m\).

The bounded anatomic coordinates within each
letter can be taken pairwise coprime by BBLT's
decomposition. Across letters they are coprime
by primitivity. Thus, in an anatomic box, any
bounded number of coordinates of degrees tending
to infinity as \(\lambda\to0\) can be left free
together. After freezing the other coordinates,
every free height is at most \(X^\lambda\) in
the limiting box, so these coordinates cost
\(o(\lambda)\) in the count exponent.
[BBLT, Lemma 2.2 and Proposition 2.1](https://arxiv.org/html/2410.12234v2#S2)

## The small-budget coefficient \(6/13\)

For sufficiently small fixed \(\lambda>0\),
\[
 N_\lambda(X)\ll_{\lambda,\varepsilon}
 X^{(6/13)\lambda+C\lambda^{49/48}+\varepsilon}
 \tag{31}
\]
for an absolute \(C\). We use the bounded-root estimate
(30) for at most five free coordinates.

Keep the heavy threshold \(\theta_T\) used in the proof
of (28), at \(K_0=\lceil\lambda^{-1/2}\rceil\).
Suppose \(D\ge v\lambda\), where
\(v=6/13+C_1\lambda^{1/48}\), for a sufficiently large
constant \(C_1\). Put \(t=v/6\), let \(s=L_A\) denote
short mass when discussing the height gap, and write
\(h_A\) for the normalized short height.

On every letter of mass at least \(t\lambda\), the
threshold is at least
\[
 ((13/12)v-1/2)\lambda-3\lambda/K_0\gg\lambda^{49/48}.
\]
Thus there are \(O(\lambda^{-1/48})\) heavy coordinates
on each such letter and the complementary height is
\(O(\sqrt\lambda)\). On any near-tall letter in this
class choose a heavy coordinate of largest weighted
height. Its weight is \(\gg\lambda^{1/48}\), so its
degree is \(\gg\lambda^{-47/48}\). A letter with
exactly one heavy coordinate of near-tall height has
heavy degree \(\gg1/\lambda\).

On a tall letter with at least two heavy coordinates,
one may retain two including the largest-weight
coordinate and freeze all others at cost
\(L_T-2\theta_T\le L-2D+O(\lambda^{3/2})\).
The proof of the real fibre (19) works for either
retained exponent of size \(\gg\lambda^{-47/48}\) and
gap \(\tau\ge\lambda^{1/4}\). In the grouped case its
spacing test uses \(1/g\ll m\lambda^{47/48}=o(\tau)\)
with \(m\asymp\sqrt{\tau/\lambda}\); all other
constants are as in (19). The two-heavy real fibre
and the one-heavy cofactor fibre then each cost
\(O(\lambda^{9/8})\) in the count exponent, beyond
any explicitly frozen complementary mass.

### All three masses at least \(t\lambda\)

If the short gap is at least \(\lambda^{1/4}\), use
the real bounds with one complete tall term fixed.
When both tall letters have one heavy coordinate,
\(D\le\min(L_B,L_C)+O(\lambda^{9/8})\) gives
\(s\le L-2D+O(\lambda^{9/8})<t\lambda\).

When both tall letters have at least two heavy
coordinates, the real bound gives
\(s\le3L-6D+O(\lambda^{9/8})\).
Fix the short term and leave two heavy roots on each
tall letter. The bounded-root estimate accounts for
the four free coordinates, including the case where
the third summand is fixed. The frozen mass bound
gives \(5D\le2L+s+O(\lambda^{49/48})\), since the
sum of the four free masses is at least
\(4D-L-s-O(\lambda^{3/2})\).
Together these give \(11D\le5L+O(\lambda^{49/48})\).

If \(F\) is the tall letter with multiple heavy
coordinates and \(J\) has one, fixing \(F\) gives
\(L_F\ge D-O(\lambda^{9/8})\) by the real fibre.
Fixing \(J\) and freeing two roots on \(F\) instead
gives \(L_J\ge3D-L-O(\lambda^{9/8})\). Thus
\(s\le2L-4D+O(\lambda^{9/8})\).
After fixing the short term, release two heavy
coordinates on \(F\) and one on \(J\) in (30).
The frozen mass estimate yields
\(8D\le5L-2L_F-L_J+O(\lambda^{49/48})\).
Using \(L_F+L_J=L-s\) and \(L_F\ge D-O(\lambda^{9/8})\)
gives \(9D\le4L+s+O(\lambda^{49/48})\).
It follows that \(13D\le6L+O(\lambda^{49/48})\).

It remains that the three heights are \(1-O(\lambda^{1/4})\).
Every letter has at least one heavy coordinate.
For three single-heavy letters use the \(2/5\) bound
following (18). For one multiple-heavy letter \(F\),
leave two of its coordinates, including its
largest-weight one of mass \(\beta_F\), and the two
single-heavy coordinates free in (30). Their masses
obey
\[
 \beta_1+\beta_2\le
 2L-3D-L_F+O(\lambda^{49/48}),\qquad
 \beta_F\le (L+L_F)/2-D+O(\lambda^{3/2}).
\]
The first inequality uses the lower bounds
\(\theta_F\) on both released masses of \(F\), and
the second uses \(\theta_F\) on the other one.
The primitive-ratio congruence counts after fixing
the cheaper single-heavy coordinate and the
modulus coordinate of mass \(\beta_F\).
Indeed its modulus has weight \(w\gg\lambda^{1/48}\);
choose \(k\asymp w/\lambda\) and let
\(\ell\to\infty\) slowly with
\(\ell^3\ll w/\sqrt\lambda\).
The unit ratio of the two low cofactors has
height \(X^{O(\sqrt\lambda)}\), and the free single
power has degree \(\gg1/\lambda\). The scalable
independent-exponent bound therefore costs
\(w/(\ell k)=o(\lambda)\).
This gives \(D\le\min(\beta_1,\beta_2)+\beta_F+o(\lambda)\).
Adding the displayed mass inequalities yields
\(7D\le3L+o(\lambda)\).

If at least two letters have multiple heavy
coordinates, leave two on each of those letters,
and one on the third letter \(J\). Estimate (30)
counts the five free coordinates at exponent
\(O(\lambda^{49/48})\), since their degrees are
at least \(K_0\). The heavy thresholds now give
\[
 6D\le (5L+L_J)/2+O(\lambda^{49/48}).
\]
Since \(L_J\le L-D\) by the pair bound,
\(13D\le6L+O(\lambda^{49/48})\).

### A mass below \(t\lambda\)

Write \(s<t\lambda\) for this mass and fix the
corresponding complete term. The other masses are
at least \((v-t)\lambda\) and at most \((1-v)\lambda\).
Both heavy thresholds are bounded below by a
positive constant times \(\lambda\).
If the fixed term is tall and the short height gap
is at least \(\lambda^{1/4}\), use the real fibre on
the other tall term. This gives \(D\le s+L-2D+
O(\lambda^{9/8})\), which is below \(v\lambda\)
for small \(\lambda\). Otherwise the two unfixed
letters have heights \(1-O(\lambda^{1/4})\).

When one has multiple heavy coordinates,
leave two of them free and one high coordinate on
the other letter. The affine estimate (22) costs
\(O(\lambda^{3/2})\) beyond the frozen mass.
Writing \(J\) for the singly released letter, we
get
\[
 4D\le s+(3L+L_J)/2+O(\lambda^{3/2}).
\]
As \(L_J\le L-D\), this gives
\(9D\le4L+2s+O(\lambda^{3/2})\).
When both remaining letters are single-heavy, the
binary count and prime-power fibre in the proof
of (28) give \(7D\le3(L+s)+O(\lambda^{5/4})\).
For this application \(\theta_T\gg\lambda\), and the
degree comparison and all fibre margins are as
there. Both inequalities contradict \(D\ge v\lambda\)
for large enough \(C_1\), since \(s<t\lambda\),
\(9v-4-2t=26v/3-4\), and
\(7v-3-3t=13v/2-3\).

In the other cases the contradictions follow by
choosing the same \(C_1\) to dominate the absolute
constants in \(O(\lambda^{49/48})\). The finite
anatomic reduction and dyadic summation are as for
(28), proving (31).


## Rational compression with controlled exact fibers

The coprime ternary estimate can also be used for rational
coordinate families. The appropriate hypothesis before
amplification is a bound on exact fibers. Namely, for any
set of positive unit rational triples of height at most \(n^{C/k}\)
satisfying the ternary congruence, suppose at most \(F_0\)
triples give any one exact rational value of the monomial.
The unamplified argument bounds the set by
\[
 O_{C,k,p,q,r,\varepsilon}(F_0 n^{2/k+\varepsilon}).
 \tag{R1}
\]
The determinant steps already use rational heights. In the
bounded primitive-character case, reconstruction fixes
the exact primitive value, and hence the original value
in at most \(n^\varepsilon\) possibilities as before.
On a constant component the assumed exact-fiber bound
takes the place of divisor counting. For the nonconstant
components, evaluate \(J_0\) at the rational points and
clear denominators, which are units modulo \(n\). Its
nonzero numerator has height \(n^{O_C(1/k)}\) with constants
depending on the fixed exponents. The same Hensel records
determine rational coordinates because the reconstruction
threshold is \(n^{2C/k}\), below their common modulus.
This proves (R1). Rational coordinates may have arbitrary
signs of their prime valuations and may share primes.

Let \(u_1,u_2,u_3,u_4\) be four positive pairwise coprime
integer roots of height at most \(n^{C/k}\), and write their
signed exponent vector as \(a=(p,q,-r,-s)\).
Suppose that some nonzero bounded integer vector
\(w=(w_1,w_2,w_3,w_4)\) satisfies \(a\cdot w=0\).
Permute the coordinates and change the sign of \(w\) to
ensure \(w_4>0\). Define
\[
 R_i=u_i^{w_4}/u_4^{w_i}\quad(1\le i\le3),\qquad
 \prod_{i=1}^3 R_i^{a_i}
   =\left(\prod_{i=1}^4 u_i^{a_i}\right)^{w_4}.
 \tag{R2}
\]
This gives a ternary congruence with a new unit target.
Invert individual rational coordinates if needed to put its
exponent signs in the form used in (R1). The height constant
is bounded in terms of the bound on \(w\). Given an exact
value of the compressed monomial, positivity determines
\(\prod_i u_i^{a_i}\) uniquely. Pairwise coprimality of the
\(u_i\) now gives divisor-many tuples: at each prime the
valuation of the ratio assigns the prime to one coordinate
of the prescribed sign, after which that coordinate's
valuation is determined. This also bounds the multiplicity
of the map to each rational triple. In particular \(F_0\)
in (R1) is \(n^{o(1)}\) for the compressed family.

Amplify before applying (R2), by coloring the primes with
four coordinate colors. As in the ternary lemma, a coloring
retains an \(n^{-o(1)}\) proportion, and multiplicative
products have \(n^{o(1)}\) input multiplicity. The amplified
integer tuples are still pairwise coprime and satisfy the
congruence with target \(K^\ell\), so (R1) applies to their
compressed rational coordinates at height constant \(O(\ell)\)
for each fixed bound on \(w\). Undoing amplification gives
\(n^{2/(\ell k)+\varepsilon}\) input tuples, for any fixed
\(\ell\) and sufficiently large \(k\) depending on \(\ell\)
and the bound on \(w\). A large same-letter gcd is a special case,
giving a relation supported on two coordinates.

## An explicit coefficient below \(5/11\)

For every sufficiently small fixed \(\lambda>0\),
\[
 N_\lambda(X)\ll_{\lambda,\varepsilon}
 X^{(5/11-10^{-7})\lambda+\varepsilon}.
 \tag{32}
\]
The numerical gap is conservative. It comes from applying the
scalable congruence with amplification order one to the only
case that can approach \(5/11\).

Suppose the bound fails on budgets tending to zero. In the
resulting boxes take limiting ratios normalized by total mass
\(L\), which is comparable to \(\lambda\) by pair counting.
Write \(d=\lim D/L\), \(l_i=\lim L_i/L\), and \(u=1-2d\).
Then \(d\ge5/11-10^{-7}\), \(u\ge0\) by the half-budget
bound, \(\sum l_i=1\), and \(l_i\le1-d\) by pair counting.
Take the usual heavy cutoff \(K_0=\lceil\lambda^{-1/2}\rceil\).
On any letter with \(l_i>u\), the stopping threshold has
positive limiting mass \((l_i-u)/2\).
The complementary weighted height is \(O(\sqrt\lambda)\),
the number of heavy coordinates is bounded along the
sequence, and each heavy mass is at least its threshold
in the limit. A near-tall letter has a heavy root carrying
positive weighted height, hence of degree \(\asymp1/\lambda\).
If two heavy roots are retained on a letter, its complementary
mass is at most \(u\) in the limit.

Several estimates below leave a single heavy root on a free
near-tall letter and keep its complement variable. Fix
another complete term or fix its single heavy root, and
lock a heavy root of positive weighted height on the modulus
letter. The residual numerical height is \(Y^{O(\sqrt\lambda)}\);
the independent-exponent scalable bound counts this fiber
at normalized exponent zero. Indeed the modulus has
height bounded below, the free degree and the fullness
threshold are of order \(1/\lambda\), and the root height
cap follows by choosing the threshold as a sufficiently
small fixed multiple of \(1/\lambda\). Use any fixed
amplification order and then let that order grow after
taking the small-budget limit. The exact-ratio recovery
is the same as for (18). When two heavy roots on a free
letter vary and its complement is fixed, use the ternary
congruence with \(W=1\), again at normalized exponent zero.
The bounded-root estimate (30) counts any bounded selection
of free heavy roots in the original equation with exponent
\(o(\lambda)\).

### A letter of mass at most \(u\)

Fix this complete letter of mass \(s\le u\). The other
masses are at least \(d-s>u\), so their stopping thresholds
are positive. If the fixed term is tall and the short
height gap is at least \(\lambda^{1/4}\), the real fiber
(19) or its single-heavy version counts the other tall
letter after retaining at most two heavy roots. This gives
\(d\le s+u\le2u\), hence \(d\le2/5\).

Otherwise the two unfixed letters are near tall, of
height \(1-O(\lambda^{1/4})\). If both are single-heavy,
lock either of the two roots. This gives
\(d\le s+\min(\beta_1,\beta_2)\). Release both roots
in the original equation to obtain
\(d\le1-\beta_1-\beta_2\). Hence \(3d\le1+2s\le1+2u\),
so \(d\le3/7\).
If both have multiple heavy roots, release two on each.
The two frozen complementary masses cost at most \(2u\)
and the fixed letter costs \(s\), so \(d\le3u\) and again
\(d\le3/7\).

If one letter has multiple heavy roots and the other is
single-heavy, lock the largest-weight heavy root of the
multiple-heavy letter, of mass \(\beta\). Let \(f\) be a
second heavy mass on that letter, and let \(\gamma\) be the
mass of the single free root on the other. The modular count
gives \(d\le s+\beta\), and releasing all three roots in
the original equation gives \(d\le1-\beta-f-\gamma\).
The thresholds give
\(f+\gamma\ge2d-\tfrac12-\tfrac{s}{2}\).
Adding yields \(4d\le\tfrac32(1+s)\le\tfrac32(1+u)\),
so \(d\le3/7\).

### Near-tall boxes with all masses above \(u\)

Suppose now that all \(l_i>u\) and that even the short
height is \(1-O(\lambda^{1/4})\). Every letter has a
heavy root of positive weighted height. If all three
are single-heavy, the count with two roots fixed and
the count with all three released give \(d\le2/5\)
as for (18).

If just one letter \(F\) has multiple heavy roots, denote
its largest-weight heavy mass by \(\beta_F\) and another
heavy mass by \(f\). The other two roots have masses
\(\beta_1,\beta_2\). The scalable count with \(F\) as
modulus gives \(2d\le\beta_1+\beta_2+2\beta_F\).
Releasing four roots gives
\(d\le1-\beta_1-\beta_2-\beta_F-f\), and
\(\beta_F\le l_F-f\). Therefore
\[
 2d\le1-d+l_F-2f\le2-3d,
\]
which gives \(d\le2/5\).

If exactly two letters \(F,J\) have multiple heavy roots,
let \(\beta_F,\beta_J\) denote their respective largest-weight
heavy masses. Fix \(F\) in full and lock the \(J\) root,
leaving the third letter single-heavy and free. Interchanging
\(F,J\) gives
\(\beta_F+\beta_J\ge2d-l_F-l_J\).
Release two roots on each of \(F,J\) and one on the third
letter. Besides \(\beta_F,\beta_J\), the three released
masses sum to at least \(3d-1\), from the stopping
thresholds. The resulting count gives \(6d\le3-l_M\).
The usual estimate releasing five roots, using the
threshold on all five, gives \(12d\le5+l_M\).
Adding gives \(d\le4/9\).
If all three letters have multiple heavy roots, release
two per letter and obtain \(d\le3u\), or \(d\le3/7\).

### A short-height gap with all masses above \(u\)

It remains that the short height gap is at least
\(\lambda^{1/4}\); write \(s\) for the short mass.
When both tall letters are single-heavy, fixing either
complete tall term and counting the other by the real
fiber gives \(s\le u\).

If one tall letter is multiple-heavy and the other is
single-heavy, the real counts give \(s\le2u\).
In detail the multiple-heavy letter has mass at least
\(d\), by counting the single-heavy fiber with that term
fixed. The single-heavy letter has mass at least \(d-u\)
by the reverse count, allowing the multiple-heavy
complement. Now fix the short term, lock a root on the
multiple-heavy letter, and free the single-heavy tall
letter. Exactly as in the small-letter mixed case, the
two counts give \(4d\le\tfrac32(1+s)\le\tfrac32(1+2u)\),
so \(d\le9/20\).

If both tall letters have multiple heavy roots, the real
counts give \(l_B,l_C\ge d-u=3d-1\), hence \(s\le3u\).
Fixing the short term and releasing two roots on each
tall letter gives \(d\le s+2u\). In particular \(d\le5/11\).
If the short letter has a heavy root, release it as well.
The five thresholds give \(12d\le5+s\le5+3u\), hence
\(d\le4/9\). The surviving short letter has no heavy
root and has weighted height \(O(\sqrt\lambda)\).

### Clustering in the remaining case

Choose two heavy roots on each tall letter, including
a largest-weight root. Let \(H\) be the sum of their
limiting masses. Releasing these four roots in the original
equation gives \(H\le 1-d\). Write
\(\nu=5/11-d\), so \(0\le\nu\le10^{-7}\).
The thresholds and \(l_B,l_C\ge3d-1\) give, for each
of the four masses \(\beta_i\),
\[
 \frac3{22}-\frac52\nu\ \le\ \beta_i\
 \le\ \frac3{22}+\frac{17}{2}\nu .
 \tag{33}
\]
No additional heavy root can occur. An extra heavy mass
would be at least \(3/22-(5/2)\nu\), whereas the total
slack between \(1-d\) and the four lower bounds is
at most \(11\nu\). Releasing five roots by (30) would
contradict that slack. Thus all the complementary
coordinates have weighted height \(O(\sqrt\lambda)\).

Choose on each tall letter its lighter root by weighted
height. Both lighter roots have limiting weight at most
\(1/2\). If either limit is less than \(1/2\), release
these two roots and all the complementary coordinates in
the geometry bound (9). The total released weight is
bounded below one, and the total mass left fixed is at most
\[
 H-\left(2d-\tfrac12-\tfrac{s}{2}\right)
 \ \le\ \tfrac32-3d+\tfrac{s}{2}
 \ \le\ 3-6d<d,
\]
a contradiction. All four heavy weighted heights are
therefore \(1/2+o(1)\).

Take \(q=\lfloor 11/(3L)\rfloor\) and use the partial-power
compression (20). Lock one heavy root on a tall letter.
On the other tall letter take multiplier one on both heavy
roots and zero on its complement. Take multiplier zero
on the short letter. In the notation of (20), (33) gives
\[
 \limsup E
 \le 2\frac{11}{3}\frac{17}{2}\nu
 \le \frac{187}{3}\,10^{-7}.
\]
For this application one can use the exact scalable
parameters rather than the implied constant in (20).
Take \(k=\lfloor2/L\rfloor\) and amplification order
\(\ell=1\). The locked root has weighted height
\(w=1/2+o(1)\) and degree at least \(k\), by (33).
The compressed root has numerical height exponent at most
\((3/11+17\nu+o(1))L\), below \(2w/k\).
The numerical residual exponent \(E+o(1)\) is strictly
below \(w/16384\). The free power \(q\) exceeds
\(\lfloor k/256\rfloor\), and \(k\ge1024\) for small
enough budgets. The scalable lemma is thus applicable.
Its cost in normalized units is at most \(1/4+o(1)\).
Adding the locked mass gives
\[
 d\le\frac3{22}+\frac{17}{2}\nu+\frac14
     <5/11-10^{-7},
\]
a contradiction in the final case.

The limiting arguments can be applied to a sequence of
finite anatomic bounds with losses \(o(\lambda)\), as in
(31). At each fixed truncation the degrees are finite,
all fiber estimates are uniform in the frozen data, and
their constants enter the implied constant before letting
\(X\) tend to infinity. A failure of (32) for arbitrarily
small budgets would give exactly the sequence above.

## Residual radical bounds for hyperbolic prime-power signatures

Let \(p\) be prime, and let \(e_A,e_B,e_C\) be fixed positive
powers of \(p\) with
\[
 \frac1{e_A}+\frac1{e_B}+\frac1{e_C}<1.
\]
For a primitive positive triple \(a+b=c\), define
\[
 R(a,b,c)=\prod_{T\in\{a,b,c\}}\ \prod_{\ell\mid T,\ e_T\nmid v_\ell(T)}\ell ,
\]
where \(e_T\) is the assigned exponent on that term. Then for
some constant \(J\) depending only on the signature, the total
number of such triples with \(R(a,b,c)\le Y\) is
\(O_{e_A,e_B,e_C}(Y(1+\log Y)^{J-1})\) for \(Y\ge1\).
The count has no height cutoff.

### A three-point cover with a \(p\)-group

Sort the desired exponents as \(p^{j_1}\le p^{j_2}\le p^{j_3}\).
Take \(n=\max(3,p^{j_3-1}+1)\), so \(p^{j_3-1}<n\le p^{j_3}\)
in the hyperbolic cases. In the group of upper unitriangular
\(n\)-by-\(n\) matrices over \(\mathbb F_p\), put
\[
 A_0=I+\sum_{i\text{ odd}}E_{i,i+1},\qquad
 B_0=I+\sum_{i\text{ even}}E_{i,i+1},
\]
with summation indices \(1\le i<n\). The two nilpotent
parts square to zero, so \(A_0,B_0\) have order \(p\).
Every first-superdiagonal entry of \(A_0B_0-I\) is \(1\),
so its nilpotence index is \(n\). Consequently \(A_0B_0\)
has order \(p^{j_3}\). In the direct product with cyclic
groups of orders \(p^{j_1},p^{j_2}\), take \(A\) to be
\(A_0\) combined with a generator of the first cyclic
group and \(B\) to be \(B_0\) combined with a generator
of the second. Then \(A,B,(AB)^{-1}\) have the three
prescribed orders. They generate a finite \(p\)-group \(G\).

The Riemann existence theorem applied to the regular
permutation action of \(G\) gives a connected Galois cover
of the projective line with these branch cycles at
\(0,1,\infty\), after permuting branch labels to fit
the assignment. Belyi's theorem gives a model of the
cover over a number field. Enlarge that fixed field
to define all deck transformations, and call the
resulting cover \(f:C\to\mathbb P^1_K\). Riemann–Hurwitz
gives \(2g(C)-2=|G|(1-1/e_A-1/e_B-1/e_C)>0\).
[Belyi, Theorem 4](https://www.mathnet.ru/eng/im1682)

### Rank over bounded-degree Galois \(p\)-extensions

Let \(J/K\) be a fixed abelian variety of dimension \(g\)
over a number field. Enlarge \(K\) once to contain
\(\mu_p\) and \(J[p]\). Fix \(M\ge1\).
If \(F/K\) is a Galois \(p\)-extension of degree at
most \(M\) unramified outside a set \(S\) of finite
primes of \(K\), then
\[
 \operatorname{rank} J(F)\ll_{J,K,p,M}1+|S|.
\]
Enlarge \(S\) by the fixed primes above \(p\) and the bad
reduction primes of \(J\). Let \(S_F\) be the primes above
\(S\), and let \(V_F\subset F^\times/F^{\times p}\) be
the vector space of classes of valuations divisible by
\(p\) at every prime outside \(S_F\).

Any \(m\) independent classes in \(V_F\) define, by Kummer
theory, an elementary abelian extension \(L/F\) of rank
\(m\) unramified outside \(S_F\). Indeed at a prime away
from \(p\), after removing a uniformizer \(p\)-th power,
a \(p\)-th root of a unit is unramified, by separability
over the residue field and Hensel lifting. Take the
compositum \(E\) of the conjugates of \(L\) over \(K\).
As \(F/K\) is Galois, \(E/F\) is again an elementary
abelian \(p\)-extension, so \(G=\operatorname{Gal}(E/K)\)
is a finite \(p\)-group unramified outside \(S\).
The \(\mathbb F_p\)-dimension of
\(\operatorname{Hom}(G,\mathbb Z/p)\) is bounded by that
of \(V_K\), since cyclic degree-\(p\) subextensions are
unramified outside \(S\) and \(K\) contains \(\mu_p\).
Here \(V_K\) uses the same valuation condition with
respect to \(S\). Mapping a class to its divisor
outside \(S\) divided by \(p\) gives an exact sequence
with the \(S\)-unit classes modulo \(p\)-th powers as
kernel and the \(p\)-torsion of the \(S\)-class group
as image. Dirichlet's unit theorem gives unit rank
\(O_K(1+|S|)\), and the \(S\)-class group is a quotient
of the fixed class group of \(K\). Thus
\(\dim_{\mathbb F_p}V_K=O_K(1+|S|)\).

The Frattini quotient of a finite \(p\)-group identifies
its minimum number of generators with the dimension of
its \(\mathbb Z/p\)-character space. Schreier's bound
then gives \(O_{K,M}(1+|S|)\) generators of
\(\operatorname{Gal}(E/F)\), a subgroup of index at
most \(M\). This group maps onto \(\operatorname{Gal}(L/F)\),
so \(m\ll_{K,M}1+|S|\) and \(\dim V_F\) has that bound.

The Kummer map for multiplication by \(p\) injects
\(J(F)/pJ(F)\) into
\(H^1(F,J[p])\cong(F^\times/F^{\times p})^{2g}\),
since \(J[p]\) and \(\mu_p\) are defined over \(F\).
Outside \(S_F\) the canonical local class is unramified.
The abelian variety has a good proper smooth model,
and division by \(p\) of a local section is finite
étale because the residue characteristic differs from
\(p\). The image therefore lies in \(V_F^{2g}\).
Mordell–Weil and the bound on \(\dim V_F\) prove the
rank estimate.

### A uniform count of fields and points

Under the same hypotheses, the number of Galois
\(p\)-extensions of \(K\) of degree at most \(M\)
inside a fixed algebraic closure and unramified outside
\(S\) is at most \(C^{1+|S|}\), for fixed \(M,K,p\).
Indeed, take any finite collection of these extensions
and form its compositum. Its finite Galois group is
again a \(p\)-group unramified outside \(S\), and the
same bound on the base character dimension shows it
has \(O_K(1+|S|)\) generators. Each of the collected
fields corresponds to a distinct normal subgroup of
index \(j\le M\), the kernel of a homomorphism into
the symmetric group on \(j\) letters via the coset
action. There are at most \((j!)^{O_K(1+|S|)}\) such
homomorphisms. Summing for \(j\le M\) gives the claimed
bound for every finite collection, hence for all fields.

Take \(J=\operatorname{Jac}(C)\), and enlarge \(K\)
once to contain the base point of an Abel–Jacobi
embedding, \(J[p]\), and \(\mu_p\).
For an abc triple the fibre at \(a/c\) is a \(G\)-torsor
over \(K\). The field \(F\) of any geometric point of
that fibre is a Galois \(p\)-extension of degree at
most \(|G|\). Indeed the absolute Galois action on
the geometric fibre commutes with the simply
transitive deck action, so after choosing a point the
Galois image identifies with a subgroup of \(G\) acting
freely on the fibre. The stabilizer of a single point
is the kernel of this action.

At primes of \(K\) outside a fixed finite exceptional
set, Beckmann's specialization theorem computes the
inertia in \(F/K\) as the branch inertia generator
raised to the local intersection multiplicity, or
gives trivial inertia if the section meets no branch.
The specialization \(a/c\) meets the branch at \(0\),
\(1\), or \(\infty\) according as the prime divides the
corresponding one of \(a,b,c\). Primitivity prevents
two of these cases at the same prime. Over a rational
prime \(\ell\) the local intersection multiplicity is
\(v_\ell(T)\) times the ramification index in the
fixed field \(K/\mathbb Q\). Thus the specialized
inertia is trivial if \(e_T\mid v_\ell(T)\). The only
possibly ramified primes of \(F/K\), outside the fixed
exceptional set, lie over \(R(a,b,c)\).
[Beckmann's specialization theorem, stated as Theorem 0.1 in Conrad](https://math.stanford.edu/~conrad/papers/inertiapaper.pdf)

Fix the squarefree residual radical \(R\) and use the
set of primes of \(K\) above \(R\) together with the
fixed primes needed for the cover and the rank lemma.
Its cardinality is \(O_K(1+\omega(R))\). The field
count bounds the number of choices of \(F\) by
\(C_1^{1+\omega(R)}\). For each such field the rank
lemma bounds \(\operatorname{rank}J(F)\) by
\(O_{K,G,C}(1+\omega(R))\). Gao–Ge–Kühne's uniform
Mordell–Lang theorem gives an exponential bound in
this rank on \(\#C(F)\), since a curve of genus at
least two in its Jacobian contains no positive-dimensional
abelian coset. Each point of \(C(F)\) gives a single
ratio \(a/c\), and each reduced ratio determines the
primitive triple. The number of triples with fixed
\(R\) is therefore \(O(C_2^{\omega(R)})\). Choose a fixed integer \(J_0\ge C_2\). Since \(R\) is squarefree,
\(C_2^{\omega(R)}\le d_{J_0}(R)\), the ordered \(J_0\)-fold divisor
function. The elementary bound
\(\sum_{n\le Y}d_{J_0}(n)\le Y(1+\log Y)^{J_0-1}\) follows by summing
over the first \(J_0-1\) factors and bounding the remaining count
by \(Y\) divided by their product. This proves the result.
[Gao–Ge–Kühne, Theorem 1.1](https://arxiv.org/pdf/2105.15085v4)

In exact prime-valuation box masses \(\alpha_{T,i}\)
the resulting exponent is
\[
 D\le\sum_{T\in\{a,b,c\}}\ \sum_{e_T\nmid i}\alpha_{T,i}
\]
up to the arbitrarily small exponent loss. The cover is
a quantitative instance of the three-point-cover method
of Darmon and Granville. The ramification filter and
the uniform rank bound replace a fixed set of allowed
primes by a bound on their product, for these signatures.
A tail coordinate can use the bound if its actual
valuation classes have been retained.
[Darmon–Granville, §3](https://www.math.mcgill.ca/darmon/pub/Articles/Research/12.Granville/pub12.pdf)

For the finite projected program one may discard the
saving from the unobserved valuation classes. In BBLT's
Lemma 2.2 only the coordinate indexed by the cutoff
\(k_0=\lceil1/\varepsilon_0\rceil\) receives contributions
from larger true valuations. The six visible masses can
therefore be taken from exact valuations by making
\(\varepsilon_0\) small. On actual primitive triples
represented by the resulting anatomic boxes the projected
row is
\[
 D+\sum_T\sum_{\substack{i\le6\\e_T\mid i}}\alpha_{T,i}\le L ,
\]
where \(L\) is the total projected mass and the error
vanishes on passing to limiting profiles. The previous
box bounds remain valid for this subset of triples.
In the six-coordinate program the useful minimal
signatures are \((4,4,4)\), permutations of \((2,4,8)\)
and \((3,3,9)\), and \((5,5,5)\). The row counts
primitive triples with their actual valuation classes.
[BBLT v2, Lemma 2.2](https://arxiv.org/html/2410.12234v2)

## Residual radical bounds for two Euclidean signatures

Let \((e_A,e_B,e_C)\) be \((3,3,3)\) or any permutation of
\((2,4,4)\). For primitive positive triples \(a+b=c\) with
these assigned exponents, put
\[
 R=\prod_{T\in\{a,b,c\}}\ \prod_{\substack{\ell\mid T\\
                                e_T\nmid v_\ell(T)}}\ell .
\]
For \(X\ge2\), \(Y\ge1\) and \(\varepsilon>0\), the number
of triples with \(c\le X\) and \(R\le Y\) is
\[
 \ll_\varepsilon YX^\varepsilon .
\]
Endomorphism descent on the associated elliptic curves bounds
the number of points in the finite height range.

### Maps to elliptic curves

For fixed assigned exponents, write each term as its unique
power-free coefficient times a full \(e_T\)-th power.
The coefficients are pairwise coprime in absolute value,
and their product has radical \(R\). At each prime of \(R\)
there are at most \(\sum_T(e_T-1)\) choices of its term and
nonzero coefficient exponent.

For the cubic signature write \(a=A x^3,\ b=B y^3,\
c=C z^3\), with positive \(A,B,C,x,y,z\). Put
\[
 u=-4AB\,\frac{xy}{z^2}, \qquad
 v=-4ABC\,\left(2\frac ac-1\right).
\]
The identity \((2a/c-1)^2=1-4ab/c^2\) shows that
\[
 v^2=u^3+(4ABC)^2.
\]
For each fixed \(A,B,C\), the coordinate \(v\) recovers
\(a/c\), whose reduced numerator and denominator give
\(a,c\) since the triple is primitive. Thus the map to
this elliptic curve is injective on the triples with
these coefficients.

For the quartic signature, rearrange the equation as
\(A x^2+B y^4=C z^4\), allowing signs on the fixed
coefficients and keeping positive \(x,y,z\). Write
\(t=y/z,\ s=A x/z^2\). Then \(s^2=AC-ABt^4\), and
\[
 u=-ABt^2,\qquad v=-ABts,\qquad
 v^2=u^3-A^2BCu .
\]
The coordinate \(u\) determines \(t>0\), hence the ratio
of the two complete quartic terms in absolute value.
Primitivity determines both terms from this ratio,
and the original equation determines the quadratic term.

In both cases the image has \(h(u),h(v)=O(\log X)\)
uniformly in the power-free coefficients. The discriminant
of the displayed Weierstrass model is supported on the
primes of \(6R\).

### Uniform rank and point bounds for the two twist families

For the curve \(E:v^2=u^3+w^2\) take \(K=\mathbb Q(\zeta_3)\)
and the automorphism \(\tau(u,v)=(\zeta_3u,v)\).
The degree-three endomorphism \(\varphi=1-\tau\) has
kernel consisting of the point at infinity and \((0,\pm w)\).
For \(E:v^2=u^3+du\) take \(K=\mathbb Q(i)\) and
\(\tau(u,v)=(-u,iv)\). Then \(\varphi=1-\tau\)
has degree two and kernel consisting of the point at
infinity and \((0,0)\). In either case the kernel
has order \(p\) and is pointwise rational over \(K\),
where \(p=3\) or \(2\), respectively. The degree
follows also from the listed fixed points since the
endomorphism is separable in characteristic zero.

The rank \(r\) of \(E(\mathbb Q)\) is the dimension
of the \(+1\)-eigenspace for conjugation on
\(V=E(K)\otimes_{\mathbb Z}\mathbb Q\). Indeed an
invariant in \(V\) is the average of a vector and
its conjugate, and the trace of each point lies in
\(E(\mathbb Q)\). Complex conjugation sends \(\tau\) to \(\tau^{-1}\).
The nonzero CM endomorphism \(\tau-\tau^{-1}\) is invertible on \(V\) and gives an
isomorphism of the two eigenspaces. Consequently \(V\)
has dimension \(r\) over the CM field \(K\), and the
determinant of \(\varphi\) on the free Mordell–Weil
lattice has absolute value
\(\lvert N_{K/\mathbb Q}(1-\tau)\rvert^r=p^r\),
where \(\tau\) is identified with the corresponding
root of unity in the CM algebra.

Kummer descent injects \(E(K)/\varphi E(K)\) into
\(H^1(K,E[\varphi])\cong K^\times/K^{\times p}\).
This can be seen by choosing a \(\varphi\)-division
point and taking its Galois differences in the kernel.
At a place of good reduction of residue characteristic
different from \(p\), the image is unramified.
For a division point \(Q\), each inertia element
\(\sigma\) satisfies \(Q^\sigma-Q\in E[\varphi]\).
The difference reduces to zero, since inertia acts
trivially on the residue field. Reduction is injective
on torsion of order prime to the residue characteristic,
so the difference vanishes. This is the reduction
argument in Silverman, Proposition VIII.1.5(b).
A local Kummer class unramified at such a place
has valuation divisible by \(p\), since valuations
in an unramified extension have the same value group.
Hence the global images lie in the space of classes
of valuations divisible by \(p\) off the primes
over \(6R\). Both fixed fields \(K\) have class
number one, so this space is the group of \(S\)-units
modulo \(p\)-th powers. Its dimension is at most
\(2\omega(R)+O(1)\), by the unit theorem and the
bound of two primes above each rational prime.
As the quotient of the free lattice has size
\(p^r\), this proves
\[
 \operatorname{rank} E(\mathbb Q)\le 2\omega(R)+O(1).
\]
[Silverman, The Arithmetic of Elliptic Curves, Proposition VIII.1.5 and the proof of Proposition VIII.1.6](https://www.math.ens.psl.eu/~obenoist/refs/Silverman.pdf)

Each model is isomorphic over a field of degree
at most six over \(\mathbb Q\) to one of the two
fixed models \(v^2=u^3+1\) and \(v^2=u^3+u\):
use a sixth or fourth root of the nonzero
coefficient. The coordinate scaling on \(u,v\)
has height \(O(\log X)\). On the fixed models,
the canonical height equals \(h(u)/2+O(1)\)
and is invariant under the isomorphism. Invariance
follows from the defining limit along doubling because
scaling \(u\) changes its logarithmic height by
a bounded amount for that single isomorphism.
Thus the image points have canonical height \(O(\log X)\).
Northcott's theorem on the two fixed models gives a
uniform positive lower bound on the canonical height
of nontorsion points of degree at most six and
a uniform finite bound on the number of such torsion
points. For the positive lower bound, take the
minimum among the nontorsion points of canonical
height at most one, if this set is nonempty.

The canonical height defines a positive definite
quadratic form on the free Mordell–Weil lattice.
The image set lies in a ball of radius \(O(\sqrt{\log X})\)
for the associated Euclidean norm, and distinct
lattice points are separated by an absolute positive
distance. Pack balls of half that separation around
them. The volume ratio, with the bounded torsion
multiplicity restored, gives at most
\((C_0\sqrt{\log(X+2)})^r O(1)\) points for an
absolute \(C_0>1\).
[Silverman, Theorems VIII.5.11 and VIII.9.3 and Proposition VIII.9.6](https://www.math.ens.psl.eu/~obenoist/refs/Silverman.pdf)

### Summation over residual radicals

Combining the coefficient count with the rank bound,
the number of triples for a fixed squarefree \(R\)
is \(\ll (\log(X+2))^{O(1)}Z^{\omega(R)}\), where
\(Z=C_1\log(X+2)\) and \(C_1\) is an absolute
constant. For fixed \(\delta>0\), put \(\sigma=1+\delta\).
Rankin's bound gives
\[
 \sum_{\substack{R\le Y\\R\ {\rm squarefree}}}
   Z^{\omega(R)}\le
 Y^\sigma\prod_{\ell}
     \left(1+\frac{Z}{\ell^\sigma}\right).
\]
Split the logarithm of this product at \(Z^{1/\sigma}\).
The smaller primes cost at most \(O(Z^{1/\sigma}\log Z)\).
For the larger ones, summing \(Z/n^\sigma\)
over the integers above the cutoff bounds the
cost by \(O_\delta(Z^{1/\sigma})\).
Thus the product is \(X^{o(1)}\) for fixed \(\delta>0\).
As \(R\le abc\le X^3\), one may take \(Y\le X^3\).
With \(\delta\) sufficiently small in terms of
\(\varepsilon\), this proves the count.

In the finite projected program this yields
\[
 D+\sum_T\ \sum_{\substack{i\le 6\\e_T\mid i}}\alpha_{T,i}
 \le L
\]
for the two Euclidean signatures, by the same
exact-valuation extraction as for the hyperbolic
signatures. These rows dominate the prior
\((3,3,9)\) and \((2,4,8)\) permutations, respectively,
and the \((2,4,4)\) permutations also dominate
\((4,4,4)\). The \((5,5,5)\) hyperbolic row remains
a separate input.

## A small-budget coefficient of \(3/7\)

For every \(\delta>0\) and all sufficiently small fixed \(\lambda>0\),
\[
N_\lambda(X)\ll_{\lambda,\delta,\varepsilon}
X^{(3/7+\delta)\lambda+\varepsilon}.
\]

Use the limiting notation of the preceding explicit coefficient
proof. If the displayed bound fails on budgets tending to zero,
there is a limiting box with \(d>3/7\), total normalized mass one,
letter masses \(l_T\), and \(u=1-2d\). Pair counting gives
\(l_T\le1-d\). On a letter with \(l_T>u\), every heavy mass is
at least \((l_T-u)/2>0\), there are only boundedly many heavy
coordinates, and the weighted height of the complement is
\(O(\sqrt\lambda)\). Put \(S_T\) for the sum of all heavy masses on
such a letter and \(\beta_T\) for the mass of any heavy root
carrying positive limiting weighted height when one is chosen.
If a letter has just one heavy root, \(S_T=\beta_T\).
The bounded-root estimate counts any bounded collection of
released heavy coordinates of degrees at least
\(\lceil\lambda^{-1/2}\rceil\) at exponent \(o(\lambda)\).

We use two modular operations. A fixed heavy root of positive
limiting weighted height gives a modulus \(n\) of positive
limiting height, with fullness threshold of order \(1/\lambda\).
If the main root of a free term has positive limiting weighted
height, and both the rest of that term and the unfixed part of
the other term have numerical heights \(X^{o(1)}\), the
independent-exponent scalable congruence gives normalized cost
zero beyond the frozen coordinates. Indeed, the residual ratio
\(V\) has height \(X^{o(1)}\). For every fixed amplification order
its height is eventually within the residual cap, and the
modulus weight and degrees permit the required root height
cap with threshold a suitable constant times \(1/\lambda\).
Taking the amplification order large gives the claimed limiting
cost. This applies with more than one heavy coordinate on the
partly fixed term: fix *all* its heavy roots, at cost \(S_T\),
leaving only its low cofactor variable.

If both unfixed cofactors have numerical heights \(X^{o(1)}\),
their reduced ratio can instead be recovered directly. The
congruence fixes the ratio modulo \(n\), by primitivity. Its
height is \(X^{o(1)}\), so clearing denominators shows that two
distinct such ratios cannot agree modulo a modulus of positive
limiting height. The recovered exact ratio, together with the
fixed factors, gives the ratio of the two complete terms.
These terms are coprime, so their complete values and then
the third term are determined. The same exact-ratio recovery
applies to the scalable count.

### A letter of mass at most \(u\)

Fix it in full at cost \(s\le u\). The other masses are at least
\(d-s>u\). If the fixed term is tall and the short height gap is
at least \(\lambda^{1/4}\), the real fibre bound (19) or its
single-heavy version counts the other tall letter. That gives
\(d\le s+u\le2u\), so \(d\le2/5\).

Otherwise both unfixed letters are near tall. If they are
single-heavy, locking either of their roots gives
\(d\le s+\min(\beta_1,\beta_2)\); releasing both gives
\(d\le1-\beta_1-\beta_2\). Hence \(3d\le1+2s\le1+2u\), so
\(d\le3/7\). If both have multiple heavy roots, release two
on each. The complement on each costs at most \(u\), so
\(d\le s+2u\le3u\), again giving \(d\le3/7\).

In the mixed case lock a positive-weight heavy root of mass
\(\beta\) on the multiple-heavy letter, allowing its second
heavy mass \(f\) to remain fixed only in the separate
three-root release. The single-heavy root on the other
letter has mass \(\gamma\). The counts give
\(d\le s+\beta\) and \(d\le1-\beta-f-\gamma\).
Since \(f+\gamma\ge 2d-(1+s)/2\), their sum gives
\(4d\le3(1+s)/2\le3(1+u)/2\). Thus \(d\le3/7\).

### All masses above \(u\), with positive limiting short height

All letters now have a positive-weight heavy root. If all
are single-heavy, fix two roots, using either as modulus, and
count the third by the scaled congruence. Adding the three
two-root estimates to twice the three-root release yields
\(5d\le2\). If exactly one letter \(F\) is multiple-heavy,
choose a positive-weight root of mass \(\beta_F\) and another
heavy root of mass \(f\) on it. Write \(\beta_1,\beta_2\) for the
singleton masses. The two modular estimates give
\(2d\le\beta_1+\beta_2+2\beta_F\). Releasing these four roots
and using \(\beta_F\le l_F-f\) gives
\(2d\le1-d+l_F-2f\le2-3d\), since
\(f\ge(l_F-u)/2\). Hence again \(d\le2/5\).

If exactly two letters \(F,J\) are multiple-heavy and \(M\)
is singleton, lock a positive-weight root on \(J\) and fix all
heavy roots on \(F\), leaving the singleton on \(M\) free.
Exchange the roles of \(F,J\). This gives
\[
2d\le S_F+S_J+\beta_F+\beta_J.
\]
Each of \(F,J\) has another heavy root, so
\(\beta_F+\beta_J\le S_F+S_J-(1-l_M-2u)/2\).
Releasing all heavy roots gives
\(S_F+S_J+\beta_M\le 1-d\), with
\(\beta_M\ge(l_M-u)/2\). Combining,
\(4d\le3/2+2u-l_M/2\), hence
\(8d\le(7-l_M)/2\le(7-u)/2=3+d\).
Thus \(d\le3/7\).
If all three letters are multiple-heavy, release two roots
on each. Their frozen complements cost at most \(3u\), so
\(d\le3u\) and \(d\le3/7\).

### All masses above \(u\), with short height tending to zero

Call a tall letter effectively single if it has only one
heavy root of positive limiting weighted height. Its entire
complement then has numerical height \(X^{o(1)}\),
because the number of heavy roots is bounded. If both tall
letters are effectively single, either main root serves as
modulus while the other main root, its cofactor, and the
short term are free. The scaled count gives
\(d\le\min(\beta_B,\beta_C)\), and releasing both main roots
in the original equation gives \(d\le1-\beta_B-\beta_C\).
Thus \(d\le1/3\). If only one tall letter is effectively
single, either of two positive-weight roots on the other
can be a modulus. The same count costs at most half the
mass of that letter, so \(3d\le1\).

Both tall letters may therefore be taken to have at least
two positive-weight heavy roots. The short height gap is
bounded below. With either whole tall term fixed, retain
two heavy roots on the other tall term, including a
positive-weight root, and use (19). The two real counts give
\(l_B,l_C\ge d-u=3d-1\). Every tall heavy mass is consequently
at least \(a=(5d-2)/2\). Write \(H=S_B+S_C\). Releasing all
tall heavy roots gives \(H\le1-d\).

Fix one positive-weight heavy root on one tall term as
modulus and freeze *all* heavy roots of the other tall term.
Its remaining low cofactor and the whole short term have
heights \(X^{o(1)}\). The direct ratio recovery shows that
the cost is only the modulus mass and the other term's
frozen heavy mass. If at least one tall letter has at least
three heavy roots, choose it as the modulus letter. The
other two of its heavy masses cost nothing, so
\(d\le H-2a\le1-d-2a\), giving \(7d\le3\).
If both tall letters have exactly two heavy roots, all four
have positive weighted height. The same count allows
any one of them to be the omitted root, by locking its
partner and freezing the two roots on the opposite
letter. Thus \(d\le H-\max\beta\le3H/4\le3(1-d)/4\),
again giving \(d\le3/7\).

This contradicts the chosen limiting exponent in every case.
The finite-truncation and dyadic passage is the same as in
(31) and (32). At each fixed truncation all chosen degrees
are fixed in the estimates; the modular constants are
uniform in the frozen data.

## Excluding the limiting coefficient \(3/7\)

There are absolute \(\eta>0\) and \(\lambda_0>0\) such that
\[
 N_\lambda(X)\ll_{\lambda,\varepsilon}
 X^{(3/7-\eta)\lambda+\varepsilon}
 \qquad (0<\lambda<\lambda_0).
\]
Indeed, if there is no such gap, the previous result produces a
limiting sequence with \(D/L\to d=3/7\), and \(u=1/7\).
It suffices to exclude equality in each of the three cases.

When a modulus of positive limiting weight is needed on a
letter of total weight \(h>0\), the part in degrees at least
\(k_0=\lfloor c/\lambda\rfloor\) has positive limiting weight
for sufficiently small fixed \(c>0\). This follows from the
bound \(L\le\lambda\) on total coordinate mass. At a modulus
of this kind the scaled independent-exponent estimates
apply at normalized exponent zero whenever the residual
height is \(X^{o(1)}\). In particular, when two heavy powers
on a letter have exponents differing by \(o(1/\lambda)\),
multiply their roots and use either exponent as the common
power; the numerical height of the residual is \(X^{o(1)}\).
These estimates still apply if the fixed factors on the
other free letter leave a variable low cofactor.

### Positive short height and all masses above \(u\)

The only equality case left in this range has three
multiple-heavy letters. Release two heavy roots per letter.
Equality in \(d\le3u\) forces each of their two masses to
be \(\beta_T=(l_T-u)/2\). A third heavy root on any
letter would add positive mass to this release by (30)
and cut equality. Thus each letter has exactly two heavy
roots, and the total high mass \(H=\sum_T2\beta_T\)
is \(4/7\).

Suppose some letter \(T\) has a lighter heavy root of
limiting weight less than \(1/2\). Release it and the
lighter heavy root of another letter \(U\), along with
all three low complements. The released weight is strictly
less than one. The geometry bound (9) therefore gives
\(d\le H-\beta_T-\beta_U\).
At equality this requires
\(\beta_T+\beta_U=(1-l_V-2u)/2\le1/7\), or \(l_V\ge3/7\),
where \(V\) is the remaining letter.
Choosing either other letter as \(U\) forces both their
masses to be at least \(3/7\), contrary to \(l_T>1/7\).
All six weighted heights must therefore be \(1/2+o(1)\).

The two heavy degrees on each letter now differ by
\(o(1/\lambda)\), since the two masses of that letter
have the same positive limit. Choose \(F,M\) to be
the letters of smallest and second smallest \(\beta_T\).
Lock one heavy root on \(M\), fix both heavy roots on \(F\),
and use the scaled count on the product of the two
heavy roots of the remaining letter \(J\).
The free exponent is the common approximate degree of the
two roots on \(J\). The residual ratio involving the
low cofactors has height \(X^{o(1)}\). This gives
\(d\le2\beta_F+\beta_M\le\sum_T\beta_T=2/7\).

### Short height tending to zero and all masses above \(u\)

As in the \(3/7\) proof, both tall letters must have at
least two positive-weight heavy roots. Each tall heavy
mass is at least \(a=1/14\), and their total \(H\le4/7\).

If one tall letter has at least three heavy coordinates,
use it as the modulus letter and freeze all the heavy
roots of the other tall term. Equality in \(d\le H-2a\)
forces \(H=4/7\), exactly three heavy coordinates on
the modulus letter, and its two nonmodulus heavy masses
to be \(a\). Switch between two different positive-weight
roots on that letter as modulus; this proves that all
three masses there are \(a\). If both tall letters
have three heavy coordinates, the same argument on the
other letter gives \(H=6a<4/7\). Otherwise the other
letter has two heavy coordinates with total mass
\(H-3a=5/14\). Use the partner of the larger of those
two as modulus. Both have positive weight. The cost is
at most \(H-5/28=11/28<3/7\).

With exactly two heavy coordinates on each tall letter,
all four masses must be \(1/7\), since equality in
\(d\le H-\max\beta\le3H/4\) forces \(H=4/7\) and all
four masses equal. If the sum of the weights of the
two lighter tall roots were strictly below one, release
them by geometry along with all the low cofactors
and the full short term. This gives \(d\le2/7\).
Thus the two lighter weights are each \(1/2+o(1)\).
The degrees within either tall pair coalesce; compress
one pair into a common power, while locking one heavy
root of the other tall term. The full short term has
height \(X^{o(1)}\), so the scalable count yields
\(d\le1/7\).

### A letter of mass at most \(u\)

In the remaining range \(s=u=1/7\) is necessary for
equality, and the two other terms are near tall,
as otherwise the real-fibre bound already saves.
If both other letters are single-heavy, their two
heavy masses are both \(2/7\) by the two modular
bounds and the exact two-root release. Their heavy
degrees therefore coalesce. If the small letter
has positive limiting height, fix it and use its
high-degree modulus part to count the two free
tall roots compressed into a common power ratio.
The low residual has height \(X^{o(1)}\), giving
\(d\le s=1/7\). If its height tends to zero, lock
one tall root and use the other with the unfixed
small term in the low residual ratio, giving \(d\le2/7\).

When one tall letter is multiple-heavy and the other
single-heavy, the previous mixed proof forces the
locked heavy mass \(\beta\) to be \(d-s=2/7\).
The second heavy of the same letter has mass
\(f=(l_F-u)/2\le3/14\), and the released singleton
has \(\gamma=(l_J-u)/2\). No other heavy root is
possible, since releasing the three selected ones
has already used the whole count margin.
If \(f\) has positive limiting weight, use that
root instead as modulus after fixing the small term.
This gives \(d\le s+f<3/7\). Otherwise the larger
root is an effective singleton with residual
numerical height \(X^{o(1)}\). Lock the singleton
of the other tall letter and reverse the modular
count with the small term fixed. Since
\(l_F\ge\beta+f\), one has \(l_F\ge3/7\),
\(l_J\le3/7\) and \(\gamma\le1/7\), so \(d\le2/7\).

It remains that both other letters are multiple-heavy.
Equality in \(d\le s+2u\) forces the two heavy masses
on each tall letter \(T\) both to equal
\(\beta_T=(l_T-u)/2\). Their total high mass is
\(H=4/7\) and no fifth heavy root on the tall pair
is possible. Also, no coordinate on the small
letter can have both positive limiting normalized
mass and degree tending to infinity: releasing
that coordinate with the four tall heavy roots
would give a strict improvement by (30).

If the small height tends to zero, the two tall
letters cannot have an effective singleton. If
one does, lock a positive-weight root on the other
tall term, leaving the whole small term and the
other residual variable. The scalable count is at
most the chosen modulus mass, which is at most
\((1-d-u)/2=3/14\) here. Thus all four tall weights
have positive limits. The same rational ratio count
as in the all-mass low-short case yields
\(d\le H-\max\beta\). Equality requires
\(\beta_B=\beta_C=1/7\). Geometry then forces all
four weights to be \(1/2+o(1)\), and compressing
one tall pair while locking a root of the other
again gives \(d\le1/7\).

Suppose finally that the small letter has positive
limiting height. Its part in degrees at least
\(\lfloor c/\lambda\rfloor\) has positively bounded
height for small fixed \(c>0\), while its maximum
individual normalized mass tends to zero by the
release just noted. For any fixed \(\zeta>0\),
partition these coordinates into a bounded number
of groups of normalized root mass at most
\(\zeta\), by consecutive greedy packing. At least
one group has positive weight. Its fixed product
gives a \(k\)-full modulus for \(k\asymp1/\lambda\)
at cost at most \(\zeta\).
Fix both tall complements and one heavy root on
the tall letter of smaller \(\beta_T\). The
three-variable congruence for the other three
coprime roots has normalized cost zero.
The resulting count costs at most
\(\zeta+2u+\min(\beta_B,\beta_C)\).
Taking \(\zeta\) arbitrarily small forces
\(\beta_B=\beta_C=1/7\).
The small letter is recovered from the two
tall terms, so no other factors of the small
letter need be fixed in this count.

If the two tall lighter weights sum to less than
one, release those two roots and both tall low
complements. There is room to release a subset of
coordinates of the small letter with positive
limiting mass in the remaining weight margin.
To see this, if it has a single coordinate of
positive limiting mass then, after passing to
a subsequence, its degree must stay bounded
and its weight tends to zero. Otherwise its
maximum normalized coordinate mass tends to
zero. Partition its full mass into \(m\) groups
of asymptotically equal mass for some sufficiently
large fixed \(m\). Since the total weight is
at most one, one group has limiting weight
at most \(1/m\), and all groups have positive
mass. Geometry now improves \(d\le3/7\) strictly:
before that extra release the fixed mass is
just \(2/7+1/7=3/7\).

Consequently all four tall weights are \(1/2+o(1)\).
Their equal limiting masses imply that the four
degrees coalesce about a single high power.
Use the arbitrarily cheap positive-weight modulus
on the small letter constructed above; compress
the two tall products at this common exponent.
The residual ratio has height \(X^{o(1)}\),
so the scaled count has normalized cost at most
\(\zeta\). Since \(\zeta\) can be arbitrarily
small, this contradicts \(d=3/7\).

Every equality case is excluded. If no absolute
gap existed, successively smaller gaps and budgets,
with decomposition losses \(o(\lambda)\), would
give a sequence with \(D/\lambda\to3/7\).
Its normalized limit would be the equality case
just ruled out.

## An explicit gap below \(3/7\)

The preceding argument can be quantified. For every sufficiently
small fixed \(\lambda>0\),
\[
 N_\lambda(X)\ll_{\lambda,\varepsilon}
 X^{(3/7-10^{-15})\lambda+\varepsilon}.
\]

Work by contradiction with a limiting sequence as in the \(3/7\)
proof. Put \(\delta=10^{-15}\), \(\rho_0=10^{-6}\), and suppose
\(3/7-\delta\le d\le3/7\). As before, \(u=1-2d\),
\(l_T\le1-d\), and the lower threshold for each heavy mass on a
letter with \(l_T>u\) is \(\theta_T=(l_T-u)/2\).
All the quantities \(d,l_T,\theta_T,S_T,\beta_T\) in this
argument are limits normalized by \(L\). In contrast,
\(h,w\) denote limiting numerical heights to base \(X\).
Thus a root with normalized mass \(\beta\) and positive
limiting weight \(w\) has degree asymptotic to \(w/(\beta L)\).

In the finite-order scaled counts below, write \(t=\limsup
\log_X H(T)/L\) and \(e=\limsup \log_X H(V)\) for the
compressed root and residual ratio. For a modulus of
positive limiting weight \(w\), the scalable independent-exponent
lemma applies with \(kL\to b>0\) and amplification order
\(\ell\) whenever
\[
 t<2w/b,\qquad e<w/(16384\ell^3),
\]
provided the modulus is \(k\)-full and its free power exceeds
\(\lfloor k/(256\ell^2)\rfloor\).
The normalized cost is then \(w/(\ell b)\).
All parameter comparisons below have strict margins.

### All letter masses above \(u+20\delta\)

If the short height has a positive limit, the preceding \(3/7\)
case estimates apply. The case of three singletons and the case of one
multiple-heavy letter each give \(d\le2/5\). With two multiple-heavy letters and a
singleton \(M\), we obtained \(8d\le(7-l_M)/2\). It follows
that \(l_M\le 7-16d\le u+14\delta\), contrary to this mass
range.

For three multiple-heavy letters, release two roots on each.
Their six lower thresholds leave slack at most \(3u-d\le7\delta\),
because their total mass is at most \(1-d\). Each threshold
exceeds \(10\delta\), so there is no seventh heavy root, and
the two heavy masses on \(T\) each lie in
\([\theta_T,\theta_T+7\delta]\).
If a lighter heavy root of a letter has limiting weight
strictly below \(1/2\), pair it with either other lighter
root and release them and all low complements in (9).
Writing \(H\le1-d\) for total heavy mass, we get
\(d\le H-\theta_T-\theta_U\) for each choice of the other
letter \(U\). Hence each remaining letter has mass at least
\(1-4u\), and \(l_T\le8u-1\le u+14\delta\). This is again
impossible. All six weights are therefore exactly \(1/2\)
in the limit.

Let \(F,M,J\) be ordered by increasing \(l_T\).
Fix the two heavy roots on \(F\) and lock a heavy
root on \(M\). On \(J\), multiply its two roots,
using the degree of one of them as the compressed
free exponent. Since \(\theta_J\ge(1/3-u)/2>9/100\),
the residual height exponent is at most \(7\delta/(2\theta_J)<39\delta\).
Take \(w=1/2,\ b=w/(\theta_J+20\delta)\), and \(\ell=1\).
Fullness follows from \(\beta_M\le\theta_J+7\delta\);
the free exponent exceeds the threshold since its
normalized degree is at least \(w/(\theta_J+7\delta)\).
The compressed root height is at most
\(2\theta_J+14\delta<2w/b\).
The resulting cost is at most
\[
 2\theta_F+\theta_M+\theta_J+41\delta
 \le \frac{8}{21}+41\delta<d .
\]

If the short height tends to zero, the effectively
single cases are cut by the preceding \(3/7\) proof.
For the two remaining multiple-heavy tall letters,
put \(a=(5d-2)/2\) and let \(H\le1-d\) be their
total heavy mass. The real counts give all heavy
masses at least \(a\). If either letter has at least
four heavy roots, use that one as modulus letter
and fix the full heavy part of the other. The direct
ratio reconstruction omits at least three masses
\(a\), giving \(d\le H-3a\), hence \(d\le8/19\).

If a tall letter has exactly three heavy roots,
any two of its positive-weight roots can be chosen
in turn as modulus. The direct ratio bound requires
the sum of the two nonmodulus masses to be at most
\(H-d\le u\) for both choices. Thus all three masses
are at most \(u-a\). If the other tall letter also has
three, we obtain \(H\le6(u-a)\), whereas
\(H\ge d+2a\), giving \(33d\le14\).
If the other has two, use the partner of the larger
root there as modulus. The resulting bound
\(d\le (H+3(u-a))/2\) again gives \(33d\le14\).
Both contradict the assumed range.

With two heavy roots on each tall letter, the
direct ratio count gives \(d\le H-\max\beta\).
Thus \(H\ge4d/3\), every \(\beta\le u\), and every
\(\beta\ge4d/3-3u\). In particular
\[
 1/7-22\delta/3\le\beta\le 1/7+2\delta .
\]
Release the two lighter tall roots along with the
tall low complements and the entire short term
in (9). If their limiting weights sum to less
than one, this gives \(d\le 2u\), impossible.
The four weights must be \(1/2\).
Lock one tall root and compress the other tall
pair at a common exponent as before. With
\(w=1/2,\ b=3w,\ \ell=2\), the root height is at most
\(2u<2/3\), the residual height exponent is at most
\((28\delta/3)/(2(1/7-22\delta/3))<34\delta\),
and fullness and the free exponent condition
follow from the displayed mass range. The cost
is at most \(u+1/6<d\).

### A mass at most \(u+20\delta\)

Write \(s\le u+20\delta\) for this mass and \(h\) for the
letter's limiting weight. The other masses are at least
\(d-s>u\). If the selected letter is tall and the short
gap stays above \(\lambda^{1/4}\), the previous real
count gives \(d\le s+u\le2u+20\delta\), impossible.
Otherwise the other two letters have limiting height
one.

If the other letters are both single-heavy, fixing the
small term and locking either tall root yields
\(\beta_j\ge d-s\). Releasing both tall roots gives
\(\beta_1+\beta_2\le1-d\), hence
\[
 2/7-23\delta\le\beta_j\le2/7+24\delta.
\]
If \(h\le\rho_0\), leave the small term free and lock
one tall root, taking \(w=1,\ b=3,\ \ell=3\).
Then \(t\le2/7+24\delta<2/3\) and
\(e\le \rho_0<1/(16384\cdot27)\).
The modulus and free degrees both have normalized
degree near \(7/2\), so the degree conditions hold.
The cost is at most \(2/7+24\delta+1/9<d\).
If \(h>\rho_0\), fix the small term and use its part
in degrees at least \(k=\lfloor2h/L\rfloor\).
This has weight \(w\ge h(1-2s)>7h/10\).
Compress both tall roots into a common power ratio
of degree \(\lfloor7/(2L)\rfloor\).
Here \(t\le 2/7+24\delta<7/10\) and
\(e\le168\delta<(7\rho_0/10)/131072\).
Thus we can take \(b=2h,\ \ell=2\); the normalized
free degree is \(7/2\), while \(b\le2\).
The cost \(s+w/(4h)\le s+1/4\) is below \(d\).

If one of the other letters is multiple-heavy and
the other singleton, call the singleton mass \(\gamma\).
If the multiple-heavy letter has two positive-weight
heavy roots, then with the small term fixed each of
these can be the modulus. Thus both masses are at
least \(d-s\), while \(\gamma\ge(d-s-u)/2\).
The three-root release gives \(7d\le2+5s+u\), hence
\(19d\le8+100\delta\), impossible.
If it has just one positive-weight root \(\beta\),
it is effectively single. With the small term fixed
we can alternately lock either main root, so both
\(\beta,\gamma\ge d-s\). An extra heavy root \(f\)
on the multiple-heavy letter \(F\) satisfies
\(f\ge(l_F-u)/2\) and \(l_F\ge\beta+f\), hence
\(f\ge d-s-u\). The three-root release gives
\(4d\le1+3s+u\), hence \(12d\le5+60\delta\).

It remains that both other letters are multiple-heavy.
Releasing two from each after fixing the small term gives
\(d\le s+2u\), so \(s\ge1/7-5\delta\).
The two thresholds on each tall letter sum, over the two
letters, to \(1-s-2u\). This leaves at most \(s+2u-d\le
27\delta\) of slack in the upper bound \(1-d\) on
released mass. Every tall heavy threshold exceeds
\(7/100\), so each tall term has exactly two heavies.
Their masses are each within \(27\delta\) of the relevant
threshold. Any coordinate of the small letter whose
degree tends to infinity has normalized mass at most
\(27\delta+o(1)\), since otherwise releasing it together
with the four tall roots in (30) contradicts this slack.

If \(h\le\rho_0\) and a tall letter is effectively single,
lock a largest-weight root on the other tall term and
use its modulus with the effective singleton and the
small letter free. Then \(w\ge1/2\),
\(\beta_{\mathrm{lock}}\le d/2+27\delta\),
\(t\le d/2+27\delta<4/7\), and \(e\le\rho_0\).
Take \(b=7w/2,\ \ell=2\) (the free degree is greater
than \(4/L\)). The cost is
\(\beta_{\mathrm{lock}}+1/7<d\).
Otherwise use either tall term's largest-weight root
as modulus, and freeze both heavy roots of the other.
The ratio recovery applies since the height exponent
of the residual is at most \(\rho_0<w/2\).
This gives \(d\le1-d-\max(\theta_B,\theta_C)\).
Thus \(\max(\theta_B,\theta_C)\le u\), and
\(\theta_B+\theta_C=(1-s-2u)/2\ge2/7-13\delta\)
gives \(\min(\theta_B,\theta_C)\ge1/7-15\delta\).

If \(h>\rho_0\), the part of the small term in degrees
at least \(\lfloor h/(2L)\rfloor\) has positive weight.
Its high coordinates have root mass at most
\(27\delta+o(1)\). Group them into a bounded number of
blocks of root mass at most \(28\delta+o(1)\).
One block has positive weight and can be fixed as
a powerful modulus of cost at most \(28\delta\).
Freeze both tall low cofactors and a root of the
tall letter with smaller threshold, costing at most
\(2u+\min(\theta_B,\theta_C)+27\delta\).
The ternary congruence charges normalized exponent zero
for the remaining three roots, so
\(\min(\theta_B,\theta_C)\ge d-2u-55\delta
\ge1/7-60\delta\).
Also \(\theta_B+\theta_C\le(1-d)/2\), using \(s\ge d-2u\).
It follows in either height range that all four tall
heavy masses differ from \(1/7\) by at most \(90\delta\).

Let \(w_B,w_C\) denote the limiting lighter tall weights.
Release these two roots and both tall low cofactors.
The fixed mass before releasing anything from the small
letter is at most \(s+1-d-\theta_B-\theta_C
\le 3/7+36\delta\).
A subset of the small letter has limiting mass at least
\(38\delta\) and weight at most \(610\delta h\).
If a bounded-degree coordinate has such mass, take it,
at weight zero in the limit. Otherwise the maximum
normalized coordinate mass is at most \(38\delta+o(1)\)
(using the \(27\delta+o(1)\) bound for unbounded degrees).
In increasing degree order stop when mass first reaches
\(38\delta\). This costs at most \(76\delta+o(1)\) in
normalized mass; its weight divided by its mass is at
most the weighted average for the whole letter.
Because \(s\ge1/7-5\delta>1/8\), the claimed weight
bound follows. Geometry (9) now forces
\(w_B+w_C\ge1-610\delta h\), since otherwise this
extra release makes the fixed mass at most
\(3/7-2\delta<d\).
All tall heavy weights are consequently within
\(610\delta h\) of \(1/2\).

Compress each tall pair using \(q=\lfloor7/(2L)\rfloor\).
Each tall heavy contributes at most
\(610\delta h+315\delta\le925\delta\) to the limiting
residual exponent. If \(h\le\rho_0\), lock a tall heavy
of weight \(w\ge1/2\) and count the other tall pair
with the small term left free. With \(b=3w,\ \ell=2\),
we have \(t\le2/7+180\delta<2/3\) and
\(e\le\rho_0+1850\delta<w/131072\);
the modulus and free degree conditions hold.
The cost is at most \(1/7+90\delta+1/6<d\).
If \(h>\rho_0\), fix the small term and use its part
in degrees at least \(k=\lfloor2h/L\rfloor\),
of weight \(w>7h/10\). With \(b=2h,\ \ell=2\), the
tall ratio has \(t\le2/7+180\delta<7/10\) and
\(e\le3700\delta<(7\rho_0/10)/131072\).
The cost is at most \(s+1/4<d\).

All limiting cases are contradicted. Taking the
anatomic loss \(o(\lambda)\) along a hypothetical
sequence of failing budgets, as in the \(3/7\)
argument, proves the stated fixed coefficient.

## A prime-root graph and a larger explicit small-budget gap

For every sufficiently small fixed \(\lambda>0\),
\[
N_\lambda(X)\ll_{\lambda,\varepsilon}
X^{(3/7-10^{-4})\lambda+\varepsilon}.
\]

### A graph count with two or three prime roots

In limiting boxes let \(n=X^{w+o(1)}\) be a \(k\)-full modulus with \(w>0\) and \(k\asymp 1/L\). Consider
\[
U\prod_{i\in I}x_i^{p_i}\equiv
K V\prod_{i\in J}x_i^{p_i}\pmod n ,
\tag{PG}
\]
where \(K\) is a fixed unit and all the variables are units. The disjoint sets \(I,J\) have total size \(m\in\{2,3\}\). Either set may be empty. The \(x_i\) are pairwise distinct primes of bounded normalized mass and positive limiting weight, so \(p_i\asymp 1/L\). The cofactors \(U,V\) have numerical heights \(X^{o(1)}\) and range over at most \(X^{(f+o(1))L+\varepsilon}\) and \(X^{(g+o(1))L+\varepsilon}\) values. The normalized count exponent is at most
\[
G(f,g)=\max\bigl(f,g,\min(f+g/2,g+f/2)\bigr).
\tag{PG1}
\]

Color primes with \(m\) colors. Some coloring retains at least \(m^{-m}\) of the tuples with the matching prime color in each coordinate. Graph the retained \(U,V\), choosing one prime-root witness per edge. At fixed \(U,V\) the amplified integer ternary congruence bounds witness multiplicity at cost \(o(L)\) in the \(X\)-exponent. For \(m=2\), pad the third coordinate by \(1\). Throughout, the root heights are at most \(n^{C/k}\) for a fixed \(C\), and the amplification order can be any fixed integer before taking the limit.

At a fixed middle vertex \(V\), call a two-edge path good if at most one prime coordinate is unchanged between its chosen witnesses. For \(m=3\), a bad path shares at least two prime coordinates. Given \(V\) and one edge, fix any two of its prime coordinates. The other edge is counted by the scalable one-high congruence with its remaining prime and \(U\) free, at cost \(o(L)\). Indeed the root exponent and the fullness threshold have order \(1/L\), the root height is \(X^{O(L)}\), and the cofactor height is \(X^{o(1)}\). One may decrease the fullness threshold by a fixed factor to fit the root-height cap and take arbitrarily large fixed amplification. For \(m=2\) a bad path shares both primes, which determines \(U\) modulo \(n\) and hence exactly. The number of bad paths is at most \(|E|X^{o(L)+\varepsilon}\) for \(E\) the edge set.

For fixed distinct endpoints \(U,U'\), eliminate \(V\) between the two congruences. This gives a ternary congruence for the rational quotients of the three corresponding prime coordinates, or for two quotients with the third coordinate \(1\). Each quotient is \(1\) or a reduced ratio of two distinct primes. To amplify, take the set of distinct quotient tuples and independently color primes by numerator or denominator. Some such coloring retains at least \(2^{-2m}\) of the tuples. In an \(\ell\)-fold coordinate product of retained tuples, the numerator and denominator have no cancellation, and the divisor bound gives at most \(n^\varepsilon\) input tuples per product, with an arbitrarily small \(\varepsilon\) for fixed \(\ell\). The original coordinate colors assign each prime to a unique coordinate in the products. Its valuation in the exact monomial therefore determines the signed valuation in that coordinate. The exact fiber bound in (R1) is \(1\). Apply (R1) with height constant \(\ell C\) and the target raised to the \(\ell\)-th power. On taking \(\ell\)-th roots, the quotient count is \(O(n^{2/(\ell k)+\varepsilon})\), which costs \(o(L)\) in the limiting \(X\)-exponent.

Each nontrivial quotient fixes both of its primes. In a good path, the remaining data consist of \(V\) and at most one common prime. A single-high scalable count costs \(o(L)\) if there is a common prime, and otherwise \(V\) is determined by its residue. Thus each distinct endpoint pair contributes \(X^{o(L)+\varepsilon}\) good paths. Summing squares of the degrees on the right gives at most \(X^{o(L)+\varepsilon}(|E|+X^{2fL})\). Interchanging the sides gives the corresponding bound with \(g\). Cauchy–Schwarz proves (PG1).

In a primitive anatomic box the fixed components of the two free terms can be absorbed into \(K\). Each pair of complete terms determines the third. The divisor bound controls the multiplicity of the original coordinates giving the low cofactor values.

We also use a consequence of (A). Write \(d=D/L\), \(u=1-2d\) and \(\theta_T=(l_T-u)/2\), where \(l_T\) is the normalized mass of the letter. Let a positive-weight heavy root on \(T\) have normalized mass \(\beta\), and suppose \(\beta<j\theta_T\) for an integer \(j\ge2\). After treating the families bounded directly by (A), the root has at most \(j-1\) prime factors counted with multiplicity, apart from a cofactor that can be fixed at cost \(o(L)\). Indeed, choose a normalized threshold strictly between \(\beta/j\) and \(\theta_T\). A divisor of a degree-\(e\) coordinate with unnormalized mass between \(2\lambda/e\) and this threshold times \(L\) gives (A) with normalized exponent strictly below \(d\). The product of the prime factors below \(2\lambda/e\), counted with multiplicity, stays below that same cutoff, or the first partial product reaching it would give such a divisor. Every remaining prime factor has normalized mass above the chosen threshold. When a later estimate uses the lower bound \(\theta_T\) for those factors, take the threshold to be \(\theta_T-c\) with arbitrarily small fixed \(c>0\). The discarded families cost at most \(d-c+o(1)\) and the surviving estimates below have strict margins. In particular, a root with \(\beta<2\theta_T\) can be treated as prime.

### Starting the limiting argument

Take the limiting notation of the \(3/7\) proof, set \(\delta=10^{-4}\), and suppose
\[
3/7-\delta\le d\le 3/7.
\]
Pair counting gives \(l_T\le 1-d\). A letter with \(l_T>u\) has boundedly many heavy coordinates, each of mass at least \(\theta_T\). Its complementary numerical height is \(X^{o(1)}\). As in the \(3/7\) proof, a fixed positive-weight root gives a positive-weight modulus with fullness threshold of order \(1/L\).

### Numerical height zero on the short letter

If a tall letter has mass at most \(u+20\delta\), fix it and use the gap estimate (19) on the other tall letter. This gives \(d\le2u+20\delta\), a contradiction. Thus the two tall letters have masses above \(u\). If either is effectively single, the modular count from the zero-short-height case of the \(3/7\) proof gives \(d\le1/3\). This estimate uses no lower mass bound on the short letter.

Both tall letters therefore have at least two positive-weight heavy roots. As in the earlier zero-short-height argument, (19) gives \(l_B,l_C\ge d-u\), so all tall heavy masses are at least
\[
a=(5d-2)/2.
\]
Their total \(H\) is at most \(1-d\) by (30). Fix a positive-weight root on one tall letter as modulus and freeze all heavy roots on the opposite tall letter. Its low cofactor and the entire short term have numerical heights \(X^{o(1)}\). Direct ratio recovery omits all other roots on the modulus letter.

With at least four roots on a tall letter, this gives \(d\le H-3a\) and \(d\le8/19\). With three roots, choose either of two positive-weight roots as modulus. Since each omitted pair has total mass at most \(H-d\le u\), each of the three roots has mass at most \(u-a\). If both tall letters have three roots, \(H\le6(u-a)\) and \(H\ge d+2a\). If the other tall letter has just two roots, omit its larger root to get \(d\le(H+3(u-a))/2\). In both cases \(33d\le14\), a contradiction.

With exactly two positive-weight roots on each tall letter, direct ratio recovery can omit any of the four. Thus \(H\ge4d/3\) and their masses satisfy
\[
1/7-22\delta/3\le\beta\le1/7+2\delta.
\tag{PG2}
\]
All four limiting weights are \(1/2\). Otherwise the lighter root on each tall letter, the tall low complements and the whole short letter would have total weight below one. Releasing them by (9) leaves mass at most \(2u<d\).

Fix the tall term of smaller mass, with value \(C\). On the other tall letter write the term as \(VT^q\), with \(T\) the product of its two roots and \(q\) the smaller exponent. The residual \(V\) is a positive integer. By (PG2) and equality of the two numerical weights, the limiting height exponent of \(V\) is at most
\[
\frac{28\delta/3}{2(1/7-22\delta/3)}<0.0034 .
\]
The height of \(T\), normalized by \(L\), is at most \(2u<0.287\). The equation gives \(VT^q=C+O(X^{o(1)})\). For fixed \(T\) at most one \(V\) is possible since \(T^q\) has positive limiting weight. Factoring \(T\) recovers the two roots in divisor-many ways.

On dyadic scales \(T_0,V_0\), put \(z=T/T_0\), \(v=V/V_0\). The graph \(v=\gamma z^{-q}\), with \(\gamma=C/(V_0T_0^q)\) in a fixed positive compact interval at each fixed truncation, approximates these points to error \(X^{-1+o(1)}\). Take the \(D_0=16(r+1)\) monomials
\[
v^i z^j,\qquad 0\le i\le15,\qquad
0\le j\le r,\qquad
r=\left\lfloor\frac{1}{20\cdot0.287 L}\right\rfloor .
\]
The powers \(j-qi\) on the exact graph are distinct since \(q>r\). After dividing the evaluation determinant by the ordinary Vandermonde in \(z\), the exact comparison determinant is bounded away from zero. Multiplying the original rows by a common power of \(z\) turns the Laurent-power quotient into a Schur polynomial with nonnegative coefficients, up to fixed units. The divided-difference row transformation amplifies the error by at most \(O(T_0^{D_0-1})\), using the spacing \(1/T_0\) of distinct normalized integers. The error in the transformed entries is consequently \(X^{-1+o(1)+0.287L(D_0-1)}\), which tends to zero because \(0.287 L D_0\to0.8\) as \(L\to0\).

In any interval for \(z\) of length \(X^{-0.037L}\), restore the dyadic column scales and the ordinary Vandermonde. The resulting nonzero integer determinant has magnitude bounded by a fixed constant times
\[
X^{(D_0/2)(15\cdot0.0034+r\cdot0.287 L-0.037 L(D_0-1))}.
\]
The first two terms in parentheses sum to at most \(0.101\), while the last tends to \(0.037\cdot0.8/0.287>0.101\). This is impossible for large \(X\) at fixed sufficiently small \(L\). The interval count gives normalized exponent at most \(0.037\), hence \(d\le\min(l_B,l_C)+0.037\). This excludes \(\min(l_B,l_C)\le0.39\).

If both tall masses exceed \(0.39\), their heavy thresholds exceed \(0.123\). Hence all four roots can be treated as prime. The short mass \(s\) is below \(0.22\). Take the tall letter with smaller low mass \(f\), lock a root on the opposite tall letter, and apply (PG1) with the free tall low cofactor and the entire short term. There are two free prime roots. Since \(f\le(1-s-H)/2\), the graph contributes at most
\[
G(f,s)\le\max(s,f+s/2)
\le\max\bigl(0.22,(1-4d/3)/2\bigr)=0.22.
\]
The modulus costs at most \(u<0.144\), again a contradiction.

### Positive limiting short height and all masses above \(u+20\delta\)

The \(3/7\) case estimates give \(d\le2/5\) for three singleton letters or exactly one multiple-heavy letter. With two multiple-heavy letters and a singleton \(M\), they give \(l_M\le7-16d\le u+14\delta\).

For three multiple-heavy letters, choose two heavy roots on each. The six lower thresholds leave at most \(3u-d\le7\delta\) of slack in (30), and each threshold exceeds \(10\delta\). Thus there are exactly six heavy roots, and each mass on \(T\) lies in \([\theta_T,\theta_T+7\delta]\). All six roots can be treated as prime. If one lighter root had weight below \(1/2\), release it with the lighter root on either other letter and all low complements in (9). This gives \(d\le H-\theta_T-\theta_U\) for each choice of \(U\), with \(H\le1-d\). Each remaining letter must have mass at least \(1-4u\), so \(l_T\le8u-1\le u+14\delta\). Thus all six limiting weights are \(1/2\).

Use the two letters with smallest thresholds. Lock a root on one as modulus and freeze a root on the other. Apply (PG1) to the three free primes and their two low cofactors. Each low mass is at most \(u\), and the two smallest thresholds sum to at most \(1/3-u\). The bound is
\[
d\le 1/3-u+14\delta+3u/2
\le17/42+15\delta<d.
\]

### Positive limiting short height and a mass at most \(u+20\delta\)

Write \(s\le u+20\delta\) for this mass and \(h>0\) for the corresponding height. The two other masses are at least \(d-s>u\). If the fixed small-mass term is tall and the short gap stays above \(\lambda^{1/4}\), (19) gives \(d\le s+u\le2u+20\delta\). Otherwise both other letters have limiting height one. The mixed singleton and multiple-heavy case is ruled out by the preceding explicit-gap estimates \(19d\le8+100\delta\) and \(12d\le5+60\delta\).

We will use a cheap modulus from the small-mass letter. If releasing specified tall roots in (30) leaves slack \(m_0\) in \(1-d\), each small-letter coordinate whose degree tends to infinity has mass at most \(m_0+o(1)\), by releasing it as well. The part of the small term in degrees at least \(\lfloor h/(2L)\rfloor\) has positive weight, since its complementary weight is at most \(hs/2+o(1)\). Divide those high coordinates into blocks of mass at most \(m_0+\delta+o(1)\) by starting a new block when the next coordinate would exceed that cap. The number of blocks is bounded, since the masses of consecutive completed blocks sum to more than the cap. Some block has positive limiting weight. It can be fixed as a modulus of cost at most \(m_0+\delta\).

If both tall letters are single-heavy, their high masses satisfy \(\beta_j\ge d-s\) by fixing the small term, locking one tall root, and counting the other by the scalable congruence. Releasing both tall roots gives \(\beta_1+\beta_2\le1-d\) and \(s\ge(3d-1)/2\). The cheap modulus then costs at most \(48\delta\), since \(1-d-\beta_1-\beta_2\le1-3d+2s\le47\delta\). With this modulus, fixing either tall term in full and counting the other singleton gives \(l_j\ge d-48\delta\). It follows that
\[
l_j\le 3/7+50.5\delta,\qquad
\theta_j\ge1/7-25.5\delta,\qquad
2/7-23\delta\le\beta_j\le2/7+24\delta<3\theta_j.
\]
Each singleton can thus have at most two prime factors of positive weight, with any tiny remainder fixed. A square uses one prime with twice the exponent. If four distinct prime roots remain, fix one, costing at most \(\beta_j-\theta_j\le1/7+49.5\delta\). The other three are counted by (PG1). Each low mass is at most \(1/7+73.5\delta\), so the total cost is at most
\[
48\delta+1/7+49.5\delta+(3/2)(1/7+73.5\delta)
=5/14+207.75\delta<0.38<d .
\]

Finally suppose both tall terms are multiple-heavy. Releasing two roots from each gives \(s\ge d-2u\ge1/7-5\delta\), and the four tall thresholds exceed \(0.07\). The slack is at most \(s+2u-d\le27\delta\). Thus both tall terms have exactly two heavy roots and each heavy mass lies within \(27\delta\) of its threshold. The cheap modulus costs at most \(28\delta\).

To verify positive weight on all four roots, let \(w_B,w_C\) be the lighter tall weights. Releasing those two roots and both tall low complements leaves fixed mass at most
\[
s+1-d-\theta_B-\theta_C\le3/7+36\delta .
\]
There is a subset of the small letter with mass at least \(38\delta\) and weight at most \(610\delta h\). If a bounded-degree coordinate has that mass, it has zero limiting weight. Otherwise the largest coordinate mass is at most \(38\delta+o(1)\), using the \(27\delta+o(1)\) bound on unbounded degrees. Take the first prefix in increasing degree order crossing \(38\delta\). Its mass is at most \(76\delta+o(1)\), and its weight per unit mass is bounded by the average for the letter. Now use \(s>1/8\). If \(w_B+w_C<1-610\delta h\), releasing this subset as well in (9) leaves fixed mass at most \(3/7-2\delta<d\). Each lighter weight is therefore at least \(1/2-610\delta>0\).

All four tall heavies can be treated as prime, since \(27\delta<\theta_T\). Freeze one on the tall letter with smaller threshold. With the cheap modulus and the graph count for the remaining three primes, the cost is at most
\[
28\delta+(1-d)/4+27\delta+3u/2<d.
\]
Here \(\theta_B+\theta_C\le(1-d)/2\) and each tall low mass is at most \(u\).

These contradictions exhaust the limiting boxes. The same fixed-truncation and dyadic passage as in the \(3/7\) argument gives the claimed bound for all sufficiently small fixed \(\lambda\).

### Improving the gap to \(10^{-3}\)

Using (PG1) with the smallest available modulus root sharpens the zero-short-height case. The resulting bound is
\[
N_\lambda(X)\ll_{\lambda,\varepsilon}
X^{(3/7-10^{-3})\lambda+\varepsilon}
\qquad\text{for all sufficiently small fixed }\lambda>0.
\]
Put \(\delta=10^{-3}\) in the limiting argument. We give the estimates that change at this value of \(\delta\).

Suppose first that the short height is zero. The earlier case reductions still give two positive-weight heavy roots on each tall letter. All four weights are \(1/2\) and the masses satisfy (PG2). In the real estimate for \(VT^q=C+O(X^{o(1)})\), suppose \(T\) has normalized height below \(t_0\) and \(V\) has limiting numerical height exponent below \(e_0\). Use the monomials \(v^i z^j\) with
\[
0\le i\le m,\qquad 0\le j\le r,\qquad
r=\left\lfloor\frac{c}{(m+1)t_0L}\right\rfloor,\qquad
D_0=(m+1)(r+1),
\]
where \(0<c<1\) and \(q>r\). The divided differences amplify the entry error by at most \(X^{t_0L(D_0-1)}\), with \(t_0LD_0\to c<1\). On intervals of normalized length \(X^{-\alpha L}\) the integer determinant vanishes for large \(X\) if
\[
\alpha c/t_0>me_0+\frac{c}{m+1}.
\tag{PG3}
\]
This gives normalized exponent at most \(\alpha\). In all the following choices \(t_0=0.29\) and \(qL\) is eventually greater than \(3\), so \(q>r\).

First (PG2) gives \(e_0=0.0345\). Take \(m=5\), \(c=0.99\) and \(\alpha=0.100\) in (PG3). Fixing the tall term of smaller mass excludes \(\min(l_B,l_C)\le0.327\), since \(0.327+0.100<3/7-\delta\). With both masses above \(0.327\), the heavy thresholds exceed \(0.091\) and all four roots can be treated as prime.

If the short mass \(s\le0.282\), apply (PG1) to the smaller tall low cofactor and the whole short term as before. The bound is at most \(u+0.282<3/7-\delta\). Suppose \(s>0.282\). Write \(f_i\) for the low mass of the tall letter \(i\in\{B,C\}\). The bound
\[
f_i\le 1-s-\min(l_B,l_C)-2(1/7-22\delta/3)<0.121<s/2
\]
gives \(G(f_i,s)=s\) with either tall letter free. Put \(\nu_j=u-\beta_j\ge0\) for the four roots and \(M=\max_j\nu_j\). We may use a root of mass \(u-M\) as modulus. If \(a_B,a_C\) are the sums of the deficits \(\nu_j\) on the two tall letters, their heavy masses are \(S_i=2u-a_i\). Applying (19) with either tall term fixed and the two roots of the other term released gives
\[
d\le1-s-\max(S_B,S_C).
\]
Combining this with the modulus estimate gives
\[
M\le s-(d-u)\le\min(a_B,a_C)\le14\delta/3.
\]
The last inequality follows from \(a_B+a_C=4u-H\le28\delta/3\). The mass difference within either tall pair is at most \(M\). Compressing at the smaller exponent therefore leaves residual height at most \(M/(2(u-M))\).

For \(s\le0.287\), use \(M\le0.287-(2/7-0.003)\) and \(u\ge1/7\). This gives \(e_0=0.0155\). Take \(m=7,c=0.995,\alpha=0.068\) in (PG3). Since \(s\ge d-u\ge2/7-0.003\), fixing the smaller tall mass and applying this count costs at most
\[
\frac{1-(2/7-0.003)}{2}+0.068<3/7-\delta .
\]
For \(s>0.287\), use \(M\le14\delta/3\) to take \(e_0=0.0169\). The same \(m,c\) permit \(\alpha=0.071\). The smaller tall mass is below \(0.3565\), and the total cost is at most \(0.4275<3/7-\delta\). This finishes the zero-short-height case.

With positive short height and all masses above \(u+20\delta\), the cost \(17/42+15\delta\) from (PG1) is still below \(3/7-\delta\). For a mass at most \(u+20\delta\), the gap bound \(2u+20\delta\) and mixed-letter bounds \(19d\le8+100\delta\) and \(12d\le5+60\delta\) remain sufficient.

If both tall letters are single-heavy, the cheap modulus still costs at most \(48\delta\), and the bounds on \(\beta_j,\theta_j\) still allow splitting into at most two prime factors. If there are at most three prime roots after splitting, fix no further factor. The low masses satisfy
\[
\max(f_B,f_C)\le1/7+73.5\delta,\qquad
f_B+f_C\le1-s-2(d-s)=u+s\le2/7+24\delta.
\]
Consequently (PG1) with the modulus costs at most
\[
48\delta+\max\bigl(1/7+73.5\delta,\tfrac34(2/7+24\delta)\bigr)
<0.281 .
\]
If each tall heavy root splits into two distinct prime factors, then \(\beta_j\ge2\theta_j\) and hence \(f_j\le u\). As in the prime-splitting argument, these inequalities permit arbitrarily small slack from the finite-box divisor thresholds. Freeze the smaller prime in either pair, at cost at most \(\beta_j/2\le1/7+12\delta\). The modulus, the frozen prime and (PG1) cost at most
\[
48\delta+1/7+12\delta+3u/2<0.421<3/7-\delta .
\]

If both tall letters are multiple-heavy, their thresholds exceed \(0.058\), so the \(27\delta\) slack still gives exactly two heavy roots on each. The small-letter subset for (9) admits a tighter bound. Traverse its coordinates in order of increasing degree until first reaching mass \(38\delta\). If this occurs below degree \(\lceil\lambda^{-1/2}\rceil\), the selected weight tends to zero. Otherwise the last coordinate has mass at most \(27\delta+o(1)\) by (30), so the total selected mass is at most \(65\delta+o(1)\). The weight per unit mass is bounded by the whole-letter average, and \(s\ge1/7-5\delta\), so the selected limiting weight is at most \(480\delta h\). Applying (9) gives \(w_B+w_C\ge1-480\delta h\), since a strict failure would leave fixed mass at most \(3/7-2\delta\). Each lighter tall root therefore has weight at least \(1/2-480\delta>0\).

The tall roots can all be treated as prime, since \(27\delta\) is less than either heavy threshold. The same cheap modulus and one additional frozen prime yield the bound
\[
28\delta+(1-d)/4+27\delta+3u/2<0.416<3/7-\delta .
\]
The limiting contradiction and the dyadic passage prove the larger gap.

## A fixed-target surface with a cubic factor and a quartic power

For fixed nonzero integers \(B,C\), consider the equation
\[
 Ux^3+B y^4=C.
\]
The variables are positive and have scales \(X^\alpha,X^\beta,X^\gamma\),
respectively. Assume \(|B|,|C|\ll X^{O(1)}\), with the polynomial growth
bound fixed. Set
\[
 H=\max(\alpha,\beta,\gamma),\qquad
 \kappa=\sqrt{\frac{\alpha\beta\gamma}{\alpha+3\beta}}
\]
when \(\alpha+3\beta>0\). The number of solutions is
\[
 \ll_\varepsilon X^{\max\{2\kappa,\kappa+H/3\}+\varepsilon}.
 \tag{SC}
\]
The constant is uniform in \(B,C\) at these fixed scales. It suffices
to establish the estimate for positive side exponents, since a
zero exponent is handled by enlargement through \(X^\delta\) and
taking \(\delta\) small. If \(\alpha=\beta=0\), choosing \(U,x\)
already gives exponent zero.

The surface is absolutely irreducible. As a polynomial in \(U\)
it is linear, and its two coefficients \(x^3\) and \(By^4-C\)
are coprime in \(\overline{\mathbb Q}[x,y]\). Apply
Browning–Verzobio, Theorem 1.1, with side modulus one and a
nonconstant side polynomial, such as \(g=x+1\).
Take a lexicographic order in which \(U\) is the first coordinate.
The highest monomial is \(Ux^3\). In formula (5) the modulus term
vanishes and the parameter is \(K\ll X^{\kappa+\varepsilon}\),
as \(\log T_{\underline m}=(\alpha+3\beta)\log X+O(1)\).
Consequently there are \(O_\varepsilon(X^{2\kappa+\varepsilon})\)
exceptional points and \(O_\varepsilon(X^{\kappa+\varepsilon})\)
proper auxiliary intersections. Their degrees, and hence the
numbers of curve components per intersection, are bounded for fixed
positive side exponents by the degree estimate in Theorem 1.1 and
Bézout. The constant in the covering theorem depends only on the
polynomial degrees and \(\varepsilon\).
[Browning–Verzobio, Theorem 1.1 and (5)](https://arxiv.org/html/2408.11453v2#S1)

A geometrically irreducible component of degree \(d\ge3\) contributes
\(O_{d,\varepsilon}(X^{H/d+\varepsilon})\) integer points by the
coefficient-uniform affine curve count
[Binyamini–Cluckers–Kato, Corollary 8](https://arxiv.org/html/2401.03982v2#S3).
If such a component is not defined over \(\mathbb Q\), its rational
points lie also on a distinct Galois conjugate and contribute only
boundedly many by Bézout. The same observation covers nonrational
degree-one and degree-two components.

It remains to control degree-one and degree-two curves. A
geometrically irreducible projective line or conic normalizes over
\(\overline{\mathbb Q}\) to \(\mathbb P^1\). Its projective coordinate
forms \(L_0,L_U,L_x,L_y\) then have degree \(d=1\) or \(2\), and satisfy
\[
 L_U L_x^3 = C L_0^4-B L_y^4.
\]
The component meets the positive affine locus, so \(L_0,L_x\)
are nonzero. Suppose that \(y=L_y/L_0\) is nonconstant.
Every zero of \(L_x\) must be a common zero of \(L_0,L_y\).
Indeed a zero of just one of \(L_0,L_y\) would give a nonzero
right side. If neither vanishes, the map \(L_y/L_0\) meets a
nonzero fourth root of \(C/B\) with multiplicity at least three,
whereas the map has degree at most two and this root is simple
in \(T^4-C/B\).

The common zeros have total degree at most \(d-1\), since \(y\)
is nonconstant. This is impossible for a line. For a conic it
forces one common linear factor and a double root of \(L_x\) at its
zero. Cancel this linear factor from \(L_0,L_y\). The remaining
forms define a degree-one map, so \(C L_0^4-B L_y^4\) vanishes
at that point to order at most five. The left side has order
at least six. Thus every low-degree component has constant \(y\).
On a constant-\(y\) component, \(Ux^3=C-B y^4\) is fixed and
nonzero on positive solutions. They have divisor-many choices.
Summing the curve and exceptional-point contributions proves (SC).

In a unit-budget box, release the degree-three coordinate of \(A\)
with mass \(\beta\). Leave all of its complementary factors
variable and absorb them into the cofactor \(U\). Its size is at
most \(X^\alpha\) with \(\alpha=1-3\beta\). Freeze \(C\) entirely
and all of \(B\) except its degree-four coordinate, of mass
\(\gamma\). The signed additive equation has the shape above,
possibly after changing the signs of \(B,C\). The frozen mass is
\(f=L_B+L_C-\gamma\) and the count for the surface gives the box
bound
\[
 D\le f+\max\left\{2W,W+\frac{\alpha}{3},
                        W+\frac{\beta}{3},W+\frac{\gamma}{3}\right\}
 +\varepsilon,\qquad W=\sqrt{\alpha\beta\gamma}.
\]
The loss from dyadic endpoints and fixed exceptional coefficients
can be absorbed in \(X^\varepsilon\). Positive \(\beta,\gamma\)
ensure that the surface freeze count counts all tuples even when
\(U\) is a product, since an exact \(U\) has divisor-many
factorizations for a fixed number of coordinates.

For fixed rational \(t,r>0\), the exact inequality
\(
 \alpha\beta\le 3t^2+(1-6t)\beta
\)
follows from \(3(\beta-t)^2\ge0\).
The arithmetic–geometric mean inequality then gives the affine
upper bound
\[
 W\le\frac{3t^2+(1-6t)\beta}{2r}+\frac{r\gamma}2.
\]
Thus any fixed \(t,r>0\) gives a maximum of four affine bounds
usable in the certificate.

At the unit-budget witness with scale denominator \(6274\),
the selected masses are \(\alpha=1189/6274\),
\(\beta=1695/6274\), \(\gamma=1342/6274\), and \(f=2048/6274\).
Taking \(t=27/100\), \(r=1/2\) gives
\(W\le 3283619/31370000<21/200\).
Each of \(\alpha/3,\beta/3,\gamma/3\) is less than \(21/200\).
In particular the bound is at most
\[
 \frac{2048}{6274}+\frac{21}{100}
 =\frac{168277}{313700}<\frac{27}{50}.
\]
These are strict material inequalities and persist on a
neighborhood of the witness.

### A diagonal surface with a cubic power and two quartic powers

Freeze every coordinate except a degree-three coordinate of
mass \(\beta\) and the degree-four coordinates of the other
two terms, of masses \(\gamma,\xi\).
With nonzero fixed coefficients the signed equation has the form
\[
 Ux^3+Vy^4+Wz^4=0
\]
on pairwise coprime positive roots.
The surface is absolutely irreducible by Eisenstein in \(x\)
at any one of the distinct linear factors of \(Vy^4+Wz^4\)
over \(\overline{\mathbb Q}\).

Use Browning–Verzobio, Theorem 1.1, with side modulus one and
leading monomial \(z^4\). With positive side exponents the
cover exponent is \(K_0=\tfrac12\sqrt{\beta\gamma}\).
A zero side exponent can be handled by enlarging the box
through \(X^\varepsilon\). The cover leaves
\(O_\varepsilon(X^{2K_0+\varepsilon})\) exceptional points and
\(O_\varepsilon(X^{K_0+\varepsilon})\) proper intersections
of bounded degree.
[Browning–Verzobio, Theorem 1.1 and (5)](https://arxiv.org/html/2408.11453v2#S1)

On a geometrically irreducible line or conic meeting the
positive affine open, write the normalized homogeneous
coordinate forms \(L_0,L_x,L_y,L_z\) of degree \(d\le2\)
on \(\mathbb P^1\). They satisfy
\[
 U L_x^3 L_0+V L_y^4+W L_z^4=0 .
\]
If \(L_y/L_z\) is nonconstant, each of the four distinct
roots of \(Vt^4+W\) has a nonempty fiber. At a point of such a
fiber outside the zeros of \(L_0\) and the common zeros of
\(L_y,L_z\), exactly one distinct quartic factor vanishes,
with order at most two. The same point must be a zero of
\(L_x\), giving order at least three on the other side.
Thus all four disjoint fibers meet only the zeros of \(L_0\)
or the common zeros of \(L_y,L_z\). There are at most
\(d+(d-1)\le3\) such points, a contradiction.
The ratio \(y/z\) is therefore constant. Primitivity fixes the
positive pair \(y,z\), and the signed equation fixes \(x\).
Curves of degree at least three contribute at most
\(X^{H/3+\varepsilon}\) apiece, where \(H=\max(\beta,\gamma,\xi)\),
by the same uniform curve bound as above. A component not
defined over \(\mathbb Q\) has only boundedly many rational
points, by intersecting it with a distinct conjugate.

Writing \(f=L-\beta-\gamma-\xi\) for the frozen mass, this gives
\[
 D\le f+\max\{2K_0,K_0+\beta/3,K_0+\gamma/3,K_0+\xi/3\}.
\]
For the exact certificate take \(K_0\le9\beta/40+5\gamma/18\)
by the arithmetic–geometric mean inequality. At the unit-budget
profile \((\beta,\gamma,\xi)=(10,8,8)/37\) with \(f=11/37\),
the bound is \(359/666<27/50\).

The cover is stronger if the frozen quartic cofactor has
positive numerical weight. Put \(n=|W|\) and shift \(y\) by
one in the side polynomial, using
\(g(x,y)=Ux^3+V(y+1)^4\). Equivalently the original \(y\) is
one plus the new box coordinate. The congruence holds
modulo \(n\) and \(g(0,0)=V\) is coprime to \(n\) by
primitivity. The leading monomial is still \(z^4\), and
the side-polynomial weight in (3) is at most
\(\log X+O(1)\), since \(3\beta,4\gamma\le1\).
If \(|W|\ge X^{1/8-o(1)}\), equation (5) bounds the
nonnegative cover exponent by \((15/16)K_0\).
Use the same affine tangent multiplied by \(15/16\)
in the four terms, with a separate branch for weight
at most \(1/8\). At the denominator-\(37\) profile either
tall quartic cofactor has weight \(5/37>1/8\).
The strengthened bound there is \(5583/10656\).

#### Using the short cofactor to count in a lattice

Suppose \(Ux^3\) is the short term, with actual weight \(h\),
and the other two terms are tall. The same frozen mass
\(f=L-\beta-\gamma-\xi\) applies. Put
\[
 s=3\beta+\gamma+\xi-1 .
\]
The stronger bound in this order is \(D\le f\) if \(s\le0\).
If \(s\ge0\) it is
\[
 D\le f+\max\left\{\frac{5s}{8},
                       \frac{7s}{24}+\frac{\beta}{3}\right\}.
 \tag{DS}
\]
In particular it uses only the retained coordinates,
even when the short term contains an aggregated tail.

For fixed coefficients, reduce modulo \(|U|\).
The pair \((y,z)\) is primitive, and both coordinates
and both quartic coefficients are units modulo \(|U|\).
It follows that \(y/z\) lies in one of \(X^\varepsilon\)
classes, since a unit has at most eight fourth roots at
any prime power. Each class \(y\equiv j z\bmod |U|\)
is a rank-two lattice of covolume \(|U|\).
The argument includes \(|U|=1\), with the full lattice.

Let \(\theta=(|W|/|V|)^{1/4}\).
In the dyadic box both \(y\) and \(\theta z\) have size
comparable to \(X^\gamma\) since both corresponding
terms have size comparable to \(X\). The equation gives
\(\big||V|y^4-|W|z^4\big|\ll X^h\), for either choice
of signs. The mean value theorem then bounds
\(\lvert y-\theta z\rvert\ll X^{\gamma+h-1}\).
Together with \(\lvert z\rvert\ll X^\xi\) this bounds
the pair in a symmetric strip of area \(O(X^{\gamma+\xi+h-1})\).
The quotient of its area by the covolume is
\(O(X^s)\) because \(\lvert U\rvert\asymp X^{h-3\beta}\).

For each class send the strip by a real linear map
to a fixed square. If its lattice points span a line,
the only possible positive primitive pair lies on its
primitive integer generator. It determines \(x\).
Otherwise take a Gauss reduced basis of the transformed
lattice. Both basis vectors have bounded length because
there are two independent vectors in the fixed square.
The angle between the basis vectors has sine at least
\(\sqrt3/2\). The integer coefficients of a lattice point
in the square therefore lie in intervals of sizes
\(O(1/\lVert v_1\rVert)\) and \(O(1/\lVert v_2\rVert)\).
Taking each box side at least two, their product is
\(O(X^s)\) by the covolume formula. Thus \(s\le0\)
gives \(O(1)\) choices of the pair per class.
For \(s>0\) the same rectangle also gives the elementary
bound \(X^{s+\varepsilon}\) for the pair, and it determines
\(x\). Thus \(D\le f+s_+\) is available alongside (DS);
in particular it keeps the bound continuous at \(s=0\).

For \(s>0\), write these two side exponents as \(u,v\)
with \(u+v\le s+o(1)\). The inverse basis substitution
in the original coordinates is integral and invertible
over \(\mathbb Q\). The new surface is still absolutely
irreducible, and its low-degree components contain at
most one admissible point by the preceding proof.
The monomial \(x^3\) is unchanged, so the largest
monomial scale in Browning–Verzobio, Theorem 1.1 with
modulus one, is at least \(X^{3\beta}\).
The cover exponent is at most
\(\sqrt{uv/3}\le s/(2\sqrt3)+o(1)\).
One can enlarge all sides by \(X^\delta\) and then
send \(\delta\) to zero. This handles zero exponents
and bounds the auxiliary degrees at each fixed
\(\delta>0\). The maximum side exponent is at most
\(\max(\beta,s)+o(1)\). Using the uniform count for
curves of degree at least three and summing the classes
yields
\[
 D\le f+\max\left\{\frac{s}{\sqrt3},
      \frac{s}{2\sqrt3}+\frac{\max(\beta,s)}{3}\right\}.
\]
Since \(1/(2\sqrt3)<7/24\), this gives (DS).
At the denominator-\(37\) profile the affine bound
is \(407/888\).

## Averaging the square–cube count over residual radicals

For a positive integer \(n\) and \(e\ge2\), let \(r_e(n)\) be
the product of its prime divisors whose valuations are not
divisible by \(e\). Fix an ordered assignment of the three
letters \((S,T,F)\), with \(S\) the square letter and \(T\) the
cube letter. Write
\[
 \mathcal R=r_2(S)r_3(T)\operatorname{rad}(F).
\]
Primitivity makes these factors pairwise coprime. Uniformly for
\(X\ge2\), \(1\le Y\le X^3\), and every \(\varepsilon>0\), the
number of primitive positive triples of height at most \(X\) with
\(\mathcal R\le Y\) is
\[
 \ll_\varepsilon X^\varepsilon Y . \tag{ACT}
\]
In an anatomic box this implies
\[
 D\le L-\sum_{2\mid i}\alpha_{S,i}
        -\sum_{3\mid j}\alpha_{T,j}. \tag{AC}
\]

#### A uniform fixed-coefficient estimate

More generally, fix signed nonzero integers \(A,B,C\) and a
positive integer \(z\), of size bounded by a fixed power of \(X\).
Count solutions with \(x,y>0\) and
\[
 A x^2+B y^3=C z^6 .
\]
For \(z>1\) suppose that the positive absolute values of the
three terms are pairwise coprime. The estimate for \(z=1\) also
holds without this supposition. Define
\[
 P=-AB\,y,\qquad Q=A^2B\,x,\qquad \Delta=A^3B^2 C.
\]
Then \(Q^2=P^3+\Delta z^6\). The replacements are injective
for fixed coefficients and respect arbitrary choices of the
signs. In a primitive \(a+b=c\) triple, put the sign originally
attached to each of the two selected terms into \(A,B\), and
the negative of the sign on the remaining term into \(C\).

Write \(\Delta=d m^2\) for the signed squarefree kernel \(d\)
and \(m>0\). First suppose \(d\ne1\) and put
\(K=\mathbb Q(\sqrt d)\). The algebraic integer
\(\eta=Q+m z^3\sqrt d\) has norm \(P^3\), so that
\[
 (\eta)(\bar\eta)=(P)^3 .
\]
At a prime ideal above \(p\nmid 2\Delta\), simultaneous
divisibility of \(\eta\) and \(\bar\eta\) would imply
\(p\mid Q\) and \(p\mid z\), by taking their sum and difference.
For \(z>1\) the latter is incompatible with \(p\mid Q\), since
\(p\nmid AB\) and \(\gcd(x,z)=1\) by primitivity. For \(z=1\),
the difference is already a unit there. Thus any valuation of
\((\eta)\) outside the primes above \(2\Delta\) is divisible
by three. At the remaining primes reduce the valuations
modulo three. This gives
\[
 (\eta)=\mathfrak b\,\mathfrak a^3
\]
with \(\mathfrak a,\mathfrak b\) integral and at most
\(3^{2\omega(2\Delta)}\) possibilities for \(\mathfrak b\).
For any such \(\mathfrak b\), the class equation
\([\mathfrak a]^3=[\mathfrak b]^{-1}\) has at most
\(h_3(K)=|\mathrm{Cl}(K)[3]|\) solutions.

Fix such a class, choose an integral ideal
\(\mathfrak a_0\) representing it, and choose a generator
\(\eta_0\) of \(\mathfrak b\mathfrak a_0^3\). Since
\(\mathfrak a=(\gamma)\mathfrak a_0\) is integral,
\(\gamma\in\mathfrak a_0^{-1}\). The principal ideal identity
gives \(\eta=\eta_0 u\gamma^3\) with a field unit \(u\).
After absorbing a cube into \(\gamma\), there are boundedly
many choices of \(u\) by the unit theorem in quadratic fields.
Write \(\gamma=s\omega_1+t\omega_2\) in a fixed
basis of the rank-two lattice \(\mathfrak a_0^{-1}\).
The difference with the conjugate gives the equation
\[
 F(s,t)=2m z^3,\qquad
 F(S,T)=\frac{\eta_0 u(S\omega_1+T\omega_2)^3
 -\overline{\eta_0 u}(S\bar\omega_1+T\bar\omega_2)^3}{\sqrt d}.
\]
Each coefficient before subtraction belongs to
\(\eta_0 \mathfrak a_0^{-3}=\mathfrak b\subseteq\mathcal O_K\).
For any \(\xi\in\mathcal O_K\), \((\xi-\bar\xi)/\sqrt d\)
is an integer, so \(F\) has integral coefficients. The two
linear forms \(S\omega_1+T\omega_2\) and
\(S\bar\omega_1+T\bar\omega_2\) are not proportional, since
\(\omega_1,\omega_2\) are independent over \(\mathbb Q\).
The difference of their nonzero scalar multiples cubed has
three distinct geometric linear factors.

For any integral binary cubic of nonzero discriminant, the
number of pairs at a fixed nonzero right side \(M\) is
\(O_\varepsilon(|M|^\varepsilon)\), uniformly in the form.
If it is irreducible, first remove \(g=\gcd(s,t)>0\), for
which \(g^3\mid |M|\). For each \(g\), Evertse's Theorem 1
with the infinite place and primes dividing \(M\) bounds the
number of \(S\)-unit cosets of solutions of the corresponding
equation. At the fixed right side \(M/g^3\), each coset
contributes at most one rational pair, since a common rational
scalar must have cube one. The divisor bound treats \(g\), and
the exponential in \(\omega(M)\) is absorbed in
\(|M|^\varepsilon\). If the cubic is reducible, it has a
primitive integral linear factor \(L\) and can be written
\(L F_2\) with \(F_2\) integral quadratic by Gauss's lemma.
For each signed divisor \(\ell\mid M\) restrict \(F_2\) to
the line \(L=\ell\). This restriction is nonconstant,
since otherwise homogeneity and \(\ell\ne0\) would give
\(F_2\propto L^2\), contradicting the distinct factors.
Thus there are at most two points per signed divisor.
[Evertse, Theorem 1](https://pub.math.leidenuniv.nl/~evertsejh/95-Thue-Mahler.pdf)

If \(d=1\), use the integer factors
\((Q+m z^3)(Q-m z^3)=P^3\) instead. They are coprime outside
\(2\Delta\), by the same sum and difference test. Their signed
cube-free parts \(u_1,u_2\) are therefore supported on \(2\Delta\)
and have \(X^{o(1)}\) choices. With these fixed one obtains the
integral cubic equation \(u_1 s^3-u_2 t^3=2m z^3\), of nonzero
discriminant, and applies the same uniform bound. We assign
the value 1 to \(h_3\) for \(d=1\).

This proves, for each fixed choice of \(A,B,C,z\), an upper
bound \(X^\varepsilon h_3(\mathbb Q(\sqrt d))\) on the number
of root pairs, with the split case understood as above. All
right sides derived here have polynomial size in \(X\).
The coefficients and basis of the Thue form need no height
bound because the estimate is uniform in that form.

#### Summing the fixed-coefficient estimate

For (ACT), in each prime valuation strip group all degrees
divisible by two on \(S\) into the square root, and all degrees
divisible by three on \(T\) into the cube root. Put the other
full prime powers of \(S,T\) into \(|A|,|B|\), respectively,
and put the complete remaining term \(F\) into \(|C|\);
use \(z=1\). This produces exactly
\(\operatorname{rad}(|ABC|)=\mathcal R\), and a triple is
determined by the coefficient choice and the two roots. Using
exact valuation residues instead of whole powers for the two
cofactors also gives the same radical.

At each fixed \(\mathcal R\), there are \(X^{o(1)}\) coefficient
assignments, uniformly for \(\mathcal R\le X^3\). Indeed for
any fixed \(\delta>0\) the number of positive integers at most
\(X^{O(1)}\) supported on primes dividing \(\mathcal R\) is
at most
\[
 X^{O(\delta)}
 \prod_{p\mid\mathcal R}(1-p^{-\delta})^{-1}
 =X^{O(\delta)+o(1)}.
\]
For the last equality, the logarithm of the product is
\(O_\delta(\sum_{j\le\omega(\mathcal R)}(j+1)^{-\delta})
=o(\log X)\) using the increasing order of its primes and
\(\omega(\mathcal R)\ll\log X/\log\log(2X)\).
Choose \(\delta\) sufficiently small for the desired loss.

The integer \(r=|d|\) is a divisor of \(\mathcal R\).
At most \(Y/r\) positive integers \(\mathcal R\le Y\) are
multiples of \(r\). The signed squarefree first moment
\[
 \sum_{\substack{1\le r\le T\\ r\ {\rm squarefree}}}
  \bigl(h_3(\mathbb Q(\sqrt r))
       +h_3(\mathbb Q(\sqrt{-r}))\bigr)\ll T
\]
is the main-term consequence of Davenport–Heilbronn,
or of Taniguchi–Thorne Theorem 1.2. Each radicand has
fundamental discriminant of absolute value at most \(4r\),
and \(d=1\) is counted by the split convention. Dyadic
summation gives
\[
 \sum_{\substack{r\le Y\\ r\ {\rm squarefree}}}
 \frac{h_3(\mathbb Q(\sqrt r))+h_3(\mathbb Q(\sqrt{-r}))}{r}
 \ll 1+\log Y .
\]
Multiply by \(Y\) and absorb the log factor into \(X^\varepsilon\).
This proves (ACT), including uniformity in \(Y\).
[Taniguchi–Thorne, Theorem 1.2](https://arxiv.org/pdf/1102.2914)

For an anatomic box, the upper exponent of \(\mathcal R\)
is precisely \(L-m_2(S)-m_3(T)\), where \(m_e(Q)\)
denotes the mass of coordinates on letter \(Q\) with degree
divisible by \(e\). Taking the standard vanishing box loss
gives (AC). With exact degree coordinates through six it
is safe to retain only \(2,4,6\) on \(S\), and \(3,6\) on \(T\).
With exact coordinates through four retain only \(2,4\) on
\(S\) and \(3\) on \(T\). Unexposed tail coordinates may be
dropped from these subtractions. For a deterministic order
the six ordered pairs \((S,T)\) are
\((a,b),(a,c),(b,a),(c,a),(b,c),(c,b)\).
Each projected row is
\(D-L+m_2(S)_{\mathrm{early}}+m_3(T)_{\mathrm{early}}\le0\),
with no height or short-letter hypothesis. In the denominator
\(37\) witness take \((S,T)=(b,a)\), giving
\(D\le1-(8+10)/37=19/37<13/25\).

For the full \((2,3,6)\) residual radical \(R\), one may
instead retain the exact square-free, cube-free and sixth-free
coefficients \(A,B,C\) and apply the same proof with \(z\)
fixed. The sixth root satisfies \(z\le X^{1/6}\). Summing
the coefficient and signed squarefree-kernel estimates gives
\[
 \#\{\text{primitive triples of height }\le X:\ R\le Y\}
   \ll_\varepsilon X^{1/6+\varepsilon}Y .
\]
In an anatomic box the extra exponent \(1/6\) can be replaced
by the mass of the degrees divisible by six on its assigned
letter. Indeed the part of \(z\) outside primes of \(R\) has
that many choices, since it is determined by the coordinates
of degree divisible by six. Its part supported on primes of
\(R\) has \(X^{o(1)}\) choices by the same Rankin estimate.
Adding this mass to the full residual exponent recovers (AC).

The argument for the full signature also counts all heights
once the canonical sixth root satisfies \(z\le Z\). Its
square-free, cube-free and sixth-free coefficients are
bounded by powers of \(Y\), since their prime supports divide
the residual radical and their valuations are bounded by
five. The right sides in the binary form equations are
therefore polynomially bounded in \(YZ\), while Evertse's
bound needs no a priori restriction on the sizes of the
solutions. The same summation gives
\(O_\varepsilon(Y^{1+\varepsilon}Z^{1+\varepsilon})\)
for \(Y,Z\ge1\).

#### Square and odd prime powers

The same descent applies with the cube replaced by a fixed odd prime
\(\ell\ge5\). For the fixed equation \(Ax^2+By^\ell=C\), use
\[
 P=-AB y,\qquad
 Q=A^{(\ell+1)/2}B^{(\ell-1)/2}x,\qquad
 \Delta=A^\ell B^{\ell-1}C.
\]
Then \(Q^2=P^\ell+\Delta\). In the nonsplit quadratic field the
ideals of \(Q\pm\sqrt\Delta\) are coprime outside \(2\Delta\) and
their product is \((P)^\ell\). Reducing the exceptional valuations
modulo \(\ell\) gives at most \(\ell^{2\omega(2\Delta)}\) ideal
choices. The class equation has at most
\(h_\ell(\mathbb Q(\sqrt{\operatorname{sf}(\Delta)}))\) choices,
where \(\operatorname{sf}\) denotes the signed odd-valuation kernel.
Taking the difference with the conjugate yields a binary form
which is a difference of nonzero multiples of two independent
linear forms raised to the \(\ell\)-th power. Its factors are
therefore geometrically distinct. The split case again uses the
two integer factors and the same binary-form estimate.

For completeness, an integral binary form of odd fixed degree
with distinct geometric factors has
\(O_{\ell,\varepsilon}(|M|^\varepsilon)\) integer points at a fixed
nonzero value \(M\). If it has an irreducible factor of degree at
least three, choose its value from the signed divisors of \(M\)
and apply Evertse as in the cubic case. If all factors have degree
at most two, odd total degree supplies a linear factor \(L\).
Its value is a signed nonzero divisor of \(M\), and the product of
the other factors restricts to a nonconstant polynomial on the
resulting line. In coordinates whose first linear form is \(L\),
homogeneity would make a constant restriction at \(L\ne0\) possible
only for a power of \(L\), contrary to distinctness.

Use \(\mathcal R=r_2(S)r_\ell(T)\operatorname{rad}(F)\).
The squarefree kernel of \(\Delta\) has absolute value \(r\)
dividing \(\mathcal R\), and only the odd valuation primes on
\(S\) and \(F\) can contribute to \(r\). In particular,
\[
 r\le\operatorname{rad}_{\mathrm{odd}}(S)
       \operatorname{rad}_{\mathrm{odd}}(F).
\]
Here \(\operatorname{rad}_{\mathrm{odd}}\) denotes the product of
the primes with odd valuation. The first moment of \(h_\ell\)
for the two signs of squarefree \(d\) up to \(T\) is
\(O_{\ell,\varepsilon}(T^{3/2-3/(2\ell+2)+\varepsilon})\).
As before, at most \(Y/r\) possible residual radicals are multiples
of \(r\), and the choices of coefficients at fixed residual
radical cost \(X^{o(1)}\). Thus the number of triples with
\(\mathcal R\le Y\) and \(r\le Z\) is
\(O_{\ell,\varepsilon}(X^\varepsilon Y Z^{1/2-3/(2\ell+2)})\),
for \(1\le Y,Z\le X^3\), by dyadic summation.
[Heath-Brown–Pierce, Theorem 1.1](https://arxiv.org/pdf/1409.3177)
· [Chan–Koymans, Theorem 1.2](https://arxiv.org/pdf/2505.00611)

For \(\ell=5\) the resulting affine inequality is
\[
 D\le L-m_2(S)-m_5(T)+\tfrac14\bigl(L_S+L_F-m_2(S)-m_2(F)\bigr).
 \tag{AQ}
\]
Both the removed mass and the upper bound for the odd mass
remain valid after dropping unexposed valuations from the
\(m_e\) sums. With six exact indices this yields six further
rows, one for each assignment of square, fifth-power and
fixed letters.

## Small-coefficient determinant-one families

For \(n>0\), put
\[
 o(n)=\prod_{\substack{p\mid n\\v_p(n)\ge3,\ v_p(n)\ {\rm odd}}}p.
\]
There is an absolute \(\eta_0>0\) such that the following
restricted part of the \(abc\) exceptional set satisfies
\[
 N^{\mathrm{det}=1}_{\lambda}(X)
 \ll_{\lambda,\varepsilon}
 X^{\lambda/3-\eta_0\lambda^3+\varepsilon}
 \qquad (0<\lambda\le1).
 \tag{DU}
\]
Here \(N^{\mathrm{det}=1}_{\lambda}\) counts primitive positive
\(a+b=c\le X\) with \(\operatorname{rad}(abc)<c^\lambda\),
\(o(a)o(b)o(c)\le c^{\eta_0\lambda^3}\), and at least one
positive pair \(u,v\le c^{\eta_0\lambda^3}\) satisfying
\(\lvert va-ub\rvert=1\).
The same determinant-one relation gives adjacent scaled
values for any pair of letters, including the pairs with \(c\).

### A support correction in a consecutive pair

Consider coprime consecutive positive integers \(m,m+1\)
of size comparable to \(M\), with radicals at most \(A,B\).
Write \(m=dx^2\), \(m+1=ey^2\) with \(d,e\) squarefree and
localize \(de\asymp W\). The Pell identity
\[
 (2m+1)^2-de(2xy)^2=1
\]
gives \(O(W\log M)\) pairs in this box. A discriminant is
nonsquare because two positive squares cannot differ by one.
Each discriminant gives only \(O(\log M)\) possible powers
of the fundamental unit, and \(2m+1\) determines the pair.

The radical identity for the square part is
\[
 \operatorname{rad}(x)=\operatorname{rad}(m)o(m)/d .
\]
Consequently, with the additional condition
\(o(m)o(m+1)\le T\), a second bound for the box when \(W\ll M\)
is
\[
 \ll_\varepsilon M^\varepsilon AB T/W.
 \tag{SP}
\]
To prove this, localize \(d,e,o(m),o(m+1)\) to dyadic
sizes \(D,E,T_1,T_2\), with \(DE\asymp W\) and
\(T_1T_2\ll T\). The ordinary radical-box estimate gives
\(O_\varepsilon(M^\varepsilon A T_1/D)\) choices of \(x\)
and \(O_\varepsilon(M^\varepsilon B T_2/E)\) choices of \(y\),
uniformly over all coefficient supports in the box. That
estimate follows by choosing the radical and applying a
Rankin bound to integers with the fixed prime support. For
fixed \(x,y\), coprimality forces coefficient solutions of
\(ey^2-dx^2=1\) into one progression of step
\((y^2,x^2)\) on \((d,e)\). Since \(y^2\asymp M/E\) and
\(x^2\asymp M/D\), the progression contributes \(O(1+W/M)\)
points. Sum over the dyadic sizes.

The fixed-support estimate in Zeng's Lemma 2.3 applies to
each fixed kernel \(d\). In the application in Lemma 4.1,
equation (4.5) needs a count over all \(d\) in the dyadic
interval. The fixed-support estimate alone leaves that
step unproved. For example, for \(d=x=p\) over primes
\(p\asymp U\), the radical outside \(d\) is one for every
\(x\), while the union of eligible \(x\) has size
\(\asymp U/\log U\) at heights \(\asymp U^3\).
This example diagnoses the use of the separate support
inequality, without giving a counterexample to the count
with the second consecutive integer also constrained.
The factor \(T\) in (SP) pays for the changing support.
[Zeng, Lemmas 2.3 and 4.1](https://hexagonmath.org/pdf/2611.00611v1)

### Counting the determinant-one family

Write \(\rho=\lambda/3\). Let \(c_D>0\) be the constant in
Diao's Theorem 1.3, and set
\(\delta=\min(1,c_D)\rho^3/4096\). Prove the claim dyadically
with \(Y/2<c\le Y\), allowing \(u,v\le H=Y^\delta\) and
\(o(a)o(b)o(c)\le Y^\delta\).

If some coordinate has radical at most \(Y^{\rho-4\delta}\),
there are \(O_\varepsilon(Y^{\rho-4\delta+\varepsilon})\)
choices for its value by the radical-box count. For a fixed
coordinate value \(n\), each of the relevant determinant-one
relations has the shape \(j n-k m=\pm1\) with \(j,k\le2H\).
For a fixed \(k\), either \(\gcd(k,n)>1\) and there are no
possibilities or \(j\) is fixed in one residue class modulo
\(k\). Summing \(O(1+H/k)\) gives \(O(H\log(2H))\) choices
of the other coordinate. The count in this case is
\(O_\varepsilon(Y^{\rho-3\delta+\varepsilon})\).

Otherwise the two coordinates with the smallest radicals
both have radicals at most \(Y^{\rho+2\delta}\), since the
three radicals multiply to at most \(Y^{3\rho}\) and the
smallest exceeds \(Y^{\rho-4\delta}\). The original relation
and its two consequences
\[
 (u+v)a-u c=va-ub,\qquad
 (u+v)b-v c=ub-va
\]
turn the selected pair into adjacent integers after
multiplying each coordinate by a positive integer at most
\(2H\). The scaled pair has height between \(Y/4\) and
\(2HY\). Each scaled pair has \(Y^\varepsilon\) preimages
because the two multipliers divide the two scaled values.
No sum over all pairs \(u,v\) is needed.

The scaled radicals are at most \(2Y^{\rho+3\delta}\).
Also \(o(jn)\le\operatorname{rad}(j)o(n)\), since a new
prime in the odd-overlap product must divide \(j\).
Thus the product of the two scaled overlap factors is at
most \(4Y^{3\delta}\). The resulting \(W\) is at most the
product of the scaled radicals, which is \(o(Y)\)
because \(2\rho+6\delta<1\).
For \(W=Y^s\) with \(s\le\rho-\delta\), use the elementary
Pell count. For \(s\ge\rho+10\delta\), use (SP), whose
exponent is at most \(2\rho+9\delta-s\le\rho-\delta\).

For the remaining range, \(\rho/2\le s\le\rho+10\delta\).
The Pell unit is \(O(Y^{1+\delta})\) and \(de\asymp Y^s\),
so it is at most \((de)^{4/s}\) for large \(Y\).
The fixed parameter \(4/s\) lies in the range of Diao's
Theorem 1.3 for large \(Y\), uniformly on this interval
at fixed \(\lambda\). The number of discriminants is
\[
 O_\lambda(Y^{s-(c_D/16)s^3}).
\]
The subtracted term is at least
\(c_D\rho^3/128\ge32\delta\), so this last exponent is
at most \(\rho-\delta\). Multiply by the logarithmic
Pell fibre count and the divisor bound for preimages.
Dyadic summation proves (DU) with
\(\eta_0=\min(1,c_D)/(27\cdot4096)\).
[Diao, Theorem 1.3](https://arxiv.org/html/2610.28001v1#S1.SS3)

## Averaged full \((2,3,6)\) count with a Mordell height lower bound

Fix an orientation of the full signature. Write the primitive
equation in canonical coefficient form
\[
 A x^2+B y^3=C z^6,\qquad \Delta=A^3B^2C,
 \qquad U=-AB\,y/z^2,\qquad V=A^2B\,x/z^3.
\]
The coefficients are signed square-free, cube-free and
sixth-free, respectively, in absolute value.
Then \(V^2=U^3+\Delta\). Pairwise coprimality of the complete
terms implies that \(\Delta\) has valuations \(1,\ldots,5\)
at the primes dividing the full residual radical \(R\),
and is a unit at the other primes. Set
\[
 R_e=\prod_{v_p(\Delta)=e}p,\qquad R=\prod_{e=1}^5R_e,
 \qquad R_2^{\rm odd}=R_1R_3R_5.
\]
Thus the signed squarefree kernel of \(\Delta\) is
\(d=\operatorname{sgn}(\Delta)R_2^{\rm odd}\).

Consider a height cutoff \(c\le X\). For any rational point
\(P=(U,V)\) on the Mordell curve away from \(U V=0\), put
\(t=U^3/\Delta\). The isomorphism to the fixed curve
\(E_1: v^2=u^3+1\) sends \((U,V)\) to
 \((U/\delta^2,V/\delta^3)\) for a compatible root \(\delta^6=\Delta\).
The scaled abscissa has cube \(t\), so by the standard
height comparison on the fixed curve,
\[
 \widehat h(P)=\tfrac16 h(t)+O(1).
 \tag{M1}
\]
[Silverman, VIII.9.1 and VIII.9.3](https://www.math.ens.psl.eu/~obenoist/refs/Silverman.pdf#page=264)

The error is absolute. Heights agree through this isomorphism:
under the limiting definition, the bounded difference between
the heights of two abscissas differing by a fixed scalar
vanishes after multiplication by \(4^{-n}\) on \([2^n]P\).
The points with \(U=0\) or \(V=0\)
are torsion: on \(E_1\), the first has order three
(by doubling or the tangent formula) and the second has
order two. The point at infinity is also torsion.
The number of rational torsion points is uniformly bounded,
since their images on \(E_1\) have bounded degree and
bounded coordinate height, hence lie in a finite set by
Northcott. A point from a primitive triple has
\(t=-B y^3/(Cz^6)\), a signed ratio of complete terms of size
at most \(X\). The signed rational ratio determines those two
coprime complete terms and hence the full primitive triple.
Therefore
\[
 \widehat h(P)\le \tfrac16\log X+O(1).
 \tag{M2}
\]

For a non-torsion point on the same curve, write \(t=P_0/Q_0\)
in lowest terms with \(Q_0>0\). At a prime \(p\mid R_e\),
\[
 v_p(t)\equiv-e\pmod3,\qquad
 v_p(t+1)\equiv-e\pmod2 .
\]
Consequently at least one of \(P_0,Q_0,P_0+Q_0\) is divisible
by \(p\). The cases are more precise.
If \(p\mid Q_0\), both valuations of \(t,t+1\) coincide
and the denominator has exponent congruent to \(e\bmod6\),
hence at least \(e\). If \(p\mid P_0\), then \(v_p(t+1)=0\),
so \(e=2\) or \(4\), and the numerator has exponent at least
\(1\) or \(2\), respectively. If \(p\mid P_0+Q_0\), then
\(v_p(t)=0\), forcing \(e=3\), with exponent of
\(P_0+Q_0\) at least one. In particular, primes in \(R_1R_5\)
all divide the denominator. Put \(\tau_e=\log R_e/\log X\)
and \(\rho=\sum_e\tau_e\). The following are lower bounds
for \(h(t)/\log X\), up to \(O(1/\log X)\):
\[
 \frac{\rho}{3},\quad \tau_1+5\tau_5,\quad
 \frac{2(\tau_2+2\tau_4)+\tau_1+5\tau_5}{3},\quad
 \frac{3\tau_3+\tau_1+5\tau_5}{4},\quad
 \frac{2(\tau_2+2\tau_4+\tau_3)+\tau_1+5\tau_5}{5}.
 \tag{M3}
\]
For the first bound use \(R\mid P_0Q_0(P_0+Q_0)\).
For the remaining bounds use, respectively, the denominator
alone, the numerator plus half of the denominator logarithm,
the sum numerator \(P_0+Q_0\) plus one third of the denominator
logarithm, and both numerator logarithms plus half of the
denominator logarithm. Each numerator and denominator
logarithm is at most \(h(t)+O(1)\). The valuation cases
above give the stated coefficients.

There is one further affine lower bound that uses all five colors
with their full coefficient weights:
\[
 h(t)\ \ge\ \frac{\log|\Delta|}{6}-O(1)
       =\frac{\log X}{6}\sum_{e=1}^5e\tau_e-O(1).
 \tag{M3a}
\]
The valuation cases imply
\(\lvert\Delta\rvert\mid Q_0P_0^2(P_0+Q_0)^3\).
Taking logarithms gives (M3a) because \(\log|P_0|\) and
\(\log Q_0\) are at most \(h(t)\), while
\(\log|P_0+Q_0|\le h(t)+\log 2\).
The Mordell coset counts can also take \(u\) from this lower bound.

Fix a real constant \(u>0\). Consider those coefficient
slices where one of (M3) is at least \(u\), and choose an
integer \(m\ge2\) with \(m^2u>2\).
The Mordell–Weil lattice modulo torsion has only boundedly
many counted vectors per class modulo \(m\).
Indeed the difference of any two distinct vectors in one
class has height at least \(d_0^2=m^2 u\log X/6-O(1)\),
by (M1) and (M3). Each vector has norm at most
\(\sqrt{H_0}\) with \(H_0=\log X/6+O(1)\) by (M2).
For \(n\) vectors in the class, the sum of squared distances
over pairs is \(n\sum_i\|P_i\|^2-\|\sum_iP_i\|^2\le
n^2 H_0\). It is at least \(n(n-1)d_0^2/2\). Thus
\(n\le d_0^2/(d_0^2-2H_0)\), which is uniformly bounded
for sufficiently large \(X\) by \(m^2u>2\).
This gives \(O_{m,u}(m^{r_\Delta})\) points per slice,
where \(r_\Delta=\operatorname{rank}E_\Delta(\mathbb Q)\).

### The 3-class group in the rank bound

Over \(K=\mathbb Q(\zeta_3,\sqrt d)\), the degree-three CM
isogeny \(1-\zeta_3\) is the difference from the automorphism
\((U,V)\mapsto(\zeta_3 U,V)\). Its kernel consists of the point
at infinity and the two points with \(U=0\). It is separable
in characteristic zero, hence has degree three and rational
kernel over \(K\), which also contains
\(\mu_3\). Its Kummer image consists of classes unramified
outside the primes of \(K\) over \(6R\), since at other
primes the elliptic curve has good reduction and the
isogeny is étale. The rank of the \(S\)-unit group is
\(O(1+\omega(R))\). The remaining Kummer contribution
is bounded by \(\mathrm{Cl}(K)[3]\), since the
3-rank of the localized class-group quotient does
not exceed the 3-rank of the full class group.

If \(K\) is biquadratic, norm and extension, using that
2 is invertible in the 3-primary group, identify its
3-torsion as the product of the 3-torsion of its three
quadratic subfields. To see this, decompose the
\(\mathbb F_3[\operatorname{Gal}(K/\mathbb Q)]\)-module into
characters. The part invariant under each index-two
subgroup identifies under norm with the 3-torsion of the
corresponding quadratic field, and the trivial character
vanishes by taking norm to \(\mathbb Q\). The fixed
\(\mathbb Q(\zeta_3)\) factor is trivial. For the other
two factors, the Scholz reflection relation gives
\[
 h_3(\mathbb Q(\sqrt{-3d}))\le
       3 h_3(\mathbb Q(\sqrt d)).
\]
Here multiplying the radicand by a rational square does
not change the field. The cases where \(K\) is only
quadratic have a fixed field and can be absorbed in the
constant. For negative \(d\), apply reflection to its positive reflected
radicand.
[Heath-Brown–Pierce, discussion after Theorem 1.2](https://arxiv.org/pdf/1409.3177)

The \(\mathbb Q(\zeta_3)\)-span of a rational basis of
\(E_\Delta(\mathbb Q)\otimes\mathbb Q\) has dimension
\(r_\Delta\). Complex conjugation fixes these basis
vectors. Subtracting the conjugate of any purported
\(\mathbb Q(\zeta_3)\)-linear relation first gives a
rational relation for the coefficients of the
non-real generator, and then a rational relation for
the other coefficients. Since \(1-\zeta_3\) has ring
norm 3, the free CM module over \(K\) contributes at
least \(3^{r_\Delta}\) to the Kummer quotient. Consequently
\[
 3^{r_\Delta}\ll C^{1+\omega(R)}h_3(\mathbb Q(\sqrt d))^2.
\]
With \(m\) fixed and \(k\ge2\log_3 m\), the bound on the
number of points in each coefficient slice is therefore
\(X^{o(1)}h_3(\mathbb Q(\sqrt d))^k\).

### Averaging and the finite affine form

Suppose \(R\le Y\), \(|d|\le Z\), all considered coefficient
slices satisfy the fixed height lower threshold \(u\)
above, and \(Y,Z\le X^3\). There are
\(X^{o(1)}\) coefficient assignments per \(R\). For fixed
\(|d|=r\) there are at most \(Y/r\) possible \(R\), since
\(r\mid R\). For real \(k\ge1\), the available bounds
for the signed squarefree-radicand moments are
\[
 \sum_{\substack{1\le r\le H\\r\ {\rm squarefree}}}
   \left(h_3(\mathbb Q(\sqrt r))^k+
         h_3(\mathbb Q(\sqrt{-r}))^k\right)
 \ll_{k,\varepsilon} H^{f(k)+\varepsilon},
\]
where one may use
\[
 f(k)=
 \begin{cases}
 (5k+13)/18, & 1\le k\le4,\\
 11/6+(1597/5000)(k-4), & k\ge4.
 \end{cases}
\]
The first range and fourth moment are Heath-Brown–Pierce.
For the second range, bound all remaining powers by the
pointwise estimate of Chan–Koymans. By dyadic partial
summation with weight \(1/r\), the full number of triples
under these conditions is at most
\[
 O_{u,m,k,\varepsilon}\bigl(X^\varepsilon Y
                            Z^{f(k)-1}\bigr).
 \tag{M4}
\]
Degenerate radicand \(d=1\) costs the same with class-group
size interpreted as one and the field fixed.
[Heath-Brown–Pierce, Theorem 1.2 and Corollary 1.4](https://arxiv.org/pdf/1409.3177)
· [Chan–Koymans, Theorem 1.1](https://arxiv.org/pdf/2505.00611)

For the projected exponent program, take a lower bound on
each \(\tau_e\) using only retained exact valuation indices.
On the square letter \(S\), every odd index contributes to
\(\tau_3\). On the cube letter \(T\), indices congruent to
\(1,2\bmod3\) contribute to \(\tau_2,\tau_4\), respectively.
On the sixth-power letter \(F\), an index congruent to
\(e\bmod6,\ e\ne0\), contributes to \(\tau_e\).
All coefficients in (M3) are nonnegative, so these are
valid lower bounds. On the sign branch where one such
affine lower bound is at least a fixed \(u\), (M4) gives
\[
 D\le \rho_{\mathrm{up}}+(f(k)-1)\sigma_{\mathrm{up}},
\]
where \(\rho_{\mathrm{up}}\) bounds the full residual mass
by dropping unexposed removable valuations and
\(\sigma_{\mathrm{up}}=L_S+L_F-E_S-E_F\) bounds the odd
residual mass, with \(E_Q\) the retained even-index mass.
The strict inequality \(m^2u>2\) absorbs the vanishing
box endpoint loss as well as the bounded fixed-curve
height error. The branch constants and \(k\) can all be
rational, with \(m\) integer.

For example, \(u=1/10\), \(m=5\), \(k=3\) are valid
since \(25>20\) and \(3^3>5^2\). The slope \(f(k)-1\)
is \(5/9\). If the retained \(e=5\) mass is at least
\(1/20\), the denominator row gives \(u=1/4\).
Then \(m=3,k=2\) give correction \(5\sigma_{\mathrm{up}}/18\).
This saves over freezing the sixth root when its mass
exceeds that correction.

### A one-class-group bound from the rational three-isogeny

The isogeny improves the same height-window method in one
additional range. If \(u>2/3\), then (M4) holds with
\(f(k)-1=0\). Thus in this range the full \((2,3,6)\)
residual radical has coefficient one, even when both parity
and cube defects have positive mass. For example,
\(\tau_5\ge7/50\) suffices, since then (M3) gives \(u=7/10\).

Suppose first that \(\Delta\) is not a square, so that
\(K_0=\mathbb Q(\sqrt{\Delta})=\mathbb Q(\sqrt d)\) is
quadratic. On \(E_\Delta(\mathbb Q)\) define
\[
 f(P)=V+\sqrt{\Delta}\quad\bmod K_0^{\times3},\qquad f(O)=1 .
\]
This is a homomorphism. For a nonvertical line \(V=sU+t\)
meeting the cubic at three points \(P_i\) with multiplicity,
the relation
\(\prod_i (V(P_i)+\sqrt\Delta)=(t+\sqrt\Delta)^3\)
follows by evaluating the monic cubic equation for the
three \(U\)-coordinates at \(U=-(t+\sqrt\Delta)/s\).
The horizontal line case is immediate. Also
\(f(-P)=f(P)^{-1}\) by the norm relation, giving the
homomorphism by the chord and tangent law including
vertical lines.

Let \(\psi:E_{-27\Delta}\to E_\Delta\) be the rational
dual three-isogeny, given off the zero abscissa by
\[
 \psi(x',y')=
 \left(\frac{x'^3-108\Delta}{9x'^2},
       \frac{y'(x'^3+216\Delta)}{27x'^3}\right).
\]
This is the standard three-isogeny followed by the change
of scale \(x\mapsto x/9,\ y\mapsto y/27\). Equivalently
the displayed rational map extends to a morphism of
elliptic curves carrying infinity to infinity. It has
degree three: its abscissa is a degree-three rational
function of \(x'\), and the degree of each abscissa
map from an elliptic curve is two.
[Bhargava–Elkies–Shnidman, Proposition 18 and Remark 19](https://arxiv.org/html/1610.05759#S3.SS1)

If \(f(P)=1\) for a finite rational point, write
\[
 V+\sqrt\Delta=(r+s\sqrt\Delta)^3,\qquad r,s\in\mathbb Q.
\]
Comparing the coefficient of \(\sqrt\Delta\) gives
\(s(3r^2+s^2\Delta)=1\), so \(s\ne0\), and taking the
norm gives \(U=r^2-s^2\Delta\). Then
\((x',y')=(3/s,9r/s)\) lies on \(E_{-27\Delta}\), and the
displayed isogeny maps it to \(P\). Hence
\(\ker f\subseteq\psi E_{-27\Delta}(\mathbb Q)\), with
infinity included.

The image of \(f\) has valuations divisible by three at
all primes of \(K_0\) above \(p\nmid2\Delta\). If \(U\) is
integral at \(p\), then so is \(V\). The algebraic integers
\(V\pm\sqrt\Delta\) locally have unit difference and cube
product, so either valuation is a multiple of three.
If \(U\) is not integral, \(2v_p(V)=3v_p(U)\) forces
\(v_p(V)\) to be a multiple of three, and the unit
\(\sqrt\Delta\) does not change that valuation.
Therefore the usual valuation map for \(S\)-units gives
\[
 [E_\Delta(\mathbb Q):\psi E_{-27\Delta}(\mathbb Q)]
 \le |\operatorname{im} f|
 \ll C^{1+\omega(R)}h_3(K_0).
 \tag{M5}
\]
Here the kernel of the map to the \(S\)-class group is
the \(S\)-units modulo cubes, and the class group
contribution is bounded by the full 3-torsion subgroup.

Canonical heights obey
\(\widehat h_{E_\Delta}(\psi Q')=3
\widehat h_{E_{-27\Delta}}(Q')\). For example, this
follows from the limiting canonical-height definition:
the displayed abscissa of \(\psi\) is a rational map
of degree three on the abscissa line, so the corresponding
naive height changes by a factor three up to a bounded
amount, which disappears in the limit. The same
valuation lower bound (M3) applies on \(E_{-27\Delta}\)
outside the fixed prime 3, since multiplying \(\Delta\)
by \(-27\) leaves its other valuations unchanged. Removing
this prime costs \(O(1)\) in the log height. Thus distinct
counted lattice vectors in one coset of the image of
\(\psi\) have mutual squared distances at least
\(3u\log X/6-O(1)\). For \(3u>2\) the same
sum-of-squared-distances estimate gives boundedly many
vectors per coset. The number of torsion points on both
curves is uniformly bounded as above.

Consequently (M5) gives \(X^{o(1)}h_3(K_0)\) counted
points per coefficient slice. If \(\Delta\) is a square,
the earlier lattice bound with \(h_3=1\) already gives
\(X^{o(1)}\) points. Averaging (M5) with the signed
first moment of Davenport–Heilbronn gives
\(O_{u,\varepsilon}(X^\varepsilon Y)\) triples when
\(R\le Y\) and the lower threshold \(u>2/3\) holds.
The upper bound on the parity radical is unnecessary
in this range.
[Taniguchi–Thorne, Theorem 1.2](https://arxiv.org/pdf/1102.2914)

More generally use the image of \(\psi\circ[m]\) for
any fixed integer \(m\ge1\) with \(3m^2u>2\). The two
curves have the same rank \(r_\Delta\). Its index on
the Mordell lattices is at most (M5) times \(m^{r_\Delta}\),
up to bounded torsion factors. The pointwise class
group weight in the average can therefore take exponent
\(k=1+2\log_3 m\), or any larger rational value, in place
of \(2\log_3 m\) for plain multiplication. The same
bound (M4) follows with this alternative choice.

## A real short-strip estimate with a free low cofactor

Fix \(\eta>0\). There is \(K_0=K_0(\eta)\) such that for \(K\ge K_0\), fixed positive integers \(p,q\), and \(H\to\infty\), the number of positive integer triples \(U,x,y\le H\) satisfying
\[
\left|\frac{Ux^py^q}{C}-1\right|\le H^{-K}
\]
is \(O_{p,q,K,\eta}(H^\eta)\), uniformly in \(C>0\). All constants in the choice of \(K_0\) are independent of \(p,q\). The precision in the small-budget application satisfies \(K\to\infty\) as \(\lambda\to0\). The same proof applies to positive integer monomials in at most three variables and, with a fixed bound on exact monomial fibers, to a signed monomial.

### The determinant at high precision

Restrict to dyadic boxes, and write \(z=(U/B_1,x/B_2,y/B_3)\in[1,2]^3\) with integer dyadic scales \(B_i\le H\). On the enlarged fixed positive box \([1/2,3]^3\), put
\[
\Phi(z)= (B_1 B_2^p B_3^q/C)z_1z_2^pz_3^q.
\]
If the solution set is nonempty, the leading factor in \(\Phi\) is bounded above and below by positive constants in terms of \(p,q\). Solving for \(z_1\) gives a point of the exact level \(\Phi=1\) within \(O_{p,q}(H^{-K})\) of each original point.

Here is the determinant estimate used twice below. Consider approximating points on a semialgebraic set of dimension at most \(d\) in a fixed box, with distance \(O(H^{-k})\) for \(k\) sufficiently large. The complexity of the semialgebraic set may depend on \(p,q\), and the parameterization count may depend on that complexity. By the Yomdin–Gromov algebraic lemma, it has finitely many charts of dimension at most \(d\) with bounded derivatives up to any fixed order. The number and derivative bounds are uniform in its defining coefficients. See [Binyamini–Novikov, Theorem 1 and its semialgebraic version](https://www.wisdom.weizmann.ac.il/~dnovikov/Papers/s40598-021-00176-w.pdf).

Choose \(N>d\) of the ambient coordinates. Use the \(D=\binom{M+N}{N}\) monomials of total degree at most \(M\) in those coordinates. Let \(b_1,\ldots,b_D\) be the least total degrees of distinct monomials in \(d\) parameters, ordered nondecreasingly. Subdivide each chart parameter box into cubes of side \(H^{-\tau}\), where \(\tau>0\) will be fixed. Taylor expansion around a cube center and multilinearity give
\[
\left|\det(z_s^{\alpha})\right|
\ll H^{-\tau\sum_{j=1}^D b_j},
\qquad k\ge\tau b_D .
\]
In the expansion, any term with repeated Taylor monomial vectors is zero. Each error column, including the error from replacing the exact chart point by the nearby rational point, has weight at least \(\tau b_D\). The determinant bound also holds for charts of smaller dimension by padding parameter coordinates. Since \(\sum b_j\asymp_d D^{1+1/d}\) and \(N>d\), we may choose \(M=M(N,d,\tau)\) with \(\tau\sum b_j>MD\). Multiplication by the dyadic monomial scales, costing at most \(H^{MD}\), gives an integer determinant of magnitude less than one. Thus the rational points assigned to each parameter cube lie on an auxiliary polynomial of degree at most \(M\) in the chosen \(N\) coordinates. It can be chosen with coefficients \(H^{O_M(1)}\) by taking a null vector built from a maximal rank minor of the original integer monomial evaluations, or of the normalized rational evaluations with denominators cleared. A fixed-degree primitive factor also has coefficient height \(H^{O_M(1)}\), for example by Kronecker substitution and the univariate factor-height bound.

### First cut and quantitative perturbation to an exact curve

Use the determinant with \(d=2,N=3\) on the exact surface \(\Phi=1\). It gives \(O_{p,q,\tau}(H^{2\tau})\) auxiliary polynomials of bounded degree. Work in the normalized coordinates and pass to the primitive irreducible rational factors \(G\). Each \(G\) has bounded degree and height \(H^{O_\tau(1)}\).

If \(\Phi\) is constant on a geometric component of \(G=0\) meeting the positive box, the original monomial \(Ux^py^q\) is constant there. If an integer point occurs, its value is an integer at most \(H^{1+p+q}\). Each factor divides that value, so this component contributes \(O_{p,q,\eta}(H^{\eta/4})\) by the divisor bound. If the component is nonconstant, the \(2\)-by-\(2\) minors of \((\nabla G,\nabla\log\Phi)\) cannot all vanish identically on it. For a rationally irreducible \(G\) the components are conjugate, so either they are all constant or all nonconstant. In the nonconstant case choose one minor that is not identically zero on \(G\), after multiplying by \(z_1z_2z_3\) to clear denominators. Points where all minors vanish belong to the bounded-degree curve cut out by \(G\) and this numerator. These two polynomials are coprime over \(\mathbb Q\), with height \(H^{O_\tau(1)}\).

At any remaining point, some minor is nonzero and bounded below by a fixed constant times \(H^{-A}\), with \(A=O_\tau(1)\) independent of \(p,q\). Indeed, the cleared logarithmic row has entries \(z_2z_3,pz_1z_3,qz_1z_2\). All minors are then integer polynomials of bounded degree evaluated at rationals of denominators at most \(H\). To keep all derivatives of \(G\) uniformly bounded, divide its coefficients by their maximum absolute value, costing another factor at most \(H^{O_\tau(1)}\). The logarithmic derivatives themselves are uniformly bounded on the enlarged box at fixed \(p,q\).

Keep the unselected coordinate fixed and solve \((G,\log\Phi)=(0,0)\) in the other two real coordinates. At the original point the residual is \(O(H^{-K})\), the inverse Jacobian has norm \(O(H^A)\), and the second derivatives are bounded independently of \(H\). For \(K>2A+1\), the quantitative inverse function theorem gives an exact point at distance \(O(H^{A-K})\). One may obtain the bound directly by the contraction mapping theorem on a ball of radius \(O(H^{A-K})\) using the original inverse Jacobian. The nearby point belongs to \(G=0,\Phi=1\) inside the enlarged box. This intersection has dimension at most one because \(G\) and \(\Phi-1\) share no geometric component in the nonconstant case.

Use the determinant again, now with \(d=1,N=2\) and \(k=K-A\), on this exact intersection. Choose two normalized coordinates that omit a variable on which \(G\) depends. The new polynomial is independent of that variable and coprime to \(G\). The \(O_{p,q,\tau}(H^\tau)\) resulting intersections therefore consist of bounded-degree algebraic curves and points cut out by two bounded-degree integer polynomials of height \(H^{O_\tau(1)}\).

### Counting on the resulting curves

For either type of curve intersection let \(F_1,F_2\) be the two coprime polynomials. An integer linear change of bounded coefficients gives a coordinate \(t\) nonconstant on every geometric curve component. Indeed \(t=z_1+jz_2+j^2z_3\) excludes only two choices of \(j\) per component, at most a bounded number altogether. Let \(v,w\) denote the other coordinates. Eliminate \(w\) and \(v\), respectively, using resultants to obtain nonzero polynomials \(A(t,v),B(t,w)\). If both input polynomials are independent of an eliminated variable, their common plane zeros are finite and the corresponding space components have constant \(t\), precluded by the choice of \(t\). Strip factors only in \(t\) and pass to primitive squarefree parts over \(\mathbb Q(t)\). The resultant and factor-height bounds leave degrees bounded and coefficients \(H^{O_\tau(1)}\). Every curve component dominates the \(t\)-line, so \(A_v B_w\) is generically nonzero on it. The removed factors and derivative zeros account for \(O_\tau(1)\) points.

On the remaining open set the rows \(\nabla A,\nabla B\) have rank two. If the logarithmic Jacobian with third row \(\nabla\log\Phi\) is identically zero on a geometric component, the monomial is constant there and the divisor argument applies. On the other components its zeros cost boundedly many points by Bézout. At all other rational points, clear the three original coordinate denominators from the logarithmic row, and normalize the sizes of the two ordinary polynomial coefficient vectors. As above, the nonzero Jacobian determinant has inverse bounded by \(H^{A'}\) for \(A'=O_\tau(1)\) independent of \(p,q\). The bounded linear change can be absorbed in the constants. Solving \((A,B,\log\Phi)=(0,0,0)\) gives an exact point of this system within \(O(H^{A'-K})\). The Jacobian remains nonzero there, so the point is isolated among the algebraic solutions of \(A=B=\Phi-1=0\). There are boundedly many such isolated points, by Bézout, with the bound allowed to depend on \(p,q\). Since distinct original points are separated in normalized coordinates by at least \(1/H\), an exact point can serve at most one original point for \(K>A'+2\). Thus the nonconstant components contribute a bounded number of points.

Take \(\tau<\eta/8\) and absorb dyadic choices and the divisor loss into the remaining exponent. Choices of degrees, height exponents and the required precision threshold \(K_0\) depend only on \(\eta\). The constants that depend on \(p,q\) occur only after the exponent inequalities and can be absorbed at fixed exponents.

### Making the low coefficient short in anatomic boxes

For a letter with \(l_T>u\), write \(\theta=(l_T-u)/2>0\). On a limiting sequence there are at most \(2/\theta\) coordinates with mass at least \(\theta L/2\). Pass to a subsequence on which this list has fixed length, its normalized masses converge, and each of its degrees is fixed or tends to infinity. Call the latter coordinates heavy. Each has normalized mass at least \(\theta\) in the limit by (A), because its unnormalized weight is eventually at least \(2\lambda\). Choose \(J\) larger than all the fixed degrees on this list. Then every remaining coordinate of degree at least \(J\) has mass below \(\theta L/2\). If any such coordinate has weight at least \(2\lambda\), (A) gives a contradiction. If a collection of them first reaches weight \(2\lambda\), its mass before the last addition is below \(2\lambda/J\), and adding the last coordinate still leaves mass below \(\theta L\) by a fixed limiting margin if \(J\) is large enough. Again this contradicts (A). A fixed \(J\) suffices since \(L\asymp\lambda\). Hence the entire heavy-free complement has weighted height exponent at most \(JL+2\lambda+o(\lambda)\). This selection has the properties of the earlier heavy classification and additionally makes the low cofactor \(X^{O_J(\lambda)}\). All selected heavy degrees tend to infinity, as required by (30).

When one tall term is fixed and the short-to-tall height gap is bounded below, the other tall term is \(D U x^p y^q\), with any additional high factors absorbed into fixed \(D\). Its relative discrepancy from the fixed term is \(X^{-c+o(1)}\) for some \(c>0\). The same holds with a gap at least \(\lambda^{1/4}\), using the low-cofactor bound just obtained. Taking \(H=X^{O_J(\lambda)}\) and \(C\) equal to the fixed tall term divided by \(D\) gives \(K\to\infty\). This counts the two high roots and the entire low cofactor at cost \(o(\lambda)\). The cofactor has only divisor-many decompositions into the fixed anatomic indices.

The order of limits is important. First fix \(\lambda\), a truncation loss, and the maximum index from BBLT Proposition 2.1. The finite set of degree tuples is fixed before \(X\to\infty\). At that stage take the largest constant in the strip lemma over the finitely many \(p,q\). The bound on precision \(K_0\), the only part needed uniformly as the degrees vary with \(\lambda\), does not involve the exponents \(p,q\). The small BBLT coefficient can be fixed at the truncation loss or included in the low cofactor bound, which takes the loss \(o(\lambda)\). [BBLT v2, Proposition 2.1](https://arxiv.org/html/2410.12234v2#S2)

A selected heavy coordinate can have zero limiting weight. It still has normalized mass at least \(\theta\), is included in \(H_{\mathrm{heavy}}\), and counts among the available omissions in the ratio bound. Thus the final \(2+2\) case has exhausted all selected heavy coordinates, including any possible zero-weight ones. A bounded-degree factor of positive normalized mass remains in the low cofactor. Unselected factors of unbounded degree are covered collectively by the same stopping argument.

## The short-height-zero bound \(d\le11/27\)

Let \(h_A=0\) be the limiting weight of the short term, with both tall weights equal to one, and suppose \(d>11/27\). Put \(u=1-2d\). We use pair counting \(l_T\le1-d\), the partial modulus lock (A), the bounded-root release (30), the scalable single-root congruence, and (PG1).
If a tall letter has mass \(s\le u\), the other tall letter has mass at least \(d-s>u\). Fix the smaller tall letter. On the other, retain up to two heavy roots and count them together with the entire low cofactor by the strip lemma. With a singleton the cost is \(s\). With two or more roots, the two retained masses sum to at least \(l_T-u\), so the frozen complementary heavy mass is at most \(u\). This gives \(d\le2u\), a contradiction.

Both tall letters now have heavy classifications. If one is effectively single, the scalable modular count uses the whole short term and the effective singleton's low complement, both of numerical height \(X^{o(1)}\). If the opposite tall letter has at least two positive-weight roots, lock either to get \(d\le l_{\mathrm{opp}}/2\le(1-d)/2\). If both are effectively single, locking either gives \(d\le\min(\beta_B,\beta_C)\), and releasing both main roots by (30) gives \(d\le1-\beta_B-\beta_C\). Both cases give \(d\le1/3\).

Each tall letter therefore has at least two positive-weight heavies. Let \(H_i\) be its total heavy mass and \(H=H_B+H_C\le1-d\), by (30). Lock a positive-weight root on one tall term and freeze all heavies on the opposite term. The modulus fixes the reduced ratio of the opposite low cofactor and the whole short term, since both have numerical height \(X^{o(1)}\). Primitivity then determines the two complete terms and hence the triple. This ratio count omits all other heavy roots on the modulus letter.

Select two positive-weight roots on each tall letter \(i\), and let \(E_i\) be the mass of all its other heavy roots. Fixing the opposite tall term and those other heavies, the real strip lemma counts the two chosen roots and the entire low cofactor at vanishing normalized cost. Thus \(l_j\ge d-E_i\) for \(j\ne i\). If both tall letters have at least three heavies, \(E_i\ge\theta_i\), so
\[
2E_B+E_C\ge d-u,\qquad 2E_C+E_B\ge d-u.
\]
The cheaper of the two selected positive roots on \(i\) costs at most \((H_i-E_i)/2\). Use it as modulus in the ratio count on each tall letter and sum:
\[
2d\le(3H-E_B-E_C)/2
\le3(1-d)/2-(d-u)/3 .
\]
It follows that \(27d\le11\), a contradiction.

If one tall letter has just two heavies and the other has more, the real count on the two-root letter gives \(l_{\mathrm{other}}\ge d\). Thus every heavy on the other letter has mass at least \((d-u)/2\). Omitting two of them in the ratio count gives \(d\le1-d-(d-u)\), hence \(d\le2/5\).

Finally suppose there are two positive-weight heavies on each tall letter. The real count gives \(l_B,l_C\ge d\), so the short mass \(s\le u\). The ratio count can omit any one heavy, giving \(H\ge4d/3\) and \(\beta\le H-d\le u\) for all four masses. Their thresholds are at least \((d-u)/2>u/2\), since \(d>2/5\). By (A), all four roots can be treated as prime. Put \(f_i=l_i-H_i\). If \(i\) is the free tall letter and \(j\) the modulus letter, locking the cheaper root on \(j\) and applying (PG1) to the two free primes on \(i\), its low cofactor and the whole short term gives
\[
d\le H_j/2+G(f_i,s)\le H_j/2+\max(s,f_i+s/2).
\]
If any \(f_i\le s/2\), this bounds \(d\) by \(u+s\le2u\). Otherwise add the two estimates:
\[
2d\le H/2+f_B+f_C+s=1-H/2\le1-2d/3.
\]
This gives \(d\le3/8\). All cases contradict \(d>11/27\).

## The full small-budget coefficient \(8/19\)

For every sufficiently small fixed \(\lambda>0\),
\[
N_\lambda(X)\ll_{\lambda,\varepsilon}
X^{(8/19)\lambda+\varepsilon}.
\]
The zero-short-height regime is proved above. We now give the positive-height proof in the limiting notation of the \(3/7\) argument. Besides pair counting, the input rows are the partial higher-moment lock (A), the bounded-root release (30), the scalable single-high congruence in “Scalable amplification and stability near the half bound,” and the signed-amplification graph (PG1) in “A prime-root graph and a larger explicit small-budget gap.” The scalable count applies to a positive-weight powerful modulus and one positive-weight root with a variable cofactor of numerical height \(X^{o(1)}\). It costs \(o(L)\) in the \(X\)-exponent. The residual variable may also be a ratio of two such low cofactors.

Suppose \(d\ge8/19\) in a limiting sequence. Write \(u=1-2d\) and \(s_0=7-16d\). We have \(d\le3/7\), \(\sum_Tl_T=1\), \(l_T\le1-d\), and \(s_0\ge u\). On a letter \(T\) of mass above \(u\), the heavy threshold is \(\theta_T=(l_T-u)/2\). Its boundedly many heavy coordinates all have normalized mass at least \(\theta_T\), and its low complement has numerical height \(X^{o(1)}\). All three terms now have positive limiting weight. A letter above \(u\) is effectively single if only one heavy root has positive limiting weight, and genuinely multiple if at least two do. For an effective singleton its other heavy roots may be included in the numerical \(X^{o(1)}\) cofactor.

We will use two forms of a prime-factor restriction. On a positive-weight heavy root, the partial lock (A) allows us to retain only prime factors of normalized mass at least \(\theta_T-c\), for any arbitrarily small fixed \(c>0\), after fixing a negligible cofactor. Indeed a divisor of degree \(e\) with unnormalized mass between \(2\lambda/e\) and \((\theta_T-c)L\) gives a normalized saving of \(c\). Sum over its dyadic scales. The product of prime factors below the lower cutoff, counted with multiplicity, is itself below that cutoff, or a first partial product reaching it would give such a divisor. In particular, a root of mass below \(2\theta_T\) can be treated as prime. We can also use a stronger threshold \(t\) if a separate count shows that a positive-weight divisor of the root costs at least \(t\). For any \(0<\zeta<(t-c)/2\), exclude a divisor of normalized mass in \([\zeta,t-c]\) by that count. All such divisors have positive limiting weight since the root has degree of order \(1/L\). The same first-crossing argument leaves a fixed remainder of mass below \(\zeta\) and prime factors of mass above \(t-c\). Take \(\zeta,c\) arbitrarily small. When selecting one prime of a root, a repeated factor of the same prime can either be frozen with the other factors or included by increasing the exponent. The anatomic coprimality makes all chosen graph variables from distinct coordinates pairwise coprime.

### All masses above \(s_0\)

If all letters are effectively single, with positive heavy masses \(\beta_T\), freezing any two and applying the scalable count to the third gives \(d\le\beta_T+\beta_U\) for each pair. Releasing the three roots in (30) gives \(\sum\beta_T\le1-d\), hence \(5d\le2\).

If exactly one letter \(T\) is genuinely multiple, select a positive heavy of mass \(\beta\) and another of mass \(f\ge\theta_T\). Lock the first as modulus and alternately freeze the main root on either effective singleton \(U,V\). Releasing these four roots in (30) yields
\[
2d\le2\beta+\beta_U+\beta_V
\le1-d+\beta-f
\le1-d+l_T-2f
\le2-3d.
\]
These two arguments only need all masses above \(u\).

If \(F,J\) are genuinely multiple and \(M\) is effectively single, let \(S_F,S_J\) be their total heavy masses, and select positive heavies \(\beta_F,\beta_J\). Freeze all heavies on \(F\), lock the selected root on \(J\), and count the effective singleton on \(M\). Interchange \(F,J\) to obtain \(2d\le S_F+S_J+\beta_F+\beta_J\). Here
\[
\beta_F+\beta_J\le S_F+S_J-(1-l_M-2u)/2,\qquad
S_F+S_J\le1-d-(l_M-u)/2.
\]
The second inequality follows by releasing all heavies on \(F,J\) and the main root on \(M\). Thus \(4d\le3/2+2u-l_M/2\) and \(l_M\le s_0\).

For three genuinely multiple letters, choose two positive roots on each. The total slack over the six thresholds in (30) is at most \(3u-d\). On each letter \(\theta_T>3u-d\), since \(l_T>s_0\). Thus there are exactly two positive-weight heavy roots per letter, both prime after the partial lock.

Whenever this profile has two prime heavies per letter, write \(m_T\) for the smaller heavy mass, \(f_T\le u\) for the low mass, \(H_{\mathrm{heavy}}\) for total heavy mass and \(F_{\mathrm{low}}=1-H_{\mathrm{heavy}}\). In three cyclic choices, use the smaller root on one letter as modulus and freeze the smaller root on the next. Apply (PG1) to the three remaining primes on the two nonmodulus letters. This gives
\[
3d\le 2\sum_Tm_T+\sum_{\{T,U\}}G(f_T,f_U).
\]
Since \(G(f,g)\le f+g-\min(f,g)/2\) and all low masses are at most \(u\), the sum of the three pair minima is at least \(2F_{\mathrm{low}}-3u\). In fact, ordering \(f_1\le f_2\le f_3\) gives \(2f_1+f_2=2F_{\mathrm{low}}-f_2-2f_3\). Also \(2\sum m_T\le H_{\mathrm{heavy}}\). Hence \(3d\le1+3u/2\), giving \(d\le5/12<8/19\).

### A mass at most \(s_0\)

Let \(S\) have smallest mass \(s\le s_0\), and let \(B,C\) be the other letters. Both are above \(u\), since otherwise two masses would be at most \(u\) and the remaining mass at most \(1-d\), giving \(d\le2u\). Each is also at least \(d-s\).

Here is a modulus construction used below. Given specified roots on \(B,C\) of total mass \(H_0\), let \(\sigma=1-d-H_0\). By (30) each coordinate on \(S\) with degree tending to infinity has mass at most \(\sigma+o(1)\), since it can be released as well. The bound is uniform across those coordinates by choosing a maximizing coordinate on a subsequence. Let \(h>0\) be the limiting numerical height of \(S\). Its part in degrees at least \(\lfloor h/(2L)\rfloor\) has positive weight, since the complementary part has weight at most \(sh/2+o(1)\). Group the high coordinates into blocks of mass at most \(\sigma+\xi+o(1)\), with any fixed \(\xi>0\), starting a new block when the next coordinate would exceed the cap. Consecutive completed blocks have masses summing to more than the cap, so there are boundedly many blocks and at least one has positive weight. This gives a modulus of limiting cost at most \(\sigma+\xi\) and fullness threshold of order \(1/L\). We may let \(\xi\) decrease to zero in the limiting bounds. When \(S\) is genuinely multiple and \(s>u\), we may instead use the cheaper of two positive heavy roots on it at cost at most \(\sigma/2\), or one positive root at cost at most \(\sigma-\theta_S\), by releasing the two chosen roots on \(S\) along with the specified \(B,C\) roots.

We will also use three direct consequences of (PG1). For \(f,g\le u\),
\[
G(f,g)\le\max\bigl(u,\tfrac34(f+g)\bigr),\qquad
G(f,g)\le(f+g+u)/2 .
\tag{G2}
\]
The first bound in fact holds with \(\max(f,g)\) in place of \(u\). For the second, when \(f\le g\) we have \(G(f,g)=\max(g,f+g/2)\). Equality in the second bound requires \((f,g)=(0,u)\) or \((u,u)\), up to exchange. For arbitrary nonnegative \(f,g\), we also have \(G(f,g)\le\max(g,f+g/2)\), since \(f\le f+g/2\) and the minimum term in (PG1) is at most \(f+g/2\).

#### Both \(B,C\) effectively single

If \(s>u\), the earlier all-single or one-genuinely-multiple calculation applies. Thus \(s\le u\). With \(S\) fixed, locking any positive-weight divisor of the main root on either \(B\) or \(C\) and counting the other effective singleton shows that the divisor costs at least \(t=d-s\). In particular \(\beta_B,\beta_C\ge t\). By (30), \(\beta_i\le1-d-t=u+s\le2u<2t\), so each main root is prime after the divisor restriction. Put \(H_0=\beta_B+\beta_C\ge2(d-s)\) and choose the small-letter modulus at cost \(\sigma=1-d-H_0\). With the modulus and either opposite term fixed in full, the scalable count gives \(l_i\ge d-\sigma=2d+H_0-1\). Thus for \(f_i=l_i-\beta_i\),
\[
f_B+f_C\le u+s,\qquad
\max(f_B,f_C)\le2-3d-H_0 .
\]
The second inequality uses \(l_i\le1-s-(2d+H_0-1)\) and \(\beta_i\ge d-s\). Apply (PG1) to the two primes without a further freeze. For \(d\ge8/19\) and \(s\le u\), we have \(\sigma\le1/19\), \(2-3d-H_0\le4/19\), and \(3(u+s)/4\le9/38\). The total normalized cost is at most \(11/38<8/19\).

#### \(B\) genuinely multiple, \(C\) effectively single

With \(S\) fixed, any positive-weight divisor of a positive heavy on \(B\) can be used as a modulus to free the effective singleton on \(C\). Its cost is thus at least \(t=d-s\). Select one prime factor from each of two positive heavies on \(B\), and denote their total mass by \(R_B\ge2t\). Let \(H_B\) be the total heavy mass on \(B\), including all factors of those two roots, and \(\beta_C\) the main heavy mass on \(C\). Retain one prime of mass at least \(\theta_C\) on \(C\). Freeze the remaining factors, at cost at most \(H_B-R_B+\beta_C-\theta_C\), and use (PG1) for the three selected primes. Its inputs are \(f_B=l_B-H_B\le u\) and \(f_C=l_C-\beta_C\). We have
\[
G(f_B,f_C)-\theta_C
\le\max\bigl(f_C-\theta_C,u+f_C/2-\theta_C\bigr)
\le\max\bigl(u,7u/4-l_C/4\bigr),
\tag{G3}
\]
using \(f_C\le l_C-\theta_C\).

If \(s\le u\), the modulus costs at most \(1-d-H_B-\beta_C\), and the heavy masses that were fixed cancel from the combined bound. Thus \(4d\le1+2s+7u/4\le1+15u/4\), giving \(d\le19/46<8/19\).

If \(s>u\) and \(S\) is effectively single, the one-genuinely-multiple calculation gives \(d\le2/5\). Otherwise \(S\) is genuinely multiple, so choose its modulus at cost at most \(1-d-H_B-\beta_C-\theta_S\). Since \(l_C\ge s\), the two terms in (G3) give
\[
4d\le1+2s-\theta_S+
\max\bigl(u,7u/4-s/4\bigr).
\]
For the first term this gives \(7d\le5/2+3s/2\le13-24d\), hence \(d\le13/31\). The second gives
\[
17d/2\le13/4+5s/4\le12-20d,
\]
hence \(d\le8/19\).

Equality in the latter calculation forces \(s=l_C=s_0=5/19\), \(\theta_C=\theta_S=1/19\), and \(f_B=u\). Equality in \(f_C\le l_C-\theta_C\) gives \(\beta_C=1/19\). Thus \(l_B=9/19\) and \(H_B=l_B-u=6/19\). The cost \(1-d-H_B-\beta_C-\theta_S\) used for the modulus is \(3/19\). But we can use the cheaper of the two positive-weight heavy roots on \(S\), costing at most \((1-d-H_B-\beta_C)/2=2/19\). This makes the count strict at \(d=8/19\).

#### Both \(B,C\) genuinely multiple

Let \(H=H_B+H_C\) be their total heavy mass. Choose two positive-weight heavy roots on each letter and retain one prime from each of the four, of total mass \(R\). Then
\[
R\ge2\theta_B+2\theta_C=1-s-2u .
\]
Freeze the remaining heavy factors, of mass \(H-R\), and the cheapest retained prime, of mass at most \(R/4\). The other three are counted by (PG1) with low masses \(f_i=l_i-H_i\le u\).

If \(s\le u\), the modulus from \(S\) costs at most \(1-d-H\). Thus
\[
d\le1-d-3R/4+3u/2
\le(1+3s)/4+3u-d,
\]
giving \(d\le8/19\). Equality forces \(s=u=3/19\), \(R=10/19\), and all four retained primes to have mass \(5/38\). It also forces \(G(f_B,f_C)=3u/2\), hence \(f_B=f_C=u\). The two thresholds coincide with the retained prime masses, so \(l_B=l_C=8/19\) and \(H=R=10/19\). The modulus bound used was \(\sigma=1/19\).

In this equality case, either \(S\) supplies a cheaper positive-weight modulus or it is effectively a single root with numerical \(X^{o(1)}\) cofactor. To verify the dichotomy, apply the high-block construction with a strict cap below \(\sigma\) if the largest high coordinate mass stays below \(\sigma\) by a fixed margin. Otherwise take a subsequence with a selected high coordinate attaining mass \(\sigma\) in the limit. By (30), releasing that coordinate and any other unbounded-degree coordinate on \(S\) alongside the \(B,C\) roots implies that each other coordinate of unbounded degree has mass \(o(1)\), uniformly. If these other factors and the remaining lower-degree factors together carry positive limiting height, apply the high-block construction to this remainder with any arbitrarily small fixed cap. This supplies a cheaper modulus. In the alternative, the \(\sigma\)-mass coordinate carries all positive height, and the remainder has height \(X^{o(1)}\). Lock a root on \(B\), fix both roots on \(C\), and use the scalable count on \(S\). The cost \(5/38+10/38=15/38\) is below \(8/19\).

Now take \(s>u\). If \(S\) is effectively single, use the smaller selected prime on \(B\) as modulus and freeze the heavy part of \(C\), then interchange the roles. The scalable count freeing \(S\) gives \(2d\le R/2+H\le3H/2\), so \(H\ge4d/3\). With the small-letter modulus costing at most \(1-d-H\), the prime graph estimate gives
\[
d\le 1-d-3R/4+
\max\bigl(u,\tfrac34(1-s-4d/3)\bigr).
\]
The first term gives \(7d\le11/4+3s/4\le8-12d\), hence \(d\le8/19\). In equality \(s=s_0=5/19\). The main root on \(S\) then has mass at least \(\theta_S=1/19\), so (30) gives \(H\le1-d-\theta_S=10/19\), contradicting \(H\ge4d/3=32/57\). The second term gives \(3d\le1+3u/2\), hence \(d\le5/12\).

If \(S\) is genuinely multiple, choose its cheaper positive modulus, of cost at most \((1-d-H)/2\). Use the second bound in (G2) to obtain
\[
\begin{aligned}
d&\le (1-d-H)/2+H-3R/4+(1-s-H+u)/2\\
 &\le 1/4-d/2+s/4+2u .
\end{aligned}
\]
This gives \(22d\le9+s\le16-16d\), hence \(d\le8/19\). Equality forces \(s=s_0=5/19\) and \(R=8/19\), with the four retained masses all \(2/19\). Therefore \(l_B=l_C=7/19\). Equality in the second bound of (G2) requires \((f_B,f_C)=(u,u)\) or \((0,u)\) up to exchange. Releasing the two positive heavies on \(S\) gives \(H\le1-d-2\theta_S=9/19\), so \(f_B+f_C\ge5/19>u\), excluding the latter. Thus \(H=R=8/19\). Equality in the cheaper-modulus bound requires that both chosen positive roots on \(S\) have mass \(3/38\) and exhaust its heavy mass, by (30). Their threshold is \(1/19\), so both can be treated as prime. Consequently all three letters have exactly two positive prime heavies and low masses at most \(u\). The cyclic graph count gives \(d\le5/12\), as already proved.

These cases exclude a limiting \(d\ge8/19\). All endpoint contradictions have a strict margin after taking the permissible divisor-threshold and modulus-block errors arbitrarily small. If the displayed fixed small-budget bound failed for arbitrarily small budgets, the fixed-truncation and dyadic passage of the \(3/7\) proof would give such a limiting sequence, since \(L\le\lambda+o(\lambda)\). This proves the bound.

## An explicit small-budget coefficient \(13/31\)

For every sufficiently small fixed \(\lambda>0\),
\[
N_\lambda(X)\ll_{\lambda,\varepsilon}
X^{(13/31)\lambda+\varepsilon}.
\]
The zero-short regime is already bounded by \(11/27\). In the positive-height regime suppose a limiting sequence has \(d\ge13/31\). Keep \(u=1-2d\), the heavy thresholds \(\theta_T=(l_T-u)/2\), and the same meaning of effectively single and genuinely multiple. We use (30), (A), and the prime-root graph (PG1) with signed amplification, as in the \(8/19\) proof.

The prime-factor restriction from (A) applies also to a heavy root of zero limiting weight. Its degree \(e\) still tends to infinity, so the lower divisor cutoff \(2\lambda/e\) is \(o(L)\). After discarding the families bounded by (A), the small product can be fixed at negligible cost and each remaining prime factor, counted with multiplicity, has normalized mass at least \(\theta_T\), with arbitrarily small slack. The number of heavy roots is bounded whenever \(l_T>u\). Below, when a prime is retained from a positive-weight root it has positive limiting weight, while an extra factor of zero weight contributes to the heavy mass accounting. If a prime occurs twice in a root, it contributes its mass twice in that accounting, and its contribution to the term can be represented by doubling its exponent.

### A modulus of half the available mass

Suppose \(S\) has positive limiting height and a chosen collection of heavy roots on the other two letters has combined mass \(H\). Set \(\sigma=1-d-H\ge0\). By (30), any fixed bounded selection of unbounded-degree coordinates on \(S\) has combined mass at most \(\sigma\). Then either \(S\) supplies a positive-weight modulus of cost at most \(\sigma/2\) with arbitrarily small slack, or \(S\) consists of one positive-weight root and a cofactor of numerical height \(X^{o(1)}\).

To prove this, restrict first to the degrees at least \(\lfloor h/(2L)\rfloor\), with \(h>0\) the limiting height of \(S\). Those degrees carry positive height. If their maximum coordinate mass is at most \(\sigma/2+o(1)\), the block construction of the \(8/19\) proof gives the modulus of cost at most \(\sigma/2+\xi\) for arbitrary fixed \(\xi>0\). Otherwise take a subsequence on which one high coordinate has limiting mass \(\rho>\sigma/2\). Every other unbounded-degree coordinate has mass at most \(\sigma-\rho+o(1)\) by releasing it alongside the chosen coordinate. If the remaining factors have positive limiting height, apply the same block construction to them, again giving cost at most \(\sigma/2+\xi\). If not, the chosen coordinate carries the positive height and its complement has numerical height \(X^{o(1)}\).

### A letter of mass at most \(u\)

Let \(S\) be a letter of smallest mass \(s\le u\). The other two letters \(B,C\) have masses above \(u\), since two masses at most \(u\) would give \(d\le2u\) by pair counting. If both are effectively single, the earlier count in the \(8/19\) proof applies. Fixing \(S\) and locking any positive-weight divisor shows that each main root is prime after negligible factors, since its factors cost at least \(d-s\) and its total mass is at most \(u+s<2(d-s)\). With \(H=\beta_B+\beta_C\ge2(d-s)\), the \(S\)-modulus costs at most \(1-d-H\). For the two low masses \(f_B,f_C\) that proof gives \(f_B+f_C\le u+s\) and \(\max(f_B,f_C)\le2-3d-H\). At \(d\ge13/31\) the modulus costs at most \(2/31\), and (PG1) costs at most \(\max(7/31,15/62)=15/62\). Thus \(d\le19/62\).

If one of \(B,C\) is genuinely multiple, the other effectively single, the mixed \(s\le u\) estimate (G3) in the \(8/19\) proof still gives \(d\le19/46\). Its proof uses no lower bound on \(d\) beyond \(d>2u\). Namely, retain two primes of mass at least \(d-s\) from the genuinely multiple letter and one prime of mass at least the threshold from the singleton. The modulus cost and all fixed heavy factors cancel as there, giving \(4d\le1+2s+7u/4\le1+15u/4\).

Now suppose both \(B,C\) are genuinely multiple, and let \(H=H_B+H_C\) be their total heavy mass. Retain one prime from each of two positive-weight roots on each letter. Their combined mass \(R\) is at least \(1-s-2u\). Freezing their cheapest prime and the remaining heavy factors costs at most \(H-3R/4\). The two low masses \(f_B,f_C\) have sum \(1-s-H\) and are individually at most \(u\). If \(S\) supplies a modulus of cost at most \((1-d-H)/2\), the bound \(G(f_B,f_C)\le(f_B+f_C+u)/2\) from (G2) gives
\[
d\le(1-d-H)/2+H-3R/4+(1-s-H+u)/2
\le1/4-d/2+s/4+2u .
\]
Thus \(22d\le9+s\le9+u=10-2d\), giving \(d\le5/12\).

In the alternate case of a numerically effective singleton \(S\), use the cheaper of the two selected primes on \(B\) as modulus and freeze all heavies on \(C\), then exchange \(B,C\). The scalable count on \(S\) gives
\[
2d\le R/2+H\le3H/2.
\]
In particular \(H\ge4d/3\). Use the \(S\)-modulus of cost at most \(1-d-H\) and the bound \(G(f_B,f_C)\le\max\bigl(u,3(1-s-4d/3)/4\bigr)\). The resulting estimate is
\[
d\le1-d-3R/4+\max\bigl(u,3(1-s-4d/3)/4\bigr).
\tag{A13}
\]
From the first term we get \(7d\le11/4+3s/4\le11/4+3u/4\), hence \(d\le7/17\). The second gives \(3d\le1+3u/2\), hence \(d\le5/12\).

### All masses above \(u\)

The all-effective-single and one-genuinely-multiple counts in the \(8/19\) proof give \(d\le2/5\). Suppose exactly two letters are genuinely multiple, and designate their effective singleton as \(S\), with mass \(s\). The construction of (A13) still applies and \(H\ge4d/3\). The main root on \(S\) has mass at least \((s-u)/2\). Releasing it along with \(H\) gives
\[
(s-u)/2\le1-d-H\le1-7d/3,\qquad s\le3-20d/3 .
\]
The first term in (A13) therefore gives
\[
7d\le11/4+3s/4\le5-5d,
\]
hence \(d\le5/12\). The second again gives \(3d\le1+3u/2\).

Finally suppose all three letters are genuinely multiple, and let \(S\) have smallest mass \(s\). The cheaper positive-weight root on \(S\) costs at most \((1-d-H)/2\) by (30), where \(H\) is the total heavy mass on the other two letters. The same graph estimate just used for this modulus gives \(22d\le9+s\). Thus
\[
\theta_{\min}=(s-u)/2\ge12d-5.
\]
After applying (A) to each heavy coordinate, there are at least two prime factors on each letter. The slack over their six thresholds is at most \(3u-d\), by (30). An additional prime factor costs at least its letter's threshold. Since
\[
2\theta_{\min}\ge24d-10\ge3u-d
\qquad(d\ge13/31),
\]
there can be at most one additional prime unless \(d=13/31\), where there can be two.

If there are no additional primes, the cyclic graph bound for two positive prime heavies per letter in the \(8/19\) proof gives \(d\le5/12\). If there is one additional prime on \(T\), choose the cheaper of two positive-weight primes from distinct heavy roots of \(T\) as modulus. It costs at most \((H_T-\theta_T)/2\), since the third prime has mass at least \(\theta_T\). Freeze the cheapest of the four primes on the other two letters, costing at most \((H-H_T)/4\) with \(H\) the total heavy mass on all three. With the graph cost at most \(3u/2\), and
\[
H\le1-d,\qquad H_T\le l_T+2u-d,
\]
we obtain \(d\le(1-d+H_T)/4-\theta_T/2+3u/2\le5u/2\), hence \(d\le5/12\).

Two additional primes at \(d=13/31\) force \(u=5/31\), \(s=7/31\), and \(\theta_{\min}=1/31\). They exhaust all slack \(3u-d=2/31\). Each additional factor lies on a letter of mass \(7/31\), has mass \(1/31\), and all six mandatory prime factors attain their letter thresholds. If the two additional factors lie on the same letter \(T\), lock there any positive-weight prime at cost \(1/31\). The other letters have combined heavy mass \(14/31\), so freeze one of their four primes at cost at most \(7/62\). The graph costs at most \(15/62\), giving total at most \(12/31<13/31\). If the additional factors lie on distinct letters, those two letters each have three factors of mass \(1/31\). Lock a positive prime on one of them and freeze two factors on the other, retaining one of its positive primes. The third letter has precisely two positive prime factors. If a frozen factor coincides with the retained factor, the freezing determines that prime and leaves two free primes. Apply (PG1) at total cost at most \(1/31+2/31+15/62=21/62<13/31\).

Thus no limiting \(d\ge13/31\) survives. The inequalities concerning prime-factor masses and block moduli can use arbitrarily small fixed threshold margins before taking limits, and the equality cases above have strict savings. The fixed-truncation and dyadic passage gives the claimed fixed coefficient.

## A limiting small-budget coefficient of \(5/12\)

For every \(\eta>0\) and all sufficiently small fixed \(\lambda>0\),
\[
N_\lambda(X)\ll_{\lambda,\eta,\varepsilon}
X^{(5/12+\eta)\lambda+\varepsilon}.
\]
In particular \(417/1000\) is an admissible fixed coefficient for all sufficiently small fixed \(\lambda\).

Use the notation and the positive-height setup of the \(13/31\) proof. The zero-short-height estimate is \(d\le11/27\). Suppose \(d>5/12\) and the short height is positive. The cases with an effective singleton or a mass at most \(u\) have bounds of at most \(5/12\).

If \(s\le u\) and both other letters are effectively single, the normalized cost from that proof is at most \(1-d-H+\max(2-3d-H,3(u+s)/4)\), with \(H\ge2(d-s)\). Here \(u\le1/6\), \(H\ge1/2\), the modulus costs at most \(1/12\), and both graph terms are at most \(1/4\). Thus \(d\le1/3\). With one of the other letters genuinely multiple and one effective singleton, the same proof gives \(d\le19/46\). If the other two are genuinely multiple, a modulus of half the available mass gives \(22d\le9+s\le9+u\), hence \(d\le5/12\). In the effective-singleton alternative, the two terms in (A13) give \(d\le7/17\) and \(d\le5/12\), respectively.

If all masses exceed \(u\), the all-single and one-genuinely-multiple arguments give \(d\le2/5\). With exactly two genuinely multiple letters, designate the effective singleton as \(S\) and use (A13). The two modular counts give \(H\ge4d/3\) for the combined heavy mass on the other two letters. Its main root has mass at least \((s-u)/2\), yielding \(s\le3-20d/3\) by release (30). The first term in (A13) gives \(7d\le11/4+3s/4\le5-5d\), and the second gives \(3d\le1+3u/2\). Both bound \(d\) by \(5/12\).

It remains that all three letters are genuinely multiple. Apply the prime-factor restriction under (A) to every heavy root, including those of zero limiting weight, as in the \(13/31\) proof. Count prime factors with multiplicity and call a letter affected if it has more than two such factors after the negligible remainders are fixed. At least two prime factors on every letter have positive weight and are from distinct coordinates. Every additional factor costs at least its letter's threshold, up to arbitrarily small slack. By (30) the available slack over the six mandatory thresholds is at most \(3u-d\). If all three letters are affected, this forces
\[
3u-d\ge\sum_T\theta_T=(1-3u)/2,
\]
so \(d\le2/5\).

If none is affected, the cyclic graph count for two prime heavies per letter in the \(8/19\) proof gives \(d\le5/12\). If exactly one letter \(T\) is affected, choose the cheaper of two positive-weight prime factors from distinct roots on \(T\) as modulus. It costs at most \((H_T-\theta_T)/2\), where \(H_T\) is its total heavy mass, since at least one other factor costs at least \(\theta_T\). Freeze the cheapest of the four primes on the other two letters, at cost at most \((H-H_T)/4\) with \(H\le1-d\) the total heavy mass. The remaining three primes and the two low cofactors are counted by (PG1) at cost at most \(3u/2\). As \(H_T\le l_T+2u-d\), this gives
\[
d\le(1-d+H_T)/4-\theta_T/2+3u/2
\le5u/2,
\]
hence \(d\le5/12\). This estimate permits any number of additional factors on \(T\).

If exactly two letters \(T_1,T_2\) are affected, the third letter \(K\) has two positive prime factors. Lock the cheaper of two positive prime factors on \(T_1\) as modulus, costing at most \((H_1-\theta_1)/2\). Retain one positive prime of mass at least \(\theta_2\) on \(T_2\) and freeze the rest of its heavy factors. A repeated copy of the retained prime can be included by multiplying its exponent. Apply (PG1) with that prime and both primes on \(K\), at cost at most \(3u/2\). Interchange \(T_1,T_2\) and sum. Writing \(s_{12}=l_1+l_2\), the result is
\[
\begin{aligned}
2d
&\le\tfrac32(H_1+H_2-\theta_1-\theta_2)+3u\\
&\le\tfrac34s_{12}+6u-\tfrac32d,
\end{aligned}
\]
where \(H_1+H_2\le s_{12}+u-d\) follows from \(H_K\ge l_K-u\). Both affected letters use at least one additional factor, so \(\theta_1+\theta_2\le3u-d\), or \(s_{12}\le8u-2d\). Substitution gives \(5d\le12u\), hence \(d\le12/29<5/12\).

We have excluded every limiting \(d>5/12\). The fixed-truncation and dyadic passage therefore gives the estimate with any positive \(\eta\).

## The unrestricted one-third endpoint

A strict uniform gap below \(1/3\) at all sufficiently small fixed
budgets is equivalent, given the limiting one-third bound, to
excluding every sequence of actual populated anatomic boxes with
\(D/L\to1/3\) and \(L\) comparable to the budget. In the converse
direction, if the gap fails take budgets \(\lambda_j\to0\) and
prospective gaps tending to zero. At each fixed \(j\) take an
unbounded height sequence witnessing failure, with arbitrarily
small fixed-truncation losses before \(j\) changes. The selected
boxes have \(D\ge(1/3-o(1))\lambda_j\) and
\(L\le(1+o(1))\lambda_j\). The limiting count forces
\(L/\lambda_j\to1\) and \(D/L\to1/3\). Their letter radical masses
are balanced, but their nonatomic total \(\Gamma\) can be anywhere
in \([1/2,1]\). A formal degree measure passing necessary tests is
not itself evidence of a populated box.

The complete-label reduction can be formulated in terms of the
actual total outside-label entropy \(h\), in units of \(L\log X\).
A family of exponent \(d\) admits a complete additive realization
after a bounded refinement of its labels and loss of edges, with
exponent at least \(10d-6h\) when this is positive. The refinement
uses anchored two-edge comparisons, including projections with
no four-cycles. Thus the range of this reduction is \(h<5/9\)
at equality and is independent of the formal upper cap \(\Gamma\).
A matching with large outside-label entropies shows that the
incidence step cannot by itself force a positive exponent there.
[Anchor reduction and matching obstruction](ABC_DIFFUSE_FACTOR_GRAPH.md#anchor-sets-for-arbitrary-nonatomic-masses)

There are three further constraints on the entire balanced endpoint.
First, every complete marginal has entropy at least
\((1/800)\log(X^L)\). A fixed dictionary of at most
\(X^{\lambda/800}\) complete values on any letter gives an
absolute saving. The proof pairs the upper-weight part of the
two free radical spectra with the actual height of the shorter
free letter. If that pair fails to reconstruct, its degree
distributions must be close enough to a shared scale for one
residue determinant to finish the count.
[Quantitative complete-value dictionaries](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-quantitative-bound-for-complete-value-dictionaries)

Second, the upper half of the pooled root mass on every pair of
letters carries uniformly more than half the combined height.
A stronger lock selection starts with one whole diffuse letter
and deletes or exchanges a small block in a bounded degree band.
It gives a **whole-coordinate pair record of cost at most
\(1/3-\eta\) and weight at least the shorter height plus \(\eta\)**
for some absolute \(\eta>0\) on at least one pair in every
actual endpoint sequence. This includes indivisible
macroscopic primes. When a particular pair fails even with
prime-resolved selection, its exceptional form is a full-mass
single prime on a letter at least as tall as its opposite one.
[Pooled degree inequality](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-pooled-degree-inequality-and-primitive-shift-fibres)
· [Whole-coordinate pair locks](ABC_WIDE_COFACTOR_CONGRUENCES.md#affordable-whole-coordinate-pair-locks-at-every-balanced-endpoint)

Third, there is an absolute exponent saving on every curve
in two complete letters, uniform in the coefficients and the
fixed degree cap in the exponent. The implicit constant depends
on the cap, which may grow arbitrarily as the outer budget
shrinks. The new input is exact divisibility of a prime valuation
by the local ramification order. In a diffuse compact band,
the few large-ramification branches account for vanishing
mass. This makes the weighted determinant work even when a
whole free letter's root mass sits in one indivisible coordinate.
An explicit support and determinant analysis also permits curve
degrees growing as a small positive power of the logarithm of
height. Both results give a strict saving after summing over an
economical curve cover.
[Absolute complete-letter curve saving](ABC_WIDE_COFACTOR_CONGRUENCES.md#an-absolute-curve-saving-with-arbitrary-fixed-budget-degrees)

The powerful-triples argument concentrates an endpoint at the
scaled order three and so contradicts the required dispersion.
The pooled inequality gives a positive total moment below
that order on every pair of letters. The height-adjusted
fullness and one-full diffuse estimates give additional
restrictions on how that low-degree mass can be arranged.
[Powerful triple endpoint transfer](POWERFUL_TRIPLES_ABC_ENDPOINT_TRANSFER.md)
· [One-full diffuse profile consequences](ABC_ONE_FULL_DIFFUSE_PROFILE_CONSEQUENCES.md)

For determinant approaches, the first rank drop alone recovers all
free radical factors outside finitely many degree positions. The
resulting count charges the atom masses at the exceptional positions
and works at arbitrary \(\Gamma\). Two separate first determinants
can instead use the exact additive relation and a function-field
height argument to force large degree on their narrow projection,
even when all three numerical heights are one. The broad formal
profile \(\nu_T=\frac13\operatorname{Unif}[2,4]\) on all three letters
exposes a limit of the direct Newton approach. One-root supports cannot meet both Taylor and chart budgets, and
direct specialization of several groups has an independent
chart-budget obstruction at total diffusion one. For rectangular
first determinants on this profile a successful chart budget would
have to use a group whose entire root dictionary costs less than
its shallow residue. Freezing such groups and rerunning the
rectangular method leaves the obstruction in place.
[Direct Newton counting and its broad-profile barrier](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#direct-newton-recovery-at-arbitrary-nonatomic-mass)
· [Rectangular multigroup obstruction](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#a-broad-diffuse-obstruction-for-rectangular-multigroup-determinants)
· [Two first determinants with the additive relation](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#two-first-residue-determinants-and-the-exact-additive-relation)

The same divisibility gives a terminal count on an auxiliary curve
when its exact odds map has degree \(O(1/L)\). An extension to
higher-dimensional full components uses a normal-crossings
model of the odds coordinates. In a diffuse band the valuation
orders force at least one arbitrarily deep local direction,
and a determinant on the narrow projection cuts dimension
after a cheap boundary record. This applies to components
of bounded full degree with a quasi-finite narrow projection.
When fixed-support removal leaves purely atomic free spectra,
the already fixed full third modulus supplies the wide-cofactor
fallback. A one-step version also handles growing full degree if the
scaled boundary orders of a smooth model have a compact
countable limiting set. A tangency example shows that narrow
degree and odds degree alone can permit a dense set of
exceptional orders on a resolution.
[Auxiliary-curve terminal](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#a-diffuse-terminal-count-on-auxiliary-curves)
· [Higher-dimensional descent](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#dimension-descent-from-diffuse-locks-on-bounded-degree-full-components)

A first determinant can also retain Taylor rank past the
least individual lifting order by capping each term at the
available prime-power precision. This yields a strict saving
on an open wholly diffuse sector missed by the former
full-grouping rectangle. Expanding the exact error as a
second local coordinate gives another strict sector beyond
that cap. On the broad uniform \([2,4]\)
example, it still fails for arbitrary sparse one-root supports.
For any fixed number of rectangular groups the obstruction
also persists with separate depths per coordinate and prime,
unless some whole group dictionary costs less than its own
shallow chart. The exact error expansion also has a broad
one-root obstruction by mixed area for sparse supports
with arbitrary nonzero signed multipliers and adaptive depths
under the primitive maximum row height bound. An exact
column-mean estimate leaves a distinct asymmetric
Archimedean optimization.
[Full-modulus Taylor precision](ABC_TAYLOR_PRECISION_SATURATION.md)

There is also an exclusion for rational projective parametrizations
of degree at least three with small coefficient height and positive
nonatomic mass on all three letters. This condition allows the
degree to grow in the outer limit. The proof compares the
ramification-weighted locks with the actual height of the
projective parameter, including the improvement when a short
letter forces that parameter near a real contact. It also excludes
a split quadratic with two ramified boundary fibres by applying
the unrestricted count to three derived linear forms.
[Rational projective parametrizations](ABC_RATIONAL_PROJECTIVE_CARRIERS.md)

### Residual obstruction

A surviving family must have large complete dictionaries on
every letter, evade economical complete-letter curve covers,
and avoid the positive-exponent anchor reduction. Where the
pooled locks have cost \(t<1/3\) and combined numerical weight
above the shorter free height, (PQ4) isolates the missing bound.
A sufficient statement is that each fixed record admits at
most \(n^{1/3-t-\sigma}\) distinct targets in the third index
box for a uniform \(\sigma>0\). Both free values then vary
through prescribed multiples of the two full divisors, with
the remaining radical caps inherited from their boxes.
Such affordable records now exist on at least one pair on
every endpoint, including the large-prime atomic sector.
On a matching most such shifts have one representation, so
even the improved uniform bounds for each shift allow all of them. There is an explicit primitive cubic
matching with saturated internal root entropy, demonstrating
that primitivity, additivity and within-letter independence
alone allow this behavior. For that family \(\operatorname{rad}(abc)\)
has height exponent two, and the family also requires many
bounded-degree complete-letter curves to cover it. A separate
primitive quartic family achieves the unsaved fixed-record
exponent with all three radical caps at a budget tending to
\(3/2\). The unresolved counting step has to
use the much smaller radical budget or the consequent high
valuation degrees in a way that couples different edges.
[Primitive matching and curve cover](ABC_DIFFUSE_FACTOR_GRAPH.md#complete-letter-curve-covers-of-the-cubic-matching)
· [Balanced large-budget two-divisor obstruction](ABC_DIFFUSE_FACTOR_GRAPH.md#primitive-two-divisor-fibres-attaining-the-unsaved-count)

### Scope of the few-prime sum-product approach

At fixed \(L\), set \(n=X^L\). For a balanced box, the number of
complete values on a letter with at least
\(k=\lceil\eta\log n/\log\log n\rceil\) distinct primes is
\(n^{1/3-\eta+o(1)}\) up to arbitrary box tolerances, for fixed
\(\eta>0\). Indeed the letter radical is at most \(n^{1/3+o(1)}\).
Multiplying this bound by
\((\sum_{p\le n}1/p)^k/k!\) bounds the number of possible squarefree
radicals with at least \(k\) prime factors. The product is
\(n^{-\eta+o(1)}\) by Stirling's formula and the prime harmonic
upper bound. At each fixed budget the fixed-support count gives
an arbitrary subpower allowance for the choices of the integer
with that radical. If a letter marginal has entropy approaching
\((1/3)\log n\), the high-prime-count values consequently have
vanishing edge probability for each fixed \(\eta\), by the entropy
bound for a distribution split into two sets. Taking \(\eta\)
to zero after the inner height limits permits a diagonal restriction to the few-prime range
in the high marginal entropy branch. Bloom's Theorem 7 supplies a lower bound for the maximum of
the sumset and product-set sizes in that range. A sparse matching
of sums is compatible with both sets being large, so a
size-disjunction of this kind does not settle the sparse branch.
The public Pell saving has a different scope, with a prescribed
one-parameter family instead of arbitrary complete dictionaries.
[Thomas Bloom, Theorem 7](https://arxiv.org/html/2512.04931v1)
· [Chunping Zeng, balanced radical triples](https://hexagonmath.org/pdf/2611.00611v1)


The bounded-degree curve count has a different range from ordinary
coefficient-uniform integral-point bounds. For example, the
general curve bound in Binyamini, Cluckers and Kato has the
factor \(X^{1/D}\) in the ambient height, whereas the argument here exploits two radical spectra at exponent
\(O(\lambda)\). Zeng's complete-modulus exponent certificate
has an unrestricted bound uniform in \(0<\lambda\le1\) and
uses a finite set of fixed degrees. Its stated exponent is independent of the small budget.
These comparisons identify the scope of the argument.
[Binyamini, Cluckers and Kato, Corollary 8](https://arxiv.org/html/2401.03982v2)
· [Zeng, Complete modulus locking](https://hexagonmath.org/pdf/2611.00621v1)
